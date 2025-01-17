from django.contrib import admin
from .models import UserAnalysis, WebSocketBackendUUID, ToolVersion, Report


@admin.register(ToolVersion)
class ToolVersionAdmin(admin.ModelAdmin):
    list_display = (
        "nextclade_version",
        "pangolin_version",
        "constellations_version",
        "scorpio_version",
        "usher_version",
        "gofasta_version",
        "minimap2_version",
        "faToVcf_version",
        "created_at",
    )
    search_fields = ("nextclade_version", "pangolin_version", "created_at")
    list_filter = ("created_at",)


@admin.register(UserAnalysis)
class UserAnalysisAdmin(admin.ModelAdmin):
    list_display = ("analysis_id", "user_id", "submission_date", "celery_task_id")
    list_filter = ("submission_date", "overall_status")
    search_fields = ("analysis_id", "user_id")
    readonly_fields = ("submission_date",)

    # Optional: Customizing the admin form to handle JSONField display
    def formfield_for_dbfield(self, db_field, **kwargs):
        if db_field.name == "analysis_status":
            from django.forms import Textarea

            kwargs["widget"] = Textarea(attrs={"rows": 10, "cols": 80})
        return super().formfield_for_dbfield(db_field, **kwargs)


@admin.register(WebSocketBackendUUID)
class WebSocketBackendUUIDAdmin(admin.ModelAdmin):
    list_display = ("uuid", "created_at")
    list_filter = ("created_at",)
    search_fields = ("uuid",)
    readonly_fields = ("created_at",)


@admin.register(Report)
class ReportAdmin(admin.ModelAdmin):
    list_display = (
        "name",
        "report_type",
        "graph_type",
        "user_analysis",
        "created_at",
    )  # Columns to display in the list view
    list_filter = ("report_type", "graph_type", "created_at")  # Filters on the sidebar
    search_fields = ("name", "user_analysis__analysis_id")  # Searchable fields
    autocomplete_fields = ("user_analysis",)  # Autocomplete for related fields
    filter_horizontal = (
        "summary_sources",
    )  # UI for selecting ManyToMany relationships
    readonly_fields = ("created_at",)  # Fields that cannot be edited

    # Customizing the admin form view
    fieldsets = (
        (None, {"fields": ("name", "report_type", "graph_type", "user_analysis")}),
        ("Content", {"fields": ("data", "text_summary")}),
        ("Relationships", {"fields": ("summary_sources",)}),
        ("Metadata", {"fields": ("created_at",)}),
    )
