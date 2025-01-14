# admin.py
from django.urls import reverse
from django.contrib import admin
from django.utils.html import format_html
from .models import (
    VirusStrain, Location, AAMutation, NucMutation,
    StructuralVariation, QualityMetrics
)

class LocationInline(admin.StackedInline):
    model = Location
    can_delete = False
    extra = 0

class AAMutationInline(admin.TabularInline):
    model = AAMutation
    extra = 0
    verbose_name = "Amino Acid Mutation"
    verbose_name_plural = "Amino Acid Mutations"

class NucMutationInline(admin.TabularInline):
    model = NucMutation
    extra = 0
    verbose_name = "Nucleotide Mutation"
    verbose_name_plural = "Nucleotide Mutations"

class StructuralVariationInline(admin.StackedInline):
    model = StructuralVariation
    can_delete = False
    extra = 0

class QualityMetricsInline(admin.StackedInline):
    model = QualityMetrics
    can_delete = False
    extra = 0

@admin.register(VirusStrain)
class VirusStrainAdmin(admin.ModelAdmin):
    list_display = [
        'strain', 'gisaid_epi_isl', 'collection_date', 
        'nextclade_clade', 'nextclade_lineage', 'pangolin_lineage',
        'coverage', 'pangolin_qc_status', 'nextclade_qc_status', 'location_info',
        'mutation_count'
    ]
    list_filter = [
        'nextclade_clade', 'nextclade_lineage', 'pangolin_lineage',
        'pangolin_qc_status', 'nextclade_qc_status', 'location__state', 'location__country'
    ]
    search_fields = [
        'strain', 'gisaid_epi_isl', 'genbank_accession',
        'location__country', 'location__state'
    ]
    date_hierarchy = 'collection_date'
    readonly_fields = ['mutation_count']
    
    inlines = [
        LocationInline,
        AAMutationInline,
        NucMutationInline,
        StructuralVariationInline,
        QualityMetricsInline
    ]

    fieldsets = (
        ('Basic Information', {
            'fields': (
                'strain', 'gisaid_epi_isl', 'genbank_accession',
                'collection_date'
            )
        }),
        ('Host Information', {
            'fields': ('host', 'age', 'sex')
        }),
        ('Classification', {
            'fields': ('nextclade_clade', 'nextclade_lineage', 'pangolin_lineage')
        }),
        ('Quality', {
            'fields': ('coverage', 'pangolin_qc_status', 'nextclade_qc_status', 'mutation_count')
        })
    )

    def location_info(self, obj):
        if hasattr(obj, 'location'):
            return f"{obj.location.country} - {obj.location.state}"
        return "No location info"
    location_info.short_description = "Location"

    def mutation_count(self, obj):
        aa_count = obj.aa_mutations.count()
        nuc_count = obj.nuc_mutations.count()
        return f"AA: {aa_count}, Nuc: {nuc_count}"
    mutation_count.short_description = "Mutation Count"

    def get_queryset(self, request):
        return super().get_queryset(request).select_related(
            'location',
            'structural_variation',
            'quality_metrics'
        ).prefetch_related(
            'aa_mutations',
            'nuc_mutations'
        )

@admin.register(AAMutation)
class AAMutationAdmin(admin.ModelAdmin):
    list_display = ['strain_link', 'protein', 'reference_aa', 'position', 'altered_aa']
    list_filter = ['protein']
    search_fields = [
        'strain__gisaid_epi_isl',
        'protein',
        'reference_aa',
        'position',
        'altered_aa'
    ]

    def strain_link(self, obj):
        url = reverse('admin:query_engine_virusstrain_change', args=[obj.strain.id])
        return format_html('<a href="{}">{}</a>', url, obj.strain)
    strain_link.short_description = 'Strain'

@admin.register(NucMutation)
class NucMutationAdmin(admin.ModelAdmin):
    list_display = ['strain_link', 'position', 'reference_base', 'altered_base']
    list_filter = ['reference_base', 'altered_base']
    search_fields = [
        'strain__gisaid_epi_isl',
        'position',
        'reference_base',
        'altered_base'
    ]

    def strain_link(self, obj):
        url = reverse('admin:query_engine_virusstrain_change', args=[obj.strain.id])
        return format_html('<a href="{}">{}</a>', url, obj.strain)
    strain_link.short_description = 'Strain'

@admin.register(Location)
class LocationAdmin(admin.ModelAdmin):
    list_display = [
        'strain_link', 'continent', 'country', 'state',
        'continent_exposure', 'country_exposure', 'state_exposure'
    ]
    list_filter = ['continent', 'country', 'state', 'continent_exposure', 'country_exposure', 'state_exposure']
    search_fields = [
        'strain__gisaid_epi_isl',
        'continent', 
        'country', 
        'state',
    ]

    def strain_link(self, obj):
        url = reverse('admin:query_engine_virusstrain_change', args=[obj.strain.id])
        return format_html('<a href="{}">{}</a>', url, obj.strain)
    strain_link.short_description = 'Strain'

@admin.register(QualityMetrics)
class QualityMetricsAdmin(admin.ModelAdmin):
    list_display = [
        'strain_link', 'total_non_acgtns', 'partially_aliased',
        'scorpio_call', 'scorpio_support'
    ]
    list_filter = ['partially_aliased', 'scorpio_call']
    search_fields = ['strain__gisaid_epi_isl']

    def strain_link(self, obj):
        url = reverse('admin:query_engine_virusstrain_change', args=[obj.strain.id])
        return format_html('<a href="{}">{}</a>', url, obj.strain)
    strain_link.short_description = 'Strain'

@admin.register(StructuralVariation)
class StructuralVariationAdmin(admin.ModelAdmin):
    list_display = [
        'strain_link', 'total_insertions',
        'total_frameshifts', 'total_missing'
    ]
    list_filter = ['total_insertions', 'total_frameshifts', 'total_missing']
    search_fields = ['strain__gisaid_epi_isl']

    def strain_link(self, obj):
        url = reverse('admin:query_engine_virusstrain_change', args=[obj.strain.id])
        return format_html('<a href="{}">{}</a>', url, obj.strain)
    strain_link.short_description = 'Strain'
