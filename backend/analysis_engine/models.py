import os
from django.db import models
from django.utils.timezone import now, timedelta, localtime


def upload_file_location(instance, filename):
    """
    Generate the file path for uploaded files.
    Path: MEDIA_ROOT/<user_id>/<analysis_id>/uploaded/<filename>
    """
    return os.path.join(
        str(instance.user_id),
        str(instance.analysis_id),
        "uploaded",
        filename,
    )


def get_expiration_date():
    return now() + timedelta(days=14)


def default_analysis_status():
    return [
        {
            "status": "start",
            "step_id": "step3",
            "step_name": "Queuing Analysis",
            "timestamp": localtime().isoformat(),
        }
    ]


# Model
class WebSocketBackendUUID(models.Model):
    uuid = models.UUIDField(unique=True, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return str(self.uuid)


class ToolVersion(models.Model):
    nextclade_version = models.CharField(max_length=255, null=True, blank=True)
    pangolin_version = models.CharField(max_length=255, null=True, blank=True)
    constellations_version = models.CharField(max_length=255, null=True, blank=True)
    scorpio_version = models.CharField(max_length=255, null=True, blank=True)
    usher_version = models.CharField(max_length=255, null=True, blank=True)
    gofasta_version = models.CharField(max_length=255, null=True, blank=True)
    minimap2_version = models.CharField(max_length=255, null=True, blank=True)
    faToVcf_version = models.CharField(max_length=255, null=True, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        ist_time = localtime(self.created_at)
        return f"Tool Versions - {ist_time.strftime('%d-%m-%Y %I:%M %p')}"


class UserAnalysis(models.Model):
    STATUS_CHOICES = [
        ("PENDING", "Pending"),
        ("ERROR", "Error"),
        ("SUCCESS", "Success"),
    ]

    user_id = models.UUIDField()
    analysis_id = models.CharField(max_length=21)
    # Overall Status will contain only these options
    overall_status = models.CharField(
        max_length=20, null=True, blank=True, choices=STATUS_CHOICES, default="PENDING"
    )
    # The JSON structure is an array of objects, where each object represents an analysis step.
    # Each object contains the following keys:
    # - "step_name": A string representing the name of the step.
    # - "step_id": Identification for the step to be used in teh frontend
    # - "status": Containing either "start" or "end"
    # - "timestamp": A string indicating the timestamp when the step was last updated.
    submission_date = models.DateTimeField(auto_now_add=True)
    metadata = models.FileField(upload_to=upload_file_location)
    sequence = models.FileField(upload_to=upload_file_location)
    total_sequences = models.IntegerField(null=True, blank=True)
    completion_date = models.DateTimeField(null=True, blank=True)
    analysis_status = models.JSONField(default=default_analysis_status)
    expiration_date = models.DateTimeField(default=get_expiration_date)
    celery_task_id = models.CharField(max_length=255, blank=True, null=True)
    tool_version = models.ForeignKey(
        ToolVersion,
        on_delete=models.CASCADE,
        related_name="analyses",
        null=True,
        blank=True,
    )

    def __str__(self):
        return f"Analysis {self.analysis_id} by User {self.user_id}"


class Report(models.Model):
    GRAPH_TYPE_CHOICES = [
        ("Bar", "Bar"),
        ("Line", "Line"),
        ("Pie", "Pie"),
        ("None", "None"),  # For non-graph reports like text summaries
    ]

    REPORT_TYPE_CHOICES = [
        ("data", "Data Report"),  # Reports with structured data
        ("text", "Text Summary"),  # Summaries or interpretations of data reports
    ]

    user_analysis = models.ForeignKey(
        "UserAnalysis", on_delete=models.CASCADE, related_name="reports"
    )
    name = models.CharField(max_length=255)
    data = models.JSONField(null=True, blank=True)  # Only for data reports
    created_at = models.DateTimeField(auto_now_add=True)
    text_summary = models.TextField(null=True, blank=True)  # Only for text summaries
    graph_type = models.CharField(
        max_length=20, choices=GRAPH_TYPE_CHOICES, default="None"
    )
    report_type = models.CharField(
        max_length=20, choices=REPORT_TYPE_CHOICES, default="data"
    )
    summary_sources = models.ManyToManyField(
        "self",  # Self-referential relationship
        blank=True,  # Allows reports without summaries
        symmetrical=False,  # Allows directional relationships
        related_name="summarized_by",  # Reverse relation for querying summaries
    )

    def __str__(self):
        return f"{self.name} ({self.name}) - {self.user_analysis.analysis_id}"
