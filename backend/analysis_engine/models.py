import os
import uuid
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

    class Meta:
        ordering = ["-submission_date"]


class Report(models.Model):
    GRAPH_TYPE_CHOICES = [
        ("Bar", "Bar"),
        ("Pie", "Pie"),
        ("Line", "Line"),
        ("None", "None"),  # For non-graph reports like text summaries
        ("Stacked Bar", "Stacked Bar"),
    ]

    REPORT_TYPE_CHOICES = [
        ("data", "Data Report"),  # Reports with structured data
        ("text", "Text Summary"),  # Summaries or interpretations of data reports
        ("table", "Table Report"),  # For combined analysis reports
    ]

    user_analysis = models.ForeignKey(
        "UserAnalysis", on_delete=models.CASCADE, related_name="reports"
    )
    data = models.JSONField(default=list)  # Only for data reports
    name = models.CharField(max_length=255)
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


class ChatMessages(models.Model):
    SENDER_CHOICES = [
        ("human", "Human"),
        ("assistant", "Assistant"),
    ]

    CONTENT_TYPE_CHOICES = [
        ("text", "Text"),
        ("rich", "Rich Content"),  # For messages containing images, charts, etc.
    ]

    content = models.JSONField()
    created_at = models.DateTimeField(auto_now_add=True)
    uuid = models.UUIDField(default=uuid.uuid4, unique=True)
    parent_message_uuid = models.UUIDField(null=True, blank=True)
    sender = models.CharField(max_length=20, choices=SENDER_CHOICES)
    content_type = models.CharField(
        max_length=20, choices=CONTENT_TYPE_CHOICES, default="text"
    )
    user_analysis = models.ForeignKey(
        "UserAnalysis", on_delete=models.CASCADE, related_name="chat_messages"
    )

    def __str__(self):
        return f"Message {self.uuid} ({self.sender}) - {self.user_analysis.analysis_id}"


# class CombinedAnalysisReport(models.Model):
#     user_analysis = models.ForeignKey(
#         "UserAnalysis", on_delete=models.CASCADE, related_name="combined_report"
#     )
#     note = models.TextField(blank=True, null=True)
#     sample_type = models.CharField(max_length=100, blank=True, null=True)
#     tag = models.CharField(max_length=100, blank=True, null=True)
#     igsl_id = models.CharField(max_length=100)
#     name = models.CharField(max_length=100)
#     virus_type = models.CharField(max_length=100)
#     passage_details = models.TextField(blank=True, null=True)
#     collection_date = models.DateField(blank=True, null=True)
#     location = models.CharField(max_length=255, blank=True, null=True)
#     country = models.CharField(max_length=100)
#     state = models.CharField(max_length=100, blank=True, null=True)
#     district = models.CharField(max_length=100, blank=True, null=True)
#     additional_location_info = models.TextField(blank=True, null=True)
#     host = models.CharField(max_length=100)
#     host_info = models.TextField(blank=True, null=True)
#     gender = models.CharField(max_length=10, blank=True, null=True)
#     patient_age = models.CharField(max_length=50, blank=True, null=True)
#     patient_status = models.CharField(max_length=100, blank=True, null=True)
#     specimen_source = models.CharField(max_length=255, blank=True, null=True)
#     outbreak = models.CharField(max_length=255, blank=True, null=True)
#     last_vaccinated = models.CharField(max_length=100, blank=True, null=True)
#     treatment = models.TextField(blank=True, null=True)
#     sequencing_technology = models.CharField(max_length=255, blank=True, null=True)
#     assembly_method = models.CharField(max_length=255, blank=True, null=True)
#     coverage = models.CharField(max_length=100, blank=True, null=True)
#     originating_lab = models.CharField(max_length=255, blank=True, null=True)
#     originating_lab_address = models.TextField(blank=True, null=True)
#     submitting_lab = models.CharField(max_length=255, blank=True, null=True)
#     submitting_lab_address = models.TextField(blank=True, null=True)
#     sample_id_submitting_lab = models.CharField(max_length=100)
#     authors = models.TextField(blank=True, null=True)
#     collection_month = models.CharField(max_length=20, blank=True, null=True)
#     collection_week = models.CharField(max_length=20, blank=True, null=True)
#     index = models.IntegerField(blank=True, null=True)

#     # Nextclade clade information
#     nextclade_clade = models.CharField(max_length=100, blank=True, null=True)
#     clade_display = models.CharField(max_length=100, blank=True, null=True)
#     clade_who = models.CharField(max_length=100, blank=True, null=True)
#     clade_nextstrain = models.CharField(max_length=100, blank=True, null=True)

#     partially_aliased = models.BooleanField(default=False)
#     nextclade_lineage = models.CharField(max_length=100, blank=True, null=True)
#     nextclade_qc_score = models.FloatField(blank=True, null=True)
#     nextclade_qc_status = models.CharField(max_length=50, blank=True, null=True)

#     # Mutation summary (simplified)
#     total_substitutions = models.IntegerField(blank=True, null=True)
#     total_deletions = models.IntegerField(blank=True, null=True)
#     total_insertions = models.IntegerField(blank=True, null=True)
#     total_frame_shifts = models.IntegerField(blank=True, null=True)
#     total_missing = models.IntegerField(blank=True, null=True)
#     total_non_acgtns = models.IntegerField(blank=True, null=True)
#     total_aa_substitutions = models.IntegerField(blank=True, null=True)
#     total_aa_deletions = models.IntegerField(blank=True, null=True)
#     total_aa_insertions = models.IntegerField(blank=True, null=True)
#     total_unknown_aa = models.IntegerField(blank=True, null=True)

#     alignment_score = models.FloatField(blank=True, null=True)
#     alignment_start = models.IntegerField(blank=True, null=True)
#     alignment_end = models.IntegerField(blank=True, null=True)
#     cds_coverage = models.TextField(blank=True, null=True)
#     is_reverse_complement = models.BooleanField(default=False)

#     # Mutation details - use TextField to hold lists/strings
#     substitutions = models.TextField(blank=True, null=True)
#     deletions = models.TextField(blank=True, null=True)
#     insertions = models.TextField(blank=True, null=True)
#     frame_shifts = models.TextField(blank=True, null=True)
#     aa_substitutions = models.TextField(blank=True, null=True)
#     aa_deletions = models.TextField(blank=True, null=True)
#     aa_insertions = models.TextField(blank=True, null=True)

#     private_nuc_reversions = models.TextField(blank=True, null=True)
#     private_nuc_labeled = models.TextField(blank=True, null=True)
#     private_nuc_unlabeled = models.TextField(blank=True, null=True)
#     total_private_reversions = models.IntegerField(blank=True, null=True)
#     total_private_labeled = models.IntegerField(blank=True, null=True)
#     total_private_unlabeled = models.IntegerField(blank=True, null=True)
#     total_private = models.IntegerField(blank=True, null=True)

#     # Founder mutations (clade, pango)
#     founder_clade_node = models.CharField(max_length=255, blank=True, null=True)
#     founder_clade_substitutions = models.TextField(blank=True, null=True)
#     founder_clade_deletions = models.TextField(blank=True, null=True)
#     founder_clade_aa_substitutions = models.TextField(blank=True, null=True)
#     founder_clade_aa_deletions = models.TextField(blank=True, null=True)

#     founder_pango_node = models.CharField(max_length=255, blank=True, null=True)
#     founder_pango_substitutions = models.TextField(blank=True, null=True)
#     founder_pango_deletions = models.TextField(blank=True, null=True)
#     founder_pango_aa_substitutions = models.TextField(blank=True, null=True)
#     founder_pango_aa_deletions = models.TextField(blank=True, null=True)

#     # Relative mutations (JN.1, XBB.1.5)
#     relmut_jn1_node = models.CharField(max_length=255, blank=True, null=True)
#     relmut_jn1_substitutions = models.TextField(blank=True, null=True)
#     relmut_jn1_deletions = models.TextField(blank=True, null=True)
#     relmut_jn1_aa_substitutions = models.TextField(blank=True, null=True)
#     relmut_jn1_aa_deletions = models.TextField(blank=True, null=True)

#     relmut_xbb15_node = models.CharField(max_length=255, blank=True, null=True)
#     relmut_xbb15_substitutions = models.TextField(blank=True, null=True)
#     relmut_xbb15_deletions = models.TextField(blank=True, null=True)
#     relmut_xbb15_aa_substitutions = models.TextField(blank=True, null=True)
#     relmut_xbb15_aa_deletions = models.TextField(blank=True, null=True)

#     # QC sections (missing data, mixed sites, etc.)
#     unknown_aa_ranges = models.TextField(blank=True, null=True)
#     non_acgtns = models.TextField(blank=True, null=True)

#     qc_missing_score = models.FloatField(blank=True, null=True)
#     qc_missing_status = models.CharField(max_length=50, blank=True, null=True)
#     qc_missing_total = models.IntegerField(blank=True, null=True)

#     qc_mixed_score = models.FloatField(blank=True, null=True)
#     qc_mixed_status = models.CharField(max_length=50, blank=True, null=True)
#     qc_mixed_total = models.IntegerField(blank=True, null=True)

#     # Continue adding qc fields as per need...

#     pangousher_lineage = models.CharField(max_length=100, blank=True, null=True)
#     scorpio_call = models.CharField(max_length=255, blank=True, null=True)
#     scorpio_support = models.FloatField(blank=True, null=True)
#     scorpio_conflict = models.FloatField(blank=True, null=True)
#     scorpio_notes = models.TextField(blank=True, null=True)

#     pangolin_version = models.CharField(max_length=50, blank=True, null=True)
#     scorpio_version = models.CharField(max_length=50, blank=True, null=True)
#     constellation_version = models.CharField(max_length=50, blank=True, null=True)

#     is_designated = models.BooleanField(default=False)
#     pangousher_qc_status = models.CharField(max_length=50, blank=True, null=True)
#     qc_notes = models.TextField(blank=True, null=True)
#     note_field = models.TextField(blank=True, null=True)
