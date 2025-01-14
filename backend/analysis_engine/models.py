import os
from django.db import models


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


# Model
class UserAnalysis(models.Model):
    user_id = models.UUIDField()
    analysis_id = models.CharField(max_length=21)
    analysis_status = models.JSONField(default=dict)
    submission_date = models.DateTimeField(auto_now_add=True)
    metadata = models.FileField(upload_to=upload_file_location)
    sequence = models.FileField(upload_to=upload_file_location)
    celery_task_id = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self):
        return f"Analysis {self.analysis_id} by User {self.user_id}"


class WebSocketBackendUUID(models.Model):
    uuid = models.UUIDField(unique=True, editable=False)
    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return str(self.uuid)
