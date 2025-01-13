<<<<<<< HEAD
from django.db import models

# Create your models here.
=======
import os
from django.db import models


def upload_file_location(instance, filename):
    """
    Generate the file path for uploaded files.
    Path: datalake/<user_id>/<analysis_id>/uploaded/<filename>
    """
    return os.path.join(
        "datalake",
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
>>>>>>> main#1
