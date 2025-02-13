from django.db import models


class LLMCache(models.Model):
    modified_query = models.TextField()
    generated_code = models.TextField()
    hit_count = models.IntegerField(default=1)
    user_query = models.TextField(unique=True)
    query_embedding = models.BinaryField(null=True)
    last_accessed = models.DateTimeField(auto_now=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        indexes = [
            models.Index(fields=["user_query"]),
        ]
