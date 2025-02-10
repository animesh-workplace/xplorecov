from analysis_engine.models import ChatMessages
from rest_framework.serializers import ModelSerializer


class ChatMessagesSerializer(ModelSerializer):
    class Meta:
        model = ChatMessages
        fields = [
            "uuid",
            "sender",
            "content",
            "created_at",
            "content_type",
            "user_analysis",
            "parent_message_uuid",
        ]
