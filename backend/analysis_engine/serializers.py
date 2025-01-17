from rest_framework.serializers import ModelSerializer, ValidationError
from .models import (
    Report,
    ToolVersion,
    UserAnalysis,
    ChatMessages,
    WebSocketBackendUUID,
)


class ToolVersionSerializer(ModelSerializer):
    class Meta:
        model = ToolVersion
        fields = [
            "nextclade_version",
            "pangolin_version",
            "constellations_version",
            "scorpio_version",
            "usher_version",
            "gofasta_version",
            "minimap2_version",
            "faToVcf_version",
        ]

    def validate(self, data):
        # Check if the key 'BACKEND_WEBSOCKET_UUID' exists in the initial data
        websocket_uuid = self.initial_data.get("BACKEND_WEBSOCKET_UUID")
        if websocket_uuid is None:
            raise ValidationError({"detail": "BACKEND_WEBSOCKET_UUID is required."})

        # Verify if the uuid exists in WebSocketBackendUUID model
        if not WebSocketBackendUUID.objects.filter(uuid=websocket_uuid).exists():
            raise ValidationError({"detail": "Invalid BACKEND_WEBSOCKET_UUID."})

        return data


class UserAnalysisSerializer(ModelSerializer):
    class Meta:
        model = UserAnalysis
        fields = [
            "user_id",
            "analysis_id",
            "metadata",
            "sequence",
            "tool_version",
            "overall_status",
            "total_sequences",
        ]

    def create(self, validated_data):
        # Fetch the most recent ToolVersion
        try:
            validated_data["tool_version"] = ToolVersion.objects.latest("created_at")
        except ToolVersion.DoesNotExist:
            raise ValidationError({"tool_version": "No tool versions are available."})

        # Create and return the UserAnalysis instance
        return super().create(validated_data)


class ReportSerializer(ModelSerializer):
    class Meta:
        model = Report
        fields = [
            "id",
            "name",
            "data",
            "text_summary",
            "graph_type",
            "report_type",
            "created_at",
        ]


class GetUserAnalysisSerializer(ModelSerializer):
    tool_version = ToolVersionSerializer()

    class Meta:
        model = UserAnalysis
        fields = [
            "user_id",
            "analysis_id",
            "metadata",
            "sequence",
            "tool_version",
            "submission_date",
            "total_sequences",
            "overall_status",
            "completion_date",
            "expiration_date",
        ]


class ChatMessagesSerializer(ModelSerializer):
    class Meta:
        model = ChatMessages
        fields = [
            "uuid",
            "content",
            "created_at",
            "parent_message_uuid",
            "sender",
            "content_type",
        ]
