from .models import UserAnalysis, ToolVersion
from rest_framework.serializers import ModelSerializer, ValidationError

class ToolVersionSerializer(ModelSerializer):
    class Meta:
        model = ToolVersion
        fields = ['nextclade_version', 'pangolin_version', 'constellations_version', 'scorpio_version', 
                  'usher_version', 'gofasta_version', 'minimap2_version', 'faToVcf_version', 'created_at']

class UserAnalysisSerializer(ModelSerializer):
    class Meta:
        model = UserAnalysis
        fields = ["user_id", "analysis_id", "metadata", "sequence", "tool_version", "overall_status", "total_sequences"]

    def create(self, validated_data):
        # Fetch the most recent ToolVersion
        try:
            validated_data['tool_version'] = ToolVersion.objects.latest('created_at')
        except ToolVersion.DoesNotExist:
            raise ValidationError({"tool_version": "No tool versions are available."})
        
        # Create and return the UserAnalysis instance
        return super().create(validated_data)
    

class GetUserAnalysisSerializer(ModelSerializer):
    tool_version = ToolVersionSerializer()

    class Meta:
        model = UserAnalysis
        fields = ["user_id", "analysis_id", "metadata", "sequence", "tool_version", "submission_date"]
