from .models import UserAnalysis, ToolVersion
from rest_framework.serializers import ModelSerializer


class UserAnalysisSerializer(ModelSerializer):
    class Meta:
        model = UserAnalysis
        fields = ["user_id", "analysis_id", "metadata", "sequence"]


class ToolVersionSerializer(ModelSerializer):
    class Meta:
        model = ToolVersion
        fields = ['nextclade_version', 'pangolin_version', 'constellations_version', 'scorpio_version', 
                  'usher_version', 'gofasta_version', 'minimap2_version', 'faToVcf_version', 'created_at']
