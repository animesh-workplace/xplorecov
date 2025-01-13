from .models import UserAnalysis
from rest_framework.serializers import ModelSerializer


class UserAnalysisSerializer(ModelSerializer):
    class Meta:
        model = UserAnalysis
        fields = ["user_id", "analysis_id", "metadata", "sequence"]
