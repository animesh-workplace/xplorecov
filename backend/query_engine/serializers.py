from rest_framework import serializers
from .models import VirusStrain, Location, AAMutation, NucMutation, StructuralVariation, QualityMetrics

class CSVUploadSerializer(serializers.Serializer):
	file = serializers.FileField(required=True)

	def validate_file(self, csv_file):
		if not csv_file.name.endswith('.csv'):
			raise serializers.ValidationError("File must be CSV format")
		return csv_file