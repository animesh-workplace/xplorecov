import io, time, pandas, tempfile
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.response import Response
from .serializers import CSVUploadSerializer
from .tasks.train_data import train_ai_model
from .tasks.upload_data import process_csv_file
from .tasks.generate_sql_command import query_using_ai_model

class VirusDataUploadView(APIView):
    def post(self, request):
        serializer = CSVUploadSerializer(data=request.data)
        if not serializer.is_valid():
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

        csv_file = serializer.validated_data['file']
        csv_content = csv_file.read().decode('utf-8')

        df = pandas.read_csv(io.StringIO(csv_content), engine="pyarrow")
        # tempdir = TemporaryDirectory()
        tempdir = tempfile.mkdtemp() 
        metadata_loc = f"{tempdir}/sars_metadata.feather"
        df.to_feather(metadata_loc)
        task = process_csv_file.delay(metadata_loc)

        return Response({"message": "Job Submitted", "task": task.id}, status=status.HTTP_200_OK)


class AIModelTrainView(APIView):
    def post(self, request):
        task = train_ai_model.delay()
        return Response({"message": "Job Submitted", "task": task.id}, status=status.HTTP_200_OK)

class AIModelQueryView(APIView):
    def post(self, request):
        task = query_using_ai_model.delay(request.data["question"])
        # return Response({"message": "Job Submitted"}, status=status.HTTP_200_OK)
        return Response({"message": "Job Submitted", "task": task.id}, status=status.HTTP_200_OK)


