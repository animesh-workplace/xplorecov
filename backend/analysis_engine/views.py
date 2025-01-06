import os
from django.conf import settings
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.response import Response
from .tasks.run_main_workflow import run_snakemake


class SnakemakeView(APIView):
    def post(self, request):
        # Validate input
        user_id = request.data.get("user_id")
        analysis_id = request.data.get("analysis_id")

        if not user_id or not analysis_id:
            return Response(
                {"error": "Both user_id and analysis_id are required"},
                status=status.HTTP_400_BAD_REQUEST,
            )

        # Construct output directory path
        output_dir = f"datalake/{user_id}/{analysis_id}"

        # Trigger Celery task
        task = run_snakemake.delay(output_dir)

        return Response(
            {
                "task_id": task.id,
                "message": "Snakemake task has been queued",
                "status": "pending",
            }
        )
