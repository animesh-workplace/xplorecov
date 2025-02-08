import ollama
from django.conf import settings
from rest_framework import status
from rest_framework.views import APIView
from rest_framework.response import Response
from .models import UserAnalysis, Report, ChatMessages
from rest_framework.parsers import MultiPartParser, FormParser
from .tasks.run_workflow import (
    ask_ai_for_code,
    run_update_workflow,
    run_analysis_workflow,
)
from .serializers import (
    AIChatSerializer,
    ReportSerializer,
    ToolVersionSerializer,
    ChatMessagesSerializer,
    UserAnalysisSerializer,
    GetUserAnalysisSerializer,
)


class ToolUpdateView(APIView):
    def post(self, request):

        # Trigger Celery task
        task = run_update_workflow.delay()

        return Response(
            {
                "task_id": task.id,
                "status": "pending",
                "message": "Update workflow task has been queued",
            }
        )


class CreateUserAnalysisView(APIView):
    parser_classes = (MultiPartParser, FormParser)

    def post(self, request, *args, **kwargs):
        serializer = UserAnalysisSerializer(data=request.data)
        if serializer.is_valid():
            analysis = serializer.save()

            user_id = request.data.get("user_id")
            analysis_id = request.data.get("analysis_id")

            if not user_id or not analysis_id:
                return Response(
                    {"error": "Both user_id and analysis_id are required"},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            # Construct output directory path
            output_dir = settings.MEDIA_ROOT / user_id / analysis_id

            # Trigger Celery task
            task = run_analysis_workflow.delay(str(output_dir), user_id, analysis_id)

            return Response(
                {"message": "Analysis submitted successfully.", "id": analysis.id},
                status=status.HTTP_201_CREATED,
            )
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class GetUserAnalysisView(APIView):
    def get(self, request, *args, **kwargs):
        user_id = request.query_params.get("user_id")

        if not user_id:
            return Response(
                {"error": "user_id is required"}, status=status.HTTP_400_BAD_REQUEST
            )

        # Query all analyses for the given user_id
        analyses = UserAnalysis.objects.filter(user_id=user_id).order_by(
            "-submission_date"
        )

        if not analyses.exists():
            return Response(
                {"message": "No analyses found for this user."},
                status=status.HTTP_404_NOT_FOUND,
            )

        # Serialize the analyses data
        serializer = GetUserAnalysisSerializer(analyses, many=True)

        return Response(serializer.data, status=status.HTTP_200_OK)


class GetAnalysisDetailView(APIView):
    def get(self, request, *args, **kwargs):
        user_id = request.query_params.get("user_id")
        analysis_id = request.query_params.get("analysis_id")

        if not analysis_id or not user_id:
            return Response(
                {"error": "Both user_id and analysis_id are required"},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            # Retrieve the specific analysis by ID
            analysis = UserAnalysis.objects.prefetch_related(
                "reports", "chat_messages"
            ).get(user_id=user_id, analysis_id=analysis_id)
        except UserAnalysis.DoesNotExist:
            return Response(
                {"error": "Analysis not found"},
                status=status.HTTP_404_NOT_FOUND,
            )

        # Serialize the analysis data
        analysis_serializer = GetUserAnalysisSerializer(analysis)
        reports_serializer = ReportSerializer(analysis.reports.all(), many=True)
        chat_messages_serializer = ChatMessagesSerializer(
            analysis.chat_messages.all(), many=True
        )

        return Response(
            {
                **analysis_serializer.data,
                "graph_reports": [
                    report
                    for report in reports_serializer.data
                    if report["graph_type"] != "None"
                ],
                "summary_reports": [
                    report
                    for report in reports_serializer.data
                    if report["graph_type"] == "None"
                ],
                "chat_messages": chat_messages_serializer.data,
            },
            status=status.HTTP_200_OK,
        )


class ToolVersionCreateView(APIView):
    def post(self, request):
        # Extract the tool version information from the request data
        tool_versions = request.data

        # Prepare the data to save in the model
        data = {
            "usher_version": tool_versions.get("usher"),
            "gofasta_version": tool_versions.get("gofasta"),
            "faToVcf_version": tool_versions.get("fatovcf"),
            "scorpio_version": tool_versions.get("scorpio"),
            "pangolin_version": tool_versions.get("pangolin"),
            "minimap2_version": tool_versions.get("minimap2"),
            "nextclade_version": tool_versions.get("nextclade"),
            "constellations_version": tool_versions.get("constellations"),
            "BACKEND_WEBSOCKET_UUID": tool_versions.get("BACKEND_WEBSOCKET_UUID"),
        }

        # Use the serializer to validate and save the data
        serializer = ToolVersionSerializer(data=data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=status.HTTP_201_CREATED)
        return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)


class AddBulkReportsView(APIView):
    def post(self, request, *args, **kwargs):
        user_id = request.data.get("user_id")
        analysis_id = request.data.get("analysis_id")
        reports = request.data.get("reports", [])

        if not user_id or not analysis_id or not reports:
            return Response(
                {"error": "user_id, analysis_id, and reports are required."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            user_analysis = UserAnalysis.objects.get(
                user_id=user_id, analysis_id=analysis_id
            )
        except UserAnalysis.DoesNotExist:
            return Response(
                {"error": "User analysis not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        # Prepare Report objects for bulk creation
        report_objects = [
            Report(
                name=report["name"],
                data=report.get("data"),
                user_analysis=user_analysis,
                text_summary=report.get("text_summary"),
                graph_type=report.get("graph_type", "None"),
                report_type=report.get("report_type", "data"),
            )
            for report in reports
        ]

        # Bulk create the reports
        Report.objects.bulk_create(report_objects)

        return Response(
            {"message": "Reports added successfully."}, status=status.HTTP_201_CREATED
        )


class AddChatMessagesView(APIView):
    def post(self, request, *args, **kwargs):
        user_id = request.data.get("user_id")
        analysis_id = request.data.get("analysis_id")
        message = request.data.get("message")

        if not user_id or not analysis_id or not message:
            return Response(
                {"error": "user_id, analysis_id, and message are required."},
                status=status.HTTP_400_BAD_REQUEST,
            )

        try:
            user_analysis = UserAnalysis.objects.get(
                user_id=user_id, analysis_id=analysis_id
            )
        except UserAnalysis.DoesNotExist:
            return Response(
                {"error": "User analysis not found."},
                status=status.HTTP_404_NOT_FOUND,
            )

        # Prepare ChatMessages objects for creation
        chat_message_object = {
            "user_analysis": user_analysis.id,
            "content": message.get("content", ""),
            "sender": message.get("sender", "human"),
            "content_type": message.get("content_type", "text"),
            "parent_message_uuid": message.get("parent_message_uuid", None),
        }

        # Use the serializer to validate and save the data
        serializer = ChatMessagesSerializer(data=chat_message_object)

        if serializer.is_valid():
            serializer.save()
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

        if message["sender"] == "human":
            ai_message = ask_ai_for_code.delay(
                message.get("content", ""), user_analysis.id, user_id, analysis_id
            )
            print(ai_message)

        return Response(
            {"message": "Messages added successfully."}, status=status.HTTP_201_CREATED
        )


# class GenerateAICodeView(APIView):
#     def post(self, request):
#         serializer = AIChatSerializer(data=request.data)
#         if serializer.is_valid():
#             try:
#                 response = ollama.chat(
#                     model="xplorecov-coder",
#                     messages=[
#                         {
#                             "role": "user",
#                             "content": serializer.validated_data["content"],
#                         }
#                     ],
#                 )
#                 return Response(
#                     {"response": response["message"]["content"]},
#                     status=status.HTTP_200_OK,
#                 )
#             except Exception as e:
#                 return Response(
#                     {"error": str(e)}, status=status.HTTP_500_INTERNAL_SERVER_ERROR
#                 )
#         return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)
