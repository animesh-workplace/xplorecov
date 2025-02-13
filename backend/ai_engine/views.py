from rest_framework import status
from .tasks.rag import ask_ai_for_code
from rest_framework.views import APIView
from .cache_manager import LLMCacheManager
from rest_framework.response import Response
from .serializers import ChatMessagesSerializer
from analysis_engine.models import UserAnalysis


# Create your views here.
class AddChatMessagesView(APIView):
    def __init__(self):
        self.cache_manager = LLMCacheManager()

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

        # Try to get from cache (including similar queries)
        cached_response = self.cache_manager.get_cached_response(
            message.get("content", "")
        )
        if cached_response:
            response_data = {
                "modified_query": cached_response["modified_query"],
                "generated_code": cached_response["generated_code"],
                "cache_hit": True,
                "match_type": cached_response["match_type"],
            }

            # Add similarity information if it was a similar match
            if cached_response["match_type"] == "similar":
                response_data.update(
                    {
                        "similarity_score": cached_response["similarity_score"],
                        "original_query": cached_response["original_query"],
                    }
                )

            print(response_data)

            return Response(
                {"message": "Messages found in cache successfully."},
                status=status.HTTP_201_CREATED,
            )
        # # If not in cache, process normally
        # modified_query = self.process_query_with_llm(message)
        # generated_code = self.generate_code(modified_query)

        # Store in cache

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
            if message["sender"] == "human":
                ask_ai_for_code.delay(
                    message.get("content", ""), user_analysis.id, user_id, analysis_id
                )
        else:
            return Response(serializer.errors, status=status.HTTP_400_BAD_REQUEST)

        return Response(
            {"message": "Messages added successfully."}, status=status.HTTP_201_CREATED
        )
