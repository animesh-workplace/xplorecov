import json
from urllib.parse import parse_qs
from asgiref.sync import sync_to_async
from django.forms.models import model_to_dict
from .models import WebSocketBackendUUID, UserAnalysis
from channels.generic.websocket import AsyncJsonWebsocketConsumer


class AnalysisConsumer(AsyncJsonWebsocketConsumer):
    async def connect(self):
        # Get user_id and analysis_id from the URL parameters
        user_id = self.scope["url_route"]["kwargs"]["user_id"]
        analysis_id = self.scope["url_route"]["kwargs"]["analysis_id"]

        existing_status = await self.get_analysis_status(user_id, analysis_id)

        task_id = f"{user_id}.{analysis_id}"
        try:
            await self.accept()
            await self.channel_layer.group_add(task_id, self.channel_name)
            data = {"message": "You have connected to Analysis Real Time Updates", **existing_status}
            await self.send_json(data)
        except Exception as e:
            raise

    async def websocket_receive(self, event):
        # This section is authenticated zone as updates to the database are done here``
        is_authenticated = await self.authenticate_backend_service()
        if not is_authenticated:
            return  # Connection closed during authentication
        
        # Validate the expected fields in the message
        user_id = self.scope["url_route"]["kwargs"]["user_id"]
        analysis_id = self.scope["url_route"]["kwargs"]["analysis_id"]
        status_update = json.loads(event["text"])
        task_id = f"{user_id}.{analysis_id}"

        if not user_id or not analysis_id or not status_update:
            print("Missing required fields in the message")
            return

        # Update the UserAnalysis model
        await self.update_analysis_status(user_id, analysis_id, status_update)

        await self.channel_layer.group_send(
            task_id, {"type": "task_message", "message": json.loads(event["text"])}
        )

    async def task_message(self, event):
        data = {"message": event["message"]}
        await self.send_json(data)

    @sync_to_async
    def get_analysis_status(self, user_id, analysis_id):
        try:
            analysis = UserAnalysis.objects.get(user_id=user_id, analysis_id=analysis_id)
            return {"current_status": analysis.analysis_status, "tools_used": {
                "tools": model_to_dict(analysis.tool_version, exclude=['id']),
                'updated_at': analysis.tool_version.created_at.strftime('%d-%m-%Y %I:%M %p')
            }}
        except UserAnalysis.DoesNotExist:
            return None


    @sync_to_async
    def update_analysis_status(self, user_id, analysis_id, status_update):
        try:
            analysis = UserAnalysis.objects.get(user_id=user_id, analysis_id=analysis_id)

            if not isinstance(analysis.analysis_status, list):
                analysis.analysis_status = [status_update]
            else:
                # Append the new status to the existing `analysis_status`
                analysis.analysis_status.append(status_update)

            # Save the changes
            analysis.save()

        except UserAnalysis.DoesNotExist:
            print(f"User Analysis with user_id = {user_id} and analysis_id = {analysis_id} not found")
        except Exception as e:
            print(f"Error updating analysis status: {e}")

    async def authenticate_backend_service(self):
        query_string = self.scope.get("query_string", b"").decode("utf-8")
        query_params = parse_qs(query_string)
        BACKEND_WEBSOCKET_UUID = query_params.get("BACKEND_WEBSOCKET_UUID", [None])[0]

        if not BACKEND_WEBSOCKET_UUID:
            # Reject the connection if the UUID is missing
            await self.send_json("Unauthorized")
            await self.close()
            return False

        # Check the database for the UUID using the ORM
        is_valid_uuid = await self.is_uuid_valid(BACKEND_WEBSOCKET_UUID)

        if not is_valid_uuid:
            # Reject the connection if the UUID is not valid
            await self.send_json("Unauthorized")
            await self.close()
            return False

        # Connection is authenticated
        return True

    @sync_to_async
    def is_uuid_valid(self, uuid):
        # Check if the UUID exists in the WebSocketBackendUUID table
        return WebSocketBackendUUID.objects.filter(uuid=uuid).exists()
