# from .tasks import *
# from channels.consumer import AsyncConsumer
# from channels.exceptions import StopConsumer
import json
from urllib.parse import parse_qs
from asgiref.sync import sync_to_async
from .models import WebSocketBackendUUID
from channels.generic.websocket import AsyncJsonWebsocketConsumer


class TestConsumer(AsyncJsonWebsocketConsumer):
    async def connect(self):
        # Get user_id and analysis_id from the URL parameters
        user_id = self.scope["url_route"]["kwargs"]["user_id"]
        analysis_id = self.scope["url_route"]["kwargs"]["analysis_id"]

        task_id = f"{user_id}.{analysis_id}"
        try:
            await self.accept()
            await self.channel_layer.group_add(task_id, self.channel_name)
            data = {"message": "You have connected to TestConsumer"}
            await self.send_json(data)
        except Exception as e:
            raise

    async def websocket_receive(self, event):
        # This section is authenticated zone as updates to the database are done here
        is_authenticated = await self.authenticate_backend_service()
        if not is_authenticated:
            return  # Connection closed during authentication
        print(f"Received message: {event}")
        await self.channel_layer.group_send(
            "Workshop", {"type": "task_message", "message": json.loads(event["text"])}
        )

    async def disconnect(self, close_code):
        print(f"Disconnect received with code: {close_code}")

    async def task_message(self, event):
        data = {"message": event["message"]}
        await self.send_json(data)

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
