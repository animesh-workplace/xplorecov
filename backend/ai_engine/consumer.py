from channels.generic.websocket import AsyncJsonWebsocketConsumer


class ChatConsumer(AsyncJsonWebsocketConsumer):
    async def connect(self):
        # Get user_id and analysis_id from the URL parameters
        user_id = self.scope["url_route"]["kwargs"]["user_id"]
        analysis_id = self.scope["url_route"]["kwargs"]["analysis_id"]

        task_id = f"{user_id}.{analysis_id}.chat.llm"
        try:
            await self.accept()
            await self.channel_layer.group_add(task_id, self.channel_name)
            data = {"message": "You have connected to LLM Real Time Updates"}
            await self.send_json(data)
        except Exception as e:
            raise

    async def websocket_receive(self, event):
        return

    async def chat_message(self, event):
        data = {"message": event["message"]}
        await self.send_json(data)
