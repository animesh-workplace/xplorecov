"""
ASGI config for nsm_server project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/2.2/howto/deployment/asgi/
"""

import json

# from .tasks import *
from time import sleep

# from asgiref.sync import async_to_sync
# from channels.consumer import AsyncConsumer
# from channels.exceptions import StopConsumer
from channels.generic.websocket import (
    AsyncJsonWebsocketConsumer,
    AsyncWebsocketConsumer,
)


class TestConsumer(AsyncJsonWebsocketConsumer):
    async def connect(self):
        task_id = "Workshop"
        try:
            await self.accept()
            await self.channel_layer.group_add(task_id, self.channel_name)
            data = {"message": "You have connected to TestConsumer"}
            await self.send_json(data)
        except Exception as e:
            raise

    async def websocket_receive(self, event):
        print(f"Received message: {event}")
        await self.channel_layer.group_send(
            "Workshop", {"type": "task_message", "message": json.loads(event["text"])}
        )

    async def disconnect(self, close_code):
        print(f"Disconnect received with code: {close_code}")

    async def task_message(self, event):
        data = {"message": event["message"]}
        await self.send_json(data)
