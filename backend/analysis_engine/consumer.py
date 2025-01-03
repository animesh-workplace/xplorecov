"""
ASGI config for nsm_server project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/2.2/howto/deployment/asgi/
"""

# import json
# from .tasks import *
from time import sleep

# from asgiref.sync import async_to_sync
# from channels.consumer import AsyncConsumer
# from channels.exceptions import StopConsumer
from channels.generic.websocket import (
    AsyncJsonWebsocketConsumer,
    AsyncWebsocketConsumer,
)


import logging

logger = logging.getLogger(__name__)


class TestConsumer(AsyncWebsocketConsumer):
    async def connect(self):
        logger.info("Connection attempt received")
        task_id = "Workshop"
        try:
            await self.accept()
            logger.info("Connection accepted")
            await self.channel_layer.group_add(task_id, self.channel_name)
            logger.info(f"Added to group {task_id}")
            data = "You have connected to TestConsumer"
            await self.send(text_data=data)  # Changed from send() to send(text_data=)
            logger.info("Initial message sent")
        except Exception as e:
            logger.error(f"Error in connect: {str(e)}")
            raise

    async def websocket_receive(self, event):
        print(f"Received message: {event}")
        # await self.channel_layer.group_send(
        #     "Workshop", {"type": "task_test_message", "message": event["text"]}
        # )

    async def disconnect(self, close_code):
        print(f"Disconnect received with code: {close_code}")
        # task_id = "Workshop"
        # await self.channel_layer.group_discard(task_id, self.channel_name)
        # await self.close()
