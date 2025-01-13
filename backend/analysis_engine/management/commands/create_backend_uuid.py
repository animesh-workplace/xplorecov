import uuid
from django.core.management.base import BaseCommand
from analysis_engine.models import WebSocketBackendUUID


class Command(BaseCommand):
    help = "Generate a random UUID for backend WebSocket verification and store it in the database."

    def handle(self, *args, **kwargs):
        # Generate a random UUID
        backend_uuid = str(uuid.uuid4())

        # Save to the database (assuming a model exists)
        WebSocketBackendUUID.objects.create(uuid=backend_uuid)

        # Display the UUID for the user
        self.stdout.write(self.style.SUCCESS(f"Generated UUID: {backend_uuid}"))
        self.stdout.write("Add this UUID to your .env file as:")
        self.stdout.write(f"BACKEND_WEBSOCKET_UUID={backend_uuid}")
