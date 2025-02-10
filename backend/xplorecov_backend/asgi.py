"""
ASGI config for xplorecov_backend project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/5.1/howto/deployment/asgi/
"""

import os
from channels.routing import (
    URLRouter,
    ProtocolTypeRouter,
)
from django.core.asgi import get_asgi_application
from django.contrib.staticfiles.handlers import ASGIStaticFilesHandler

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "xplorecov_backend.settings")

asgi_app = get_asgi_application()


from django.urls import path
from dotenv import load_dotenv
from django.conf import settings
from channels.auth import AuthMiddlewareStack
from ai_engine.consumer import ChatConsumer
from analysis_engine.consumer import AnalysisConsumer
from channels.security.websocket import AllowedHostsOriginValidator

load_dotenv(settings.BASE_DIR / ".env")

application = ProtocolTypeRouter(
    {
        "http": asgi_app,
        "websocket": AllowedHostsOriginValidator(
            AuthMiddlewareStack(
                URLRouter(
                    [
                        path(
                            os.getenv("BASE_URL"),
                            URLRouter(
                                [
                                    path(
                                        "ws/analysis/<user_id>/<analysis_id>/",
                                        AnalysisConsumer.as_asgi(),
                                        name="analysis-consumer",
                                    ),
                                    path(
                                        "ws/chat/<user_id>/<analysis_id>/",
                                        ChatConsumer.as_asgi(),
                                        name="chat-consumer",
                                    ),
                                ]
                            ),
                        )
                    ]
                )
            )
        ),
    }
)
