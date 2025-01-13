"""
ASGI config for xplorecov_backend project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/5.1/howto/deployment/asgi/
"""

import os
<<<<<<< HEAD

from django.core.asgi import get_asgi_application

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'xplorecov_backend.settings')

application = get_asgi_application()
=======
from channels.routing import (
    URLRouter,
    ProtocolTypeRouter,
)
from django.core.asgi import get_asgi_application
from django.contrib.staticfiles.handlers import ASGIStaticFilesHandler

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "xplorecov_backend.settings")

asgi_app = get_asgi_application()


from dotenv import load_dotenv
from django.urls import re_path
from django.conf import settings
from channels.auth import AuthMiddlewareStack
from channels.security.websocket import AllowedHostsOriginValidator
from analysis_engine.consumer import TestConsumer

load_dotenv(settings.BASE_DIR / ".env")

application = ProtocolTypeRouter(
    {
        "http": asgi_app,
        "websocket": AllowedHostsOriginValidator(
            AuthMiddlewareStack(
                URLRouter(
                    [
                        re_path(
                            os.getenv("BASE_URL"),
                            URLRouter(
                                [
                                    re_path(
                                        r"^wsa/backend/$",
                                        TestConsumer.as_asgi(),
                                        name="backend-consumer",
                                    ),
                                    # re_path(r'^wsa/frontend/$', FrontendConsumer.as_asgi(), name='frontend-consumer'),
                                ]
                            ),
                        )
                    ]
                )
            )
        ),
    }
)
>>>>>>> main#1
