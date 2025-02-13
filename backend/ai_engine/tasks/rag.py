import json, re
from django.apps import apps
from celery import shared_task
from ..llm_chain import AgentChain
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from ..serializers import ChatMessagesSerializer

# from openai import OpenAI
# from guardrails import Guard
# from guardrails.errors import ValidationError
# from guardrails.hub import NSFWText, UnusualPrompt, WebSanitization


def get_schema_dict():
    schema = {}
    for model in apps.get_app_config("query_engine").get_models():
        model_name = model._meta.model_name
        fields = {}
        for field in model._meta.get_fields():
            field_info = {
                "type": field.get_internal_type(),
                "null": getattr(field, "null", False),
                "unique": getattr(field, "unique", False),
            }
            if field.is_relation:
                field_info["related_model"] = (
                    f"{field.related_model._meta.app_label}.{field.related_model._meta.model_name}"
                    if field.related_model
                    else None
                )
            fields[field.name] = field_info
        schema[model_name] = fields
    return schema


@shared_task
def ask_ai_for_code(content, user_analysis_id, user_id, analysis_id):
    db_schema = get_schema_dict()

    chain = AgentChain()
    result = chain.process_request(content, db_schema)

    if result.success:
        result_content = result.data["summary"]
    else:
        result_content = result.message

    chat_message_object = {
        "sender": "assistant",
        "content_type": "text",
        "parent_message_uuid": None,
        "user_analysis": user_analysis_id,
        "content": json.dumps(result_content),
    }
    serializer = ChatMessagesSerializer(data=chat_message_object)
    if serializer.is_valid():
        serializer.save()
        channel_layer = get_channel_layer()
        async_to_sync(channel_layer.group_send)(
            f"{user_id}.{analysis_id}.chat.llm",
            {
                "type": "chat_message",
                "message": json.dumps(result_content),
            },
        )
    else:
        print(serializer.errors)
