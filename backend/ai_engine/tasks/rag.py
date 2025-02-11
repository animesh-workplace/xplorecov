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
    print("Resultt", result)

    # if result.success:
    #     print(result.data["summary"])
    # else:
    #     print(f"Error: {result.message}")

    # client = OpenAI(
    #     base_url="http://10.10.6.80/xplorecov/ai/code/v1",
    #     api_key="sk-no-key-required",
    # )
    # completion = client.chat.completions.create(
    #     temperature=1,
    #     model="LLaMA_CPP",
    #     messages=[
    #         {
    #             "role": "system",
    #             "content": """
    #                 You are an AI assistant specifically validates if the request is safe and feasible to process. Returns TRUE/FALSE ONLY nothing else.
    #             """,
    #         },
    #         {"role": "user", "content": content},
    #     ],
    # )

    # print(completion.choices[0].message.content)

    # santized_message_content = re.sub(
    #     r"\n+", " ", completion.choices[0].message.content
    # ).strip()

    # chat_message_object = {
    #     "sender": "assistant",
    #     "content_type": "text",
    #     "parent_message_uuid": None,
    #     "user_analysis": user_analysis_id,
    #     "content": json.dumps(santized_message_content),
    # }
    # serializer = ChatMessagesSerializer(data=chat_message_object)
    # if serializer.is_valid():
    #     serializer.save()
    #     channel_layer = get_channel_layer()
    #     async_to_sync(channel_layer.group_send)(
    #         f"{user_id}.{analysis_id}.chat.llm",
    #         {
    #             "type": "chat_message",
    #             "message": santized_message_content,
    #         },
    #     )
    # else:
    #     print(serializer.errors)
