import json, re
from openai import OpenAI
from celery import shared_task
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from ..serializers import ChatMessagesSerializer


@shared_task
def ask_ai_for_code(content, user_analysis_id, user_id, analysis_id):
    client = OpenAI(
        base_url="http://10.10.6.80/xplorecov/ai/code/v1",
        api_key="sk-no-key-required",
    )
    completion = client.chat.completions.create(
        temperature=1,
        model="LLaMA_CPP",
        messages=[
            {
                "role": "system",
                "content": """
                    You are an AI assistant specifically validates if the request is safe and feasible to process. Returns TRUE/FALSE ONLY nothing else.
                """,
            },
            {"role": "user", "content": content},
        ],
    )

    print(completion.choices[0].message.content)

    santized_message_content = re.sub(
        r"\n+", " ", completion.choices[0].message.content
    ).strip()

    chat_message_object = {
        "sender": "assistant",
        "content_type": "text",
        "parent_message_uuid": None,
        "user_analysis": user_analysis_id,
        "content": json.dumps(santized_message_content),
    }
    serializer = ChatMessagesSerializer(data=chat_message_object)
    if serializer.is_valid():
        serializer.save()
        channel_layer = get_channel_layer()
        async_to_sync(channel_layer.group_send)(
            f"{user_id}.{analysis_id}.chat.llm",
            {
                "type": "chat_message",
                "message": santized_message_content,
            },
        )
    else:
        print(serializer.errors)
