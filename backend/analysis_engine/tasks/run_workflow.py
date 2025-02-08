import subprocess, json, re
from openai import OpenAI
from datetime import datetime
from celery import shared_task
from django.conf import settings
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer
from ..serializers import ChatMessagesSerializer


@shared_task
def run_analysis_workflow(output_dir, user_id, analysis_id):
    try:
        cmd = [
            "snakemake",
            "--snakefile",
            "workflow/Main.Snakefile",
            "--config",
            f"UserID={user_id}",
            f"OutputDir={output_dir}",
            f"AnalysisID={analysis_id}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stderr, result.stdout)
        return {
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


@shared_task
def run_update_workflow():
    try:
        cmd = [
            "snakemake",
            "--snakefile",
            "workflow/Update.Snakefile",
            "--config",
            f"Date={datetime.now().strftime('%d-%m-%Y')}",
            "UpdateDir=workflow",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        print(result.stderr, result.stdout)
        return {
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except Exception as e:
        return {"success": False, "error": str(e)}


@shared_task
def ask_ai_for_code(content, user_analysis_id, user_id, analysis_id):
    client = OpenAI(
        base_url="http://10.10.6.80/xplorecov/ai/content/v1",
        api_key="sk-no-key-required",
    )
    completion = client.chat.completions.create(
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
