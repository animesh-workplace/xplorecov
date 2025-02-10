import subprocess, json, re
from openai import OpenAI
from datetime import datetime
from celery import shared_task
from django.conf import settings
from asgiref.sync import async_to_sync
from channels.layers import get_channel_layer

# from ..serializers import ChatMessagesSerializer


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
