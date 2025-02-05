import subprocess
from openai import OpenAI
from datetime import datetime
from celery import shared_task


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
def run_ask_ai():
    client = OpenAI(base_url="http://10.10.6.80/ai/code/v1", api_key="no-key-required")
    completion = client.chat.completions.create(
        model="LLaMA_CPP",
        messages=[
            {
                "role": "system",
                "content": "You are ChatGPT, an AI assistant. Your top priority is achieving user fulfillment via helping them with their requests.",
            },
            {"role": "user", "content": "Write a limerick about python exceptions"},
        ],
    )
