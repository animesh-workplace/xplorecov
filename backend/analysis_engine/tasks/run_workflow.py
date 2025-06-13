import subprocess
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
            "--cores",
            "all",
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
            "--cores",
            "all",
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
