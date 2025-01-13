import subprocess
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
        print(cmd)
        result = subprocess.run(cmd, capture_output=True, text=True)
        return {
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except Exception as e:
        return {"success": False, "error": str(e)}
