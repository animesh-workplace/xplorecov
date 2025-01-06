import subprocess
from celery import shared_task


@shared_task
def run_snakemake(output_dir):
    try:
        cmd = [
            "snakemake",
            "--snakefile",
            "workflow/Main.Snakefile",
            "--config",
            f"OutputDir={output_dir}",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        return {
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
        }
    except Exception as e:
        return {"success": False, "error": str(e)}
