import requests, time, os
from datetime import datetime


rule pangolin_usher:
    input:
        sequences=f'{config["OutputDir"]}/uploaded/sequences.fasta',
    output:
        lineage_report=f"{config['OutputDir']}/result/pangolin-usher/lineage.tsv",
    resources:
        tempdir=f"{config['OutputDir']}/result/pangolin-usher/temp",
    log:
        f'{config["OutputDir"]}/log/pangolin-usher.log',
    threads: 10
    run:
        run_websocket_message("Pangolin Analysis Execution", "step5", "start")
        output_dir = os.path.dirname(output.lineage_report)
        os.makedirs(output_dir, exist_ok=True)
        shell(
            """
                time micromamba run -p ".workflow-venv/envs/xplorecov" pangolin {input.sequences} --outfile {output.lineage_report} -t {threads} > {log} 2>&1
            """
        )
        run_websocket_message("Pangolin Analysis Execution", "step5", "end")
