import requests, time
from datetime import datetime


rule nextclade:
    input:
        dataset="workflow/resources/nextclade/data",
        sequences=f'{config["OutputDir"]}/uploaded/sequences.fasta',
    output:
        clade_report=f"{config['OutputDir']}/result/nextclade/clade_report.tsv",
        clade_folder=directory(f'{config["OutputDir"]}/result/nextclade/others'),
    log:
        f'{config["OutputDir"]}/log/nextclade.log',
    threads: 10
    run:
        run_websocket_message("Nextclade Analysis Execution", "step4", "start")
        shell(
            """
            time micromamba run -p "/home/nsm/Desktop/All_Development/Manuscript_Work/xplorecov/backend/.workflow-venv/envs/xplorecov" nextclade run \
            -D "{input.dataset}" -j {threads} \
            --output-tsv "{output.clade_report}" --output-all "{output.clade_folder}" \
            "{input.sequences}" > {log} 2>&1
            """
        )
        run_websocket_message("Nextclade Analysis Execution", "step4", "end")

