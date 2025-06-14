import requests, time, os
from datetime import datetime


rule nextclade:
    input:
        dataset="workflow/resources/nextclade/data",
        sequences=f'{config["OutputDir"]}/uploaded/sequences.fasta',
    output:
        clade_report=f"{config['OutputDir']}/result/nextclade/clade_report.tsv",
        clade_folder=directory(f'{config["OutputDir"]}/result/nextclade/others'),
        aligned_fasta=f'{config["OutputDir"]}/result/snpeff/aligned.fasta',
    log:
        f'{config["OutputDir"]}/log/nextclade.log',
    threads: 10
    run:
        run_websocket_message("Nextclade Analysis Execution", "step4", "start")
        output_dir = os.path.dirname(output.aligned_fasta)
        os.makedirs(output_dir, exist_ok=True)
        shell(
            """
            time micromamba run -p ".workflow-venv/envs/xplorecov" nextclade run \
            -D "{input.dataset}" -j {threads} \
            --output-tsv "{output.clade_report}" --output-all "{output.clade_folder}" \
            "{input.sequences}" > {log} 2>&1
            cat "{input.dataset}/reference.fasta" "{output.clade_folder}/nextclade.aligned.fasta" > {output.aligned_fasta}
            """
        )
        run_websocket_message("Nextclade Analysis Execution", "step4", "end")
