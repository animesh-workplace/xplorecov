import requests
from datetime import datetime

shell.executable("bash")


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
        print("Started Nextclade")
        shell(
            """
            time micromamba run -p "~/micromamba-env/.workflow-venv/envs/nibmg_tool" nextclade run \
            -D "{input.dataset}" -j {threads} \
            --output-tsv "{output.clade_report}" --output-all "{output.clade_folder}" \
            "{input.sequences}" > {log} 2>&1
            """
        )
        print("Finished Nextclade")
        # requests.post(
        #     "http://localhost:5000/print",
        #     data={
        #         "task": "Qualimap BamQC",
        #         "status": "Started",
        #         "time": str(datetime.now()),
        #         "db_loc": snakemake.config["db_loc"],
        #         "sample": f"{snakemake.wildcards.sample}",
        #     },
        # )
        # requests.post(
        #     "http://localhost:5000/print",
        #     data={
        #         "task": "Qualimap BamQC",
        #         "status": "Finished",
        #         "time": str(datetime.now()),
        #         "db_loc": snakemake.config["db_loc"],
        #         "sample": f"{snakemake.wildcards.sample}",
        #     },
        # )

