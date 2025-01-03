import requests
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
        print("Started Pangolin: Usher")
        run_websocket_message("start")
        shell(
            """
                time micromamba run -p ".workflow-venv/envs/xplorecov" pangolin {input.sequences} --outfile {output.lineage_report} -t {threads} > {log} 2>&1
            """
        )
        run_websocket_message("end")
        print("Finished Pangolin: Usher")
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
