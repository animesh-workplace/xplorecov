import requests
from datetime import datetime
from snakemake.shell import shell

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

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
print("Started Nextclade")
shell(
    """
        time nextclade run -D {snakemake.input.dataset} -j {snakemake.threads} \
        --output-tsv {snakemake.output.clade_report} --output-all {snakemake.output.clade_folder} \
        {snakemake.input.sequences} {log}
    """
)
print("Finished Nextclade")
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
