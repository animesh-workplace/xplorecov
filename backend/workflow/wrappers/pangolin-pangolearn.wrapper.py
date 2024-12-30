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
print("Started Pangolin: Pangolearn")
shell(
    """
        time pangolin {snakemake.input.sequences} --outfile {snakemake.output.lineage_report} \
        -t {snakemake.threads} --analysis-mode fast {log}
    """
)
print("Finished Pangolin: Pangolearn")
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
