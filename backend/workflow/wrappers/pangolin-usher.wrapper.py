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
print("Started Pangolin: Usher")
shell(
    """
        time python scripts/parallel_pango_usher.py --sequences {snakemake.input.sequences} \
        --threads {snakemake.threads} --sequences-per-group 100 --output {snakemake.output.lineage_report} {log}
    """
)
print("Finished Pangolin: Usher")
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
