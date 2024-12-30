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
print("Started Usher")
shell(
    """
        time python scripts/parallel_usher.py --sequences {snakemake.input.sequence_folder}/nextclade.aligned.fasta \
        --maskSites {snakemake.input.problem_vcf} --mat-tree {snakemake.input.dataset} \
        --reference {snakemake.input.reference} --sequences-per-group 100 \
        --output {snakemake.output.usher_report} --threads {snakemake.threads} {log}
    """
)
print("Finished Usher")
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
