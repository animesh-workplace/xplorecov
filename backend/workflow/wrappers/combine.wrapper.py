import pandas
from functools import reduce

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

print("Started Combined")
nextclade = pandas.read_csv(
    snakemake.input.nextclade, delimiter="\t", encoding="utf-8", low_memory=False
)
pangousher = pandas.read_csv(
    snakemake.input.pangolin_usher, delimiter="\t", encoding="utf-8", low_memory=False
)

nextclade.rename(
    columns={
        "seqName": "Name",
        "clade": "Nextclade-Clade",
        "Nextclade_pango": "Nextclade-Lineage",
        "qc.overallScore": "Nextclade-QC-Score",
        "qc.overallStatus": "Nextclade-QC-Status",
    },
    inplace=True,
)

pangousher.rename(
    columns={
        "taxon": "Name",
        "lineage": "Pangousher-Lineage",
        "qc_status": "Pangousher-QC-Status",
    },
    inplace=True,
)

combined_list = [
    nextclade[
        [
            "Name",
            "Nextclade-Clade",
            "Nextclade-Lineage",
            "Nextclade-QC-Score",
            "Nextclade-QC-Status",
        ]
    ],
    pangousher[["Name", "Pangousher-Lineage", "Pangousher-QC-Status"]],
]

combined_report = reduce(
    lambda left, right: pandas.merge(left, right, on="Name", how="inner"), combined_list
)
combined_report.to_csv(snakemake.output.report, sep="\t", index=False)

print("Finished Combined")
