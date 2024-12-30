rule all:
    input:
        config["OutputDir"] / "result" / "combined_report.tsv",


include: "rules/annotation/nextclade.smk"
include: "rules/annotation/pangolin-usher.smk"
include: "rules/annotation/pangolin-pangolearn.smk"
include: "rules/combine/index.smk"
