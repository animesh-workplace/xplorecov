rule all:
    input:
        # f"{config['OutputDir']}/result/nextclade/clade_report.tsv",
        # f"{config['OutputDir']}/result/pangolin-usher/lineage.tsv",
        f'{config["OutputDir"]}/result/combined_report.tsv',


include: "rules/annotation/nextclade.smk"
include: "rules/annotation/pangolin-usher.smk"
include: "rules/combine/index.smk"
