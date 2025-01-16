import pandas
from functools import reduce


rule combine:
    input:
        nextclade=rules.nextclade.output.clade_report,
        metadata=f'{config["OutputDir"]}/uploaded/metadata.tsv',
        pangolin_usher=rules.pangolin_usher.output.lineage_report,
    output:
        report=f'{config["OutputDir"]}/result/combined_report.tsv',
    log:
        f'{config["OutputDir"]}/log/combined.log',
    run:
        run_websocket_message("Results Summarization", "step6", "start")
        nextclade = pandas.read_csv(
            input.nextclade,
            delimiter="\t",
            encoding="utf-8",
            low_memory=False,
        )
        pangousher = pandas.read_csv(
            input.pangolin_usher,
            delimiter=",",
            encoding="utf-8",
            low_memory=False,
        )
        metadata = pandas.read_csv(
            input.metadata,
            delimiter="\t",
            encoding="utf-8",
            low_memory=False,
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

        metadata.rename(columns={"Virus name": "Name"}, inplace=True)

        combined_report = reduce(
            lambda left, right: pandas.merge(left, right, on="Name", how="inner"),
            [metadata, nextclade, pangousher],
        )

        # Generate the data for top voc/vui/voi
        reports = [
            {
                "graph_type": "Bar",
                "name": "Top Nextclade's Clade",
                "data": combined_report["Nextclade-Clade"].value_counts().to_dict(),
            },
            {
                "graph_type": "Bar",
                "name": "Top Nextclade's Pangolineage",
                "data": combined_report["Nextclade-Lineage"].value_counts().to_dict(),
            },
            {
                "graph_type": "Bar",
                "name": "Top Pangolin's Lineage",
                "data": combined_report["Pangousher-Lineage"].value_counts().to_dict(),
            },
            {
                "graph_type": "Bar",
                "name": "Top Pangolin's Scorpio Call",
                "data": combined_report["scorpio_call"].value_counts().to_dict(),
            },
            {
                "graph_type": "Bar",
                "name": "Top WHO's Clade",
                "data": combined_report["clade_who"].value_counts().to_dict(),
            },
        ]
        # sequence distribution in country and state
        # sequence count month wise and week wise
        # Update the output to feather to be used
        combined_report.to_csv(output.report, sep="\t", index=False)
        run_websocket_message("Results Summarization", "step6", "end")
