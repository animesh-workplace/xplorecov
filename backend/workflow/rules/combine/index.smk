import pandas, time
from functools import reduce


rule combine:
    input:
        nextclade=rules.nextclade.output.clade_report,
        pangolin_usher=rules.pangolin_usher.output.lineage_report,
    output:
        report=f'{config["OutputDir"]}/result/combined_report.tsv',
    log:
        f'{config["OutputDir"]}/log/combined.log',
    threads: 1
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

        combined_report = reduce(
            lambda left, right: pandas.merge(left, right, on="Name", how="inner"),
            [nextclade, pangousher],
        )

        # Generate the data for top voc/vui/voi
        # sequence distribution in country and state
        # sequence count month wise and week wise
        # Update the output to feather to be used 
        combined_report.to_csv(output.report, sep="\t", index=False)
        run_websocket_message("Results Summarization", "step6", "end")
