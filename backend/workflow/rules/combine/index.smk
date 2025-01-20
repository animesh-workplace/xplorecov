import pandas
import requests, os, json
from functools import reduce


def generate_reports(combined_report):
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
        {
            "graph_type": "None",
            "report_type": "text",
            "name": "Sequences failed to annotate",
            "data": combined_report.loc[
                combined_report["Nextclade-QC-Status"].isna()
                & (combined_report["qc_notes"] == "failed to map"),
                "Name",
            ].to_list(),
        },
        {
            "graph_type": "Bar",
            "name": "Country wise distribution",
            "data": combined_report["Country"].value_counts().to_dict(),
        },
        {
            "graph_type": "Bar",
            "name": "State wise distribution",
            "data": combined_report["State"].value_counts().to_dict(),
        },
        {
            "graph_type": "Bar",
            "name": "Gender distribution",
            "data": combined_report["Gender"].value_counts().to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "State wise distribution of Nextclade-Clade",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Clade"],
                columns=combined_report["State"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Country wise distribution of Nextclade-Clade",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Clade"],
                columns=combined_report["Country"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "State wise distribution of Nextclade-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Lineage"],
                columns=combined_report["State"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Country wise distribution of Nextclade-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Lineage"],
                columns=combined_report["Country"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "State wise distribution of Pangousher-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Pangousher-Lineage"],
                columns=combined_report["State"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Country wise distribution of Pangousher-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Pangousher-Lineage"],
                columns=combined_report["Country"],
            ).to_dict(),
        },
        {
            "graph_type": "Bar",
            "name": "Month wise distribution",
            "data": combined_report["Collection month"]
            .value_counts()
            .sort_index(key=lambda x: pandas.to_datetime(x, format="%b-%Y"))
            .to_dict(),
        },
        {
            "graph_type": "Bar",
            "name": "Week wise distribution",
            "data": (
                combined_report["Collection week"]
                .value_counts()
                .sort_index(
                    key=lambda x: pandas.to_datetime(
                        [i.replace("W", "") + " 1" for i in x], format="%V-%G %w"
                    )
                )
                .to_dict()
            ),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Month wise distribution of Nextclade-Clade",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Clade"],
                columns=combined_report["Collection month"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Week wise distribution of Nextclade-Clade",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Clade"],
                columns=combined_report["Collection week"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Month wise distribution of Nextclade-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Lineage"],
                columns=combined_report["Collection month"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Week wise distribution of Nextclade-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Nextclade-Lineage"],
                columns=combined_report["Collection week"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Month wise distribution of Pangousher-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Pangousher-Lineage"],
                columns=combined_report["Collection month"],
            ).to_dict(),
        },
        {
            "graph_type": "Stacked Bar",
            "name": "Week wise distribution of Pangousher-Lineage",
            "data": pandas.crosstab(
                index=combined_report["Pangousher-Lineage"],
                columns=combined_report["Collection week"],
            ).to_dict(),
        },
        {
            "graph_type": "None",
            "report_type": "text",
            "name": "Conflicts between Nextclade and Pangolin lineage annotaters",
            "data": combined_report[
                combined_report["Nextclade-Lineage"]
                != combined_report["Pangousher-Lineage"]
            ]["Name"].to_list(),
        },
    ]

    return reports


rule combine:
    input:
        nextclade=rules.nextclade.output.clade_report,
        metadata=f'{config["OutputDir"]}/uploaded/metadata.tsv',
        pangolin_usher=rules.pangolin_usher.output.lineage_report,
    output:
        report=f'{config["OutputDir"]}/result/combined_report.tsv',
        report_feather=f'{config["OutputDir"]}/result/combined_report.feather',
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
            input.metadata, delimiter="\t", encoding="utf-8", low_memory=False
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
        metadata["Collection date"] = pandas.to_datetime(metadata["Collection date"])
        metadata["Collection month"] = metadata["Collection date"].dt.strftime("%b-%Y")
        metadata["Collection week"] = metadata["Collection date"].dt.strftime("W%V-%Y")

        combined_report = reduce(
            lambda left, right: pandas.merge(left, right, on="Name", how="inner"),
            [metadata, nextclade, pangousher],
        )

        # Generate the data for top voc/vui/voi
        reports = generate_reports(combined_report)

        try:
            response = requests.post(
                "http://10.10.6.80/xplorecov/api/job/reports/create/",
                data=json.dumps(
                    {
                        "reports": reports,
                        "user_id": config["UserID"],
                        "analysis_id": config["AnalysisID"],
                        "BACKEND_WEBSOCKET_UUID": os.getenv("BACKEND_WEBSOCKET_UUID"),
                    }
                ),
                headers={"Content-Type": "application/json"},
            )
            if response.status_code == 200 or response.status_code == 201:
                print("Reports sent successfully!")
            else:
                print(f"Failed to send reports. Status code: {response.status_code}")
        except Exception as e:
            print(e)

            # Update the output to feather to be used
        combined_report.to_csv(output.report, sep="\t", index=False)
        combined_report.to_feather(output.report_feather)
        run_websocket_message("Results Summarization", "step6", "end")
