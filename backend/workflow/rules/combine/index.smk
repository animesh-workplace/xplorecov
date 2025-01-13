<<<<<<< HEAD
import pandas
=======
import pandas, time
>>>>>>> main#1
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
        print("Started Combined")
<<<<<<< HEAD
=======
        run_websocket_message("combine", "start")
>>>>>>> main#1
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
<<<<<<< HEAD
        print(nextclade)
=======
>>>>>>> main#1
        pangousher.rename(
            columns={
                "taxon": "Name",
                "lineage": "Pangousher-Lineage",
                "qc_status": "Pangousher-QC-Status",
            },
            inplace=True,
        )
<<<<<<< HEAD
        print(pangousher)
=======
>>>>>>> main#1

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
            lambda left, right: pandas.merge(left, right, on="Name", how="inner"),
            combined_list,
        )
        combined_report.to_csv(output.report, sep="\t", index=False)
<<<<<<< HEAD
=======
        run_websocket_message("combine", "end")
>>>>>>> main#1
        print("Finished Combined")
