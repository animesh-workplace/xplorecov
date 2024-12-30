rule pangolin_usher:
    input:
        sequences = config["OutputDir"] / "uploaded" / "sequences.fasta",
    output:
        lineage_report = config["OutputDir"] / "result" / "pangolin-usher" / "lineage.tsv",
    resources:
        tempdir = str(config["OutputDir"] / "result" / "pangolin-usher" / "temp")
    log:
        config["OutputDir"] / "log" / "pangolin-usher.log"
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    threads: 24
    wrapper:
        "file:wrappers/pangolin-usher.wrapper.py"