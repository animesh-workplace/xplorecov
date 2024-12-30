rule pangolin_pangolearn:
    input:
        sequences = config["OutputDir"] / "uploaded" / "sequences.fasta",
    output:
        lineage_report = config["OutputDir"] / "result" / "pangolin-pangolearn" / "lineage.csv",
    log:
        config["OutputDir"] / "log" / "pangolin-pangolearn.log"
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    threads: 10
    wrapper:
        "file:wrappers/pangolin-pangolearn.wrapper.py"