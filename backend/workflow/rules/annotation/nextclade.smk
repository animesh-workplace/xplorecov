rule nextclade:
    input:
        dataset = "resources/nextclade/data",
        sequences = config["OutputDir"] / "uploaded" / "sequences.fasta",
    output:
        clade_report = config["OutputDir"] / "result" / "nextclade" / "clade_report.tsv",
        clade_folder = directory(config["OutputDir"] / "result" / "nextclade" / "others"),
    log:
        config["OutputDir"] / "log" / "nextclade.log"
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    threads: 10
    wrapper: 
        "file:wrappers/nextclade.wrapper.py"