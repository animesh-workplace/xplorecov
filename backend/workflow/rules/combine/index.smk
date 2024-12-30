rule combine:
    input:
        usher = rules.usher.output.usher_report,
        nextclade = rules.nextclade.output.clade_report,
        pangolin_usher = rules.pangolin_usher.output.lineage_report,
        pangolin_pangolearn = rules.pangolin_pangolearn.output.lineage_report,
    output: 
        report = config["OutputDir"] / "result" / "combined_report.tsv"
    log:
        config["OutputDir"] / "log" / "combined.log"        
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    threads: 1
    wrapper: 
        "file:wrappers/combine.wrapper.py"