rule usher:
    input:
        nextclade = rules.nextclade.output.clade_report,
        sequence_folder = rules.nextclade.output.clade_folder,
        reference = "resources/nextclade/data/reference.fasta",
        dataset = "resources/usher/data/public-latest.all.masked.pb.gz",
        problem_vcf = "resources/usher/data/problematic_sites_sarsCov2.vcf",
    output: 
        usher_report = config["OutputDir"] / "result" / "usher" / "usher_clades.tsv"
    log:
        config["OutputDir"] / "log" / "usher.log"        
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    threads: 28
    wrapper: 
        "file:wrappers/usher.wrapper.py"