rule update_usher:
    output:
        dataset = "resources/usher/data/public-latest.all.masked.pb.gz",
        problem_sites = "resources/usher/data/problematic_sites_sarsCov2.vcf",
    threads: 1
    resources: 
        usher_tree = "http://hgdownload.soe.ucsc.edu/goldenPath/wuhCor1/UShER_SARS-CoV-2/public-latest.all.masked.pb.gz",
        problem_sites = "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf"
    log:
        config["UpdateDir"] / "update_log" / config["Date"] / "update_usher.log"
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    shell:
        """
            curl -fsSL --retry 5 --retry-all-errors --retry-delay 1 {resources.usher_tree} \
            -o {output.dataset} > {log} 2>&1
            curl -fsSL --retry 5 --retry-all-errors --retry-delay 1 {resources.problem_sites} \
            -o {output.problem_sites} >> {log} 2>&1
            conda update -y usher >> {log} 2>&1
        """