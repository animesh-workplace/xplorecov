rule update_nextclade_pangolin:
    output:
        resources = directory("resources/nextclade/data"),
    threads: 1
    log:
        pangolin = config["UpdateDir"] / "update_log" / config["Date"] / "update_pangolin.log",
        nextclade = config["UpdateDir"] / "update_log" / config["Date"] / "update_nextclade.log",
    conda:
        config["UpdateDir"] / "envs" / "tool.yaml"
    shell:
        """
            pangolin --update > {log.pangolin} 2>&1
            pangolin --update-data >> {log.pangolin} 2>&1
            conda update -y nextclade > {log.nextclade} 2>&1
            nextclade dataset get --name 'sars-cov-2' --output-dir {output.resources}
        """