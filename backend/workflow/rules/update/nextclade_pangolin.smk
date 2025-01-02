rule update_nextclade_pangolin:
    output:
        resources=directory("resources/nextclade/data"),
    threads: 1
    log:
        pangolin=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_pangolin.log',
        nextclade=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_nextclade.log',
    shell:
        """
            micromamba run -p ".workflow-venv/envs/xplorecov" pangolin --update >> {log.pangolin} 2>&1
            micromamba run -p ".workflow-venv/envs/xplorecov" pangolin --update-data >> {log.pangolin} 2>&1
            micromamba update -y nextclade -r ".workflow-venv/" -n xplorecov >> {log.nextclade} 2>&1
            micromamba run -p ".workflow-venv/envs/xplorecov" nextclade dataset get --name 'sars-cov-2' --output-dir {output.resources}
        """
