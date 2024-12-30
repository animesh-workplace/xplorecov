rule update_nextclade_pangolin:
    output:
        resources=directory("resources/nextclade/data"),
    threads: 1
    log:
        pangolin=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_pangolin.log',
        nextclade=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_nextclade.log',
    shell:
        """
            "/mnt/d/Development Projects/PhD Work - NIBMG/xplorecov/backend/micromamba" run -p "/home/brainiac_iii/micromamba-env/.workflow-venv/envs/nibmg_tool" pangolin --update >> {log.pangolin} 2>&1
            "/mnt/d/Development Projects/PhD Work - NIBMG/xplorecov/backend/micromamba" run -p "/home/brainiac_iii/micromamba-env/.workflow-venv/envs/nibmg_tool" pangolin --update-data >> {log.pangolin} 2>&1
            "/mnt/d/Development Projects/PhD Work - NIBMG/xplorecov/backend/micromamba" update -y nextclade -r ~/micromamba-env/.workflow-venv/ -n nibmg_tool >> {log.nextclade} 2>&1
            "/mnt/d/Development Projects/PhD Work - NIBMG/xplorecov/backend/micromamba" run -p "/home/brainiac_iii/micromamba-env/.workflow-venv/envs/nibmg_tool" nextclade dataset get --name 'sars-cov-2' --output-dir {output.resources}
        """
