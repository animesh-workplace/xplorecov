rule nextclade:
    input:
        dataset="workflow/resources/nextclade/data",
        sequences=f'{config["OutputDir"]}/uploaded/sequences.fasta',
    output:
        clade_report=f"{config['OutputDir']}/result/nextclade/clade_report.tsv",
        clade_folder=directory(f'{config["OutputDir"]}/result/nextclade/others'),
    log:
        f'{config["OutputDir"]}/log/nextclade.log',
    threads: 10
    conda:
        "~/micromamba-env/.workflow-venv/envs/nibmg_tool"
    wrapper:
        "file:workflow/wrappers/nextclade.wrapper.py"
