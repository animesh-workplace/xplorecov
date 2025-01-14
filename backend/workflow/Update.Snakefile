import shutil
from pathlib import Path

# If the folder exists remove it to allow update of tools
if Path("workflow/resources").is_dir():
    shutil.rmtree("workflow/resources")


rule all:
    input:
        "workflow/resources/nextclade/data",


include: "rules/update/nextclade_pangolin.smk"
