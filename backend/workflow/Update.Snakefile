import shutil
from pathlib import Path

# If the folder exists remove it to allow update of tools 
if(Path('resources').is_dir()):
    shutil.rmtree("resources")

rule all:
    input: 
        "resources/nextclade/data", 
        "resources/usher/data/public-latest.all.masked.pb.gz", 
        "resources/usher/data/problematic_sites_sarsCov2.vcf",

include: "rules/update/nextclade_pangolin.smk"
include: "rules/update/usher.smk"