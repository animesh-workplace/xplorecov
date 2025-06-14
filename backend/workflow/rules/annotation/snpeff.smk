import requests, time
from datetime import datetime


rule snpeff:
    input:
        aligned_fasta=rules.nextclade.output.aligned_fasta,
        problematic_sites="workflow/resources/snpeff/problematic_sites_sarsCov2.vcf",
    output:
        vcf_report=f"{config['OutputDir']}/result/snpeff/output.vcf",
        annotated_vcf_report=f"{config['OutputDir']}/result/snpeff/output.annotated.vcf",
    log:
        f'{config["OutputDir"]}/log/snpeff.log',
    threads: 10
    run:
        run_websocket_message("Deep Analysis Execution", "step6", "start")
        shell(
            """
            time micromamba run -p ".workflow-venv/envs/xplorecov" faToVcf -maskSites={input.problematic_sites} {input.aligned_fasta} {output.vcf_report} > {log} 2>&1
            time sed -i 's/^MN908947/MN908947.3/' {output.vcf_report}
            time micromamba run -p ".workflow-venv/envs/xplorecov" snpEff MN908947.3 {output.vcf_report} > {output.annotated_vcf_report}
            """
        )
        run_websocket_message("Deep Analysis Execution", "step6", "end")
