import json, requests, os


def parse_version_output(output):
    version_dict = {}
    for line in output.splitlines():
        # Split each line into the name and version
        name, version = line.split(":") if ":" in line else line.split()
        name = name.strip().lower().replace(" ", "_")  # Normalize the key
        version = (
            version.strip() if ":" in line else version.strip()
        )  # Normalize version string
        version_dict[name] = version
    return version_dict


rule update_nextclade_pangolin:
    output:
        resources=directory("workflow/resources/nextclade/data"),
        snpeff_resources=directory("workflow/resources/snpeff"),
    threads: 1
    log:
        pangolin=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_pangolin.log',
        nextclade=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_nextclade.log',
        snpeff=f'{config["UpdateDir"]}/update_log/{config["Date"]}/update_snpeff.log',
    run:
        shell(
            """
                micromamba run -p ".workflow-venv/envs/xplorecov" pangolin --update >> {log.pangolin} 2>&1
                micromamba run -p ".workflow-venv/envs/xplorecov" pangolin --update-data >> {log.pangolin} 2>&1
                micromamba update -y nextclade -r ".workflow-venv/" -n xplorecov >> {log.nextclade} 2>&1
                micromamba run -p ".workflow-venv/envs/xplorecov" nextclade dataset get --name 'sars-cov-2' --output-dir {output.resources}
                mkdir -p {output.snpeff_resources}
                wget -P {output.snpeff_resources} "https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf" >> {log.snpeff} 2>&1
            """
        )
        nextclade_version = parse_version_output(
            shell(
                """micromamba run -p ".workflow-venv/envs/xplorecov" nextclade --version""",
                read=True,
            ).strip()
        )
        pangolin_version = parse_version_output(
            shell(
                """micromamba run -p ".workflow-venv/envs/xplorecov" pangolin --all-version""",
                read=True,
            ).strip()
        )
        tool_versions = {
            **nextclade_version,
            **pangolin_version,
            "BACKEND_WEBSOCKET_UUID": os.getenv("BACKEND_WEBSOCKET_UUID"),
        }

        try:
            response = requests.post(
                "http://10.10.6.80/xplorecov/api/job/update-tool-version/",
                data=json.dumps(tool_versions),
                headers={"Content-Type": "application/json"},
            )
            if response.status_code == 200 or response.status_code == 201:
                print("Versions sent successfully!")
            else:
                print(f"Failed to send versions. Status code: {response.status_code}")
        except Exception as e:
            print(e)
