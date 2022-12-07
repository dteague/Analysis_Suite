configfile: "data/snake_config.yaml"

UL_years = ["2016preVFP", "2016postVFP", "2017", "2018"]
years = ["2016", "2017", "2018"]

fake_dir="workspace/fake_rate"
misId_dir="workspace/charge_misId"
fake_config = config["fake_rate"]

# Stuff related to nonprompt fake rate
rule nonprompt_sideband:
    output:
        fake_dir+"/{area}/mc_scales_{year}.pickle"
    shell:
        "calculate_fr.py -y {wildcards.year} -r sideband -d {wildcards.area} > {output}"

rule nonprompt_measurement:
    input:
        fake_dir+"/{area}/mc_scales_{year}.pickle"
    output:
        fake_dir+"/{area}/fr_{year}.pickle"
    shell:
        "calculate_fr.py -y {wildcards.year} -r measurement -d {wildcards.area} > {output}"

# Stuff related to charge misId
rule misId_measurement:
    output:
        misId_dir+"/{area}/charge_misid_rate_{year}.pickle"
    shell:
        "calculate_mis.py -y {wildcards.year} -r measurement -d {wildcards.area} > {output}"

rule misId_closure:
    input:
        misId_dir+"/{area}/charge_misid_rate_{year}.pickle"
    output:
        misId_dir+"/{area}/fr_{year}.pickle"
    shell:
        "calculate_fr.py -y {wildcards.year} -r measurement -d {wildcards.area} > {output}"

# Final Fake Rates
rule fake_rates:
    input:
        expand(misId_dir+"/{area}/fr_{{year}}.pickle", area=fake_config["misId_directory"]),
        expand(fake_dir+"/{area}/fr_{{year}}.pickle", area=fake_config["nonprompt_directory"])
    output:
        "data/POG/USER/{year}_UL/fake_rates.json.gz"
    shell:
        expand(
            "convert_fakerates.py -y {{wildcards.year}} --nonprompt {fake_area} --charge {misId_area} > {{output}}",
            fake_area=fake_config["nonprompt_directory"],
            misId_area=fake_config["misId_directory"]
        )[0]

rule fake_rates_2016:
    input:
        "data/POG/USER/2016_UL/fake_rates.json.gz"

rule fake_rates_all:
    input:
        expand("data/POG/USER/{years}_UL/fake_rates.json.gz", years=years)

# BTagging efficiencies
rule befficiency:
    output:
        "data/POG/USER/{year}_UL/beff.json.gz"
    shell:
        "beff.py -y {wildcards.year} > {output}"

rule befficiency_all:
    input:
        expand("data/POG/USER/{years}_UL/beff.json.gz", year=UL_years)
