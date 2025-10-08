# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

from os.path import normpath, exists
from shutil import copyfile, move

if not exists("config.yaml"):
    copyfile("config.default.yaml", "config.yaml")


configfile: "config.yaml"


wildcard_constraints:
    sector=r"(p|p-i|p-t|p-i-t)",
    years=r"\d{4}",


RDIR = (
    config["path_configs"]["results_dir"]
    + config["path_configs"]["project_name"]
    + "/"
    + config["path_configs"]["output_scenario_name"]
    + "/"
)
PDIR = config["path_configs"]["input_dir"] + config["path_configs"]["project_name"]
SDIR = PDIR + "/" + config["path_configs"]["input_scenario_name"]
GPDIR = config["path_configs"]["input_dir"] + "global_csv_templates"

PP_COSTS = PDIR + "/" + "PowerPlant_costs.csv"
TECHNOLOGIES = PDIR + "/" + "Technologies.csv"
STORAGE_COSTS = PDIR + "/" + "Storage_costs.csv"
EV_PARAMETERS = PDIR + "/" + "EV_Parameters.csv"
DMD_PROFILES = PDIR + "/" + "Demand_Profile.csv"
AVAILABILITY = PDIR + "/" + "Availability.csv"
INFLOWS = PDIR + "/" + "Storage-Inflows.csv"
RE_TECH_CAP = PDIR + "/" + "Renewables_technical_potential.csv"


rule add_baseyear:
    input:
        fuel_supplies=SDIR + "/Power/fuel_supplies.csv",
        interconnection=SDIR + "/Power/interconnector.csv",
        elec_buses=SDIR + "/Power/buses.csv",
        elec_loads=SDIR + "/Power/loads.csv",
        elec_power_generators=SDIR + "/Power/power_generators.csv",
        elec_storage_capacity=SDIR + "/Power/storage_capacity.csv",
        elec_storage_energy=SDIR + "/Power/storage_energy.csv",
        elec_power_links=SDIR + "/Power/power_links.csv",
        ind_buses=SDIR + "/Industry/buses.csv",
        tra_buses=SDIR + "/Transport/buses.csv",
        ind_loads=SDIR + "/Industry/loads.csv",
        tra_loads=SDIR + "/Transport/loads.csv",
        ind_generators=SDIR + "/Industry/heat_generators.csv",
        ind_storage_capacity=SDIR + "/Industry/storage_capacity.csv",
        ind_storage_energy=SDIR + "/Industry/storage_energy.csv",
        ind_heat_links=SDIR + "/Industry/heat_links.csv",
        ind_fuel_conversion_links=SDIR + "/Industry/fuel_conversion.csv",
        ind_dac_links=SDIR + "/Industry/direct_air_capture.csv",
        tra_pev_chargers=SDIR + "/Transport/pev_chargers.csv",
        tra_pev_storages=SDIR + "/Transport/pev_storages.csv",
        dmd_profiles=DMD_PROFILES,
        pp_availability=AVAILABILITY,
        stor_inflows=INFLOWS,
        powerplant_costs=PP_COSTS,
        powerplant_type=TECHNOLOGIES,
        storage_costs=STORAGE_COSTS,
        ev_parameters=EV_PARAMETERS,
    output:
        network=RDIR + "pre-solve-brownfield/network_{sector}_{years}.nc",
    params:
        country_region=config["base_configs"]["regions"],
        interest=config["scenario_configs"]["interest"],
        currency=config["base_configs"]["currency"],
        snapshots=config["scenario_configs"]["snapshots"],
        method=config["scenario_configs"]["resolution"]["method"],
        numDays=config["scenario_configs"]["resolution"]["number_of_days"],
        stepsize=config["scenario_configs"]["resolution"]["stepsize"],
        solve_name=config["solving"]["solver"]["name"],
        co2_management=config["co2_management"]
    wildcard_constraints:
        years=config["base_configs"]["years"][0],  #only applies to baseyear
    log:
        "logs/add_electricity_{sector}_{years}.log",
    benchmark:
        "benchmarks/add_electricity_{sector}_{years}"
    threads: 1
    resources:
        mem_mb=5000,
    script:
        "scripts/add_baseyear.py"


def solved_previous_year(wildcards):
    years = config["base_configs"]["years"]
    i = years.index(int(wildcards.years))
    years_p = str(years[i - 1])
    return RDIR + "post-solve/network_{sector}_" + years_p + ".nc"


def previous_year_outputs(wildcards):
    years = config["base_configs"]["years"]
    i = years.index(int(wildcards.years))
    years_p = str(years[i - 1])
    return RDIR + "csvs/{sector}/" + years_p + "/summary.txt"


rule add_brownfield:
    input:
        previous_year_output_folder=previous_year_outputs,
        previous_year_network=solved_previous_year,
        fuel_supplies=SDIR + "/Power/fuel_supplies.csv",
        interconnection=SDIR + "/Power/interconnector.csv",
        elec_buses=SDIR + "/Power/buses.csv",
        elec_loads=SDIR + "/Power/loads.csv",
        elec_power_generators=SDIR + "/Power/power_generators.csv",
        elec_storage_capacity=SDIR + "/Power/storage_capacity.csv",
        elec_storage_energy=SDIR + "/Power/storage_energy.csv",
        elec_power_links=SDIR + "/Power/power_links.csv",
        ind_buses=SDIR + "/Industry/buses.csv",
        tra_buses=SDIR + "/Transport/buses.csv",
        ind_loads=SDIR + "/Industry/loads.csv",
        tra_loads=SDIR + "/Transport/loads.csv",
        ind_generators=SDIR + "/Industry/heat_generators.csv",
        ind_storage_capacity=SDIR + "/Industry/storage_capacity.csv",
        ind_storage_energy=SDIR + "/Industry/storage_energy.csv",
        ind_heat_links=SDIR + "/Industry/heat_links.csv",
        ind_fuel_conversion_links=SDIR + "/Industry/fuel_conversion.csv",
        ind_dac_links=SDIR + "/Industry/direct_air_capture.csv",
        tra_pev_chargers=SDIR + "/Transport/pev_chargers.csv",
        tra_pev_storages=SDIR + "/Transport/pev_storages.csv",
        pow_decom=SDIR + "/Power/decomission_capacity.csv",
        ind_decom=SDIR + "/Industry/decomission_capacity.csv",
        dmd_profiles=DMD_PROFILES,
        pp_availability=AVAILABILITY,
        stor_inflows=INFLOWS,
        powerplant_costs=PP_COSTS,
        powerplant_type=TECHNOLOGIES,
        storage_costs=STORAGE_COSTS,
        ev_parameters=EV_PARAMETERS,
    output:
        brownfield_network=RDIR + "pre-solve-brownfield/network_{sector}_{years}.nc",
    params:
        country_region=config["base_configs"]["regions"],
        years=config["base_configs"]["years"],
        interest=config["scenario_configs"]["interest"],
        currency=config["base_configs"]["currency"],
        remove_threshold=config["scenario_configs"]["remove_threshold"],
        co2_management=config["co2_management"],
    log:
        "logs/add_brownfield_{sector}_{years}.log",
    benchmark:
        "benchmarks/add_brownfield_{sector}_{years}"
    threads: 1
    resources:
        mem_mb=5000,
    script:
        "scripts/add_brownfield.py"


ruleorder: add_baseyear > add_brownfield


rule solve_network:
    input:
        re_technical_potential=RE_TECH_CAP,
        fuel_limits= SDIR + "/Power/fuel_supplies.csv",
        network=RDIR + "pre-solve-brownfield/network_{sector}_{years}.nc",
    output:
        final_network=RDIR + "post-solve/network_{sector}_{years}.nc",
        pre_solved=RDIR + "pre-solve/network_{sector}_{years}.nc",
    params:
        years=config["base_configs"]["years"],
        resolution=config["scenario_configs"]["resolution"],
        currency=config["base_configs"]["currency"]
    threads: config["solving"]["solver"].get("threads", 1)
    script:
        "scripts/solve_network.py"


rule solve_all_networks:
    input:
        expand(
            RDIR + "csvs/{sector}/all_years/combined_summary.txt",
            **config["base_configs"],
        ),
    default_target: True


rule make_summary:
    input:
        network=RDIR + "post-solve/network_{sector}_{years}.nc",
    output:
        summary=RDIR + "csvs/{sector}/{years}/summary.txt",
    params:
        results_dir=config["path_configs"]["results_dir"],
        project_name=config["path_configs"]["project_name"],
        scenario_name=config["path_configs"]["output_scenario_name"],
        resolution=config["scenario_configs"]["resolution"],
        currency=config["base_configs"]["currency"]
    script:
        "scripts/make_summary.py"

rule combine_summaries:
    input:
        networks=expand(
            RDIR + "post-solve/network_{sector}_{years}.nc", **config["base_configs"]
        ),
        summary=expand(RDIR + "csvs/{sector}/{years}/summary.txt", **config["base_configs"]),
    output:
        combined_summaries=RDIR + "csvs/{sector}/all_years/combined_summary.txt",
    script:
        "scripts/combine_summaries.py"
