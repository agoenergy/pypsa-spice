# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Generate scenario config file for the assumption of the energy system modelling.

It creates configuration of different scenario settings.
"""

import glob
import os

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

HEADER = """
# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

# Scenario configuration includes model constraints, scenario and solver settings.

"""

SCENAIOR_DEFAULT_CONFIGS = """
version: 0.1.0

logging:
  level: INFO
  format: "%(levelname)s:%(name)s:%(message)s"

#####Scenario settings#####
scenario_configs:
  snapshots:
  resolution:
    method: "nth_hour" # options: ["nth_hour", "clustered"]
    number_of_days: 3 # the number of representative days (used fo "clustered" method)
    stepsize: 25 # used for "nth_hour" method
  interest:
  remove_threshold: 0.1 # unit: MW (remove assests with capacity below this threshold)


#####Model constraints#####
co2_management:

custom_constraints:

#####Solver settings#####
solving:
  solver:
    name: highs
    options: highs-default
  oetc:
    activate: false
    name: test-agora-job
    cpu_cores: 4
    disk_space_gb: 20
    delete_worker_on_error: false
  solver_options:
    default: {} # used if no solver option is available
    cbc-default:
      threads: 8
      cuts: 0 # no cuts
      maxsol: 1 # only one solution
      ratio: 0.1 # ratio of the solution to the best solution
      presolve: 1 # presolve on
      time_limit: 3600 # time limit in seconds
    gurobi-default:
      threads: 8
      method: 2 # barrier
      crossover: 0
      BarConvTol: 1.e-5
      AggFill: 0
      PreDual: 0
      GURO_PAR_BARDENSETHRESH: 200
    gurobi-numeric-focus:
      name: gurobi
      NumericFocus: 3 # Favour numeric stability over speed
      method: 2 # barrier
      crossover: 0 # do not use crossover
      BarHomogeneous: 1 # Use homogeneous barrier if standard does not converge
      BarConvTol: 1.e-5
      FeasibilityTol: 1.e-4
      OptimalityTol: 1.e-4
      ObjScale: -0.5
      threads: 8
      Seed: 123
    cplex-default:
      threads: 4
      lpmethod: 4 # barrier
      solutiontype: 2 # non basic solution, ie no crossover
      barrier_convergetol: 1.e-5
      feasopt_tolerance: 1.e-6
    highs-default:
      # refer to https://ergo-code.github.io/HiGHS/dev/options/definitions/
      threads: 4
      solver: "ipm"
      run_crossover: "off"
      small_matrix_value: 1e-6
      large_matrix_value: 1e9
      primal_feasibility_tolerance: 1e-5
      dual_feasibility_tolerance: 1e-5
      ipm_optimality_tolerance: 1e-4
      parallel: "on"
      random_seed: 123
    highs-simplex:
      solver: "simplex"
      parallel: "on"
      primal_feasibility_tolerance: 1e-5
      dual_feasibility_tolerance: 1e-5
      random_seed: 123
"""


def add_activate_comments(list_of_targeted_dict: dict):
    """Add comments to all targeted 'activate' values in the given dictionary.

    Parameters
    ----------
    list_of_targeted_dict : dict
        targeted dictionary to add comments to
    """
    for targeted_value in list_of_targeted_dict.keys():
        if targeted_value == "activate":
            list_of_targeted_dict.yaml_add_eol_comment(
                "Change this to true for activating this constraint", targeted_value
            )


def add_please_fill_here_comments_for_dict(
    list_of_targeted_dict: dict, exception_list: list = None
):
    """Add "Please fill here" comments to all targeted values in the given dictionary.

    Parameters
    ----------
    list_of_targeted_dict : dict
        targeted dictionary to add comments to
    exception_list : list, optional
        list of keys to exclude from adding comments, by default []
    """
    if exception_list:
        final_targeted_keys_list = [
            item for item in list_of_targeted_dict.keys() if item not in exception_list
        ]
    else:
        final_targeted_keys_list = list_of_targeted_dict.keys()

    for targeted_value in final_targeted_keys_list:
        list_of_targeted_dict.yaml_add_eol_comment("Please fill here", targeted_value)


def add_country_specific_parameters(input_scenario_data: YAML, config_dict: dict):
    """Add country specific parameters to the scenario configuration data.

    Parameters
    ----------
    scenario_data : YAML
        YAML dictionary to add country specific constraints to
    config_dict : dict
        configurations containing country/region information
    """
    # Add country specific inputs
    regions = config_dict["base_configs"]["regions"]
    years = config_dict["base_configs"]["years"]

    # Iterate through all countries/regions in base_config.yaml
    for country in regions.keys():
        # 1) Setting interest rate block (per country)
        interest_by_country = input_scenario_data["scenario_configs"]["interest"]
        interest_by_country[country] = None
        interest_by_country.yaml_add_eol_comment("Please fill here", country)
        # 2) Setting CO2 management block (per country)
        co2_country_block = CommentedMap(
            {
                "option": DoubleQuotedScalarString("co2_cap"),
                "value": CommentedMap({}),
            }
        )
        input_scenario_data["co2_management"][country] = co2_country_block
        # Populate yearly values with placeholders
        for year in years:
            co2_country_block["value"][year] = None
        # Add "Please fill here" comments for all year entries
        add_please_fill_here_comments_for_dict(co2_country_block["value"])
        # 3) Setting custom constraints block (per country)
        fuels_list = CommentedSeq([])
        fuels_list.fa.set_flow_style()

        custom_constraints_country = CommentedMap(
            {
                # ---- Energy independence constraint ----
                "energy_independence": CommentedMap(
                    {
                        "activate": False,
                        # Primary energy conversion fractions (technology -> fraction)
                        "pe_conv_fraction": CommentedMap(
                            {
                                "Solar": None,
                                "Wind": None,
                                "Geothermal": None,
                                "Water": None,
                            }
                        ),
                        "ei_fraction": CommentedMap({}),
                    }
                ),
                # ---- Production constraint by fuels ----
                "production_constraint_fuels": CommentedMap(
                    {
                        "activate": False,
                        "fuels": fuels_list,
                    }
                ),
                # ---- Reserve margin constraint ----
                "reserve_margin": CommentedMap(
                    {
                        "activate": False,
                        "epsilon_load": None,
                        "epsilon_vre": None,
                        "contingency": None,
                        "method": DoubleQuotedScalarString("static"),
                    }
                ),
                # ---- Renewable generation share constraint ----
                "res_generation": CommentedMap(
                    {
                        "activate": False,
                        "math_symbol": DoubleQuotedScalarString("<="),
                        "res_generation_share": CommentedMap({}),
                    }
                ),
                # ---- Thermal must-run constraint ----
                "thermal_must_run": CommentedMap(
                    {
                        "activate": False,
                        "min_must_run_ratio": None,
                    }
                ),
                # ---- Capacity factor constraint ----
                "capacity_factor_constraint": CommentedMap(
                    {
                        "activate": False,
                        "value": CommentedMap({}),
                    }
                ),
            }
        )

        input_scenario_data["custom_constraints"][country] = custom_constraints_country

        # Fill in year-indexed placeholders inside custom constraints
        energy_independence = custom_constraints_country["energy_independence"]
        res_generation = custom_constraints_country["res_generation"]
        capacity_factor = custom_constraints_country["capacity_factor_constraint"]

        for year in years:
            energy_independence["ei_fraction"][year] = None
            res_generation["res_generation_share"][year] = None
            capacity_factor["value"][year] = None

        # Add comments per constraint (activate toggles, please-fill hints, etc.)
        # Energy independence comments
        add_activate_comments(energy_independence)
        add_please_fill_here_comments_for_dict(energy_independence["pe_conv_fraction"])
        add_please_fill_here_comments_for_dict(energy_independence["ei_fraction"])

        # Production constraint fuels comments (skip activate field)
        prod_fuels = custom_constraints_country["production_constraint_fuels"]
        add_activate_comments(prod_fuels)
        add_please_fill_here_comments_for_dict(prod_fuels, exception_list=["activate"])

        # Reserve margin comments (skip activate + method because method has its own hint)
        reserve_margin = custom_constraints_country["reserve_margin"]
        add_activate_comments(reserve_margin)
        add_please_fill_here_comments_for_dict(
            reserve_margin, exception_list=["activate", "method"]
        )
        reserve_margin.yaml_add_eol_comment(
            "static or dynamic. If static, ignore epsilon_vre", "method"
        )

        # RES generation comments
        add_activate_comments(res_generation)
        add_please_fill_here_comments_for_dict(res_generation["res_generation_share"])

        # Thermal must-run comments
        thermal_must_run = custom_constraints_country["thermal_must_run"]
        add_activate_comments(thermal_must_run)
        add_please_fill_here_comments_for_dict(
            thermal_must_run, exception_list=["activate"]
        )

        # Capacity factor comments
        add_activate_comments(capacity_factor)
        add_please_fill_here_comments_for_dict(capacity_factor["value"])


def build_scenario_config_file(configurations: dict):
    """Combine all the functions of building scenario configuration file.

    Parameters
    ----------
    configurations :
        configurations from the snakemake workflow
    """
    # Skeleton porject folder path
    project_folder_path = (
        "data/"
        + configurations["path_configs"]["data_folder_name"]
        + "/"
        + configurations["path_configs"]["project_name"]
    )

    input_scenario_name = configurations["path_configs"]["input_scenario_name"]

    # Skeleton input scenario folder path
    input_folder_path = project_folder_path + "/input"

    if not os.path.exists(input_folder_path):
        raise FileNotFoundError(
            f"Input folder not found at: {input_folder_path}. "
            "Please ensure that the 'build_skeleton' rule has been executed or that "
            "the project folder exists."
        )

    # Initialize YAML instance
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.default_flow_style = False
    yaml.width = 4096  # Prevent line wrapping

    scenario_data = yaml.load(SCENAIOR_DEFAULT_CONFIGS)

    # Initialize certain keys as a dictionary if it is None
    if scenario_data["scenario_configs"].get("snapshots") is None:
        scenario_data["scenario_configs"]["snapshots"] = {}
    if scenario_data["scenario_configs"].get("interest") is None:
        scenario_data["scenario_configs"]["interest"] = CommentedMap()

    if scenario_data.get("co2_management") is None:
        scenario_data["co2_management"] = CommentedMap()

    if scenario_data.get("custom_constraints") is None:
        scenario_data["custom_constraints"] = CommentedMap()

    start_year = configurations["base_configs"]["years"][0]
    end_year = start_year + 1
    scenario_data["scenario_configs"]["snapshots"] = CommentedMap(
        {
            "start": DoubleQuotedScalarString(f"{start_year}-01-01"),
            "end": DoubleQuotedScalarString(f"{end_year}-01-01"),
            "inclusive": DoubleQuotedScalarString("left"),
        }
    )

    # Add comments before or after keys
    SNAPSHOTS_COMMENT_TEXT = 'shall be should be "12-31" if base year is a leap year'
    scenario_data["scenario_configs"]["snapshots"].yaml_add_eol_comment(
        SNAPSHOTS_COMMENT_TEXT, "end"
    )
    scenario_data["scenario_configs"].yaml_set_comment_before_after_key(
        "interest", after="interest rate in decimals (e.g. 0.05 represents 5%)"
    )
    CO2_MANAGEMENT_COMMENT_TEXT = (
        "please indicate country/region specific CO2 management options\n"
        'options: ["co2_cap", "co2_price"]\n'
        '"co2_cap" is a cap on the total CO2 emissions in the system. '
        "Values are in mtCO2.\n"
        '"co2_price" is a price on CO2 emissions in the system.'
        "Values are in currency/tonne of CO2."
    )
    scenario_data.yaml_set_comment_before_after_key(
        "co2_management", after=CO2_MANAGEMENT_COMMENT_TEXT
    )
    CUSTOM_CONSTRAINTS_COMMENT_TEXT = (
        "please indicate country/region and their specific custom constraints"
    )
    scenario_data.yaml_set_comment_before_after_key(
        "custom_constraints", after=CUSTOM_CONSTRAINTS_COMMENT_TEXT
    )

    add_country_specific_parameters(scenario_data, configurations)

    # Create a scenario_config.yaml template file
    yaml_files = glob.glob(f"{input_folder_path}/{input_scenario_name}/*.yaml")
    if len(yaml_files) == 0:
        with open(
            f"{input_folder_path}/{input_scenario_name}/scenario_config.yaml",
            "w",
            encoding="utf-8",
        ) as f:
            f.write(HEADER)
            yaml.dump(
                scenario_data,
                f,
            )
    else:
        raise (
            FileExistsError(
                "scenario_config.yaml file already exists in the "
                f"{input_folder_path}/{input_scenario_name}. "
                "Please check the folder or adjust path settings in the "
                "based_config.yaml if you want to generate a "
                "new scenario_config.yaml file."
            )
        )
