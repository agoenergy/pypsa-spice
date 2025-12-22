# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Generate scenario config file for the assumption of the energy system modelling.

It creates configuration of different scenario settings.
"""

import glob
import os
from typing import Any

from _helpers import configure_logging
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

HEADER = """
# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

# Scenario configuration includes:
# scenario settings,
# model constraints,
# and solver settings.

"""

SCENAIOR_CONFIGS = """
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

if __name__ == "__main__":
    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake  # pylint: disable=ungrouped-imports

        snakemake = mock_snakemake("build_scenario_config")
    # Getting global config params
    configure_logging(snakemake)
    configurations = snakemake.params.config

    # Skeleton porject folder path
    project_folder_path = (
        "data/"
        + configurations["path_configs"]["data_folder_name"]
        + "/"
        + configurations["path_configs"]["project_name"]
    )

    # Skeleton input scenario folder path
    input_folder_path = project_folder_path + "/input/"

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

    data = yaml.load(SCENAIOR_CONFIGS)

    # Initialize certain keys as a dictionary if it is None
    if data["scenario_configs"].get("snapshots") is None:
        data["scenario_configs"]["snapshots"] = {}
    if data["scenario_configs"].get("interest") is None:
        data["scenario_configs"]["interest"] = CommentedMap()

    if data.get("co2_management") is None:
        data["co2_management"] = CommentedMap()

    if data.get("custom_constraints") is None:
        data["custom_constraints"] = CommentedMap()

    start_year = configurations["base_configs"]["years"][0]
    end_year = start_year + 1
    data["scenario_configs"]["snapshots"] = {
        "start": DoubleQuotedScalarString(f"{start_year}-01-01"),
        "end": DoubleQuotedScalarString(f"{end_year}-01-01"),
        "inclusive": DoubleQuotedScalarString("left"),
    }

    # Add comments before or after keys
    SNAPSHOTS_COMMENT_TEXT = (
        '"end" shall be should be "12-31" if base year is a leap year'
    )
    data["scenario_configs"].yaml_set_comment_before_after_key(
        "snapshots", after=SNAPSHOTS_COMMENT_TEXT
    )
    data["scenario_configs"].yaml_set_comment_before_after_key(
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
    data.yaml_set_comment_before_after_key(
        "co2_management", after=CO2_MANAGEMENT_COMMENT_TEXT
    )
    CUSTOM_CONSTRAINTS_COMMENT_TEXT = (
        "please indicate country/region and their specific custom constraints"
    )
    data.yaml_set_comment_before_after_key(
        "custom_constraints", after=CUSTOM_CONSTRAINTS_COMMENT_TEXT
    )

    # Add country specific inputs
    for country in configurations["base_configs"]["regions"].keys():
        # Create a new list for each loop to solve YAML anchors (&id001) issue
        fuels_list = CommentedSeq(
            [
                DoubleQuotedScalarString("Bio"),
                DoubleQuotedScalarString("Bit"),
                DoubleQuotedScalarString("Gas"),
                DoubleQuotedScalarString("Oil"),
            ]
        )
        fuels_list.fa.set_flow_style()  # Force inline style [a, b, c]

        data["scenario_configs"]["interest"][country] = 0.05
        data["co2_management"][country] = {
            "option": DoubleQuotedScalarString("co2_cap"),
            "value": {
                2025: 100,
                2030: 90,
                2035: 80,
                2040: 70,
                2045: 60,
                2050: 50,
            },
        }
        data["custom_constraints"][country] = {
            "energy_independence": {
                "activate": False,
                "pe_conv_fraction": {
                    "Solar": 1,
                    "Wind": 1,
                    "Geothermal": 1,
                    "Water": 1,
                },
                "ei_fraction": {
                    2025: 0.3,
                    2030: 0.4,
                    2035: 0.5,
                    2040: 0.6,
                    2045: 0.7,
                    2050: 0.8,
                },
            },
            "production_constraint_fuels": {
                "activate": False,
                "fuels": fuels_list,
            },
            "reserve_margin": {
                "activate": False,
                "epsilon_load": 0.1,
                "epsilon_vre": 0.1,
                "contingency": 1000,
                "method": DoubleQuotedScalarString("static"),
            },
            "res_generation": {
                "activate": False,
                "math_symbol": DoubleQuotedScalarString("<="),
                "res_generation_share": {
                    2030: 0.25,
                    2035: 0.35,
                    2040: 0.4,
                    2045: 0.45,
                    2050: 0.5,
                },
            },
            "thermal_must_run": {
                "activate": False,
                "min_must_run_ratio": 0.2,
            },
            "capacity_factor_constraint": {
                "activate": False,
                "values": {
                    DoubleQuotedScalarString("SubC"): 0.6,
                    DoubleQuotedScalarString("SupC"): 0.6,
                    DoubleQuotedScalarString("HDAM"): 0.4,
                },
            },
        }

    # create a scenario_config.yaml template file
    yaml_files = glob.glob(f"{input_folder_path}/*.yaml")
    if len(yaml_files) == 0:
        with open(
            f"{input_folder_path}/scenario_config.yaml", "w", encoding="utf-8"
        ) as f:
            f.write(HEADER)
            yaml.dump(
                data,
                f,
                # sort_keys=False, default_flow_style=False, indent=2
            )

        print("YAML file generated successfully.")
