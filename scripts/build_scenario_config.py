# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Generate scenario config file for the assumption of the energy system modelling.

It creates configuration of different scenario settings.
"""

import glob
import os
from typing import Any

import yaml
from _helpers import configure_logging

HEADER = """
# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

# Scenario configuration includes:
# scenario settings,
# model constraints,
# and solver settings.

"""

SCENAIOR_CONFIGS = {
    "version": "0.1.0",
    "logging": {"level": "INFO", "format": "%(levelname)s:%(name)s:%(message)s"},
    "scenario_configs": {
        "snapshots": {
            "start": "",
            "end": "",  # should be "12-31" if base year is a leap year
            "inclusive": "left",
        },
        "resolution": {
            "method": "nth_hour",  # options: ["nth_hour", "clustered"]
            "number_of_days": 3,  # the number of representative days (used fo "clustered" method) # noqa:E501
            "stepsize": 25,  # used for "nth_hour" method
        },
        "interest": {},
        "remove_threshold": 0.1,  # unit: MW (remove assests with capacity below this threshold) # noqa:E501
    },
    "co2_management": {},
    "custom_constraints": {},
    "solving": {
        "solver": {"name": "highs", "options": "highs-default"},
        "oetc": {
            "activate": False,
            "name": "test-agora-job",
            "cpu_cores": 4,
            "disk_space_gb": 20,
            "delete_worker_on_error": False,
        },
        "solver_options": {
            "default": {},
            "cbc-default": {
                "threads": 8,
                "cuts": 0,
                "maxsol": 1,
                "ratio": 0.1,
                "presolve": 1,
                "time_limit": 3600,
            },
            "gurobi-default": {
                "threads": 8,
                "method": 2,
                "crossover": 0,
                "BarConvTol": 1.0e-5,
                "AggFill": 0,
                "PreDual": 0,
                "GURO_PAR_BARDENSETHRESH": 200,
            },
            "gurobi-numeric-focus": {
                "name": "gurobi",
                "NumericFocus": 3,
                "method": 2,
                "crossover": 0,
                "BarHomogeneous": 1,
                "BarConvTol": 1.0e-5,
                "FeasibilityTol": 1.0e-4,
                "OptimalityTol": 1.0e-4,
                "ObjScale": -0.5,
                "threads": 8,
                "Seed": 123,
            },
            "cplex-default": {
                "threads": 4,
                "lpmethod": 4,
                "solutiontype": 2,
                "barrier_convergetol": 1.0e-5,
                "feasopt_tolerance": 1.0e-6,
            },
            "highs-default": {
                "threads": 4,
                "solver": "ipm",
                "run_crossover": "off",
                "small_matrix_value": 1e-6,
                "large_matrix_value": 1e9,
                "primal_feasibility_tolerance": 1e-5,
                "dual_feasibility_tolerance": 1e-5,
                "ipm_optimality_tolerance": 1e-4,
                "parallel": "on",
                "random_seed": 123,
            },
            "highs-simplex": {
                "solver": "simplex",
                "parallel": "on",
                "primal_feasibility_tolerance": 1e-5,
                "dual_feasibility_tolerance": 1e-5,
                "random_seed": 123,
            },
        },
    },
}


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

    # create a scenario_config.yaml template file
    yaml_files = glob.glob(f"{input_folder_path}/*.yaml")
    if len(yaml_files) == 0:
        with open(f"{input_folder_path}/scenario_config.yaml", "w") as f:
            f.write(HEADER)
            yaml.dump(
                SCENAIOR_CONFIGS, f, sort_keys=False, default_flow_style=False, indent=2
            )

        print("YAML file generated successfully.")
