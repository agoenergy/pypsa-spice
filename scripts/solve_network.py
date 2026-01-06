# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Solve an energy system optimization problem using the PyPSA library.

The script configures solver settings, and solves the optimization problem for a
specified year and sector. It includes functionality for adding custom constraints to
the model. It supports different solvers and can handle scenarios with multiple years
and countries.
"""

import logging
import os
import pathlib
from functools import partial
from typing import Any

import pandas as pd
import pypsa
from _helpers import configure_logging, load_scenario_config
from custom_constraints import (
    add_energy_independence_constraint,
    add_generation_constraint,
    add_reserve_margin,
    add_storage_constraints,
    capacity_factor_constraint,
    co2_cap_constraint,
    fuel_supply_constraint,
    re_pow_generation_constraint,
    renewable_potential_constraint,
    thermal_must_run_constraint,
)
from dotenv import load_dotenv
from linopy.remote.oetc import OetcCredentials, OetcHandler, OetcSettings

logger = logging.getLogger(__name__)


def extra_functionality_linopt(
    network: pypsa.Network, snapshots: pd.Series, scenario_configs: dict
):
    """Add all custom constraints to the network before solving.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network object containing all components and functions.
    snapshots (pd.Series):
        hours being optimised in the model.
    scenario_configs : dict
        Configuration dictionary containing scenario settings.
    """
    add_storage_constraints(network)
    constraint_added = False
    # pylint: disable=possibly-used-before-assignment
    config = snakemake.config
    year = int(snakemake.wildcards.years)
    base_year = config["base_configs"]["years"][0]
    country = config
    for country in config["base_configs"]["regions"].keys():
        country_emission_settings = scenario_configs["co2_management"][country]
        if country_emission_settings.get("option") == "co2_cap":
            co2_cap_constraint(
                network,
                country=country,
                co2_cap=country_emission_settings["value"][year],
            )
        custom_constraint_set = scenario_configs.get("custom_constraints") or {}
        if country not in custom_constraint_set:
            continue

        country_constraints = scenario_configs["custom_constraints"][country]

        # Capacity factor constraint
        if (
            country_constraints["capacity_factor_constraint"].get("activate", False)
            and year > base_year
        ):
            capacity_factor_constraint(
                network,
                country=country,
                cf_dict=country_constraints["capacity_factor_constraint"]["values"],
            )
            constraint_added = True

        # Power generation constraint
        if country_constraints["generation_constraint"].get("activate", False):
            add_generation_constraint(
                network,
                country=country,
                year=year,
                gen_dict=country_constraints["generation_constraint"]["values"],
            )
            constraint_added = True

        # Energy independence constraint
        if (
            country_constraints["energy_independence"].get("activate", False)
            and year > base_year
        ):
            energy_independence = country_constraints["energy_independence"]
            add_energy_independence_constraint(
                network,
                ei_frac=energy_independence["ei_fraction"][year],
                pe_con_frac=energy_independence["pe_conv_fraction"],
                country=country,
            )
            constraint_added = True

        # Fuel supply/production limit constraint
        if country_constraints["production_constraint_fuels"].get("activate", False):
            fuel_supply_limits = pd.read_csv(snakemake.input.fuel_limits)
            fuel_supply_limits = fuel_supply_limits[
                (fuel_supply_limits.year == year)
                & (fuel_supply_limits.country == country)
            ].set_index("carrier")
            restricted_carriers = country_constraints["production_constraint_fuels"][
                "fuels"
            ]
            fuel_supply_limits = fuel_supply_limits.loc[
                restricted_carriers, "max_supply__mwh_year"
            ].to_dict()
            fuel_supply_constraint(
                network, country=country, supply_limits=fuel_supply_limits
            )
            constraint_added = True

        # Reserve margin constraint
        if country_constraints["reserve_margin"].get("activate", False):
            if year > base_year:
                reserve_params = country_constraints["reserve_margin"]
                add_reserve_margin(
                    network,
                    EP_LOAD=reserve_params["epsilon_load"],
                    EP_VRE=reserve_params["epsilon_vre"],
                    CONT=reserve_params["contingency"],
                    method=reserve_params["method"],
                    country=country,
                )
                constraint_added = True
            elif year == base_year:
                print(
                    "Reserve margin constraint is not added in the first year of the "
                    "simulation"
                )

        # Renewable generation constraint
        if (
            country_constraints["res_generation"].get("activate", False)
            and year > base_year
        ):
            res_generation = country_constraints["res_generation"]
            re_pow_generation_constraint(
                network,
                res_generation_share=res_generation["res_generation_share"][year],
                math_symbol=res_generation["math_symbol"],
                country=country,
            )
            constraint_added = True

        # Thermal must-run constraint
        if country_constraints["thermal_must_run"].get("activate", False):
            thermal_must_run = country_constraints["thermal_must_run"]
            thermal_must_run_constraint(
                network,
                min_must_run_ratio=thermal_must_run["min_must_run_ratio"],
                country=country,
            )
            constraint_added = True

    if not constraint_added:
        print("No custom constraint was added to the model")


def solve_network(
    network: pypsa.Network, year: int, scenario_configs: dict
) -> pypsa.Network:
    """Solve the optimization problem for the given network and year.

    Parameters
    ----------
    network : pypsa.Network
        PyPSA network object containing all components and functions.
    year : int
        Year for which the optimization is being solved.
    scenario_configs : dict
        Configuration dictionary containing solver settings and other parameters.

    Returns
    -------
    pypsa.Network
        The solved PyPSA network object.
    """
    # solver settings
    set_of_options = scenario_configs["solving"]["solver"]["options"]
    solver_options = (
        scenario_configs["solving"]["solver_options"][set_of_options]
        if set_of_options
        else {}
    )
    solver_name = (scenario_configs["solving"]["solver"]["name"]).lower()
    oetc = scenario_configs["solving"]["oetc"]

    if oetc["activate"]:
        # remove activate key from oetc dict
        oetc_setup = {k: v for k, v in oetc.items() if k != "activate"}

        ROOT_FOLDER = pathlib.Path(__file__).parent.parent
        load_dotenv()
        env_path = ROOT_FOLDER / "envs" / ".env"
        load_dotenv(dotenv_path=env_path)
        oetc_setup["authentication_server_url"] = os.getenv("AUTH_SERVER")
        oetc_setup["orchestrator_server_url"] = os.getenv("ORCH_SERVER")
        oetc_setup["credentials"] = OetcCredentials(
            email=os.getenv("OETC_EMAIL"), password=os.getenv("OETC_PASSWORD")
        )
        oetc_setup["solver"] = solver_name
        oetc_setup["solver_options"] = solver_options
        oetc_settings = OetcSettings(**oetc_setup)
        oetc_handler = OetcHandler(oetc_settings)
    else:
        oetc_settings = None
        oetc_handler = None

    print(f"######## 🌶 Solving model with {solver_name.capitalize()}")
    extra_functionality_linopt_config = partial(
        extra_functionality_linopt, scenario_configs=scenario_configs
    )
    # ================================ Solve with Gurobi ===============================
    if solver_name.lower() == "gurobi":
        # solve network
        network.optimize(
            solver_name=solver_name,
            extra_functionality=extra_functionality_linopt_config,
            **({"remote": oetc_handler} if oetc_handler else {}),
            **solver_options,
        )
        if not hasattr(network, "objective"):
            print("####################################################")
            print(f"warning network not solved with barrier method in year {year}")

            solver_options_numerical = {
                "threads": 16,
                "method": -1,
                "BarHomogeneous": 0,
                "crossover": 0,
                "BarConvTol": 1e-5,
                "AggFill": 0,
                "PreDual": 0,
                "FeasibilityTol": 1e-5,
                "GURO_PAR_BARDENSETHRESH": 200,
                "Seed": 123,
                "MarkowitzTol": 0.5,
                "Quad": 1,
            }

            network.optimize(
                solver_name=solver_name,
                solver_options=solver_options_numerical,
                extra_functionality=extra_functionality_linopt_config,
            )
    # ======= Solve with solver with solver-specific options (e.g., CPLEX, HiGHS) ======
    elif solver_name.lower() in ["cplex", "highs"]:
        network.optimize(
            solver_name=solver_name,
            extra_functionality=extra_functionality_linopt_config,
            **({"remote": oetc_handler} if oetc_handler else {}),
            **solver_options,
        )
    # =========== Solve with solver without solver options (e.g., GLPK, CBC) ===========
    else:
        network.optimize(
            solver_name=solver_name,
            keep_references=True,
            extra_functionality=extra_functionality_linopt_config,
        )
    if n.model.status != "warning":
        print("model feasible! 🌶 ")
        return network
    else:
        print("model infeasible compute infeasibilites")
        if solver_name.lower() in ["gurobi", "cplex"]:
            n.model.print_infeasibilities()


if __name__ == "__main__":
    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake  # pylint: disable=ungrouped-imports

        snakemake = mock_snakemake("solve_network", sector="p", years=2023)
    configure_logging(snakemake)
    scenario_configs = load_scenario_config(
        "data/"
        + snakemake.config["path_configs"]["data_folder_name"]
        + "/"
        + snakemake.config["path_configs"]["project_name"]
        + "/input/"
        + snakemake.config["path_configs"]["input_scenario_name"]
    )
    y = int(snakemake.wildcards.years)
    n = pypsa.Network(
        snakemake.input.network,
    )
    if y > snakemake.params.years[0]:  # Only do for year which is not baseyear
        renewable_potential_constraint(
            n, snakemake.input.re_technical_potential, year=y
        )
    n.export_to_netcdf(snakemake.output.pre_solved)
    n = solve_network(n, y, scenario_configs=scenario_configs)
    n.export_to_netcdf(snakemake.output.final_network)
