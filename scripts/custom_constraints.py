# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Define custom constraints in a PyPSA network.

These constraints include:
- limits on renewable potential,
- storage behavior,
- CO2 emissions,
- capacity factors,
- thermal generation,
- renewable generation share,
- fuel supply,
- energy independence,
- reserve margins,
- ramp costs,
- demand side management (DSM),
- etc.

DSM is set per country and per load type (``profile_type``) in the
``demand_side_management`` section of ``scenario_config.yaml``. Two types exist:

- ``shed``: one Generator per bus (carrier ``DSM_SHED``, type ``DSMD``). Load can be
  left unserved at a marginal cost. The power is capped at a share of the peak of the
  load on the bus and the weighted energy at a share of the annual demand.
- ``shift``: one lossless StorageUnit per bus (carrier ``DSM_SHIFT``, type ``DSMF``).
  Charging means consuming earlier, discharging means consuming less. The state of
  charge is fixed at half of the energy capacity at every window boundary, so load can
  move earlier or later and the energy balances to zero in each window.

DSM is only applied with hourly resolution (``nth_hour`` with ``stepsize: 1``). With
any other resolution, a warning is logged and DSM is skipped.
"""
import logging
import math
import numbers
import os

import colorama
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import FilePath
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

logger = logging.getLogger(__name__)

RAMP_COST_COLUMNS = ["ramp_up_cost", "ramp_down_cost"]

DSM_SECTION = "demand_side_management"
# Developer-only escape hatch for fast test runs: setting this environment variable to
# 1 applies DSM with any resolution. Results of such runs are not valid.
DSM_TEST_OVERRIDE_ENV = "PYPSA_SPICE_DSM_ANY_RESOLUTION"
DSM_TYPES = {
    "shed": {"carrier": "DSM_SHED", "type": "DSMD"},
    "shift": {"carrier": "DSM_SHIFT", "type": "DSMF"},
}
DSM_CARRIERS = [v["carrier"] for v in DSM_TYPES.values()]
DSM_COMPONENT_TYPES = [v["type"] for v in DSM_TYPES.values()]
DSM_REQUIRED_KEYS = {
    "shed": {"type", "max_peak_share", "max_energy_share", "marginal_cost"},
    "shift": {"type", "max_peak_share", "window_hours", "marginal_cost"},
}
DSM_OPTIONAL_KEYS = {"shed": set(), "shift": {"max_hours"}}
DSM_YEAR_KEYS = ["max_peak_share", "max_energy_share", "marginal_cost", "max_hours"]


def renewable_potential_constraint(
    n: pypsa.Network, technical_potential: FilePath, year: int
):
    """Add constraint on installed capacity of renewables.

    How much installed capacity of renewables per optimization year per bus cannot be
    higher than network total technical potential?

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    technical_potential : FilePath
        Path to a CSV file containing the technical potential data.
    year : int
        The year of optimization for which the constraint is applied.

    Raises
    ------
    ValueError
        If existing capacities exceed the technical potential for any asset.
    """
    technical_potential_df = pd.read_csv(technical_potential)

    for c in ["Generator", "StorageUnit", "Link"]:
        comp_potential_df = technical_potential_df[
            technical_potential_df["class"] == c
        ]  # filter technical technical potential sheet by type
        bus = "bus" if c != "Link" else "bus1"
        if not comp_potential_df.empty:
            print(f"Constraining technical potential for {c}")
            comp_potential_df = comp_potential_df.to_dict("records")
            for row in comp_potential_df:
                new_asset_index = (
                    row["name"] + "_" + str(year)
                )  # construct name of current year asset index
                avail_potential = row[
                    "technical_potential__mw"
                ]  # This is the maximum capacity possible per type per carrier
                if new_asset_index in n.df(c).index:
                    if c == "Link":
                        avail_potential = avail_potential / (
                            n.df(c).loc[new_asset_index, "efficiency"]
                        )  # if component is link,
                        # converting technical power potential to thermal potential
                    existing_assets = [x for x in n.df(c).index if row["name"] in x]
                    existing = (
                        n.df(c)
                        .loc[existing_assets, "p_nom"]
                        .groupby(n.df(c)[bus])
                        .sum()
                        .sum()
                    )  # group all existing capacity per type per carrier
                    remain_tech_potential = (avail_potential - existing).round()
                    if remain_tech_potential < 0:
                        raise ValueError(
                            "existing capacities are already "
                            + "larger than technical potential "
                            + f"for {existing_assets}. Please check the input again"
                        )
                    if (
                        n.df(c).loc[new_asset_index, "p_nom_max"]
                        > remain_tech_potential
                    ):
                        n.df(c).loc[
                            new_asset_index, "p_nom_max"
                        ] = remain_tech_potential
                        # if remaining potential is small than new p_max_pu,
                        # then set the remaining potential as new p_max_pu
        check_assets = (
            n.df(c)
            .loc[np.where(n.df(c)["p_nom_min"] > n.df(c)["p_nom_max"], True, False)]
            .index
        )
        if len(check_assets) > 0:
            print(
                colorama.Fore.RED
                + f"Warning!!! These generators {check_assets} "
                + "have p_nom_min larger than p_nom_max"
            )


def add_storage_constraints(n: pypsa.Network):
    """Add constraint to ensure that charger = discharger for storage links.

    Formulation: 1 * charger_size - efficiency * discharger_size = 0

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    """
    if not n.links.p_nom_extendable.any():
        return

    discharger_bool = n.links.type.str.contains("DISCHARGE")
    charger_bool = n.links.type.str.contains("STORE")

    dischargers_ext = n.links[discharger_bool].query("p_nom_extendable").index
    chargers_ext = n.links[charger_bool].query("p_nom_extendable").index

    eff = n.links.efficiency[dischargers_ext].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers_ext]
        - n.model["Link-p_nom"].loc[dischargers_ext] * eff
    )

    n.model.add_constraints(lhs, "==", 0, name="Link-charger_ratio")


def co2_cap_constraint(n: pypsa.Network, country: str, co2_cap: float):
    """Add constraint on CO2 emissions for a specific country.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the CO2 cap is applied.
    co2_cap : float
        Limit on CO2 emissions in million tonnes (MtCO2).
    """
    bus_carrier = n.stores.bus.map(n.buses.carrier)
    stores = n.stores[
        (n.stores.bus == f"{country}_ATMP")
        & (bus_carrier.isin(["CO2"]))
        & (~n.stores.e_cyclic)
        & (n.stores.country == country)
    ].index
    rhs = co2_cap * 1e6
    lhs = n.model["Store-e"].loc[n.snapshots[-1], stores]
    print(f"Adding CO2 Cap of {co2_cap}mtCO2 for {country}")
    n.model.add_constraints(lhs, "<=", rhs, name=f"co2_cap_{country}")


def capacity_factor_constraint(n: pypsa.Network, country: str, cf_dict: dict):
    """Add constraint on capacity factor for specific generator types.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    cf_dict : dict
        Dictionary with generator types as keys and their corresponding capacity
        factor limits as values.
    """
    lhs = 0
    for gen_type in cf_dict:
        print(
            "....adding maximum capacity factor constraint: "
            f"{gen_type} - {(cf_dict[gen_type]*100):.2f}%"
        )
        for c in ["Generator", "StorageUnit", "Link"]:
            df = n.df(c)
            p_gen = "p" if c != "StorageUnit" else "p_dispatch"
            bus_name = "bus1" if c == "Link" else "bus"
            # Get p_max_pu
            p_max_pu = get_as_dense(n, c, "p_max_pu")
            if not df.empty:
                gen = df[
                    (df.type == gen_type)
                    & (df[bus_name].str.contains("HVELEC"))
                    & (df.country == country)
                ].index

                if not gen.empty:
                    ext_i, fix_i = __get_asset_indices(df=df.loc[gen])
                    # Get var and cf
                    cf = cf_dict.get(gen_type)
                    dispatch = n.model[f"{c}-{p_gen}"]
                    if not ext_i.empty:
                        capacity_variable = (
                            n.model[f"{c}-p_nom"].sel(name=ext_i).rename({"name": c})
                        )
                        # LHS
                        if c == "Link":
                            eff = xr.DataArray(df.loc[gen, "efficiency"])
                            lhs = (
                                dispatch.loc[:, ext_i].mul(eff.loc[ext_i]).sum().sum()
                                - capacity_variable.mul(p_max_pu.loc[:, ext_i])
                                .mul(eff.loc[ext_i])
                                .mul(cf)
                                .sum()
                                .sum()
                            )
                        else:
                            lhs = (
                                dispatch.loc[:, ext_i].sum().sum()
                                - capacity_variable.mul(p_max_pu.loc[:, ext_i])
                                .mul(cf)
                                .sum()
                                .sum()
                            )
                        n.model.add_constraints(
                            lhs,
                            "<=",
                            0,
                            name=f"updated_capacity_factor_{gen_type}_{country}"
                            + "_constraint_ext",
                        )
                    if not fix_i.empty:
                        # Get var and cf
                        capacity_fixed = df["p_nom"].loc[fix_i]
                        if c == "Link":
                            eff = xr.DataArray(df.loc[gen, "efficiency"])
                            capacity_fixed = df["p_nom"].loc[fix_i]
                            lhs = dispatch.loc[:, fix_i].mul(eff.loc[fix_i]).sum().sum()
                            rhs = (
                                (
                                    p_max_pu.loc[:, fix_i]
                                    .mul(capacity_fixed)
                                    .mul(eff.loc[fix_i])
                                    .mul(cf)
                                )
                                .sum()
                                .sum()
                            )
                        else:
                            lhs = dispatch.loc[:, fix_i].sum().sum()
                            rhs = (
                                p_max_pu.loc[:, fix_i]
                                .mul(capacity_fixed)
                                .mul(cf)
                                .sum()
                                .sum()
                            )
                        n.model.add_constraints(
                            lhs,
                            "<=",
                            rhs,
                            name=f"updated_capacity_factor_{gen_type}_{country}"
                            + "_constraint_fix",
                        )


def thermal_must_run_constraint(
    n: pypsa.Network, min_must_run_ratio: float, country: str
):
    """Add constraint on minimum thermal generation as a fraction of total load.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    min_must_run_ratio : float
        Share of thermal generation / total load per snapshot
    country : str
        Country for which the constraint is applied.
    """
    # RHS (Right Hand Side):
    # min_must_run_ratio * total load per snapshot
    rhs = (
        min_must_run_ratio
        * n.loads_t.p_set.T.groupby(n.loads.carrier).sum().T["Electricity"]
    )
    # LHS (Left Hand Side):
    # sum of all thermal generation dispatch variables per snapshot
    lhs = 0
    for c in ["Generator", "StorageUnit", "Link"]:
        df = n.df(c)
        p_gen = "p" if c != "StorageUnit" else "p_dispatch"
        # thermal assets index --------------------------------------------------------
        thermal_asset_index = df[
            (df.country == country)
            & (
                df.type.isin(
                    [
                        "SubC",
                        "SupC",
                        "OILT",
                        "OCGT",
                        "BIOT",
                        "CCGT",
                        "OCHT",
                        "WSTT",
                        "NUCL",
                    ]
                )
            )
        ].index
        if not thermal_asset_index.empty:
            # get dispatch variables
            gen_var = n.model[f"{c}-{p_gen}"].loc[:, thermal_asset_index]
            eff = xr.DataArray(df.loc[thermal_asset_index, "efficiency"])
            if c == "Link":
                gen_var = gen_var.mul(eff)
            # get dispatch variable sum
            lhs += gen_var.sum("name")
    # define constraint ------------------------------------------------------
    print(
        f"....add minimum thermal must run as {round(min_must_run_ratio*100)}% "
        f"of power load per snapshot for {country}"
    )
    n.model.add_constraints(lhs, ">=", rhs, name=f"thermal_must_run_limit_{country}")


def re_pow_generation_constraint(
    n: pypsa.Network, res_generation_share: float, country: str, math_symbol: str = ">="
):
    """Add renewable generation lower limit as fraction of total electrical load.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    res_generation_share : float
        Share of renewable / total power load
    country : str
        Country for which the constraint is applied.
    math_symbol : str, optional
        Math symbol of the constraint formula, by default ">="
    """
    # LHS (Left Hand Side): objectives of the renewable generation
    lhs = 0
    for c in ["Generator", "StorageUnit", "Link"]:
        df = n.df(c)
        if not df.empty:
            p_gen = "p" if c != "StorageUnit" else "p_dispatch"
            bus_name = "bus1" if c == "Link" else "bus"
            res_gen = df[
                (df.country == country)
                & (
                    df.type.isin(
                        [
                            "PHOT",
                            "FLOT",
                            "CSP",
                            "RTPV",
                            "GEOT",
                            "GEOX",
                            "WTON",
                            "WTOF",
                            "HDAM",
                            "HROR",
                            "HPHS",
                            "BATS",
                            "BIOT",
                            "WSTT",
                        ]
                    )
                )
                & (df[bus_name].str.contains("HVELEC"))
            ].index
            if not res_gen.empty:
                # Get var and weights
                gen_var = n.model[f"{c}-{p_gen}"].loc[:, res_gen]
                weight_gen = (
                    n.snapshot_weightings["objective"].reindex(df.index).to_xarray()
                )
                if c == "Link":
                    eff = xr.DataArray(df.loc[res_gen, "efficiency"])
                    lhs += (weight_gen * gen_var.mul(eff)).sum().sum()
                else:
                    lhs += (weight_gen * gen_var).sum().sum()
    # RHS (Right Hand Side): targeted renewable generation share
    rhs = res_generation_share * (
        n.loads_t.p_set.multiply(n.snapshot_weightings.objective, axis=0)
        .T.groupby([n.loads.carrier, n.loads.country])
        .sum()
        .loc["Electricity", country]
        .sum()
    )
    print(
        f"....add minimum renewable generation as: {round(res_generation_share*100)}% "
        f"of total power load in {country}"
    )

    n.model.add_constraints(
        lhs,
        math_symbol,
        rhs,
        name=f"RE share of {round(res_generation_share*100)}% in {country}",
    )


def fuel_supply_constraint(n: pypsa.Network, country: str, supply_limits: dict):
    """Add constraint on fuel supply limits for specific carriers in a country.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    supply_limits : dict
        Dictionary with fuel carriers as keys and their corresponding supply limits.
    """
    for carrier in supply_limits:
        supply_gen = n.generators[
            (n.generators.country == country) & (n.generators.carrier == carrier)
        ].index
        rhs = supply_limits.get(carrier)
        lhs = (
            n.model["Generator-p"]
            .loc[:, supply_gen]
            .mul(n.snapshot_weightings.generators)
            .sum("name")
            .sum()
        )
        print(f"....add a supply limit for {carrier} in {country}")
        n.model.add_constraints(lhs <= rhs, name=f"supply_limit_{carrier}_{country}")


def add_energy_independence_constraint(
    n: pypsa.Network, country: str, ei_frac: float, pe_con_frac: dict
):
    """Add constraint to ensure a certain level of energy independence.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    ei_frac : float
        Fraction of energy independence to be enforced.
    pe_con_frac : dict
        Dictionary with energy carriers as keys and their corresponding primary energy
        conversion factors.
    """
    print(f"....adding energy independence constraint for {country}: {ei_frac}")

    # Theoretical fuel supply plants (e.g., Oil, Gas, Biomass)
    # Get the indices of generators that have "SUPPLY" in their type
    fs_gen = n.generators[
        (n.generators.type.str.contains("SUPPLY")) & (n.generators.country == country)
    ].index
    # Filter out the import fuel supply generators (those containing "-IMP")
    fs_imp_gen = fs_gen[fs_gen.str.contains("-IMP")]  # Imported fuel generators
    # Get the local generators by excluding the import fuel supply generators
    fs_loc_gen = [ele for ele in fs_gen if ele not in fs_imp_gen]  # local generators

    # Get var and weights
    gen_var = n.model["Generator-p"]
    weight_gen = (
        n.snapshot_weightings["objective"].reindex(n.generators.index).to_xarray()
    )
    # Primary energy import (for imported fossil generators)
    gen_var_imp = gen_var.loc[:, fs_imp_gen]
    pe_imp = (weight_gen * gen_var_imp).sum().sum()

    # Primary energy local (for local fossil generators)
    gen_var_loc = gen_var.loc[:, fs_loc_gen]
    pe_loc = (weight_gen * gen_var_loc).sum().sum()

    # Fossil-free Generators (e.g., Solar, Wind, Hydro, Geothermal)
    for c in ["Generator", "StorageUnit"]:
        df = n.df(c)
        p_gen = "p" if c != "StorageUnit" else "p_dispatch"
        # Get var and weights
        gen_var = n.model[f"{c}-{p_gen}"]
        weight_gen = n.snapshot_weightings["objective"].reindex(df.index).to_xarray()
        for res in ["Solar", "Wind", "Geothermal", "Water"]:
            # Local generators
            res_gen_loc = df[
                (df.country == country)
                & (df.carrier == res)
                & ~(df.index.str.contains("-IMP"))
            ].index
            if not res_gen_loc.empty:
                pe_loc += (
                    (weight_gen * gen_var.loc[:, res_gen_loc] * pe_con_frac[res])
                    .sum()
                    .sum()
                )
            # Imported generators
            res_gen_imp = df[
                (df.country == country)
                & (df.carrier == res)
                & (df.index.str.contains("-IMP"))
            ].index
            if not res_gen_imp.empty:
                pe_imp += (
                    (weight_gen * gen_var.loc[:, res_gen_imp] * pe_con_frac[res])
                    .sum()
                    .sum()
                )
            # Add PE conversion faction to pypsa components
            combined_gen = res_gen_loc.union(res_gen_imp)
            df.loc[combined_gen, "pe_factor"] = pe_con_frac.get(res)

    # Final equation
    lhs = (1 - ei_frac) * pe_loc - ei_frac * pe_imp
    rhs = 0

    n.model.add_constraints(
        lhs,
        ">=",
        rhs,
        f"Energy Independence Constraint of >= {ei_frac} for {country}",
    )


def add_maximum_power_generation_constraint(
    n: pypsa.Network, country: str, year: int, gen_dict: dict
):
    """Add constraint on maximum generation per techology, country and year.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    year: int
        The year of optimization for which the constraint is applied.
    gen_dict : dictionary
        Maximum generation (TWh) of each technology per each modelling year.
    """
    weight_da = xr.DataArray(
        n.snapshot_weightings["objective"],
        dims="snapshot",
        coords={"snapshot": n.snapshots},
    )

    for gen_type in gen_dict.keys():
        if year in gen_dict[gen_type].keys():
            lhs = 0
            for c in ["Generator", "StorageUnit", "Link"]:
                df = n.df(c)
                p_gen = "p" if c != "StorageUnit" else "p_dispatch"
                bus_name = "bus1" if c == "Link" else "bus"

                if not df.empty:

                    gen_index = df[
                        (df.country == country)
                        & (df.type == gen_type)
                        & (df[bus_name].str.contains("HVELEC"))
                    ].index

                    if gen_index.empty:
                        continue

                    # Get var and weights
                    gen_var = n.model[f"{c}-{p_gen}"].loc[:, gen_index]

                    if c == "Link":
                        eff = xr.DataArray(
                            df.loc[gen_index, "efficiency"],
                            dims="name",
                            coords={"name": gen_index},
                        )
                        lhs += (gen_var * eff * weight_da).sum("snapshot").sum("name")
                    else:
                        lhs += (gen_var * weight_da).sum("snapshot").sum("name")

            # rhs = maximum generation from gen_dict in converted from TWh to MWh
            total_gen = gen_dict[gen_type][year]
            rhs = total_gen * 1e6

            print(
                "....adding maximum power generation constraint: "
                f"{gen_type} - {total_gen} TWh "
                f"in {country} in {year}"
            )

            n.model.add_constraints(
                lhs, "<=", rhs, name=f"max_generation_{gen_type}_{country}_{year}"
            )


def add_ramp_costs(n: pypsa.Network, country: str, technologies: list):
    """Add ramp-up and ramp-down costs of selected technologies to the objective.

    Ramping of each asset is split into non-negative ramp-up and ramp-down
    variables: ramp_up(t) - ramp_down(t) = P(t) - P(t-1), for t = 1, ..., T-1.
    P is the power output, i.e. p for Generators and p * efficiency (output at bus1)
    for Links. Ramp costs (currency/MW) are the ramp_up_cost and ramp_down_cost
    attributes of the assets (from technologies.csv), weighted by the objective
    snapshot weightings.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the ramp costs will be applied.
    country : str
        Country for which the ramp costs are applied.
    technologies : list
        Technology types (e.g. ["SubC", "SupC", "CCGT"]) to which ramp costs apply.
    """
    # ramping snapshots: t = 1, ..., T-1
    snapshots = n.snapshots[1:]
    weightings = xr.DataArray(
        n.snapshot_weightings["objective"].loc[snapshots],
        dims="snapshot",
        coords={"snapshot": snapshots},
    )

    for c in ["Generator", "Link"]:
        df = n.df(c)
        missing_cols = set(RAMP_COST_COLUMNS).difference(df.columns)
        if missing_cols:
            raise KeyError(
                f"{c}s have no {sorted(missing_cols)} attributes. Please add them to "
                "technologies.csv and rebuild the network to use ramp costs."
            )
        df = df[
            (df.country == country)
            & (df.type.isin(technologies))
            & ((df.p_nom > 0) | df.p_nom_extendable)
        ]
        costs = df[RAMP_COST_COLUMNS].fillna(0)
        costs = costs[(costs > 0).any(axis=1)].rename_axis("name")
        if costs.empty:
            continue
        assets = costs.index

        # Power output: p for Generators and p * efficiency for Links
        p = n.model[f"{c}-p"].sel(name=assets)
        if c == "Link":
            p = p * xr.DataArray(df.loc[assets, "efficiency"])

        ramp_up = n.model.add_variables(
            lower=0,
            coords=[snapshots, assets],
            name=f"{c}-ramp_up-{country}",
        )
        ramp_down = n.model.add_variables(
            lower=0,
            coords=[snapshots, assets],
            name=f"{c}-ramp_down-{country}",
        )
        # diff keeps the first snapshot, so restrict it to t = 1, ..., T-1
        n.model.add_constraints(
            ramp_up - ramp_down == p.diff("snapshot").sel(snapshot=snapshots),
            name=f"{c}-ramp_definition-{country}",
        )

        n.model.objective += (
            ramp_up * xr.DataArray(costs["ramp_up_cost"]) * weightings
        ).sum() + (ramp_down * xr.DataArray(costs["ramp_down_cost"]) * weightings).sum()
        print(
            f"....adding ramp costs for {len(assets)} {c}s "
            f"of types {sorted(df.loc[assets, 'type'].unique())} in {country}"
        )


def add_reserve_margin(
    n: pypsa.Network,
    EP_LOAD: float,
    EP_VRE: float,
    CONT: int,
    country: str,
    method: str = "dynamic",
):
    """Add constraint on the reserve margin.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    EP_LOAD : float
        Fraction of the total load considered as reserve.
    EP_VRE : float
        Contribution of variable renewable energy to the reserve.
    CONT : int
        Fixed contingency in MW.
    country : str
        Country for which the constraint is applied.
    method : str, optional
        "static" means no VRE and "dynamic" includes VRE, by default "dynamic"

    Raises
    ------
    ValueError
        If the method is not "static" or "dynamic".
    """
    if method not in ["static", "dynamic"]:
        raise ValueError(
            "Please choose either 'static' or 'dynamic' method for reserve constraint"
        )

    if method == "static":
        _add_reserve_margin_static(n, EP_LOAD, country)
    elif method == "dynamic":

        for c in ["Generator", "StorageUnit", "Link"]:
            df = n.df(c)
            # select index of the components
            reserve_assets = df.query("r_rating > 0").index
            # define linopy variables for reserve
            try:
                n.model.add_variables(
                    0, np.inf, coords=[n.snapshots, reserve_assets], name=f"{c}-r"
                )
            except ValueError:
                # Skip if component reserve variables are already added
                continue

        # add constraints
        _add_reserve_margin_dynamic(n, EP_LOAD, EP_VRE, CONT, country)
        _update_capacity_constraint_non_link(n, country)
        _update_capacity_constraint_links(n, country)
        _update_storage_reserve_constraint(n, country)


def _add_reserve_margin_static(n: pypsa.Network, ep_load: float, country: str):
    """Add operational static reserve margin as fraction of peak load.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    ep_load : float
        Fraction of the total load considered as reserve.
    country : str
        Country for which the constraint is applied.
    """
    print(f"....adding static reserve margin condition for {country}")
    fix_cap = 0
    lhs = 0
    for c in ["Generator", "StorageUnit", "Link"]:
        fix_i = (
            n.df(c)
            .query("not p_nom_extendable & r_rating > 0")
            .loc[lambda df: df["country"] == country]
            .index
        )
        ext_i = (
            n.df(c)
            .query("p_nom_extendable & r_rating > 0")
            .loc[lambda df: df["country"] == country]
            .index
        )
        r_rating = xr.DataArray(
            n.df(c).loc[ext_i, "r_rating"].rename({f"{c}": f"{c}-ext"})
        )
        if not fix_i.empty:
            if c == "Link":
                fix_cap += (
                    n.df(c)
                    .loc[fix_i, ["r_rating", "p_nom", "efficiency"]]
                    .prod(axis=1)
                    .sum()
                )
            else:
                fix_cap += (
                    n.df(c).loc[fix_i, "r_rating"].mul(n.df(c).loc[fix_i, "p_nom"])
                ).sum()

        if not ext_i.empty:
            eff = None
            if c == "Link":
                eff = xr.DataArray(
                    n.df(c).loc[ext_i, "efficiency"],
                    coords={"name": ext_i},
                    dims=["name"],
                )
            lhs += (
                n.model.variables[f"{c}-p_nom"].sel(name=ext_i)
                * r_rating
                * (eff if eff is not None else 1)
            ).sum("name")

    # rhs = (1+ep_load) * peak load - fix_cap
    load_id = n.loads[
        (n.loads.carrier == "Electricity") & (n.loads.country == country)
    ].index
    peak_load = get_as_dense(n, "Load", "p_set").loc[:, load_id].sum(axis=1).max()
    rhs = (1 + ep_load) * peak_load - fix_cap

    n.model.add_constraints(lhs, ">=", rhs, name=f"reserve_margin_static_{country}")


def _add_reserve_margin_dynamic(
    n: pypsa.Network,
    ep_load: float,
    ep_vre: float,
    cont: int,
    country: str,
):
    """Add operational dynamic reserve margin condition.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    ep_load : float
        Fraction of the total load considered as reserve.
    ep_vre : float
        Contribution of variable renewable energy to the reserve.
    cont : int
        Fixed contingency in MW.
    country : str
        Country for which the constraint is applied.
    """
    print(f"....adding dynamic reserve margin condition for {country}")

    # select index of generators and links
    # ------ LHS----------------------------------------------
    # First get all variable parts of the equation to the left-hand-side (lhs)
    # lhs = reserve * dr factors - ep_vres * sum_g (CF*G)
    # reserve margin of asset g at time t
    lhs = 0
    for c in ["Generator", "StorageUnit", "Link"]:
        df = n.df(c)
        # select index of the components
        reserve_asset_indexes = (
            df.query("r_rating > 0").loc[lambda df: df["country"] == country].index
        )
        assets_reserve = n.model[f"{c}-r"].loc[:, reserve_asset_indexes]
        r_rating = xr.DataArray(df.loc[reserve_asset_indexes, "r_rating"])
        assets_sumed_reserve = (assets_reserve * r_rating).sum("name")
        # add total reserve of the component to the lhs
        lhs += assets_sumed_reserve
    # subtract share of extendable renewable capacities
    ext_i = (
        n.generators[n.generators.country == country].query("p_nom_extendable").index
    )
    vres_i = n.generators[n.generators.country == country].index.intersection(
        n.generators_t.p_max_pu.columns
    )
    if not ext_i.empty and not vres_i.empty:
        # intersection of VRES and extendable generators
        vres_ext = vres_i.intersection(ext_i)

        if not vres_ext.empty:
            # capacity factor CF: DataFrame (snapshots x generators)
            # -> DataArray (snapshot, name)
            cf_df = n.generators_t.p_max_pu[vres_ext]  # pandas DataFrame
            capacity_factor_ext = xr.DataArray(
                cf_df.values,
                coords={"snapshot": cf_df.index, "name": cf_df.columns},
                dims=("snapshot", "name"),
            )

            # capacities of VRES: linopy/xarray variable with dim "name"
            p_nom_vres = n.model["Generator-p_nom"].sel(name=vres_ext)

            # broadcast over snapshot, sum over generators
            lhs += (p_nom_vres * (-ep_vre * capacity_factor_ext)).sum("name")

    # minus the electricity going to storage
    store_links = n.links[
        (n.links.country == country)
        & (n.links.type.str.contains("STORE"))
        & (n.links.carrier == "Electricity")
    ].index
    storageunits_id = n.storage_units[
        (n.storage_units.country == country)
        & (n.storage_units.carrier == "Electricity")
    ].index
    if not store_links.empty:
        storage_load = n.model["Link-p"].loc[:, store_links]
        store_eff = xr.DataArray(n.links.loc[store_links, "efficiency"])
        lhs += (-storage_load.mul(store_eff)).sum("name")
    if not storageunits_id.empty:
        su_store = n.model["StorageUnit-p_store"].loc[:, storageunits_id]
        lhs += (-su_store).sum("name")

    # ------------- RHS -------------------------------------------------
    # ------------- Generators
    # Right-hand-side (rhs) all terms which do not depend on variables

    # Total electricity demand per t summed over all buses
    load_id = n.loads[
        (n.loads.carrier == "Electricity") & (n.loads.country == country)
    ].index
    demand = get_as_dense(n, "Load", "p_set").loc[:, load_id].sum(axis=1)

    # VRES potential of non extendable generators
    capacity_factor_fix = n.generators_t.p_max_pu[vres_i.difference(ext_i)]
    renewable_capacity_fix = n.generators.p_nom[vres_i.difference(ext_i)]
    potential = (capacity_factor_fix * renewable_capacity_fix).sum(axis=1)
    rhs = ep_load * demand + ep_vre * potential + cont
    n.model.add_constraints(lhs, ">=", rhs, name=f"dynamic_reserve_margin_{country}")


def _update_capacity_constraint_non_link(n: pypsa.Network, country: str):
    """Update capacity constraints for non-link components.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    """
    for c in ["Generator", "StorageUnit"]:
        print(f"....updating capacity constraint for {c.lower()}s in {country}")
        df = n.df(c)
        p_gen = "p" if c != "StorageUnit" else "p_dispatch"
        # select index of the components
        reserve_asset_indices, ext_i, fix_i = __get_reserve_asset_indices(df, country)
        dispatch = n.model[f"{c}-{p_gen}"].loc[:, reserve_asset_indices]
        reserve = n.model[f"{c}-r"].loc[:, reserve_asset_indices]
        # Get p_max_pu
        p_max_pu = get_as_dense(n, c, "p_max_pu")

        if not ext_i.empty:
            capacity_variable = n.model[f"{c}-p_nom"].sel(name=ext_i)

            # align everything on axis=1 (asset name)
            lhs = (
                dispatch.loc[:, ext_i]
                + reserve.loc[:, ext_i]
                - (p_max_pu.loc[:, ext_i] * capacity_variable)
            )
            # variable assets constraint
            n.model.add_constraints(
                lhs,
                "<=",
                0,
                name=f"updated_capacity_constraint_{c.lower()}s_ext_{country}",
            )
        if not fix_i.empty:
            # fix assets
            capacity_fixed = df["p_nom"].loc[fix_i]
            lhs = dispatch.loc[:, fix_i] + reserve.loc[:, fix_i]
            rhs = p_max_pu.loc[:, fix_i].mul(capacity_fixed)
            # fix assets constraint
            n.model.add_constraints(
                lhs,
                "<=",
                rhs,
                name=f"updated_capacity_constraint_{c.lower()}s_fix_{country}",
            )


def _update_capacity_constraint_links(n: pypsa.Network, country: str):
    """Update capacity constraints for link components.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    """
    print(f"....updating capacity constraint for links in {country}")
    df = n.df("Link")
    # select index of the components
    reserve_asset_indices, ext_i, fix_i = __get_reserve_asset_indices(df, country)
    dispatch = n.model["Link-p"].loc[:, reserve_asset_indices]
    reserve = n.model["Link-r"].loc[:, reserve_asset_indices]
    # Get r_rating and efficiencies of links
    eff = xr.DataArray(
        df.loc[reserve_asset_indices, "efficiency"],
        coords={"name": reserve_asset_indices},
        dims=["name"],
    )
    # Get p_max_pu
    p_max_pu = get_as_dense(n, "Link", "p_max_pu")
    if not ext_i.empty:
        # variable assets
        capacity_variable = n.model["Link-p_nom"].sel(name=ext_i)
        lhs = (
            dispatch.loc[:, ext_i].mul(eff.loc[ext_i])
            + reserve.loc[:, ext_i]
            - capacity_variable.mul(p_max_pu.loc[:, ext_i]).mul(eff.loc[ext_i])
        )
        # variable assets constraint
        n.model.add_constraints(
            lhs, "<=", 0, name=f"updated_capacity_constraint_links_ext_{country}"
        )
    if not fix_i.empty:  # fix assets
        capacity_fixed = df["p_nom"].loc[fix_i]
        lhs = dispatch.loc[:, fix_i].mul(eff.loc[fix_i]) + reserve.loc[:, fix_i]
        rhs = p_max_pu.loc[:, fix_i].mul(capacity_fixed).mul(eff.loc[fix_i])
        # fix assets constraint
        n.model.add_constraints(
            lhs, "<=", rhs, name=f"updated_capacity_constraint_links_fix_{country}"
        )


def _update_storage_reserve_constraint(n: pypsa.Network, country: str):
    """Update state of charge constraint for storage assets.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the constraint is applied.
    """
    print(f"....updating state of charge constraint for storage assets in {country}")
    # constraining soc for StorageUnit components
    df = n.df("StorageUnit")
    reserve_asset_indexes = (
        df.query("r_rating > 0").loc[lambda df: df["country"] == country].index
    )
    if not reserve_asset_indexes.empty:
        assets_reserve = n.model["StorageUnit-r"].loc[:, reserve_asset_indexes]
        assets_soc = n.model["StorageUnit-state_of_charge"].loc[
            :, reserve_asset_indexes
        ]
        # lhs = (reserve of asset g at time t * r_rating)
        # - (0.25 * state of charge of asset g at time t)
        lhs = assets_reserve - assets_soc.mul(0.25)
        n.model.add_constraints(
            lhs, "<=", 0, name=f"updated_soc_constraint_storage_units_{country}"
        )

    # constraining soc for links plus store buses
    df = n.df("Link")
    reserve_asset_indexes = (
        df.query('r_rating > 0 & carrier == "Electricity"')
        .loc[lambda df: df["country"] == country]
        .index
    )
    if not reserve_asset_indexes.empty:
        assets_reserve = n.model["Link-r"].loc[:, reserve_asset_indexes]
        # Calculate storage bus cumulative reserve variables
        storage_buses = df["bus0"].loc[reserve_asset_indexes]
        bus_reserve = (
            assets_reserve.groupby(storage_buses.to_xarray())
            .sum()
            .rename({"bus0": "bus"})
        )
        # Get state of charge variables
        store_asset_indexes = n.stores[
            (n.stores.bus.isin(storage_buses)) & (n.stores.country == country)
        ].index
        store_asset_soc = (
            n.model["Store-e"]
            .sel(name=store_asset_indexes)
            .groupby(n.stores.loc[store_asset_indexes, "bus"].to_xarray())
            .sum()
        )
        # lhs = (sum_g,h reserve(g,h,t) * r_rating of asset g at bus h)
        # - (0.25 * sum_s,h soc(s,h,t) of store s at bus h)
        lhs = bus_reserve - store_asset_soc.mul(0.25)
        n.model.add_constraints(
            lhs, "<=", 0, name=f"updated_soc_constraint_storage_links_{country}"
        )


def __get_reserve_asset_indices(df: pd.DataFrame, country: str) -> pd.Index:
    """Filter reserve assets from n.df(c).

    Parameters
    ----------
    df : pd.DataFrame
        n.df(c) with c is pypsa component
    country : str
        Country for which the constraint is applied.

    Returns
    -------
    pd.Index
        Indices for all reserve assets, the extendable ones and the fixed ones.
    """
    reserve_indices = (
        df.query("r_rating > 0").loc[lambda df: df["country"] == country].index
    )
    ext_indices = (
        df.query("p_nom_extendable").loc[lambda df: df["country"] == country].index
    )
    ext_i = reserve_indices.intersection(ext_indices)
    fix_i = reserve_indices.difference(ext_i)
    return reserve_indices, ext_i, fix_i


def __get_asset_indices(df: pd.DataFrame) -> pd.Index:
    """Filter reserve assets from n.df(c).

    Parameters
    ----------
    df : pd.DataFrame
        n.df(c) with c is pypsa component

    Returns
    -------
    pd.Index
        Indices for the extendable ones and the fixed ones.
    """
    ext_i = df.query("p_nom_extendable").index
    fix_i = df.index.difference(ext_i)
    return ext_i, fix_i


# ---- Demand side management (DSM): config, components and constraints ----


def _is_number(value) -> bool:
    """Return True for finite real numbers (not bool, NaN or inf)."""
    return (
        isinstance(value, numbers.Real)
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def value_for_year(value, year: int) -> float:
    """Return the value of a DSM parameter for a year (scalar or dict keyed by year)."""
    if isinstance(value, dict):
        return float(value[year])
    return float(value)


def validate_dsm_config(scenario_configs: dict, years: list, countries: list):
    """Validate the ``demand_side_management`` section of the scenario config.

    Validation always runs, also when the resolution does not support DSM, so that a
    broken config is found early.

    Parameters
    ----------
    scenario_configs : dict
        Scenario configurations.
    years : list
        All model years.
    countries : list
        Countries in the model.

    Raises
    ------
    ValueError
        If the section is malformed.
    """
    section = scenario_configs.get(DSM_SECTION) or {}
    errors = []
    if not isinstance(section, dict):
        raise ValueError(f"'{DSM_SECTION}' must be a mapping of countries.")
    for country, country_cfg in section.items():
        prefix = f"{DSM_SECTION}.{country}"
        if country not in countries:
            errors.append(f"{prefix}: country is not in base_configs.regions.")
            continue
        if not isinstance(country_cfg, dict):
            errors.append(f"{prefix}: must be a mapping with 'activate' and 'loads'.")
            continue
        unknown = set(country_cfg) - {"activate", "loads"}
        if unknown:
            errors.append(f"{prefix}: unknown keys {sorted(unknown)}.")
        if not isinstance(country_cfg.get("activate", False), bool):
            errors.append(f"{prefix}.activate: must be true or false.")
        loads = country_cfg.get("loads")
        # Inactive countries may leave 'loads' empty, but any settings given are
        # still validated so that errors show up before DSM is switched on.
        if not country_cfg.get("activate", False) and not loads:
            continue
        if not isinstance(loads, dict) or not loads:
            errors.append(f"{prefix}.loads: must list at least one load type.")
            continue
        for profile_type, load_cfg in loads.items():
            errors.extend(
                _validate_load_cfg(f"{prefix}.loads.{profile_type}", load_cfg, years)
            )
    if errors:
        raise ValueError(
            "Invalid demand_side_management configuration:\n  - "
            + "\n  - ".join(errors)
        )


def _validate_load_cfg(prefix: str, load_cfg, years: list) -> list:
    """Validate the DSM settings of one load type and return a list of errors."""
    if not isinstance(load_cfg, dict):
        return [f"{prefix}: must be a mapping."]
    dsm_type = load_cfg.get("type")
    if dsm_type not in DSM_TYPES:
        return [f"{prefix}.type: must be one of {sorted(DSM_TYPES)}, got {dsm_type!r}."]
    errors = []
    missing = DSM_REQUIRED_KEYS[dsm_type] - set(load_cfg)
    if missing:
        errors.append(f"{prefix}: missing keys {sorted(missing)} for type {dsm_type}.")
    unknown = set(load_cfg) - DSM_REQUIRED_KEYS[dsm_type] - DSM_OPTIONAL_KEYS[dsm_type]
    if unknown:
        errors.append(f"{prefix}: unknown keys {sorted(unknown)} for type {dsm_type}.")
    for key in DSM_YEAR_KEYS:
        if key not in load_cfg:
            continue
        value = load_cfg[key]
        values = value if isinstance(value, dict) else {None: value}
        if isinstance(value, dict):
            missing_years = [y for y in years if y not in value]
            if missing_years:
                errors.append(f"{prefix}.{key}: missing years {missing_years}.")
        for year, v in values.items():
            label = f"{prefix}.{key}" + (f"[{year}]" if year is not None else "")
            if not _is_number(v):
                errors.append(f"{label}: must be a finite number, got {v!r}.")
            elif key in ("max_peak_share", "max_energy_share") and not 0 <= v <= 1:
                errors.append(f"{label}: must be between 0 and 1, got {v}.")
            elif key == "marginal_cost" and v < 0:
                errors.append(f"{label}: must not be negative, got {v}.")
            elif key == "max_hours" and v <= 0:
                errors.append(f"{label}: must be positive, got {v}.")
    if "window_hours" in load_cfg:
        window = load_cfg["window_hours"]
        if not isinstance(window, int) or isinstance(window, bool) or window < 1:
            errors.append(
                f"{prefix}.window_hours: must be a positive integer, got {window!r}."
            )
    return errors


def dsm_resolution_supported(scenario_configs: dict) -> bool:
    """Return True if the temporal resolution supports DSM (hourly ``nth_hour``).

    For developers only: if the environment variable ``PYPSA_SPICE_DSM_ANY_RESOLUTION``
    is set to ``1``, DSM is applied with any resolution so that the code can be tested
    on fast, coarse runs. A warning is logged, and the results of such runs are not
    valid.
    """
    resolution = scenario_configs["scenario_configs"]["resolution"]
    hourly = (
        str(resolution.get("method", "")).lower() == "nth_hour"
        and resolution.get("stepsize") == 1
    )
    if not hourly and os.environ.get(DSM_TEST_OVERRIDE_ENV) == "1":
        logger.warning(
            "%s=1: DSM is applied with resolution method=%s, stepsize=%s for "
            "testing only. DSM results of this run are not valid.",
            DSM_TEST_OVERRIDE_ENV,
            resolution.get("method"),
            resolution.get("stepsize"),
        )
        return True
    return hourly


def get_active_dsm_countries(
    scenario_configs: dict, countries: list, year: int
) -> dict:
    """Return the DSM settings of every country where DSM is active and supported.

    Logs a warning for every active country if the resolution does not support DSM.

    Returns
    -------
    dict
        ``{country: {profile_type: settings}}``, empty if DSM is skipped.
    """
    section = scenario_configs.get(DSM_SECTION) or {}
    active = {
        country: cfg["loads"]
        for country, cfg in section.items()
        if country in countries and cfg.get("activate", False)
    }
    if active and not dsm_resolution_supported(scenario_configs):
        resolution = scenario_configs["scenario_configs"]["resolution"]
        for country in active:
            logger.warning(
                "Demand side management is not implemented for resolution "
                "method=%s, stepsize=%s. Only nth_hour with stepsize 1 (hourly) is "
                "supported. DSM is skipped for %s in %s.",
                resolution.get("method"),
                resolution.get("stepsize"),
                country,
                year,
            )
        return {}
    return active


def remove_dsm_components(n: pypsa.Network):
    """Remove all DSM components, e.g. those carried over from the previous year."""
    for comp in ["Generator", "StorageUnit"]:
        df = n.components[comp].static
        if df.empty or "type" not in df:
            continue
        n.remove(comp, df.index[df["type"].isin(DSM_COMPONENT_TYPES)])


def add_dsm_components(
    n: pypsa.Network,
    full_loads: pd.DataFrame,
    dsm_settings: dict,
    year: int,
):
    """Add DSM shed generators and shift storage units for one year.

    Parameters
    ----------
    n : pypsa.Network
        Network to add the components to.
    full_loads : pd.DataFrame
        Loads with the full 8760 h profile before temporal aggregation. Index is the
        load name. Columns are ``bus``, ``country``, ``profile_type`` and one column
        per hour.
    dsm_settings : dict
        Output of :func:`get_active_dsm_countries`.
    year : int
        Model year.
    """
    meta_cols = ["bus", "country", "profile_type", "carrier", "node"]
    hour_cols = [c for c in full_loads.columns if c not in meta_cols]
    for country, loads_cfg in dsm_settings.items():
        for profile_type, cfg in loads_cfg.items():
            dsm_type = cfg["type"]
            sel = full_loads[
                (full_loads["country"] == country)
                & (full_loads["profile_type"] == profile_type)
            ]
            if sel.empty:
                logger.warning(
                    "DSM: no %s loads found in %s for %s. Nothing added.",
                    profile_type,
                    country,
                    year,
                )
                continue
            # peak of the load on each bus, from the full hourly profile
            peak = sel[hour_cols].astype(float).groupby(sel["bus"]).sum().max(axis=1)
            share = value_for_year(cfg["max_peak_share"], year)
            p_nom = (peak * share).clip(lower=0)
            p_nom = p_nom[p_nom > 0]
            if p_nom.empty:
                logger.warning(
                    "DSM: %s potential of %s in %s is zero in %s. Nothing added.",
                    dsm_type,
                    profile_type,
                    country,
                    year,
                )
                continue
            buses = p_nom.index
            carrier = DSM_TYPES[dsm_type]["carrier"]
            names = [f"{carrier} - {profile_type} - {bus}" for bus in buses]
            # Shed power and shift discharge (load reduction) cannot exceed the load
            # of this type on the bus in each snapshot, so DSM never supplies energy
            # to other loads.
            demand = (
                n.loads_t.p_set.reindex(columns=sel.index, fill_value=0.0)
                .T.groupby(sel["bus"])
                .sum()
                .T.reindex(columns=buses, fill_value=0.0)
            )
            p_max_pu = demand.div(p_nom).clip(lower=0.0, upper=1.0)
            p_max_pu.columns = names
            common = dict(
                bus=buses,
                country=country,
                carrier=carrier,
                type=DSM_TYPES[dsm_type]["type"],
                p_nom=p_nom.to_numpy(),
                p_nom_extendable=False,
                marginal_cost=value_for_year(cfg["marginal_cost"], year),
                build_year=year,
                lifetime=np.inf,
                p_max_pu=p_max_pu,
                dsm_profile_type=profile_type,
            )
            if dsm_type == "shed":
                n.add(
                    "Generator",
                    names,
                    dsm_max_energy_share=value_for_year(cfg["max_energy_share"], year),
                    **common,
                )
            else:
                window = int(cfg["window_hours"])
                max_hours = value_for_year(cfg.get("max_hours", window), year)
                if max_hours > window:
                    logger.warning(
                        "DSM: max_hours=%s of %s in %s is larger than "
                        "window_hours=%s and is clipped to %s.",
                        max_hours,
                        profile_type,
                        country,
                        window,
                        window,
                    )
                    max_hours = float(window)
                n.add(
                    "StorageUnit",
                    names,
                    max_hours=max_hours,
                    efficiency_store=1.0,
                    efficiency_dispatch=1.0,
                    standing_loss=0.0,
                    cyclic_state_of_charge=False,
                    state_of_charge_initial=p_nom.to_numpy() * max_hours / 2,
                    dsm_window_hours=float(window),
                    **common,
                )
            logger.info(
                "DSM: added %s %s on %d bus(es) for %s in %s (%.1f MW in total).",
                dsm_type,
                profile_type,
                len(buses),
                country,
                year,
                p_nom.sum(),
            )


def add_dsm_shed_energy_constraint(n: pypsa.Network, country: str):
    """Cap the weighted energy shed by DSM per bus.

    For every DSM shed generator g on bus b for load type k:
    sum_t w_t * p_g,t <= max_energy_share_g * sum_t w_t * sum_(l on b of type k) d_l,t

    The annual demand on the right hand side is the weighted load of the network, so
    both sides use the same snapshot weightings.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network with DSM shed generators (type DSMD).
    country : str
        Country for which the constraint is applied.
    """
    gens = n.generators
    if gens.empty or "type" not in gens:
        return
    shed = gens[(gens["type"] == "DSMD") & (gens["country"] == country)]
    if shed.empty:
        return
    weights = n.snapshot_weightings["generators"]
    loads = n.loads
    demand = n.loads_t.p_set.mul(weights, axis=0).sum()
    rhs = pd.Series(0.0, index=shed.index)
    for gen, row in shed.iterrows():
        mask = (loads["bus"] == row["bus"]) & (
            loads["profile_type"] == row["dsm_profile_type"]
        )
        rhs[gen] = row["dsm_max_energy_share"] * demand.reindex(loads.index[mask]).sum()
    weight_da = xr.DataArray(weights, dims="snapshot", coords={"snapshot": n.snapshots})
    p_shed = n.model["Generator-p"].loc[:, shed.index]
    lhs = (p_shed * weight_da).sum("snapshot")
    rhs_da = xr.DataArray(rhs.to_numpy(), dims="name", coords={"name": shed.index})
    print(f"....adding DSM shed energy constraint for {len(shed)} bus(es) in {country}")
    n.model.add_constraints(lhs <= rhs_da, name=f"dsm_shed_energy_{country}")


def add_dsm_shift_window_constraint(n: pypsa.Network, country: str):
    """Fix the state of charge of DSM shift storage units at every window boundary.

    Windows start at the first snapshot and are counted in modelled hours (the
    cumulative snapshot weightings), so a window closes at the snapshot where the
    cumulative hours reach a multiple of ``window_hours``. With hourly resolution this
    is every ``window_hours`` snapshots. At the last snapshot of every window (the last
    window may be shorter), the state
    of charge equals half of the energy capacity, which is also the initial state of
    charge. Load can therefore move earlier or later, and the energy balances to zero
    in each window.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network with DSM shift storage units (type DSMF).
    country : str
        Country for which the constraint is applied.
    """
    sus = n.storage_units
    if sus.empty or "type" not in sus:
        return
    shift = sus[(sus["type"] == "DSMF") & (sus["country"] == country)]
    if shift.empty:
        return
    n_snapshots = len(n.snapshots)
    soc = n.model["StorageUnit-state_of_charge"]
    # modelled hours at the end of each snapshot (state of charge uses store weights)
    cum_hours = n.snapshot_weightings["stores"].cumsum().to_numpy()
    for window, group in shift.groupby("dsm_window_hours"):
        window = int(window)
        # index of the window each snapshot ends in; a window closes where it changes
        window_idx = np.floor(cum_hours / window + 1e-9)
        ends = list(np.flatnonzero(np.diff(window_idx, prepend=0) > 0))
        if not ends or ends[-1] != n_snapshots - 1:
            ends.append(n_snapshots - 1)
        end_snapshots = n.snapshots[ends]
        half_energy = group["p_nom"] * group["max_hours"] / 2
        rhs = xr.DataArray(
            half_energy.to_numpy(), dims="name", coords={"name": group.index}
        )
        lhs = soc.loc[end_snapshots, group.index]
        print(
            f"....adding DSM shift window constraint ({window} h, {len(ends)} "
            f"windows) for {len(group)} bus(es) in {country}"
        )
        n.model.add_constraints(
            lhs == rhs, name=f"dsm_shift_window_{window}h_{country}"
        )
