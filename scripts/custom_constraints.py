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
- etc.
"""
import colorama
import numpy as np
import pandas as pd
import pypsa
import xarray as xr
from _helpers import FilePath
from pypsa.descriptors import expand_series
from pypsa.descriptors import get_switchable_as_dense as get_as_dense
from pypsa.optimization.compat import define_constraints, get_var


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

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


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
    n.model.add_constraints(lhs <= rhs, name=f"co2_cap_{country}")


def capacity_factor_constraint(n: pypsa.Network, country: str, cf_dict: dict):
    """Add constraint on capacity factor for specific generator types.

    Parameters
    ----------
    n : pypsa.Network
        PyPSA network object to which the constraint will be applied.
    country : str
        Country for which the CO2 cap is applied.
    cf_dict : dict
        Dictionary with generator types as keys and their corresponding capacity
        factor limits as values.
    """
    lhs = 0
    for gen_type in cf_dict:
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
                            n.model[f"{c}-p_nom"].rename({f"{c}-ext": c}).loc[ext_i]
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
                            lhs <= 0,
                            name=(
                                f"updated_capacity_factor_{gen_type}_{country}"
                                + "_constraint_ext"
                            ),
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
                            lhs <= rhs,
                            name=(
                                f"updated_capacity_factor_{gen_type}_{country}"
                                + "_constraint_fix"
                            ),
                        )


def thermal_must_run_constraint(
    n: pypsa.Network, min_must_run_ratio: float, country: str, sense: str = ">="
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
    sense : str, optional
        , by default ">=" TODO: REMOVE THIS?
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
            (df.index.astype(str).str.contains(country + "_"))
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
            lhs += gen_var.sum(c)
    # define constraint ------------------------------------------------------
    print(
        f"....add minimum thermal must run as {round(min_must_run_ratio*100)}% "
        f"of power load per snapshot for {country}"
    )
    n.model.add_constraints(lhs >= rhs, name=f"thermal_must_run_limit_{country}")


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
                gen_var = get_var(n, c, p_gen).loc[:, res_gen]
                weight_gen = xr.DataArray(
                    expand_series(n.snapshot_weightings.objective, df.index)
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
    define_constraints(
        n,
        lhs,
        math_symbol,
        rhs,
        f"RE share of {round(res_generation_share*100)}% in {country}",
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
            .sum("Generator")
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
    gen_var = get_var(n, "Generator", "p")
    weight_gen = xr.DataArray(
        expand_series(n.snapshot_weightings.objective, n.generators.index)
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
        gen_var = get_var(n, c, p_gen)
        weight_gen = xr.DataArray(
            expand_series(n.snapshot_weightings.objective, df.index)
        )
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

    define_constraints(
        n,
        lhs,
        ">=",
        rhs,
        f"Energy Independence Constraint of >= {ei_frac} for {country}",
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
                    n.df(c).loc[ext_i, "efficiency"].rename({f"{c}": f"{c}-ext"})
                )
            lhs += (
                n.model.variables[f"{c}-p_nom"].sel({f"{c}-ext": ext_i})
                * r_rating
                * (eff if eff is not None else 1)
            ).sum(f"{c}-ext")

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
        assets_sumed_reserve = (assets_reserve * r_rating).sum(f"{c}")
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
        # capacity factor CF
        capacity_factor_ext = n.generators_t.p_max_pu[vres_i.intersection(ext_i)]
        # capacities of VRES
        p_nom_vres = (
            n.model["Generator-p_nom"]
            .loc[vres_i.intersection(ext_i)]
            .rename({"Generator-ext": "Generator"})
        )
        lhs += (p_nom_vres * (-ep_vre * capacity_factor_ext)).sum("Generator")

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
        lhs += (-storage_load.mul(store_eff)).sum("Link")
    if not storageunits_id.empty:
        SU_store = n.model["StorageUnit-p_store"].loc[:, storageunits_id]
        lhs += (-SU_store).sum("StorageUnit")

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
    n.model.add_constraints(lhs >= rhs, name=f"dynamic_reserve_margin_{country}")


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
            # variable assets
            capacity_variable = n.model[f"{c}-p_nom"].rename({f"{c}-ext": c}).loc[ext_i]
            lhs = (
                dispatch.loc[:, ext_i]
                + reserve.loc[:, ext_i]
                - capacity_variable.mul(p_max_pu.loc[:, ext_i])
            )
            # variable assets constraint
            n.model.add_constraints(
                lhs <= 0, name=f"updated_capacity_constraint_{c.lower()}s_ext_{country}"
            )
        if not fix_i.empty:
            # fix assets
            capacity_fixed = df["p_nom"].loc[fix_i]
            lhs = dispatch.loc[:, fix_i] + reserve.loc[:, fix_i]
            rhs = p_max_pu.loc[:, fix_i].mul(capacity_fixed)
            # fix assets constraint
            n.model.add_constraints(
                lhs <= rhs,
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
    eff = xr.DataArray(df.loc[reserve_asset_indices, "efficiency"])
    # Get p_max_pu
    p_max_pu = get_as_dense(n, "Link", "p_max_pu")
    if not ext_i.empty:
        # variable assets
        capacity_variable = (
            n.model["Link-p_nom"].rename({"Link-ext": "Link"}).loc[ext_i]
        )
        lhs = (
            dispatch.loc[:, ext_i].mul(eff.loc[ext_i])
            + reserve.loc[:, ext_i]
            - capacity_variable.mul(p_max_pu.loc[:, ext_i]).mul(eff.loc[ext_i])
        )
        # variable assets constraint
        n.model.add_constraints(
            lhs <= 0, name=f"updated_capacity_constraint_links_ext_{country}"
        )
    if not fix_i.empty:  # fix assets
        capacity_fixed = df["p_nom"].loc[fix_i]
        lhs = dispatch.loc[:, fix_i].mul(eff.loc[fix_i]) + reserve.loc[:, fix_i]
        rhs = p_max_pu.loc[:, fix_i].mul(capacity_fixed).mul(eff.loc[fix_i])
        # fix assets constraint
        n.model.add_constraints(
            lhs <= rhs, name=f"updated_capacity_constraint_links_fix_{country}"
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
            lhs <= 0, name=f"updated_soc_constraint_storage_units_{country}"
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
            .sel({"Store": store_asset_indexes})
            .groupby(n.stores.loc[store_asset_indexes, "bus"].to_xarray())
            .sum()
        )
        # lhs = (sum_g,h reserve(g,h,t) * r_rating of asset g at bus h)
        # - (0.25 * sum_s,h soc(s,h,t) of store s at bus h)
        lhs = bus_reserve - store_asset_soc.mul(0.25)
        n.model.add_constraints(
            lhs <= 0, name=f"updated_soc_constraint_storage_links_{country}"
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
