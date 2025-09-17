# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Add brownfield capacities and future assets to a PyPSA network for the current year.

This script performs the following tasks:
- Imports a solved PyPSA network from previous year and updates brownfield capacities.
- Removes outdated assets and capacities below a specified threshold.
- Updates fuel prices and decommissions base assets according to input data.
- Adds new assets for the current year, including interconnectors, storage units,
  stores, loads, generators, links, and EV infrastructure.
- Handles temporal aggregation and snapshot management.
- Applies CO2 management options (e.g., CO2 price adjustments).
- Exports the updated network to a NetCDF file.

The script is designed to be run as part of a Snakemake workflow and supports
sector-specific asset addition for power, industry, and transport sectors.
"""

import logging
from typing import Any

import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    FilePath,
    configure_logging,
    get_capital_cost,
    get_link_availabilities,
    get_plant_availabilities,
    get_storage_units_inflows,
    get_store_min_availabilities,
    get_time_series_demands,
    update_ev_char_parameters,
    update_ev_store_parameters,
    update_storage_costs,
    update_tech_fact_table,
)

idx = pd.IndexSlice

logger = logging.getLogger(__name__)


def add_brownfield(n: pypsa.Network, year: int, threshold: float):
    """Import the solved network from previous year and update brownfield capacities.

    Parameters
    ----------
    n : pypsa.Network
        Solved Pypsa network from previous year
    year : int
        Current model year
    threshold : float
        Capacity threshold below which assets are removed to avoid numerical issues
    """
    # remove all loads
    n.remove("Load", n.loads.index)

    # remove co2 stores in CO2STORN
    co2_store = n.stores[n.stores.bus.str.contains("CO2STORN")].index
    n.remove("Store", co2_store)

    for c in n.iterate_components(["Link", "Generator", "Store", "StorageUnit"]):
        attr = "e" if c.name == "Store" else "p"
        # fix capacities from older assets
        c.df[attr + "_nom"] = c.df[attr + "_nom_opt"]
        c.df[attr + "_nom_extendable"] = False

        # re-enable theoretical store
        if c.name == "Store":
            tstor = c.df[
                (c.df.carrier.isin(["CO2"])) | (c.df.type.str.contains("_STORE"))
            ].index
            c.df.loc[tstor, attr + "_nom_extendable"] = True

        # re-enable theoretical generators
        if c.name == "Generator":
            tgen = c.df[c.df.type.str.contains("SUPPLY")].index
            c.df.loc[tgen, attr + "_nom_extendable"] = True

        for country in snakemake.params.country:
            # subtract co2 price from marginal cost of previous year fossil links
            if (
                c.name == "Link"
                and snakemake.params.co2_management[country]["option"] == "co2_price"
            ):
                years = snakemake.params.years
                i = years.index(year)
                year_previous = years[i - 1]
                co2_price_previous = snakemake.params.co2_management[country]["value"][
                    year_previous
                ]
                emit_links = c.df[c.df.bus2 == f"{country}_ATMP"].index
                c.df.loc[emit_links, "marginal_cost"] = c.df.loc[
                    emit_links, "marginal_cost"
                ] - c.df.loc[emit_links, "efficiency2"].mul(co2_price_previous)

        # remove assets whose build_year + lifetime < year
        assets_phased_out = c.df.index[c.df.build_year + c.df.lifetime < year]

        n.remove(c.name, assets_phased_out)

        if not assets_phased_out.empty:
            print(
                "Remove: \n"
                + "\n".join(assets_phased_out)
                + "\nsince assets reached end of their lifetime."
            )

        # remove capacities below a certain threshold (this is mainly to avoid
        # numerical issues)
        assets_below_threshold = c.df.index[
            (c.df[attr + "_nom_opt"] < threshold) & ~(c.df[attr + "_nom_extendable"])
        ]
        n.remove(c.name, assets_below_threshold)

    print("Finish adding brownfield capacities")


def update_decommission_base_assets(
    n: pypsa.Network, year: int, pow_decom: FilePath, ind_decom: FilePath
):
    """Update the brownfield capacities by decommissioning base assets.

    How?
    By subtracting base year assets' capacities (per year) by decommission fleet
    specified in the `decommission_capacity.csv` files.

    Parameters
    ----------
    n : pypsa.Network
        Pypsa network after adding brownfield
    year : int
        current model year
    pow_decom : FilePath
        Path to the `decommission_capacity.csv` in the power folder
    ind_decom : FilePath
        Path to the `decommission_capacity.csv` in the industry folder
    """
    pow_decom_df = pd.read_csv(pow_decom, index_col=["name"])
    ind_decom_df = pd.read_csv(ind_decom, index_col=["name"])
    decap_df = pd.concat([pow_decom_df, ind_decom_df], axis=0)
    for col in decap_df.columns:
        if col in ["country", "name", "class"]:
            decap_df[col] = decap_df[col].astype(str)
        else:
            decap_df[col] = decap_df[col].astype(float).fillna(0)

    for c in n.iterate_components(
        [
            "Store",
            "StorageUnit",
            "Link",
            "Generator",
        ]
    ):
        decap = decap_df[decap_df["class"] == c.name][str(year)]
        if len(decap) == 0:
            print(f"No decommissioning for {c.name}")
            continue
        attr = "e" if c.name == "Store" else "p"
        for plant in decap.index:
            if plant in c.df.index:
                if c.name == "Link":
                    c.df.loc[plant, f"{attr}_nom"] = c.df.loc[plant, f"{attr}_nom"] - (
                        decap.loc[plant] / c.df.loc[plant, "efficiency"]
                    )
                else:
                    c.df.loc[plant, f"{attr}_nom"] = (
                        c.df.loc[plant, f"{attr}_nom"] - decap.loc[plant]
                    )


def update_fuel_data(n: pypsa.Network, hubs: FilePath, year: int, currency: str):
    """Update fuel prices for supply generators based on the current year.

    Parameters
    ----------
    n : pypsa.Network
        Solved Pypsa network from previous year
    hubs : FilePath
        CSV file where all supply hubs in the network can be found
    year : int
        Current model year
    currency : str
        Currency of the model
    """
    print(f"Updating fuel prices for {year}")
    hubs_raw = pd.read_csv(hubs)
    hubs_df = hubs_raw[hubs_raw["year"] == year]
    hubs_df.index = hubs_df["country"] + "_" + hubs_df["supply_plant"]
    tgen = n.generators[n.generators.type.str.contains("SUPPLY")].index
    n.generators.loc[tgen, "marginal_cost"] = hubs_df.loc[
        tgen, f"fuel_cost [{currency}/MWh]"
    ].reindex(tgen)


def get_previous_year_red_hours(n: pypsa.Network) -> tuple[list, pd.DataFrame]:
    """Get the reduced hours and weights from previous year network.

    Parameters
    ----------
    n : pypsa.Network
        Pypsa network from previous optimization

    Returns
    -------
    tuple[list, pd.DataFrame]
        red_hours : list
            List of reduced hours from previous year network
        weights : pd.DataFrame
            DataFrame of snapshot weightings from previous year network
    """
    previous_year = int(n.snapshots.year.unique()[0])
    year_hours = pd.date_range(
        start=f"{previous_year}-01-01",
        end=f"{previous_year+1}-01-01",
        freq="h",
        inclusive="left",
    )
    red_hours = [year_hours.get_loc(x) for x in n.snapshots]
    weights = n.snapshot_weightings.reset_index(drop=True)
    return red_hours, weights


def update_brownfield_snapshots(
    n: pypsa.Network, red_hour: list, weights: pd.DataFrame, year: int
):
    """Update snapshots of the brownfield network to the current year's reduced hours.

    Parameters
    ----------
    n : pypsa.Network
        Pypsa network after adding brownfield
    red_hour : list
        List of reduced hours from previous year network
    weights : pd.DataFrame
        DataFrame of snapshot weightings from previous year network
    year : int
        Current model year
    """
    year_hours = pd.date_range(
        start=f"{year}-01-01",
        end=f"{year+1}-01-01",
        freq="h",
        inclusive="left",
    )
    reduced_year_hours = year_hours[red_hour]
    n_p = pypsa.Network(
        snakemake.input.previous_year_network,
    )
    n.set_snapshots(reduced_year_hours)
    n.snapshot_weightings = weights.set_index(n.snapshots)
    for c in n_p.iterate_components(["Link", "Generator", "Store", "StorageUnit"]):
        # copy time-dependent
        selection = n.component_attrs[c.name].type.str.contains(
            "series"
        ) & n.component_attrs[c.name].status.str.contains("Input")
        for tattr in n.component_attrs[c.name].index[selection]:
            n.import_series_from_dataframe(
                c.pnl[tattr].set_index(n.snapshots), c.name, tattr
            )


class AddFutureAssets:
    """Add future assets to the brownfield network for the current year."""

    def __init__(self, network: pypsa.Network, year: int, red_hours: list):
        """Initialize the AddFutureAssets class.

        Parameters
        ----------
        network : pypsa.Network
            Pypsa network after adding brownfield
        year : int
            Current model year
        red_hours : list
            List of reduced hours from previous year network
        """
        self.network = network
        self.year = year
        self.country = snakemake.params.country
        self.red_hours = red_hours
        # Getting path for all database files
        self.technologies_dir = snakemake.input.powerplant_type
        self.tech_cost_dir = snakemake.input.powerplant_costs
        self.storage_cost_path = snakemake.input.storage_costs
        self.dmd_profile_path = snakemake.input.dmd_profiles
        self.availability_dir = snakemake.input.pp_availability
        self.inflows_path = snakemake.input.stor_inflows
        self.interest = snakemake.params.interest
        self.currency = snakemake.params.currency

    # ========================== ADDING future interconnectors =========================

    def add_interconnectors(self):
        """Add interconnector capacity to the pyPSA network."""
        interconnectors = pd.read_csv(snakemake.input.interconnection).set_index("link")
        p_nom_max_min_columns = [
            x for x in interconnectors.columns if "p_nom_max" in x or "p_nom_min" in x
        ]
        interconnectors[p_nom_max_min_columns] = interconnectors[
            p_nom_max_min_columns
        ].astype("float64")
        interconnectors["LIFE"] = (
            50  # assuming lifetime of 50 years for interconnectors
        )
        interconnectors[f"CAP[{self.currency}/MW]"] = get_capital_cost(
            plant_type="ITCN",
            tech_costs=interconnectors.set_index("type"),
            interest=self.interest,
            currency=self.currency,
        ).values  # annualized capital cost
        # Make distribution grid expandable
        distribution_grid = interconnectors[
            (interconnectors.type == "ITCN")
            & (interconnectors.bus1.str.contains("LVELEC"))
        ].index
        interconnectors.loc[distribution_grid, "p_nom_extendable"] = True
        interconnectors.loc[distribution_grid, f"p_nom_min_{self.year}"] = 0
        interconnectors.loc[distribution_grid, f"p_nom_max_{self.year}"] = np.inf

        # Filter only expandable interconnectors
        interconnectors = interconnectors[interconnectors["p_nom_extendable"]]
        interconnectors.index = [
            s + f"_{str(self.year)}" for s in interconnectors.index
        ]

        self.network.add(
            class_name="Link",
            name=interconnectors.index,
            bus0=interconnectors["bus0"],
            bus1=interconnectors["bus1"],
            carrier=interconnectors["carrier"],
            type=interconnectors["type"],
            efficiency=interconnectors["efficiency"],
            p_max_pu=interconnectors["p_max_pu"],
            p_min_pu=interconnectors["p_min_pu"],
            p_nom_max=interconnectors[f"p_nom_max_{self.year}"],
            p_nom_min=interconnectors[f"p_nom_min_{self.year}"],
            p_nom_extendable=True,
            capital_cost=interconnectors[f"CAP[{self.currency}/MW]"],
            build_year=self.year,
            lifetime=interconnectors["LIFE"],
            country=interconnectors["country"],
        )

    # ================================= ADDING storage =================================

    def add_storage_capacity(self, storage_capacity_dir: FilePath):
        """Add storage capacity (StorageUnit) assets to the PyPSA network.

        Parameters
        ----------
        storage_capacity_dir : FilePath
            Path to the `storage_capacity.csv` file
        """
        # add StorageUnit (storage_capacity)
        storage_capacity = pd.read_csv(storage_capacity_dir)
        storage_capacity = update_tech_fact_table(
            tech_table=storage_capacity,
            technologies_dir=self.technologies_dir,
            tech_costs_dir=self.tech_cost_dir,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )
        storage_capacity = storage_capacity[(storage_capacity["p_nom_extendable"])]

        # add StorageUnit (storage_capacity)
        storage_unit = storage_capacity[(storage_capacity["class"] == "StorageUnit")]

        # StorageUnit with inflow (e.g. hydro dam, hydro reservoir)
        indexed_inflow_storage_unit = storage_unit[
            storage_unit["type"].isin(["HDAM", "HPHS"])
        ].set_index("name")

        if len(indexed_inflow_storage_unit) > 0:
            store_inflow = get_storage_units_inflows(
                storage_capacity=indexed_inflow_storage_unit.reset_index(),
                inflows_path=self.inflows_path,
            ).reindex(indexed_inflow_storage_unit.index)
            indexed_inflow_storage_unit.index = [
                s + f"_{str(self.year)}" for s in indexed_inflow_storage_unit.index
            ]
            store_inflow.index = indexed_inflow_storage_unit.index
            inflow = store_inflow.T.iloc[self.red_hours].set_index(
                self.network.snapshots
            )
            inflow.columns.name = "StorageUnit"

            self.network.add(
                class_name="StorageUnit",
                name=indexed_inflow_storage_unit.index,
                bus=indexed_inflow_storage_unit["bus"],
                carrier=indexed_inflow_storage_unit["carrier"],
                p_nom_extendable=True,
                p_nom_max=indexed_inflow_storage_unit["p_nom_max"],
                p_nom_min=indexed_inflow_storage_unit["p_nom_min"],
                inflow=inflow,
                p_max_pu=indexed_inflow_storage_unit["p_max_pu"],
                p_min_pu=indexed_inflow_storage_unit["p_min_pu"],
                efficiency_dispatch=indexed_inflow_storage_unit["efficiency"],
                efficiency_store=indexed_inflow_storage_unit["efficiency_store"],
                capital_cost=indexed_inflow_storage_unit["capital_cost"],
                marginal_cost=indexed_inflow_storage_unit["marginal_cost"],
                type=indexed_inflow_storage_unit["type"],
                max_hours=indexed_inflow_storage_unit["max_hours"],
                standing_loss=indexed_inflow_storage_unit["standing_loss"],
                build_year=self.year,
                lifetime=indexed_inflow_storage_unit["lifetime"],
                node=indexed_inflow_storage_unit["node"],
                inv_cost=indexed_inflow_storage_unit["inv_cost"],
                fom_cost=indexed_inflow_storage_unit["fom_cost"],
                r_rating=indexed_inflow_storage_unit["r_rating"],
                cyclic_state_of_charge=True,
                country=indexed_inflow_storage_unit["country"],
            )

        # StorageUnit without inflow (e.g. battery, industrial heat storage)
        indexed_storage_unit = storage_unit[
            ~storage_unit["type"].isin(["HDAM", "HPHS"])
        ].set_index("name")

        if len(indexed_inflow_storage_unit) > 0:
            indexed_storage_unit.index = [
                s + f"_{str(self.year)}" for s in indexed_storage_unit.index
            ]
            inflow = 0

            self.network.add(
                class_name="StorageUnit",
                name=indexed_storage_unit.index,
                bus=indexed_storage_unit["bus"],
                carrier=indexed_storage_unit["carrier"],
                p_nom_extendable=True,
                p_nom_max=indexed_storage_unit["p_nom_max"],
                p_nom_min=indexed_storage_unit["p_nom_min"],
                inflow=inflow,
                p_max_pu=indexed_storage_unit["p_max_pu"],
                p_min_pu=indexed_storage_unit["p_min_pu"],
                efficiency_dispatch=indexed_storage_unit["efficiency"],
                efficiency_store=indexed_storage_unit["efficiency_store"],
                capital_cost=indexed_storage_unit["capital_cost"],
                marginal_cost=indexed_storage_unit["marginal_cost"],
                type=indexed_storage_unit["type"],
                max_hours=indexed_storage_unit["max_hours"],
                standing_loss=indexed_storage_unit["standing_loss"],
                build_year=self.year,
                lifetime=indexed_storage_unit["lifetime"],
                node=indexed_storage_unit["node"],
                inv_cost=indexed_storage_unit["inv_cost"],
                fom_cost=indexed_storage_unit["fom_cost"],
                r_rating=indexed_storage_unit["r_rating"],
                cyclic_state_of_charge=True,
                country=indexed_storage_unit["country"],
            )

        # add charging and discharging links for stores (storage_energy)
        store_links = storage_capacity[(storage_capacity["class"] != "StorageUnit")]
        store_links = store_links.set_index("name")
        store_links = store_links[~store_links.efficiency_store.isna()]
        store_links.index = [s + f"_{str(self.year)}" for s in store_links.index]
        store_links["bus0"] = store_links["node"] + "_" + store_links["type"] + "N"
        store_links = store_links.rename(columns={"bus": "bus1"})

        # add discharging links
        self.network.add(
            class_name="Link",
            name=[i + "_DISCHARGE" for i in store_links.index],
            bus0=store_links["bus0"].values,
            bus1=store_links["bus1"].values,
            carrier=store_links["carrier"].values,
            efficiency=store_links["efficiency"].values,
            type=list(store_links["type"].str.upper() + "_DISCHARGE"),
            capital_cost=store_links["capital_cost"].values,
            marginal_cost=store_links["marginal_cost"].values,
            p_nom_extendable=True,
            p_nom_max=store_links["p_nom_max"].values,
            p_nom_min=store_links["p_nom_min"].values,
            build_year=self.year,
            lifetime=store_links["lifetime"].values,
            inv_cost=store_links["inv_cost"].values,
            fom_cost=store_links["fom_cost"].values,
            r_rating=store_links["r_rating"].values,
            country=store_links["country"].values,
        )

        # add store links
        self.network.add(
            class_name="Link",
            name=[i + "_STORE" for i in store_links.index],
            bus0=store_links["bus1"].values,
            bus1=store_links["bus0"].values,
            carrier=store_links["carrier"].values,
            efficiency=store_links["efficiency_store"].values,
            type=list(store_links["type"].str.upper() + "_STORE"),
            p_nom_extendable=True,
            marginal_cost=0
            + 1e-2
            + 2e-3 * (np.random.random(len(store_links["marginal_cost"].values)) - 0.5),
            build_year=self.year,
            lifetime=store_links["lifetime"].values,
            country=store_links["country"].values,
        )

    def add_storage_energy(self, storage_energy_dir: FilePath):
        """Add storage energy (store+links) assets to the PyPSA network.

        Parameters
        ----------
        storage_energy_dir : FilePath
            Path to the `storage_energy.csv` file
        """
        # Add store (storage_energy)
        storage_energy_raw = pd.read_csv(storage_energy_dir).set_index("store")
        storage_energy_raw["cyclic"] = storage_energy_raw["type"].apply(
            lambda x: x != "CO2STOR"
        )
        storage_energy = update_storage_costs(
            storage_energy_raw,
            storage_costs=self.storage_cost_path,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )

        storage_energy = storage_energy[
            (storage_energy["e_nom_extendable [True/False]"])
            | (storage_energy_raw["type"] == "CO2STOR")
        ]
        storage_energy.index = [s + f"_{str(self.year)}" for s in storage_energy.index]

        self.network.add(
            class_name="Store",
            bus=storage_energy["bus"],
            name=storage_energy.index,
            type=storage_energy["type"],
            carrier=storage_energy["carrier"],
            capital_cost=storage_energy["capital_cost"],
            marginal_cost=storage_energy["marginal_cost"],
            e_nom=storage_energy["e_nom [MWh]"],
            e_nom_extendable=False,
            standing_loss=storage_energy["standing_loss [%/hour]"],
            e_cyclic=storage_energy["cyclic"],
            build_year=self.year,
            lifetime=np.inf,
            inv_cost=storage_energy["inv_cost"],
            fom_cost=storage_energy["fom_cost"],
            country=storage_energy["country"],
        )

    # ================================== ADDING loads ==================================

    def add_load(self, loads: FilePath):
        """Add all kinds of load to the PyPSA network.

        Parameters
        ----------
        loads : FilePath
            Path to the `demand_profiles.csv` file
        """
        final_load = get_time_series_demands(loads, self.dmd_profile_path, self.year)
        final_load.reset_index(["country", "bus", "carrier", "node"], inplace=True)
        p_set = (
            final_load.T.loc[
                ~final_load.columns.isin(["country", "bus", "carrier", "node"])
            ]
            .iloc[self.red_hours]
            .astype(float)
            .set_index(self.network.snapshots)
        )
        p_set.columns.name = "Load"

        self.network.add(
            class_name="Load",
            name=final_load.index,
            bus=final_load["bus"],
            carrier=final_load["carrier"],
            p_set=p_set,
            country=final_load["country"],
        )

    # ================================ ADDING generators ===============================

    def add_generators(self, plants: FilePath):
        """Add power/industry/transport generators to the PyPSA network.

        Parameters
        ----------
        plants : FilePath
            Path to the `powerplants.csv` or `industrial_plants.csv` or
            `transport_plants.csv` file TODO: clarify
        """
        clean_pps = pd.read_csv(plants)
        clean_pps = update_tech_fact_table(
            tech_table=clean_pps,
            technologies_dir=self.technologies_dir,
            tech_costs_dir=self.tech_cost_dir,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )

        clean_pps = clean_pps[clean_pps["p_nom_extendable"]]
        clean_gen = clean_pps[clean_pps["class"] == "Generator"]
        if len(clean_gen) > 0:
            clean_gen = clean_gen.set_index("name")
            gen_avail = get_plant_availabilities(
                plants,
                self.availability_dir,
                self.technologies_dir,
            ).reindex(clean_gen.index)
            clean_gen.index = [s + f"_{str(self.year)}" for s in clean_gen.index]
            gen_avail.index = clean_gen.index
            p_max_pu = gen_avail.T.iloc[self.red_hours].set_index(
                self.network.snapshots
            )
            p_max_pu.columns.name = "Generator"
            self.network.add(
                class_name="Generator",
                name=clean_gen.index,
                bus=clean_gen["bus"],
                carrier=clean_gen["carrier"],
                p_max_pu=p_max_pu,
                p_nom_extendable=True,
                p_nom_max=clean_gen["p_nom_max"],
                p_nom_min=clean_gen["p_nom_min"],
                p_min_pu=clean_gen["p_min_pu"],
                ramp_limit_up=clean_gen["ramp_limit_up"],
                ramp_limit_down=clean_gen["ramp_limit_down"],
                efficiency=clean_gen["efficiency"],
                capital_cost=clean_gen["capital_cost"],
                marginal_cost=clean_gen["marginal_cost"],
                type=clean_gen["type"],
                build_year=self.year,
                lifetime=clean_gen["lifetime"],
                node=clean_gen["node"],
                inv_cost=clean_gen["inv_cost"],
                fom_cost=clean_gen["fom_cost"],
                r_rating=clean_gen["r_rating"],
                country=clean_gen["country"],
            )

    # ==================== ADDING power/heat or fuel converter links ===================

    def add_links(
        self,
        links: FilePath,
    ):
        """Add power/heat or fuel converter links to the PyPSA network.

        Parameters
        ----------
        links : FilePath
            Path to the `links.csv` file TODO: clarify
        """
        links_df = pd.read_csv(links)
        links_df = update_tech_fact_table(
            tech_table=links_df,
            technologies_dir=self.technologies_dir,
            tech_costs_dir=self.tech_cost_dir,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )
        links_df = links_df[links_df["p_nom_extendable"]]
        # remove storage links from the list of links
        links_df = links_df[links_df.efficiency_store.isna()]
        # fillna for all efficiency columns
        eff_columns = [x for x in links_df.columns if "efficiency" in x]
        links_df[eff_columns] = links_df[eff_columns].fillna(value=0)
        # add info if bus3 is missing in input dataframe
        if "bus3" not in links_df.columns:
            links_df["bus3"] = np.NAN
            links_df["efficiency3"] = 0
        links_df.set_index("link", inplace=True)
        links_avail = get_link_availabilities(
            links, self.availability_dir, self.technologies_dir
        ).reindex(links_df.index)
        # add year information to index
        links_df.index = [s + f"_{str(self.year)}" for s in links_df.index]
        links_avail.index = links_df.index
        p_max_pu = links_avail.T.iloc[self.red_hours].set_index(self.network.snapshots)
        p_max_pu.columns.name = "Link"

        self.network.add(
            class_name="Link",
            name=links_df.index,
            bus0=links_df["bus0"],
            bus1=links_df["bus1"],
            bus2=links_df["bus2"],
            bus3=links_df["bus3"],
            carrier=links_df["carrier"],
            efficiency=links_df["efficiency"],
            efficiency2=links_df["efficiency2"],
            efficiency3=links_df["efficiency3"],
            type=links_df["type"],
            capital_cost=links_df["capital_cost"],
            marginal_cost=links_df["marginal_cost"],
            p_max_pu=p_max_pu,
            p_min_pu=links_df["p_min_pu"],
            p_nom_max=links_df["p_nom_max"],
            p_nom_min=links_df["p_nom_min"],
            ramp_limit_up=links_df["ramp_limit_up"],
            ramp_limit_down=links_df["ramp_limit_down"],
            p_nom_extendable=True,
            build_year=self.year,
            lifetime=links_df["lifetime"],
            inv_cost=links_df["inv_cost"],
            fom_cost=links_df["fom_cost"],
            r_rating=links_df["r_rating"],
            country=links_df["country"],
        )

    # ========================== ADDING PEV links and storages =========================

    def add_ev_chargers(self):
        """Add new PEV chargers to the PyPSA network."""
        ev_links_df = pd.read_csv(snakemake.input.tra_pev_chargers).set_index("link")
        tech_cost_df = pd.read_csv(snakemake.input.powerplant_costs).set_index(
            ["powerplant_type", "country"]
        )
        ev_links = update_ev_char_parameters(
            tech_df=ev_links_df,
            year=self.year,
            ev_param_dir=snakemake.input.ev_parameters,
            cost_df=tech_cost_df,
            interest_rate=self.interest,
            currency=self.currency,
        )
        link_p_max_pu = get_link_availabilities(
            snakemake.input.tra_pev_chargers,
            self.availability_dir,
            self.technologies_dir,
        ).reindex(ev_links.index)
        ev_links.index = [s + f"_{str(self.year)}" for s in ev_links.index]
        link_p_max_pu.index = ev_links.index
        link_p_max_pu = link_p_max_pu.T.iloc[self.red_hours].set_index(
            self.network.snapshots
        )
        link_p_max_pu.columns.name = "Link"

        self.network.add(
            class_name="Link",
            name=ev_links.index,
            bus0=ev_links["bus0"],
            bus1=ev_links["bus1"],
            carrier=ev_links["carrier"],
            efficiency=ev_links["efficiency"],
            type=ev_links["type"],
            p_nom_max=ev_links["p_nom_char"],
            p_nom_min=ev_links["p_nom_char"],
            p_max_pu=link_p_max_pu,
            p_min_pu=0,
            p_nom_extendable=True,
            build_year=self.year,
            capital_cost=ev_links["capital_cost"],
            marginal_cost=ev_links["marginal_cost"],
            lifetime=ev_links["lifetime"],
            inv_cost=ev_links["inv_cost"],
            fom_cost=ev_links["fom_cost"],
            country=ev_links["country"],
        )

    def add_ev_storage(self):
        """Add PEV storages to the network."""
        ev_storages = pd.read_csv(snakemake.input.tra_pev_storages).set_index("name")
        ev_storages = update_ev_store_parameters(
            tech_table=ev_storages,
            year=self.year,
            ev_param_dir=snakemake.input.ev_parameters,
        )
        ev_storages.index = [s + f"_{str(self.year)}" for s in ev_storages.index]

        ev_store_e_min_pu = get_store_min_availabilities(
            store_dir=snakemake.input.tra_pev_storages, avail_dir=self.availability_dir
        )
        ev_store_e_min_pu.index = ev_storages.index
        ev_store_e_min_pu = ev_store_e_min_pu.T.iloc[self.red_hours].set_index(
            self.network.snapshots
        )
        ev_store_e_min_pu.columns.name = "Store"

        self.network.add(
            class_name="Store",
            bus=ev_storages["bus"],
            name=ev_storages.index,
            type=ev_storages["type"],
            carrier=ev_storages["carrier"],
            e_nom=ev_storages["ev_e_nom"],
            e_nom_extendable=False,
            e_min_pu=ev_store_e_min_pu,
            e_cyclic=True,
            build_year=self.year,
            lifetime=np.inf,
            country=ev_storages["country"],
        )

    def add_co2_option(self):
        """Add CO2 management options."""
        for country in self.country:
            co2_option = snakemake.params.co2_management[country]["option"]
            if co2_option == "co2_price":
                # assign co2 price to marginal cost of all links with bus 2
                # attached to ATMP bus (meaning emitting co2)
                co2price = snakemake.params.co2_management[country]["value"][self.year]
                print(
                    f"Adding CO2 Price as {co2price}{self.currency}/tCO2 for {country}"
                )
                emit_links = self.network.links[
                    self.network.links.bus2 == f"{country}_ATMP"
                ].index
                # Marginal cost of Links (CURRENCY/MWh_th):
                # CO2 price (CURRENCY/t_CO2) * emission factor (t_CO2/MWh_th)
                self.network.links.loc[
                    emit_links, "marginal_cost"
                ] += self.network.links.loc[emit_links, "efficiency2"].mul(co2price)


if __name__ == "__main__":
    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake  # pylint: disable=ungrouped-imports

        snakemake = mock_snakemake("add_brownfield", sector="p-i-t", years=2030)
    configure_logging(snakemake)
    sm_threshold = snakemake.params.remove_threshold
    sm_currency = snakemake.params.currency
    sm_year = int(snakemake.wildcards.years)
    # Define valid sector options
    VALID_SECTORS = {"p", "p-i", "p-t", "p-i-t"}
    # Validate sector value
    if snakemake.wildcards.sector not in VALID_SECTORS:
        raise ValueError(
            f"Wrong value for sector: '{snakemake.wildcards.sector}'. "
            f"Please choose one from following options: {', '.join(VALID_SECTORS)}"
        )
    sm_n = pypsa.Network(
        snakemake.input.previous_year_network,
    )
    sm_n.name = f"network_{snakemake.wildcards.sector}_{sm_year}"
    sm_red_hours, sm_weights = get_previous_year_red_hours(sm_n)
    add_brownfield(sm_n, sm_year, sm_threshold)
    update_decommission_base_assets(
        sm_n, sm_year, snakemake.input.pow_decom, snakemake.input.ind_decom
    )
    update_fuel_data(sm_n, snakemake.input.fuel_supplies, sm_year, sm_currency)
    # Executing each add functions
    sm_c = AddFutureAssets(network=sm_n, year=sm_year, red_hours=sm_red_hours)
    sm_c.add_interconnectors()
    sm_c.add_storage_capacity(
        storage_capacity_dir=snakemake.input.elec_storage_capacity
    )
    sm_c.add_storage_energy(storage_energy_dir=snakemake.input.elec_storage_energy)
    sm_c.add_load(loads=snakemake.input.elec_loads)
    sm_c.add_generators(
        plants=snakemake.input.elec_power_generators,
    )
    sm_c.add_links(
        links=snakemake.input.elec_power_links,
    )
    print(f"Finish adding {sm_year} assets for power sector")
    if "i" in snakemake.wildcards.sector:
        sm_c.add_load(loads=snakemake.input.ind_loads)
        sm_c.add_storage_capacity(
            storage_capacity_dir=snakemake.input.ind_storage_capacity
        )
        sm_c.add_storage_energy(storage_energy_dir=snakemake.input.ind_storage_energy)
        sm_c.add_generators(
            plants=snakemake.input.ind_generators,
        )
        sm_c.add_links(
            links=snakemake.input.ind_heat_links,
        )
        sm_c.add_links(
            links=snakemake.input.ind_fuel_conversion_links,
        )
        sm_c.add_links(
            links=snakemake.input.ind_dac_links,
        )
        print(f"Finish adding {sm_year} assets for industry sector")

    if "t" in snakemake.wildcards.sector:
        sm_c.add_load(loads=snakemake.input.tra_loads)
        sm_c.add_ev_chargers()
        sm_c.add_ev_storage()
        print(f"Finish adding {sm_year} assets for EV sector")
    sm_c.add_co2_option()
    sm_c.network.export_to_netcdf(path=snakemake.output.brownfield_network)
