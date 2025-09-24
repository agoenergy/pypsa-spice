# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Add all necessary components to the pyPSA network for the base year.

This including buses, generators, loads, storage units, links, and interconnectors.
It also handles temporal aggregation and CO2 management options.
"""

import logging
from typing import Any

import numpy as np
import pandas as pd
import pypsa
import tsam.timeseriesaggregation as tsam
from _helpers import (
    FilePath,
    configure_logging,
    filter_selected_countries_and_regions,
    get_buses,
    get_link_availabilities,
    get_plant_availabilities,
    get_reduced_hours,
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


class AddBaseNetwork:
    """Add base year assets to pypsa network."""

    def __init__(self, network: pypsa.Network, year: int):
        self.network = network
        self.year = year
        self.country_region = snakemake.params.country_region
        self.country = list(self.country_region.keys())
        # Getting path for all database files
        self.technologies_dir = snakemake.input.powerplant_type
        self.tech_cost_dir = snakemake.input.powerplant_costs
        self.storage_cost_path = snakemake.input.storage_costs
        self.dmd_profile_path = snakemake.input.dmd_profiles
        self.availability_dir = snakemake.input.pp_availability
        self.inflows_path = snakemake.input.stor_inflows
        self.interest = snakemake.params.interest
        self.currency = snakemake.params.currency

    # ================================== ADDING buses ==================================
    def add_buses(self, bus_dir: FilePath):
        """Add buses to the pyPSA network.

        Parameters
        ----------
        bus_dir : FilePath
            CSV file where all buses can be found
        """
        bus_df = pd.read_csv(bus_dir)
        bus_df = filter_selected_countries_and_regions(
            df=bus_df, column="node", country_region=self.country_region, buses_csv=True
        )
        bus_df = get_buses(bus_df=bus_df)
        self.network.add(
            class_name="Bus",
            name=bus_df.index,
            carrier=bus_df["carrier"],
            country=bus_df["country"],
            x=bus_df["x"], 
            y=bus_df["y"]
        )

    # ============================= ADDING atmosphere store ============================
    def add_atmosphere(self):
        """Add a theoretical atmosphere store to PyPSA network."""
        for country in self.country:
            # add atmosphere store
            self.network.add(
                class_name="Store",
                name=f"{country}_atmosphere",
                e_nom_extendable=True,
                # e_min_pu=-1, # activate to allow negative emission # noqa: E800
                carrier="CO2",
                bus=f"{country}_ATMP",
                country=country,
            )

    # ============================= ADDING network carriers ============================
    def add_network_carriers(self):
        """Add carriers to network and set emission factor for CO2 carrier."""
        bus_df = pd.DataFrame(self.network.buses.carrier)
        bus_df["emission factor"] = bus_df["carrier"].apply(
            lambda x: -1 if x == "CO2" else 0
        )
        carrier_df = bus_df[["carrier", "emission factor"]].reset_index(drop=True)
        # add extra missing carriers, which are normally not in bus_df
        missing_carriers = [
            "Solar",
            "Wind",
            "Water",
            "Geothermal",
            "Low_Heat",
            "High_Heat",
            "Bit-imp",
            "Gas-imp",
            "Oil-imp",
            "Uranium-imp",
            "ENS",
            "Ammonia"
        ]
        missing_carriers_df = pd.DataFrame(missing_carriers, columns=["carrier"])
        missing_carriers_df["emission factor"] = 0
        carrier_df = pd.concat([carrier_df, missing_carriers_df], axis=0)
        carrier_df = carrier_df.drop_duplicates().set_index("carrier")
        # add carriers to network
        self.network.add(
            class_name="Carrier",
            name=carrier_df.index,
            co2_emissions=carrier_df["emission factor"],
        )

    # =========================== ADDING Base interconnector ===========================
    def add_interconnectors(self):
        """Add interconnector capacity to the pyPSA network."""
        interc = pd.read_csv(snakemake.input.interconnection)
        interc = filter_selected_countries_and_regions(
            df=interc,
            column="link",
            country_region=self.country_region,
        ).set_index("link")

        self.network.add(
            class_name="Link",
            name=interc.index,
            bus0=interc["bus0"],
            bus1=interc["bus1"],
            carrier=interc["carrier"],
            type=interc["type"],
            efficiency=interc["efficiency"],
            p_max_pu=interc["p_max_pu"],
            p_min_pu=interc["p_min_pu"],
            p_nom=interc["p_nom"],
            p_nom_extendable=False,
            build_year=self.year,
            lifetime=np.inf,
            country=interc["country"],
        )

        # Distribution grid activation
        distribution_grid = self.network.links[
            (self.network.links.type == "ITCN")
            & (self.network.links.bus1.str.contains("LVELEC"))
        ].index
        # Allow extendable distribution grid capacities as these number are often
        # endogenously optimized
        self.network.links.loc[distribution_grid, "p_nom_extendable"] = True

    # ========================= ADDING fuel_supplies generator =========================
    def add_fuel_supplies(self):
        """Add fuel_supplies to the pyPSA network."""
        hubs_raw = pd.read_csv(snakemake.input.fuel_supplies)
        hubs_raw = hubs_raw[(hubs_raw["country"].str.contains("|".join(self.country)))]
        hubs_df = hubs_raw[hubs_raw["year"] == self.year]
        hubs_df.index = hubs_df["country"] + "_" + hubs_df["supply_plant"]
        self.network.add(
            class_name="Generator",
            name=hubs_df.index,
            bus=hubs_df["bus"],
            carrier=hubs_df["carrier"],
            marginal_cost=hubs_df[f"fuel_cost [{self.currency}/MWh]"],
            type=hubs_df["carrier"].str.upper() + "_SUPPLY",
            p_nom_extendable=True,
            country=hubs_df["country"],
        )

    # ================================= ADDING storage =================================
    def add_storage_capacity(self, storage_capacity_dir: FilePath):
        """Add storage_capacity (StorageUnit) assets to PyPSA network."""
        # add StorageUnit (storage_capacity)
        storage_capacity = pd.read_csv(storage_capacity_dir)
        storage_capacity = filter_selected_countries_and_regions(
            df=storage_capacity,
            column="node",
            country_region=self.country_region,
        )
        storage_capacity = update_tech_fact_table(
            tech_table=storage_capacity,
            technologies_dir=self.technologies_dir,
            tech_costs_dir=self.tech_cost_dir,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )

        # add StorageUnit (storage_capacity)
        storage_unit = storage_capacity[storage_capacity["class"] == "StorageUnit"]

        # StorageUnit with inflow (e.g. hydro dam, hydro reservoir)
        indexed_inflow_storage_unit = storage_unit[
            storage_unit["type"].isin(["HDAM", "HPHS"])
        ].set_index("name")

        if len(indexed_inflow_storage_unit) > 0:
            inflow = (
                get_storage_units_inflows(
                    storage_capacity=indexed_inflow_storage_unit.reset_index(),
                    inflows_path=self.inflows_path,
                )
                .reindex(indexed_inflow_storage_unit.index)
                .T.set_index(self.network.snapshots)
            )

            self.network.add(
                class_name="StorageUnit",
                name=indexed_inflow_storage_unit.index,
                bus=indexed_inflow_storage_unit["bus"],
                carrier=indexed_inflow_storage_unit["carrier"],
                p_nom=indexed_inflow_storage_unit["p_nom"],
                p_nom_extendable=False,
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
                lifetime=np.inf,
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
        inflow = 0

        if len(indexed_inflow_storage_unit) > 0:
            self.network.add(
                class_name="StorageUnit",
                name=indexed_storage_unit.index,
                bus=indexed_storage_unit["bus"],
                carrier=indexed_storage_unit["carrier"],
                p_nom=indexed_storage_unit["p_nom"],
                p_nom_extendable=False,
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
                lifetime=np.inf,
                node=indexed_storage_unit["node"],
                inv_cost=indexed_storage_unit["inv_cost"],
                fom_cost=indexed_storage_unit["fom_cost"],
                r_rating=indexed_storage_unit["r_rating"],
                cyclic_state_of_charge=True,
                country=indexed_storage_unit["country"],
            )

        # add charging and discharging links for stores (storage_energy)
        store_links = storage_capacity[storage_capacity["class"] != "StorageUnit"]
        store_links = store_links.set_index("name")
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
            p_nom=store_links["p_nom"].values,
            p_nom_extendable=False,
            build_year=self.year,
            lifetime=np.inf,
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
            p_nom=store_links["p_nom"].values,
            p_nom_extendable=False,
            marginal_cost=0
            + 1e-2
            + 2e-3 * (np.random.random(len(store_links["marginal_cost"].values)) - 0.5),
            build_year=self.year,
            lifetime=np.inf,
            country=store_links["country"].values,
        )

    def add_storage_energy(self, storage_energy_dir: FilePath):
        """Add storage_energy (store+links) assets to PyPSA network."""
        # Add store (storage_energy)
        storage_energy_raw = pd.read_csv(storage_energy_dir).set_index("store")
        storage_energy_raw = filter_selected_countries_and_regions(
            df=storage_energy_raw.reset_index(),
            column="store",
            country_region=self.country_region,
        ).set_index("store")
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

    # =================================== ADDING load ==================================
    def add_load(self, load_dir: FilePath):
        """Add loads of all kind to the PyPSA network."""
        load_df = filter_selected_countries_and_regions(
            df=pd.read_csv(load_dir),
            column="node",
            country_region=self.country_region,
        )
        final_load = get_time_series_demands(load_df, self.dmd_profile_path, self.year)
        final_load.reset_index(["country", "bus", "carrier", "node"], inplace=True)
        self.network.add(
            class_name="Load",
            name=final_load.index,
            bus=final_load["bus"],
            carrier=final_load["carrier"],
            p_set=final_load.T[
                ~final_load.columns.isin(["country", "bus", "carrier", "node"])
            ]
            .astype(float)
            .set_index(self.network.snapshots),
            country=final_load["country"],
        )

    # ================================ ADDING generators ===============================
    def add_generators(self, plant_dir: FilePath):
        """Add power/industry/transport generators to the pyPSA network.

        Parameters
        ----------
        plant_dir : FilePath
            CSV file containing plants' data.
        """
        clean_pps = pd.read_csv(plant_dir)
        clean_pps = filter_selected_countries_and_regions(
            df=clean_pps,
            column="node",
            country_region=self.country_region,
        )
        clean_pps = update_tech_fact_table(
            tech_table=clean_pps,
            technologies_dir=self.technologies_dir,
            tech_costs_dir=self.tech_cost_dir,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )
        clean_gen = clean_pps[clean_pps["class"] == "Generator"]
        if len(clean_gen) > 0:
            clean_gen = clean_gen.set_index("name")
            self.network.add(
                class_name="Generator",
                name=clean_gen.index,
                bus=clean_gen["bus"],
                carrier=clean_gen["carrier"],
                p_nom=clean_gen["p_nom"],
                p_max_pu=get_plant_availabilities(
                    clean_pps.reset_index(),
                    self.availability_dir,
                    self.technologies_dir,
                )
                .reindex(clean_gen.index)
                .T.set_index(self.network.snapshots),
                p_nom_extendable=False,
                p_min_pu=clean_gen["p_min_pu"],
                ramp_limit_up=clean_gen["ramp_limit_up"],
                ramp_limit_down=clean_gen["ramp_limit_down"],
                efficiency=clean_gen["efficiency"],
                capital_cost=clean_gen["capital_cost"],
                marginal_cost=clean_gen["marginal_cost"],
                type=clean_gen["type"],
                build_year=self.year,
                lifetime=np.inf,
                node=clean_gen["node"],
                inv_cost=clean_gen["inv_cost"],
                fom_cost=clean_gen["fom_cost"],
                r_rating=clean_gen["r_rating"],
                country=clean_gen["country"],
            )

    # ==================== ADDING power/heat or fuel converter links ===================
    def add_links(self, links_dir: FilePath):
        """Add power/heat or fuel converter links to the PyPSA network."""
        links_df = pd.read_csv(links_dir)
        links_df = filter_selected_countries_and_regions(
            df=links_df,
            column="link",
            country_region=self.country_region,
        )
        links_df = update_tech_fact_table(
            tech_table=links_df,
            technologies_dir=self.technologies_dir,
            tech_costs_dir=self.tech_cost_dir,
            year=self.year,
            interest=self.interest,
            currency=self.currency,
        )
        # ensure all efficiency columns are filled
        eff_columns = [x for x in links_df.columns if "efficiency" in x]
        links_df[eff_columns] = links_df[eff_columns].fillna(value=0)
        if "bus3" not in links_df.columns:
            links_df["bus3"] = np.NAN
            links_df["efficiency3"] = 0
        links_df = links_df.set_index("link")
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
            p_nom=links_df["p_nom"],
            p_max_pu=get_link_availabilities(
                links_df.reset_index(), self.availability_dir, self.technologies_dir
            )
            .reindex(links_df.index)
            .T.set_index(self.network.snapshots),
            p_min_pu=links_df["p_min_pu"],
            ramp_limit_up=links_df["ramp_limit_up"],
            ramp_limit_down=links_df["ramp_limit_down"],
            p_nom_extendable=False,
            build_year=self.year,
            lifetime=np.inf,
            inv_cost=links_df["inv_cost"],
            fom_cost=links_df["fom_cost"],
            r_rating=links_df["r_rating"],
            country=links_df["country"],
        )

    # ========================== ADDING pev links and storages =========================
    def add_ev_chargers(self):
        """Add PEV chargers to the PyPSA network."""
        ev_links_df = pd.read_csv(snakemake.input.tra_pev_chargers).set_index("link")
        tech_costs_df = pd.read_csv(snakemake.input.powerplant_costs)
        tech_costs_df = tech_costs_df[
            (tech_costs_df["country"].str.contains("|".join(self.country)))
        ].set_index(["powerplant_type", "country"])
        ev_links_df = filter_selected_countries_and_regions(
            df=ev_links_df.reset_index(),
            column="link",
            country_region=self.country_region,
        ).set_index("link")
        ev_links_df = update_ev_char_parameters(
            tech_df=ev_links_df,
            year=self.year,
            ev_param_dir=snakemake.input.ev_parameters,
            cost_df=tech_costs_df,
            interest_rate=self.interest,
            currency=self.currency,
        )
        link_p_max_pu = (
            get_link_availabilities(
                ev_links_df.reset_index(),
                self.availability_dir,
                self.technologies_dir,
            )
            .reindex(ev_links_df.index)
            .T.set_index(self.network.snapshots)
        )

        self.network.add(
            class_name="Link",
            name=ev_links_df.index,
            bus0=ev_links_df["bus0"],
            bus1=ev_links_df["bus1"],
            carrier=ev_links_df["carrier"],
            efficiency=ev_links_df["efficiency"],
            type=ev_links_df["type"],
            p_nom=ev_links_df["p_nom_char"],
            p_max_pu=link_p_max_pu,
            p_min_pu=0,
            p_nom_extendable=False,
            build_year=self.year,
            capital_cost=ev_links_df["capital_cost"],
            marginal_cost=ev_links_df["marginal_cost"],
            lifetime=ev_links_df["lifetime"],
            inv_cost=ev_links_df["inv_cost"],
            fom_cost=ev_links_df["fom_cost"],
            country=ev_links_df["country"],
        )

    def add_ev_storage(self):
        """Add PEV storages to network."""
        ev_storages = pd.read_csv(snakemake.input.tra_pev_storages).set_index("name")
        ev_storages = filter_selected_countries_and_regions(
            df=ev_storages,
            column="node",
            country_region=self.country_region,
        )
        ev_storages = update_ev_store_parameters(
            tech_table=ev_storages,
            year=self.year,
            ev_param_dir=snakemake.input.ev_parameters,
        )

        ev_store_e_min_pu = (
            get_store_min_availabilities(
                ev_storages.reset_index(),
                avail_dir=self.availability_dir,
            )
            .reindex(ev_storages.index)
            .T.set_index(self.network.snapshots)
        )

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

    def add_load_shedding(self):
        """Add load shedding to the network to ensure successful optimisation."""
        self.network.add(
            class_name="Generator",
            name=[f"LOSTLOAD - {x}" for x in self.network.buses.index],
            bus=self.network.buses.index,
            country=[x.split("_")[0] for x in self.network.buses.index],
            carrier="ENS",
            # sign=1e-3,  # Adjust sign to measure p or p_nom in kW or MW  # noqa: E800
            marginal_cost=1e5,  # Eur/MWh
            # intersect between macroeconomic and survey-based willingness to pay
            # http://journal.frontiersin.org/article/10.3389/fenrg.2015.00055/full
            type="LSLO",
            p_nom=1e9,  # MW
        )

    def add_noisy_cost(self):
        """Add noisy cost for faster solving time."""
        for t in self.network.iterate_components(self.network.one_port_components):
            non_co2_assets = t.df[t.df.carrier != "CO2"].index
            if "marginal_cost" in t.df:
                t.df.loc[non_co2_assets, "marginal_cost"] += 1e-2 + 2e-3 * (
                    np.random.random(len(t.df.loc[non_co2_assets])) - 0.5
                )
            if "capital_cost" in t.df:
                t.df.loc[non_co2_assets, "capital_cost"] += 1e-1 + 2e-2 * (
                    np.random.random(len(t.df.loc[non_co2_assets])) - 0.5
                )

    def remove_abundant_components(self):
        """Remove components with 0 capacities and non extendable status."""
        for comp in self.network.iterate_components(
            ["Link", "Generator", "Store", "StorageUnit"]
        ):
            attr = "e" if comp.name == "Store" else "p"
            # remove asset with capacity = 0 and non extendable
            assets_below_threshold = comp.df.index[
                (comp.df[attr + "_nom"] < 0.1) & ~(comp.df[attr + "_nom_extendable"])
            ]
            self.network.remove(comp.name, assets_below_threshold)

    def add_temporal_aggregation(self, method: str = "nth_hour"):
        """Aggregate the timeseries by stepsize.

        Parameters
        ----------
        method : str, optional
            'nth_hour' or 'clustered' from the config file, by default "nth_hour"
        """
        if (method).lower() == "clustered":
            solve_name = (snakemake.params.solve_name).lower()
            num_days = snakemake.params.numDays
            period_length = 24
            # Get all time-dependent data
            dfs = [
                pnl
                for c in self.network.iterate_components()
                for attr, pnl in c.pnl.items()
                if not pnl.empty and attr != "e_min_pu"
            ]
            df = pd.concat(dfs, axis=1)
            # Reset columns to flat index
            df = df.T.reset_index(drop=True).T

            # Normalise all time-dependent data
            annual_max = df.max().replace(0, 1)
            df = df.div(annual_max, level=0)
            # Get representative segments
            agg = tsam.TimeSeriesAggregation(
                df,
                hoursPerPeriod=period_length,
                noTypicalPeriods=num_days,
                rescaleClusterPeriods=False,
                solver=solve_name,
            )
            clustered = agg.createTypicalPeriods()
            map_snapshots_to_periods = agg.indexMatching()
            map_snapshots_to_periods["day_of_year"] = (
                map_snapshots_to_periods.index.day_of_year
            )
            cluster_weights = agg.clusterPeriodNoOccur
            cluster_center_indices = agg.clusterCenterIndices
            # pandas Day of year starts at 1, clusterCenterIndices at 0

            new_snapshots = map_snapshots_to_periods[
                (map_snapshots_to_periods.day_of_year - 1).isin(cluster_center_indices)
            ]
            new_snapshots["weightings"] = (
                new_snapshots["PeriodNum"].map(cluster_weights).astype(float)
            )
            clustered.set_index(new_snapshots.index, inplace=True)
            # last hour of typical period
            last_hour = new_snapshots[new_snapshots["TimeStep"] == period_length - 1]
            # first hour
            first_hour = new_snapshots[new_snapshots["TimeStep"] == 0]

            # add typical period name and last hour to mapping original snapshot
            map_snapshots_to_periods["RepresentativeDay"] = map_snapshots_to_periods[
                "PeriodNum"
            ].map(last_hour.set_index(["PeriodNum"])["day_of_year"].to_dict())
            map_snapshots_to_periods[
                "last_hour_RepresentativeDay"
            ] = map_snapshots_to_periods["PeriodNum"].map(
                last_hour.reset_index()
                .set_index(["PeriodNum"])["day_of_year"]
                .to_dict()
            )
            map_snapshots_to_periods[
                "first_hour_RepresentativeDay"
            ] = map_snapshots_to_periods["PeriodNum"].map(
                first_hour.reset_index()
                .set_index(["PeriodNum"])["day_of_year"]
                .to_dict()
            )

            snapshot_weightings = self.network.snapshot_weightings.mul(
                new_snapshots.weightings, axis=0
            )
            self.network.snapshot_weightings = snapshot_weightings

            # Define a series used for aggregation, mapping each hour in
            # n.snapshots to the closest previous time step in
            # snapshot_weightings.index
            aggregation_map = (
                pd.Series(
                    snapshot_weightings.index.get_indexer(self.network.snapshots),
                    index=self.network.snapshots,
                )
                .replace(-1, np.nan)
                .ffill()
                .astype(int)
                .map(lambda i: snapshot_weightings.index[i])
            )

            # save map original snapshots to typical periods
            self.network.cluster = map_snapshots_to_periods
            self.network.set_snapshots(new_snapshots.index)
            # Aggregation all time-varying data.
            for comp in self.network.iterate_components():
                pnl = getattr(self.network, comp.list_name + "_t")
                for k, df in comp.pnl.items():
                    if not df.empty:
                        if comp.list_name == "stores" and k == "e_max_pu":
                            pnl[k] = df.groupby(aggregation_map).min()
                        elif comp.list_name == "stores" and k == "e_min_pu":
                            pnl[k] = df.groupby(aggregation_map).max()
                        else:
                            pnl[k] = df.groupby(aggregation_map).mean()
        else:
            stepsize = snakemake.params.stepsize
            if stepsize > 1:
                print(f"Aggregating the network at every {stepsize}th hour")
            elec_load = (
                self.network.loads_t.p_set.T.groupby(self.network.loads.carrier).sum().T
            )
            elec_load = elec_load["Electricity"]
            red_hours, weights = get_reduced_hours(step_size=stepsize, load=elec_load)
            self.network.set_snapshots(self.network.snapshots[red_hours])
            self.network.snapshot_weightings = weights.set_index(self.network.snapshots)

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
                # marginal_cost of Links (CURRENCY/MWh_th) = co2 price (CURRENCY/tCO2) *
                # emission factors (tCO2/MWh_th)
                self.network.links.loc[
                    emit_links, "marginal_cost"
                ] += self.network.links.loc[emit_links, "efficiency2"].mul(co2price)


if __name__ == "__main__":
    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake  # pylint: disable=ungrouped-imports

        snakemake = mock_snakemake("add_baseyear", sector="p-i-t", years=2030)
    # Getting global config params
    configure_logging(snakemake)
    selected_year = int(snakemake.wildcards.years)
    # Define valid sector options
    VALID_SECTORS = {"p", "p-i", "p-t", "p-i-t"}
    # Validate sector value
    if snakemake.wildcards.sector not in VALID_SECTORS:
        raise ValueError(
            f"Wrong value for sector: '{snakemake.wildcards.sector}'. "
            f"Please choose one from following options: {', '.join(VALID_SECTORS)}"
        )
    # Establishing network
    n = pypsa.Network(name=f"network_{snakemake.wildcards.sector}_{selected_year}")

    snapshots = snakemake.params.snapshots
    resolution = pd.date_range(
        start=snapshots["start"],
        end=snapshots["end"],
        freq="h",
        inclusive=snapshots["inclusive"],
    )
    n.set_snapshots(resolution)
    # executing each add functions
    c = AddBaseNetwork(network=n, year=selected_year)
    c.add_buses(bus_dir=snakemake.input.elec_buses)
    c.add_atmosphere()
    c.add_interconnectors()
    c.add_fuel_supplies()
    c.add_storage_capacity(storage_capacity_dir=snakemake.input.elec_storage_capacity)
    c.add_storage_energy(storage_energy_dir=snakemake.input.elec_storage_energy)
    c.add_load(load_dir=snakemake.input.elec_loads)
    c.add_generators(plant_dir=snakemake.input.elec_power_generators)
    c.add_links(links_dir=snakemake.input.elec_power_links)
    print("Finish adding power sector for base year")
    if "i" in snakemake.wildcards.sector:
        c.add_buses(bus_dir=snakemake.input.ind_buses)
        c.add_load(load_dir=snakemake.input.ind_loads)
        c.add_storage_capacity(
            storage_capacity_dir=snakemake.input.ind_storage_capacity
        )
        c.add_storage_energy(storage_energy_dir=snakemake.input.ind_storage_energy)
        c.add_generators(plant_dir=snakemake.input.ind_generators)
        c.add_links(links_dir=snakemake.input.ind_heat_links)
        c.add_links(links_dir=snakemake.input.ind_fuel_conversion_links)
        c.add_links(links_dir=snakemake.input.ind_dac_links)
        print("Finish adding industry sector for base year")
    if "t" in snakemake.wildcards.sector:
        c.add_buses(bus_dir=snakemake.input.tra_buses)
        c.add_load(load_dir=snakemake.input.tra_loads)
        c.add_ev_chargers()
        c.add_ev_storage()
        print("Finish adding transport sector for base year")
    c.add_load_shedding()
    # c.add_noisy_cost() # noqa: E800
    c.remove_abundant_components()
    c.add_temporal_aggregation(method=snakemake.params.method)
    c.add_co2_option()
    c.add_network_carriers()
    c.network.export_to_netcdf(snakemake.output[0])
