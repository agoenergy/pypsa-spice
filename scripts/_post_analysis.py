# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Post analysis functions and classes for processing model results."""

import numpy as np
import pandas as pd
import pypsa
from _helpers import (
    generation_type_mapping,
    plot_table,
    scaling_conversion,
)
from pypsa.descriptors import get_switchable_as_dense as get_as_dense

# ================================= NAMING CONVENTIONS =================================
#
# sector_KPI_geo_resolution
#
# Example --> pow_genbyfuel_nat_hourly
#             Power sector generation by fuel at national level and hourly resolution
#
# sector            --> [pow, ind, bdg, tra]
# geo               --> [tot, nat, zone]
# time_resolution   --> [y, m, d, h]
#
# ======================================= CLASSES ======================================


def standard_plot(
    file_name: str,
    save_path: str = None,
    y: str = "value",
    table: pd.DataFrame = None,
    number_of_days: int = 0,
):
    """Generate and save standard plots from a PyPSA networks.

    Parameters
    ----------
    file_name : str
        Name of the file to be saved.
    SavePath : str, optional
        Path to save the charts, by default None
    y : str, optional
        Column name used for y-axis, by default "value"
    Table : pd.DataFrame, optional
        _description_, by default None
    number_of_days : int, optional
        selected days range for hourly plot, by default 0
    """
    x = ["snapshot" if "hourly" in file_name else "year"][0]
    plot_table(
        file_dir=file_name,
        save_dir=save_path,
        x=x,
        y=y,
        table=table,
        days_range=number_of_days,
    )


class Plots:
    """Generate and save standard plots from a PyPSA networks."""

    def __init__(self, network_dict):
        self.network_dict = network_dict
        self.output_tables = OutputTables()

    # ===================================== ENERGY =====================================
    def ene_dmd_by_carrier_yearly_plot(self):
        """Plot annual energy demand by country and carrier."""
        file_name = self.output_tables.ene_dmd_by_carrier_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.ene_dmd_by_carrier_yearly()
        )

    # ===================================== POWER ======================================
    def pow_cap_by_type_yearly_plot(self):
        """Plot annual power generation capacity by country and type."""
        file_name = self.output_tables.pow_cap_by_type_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.pow_cap_by_type_yearly()
        )

    def pow_gen_by_type_yearly_plot(self):
        """Plot annual power generation by country and year."""
        file_name = self.output_tables.pow_gen_by_type_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.pow_gen_by_type_yearly()
        )

    def pow_gen_by_category_share_yearly_plot(self):
        """Plot annual power generation share by country and carrier."""
        file_name = self.output_tables.pow_gen_by_category_share_yearly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_gen_by_category_share_yearly(),
        )

    def pow_capex_by_type_yearly_plot(self):
        """Plot annual CAPEX in the power sector by country and type."""
        file_name = self.output_tables.pow_capex_by_type_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.pow_capex_by_type_yearly()
        )

    def pow_gen_by_type_hourly_plot(self, year: int, number_of_days: int = 0):
        """Plot hourly power generation (i.e., dispatch) by country and type."""
        file_name = self.output_tables.pow_gen_by_type_hourly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_gen_by_type_hourly(year=year),
            number_of_days=number_of_days,
        )

    def pow_nodal_flow_hourly_plot(self, year: int, number_of_days: int = 0):
        """Plot hourly nodal power flow between regions."""
        file_name = self.output_tables.pow_nodal_flow_hourly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_nodal_flow_hourly(year=year),
            number_of_days=number_of_days,
        )

    def pow_battery_flows_by_region_hourly_plot(
        self, year: int, number_of_days: int = 0
    ):
        """Plot hourly battery flow by country."""
        file_name = self.output_tables.pow_battery_flows_by_region_hourly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_battery_flows_by_region_hourly(year=year),
            number_of_days=number_of_days,
        )

    def pow_hphs_flows_by_region_hourly_plot(self, year: int, number_of_days: int = 0):
        """Plot hourly hydro pumped-storage flow by country."""
        file_name = self.output_tables.pow_hphs_flows_by_region_hourly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_hphs_flows_by_region_hourly(year=year),
            number_of_days=number_of_days,
        )

    def pow_marginal_price_by_region_hourly_plot(
        self, year: int, number_of_days: int = 0
    ):
        """Plot hourly marginal price of electricity by country and region."""
        file_name = self.output_tables.pow_marginal_price_by_region_hourly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_marginal_price_by_region_hourly(year=year),
            number_of_days=number_of_days,
        )

    def pow_emi_by_carrier_yearly_plot(self):
        """Plot total annual emissions by country, sector, and carrier."""
        file_name = self.output_tables.pow_emi_by_carrier_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.pow_emi_by_carrier_yearly()
        )

    def pow_intercap_by_region_yearly_plot(self):
        """Plot interconnection capacities across regions by country and year."""
        file_name = self.output_tables.pow_intercap_by_region_yearly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.pow_intercap_by_region_yearly(),
        )

    def pow_cap_by_region_yearly_plot(self):
        """Plot annual power generation capacity by country and region."""
        file_name = self.output_tables.pow_cap_by_region_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.pow_cap_by_region_yearly()
        )

    # ==================================== INDUSTRY ====================================
    def ind_cap_by_carrier_by_region_yearly_plot(self):
        """Plot annual industrial cap. by country, carrier, region and heat type."""
        file_name = self.output_tables.ind_cap_by_carrier_by_region_yearly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.ind_cap_by_carrier_by_region_yearly(),
        )

    def ind_gen_by_carrier_by_region_yearly_plot(self):
        """Plot annual industrial gen. by country, carrier, region and heat type."""
        file_name = self.output_tables.ind_gen_by_carrier_by_region_yearly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.ind_gen_by_carrier_by_region_yearly(),
        )

    def ind_emi_by_carrier_yearly_plot(self):
        """Plot annual industrial CO2 emissions by country and carrier."""
        file_name = self.output_tables.ind_emi_by_carrier_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.ind_emi_by_carrier_yearly()
        )

    def ind_cap_by_type_by_carrier_yearly_plot(self):
        """Plot annual industrial cap. by country, type, carrier and heat type."""
        file_name = self.output_tables.ind_cap_by_type_by_carrier_yearly.__name__
        standard_plot(
            file_name=file_name,
            table=self.output_tables.ind_cap_by_type_by_carrier_yearly(),
        )

    # =================================== TRANSPORT ====================================
    def tran_capex_by_type_yearly_plot(self):
        """Plot annual transport CAPEX by country and type."""
        file_name = self.output_tables.tran_capex_by_type_yearly.__name__
        standard_plot(
            file_name=file_name, table=self.output_tables.tran_capex_by_type_yearly()
        )


class OutputTables(Plots):
    """Create output tables for post analysis."""

    def __init__(
        self,
        network_list: list,
        config: dict,
    ):
        """Initialize OutputTables class.

        Parameters
        ----------
        network_list : list
            list of networks
        config : dict
            dictionary containing network elements
        """
        self.networks = network_list
        self.network_dict = self.get_network_dict()
        self.config = config
        self.countries = self.config["base_configs"]["countries"]

    def get_network_dict(self) -> dict:
        """Get a dictionary of networks from a list of file paths.

        Returns
        -------
        dict
            A dictionary where keys are years (int) and values are PyPSA Network objects
        """
        network_dict = {}
        for path in self.networks:
            network = pypsa.Network(
                path,
            )
            # Extract year from the file path assuming the format contains '_YEAR.nc'
            year = path.split("_")[-1].replace(".nc", "")
            network_dict[int(year)] = network
        return network_dict

    # ===================================== ENERGY =====================================
    def ene_dmd_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual energy demand by country and carrier.

        << Unit: TWh >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            demand = (
                n.loads_t.p_set.multiply(n.snapshot_weightings.generators, axis=0)
                .sum()
                .groupby([n.loads.country, n.loads.carrier])
                .sum()
                .to_frame()
            )
            demand.columns = ["value"]
            demand["year"] = year
            final_df = pd.concat(
                [final_df, demand], axis=0
            )  # return demand for all year
        final_df.index.names = ["country", "carrier"]
        final_df = scaling_conversion(df=final_df, scaling_number=1e6, decimals=1)

        return final_df

    def ene_fom_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual fixed O&M costs by type, carrier, and country.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology(=type), carrier, year)
            and column (value)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            # Filter out the rows where the index contains specific substrings
            excluded_substrings = ["_SUPPLY", "LSLO", "_STORE"]
            capex_val = (
                n.statistics.capex(
                    cost_attribute="fom_cost", groupby=["type", "carrier", "country"]
                )
                .reset_index(drop=True, level=0)
                .to_frame(name="value")
            )
            keep_types = [
                x
                for x in capex_val.index.levels[0]
                if not any(substring in x for substring in excluded_substrings)
            ]
            capex_val = capex_val.loc[
                capex_val.index.get_level_values("type").isin(keep_types)
            ]
            capex_val["year"] = year
            final_df = pd.concat(
                [final_df, capex_val], axis=0
            )  # return cost mix for a all years
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.reset_index().rename(columns={"type": "technology"})
        final_df = (
            final_df.groupby(["country", "technology", "carrier", "year"])
            .sum()
            .sort_values(by=["year", "technology", "carrier"])[["value"]]
        )  # groupby type and ignore regional information

        return final_df

    def ene_capex_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual CAPEX by type, carrier, and country.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology(=type), carrier, year)
            and column (value)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            # Filter out the rows where the index contains specific substrings
            excluded_substrings = ["_SUPPLY", "LSLO", "_STORE"]
            capex_val = (
                n.statistics.capex(groupby=["type", "carrier", "country"])
                .sub(
                    n.statistics.capex(
                        cost_attribute="fom_cost",
                        groupby=["type", "carrier", "country"],
                    ),
                    fill_value=0,
                )
                .reset_index(drop=True, level=0)
                .to_frame(name="value")
            )
            keep_types = [
                x
                for x in capex_val.index.levels[0]
                if not any(substring in x for substring in excluded_substrings)
            ]
            capex_val = capex_val.loc[
                capex_val.index.get_level_values("type").isin(keep_types)
            ]  # return cost mix for a single year in Million CURRENCY
            capex_val["year"] = year
            final_df = pd.concat(
                [final_df, capex_val], axis=0
            )  # return cost mix for a all years
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.reset_index().rename(columns={"type": "technology"})
        final_df = (
            final_df.groupby(["country", "technology", "carrier", "year"])
            .sum()
            .sort_values(by=["year", "technology", "carrier"])[["value"]]
        )  # groupby type and ignore regional information

        return final_df

    def ene_opex_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual OPEX by type, carrier, and country.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology(=type), carrier, year)
            and column (value)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            # Filter out the rows where the index contains specific substrings
            excluded_substrings = ["LSLO", "_STORE"]
            opex_val = (
                n.statistics.opex(groupby=["type", "carrier", "country"])
                .add(
                    n.statistics.capex(
                        cost_attribute="fom_cost",
                        groupby=["type", "carrier", "country"],
                    ),
                    fill_value=0,
                )
                .reset_index(drop=True, level=0)
                .to_frame(name="value")
            )
            keep_types = [
                x
                for x in opex_val.index.levels[0]
                if not any(substring in x for substring in excluded_substrings)
            ]
            opex_val = opex_val.loc[
                opex_val.index.get_level_values("type").isin(keep_types)
            ]  # return cost mix for a single year in Million CURRENCY
            opex_val["year"] = year
            final_df = pd.concat(
                [final_df, opex_val], axis=0
            )  # return cost mix for a all years
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.reset_index().rename(columns={"type": "technology"})
        final_df = (
            final_df.groupby(["country", "technology", "carrier", "year"])
            .sum()
            .sort_values(by=["year", "technology", "carrier"])[["value"]]
        )  # groupby type and ignore regional information

        return final_df

    def ene_emi_by_carrier_by_sector_yearly(self) -> pd.DataFrame:
        """Calculate total annual emissions by country, sector, and carrier.

        << Unit: Mt_CO2 >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, sector, carrier, year) and
            column (value)
        """
        df_all = pd.DataFrame()
        for country in self.countries:
            pow_emi_df = (
                self.pow_emi_by_carrier_yearly()
                .loc[country]
                .groupby(["carrier", "year"])
                .sum()
                .reset_index()
                .assign(sector="power")
            )
            if "i" in self.config["base_configs"]["sector"][0]:
                ind_emi_df = (
                    self.ind_emi_by_carrier_yearly()
                    .loc[country]
                    .groupby(["carrier", "year"])
                    .sum()
                    .reset_index()
                    .assign(sector="industry")
                )
            else:
                ind_emi_df = pd.DataFrame()  # empty dataframe if no industry
            country_df = pd.concat([pow_emi_df, ind_emi_df], axis=0)
            country_df["country"] = country  # add country column
            df_all = pd.concat([df_all, country_df], axis=0)

        if df_all.empty:
            # return an empty dataframe
            return pd.DataFrame(
                columns=["country", "sector", "carrier", "year", "value"]
            )

        # Set index and return the dataframe when not empty
        df_all = df_all.set_index(["country", "sector", "carrier", "year"])
        return df_all

    def ene_emi_cost_by_sector_yearly(self) -> pd.DataFrame:
        """Calculate total annual emission costs by country and sector.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, sector, year) and column (value)
        """
        df_all = pd.DataFrame()
        for country in self.countries:
            if self.config["co2_management"][country]["option"] != "co2_price":
                continue  # skip if no CO2 price is applied
            pow_emi_cost_df = (
                self.pow_emi_cost_by_carrier_yearly()
                .loc[country]
                .groupby("year")
                .sum()
                .reset_index()
                .assign(sector="power")
            )
            if "i" in self.config["base_configs"]["sector"][0]:
                ind_emi_cost_df = (
                    self.ind_emi_cost_by_carrier_yearly()
                    .loc[country]
                    .groupby("year")
                    .sum()
                    .reset_index()
                    .assign(sector="industry")
                )
            else:
                ind_emi_cost_df = pd.DataFrame()  # empty dataframe if no industry
            country_df = pd.concat([pow_emi_cost_df, ind_emi_cost_df], axis=0)
            country_df["country"] = country  # add country column
            df_all = pd.concat([df_all, country_df], axis=0)

        # return an empty dataframe
        if df_all.empty:
            return pd.DataFrame(columns=["country", "sector", "year", "value"])

        # Set index and return the dataframe when not empty
        df_all = df_all.set_index(["country", "sector", "year"])
        return df_all

    def ene_costs_opex_capex_yearly(self) -> pd.DataFrame:
        """Calculate total annual system costs broken down by CAPEX and OPEX.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, finance, year) and column (value)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            opex_val = (
                self.ene_opex_by_type_by_carrier_yearly()
                .loc[:, :, :, year, :]
                .groupby("country")
                .sum()
            )
            opex_val["finance"] = "opex"
            capex_val = (
                self.ene_capex_by_type_by_carrier_yearly()
                .loc[:, :, :, year, :]
                .groupby("country")
                .sum()
            )
            capex_val["finance"] = "capex"
            cost_yr = pd.concat([opex_val, capex_val], axis=0)
            cost_yr = cost_yr.set_index(["finance"], append=True)
            cost_yr["year"] = year
            final_df = pd.concat([final_df, cost_yr], axis=0)

        return final_df

    def ene_avg_fuel_costs_fuel_yearly(self) -> pd.DataFrame:
        """Calculate average fuel costs incl. synthetic fuels by country and carrier.

        << Unit: CURRENCY/MWh_th >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier, year) and column (value)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            # marginal prices of the fuel buses
            bus_id = n.buses[
                ~n.buses.carrier.isin(
                    ["Electricity", "Co2stor", "CO2", "Low_Heat", "High_Heat"]
                )
            ].index
            mar_prices = (
                n.buses_t.marginal_price[bus_id]
                .T.groupby([n.buses.country, n.buses.carrier])
                .mean()
                .T
            )
            # energy balance of the fuel buses
            eb = n.statistics.energy_balance(
                groupby=n.statistics.groupers.get_bus_and_carrier_and_bus_carrier,
                aggregate_time=False,
            ).reset_index()
            eb["country"] = eb["bus"].apply(lambda x: n.buses.country[x])
            eb = eb.loc[eb.bus.isin(bus_id)].round()
            # taking only the links's energy balance
            eb = eb.loc[eb.component == "Link"].drop(
                columns=["component", "bus", "bus_carrier"]
            )
            eb = eb.set_index(["country", "carrier"]).T
            # selecting only the exports from the link
            eb_links_pos = pd.DataFrame()
            for col in eb.columns:
                eb_links_pos[col] = eb[col].apply(lambda x: 0 if x > 0 else x)
            average_fuel = (
                (
                    (mar_prices[eb_links_pos.columns]).mul(eb_links_pos).sum()
                    / eb_links_pos.sum()
                )
                .round(2)
                .to_frame(name="value")
            )
            average_fuel["year"] = year
            final_df = pd.concat([final_df, average_fuel], axis=0)
        final_df = final_df.reset_index().set_index(["country", "carrier", "year"])

        return final_df

    def ene_pe_import_local_ei_yearly(self) -> pd.DataFrame:
        """Calculate energy independence fraction.

        Calculation is based on primary energy imports and local primary energy
        generation.

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (year, type, value, country)
        """
        df_all = pd.DataFrame()
        for country in self.countries:
            if country not in self.config.get("custom_constraints", {}):
                continue
            if not self.config["custom_constraints"][country].get(
                "energy_independence", False
            ):
                continue
            final_df = pd.DataFrame()
            for year in self.network_dict:
                df_year = pd.DataFrame()
                n = self.network_dict[year]
                ff_gen = n.generators[
                    (n.generators.country == country)
                    & (n.generators.type.str.contains("SUPPLY"))
                ].index
                ff_imp_gen = ff_gen[ff_gen.str.contains("-IMP")]
                ff_loc_gen = [ele for ele in ff_gen if ele not in ff_imp_gen]

                # primary energy import
                pe_imp = (
                    n.generators_t.p[ff_imp_gen]
                    .multiply(n.snapshot_weightings.generators, axis=0)
                    .sum()
                    .sum()
                )
                # local energy generators
                pe_loc = (
                    n.generators_t.p[ff_loc_gen]
                    .multiply(n.snapshot_weightings.generators, axis=0)
                    .sum()
                    .sum()
                )

                # electricity generators from renewable sources (res)
                pe_conv_frac = self.config["custom_constraints"][country][
                    "energy_independence"
                ]["pe_conv_fraction"]
                for res in ["Solar", "Wind", "Geothermal", "Water"]:
                    res_gen_loc = n.generators[
                        (n.generators.carrier == res)
                        & (n.generators.country == country)
                        & ~(n.generators.index.str.contains("-IMP"))
                    ].index
                    if not res_gen_loc.empty:
                        res_gen_pe = (
                            n.generators_t.p[res_gen_loc]
                            .multiply(n.snapshot_weightings.generators, axis=0)
                            .sum()
                            .sum()
                            * pe_conv_frac[res]
                        )
                        pe_loc = pe_loc + res_gen_pe
                # If HDAM implemented as storage unit

                hdam_str_loc = n.storage_units[
                    (n.storage_units.type == "HDAM")
                    & (n.storage_units.country == country)
                    & ~(n.storage_units.index.str.contains("-IMP"))
                ].index
                hdam_str_imp = n.storage_units[
                    (n.storage_units.type == "HDAM")
                    & (n.storage_units.country == country)
                    & (n.storage_units.index.str.contains("-IMP"))
                ].index

                if not hdam_str_loc.empty:  # local HDAM
                    hdam_loc_gen_pe = (
                        n.storage_units_t.p_dispatch[hdam_str_loc]
                        .multiply(n.snapshot_weightings.stores, axis=0)
                        .sum()
                        .sum()
                        * pe_conv_frac["Water"]
                    )
                    pe_loc = pe_loc + hdam_loc_gen_pe

                if not hdam_str_imp.empty:  # imported HDAM
                    hdam_imp_gen_pe = (
                        n.storage_units_t.p_dispatch[hdam_str_loc]
                        .multiply(n.snapshot_weightings.stores, axis=0)
                        .sum()
                        .sum()
                        * pe_conv_frac["Water"]
                    )
                    pe_imp = pe_imp + hdam_imp_gen_pe
                df_year = pd.concat(
                    [df_year, pd.DataFrame([year, "import", (pe_imp / 1e6).round(2)])],
                    axis=1,
                )
                df_year = pd.concat(
                    [df_year, pd.DataFrame([year, "local", (pe_loc / 1e6).round(2)])],
                    axis=1,
                )
                df_year = pd.concat(
                    [
                        df_year,
                        pd.DataFrame(
                            [year, "ei_frac", (1 - pe_imp / (pe_imp + pe_loc)).round(3)]
                        ),
                    ],
                    axis=1,
                )
                df_year = df_year.T
                final_df = pd.concat([final_df, df_year], axis=0)
            final_df.columns = ["year", "type", "value"]
            final_df["country"] = country
            df_all = pd.concat([df_all, final_df], axis=0)
            final_df.set_index(["country", "type"], inplace=True)
            return final_df

    # ====================================== POWER =====================================
    def pow_cap_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual power generation capacity by country and type.

        << Unit: GW >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cap_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                pow_supply = df[
                    (df[bus_name].str.contains("HVELEC"))
                    | (df[bus_name].str.contains("LVELEC"))
                ]
                if c.name == "Link":
                    pow_supply_cap = (
                        (pow_supply.p_nom_opt * pow_supply.efficiency)
                        .groupby([pow_supply.country, pow_supply.type])
                        .sum()
                        .to_frame()
                    )
                else:
                    pow_supply_cap = (
                        pow_supply.p_nom_opt.groupby(
                            [pow_supply.country, pow_supply.type]
                        )
                        .sum()
                        .to_frame()
                    )
                pow_supply_cap.columns = ["value"]
                cap_mix = pd.concat(
                    [cap_mix, pow_supply_cap], axis=0
                )  # return cap mix for a single year
            cap_mix["year"] = year
            final_df = pd.concat(
                [final_df, cap_mix], axis=0
            )  # return cap mix for all years
        final_df.index.names = ["country", "technology"]
        final_df = final_df.fillna(0)
        final_df = final_df[
            ~final_df.index.get_level_values("technology").isin(["ITCN", "LSLO"])
        ]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e3,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    def pow_gen_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual power generation by country and type.

        << Unit: TWh >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            gen_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                pow_supply = c.df[
                    (c.df[bus_name].str.contains("HVELEC"))
                    | (c.df[bus_name].str.contains("LVELEC"))
                ].index
                if c.name == "Generator":
                    c_gen = (
                        c.pnl.p[pow_supply]
                        .multiply(n.snapshot_weightings.generators, axis=0)
                        .sum()
                        .multiply(c.df.sign)
                        .groupby([c.df.country, c.df.type])
                        .sum()
                        .to_frame()
                    )
                elif c.name == "StorageUnit":
                    c_gen = (
                        c.pnl.p[[x for x in c.pnl.p.columns if x in pow_supply]]
                        .multiply(n.snapshot_weightings.stores, axis=0)
                        .sum()
                        .multiply(c.df.sign)
                        .groupby([c.df.country, c.df.type])
                        .sum()
                        .to_frame()
                    )
                else:
                    c_gen = (
                        (
                            -c.pnl.p1[pow_supply]
                            .multiply(n.snapshot_weightings.generators, axis=0)
                            .sum()
                        )
                        .groupby([c.df.country, c.df.type])
                        .sum()
                        .to_frame()
                    )
                c_gen.columns = ["value"]
                gen_mix = pd.concat(
                    [gen_mix, c_gen], axis=0
                )  # return gen mix for a single year
            gen_mix["year"] = year
            final_df = pd.concat(
                [final_df, gen_mix], axis=0
            )  # return gen mix for all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology"]
        final_df = final_df[
            ~final_df.index.get_level_values("technology").isin(
                [
                    "ITCN",
                    "LSLO",
                    "EVST_PRV",
                    "EVST_PUB",
                    "BATS_DISCHARGE",
                    "HHBS_DISCHARGE",
                    "HPHS_DISCHARGE",
                ]
            )
        ]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    def pow_gen_by_category_share_yearly(self) -> pd.DataFrame:
        """Calculate annual power generation share by country and carrier.

        << Unit: % >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, year) and column (value)
        """
        final_df = self.pow_gen_by_type_yearly().reset_index()
        # import and use functions from _helpers.py
        # for easier adjustment in the long run
        final_df["carrier"] = final_df["technology"].map(generation_type_mapping)
        final_df = (
            final_df.drop("technology", axis=1)
            .groupby(["country", "year", "carrier"])
            .sum()["value"]
            .unstack()
        )
        total = total = final_df.sum(axis=1)
        for columns in final_df.columns:
            final_df[columns] = final_df[columns] / total * 100
        final_df = final_df.fillna(0)
        final_df = pd.melt(
            final_df.reset_index(),
            id_vars=["country", "year"],
            var_name="technology",
            value_name="value",
        )
        final_df = final_df.set_index(["country", "technology", "year"])
        return final_df

    def pow_capex_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual CAPEX in the power sector by country and type.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in ["Link", "Generator", "StorageUnit", "Store"]:
                bus_name = "bus1" if c == "Link" else "bus"
                included_substrings = ["HVELEC", "LVELEC", "HHBSN", "BATSN", "HPHSN"]
                if c == "Store":
                    pow_cost = n.statistics.capex(
                        comps=c, groupby=[bus_name, "country", "type"]
                    )
                else:
                    pow_cost = n.statistics.capex(
                        comps=c, groupby=[bus_name, "country", "type"]
                    ).sub(
                        n.statistics.capex(
                            cost_attribute="fom_cost",
                            comps=c,
                            groupby=[bus_name, "country", "type"],
                        ),
                        fill_value=0,
                    )
                if not pow_cost.empty:
                    list_of_bus = (
                        pow_cost.index.get_level_values(bus_name).unique().tolist()
                    )
                    elec_bus = [
                        x
                        for x in list_of_bus
                        if any(substring in x for substring in included_substrings)
                    ]
                    if len(elec_bus) > 0:
                        pow_cost = (
                            pow_cost.loc[elec_bus]
                            .groupby(["country", "type"])
                            .sum()
                            .to_frame(name="value")
                        )
                        cost_mix = pd.concat(
                            [cost_mix, pow_cost], axis=0
                        )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology"]
        final_df = final_df.loc[
            ~final_df.index.get_level_values("technology").isin(["LSLO"])
        ]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "year"]).sum()

        return final_df

    def pow_fom_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual fixed O&M costs in the power sector by country and type.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in ["Link", "Generator", "StorageUnit"]:
                bus_name = "bus1" if c == "Link" else "bus"
                included_substrings = ["HVELEC", "LVELEC", "HHBSN", "BATSN", "HPHSN"]
                pow_cost = n.statistics.capex(
                    cost_attribute="fom_cost",
                    comps=c,
                    groupby=[bus_name, "country", "type"],
                )
                if not pow_cost.empty:
                    list_of_bus = (
                        pow_cost.index.get_level_values(bus_name).unique().tolist()
                    )
                    elec_bus = [
                        x
                        for x in list_of_bus
                        if any(substring in x for substring in included_substrings)
                    ]
                    if len(elec_bus) > 0:
                        pow_cost = (
                            pow_cost.loc[elec_bus]
                            .groupby(["country", "type"])
                            .sum()
                            .to_frame(name="value")
                        )
                        cost_mix = pd.concat(
                            [cost_mix, pow_cost], axis=0
                        )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology"]
        final_df = final_df.loc[
            ~final_df.index.get_level_values("technology").isin(["LSLO"])
        ]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "year"]).sum()

        return final_df

    def pow_opex_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual OPEX in the power sector by country and type.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in ["Link", "Generator", "StorageUnit", "Store"]:
                bus_name = "bus1" if c == "Link" else "bus"
                included_substrings = ["HVELEC", "LVELEC"]
                # new opex = opex + fom_cost + fuel_cost
                if c == "Store":
                    pow_cost = n.statistics.opex(
                        comps=c, groupby=[bus_name, "country", "type"]
                    )
                else:
                    pow_cost = n.statistics.opex(
                        comps=c, groupby=[bus_name, "country", "type"]
                    ).add(
                        n.statistics.capex(
                            cost_attribute="fom_cost",
                            comps=c,
                            groupby=[bus_name, "country", "type"],
                        ),
                        fill_value=0,
                    )
                if c == "Link":
                    # Create a dictionary comprehension
                    # to map links to their prices from primary bus0
                    fuel_price_dict = {
                        link: n.buses_t.marginal_price[n.links.loc[link, "bus0"]]
                        for link in n.links_t.p0.columns
                    }
                    # Convert the dictionary into a DataFrame
                    fuel_price_df = pd.DataFrame.from_dict(fuel_price_dict)
                    fuel_cost = (
                        n.links_t.p0.multiply(n.snapshot_weightings.objective, axis=0)
                        .multiply(fuel_price_df)
                        .T.groupby([n.links["bus1"], n.links.country, n.links.type])
                        .sum()
                        .T.sum()
                    )
                    # Set fuel cost of storage and interconnection link
                    # as 0 to avoid double counting
                    target_substrings = {"_STORE", "_DISCHARGE", "ITCN"}
                    targeted_type = [
                        any(substring in level for substring in target_substrings)
                        for level in fuel_cost.index.get_level_values("type")
                    ]
                    fuel_cost.loc[targeted_type] = 0
                    # add fuel cost to opex
                    pow_cost = pow_cost.add(fuel_cost)
                if not pow_cost.empty:
                    list_of_bus = (
                        pow_cost.index.get_level_values(bus_name).unique().tolist()
                    )
                    elec_bus = [
                        x
                        for x in list_of_bus
                        if any(substring in x for substring in included_substrings)
                    ]
                    if len(elec_bus) > 0:
                        pow_cost = (
                            pow_cost.loc[elec_bus]
                            .groupby(["country", "type"])
                            .sum()
                            .to_frame(name="value")
                        )
                        cost_mix = pd.concat(
                            [cost_mix, pow_cost], axis=0
                        )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology"]
        final_df = final_df.loc[
            ~final_df.index.get_level_values("technology").isin(["LSLO"])
        ]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "year"]).sum()

        return final_df

    def pow_overnight_inv_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual overnight investment in power sector by country and type.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in n.iterate_components(
                ["Link", "Generator", "StorageUnit", "Store"]
            ):
                attr = "e" if c.name == "Store" else "p"
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                df = df.loc[df[attr + "_nom_extendable"]]
                pow_supply = df[
                    (df[bus_name].str.contains("HVELEC"))
                    | (df[bus_name].str.contains("LVELEC"))
                ]
                inv_cost = (
                    (pow_supply[attr + "_nom_opt"] * pow_supply.inv_cost)
                    .groupby([pow_supply.country, pow_supply.type])
                    .sum()
                ).to_frame()
                inv_cost.columns = ["value"]
                cost_mix = pd.concat(
                    [cost_mix, inv_cost], axis=0
                )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df.index.names = ["country", "technology"]
        final_df = final_df.fillna(0)
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "year"]).sum()

        return final_df

    def pow_gen_by_type_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly power generation (i.e., dispatch) by country and type.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, snapshot) and
            column (value)
        """
        n = self.network_dict[year]
        all_countries_df = pd.DataFrame()
        for country in self.countries:
            country_df = pd.DataFrame()
            for comp in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_col = "bus1" if comp.name == "Link" else "bus"
                # Pre-compute the mask for bus_col contains any of the relevant
                # substrings
                bus_mask = comp.df[bus_col].str.contains(
                    "HVELEC|LVELEC|BATSN|HPHSN", regex=True
                )
                supply_idx = comp.df[(comp.df.country == country) & bus_mask].index
                if comp.name in ("Generator", "StorageUnit"):
                    gen_df = (
                        comp.pnl.p[[x for x in comp.pnl.p.columns if x in supply_idx]]
                        .T.groupby([comp.df.type])
                        .sum()
                        .T
                    )
                else:
                    gen_df = (
                        (
                            -comp.pnl.p1[
                                [x for x in comp.pnl.p1.columns if x in supply_idx]
                            ]
                        )
                        .T.groupby(comp.df.type)
                        .sum()
                        .T
                    )
                country_df = pd.concat([country_df, gen_df], axis=1)
            country_df = country_df.loc[
                :, ~country_df.columns.isin(["ITCN", "EVST_PRV", "EVST_PUB"])
            ]
            country_df = (country_df.loc[:, ~(country_df == 0).all(axis=0)]).round(1)
            country_df = country_df.reset_index()
            country_df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(country_df), freq=f"{nth_hour}h"
            )

            country_df = pd.melt(
                country_df, id_vars=["snapshot"], var_name="technology"
            ).set_index("snapshot")
            country_df["country"] = country
            # reverse value for all store links
            country_df.loc[country_df.technology.str.contains("_STORE"), "value"] *= -1
            country_df = country_df.set_index(["country", "technology"], append=True)
            all_countries_df = pd.concat([all_countries_df, country_df], axis=0)
        return all_countries_df

    def pow_marginal_price_by_region_hourly(
        self, year: int, nth_hour: int
    ) -> pd.DataFrame:
        """Calculate hourly marginal price of electricity by country and region.

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, region, snapshot) and column (value)
        """
        all_countries_df = pd.DataFrame()
        for country in self.countries:
            n = self.network_dict[year]
            bus_id = n.buses[
                (n.buses.country == country) & (n.buses.index.str.contains("HVELEC"))
            ].index
            marginal_price = n.buses_t.marginal_price[bus_id]
            marginal_price.columns = marginal_price.columns.str.split("_").str[1]
            final_df = marginal_price.reset_index()
            final_df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(final_df), freq=f"{nth_hour}h"
            )
            final_df = pd.melt(
                final_df, id_vars=["snapshot"], var_name="region"
            ).set_index("snapshot")
            final_df["country"] = country
            final_df = final_df.set_index(["country", "region"], append=True)
            all_countries_df = pd.concat([all_countries_df, final_df], axis=0)
        return all_countries_df

    def pow_nodal_flow_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly nodal power flow between regions.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (flow, snapshot) and column (value)
        """
        n = self.network_dict[year]
        interlinks = n.links[
            (n.links.type == "ITCN") & (n.links.bus1.str.contains("HVELEC"))
        ].index
        flows = (
            n.links_t.p0[interlinks]
            .T.groupby(
                [
                    n.links.bus0,
                    n.links.bus1,
                ],
            )
            .sum()
        )
        number_of_flows = len(flows)
        flows = flows.reset_index()
        flows[["bus0", "bus1"]] = flows[["bus0", "bus1"]].replace(
            "_HVELEC", "", regex=True
        )
        flows["flow"] = flows["bus0"] + "_to_" + flows["bus1"]
        flows = flows.drop(["bus0", "bus1"], axis=1).set_index("flow")
        flows = pd.melt(
            flows, var_name="snapshot", value_name="value", ignore_index=False
        ).reset_index()
        flows = flows.sort_values(by=["flow", "snapshot"]).set_index("flow")
        time_range = pd.date_range(
            start=f"{year}-01-01",
            periods=int(len(flows) / number_of_flows),
            freq=f"{nth_hour}h",
        )
        flows["snapshot"] = np.tile(time_range, number_of_flows)

        return flows

    def pow_gen_by_type_region_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly power dispatch by country and region.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, region, type_and_flow, snapshot)
            and column (value)
        """
        n = self.network_dict[year]
        all_countries_dispatch = pd.DataFrame()
        for country in self.countries:
            region_generation_dispatch = pd.DataFrame()
            for component in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_column = "bus1" if component.name == "Link" else "bus"
                # Pre-compute the mask for bus_column contains any of the relevant
                # substrings
                is_electric_bus = component.df[bus_column].str.contains(
                    "HVELEC|LVELEC", regex=True
                )
                supply_indices = component.df[
                    (component.df.country == country) & is_electric_bus
                ].index
                if component.name in ("Generator", "StorageUnit"):
                    generation = (
                        component.pnl.p[
                            [x for x in component.pnl.p.columns if x in supply_indices]
                        ]
                        .T.groupby(
                            [
                                component.df[bus_column].str.split("_", n=2).str[1],
                                component.df.type,
                            ]
                        )
                        .sum()
                        .T
                    )
                else:
                    generation = (
                        (
                            -component.pnl.p1[
                                [
                                    x
                                    for x in component.pnl.p1.columns
                                    if x in supply_indices
                                ]
                            ]
                        )
                        .T.groupby(
                            [
                                component.df[bus_column].str.split("_", n=2).str[1],
                                component.df.type,
                            ]
                        )
                        .sum()
                        .T
                    )
                region_generation_dispatch = pd.concat(
                    [region_generation_dispatch, generation], axis=1
                )

            region_generation_dispatch = region_generation_dispatch.loc[
                :,
                ~region_generation_dispatch.columns.get_level_values("type").isin(
                    ["ITCN", "EVST-PRV", "EVST-PUB"]
                ),
            ]

            region_generation_dispatch = (
                region_generation_dispatch.loc[
                    :, ~(region_generation_dispatch == 0).all(axis=0)
                ]
            ).round(1)
            region_generation_dispatch = region_generation_dispatch.reset_index()
            region_generation_dispatch["snapshot"] = pd.date_range(
                start=f"{year}-01-01",
                periods=len(region_generation_dispatch),
                freq=f"{nth_hour}h",
            )
            region_generation_dispatch = region_generation_dispatch.set_index(
                "snapshot"
            )

            # use function from pow_nodal_flow_hourly, with edits in groupby
            interconnector_links = n.links[
                (n.links.index.str.contains(country))
                & (n.links.type == "ITCN")
                & (n.links.bus1.str.contains("HVELEC"))
            ].index
            interconnector_flows = (
                n.links_t.p0[interconnector_links]
                .T.groupby(
                    [
                        n.links.bus0.str.replace("_HVELEC", "", regex=True),
                        n.links.bus1.str.replace("_HVELEC", "", regex=True),
                    ],
                )
                .sum()
                .T
            )
            interconnector_flows = interconnector_flows.reset_index()
            interconnector_flows["snapshot"] = pd.date_range(
                start=f"{year}-01-01",
                periods=len(interconnector_flows),
                freq=f"{nth_hour}h",
            )
            interconnector_flows = interconnector_flows.set_index("snapshot")

            # flip the flows dataframe
            # this is needed since if we filter and visualize the df by region,
            # we want the ITCN to reflect on both regions, but with opposite values
            flipped_flows = -1 * interconnector_flows.swaplevel(axis=1)

            # concat flipped dataframe with the original dataframe
            all_flows = pd.concat([interconnector_flows, flipped_flows], axis=1)

            # `first` is the region where the ITCN is connected
            # `second` is the ITCN flow, and is values are flipped if viewed from the
            # other side of the ITCN
            all_flows.columns = pd.MultiIndex.from_tuples(
                [
                    (to_region, f"{from_region}_to_{to_region}")
                    for from_region, to_region in all_flows.columns
                ],
                names=["region", "flow"],
            )

            # concat generators and flows, sort index on level=0 (regions)
            combined_dispatch = pd.concat(
                [region_generation_dispatch, all_flows], axis=1
            ).sort_index(level=0, axis=1)

            # pd.melt cannot be used here since there are multiple levels. thus, we use
            # stack twice.
            combined_dispatch = (
                combined_dispatch.stack(
                    level=0, future_stack=True
                )  # stacks the outermost column level using the new implementation
                .stack(
                    future_stack=True
                )  # stacks the outermost column level using the new implementation
                .fillna(0)
            ).reset_index()
            combined_dispatch = combined_dispatch.rename(
                columns={
                    "level_1": "region",
                    "level_2": "type_and_flow",
                    0: "value",
                }
            )
            combined_dispatch["country"] = country
            combined_dispatch["region"] = combined_dispatch["region"].str.replace(
                f"{country}_", "", regex=True
            )
            combined_dispatch = combined_dispatch.set_index("snapshot")
            all_countries_dispatch = pd.concat(
                [all_countries_dispatch, combined_dispatch], axis=0
            )
        return all_countries_dispatch

    def pow_cap_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual power generation capacity by country and region.

        << Unit: GW >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (country, region, value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cap_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                pow_supply = df[
                    (
                        (df[bus_name].str.contains("HVELEC"))
                        | (df[bus_name].str.contains("LVELEC"))
                    )
                    & ~(
                        (df.type == "LSLO")
                        | (df.type == "ITCN")
                        | (df.type.str.contains("SUPPLY"))
                    )  # removes the non-capacity items
                ]
                if c.name == "Link":
                    pow_supply_cap = (
                        (pow_supply.p_nom_opt * pow_supply.efficiency)
                        .groupby(
                            [
                                pow_supply.country,
                                pow_supply[bus_name].str.split("_", n=2).str[1],
                            ]
                        )  # groupby using the region_name of the bus
                        .sum()
                        .to_frame()
                    )
                else:
                    pow_supply_cap = (
                        pow_supply.p_nom_opt.groupby(
                            [
                                pow_supply.country,
                                pow_supply[bus_name].str.split("_", n=2).str[1],
                            ]
                        )
                        .sum()
                        .to_frame()  # groupby using the region_name of the bus
                    )
                pow_supply_cap.columns = ["value"]
                cap_mix = pd.concat(
                    [cap_mix, pow_supply_cap], axis=0
                )  # return cap mix for a single year
            cap_mix = cap_mix.reset_index()
            cap_mix.columns = ["country", "region", "value"]
            cap_mix = cap_mix.groupby(
                ["country", "region"]
            ).sum()  # need to groupby-sum per region
            cap_mix["year"] = year
            final_df = pd.concat(
                [final_df, cap_mix], axis=0
            )  # return cap mix for all years
        final_df = final_df.fillna(0)
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e3,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    def pow_flh_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual full load hours by country and type.

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        all_countries_flh = pd.DataFrame()
        for country in self.countries:
            generation_df = (
                self.pow_gen_by_type_yearly().loc[country].reset_index(drop=True)
            )
            capacity_df = (
                self.pow_cap_by_type_yearly().loc[country].reset_index(drop=True)
            )
            generation_by_type_year = generation_df.set_index(["technology", "year"])
            capacity_by_type_year = capacity_df.set_index(["technology", "year"])
            flh_df = (generation_by_type_year * 1e3) / capacity_by_type_year
            flh_df = flh_df.reset_index().dropna()
            flh_df["country"] = country
            all_countries_flh = pd.concat([all_countries_flh, flh_df], axis=0)
        return all_countries_flh.set_index(["country", "technology"])

    def pow_reserve_by_type_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly reserve by country and type.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, snapshot) and columns (type, value)
        """
        network = self.network_dict[year]
        all_countries_reserve = pd.DataFrame()
        for country in self.countries:
            country_reserve = pd.DataFrame()
            country_generators = network.generators[
                network.generators.country == country
            ].index
            try:
                generator_reserve = network.generators_t.r[
                    [
                        x
                        for x in network.generators_t.r.columns
                        if x in country_generators
                    ]
                ]  # for handling when there is no reserve
            except AttributeError:
                return pd.DataFrame()
            # drop generators without defined reserve variable
            generator_reserve = (
                generator_reserve.loc[:, (~(generator_reserve == -1).all())]
                .T.groupby([network.generators.type])
                .sum()
                .T
            )

            country_links = network.links[network.links.country == country].index
            link_reserve = network.links_t.r[
                [x for x in network.links_t.r.columns if x in country_links]
            ]
            link_reserve = (
                link_reserve.loc[:, ~(link_reserve == -1).all()]
                .T.groupby([network.links.type])
                .sum()
                .T
            )
            reserve_time_series = generator_reserve.join(link_reserve)
            for col in reserve_time_series:
                if reserve_time_series[col].max() <= 1:
                    reserve_time_series = reserve_time_series.drop([col], axis=1)
            reserve_time_series["reserve_total"] = reserve_time_series.sum(axis=1)
            load_indices = network.loads[
                (network.loads.carrier == "Electricity")
                & (network.loads.country == country)
            ].index
            demand_reserve = (
                get_as_dense(network, "Load", "p_set").loc[:, load_indices].sum(axis=1)
            )
            demand_reserve = demand_reserve.to_frame(name="demand_reserve")
            all_generators_indices = network.generators[
                (network.generators.country == country)
                & (network.generators.bus.map(network.buses.carrier) == "Electricity")
                & (network.generators.carrier != "ENS")
            ].index
            vre_reserve = (
                network.generators.p_nom[all_generators_indices]
                .multiply(network.generators_t.p_max_pu)
                .sum(axis=1)
            )
            vre_reserve = vre_reserve.to_frame(name="vre_reserve")
            country_reserve = pd.concat(
                [reserve_time_series, demand_reserve, vre_reserve], axis=1
            )
            country_reserve.index = pd.date_range(
                start=f"{year}-01-01", periods=len(country_reserve), freq=f"{nth_hour}h"
            )
            country_reserve["country"] = country
            country_reserve = country_reserve.set_index("country", append=True)
            all_countries_reserve = pd.concat(
                [all_countries_reserve, country_reserve], axis=0
            )
        return all_countries_reserve

    def pow_bats_ep_ratio(self) -> pd.DataFrame:
        """Calculate battery energy to power ratio by country and year.

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, year) and column (value)
        """
        data = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            e = (
                n.stores[n.stores.type == "BATS"]["e_nom_opt"]
                .groupby([n.stores.country, n.stores.bus])
                .sum()
            )
            p = (
                n.links[n.links.type == "BATS_DISCHARGE"]["p_nom_opt"]
                .groupby([n.links.country, n.links.bus0])
                .sum()
            )
            ratio = e.div(p).to_frame(name="value")
            ratio["year"] = year
            ratio = ratio.set_index(["year"], append=True)
            data = pd.concat([data, ratio])
        return data

    def pow_bats_charging_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly battery charging profile and state of charge by country.

        << Unit: MWh >>

        Parameters
        ----------
        year : int
            _description_
        nth_hour : int
            _description_

        Returns
        -------
        pd.DataFrame
            _description_
        """
        n = self.network_dict[year]
        all_countries_df = pd.DataFrame()
        for country in self.countries:
            discharge_link_indices = (
                n.links[
                    (n.links.country == country)
                    & (n.links.type.str.contains("BATS_DISCHARGE"))
                ].index
                if "BATS_DISCHARGE" in n.links.type.unique()
                else None
            )
            store_link_indices = (
                n.links[
                    (n.links.country == country)
                    & (n.links.type.str.contains("BATS_STORE"))
                ].index
                if "BATS_STORE" in n.links.type.unique()
                else None
            )
            storageunit_indices = (
                n.storage_units[
                    (n.storage_units.country == country)
                    & (n.storage_units.type.str.contains("BATS"))
                ].index
                if "BATS" in n.storage_units.type.unique()
                else None
            )
            df = pd.DataFrame()
            if discharge_link_indices is not None and store_link_indices is not None:
                battery_power_flow = (
                    pd.concat(
                        [
                            -1 * n.links_t.p0[discharge_link_indices],
                            n.links_t.p0[store_link_indices],
                        ],
                        axis=1,
                    )
                    .T.groupby(n.links.bus0)
                    .sum()
                    .T.sum(axis=1)
                )
                df = pd.DataFrame(battery_power_flow, columns=["powerFlow"])

                store_indices = n.stores[
                    (n.stores.country == country) & (n.stores.type.str.contains("BATS"))
                ].index
                battery_state_of_charge = (
                    n.stores_t.e[store_indices]
                    .T.groupby(n.stores.bus)
                    .sum()
                    .T.sum(axis=1)
                )

                df["stateOfCharge"] = pd.DataFrame(battery_state_of_charge)
                df = df.reset_index()
            elif storageunit_indices is not None:
                storage_unit_indices = n.storage_units[
                    (n.storage_units.country == country)
                    & (n.storage_units.type.str.contains("BATS"))
                ].index
                battery_power_flow = (
                    n.storage_units_t.p[storage_unit_indices]
                    .T.groupby([n.storage_units.bus])
                    .sum()
                    .T.sum(axis=1)
                )
                df = pd.DataFrame(battery_power_flow, columns=["powerFlow"])
                battery_state_of_charge = (
                    n.storage_units_t.state_of_charge[storage_unit_indices]
                    .T.groupby([n.storage_units.bus])
                    .sum()
                    .T.sum(axis=1)
                )
                df["state_of_charge"] = pd.DataFrame(battery_state_of_charge)
                df = df.reset_index()

            df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(df), freq=f"{nth_hour}h"
            )
            df = pd.melt(df, id_vars=["snapshot"], var_name="status").set_index(
                "snapshot"
            )
            df["country"] = country
            df = df.set_index("country", append=True)
            all_countries_df = pd.concat([all_countries_df, df])

        return all_countries_df

    def pow_battery_flows_by_region_hourly(
        self, year: int, nth_hour: int
    ) -> pd.DataFrame:
        """Calculate hourly battery flow by country.

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, flow, snapshot) and column (value)
        """
        all_country_df = pd.DataFrame()
        for country in self.countries:
            n = self.network_dict[year]
            store_links = n.links[
                (n.links.country == country) & (n.links.type == ("BATS_STORE"))
            ].index
            discharge_links = n.links[
                (n.links.country == country) & (n.links.type == ("BATS_DISCHARGE"))
            ].index
            charging = n.links_t.p0[store_links].T.groupby(n.links.bus0).sum()
            charging.index.name = "bus"
            discharging = (
                n.links_t.p1[discharge_links].mul(-1).T.groupby(n.links.bus1).sum()
            )
            discharging.index.name = "bus"
            flow = pd.concat([charging, discharging], axis=0).groupby("bus").sum().T
            flow = flow.reset_index()
            flow["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(flow), freq=f"{nth_hour}h"
            )
            flow = (
                pd.melt(flow, id_vars=["snapshot"], var_name="flow")
                .replace("_HVELEC", "", regex=True)
                .set_index("snapshot")
            )
            flow["country"] = country
            flow = flow.set_index("country", append=True)
            all_country_df = pd.concat([all_country_df, flow], axis=1)
        return all_country_df

    def pow_hphs_flows_by_region_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly hydro pumped-storage flow by country.

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, flow, snapshot) and column (value)
        """
        all_country_df = pd.DataFrame()
        for country in self.countries:
            n = self.network_dict[year]
            store_links = n.links[
                (n.links.country == country) & (n.links.type == ("HPHS_STORE"))
            ].index
            discharge_links = n.links[
                (n.links.country == country) & (n.links.type == ("HPHS_DISCHARGE"))
            ].index
            charging = n.links_t.p0[store_links].T.groupby(n.links.bus0).sum()
            charging.index.name = "bus"
            discharging = (
                n.links_t.p1[discharge_links].mul(-1).T.groupby(n.links.bus1).sum()
            )
            discharging.index.name = "bus"
            flow = pd.concat([charging, discharging], axis=0).groupby("bus").sum().T
            flow = flow.reset_index()
            flow["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(flow), freq=f"{nth_hour}h"
            )
            flow = (
                pd.melt(flow, id_vars=["snapshot"], var_name="flow")
                .replace("_HVELEC", "", regex=True)
                .set_index("snapshot")
            )
            flow["country"] = country
            flow = flow.set_index("country", append=True)
            all_country_df = pd.concat([all_country_df, flow], axis=1)
        return all_country_df

    # Electricity hourly load by sector (MW)
    def pow_elec_load_by_sector_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly electricity load by sector and country.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, sector, snapshot) and column (value)
        """
        n = self.network_dict[year]
        all_countries_load = pd.DataFrame()
        for country in self.countries:
            country_load_df = pd.DataFrame()
            # Inflexible electrical load attached to HV/LV buses
            hv_lv_load_indices = n.loads[
                (n.loads.country == country)
                & (n.loads.bus.str.contains("LVELEC|HVELEC", regex=True))
            ].index

            inflexible_power_load = n.loads_t.p
            if len(inflexible_power_load.columns) > 0:
                inflexible_power_load = (
                    n.loads_t.p[hv_lv_load_indices].T.groupby(n.loads.carrier).sum().T
                )

                inflexible_power_load.columns = ["Inflexible_Power"]
            country_load_df = pd.concat(
                [country_load_df, inflexible_power_load], axis=1
            )
            # Flexible industry electrical load for heat supply
            industry_heat_load_indices = n.links[
                (n.links.country == country)
                & (n.links.bus1.str.contains("IND"))
                & (n.links.carrier == "Electricity")
            ].index
            flexible_industry_heat_load = pd.DataFrame(
                (-n.links_t.p1[industry_heat_load_indices]).sum(axis=1),
                columns=["Flexible_Industry_Heat"],
            )
            country_load_df = pd.concat(
                [country_load_df, flexible_industry_heat_load], axis=1
            ).reset_index()

            country_load_df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(country_load_df), freq=f"{nth_hour}h"
            )

            # # Flexible EV charging load
            ev_charging_load = self.tran_load_charging_all_hourly(year, nth_hour)
            ev_charging_load = ev_charging_load[
                ev_charging_load["country"] == country
            ].drop("country", axis=1)
            ev_charging_load = ev_charging_load.loc[
                ev_charging_load.status == "charging"
            ][["value"]]
            ev_charging_load = ev_charging_load.rename(
                columns={"value": "Flexible_EV_Charging"}
            )

            country_load_df = pd.concat(
                [country_load_df.set_index("snapshot"), ev_charging_load], axis=1
            ).reset_index()

            country_load_df = pd.melt(
                country_load_df, id_vars=["snapshot"], var_name="sector"
            ).set_index("snapshot")
            country_load_df["country"] = country
            all_countries_load = pd.concat(
                [all_countries_load, country_load_df], axis=0
            )
        return all_countries_load.set_index("country", append=True)

    def pow_emi_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual power sector CO2 emissions by country and carrier.

        << Unit: Mt_CO2 >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            pow_emit = n.links[
                (n.links.bus2.str.contains("ATMP")) & (n.links.type != "IND_BOILER")
            ].index
            emission = (
                (
                    n.links_t.p0[pow_emit].multiply(
                        n.snapshot_weightings.objective, axis=0
                    )
                )
                .multiply(n.links.reindex(pow_emit).efficiency2, axis=1)
                .T.groupby([n.links.country, n.links.carrier])
                .sum()
                .T.sum()
            ).to_frame()
            emission.columns = ["value"]  # return emission mix for a single year
            emission["year"] = year
            final_df = pd.concat(
                [final_df, emission], axis=0
            )  # return emission mix for all years
        final_df.index.names = ["country", "carrier"]
        final_df = final_df.fillna(0)
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        return final_df

    def pow_emi_cost_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual power sector CO2 emission costs by country and carrier.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier) and columns (value, year)
        """
        df_all = pd.DataFrame()
        for country in self.countries:
            if self.config["co2_management"][country]["option"] != "co2_price":
                continue
            df_country = pd.DataFrame()
            for year in self.network_dict:
                co2_price = self.config["co2_management"][country]["value"][year]
                emi = (
                    self.pow_emi_by_carrier_yearly()
                    .loc[country]
                    .set_index(["carrier", "year"], append=True)["value"]
                    .unstack()[year]
                    .fillna(0)
                )
                emi_cost = emi.mul(co2_price).to_frame("value")
                emi_cost["year"] = year
                df_country = pd.concat([df_country, emi_cost], axis=0)
                df_country.index.names = ["country", "carrier"]
            df_all = pd.concat([df_all, df_country], axis=0)
        if df_all.empty:
            return pd.DataFrame(columns=["year", "value"])  # return an empty dataframe

        df_all = df_all.loc[df_all.value != 0]
        return df_all[["year", "value"]]

    def pow_intercap_by_region_yearly(self) -> pd.DataFrame:
        """Calculate interconnection capacities across regions by country and year.

        << Unit: GW >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, from, to) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            inter_cap = (
                n.links[
                    (n.links.type == "ITCN") & (n.links.bus1.str.contains("HVELEC"))
                ]["p_nom_opt"]
                .groupby([n.links.country, n.links.bus0, n.links.bus1])
                .sum()
                .to_frame()
            )
            inter_cap.columns = ["value"]
            inter_cap["year"] = year
            final_df = pd.concat(
                [final_df, inter_cap], axis=0
            )  # return interconnection cap for all years
        final_df.index.names = ["country", "from", "to"]
        final_df = final_df.fillna(0)
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e3,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.reset_index()
        final_df = final_df.replace("_HVELEC", "", regex=True).set_index(
            ["country", "from", "to"]
        )
        return final_df

    def pow_inter_cf_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual interconnection capacity factors across regions by country.

        << Unit: per unit >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, from, to) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            inter_cf = (
                n.statistics.capacity_factor(
                    comps="Link", groupby=["type", "country", "bus0", "bus1"]
                )
                .loc["ITCN"]
                .reset_index()
            )
            inter_cf = inter_cf[inter_cf.bus1.str.contains("HVELEC")]
            inter_cf["year"] = year
            final_df = pd.concat(
                [final_df, inter_cf], axis=0
            )  # return interconnection flh for all years
        final_df = final_df.rename(
            columns={
                0: "value",
                "bus0": "from",
                "bus1": "to",
            }
        )
        final_df["value"] = (final_df["value"]).round(2)
        final_df = final_df.fillna(0)
        final_df = (
            final_df.replace("_HVELEC", "", regex=True)
            .sort_index()
            .set_index(["country", "from", "to"])
        )
        return final_df

    def pow_cap_by_type_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual power generation capacity by country, region, and type.

        << Unit: GW >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, region, technology) and
            columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cap_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                pow_supply = df[
                    (
                        (df[bus_name].str.contains("HVELEC"))
                        | (df[bus_name].str.contains("LVELEC"))
                    )
                    & ~(
                        (df.type == "LSLO")
                        | (df.type == "ITCN")
                        | (df.type.str.contains("SUPPLY"))
                    )  # removes the non-capacity items
                ]
                if c.name == "Link":
                    pow_supply_cap = (
                        (pow_supply.p_nom_opt * pow_supply.efficiency)
                        .groupby(
                            [
                                pow_supply.country,
                                pow_supply[bus_name].str.split("_", n=2).str[1],
                                pow_supply.type,
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                else:
                    pow_supply_cap = (
                        pow_supply.p_nom_opt.groupby(
                            [
                                pow_supply.country,
                                pow_supply[bus_name].str.split("_", n=2).str[1],
                                pow_supply.type,
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                pow_supply_cap.columns = ["value"]
                cap_mix = pd.concat(
                    [cap_mix, pow_supply_cap], axis=0
                )  # return cap mix for a single year
            cap_mix["year"] = year
            final_df = pd.concat(
                [final_df, cap_mix], axis=0
            )  # return cap mix for all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "region", "technology"]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e3,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        return final_df

    # ==================================== INDUSTRY ====================================
    def ind_emi_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial CO2 emissions by country and carrier.

        << Unit: Mt_CO2 >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            ind_emit = n.links[
                (n.links.bus2.str.contains("ATMP"))
                & (n.links.type.str.contains("IND_BOILER"))
            ].index
            emission = (
                (
                    n.links_t.p0[ind_emit].multiply(
                        n.snapshot_weightings.objective, axis=0
                    )
                )
                .multiply(n.links.reindex(ind_emit).efficiency2, axis=1)
                .T.groupby([n.links.country, n.links.carrier])
                .sum()
                .T.sum()
            ).to_frame()
            emission.columns = ["value"]  # return emission mix for a single year
            emission["year"] = year
            final_df = pd.concat(
                [final_df, emission], axis=0
            )  # return emission mix for all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "carrier"]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        return final_df

    def ind_emi_cost_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial CO2 emission costs by country and carrier.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier) and columns (value, year)
        """
        df_all = pd.DataFrame()
        for country in self.countries:
            if self.config["co2_management"][country]["option"] != "co2_price":
                continue
            df_country = pd.DataFrame()
            for year in self.network_dict:
                co2_price = self.config["co2_management"][country]["value"][year]
                emi = (
                    self.ind_emi_by_carrier_yearly()
                    .loc[country]
                    .set_index(["carrier", "year"], append=True)["value"]
                    .unstack()[year]
                    .fillna(0)
                )
                emi_cost = emi.mul(co2_price).to_frame("value")
                emi_cost["year"] = year
                df_country = pd.concat([df_country, emi_cost], axis=0)
                df_country.index.names = ["country", "carrier"]
            df_all = pd.concat([df_all, df_country], axis=0)

        # return an empty dataframe
        if df_all.empty:
            return pd.DataFrame(columns=["year", "value"])

        # Filter out zero values and return the dataframe when not empty
        df_all = df_all.loc[df_all.value != 0]
        return df_all[["year", "value"]]

    def ind_gen_by_type_by_carrier_by_heatgroup_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial gen. by country, type, carrier, and heat group.

        << Unit: TWh >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, carrier, heat_type) and
            columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            gen_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                pow_supply = c.df[(c.df[bus_name].str.contains("IND"))].index
                if c.name == "Generator":
                    c_gen = (
                        c.pnl.p[pow_supply]
                        .multiply(n.snapshot_weightings.generators, axis=0)
                        .sum()
                        .multiply(c.df.sign)
                        .groupby(
                            [
                                c.df.country,
                                c.df.type,
                                c.df.carrier,
                                c.df[bus_name].str.split("_", n=2).str[2],
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                elif c.name == "StorageUnit":
                    c_gen = (
                        c.pnl.p[[x for x in c.pnl.p.columns if x in pow_supply]]
                        .multiply(n.snapshot_weightings.stores, axis=0)
                        .sum()
                        .multiply(c.df.sign)
                        .groupby(
                            [
                                c.df.country,
                                c.df.type,
                                c.df.carrier,
                                c.df[bus_name].str.split("_", n=2).str[2],
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                else:
                    c_gen = (
                        (
                            -c.pnl.p1[pow_supply]
                            .multiply(n.snapshot_weightings.generators, axis=0)
                            .sum()
                        )
                        .groupby(
                            [
                                c.df.country,
                                c.df.type,
                                c.df.carrier,
                                c.df[bus_name].str.split("_", n=2).str[2],
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                c_gen.columns = ["value"]
                gen_mix = pd.concat(
                    [gen_mix, c_gen], axis=0
                )  # return gen mix for a single year
            gen_mix["year"] = year
            final_df = pd.concat(
                [final_df, gen_mix], axis=0
            )  # return gen mix for all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology", "carrier", "heat_type"]
        final_df = final_df[
            ~final_df.index.get_level_values("technology").str.contains(
                "LSLO|DISCHARGE|STORE"
            )
        ]  # remove store and discharge, since we only consider generation
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=4,
        )
        final_df = final_df.loc[final_df.value != 0]
        return final_df

    def ind_cap_by_carrier_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial cap. by country, carrier, region and heat type.

        << Unit: GW >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier, region, heat_type) and
            columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cap_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                ind_supply = df[
                    ((df[bus_name].str.contains("IND")))
                    & ~(
                        (df.type == "LSLO")
                        | (df.type == "ITCN")
                        | (df.type.str.contains("SUPPLY"))
                    )
                    & ~(df.carrier == "ENS")
                ]
                if c.name == "Link":

                    ind_supply_cap = (
                        (ind_supply.p_nom_opt * ind_supply.efficiency)
                        .groupby(
                            [
                                ind_supply.country,
                                ind_supply.carrier,
                                ind_supply[bus_name].str.split("_", n=2).str[1],
                                ind_supply[bus_name].str.split("_", n=2).str[2],
                            ]
                        )  # groupby using the region_name of the bus
                        .sum()
                        .to_frame()
                    )
                else:
                    ind_supply_cap = (
                        ind_supply.p_nom_opt.groupby(
                            [
                                ind_supply.country,
                                ind_supply.carrier,
                                ind_supply[bus_name].str.split("_", n=2).str[1],
                                ind_supply[bus_name].str.split("_", n=2).str[2],
                            ]
                        )
                        .sum()
                        .to_frame()  # groupby using the region_name of the bus
                    )
                ind_supply_cap.columns = ["value"]
                cap_mix = pd.concat(
                    [cap_mix, ind_supply_cap], axis=0
                )  # return cap mix for a single year

            cap_mix["year"] = year
            final_df = pd.concat(
                [final_df, cap_mix], axis=0
            )  # return cap mix for all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "carrier", "region", "heat_type"]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e3,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    # Regional Ind generation by year (in TWh)
    def ind_gen_by_carrier_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial gen. by country, carrier, region and heat type.

        << Unit: TWh >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier, region, heat_type) and
            columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            gen_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                ind_supply = c.df[
                    ((c.df[bus_name].str.contains("IND")))
                    & ~(
                        (c.df.type == "LSLO")
                        | (c.df.type == "ITCN")
                        | (c.df.type.str.contains("SUPPLY"))
                        | (
                            c.df.type.isin(
                                ["EVST_PRV", "EVST_PUB", "BATS", "HHBS", "HPHS"]
                            )
                        )
                    )  # removes the non-capacity items
                ].index
                if c.name == "Generator":
                    c_gen = (
                        c.pnl.p[ind_supply]
                        .multiply(n.snapshot_weightings.generators, axis=0)
                        .sum()
                        .multiply(c.df.sign)
                        .groupby(
                            [
                                c.df.country,
                                c.df.carrier,
                                c.df[bus_name].str.split("_", n=2).str[1],
                                c.df[bus_name].str.split("_", n=2).str[2],
                            ]
                        )  # replace with region
                        .sum()
                        .to_frame()
                    )
                elif c.name == "StorageUnit":
                    c_gen = (
                        c.pnl.p[[x for x in c.pnl.p.columns if x in ind_supply]]
                        .multiply(n.snapshot_weightings.stores, axis=0)
                        .sum()
                        .multiply(c.df.sign)
                        .groupby(
                            [
                                c.df.country,
                                c.df.carrier,
                                c.df[bus_name].str.split("_", n=2).str[1],
                                c.df[bus_name].str.split("_", n=2).str[2],
                            ]
                        )  # replace with region
                        .sum()
                        .to_frame()
                    )
                else:

                    c_gen = (
                        (
                            -c.pnl.p1[ind_supply]
                            .multiply(n.snapshot_weightings.generators, axis=0)
                            .sum()
                        )
                        .groupby(
                            [
                                c.df.country,
                                c.df.carrier,
                                c.df[bus_name].str.split("_", n=2).str[1],
                                c.df[bus_name].str.split("_", n=2).str[2],
                            ]
                        )  # replace with region
                        .sum()
                        .to_frame()
                    )
                c_gen.columns = ["value"]
                c_gen.index.names = ["country", "carrier", "region", "heat_type"]
                gen_mix = pd.concat(
                    [gen_mix, c_gen], axis=0
                )  # return gen mix for a single year

                gen_mix = gen_mix.groupby(
                    ["country", "carrier", "region", "heat_type"]
                ).sum()  # need to groupby-sum per region

            gen_mix["year"] = year
            final_df = pd.concat(
                [final_df, gen_mix], axis=0
            )  # return gen mix for all years
        final_df = final_df.fillna(0)
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    def ind_cap_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial cap. by country, type, carrier and heat type.

        << Unit: GW >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, carrier, heat_type) and
            columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cap_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                ind_supply = df[(df[bus_name].str.contains("IND"))]
                if c.name == "Link":
                    ind_supply_cap = (
                        (ind_supply.p_nom_opt * ind_supply.efficiency)
                        .groupby(
                            [
                                ind_supply.country,
                                ind_supply.type,
                                ind_supply.carrier,
                                ind_supply[bus_name].str.split("_", n=2).str[2],
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                else:
                    ind_supply_cap = (
                        ind_supply.p_nom_opt.groupby(
                            [
                                ind_supply.country,
                                ind_supply.type,
                                ind_supply.carrier,
                                ind_supply[bus_name].str.split("_", n=2).str[2],
                            ]
                        )
                        .sum()
                        .to_frame()
                    )
                ind_supply_cap.columns = ["value"]
                cap_mix = pd.concat(
                    [cap_mix, ind_supply_cap], axis=0
                )  # return cap mix for a single year
            cap_mix["year"] = year
            final_df = pd.concat(
                [final_df, cap_mix], axis=0
            )  # return cap mix for all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology", "carrier", "heat_type"]
        final_df = final_df[
            ~final_df.index.get_level_values("technology").str.contains("ITCN|LSLO")
        ]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e3,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    def ind_gen_by_type_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly industrial generation by country and type.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, region, technology, snapshot) and
            column (value)
        """
        network = self.network_dict[year]
        all_countries_df = pd.DataFrame()
        for country in self.countries:
            country_dispatch_df = pd.DataFrame()
            for component in network.iterate_components(
                ["Link", "Generator", "StorageUnit"]
            ):
                bus_col = "bus1" if component.name == "Link" else "bus"
                heat_supply_indices = component.df[
                    (component.df.country == country)
                    & (component.df[bus_col].str.contains("IND"))
                ].index
                if component.name in ["Generator", "StorageUnit"]:
                    generation_df = (
                        component.pnl.p[
                            [
                                idx
                                for idx in component.pnl.p.columns
                                if idx in heat_supply_indices
                            ]
                        ]
                        .T.groupby(
                            [
                                component.df[bus_col].str.split("_", n=2).str[1],
                                component.df.type,
                            ]
                        )
                        .sum()
                        .T
                    )
                else:
                    generation_df = (
                        (
                            -component.pnl.p1[
                                [
                                    idx
                                    for idx in component.pnl.p1.columns
                                    if idx in heat_supply_indices
                                ]
                            ]
                        )
                        .T.groupby(
                            [
                                component.df[bus_col].str.split("_", n=2).str[1],
                                component.df.type,
                            ]
                        )
                        .sum()
                        .T
                    )
                country_dispatch_df = pd.concat(
                    [country_dispatch_df, generation_df], axis=1
                )

            country_dispatch_df = (
                country_dispatch_df.loc[:, (country_dispatch_df != 0).all(axis=0)]
            ).round(1)
            country_dispatch_df = country_dispatch_df.reset_index()
            country_dispatch_df["snapshot"] = pd.date_range(
                start=f"{year}-01-01",
                periods=len(country_dispatch_df),
                freq=f"{nth_hour}h",
            )

            storage_cols = list(
                country_dispatch_df.columns.get_level_values("type").str.contains(
                    "_STORE"
                )
            )
            country_dispatch_df[country_dispatch_df.columns[storage_cols]] *= -1

            # Flatten columns if the df is multiIndex
            country_dispatch_df.columns = [
                "!".join(filter(None, map(str, col))) if isinstance(col, tuple) else col
                for col in country_dispatch_df.columns
            ]

            melted_df = pd.melt(
                country_dispatch_df,
                id_vars="snapshot",
                var_name="region!technology",
            ).set_index("snapshot")

            melted_df[["region", "technology"]] = melted_df[
                "region!technology"
            ].str.split("!", expand=True)
            melted_df = melted_df.drop(columns=["region!technology"])
            melted_df["country"] = country
            all_countries_df = pd.concat([all_countries_df, melted_df])
        return all_countries_df.set_index("country", append=True)

    def ind_capex_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial CAPEX by country, type, and carrier.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, carrier) and columns
            (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in ["Link", "Generator", "StorageUnit", "Store"]:
                bus_name = "bus1" if c == "Link" else "bus"
                included_substrings = ["IND", "INLHSTORN"]
                if c == "Store":
                    ind_cost = n.statistics.capex(
                        comps=c, groupby=[bus_name, "country", "type", "carrier"]
                    )
                else:
                    ind_cost = n.statistics.capex(
                        comps=c, groupby=[bus_name, "country", "type", "carrier"]
                    ).sub(
                        n.statistics.capex(
                            cost_attribute="fom_cost",
                            comps=c,
                            groupby=[bus_name, "country", "type", "carrier"],
                        ),
                        fill_value=0,
                    )
                if not ind_cost.empty:
                    list_of_bus = (
                        ind_cost.index.get_level_values(bus_name).unique().tolist()
                    )
                    ind_bus = [
                        x
                        for x in list_of_bus
                        if any(substring in x for substring in included_substrings)
                    ]
                    if len(ind_bus) > 0:
                        ind_cost = (
                            ind_cost.loc[ind_bus]
                            .groupby(["country", "type", "carrier"])
                            .sum()
                            .to_frame(name="value")
                        )
                        cost_mix = pd.concat(
                            [cost_mix, ind_cost], axis=0
                        )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology", "carrier"]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "carrier", "year"]).sum()

        return final_df

    def ind_fom_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial fixed O&M costs by country, type, and carrier.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, carrier) and columns
            (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in ["Link", "Generator", "StorageUnit"]:
                bus_name = "bus1" if c == "Link" else "bus"
                included_substrings = ["IND", "INLHSTORN"]
                pow_cost = n.statistics.capex(
                    cost_attribute="fom_cost",
                    comps=c,
                    groupby=[bus_name, "country", "type", "carrier"],
                )
                if not pow_cost.empty:
                    list_of_bus = (
                        pow_cost.index.get_level_values(bus_name).unique().tolist()
                    )
                    elec_bus = [
                        x
                        for x in list_of_bus
                        if any(substring in x for substring in included_substrings)
                    ]
                    if len(elec_bus) > 0:
                        pow_cost = (
                            pow_cost.loc[elec_bus]
                            .groupby(["country", "type", "carrier"])
                            .sum()
                            .to_frame(name="value")
                        )
                        cost_mix = pd.concat(
                            [cost_mix, pow_cost], axis=0
                        )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology", "carrier"]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "carrier", "year"]).sum()

        return final_df

    def ind_opex_by_type_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual industrial OPEX by country, type, and carrier.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, carrier) and columns
            (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in ["Link", "Generator", "StorageUnit", "Store"]:
                bus_name = "bus1" if c == "Link" else "bus"
                included_substrings = ["IND"]
                # new opex = opex + fom_cost + fuel_cost
                if c == "Store":
                    ind_cost = n.statistics.opex(
                        comps=c, groupby=[bus_name, "country", "type", "carrier"]
                    )
                else:
                    ind_cost = n.statistics.opex(
                        comps=c, groupby=[bus_name, "country", "type", "carrier"]
                    ).add(
                        n.statistics.capex(
                            cost_attribute="fom_cost",
                            comps=c,
                            groupby=[bus_name, "country", "type", "carrier"],
                        ),
                        fill_value=0,
                    )
                if c == "Link":
                    # Create a dictionary comprehension
                    # to map links to their prices from primary bus0
                    fuel_price_dict = {
                        link: n.buses_t.marginal_price[n.links.loc[link, "bus0"]]
                        for link in n.links_t.p0.columns
                    }
                    # Convert the dictionary into a DataFrame
                    fuel_price_df = pd.DataFrame.from_dict(fuel_price_dict)
                    fuel_cost = (
                        n.links_t.p0.multiply(n.snapshot_weightings.objective, axis=0)
                        .multiply(fuel_price_df)
                        .T.groupby(
                            [
                                n.links["bus1"],
                                n.links.country,
                                n.links.type,
                                n.links.carrier,
                            ]
                        )
                        .sum()
                        .T.sum()
                    )
                    # Set fuel cost of storage and interconnection link
                    # as 0 to avoid double counting
                    target_substrings = {"_STORE", "_DISCHARGE", "ITCN"}
                    targeted_type = [
                        any(substring in level for substring in target_substrings)
                        for level in fuel_cost.index.get_level_values("type")
                    ]
                    fuel_cost.loc[targeted_type] = 0
                    # add fuel cost to opex
                    ind_cost = ind_cost.add(fuel_cost)
                if not ind_cost.empty:
                    list_of_bus = (
                        ind_cost.index.get_level_values(bus_name).unique().tolist()
                    )
                    ind_bus = [
                        x
                        for x in list_of_bus
                        if any(substring in x for substring in included_substrings)
                    ]
                    if len(ind_bus) > 0:
                        ind_cost = (
                            ind_cost.loc[ind_bus]
                            .groupby(["country", "type", "carrier"])
                            .sum()
                            .to_frame(name="value")
                        )
                        cost_mix = pd.concat(
                            [cost_mix, ind_cost], axis=0
                        )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology", "carrier"]

        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]
        final_df = final_df.groupby(["country", "technology", "carrier", "year"]).sum()

        return final_df

    def ind_lh_gen_by_type_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly low-heat industrial generation by country and type.

        << Unit: MWh >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, snapshot) and
            column (value)
        """
        all_countries_df = pd.DataFrame()
        for country in self.countries:
            network = self.network_dict[year]
            country_df = pd.DataFrame()
            for component in network.iterate_components(
                ["Link", "Generator", "StorageUnit"]
            ):
                bus_col = "bus1" if component.name == "Link" else "bus"
                low_heat_indices = component.df[
                    (component.df.country == country)
                    & (component.df[bus_col].str.contains("IND-LH"))
                ].index
                if component.name in ["Generator", "StorageUnit"]:
                    generation_df = (
                        component.pnl.p[
                            [
                                idx
                                for idx in component.pnl.p.columns
                                if idx in low_heat_indices
                            ]
                        ]
                        .T.groupby(component.df.type)
                        .sum()
                        .T
                    )
                else:
                    generation_df = (
                        (
                            -component.pnl.p1[
                                [
                                    idx
                                    for idx in component.pnl.p1.columns
                                    if idx in low_heat_indices
                                ]
                            ]
                        )
                        .T.groupby(component.df.type)
                        .sum()
                        .T
                    )
                country_df = pd.concat([country_df, generation_df], axis=1)

            country_df = (country_df.loc[:, (country_df != 0).all(axis=0)]).round(1)
            country_df = country_df.reset_index()
            country_df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(country_df), freq=f"{nth_hour}h"
            )
            storage_columns = list(country_df.columns.str.contains("_STORE"))
            country_df[country_df.columns[storage_columns]] *= -1

            melted_df = pd.melt(
                country_df, id_vars=["snapshot"], var_name="technology"
            ).set_index("snapshot")
            melted_df["country"] = country
            all_countries_df = pd.concat([all_countries_df, melted_df])
        return all_countries_df.set_index("country", append=True)

    def ind_hh_gen_by_type_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate hourly high-heat industrial generation by country and type.

        << Unit: MWh >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology, snapshot) and
            column (value)
        """
        all_countries_df = pd.DataFrame()
        for country in self.countries:
            network = self.network_dict[year]
            country_df = pd.DataFrame()
            for component in network.iterate_components(
                ["Link", "Generator", "StorageUnit"]
            ):
                bus_col = "bus1" if component.name == "Link" else "bus"
                low_heat_indices = component.df[
                    (component.df.country == country)
                    & (component.df[bus_col].str.contains("IND-HH"))
                ].index
                if component.name in ["Generator", "StorageUnit"]:
                    generation_df = (
                        component.pnl.p[
                            [
                                idx
                                for idx in component.pnl.p.columns
                                if idx in low_heat_indices
                            ]
                        ]
                        .T.groupby(component.df.type)
                        .sum()
                        .T
                    )
                else:
                    generation_df = (
                        (
                            -component.pnl.p1[
                                [
                                    idx
                                    for idx in component.pnl.p1.columns
                                    if idx in low_heat_indices
                                ]
                            ]
                        )
                        .T.groupby(component.df.type)
                        .sum()
                        .T
                    )
                country_df = pd.concat([country_df, generation_df], axis=1)

            country_df = (country_df.loc[:, (country_df != 0).all(axis=0)]).round(1)
            country_df = country_df.reset_index()
            country_df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(country_df), freq=f"{nth_hour}h"
            )

            storage_columns = list(country_df.columns.str.contains("_STORE"))
            country_df[country_df.columns[storage_columns]] *= -1

            melted_df = pd.melt(
                country_df, id_vars=["snapshot"], var_name="technology"
            ).set_index("snapshot")
            melted_df["country"] = country
            all_countries_df = pd.concat([all_countries_df, melted_df])
        return all_countries_df.set_index("country", append=True)

    # ==================================== TRANSPORT ===================================

    def tran_capex_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual transport CAPEX by country and type.

        << Unit: Million CURRENCY >>

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, technology) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            cost_mix = pd.DataFrame()
            for c in n.iterate_components(["Link", "Generator", "StorageUnit"]):
                bus_name = "bus1" if c.name == "Link" else "bus"
                df = c.df
                df["capex"] = df["capital_cost"] - df["fom_cost"]
                tran_supply = df[(df[bus_name].str.contains("TRAN"))]
                if c.name == "Link":
                    tran_cost = (
                        (
                            tran_supply.p_nom_opt
                            * tran_supply.capex
                            * tran_supply.efficiency
                        )
                        .groupby([tran_supply.country, tran_supply.type])
                        .sum()
                    ).to_frame()
                else:
                    tran_cost = (
                        (tran_supply.p_nom_opt * tran_supply.capex)
                        .groupby([tran_supply.country, tran_supply.type])
                        .sum()
                    ).to_frame()
                tran_cost.columns = ["value"]
                cost_mix = pd.concat(
                    [cost_mix, tran_cost], axis=0
                )  # return cost mix for a single year
            cost_mix["year"] = year
            final_df = pd.concat(
                [final_df, cost_mix], axis=0
            )  # return cost mix for a all years
        final_df = final_df.fillna(0)
        final_df.index.names = ["country", "technology"]
        final_df = scaling_conversion(
            df=final_df.loc[~(final_df == 0).all(axis=1), :],
            scaling_number=1e6,
            decimals=2,
        )
        final_df = final_df.loc[final_df.value != 0]

        return final_df

    def tran_charger_capacity_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual EV charger capacity by country and type.

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            df = (
                n.links[n.links.type.str.contains("EVCH")]
                .groupby([n.links.country, n.links.type.str.split("-", n=1).str[1]])
                .sum()
                .p_nom.to_frame()
            )
            df.columns = ["value"]
            df["year"] = year
            final_df = pd.concat([final_df, df], axis=0)
        return final_df

    def tran_charger_capacity_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual EV charger capacity by country and region.

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, region) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            df = (
                n.links[n.links.type.str.contains("EVCH")]
                .groupby([n.links.country, n.links.bus0.str.split("_", n=2).str[1]])
                .sum()
                .p_nom.to_frame()
            )
            df = df.reset_index()
            df.columns = ["country", "region", "value"]
            df = df.set_index(["country", "region"])
            df["year"] = year
            final_df = pd.concat([final_df, df], axis=0)
        return final_df

    def tran_storage_capacity_by_type_yearly(self) -> pd.DataFrame:
        """Calculate annual EV storage capacity by country and type.

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            df = (
                n.stores[
                    (n.stores.carrier == "Electricity")
                    & (n.stores.index.str.contains("TRAN"))
                ]
                .groupby([n.stores.country, n.stores.type.str.split("-", n=1).str[1]])
                .sum()
                .e_nom.to_frame()
            )
            df.columns = ["value"]
            df["year"] = year
            final_df = pd.concat([final_df, df], axis=0)
        return final_df

    def tran_storage_capacity_by_region_yearly(self) -> pd.DataFrame:
        """Calculate annual EV storage capacity by country and region.

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, region) and columns (value, year)
        """
        final_df = pd.DataFrame()
        for year in self.network_dict:
            n = self.network_dict[year]
            df = (
                n.stores[
                    (n.stores.carrier == "Electricity")
                    & (n.stores.index.str.contains("TRAN"))
                ]
                .groupby([n.stores.country, n.stores.bus.str.split("_", n=2).str[1]])
                .sum()
                .e_nom.to_frame()
            )
            df = df.reset_index()
            df.columns = ["country", "region", "value"]
            df = df.set_index(["country", "region"])
            df["year"] = year
            final_df = pd.concat([final_df, df], axis=0)
        return final_df

    def tran_load_charging_all_hourly(self, year: int, nth_hour: int) -> pd.DataFrame:
        """Calculate total EV load and total charging for all regions combined.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (snapshot, status, country) and index (snapshot)
        """
        df_all = pd.DataFrame()
        for country in self.countries:
            n = self.network_dict[year]
            evlinks = n.links[
                (n.links.country == country) & (n.links.type.str.contains("EVCH"))
            ].index

            ev_charging = (
                n.links_t.p0[evlinks]
                .T.groupby(
                    [n.links.bus0]
                )  # does not differentiate between PRV and PUB EVs
                .sum()
                .T.sum(axis=1)
            )  # sums up different regions
            df = pd.DataFrame(ev_charging, columns=["charging"])

            loadsidx = n.loads[
                (n.loads.country == country)
                & (n.loads.carrier == "Electricity")
                & (n.loads.index.str.contains("TRAN"))
            ].index
            evload = n.loads_t.p_set[
                n.loads_t.p_set.columns.intersection(loadsidx)
            ].sum(axis=1)
            storeidx = n.stores[
                (n.stores.country == country)
                & (n.stores.carrier == "Electricity")
                & (n.stores.index.str.contains("TRAN"))
            ].index
            store_discharge = (
                n.stores_t.p[storeidx].T.groupby([n.stores.bus]).sum().T.sum(axis=1)
            )
            store_soc = (
                n.stores_t.e[storeidx].T.groupby([n.stores.bus]).sum().T.sum(axis=1)
            )

            # use of intersection as in some cases loads_t does not have all columns as
            # loadidx
            df["load"] = pd.DataFrame(evload)
            df["store_discharge"] = pd.DataFrame(store_discharge)
            df["state_of_charge"] = pd.DataFrame(store_soc)
            df = df.reset_index()
            df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(df), freq=f"{nth_hour}h"
            )
            df = pd.melt(df, id_vars=["snapshot"], var_name="status").set_index(
                "snapshot"
            )
            df["country"] = country
            df_all = pd.concat([df_all, df], axis=0)
        return df_all

    # EV charging profile and load all regions combined (in MW)
    def tran_load_charging_all_hourly_by_type(
        self, year: int, nth_hour: int
    ) -> pd.DataFrame:
        """Calculate total EV load and total charging for all regions combined by type.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (snapshot, type, status, country) and
            index (snapshot)
        """
        n = self.network_dict[year]
        df_all = pd.DataFrame()
        for country in self.countries:
            evlinks = n.links[
                (n.links.country == country) & (n.links.type.str.contains("EVCH"))
            ].index
            ev_charging = (
                n.links_t.p0[evlinks]
                .T.groupby([n.links.bus1.str.split("_", n=2).str[2]])
                .sum()
                .T
            )
            ev_charging = ev_charging.stack().to_frame().rename({0: "charging"}, axis=1)

            loadsidx = n.loads[
                (n.loads.country == country)
                & (n.loads.carrier == "Electricity")
                & (n.loads.index.str.contains("TRAN"))
                & ~(n.loads.index.str.contains("Rail"))
            ].index

            evload = (
                n.loads_t.p_set[n.loads_t.p_set.columns.intersection(loadsidx)]
                .T.groupby([n.loads.bus.str.split("_", n=2).str[2]])
                .sum()
                .T
            )
            evload = evload.stack().to_frame().rename({0: "load"}, axis=1)

            storeidx = n.stores[
                (n.stores.country == country)
                & (n.stores.carrier == "Electricity")
                & (n.stores.index.str.contains("TRAN"))
            ].index
            store_discharge = (
                n.stores_t.p[storeidx]
                .T.groupby([n.stores.bus.str.split("_", n=2).str[2]])
                .sum()
                .T
            )
            store_discharge = (
                store_discharge.stack()
                .to_frame()
                .rename({0: "store_discharge"}, axis=1)
            )

            store_soc = (
                n.stores_t.e[storeidx]
                .T.groupby([n.stores.bus.str.split("_", n=2).str[2]])
                .sum()
                .T
            )
            store_soc = (
                store_soc.stack().to_frame().rename({0: "state_of_charge"}, axis=1)
            )

            df = (
                pd.concat([ev_charging, evload, store_discharge, store_soc], axis=1)
                .reset_index()
                .rename({"level_1": "type"}, axis=1)
            )
            df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(df), freq=f"{nth_hour}h"
            )
            df = pd.melt(df, id_vars=["snapshot", "type"], var_name="status").set_index(
                "snapshot"
            )
            df["country"] = country
            df_all = pd.concat([df_all, df], axis=0)
        return df_all

    def tran_load_charging_all_hourly_by_region(
        self, year: int, nth_hour: int
    ) -> pd.DataFrame:
        """Calculate total EV load and charging for all regions combined by region.

        << Unit: MW >>

        Parameters
        ----------
        year : int
            year of interest
        nth_hour : int
            frequency of the time series in hours (e.g., 1 for hourly, 3 for every
            3 hours)

        Returns
        -------
        pd.DataFrame
            A DataFrame with columns (snapshot, region, status, country) and
            index (snapshot)
        """
        n = self.network_dict[year]
        df_all = pd.DataFrame()
        for country in self.countries:
            evlinks = n.links[
                (n.links.country == country) & (n.links.type.str.contains("EVCH"))
            ].index
            ev_charging = (
                n.links_t.p0[evlinks]
                .T.groupby([n.links.bus1.str.split("_", n=2).str[1]])
                .sum()
                .T
            )
            ev_charging = ev_charging.stack().to_frame().rename({0: "charging"}, axis=1)

            loadsidx = n.loads[
                (n.loads.country == country)
                & (n.loads.carrier == "Electricity")
                & (n.loads.index.str.contains("TRAN"))
                & ~(n.loads.index.str.contains("Rail"))
            ].index

            evload = (
                n.loads_t.p_set[n.loads_t.p_set.columns.intersection(loadsidx)]
                .T.groupby([n.loads.bus.str.split("_", n=2).str[1]])
                .sum()
                .T
            )
            evload = evload.stack().to_frame().rename({0: "load"}, axis=1)

            storeidx = n.stores[
                (n.stores.country == country)
                & (n.stores.carrier == "Electricity")
                & (n.stores.index.str.contains("TRAN"))
            ].index
            store_discharge = (
                n.stores_t.p[storeidx]
                .T.groupby([n.stores.bus.str.split("_", n=2).str[1]])
                .sum()
                .T
            )
            store_discharge = (
                store_discharge.stack()
                .to_frame()
                .rename({0: "store_discharge"}, axis=1)
            )

            store_soc = (
                n.stores_t.e[storeidx]
                .T.groupby([n.stores.bus.str.split("_", n=2).str[1]])
                .sum()
                .T
            )
            store_soc = (
                store_soc.stack().to_frame().rename({0: "state_of_charge"}, axis=1)
            )

            df = (
                pd.concat([ev_charging, evload, store_discharge, store_soc], axis=1)
                .reset_index()
                .rename({"level_1": "region"}, axis=1)
            )
            df["snapshot"] = pd.date_range(
                start=f"{year}-01-01", periods=len(df), freq=f"{nth_hour}h"
            )
            df = pd.melt(
                df, id_vars=["snapshot", "region"], var_name="status"
            ).set_index("snapshot")
            df["country"] = country
            df_all = pd.concat([df_all, df], axis=0)
        return df_all

    # =================================== TEST/OTHERS ==================================
    def test_energy_not_served_warning(self):
        """Check if there is any energy not served (ENS) in the model."""
        for year in self.network_dict:
            n = self.network_dict[year]
            lslo = (
                n.generators_t.p[n.generators[n.generators.carrier == "ENS"].index]
                .sum()
                .sum()
            )
            if lslo < 100:
                pass
                print("No energy not served observed")
            else:
                print("WARNING: Loss of load observed")
                gen_loss_of_load = n.generators_t.p[
                    n.generators[n.generators.carrier == "ENS"].index
                ].sum()
                print(
                    "Generators with loss of load (Unit: GWh)",
                    gen_loss_of_load[gen_loss_of_load > 10].div(1e3).round(),
                )

    def ene_gen_by_carrier_yearly(self) -> pd.DataFrame:
        """Calculate annual electricity generation by country and carrier.

        Notes
        -----
        Fossil resources are already primary side, but renewable resources still
        need transformation, subject to the RE transformation ratio used by country.

        Returns
        -------
        pd.DataFrame
            A DataFrame with multi-index (country, carrier) and columns (value, year)
        """
        final_df = pd.DataFrame()

        for year in self.network_dict:
            n = self.network_dict[year]

            # Components modelled as generators
            re_generators = n.generators[
                (n.generators["bus"].str.contains("HVELEC"))
                | (n.generators["bus"].str.contains("LVELEC"))
            ].index

            re_generators = (
                n.generators_t.p[re_generators]
                .multiply(n.snapshot_weightings.generators, axis=0)
                .sum()
                .multiply(n.generators.sign)
                .groupby([n.generators.country, n.generators.carrier])
                .sum()
                .to_frame()
            )

            # Components modelled as links + generators in the primary-side
            fuel_supply = n.generators[(n.generators.type.str.contains("SUPPLY"))].index

            fuel_supply = (
                n.generators_t.p[fuel_supply]
                .multiply(n.snapshot_weightings.generators, axis=0)
                .sum()
                .multiply(n.generators.sign)
                .groupby([n.generators.country, n.generators.carrier])
                .sum()
                .to_frame()
            )

            gen_mix = pd.concat([re_generators, fuel_supply], axis=0)

            gen_mix.columns = ["value"]
            gen_mix["year"] = year

            final_df = pd.concat([final_df, gen_mix], axis=0)

        final_df = final_df.fillna(0)

        # Filter out unnecessary technology
        final_df = final_df[~final_df.index.get_level_values("carrier").isin(["ENS"])]

        final_df = final_df.loc[final_df.value != 0]
        return final_df
