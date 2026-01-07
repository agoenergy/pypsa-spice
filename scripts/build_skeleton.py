# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Generate skeleton CSV files for setting up an energy system modelling.

It creates templates for buses, generators, storage, loads, interconnectors,
and converter links for power, industry and transport sectors.
"""

import itertools
import os
import shutil
import warnings
from itertools import combinations
from typing import Any

import numpy as np
import pandas as pd
from _helpers import (
    FilePath,
    configure_logging,
)
from build_scenario_config import build_scenario_config_file
from pandas.errors import SettingWithCopyWarning

warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
warnings.simplefilter(action="ignore", category=FutureWarning)

# Data definitions and helper functions

# List of power bus types
POWER_BUS_TYPES = [
    "HVELEC",
    "LVELEC",
    "BATSN",
    "HHBSN",
    "HPHSN",
    "LIGN",
    "BITN",
    "HRDCN",
    "OILN",
    "GASN",
    "LNGN",
    "BION",
    "ATMP",
    "HYDN",
    "CO2STORN",
    "WSTN",
    "NUCLN",
]

# List of CO2-free power plant technologies
POWER_TECHNOLOGIES = ["PHOT", "CSP", "HROR", "WTON", "WTOF", "GEOT", "RTPV"]

# Industry bus types and technologies
INDUSTRY_BUS_TYPES = ["IND-LH", "IND-HH", "INLHSTORN"]
INDUSTRY_INDICATORS = ["coal_cap", "oil_cap", "gas_cap", "bio_cap", "elec_cap"]
INDUSTRY_TECHNOLOGIES = ["SWHT"]
INDUSTRY_STORAGE_TECHNOLOGIES = ["INLHSTOR"]

# Transport bus types and technologies
TRANSPORT_BUS_TYPES = ["TRAN-PUB", "TRAN-PRV"]
TRANSPORT_TECHNOLOGIES = ["EVST-PUB", "EVST-PRV"]

# Demand profile types mapping
PROFILE_MAPPING = {
    "HVELEC": "HV_LOAD",
    "LVELEC": "LV_LOAD",
    "IND-LH": "IND_LOAD",
    "IND-HH": "IND_LOAD",
    "TRAN-PUB": "HPV_LOAD",
    "TRAN-PRV": "LPV_LOAD",
}


def get_carrier(bus_type: str) -> str:
    """Map a bus type to its corresponding energy carrier.

    Parameters
    ----------
    bus_type : str
        Type of bus.

    Returns
    -------
    str
        Corresponding carrier name.
    """
    carrier = None

    if bus_type in [
        "HVELEC",
        "LVELEC",
        "TRAN-PUB",
        "TRAN-PRV",
        "BATSN",
        "HHBSN",
        "HPHSN",
    ]:
        carrier = "Electricity"
    elif bus_type in [
        "LIGN",
        "BITN",
        "HRDCN",
        "OILN",
        "GASN",
        "LNGN",
        "BION",
        "HYDN",
        "CO2STORN",
    ]:
        carrier = bus_type[:-1].capitalize()
    elif bus_type == "NUCLN":
        carrier = "Uranium"
    elif bus_type in ["IND-LH", "INLHSTORN", "INLHSTOR"]:
        carrier = "Low_Heat"
    elif bus_type == "IND-HH":
        carrier = "High_Heat"
    elif bus_type == "ATMP":
        carrier = "CO2"
    elif bus_type in ["PHOT", "CSP", "RTPV", "SWHT"]:
        carrier = "Solar"
    elif bus_type in ["HROR", "HDAM"]:
        carrier = "Water"
    elif bus_type in ["WTON", "WTOF"]:
        carrier = "Wind"
    elif bus_type == "GEOT":
        carrier = "Geothermal"
    elif bus_type == "WSTN":
        carrier = "Waste"
    elif "EVST" in bus_type:
        carrier = "Electricity"
    else:
        carrier = f"{bus_type}'s carrier is not assigned. Please add manually."

    return carrier


def get_bus_type(technology: str) -> str:
    """Map a generator or storage technology to its bus type.

    Parameters
    ----------
    technology : str
        Generator or storage technology.

    Returns
    -------
    str
        Corresponding bus type.
    """
    bus_type = None

    if technology in ["PHOT", "CSP", "HROR", "HDAM", "WTON", "WTOF", "GEOT", "WSTT"]:
        bus_type = "HVELEC"
    elif technology == "RTPV":
        bus_type = "LVELEC"
    elif technology in ["INLHSTORN", "SWHT"]:
        bus_type = "IND-LH"
    elif "EVST" in technology:
        bus_type = f"TRAN-{technology.split('-')[1]}"
    else:
        bus_type = f"{technology}'s bus_type is not assigned. Please add manually"

    return bus_type


def create_buses(
    countries: list,
    nodes: list,
    bus_types: list,
    file_path: FilePath = None,
    return_df: bool = False,
) -> pd.DataFrame:
    """Generate a DataFrame of energy system buses for all nodes and countries.

    Parameters
    ----------
    countries : list
        List of country names
    node_list : list
        List of node names with this format: "country_region"
    bus_types : list
        List of bus types to include.
    file_path : FilePath (str or path object), optional
        If provided, saves the resulting DataFrame as a CSV fiel to this given path,
        by default None
    return_df : bool, optional
        If set to True, returns the DataFrame instead of saving it, by default False

    Returns
    -------
    pd.DataFrame
        DataFrame of buses (only if `return_df` is True).
    """
    # Define which bus types are global (country-level, e.g., fuels, CO2, hydrogen)
    global_bus_types = [
        x
        for x in bus_types
        if x
        in [
            "LIGN",
            "BITN",
            "HRDCN",
            "OILN",
            "GASN",
            "LNGN",
            "BION",
            "ATMP",
            "HYDN",
            "CO2STORN",
            "WSTN",
            "NUCLN",
        ]
    ]
    # Build DataFrame for global buses (one per country)
    global_buses = pd.DataFrame()
    for country in countries:
        df = pd.DataFrame({"bus_type": global_bus_types})
        df["node"] = country
        global_buses = pd.concat([global_buses, df], axis=0)

    # Define which bus types are regional (node-level, e.g., electricity, storage,
    # transport, industry)
    regional_bus_types = [
        x
        for x in bus_types
        if x
        in [
            "HVELEC",
            "LVELEC",
            "BATSN",
            "HHBSN",
            "HPHSN",
            "TRAN-PUB",
            "TRAN-PRV",
            "IND-LH",
            "IND-HH",
            "INLHSTORN",
        ]
    ]
    # Build DataFrame for regional buses (one per node)
    regional_buses = pd.DataFrame()
    regional_buses["node"] = nodes * len(regional_bus_types)
    regional_buses["bus_type"] = np.repeat(regional_bus_types, len(nodes))

    # Combine global and regional buses
    bus_df = pd.concat([global_buses, regional_buses], axis=0)
    bus_df["bus"] = bus_df["node"] + "_" + bus_df["bus_type"]
    bus_df["carrier"] = bus_df["bus_type"].apply(get_carrier)
    bus_df["country"] = bus_df["bus"].apply(lambda x: x.split("_")[0])
    bus_df = bus_df.set_index(["country", "node"]).sort_values(by=["country", "node"])

    result = None
    if return_df:
        result = bus_df
    else:
        # Save to CSV, dropping the bus_type column
        bus_df.drop(columns=["bus_type"]).to_csv(path_or_buf=file_path)

    return result


def create_fuel_supplies(
    countries: list, years: list, file_path: FilePath, currency: str
) -> pd.DataFrame:
    """Generate a DataFrame of fuel supply plants for all countries and years.

    Each row in this DataFrame represents a fuel supply bus for a country and year,
    with columns for:
    - bus: Bus name (country + fuel type)
    - supply_plant: Plant name (e.g., TGEN_LIGN)
    - carrier: Energy carrier (from get_carrier)
    - fuel_cost [CURRENCY/MWh]: Placeholder for fuel cost (to be filled)
    - year: Year
    - country: Country

    Parameters
    ----------
    countries : list
        List of country names
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])

    file_path : FilePath (str or path object)
        Path to save the resulting CSV file.
    currency : str
        Currency for fuel cost (e.g., "EUR", "USD").

    Returns
    -------
    pd.DataFrame
        DataFrame of fuel supply plants (also written to CSV).
    """
    # List of fuel types for which supply plants are created
    fuel_types = [
        "LIGN",
        "BITN",
        "HRDCN",
        "GASN",
        "LNGN",
        "BION",
        "OILN",
        "WSTN",
        "NUCLN",
    ]
    # Build list of all fuel supply buses (one per country and fuel type)
    fuel_supply_buses = [
        f"{country}_{ftype}" for country in countries for ftype in fuel_types
    ]

    # Build DataFrame for all years
    fuel_supply_df = pd.DataFrame()
    for year in years:
        temp = pd.DataFrame()
        temp["bus"] = fuel_supply_buses
        temp["bus_type"] = fuel_types * len(countries)
        temp["supply_plant"] = temp["bus_type"].apply(lambda x: f"TGEN_{x}")
        temp["carrier"] = temp["bus_type"].apply(get_carrier)
        temp["max_supply__mwh_year"] = "Please fill here"
        temp[f"fuel_cost__{str(currency).lower()}_mwh"] = "Please fill here"
        temp["year"] = year
        temp["country"] = temp["bus"].apply(lambda x: x.split("_")[0])
        temp = temp.drop("bus_type", axis=1)
        fuel_supply_df = pd.concat([fuel_supply_df, temp], axis=0)

    # Set index and sort for clarity
    fuel_supply_df = fuel_supply_df.set_index(["country", "bus"]).sort_values(
        by=["country", "bus"]
    )
    fuel_supply_df.to_csv(file_path)


def create_storage_energy(
    countries: list, nodes: list, years: list, sector: str, file_path: FilePath
) -> pd.DataFrame:
    """Create a template DataFrame for storage energy assets in the specified sector.

    For the given sector ("power" or "industry"), this function generates a table of
    storage assets (buses), with columns for storage type, carrier, standing loss,
    extendability, nominal energy, and placeholders for maximum storage per year.
    The resulting DataFrame is saved in the provided file_path.


    Parameters
    ----------
    countries : list
        List of country names
    node_list : list
        List of node names with this format: "country_region"
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    sector : str
        Sector for which to create storage energy ("power" or "industry").
    file_path : FilePath (str or path object)
        Path to save the CSV.

    Returns
    -------
    pd.DataFrame
        DataFrame of storage energy (also written to CSV).

    Raises
    ------
    ValueError
        If sector is not "power" or "industry".
    """
    # Determine storage buses by sector
    if sector.lower() == "power":
        # Global storage buses (country-level, e.g., hydrogen, CO2 storage)
        country_storage_buses = []
        for country in countries:
            country_storage_buses += [
                f"{country}_{bus_type}" for bus_type in ["HYDN", "CO2STORN"]
            ]
        # Regional storage buses (node-level, e.g., batteries, pumped hydro, household
        # batteries)
        node_storage_buses = [
            f"{node}_{bus_type}"
            for node, bus_type in itertools.product(nodes, ["BATSN", "HHBSN", "HPHSN"])
        ]
        all_storage_buses = country_storage_buses + node_storage_buses
    elif sector.lower() == "industry":
        # Only low-heat industry storage (node-level)
        all_storage_buses = [f"{node}_INLHSTORN" for node in nodes]
    else:
        raise ValueError("Sector must be 'power' or 'industry'.")

    # Build DataFrame
    storage_df = pd.DataFrame()
    storage_df["bus"] = all_storage_buses
    storage_df["bus_type"] = storage_df["bus"].apply(lambda x: x.split("_")[-1])
    storage_df["store"] = storage_df["bus"] + "_STOR"
    # Remove trailing 'N' from bus_type for type (e.g., BATSN -> BATS)
    storage_df["type"] = storage_df["bus_type"].apply(lambda x: x[:-1])
    storage_df["carrier"] = storage_df["bus_type"].apply(get_carrier)
    storage_df["standing_loss"] = "Please fill here"
    storage_df["e_nom_extendable"] = "Please fill here"
    storage_df["e_nom"] = "Please fill here"
    storage_df["country"] = storage_df["bus"].apply(lambda x: x.split("_")[0])
    # Add columns for each year
    for year in years:
        storage_df[f"max_store_{year}"] = "Please fill here"
    # Finalize DataFrame: drop bus_type, set index, sort
    storage_df = (
        storage_df.drop("bus_type", axis=1)
        .set_index(["country", "bus"])
        .sort_values(by=["country", "bus"])
    )
    storage_df.to_csv(file_path)
    return storage_df


def create_storage_capacity(
    nodes: list,
    years: list,
    sector: str,
    file_path: FilePath,
    elec_buses: pd.DataFrame,
    ind_buses: pd.DataFrame,
) -> pd.DataFrame:
    """Generate a DataFrame of storage capacities for the given sector and save as CSV.

    For the "Power" sector, this includes:
    - Hydro dam storage (HDAM) on high-voltage electricity buses
    - Battery (BATS) and pumped hydro (HPHS) storage on high-voltage electricity buses
    - Household battery storage (HHBS) on low-voltage electricity buses

    For the "industry" sector, this includes:
    - Low-heat industry storage (INLHSTOR) on low-heat industry buses

    Parameters
    ----------
    node_list : list
        List of node names with this format: "country_region"
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    sector : str
        Either "power" or "industry".
    file_path : FilePath (str or path object)
        Path to save the resulting CSV file.
    elec_buses : pd.DataFrame
        DataFrame of electricity sector buses.
    ind_buses : pd.DataFrame
        DataFrame of industry sector buses.

    Returns
    -------
    pd.DataFrame
        DataFrame of storage capacities (also written to CSV).

    Raises
    ------
    ValueError
        If sector is not "power" or "industry".
    """
    # Combine all buses for lookup
    all_buses = pd.concat([elec_buses, ind_buses], axis=0)

    if sector.lower() == "power":
        # Get high-voltage and low-voltage electricity buses
        hv_buses = all_buses[all_buses.bus_type.str.contains("HVELEC")]["bus"].to_list()
        lv_buses = all_buses[all_buses.bus_type.str.contains("LVELEC")]["bus"].to_list()

        # Hydro dam storage (HDAM) on HV buses
        hydro_df = pd.DataFrame(
            {
                "type": ["HDAM"] * len(hv_buses),
                "node": nodes,
                "carrier": ["Water"] * len(hv_buses),
                "bus": hv_buses,
            }
        )

        # Battery (BATS) and pumped hydro (HPHS) storage on HV buses
        storage_types = ["BATS", "HPHS"]
        storage_rows = []
        for stype in storage_types:
            storage_rows.append(
                pd.DataFrame(
                    {
                        "type": [stype] * len(hv_buses),
                        "node": nodes,
                        "carrier": ["Electricity"] * len(hv_buses),
                        "bus": hv_buses,
                    }
                )
            )
        storage_df = pd.concat(storage_rows, axis=0)

        # Household battery storage (HHBS) on LV buses
        hhbs_df = pd.DataFrame(
            {
                "type": ["HHBS"] * len(lv_buses),
                "node": nodes,
                "carrier": ["Electricity"] * len(lv_buses),
                "bus": lv_buses,
            }
        )

        # Combine all power sector storage
        storage_capacity_df = pd.concat([hydro_df, storage_df, hhbs_df], axis=0)[
            ["type", "node", "carrier", "bus"]
        ]

    elif sector.lower() == "industry":
        # Low-heat industry storage (INLHSTOR) on low-heat industry buses
        lh_buses = all_buses[all_buses.bus_type.str.contains("IND-LH")]["bus"].to_list()
        storage_capacity_df = pd.DataFrame(
            {
                "type": ["INLHSTOR"] * len(lh_buses),
                "node": nodes,
                "carrier": ["Low_Heat"] * len(lh_buses),
                "bus": lh_buses,
            }
        )[["type", "node", "carrier", "bus"]]

    else:
        raise ValueError("Sector must be 'power' or 'industry'.")

    # Add name and country columns
    storage_capacity_df["name"] = (
        storage_capacity_df["bus"] + "_" + storage_capacity_df["type"]
    )
    storage_capacity_df["country"] = storage_capacity_df["bus"].apply(
        lambda x: x.split("_")[0]
    )

    # Add placeholder columns for user input
    fill_columns = (
        ["p_nom", "p_nom_extendable"]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        storage_capacity_df[column] = "Please fill here"

    # Sort and set index for clarity
    storage_capacity_df = storage_capacity_df.sort_values(
        by=["country", "node", "type"]
    ).set_index(["country", "node", "type"])

    storage_capacity_df.to_csv(file_path)
    return storage_capacity_df


def create_loads(
    years: list, bus_df: pd.DataFrame, file_path: FilePath
) -> pd.DataFrame:
    """Create a DataFrame of loads for each bus and year.

    Parameters
    ----------
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    bus_df : pd.DataFrame
        DataFrame of buses.
    file_path : FilePath (str or path object)
        Path to save the CSV.

    Returns
    -------
    pd.DataFrame
        DataFrame of loads.
    """
    load_df = bus_df[
        bus_df["bus_type"].isin(
            ["HVELEC", "LVELEC", "IND-LH", "IND-HH", "TRAN-PUB", "TRAN-PRV"]
        )
    ].copy()
    load_df["profile_type"] = load_df["bus_type"].apply(
        lambda x: PROFILE_MAPPING.get(
            x, f"{x}'s profile is not assigned. Please add manually"
        )
    )
    load_df["total_load__mwh"] = "Please fill here"
    temp = pd.DataFrame()
    for year in years:
        load_df["year"] = year
        temp = pd.concat([temp, load_df], axis=0)
    load_df = temp
    load_df["name"] = load_df["bus"] + "_" + load_df["profile_type"]
    load_df = load_df[
        ["bus", "profile_type", "name", "total_load__mwh", "carrier", "year"]
    ]
    load_df = load_df.sort_values(by=["node", "year"])
    load_df.to_csv(file_path)


def create_generators(
    nodes: list, years: list, technologies: list, file_path: FilePath
) -> pd.DataFrame:
    """Generate a template DataFrame of generators for each technology and node.

    For each node and generator technology, this function creates a row with:
    - type: Generator technology (e.g., PHOT, CSP, HROR, etc.)
    - node: Node name (region)
    - carrier: Energy carrier (from get_carrier)
    - bus: Bus name (node + bus_type)
    - country: Country name (parsed from node)
    - name: Unique generator name (bus + type)
    - p_nom, p_nom_extendable, p_nom_max_YEAR, p_nom_min_YEAR: Placeholders for user
      input

    Parameters
    ----------
    node_list : list
        List of node names with this format: "country_region"
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    technologies : list
        List of generator technologies
    file_path : FilePath (str or path object)
        Path to save the resulting CSV file

    Returns
    -------
    pd.DataFrame
        DataFrame of generators (also written to CSV).
    """
    # Build DataFrame for all combinations of node and technology
    gen_df = pd.DataFrame()
    gen_df["type"] = technologies * len(nodes)
    gen_df["node"] = np.repeat(nodes, len(technologies))
    gen_df["bus_type"] = gen_df["type"].apply(get_bus_type)
    gen_df["carrier"] = gen_df["type"].apply(get_carrier)
    gen_df["bus"] = gen_df["node"] + "_" + gen_df["bus_type"]
    gen_df["country"] = gen_df["node"].apply(lambda x: x.split("_")[0])
    gen_df["name"] = gen_df["bus"] + "_" + gen_df["type"]

    # Add placeholder columns for user input
    fill_columns = (
        ["p_nom", "p_nom_extendable"]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        gen_df[column] = "Please fill here"

    # Finalize DataFrame: drop bus_type, set index, sort
    gen_df = gen_df.drop("bus_type", axis=1)
    gen_df = gen_df.set_index(["country", "node", "type"]).sort_values(
        by=["country", "node", "type"]
    )
    gen_df.to_csv(file_path)
    return gen_df


def create_interconnector(
    nodes: list, years: list, file_path: FilePath, currency: str
) -> pd.DataFrame:
    """Generate transmission and distribution interconnectors for the power sector.

    This function creates:
    - High-voltage (HV) interconnectors between all pairs of nodes(representing
      transmission lines).
    - Low-voltage (LV) interconnectors within each node (representing distribution lines
      from HV to LV).

    The resulting DataFrame includes columns for bus connections, country, link name,
    carrier, type, efficiency, power limits, and placeholders for user input (capacity,
    costs, etc.), and is saved to the specified CSV file.

    Parameters
    ----------
    node_list : list
        List of node names with this format: "country_region"
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Path to save the resulting CSV file.
    currency : str
        Currency for costs (e.g., "EUR", "USD").

    Returns
    -------
    pd.DataFrame
        DataFrame of interconnectors (also written to CSV).
    """
    # =========== High-voltage transmission lines: all unique pairs of nodes ===========
    hv_bus_names = [f"{node}_HVELEC" for node in nodes]
    hv_pairs = list(combinations(hv_bus_names, 2))
    hv_df = pd.DataFrame(hv_pairs, columns=["bus0", "bus1"])

    # =========== Low-voltage distribution lines: within each node (HV to LV) ==========
    lv_df = pd.DataFrame(
        {"bus0": hv_bus_names, "bus1": [f"{node}_LVELEC" for node in nodes]}
    )

    # ======================== Combine HV and LV interconnectors =======================
    intercon_df = pd.concat([hv_df, lv_df], axis=0, ignore_index=True)

    # ============================== Add metadata columns ==============================
    intercon_df["country"] = intercon_df["bus0"].apply(lambda x: x.split("_")[0])
    intercon_df["link"] = intercon_df["bus0"] + "_to_" + intercon_df["bus1"]
    intercon_df["carrier"] = "Electricity"
    intercon_df["type"] = "ITCN"
    intercon_df["efficiency"] = 1
    intercon_df["p_max_pu"] = 1
    # HV lines: allow negative flow; LV lines: only positive flow
    intercon_df["p_min_pu"] = intercon_df["bus1"].apply(
        lambda x: -1 if "HVELEC" in x else 0
    )

    # ===================== Add placeholder columns for user input =====================
    fill_columns = (
        [
            "p_nom",
            "p_nom_extendable",
            f"cap__{str(currency).lower()}_mw",
            f"fom__{str(currency).lower()}_mwa",
            "marginal_cost",
        ]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        intercon_df[column] = "Please fill here"

    # ======================== 🌶 Finalize and save 🌶 ===================
    intercon_df = intercon_df.set_index(["country", "link"]).sort_values(
        by=["bus0", "bus1"]
    )
    intercon_df.to_csv(file_path)


def create_power_links(
    countries: list,
    years: list,
    file_path: FilePath,
    elec_buses: pd.DataFrame,
    ind_buses: pd.DataFrame,
) -> pd.DataFrame:
    """Create power sector converter links (thermal power plants) for all countries.

    For each country, this function generates links representing thermal power plant
    technologies (coal, oil, gas, biomass, hydrogen, waste, nuclear) connecting fuel
    supply buses to high-voltage electricity buses. Emitting technologies are also
    linked to atmosphere and CO2 storage buses.

    Parameters
    ----------
    countries : list
        List of country names
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Output CSV file path.
    elec_buses : pd.DataFrame
        DataFrame of electricity sector buses.
    ind_buses : pd.DataFrame
        DataFrame of industry sector buses.

    Returns
    -------
    pd.DataFrame
        DataFrame of power sector converter links (also written to CSV).
    """
    # Mapping: fuel bus types -> power plant technology types
    fuel_to_tech_map = [
        (["LIGN", "BITN", "HRDCN"], ["SubC", "SupC"]),  # Coal
        (["OILN"], ["OILT"]),  # Oil
        (["GASN", "LNGN"], ["OCGT", "CCGT"]),  # Gas
        (["BION"], ["BIOT"]),  # Biomass
        (["HYDN"], ["OCHT"]),  # Hydrogen
        (["WSTN"], ["WSTT"]),  # Waste
        (["NUCLN"], ["NUCL"]),  # Nuclear
    ]

    # Combine bus DataFrames for lookup
    bus_df = pd.concat([elec_buses, ind_buses], axis=0)

    all_power_links = pd.DataFrame()

    for country in countries:
        # All HV electricity buses for this country
        hv_buses = bus_df[
            (bus_df.bus_type.str.contains("HVELEC"))
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()

        country_links = []

        for fuel_types, tech_types in fuel_to_tech_map:
            # All fuel supply buses of these types for this country
            fuel_buses = bus_df[
                (bus_df.bus_type.isin(fuel_types)) & (bus_df.bus.str.contains(country))
            ]["bus"].to_list()

            # For each HV bus, create all (fuel_bus, tech_type) combinations
            for hv_bus in hv_buses:
                for fuel_bus in fuel_buses:
                    for tech in tech_types:
                        link = {
                            "country": country,
                            "bus0": fuel_bus,
                            "bus1": hv_bus,
                            "type": tech,
                        }
                        country_links.append(link)

        country_links_df = pd.DataFrame(country_links)
        if country_links_df.empty:
            continue

        # Add emission management buses for emitting technologies (all except hydrogen)
        emitting_mask = country_links_df["type"] != "OCHT"
        if emitting_mask.any():
            country_links_df.loc[emitting_mask, "bus2"] = bus_df.loc[
                (bus_df.bus == f"{country}_ATMP"), "bus"
            ].values[0]
            country_links_df.loc[emitting_mask, "bus3"] = bus_df.loc[
                (bus_df.bus == f"{country}_CO2STORN"), "bus"
            ].values[0]

        # Add unique link name
        country_links_df["link"] = (
            country_links_df["bus0"]
            + "_to_"
            + country_links_df["bus1"]
            + "_by_"
            + country_links_df["type"]
        )

        all_power_links = pd.concat([all_power_links, country_links_df], axis=0)

    # Assign carrier based on fuel bus type
    all_power_links["carrier"] = all_power_links["bus0"].apply(
        lambda x: get_carrier(x.split("_")[-1])
    )

    # Add placeholder columns for user input
    fill_columns = (
        ["p_nom", "p_nom_extendable"]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        all_power_links[column] = "Please fill here"

    # Set index and sort for clarity
    index_cols = ["country"] + sorted(
        [col for col in all_power_links.columns if "bus" in col]
    )
    all_power_links = all_power_links.set_index(index_cols, append=False)
    all_power_links = all_power_links.sort_values(
        by=["country", "bus0", "bus1", "type"]
    )

    all_power_links.to_csv(file_path)
    return all_power_links


def create_ind_heat_links(
    countries: list,
    years: list,
    file_path: FilePath,
    elec_buses: pd.DataFrame,
    ind_buses: pd.DataFrame,
) -> pd.DataFrame:
    """Create industry sector heat production/conversion links for each country/region.

    This includes:
    - Fossil, bio, and hydrogen boilers (fuel -> industry heat bus, with CO2 management
      for non-hydrogen)
    - Low-temperature electric heaters (LV electricity -> low heat bus)
    - High-temperature electric heaters (HV electricity -> high heat bus)
    - Combined Heat and Power (CHP) links (bio -> HV electricity + low heat bus, with
      CO2 management)

    Parameters
    ----------
    countries : list
        List of country names
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Output CSV file path.
    elec_buses : pd.DataFrame
        DataFrame of electricity sector buses.
    ind_buses : pd.DataFrame
        DataFrame of industry sector buses.

    Returns
    -------
    pd.DataFrame
        DataFrame of industry sector converter links with columns for capacity (p_nom),
        investment status (p_nom_extendable), year-specific maximum/minimum installed
        capacities (p_nom_max_year, p_nom_min_year), and efficiency placeholders to be
        filled.
    """
    # Combine bus DataFrames for lookup
    bus_df = pd.concat([elec_buses, ind_buses], axis=0)
    all_links = pd.DataFrame()

    for country in countries:
        # Get all high and low heat industry buses for this country
        ind_hh_buses = bus_df[
            (bus_df.bus_type.str.contains("IND-HH"))
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()
        ind_lh_buses = bus_df[
            (bus_df.bus_type.str.contains("IND-LH"))
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()
        all_ind_buses = ind_hh_buses + ind_lh_buses

        # ======== Boiler links: fuel (fossil/bio/hydrogen) to industry heat bus =======
        boiler_fuel_types = ["HRDCN", "GASN", "OILN", "BION", "HYDN"]
        boiler_fuel_buses = bus_df[
            (bus_df.bus_type.isin(boiler_fuel_types))
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()
        boiler_links = pd.DataFrame()
        boiler_links["bus1"] = all_ind_buses * len(boiler_fuel_buses)
        boiler_links["bus0"] = np.repeat(boiler_fuel_buses, len(all_ind_buses))
        # Add emission management buses for non-hydrogen boilers
        non_hydrogen = ~boiler_links.bus0.str.contains("HYDN")
        boiler_links.loc[non_hydrogen, "bus2"] = bus_df.loc[
            (bus_df.bus == f"{country}_ATMP"), "bus"
        ].values[0]
        boiler_links.loc[non_hydrogen, "bus3"] = bus_df.loc[
            (bus_df.bus == f"{country}_CO2STORN"), "bus"
        ].values[0]
        boiler_links["type"] = "IND_BOILER"
        boiler_links["carrier"] = boiler_links["bus0"].apply(
            lambda x: get_carrier(x.split("_")[-1])
        )

        # ========== Low-temp electric heaters: LV electricity to low heat bus =========
        low_elec_types = ["EDLH", "EHPP"]
        low_elec_links = pd.DataFrame()
        low_elec_links["type"] = low_elec_types * len(ind_lh_buses)
        low_elec_links["bus1"] = np.repeat(ind_lh_buses, len(low_elec_types))
        # LV electricity bus for each region
        low_elec_links["bus0"] = (
            low_elec_links["bus1"].str.replace("IND-LH", "", regex=False) + "LVELEC"
        )
        low_elec_links["carrier"] = "Electricity"

        # ========= High-temp electric heaters: HV electricity to high heat bus ========
        high_elec_types = ["EIDT", "EERH"]
        high_elec_links = pd.DataFrame()
        high_elec_links["type"] = high_elec_types * len(ind_hh_buses)
        high_elec_links["bus1"] = np.repeat(ind_hh_buses, len(high_elec_types))
        # HV electricity bus for each region
        high_elec_links["bus0"] = (
            high_elec_links["bus1"].str.replace("IND-HH", "", regex=False) + "HVELEC"
        )
        high_elec_links["carrier"] = "Electricity"

        # == CHP links: bio bus to HV elec and low heat bus, with emission management ==
        chp_links = pd.DataFrame()
        chp_links["type"] = ["CHP"] * len(ind_lh_buses)
        chp_links["bus3"] = np.repeat(ind_lh_buses, 1)
        chp_links["bus1"] = (
            chp_links["bus3"].str.replace("IND-LH", "", regex=False) + "HVELEC"
        )
        chp_links["bus2"] = bus_df.loc[(bus_df.bus == f"{country}_ATMP"), "bus"].values[
            0
        ]
        chp_links["bus0"] = bus_df.loc[(bus_df.bus == f"{country}_BION"), "bus"].values[
            0
        ]
        chp_links["carrier"] = "Bio"

        # Combine all link types for this country
        country_links = pd.concat(
            [boiler_links, low_elec_links, high_elec_links, chp_links], axis=0
        )
        country_links["country"] = country
        all_links = pd.concat([all_links, country_links], axis=0)

    # Create unique link name and add placeholder columns for capacity and investment
    # status
    all_links["link"] = (
        all_links["bus0"] + "_to_" + all_links["bus1"] + "_by_" + all_links["type"]
    )
    fill_columns = (
        ["p_nom", "p_nom_extendable"]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        all_links[column] = "Please fill here"

    # Set index for easier lookup and save to CSV
    all_links = all_links.set_index(
        ["country", "link"] + sorted([col for col in all_links.columns if "bus" in col])
    )
    all_links.to_csv(file_path)


def create_fuel_conversion_links(
    countries: list,
    years: list,
    file_path: FilePath,
    elec_buses: pd.DataFrame,
    ind_buses: pd.DataFrame,
) -> pd.DataFrame:
    """Create synthetic fuel production links for each country.

    This includes:
    - electrolysers (electricity to hydrogen)
    - Fischer-Tropsch (hydrogen + CO2 to synthetic oil), and
    - methanation (hydrogen + CO2 to synthetic gas)

    Parameters
    ----------
    countries : list
        List of country names
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Output CSV file path.
    elec_buses : pd.DataFrame
        Electricity sector buses.
    ind_buses : pd.DataFrame
        industry sector buses.

    Returns
    -------
    pd.DataFrame
        DataFrame of synthetic fuel conversion links (also written to CSV).
    """
    bus_df = pd.concat([elec_buses, ind_buses], axis=0)
    all_links = pd.DataFrame()

    for country in countries:
        # Hydrogen bus for this country
        hyd_bus = bus_df[
            (bus_df.bus_type.str.contains("HYDN")) & (bus_df.bus.str.contains(country))
        ]["bus"].values[0]

        # Electrolyser links: HV electricity buses to hydrogen bus
        hv_elec_buses = bus_df[
            (bus_df.bus_type.str.contains("HVELEC"))
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()
        eltz_links = pd.DataFrame(
            {
                "bus0": hv_elec_buses,
                "bus1": hyd_bus,
                "type": "ELTZ",
                "carrier": "Hyd",
                # Electrolyser efficiency (electricity to hydrogen)
                # "efficiency": 0.7, # noqa E800
            }
        )

        # Fischer-Tropsch link: hydrogen + CO2 to synthetic oil
        fitr_link = pd.Series(
            {
                "bus0": hyd_bus,
                "bus1": f"{country}_OILN",
                "bus2": f"{country}_CO2STORN",
                "type": "FITR",
                "carrier": "Oil",
                # Fischer-Tropsch efficiency (hydrogen + CO2 to oil)
                # "efficiency": 0.5, # noqa E800
            }
        )

        # Methanation link: hydrogen + CO2 to synthetic gas
        meth_link = pd.Series(
            {
                "bus0": hyd_bus,
                "bus1": f"{country}_GASN",
                "bus2": f"{country}_CO2STORN",
                "bus3": np.NAN,
                "type": "METH",
                "carrier": "Gas",
                # Methanation efficiency (hydrogen + CO2 to gas)
                # "efficiency": 0.6, # noqa E800
            }
        )

        country_links = pd.concat(
            [eltz_links, pd.DataFrame([fitr_link, meth_link])], axis=0
        )
        country_links["country"] = country
        all_links = pd.concat([all_links, country_links], axis=0)

    all_links["link"] = (
        all_links["bus0"] + "_to_" + all_links["bus1"] + "_by_" + all_links["type"]
    )

    fill_columns = (
        ["p_nom", "p_nom_extendable"]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        all_links[column] = "Please fill here"

    all_links = all_links.set_index(
        ["country", "link"] + sorted([col for col in all_links.columns if "bus" in col])
    )
    all_links.to_csv(file_path)


def create_dac_links(
    countries: list,
    years: list,
    file_path: FilePath,
    elec_buses: pd.DataFrame,
    ind_buses: pd.DataFrame,
) -> pd.DataFrame:
    """Create Direct Air Capture link for each country and high-voltage electricity bus.

    Each Direct Air Capture (DAC) link connects:
    - bus0: Atmosphere bus (ATMP)
    - bus1: CO2 storage bus (CO2STORN)
    - bus2: High-voltage electricity bus (HVELEC) (electricity input)
    - bus3: (unused, set to NaN)

    The resulting DataFrame is saved to file_path, with columns for capacity (p_nom),
    investment status (p_nom_extendable), and year-specific maximum/minimum installed
    capacities (p_nom_max_year, p_nom_min_year) left to be filled.

    Parameters
    ----------
    countries : list
        List of country names
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Output CSV file path.
    elec_buses : pd.DataFrame
        DataFrame of electricity sector buses.
    ind_buses : pd.DataFrame
        DataFrame of industry sector buses.

    Returns
    -------
    pd.DataFrame
        DataFrame of DAC links (also written to CSV).
    """
    # Combine bus DataFrames for lookup
    bus_df = pd.concat([elec_buses, ind_buses], axis=0)

    dac_links = pd.DataFrame()
    for country in countries:
        # Get all HV electricity buses for this country
        hv_elec_buses = bus_df[
            bus_df.bus_type.str.contains("HVELEC") & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()

        # Atmosphere and CO2 storage buses for this country
        atmp_bus = bus_df.loc[(bus_df.bus == f"{country}_ATMP"), "bus"].values[0]
        co2stor_bus = bus_df.loc[(bus_df.bus == f"{country}_CO2STORN"), "bus"].values[0]

        # Create DAC links for each HV bus
        country_dac_df = pd.DataFrame(
            {
                "bus0": atmp_bus,  # Atmosphere
                "bus1": co2stor_bus,  # CO2 storage
                "bus2": hv_elec_buses,  # Electricity input
                "bus3": np.NAN,  # Unused
                "type": "DAC",
                "carrier": "CO2",
                # "efficiency": 1,          # Efficiency of CO2 capture/storage
                # "efficiency2": -0.5,      # Electricity consumed per unit CO2 captured
            }
        )
        country_dac_df["country"] = country
        dac_links = pd.concat([dac_links, country_dac_df], axis=0)

    # Create unique link name
    dac_links["link"] = dac_links["type"] + "_at_" + dac_links["bus2"]

    # Optionally, set defaulted settings for power limits
    # dac_links['p_max_pu'] = 1 # noqa E800
    # dac_links['p_min_pu'] = 0 # noqa E800

    # Add columns to be filled by user
    fill_columns = (
        ["p_nom", "p_nom_extendable"]
        + [f"p_nom_max_{year}" for year in years]
        + [f"p_nom_min_{year}" for year in years]
    )
    for column in fill_columns:
        dac_links[column] = "Please fill here"

    dac_links = dac_links.set_index(
        ["country", "link"] + sorted([col for col in dac_links.columns if "bus" in col])
    )
    dac_links.to_csv(file_path)


def create_ev_chargers(
    countries: list,
    years: list,
    file_path: FilePath,
    elec_buses: pd.DataFrame,
    tra_buses: pd.DataFrame,
) -> pd.DataFrame:
    """Create EV charger links between electricity and transport buses for each country.

    For each country:
    - Public EV chargers: Connect low voltage electricity buses (LVELEC) to public
      transport EV buses (TRAN-PUB).
    - Private EV chargers: Connect low voltage electricity buses (LVELEC) to private
      transport EV buses (TRAN-PRV).

    The resulting DataFrame is saved to file_path, with columns for charger numbers
    (num_ch_YEAR) left to be filled.

    Parameters
    ----------
    countries : list
        List of country names
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Output CSV file path.
    elec_buses : pd.DataFrame
        DataFrame of electricity sector buses.
    tra_buses : pd.DataFrame
        DataFrame of transport sector buses.

    Returns
    -------
    pd.DataFrame
        DataFrame of EV charger links (also written to CSV).
    """
    bus_df = pd.concat([elec_buses, tra_buses], axis=0)
    all_ev_links = pd.DataFrame()

    for country in countries:
        # Get all low voltage electricity buses for this country
        lv_elec_buses = bus_df[
            bus_df.bus_type.str.contains("LVELEC") & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()

        # Get all public EV buses for this country
        public_ev_buses = bus_df[
            bus_df.bus_type.str.contains("TRAN-PUB")
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()

        # Create links for public EV chargers
        public_ev_links = pd.DataFrame(
            {
                "bus0": sorted(lv_elec_buses),
                "bus1": sorted(public_ev_buses),
                "type": "EVCH-PUB",
                "carrier": "Electricity",
                "p_max_pu": "EVCH-PUB",
            }
        )

        # Get all private EV buses for this country
        private_ev_buses = bus_df[
            bus_df.bus_type.str.contains("TRAN-PRV")
            & (bus_df.bus.str.contains(country))
        ]["bus"].to_list()

        # Create links for private EV chargers
        private_ev_links = pd.DataFrame(
            {
                "bus0": sorted(lv_elec_buses),
                "bus1": sorted(private_ev_buses),
                "type": "EVCH-PRV",
                "carrier": "Electricity",
                "p_max_pu": "EVCH-PRV",
            }
        )

        # Combine public and private EV charger links for this country
        country_ev_links = pd.concat([public_ev_links, private_ev_links], axis=0)
        country_ev_links["country"] = country
        all_ev_links = pd.concat([all_ev_links, country_ev_links], axis=0)

    # Create unique link name and add columns to be filled
    all_ev_links["link"] = (
        all_ev_links["bus0"]
        + "_to_"
        + all_ev_links["bus1"]
        + "_by_"
        + all_ev_links["type"]
    )
    fill_columns = [f"num_ch_{year}" for year in years]
    for column in fill_columns:
        all_ev_links[column] = "Please fill here"

    all_ev_links = all_ev_links.set_index(
        ["country", "link"]
        + sorted([col for col in all_ev_links.columns if "bus" in col])
    )
    all_ev_links.to_csv(file_path)


def create_ev_storages(
    nodes: list, years: list, technologies_list: list, file_path: FilePath
) -> pd.DataFrame:
    """Create a DataFrame for plug-in electric vehicle (PEV) storage units.

    For each node and transport storage technologies, this function generates a row with
    columns for:
    - type: Storage technologies (e.g., EVST-PUB, EVST-PRV)
    - node: Node name (region)
    - carrier: Energy carrier (from get_carrier)
    - bus: Bus name (node + bus_type)
    - name: Unique storage name (bus + type)
    - num_ev_YEAR: Placeholder columns for number of EVs per year (to be filled)

    Parameters
    ----------
    node_list : list
        List of node names with this format: "country_region"
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    technologies_list : list
        List of transport storage technologies (e.g., ["EVST-PUB", "EVST-PRV"])
    file_path : FilePath (str or path object)
        Path to save the resulting CSV file.

    Returns
    -------
    pd.DataFrame
        DataFrame of PEV storages (also written to CSV).
    """
    # Build DataFrame for all combinations of node and technologies
    storage_df = pd.DataFrame()
    storage_df["type"] = technologies_list * len(nodes)
    storage_df["node"] = np.repeat(nodes, len(technologies_list))
    storage_df["bus_type"] = storage_df["type"].apply(get_bus_type)
    storage_df["carrier"] = storage_df["type"].apply(get_carrier)
    storage_df["bus"] = storage_df["node"] + "_" + storage_df["bus_type"]
    storage_df["name"] = storage_df["bus"] + "_" + storage_df["type"]
    storage_df["country"] = storage_df["node"].apply(lambda x: x.split("_")[0])

    # Add placeholder columns for number of EVs per year
    for year in years:
        storage_df[f"num_ev_{year}"] = "Please fill here"

    storage_df = storage_df.drop("bus_type", axis=1)
    storage_df = storage_df.set_index(["country", "node", "type"]).sort_values(
        by=["country", "node", "type"]
    )
    # Save to CSV
    storage_df.to_csv(file_path)

    return storage_df


def create_decom_csv(years: list, file_path: FilePath) -> pd.DataFrame:
    """Create a template CSV for decommissioned capacity by year.

    The CSV will have the following columns:
    - name: Asset name
    - class: Asset class/type
    - One column per year (excluding the first year in year_list)

    Parameters
    ----------
    year_list : list
        List of years (e.g., [2020, 2030, 2040, 2050])
    file_path : FilePath (str or path object)
        Path to save the CSV file.

    Returns
    -------
    pd.DataFrame
        The empty decommissioning DataFrame (also written to CSV).
    """
    # Only include years after the first year for decommissioning columns
    year_columns = [year for year in years if year != years[0]]
    columns = ["country", "name", "class"] + year_columns
    decom_df = pd.DataFrame(columns=columns)
    decom_df.to_csv(file_path, index=False)
    return decom_df


def update_ev_param(
    save_path: FilePath,
    countries: list,
):
    """Add countries into the ev_parameters csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    countries : list
        List of country names
    """
    ev_param = pd.read_csv(save_path)
    final_ev_df = pd.DataFrame()
    for country in countries:
        country_ev_df = ev_param.copy()
        country_ev_df["country"] = [country] * len(ev_param)
        final_ev_df = pd.concat([final_ev_df, country_ev_df], axis=0)

    final_ev_df.to_csv(save_path, index=False)


def update_availability(
    save_path: FilePath,
    countries: list,
    nodes: list,
):
    """Add countries and nodes into the availability csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    countries : list
        List of country names
    nodes : list
        List of node names with this format: "country_region"
    """
    availability_df = pd.read_csv(save_path)

    for country in countries:
        node_tech_list = [
            "EVST-PRV",
            "EVST-PUB",
            "EVCH-PUB",
            "EVCH-PRV",
        ]
        country_availability_df = pd.DataFrame(columns=availability_df.columns)
        country_availability_df["country"] = [country] * len(node_tech_list)
        country_availability_df["node"] = [country] * len(node_tech_list)
        country_availability_df["technology"] = node_tech_list
        availability_df = pd.concat([country_availability_df, availability_df], axis=0)

    for node in nodes:
        node_tech_list = ["PHOT", "WTON", "RTPV", "WTOF", "CSP", "SWHT", "HROR"]
        country_availability_df = pd.DataFrame(columns=availability_df.columns)
        country_availability_df["node"] = [node] * len(node_tech_list)
        country_availability_df["country"] = country_availability_df["node"].apply(
            lambda x: x.split("_")[0]
        )
        country_availability_df["technology"] = node_tech_list
        availability_df = pd.concat([country_availability_df, availability_df], axis=0)

    availability_df = availability_df.fillna("Please fill here")
    availability_df.to_csv(save_path, index=False)


def update_demand_profiles(
    save_path: FilePath,
    countries: list,
    nodes: list,
):
    """Add countries and nodes into the demand_profile csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    countries : list
        List of country names
    nodes : list
        List of node names with this format: "country_region"
    """
    demand_df = pd.read_csv(save_path)
    final_demand_df = pd.DataFrame()
    for country in countries:
        profile_list = ["IND_LOAD", "HPV_LOAD", "LPV_LOAD"]
        country_demand_df = demand_df.copy()
        country_demand_df["country"] = [country] * len(profile_list)
        country_demand_df["node"] = [country] * len(profile_list)
        country_demand_df["profile_type"] = profile_list
        final_demand_df = pd.concat([final_demand_df, country_demand_df], axis=0)

    for node in nodes:
        profile_list = ["HV_LOAD", "LV_LOAD"]
        node_demand_df = demand_df.copy()
        node_demand_df["node"] = [node] * len(profile_list)
        node_demand_df["country"] = node_demand_df["node"].apply(
            lambda x: x.split("_")[0]
        )
        node_demand_df["profile_type"] = profile_list
        final_demand_df = pd.concat([final_demand_df, node_demand_df], axis=0)

    final_demand_df = final_demand_df.fillna("Please fill here")
    final_demand_df.to_csv(save_path, index=False)


def update_power_plant_costs(save_path: FilePath, countries: list, currency: str):
    """Add countries into the power_plant_costs csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    countries : list
        List of country names
    currency : str
        Currency for costs (e.g., "EUR", "USD").
    """
    power_plant_df = pd.read_csv(save_path)
    final_power_plant_df = pd.DataFrame()
    currency = str(currency).lower()
    for country in countries:
        country_power_plant_df = power_plant_df.copy()
        country_power_plant_df["country"] = country
        country_power_plant_df[f"cap__{currency}_mw"] = "Please fill here"
        country_power_plant_df[f"fom__{currency}_mwa"] = "Please fill here"
        country_power_plant_df[f"vom__{currency}_mwh"] = "Please fill here"
        country_power_plant_df["life__years"] = "Please fill here"
        final_power_plant_df = pd.concat(
            [final_power_plant_df, country_power_plant_df], axis=0
        )

    final_power_plant_df.to_csv(save_path, index=False)


def update_technical_potentials(
    save_path: FilePath,
    nodes: list,
):
    """Add nodes into the renewables_technical_potential csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    nodes : list
        List of node names with this format: "country_region"
    """
    potentials_df = pd.read_csv(save_path)

    for node in nodes:
        node_tech_list = ["PHOT", "WTON", "RTPV", "WTOF", "HROR"]
        country_potentials_df = pd.DataFrame(columns=potentials_df.columns)
        country_potentials_df["class"] = ["Generator"] * len(node_tech_list)
        country_potentials_df["node"] = [node] * len(node_tech_list)
        country_potentials_df["country"] = country_potentials_df["node"].apply(
            lambda x: x.split("_")[0]
        )
        country_potentials_df["type"] = node_tech_list
        country_potentials_df["bus"] = country_potentials_df["node"] + "_HVELEC"
        country_potentials_df["name"] = (
            country_potentials_df["bus"] + "_" + country_potentials_df["type"]
        )
        potentials_df = pd.concat([potentials_df, country_potentials_df], axis=0)

    potentials_df = potentials_df.fillna("Please fill here")
    potentials_df.to_csv(save_path, index=False)


def update_storage_costs(save_path: FilePath, countries: list, currency: str):
    """Add countries into the storage_costs csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    countries : list
        List of country names
    currency : str
        Currency for costs (e.g., "EUR", "USD").
    """
    storage_cost_df = pd.read_csv(save_path)
    final_storage_cost_df = pd.DataFrame()
    currency = str(currency).lower()
    for country in countries:
        country_storage_cost_df = storage_cost_df.copy()
        country_storage_cost_df["country"] = country
        country_storage_cost_df[f"cap__{currency}_mwh"] = "Please fill here"
        country_storage_cost_df[f"fom__{currency}_mwh"] = "Please fill here"
        country_storage_cost_df[f"vom__{currency}_mwh"] = "Please fill here"
        country_storage_cost_df["life__years"] = "Please fill here"
        final_storage_cost_df = pd.concat(
            [final_storage_cost_df, country_storage_cost_df], axis=0
        )

    final_storage_cost_df.to_csv(save_path, index=False)


def update_storage_inflows(
    save_path: FilePath,
    nodes: list,
):
    """Add countries and nodes into the storage_inflows csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    nodes : list
        List of node names with this format: "country_region"
    """
    inflows_df = pd.read_csv(save_path)

    for node in nodes:
        country_inflows_df = pd.DataFrame(columns=inflows_df.columns)
        country_inflows_df["node"] = [node]
        country_inflows_df["country"] = country_inflows_df["node"].apply(
            lambda x: x.split("_")[0]
        )
        country_inflows_df["technology"] = ["HDAM"]
        inflows_df = pd.concat([inflows_df, country_inflows_df], axis=0)

    inflows_df = inflows_df.fillna("Please fill here")
    inflows_df.to_csv(save_path, index=False)


def update_technologies(
    save_path: FilePath,
    countries: list,
):
    """Add countries into the technologies csv.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    countries : list
        List of country names
    """
    technologies = pd.read_csv(save_path)
    final_technologies = pd.DataFrame()
    for country in countries:
        country_technologies = technologies.copy()
        country_technologies["country"] = country
        final_technologies = pd.concat(
            [final_technologies, country_technologies], axis=0
        )

    final_technologies.to_csv(save_path, index=False)


def create_folders(save_path: FilePath, sector_folder: bool = False):
    """Create the main output directory and subdirectories for each sector.

    If the directory already exists, it will be deleted and recreated.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to create.
    sector_folder : bool, optional
        Whether to create subfolders for each sector. Default is False.
    """
    # Create subfolders for each sector
    if sector_folder:
        for sector in ["power", "industry", "transport"]:
            os.makedirs(os.path.join(save_path, sector), exist_ok=True)
    else:  # general input folder structure
        # Remove existing directory if present, then create a fresh one
        if os.path.exists(save_path):
            shutil.rmtree(save_path)
        os.makedirs(save_path, exist_ok=True)
        os.makedirs(os.path.join(save_path, "input"), exist_ok=True)
        os.makedirs(os.path.join(save_path, "input", "global_input"), exist_ok=True)


def check_and_create_global_template_csv(
    save_path: FilePath, global_template_path: FilePath
):
    """Compare or clone each global template csv in the targeted folder.

    Parameters
    ----------
    save_path : FilePath (str or path object)
        Root directory to check or create.
    global_template_path : FilePath (str or path object)
        Root directory of global csv templates.
    """
    for filename in os.listdir(global_template_path):
        if filename.endswith(".csv"):
            source_path = os.path.join(global_template_path, filename)
            target_path = os.path.join(save_path, filename)
            if not os.path.exists(target_path):
                shutil.copy2(source_path, target_path)


if __name__ == "__main__":
    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake  # pylint: disable=ungrouped-imports

        snakemake = mock_snakemake("build_skeleton")
    # Getting global config params
    configure_logging(snakemake)
    configurations = snakemake.params.config

    # Country names from the configuration file
    cfg_countries = configurations["base_configs"]["regions"].keys()
    cfg_nodes = []
    for cfg_country in cfg_countries:
        # List of regions based on the given inputs in the configuration file
        cfg_regions = configurations["base_configs"]["regions"][cfg_country]
        # List of nodes based on the configuration file
        cfg_country_node_list = [f"{cfg_country}_{x}" for x in cfg_regions]
        cfg_nodes += cfg_country_node_list
    # List of years
    cfg_years = configurations["base_configs"]["years"]
    global_csv_skeletons_path = "data/global_input_template"

    # Skeleton porject folder path
    project_folder_path = (
        "data/"
        + configurations["path_configs"]["data_folder_name"]
        + "/"
        + configurations["path_configs"]["project_name"]
    )

    # create_skeleton_inputs
    global_input_path = project_folder_path + "/input/global_input/"
    create_folders(save_path=project_folder_path)
    check_and_create_global_template_csv(
        save_path=global_input_path,
        global_template_path=global_csv_skeletons_path,
    )

    # Skeleton input scenario folder path
    input_scenario_folder_path = (
        project_folder_path
        + "/input/"
        + configurations["path_configs"]["input_scenario_name"]
    )
    cfg_currency = configurations["base_configs"]["currency"]

    # output paths for global templates csvs
    path_availability = global_input_path + "/availability.csv"
    path_demand_profiles = global_input_path + "/demand_profile.csv"
    path_ev_params = global_input_path + "/ev_parameters.csv"
    path_power_plant_cost = global_input_path + "/power_plant_costs.csv"
    path_potentials = global_input_path + "/renewables_technical_potential.csv"
    path_storage_costs = global_input_path + "/storage_costs.csv"
    path_storage_inflows = global_input_path + "/storage_inflows.csv"
    path_technologies = global_input_path + "/technologies.csv"

    # # Creating global template CSVs and save to defined output paths
    update_availability(
        save_path=path_availability,
        countries=cfg_countries,
        nodes=cfg_nodes,
    )
    update_demand_profiles(
        save_path=path_demand_profiles,
        countries=cfg_countries,
        nodes=cfg_nodes,
    )
    update_ev_param(
        save_path=path_ev_params,
        countries=cfg_countries,
    )
    update_power_plant_costs(
        save_path=path_power_plant_cost,
        countries=cfg_countries,
        currency=cfg_currency,
    )
    update_technical_potentials(
        save_path=path_potentials,
        nodes=cfg_nodes,
    )
    update_storage_costs(
        save_path=path_storage_costs,
        countries=cfg_countries,
        currency=cfg_currency,
    )
    update_storage_inflows(
        save_path=path_storage_inflows,
        nodes=cfg_nodes,
    )
    update_technologies(
        save_path=path_technologies,
        countries=cfg_countries,
    )

    # create_skeleton_inputs
    create_folders(save_path=input_scenario_folder_path, sector_folder=True)

    # output paths for regional templates csvs
    path_p_buses = input_scenario_folder_path + "/power/buses.csv"
    path_p_fuel_supplies = input_scenario_folder_path + "/power/fuel_supplies.csv"
    path_p_storage_energy = input_scenario_folder_path + "/power/storage_energy.csv"
    path_p_loads = input_scenario_folder_path + "/power/loads.csv"
    path_p_generators = input_scenario_folder_path + "/power/power_generators.csv"
    path_p_storage_capacity = input_scenario_folder_path + "/power/storage_capacity.csv"
    path_p_interconnector = input_scenario_folder_path + "/power/interconnector.csv"
    path_p_links = input_scenario_folder_path + "/power/power_links.csv"
    path_p_decom_capacity = (
        input_scenario_folder_path + "/power/decomission_capacity.csv"
    )

    path_i_buses = input_scenario_folder_path + "/industry/buses.csv"
    path_i_storage_energy = input_scenario_folder_path + "/industry/storage_energy.csv"
    path_i_loads = input_scenario_folder_path + "/industry/loads.csv"
    path_i_heat_generators = (
        input_scenario_folder_path + "/industry/heat_generators.csv"
    )
    path_i_storage_capacity = (
        input_scenario_folder_path + "/industry/storage_capacity.csv"
    )
    path_i_heat_links = input_scenario_folder_path + "/industry/heat_links.csv"
    path_i_fuel_conversion = (
        input_scenario_folder_path + "/industry/fuel_conversion.csv"
    )
    path_i_dac = input_scenario_folder_path + "/industry/direct_air_capture.csv"
    path_i_decom_capacity = (
        input_scenario_folder_path + "/industry/decomission_capacity.csv"
    )

    path_t_buses = input_scenario_folder_path + "/transport/buses.csv"
    path_t_loads = input_scenario_folder_path + "/transport/loads.csv"
    path_t_pev_storages = input_scenario_folder_path + "/transport/pev_storages.csv"
    path_t_pev_chargers = input_scenario_folder_path + "/transport/pev_chargers.csv"

    # Creating regional template CSVs and save to defined output paths
    # ====================================== Power =====================================
    elec_buses_df = create_buses(
        countries=cfg_countries,
        nodes=cfg_nodes,
        bus_types=POWER_BUS_TYPES,
        return_df=True,
    )
    ind_buses_df = create_buses(
        countries=cfg_countries,
        nodes=cfg_nodes,
        bus_types=INDUSTRY_BUS_TYPES,
        return_df=True,
    )
    tra_buses_df = create_buses(
        countries=cfg_countries,
        nodes=cfg_nodes,
        bus_types=TRANSPORT_BUS_TYPES,
        return_df=True,
    )
    create_buses(
        countries=cfg_countries,
        nodes=cfg_nodes,
        bus_types=POWER_BUS_TYPES,
        file_path=path_p_buses,
    )
    create_fuel_supplies(
        countries=cfg_countries,
        years=cfg_years,
        file_path=path_p_fuel_supplies,
        currency=cfg_currency,
    )
    create_storage_energy(
        countries=cfg_countries,
        nodes=cfg_nodes,
        years=cfg_years,
        sector="Power",
        file_path=path_p_storage_energy,
    )
    create_storage_capacity(
        nodes=cfg_nodes,
        years=cfg_years,
        sector="Power",
        file_path=path_p_storage_capacity,
        elec_buses=elec_buses_df,
        ind_buses=ind_buses_df,
    )
    create_power_links(
        countries=cfg_countries,
        years=cfg_years,
        file_path=path_p_links,
        elec_buses=elec_buses_df,
        ind_buses=ind_buses_df,
    )
    create_loads(
        years=cfg_years,
        bus_df=elec_buses_df,
        file_path=path_p_loads,
    )
    create_generators(
        nodes=cfg_nodes,
        years=cfg_years,
        technologies=POWER_TECHNOLOGIES,
        file_path=path_p_generators,
    )
    create_interconnector(
        nodes=cfg_nodes,
        years=cfg_years,
        file_path=path_p_interconnector,
        currency=cfg_currency,
    )
    create_decom_csv(years=cfg_years, file_path=path_p_decom_capacity)

    # ==================================== Industry ====================================
    create_buses(
        countries=cfg_countries,
        nodes=cfg_nodes,
        bus_types=INDUSTRY_BUS_TYPES,
        file_path=path_i_buses,
    )
    create_storage_energy(
        countries=cfg_countries,
        nodes=cfg_nodes,
        years=cfg_years,
        sector="industry",
        file_path=path_i_storage_energy,
    )
    create_storage_capacity(
        nodes=cfg_nodes,
        years=cfg_years,
        sector="industry",
        file_path=path_i_storage_capacity,
        elec_buses=elec_buses_df,
        ind_buses=ind_buses_df,
        # technologies_list=ind_storage_technologies_list, # noqa E800
    )
    create_loads(
        years=cfg_years,
        bus_df=create_buses(
            countries=cfg_countries,
            nodes=cfg_nodes,
            bus_types=INDUSTRY_BUS_TYPES,
            return_df=True,
        ),
        file_path=path_i_loads,
    )
    create_generators(
        nodes=cfg_nodes,
        years=cfg_years,
        technologies=INDUSTRY_TECHNOLOGIES,
        file_path=path_i_heat_generators,
    )
    create_ind_heat_links(
        countries=cfg_countries,
        years=cfg_years,
        file_path=path_i_heat_links,
        elec_buses=elec_buses_df,
        ind_buses=ind_buses_df,
    )
    create_fuel_conversion_links(
        countries=cfg_countries,
        years=cfg_years,
        file_path=path_i_fuel_conversion,
        elec_buses=elec_buses_df,
        ind_buses=ind_buses_df,
    )
    create_dac_links(
        countries=cfg_countries,
        years=cfg_years,
        file_path=path_i_dac,
        elec_buses=elec_buses_df,
        ind_buses=ind_buses_df,
    )
    create_decom_csv(years=cfg_years, file_path=path_i_decom_capacity)

    # ==================================== Transport ===================================

    bus_df = create_buses(
        countries=cfg_countries,
        nodes=cfg_nodes,
        bus_types=TRANSPORT_BUS_TYPES,
        file_path=path_t_buses,
    )
    create_loads(
        years=cfg_years,
        bus_df=create_buses(
            countries=cfg_countries,
            nodes=cfg_nodes,
            bus_types=TRANSPORT_BUS_TYPES,
            return_df=True,
        ),
        file_path=path_t_loads,
    )
    create_ev_storages(
        nodes=cfg_nodes,
        years=cfg_years,
        technologies_list=TRANSPORT_TECHNOLOGIES,
        file_path=path_t_pev_storages,
    )
    create_ev_chargers(
        countries=cfg_countries,
        years=cfg_years,
        file_path=path_t_pev_chargers,
        elec_buses=elec_buses_df,
        tra_buses=tra_buses_df,
    )

    # =============== Create a scenario config file in the scenario folder =============
    build_scenario_config_file(configurations)
