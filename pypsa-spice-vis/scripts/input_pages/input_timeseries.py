# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Timeseries page under Input section."""

import os

import pandas as pd
import streamlit as st

from scripts.data_utils import render_countries_n_scenario_pills
from scripts.input_st_handler import DataFrameWidgetsHandler


def get_input_scenario_list(base_config: dict) -> list[str]:
    """Return available input scenarios with the default scenario first."""
    data_folder_path = base_config["input_folder_path"]

    if not os.path.exists(data_folder_path):
        raise FileNotFoundError(f"folder not found: {data_folder_path}")

    scenario_list = [
        scenario
        for scenario in os.listdir(data_folder_path)
        if scenario not in [".DS_Store"] and scenario != "global_input"
    ]

    for sce in (base_config["path_configs"]["input_scenario_name"], ""):
        if sce in scenario_list:
            scenario_list.insert(0, scenario_list.pop(scenario_list.index(sce)))

    return scenario_list


def get_mapping_list(*dfs: pd.DataFrame) -> list[str]:
    """Get sorted technology/profile types from one or more dataframes."""
    type_set = set()

    for df in dfs:
        if "technology" in df.columns:
            type_set |= set(df["technology"].unique())
        if "profile_type" in df.columns:
            type_set |= set(df["profile_type"].unique())

    return sorted(type_set)


if __name__ == "__main__":
    base_config = st.session_state.base_config
    selected_sector = st.session_state.get("selected_input_sector", "Power")
    sector_lower = selected_sector.lower()

    tech_config = st.session_state.input_config["Static"]
    time_series_config = st.session_state.input_config["Time_series"]

    df_widgets_handler = DataFrameWidgetsHandler(input_config=time_series_config)

    all_countries = base_config["base_configs"]["regions"].keys()

    st.subheader(f":material/timeline: Timeseries | {selected_sector}")

    dfs = df_widgets_handler.load_all_dfs()

    selected_countries, selected_scenario = render_countries_n_scenario_pills(
        scenario_options=get_input_scenario_list(base_config),
        all_countries=all_countries,
        key=f"{sector_lower}_timeseries_pills",
    )

    # Reload dfs after the scenario is selected
    dfs = df_widgets_handler.load_all_dfs(selected_scenario=selected_scenario)

    tech_path = os.path.join(
        df_widgets_handler.base_input_path,
        "global_input",
        tech_config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)

    types = get_mapping_list(tech_df)
    tech_mapping = dict(zip(tech_df["technology"], tech_df["technology_nomenclature"]))
    types_full_names = [tech_mapping.get(t, t) for t in types]
    types_full_names.sort()
    reverse_mapping = {v: k for k, v in tech_mapping.items()}

    default_supply_selection = [types_full_names[0]] if types_full_names else []
    selected_supply_types_full = st.multiselect(
        f"Select technology types for {selected_sector.lower()} supply timeseries:",
        types_full_names,
        default=default_supply_selection,
        key=f"timeseries_supply_types_{sector_lower}",
    )
    if not selected_supply_types_full and default_supply_selection:
        selected_supply_types_full = default_supply_selection
    selected_supply_types = [
        reverse_mapping.get(v, v) for v in selected_supply_types_full
    ]

    # Filter availability/storage inflow technologies by selected sector where possible.
    if "class" in tech_df.columns:
        class_mapping = dict(zip(tech_df["technology"], tech_df["class"]))
        if selected_sector == "Power":
            allowed_classes = {"Generator", "StorageUnit", "Store", "Link"}
        elif selected_sector == "Industry":
            allowed_classes = {"Link", "StorageUnit", "Store"}
        else:
            allowed_classes = {"Link", "StorageUnit", "Store"}

        filtered_supply_types = [
            t for t in selected_supply_types if class_mapping.get(t) in allowed_classes
        ]
        if filtered_supply_types:
            selected_supply_types = filtered_supply_types

    st.markdown("#### Supply Profiles")
    for title in time_series_config["Global_input"].keys():
        if "demand" not in time_series_config["Global_input"][title]["tag_name"]:
            df_widgets_handler.set_up_df_with_charts(
                sector="Global_input",
                title=title,
                input_df=dfs[title + "_df"],
                selected_types=selected_supply_types,
                selected_countries=selected_countries,
            )

    st.markdown("#### Demand Profiles")
    load_mapping = {
        "HV_LOAD": "Wholesale market load (High voltage level)",
        "LV_LOAD": "Building load (low/medium voltage level)",
    }
    profile_full_names = list(load_mapping.values())
    reverse_load_mapping = {v: k for k, v in load_mapping.items()}
    default_by_sector = {
        "Power": ["HV_LOAD", "LV_LOAD"],
        "Industry": ["LV_LOAD"],
        "Transport": ["LV_LOAD"],
    }
    default_profile_keys = [
        p
        for p in default_by_sector.get(selected_sector, ["HV_LOAD"])
        if p in load_mapping
    ]
    default_profile_selection = [load_mapping[p] for p in default_profile_keys]

    selected_profile_types_full = st.multiselect(
        f"Select demand/load profiles for {selected_sector.lower()}:",
        profile_full_names,
        default=default_profile_selection,
        key=f"timeseries_demand_types_{sector_lower}",
    )
    if not selected_profile_types_full and default_profile_selection:
        selected_profile_types_full = default_profile_selection
    selected_profile_types = [
        reverse_load_mapping.get(v, v) for v in selected_profile_types_full
    ]

    df_widgets_handler.set_up_df_with_charts(
        sector="Global_input",
        title="Demand_Profiles",
        input_df=dfs["Demand_Profiles_df"],
        selected_types=selected_profile_types,
        selected_countries=selected_countries,
    )
