# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Static page under Input section."""

import os

import pandas as pd
import streamlit as st

from scripts.data_utils import (
    get_input_scenario_list,
    render_countries_n_scenario_pills,
    render_type_and_class_filters,
)
from scripts.input_st_handler import DataFrameWidgetsHandler


def render_demand_profiles_widget(
    sector_selected_countries: list,
    df_widgets_handler: DataFrameWidgetsHandler,
    demand_profile_df: pd.DataFrame,
    specify_sector: str,
):
    """Render demand profile section."""
    load_mapping = {
        "HV_LOAD": "Transmission/Wholesale market load (High voltage level)",
        "LV_LOAD": "Distribution/Building load (low/medium voltage level)",
        "IND_LOAD": "Industrial load (both high- and low-temperature heat)",
        "HPV_LOAD": "Transport load (high voltage level)",
        "LPV_LOAD": "Transport load (low/medium voltage level)",
    }
    reverse_load_mapping = {v: k for k, v in load_mapping.items()}
    default_by_sector = {
        "Power": ["HV_LOAD", "LV_LOAD"],
        "Industry": ["IND_LOAD"],
        "Transport": ["HPV_LOAD", "LPV_LOAD"],
    }
    default_profile_keys = [
        p for p in default_by_sector.get(specify_sector) if p in load_mapping
    ]
    default_profile_selection = [load_mapping[p] for p in default_profile_keys]

    st.selectbox(
        "Select type of demand/load profiles:",
        options=default_profile_selection,
        index=0,
        key="timeseries_demand_type",
    )
    st.markdown("Class: **Load**")
    selected_profile_label = st.session_state.get("timeseries_demand_type")
    selected_profile_key = reverse_load_mapping.get(selected_profile_label)

    df_widgets_handler.set_up_df_with_charts(
        sector="Global_input",
        title="Demand_Profiles",
        input_df=demand_profile_df,
        selected_types=[selected_profile_key],
        selected_countries=sector_selected_countries,
    )


if __name__ == "__main__":
    base_config = st.session_state.base_config
    input_config = st.session_state.input_config
    selected_sector = st.session_state.get("selected_input_sector", "Power")

    df_widgets_handler = DataFrameWidgetsHandler(input_config=input_config)

    all_countries = list(base_config["base_configs"]["regions"].keys())

    if selected_sector == "Power":
        SECTOR_TITLE = ":material/bolt: Power"
    elif selected_sector == "Industry":
        SECTOR_TITLE = ":material/construction: Industry"
    else:
        SECTOR_TITLE = ":material/directions_car: Transport"

    st.subheader(":globe_with_meridians: Global input  | " + SECTOR_TITLE)

    sector_lower = selected_sector.lower()
    sector_selected_countries, sector_selected_scenario = (
        render_countries_n_scenario_pills(
            scenario_options=get_input_scenario_list(base_config),
            all_countries=all_countries,
            key=f"{sector_lower}_global_pills",
        )
    )

    # Reload dfs
    dfs = df_widgets_handler.load_all_dfs(
        selected_scenario=sector_selected_scenario, specify_sector="Global_input"
    )

    tech_path = os.path.join(
        st.session_state.input_path,
        "global_input",
        input_config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)

    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df, key=f"{sector_lower}_global"
    )

    for title in input_config["Global_input"].keys():
        if not input_config["Global_input"][title].get("timeseries"):
            df_widgets_handler.set_up_df_with_charts(
                title=title,
                input_df=dfs[title + "_df"],
                selected_types=sector_selected_types,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
            )

    st.subheader(":material/timer: Timeseries Profiles  | " + SECTOR_TITLE)

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
            t for t in sector_selected_types if class_mapping.get(t) in allowed_classes
        ]
        if filtered_supply_types:
            sector_selected_types = filtered_supply_types

    for title in input_config["Global_input"].keys():
        if (
            input_config["Global_input"][title].get("timeseries")
            and "demand" not in input_config["Global_input"][title]["tag_name"]
        ):
            df_widgets_handler.set_up_df_with_charts(
                sector="Global_input",
                title=title,
                input_df=dfs[title + "_df"],
                selected_types=sector_selected_types,
                selected_countries=sector_selected_countries,
            )

    st.subheader(":material/timer: Demand Profiles  | " + SECTOR_TITLE)

    demand_selected_countries, demand_selected_scenario = (
        render_countries_n_scenario_pills(
            scenario_options=get_input_scenario_list(base_config),
            all_countries=all_countries,
            key="demand_global_pills",
        )
    )

    demand_profile_path = os.path.join(
        df_widgets_handler.base_input_path,
        "global_input",
        input_config["Global_input"]["Demand_Profiles"]["csv_name"],
    )
    demand_profile_df = pd.read_csv(demand_profile_path)

    render_demand_profiles_widget(
        sector_selected_countries=sector_selected_countries,
        df_widgets_handler=df_widgets_handler,
        demand_profile_df=demand_profile_df,
        specify_sector=selected_sector,
    )
