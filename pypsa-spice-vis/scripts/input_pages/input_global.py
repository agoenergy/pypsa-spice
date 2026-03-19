# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Static page under Input section."""

import os

import pandas as pd
import streamlit as st

from scripts.data_utils import (
    load_tech_info_mapping_df,
    render_countries_pills,
    render_type_and_class_filters,
)
from scripts.input_st_handler import DataFrameWidgetsHandler


def render_demand_profiles_widget(
    sector_selected_countries: list | None,
    df_widgets_handler: DataFrameWidgetsHandler,
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
        p for p in default_by_sector.get(specify_sector, []) if p in load_mapping
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

    if selected_profile_key is None:
        st.info("Select a demand profile type to load data.")
        return

    df_widgets_handler.set_up_df_with_charts(
        sector="Global_input",
        title="Demand_Profiles",
        selected_types=[selected_profile_key],
        selected_countries=sector_selected_countries,
    )


if __name__ == "__main__":
    base_config = st.session_state.base_config
    input_config = st.session_state.input_config
    selected_sector = st.session_state.get("selected_input_sector", "Power")

    df_widgets_handler = DataFrameWidgetsHandler(input_config=input_config)

    all_countries = list(base_config["base_configs"]["regions"].keys())

    # df loaded from setting/tech_mapping.csv
    tech_info_df = load_tech_info_mapping_df()

    if selected_sector == "Power":
        SECTOR_TITLE = ":material/bolt: Power"
        tech_info_df = tech_info_df[tech_info_df["sector"] == "Power"]
    elif selected_sector == "Industry":
        SECTOR_TITLE = ":material/construction: Industry"
        tech_info_df = tech_info_df[tech_info_df["sector"] == "Industry"]
    else:
        SECTOR_TITLE = ":material/directions_car: Transport"
        tech_info_df = tech_info_df[tech_info_df["sector"] == "Transport"]

    # =============================================================================
    # Global input
    # =============================================================================

    st.subheader(":globe_with_meridians: Global input  | " + SECTOR_TITLE)

    st.markdown(
        """
        <p style="color:orange; font-weight:bold;">
        ⚠️ Changes made to the global input files will be automatically applied
        across all scenarios.
        </p>
        """,
        unsafe_allow_html=True,
    )

    sector_lower = selected_sector.lower()
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key=f"{sector_lower}_global_pills",
    )

    # df from global input technology csv
    tech_path = os.path.join(
        st.session_state.input_path,
        "global_input",
        input_config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)

    # Keep only technologies defined in the technical mapping/info table.
    allowed_technologies = set()
    if isinstance(tech_info_df, pd.DataFrame):
        if not tech_info_df.index.empty:
            allowed_technologies.update(
                tech_info_df.index.astype(str).str.strip().tolist()
            )

        for col in ["original_names", "nice_names"]:
            if col in tech_info_df.columns:
                allowed_technologies.update(
                    tech_info_df[col].dropna().astype(str).str.strip().tolist()
                )

    if "technology" in tech_df.columns and allowed_technologies:
        tech_df["technology"] = tech_df["technology"].astype(str).str.strip()
        tech_df = tech_df[tech_df["technology"].isin(allowed_technologies)]

    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df, key=f"{sector_lower}_global"
    )

    for title in input_config["Global_input"].keys():
        if not input_config["Global_input"][title].get("timeseries"):
            df_widgets_handler.set_up_df_with_charts(
                sector="Global_input",
                title=title,
                selected_types=sector_selected_types,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
            )

    # =============================================================================
    # Timeseries profiles
    # =============================================================================

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
                selected_types=sector_selected_types,
                selected_countries=sector_selected_countries,
            )

    st.subheader(":material/timer: Demand Profiles  | " + SECTOR_TITLE)

    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_global_pills",
    )

    render_demand_profiles_widget(
        sector_selected_countries=demand_selected_countries,
        df_widgets_handler=df_widgets_handler,
        specify_sector=selected_sector,
    )
