# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Scenario specific page under Input section."""

import os

import pandas as pd
import streamlit as st

from scripts.data_utils import (
    load_tech_info_mapping_df,
    render_countries_pills,
    render_type_and_class_filters,
)
from scripts.input_st_handler import DataFrameWidgetsHandler


def get_mapping_list(*dfs: pd.DataFrame) -> list[str]:
    """Get sorted technology/profile types from one or more dataframes."""
    type_set = set()

    for df in dfs:
        if "technology" in df.columns:
            type_set |= set(df["technology"].unique())
        if "profile_type" in df.columns:
            type_set |= set(df["profile_type"].unique())

    return sorted(type_set)


def render_interconnections_widget(
    grid_selected_countries: list,
    df_widgets_handler: DataFrameWidgetsHandler,
    grid_df: pd.DataFrame,
):
    """Render interconnections section."""
    st.selectbox(
        "Select type of interconnection:",
        options=["Transmission", "Distribution"],
        index=0,
        key="static_intercon_type",
    )
    st.markdown("Class: **Link**")
    grid_type = st.session_state.get("static_intercon_type", None)

    if grid_type == "Transmission":
        intercon_df = grid_df[
            (grid_df["bus0"].str.contains("HVELEC"))
            & (grid_df["bus1"].str.contains("HVELEC"))
        ]
    else:
        intercon_df = grid_df[
            (grid_df["bus0"].str.contains("LVELEC"))
            & (grid_df["bus1"].str.contains("LVELEC"))
        ]

    df_widgets_handler.set_up_df_with_charts(
        sector="Grids",
        title="Interconnectors",
        selected_types=["ITCN"],
        input_df=intercon_df,
        selected_classes=["Link"],
        selected_countries=grid_selected_countries,
    )


if __name__ == "__main__":
    base_config = st.session_state.base_config
    input_config = st.session_state.input_config
    sector_selected_scenario = st.session_state.input_sce1
    selected_sector = st.session_state.get("selected_input_sector", "Power")
    sector_lower = selected_sector.lower()

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

    st.subheader(":material/timeline: Scenario specific input | " + SECTOR_TITLE)

    sector_lower = selected_sector.lower()
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key=f"{sector_lower}_global_pills",
    )

    # df from global input technology csv
    tech_path = os.path.join(
        df_widgets_handler.base_input_path,
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
        tech_df, key=f"{sector_lower}_scenario"
    )

    for title in input_config[selected_sector].keys():
        df_widgets_handler.set_up_df_with_charts(
            title=title,
            sector=selected_sector,
            selected_types=sector_selected_types,
            selected_classes=sector_selected_classes,
            selected_countries=sector_selected_countries,
            selected_scenario=sector_selected_scenario,
        )

    if selected_sector == "Power":
        st.subheader(":material/diagonal_line: Interconnections")
        grid_selected_countries = (
            render_countries_pills(
                all_countries=all_countries,
                key="grid_scenario_pills",
            )
            or []
        )
        grid_path = os.path.join(
            df_widgets_handler.base_input_path,
            sector_selected_scenario,
            "power",
            input_config["Grids"]["Interconnectors"]["csv_name"],
        )
        grid_df = pd.read_csv(grid_path)

        render_interconnections_widget(
            sector_selected_countries, df_widgets_handler, grid_df
        )
