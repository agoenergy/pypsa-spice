# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Scenario-specific power input page."""

import os

import pandas as pd
import streamlit as st

from scripts.data_utils import render_countries_pills, render_type_and_class_filters
from scripts.input_st_handler import (
    generate_sector_title,
    get_all_countries,
    render_input_table_section,
    set_available_technology_df,
)


def render_interconnections_widget(
    grid_selected_countries: list[str],
    grid_df: pd.DataFrame,
    input_config: dict,
    selected_scenario: str,
) -> None:
    """Render the interconnection editor for the power sector."""
    st.selectbox(
        "Select type of interconnection:",
        options=["Transmission", "Distribution"],
        index=0,
        key="power_static_intercon_type",
        help="Choose the grid voltage level to inspect and edit interconnectors.",
    )
    st.markdown("Class: **Link**")
    grid_type = st.session_state.get("power_static_intercon_type")

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

    render_input_table_section(
        sector="Grids",
        title="Interconnectors",
        selected_types=["ITCN"],
        input_df=intercon_df,
        selected_classes=["Link"],
        selected_countries=grid_selected_countries,
        selected_scenario=selected_scenario,
        input_config=input_config,
    )


def render_page() -> None:
    """Render the scenario-specific power input page."""
    input_config = st.session_state.input_config
    selected_sector = "Power"
    sector_title = generate_sector_title(selected_sector)
    selected_scenario = st.session_state.input_sce1
    all_countries = get_all_countries()
    tech_df = set_available_technology_df(selected_sector, input_config)

    st.subheader(f":material/timeline: Scenario specific input | {sector_title}")
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="power_scenario_pills",
    )
    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df,
        key="power_scenario",
    )

    for title in input_config[selected_sector]:
        render_input_table_section(
            title=title,
            sector=selected_sector,
            selected_types=sector_selected_types,
            selected_classes=sector_selected_classes,
            selected_countries=sector_selected_countries,
            selected_scenario=selected_scenario,
            input_config=input_config,
        )

    st.subheader(":material/diagonal_line: Interconnections")
    grid_selected_countries = (
        render_countries_pills(
            all_countries=all_countries,
            key="grid_scenario_pills",
        )
        or []
    )
    grid_path = os.path.join(
        st.session_state.input_path,
        selected_scenario,
        "power",
        input_config["Grids"]["Interconnectors"]["csv_name"],
    )
    grid_df = pd.read_csv(grid_path)
    render_interconnections_widget(
        grid_selected_countries,
        grid_df,
        input_config,
        selected_scenario,
    )


if __name__ == "__main__":
    render_page()
