# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Scenario-specific industry input page."""

import streamlit as st

from scripts.data_utils import render_countries_pills, render_type_and_class_filters
from scripts.input_st_handler import (
    generate_sector_title,
    get_all_countries,
    render_input_table_section,
    set_available_technology_df,
)


def render_page() -> None:
    """Render the scenario-specific industry input page."""
    input_config = st.session_state.input_config
    selected_sector = "Industry"
    sector_title = generate_sector_title(selected_sector)
    selected_scenario = st.session_state.input_sce1
    all_countries = get_all_countries()
    tech_df = set_available_technology_df(selected_sector, input_config)

    st.subheader(f":material/timeline: Scenario specific input | {sector_title}")
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="industry_scenario_pills",
    )
    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df,
        key="industry_scenario",
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


if __name__ == "__main__":
    render_page()
