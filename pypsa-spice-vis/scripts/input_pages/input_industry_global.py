# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Global industry input page."""

import streamlit as st

from scripts.data_utils import render_countries_pills, render_type_and_class_filters
from scripts.input_st_handler import (
    generate_global_markdown_message,
    generate_sector_title,
    get_all_countries,
    get_class_type_for_timeseries_data,
    render_demand_profiles_widget,
    render_input_table_section,
    set_available_technology_df,
)

if __name__ == "__main__":
    input_config = st.session_state.input_config
    selected_sector = "Industry"
    sector_title = generate_sector_title(selected_sector)
    all_countries = get_all_countries()
    tech_df = set_available_technology_df(selected_sector, input_config)

    st.subheader(f":globe_with_meridians: Global input  | {sector_title}")
    generate_global_markdown_message()

    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="industry_global_pills",
    )
    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df,
        key="industry_global",
    )

    for title, table_config in input_config["Global_input"].items():
        if not table_config.get("timeseries"):
            render_input_table_section(
                sector="Global_input",
                title=title,
                selected_types=sector_selected_types,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
                input_config=input_config,
            )

    st.subheader(f":material/timer: Timeseries Profiles  | {sector_title}")
    filtered_supply_types = get_class_type_for_timeseries_data(
        tech_df,
        selected_sector,
        sector_selected_types,
    )
    for title, table_config in input_config["Global_input"].items():
        if table_config.get("timeseries") and "demand" not in table_config["tag_name"]:
            render_input_table_section(
                sector="Global_input",
                title=title,
                selected_types=filtered_supply_types,
                selected_countries=sector_selected_countries,
                input_config=input_config,
            )

    st.subheader(f":material/timer: Demand Profiles  | {sector_title}")
    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_industry_global_pills",
    )
    render_demand_profiles_widget(
        selected_sector=selected_sector,
        sector_selected_countries=demand_selected_countries,
        input_config=input_config,
    )
