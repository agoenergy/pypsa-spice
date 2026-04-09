# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create scenario-specific industry input page.

Page displays and allows editing of scenario-specific input tables for the industry
sector.
"""

# pylint: disable=too-many-arguments,too-many-locals, too-many-positional-arguments
import os

import pandas as pd
import streamlit as st

from scripts.data_utils import render_countries_pills, render_type_and_class_filters
from scripts.input_st_handler import (
    create_editable_df,
    generate_sector_title,
    get_all_countries,
    get_empty_df_notice_message,
    get_fuel_mapping,
    get_table_config_and_path,
    get_unique_type_key,
    render_demand_profiles_selectbox,
    render_line_chart,
    render_save_button_for_input_df,
    set_available_technology_df,
    set_decommission_filter_df,
    set_general_filter_df,
)


def render_input_industry_scenario_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
    selected_scenario: str | None = None,
) -> None:
    """Render scenario specific industry input section with editable tables and charts.

    Parameters
    ----------
    title : str
        key or heading of the indicator from input_settings.yaml
    selected_types : list[str]
        list of selected technology types for filtering the input table
    input_config : dict
        configuration dictionary for input_config
    input_df : pd.DataFrame | None, optional
        default input dataframe, by default None
    selected_classes : list[str] | None, optional
        list of selected PyPSA classes for filtering the input table, by default None
    selected_countries : list[str] | None, optional
        list of selected countries for filtering the input table, by default None
    selected_scenario : str | None, optional
        selected scenario for filtering the input table, by default None
    """
    # Get table configuration and input CSV path based on the title and sector
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Industry",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Industry",
        title,
        selected_types,
        selected_classes,
        selected_countries,
        selected_scenario,
    )

    edited_df_key = f"{title}_{csv_identifier}_editor_{unique_type_key}"
    has_changes_key = f"has_changes_{title}_{csv_identifier}_{unique_type_key}"
    save_button_key = f"save_{title}_{csv_identifier}_{unique_type_key}"

    # Render the section in an expander with the title and CSV path
    with st.expander(title):
        st.write(f"### {title}")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )

        if input_df is None:
            if not os.path.exists(input_csv_path):
                st.error(f"File not found: {input_csv_path}")
                return
            input_df = pd.read_csv(input_csv_path)

        if "decommission" in title.lower():
            # For decommissioning, filter column is different defined in table_config
            filtered_df = set_decommission_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=selected_types,
            )
        elif table_config["filter_col"] == "carrier":
            # Tables filtered by carriers
            fuels = (
                get_fuel_mapping(selected_types, input_config)
                if selected_classes and "Link" in selected_classes
                else {}
            )
            filtered_df = set_general_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=list(fuels.values()),
            )
        elif "direct" in title.lower():
            # For direct air capture table, all rows will be shown
            filtered_df = input_df
        else:
            # For other tables, apply general filtering based on selected types
            filtered_df = set_general_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=selected_types,
            )

        # Further filter df based on selected countries if the 'country' column exists
        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        if table_config["with_charts"]:
            # If the table is with line charts, render them in a separate tab
            table_tab, visualisation_tab = st.tabs(["Table", "Visualisation"])
            with table_tab:
                edited_df, to_save = create_editable_df(
                    filtered_df,
                    edited_df_key,
                    has_changes_key,
                )
            with visualisation_tab:
                render_line_chart(filtered_df, table_config, unique_type_key)
        else:
            # Otherwise, just show the editable table
            edited_df, to_save = create_editable_df(
                filtered_df,
                edited_df_key,
                has_changes_key,
            )

        # Check if there are changes to save and render the save button
        has_changes = st.session_state.get(has_changes_key, False)
        if to_save:
            render_save_button_for_input_df(
                filtered_df,
                edited_df,
                has_changes,
                has_changes_key,
                save_button_key,
                input_csv_path,
                message_delay=1,
            )


def render_input_industry_annual_demand_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
    selected_scenario: str | None = None,
) -> None:
    """Render scenario specific load section with editable tables and charts.

    Parameters
    ----------
    title : str
        key or heading of the indicator from input_settings.yaml
    selected_types : list[str]
        list of selected technology types for filtering the input table
    input_config : dict
        configuration dictionary for input_config
    input_df : pd.DataFrame | None, optional
        default input dataframe, by default None
    selected_classes : list[str] | None, optional
        list of selected PyPSA classes for filtering the input table, by default None
    selected_countries : list[str] | None, optional
        list of selected countries for filtering the input table, by default None
    selected_scenario : str | None, optional
        selected scenario for filtering the input table, by default None
    """
    # Get table configuration and input CSV path based on the title and sector
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Industry",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Industry",
        title,
        selected_types,
        selected_classes,
        selected_countries,
        selected_scenario,
    )

    edited_df_key = f"{title}_{csv_identifier}_editor_{unique_type_key}"
    has_changes_key = f"has_changes_{title}_{csv_identifier}_{unique_type_key}"
    save_button_key = f"save_{title}_{csv_identifier}_{unique_type_key}"

    # Render the section in an expander with the title and CSV path
    with st.expander(title):
        st.write(f"### {title}")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )

        if input_df is None:
            if not os.path.exists(input_csv_path):
                st.error(f"File not found: {input_csv_path}")
                return
            input_df = pd.read_csv(input_csv_path)

        # Apply general filtering based on selected types
        filtered_df = set_general_filter_df(
            df=input_df,
            filter_col=table_config["filter_col"],
            selected_types=selected_types,
        )

        # Further filter df based on selected countries if the 'country' column exists
        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        if table_config["with_charts"]:
            # If the table is with line charts, render them in a separate tab
            table_tab, visualisation_tab = st.tabs(["Table", "Visualisation"])
            with table_tab:
                edited_df, to_save = create_editable_df(
                    filtered_df,
                    edited_df_key,
                    has_changes_key,
                )
            with visualisation_tab:
                render_line_chart(filtered_df, table_config, unique_type_key)
        else:
            # Otherwise, just show the editable table
            edited_df, to_save = create_editable_df(
                filtered_df,
                edited_df_key,
                has_changes_key,
            )

        # Check if there are changes to save and render the save button
        has_changes = st.session_state.get(has_changes_key, False)
        if to_save:
            render_save_button_for_input_df(
                filtered_df,
                edited_df,
                has_changes,
                has_changes_key,
                save_button_key,
                input_csv_path,
                message_delay=1,
            )


if __name__ == "__main__":
    app_input_config = st.session_state.input_config
    SELECTED_SECTOR = "Industry"
    sector_title = generate_sector_title(SELECTED_SECTOR)
    DOCS_PATH = "getting-started/input-data/regional_csv_template"
    st.markdown(
        "Detailed explanation can be found in: "
        f"[scenario input guides](https://agoenergy.github.io/pypsa-spice/{DOCS_PATH})"
    )

    scenario_name = st.session_state.input_sce1
    all_countries = get_all_countries()
    tech_df = set_available_technology_df(SELECTED_SECTOR, app_input_config)

    st.subheader(f":material/timeline: Scenario specific input  | {sector_title}")

    # Render country selection pills for the industry scenario input section
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="industry_scenario_pills",
    )

    # Render type and PyPSA class filters for the industry scenario input section
    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df,
        key="industry_scenario",
    )

    # Render scenario industry input relevant sections for non-timeseries tables
    for table_name in app_input_config[SELECTED_SECTOR]:
        if table_name != "Heat_loads":
            render_input_industry_scenario_section(
                title=table_name,
                selected_types=sector_selected_types,
                input_config=app_input_config,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
                selected_scenario=scenario_name,
            )

    st.subheader(f":material/timeline: Loads  | {sector_title}")

    # Render country selection pills for the industry load input section
    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_industry_scenario_pills",
    )

    # Render demand profiles selectbox for the industry demand profiles section
    demand_profile_types = render_demand_profiles_selectbox(selected_sector="Industry")

    # Render scenario industry demand profiles section
    render_input_industry_annual_demand_section(
        title="Heat_loads",
        selected_types=demand_profile_types,
        input_config=app_input_config,
        selected_classes=["Load"],
        selected_countries=demand_selected_countries,
        selected_scenario=scenario_name,
    )
