# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create scenario-specific power input page.

Page displays and allows editing of scenario-specific input tables for the power
sector.
"""

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


def render_input_power_scenario_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
    selected_scenario: str | None = None,
) -> None:
    """Render scenario specific power input section with editable tables and charts.

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
        sector="Power",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Power",
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
                filter_col="carrier",
                selected_types=list(fuels.values()),
            )
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


def render_input_load_section(
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
        sector="Power",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Power",
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


def render_interconnections_section(
    grid_selected_countries: list[str],
    grid_df: pd.DataFrame,
    input_config: dict,
    selected_scenario: str,
) -> None:
    """Render interconnections section with editable tables.

    Parameters
    ----------
    grid_selected_countries : list[str]
        list of selected countries for filtering the interconnections table
    grid_df : pd.DataFrame
        dataframe containing the grid data
    input_config : dict
        configuration dictionary for input_config
    selected_scenario : str
        selected scenario for filtering the interconnections table
    """
    # Render interconnection type selectbox between Transmission and Distribution
    st.selectbox(
        "Select type of interconnection:",
        options=["Transmission", "Distribution"],
        index=0,
        key="power_static_intercon_type",
        help="Choose the grid voltage level to inspect and edit interconnectors.",
    )
    st.markdown("Class: **Link**")
    grid_type = st.session_state.get("power_static_intercon_type")

    # Filter the grid dataframe based on the selected interconnection type
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

    # Get table configuration and input CSV path based on the title and sector
    table_config, input_csv_path = get_table_config_and_path(
        title="Interconnectors",
        sector="Grids",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]
    selected_types = ["ITCN"]
    selected_classes = ["Link"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Power",
        "Interconnectors",
        selected_types,
        selected_classes,
        grid_selected_countries,
        selected_scenario,
    )

    edited_df_key = f"Interconnectors_{csv_identifier}_editor_{unique_type_key}"
    has_changes_key = f"has_changes_Interconnectors_{csv_identifier}_{unique_type_key}"
    save_button_key = f"save_Interconnectors_{csv_identifier}_{unique_type_key}"

    # Render the section in an expander with the title and CSV path
    with st.expander("Interconnectors"):
        st.write("### Interconnectors")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )

        if intercon_df is None:
            if not os.path.exists(input_csv_path):
                st.error(f"File not found: {input_csv_path}")
                return
            intercon_df = pd.read_csv(input_csv_path)

        # Apply general filtering based on selected types
        filtered_df = set_general_filter_df(
            df=intercon_df,
            filter_col=table_config["filter_col"],
            selected_types=selected_types,
        )

        # Further filter df based on selected countries if the 'country' column exists
        if grid_selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[
                filtered_df["country"].isin(grid_selected_countries)
            ]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        # Render the editable table
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
    input_config = st.session_state.input_config
    selected_sector = "Power"
    sector_title = generate_sector_title(selected_sector)
    DOCS_PATH = "getting-started/input-data/regional_csv_template"
    st.markdown(
        "Detailed explanation can be found in: "
        f"[scenario input guides](https://agoenergy.github.io/pypsa-spice/{DOCS_PATH})"
    )

    selected_scenario = st.session_state.input_sce1
    all_countries = get_all_countries()
    tech_df = set_available_technology_df(selected_sector, input_config)

    st.subheader(f":material/timeline: Scenario specific input  | {sector_title}")

    # Render country selection pills for the power scenario input section
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="power_scenario_pills",
    )

    # Render type and PyPSA class filters for the power scenario input section
    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df,
        key="power_scenario",
    )

    # Render scenario power input relevant sections for non-timeseries tables
    for title in input_config[selected_sector]:
        if title != "Power_loads":
            render_input_power_scenario_section(
                title=title,
                selected_types=sector_selected_types,
                input_config=input_config,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
                selected_scenario=selected_scenario,
            )

    st.subheader(f":material/timeline: Loads  | {sector_title}")

    # Render country selection pills for the power load input section
    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_power_scenario_pills",
    )

    # Render demand profiles selectbox for the power demand profiles section
    selected_types = render_demand_profiles_selectbox(selected_sector="Power")

    # Render scenario power demand profiles section
    render_input_load_section(
        title="Power_loads",
        selected_types=selected_types,
        input_config=input_config,
        selected_classes=["Load"],
        selected_countries=demand_selected_countries,
        selected_scenario=selected_scenario,
    )

    st.subheader(":material/diagonal_line: Interconnections")

    # Render country selection pills for the interconnections section
    grid_selected_countries = (
        render_countries_pills(
            all_countries=all_countries,
            key="grid_scenario_pills",
        )
        or []
    )

    # Load grid data for the selected scenario
    grid_path = os.path.join(
        st.session_state.input_path,
        selected_scenario,
        "power",
        input_config["Grids"]["Interconnectors"]["csv_name"],
    )
    grid_df = pd.read_csv(grid_path)

    # Render interconnections section
    render_interconnections_section(
        grid_selected_countries,
        grid_df,
        input_config,
        selected_scenario,
    )
