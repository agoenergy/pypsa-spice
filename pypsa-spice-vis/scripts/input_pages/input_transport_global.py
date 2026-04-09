# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Global transport input page.

Page displays and allows editing of global input tables for the transport sector.
"""


import os

import pandas as pd
import streamlit as st

from scripts.data_utils import render_countries_pills, render_type_and_class_filters
from scripts.input_st_handler import (
    create_editable_df,
    generate_global_markdown_message,
    generate_sector_title,
    get_all_countries,
    get_empty_df_notice_message,
    get_table_config_and_path,
    get_unique_type_key,
    render_demand_profiles_selectbox,
    render_line_chart,
    render_save_button_for_input_df,
    set_available_technology_df,
    set_general_filter_df,
)


def render_input_transport_global_section(
    title_string: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
) -> None:
    """Render global transport input section with editable tables and charts.

    Parameters
    ----------
    title_string : str
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
    """
    # Get table configuration and input CSV path based on the title_string and sector
    table_config, input_csv_path = get_table_config_and_path(
        title_string=title_string,
        sector="Global_input",
        input_config=input_config,
    )
    csv_identifier = table_config["identifier"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Global_input",
        title_string,
        selected_types,
        selected_classes,
        selected_countries,
    )

    edited_df_key = f"{title_string}_{csv_identifier}_editor_{unique_type_key}"
    has_changes_key = f"has_changes_{title_string}_{csv_identifier}_{unique_type_key}"
    save_button_key = f"save_{title_string}_{csv_identifier}_{unique_type_key}"

    # Render the section in an expander with the title_string and CSV path
    with st.expander(title_string):
        st.write(f"### {title_string}")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )

        if input_df is None:
            if not os.path.exists(input_csv_path):
                st.error(f"File not found: {input_csv_path}")
                return
            input_df = pd.read_csv(input_csv_path)

        if title_string == "EV_parameters":
            # For ev parameters table, all rows will be shown
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


def render_input_transport_timeseries_section(
    title_string: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
) -> None:
    """Render timeseries transport input section with charts.

    Parameters
    ----------
    title_string : str
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
    """
    # Get table configuration and input CSV path based on the title_string and sector
    table_config, input_csv_path = get_table_config_and_path(
        title_string=title_string,
        sector="Global_input",
        input_config=input_config,
    )

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Global_input",
        title_string,
        selected_types,
        selected_classes,
        selected_countries,
    )

    # Render the section in an expander with the title_string and CSV path
    with st.expander(title_string):
        st.write(f"### {title_string}")
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

        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        # Render line charts for the timeseries data
        render_line_chart(filtered_df, table_config, unique_type_key)


def render_input_transport_demand_profile_section(
    title_string: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
) -> None:
    """Render global transport demand input section with editable charts.

    Parameters
    ----------
    title_string : str
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
    """
    # Get table configuration and input CSV path based on the title_string and sector
    table_config, input_csv_path = get_table_config_and_path(
        title_string=title_string,
        sector="Global_input",
        input_config=input_config,
    )

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Global_input",
        title_string,
        selected_types,
        selected_classes,
        selected_countries,
    )

    # Render the section in an expander with the title_string and CSV path
    with st.expander(title_string):
        st.write(f"### {title_string}")
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

        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        # Render line charts for the timeseries data
        render_line_chart(filtered_df, table_config, unique_type_key)


if __name__ == "__main__":
    input_config = st.session_state.input_config
    selected_sector = "Transport"
    sector_title = generate_sector_title(selected_sector)
    DOCS_PATH = "getting-started/input-data/global_csv_template"
    st.markdown(
        "Detailed explanation can be found in: "
        f"[global input guides](https://agoenergy.github.io/pypsa-spice/{DOCS_PATH})"
    )

    all_countries = get_all_countries()
    tech_df = set_available_technology_df(selected_sector, input_config)

    st.subheader(f":globe_with_meridians: Global input  | {sector_title}")
    generate_global_markdown_message()

    # Render country selection pills for the transport global input section
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="transport_global_pills",
    )

    # Render type and PyPSA class filters for the transport global input section
    sector_selected_types = ["EVCH", "EVST"]
    sector_selected_classes = ["Link", "Store"]
    st.markdown(f"Tech: **{', '.join(sector_selected_types)}**")
    st.markdown(f"Class: **{', '.join(sector_selected_classes)}**")

    # Render global transport input relevant sections for non-timeseries tables
    for title_string, table_config in input_config["Global_input"].items():
        if (
            not table_config.get("timeseries")
            and "storage_costs" not in title_string.lower()
            and "inflows" not in title_string.lower()
        ):
            render_input_transport_global_section(
                title_string=title_string,
                selected_types=sector_selected_types,
                input_config=input_config,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
            )

    st.subheader(f":material/timer: Demand Profiles  | {sector_title}")

    # Render country selection pills for the transport demand profiles section
    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_transport_global_pills",
    )

    # Render type and PyPSA class filters for the transport demand profiles section
    demand_profile_types, selected_classes = render_type_and_class_filters(
        tech_df,
        key="transport_demand_global",
    )

    # Render demand profiles selectbox for the transport demand profiles section
    demand_profile_types = render_demand_profiles_selectbox(selected_sector="Transport")

    # Render global transport demand profiles section
    render_input_transport_demand_profile_section(
        title_string="Demand_Profiles",
        selected_types=demand_profile_types,
        input_config=input_config,
        selected_countries=demand_selected_countries,
        selected_classes=selected_classes,
    )
