# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Global industry input page.

Page displays and allows editing of global input tables for the industry sector.
"""

# pylint: disable=too-many-arguments,too-many-locals, too-many-positional-arguments
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


def render_input_industry_global_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
) -> None:
    """Render global industry input section with editable tables and charts.

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
    """
    # Get table configuration and input CSV path based on the title and sector
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Global_input",
        input_config=input_config,
    )
    csv_identifier = table_config["identifier"]

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Global_input",
        title,
        selected_types,
        selected_classes,
        selected_countries,
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

        if "direct_air_capture" in title.lower():
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


def render_input_industry_timeseries_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
) -> None:
    """Render timeseries industry input section with charts.

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
    """
    # Get table configuration and input CSV path based on the title and sector
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Global_input",
        input_config=input_config,
    )

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Global_input",
        title,
        selected_types,
        selected_classes,
        selected_countries,
    )

    # Render the section in an expander with the title and CSV path
    with st.expander(title):
        st.write(f"### {title}")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )
        st.markdown(
            "⚠️ To modify the time series data, please make changes locally. "
            "This section in this app contains visualisation only.",
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


def render_input_industry_demand_profile_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
) -> None:
    """Render global industry demand input section with editable charts.

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
    """
    # Get table configuration and input CSV path based on the title and sector
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Global_input",
        input_config=input_config,
    )

    # Generate a unique key for Streamlit session state
    unique_type_key = get_unique_type_key(
        "Global_input",
        title,
        selected_types,
        selected_classes,
        selected_countries,
    )

    # Render the section in an expander with the title and CSV path
    with st.expander(title):
        st.write(f"### {title}")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )
        st.markdown(
            "⚠️ To modify the time series data, please make changes locally. "
            "This section in this app contains visualisation only.",
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
    app_input_config = st.session_state.input_config
    SELECTED_SECTOR = "Industry"
    sector_title = generate_sector_title(SELECTED_SECTOR)
    DOCS_PATH = "getting-started/input-data/global_csv_template"
    st.markdown(
        "Detailed explanation can be found in: "
        f"[global input guides](https://agoenergy.github.io/pypsa-spice/{DOCS_PATH})"
    )

    all_countries = get_all_countries()
    tech_df = set_available_technology_df(SELECTED_SECTOR, app_input_config)

    st.subheader(f":globe_with_meridians: Global input  | {sector_title}")
    generate_global_markdown_message()

    # Render country selection pills for the industry global input section
    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="industry_global_pills",
    )

    # Render type and PyPSA class filters for the industry global input section
    with st.sidebar:
        st.divider()
        st.markdown("#### Technology filter | :material/construction: Industry")
        sector_selected_types, sector_selected_classes = render_type_and_class_filters(
            tech_df,
            key="industry_global",
        )

    # Render global industry input relevant sections for non-timeseries tables
    for table_name, section_config in app_input_config["Global_input"].items():
        if (
            not section_config.get("timeseries")
            and "inflows" not in table_name.lower()
            and "ev_parameters" not in table_name.lower()
        ):
            render_input_industry_global_section(
                title=table_name,
                selected_types=sector_selected_types,
                input_config=app_input_config,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
            )

    st.subheader(f":material/timer: Timeseries Profiles  | {sector_title}")

    # Render global industry input relevant sections for timeseries tables
    for table_name, section_config in app_input_config["Global_input"].items():
        if section_config.get("timeseries") and "demand" not in table_name.lower():

            render_input_industry_timeseries_section(
                title=table_name,
                selected_types=sector_selected_types,
                input_config=app_input_config,
                selected_classes=sector_selected_classes,
                selected_countries=sector_selected_countries,
            )

    st.subheader(f":material/timer: Demand Profiles  | {sector_title}")

    # Render country selection pills for the industry demand profiles section
    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_industry_global_pills",
    )

    # Render demand profiles selectbox for the industry demand profiles section
    demand_profile_types = render_demand_profiles_selectbox(
        selected_sector=SELECTED_SECTOR
    )

    # Render global industry demand profiles section
    render_input_industry_demand_profile_section(
        title="Demand_Profiles",
        selected_types=demand_profile_types,
        input_config=app_input_config,
        selected_countries=demand_selected_countries,
        selected_classes=["Loads"],
    )
