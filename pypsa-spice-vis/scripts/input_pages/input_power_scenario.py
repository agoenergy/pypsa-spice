# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Scenario-specific power input page."""

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
    render_save_button,
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
    selected_scenario: list[str] | None = None,
) -> None:
    """Render one input table section with editing and charts."""
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Power",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]
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
            filtered_df = set_decommission_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=selected_types,
            )
        elif table_config["filter_col"] == "carrier":
            fuels = (
                get_fuel_mapping(selected_types, input_config)
                if "Link" in selected_classes
                else {}
            )
            filtered_df = set_general_filter_df(
                df=input_df,
                filter_col="carrier",
                selected_types=list(fuels.values()),
            )
        else:
            filtered_df = set_general_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=selected_types,
            )
        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        if table_config["with_charts"]:
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
            edited_df, to_save = create_editable_df(
                filtered_df,
                edited_df_key,
                has_changes_key,
            )

        has_changes = st.session_state.get(has_changes_key, False)
        if to_save:
            render_save_button(
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
    selected_scenario: list[str] | None = None,
) -> None:
    """Render one input table section with editing and charts."""
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector="Power",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]
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

        filtered_df = set_general_filter_df(
            df=input_df,
            filter_col=table_config["filter_col"],
            selected_types=selected_types,
        )

        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        if table_config["with_charts"]:
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
            edited_df, to_save = create_editable_df(
                filtered_df,
                edited_df_key,
                has_changes_key,
            )

        has_changes = st.session_state.get(has_changes_key, False)
        if to_save:
            render_save_button(
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

    table_config, input_csv_path = get_table_config_and_path(
        title="Interconnectors",
        sector="Grids",
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]
    selected_types = ["ITCN"]
    selected_classes = ["Link"]

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

        filtered_df = set_general_filter_df(
            df=intercon_df,
            filter_col=table_config["filter_col"],
            selected_types=selected_types,
        )

        if grid_selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[
                filtered_df["country"].isin(grid_selected_countries)
            ]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        edited_df, to_save = create_editable_df(
            filtered_df,
            edited_df_key,
            has_changes_key,
        )

        has_changes = st.session_state.get(has_changes_key, False)
        if to_save:
            render_save_button(
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
    selected_scenario = st.session_state.input_sce1
    all_countries = get_all_countries()
    tech_df = set_available_technology_df(selected_sector, input_config)

    st.subheader(f":material/timeline: Scenario specific input  | {sector_title}")

    sector_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="power_scenario_pills",
    )

    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df,
        key="power_scenario",
    )

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

    st.subheader(":material/timeline: Load")
    demand_selected_countries = render_countries_pills(
        all_countries=all_countries,
        key="demand_power_scenario_pills",
    )

    selected_types = render_demand_profiles_selectbox(selected_sector="Power")

    render_input_load_section(
        title="Power_loads",
        selected_types=selected_types,
        input_config=input_config,
        selected_classes=["Load"],
        selected_countries=demand_selected_countries,
        selected_scenario=selected_scenario,
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
    render_interconnections_section(
        grid_selected_countries,
        grid_df,
        input_config,
        selected_scenario,
    )
