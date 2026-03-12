# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Timeseries page under Input section."""

import pandas as pd
import streamlit as st

from scripts.input_st_handler import DFWidgetsHandler

pd.set_option("future.no_silent_downcasting", True)


TIMESERIES_SUPPLY_WIDGETS = [
    {
        "csv_key": "availability",
        "df_key": "avail_df",
        "widget": "double",
        "extra_kwargs": lambda dfs: {"secondary_df": dfs["tech_df"]},
    },
    {
        "csv_key": "storage_inflows",
        "df_key": "storage_inflows_df",
        "widget": "single",
    },
]

TIMESERIES_DEMAND_WIDGETS = [
    {"csv_key": "demand", "df_key": "demand_df", "widget": "double"},
    {"csv_key": "load", "df_key": "load_df", "widget": "double"},
]


def _render_timeseries_widgets(
    input_ui_handler,
    widget_configs,
    dfs,
    csvs_dict,
    selected_types,
    selected_countries,
):
    """Render one configured timeseries widget group."""
    for widget_config in widget_configs:
        render_fn = (
            input_ui_handler.set_up_double_tab_widget
            if widget_config["widget"] == "double"
            else input_ui_handler.set_up_single_tab_widget
        )

        extra_kwargs = {}
        if "extra_kwargs" in widget_config:
            extra_kwargs = widget_config["extra_kwargs"](dfs)

        csv_key = widget_config["csv_key"]
        df_key = widget_config["df_key"]

        render_fn(
            csv_key,
            dfs[df_key],
            selected_types,
            csvs_dict[csv_key].path,
            selected_countries,
            **extra_kwargs,
        )


def main(getters):
    """Render the Timeseries input page."""
    df_widgets_handler = DFWidgetsHandler()
    input_ui_handler = df_widgets_handler.input_ui_handler

    dfs = df_widgets_handler.load_all_dfs()
    if not dfs:
        return

    csvs_dict = df_widgets_handler.csvs_dict

    all_countries = set()
    for df in [
        dfs["avail_df"],
        dfs["demand_df"],
        dfs["load_df"],
        dfs["storage_inflows_df"],
    ]:
        all_countries.update(getters.get_country_list(df))

    st.header(":material/timeline: Timeseries")

    col11, col12 = st.columns([1, 1])

    with col11:
        if all_countries:
            selected_countries = st.pills(
                "Select Countries:",
                options=sorted(all_countries),
                default=sorted(all_countries),
                help="Select countries to filter the data.",
                selection_mode="multi",
                key="timeseries_selection_pills",
            )
        else:
            selected_countries = None
            st.info("No countries found")

    with col12:
        scenario_options = getters.get_input_scenario_list()
        if "scenario" not in st.session_state:
            st.session_state.scenario = (
                scenario_options[0] if scenario_options else None
            )

        selected_scenario = st.pills(
            "Select Scenario:",
            options=scenario_options,
            default=(
                st.session_state.scenario
                if st.session_state.scenario in scenario_options
                else scenario_options[0]
            ),
            help="Select scenario to view/edit data.",
            selection_mode="single",
            key="timeseries_scenario_pills",
        )

        if selected_scenario:
            st.session_state.scenario = selected_scenario

    df_widgets_handler.reload_scenario_dfs(dfs, selected_scenario)

    st.subheader("Supply Profiles")

    types = getters.get_mapping_list(dfs["tech_df"])
    tech_mapping = dict(
        zip(dfs["tech_df"]["technology"], dfs["tech_df"]["technology_nomenclature"])
    )
    types_full_names = [tech_mapping.get(t, t) for t in types]
    types_full_names.sort()
    reverse_mapping = {v: k for k, v in tech_mapping.items()}

    default_supply_selection = [types_full_names[0]] if types_full_names else []
    selected_supply_types_full = st.multiselect(
        "Select technology types for supply timeseries:",
        types_full_names,
        default=default_supply_selection,
        key="timeseries_supply_types",
    )
    if not selected_supply_types_full and default_supply_selection:
        selected_supply_types_full = default_supply_selection
    selected_supply_types = [
        reverse_mapping.get(v, v) for v in selected_supply_types_full
    ]

    _render_timeseries_widgets(
        input_ui_handler=input_ui_handler,
        widget_configs=TIMESERIES_SUPPLY_WIDGETS,
        dfs=dfs,
        csvs_dict=csvs_dict,
        selected_types=selected_supply_types,
        selected_countries=selected_countries,
    )

    st.subheader("Demand Profiles")

    load_mapping = {
        "HV_LOAD": "Wholesale market load (High voltage level)",
        "LV_LOAD": "Building load (low/medium voltage level)",
    }
    profile_full_names = list(load_mapping.values())
    reverse_load_mapping = {v: k for k, v in load_mapping.items()}
    default_profile_selection = [profile_full_names[0]] if profile_full_names else []

    selected_profile_types_full = st.multiselect(
        "Select demand/load profiles:",
        profile_full_names,
        default=default_profile_selection,
        key="timeseries_demand_types",
    )
    if not selected_profile_types_full and default_profile_selection:
        selected_profile_types_full = default_profile_selection
    selected_profile_types = [
        reverse_load_mapping.get(v, v) for v in selected_profile_types_full
    ]

    _render_timeseries_widgets(
        input_ui_handler=input_ui_handler,
        widget_configs=TIMESERIES_DEMAND_WIDGETS,
        dfs=dfs,
        csvs_dict=csvs_dict,
        selected_types=selected_profile_types,
        selected_countries=selected_countries,
    )
