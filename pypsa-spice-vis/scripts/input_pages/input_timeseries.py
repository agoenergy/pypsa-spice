# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Timeseries page under Input section."""

import pandas as pd
import streamlit as st

from scripts.data_utils import (
    render_countries_n_scenario_pills,
    render_widgets_from_config,
)
from scripts.input_st_handler import DataFrameWidgetsHandler

pd.set_option("future.no_silent_downcasting", True)


TIMESERIES_SUPPLY_WIDGETS = [
    {
        "csv_key": "availability",
        "widget": "double",
        "extra_kwargs": lambda render_context: {
            "secondary_df": render_context["dfs"]["technologies_df"]
        },
    },
    {"csv_key": "storage_inflows", "widget": "double"},
]

TIMESERIES_DEMAND_WIDGETS = [
    {"csv_key": "demand", "widget": "double"},
]


def main(get_params):
    """Render the Timeseries input page."""
    df_widgets_handler = DataFrameWidgetsHandler()

    dfs = df_widgets_handler.load_all_dfs()
    if not dfs:
        return

    all_countries = get_params.init_config["base_configs"]["regions"].keys()

    st.header("Supply Profiles")
    selected_countries, selected_scenario = render_countries_n_scenario_pills(
        get_params=get_params,
        all_countries=all_countries,
        key="timeseries_scenario_pills",
    )

    df_widgets_handler.reload_scenario_dfs(dfs, selected_scenario)

    tech_df = dfs["technologies_df"]

    types = get_params.get_mapping_list(tech_df)
    tech_mapping = dict(zip(tech_df["technology"], tech_df["technology_nomenclature"]))
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

    supply_render_context = {
        "dfs": dfs,
        "selected_types": selected_supply_types,
        "selected_countries": selected_countries,
    }

    render_widgets_from_config(
        input_ui_handler=df_widgets_handler.input_ui_handler,
        csvs_dict=df_widgets_handler.csvs_dict,
        widget_configs=TIMESERIES_SUPPLY_WIDGETS,
        render_context=supply_render_context,
    )

    st.header("Demand Profiles")

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

    demand_render_context = {
        "dfs": dfs,
        "selected_types": selected_profile_types,
        "selected_countries": selected_countries,
    }

    render_widgets_from_config(
        input_ui_handler=df_widgets_handler.input_ui_handler,
        csvs_dict=df_widgets_handler.csvs_dict,
        widget_configs=TIMESERIES_DEMAND_WIDGETS,
        render_context=demand_render_context,
    )
