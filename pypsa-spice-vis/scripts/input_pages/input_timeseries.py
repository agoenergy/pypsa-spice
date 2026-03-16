# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Timeseries page under Input section."""

import os

import pandas as pd
import streamlit as st
import yaml

from scripts.data_utils import render_countries_n_scenario_pills
from scripts.input_st_handler import DataFrameWidgetsHandler


def main(get_params):
    """Render the Timeseries input page."""
    with open(
        os.path.join(st.session_state.current_dir, "setting/input_settings.yaml"),
        encoding="utf-8",
    ) as file:
        config = yaml.safe_load(file)

    tech_config = config["Static"]
    time_series_config = config["Time_series"]

    df_widgets_handler = DataFrameWidgetsHandler(config=time_series_config)
    all_countries = get_params.init_config["base_configs"]["regions"].keys()

    # ===================== Supply Profiles =====================
    st.header("Supply Profiles")
    dfs = df_widgets_handler.load_all_dfs()
    if not dfs:
        return

    selected_countries, selected_scenario = render_countries_n_scenario_pills(
        get_params=get_params,
        all_countries=all_countries,
        key="supply_timeseries_pills",
    )

    # Reload dfs after the scenario is selected
    dfs = df_widgets_handler.load_all_dfs(selected_scenario=selected_scenario)

    tech_path = os.path.join(
        df_widgets_handler.base_input_path,
        "global_input",
        tech_config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)

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

    for title in time_series_config["Global_input"].keys():
        if "demand" not in time_series_config["Global_input"][title]["tag_name"]:
            df_widgets_handler.set_up_df_with_charts(
                sector="Global_input",
                title=title,
                input_df=dfs[title + "_df"],
                selected_types=selected_supply_types,
                selected_countries=selected_countries,
            )

    # ===================== Demand Profiles =====================
    st.header("Demand Profiles")
    demand_dfs = df_widgets_handler.load_all_dfs()
    if not demand_dfs:
        return

    demand_selected_countries, demand_selected_scenario = (
        render_countries_n_scenario_pills(
            get_params=get_params,
            all_countries=all_countries,
            key="demand_timeseries_pills",
        )
    )

    # Reload demand_dfs after the scenario is selected
    dfs = df_widgets_handler.load_all_dfs(selected_scenario=demand_selected_scenario)

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

    df_widgets_handler.set_up_df_with_charts(
        sector="Global_input",
        title="Demand_Profiles",
        input_df=dfs["Demand_Profiles_df"],
        selected_types=selected_profile_types,
        selected_countries=demand_selected_countries,
    )
