# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import streamlit as st
import pandas as pd
from scripts.getters import Getters
from scripts.input_helpers import dfWidgetsHandler

pd.set_option("future.no_silent_downcasting", True)


def main(getters):
    df_widgets_handler = dfWidgetsHandler()

    input_ui_handler = df_widgets_handler.input_ui_handler

    dfs = df_widgets_handler.load_all_dfs()
    csvs_dict = df_widgets_handler.csvs_dict

    # Get all unique countries from all dataframes and remove "world"
    all_countries = set()
    for df in [
        dfs["decomission_capacity_df"],
        dfs["fuel_costs_df"],
        dfs["intercon_df"],
        dfs["load_df"],
        dfs["generator_df"],
        dfs["links_df"],
        dfs["storageunit_df"],
        dfs["store_df"]
        ]:
        all_countries.update(getters.get_country_list(df))

    scenario_options = getters.get_project_folder_list(df_widgets_handler.base_input_path)

########################## Demand Side #################################################
########################## Render the UI part ##########################################

    st.header(":material/Lightbulb: Demand")

    col11, col12 = st.columns([1, 1])
    col21, col22 = st.columns([1, 1])
    
    with col11:
        if all_countries:
            selected_countries = st.pills(
            "Select Countries:",
            options=all_countries,
            default=all_countries,
            help="Select countries to filter the data.",
            selection_mode="multi"
            )
        else:
            selected_countries = None
            st.info("No countries found")
            
    with col12:
        scenario_options = getters.get_project_folder_list(df_widgets_handler.base_input_path)

        # Use the first scenario as default if no scenario is set
        if "scenario" not in st.session_state:
            st.session_state.scenario = scenario_options[0] if scenario_options else None

        selected_scenario = st.pills(
            "Select Scenario:",
            options=scenario_options,
            default=st.session_state.scenario 
            if st.session_state.scenario in scenario_options 
            else scenario_options[0],
            help="Select scenario to view/edit data.",
            selection_mode="single"
        )

        if selected_scenario:
            st.session_state.scenario = selected_scenario

    df_widgets_handler.reload_scenario_dfs(dfs, selected_scenario)

    with col21:
        load_mapping = {
            "HV_LOAD": "Wholesale market load (High voltage level)",
            "LV_LOAD": "Building load (low/medium voltage level)",
        }
        types_full_names = list(load_mapping.values())
        reverse_mapping = {v: k for k, v in load_mapping.items()}
        default_profile_types_selection = [types_full_names[0]]

        selected_profile_types_full = st.multiselect(
            "Select load profiles:",
            types_full_names,
            default=default_profile_types_selection,
        )

        if not selected_profile_types_full:
            st.warning("At least one load profile must be selected. Resetting to default.")
            selected_profile_types_full = default_profile_types_selection

        selected_profile_types = [reverse_mapping.get(v, v) for v in selected_profile_types_full]
        selected_profile_types_str = ", ".join(selected_profile_types)
    with col22:
        selecte_types_str = ", ".join(selected_profile_types)
        st.markdown(f"Tech: **{selecte_types_str}**")
        st.markdown(f"Class: **Load**")

    ## Demand Profiles
    input_ui_handler.set_up_double_tab_widget(
        "demand",
        dfs["demand_df"],
        selected_profile_types,
        csvs_dict["demand"].path,
        selected_countries,
        # secondary_df=dfs["tech_df"]
    )

    ## Load
    input_ui_handler.set_up_double_tab_widget(
        "load",
        dfs["load_df"],
        selected_profile_types,
        csvs_dict["load"].path,
        selected_countries
    )
    

if __name__ == "__main__":
    getters = Getters()
    main(getters)