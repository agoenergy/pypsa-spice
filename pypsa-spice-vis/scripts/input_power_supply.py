# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Power - Supply page under Input section.

Page shows editable power supply related
dataframes and visualisations from the modelling inputs.
"""

import pandas as pd
import streamlit as st

from scripts.getters import Getters
from scripts.input_st_handler import DFWidgetsHandler

pd.set_option("future.no_silent_downcasting", True)


def main(getters):
    """Render the Power - Supply page."""
    df_widgets_handler = DFWidgetsHandler()

    input_ui_handler = df_widgets_handler.input_ui_handler

    dfs = df_widgets_handler.load_all_dfs()
    csvs_dict = df_widgets_handler.csvs_dict

    # Get all unique countries from all dataframes and remove "world"
    all_countries = set()
    for df in [
        dfs["decomission_capacity_df"],
        dfs["fuel_costs_df"],
        dfs["intercon_df"],
        dfs["generator_df"],
        dfs["links_df"],
        dfs["storageunit_df"],
        dfs["store_df"],
    ]:
        all_countries.update(getters.get_country_list(df))

    # =============================== Supply Side ====================================
    # ============================== Render the UI part ==============================

    st.header(":material/bolt: Supply")

    col11, col12 = st.columns([1, 1])
    col21, col22 = st.columns([1, 1])
    with col11:
        if all_countries:
            selected_countries = st.pills(
                "Select Countries (can be multiple countries):",
                options=sorted(all_countries),
                default=all_countries,
                help="Select countries to filter the data.",
                selection_mode="multi",
                key="supply_selection_pills",
            )
        else:
            selected_countries = None
            st.info("No countries found")

    with col12:
        scenario_options = getters.get_input_scenario_list()
        # Use the first scenario as default if no scenario is set
        if "scenario" not in st.session_state:
            st.session_state.scenario = scenario_options if scenario_options else None
        selected_scenario = st.pills(
            "Select Scenario:",
            options=scenario_options,
            default=scenario_options[0],
            help="Select scenario to view/edit data.",
            selection_mode="single",
            key="supply_scenario_pills",
        )

        if selected_scenario:
            st.session_state.scenario = selected_scenario

    df_widgets_handler.reload_scenario_dfs(dfs, selected_scenario)

    types = getters.get_mapping_list(dfs["tech_df"])
    tech_mapping = dict(
        zip(dfs["tech_df"]["technology"], dfs["tech_df"]["technology_nomenclature"])
    )
    types_full_names = [tech_mapping.get(t, t) for t in types]
    types_full_names.sort()

    with col21:
        reverse_mapping = {v: k for k, v in tech_mapping.items()}
        default_type_selection = [types_full_names[0]]

        selected_type_full = st.multiselect(
            "Select Technology types:",
            types_full_names,
            default=default_type_selection,
        )

        if not selected_type_full:
            st.warning(
                "At least one technology type must be selected. Resetting to default."
            )
            selected_type_full = default_type_selection

        selected_types = [reverse_mapping.get(v, v) for v in selected_type_full]
        selected_types_str = ", ".join(selected_types)
    with col22:
        selected_classes = (
            dfs["tech_df"]
            .loc[dfs["tech_df"]["technology"].isin(selected_types), "class"]
            .unique()
            .tolist()
        )
        selected_classes_str = ", ".join(selected_classes)
        st.markdown(f"Tech: **{selected_types_str}**")
        st.markdown(f"Class: **{selected_classes_str}**")

    # Technologies
    input_ui_handler.set_up_single_tab_widget(
        "technologies",
        dfs["tech_df"],
        selected_types,
        csvs_dict["technologies"].path,
        selected_countries,
    )

    # Availability
    input_ui_handler.set_up_double_tab_widget(
        "availability",
        dfs["avail_df"],
        selected_types,
        csvs_dict["availability"].path,
        selected_countries,
        secondary_df=dfs["tech_df"],
    )

    # Power plant costs
    input_ui_handler.set_up_double_tab_widget(
        "pp_costs",
        dfs["pp_costs_df"],
        selected_types,
        csvs_dict["pp_costs"].path,
        selected_countries,
    )

    # Potentials
    input_ui_handler.set_up_single_tab_widget(
        "potentials",
        dfs["potentials_df"],
        selected_types,
        csvs_dict["potentials"].path,
        selected_countries,
    )

    # Storage Costs
    input_ui_handler.set_up_single_tab_widget(
        "storage_cost",
        dfs["storage_cost_df"],
        selected_types,
        csvs_dict["storage_cost"].path,
        selected_countries,
    )

    # Storage Inflows
    input_ui_handler.set_up_single_tab_widget(
        "storage_inflows",
        dfs["storage_inflows_df"],
        selected_types,
        csvs_dict["storage_inflows"].path,
        selected_countries,
    )

    # Fuel costs
    input_ui_handler.set_up_single_tab_widget(
        "fuel_costs",
        dfs["fuel_costs_df"],
        selected_types,
        csvs_dict["fuel_costs"].path,
        selected_countries,
        selected_classes=selected_classes,
        secondary_df=dfs["tech_df"],
    )

    # Asset - Generators
    if "Generator" in selected_classes:
        input_ui_handler.set_up_single_tab_widget(
            "generator",
            dfs["generator_df"],
            selected_types,
            csvs_dict["generator"].path,
            selected_countries,
        )

    if "Storage Unit" in selected_classes:
        input_ui_handler.set_up_single_tab_widget(
            "storageunit",
            dfs["storageunit_df"],
            selected_types,
            csvs_dict["storageunit"].path,
            selected_countries,
        )

    if "Store" in selected_classes:
        input_ui_handler.set_up_single_tab_widget(
            "store",
            dfs["store_df"],
            selected_types,
            csvs_dict["store"].path,
            selected_countries,
        )

    input_ui_handler.set_up_single_tab_widget(
        "links",
        dfs["links_df"],
        selected_types,
        csvs_dict["links"].path,
        selected_countries,
    )

    # Asset - Decomissioning
    input_ui_handler.set_up_single_tab_widget(
        "decomission",
        dfs["decomission_capacity_df"],
        selected_types,
        csvs_dict["decomission"].path,
        selected_countries,
    )

    # ============================== Interconnection =================================
    # ============================== Render the UI part ==============================

    st.header(":material/diagonal_line: Interconnections")

    col11, col12, col13 = st.columns([1, 1, 1])
    with col11:
        if all_countries:
            base_country = st.pills(
                "Select Base Country (only one country):",
                options=sorted(all_countries),
                default=sorted(all_countries)[0],
                help="Select a base country to filter the data.",
                selection_mode="single",
                key="intercon_selection_pills",
            )
        else:
            base_country = None
            st.info("No countries found")

    with col12:
        scenario_options = getters.get_input_scenario_list()

        # Use the first scenario as default if no scenario is set
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
            key="intercon_scenario_pills",
        )

        if selected_scenario:
            st.session_state.scenario = selected_scenario
    with col13:
        st.markdown("Tech: **ITCN**")
        st.markdown("Class: **Link**")

    df_widgets_handler.reload_scenario_dfs(dfs, selected_scenario)

    # Interconnectors
    input_ui_handler.set_up_single_tab_widget(
        csv_dict_key="interconnector",
        input_df=dfs["intercon_df"],
        selected_types=["ITCN"],
        input_csv_path=csvs_dict["interconnector"].path,
        selected_countries=[base_country],
    )


if __name__ == "__main__":
    getters = Getters()
    main(getters)
