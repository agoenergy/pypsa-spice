# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Static page under Input section."""

import pandas as pd
import streamlit as st

from scripts.input_st_handler import DFWidgetsHandler

pd.set_option("future.no_silent_downcasting", True)


STATIC_WIDGETS = [
    {"csv_key": "technologies", "df_key": "tech_df", "widget": "single"},
    {"csv_key": "pp_costs", "df_key": "pp_costs_df", "widget": "double"},
    {"csv_key": "potentials", "df_key": "potentials_df", "widget": "single"},
    {
        "csv_key": "storage_cost",
        "df_key": "storage_cost_df",
        "widget": "single",
    },
    {
        "csv_key": "fuel_costs",
        "df_key": "fuel_costs_df",
        "widget": "single",
        "extra_kwargs": lambda dfs, selected_classes: {
            "selected_classes": selected_classes,
            "secondary_df": dfs["tech_df"],
        },
    },
    {"csv_key": "links", "df_key": "links_df", "widget": "single"},
    {
        "csv_key": "decomission",
        "df_key": "decomission_capacity_df",
        "widget": "single",
    },
]


CLASS_DEPENDENT_WIDGETS = {
    "Generator": {"csv_key": "generator", "df_key": "generator_df"},
    "Storage Unit": {"csv_key": "storageunit", "df_key": "storageunit_df"},
    "Store": {"csv_key": "store", "df_key": "store_df"},
}


def _render_widget(
    widget_config,
    render_context,
):
    """Render one widget from config."""
    input_ui_handler = render_context["input_ui_handler"]
    dfs = render_context["dfs"]
    csvs_dict = render_context["csvs_dict"]
    selected_types = render_context["selected_types"]
    selected_countries = render_context["selected_countries"]
    selected_classes = render_context["selected_classes"]

    render_fn = (
        input_ui_handler.set_up_double_tab_widget
        if widget_config["widget"] == "double"
        else input_ui_handler.set_up_single_tab_widget
    )

    extra_kwargs = {}
    if "extra_kwargs" in widget_config:
        extra_kwargs = widget_config["extra_kwargs"](dfs, selected_classes)

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
    """Render the Static input page."""
    df_widgets_handler = DFWidgetsHandler()
    input_ui_handler = df_widgets_handler.input_ui_handler

    dfs = df_widgets_handler.load_all_dfs()
    if not dfs:
        return

    csvs_dict = df_widgets_handler.csvs_dict

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

    st.header(":material/view_list: Static")

    col11, col12 = st.columns([1, 1])
    col21, col22 = st.columns([1, 1])

    with col11:
        if all_countries:
            selected_countries = st.pills(
                "Select Countries (can be multiple countries):",
                options=sorted(all_countries),
                default=sorted(all_countries),
                help="Select countries to filter the data.",
                selection_mode="multi",
                key="static_selection_pills",
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
            key="static_scenario_pills",
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
        default_type_selection = [types_full_names[0]] if types_full_names else []

        selected_type_full = st.multiselect(
            "Select Technology types:",
            types_full_names,
            default=default_type_selection,
        )

        if not selected_type_full and default_type_selection:
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

    render_context = {
        "input_ui_handler": input_ui_handler,
        "dfs": dfs,
        "csvs_dict": csvs_dict,
        "selected_types": selected_types,
        "selected_countries": selected_countries,
        "selected_classes": selected_classes,
    }

    for widget_config in STATIC_WIDGETS:
        _render_widget(
            widget_config=widget_config,
            render_context=render_context,
        )

    for class_name, class_widget in CLASS_DEPENDENT_WIDGETS.items():
        if class_name in selected_classes:
            input_ui_handler.set_up_single_tab_widget(
                class_widget["csv_key"],
                dfs[class_widget["df_key"]],
                selected_types,
                csvs_dict[class_widget["csv_key"]].path,
                selected_countries,
            )

    st.header(":material/diagonal_line: Interconnections")

    if all_countries:
        default_base = sorted(all_countries)[0]
        base_country = st.selectbox(
            "Select Base Country (single country):",
            options=sorted(all_countries),
            index=0,
            key="static_intercon_country",
        )
    else:
        base_country = None
        default_base = None
        st.info("No countries found")

    st.markdown("Tech: **ITCN**")
    st.markdown("Class: **Link**")

    if base_country or default_base:
        input_ui_handler.set_up_single_tab_widget(
            csv_dict_key="interconnector",
            input_df=dfs["intercon_df"],
            selected_types=["ITCN"],
            input_csv_path=csvs_dict["interconnector"].path,
            selected_countries=[base_country or default_base],
        )
