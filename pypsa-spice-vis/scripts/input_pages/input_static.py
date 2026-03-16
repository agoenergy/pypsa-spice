# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Static page under Input section."""

import streamlit as st

from scripts.data_utils import (
    render_countries_n_scenario_pills,
    render_widgets_from_config,
)
from scripts.input_st_handler import DataFrameWidgetsHandler

STATIC_WIDGETS = [
    {"csv_key": "technologies", "widget": "single"},
    {"csv_key": "pp_costs", "widget": "double"},
    {"csv_key": "potentials", "widget": "single"},
    {"csv_key": "storage_cost", "widget": "double"},
    {"csv_key": "links", "widget": "single"},
    {"csv_key": "decomission", "widget": "single"},
    {"csv_key": "load", "widget": "double"},
]

POWER_CLASS_DEPENDENT_WIDGETS = {
    "Generator": {"csv_key": "generator"},
    "Storage Unit": {"csv_key": "storageunit"},
    "Store": {"csv_key": "store"},
}


def render_type_and_class_filters(dfs: dict, get_params) -> tuple[list, list]:
    """Render type and class controls, and return selected values."""
    col21, col22 = st.columns([1, 1])
    tech_df = dfs["technologies_df"]

    types = get_params.get_mapping_list(tech_df)
    tech_mapping = dict(zip(tech_df["technology"], tech_df["technology_nomenclature"]))
    types_full_names = sorted([tech_mapping.get(t, t) for t in types])

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

    with col22:
        selected_classes = (
            tech_df.loc[tech_df["technology"].isin(selected_types), "class"]
            .unique()
            .tolist()
        )
        st.markdown(f"Tech: **{', '.join(selected_types)}**")
        st.markdown(f"Class: **{', '.join(selected_classes)}**")

    return selected_types, selected_classes


def render_class_dependent_widgets(
    df_widgets_handler: DataFrameWidgetsHandler,
    dfs: dict,
    selected_types: list,
    selected_classes: list,
    selected_countries: list,
):
    """Render widgets that depend on selected technology classes."""
    for class_name, class_widget in POWER_CLASS_DEPENDENT_WIDGETS.items():
        key = class_widget["csv_key"]
        if class_name in selected_classes:
            df_widgets_handler.input_ui_handler.set_up_single_tab_widget(
                key,
                dfs[key + "_df"],
                selected_types,
                df_widgets_handler.csvs_dict[key].path,
                selected_countries,
            )


def render_interconnections_widget(
    all_countries: set,
    df_widgets_handler: DataFrameWidgetsHandler,
    dfs: dict,
):
    """Render interconnections section."""
    col21, col22 = st.columns([1, 1])

    with col21:
        if all_countries:
            st.selectbox(
                "Select Base Country (single country):",
                options=sorted(all_countries),
                index=0,
                key="static_intercon_country",
            )
        else:
            st.info("No countries found")

        st.markdown("Tech: **ITCN**")

    with col22:
        st.selectbox(
            "Select type of interconnection:",
            options=["Transmission", "Distribution"],
            index=0,
            key="static_intercon_type",
        )
        st.markdown("Class: **Link**")

    base_country = st.session_state.get("static_intercon_country", None)
    grid_type = st.session_state.get("static_intercon_type", None)

    if grid_type == "Transmission":
        intercon_df = dfs["interconnector_df"][
            (dfs["interconnector_df"]["bus0"].str.contains("HVELEC"))
            & (dfs["interconnector_df"]["bus1"].str.contains("HVELEC"))
        ]
    else:
        intercon_df = dfs["interconnector_df"][
            (dfs["interconnector_df"]["bus0"].str.contains("LVELEC"))
            & (dfs["interconnector_df"]["bus1"].str.contains("LVELEC"))
        ]

    df_widgets_handler.input_ui_handler.set_up_single_tab_widget(
        csv_dict_key="interconnector",
        input_df=intercon_df,
        selected_types=["ITCN"],
        input_csv_path=df_widgets_handler.csvs_dict["interconnector"].path,
        selected_countries=[base_country],
    )


def main(get_params):
    """Render the Static input page."""
    df_widgets_handler = DataFrameWidgetsHandler()

    dfs = df_widgets_handler.load_all_dfs()
    if not dfs:
        return

    all_countries = get_params.init_config["base_configs"]["regions"].keys()

    st.header(":globe_with_meridians: Global input")
    selected_countries, selected_scenario = render_countries_n_scenario_pills(
        get_params=get_params,
        all_countries=all_countries,
        key="static_scenario_pills",
    )

    df_widgets_handler.reload_scenario_dfs(dfs, selected_scenario)
    selected_types, selected_classes = render_type_and_class_filters(dfs, get_params)

    render_context = {
        "input_ui_handler": df_widgets_handler.input_ui_handler,
        "dfs": dfs,
        "csvs_dict": df_widgets_handler.csvs_dict,
        "selected_types": selected_types,
        "selected_countries": selected_countries,
        "selected_classes": selected_classes,
    }

    render_widgets_from_config(
        input_ui_handler=df_widgets_handler.input_ui_handler,
        csvs_dict=df_widgets_handler.csvs_dict,
        widget_configs=STATIC_WIDGETS,
        render_context=render_context,
    )

    st.header(":material/bolt: Power")
    render_class_dependent_widgets(
        df_widgets_handler,
        dfs,
        selected_types,
        selected_classes,
        selected_countries,
    )

    st.header(":material/diagonal_line: Interconnections")
    render_interconnections_widget(all_countries, df_widgets_handler, dfs)
