# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Static page under Input section."""

import os

import pandas as pd
import streamlit as st
import yaml

from scripts.data_utils import render_countries_n_scenario_pills
from scripts.input_st_handler import DataFrameWidgetsHandler


def render_type_and_class_filters(
    tech_df: pd.DataFrame,
    get_params,
    key: str = "default",
) -> tuple[list, list]:
    """Render type and class controls, and return selected values."""
    col21, col22 = st.columns([1, 1])
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
            key=f"type_filter_multiselect_{key}",
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


def render_interconnections_widget(
    grid_selected_countries: list,
    df_widgets_handler: DataFrameWidgetsHandler,
    grid_df: pd.DataFrame,
):
    """Render interconnections section."""
    st.selectbox(
        "Select type of interconnection:",
        options=["Transmission", "Distribution"],
        index=0,
        key="static_intercon_type",
    )
    st.markdown("Class: **Link**")
    grid_type = st.session_state.get("static_intercon_type", None)

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

    df_widgets_handler.set_up_df_with_charts(
        sector="Grids",
        title="Interconnectors",
        input_df=intercon_df,
        selected_types=["ITCN"],
        selected_classes=["Link"],
        selected_countries=grid_selected_countries,
    )


def main(get_params):
    """Render the Static input page."""
    with open(
        os.path.join(st.session_state.current_dir, "setting/input_settings.yaml"),
        encoding="utf-8",
    ) as file:
        config = yaml.safe_load(file)["Static"]

    df_widgets_handler = DataFrameWidgetsHandler(config=config)

    tech_path = os.path.join(
        df_widgets_handler.base_input_path,
        "global_input",
        config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)

    all_countries = get_params.init_config["base_configs"]["regions"].keys()

    # ============================== Global input ==============================
    st.header(":globe_with_meridians: Global input")
    dfs = df_widgets_handler.load_all_dfs()
    if not dfs:
        return

    selected_countries, selected_scenario = render_countries_n_scenario_pills(
        get_params=get_params,
        all_countries=all_countries,
        key="global_static_pills",
    )

    # Reload dfs after the scenario is selected
    dfs = df_widgets_handler.load_all_dfs(selected_scenario=selected_scenario)

    selected_types, selected_classes = render_type_and_class_filters(
        tech_df, get_params, key="global"
    )

    for title in config["Global_input"].keys():
        df_widgets_handler.set_up_df_with_charts(
            sector="Global_input",
            title=title,
            input_df=dfs[title + "_df"],
            selected_types=selected_types,
            selected_classes=selected_classes,
            selected_countries=selected_countries,
        )

    # ============================== Power ==============================
    st.header(":material/bolt: Power")

    power_dfs = df_widgets_handler.load_all_dfs(specific_sector="Power")

    if not power_dfs:
        return
    power_selected_countries, power_selected_scenario = (
        render_countries_n_scenario_pills(
            get_params=get_params,
            all_countries=all_countries,
            key="power_static_pills",
        )
    )

    # Reload power_dfs after the scenario is selected
    power_dfs = df_widgets_handler.load_all_dfs(
        specific_sector="Power", selected_scenario=power_selected_scenario
    )
    power_selected_types, power_selected_classes = render_type_and_class_filters(
        tech_df, get_params, key="power"
    )

    for title in config["Power"].keys():
        df_widgets_handler.set_up_df_with_charts(
            sector="Power",
            title=title,
            input_df=power_dfs[title + "_df"],
            selected_types=power_selected_types,
            selected_classes=power_selected_classes,
            selected_countries=power_selected_countries,
        )

    # ===================== Industry (if exists in base_config) =====================
    if df_widgets_handler.has_industry_sector:
        st.header(":material/factory: Industry sector")
        industry_selected_countries, industry_selected_scenario = (
            render_countries_n_scenario_pills(
                get_params=get_params,
                all_countries=all_countries,
                key="industry_static_pills",
            )
        )

        industry_dfs = df_widgets_handler.load_all_dfs(specific_sector="Industry")
        if not industry_dfs:
            return

        # Reload dfs after the scenario is selected
        industry_dfs = df_widgets_handler.load_all_dfs(
            specific_sector="Industry", selected_scenario=industry_selected_scenario
        )
        industry_selected_types, industry_selected_classes = (
            render_type_and_class_filters(tech_df, get_params, key="industry")
        )

        for title in config["Industry"].keys():
            df_widgets_handler.set_up_df_with_charts(
                sector="Industry",
                title=title,
                input_df=industry_dfs[title + "_df"],
                selected_types=industry_selected_types,
                selected_classes=industry_selected_classes,
                selected_countries=industry_selected_countries,
            )

    # ===================== Transport (if exists in base_config) =====================
    if df_widgets_handler.has_transport_sector:
        st.header(":material/factory: Transport sector")
        transport_selected_countries, transport_selected_scenario = (
            render_countries_n_scenario_pills(
                get_params=get_params,
                all_countries=all_countries,
                key="transport_static_pills",
            )
        )

        transport_dfs = df_widgets_handler.load_all_dfs(specific_sector="Transport")
        if not transport_dfs:
            return

        # Reload dfs after the scenario is selected
        transport_dfs = df_widgets_handler.load_all_dfs(
            specific_sector="Transport", selected_scenario=transport_selected_scenario
        )
        transport_selected_types, transport_selected_classes = (
            render_type_and_class_filters(tech_df, get_params, key="transport")
        )

        for title in config["Transport"].keys():
            df_widgets_handler.set_up_df_with_charts(
                sector="Transport",
                title=title,
                input_df=transport_dfs[title + "_df"],
                selected_types=transport_selected_types,
                selected_classes=transport_selected_classes,
                selected_countries=transport_selected_countries,
            )

    # ===================== Interconnections =====================
    st.header(":material/diagonal_line: Interconnections")
    grid_selected_countries, grid_selected_scenario = render_countries_n_scenario_pills(
        get_params=get_params,
        all_countries=all_countries,
        key="grid_static_pills",
    )
    grid_path = os.path.join(
        df_widgets_handler.base_input_path,
        grid_selected_scenario,
        "power",
        config["Grids"]["Interconnectors"]["csv_name"],
    )
    grid_df = pd.read_csv(grid_path)

    render_interconnections_widget(grid_selected_countries, df_widgets_handler, grid_df)
