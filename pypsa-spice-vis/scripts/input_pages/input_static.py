# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create Static page under Input section."""

import os

import pandas as pd
import streamlit as st

from scripts.data_utils import render_countries_n_scenario_pills
from scripts.input_st_handler import DataFrameWidgetsHandler


def get_input_scenario_list(base_config: dict) -> list[str]:
    """Return available input scenarios with the default scenario first."""
    data_folder_path = base_config["input_folder_path"]

    if not os.path.exists(data_folder_path):
        raise FileNotFoundError(f"folder not found: {data_folder_path}")

    scenario_list = [
        scenario
        for scenario in os.listdir(data_folder_path)
        if scenario not in [".DS_Store"] and scenario != "global_input"
    ]

    for sce in (base_config["path_configs"]["input_scenario_name"], ""):
        if sce in scenario_list:
            scenario_list.insert(0, scenario_list.pop(scenario_list.index(sce)))

    return scenario_list


def get_mapping_list(*dfs: pd.DataFrame) -> list[str]:
    """Get sorted technology/profile types from one or more dataframes."""
    type_set = set()

    for df in dfs:
        if "technology" in df.columns:
            type_set |= set(df["technology"].unique())
        if "profile_type" in df.columns:
            type_set |= set(df["profile_type"].unique())

    return sorted(type_set)


def render_type_and_class_filters(
    tech_df: pd.DataFrame,
    key: str = "default",
) -> tuple[list, list]:
    """Render type and class controls, and return selected values."""
    col21, col22 = st.columns([1, 1])
    types = get_mapping_list(tech_df)
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


if __name__ == "__main__":
    base_config = st.session_state.base_config
    input_config = st.session_state.input_config["Static"]
    selected_sector = st.session_state.get("selected_input_sector", "Power")

    df_widgets_handler = DataFrameWidgetsHandler(input_config=input_config)

    tech_path = os.path.join(
        st.session_state.input_path,
        "global_input",
        input_config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)

    all_countries = list(base_config["base_configs"]["regions"].keys())

    st.subheader(f":material/database: Static | {selected_sector}")

    sector_lower = selected_sector.lower()
    sector_selected_countries, sector_selected_scenario = (
        render_countries_n_scenario_pills(
            scenario_options=get_input_scenario_list(base_config),
            all_countries=all_countries,
            key=f"{sector_lower}_static_pills",
        )
    )

    sector_dfs = df_widgets_handler.load_all_dfs(
        specific_sector=selected_sector, selected_scenario=sector_selected_scenario
    )
    sector_selected_types, sector_selected_classes = render_type_and_class_filters(
        tech_df, key=f"{sector_lower}_static"
    )

    for title in input_config[selected_sector].keys():
        df_widgets_handler.set_up_df_with_charts(
            sector=selected_sector,
            title=title,
            input_df=sector_dfs[title + "_df"],
            selected_types=sector_selected_types,
            selected_classes=sector_selected_classes,
            selected_countries=sector_selected_countries,
        )

    if selected_sector == "Power":
        st.header(":material/diagonal_line: Interconnections")
        grid_selected_countries, grid_selected_scenario = (
            render_countries_n_scenario_pills(
                scenario_options=get_input_scenario_list(base_config),
                all_countries=all_countries,
                key="grid_static_pills",
            )
        )
        grid_path = os.path.join(
            df_widgets_handler.base_input_path,
            grid_selected_scenario,
            "power",
            input_config["Grids"]["Interconnectors"]["csv_name"],
        )
        grid_df = pd.read_csv(grid_path)

        render_interconnections_widget(
            grid_selected_countries, df_widgets_handler, grid_df
        )
