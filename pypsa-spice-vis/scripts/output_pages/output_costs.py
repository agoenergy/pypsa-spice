# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Costs page under Results section.

Page shows costs related dataframes and visualisations from the modelling results.
"""

import streamlit as st

from scripts.data_utils import (
    add_nice_names,
    clean_df_for_plotting,
    normalize_dataframe,
    prepare_y_range,
    read_result_csv,
)
from scripts.output_st_handler import (
    generate_sidebar,
    render_chart_layout,
    render_download_with_data_table,
    render_section_header,
    setup_country_filter,
)
from scripts.plot_functions import (
    create_nice_names_and_color_mapping,
    plot_simple_bar_yearly,
    plot_simple_line_yearly,
)

# =========================== Render functions for each chart =========================


def render_c1_pow_capex_by_type(graph_config: dict) -> None:
    """Render power sector capex by technology type (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c2_pow_opex_cost_by_type(graph_config: dict) -> None:
    """Render power sector opex by technology type (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c3_pow_inv_cost_by_type(graph_config: dict) -> None:
    """Render power sector overnight investment by technology type (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c4_ind_capex_by_type(graph_config: dict) -> None:
    """Render industry sector capex by technology type (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c5_tra_capex_by_type(graph_config: dict) -> None:
    """Render transport sector capex by technology type (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c6_ene_capex_by_type(graph_config: dict) -> None:
    """Render all sectors capex by technology type (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c7_opex_capex_yearly(graph_config: dict) -> None:
    """Render all sectors opex_capex_yearly (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_c8_ene_avg_fuel_costs_fuel(graph_config: dict) -> None:
    """Render all sectors ene_avg_fuel_costs_fuel_yearly (yearly line chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_line_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


if __name__ == "__main__":
    st.title(":material/attach_money: Costs")
    output_config = st.session_state.output_config["costs"]
    sector = str(st.session_state.get("sector", "")).lower()
    show_industry = "i" in sector
    show_transport = "t" in sector

    # Decide which charts should appear (and therefore be in the sidebar)
    COST_CHART_KEYS = ["c1", "c2", "c3"]
    if show_industry:
        COST_CHART_KEYS.append("c4")
    if show_transport:
        COST_CHART_KEYS.append("c5")
    COST_CHART_KEYS += ["c6", "c7", "c8"]

    # Render (power always)
    render_c1_pow_capex_by_type(output_config["c1"])
    render_c2_pow_opex_cost_by_type(output_config["c2"])
    render_c3_pow_inv_cost_by_type(output_config["c3"])

    if show_industry:
        render_c4_ind_capex_by_type(output_config["c4"])
    if show_transport:
        render_c5_tra_capex_by_type(output_config["c5"])

    # Always show energy & totals
    render_c6_ene_capex_by_type(output_config["c6"])
    render_c7_opex_capex_yearly(output_config["c7"])
    render_c8_ene_avg_fuel_costs_fuel(output_config["c8"])

    # Sidebar (only include charts present in config)
    cost_charts = [output_config[k] for k in COST_CHART_KEYS if k in output_config]
    generate_sidebar([c["name"] for c in cost_charts])
