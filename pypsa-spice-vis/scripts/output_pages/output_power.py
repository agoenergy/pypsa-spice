# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Power page under Results section.

Page shows editable power related
dataframes and visualisations from the modelling results.
"""

import os

import streamlit as st
import yaml

from scripts.data_utils import (
    add_nice_names,
    clean_df_for_plotting,
    filter_and_prepare_hourly_data,
    filter_dataframe_by_date_range,
    get_filtered_df_and_date_range,
    handle_small_values,
    load_and_validate_hourly_data,
    normalize_dataframe,
    prepare_y_range,
    prettify_label,
    read_result_csv,
)
from scripts.output_st_handler import (
    generate_colour_mapping_dict,
    generate_sidebar,
    render_chart_layout,
    render_download_with_data_table,
    render_download_without_data_table,
    render_section_header,
    setup_country_filter,
    setup_hourly_filters,
    setup_radio_button_filter,
    setup_year_filter,
)
from scripts.plot_functions import (
    create_nice_names_and_color_mapping,
    plot_area_share_yearly,
    plot_bar_with_filter,
    plot_diff_bar_yearly,
    plot_filtered_bar_hourly,
    plot_line_with_secondary_y_hourly,
    plot_simple_bar_hourly,
    plot_simple_bar_yearly,
    plot_simple_line_hourly,
)

st.title(":material/bolt: Power")

with open(
    os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"),
    encoding="utf-8",
) as file:
    config = yaml.safe_load(file)["power"]

POWER_CHART_KEYS = [
    "p1",
    "p2",
    "p3",
    "p4",
    "p6",
    "p7",
    "p8",
    "p9",
    "p10",
    "p11",
    "p12",
    "p13",
    "p14",
]

# Only include charts in the sidebar if they are present in the config
power_charts = [config[key] for key in POWER_CHART_KEYS if key in config]
table_of_content = [chart["name"] for chart in power_charts]
generate_sidebar(table_of_content)


# =========================== Shared helpers =========================


def _render_scenario_diff_chart(
    df1: "pd.DataFrame",
    df2: "pd.DataFrame",
    graph_config: dict,
    mapping_df: "pd.DataFrame | None",
) -> None:
    """Compute (df2 - df1) per year × legend and render a grouped diff bar chart."""
    import pandas as pd

    s1_pivot = df1.pivot_table(
        index="year", columns="legend", values="value", aggfunc="sum"
    ).fillna(0)
    s2_pivot = df2.pivot_table(
        index="year", columns="legend", values="value", aggfunc="sum"
    ).fillna(0)
    all_cols = s1_pivot.columns.union(s2_pivot.columns)
    diff_df = (
        s2_pivot.reindex(columns=all_cols, fill_value=0)
        - s1_pivot.reindex(columns=all_cols, fill_value=0)
    )
    diff_df = (
        diff_df.reset_index()
        .melt(id_vars="year", var_name="legend", value_name="value")
        .query("value != 0")
    )
    if diff_df.empty:
        return

    colour_mapping_diff = generate_colour_mapping_dict(
        graph_config["table_name"], mapping_df, diff_df, graph_config["leg_col"]
    )
    st.caption(
        f"Difference ({st.session_state.sce2} \u2212 {st.session_state.sce1})"
    )
    plot_diff_bar_yearly(
        diff_df,
        graph_config,
        colour_mapping_diff,
        key=(
            "plotly_chart_diff_"
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
    )


# =========================== Render functions for each chart =========================
def render_p1_capacity_by_type(graph_config: dict) -> None:
    """Render installed capacity by technology type (yearly bar chart).

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
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
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
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
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
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
    )

    if has_dual_data:
        _render_scenario_diff_chart(
            scenario_1_grouped, scenario_2_grouped, graph_config, mapping_df
        )

    st.divider()


def render_p2_capacity_by_region(graph_config: dict) -> None:
    """Render installed capacity by region with filter (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "bar_with_filter"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_filter = setup_radio_button_filter(graph_config, is_dual)

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    filter_col = graph_config["fil_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    filter_value_s1 = shared_filter or st.radio(
        f"{graph_config['slider_id'].format(st.session_state.sce1)} "
        f"Select {filter_col} ({st.session_state.sce1}):",
        options=scenario_1_raw[filter_col].unique(),
        format_func=prettify_label,
        horizontal=True,
        label_visibility="collapsed",
    )
    scenario_1_filtered = scenario_1_raw[
        scenario_1_raw[filter_col] == filter_value_s1
    ].copy()
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)
    scenario_1_filtered = clean_df_for_plotting(legend_col, scenario_1_filtered)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_filtered = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
        )
        if scenario_2_raw is not None:
            filter_value_s2 = shared_filter or st.radio(
                f"{graph_config['slider_id'].format(st.session_state.sce2)} "
                f"Select {filter_col} ({st.session_state.sce2}):",
                options=scenario_2_raw[filter_col].unique(),
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
            scenario_2_filtered = scenario_2_raw[
                scenario_2_raw[filter_col] == filter_value_s2
            ].copy()
            scenario_2_filtered = add_nice_names(
                scenario_2_filtered, legend_col, mapping_df
            )
            scenario_2_filtered = clean_df_for_plotting(legend_col, scenario_2_filtered)

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_filtered is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_bar_with_filter,
        render_download_function=render_download_with_data_table,
    )

    if has_dual_data:
        _render_scenario_diff_chart(
            scenario_1_filtered, scenario_2_filtered, graph_config, mapping_df
        )

    st.divider()


def render_p3_generation_by_type(graph_config: dict) -> None:
    """Render electricity generation by technology type (yearly bar chart).

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
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
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
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
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
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
    )

    if has_dual_data:
        _render_scenario_diff_chart(
            scenario_1_grouped, scenario_2_grouped, graph_config, mapping_df
        )

    st.divider()


def render_p4_share_category(graph_config: dict) -> None:
    """Render generation share by category (yearly area chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "area_share_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_plot = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_plot = add_nice_names(scenario_1_plot, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_plot = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
        )
        if scenario_2_raw is not None:
            scenario_2_plot = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_plot = add_nice_names(scenario_2_plot, legend_col, mapping_df)

    has_dual_data = (
        is_dual and scenario_2_plot is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_plot,
        scenario_2_vis_display_data=scenario_2_plot,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range={},  # No y_range for area charts
        plot_function=plot_area_share_yearly,
        render_download_function=render_download_with_data_table,
    )

    st.divider()


def render_p6_transmission_capacity_between_regions(graph_config: dict) -> None:
    """Render transmission capacity between regions with filter (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "bar_with_filter"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_filter = setup_radio_button_filter(graph_config, is_dual)

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    filter_col = graph_config["fil_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    filter_value_s1 = shared_filter or st.radio(
        f"{graph_config['slider_id'].format(st.session_state.sce1)} "
        f"Select {filter_col} ({st.session_state.sce1}):",
        options=scenario_1_raw[filter_col].unique(),
        format_func=prettify_label,
        horizontal=True,
        label_visibility="collapsed",
    )
    scenario_1_filtered = scenario_1_raw[
        scenario_1_raw[filter_col] == filter_value_s1
    ].copy()
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)
    scenario_1_filtered = clean_df_for_plotting(legend_col, scenario_1_filtered)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_filtered = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
        )
        if scenario_2_raw is not None:
            filter_value_s2 = shared_filter or st.radio(
                f"{graph_config['slider_id'].format(st.session_state.sce2)} "
                f"Select {filter_col} ({st.session_state.sce2}):",
                options=scenario_2_raw[filter_col].unique(),
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
            scenario_2_filtered = scenario_2_raw[
                scenario_2_raw[filter_col] == filter_value_s2
            ].copy()
            scenario_2_filtered = add_nice_names(
                scenario_2_filtered, legend_col, mapping_df
            )
            scenario_2_filtered = clean_df_for_plotting(legend_col, scenario_2_filtered)

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, "year")

    has_dual_data = (
        is_dual and scenario_2_filtered is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_bar_with_filter,
        render_download_function=render_download_with_data_table,
    )

    st.divider()


def render_p7_hourly_generation(graph_config: dict) -> None:
    """Render hourly electricity generation (hourly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )
    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario 1 data
    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    # Load and validate scenario 2 data (if dual mode)
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare scenario 1 data
    scenario_1_filtered, start_date, end_date, is_complete = (
        filter_and_prepare_hourly_data(
            scenario_1_raw, graph_config, legend_col, mapping_df
        )
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    # Filter and prepare scenario 2 data (if dual mode)
    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = filter_and_prepare_hourly_data(
            scenario_2_raw, graph_config, legend_col, mapping_df
        )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, "snapshot")

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    render_chart_layout(
        is_dual_scenario=has_dual_filters and scenario_2_filtered is not None,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        # Use filtered data for download in hourly charts
        table_1_display_data=scenario_1_filtered,
        table_2_display_data=scenario_2_filtered,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_hourly,
        render_download_function=render_download_without_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
        **plot_kwargs,
    )

    st.divider()


def render_p8_regional_hourly_generation(graph_config: dict) -> None:
    """Render regional hourly generation with filter (hourly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "filtered_bar_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data
    scenario_1_filtered, start_date, end_date, is_complete = (
        filter_and_prepare_hourly_data(
            scenario_1_raw, graph_config, legend_col, mapping_df
        )
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = filter_and_prepare_hourly_data(
            scenario_2_raw, graph_config, legend_col, mapping_df
        )

    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, "snapshot")

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    render_chart_layout(
        is_dual_scenario=has_dual_filters and scenario_2_filtered is not None,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        table_1_display_data=scenario_1_filtered,
        table_2_display_data=scenario_2_filtered,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_filtered_bar_hourly,
        render_download_function=render_download_without_data_table,
        **plot_kwargs,
    )

    st.divider()


def render_p9_energy_demand_by_carrier(graph_config: dict) -> None:
    """Render energy demand by carrier type (yearly bar chart).

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
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
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
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
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

    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

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
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
    )

    if has_dual_data:
        _render_scenario_diff_chart(
            scenario_1_grouped, scenario_2_grouped, graph_config, mapping_df
        )

    st.divider()


def render_p10_hourly_demand(graph_config: dict) -> None:
    """Render hourly energy demand (hourly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    scenario_1_filtered, start_date, end_date, is_complete = (
        filter_and_prepare_hourly_data(
            scenario_1_raw, graph_config, legend_col, mapping_df
        )
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = filter_and_prepare_hourly_data(
            scenario_2_raw, graph_config, legend_col, mapping_df
        )

    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, "snapshot")

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    render_chart_layout(
        is_dual_scenario=has_dual_filters and scenario_2_filtered is not None,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        table_1_display_data=scenario_1_filtered,
        table_2_display_data=scenario_2_filtered,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_hourly,
        render_download_function=render_download_without_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
        **plot_kwargs,
    )

    st.divider()


def render_p11_hourly_elec_price(graph_config: dict) -> None:
    """Render hourly electricity price (hourly line chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_line_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    scenario_1_filtered, start_date, end_date, is_complete = (
        filter_and_prepare_hourly_data(
            scenario_1_raw, graph_config, legend_col, mapping_df
        )
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = filter_and_prepare_hourly_data(
            scenario_2_raw, graph_config, legend_col, mapping_df
        )

    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    render_chart_layout(
        is_dual_scenario=has_dual_filters and scenario_2_filtered is not None,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        table_1_display_data=scenario_1_filtered,
        table_2_display_data=scenario_2_filtered,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_line_hourly,
        render_download_function=render_download_without_data_table,
        **plot_kwargs,
    )

    st.divider()


def render_p12_nodal_flow_between_regions(graph_config: dict) -> None:
    """Render nodal flow between regions (hourly line chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_line_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    scenario_1_filtered, start_date, end_date, is_complete = (
        filter_and_prepare_hourly_data(
            scenario_1_raw, graph_config, legend_col, mapping_df
        )
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = filter_and_prepare_hourly_data(
            scenario_2_raw, graph_config, legend_col, mapping_df
        )

    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    render_chart_layout(
        is_dual_scenario=has_dual_filters and scenario_2_filtered is not None,
        scenario_1_vis_display_data=scenario_1_filtered,
        scenario_2_vis_display_data=scenario_2_filtered,
        table_1_display_data=scenario_1_filtered,
        table_2_display_data=scenario_2_filtered,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_line_hourly,
        render_download_function=render_download_without_data_table,
        **plot_kwargs,
    )

    st.divider()


def render_p13_battery_ep_ratio(graph_config: dict) -> None:
    """Render battery energy-to-power ratio (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=graph_config["shared_country"]
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

    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=graph_config["shared_country"]
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

    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

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
            f"{graph_config['download_id'].format(st.session_state.sce1)}"
        ),
    )

    st.divider()


def render_p14_battery_charging_profile(graph_config: dict) -> None:
    """Render battery charging profile with dual y-axes (hourly line chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "line_with_secondary_y_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data (without nice names for p14 special handling)
    monthly_1, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        scenario_1_raw, graph_config
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        monthly_1, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)

    scenario_2_filtered = None
    if has_dual_filters and scenario_2_raw is not None:
        monthly_2, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, graph_config
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            monthly_2, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    # Create label mapping for primary and secondary y-axes
    label_map = {
        label: (
            mapping_df.loc[label, "nice_names"]
            if (mapping_df is not None and label in mapping_df.index)
            else prettify_label(label)
        )
        for label in graph_config["primary_y_lab"] + graph_config["secondary_y_lab"]
    }
    graph_config["label_map"] = label_map

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    if not is_dual or scenario_2_filtered is None:
        st.caption(f"{st.session_state.sce1}")
        colour_mapping = generate_colour_mapping_dict(
            table_name,
            mapping_df,
            add_nice_names(scenario_1_filtered, legend_col, mapping_df),
            legend_col,
        )
        plot_line_with_secondary_y_hourly(
            scenario_1_filtered,
            graph_config,
            colour_mapping,
            y_range,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            **plot_kwargs,
        )
        render_download_without_data_table(
            scenario_1_filtered, graph_config, st.session_state.sce1
        )
    else:
        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1}")
            colour_mapping_1 = generate_colour_mapping_dict(
                table_name,
                mapping_df,
                add_nice_names(scenario_1_filtered, legend_col, mapping_df),
                legend_col,
            )
            plot_line_with_secondary_y_hourly(
                scenario_1_filtered,
                graph_config,
                colour_mapping_1,
                y_range,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
                **plot_kwargs,
            )

        with col3:
            st.caption(f"{st.session_state.sce2}")
            colour_mapping_2 = generate_colour_mapping_dict(
                table_name,
                mapping_df,
                add_nice_names(scenario_2_filtered, legend_col, mapping_df),
                legend_col,
            )
            plot_line_with_secondary_y_hourly(
                scenario_2_filtered,
                graph_config,
                colour_mapping_2,
                y_range,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
                **plot_kwargs,
            )

        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data_table(
                scenario_1_filtered, graph_config, st.session_state.sce1
            )
        with col3:
            render_download_without_data_table(
                scenario_2_filtered, graph_config, st.session_state.sce2
            )

    st.divider()


render_p1_capacity_by_type(config["p1"])
render_p2_capacity_by_region(config["p2"])
render_p3_generation_by_type(config["p3"])
render_p4_share_category(config["p4"])
render_p6_transmission_capacity_between_regions(config["p6"])
render_p7_hourly_generation(config["p7"])
render_p8_regional_hourly_generation(config["p8"])
render_p9_energy_demand_by_carrier(config["p9"])
render_p10_hourly_demand(config["p10"])
render_p11_hourly_elec_price(config["p11"])
render_p12_nodal_flow_between_regions(config["p12"])
render_p13_battery_ep_ratio(config["p13"])
render_p14_battery_charging_profile(config["p14"])
