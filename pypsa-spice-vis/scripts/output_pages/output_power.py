# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Power page under Results section.

Page shows editable power related
dataframes and visualisations from the modelling results.
"""

import os
from typing import Optional

import pandas as pd
import streamlit as st
import yaml
from scripts.data_utils import (
    calculate_min_max_y_scale,
    clean_df_for_plotting,
    filter_dataframe_by_date_range,
    get_filtered_df_and_date_range,
    handle_small_values,
    prettify_label,
    read_result_csv,
)
from scripts.output_st_handler import (
    add_nice_names,
    generate_sidebar,
    get_colour_mapping,
    render_download_with_table,
    render_download_without_data,
    render_section_header,
    setup_country_filter,
    setup_hourly_data_filters,
    setup_radio_button_filter,
    setup_region_filter,
    setup_year_filter,
)
from scripts.plot_settings import create_nice_names_and_color_mapping
from scripts.plot_functions import (
    plot_area_share_yearly,
    plot_bar_with_filter,
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

power_charts = [config[key] for key in POWER_CHART_KEYS if key in config]
table_of_content = [chart["name"] for chart in power_charts]


# ============================================================================
# Helper Functions
# ============================================================================


def _prepare_y_range(
    scenario_1: pd.DataFrame,
    scenario_2: Optional[pd.DataFrame],
    x_col: Optional[str],
) -> dict:
    """Calculate y-axis range for consistent scaling across scenarios."""
    y_range = calculate_min_max_y_scale(scenario_1, scenario_2, x_col)  # type: ignore
    return {"max_scale": y_range["max"], "min_scale": y_range["min"]}


def _normalize_dataframe(df: pd.DataFrame | pd.Series) -> pd.DataFrame:
    """
    Ensure DataFrame is properly formatted (not a Series).

    Handles the case where groupby operations return Series instead of DataFrame.
    """
    return df.reset_index() if isinstance(df, pd.Series) else df


def _render_single_chart_layout(
    scenario_data: pd.DataFrame,
    raw_data: pd.DataFrame,
    config_dict: dict,
    table_name: str,
    legend_col: str,
    mapping_df: pd.DataFrame,
    y_range: dict,
    plot_func,
    render_download_func,
    **plot_kwargs,
):
    """Render single scenario chart with download button."""
    scenario_name = st.session_state.sce1
    st.markdown(f"{scenario_name}")

    colour_mapping = get_colour_mapping(table_name, mapping_df, scenario_data, legend_col)

    plot_key = plot_kwargs.pop("key", f"plotly_chart_{scenario_name}_{table_name}")
    plot_func(scenario_data, config_dict, colour_mapping, y_range, key=plot_key, **plot_kwargs)

    render_download_func(raw_data, config_dict, scenario_name)


def _render_dual_chart_layout(
    scenario_1_data: pd.DataFrame,
    scenario_2_data: pd.DataFrame | None,
    raw_data_1: pd.DataFrame,
    raw_data_2: pd.DataFrame | None,
    config_dict: dict,
    table_name: str,
    legend_col: str,
    mapping_df: pd.DataFrame,
    y_range: dict,
    plot_func,
    render_download_func,
    **plot_kwargs,
):
    """
    Render dual scenario charts side-by-side with download buttons.

    Note: This function should only be called when scenario_2 data is confirmed
    to be non-None via has_dual_filters check.
    """
    # Type guard: function should only be called with valid data
    if scenario_2_data is None or raw_data_2 is None:
        raise ValueError("_render_dual_chart_layout requires non-None scenario 2 data")

    scenario_1_name = st.session_state.sce1
    scenario_2_name = st.session_state.sce2

    # Remove 'key' from plot_kwargs if present (we generate our own keys)
    plot_kwargs.pop("key", None)

    # Charts
    col1, col2 = st.columns([6, 6])
    with col1:
        st.caption(f"{scenario_1_name}")
        colour_mapping_1 = get_colour_mapping(
            table_name, mapping_df, scenario_1_data, legend_col
        )
        plot_key_1 = f"plotly_chart_{scenario_1_name}_{table_name}"
        plot_func(
            scenario_1_data, config_dict, colour_mapping_1, y_range,
            key=plot_key_1, **plot_kwargs
        )

    with col2:
        st.caption(f"{scenario_2_name}")
        colour_mapping_2 = get_colour_mapping(
            table_name, mapping_df, scenario_2_data, legend_col
        )
        plot_key_2 = f"plotly_chart_{scenario_2_name}_{table_name}"
        plot_func(
            scenario_2_data, config_dict, colour_mapping_2, y_range,
            key=plot_key_2, **plot_kwargs
        )

    # Download buttons
    col1, col2 = st.columns([6, 6])
    with col1:
        render_download_func(raw_data_1, config_dict, scenario_1_name)
    with col2:
        render_download_func(raw_data_2, config_dict, scenario_2_name)


def _load_and_validate_hourly_data(
    scenario_name: str, table_name: str, year: str, country: str
) -> Optional[pd.DataFrame]:
    """
    Load hourly data and convert snapshot to datetime.

    Returns None if data is empty or doesn't exist.
    """
    raw_data = read_result_csv(scenario_name, table_name, year=year, country=country)
    if raw_data is None or raw_data.empty:
        return None
    raw_data["snapshot"] = pd.to_datetime(raw_data["snapshot"])
    return raw_data


def _setup_hourly_filters(
    config_dict: dict,
    scenario_1_raw: pd.DataFrame,
    scenario_2_raw: Optional[pd.DataFrame],
    has_dual: bool,
):
    """
    Setup region and date filters for hourly data.

    Updates config_dict in place with filter settings.
    """
    config_dict["shared_region"] = setup_region_filter(
        config_dict, scenario_1_raw, scenario_2_raw, has_dual  # type: ignore
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw, scenario_2_raw, config_dict, has_dual
    )
    config_dict.update(filter_results)


def _filter_and_prepare_hourly_data(
    raw_data: pd.DataFrame | None,
    config_dict: dict,
    legend_col: str,
    mapping_df: pd.DataFrame,
) -> tuple:
    """
    Filter hourly data by date range and prepare for plotting.

    Returns tuple of (filtered_data, start_date, end_date, is_complete).
    Returns (None, None, None, False) if raw_data is None.
    """
    if raw_data is None:
        return None, None, None, False

    monthly_data, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        raw_data, config_dict
    )
    filtered_data = filter_dataframe_by_date_range(
        monthly_data, start_date=start_date, end_date=end_date
    )
    filtered_data = handle_small_values(filtered_data)
    filtered_data = add_nice_names(filtered_data, legend_col, mapping_df)
    return filtered_data, start_date, end_date, is_complete


def render_p1_capacity_by_type(config_p1: dict) -> None:
    """Render installed capacity by technology type (yearly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p1 = {**config_p1, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p1["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p1["shared_country"] = setup_country_filter(
        config_p1, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = config_p1["table_name"]
    legend_col = config_p1["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p1["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = _normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p1["shared_country"]
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = _normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_grouped,
            scenario_1_raw,
            config_p1,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
            key=f"plotly_chart_{config_p1['download_id'].format(st.session_state.sce1)}",
        )
    else:
        _render_dual_chart_layout(
            scenario_1_grouped,
            scenario_2_grouped,
            scenario_1_raw,
            scenario_2_raw,
            config_p1,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
        )

    st.divider()


def render_p2_capacity_by_region(config_p2: dict) -> None:
    """Render installed capacity by region with filter (yearly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p2 = {**config_p2, "graph_type": "bar_with_filter"}
    render_section_header(config_p2["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p2["shared_country"] = setup_country_filter(
        config_p2, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_filter = setup_radio_button_filter(config_p2, is_dual)

    # Extract config values
    table_name = config_p2["table_name"]
    legend_col = config_p2["leg_col"]
    filter_col = config_p2["fil_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p2["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    filter_value_s1 = shared_filter or st.radio(
        f"{config_p2['slider_id'].format(st.session_state.sce1)} "
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
            st.session_state.sce2, table_name, country=config_p2["shared_country"]
        )
        if scenario_2_raw is not None:
            filter_value_s2 = shared_filter or st.radio(
                f"{config_p2['slider_id'].format(st.session_state.sce2)} "
                f"Select {filter_col} ({st.session_state.sce2}):",
                options=scenario_2_raw[filter_col].unique(),
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
            scenario_2_filtered = scenario_2_raw[
                scenario_2_raw[filter_col] == filter_value_s2
            ].copy()
            scenario_2_filtered = add_nice_names(scenario_2_filtered, legend_col, mapping_df)
            scenario_2_filtered = clean_df_for_plotting(legend_col, scenario_2_filtered)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_filtered is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_raw,
            config_p2,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_bar_with_filter,
            render_download_with_table,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_raw,
            scenario_2_raw,
            config_p2,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_bar_with_filter,
            render_download_with_table,
        )

    st.divider()


def render_p3_generation_by_type(config_p3: dict) -> None:
    """Render electricity generation by technology type (yearly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p3 = {**config_p3, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p3["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p3["shared_country"] = setup_country_filter(
        config_p3, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = config_p3["table_name"]
    legend_col = config_p3["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p3["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = _normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p3["shared_country"]
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = _normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_grouped,
            scenario_1_raw,
            config_p3,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
            key=f"plotly_chart_{config_p3['download_id'].format(st.session_state.sce1)}",
        )
    else:
        _render_dual_chart_layout(
            scenario_1_grouped,
            scenario_2_grouped,
            scenario_1_raw,
            scenario_2_raw,
            config_p3,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
        )

    st.divider()


def render_p4_share_category(config_p4: dict) -> None:
    """Render generation share by category (yearly area chart)."""
    # Inject graph_type for legacy compatibility
    config_p4 = {**config_p4, "graph_type": "area_share_yearly"}
    render_section_header(config_p4["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p4["shared_country"] = setup_country_filter(
        config_p4, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = config_p4["table_name"]
    legend_col = config_p4["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p4["shared_country"]
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
            st.session_state.sce2, table_name, country=config_p4["shared_country"]
        )
        if scenario_2_raw is not None:
            scenario_2_plot = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_plot = add_nice_names(scenario_2_plot, legend_col, mapping_df)

    # Render charts (area charts don't use y_range, only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_plot is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_plot,
            scenario_1_raw,
            config_p4,
            table_name,
            legend_col,
            mapping_df,
            {},  # No y_range for area charts
            plot_area_share_yearly,
            render_download_with_table,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_plot,
            scenario_2_plot,
            scenario_1_raw,
            scenario_2_raw,
            config_p4,
            table_name,
            legend_col,
            mapping_df,
            {},  # No y_range for area charts
            plot_area_share_yearly,
            render_download_with_table,
        )

    st.divider()


def render_p6_transmission_capacity_between_regions(config_p6: dict) -> None:
    """Render transmission capacity between regions with filter (yearly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p6 = {**config_p6, "graph_type": "bar_with_filter"}
    render_section_header(config_p6["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p6["shared_country"] = setup_country_filter(
        config_p6, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_filter = setup_radio_button_filter(config_p6, is_dual)

    # Extract config values
    table_name = config_p6["table_name"]
    legend_col = config_p6["leg_col"]
    filter_col = config_p6["fil_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p6["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    filter_value_s1 = shared_filter or st.radio(
        f"{config_p6['slider_id'].format(st.session_state.sce1)} "
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
            st.session_state.sce2, table_name, country=config_p6["shared_country"]
        )
        if scenario_2_raw is not None:
            filter_value_s2 = shared_filter or st.radio(
                f"{config_p6['slider_id'].format(st.session_state.sce2)} "
                f"Select {filter_col} ({st.session_state.sce2}):",
                options=scenario_2_raw[filter_col].unique(),
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
            scenario_2_filtered = scenario_2_raw[
                scenario_2_raw[filter_col] == filter_value_s2
            ].copy()
            scenario_2_filtered = add_nice_names(scenario_2_filtered, legend_col, mapping_df)
            scenario_2_filtered = clean_df_for_plotting(legend_col, scenario_2_filtered)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_filtered is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_raw,
            config_p6,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_bar_with_filter,
            render_download_with_table,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_raw,
            scenario_2_raw,
            config_p6,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_bar_with_filter,
            render_download_with_table,
        )

    st.divider()


def render_p7_hourly_generation(config_p7: dict) -> None:
    """Render hourly electricity generation (hourly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p7 = {**config_p7, "graph_type": "simple_bar_hourly"}
    render_section_header(config_p7["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p7["shared_country"] = setup_country_filter(
        config_p7, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p7, is_dual)
    config_p7["shared_year"] = str(shared_year)

    # Extract config values
    table_name = config_p7["table_name"]
    legend_col = config_p7["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario 1 data
    scenario_1_raw = _load_and_validate_hourly_data(
        st.session_state.sce1, table_name, str(shared_year), config_p7["shared_country"]
    )
    if scenario_1_raw is None:
        st.divider()
        return

    # Load and validate scenario 2 data (if dual mode)
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = _load_and_validate_hourly_data(
            st.session_state.sce2, table_name, str(shared_year), config_p7["shared_country"]
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    _setup_hourly_filters(config_p7, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare scenario 1 data
    scenario_1_filtered, start_date, end_date, is_complete = _filter_and_prepare_hourly_data(
        scenario_1_raw, config_p7, legend_col, mapping_df
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    # Filter and prepare scenario 2 data (if dual mode)
    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = _filter_and_prepare_hourly_data(
            scenario_2_raw, config_p7, legend_col, mapping_df
        )

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, "snapshot")

    # Render charts
    plot_kwargs = {"start_date": start_date, "end_date": end_date, "is_complete": is_complete}
    if not is_dual:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_filtered,  # Use filtered data for download in hourly charts
            config_p7,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_hourly,
            render_download_without_data,
            key=f"plotly_chart_{config_p7['download_id'].format(st.session_state.sce1)}",
            **plot_kwargs,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_filtered,  # Use filtered data for download in hourly charts
            scenario_2_filtered,
            config_p7,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_hourly,
            render_download_without_data,
            **plot_kwargs,
        )

    st.divider()


def render_p8_regional_hourly_generation(config_p8: dict) -> None:
    """Render regional hourly generation with filter (hourly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p8 = {**config_p8, "graph_type": "filtered_bar_hourly"}
    render_section_header(config_p8["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p8["shared_country"] = setup_country_filter(
        config_p8, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p8, is_dual)
    config_p8["shared_year"] = str(shared_year)

    # Extract config values
    table_name = config_p8["table_name"]
    legend_col = config_p8["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = _load_and_validate_hourly_data(
        st.session_state.sce1, table_name, str(shared_year), config_p8["shared_country"]
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = _load_and_validate_hourly_data(
            st.session_state.sce2, table_name, str(shared_year), config_p8["shared_country"]
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    _setup_hourly_filters(config_p8, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data
    scenario_1_filtered, start_date, end_date, is_complete = _filter_and_prepare_hourly_data(
        scenario_1_raw, config_p8, legend_col, mapping_df
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = _filter_and_prepare_hourly_data(
            scenario_2_raw, config_p8, legend_col, mapping_df
        )

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, "snapshot")

    # Render charts
    plot_kwargs = {"start_date": start_date, "end_date": end_date, "is_complete": is_complete}
    if not is_dual:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_filtered,
            config_p8,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_filtered_bar_hourly,
            render_download_without_data,
            **plot_kwargs,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_filtered,
            scenario_2_filtered,
            config_p8,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_filtered_bar_hourly,
            render_download_without_data,
            **plot_kwargs,
        )

    st.divider()


def render_p9_energy_demand_by_carrier(config_p9: dict) -> None:
    """Render energy demand by carrier type (yearly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p9 = {**config_p9, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p9["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p9["shared_country"] = setup_country_filter(
        config_p9, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = config_p9["table_name"]
    legend_col = config_p9["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p9["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = _normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p9["shared_country"]
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = _normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_grouped,
            scenario_1_raw,
            config_p9,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
            key=f"plotly_chart_{config_p9['download_id'].format(st.session_state.sce1)}",
        )
    else:
        _render_dual_chart_layout(
            scenario_1_grouped,
            scenario_2_grouped,
            scenario_1_raw,
            scenario_2_raw,
            config_p9,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
        )

    st.divider()


def render_p10_hourly_demand(config_p10: dict) -> None:
    """Render hourly energy demand (hourly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p10 = {**config_p10, "graph_type": "simple_bar_hourly"}
    render_section_header(config_p10["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p10["shared_country"] = setup_country_filter(
        config_p10, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p10, is_dual)
    config_p10["shared_year"] = str(shared_year)

    # Extract config values
    table_name = config_p10["table_name"]
    legend_col = config_p10["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = _load_and_validate_hourly_data(
        st.session_state.sce1, table_name, str(shared_year), config_p10["shared_country"]
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = _load_and_validate_hourly_data(
            st.session_state.sce2, table_name, str(shared_year), config_p10["shared_country"]
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    _setup_hourly_filters(config_p10, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data
    scenario_1_filtered, start_date, end_date, is_complete = _filter_and_prepare_hourly_data(
        scenario_1_raw, config_p10, legend_col, mapping_df
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = _filter_and_prepare_hourly_data(
            scenario_2_raw, config_p10, legend_col, mapping_df
        )

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, "snapshot")

    # Render charts
    plot_kwargs = {"start_date": start_date, "end_date": end_date, "is_complete": is_complete}
    if not is_dual:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_filtered,
            config_p10,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_hourly,
            render_download_without_data,
            key=f"plotly_chart_{config_p10['download_id'].format(st.session_state.sce1)}",
            **plot_kwargs,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_filtered,
            scenario_2_filtered,
            config_p10,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_hourly,
            render_download_without_data,
            **plot_kwargs,
        )

    st.divider()


def render_p11_hourly_elec_price(config_p11: dict) -> None:
    """Render hourly electricity price (hourly line chart)."""
    # Inject graph_type for legacy compatibility
    config_p11 = {**config_p11, "graph_type": "simple_line_hourly"}
    render_section_header(config_p11["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p11["shared_country"] = setup_country_filter(
        config_p11, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p11, is_dual)
    config_p11["shared_year"] = str(shared_year)

    # Extract config values
    table_name = config_p11["table_name"]
    legend_col = config_p11["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = _load_and_validate_hourly_data(
        st.session_state.sce1, table_name, str(shared_year), config_p11["shared_country"]
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = _load_and_validate_hourly_data(
            st.session_state.sce2, table_name, str(shared_year), config_p11["shared_country"]
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    _setup_hourly_filters(config_p11, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data
    scenario_1_filtered, start_date, end_date, is_complete = _filter_and_prepare_hourly_data(
        scenario_1_raw, config_p11, legend_col, mapping_df
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = _filter_and_prepare_hourly_data(
            scenario_2_raw, config_p11, legend_col, mapping_df
        )

    # Calculate common y-axis range (None for x_col means no grouping)
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    # Render charts
    plot_kwargs = {"start_date": start_date, "end_date": end_date, "is_complete": is_complete}
    if not is_dual:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_filtered,
            config_p11,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_line_hourly,
            render_download_without_data,
            **plot_kwargs,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_filtered,
            scenario_2_filtered,
            config_p11,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_line_hourly,
            render_download_without_data,
            **plot_kwargs,
        )

    st.divider()


def render_p12_nodal_flow_between_regions(config_p12: dict) -> None:
    """Render nodal flow between regions (hourly line chart)."""
    # Inject graph_type for legacy compatibility
    config_p12 = {**config_p12, "graph_type": "simple_line_hourly"}
    render_section_header(config_p12["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p12["shared_country"] = setup_country_filter(
        config_p12, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p12, is_dual)
    config_p12["shared_year"] = str(shared_year)

    # Extract config values
    table_name = config_p12["table_name"]
    legend_col = config_p12["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = _load_and_validate_hourly_data(
        st.session_state.sce1, table_name, str(shared_year), config_p12["shared_country"]
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = _load_and_validate_hourly_data(
            st.session_state.sce2, table_name, str(shared_year), config_p12["shared_country"]
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    _setup_hourly_filters(config_p12, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data
    scenario_1_filtered, start_date, end_date, is_complete = _filter_and_prepare_hourly_data(
        scenario_1_raw, config_p12, legend_col, mapping_df
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_filtered, _, _, _ = _filter_and_prepare_hourly_data(
            scenario_2_raw, config_p12, legend_col, mapping_df
        )

    # Calculate common y-axis range (None for x_col means no grouping)
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    # Render charts
    plot_kwargs = {"start_date": start_date, "end_date": end_date, "is_complete": is_complete}
    if not is_dual:
        _render_single_chart_layout(
            scenario_1_filtered,
            scenario_1_filtered,
            config_p12,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_line_hourly,
            render_download_without_data,
            **plot_kwargs,
        )
    else:
        _render_dual_chart_layout(
            scenario_1_filtered,
            scenario_2_filtered,
            scenario_1_filtered,
            scenario_2_filtered,
            config_p12,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_line_hourly,
            render_download_without_data,
            **plot_kwargs,
        )

    st.divider()


def render_p13_battery_ep_ratio(config_p13: dict) -> None:
    """Render battery energy-to-power ratio (yearly bar chart)."""
    # Inject graph_type for legacy compatibility
    config_p13 = {**config_p13, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p13["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p13["shared_country"] = setup_country_filter(
        config_p13, is_dual, scenario_tag=st.session_state.sce1
    )

    # Extract config values
    table_name = config_p13["table_name"]
    legend_col = config_p13["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p13["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = _normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p13["shared_country"]
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = _normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    if not has_dual_data:
        _render_single_chart_layout(
            scenario_1_grouped,
            scenario_1_raw,
            config_p13,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
            key=f"plotly_chart_{config_p13['download_id'].format(st.session_state.sce1)}",
        )
    else:
        _render_dual_chart_layout(
            scenario_1_grouped,
            scenario_2_grouped,
            scenario_1_raw,
            scenario_2_raw,
            config_p13,
            table_name,
            legend_col,
            mapping_df,
            y_range,
            plot_simple_bar_yearly,
            render_download_with_table,
        )

    st.divider()


def render_p14_battery_charging_profile(config_p14: dict) -> None:
    """Render battery charging profile with dual y-axes (hourly line chart)."""
    # Inject graph_type for legacy compatibility
    config_p14 = {**config_p14, "graph_type": "line_with_secondary_y_hourly"}
    render_section_header(config_p14["name"])

    # Setup filters
    is_dual = bool(st.session_state.sce2 and st.session_state.sce2 != "")
    config_p14["shared_country"] = setup_country_filter(
        config_p14, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p14, is_dual)
    config_p14["shared_year"] = str(shared_year)

    # Extract config values
    table_name = config_p14["table_name"]
    legend_col = config_p14["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = _load_and_validate_hourly_data(
        st.session_state.sce1, table_name, str(shared_year), config_p14["shared_country"]
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = _load_and_validate_hourly_data(
            st.session_state.sce2, table_name, str(shared_year), config_p14["shared_country"]
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    _setup_hourly_filters(config_p14, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data (without nice names for p14 special handling)
    monthly_1, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        scenario_1_raw, config_p14
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
        monthly_2, _, _, _ = get_filtered_df_and_date_range(scenario_2_raw, config_p14)
        scenario_2_filtered = filter_dataframe_by_date_range(
            monthly_2, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)

    # Calculate common y-axis range
    y_range = _prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    # Create label mapping for primary and secondary y-axes
    label_map = {
        label: (
            mapping_df.loc[label, "nice_names"]
            if (mapping_df is not None and label in mapping_df.index)
            else prettify_label(label)
        )
        for label in config_p14["primary_y_lab"] + config_p14["secondary_y_lab"]
    }
    config_p14["label_map"] = label_map

    # Render charts (special handling for dual-axis plot)
    plot_kwargs = {"start_date": start_date, "end_date": end_date, "is_complete": is_complete}

    if not is_dual or scenario_2_filtered is None:
        st.caption(f"{st.session_state.sce1}")
        colour_mapping = get_colour_mapping(
            table_name,
            mapping_df,
            add_nice_names(scenario_1_filtered, legend_col, mapping_df),
            legend_col,
        )
        plot_line_with_secondary_y_hourly(
            scenario_1_filtered,
            config_p14,
            colour_mapping,
            y_range,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            **plot_kwargs,
        )
        render_download_without_data(scenario_1_filtered, config_p14, st.session_state.sce1)
    else:
        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1}")
            colour_mapping_1 = get_colour_mapping(
                table_name,
                mapping_df,
                add_nice_names(scenario_1_filtered, legend_col, mapping_df),
                legend_col,
            )
            plot_line_with_secondary_y_hourly(
                scenario_1_filtered,
                config_p14,
                colour_mapping_1,
                y_range,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
                **plot_kwargs,
            )

        with col3:
            st.caption(f"{st.session_state.sce2}")
            colour_mapping_2 = get_colour_mapping(
                table_name,
                mapping_df,
                add_nice_names(scenario_2_filtered, legend_col, mapping_df),
                legend_col,
            )
            plot_line_with_secondary_y_hourly(
                scenario_2_filtered,
                config_p14,
                colour_mapping_2,
                y_range,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
                **plot_kwargs,
            )

        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(scenario_1_filtered, config_p14, st.session_state.sce1)
        with col3:
            render_download_without_data(scenario_2_filtered, config_p14, st.session_state.sce2)

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

generate_sidebar(table_of_content)
