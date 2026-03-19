# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Industry page under Results section.

Page shows editable power related
dataframes and visualisations from the modelling results.
"""

import streamlit as st

from scripts.data_utils import (
    add_nice_names,
    clean_df_for_plotting,
    normalize_dataframe,
    prepare_y_range,
    prettify_label,
    read_result_csv,
)
from scripts.output_st_handler import (
    generate_sidebar,
    render_chart_layout,
    render_download_with_data_table,
    render_section_header,
    setup_country_filter,
    setup_radio_button_filter,
)
from scripts.plot_functions import (
    create_nice_names_and_color_mapping,
    plot_bar_with_filter,
    plot_simple_bar_yearly,
)

# =========================== Render functions for each chart =========================


def render_i1_capacity_by_carrier(graph_config: dict) -> None:
    """Render installed capacity by carrier (yearly bar chart).

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

    st.divider()


def render_i2_capacity_by_region(graph_config: dict) -> None:
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
    # Aggregate low and heat heat
    scenario_1_filtered = (
        scenario_1_filtered.groupby(["region", "legend", "year"])
        .sum()["value"]
        .reset_index()
    )
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
            # Aggregate low and heat heat
            scenario_2_filtered = (
                scenario_2_filtered.groupby(["region", "legend", "year"])
                .sum()["value"]
                .reset_index()
            )
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

    st.divider()


def render_i3_generation_by_region(graph_config: dict) -> None:
    """Render generation by region with filter (yearly bar chart).

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
    # Aggregate low and heat heat
    scenario_1_filtered = (
        scenario_1_filtered.groupby(["region", "legend", "year"])
        .sum()["value"]
        .reset_index()
    )
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
            # Aggregate low and heat heat
            scenario_2_filtered = (
                scenario_2_filtered.groupby(["region", "legend", "year"])
                .sum()["value"]
                .reset_index()
            )
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

    st.divider()


def render_i4_generation_by_type(graph_config: dict) -> None:
    """Render generation by technology type (yearly bar chart).

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

    st.divider()


if __name__ == "__main__":
    st.title(":material/construction: Industry")
    output_config = st.session_state.output_config["industry"]
    INDUSTRY_CHART_KEYS = [
        "i1",
        "i2",
        "i3",
        "i4",
    ]
    show_industry = "i" in str(st.session_state.get("sector", "")).lower()
    if show_industry:
        render_i1_capacity_by_carrier(output_config["i1"])
        render_i2_capacity_by_region(output_config["i2"])
        render_i3_generation_by_region(output_config["i3"])
        render_i4_generation_by_type(output_config["i4"])

        # Only include charts in the sidebar if they are present in the config
        industry_charts = [
            output_config[key] for key in INDUSTRY_CHART_KEYS if key in output_config
        ]
        table_of_content = [chart["name"] for chart in industry_charts]
        generate_sidebar(table_of_content)
    else:
        st.warning("Industry sector is not selected. No chart is displayed.")
