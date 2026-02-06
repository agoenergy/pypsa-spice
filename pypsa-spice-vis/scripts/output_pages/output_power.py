# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Power page under Results section.

Page shows editable power related
dataframes and visualisations from the modelling results.
"""

import os

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
from scripts.power_plotting import (
    plot_area_share_yearly,
    plot_bar_with_filter,
    plot_filtered_bar_hourly,
    plot_line_with_secondary_y_hourly,
    plot_simple_bar_hourly,
    plot_simple_bar_yearly,
    plot_simple_line_hourly,
    plot_simple_line_yearly,
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


def render_p1_capacity_by_type(config_p1: dict) -> None:
    config_p1 = {**config_p1, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p1["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p1["shared_country"] = setup_country_filter(
        config_p1, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = config_p1["table_name"]
    legend_col = config_p1["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p1["shared_country"]
    )
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p1["shared_country"]
        )
        scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
        scenario_2_grouped = scenario_2_grouped.groupby(
            ["year", legend_col], as_index=False
        )["value"].sum()
        scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    y_range = calculate_min_max_y_scale(scenario_1_grouped, scenario_2_grouped, "year")
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.markdown(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_grouped, legend_col
        )
        plot_simple_bar_yearly(
            scenario_1_grouped,
            config_p1,
            colour_mapping,
            y_range,
            key=(
                f"plotly_chart_{config_p1['download_id'].format(st.session_state.sce1)}"
            ),
        )
        render_download_with_table(scenario_1_raw, config_p1, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_1_grouped,
                config_p1,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p1['download_id'].format(st.session_state.sce1)}"
                ),
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_2_grouped,
                config_p1,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p1['download_id'].format(st.session_state.sce2)}"
                ),
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(scenario_1_raw, config_p1, st.session_state.sce1)
        with col3:
            render_download_with_table(scenario_2_raw, config_p1, st.session_state.sce2)

    st.divider()


def render_p2_capacity_by_region(config_p2: dict) -> None:
    config_p2 = {**config_p2, "graph_type": "bar_with_filter"}
    render_section_header(config_p2["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p2["shared_country"] = setup_country_filter(
        config_p2, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_filter = setup_radio_button_filter(config_p2, is_dual)

    table_name = config_p2["table_name"]
    legend_col = config_p2["leg_col"]
    filter_col = config_p2["fil_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p2["shared_country"]
    )
    filter_value_s1 = (
        shared_filter
        if shared_filter
        else st.radio(
            (
                f"{config_p2['slider_id'].format(st.session_state.sce1)} "
                f"Select {filter_col} ({st.session_state.sce1}):"
            ),
            options=scenario_1_raw[filter_col].unique(),
            format_func=prettify_label,
            horizontal=True,
            label_visibility="collapsed",
        )
    )
    scenario_1_filtered = scenario_1_raw.loc[
        scenario_1_raw[filter_col] == filter_value_s1
    ].copy()
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)
    scenario_1_filtered = clean_df_for_plotting(legend_col, scenario_1_filtered)

    scenario_2_filtered = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p2["shared_country"]
        )
        filter_value_s2 = shared_filter
        if not filter_value_s2:
            filter_value_s2 = st.radio(
                (
                    f"{config_p2['slider_id'].format(st.session_state.sce2)} "
                    f"Select {filter_col} ({st.session_state.sce2}):"
                ),
                options=scenario_2_raw[filter_col].unique(),
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
        scenario_2_filtered = scenario_2_raw.loc[
            scenario_2_raw[filter_col] == filter_value_s2
        ].copy()
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )
        scenario_2_filtered = clean_df_for_plotting(legend_col, scenario_2_filtered)

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered, scenario_2_filtered, "year"
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1}")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_bar_with_filter(
            scenario_1_filtered,
            config_p2,
            colour_mapping,
            y_range,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_with_table(scenario_1_raw, config_p2, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_bar_with_filter(
                scenario_1_filtered,
                config_p2,
                colour_mapping,
                y_range,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_bar_with_filter(
                scenario_2_filtered,
                config_p2,
                colour_mapping,
                y_range,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(scenario_1_raw, config_p2, st.session_state.sce1)
        with col3:
            render_download_with_table(scenario_2_raw, config_p2, st.session_state.sce2)

    st.divider()


def render_p3_generation_by_type(config_p3: dict) -> None:
    config_p3 = {**config_p3, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p3["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p3["shared_country"] = setup_country_filter(
        config_p3, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = config_p3["table_name"]
    legend_col = config_p3["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p3["shared_country"]
    )
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p3["shared_country"]
        )
        scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
        scenario_2_grouped = scenario_2_grouped.groupby(
            ["year", legend_col], as_index=False
        )["value"].sum()
        scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    y_range = calculate_min_max_y_scale(scenario_1_grouped, scenario_2_grouped, "year")
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_grouped, legend_col
        )
        plot_simple_bar_yearly(
            scenario_1_grouped,
            config_p3,
            colour_mapping,
            y_range,
            key=(
                f"plotly_chart_{config_p3['download_id'].format(st.session_state.sce1)}"
            ),
        )
        render_download_with_table(scenario_1_raw, config_p3, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_1_grouped,
                config_p3,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p3['download_id'].format(st.session_state.sce1)}"
                ),
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_2_grouped,
                config_p3,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p3['download_id'].format(st.session_state.sce2)}"
                ),
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(scenario_1_raw, config_p3, st.session_state.sce1)
        with col3:
            render_download_with_table(scenario_2_raw, config_p3, st.session_state.sce2)

    st.divider()


def render_p4_share_category(config_p4: dict) -> None:
    config_p4 = {**config_p4, "graph_type": "area_share_yearly"}
    render_section_header(config_p4["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p4["shared_country"] = setup_country_filter(
        config_p4, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = config_p4["table_name"]
    legend_col = config_p4["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p4["shared_country"]
    )
    scenario_1_plot = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_plot = add_nice_names(scenario_1_plot, legend_col, mapping_df)

    scenario_2_plot = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p4["shared_country"]
        )
        scenario_2_plot = clean_df_for_plotting(legend_col, scenario_2_raw)
        scenario_2_plot = add_nice_names(scenario_2_plot, legend_col, mapping_df)

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_plot, legend_col
        )
        plot_area_share_yearly(
            scenario_1_plot,
            config_p4,
            colour_mapping,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_with_table(scenario_1_raw, config_p4, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_plot, legend_col
            )
            plot_area_share_yearly(
                scenario_1_plot,
                config_p4,
                colour_mapping,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_plot, legend_col
            )
            plot_area_share_yearly(
                scenario_2_plot,
                config_p4,
                colour_mapping,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(scenario_1_raw, config_p4, st.session_state.sce1)
        with col3:
            render_download_with_table(scenario_2_raw, config_p4, st.session_state.sce2)

    st.divider()


def render_p6_transmission_capacity_between_regions(config_p6: dict) -> None:
    config_p6 = {**config_p6, "graph_type": "bar_with_filter"}
    render_section_header(config_p6["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p6["shared_country"] = setup_country_filter(
        config_p6, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_filter = setup_radio_button_filter(config_p6, is_dual)

    table_name = config_p6["table_name"]
    legend_col = config_p6["leg_col"]
    filter_col = config_p6["fil_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p6["shared_country"]
    )
    filter_value_s1 = (
        shared_filter
        if shared_filter
        else st.radio(
            (
                f"{config_p6['slider_id'].format(st.session_state.sce1)} "
                f"Select {filter_col} ({st.session_state.sce1}):"
            ),
            options=scenario_1_raw[filter_col].unique(),
            format_func=prettify_label,
            horizontal=True,
            label_visibility="collapsed",
        )
    )
    scenario_1_filtered = scenario_1_raw.loc[
        scenario_1_raw[filter_col] == filter_value_s1
    ].copy()
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)
    scenario_1_filtered = clean_df_for_plotting(legend_col, scenario_1_filtered)

    scenario_2_filtered = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p6["shared_country"]
        )
        filter_value_s2 = shared_filter
        if not filter_value_s2:
            filter_value_s2 = st.radio(
                (
                    f"{config_p6['slider_id'].format(st.session_state.sce2)} "
                    f"Select {filter_col} ({st.session_state.sce2}):"
                ),
                options=scenario_2_raw[filter_col].unique(),
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
        scenario_2_filtered = scenario_2_raw.loc[
            scenario_2_raw[filter_col] == filter_value_s2
        ].copy()
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )
        scenario_2_filtered = clean_df_for_plotting(legend_col, scenario_2_filtered)

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered, scenario_2_filtered, "year"
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_bar_with_filter(
            scenario_1_filtered,
            config_p6,
            colour_mapping,
            y_range,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_with_table(scenario_1_raw, config_p6, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_bar_with_filter(
                scenario_1_filtered,
                config_p6,
                colour_mapping,
                y_range,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_bar_with_filter(
                scenario_2_filtered,
                config_p6,
                colour_mapping,
                y_range,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(scenario_1_raw, config_p6, st.session_state.sce1)
        with col3:
            render_download_with_table(scenario_2_raw, config_p6, st.session_state.sce2)

    st.divider()


def render_p7_hourly_generation(config_p7: dict) -> None:
    config_p7 = {**config_p7, "graph_type": "simple_bar_hourly"}
    render_section_header(config_p7["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p7["shared_country"] = setup_country_filter(
        config_p7, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p7, is_dual)
    config_p7["shared_year"] = str(shared_year)

    table_name = config_p7["table_name"]
    legend_col = config_p7["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1,
        table_name,
        year=str(shared_year),
        country=config_p7["shared_country"],
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return

    scenario_1_raw["snapshot"] = pd.to_datetime(scenario_1_raw["snapshot"])

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2,
            table_name,
            year=str(shared_year),
            country=config_p7["shared_country"],
        )
        if scenario_2_raw is not None and not scenario_2_raw.empty:
            scenario_2_raw["snapshot"] = pd.to_datetime(scenario_2_raw["snapshot"])

    has_dual_filters = (
        is_dual and scenario_2_raw is not None and not scenario_2_raw.empty
    )
    config_p7["shared_region"] = setup_region_filter(
        config_p7,
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        has_dual_filters,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        config_p7,
        has_dual_filters,
    )
    config_p7.update(filter_results)

    scenario_1_month, start_date, end_date, is_complete = (
        get_filtered_df_and_date_range(scenario_1_raw, config_p7)
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        scenario_1_month, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_month, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, config_p7
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            scenario_2_month, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered,
        scenario_2_filtered if has_dual_filters else None,
        "snapshot",
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_simple_bar_hourly(
            scenario_1_filtered,
            config_p7,
            colour_mapping,
            y_range,
            start_date,
            end_date,
            is_complete,
            key=(
                f"plotly_chart_{config_p7['download_id'].format(st.session_state.sce1)}"
            ),
        )
        render_download_without_data(
            scenario_1_filtered, config_p7, st.session_state.sce1
        )
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_simple_bar_hourly(
                scenario_1_filtered,
                config_p7,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=(
                    "plotly_chart_"
                    f"{config_p7['download_id'].format(st.session_state.sce1)}"
                ),
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_simple_bar_hourly(
                scenario_2_filtered,
                config_p7,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=(
                    "plotly_chart_"
                    f"{config_p7['download_id'].format(st.session_state.sce2)}"
                ),
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(
                scenario_1_filtered, config_p7, st.session_state.sce1
            )
        with col3:
            render_download_without_data(
                scenario_2_filtered, config_p7, st.session_state.sce2
            )

    st.divider()


def render_p8_regional_hourly_generation(config_p8: dict) -> None:
    config_p8 = {**config_p8, "graph_type": "filtered_bar_hourly"}
    render_section_header(config_p8["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p8["shared_country"] = setup_country_filter(
        config_p8, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p8, is_dual)
    config_p8["shared_year"] = str(shared_year)

    table_name = config_p8["table_name"]
    legend_col = config_p8["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1,
        table_name,
        year=str(shared_year),
        country=config_p8["shared_country"],
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return

    scenario_1_raw["snapshot"] = pd.to_datetime(scenario_1_raw["snapshot"])

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2,
            table_name,
            year=str(shared_year),
            country=config_p8["shared_country"],
        )
        if scenario_2_raw is not None and not scenario_2_raw.empty:
            scenario_2_raw["snapshot"] = pd.to_datetime(scenario_2_raw["snapshot"])

    has_dual_filters = (
        is_dual and scenario_2_raw is not None and not scenario_2_raw.empty
    )
    config_p8["shared_region"] = setup_region_filter(
        config_p8,
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        has_dual_filters,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        config_p8,
        has_dual_filters,
    )
    config_p8.update(filter_results)

    scenario_1_month, start_date, end_date, is_complete = (
        get_filtered_df_and_date_range(scenario_1_raw, config_p8)
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        scenario_1_month, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_month, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, config_p8
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            scenario_2_month, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered,
        scenario_2_filtered if has_dual_filters else None,
        "snapshot",
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_filtered_bar_hourly(
            scenario_1_filtered,
            config_p8,
            colour_mapping,
            y_range,
            start_date,
            end_date,
            is_complete,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_without_data(
            scenario_1_filtered, config_p8, st.session_state.sce1
        )
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_filtered_bar_hourly(
                scenario_1_filtered,
                config_p8,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_filtered_bar_hourly(
                scenario_2_filtered,
                config_p8,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(
                scenario_1_filtered, config_p8, st.session_state.sce1
            )
        with col3:
            render_download_without_data(
                scenario_2_filtered, config_p8, st.session_state.sce2
            )

    st.divider()


def render_p9_energy_demand_by_carrier(config_p9: dict) -> None:
    config_p9 = {**config_p9, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p9["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p9["shared_country"] = setup_country_filter(
        config_p9, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = config_p9["table_name"]
    legend_col = config_p9["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p9["shared_country"]
    )
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p9["shared_country"]
        )
        scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
        scenario_2_grouped = scenario_2_grouped.groupby(
            ["year", legend_col], as_index=False
        )["value"].sum()
        scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    y_range = calculate_min_max_y_scale(scenario_1_grouped, scenario_2_grouped, "year")
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_grouped, legend_col
        )
        plot_simple_bar_yearly(
            scenario_1_grouped,
            config_p9,
            colour_mapping,
            y_range,
            key=(
                f"plotly_chart_{config_p9['download_id'].format(st.session_state.sce1)}"
            ),
        )
        render_download_with_table(scenario_1_raw, config_p9, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_1_grouped,
                config_p9,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p9['download_id'].format(st.session_state.sce1)}"
                ),
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_2_grouped,
                config_p9,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p9['download_id'].format(st.session_state.sce2)}"
                ),
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(scenario_1_raw, config_p9, st.session_state.sce1)
        with col3:
            render_download_with_table(scenario_2_raw, config_p9, st.session_state.sce2)

    st.divider()


def render_p10_hourly_demand(config_p10: dict) -> None:
    config_p10 = {**config_p10, "graph_type": "simple_bar_hourly"}
    render_section_header(config_p10["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p10["shared_country"] = setup_country_filter(
        config_p10, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p10, is_dual)
    config_p10["shared_year"] = str(shared_year)

    table_name = config_p10["table_name"]
    legend_col = config_p10["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1,
        table_name,
        year=str(shared_year),
        country=config_p10["shared_country"],
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return

    scenario_1_raw["snapshot"] = pd.to_datetime(scenario_1_raw["snapshot"])

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2,
            table_name,
            year=str(shared_year),
            country=config_p10["shared_country"],
        )
        if scenario_2_raw is not None and not scenario_2_raw.empty:
            scenario_2_raw["snapshot"] = pd.to_datetime(scenario_2_raw["snapshot"])

    has_dual_filters = (
        is_dual and scenario_2_raw is not None and not scenario_2_raw.empty
    )
    config_p10["shared_region"] = setup_region_filter(
        config_p10,
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        has_dual_filters,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        config_p10,
        has_dual_filters,
    )
    config_p10.update(filter_results)

    scenario_1_month, start_date, end_date, is_complete = (
        get_filtered_df_and_date_range(scenario_1_raw, config_p10)
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        scenario_1_month, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_month, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, config_p10
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            scenario_2_month, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered,
        scenario_2_filtered if has_dual_filters else None,
        "snapshot",
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_simple_bar_hourly(
            scenario_1_filtered,
            config_p10,
            colour_mapping,
            y_range,
            start_date,
            end_date,
            is_complete,
            key=(
                "plotly_chart_"
                f"{config_p10['download_id'].format(st.session_state.sce1)}"
            ),
        )
        render_download_without_data(
            scenario_1_filtered, config_p10, st.session_state.sce1
        )
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_simple_bar_hourly(
                scenario_1_filtered,
                config_p10,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=(
                    "plotly_chart_"
                    f"{config_p10['download_id'].format(st.session_state.sce1)}"
                ),
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_simple_bar_hourly(
                scenario_2_filtered,
                config_p10,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=(
                    "plotly_chart_"
                    f"{config_p10['download_id'].format(st.session_state.sce2)}"
                ),
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(
                scenario_1_filtered, config_p10, st.session_state.sce1
            )
        with col3:
            render_download_without_data(
                scenario_2_filtered, config_p10, st.session_state.sce2
            )

    st.divider()


def render_p11_hourly_elec_price(config_p11: dict) -> None:
    config_p11 = {**config_p11, "graph_type": "simple_line_hourly"}
    render_section_header(config_p11["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p11["shared_country"] = setup_country_filter(
        config_p11, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p11, is_dual)
    config_p11["shared_year"] = str(shared_year)

    table_name = config_p11["table_name"]
    legend_col = config_p11["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1,
        table_name,
        year=str(shared_year),
        country=config_p11["shared_country"],
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return

    scenario_1_raw["snapshot"] = pd.to_datetime(scenario_1_raw["snapshot"])

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2,
            table_name,
            year=str(shared_year),
            country=config_p11["shared_country"],
        )
        if scenario_2_raw is not None and not scenario_2_raw.empty:
            scenario_2_raw["snapshot"] = pd.to_datetime(scenario_2_raw["snapshot"])

    has_dual_filters = (
        is_dual and scenario_2_raw is not None and not scenario_2_raw.empty
    )
    config_p11["shared_region"] = setup_region_filter(
        config_p11,
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        has_dual_filters,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        config_p11,
        has_dual_filters,
    )
    config_p11.update(filter_results)

    scenario_1_month, start_date, end_date, is_complete = (
        get_filtered_df_and_date_range(scenario_1_raw, config_p11)
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        scenario_1_month, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_month, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, config_p11
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            scenario_2_month, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered, scenario_2_filtered if has_dual_filters else None, None
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_simple_line_hourly(
            scenario_1_filtered,
            config_p11,
            colour_mapping,
            y_range,
            start_date,
            end_date,
            is_complete,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_without_data(
            scenario_1_filtered, config_p11, st.session_state.sce1
        )
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_simple_line_hourly(
                scenario_1_filtered,
                config_p11,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_simple_line_hourly(
                scenario_2_filtered,
                config_p11,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(
                scenario_1_filtered, config_p11, st.session_state.sce1
            )
        with col3:
            render_download_without_data(
                scenario_2_filtered, config_p11, st.session_state.sce2
            )

    st.divider()


def render_p12_nodal_flow_between_regions(config_p12: dict) -> None:
    config_p12 = {**config_p12, "graph_type": "simple_line_hourly"}
    render_section_header(config_p12["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p12["shared_country"] = setup_country_filter(
        config_p12, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p12, is_dual)
    config_p12["shared_year"] = str(shared_year)

    table_name = config_p12["table_name"]
    legend_col = config_p12["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1,
        table_name,
        year=str(shared_year),
        country=config_p12["shared_country"],
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return

    scenario_1_raw["snapshot"] = pd.to_datetime(scenario_1_raw["snapshot"])

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2,
            table_name,
            year=str(shared_year),
            country=config_p12["shared_country"],
        )
        if scenario_2_raw is not None and not scenario_2_raw.empty:
            scenario_2_raw["snapshot"] = pd.to_datetime(scenario_2_raw["snapshot"])

    has_dual_filters = (
        is_dual and scenario_2_raw is not None and not scenario_2_raw.empty
    )
    config_p12["shared_region"] = setup_region_filter(
        config_p12,
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        has_dual_filters,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        config_p12,
        has_dual_filters,
    )
    config_p12.update(filter_results)

    scenario_1_month, start_date, end_date, is_complete = (
        get_filtered_df_and_date_range(scenario_1_raw, config_p12)
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        scenario_1_month, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)
    scenario_1_filtered = add_nice_names(scenario_1_filtered, legend_col, mapping_df)

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_month, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, config_p12
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            scenario_2_month, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)
        scenario_2_filtered = add_nice_names(
            scenario_2_filtered, legend_col, mapping_df
        )

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered, scenario_2_filtered if has_dual_filters else None, None
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_filtered, legend_col
        )
        plot_simple_line_hourly(
            scenario_1_filtered,
            config_p12,
            colour_mapping,
            y_range,
            start_date,
            end_date,
            is_complete,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_without_data(
            scenario_1_filtered, config_p12, st.session_state.sce1
        )
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_filtered, legend_col
            )
            plot_simple_line_hourly(
                scenario_1_filtered,
                config_p12,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_filtered, legend_col
            )
            plot_simple_line_hourly(
                scenario_2_filtered,
                config_p12,
                colour_mapping,
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(
                scenario_1_filtered, config_p12, st.session_state.sce1
            )
        with col3:
            render_download_without_data(
                scenario_2_filtered, config_p12, st.session_state.sce2
            )

    st.divider()


def render_p13_battery_ep_ratio(config_p13: dict) -> None:
    config_p13 = {**config_p13, "graph_type": "simple_bar_yearly"}
    render_section_header(config_p13["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p13["shared_country"] = setup_country_filter(
        config_p13, is_dual, scenario_tag=st.session_state.sce1
    )

    table_name = config_p13["table_name"]
    legend_col = config_p13["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1, table_name, country=config_p13["shared_country"]
    )
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2, table_name, country=config_p13["shared_country"]
        )
        scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
        scenario_2_grouped = scenario_2_grouped.groupby(
            ["year", legend_col], as_index=False
        )["value"].sum()
        scenario_2_grouped = add_nice_names(scenario_2_grouped, legend_col, mapping_df)

    y_range = calculate_min_max_y_scale(scenario_1_grouped, scenario_2_grouped, "year")
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        colour_mapping = get_colour_mapping(
            table_name, mapping_df, scenario_1_grouped, legend_col
        )
        plot_simple_bar_yearly(
            scenario_1_grouped,
            config_p13,
            colour_mapping,
            y_range,
            key=(
                "plotly_chart_"
                f"{config_p13['download_id'].format(st.session_state.sce1)}"
            ),
        )
        render_download_with_table(scenario_1_raw, config_p13, st.session_state.sce1)
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_1_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_1_grouped,
                config_p13,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p13['download_id'].format(st.session_state.sce1)}"
                ),
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            colour_mapping = get_colour_mapping(
                table_name, mapping_df, scenario_2_grouped, legend_col
            )
            plot_simple_bar_yearly(
                scenario_2_grouped,
                config_p13,
                colour_mapping,
                y_range,
                key=(
                    "plotly_chart_"
                    f"{config_p13['download_id'].format(st.session_state.sce2)}"
                ),
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_with_table(
                scenario_1_raw, config_p13, st.session_state.sce1
            )
        with col3:
            render_download_with_table(
                scenario_2_raw, config_p13, st.session_state.sce2
            )

    st.divider()


def render_p14_battery_charging_profile(config_p14: dict) -> None:
    config_p14 = {**config_p14, "graph_type": "line_with_secondary_y_hourly"}
    render_section_header(config_p14["name"])

    is_dual = st.session_state.sce2 and st.session_state.sce2 != ""

    config_p14["shared_country"] = setup_country_filter(
        config_p14, is_dual, scenario_tag=st.session_state.sce1
    )
    shared_year = setup_year_filter(config_p14, is_dual)
    config_p14["shared_year"] = str(shared_year)

    table_name = config_p14["table_name"]
    legend_col = config_p14["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    scenario_1_raw = read_result_csv(
        st.session_state.sce1,
        table_name,
        year=str(shared_year),
        country=config_p14["shared_country"],
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return

    scenario_1_raw["snapshot"] = pd.to_datetime(scenario_1_raw["snapshot"])

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.sce2,
            table_name,
            year=str(shared_year),
            country=config_p14["shared_country"],
        )
        if scenario_2_raw is not None and not scenario_2_raw.empty:
            scenario_2_raw["snapshot"] = pd.to_datetime(scenario_2_raw["snapshot"])

    has_dual_filters = (
        is_dual and scenario_2_raw is not None and not scenario_2_raw.empty
    )
    config_p14["shared_region"] = setup_region_filter(
        config_p14,
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        has_dual_filters,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw if has_dual_filters else None,
        config_p14,
        has_dual_filters,
    )
    config_p14.update(filter_results)

    scenario_1_month, start_date, end_date, is_complete = (
        get_filtered_df_and_date_range(scenario_1_raw, config_p14)
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        scenario_1_month, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)

    scenario_2_filtered = None
    if has_dual_filters:
        scenario_2_month, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, config_p14
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            scenario_2_month, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)

    y_range = calculate_min_max_y_scale(
        scenario_1_filtered,
        scenario_2_filtered if has_dual_filters else None,
        None,
    )
    y_range = {"max_scale": y_range["max"], "min_scale": y_range["min"]}

    label_map = {
        label: (
            mapping_df.loc[label, "nice_names"]
            if (mapping_df is not None and label in mapping_df.index)
            else prettify_label(label)
        )
        for label in config_p14["primary_y_lab"] + config_p14["secondary_y_lab"]
    }
    config_p14["label_map"] = label_map

    if not is_dual:
        st.caption(f"{st.session_state.sce1} ")
        plot_line_with_secondary_y_hourly(
            scenario_1_filtered,
            config_p14,
            get_colour_mapping(
                table_name,
                mapping_df,
                add_nice_names(scenario_1_filtered, legend_col, mapping_df),
                legend_col,
            ),
            y_range,
            start_date,
            end_date,
            is_complete,
            key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
        )
        render_download_without_data(
            scenario_1_filtered, config_p14, st.session_state.sce1
        )
    else:
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.sce1} ")
            plot_line_with_secondary_y_hourly(
                scenario_1_filtered,
                config_p14,
                get_colour_mapping(
                    table_name,
                    mapping_df,
                    add_nice_names(scenario_1_filtered, legend_col, mapping_df),
                    legend_col,
                ),
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce1}_{table_name}",
            )
        with col3:
            st.caption(f"{st.session_state.sce2} ")
            plot_line_with_secondary_y_hourly(
                scenario_2_filtered,
                config_p14,
                get_colour_mapping(
                    table_name,
                    mapping_df,
                    add_nice_names(scenario_2_filtered, legend_col, mapping_df),
                    legend_col,
                ),
                y_range,
                start_date,
                end_date,
                is_complete,
                key=f"plotly_chart_{st.session_state.sce2}_{table_name}",
            )

        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data(
                scenario_1_filtered, config_p14, st.session_state.sce1
            )
        with col3:
            render_download_without_data(
                scenario_2_filtered, config_p14, st.session_state.sce2
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

generate_sidebar(table_of_content)
