# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Utility functions for plot configuration and shared settings."""

import datetime as dt
import os
import re
from collections.abc import Callable
from itertools import cycle
from typing import Any

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.graph_objs._figure import Figure
from plotly.subplots import make_subplots

from scripts.data_utils import (
    calculate_min_max_y_scale,
    clean_df_for_plotting,
    filter_dataframe_by_date_range,
    get_filtered_df_and_date_range,
    get_hourly_dfs_for_both_scenarios,
    handle_small_values,
    prettify_label,
    read_result_csv,
)

# pylint: disable=too-many-locals

# =============================================================================
# Colour + label helpers
# =============================================================================


def get_default_colour_list() -> list[str]:
    """Return default colours used when no mapping CSV provides hex codes."""
    return ["#64B9E4", "#48A8AE", "#AD86B0", "#1E83B3", "#8393BE", "#637596"]


def create_nice_names_and_color_mapping(table_name: str) -> pd.DataFrame | None:
    """Return the mapping DataFrame (nice names + hex codes) for a given table."""
    pattern_name_map = [
        (r"^ene_avg_fuel_costs_fuel_yearly$", "carrier_mapping.csv"),
        (r"_by_heatgroup", "tech_mapping.csv"),
        (r"_capex", "tech_mapping.csv"),
        (r"_opex", "tech_mapping.csv"),
        (r"_by_type_by_carrier", "carrier_mapping.csv"),
        (r"_by_carrier", "carrier_mapping.csv"),
        (r"_by_type", "tech_mapping.csv"),
        (r"_by_category", "tech_mapping.csv"),
    ]

    file_name = None
    for pattern, name in pattern_name_map:
        if re.search(pattern, table_name):
            file_name = name
            break

    if file_name is None:
        return None

    file_path = os.path.join(st.session_state.current_dir, f"setting/{file_name}")
    return pd.read_csv(file_path, index_col="original_names")


def handle_y_axis_list(title_list: list[str]) -> str:
    """Convert a list of y-axis labels into a single prettified string."""
    prettified_list = [
        re.sub(r"([a-z])([A-Z])", r"\1 \2", label).capitalize() for label in title_list
    ]
    return ", ".join(prettified_list)


def generate_default_colour_mapping_dict_for_chart(
    df: pd.DataFrame, leg_col: str
) -> dict[str, str]:
    """Generate a default colour mapping dict based on prettified legend values."""
    unique_legends = df[leg_col].unique().tolist()
    prettified_legends = [prettify_label(label) for label in unique_legends]
    return dict(zip(prettified_legends, cycle(get_default_colour_list())))


def generate_color_mapping_dict_for_chart(
    table_name: str, legend_labels: list[str] | None = None
) -> dict[str, str]:
    """Return a mapping dict of legend label -> hex colour for a given chart."""
    mapping_df = create_nice_names_and_color_mapping(table_name)
    if mapping_df is None:
        return {}

    default_colours = get_default_colour_list()

    # Remove entries missing hex codes
    mapping_df = mapping_df.dropna(subset=["hex_codes"])

    colour_dict = mapping_df["hex_codes"].to_dict()
    mapping_dict = {
        mapping_df.loc[k, "nice_names"] if k in mapping_df.index else k: v
        for k, v in colour_dict.items()
    }

    # Assign default colours for legends not present in mapping_df
    default_colour_index = 0
    if legend_labels:
        for label in legend_labels:
            if label not in mapping_dict:
                mapping_dict[label] = default_colours[
                    default_colour_index % len(default_colours)
                ]
                default_colour_index += 1

    return mapping_dict


# =============================================================================
# Plot layout + axis helpers
# =============================================================================


def add_stackedbar_total(fig: Figure, df: pd.DataFrame) -> Figure:
    """Add total values on top of each stacked bar in a Plotly bar chart."""
    if "year" not in df.columns or "value" not in df.columns:
        raise ValueError("Dataframe must contain 'year' and 'value' columns.")

    totals = df.groupby("year")["value"].sum()
    for year, total in totals.items():
        fig.add_trace(
            go.Scatter(
                x=[year],
                y=[total],
                mode="text",
                text=f"{round(total, 1)}",
                textposition="top center",
                showlegend=False,
            )
        )
    return fig


def update_hourly_plot_x_axis(
    fig: Figure,
    filtered_df: pd.DataFrame,
    start_date: dt.datetime,
    end_date: dt.datetime,
    is_complete: bool,
) -> Figure:
    """Set the x-axis ticks for hourly graphs based on completeness and range."""
    if is_complete:
        time_diff = end_date - start_date
        tick_spacing = (
            1
            if time_diff <= dt.timedelta(hours=24)
            else max(2, round(time_diff.total_seconds() / 3600 / 10))
        )
        fig.update_xaxes(
            type="date",
            tickformat="%d/%m %H:%M",
            dtick=3600000 * tick_spacing,
            tickangle=40,
            tickfont={"size": 11},
        )
        return fig

    unique_snapshots = filtered_df["snapshot"].unique()
    num_points = len(unique_snapshots)

    max_num_points = (
        16
        if (
            st.session_state.window_width < 1130
            and st.session_state.window_width > 608
            and st.session_state.sce2 != ""
        )
        else 24
    )

    if num_points <= max_num_points:
        tick_positions = unique_snapshots
        tick_labels = [str(i + 1) for i in range(num_points)]
    else:
        step = max(1, num_points // 10)
        tick_positions = unique_snapshots[::step]
        tick_labels = ["1"] + [str(i) for i in range(step + 1, num_points + 1, step)]

    fig.update_xaxes(
        tickmode="array",
        ticktext=tick_labels,
        tickvals=tick_positions,
        tickangle=0,
        tickfont={"size": 13},
    )
    return fig


# pylint: disable=too-many-branches
def configure_plot_layout(
    fig: Figure,
    df: pd.DataFrame,
    yaxis_scales: dict[str, float] | None = None,
    graph_config: dict[str, Any] | None = None,
) -> Figure:
    """Apply consistent styling and optional y-axis scaling to a Plotly figure."""
    legend_orientation = "v"
    legend_x_pos = 1.05
    legend_y_pos = 0
    legend_x_anchor = "left"
    legend_y_anchor = "bottom"
    margin_b = 0

    if st.session_state.sce2 != "":
        legend_orientation = "h"
        legend_x_pos = 0.5
        legend_y_pos = -0.25
        legend_x_anchor = "center"
        legend_y_anchor = "top"
        margin_b = 150

    xaxis_tickfont_size = getattr(fig.layout.xaxis.tickfont, "size", None)
    x_tick_font_size = 15 if xaxis_tickfont_size is None else xaxis_tickfont_size

    units = ""
    if graph_config and "units" in graph_config:
        units = f" [{graph_config['units']}]"

    if yaxis_scales is None and graph_config and "yaxis_scales" in graph_config:
        yaxis_scales = graph_config["yaxis_scales"]

    layout_dict: dict[str, Any] = {
        "showlegend": True,
        "font": {"family": "Flexo, sans-serif", "size": 15},
        "legend": {
            "orientation": legend_orientation,
            "y": legend_y_pos,
            "x": legend_x_pos,
            "xanchor": legend_x_anchor,
            "yanchor": legend_y_anchor,
            "title_text": graph_config["leg_col"].capitalize() if graph_config else "",
            "title_font_size": 16,
        },
        "margin": {"t": 50, "b": margin_b},
        "xaxis": {"tickfont": {"size": x_tick_font_size}},
        "yaxis": {"tickfont": {"size": 15}},
        "xaxis_title": "",
        "yaxis_title": "",
        "annotations": [
            {
                "text": f"{units}",
                "x": 0,
                "y": 1.02,
                "xref": "paper",
                "yref": "paper",
                "xanchor": "center",
                "yanchor": "bottom",
                "xshift": -15,
                "yshift": 20,
                "showarrow": False,
                "font": {"size": 15},
            }
        ],
    }

    if yaxis_scales is not None:
        fig.update_yaxes(range=[yaxis_scales["min_scale"], yaxis_scales["max_scale"]])

    if "year" in df.columns:
        years_with_data = df["year"].unique()
        layout_dict["xaxis"].update(
            {"tickvals": list(years_with_data), "tickmode": "array"}
        )

    fig.update_layout(**layout_dict)
    return fig


# =============================================================================
# Scenario data helpers used by chart builders
# =============================================================================


def get_yearly_dfs_for_both_scenarios(
    graph_config: dict[str, Any],
    func: Callable[[pd.DataFrame], pd.DataFrame] | None = None,
) -> tuple[pd.DataFrame | None, pd.DataFrame | None]:
    """Return yearly DataFrames for both scenarios, optionally processed with func."""
    dfs: list[pd.DataFrame | None] = []
    table_name = graph_config["table_name"]
    shared_country = graph_config["shared_country"]

    for scenario in [st.session_state.sce1, st.session_state.sce2]:
        df = read_result_csv(scenario, table_name, shared_country)
        if df is not None and func:
            df = func(df)
        dfs.append(df)

    return dfs[0], dfs[1]


# =============================================================================
# Yearly charts
# =============================================================================


@st.fragment
def simple_bar_yearly(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate a yearly stacked bar chart."""
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    table_name = graph_config["table_name"]

    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )
    if df is None or df.empty:
        return

    df = clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(table_name)
    df_grouped = df.groupby(["year", leg_col], as_index=False)["value"].sum()

    df_grouped["nice_names"] = df_grouped[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = df_grouped["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(df, leg_col)
    )

    df1_grouped: pd.DataFrame | None = df_grouped
    df2_grouped: pd.DataFrame | None = None
    if graph_config.get("shared_country") is not None and st.session_state.sce2:
        df1_grouped, df2_grouped = get_yearly_dfs_for_both_scenarios(
            graph_config,
            func=lambda d: d.groupby(["year", leg_col], as_index=False)["value"].sum(),
        )

    y_range = calculate_min_max_y_scale(df1_grouped, df2_grouped, "year")
    fig = px.bar(
        df_grouped,
        x="year",
        y="value",
        color="nice_names",
        barmode="group" if table_name == "pow_bats_ep_ratio" else "stack",
        color_discrete_map=colour_mapping,
    )
    fig = add_stackedbar_total(fig, df_grouped)
    fig = configure_plot_layout(
        fig,
        df_grouped,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )
    st.plotly_chart(fig, use_container_width=True, key=f"plotly_chart_{download_id}")


@st.fragment
def simple_line_yearly(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate a yearly line chart."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]

    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )
    if df is None or df.empty:
        return

    df = clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(table_name)
    df = df.groupby(["year", leg_col], as_index=False)["value"].mean()
    df["nice_names"] = df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = df["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(df, leg_col)
    )

    df1: pd.DataFrame | None = df
    df2: pd.DataFrame | None = None
    if graph_config.get("shared_country") is not None and st.session_state.sce2:
        df1, df2 = get_yearly_dfs_for_both_scenarios(graph_config, None)

    y_range = calculate_min_max_y_scale(df1, df2, None)
    fig = px.line(
        df, x="year", y="value", color="nice_names", color_discrete_map=colour_mapping
    )
    fig = configure_plot_layout(
        fig,
        df,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def bar_with_filter(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate a yearly stacked bar chart with a filter."""
    leg_col = graph_config["leg_col"]
    fil_col = graph_config["fil_col"]
    slider_id = graph_config["slider_id"].format(scenario_name)
    table_name = graph_config["table_name"]

    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )
    if df is None or df.empty:
        return

    mapping_df = create_nice_names_and_color_mapping(table_name)

    shared_filter = graph_config.get("shared_filter")
    if shared_filter is None:
        shared_filter = st.radio(
            f"{slider_id} Select {fil_col} ({scenario_name}):",
            options=df[fil_col].unique(),
            format_func=prettify_label,
            horizontal=True,
            label_visibility="collapsed",
        )

    df_reg = df.loc[df[fil_col] == shared_filter].copy()
    df_reg["nice_names"] = df_reg[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = df_reg["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(df, leg_col)
    )

    df1_reg: pd.DataFrame | None = df_reg
    df2_reg: pd.DataFrame | None = None
    if graph_config.get("shared_filter") is not None and st.session_state.sce2:
        df1_reg, df2_reg = get_yearly_dfs_for_both_scenarios(
            graph_config, func=lambda d: d.loc[d[fil_col] == shared_filter]
        )

    y_range = calculate_min_max_y_scale(df1_reg, df2_reg, "year")
    df_reg = clean_df_for_plotting(leg_col, df_reg)

    fig = px.bar(
        df_reg,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = add_stackedbar_total(fig, df_reg)
    fig = configure_plot_layout(
        fig,
        df_reg,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def area_share_yearly(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate yearly area chart (percentage share)."""
    leg_col = graph_config["leg_col"]
    table_name = graph_config["table_name"]

    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )
    if df is None or df.empty:
        return

    df = clean_df_for_plotting(leg_col, df)
    mapping_df = create_nice_names_and_color_mapping(table_name)

    df = df.groupby(["year", leg_col], as_index=False)["value"].mean()
    df["nice_names"] = df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = df["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(df, leg_col)
    )

    fig = px.area(
        df, x="year", y="value", color="nice_names", color_discrete_map=colour_mapping
    )
    fig = configure_plot_layout(fig, df, None, graph_config)
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


# =============================================================================
# Hourly charts
# =============================================================================


def _validate_date_range(start_date: dt.datetime, end_date: dt.datetime) -> bool:
    """Return True if the range is valid, else show Streamlit error and return False."""
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        return False
    return True


@st.fragment
def simple_bar_hourly(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate hourly stacked bar chart for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)

    df = read_result_csv(
        scenario_name,
        table_name,
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )
    if df is None or df.empty:
        return

    df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        df, graph_config
    )
    if not _validate_date_range(start_date, end_date):
        return

    mapping_df = create_nice_names_and_color_mapping(table_name)

    filtered_df = filter_dataframe_by_date_range(
        df_m, start_date=start_date, end_date=end_date
    ).copy()
    filtered_df["nice_names"] = filtered_df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = filtered_df["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(filtered_df, leg_col)
    )

    filtered_df = handle_small_values(filtered_df)

    fig = px.bar(
        filtered_df,
        x="snapshot",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, filtered_df, start_date, end_date, is_complete)

    if st.session_state.sce2:
        df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
        y_range = calculate_min_max_y_scale(df1_sel, df2_sel, "snapshot")
    else:
        y_range = calculate_min_max_y_scale(filtered_df, None, "snapshot")

    fig = configure_plot_layout(
        fig,
        filtered_df,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )
    st.plotly_chart(fig, use_container_width=True, key=f"plotly_chart_{download_id}")


@st.fragment
def simple_line_hourly(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate hourly line chart with filters for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]

    df = read_result_csv(
        scenario_name,
        table_name,
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )
    if df is None or df.empty:
        return

    df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        df, graph_config
    )
    if not _validate_date_range(start_date, end_date):
        return

    mapping_df = create_nice_names_and_color_mapping(table_name)

    filtered_df = filter_dataframe_by_date_range(
        df_m, start_date=start_date, end_date=end_date
    ).copy()
    filtered_df["nice_names"] = filtered_df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = filtered_df["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(filtered_df, leg_col)
    )

    filtered_df = handle_small_values(filtered_df)

    fig = px.line(
        filtered_df,
        x="snapshot",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, filtered_df, start_date, end_date, is_complete)

    if st.session_state.sce2:
        df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
        y_range = calculate_min_max_y_scale(df1_sel, df2_sel, None)
    else:
        y_range = calculate_min_max_y_scale(filtered_df, None, None)

    fig = configure_plot_layout(
        fig,
        filtered_df,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def filtered_bar_hourly(scenario_name: str, graph_config: dict[str, Any]) -> None:
    """Generate hourly stacked bar chart with optional line overlay."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    fil_col = graph_config["fil_col"]

    df = read_result_csv(
        scenario_name,
        table_name,
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )
    if df is None or df.empty:
        return

    df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        df, graph_config
    )
    if not _validate_date_range(start_date, end_date):
        return

    mapping_df = create_nice_names_and_color_mapping(table_name)

    filtered_df = filter_dataframe_by_date_range(
        df_m, start_date=start_date, end_date=end_date
    ).copy()
    filtered_df["nice_names"] = filtered_df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )
    filtered_df = handle_small_values(filtered_df)

    unique_legends = filtered_df["nice_names"].unique().tolist()
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(filtered_df, leg_col)
    )

    fig = px.bar(
        filtered_df[filtered_df["value"] != 0],
        x="snapshot",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, filtered_df, start_date, end_date, is_complete)

    if st.session_state.sce2:
        df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
        y_range = calculate_min_max_y_scale(df1_sel, df2_sel, "snapshot")
    else:
        y_range = calculate_min_max_y_scale(filtered_df, None, "snapshot")

    fig = configure_plot_layout(
        fig,
        filtered_df,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )

    if fil_col in filtered_df.columns:
        df_line = (
            filtered_df.groupby(["snapshot", fil_col])["value"].sum().reset_index()
        )
        line_chart_trace = px.line(df_line, x="snapshot", y="value")
        line_chart_trace.update_traces(line={"color": "blue", "width": 3})
        for trace in line_chart_trace.data:
            fig.add_trace(trace)

    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def line_with_secondary_y_hourly(
    scenario_name: str, graph_config: dict[str, Any]
) -> None:
    """Generate hourly line chart with secondary y-axis for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    primary_y_lab = graph_config["primary_y_lab"]
    secondary_y_lab = graph_config["secondary_y_lab"]

    df = read_result_csv(
        scenario_name,
        table_name,
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )
    if df is None or df.empty:
        return

    df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        df, graph_config
    )
    if not _validate_date_range(start_date, end_date):
        return

    mapping_df = create_nice_names_and_color_mapping(table_name)
    colour_mapping = (
        generate_color_mapping_dict_for_chart(table_name)
        if mapping_df is not None
        else generate_default_colour_mapping_dict_for_chart(df_m, leg_col)
    )

    filtered_df = filter_dataframe_by_date_range(
        df_m, start_date=start_date, end_date=end_date
    )

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    for prim_y in primary_y_lab:
        nice_name = (
            mapping_df.loc[prim_y, "nice_names"]
            if (mapping_df is not None and prim_y in mapping_df.index)
            else prettify_label(prim_y)
        )
        fig.add_trace(
            go.Line(
                y=filtered_df[filtered_df[leg_col] == prim_y]["value"],
                x=filtered_df["snapshot"],
                name=nice_name,
                line={"color": colour_mapping.get(nice_name, "#a0a0a0")},
            ),
            secondary_y=False,
        )

    for secd_y in secondary_y_lab:
        nice_name = (
            mapping_df.loc[secd_y, "nice_names"]
            if (mapping_df is not None and secd_y in mapping_df.index)
            else prettify_label(secd_y)
        )
        fig.add_trace(
            go.Line(
                y=filtered_df[filtered_df[leg_col] == secd_y]["value"],
                x=filtered_df["snapshot"],
                name=nice_name,
                line={"color": colour_mapping.get(nice_name, "#a0a0a0")},
            ),
            secondary_y=True,
        )

    fig = update_hourly_plot_x_axis(fig, filtered_df, start_date, end_date, is_complete)

    if st.session_state.sce2:
        df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
        y_range = calculate_min_max_y_scale(df1_sel, df2_sel, None)
    else:
        y_range = calculate_min_max_y_scale(filtered_df, None, None)

    fig = configure_plot_layout(
        fig,
        filtered_df,
        {"max_scale": y_range["max"], "min_scale": y_range["min"]},
        graph_config,
    )
    fig.update_yaxes(title_text=handle_y_axis_list(primary_y_lab), secondary_y=False)
    fig.update_yaxes(title_text=handle_y_axis_list(secondary_y_lab), secondary_y=True)

    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )
