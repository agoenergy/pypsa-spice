# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Plot-only helpers for Power page charts."""

import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots
import os
import re
import pandas as pd

from scripts.plot_settings import (
    add_stackedbar_total,
    handle_y_axis_list,
    update_hourly_plot_x_axis,
    configure_plot_layout,
)


def create_nice_names_and_color_mapping(
    table_name: str,
) -> pd.DataFrame:
    """Get the names to hex codes mapping df for a given graph.

    Parameters
    ----------
    table_name : str
        Tab name of the graph as per the config file

    Returns
    -------
    pd.DataFrame
        Dataframe of tech/carrier csv or None if no mapping file matches
    """
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
    df = pd.read_csv(file_path, index_col="original_names")

    return df


@st.fragment
def plot_simple_bar_yearly(df, graph_config, colour_mapping, y_range, key):
    """Plot yearly stacked bar chart from pre-processed data."""
    fig = px.bar(
        df,
        x="year",
        y="value",
        color="nice_names",
        barmode=(
            "group"
            if graph_config.get("table_name") == "pow_bats_ep_ratio"
            else "stack"
        ),
        color_discrete_map=colour_mapping,
    )
    fig = add_stackedbar_total(fig, df)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    # For the yearly bar charts, adjust the space between bars
    fig.update_layout(bargap=0.4)
    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_simple_line_yearly(df, graph_config, colour_mapping, y_range, key):
    """Plot yearly line chart from pre-processed data."""
    fig = px.line(
        df,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_area_share_yearly(df, graph_config, colour_mapping, y_range, key):
    """Plot yearly area chart from pre-processed data.

    Note: y_range parameter is accepted for interface consistency but not used.
    """
    fig = px.area(
        df,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = configure_plot_layout(fig, df, None, graph_config)
    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_bar_with_filter(df, graph_config, colour_mapping, y_range, key):
    """Plot yearly stacked bar chart with pre-applied filter."""
    fig = px.bar(
        df,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = add_stackedbar_total(fig, df)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    # For the yearly bar charts, adjust the space between bars
    fig.update_layout(bargap=0.4)
    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_simple_bar_hourly(
    df,
    graph_config,
    colour_mapping,
    y_range,
    start_date,
    end_date,
    is_complete,
    key,
):
    """Plot hourly stacked bar chart from pre-filtered data."""
    fig = px.bar(
        df,
        x="snapshot",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, df, start_date, end_date, is_complete)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_simple_line_hourly(
    df,
    graph_config,
    colour_mapping,
    y_range,
    start_date,
    end_date,
    is_complete,
    key,
):
    """Plot hourly line chart from pre-filtered data."""
    fig = px.line(
        df,
        x="snapshot",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, df, start_date, end_date, is_complete)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_filtered_bar_hourly(
    filtered_df,
    graph_config,
    colour_mapping,
    y_range,
    start_date,
    end_date,
    is_complete,
    key,
):
    """Plot hourly stacked bar chart with line overlay from pre-filtered data."""
    fig = px.bar(
        filtered_df[filtered_df["value"] != 0],
        x="snapshot",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, filtered_df, start_date, end_date, is_complete)
    fig = configure_plot_layout(fig, filtered_df, y_range, graph_config)

    fil_col = graph_config.get("fil_col")
    if fil_col and fil_col in filtered_df.columns:
        df_line = (
            filtered_df.groupby(["snapshot", fil_col])["value"].sum().reset_index()
        )
        line_chart_trace = px.line(df_line, x="snapshot", y="value")
        line_chart_trace.update_traces(line={"color": "blue", "width": 3})
        for trace in line_chart_trace.data:
            fig.add_trace(trace)

    st.plotly_chart(fig, use_container_width=True, key=key)


@st.fragment
def plot_line_with_secondary_y_hourly(
    df,
    graph_config,
    colour_mapping,
    y_range,
    start_date,
    end_date,
    is_complete,
    key,
):
    """Plot hourly line chart with secondary y-axis from pre-filtered data."""
    leg_col = graph_config["leg_col"]
    primary_y_lab = graph_config["primary_y_lab"]
    secondary_y_lab = graph_config["secondary_y_lab"]
    label_map = graph_config.get("label_map", {})

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    for prim_y in primary_y_lab:
        nice_name = label_map.get(prim_y, prim_y)
        fig.add_trace(
            go.Line(
                y=df[df[leg_col] == prim_y]["value"],
                x=df["snapshot"],
                name=nice_name,
                line={"color": colour_mapping.get(nice_name, "#a0a0a0")},
            ),
            secondary_y=False,
        )

    for secd_y in secondary_y_lab:
        nice_name = label_map.get(secd_y, secd_y)
        fig.add_trace(
            go.Line(
                y=df[df[leg_col] == secd_y]["value"],
                x=df["snapshot"],
                name=nice_name,
                line={"color": colour_mapping.get(nice_name, "#a0a0a0")},
            ),
            secondary_y=True,
        )

    fig = update_hourly_plot_x_axis(fig, df, start_date, end_date, is_complete)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    fig.update_yaxes(title_text=handle_y_axis_list(primary_y_lab), secondary_y=False)
    fig.update_yaxes(title_text=handle_y_axis_list(secondary_y_lab), secondary_y=True)
    st.plotly_chart(fig, use_container_width=True, key=key)
