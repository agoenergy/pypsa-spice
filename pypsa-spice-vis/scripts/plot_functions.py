# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Plotly helper functions for pypsa-spice-vis charts."""

# pylint: disable=too-many-arguments,too-many-locals, too-many-positional-arguments
import os
import re
from collections.abc import Mapping
from typing import Any, cast

import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from plotly.subplots import make_subplots
from streamlit_sortables import sort_items

from scripts.plot_settings import (
    add_stackedbar_total,
    configure_plot_layout,
    handle_y_axis_list,
    update_hourly_plot_x_axis,
)


def create_nice_names_and_color_mapping(table_name: str) -> pd.DataFrame | None:
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

    file_path = os.path.join(
        st.session_state.streamlit_base_dir, f"setting/{file_name}"
    )
    df = pd.read_csv(file_path, index_col="original_names")

    return df


def get_legend_order_by_first_years_value(df: pd.DataFrame) -> list[str]:
    """Return legend order sorted by value in the base year (descending).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the data with columns "year", "legend", and "value".

    Returns
    -------
    list[str]
        List of legend names sorted by their value in the base year (descending).
    """
    if df.empty or "year" not in df.columns:
        return []

    first_x = df["year"].sort_values().iloc[0]
    all_legends = pd.Index(df["legend"].unique())

    first_x_values = (
        df[df["year"] == first_x]
        .groupby("legend", sort=False)["value"]
        .sum()
        .reindex(all_legends, fill_value=0)
        .sort_values(ascending=False, kind="stable")
    )

    return first_x_values.index.tolist()


def get_current_legend_order(df: pd.DataFrame) -> list[str]:
    """Return the current legend order based on the dataframe rows."""
    if df.empty or "legend" not in df.columns:
        return []

    return cast(list[str], pd.Index(df["legend"].dropna().unique()).tolist())


def render_hourly_legend_order_control(
    legends: list[str],
    key: str,
) -> list[str]:
    """Render a drag-and-drop control for ordering hourly chart legends."""
    if len(legends) < 2:
        return legends

    state_key = f"{key}_legend_order"
    stored_legends = st.session_state.get(state_key, [])
    ordered_legends = [legend for legend in stored_legends if legend in legends]
    ordered_legends.extend(
        legend for legend in legends if legend not in ordered_legends
    )

    with st.popover(
        "Legend order",
        icon=":material/reorder:",
        key=f"{key}_legend_order_popover",
    ):
        st.caption("Drag the items to reorder the current chart legend.")
        ordered_legends = sort_items(
            ordered_legends,
            direction="vertical",
            key=f"{key}_legend_order_sortable",
            custom_style="""
                .sortable-component.vertical {
                    gap: 0.35rem;
                }
                .sortable-item {
                    background: rgba(248, 249, 251, 0.95);
                    border: 1px solid rgba(49, 51, 63, 0.16);
                    border-radius: 0.5rem;
                    color: rgb(49, 51, 63);
                    cursor: grab;
                    font-size: 0.9rem;
                    margin: 0;
                    padding: 0.45rem 0.7rem;
                }
            """,
        )
        if st.button("Reset", key=f"{key}_legend_order_reset", width="stretch"):
            ordered_legends = legends

    st.session_state[state_key] = ordered_legends
    return ordered_legends


@st.fragment
def plot_diff_bar_yearly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    key: str,
) -> None:
    """Plot yearly stacked bar chart showing the difference between two scenarios.

    Positive values indicate scenario 2 is higher; negative values indicate
    scenario 1 is higher. A reference line is drawn at y=0.
    """
    fig = px.bar(
        df,
        x="year",
        y="value",
        color="legend",
        barmode="group",
        color_discrete_map=colour_mapping,
    )
    fig.add_hline(y=0, line_width=1.5, line_color="black")
    units = f" [{graph_config['units']}]" if "units" in graph_config else ""
    fig.update_layout(
        font={"family": "Flexo, sans-serif", "size": 15},
        legend={
            "orientation": "h",
            "x": 0.5,
            "y": -0.25,
            "xanchor": "center",
            "yanchor": "top",
            "title_text": graph_config["leg_col"].capitalize(),
            "title_font_size": 16,
        },
        margin={"t": 50, "b": 150},
        xaxis_title="",
        yaxis_title="",
        bargap=0.6,
        bargroupgap=0.1,
        annotations=[
            {
                "text": units,
                "x": 0,
                "y": 1.02,
                "xref": "paper",
                "yref": "paper",
                "xanchor": "center",
                "yanchor": "bottom",
                "xshift": -15,
                "yshift": 20,
                "showarrow": False,
                "font": {"size": 13},
            }
        ],
    )
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_share_comparison_lines(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    key: str,
) -> None:
    """Plot two scenario lines for a share metric (e.g. renewable share) over years.

    Parameters
    ----------
    df : pd.DataFrame
        Combined dataframe with columns ``year``, ``value``, ``scenario``.
    graph_config : dict
        Chart config — only ``units`` is used for the y-axis annotation.
    key : str
        Unique Streamlit/Plotly chart key.
    """
    fig = px.line(
        df,
        x="year",
        y="value",
        color="scenario",
        markers=True,
    )
    units = f" [{graph_config['units']}]" if "units" in graph_config else ""
    fig.update_layout(
        font={"family": "Flexo, sans-serif", "size": 15},
        legend={
            "orientation": "h",
            "x": 0.5,
            "y": -0.25,
            "xanchor": "center",
            "yanchor": "top",
            "title_text": "",
            "title_font_size": 16,
        },
        margin={"t": 50, "b": 150},
        xaxis_title="",
        yaxis_title="",
        yaxis={"range": [0, 100]},
        annotations=[
            {
                "text": units,
                "x": 0,
                "y": 1.02,
                "xref": "paper",
                "yref": "paper",
                "xanchor": "center",
                "yanchor": "bottom",
                "xshift": -15,
                "yshift": 20,
                "showarrow": False,
                "font": {"size": 13},
            }
        ],
    )
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_simple_bar_yearly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    key: str,
) -> None:
    """Plot yearly stacked bar chart from pre-processed data."""
    legend_order = get_legend_order_by_first_years_value(df)
    fig = px.bar(
        df,
        x="year",
        y="value",
        color="legend",
        barmode=(
            "group"
            if graph_config.get("table_name") == "pow_bats_ep_ratio"
            else "stack"
        ),
        color_discrete_map=colour_mapping,
        category_orders={"legend": legend_order},
    )
    fig = add_stackedbar_total(fig, df)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    # For the yearly bar charts, adjust the space between bars
    fig.update_layout(bargap=0.4)
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_simple_line_yearly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    key: str,
) -> None:
    """Plot yearly line chart from pre-processed data."""
    fig = px.line(
        df,
        x="year",
        y="value",
        color="legend",
        color_discrete_map=colour_mapping,
    )
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_area_share_yearly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    key: str,
) -> None:
    """Plot yearly area chart from pre-processed data.

    Note: y_range parameter is accepted for interface consistency but not used.
    """
    fig = px.area(
        df,
        x="year",
        y="value",
        color="legend",
        color_discrete_map=colour_mapping,
    )
    fig = configure_plot_layout(fig, df, None, graph_config)
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_bar_with_filter(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    key: str,
) -> None:
    """Plot yearly stacked bar chart with pre-applied filter."""
    legend_order = get_legend_order_by_first_years_value(df)
    fig = px.bar(
        df,
        x="year",
        y="value",
        color="legend",
        color_discrete_map=colour_mapping,
        category_orders={"legend": legend_order},
    )
    fig = add_stackedbar_total(fig, df)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    # For the yearly bar charts, adjust the space between bars
    fig.update_layout(bargap=0.4)
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_simple_bar_hourly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    start_date: Any,
    end_date: Any,
    is_complete: bool,
    key: str,
) -> None:
    """Plot hourly stacked bar chart from pre-filtered data."""
    legend_order = render_hourly_legend_order_control(
        get_current_legend_order(df),
        key,
    )
    fig = px.bar(
        df,
        x="snapshot",
        y="value",
        color="legend",
        color_discrete_map=colour_mapping,
        category_orders={"legend": legend_order},
    )
    fig = update_hourly_plot_x_axis(fig, df, start_date, end_date, is_complete)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_simple_line_hourly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    start_date: Any,
    end_date: Any,
    is_complete: bool,
    key: str,
) -> None:
    """Plot hourly line chart from pre-filtered data."""
    fig = px.line(
        df,
        x="snapshot",
        y="value",
        color="legend",
        color_discrete_map=colour_mapping,
    )
    fig = update_hourly_plot_x_axis(fig, df, start_date, end_date, is_complete)
    fig = configure_plot_layout(fig, df, y_range, graph_config)
    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_filtered_bar_hourly(
    filtered_df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    start_date: Any,
    end_date: Any,
    is_complete: bool,
    key: str,
) -> None:
    """Plot hourly stacked bar chart with line overlay from pre-filtered data."""
    bar_df = filtered_df[filtered_df["value"] != 0]
    legend_order = render_hourly_legend_order_control(
        get_current_legend_order(bar_df),
        key,
    )
    fig = px.bar(
        bar_df,
        x="snapshot",
        y="value",
        color="legend",
        color_discrete_map=colour_mapping,
        category_orders={"legend": legend_order},
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

    st.plotly_chart(fig, width="stretch", key=key)


@st.fragment
def plot_line_with_secondary_y_hourly(
    df: pd.DataFrame,
    graph_config: Mapping[str, Any],
    colour_mapping: Mapping[str, str],
    y_range: Mapping[str, Any],
    start_date: Any,
    end_date: Any,
    is_complete: bool,
    key: str,
) -> None:
    """Plot hourly line chart with secondary y-axis from pre-filtered data."""
    leg_col = graph_config["leg_col"]
    primary_y_lab = graph_config["primary_y_lab"]
    secondary_y_lab = graph_config["secondary_y_lab"]
    label_map = graph_config.get("label_map", {})

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    for prim_y in primary_y_lab:
        nice_name = label_map.get(prim_y, prim_y)
        fig.add_trace(
            go.Scatter(
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
            go.Scatter(
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
    st.plotly_chart(fig, width="stretch", key=key)
