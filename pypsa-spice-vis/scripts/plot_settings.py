# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Utility functions for plot configuration and shared settings."""

import datetime as dt
import os
import re
from itertools import cycle
from typing import Any

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from plotly.graph_objs._figure import Figure

from scripts.data_utils import (
    prettify_label,
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
