# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Plot setting functions that are used across multiple plot types."""

import datetime as dt
import os
import re
from collections.abc import Callable
from itertools import cycle

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

# pylint: disable=too-many-locals, broad-exception-caught

# =========================== General functions for plotting =========================


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
        Dataframe of tech/carrier csv
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


def handle_y_axis_list(title_list: list[str]) -> str:
    """Convert a list of y-axis labels into a single, prettified string.

    This function will split the input list into separate words and convert any
    camelCase strings into capitalised plain text.

    Parameters
    ----------
    title_list : list
        _description_

    Returns
    -------
    str
        _description_
    """
    prettified_list = [
        re.sub(r"([a-z])([A-Z])", r"\1 \2", label).capitalize() for label in title_list
    ]
    y_axis_title = ", ".join(prettified_list)
    return y_axis_title


def handle_color_mapping_for_chart(
    table_name: str, legend_labels: list[str] | None = None
) -> dict[str, str]:
    """Get the colour mapping dict for a given graph.

    Parameters
    ----------
    table_name : str
        Tab name of the graph as per the config file
    legend_labels : list, optional
        List of unique legends for the current graph

    Returns
    -------
    Dict[str, str]
        Dict of nice names to hex code mapping
    """
    df = create_nice_names_and_color_mapping(table_name)

    default_colours = get_default_colour_list()

    if df is None:
        return {}

    # Remove entries from the mapping df that are missing a hex code
    df = df.dropna(subset=["hex_codes"])

    colour_dict = df["hex_codes"].to_dict()
    nice_mapping = {
        df.loc[k, "nice_names"] if k in df.index else k: v
        for k, v in colour_dict.items()
    }

    # For legends that are not present in the mapping df (either because they were
    # dropped earlier, or because they don't have an assigned colour), cycle through
    # the default colours and assign a hex code
    default_colour_index = 0
    if legend_labels:
        for label in legend_labels:
            if label not in nice_mapping:
                nice_mapping[label] = default_colours[
                    default_colour_index % len(default_colours)
                ]
                default_colour_index += 1

    return nice_mapping


def generate_default_colour_mapping(df: pd.DataFrame, leg_col: str) -> dict[str, str]:
    """Generate a default colour mapping dictionary.

    Generate a colour mapping dict for the legend series in a graph using a default
    colour scheme. This function is called for charts that don't use the tech or
    carrier mapping csvs.

    Parameters
    ----------
    df : pd.DataFrame
        The dataframe to extract the legend series from.
    leg_col : str
        The name of the legend column.

    Returns
    -------
    Dict
        The dictionary of legend names to hex code mapping.
    """
    unique_legends = df[leg_col].unique().tolist()
    prettified_legends = [prettify_label(label) for label in unique_legends]

    default_colours = get_default_colour_list()

    mapping = dict(zip(prettified_legends, cycle(default_colours)))

    return mapping


def get_default_colour_list() -> list:
    """Get the list of default colours to use.

    Colours are assigned to legends that do not have a specified hex_code in
    tech_mapping or carrier_mapping. Currently using Agora EW default colours.

    Returns
    -------
    list
        The list of default colours
    """
    return ["#64B9E4", "#48A8AE", "#AD86B0", "#1E83B3", "#8393BE", "#637596"]


def get_yearly_dfs_for_both_scenarios(
    graph_config: dict, func: Callable[[pd.DataFrame], pd.DataFrame] | None = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the yearly dfs for two scenarios.

    optionally process (before min/max y calculation takes place).

    Parameters
    ----------
    graph_config : Dict
        Configuration dictionary for the current graph
    func: Callable
        Optional function to further process each dataframe (e.g., groupby or filter)

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        The dataframes (optionally processed) for the two scenarios
    """
    dfs = []
    table_name = graph_config["table_name"]
    shared_country = graph_config["shared_country"]

    for scenario in [st.session_state.sce1, st.session_state.sce2]:
        df = read_result_csv(scenario, table_name, shared_country)
        if func:
            df = func(df)
        dfs.append(df)

    return dfs


def update_hourly_plot_x_axis(
    fig: Figure,
    filtered_df: pd.DataFrame,
    start_date: dt.datetime,
    end_date: dt.datetime,
    is_complete: bool,
) -> Figure:
    """Set the x axis values for hourly graphs.

    The values are set based on a) whether data is complete (no discontinuous hours)
    and b) the number of selected hours - for fewer than 36 hours, all hours are shown,
    otherwise ~10 ticks are shown.

    Parameters
    ----------
    fig : Figure
        The plotly figure to update
    filtered_df : pd.DataFrame
        Dataframe filtered by the selected date range
    start_date : dt.datetime
        Selected start date
    end_date : dt.datetime
        Selected end date
    is_complete : bool
        Whether the data is complete

    Returns
    -------
    Figure
        The updated plotly figure with the appropriate x-axis settings applied
    """
    if is_complete:
        time_diff = end_date - start_date
        # Set x axis ticks based on no. of hours selected - for fewer than 24 hours,
        # show all hours, otherwise show ~10 ticks
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
    else:
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
        # Set x axis ticks based on no. of hours selected - for fewer than 24 hours,
        # show increment of 1, otherwise show ~10 evenly spaced ticks
        if num_points <= max_num_points:
            tick_positions = unique_snapshots
            tick_labels = [str(i + 1) for i in range(num_points)]
        else:
            step = max(1, num_points // 10)
            tick_positions = unique_snapshots[::step]
            tick_labels = ["1"] + [
                str(i) for i in range(step + 1, num_points + 1, step)
            ]

        fig.update_xaxes(
            tickmode="array",
            ticktext=tick_labels,
            tickvals=tick_positions,
            tickangle=0,
            tickfont={"size": 13},
        )

    return fig


def update_layout(
    fig: Figure, df: pd.DataFrame, yaxis_scales: dict = None, graph_config: dict = None
) -> Figure:
    """Update the layout of a Plotly figure to improve readability and aesthetics.

    This applies consistent styling (set in layout_dict) and y axis scaling across all
    graphs.

    Parameters
    ----------
    fig : Figure
        The plotly figure object to update
    df : pd.DataFrame
        DataFrame containing the data to plot
    yaxis_scales : dict, optional
        Dictionary containing y-axis scale settings with keys 'min_scale' and
            'max_scale', by default None
    graph_config : dict, optional
        Configuration dictionary that may contain:
            - 'units': str, which will be displayed at the top of the y-axis
            - 'yaxis_scales': dict, fallback y-axis scales if yaxis_scales is None

    Returns
    -------
    Figure
        Updated The updated plotly figure with modified layout settings
    """
    # Default legend orientation and position: vertical, to the left of graph
    legend_orientation = "v"
    legend_x_pos = 1.05
    legend_y_pos = 0
    legend_x_anchor = "left"
    legend_y_anchor = "bottom"
    margin_b = 0
    # In narrow widths (but before graphs have stacked), in two scenario cases adjust
    # the legend orientation and position: horizontal and below the graph
    if (
        st.session_state.sce2 != ""
        and st.session_state.window_width < 1130
        and st.session_state.window_width > 608
    ):
        legend_orientation = "h"
        legend_x_pos = 0.2
        legend_y_pos = -0.4
        legend_x_anchor = "center"
        legend_y_anchor = "top"
        margin_b = 100

    # Check if the x tick font size has been set already (in the case of hourly graphs
    # this is set in _update_hourly_plot_x_axis, and we do not want to overwrite it)
    xaxis_tickfont_size = getattr(fig.layout.xaxis.tickfont, "size", None)
    x_tick_font_size = 15 if xaxis_tickfont_size is None else xaxis_tickfont_size

    # Set units for the graph
    units = ""
    if graph_config and "units" in graph_config:
        units = f" [{graph_config['units']}]"

    # Handle y-axis scaling
    if yaxis_scales is None and graph_config and "yaxis_scales" in graph_config:
        yaxis_scales = graph_config["yaxis_scales"]

    layout_dict = {
        "showlegend": True,
        "font": {"family": "Flexo, sans-serif", "size": 15},
        "legend": dict(
            orientation=legend_orientation,
            y=legend_y_pos,
            x=legend_x_pos,
            xanchor=legend_x_anchor,
            yanchor=legend_y_anchor,
            title_text=graph_config["leg_col"].capitalize(),
            title_font_size=16,
        ),  # pylint: disable=use-dict-literal
        "margin": {"t": 50, "b": margin_b},
        "xaxis": {"tickfont": {"size": x_tick_font_size}},
        "yaxis": {"tickfont": {"size": 15}},
        "xaxis_title": "",
        "yaxis_title": "",
        "annotations": [
            dict(
                text=f"{units}",
                x=0,
                y=1.02,
                xref="paper",
                yref="paper",
                xanchor="center",
                yanchor="bottom",
                xshift=-15,
                yshift=20,
                showarrow=False,
                font=dict(size=15),
            )  # pylint: disable=use-dict-literal
        ],
    }

    # For the yearly bar charts, adjust the space between bars
    if graph_config["graph_type"] in [
        "simple_bar_yearly",
        "simple_bar_yearly_2",
        "bar_with_filter",
    ]:
        fig.update_layout(bargap=0.4)

    if yaxis_scales is not None:
        fig.update_yaxes(range=[yaxis_scales["min_scale"], yaxis_scales["max_scale"]])

    if "year" in df.columns:
        years_with_data = df["year"].unique()
        layout_dict["xaxis"].update(
            {"tickvals": list(years_with_data), "tickmode": "array"}
        )

    fig.update_layout(**layout_dict)

    return fig


# =========================== Add different charts =========================


@st.fragment
def simple_bar_yearly(scenario_name: str, graph_config: dict) -> None:
    """Generate a yearly stacked bar chart with a downloadable data table."""
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    table_name = graph_config["table_name"]

    # Construct file path
    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        country=graph_config["shared_country"],
    )
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
        handle_color_mapping_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    df1_grouped = df_grouped
    df2_grouped = None
    if graph_config.get("shared_country") is not None and st.session_state.sce2:
        df1_grouped, df2_grouped = get_yearly_dfs_for_both_scenarios(
            graph_config,
            func=lambda df: df.groupby(["year", leg_col], as_index=False)[
                "value"
            ].sum(),
        )

    y_range = calculate_min_max_y_scale(df1_grouped, df2_grouped, "year")
    min_y = y_range["min"]
    max_y = y_range["max"]

    # Generate and display the chart
    try:
        fig = px.bar(
            df_grouped,
            x="year",
            y="value",
            color="nice_names",
            barmode="group" if table_name == "pow_bats_ep_ratio" else "stack",
            color_discrete_map=colour_mapping,
        )

        # Add stacked bar total and update layout
        fig = add_stackedbar_total(fig, df_grouped)
        fig = update_layout(
            fig, df_grouped, {"max_scale": max_y, "min_scale": min_y}, graph_config
        )

        # Display the chart with a unique key
        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{download_id}"
        )

    except ValueError as e:
        st.error(f"ValueError encountered: {e}")
        st.dataframe(df)


@st.fragment
def simple_line_yearly(scenario_name: str, graph_config: dict):
    """Generate a yearly line chart with a downloadable data table."""
    # Read data from file
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]

    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )

    df = clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(table_name)

    df["nice_names"] = df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = df["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    df1 = df
    df2 = None
    if graph_config.get("shared_country") is not None and st.session_state.sce2:
        df1, df2 = get_yearly_dfs_for_both_scenarios(graph_config, None)

    y_range = calculate_min_max_y_scale(df1, df2, None)
    min_y = y_range["min"]
    max_y = y_range["max"]

    fig = px.line(
        df, x="year", y="value", color="nice_names", color_discrete_map=colour_mapping
    )

    fig = update_layout(fig, df, {"max_scale": max_y, "min_scale": min_y}, graph_config)
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def bar_with_filter(scenario_name: str, graph_config: dict):
    """
    Generate a yearly stacked bar chart.

    Chart with a filter and  with a downloadable data table.
    """
    leg_col = graph_config["leg_col"]
    fil_col = graph_config["fil_col"]
    slider_id = graph_config["slider_id"].format(scenario_name)
    table_name = graph_config["table_name"]
    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )

    mapping_df = create_nice_names_and_color_mapping(table_name)

    # This creates a shared filter that applies to both graphs if the entry is found in
    # the config dict, otherwise a local filter is generated for the single graph.
    if "shared_filter" in graph_config:
        shared_filter = graph_config["shared_filter"]
    else:
        shared_filter = st.radio(
            f"{slider_id} Select {fil_col} ({scenario_name})" + ":",
            options=df[fil_col].unique(),
            format_func=prettify_label,
            horizontal=True,
            label_visibility="collapsed",
        )

    df_reg = df.copy()
    df_reg = df_reg.loc[df_reg[fil_col] == shared_filter]
    df_reg["nice_names"] = df_reg[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )

    unique_legends = df_reg["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    # If scenario2 exists, recalculate the maximum y tick value considering both
    # datasets, else calculate the maximum y for just the single dataset
    df1_reg = df_reg
    df2_reg = None
    if graph_config.get("shared_filter") is not None and st.session_state.sce2:
        df1_reg, df2_reg = get_yearly_dfs_for_both_scenarios(
            graph_config, func=lambda df: df.loc[df[fil_col] == shared_filter]
        )

    y_range = calculate_min_max_y_scale(df1_reg, df2_reg, "year")
    min_y = y_range["min"]
    max_y = y_range["max"]

    df_reg = clean_df_for_plotting(leg_col, df_reg)

    fig = px.bar(
        df_reg,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )

    fig = add_stackedbar_total(fig, df_reg)
    fig = update_layout(
        fig, df_reg, {"max_scale": max_y, "min_scale": min_y}, graph_config
    )
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def area_share_yearly(scenario_name: str, graph_config: dict):
    """Generate yearly area chart (percentage share) with a downloadable data table."""
    leg_col = graph_config["leg_col"]
    table_name = graph_config["table_name"]

    df = read_result_csv(
        scenario_name, table_name, country=graph_config["shared_country"]
    )
    df = clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(table_name)

    df["nice_names"] = df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else prettify_label(x)
        )
    )
    unique_legends = df["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    # Plot stacked bar chart using Plotly
    fig = px.area(
        df,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )

    fig = update_layout(fig, df, None, graph_config)
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{scenario_name}_{table_name}"
    )


@st.fragment
def simple_bar_hourly(scenario_name: str, graph_config: dict[str, str]) -> None:
    """Generate hourly stacked bar chart for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
            df, graph_config
        )

        # Validate date range
        if start_date > end_date:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        mapping_df = create_nice_names_and_color_mapping(table_name)

        # Filter data by date range
        filtered_df = filter_dataframe_by_date_range(
            df_m, start_date=start_date, end_date=end_date
        )
        filtered_df = filtered_df.copy()
        filtered_df["nice_names"] = filtered_df[leg_col].map(
            lambda x: (
                mapping_df.loc[x, "nice_names"]
                if (mapping_df is not None and x in mapping_df.index)
                else prettify_label(x)
            )
        )
        unique_legends = filtered_df["nice_names"].unique().tolist()
        colour_mapping = (
            handle_color_mapping_for_chart(table_name, unique_legends)
            if mapping_df is not None
            else generate_default_colour_mapping(filtered_df, leg_col)
        )

        filtered_df = handle_small_values(filtered_df)

        fig = px.bar(
            filtered_df,
            x="snapshot",
            y="value",
            color="nice_names",
            color_discrete_map=colour_mapping,
        )

        fig = update_hourly_plot_x_axis(
            fig, filtered_df, start_date, end_date, is_complete
        )

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
            y_range = calculate_min_max_y_scale(df1_sel, df2_sel, "snapshot")
        else:
            y_range = calculate_min_max_y_scale(filtered_df, None, "snapshot")

        fig = update_layout(
            fig,
            filtered_df,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            graph_config,
        )
        # Display the chart with a unique key
        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{download_id}"
        )


@st.fragment
def simple_line_hourly(scenario_name: str, graph_config: dict):
    """Generate hourly line chart with filters for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
            df, graph_config
        )

        mapping_df = create_nice_names_and_color_mapping(table_name)

        if start_date <= end_date:
            pass
        else:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        filtered_df = filter_dataframe_by_date_range(
            df_m, start_date=start_date, end_date=end_date
        )
        filtered_df = filtered_df.copy()
        filtered_df["nice_names"] = filtered_df[leg_col].map(
            lambda x: (
                mapping_df.loc[x, "nice_names"]
                if (mapping_df is not None and x in mapping_df.index)
                else prettify_label(x)
            )
        )
        unique_legends = filtered_df["nice_names"].unique().tolist()
        colour_mapping = (
            handle_color_mapping_for_chart(table_name, unique_legends)
            if mapping_df is not None
            else generate_default_colour_mapping(filtered_df, leg_col)
        )

        filtered_df = handle_small_values(filtered_df)

        fig = px.line(
            filtered_df,
            x="snapshot",
            y="value",
            color="nice_names",
            color_discrete_map=colour_mapping,
        )

        fig = update_hourly_plot_x_axis(
            fig, filtered_df, start_date, end_date, is_complete
        )

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
            y_range = calculate_min_max_y_scale(df1_sel, df2_sel, None)
        else:
            y_range = calculate_min_max_y_scale(filtered_df, None, None)

        fig = update_layout(
            fig,
            filtered_df,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            graph_config,
        )
        st.plotly_chart(
            fig,
            use_container_width=True,
            key=f"plotly_chart_{scenario_name}_{table_name}",
        )


@st.fragment
def filtered_bar_hourly(scenario_name: str, graph_config: dict):
    """Generate hourly stacked bar chart with filters for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    fil_col = graph_config["fil_col"]

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
            df, graph_config
        )

        mapping_df = create_nice_names_and_color_mapping(table_name)

        if start_date <= end_date:
            pass
        else:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        filtered_df = filter_dataframe_by_date_range(
            df_m, start_date=start_date, end_date=end_date
        )
        filtered_df = filtered_df.copy()
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
            handle_color_mapping_for_chart(table_name, unique_legends)
            if mapping_df is not None
            else generate_default_colour_mapping(filtered_df, leg_col)
        )

        fig = px.bar(
            filtered_df[filtered_df["value"] != 0],
            x="snapshot",
            y="value",
            color="nice_names",
            color_discrete_map=colour_mapping,
        )

        fig = update_hourly_plot_x_axis(
            fig, filtered_df, start_date, end_date, is_complete
        )

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
            y_range = calculate_min_max_y_scale(df1_sel, df2_sel, "snapshot")
        else:
            y_range = calculate_min_max_y_scale(filtered_df, None, "snapshot")

        fig = update_layout(
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

        df_line = (
            filtered_df.groupby(["snapshot", fil_col])["value"].sum().reset_index()
        )
        line_chart_trace = px.line(df_line, x="snapshot", y="value")
        line_chart_trace.update_traces(line={"color": "blue", "width": 3})

        for trace in line_chart_trace.data:
            fig.add_trace(trace)

        st.plotly_chart(
            fig,
            use_container_width=True,
            key=f"plotly_chart_{scenario_name}_{table_name}",
        )


@st.fragment
def line_with_secondary_y_hourly(scenario_name: str, graph_config: dict):
    """Generate hourly line chart with secondary y-axis for datetime."""
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    primary_y_lab = graph_config["primary_y_lab"]
    secondary_y_lab = graph_config["secondary_y_lab"]

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = get_filtered_df_and_date_range(
            df, graph_config
        )

        mapping_df = create_nice_names_and_color_mapping(table_name)
        colour_mapping = (
            handle_color_mapping_for_chart(table_name)
            if mapping_df is not None
            else generate_default_colour_mapping(df_m, leg_col)
        )

        if start_date <= end_date:
            pass
        else:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        filtered_df = filter_dataframe_by_date_range(
            df_m, start_date=start_date, end_date=end_date
        )

        # Create figure with secondary y-axis
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

        fig = update_hourly_plot_x_axis(
            fig, filtered_df, start_date, end_date, is_complete
        )

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = get_hourly_dfs_for_both_scenarios(graph_config)
            y_range = calculate_min_max_y_scale(df1_sel, df2_sel, None)
        else:
            y_range = calculate_min_max_y_scale(filtered_df, None, None)

        fig = update_layout(
            fig,
            filtered_df,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            graph_config,
        )
        fig.update_yaxes(
            title_text=handle_y_axis_list(primary_y_lab), secondary_y=False
        )
        fig.update_yaxes(
            title_text=handle_y_axis_list(secondary_y_lab), secondary_y=True
        )
        st.plotly_chart(
            fig,
            use_container_width=True,
            key=f"plotly_chart_{scenario_name}_{table_name}",
        )


@st.fragment
def sankey_diagram(
    scenario_name: str,
    graph_config: dict
):

    table_name = graph_config['table_name']
    years = graph_config['years']

    
    # Read data from file
    st.text("year:")
    year = st.radio(
        "Select Year ({})".format(scenario_name) + ":",
        options=years,
        horizontal=True,
        label_visibility="collapsed",
        key=table_name + "_" + scenario_name + "_year",
    )

    file_path = os.path.abspath(
        st.session_state.result_path
        + "/"
        + scenario_name
        + "/csvs/"
        + st.session_state.sector
        + "/"
        + str(year)
        + "/"
        + table_name
        + ".csv"
    )

    try:
        df = pd.read_csv(os.path.abspath(file_path))
        df = df[df['year']==year] 
    except FileNotFoundError:
        with st.container(height=450, border=True):
            st.write(":material/warning: File not found: {}".format(file_path))
            return None

    # Prepare the unique source and target mapping
    unique_source_target = list(pd.unique(df[['source', 'target']].values.ravel('K')))
    mapping_dict = {k: v for v, k in enumerate(unique_source_target)}

    # Map the source and target to unique numbers
    df['source'] = df['source'].map(mapping_dict)
    df['target'] = df['target'].map(mapping_dict)

    # Convert the DataFrame to a dictionary for Plotly
    links_dict = df.to_dict(orient='list')

    # Function to convert HEX color to RGBA with a given opacity
    def hex_to_rgba(hex_color, opacity=0.9):
        hex_color = hex_color.lstrip('#')
        rgb = tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
        return f'rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity})'

    # Apply the opacity change to all link colors
    link_colors_with_opacity = [hex_to_rgba(color) for color in links_dict.get('color', [])]

    color_map = {
        "Indigenous": "#F2D7A6",
        "Imported": "#2C3E50",
        "Bio": "#8b737f",
        "Bit": "#163c47",
        "ENS": "#A9A9A9",
        "Gas": "#ff7967",
        "Geothermal": "#c0d88d",
        "Oil": "#b45340",
        "Solar": "#ffae63",
        "Uranium": "#967bb6",
        "Waste": "#8b737f",
        "Water": "#1d6897",
        "Wind": "#42b2b7",
        "Low_Heat": "#F2D7A6",
        "Electricity": "#4682b4",
        "Agriculture": "#4682b4",
        "Commercial": "#4682b4",
        "Industries": "#4682b4",
        "Non-Energy Use": "#4682b4",
        "Residential": "#4682b4",
        "Transportation": "#4682b4",
        "Hyd": "#5CC9F5",
        "Other": "#FFFFFF",
        }
    
    # Create the Sankey diagram figure
    fig = go.Figure(data=[go.Sankey(
        arrangement='snap',
        valuesuffix = " TWh",
        node=dict(
            label=unique_source_target,
            pad=20,  # Increase padding to reduce overlap and improve readability
            thickness=30,
            color=[color_map[x] for x in unique_source_target],
            align='justify',
            line=dict(color='rgba(0,0,0,0.4)', width=2),  # Darker border for better contrast
        ),
        link=dict(
            source=links_dict["source"],
            target=links_dict["target"],
            value=links_dict["value"],
            color=link_colors_with_opacity,  # Set the new colors with opacity
        ),
        textfont=dict(size=16, color='black', family='Arial, sans-serif', weight=600),  # Larger, bold text
    )])

    # Update layout settings
    fig.update_layout(
        font=dict(size=16, color='black', family='Arial, sans-serif', weight=600),
        width=900,
        height=650,
        paper_bgcolor='white',
        plot_bgcolor='white',
        margin=dict(l=20, r=20, t=40, b=40)
    )

    # Streamlit app
    st.plotly_chart(fig, use_container_width=True, key=f"sankey_diagram_{scenario_name}_{table_name}")


# =========================== Interactive Map Chart =========================


@st.fragment
def interactive_map(scenario_name: str, graph_config: dict) -> None:
    """Generate a network map showing generation capacity and flows.
    
    This function creates a two-panel visualization:
    - Left: Power generation flow with colored buses by technology
    - Right: Load demand requirements
    
    Parameters
    ----------
    scenario_name : str
        The name of the scenario to visualize.
    graph_config : dict
        Configuration dictionary containing:
        - shared_year: str - Year to plot
        - title: str - Optional title for the map
    """
    import matplotlib.pyplot as plt
    import pypsa
    import numpy as np
    import cartopy.crs as ccrs
    
    shared_year = graph_config.get("shared_year")
    title = graph_config.get("title", "Network Map")
    
    if not shared_year:
        st.warning("Please select a year to display the map.")
        return
    
    # Construct path to network file
    network_path = os.path.join(
        st.session_state.result_path,
        scenario_name,
        "post-solve",
        f"network_{st.session_state.sector}_{shared_year}.nc"
    )
    
    # Check if file exists
    if not os.path.exists(network_path):
        st.error(f"Network file not found: {network_path}")
        return
    
    try:
        # Load PyPSA network
        with st.spinner("Loading network..."):
            n = pypsa.Network(network_path)
        
        # Define color mappings for technologies and regions
        tech_colors = {
            'BIOT': "#A5CCA9",   'CCGT': "#A6A3B9", 
            'OCGT': "#E1D8EB",   'OILT': "#F7946A", 
            'SubC': "#4A122B",   'SupC': "#6E2C57", 
            'GEOT': "#3A683E",   'GEOX': "#3A683E", 
            'HROR': "#3399CC",   'PHOT': "#FFC000",
            'FLOT': "#FFA500",   'WTON': "#66BBC5",
            'WTOF': "#00A3B8",   'HDAM': "#21729B", 
            'HPHS': "#196D83",   'RTPV': "#FF8C00",
            'LSLO': "#0BDA51",   'BATS': "#F5A1AD",
            'BES1': '#008080',   'BES2': '#8A2BE2',
            'BES4': '#004080',   'NSMR': '#FF10F0',
            'NUCL': '#FF10F0', 
            'LUZ-N_to_LUZ-MM': '#1F77B4', 
            'LUZ-MM_to_LUZ-S': '#FF7F0E',
            'LUZ-S_to_V-LEY':  '#2CA02C', 
            'V-LEY_to_V-BHL': '#D62728',
            'V-LEY_to_V-CEB': '#9467BD', 
            'V-CEB_to_V-NEG': '#8C564B',
            'V-NEG_to_V-PAN': '#E377C2', 
            'V-CEB_to_MIN-NW': '#7F7F7F',
            'MIN-NW_to_MIN-NE': '#BCBD22', 
            'MIN-NE_to_MIN-S': '#17BECF',   
            'LUZ-MM_to_LUZ-N': '#1F77B4',
            'LUZ-S_to_LUZ-MM': '#FF7F0E',
            'V-LEY_to_LUZ-S': '#2CA02C',
            'V-BHL_to_V-LEY': '#D62728',
            'V-CEB_to_V-LEY': '#9467BD',
            'V-NEG_to_V-CEB': '#8C564B',
            'V-PAN_to_V-NEG': '#E377C2',
            'MIN-NW_to_V-CEB': '#7F7F7F',
            'MIN-NE_to_MIN-NW': '#BCBD22',
            'MIN-S_to_MIN-NE': '#17BECF',
        }
        
        region_colors = {
            'LUZ-N':  '#1F77B4',   
            'LUZ-MM': '#FF7F0E',
            'LUZ-S': '#2CA02C',
            'V-LEY': '#D62728',   
            'V-PAN': '#9467BD', 
            'V-CEB': '#8C564B',   
            'V-BHL': '#E377C2', 
            'V-NEG': '#7F7F7F',
            'MIN-NW': '#BCBD22',
            'MIN-NE': '#17BECF',
            'MIN-S': '#9EDAE5'
        }
        
        # Get all plants (exclude SUPPLY and LSLO)
        plants = n.generators[
            (~n.generators.type.str.contains('_SUPPLY', na=False)) &
            (~n.generators.type.str.contains('LSLO', na=False))
        ]
        
        plants = plants.assign(g=n.generators_t.p.sum()).groupby(["bus", "type"]).g.sum()
        
        # Get all power converters (exclude ITCN and BATS)
        pow_gens = n.links[
            (~n.links.type.str.contains('ITCN', na=False)) &
            (~n.links.type.str.contains('BATS', na=False)) &
            (~n.links.type.str.contains('IND-BOILER', na=False)) &
            (~n.links.type.str.contains('EVCH', na=False)) &
            (~n.links.type.str.contains('EDLH|EHPP|EIDT|EERH|ELTZ', na=False))
        ]
        
        pow_gens = -pow_gens.assign(g=n.links_t.p1.sum()).groupby(["bus1", "type"]).g.sum()
        
        # Concatenate generation data
        generation = pd.concat([plants, pow_gens]).sort_index(level=0)
        
        # Ensure all technology types in generation data have colors
        unique_types = generation.index.get_level_values(1).unique()
        default_colors = get_default_colour_list()
        for i, tech_type in enumerate(unique_types):
            if tech_type not in tech_colors:
                tech_colors[tech_type] = default_colors[i % len(default_colors)]
        
        # Remove non-ITCN links and LVELEC links
        n.mremove('Link', n.links[~n.links.type.str.contains('ITCN', na=False)].index)
        n.mremove('Link', n.links[n.links.index.str.contains('LVELEC', na=False)].index)
        
        # Remove duplicate buses, but groupby sum first
        n.links = n.links.reset_index().set_index(['bus0', 'bus1'])
        n.links.p_nom = n.links.groupby([n.links.index]).p_nom.sum()
        n.links = n.links.reset_index().set_index("Link")
        n.links = n.links.sort_index()
        
        duplicated = n.links[n.links.duplicated(['bus0', 'bus1'])].index
        n.mremove('Link', duplicated)
        
        # Calculate flow
        flow = n.links_t.p0.sum().to_frame()
        flow['component'] = 'Link'
        flow = flow.groupby(['component', flow.index]).sum().mean(axis=1).apply(lambda x: 5 * np.sign(x))
        
        link_color = n.links_t.p0.mean().abs()
        
        # Calculate loads - aggregate by bus and filter to only existing buses
        loads = n.loads_t.p.T
        loads = loads.groupby(loads.index.str[:-9]).sum().sum(axis=1)
        
        # Filter loads to only include buses that exist in the network
        loads = loads[loads.index.isin(n.buses.index)]
        
        # Create single plot figure
        fig, ax = plt.subplots(1, 1, figsize=(18, 10), subplot_kw={"projection": ccrs.PlateCarree()})
        
        # Plot the network - generation flow
        collection = n.plot(
            bus_sizes=generation / 15e7,
            bus_colors=tech_colors,
            margin=0.4,
            flow=flow,
            color_geomap=True,
            link_widths=1.5,
            link_colors=link_color,
            ax=ax
        )
        
        # Add colorbar if available
        try:
            if isinstance(collection, (list, tuple)) and len(collection) > 2:
                fig.colorbar(collection[2], ax=ax, fraction=0.04, pad=0.001, label="Flow in MW")
        except (KeyError, IndexError, AttributeError):
            pass  # Skip colorbar if not available
        
        # Add legend for technology colors
        from matplotlib.patches import Patch
        unique_techs = generation.index.get_level_values(1).unique()
        legend_elements = [
            Patch(facecolor=tech_colors.get(tech, '#888888'), edgecolor='black', label=tech)
            for tech in sorted(unique_techs)
        ]
        ax.legend(
            handles=legend_elements,
            loc='upper left',
            bbox_to_anchor=(1.05, 1),
            title='Technology',
            fontsize=10,
            title_fontsize=12,
            frameon=True,
            fancybox=True,
            shadow=True
        )
        
        plt.tight_layout()
        
        # Display in Streamlit
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)
        
    except Exception as e:
        st.error(f"Error loading or plotting network: {str(e)}")
        st.exception(e)


@st.fragment
def interactive_map_plotly(scenario_name: str, graph_config: dict) -> None:
    """Generate an interactive Plotly map with values plotted as circles or pie charts.
    
    This is an alternative to the PyPSA native plotting that provides interactivity.
    
    Parameters
    ----------
    scenario_name : str
        The name of the scenario to visualize.
    graph_config : dict
        Configuration dictionary containing:
        - table_name: str - Name of the data table to read
        - leg_col: str - Column name for legend/categories
        - download_id: str - Unique identifier for download button
        - map_type: str - Type of map ("circle" or "pie")
        - shared_country: str - Optional country filter
        - shared_year: str - Optional year filter for time-series data
        - title: str - Optional title for the map
        - unit: str - Optional unit label (e.g., "MW", "TWh")
    
    Notes
    -----
    - Only plots electricity carriers (filters out non-electricity buses)
    - Requires buses.csv with columns: bus, carrier, x (longitude), y (latitude)
    - Interactive tooltips show location name and values on hover
    """
    table_name = graph_config["table_name"]
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    map_type = graph_config.get("map_type", "pie")  # "circle" or "pie" - default is pie
    unit = graph_config.get("unit", "")
    
    # Read the data from results
    df = read_result_csv(
        scenario_name,
        table_name,
        country=graph_config.get("shared_country"),
        year=graph_config.get("shared_year"),
    )
    
    if df is None or df.empty:
        st.warning("No data available for the selected scenario and filters.")
        return
    
    # Read bus coordinates from input directory
    # Path structure: data/{data_folder}/{project}/input/{scenario}/power/buses.csv
    # Get the base data folder path from init config
    input_folder_path = st.session_state.input_data_folder_path
    # Remove "/input" suffix to get to project level, then reconstruct full path
    data_folder_base = os.path.dirname(input_folder_path)
    
    # Clean scenario name (remove "(FINAL)" prefix if present)
    clean_scenario = scenario_name.split("(FINAL) ")[-1].strip()
    
    input_base_path = os.path.join(
        data_folder_base,
        "input",
        clean_scenario,
        "power",
        "buses.csv"
    )
    
    try:
        buses_df = pd.read_csv(input_base_path)
    except FileNotFoundError:
        st.error(f"Could not find buses.csv at: {input_base_path}")
        return
    
    # Filter for electricity carriers only
    buses_df = buses_df[buses_df["carrier"] == "Electricity"].copy()
    
    if buses_df.empty:
        st.warning("No electricity buses found in buses.csv")
        return
    
    # Check if data has region column
    if "region" not in df.columns:
        st.error(
            f"Data must contain a 'region' column to map to geographic locations. "
            f"Available columns: {', '.join(df.columns)}"
        )
        return
    
    # Clean data for plotting
    df = clean_df_for_plotting(leg_col, df)
    df = handle_small_values(df)
    
    # Get one representative bus per region (prefer HVELEC buses for electricity data)
    buses_per_region = buses_df[buses_df["bus"].str.contains("HVELEC")].copy()
    
    # Extract region names from node column (remove country prefix if present)
    # e.g., "PH_LUZ-N" -> "LUZ-N"
    buses_per_region["region_clean"] = buses_per_region["node"].str.split("_").str[-1]
    
    # Get unique regions from data and buses for debugging
    data_regions = set(df["region"].unique())
    bus_regions = set(buses_per_region["region_clean"].unique())
    
    # Merge with bus locations using region -> region_clean mapping
    merged_df = pd.merge(
        df,
        buses_per_region[["node", "region_clean", "x", "y"]],
        left_on="region",
        right_on="region_clean",
        how="inner"
    )
    
    if merged_df.empty:
        # Provide detailed debugging information
        st.error("❌ No matching data found between results and bus locations.")
        
        with st.expander("🔍 Debug Information - Click to expand"):
            st.write("**Regions in your data:**")
            st.write(sorted(data_regions))
            
            st.write("**Regions found in buses.csv:**")
            st.write(sorted(bus_regions))
            
            st.write("**Regions that don't match:**")
            unmatched = data_regions - bus_regions
            if unmatched:
                st.write(f"In data but not in buses: {sorted(unmatched)}")
            
            unmatched_buses = bus_regions - data_regions
            if unmatched_buses:
                st.write(f"In buses but not in data: {sorted(unmatched_buses)}")
            
            st.write("**Suggestions:**")
            st.write("1. Check that region names match exactly (case-sensitive)")
            st.write("2. Verify buses.csv has HVELEC buses for all regions")
            st.write("3. Check for extra spaces or special characters")
        
        return
    
    # Get color mapping
    mapping_df = create_nice_names_and_color_mapping(table_name)
    
    if mapping_df is not None:
        merged_df["nice_names"] = merged_df[leg_col].map(
            lambda x: (
                mapping_df.loc[x, "nice_names"]
                if x in mapping_df.index
                else prettify_label(x)
            )
        )
    else:
        merged_df["nice_names"] = merged_df[leg_col].apply(prettify_label)
    
    unique_legends = merged_df["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(table_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(merged_df, "nice_names")
    )
    
    # Aggregate data by location and category
    grouped_df = merged_df.groupby(["node", "region_clean", "x", "y", leg_col, "nice_names"], as_index=False)["value"].sum()
    
    # Check for overlapping coordinates (multiple categories at same location)
    coord_counts = grouped_df.groupby(["x", "y"]).size()
    has_overlapping = (coord_counts > 1).any()
    
    # If coordinates overlap, force pie chart mode to make all data visible
    if has_overlapping and map_type == "circle":
        map_type = "pie"
        st.info("📊 Automatically using pie chart mode because multiple data points share the same coordinates.")
    
    # Create the map based on map_type
    if map_type == "pie":
        # Create pie chart markers for each location
        fig = _create_pie_chart_map(grouped_df, colour_mapping, unit, graph_config)
    else:
        # Create circle/bubble markers (default)
        fig = _create_circle_map(grouped_df, colour_mapping, unit, graph_config)
    
    # Update map layout
    title = graph_config.get("title", "Geographic Distribution")
    fig.update_layout(
        title=title,
        geo=dict(
            scope="asia",
            projection_type="mercator",
            showland=True,
            landcolor="rgb(229, 229, 229)",
            showlakes=True,
            lakecolor="rgb(209, 230, 245)",
            showocean=True,
            oceancolor="rgb(230, 245, 255)",
            showcountries=True,
            countrycolor="rgb(150, 150, 150)",
            countrywidth=1.5,
            showcoastlines=True,
            coastlinecolor="rgb(100, 100, 100)",
            coastlinewidth=1,
            showframe=True,
            framecolor="rgb(100, 100, 100)",
            framewidth=2,
            lonaxis=dict(
                range=[merged_df["x"].min() - 1, merged_df["x"].max() + 1],
                showgrid=True,
                gridwidth=0.5,
                gridcolor="rgb(200, 200, 200)"
            ),
            lataxis=dict(
                range=[merged_df["y"].min() - 1, merged_df["y"].max() + 1],
                showgrid=True,
                gridwidth=0.5,
                gridcolor="rgb(200, 200, 200)"
            ),
            resolution=50,  # Higher resolution for more detail
        ),
        showlegend=True,
        height=900,
        width=1400,
        margin=dict(l=20, r=20, t=60, b=20),
    )
    
    # Display the map
    st.plotly_chart(fig, use_container_width=False, key=f"map_{download_id}")


def _create_circle_map(
    df: pd.DataFrame,
    colour_mapping: dict,
    unit: str,
    graph_config: dict
) -> go.Figure:
    """Create a map with circle markers sized by value.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with columns: node, x, y, nice_names, value
    colour_mapping : dict
        Mapping of category names to colors
    unit : str
        Unit label for values
    graph_config : dict
        Configuration dictionary
        
    Returns
    -------
    go.Figure
        Plotly figure object
    """
    leg_col = graph_config["leg_col"]
    
    # Aggregate total by location for sizing
    location_totals = df.groupby(["node", "x", "y"], as_index=False)["value"].sum()
    location_totals = location_totals.rename(columns={"value": "total_value"})
    
    # Merge totals back
    df_with_totals = df.merge(location_totals, on=["node", "x", "y"])
    
    # Calculate marker sizes (scale appropriately)
    max_val = df_with_totals["total_value"].max()
    min_size, max_size = 10, 50
    df_with_totals["marker_size"] = (
        (df_with_totals["total_value"] / max_val) * (max_size - min_size) + min_size
    )
    
    # Create figure
    fig = go.Figure()
    
    # Add a trace for each category
    for category in df_with_totals["nice_names"].unique():
        category_df = df_with_totals[df_with_totals["nice_names"] == category]
        
        # Create hover text
        hover_text = [
            f"<b>{row['node']}</b><br>"
            f"{row['nice_names']}: {row['value']:.2f} {unit}<br>"
            f"Total at location: {row['total_value']:.2f} {unit}"
            for _, row in category_df.iterrows()
        ]
        
        fig.add_trace(go.Scattergeo(
            lon=category_df["x"],
            lat=category_df["y"],
            mode="markers",
            marker=dict(
                size=category_df["marker_size"],
                color=colour_mapping.get(category, "#888888"),
                line=dict(width=1, color="white"),
                sizemode="diameter",
            ),
            name=category,
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
        ))
    
    return fig


def _create_pie_chart_map(
    df: pd.DataFrame,
    colour_mapping: dict,
    unit: str,
    graph_config: dict
) -> go.Figure:
    """Create a map with pie chart markers at each location.
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with columns: node, x, y, nice_names, value
    colour_mapping : dict
        Mapping of category names to colors
    unit : str
        Unit label for values
    graph_config : dict
        Configuration dictionary
        
    Returns
    -------
    go.Figure
        Plotly figure object
    """
    from plotly import graph_objects as go
    
    fig = go.Figure()
    
    # Group by location
    locations = df.groupby(["node", "x", "y"])
    
    for (node, lon, lat), location_df in locations:
        # Calculate total for this location
        total_value = location_df["value"].sum()
        
        if total_value == 0:
            continue
        
        # Get categories and values for this location
        categories = location_df["nice_names"].tolist()
        values = location_df["value"].tolist()
        colors = [colour_mapping.get(cat, "#888888") for cat in categories]
        
        # Create hover text
        hover_text = f"<b>{node}</b><br>Total: {total_value:.2f} {unit}<br><br>"
        for cat, val in zip(categories, values):
            pct = (val / total_value) * 100
            hover_text += f"{cat}: {val:.2f} {unit} ({pct:.1f}%)<br>"
        
        # Add pie chart trace for this location
        # Size pie chart based on total value
        max_total = df.groupby(["node"])["value"].sum().max()
        size_scale = (total_value / max_total) * 0.3 + 0.05  # Scale between 0.05 and 0.35 degrees
        
        # Add a scatter point for interactivity
        fig.add_trace(go.Scattergeo(
            lon=[lon],
            lat=[lat],
            mode="markers",
            marker=dict(
                size=20,
                color="rgba(0,0,0,0)",  # Transparent
            ),
            name=node,
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
            showlegend=False,
        ))
        
        # Note: Plotly doesn't natively support pie charts on geo maps
        # We'll use the Scattergeo approach with sized circles colored by dominant category
        dominant_idx = values.index(max(values))
        dominant_color = colors[dominant_idx]
        
        # Calculate marker size based on total value
        min_size, max_size = 15, 60
        marker_size = (total_value / df.groupby("node")["value"].sum().max()) * (max_size - min_size) + min_size
        
        fig.add_trace(go.Scattergeo(
            lon=[lon],
            lat=[lat],
            mode="markers",
            marker=dict(
                size=marker_size,
                color=dominant_color,
                line=dict(width=2, color="white"),
                opacity=0.7,
            ),
            name=f"{node}",
            text=hover_text,
            hovertemplate="%{text}<extra></extra>",
            showlegend=False,
        ))
    
    # Add legend entries for categories
    for category, color in colour_mapping.items():
        if category in df["nice_names"].values:
            fig.add_trace(go.Scattergeo(
                lon=[None],
                lat=[None],
                mode="markers",
                marker=dict(size=10, color=color),
                name=category,
                showlegend=True,
            ))
    
    return fig