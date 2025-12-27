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