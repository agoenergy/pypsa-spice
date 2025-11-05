# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import os
import datetime as dt
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import streamlit as st
import re
from plotly.graph_objs._figure import Figure
from plotly.subplots import make_subplots
from typing import Dict, Optional, Callable
from itertools import cycle
from styles import use_flexo, apply_sidebar_chart_nav_styles, apply_radio_menu_styles

use_flexo()

pd.options.mode.chained_assignment = None


##TODO `Move to graph_settings.py`
def create_nice_names_and_color_mapping(
    tab_name: str,
) -> pd.DataFrame:
    """Get the names to hex codes mapping df for a given graph.

    Parameters
    ----------
    tab_name : str
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
        if re.search(pattern, tab_name):
            file_name = name
            break
    if file_name is None:
        return None

    file_path = os.path.join(st.session_state.current_dir, f"setting/{file_name}")
    df = pd.read_csv(file_path, index_col="original_names")

    return df


##TODO `Move to graph_settings.py`
def handle_color_mapping_for_chart(
    tab_name: str, legend_labels: Optional[list[str]] = None
) -> Dict[str, str]:
    """Get the colour mapping dict for a given graph.

    Parameters
    ----------
    tab_name : str
        Tab name of the graph as per the config file
    legend_labels : list, optional
        List of unique legends for the current graph

    Returns
    -------
    Dict[str, str]
        Dict of nice names to hex code mapping
    """
    df = create_nice_names_and_color_mapping(tab_name)

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


##TODO `Move to graph_settings.py`
def generate_default_colour_mapping(df: pd.DataFrame, leg_col: str) -> Dict[str, str]:
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
    prettified_legends = [_prettify_label(label) for label in unique_legends]

    default_colours = get_default_colour_list()

    mapping = dict(zip(prettified_legends, cycle(default_colours)))

    return mapping


##TODO `Move to graph_settings.py`
def get_default_colour_list() -> list:
    """Gets the list of default colours to use.

    Colours are assigned to legends that do not have a specified hex_code in
    tech_mapping or carrier_mapping. Currently using Agora EW default colours.

    Returns
    -------
    list
        The list of default colours
    """
    return ["#64B9E4", "#48A8AE", "#AD86B0", "#1E83B3", "#8393BE", "#637596"]


##TODO `Move to utils.py`
def read_result_csv(
    scenario_name: str,
    table_name: str,
    country: str = None,
    year: str = None,
) -> pd.DataFrame:
    """Read model ouput csv files for a given scenario and tab name.

    Parameters
    ----------
    scenario_name : str
        Selected scenario in streamlit UI
    table_name : str
        Output table name
    country : str, optional
        If not None, filter the csv by inputted country, by default None
    year : str, optional
        If not None, read csv from year specific folder else all_years folder,
        by default None

    Returns
    -------
    pd.DataFrame
        _description_
    """
    if year:
        file_path = os.path.abspath(
            st.session_state.result_path
            + "/"
            + scenario_name
            + "/csvs/"
            + st.session_state.sector
            + "/"
            + year
            + "/"
            + table_name
            + ".csv"
        )
    else:
        file_path = os.path.abspath(
            st.session_state.result_path
            + "/"
            + scenario_name
            + "/csvs/"
            + st.session_state.sector
            + "/"
            + "/all_years/"
            + table_name
            + ".csv"
        )

    try:
        df = pd.read_csv(os.path.abspath(file_path))
    except FileNotFoundError:
        with st.container(height=450, border=True):
            st.write(
                ":material/warning: File dose not exist or is empty: {}".format(
                    file_path
                )
            )
        return None
    if "country" in df.columns and country != None:
        df = df[df["country"] == country]

    df = df.fillna(0)

    return df


##TODO `Move to utils.py`
def slugify_text(text: str):
    """Helper function to slugify text string.

    This is used to generate safe anchor IDs (URL fragments) in the sidebar.

    Parameters
    ----------
    text : str
        The input text to slugify.

    Returns
    -------
    str
        The output (slugified) text.
    """
    text = text.lower()  # Lowercase text
    text = re.sub(r"[^a-z0-9]+", "-", text)  # Replace special characters with hyphens
    text = text.strip("-")  # Remove trailing/leading hyphens
    return text


def generate_sidebar(table_of_content):
    """generate sidebar with navigation links."""
    with st.sidebar:
        st.divider()

        apply_sidebar_chart_nav_styles()

        for section in table_of_content:
            anchor_id = slugify_text(section)
            st.markdown(
                f'<a href="#{anchor_id}" class="nav-link">{section}</a>',
                unsafe_allow_html=True,
            )


def setup_year_filter(config_plot: Dict, is_dual_scenario: bool) -> str:
    """Setup the year filter that appears in graphs with hourly data.

    Parameters
    ----------
    config_plot : Dict
        Configuration dictionary for the plot.
    is_dual_scenario : bool
        True if sce2 is selected, False otherwise.

    Returns
    -------
    str
        The selected year in the widget. Defaults to the first year.

    """
    # Set widget configuration params based on one or two scenarios
    if is_dual_scenario:
        years = sorted(
            list(set(st.session_state.sce1_years + st.session_state.sce2_years))
        )
        scenario_text = "both"
        key_prefix = "shared"
    else:
        years = st.session_state.sce1_years
        scenario_text = st.session_state.sce1
        key_prefix = "single"

    slider_id = config_plot["slider_id"].format(scenario_text)
    key = f"{key_prefix}_year_{config_plot['tab_name']}"
    label = f"{slider_id} Select Year:"

    # Pills widget for the year filter
    pills_widget = lambda: st.pills(
        label,
        options=years,
        key=key,
        default=years[0],
        label_visibility="collapsed",
    )

    if is_dual_scenario:
        spacer1, filter_col, spacer2 = st.columns([1, 3, 1])
        with filter_col:
            shared_year = pills_widget()
    else:
        shared_year = pills_widget()

    return shared_year


def setup_country_filter(config_plot, is_dual_scenario=False, scenario_tag=None) -> str:
    """Setup the country filter that appears in filtered_bar_hourly graphs.

    Parameters
    ----------
    config_plot : _type_
        Configuration dictionary for the plot
    is_dual_scenario : bool, optional
        True if sce2 is selected, by default False
    scenario_tag : _type_, optional
        Optional tag for the scenario, used in the label, by default None

    Returns
    -------
    str
        ountry selected in the widget. Defaults to the first country.
    """
    df = read_result_csv(scenario_tag, "pow_cap_by_type_yearly")

    # Setup country filters if needed
    if "country" in df.columns:
        # Set widget configuration params based on one or two scenarios
        if is_dual_scenario:
            country_options = sorted(list(set(df["country"].unique().tolist())))
            scenario_text = "both"
        else:
            country_options = df["country"].unique()
            scenario_text = st.session_state.sce1

        slider_id = config_plot["tab_name"]
        if "shared_country" in config_plot:
            country_id = config_plot["shared_country"]
        else:
            country_id = "all"
        key = f"shared_country_{country_id}_{scenario_tag}_{slider_id}"
        label = f"{slider_id} Select country:"

        # Pills widget for the country selection element
        shared_country = st.pills(
            label,
            options=country_options,
            key=key,
            default=country_options[0],
            label_visibility="collapsed",
        )
    else:
        shared_country = None

    return shared_country


def setup_region_filter(
    config_plot: Dict,
    df1: pd.DataFrame,
    df2: pd.DataFrame = None,
    is_dual_scenario: bool = False,
) -> Optional[str]:
    """Setup the region filter that appears in filtered_bar_hourly graphs.

    Parameters
    ----------
    config_plot : Dict
        Configuration dictionary for the plot.
    df1 : pd.DataFrame
        Data for sce1.
    df2 : pd.DataFrame, optional
        Data for sce2 (only if sce2 is selected). Defaults to None.
    is_dual_scenario : bool, optional
        True if sce2 is selected, False otherwise.

    Returns
    -------
    Optional[str]
        The selected region as a string, or None if there is no region filter.
    """
    if "fil_col" not in config_plot:
        return None
    fil_col = config_plot["fil_col"]
    # Set widget configuration params based on one or two scenarios
    if is_dual_scenario:
        region_options = sorted(
            list(set(df1[fil_col].unique().tolist() + df2[fil_col].unique().tolist()))
        )
        scenario_text = "both"
    else:
        region_options = df1[fil_col].unique()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    key = f"shared_region_{config_plot['tab_name']}"
    label = f"{slider_id} Select {fil_col}:"

    # Pills widget for the region selection element
    shared_region = st.pills(
        label,
        options=region_options,
        key=key,
        default=region_options[0],
        label_visibility="collapsed",
    )

    return shared_region


def setup_month_filter(
    config_plot: Dict,
    df1: pd.DataFrame,
    df2: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> int:
    """Setup month selection filter.

    This is called when data is complete (no discontinuous hours).

    Parameters
    ----------
    config_plot : Dict
        Configuration dictionary for the plot.
    df1 : pd.DataFrame
        Data for sce1.
    df2 : pd.DataFrame, optional
        Data for sce2 (only if sce2 is selected). Defaults to None.
    is_dual_scenario : bool, optional
        True if sce2 is selected, False otherwise. Defaults to False.

    Returns
    -------
    int
        The selected month as an integer (1-12).
    """
    # Set widget configuration params based on one or two scenarios
    if is_dual_scenario:
        months_sce1 = set(df1["snapshot"].dt.month.unique())
        months_sce2 = set(df2["snapshot"].dt.month.unique())
        months_all = sorted(list(months_sce1.union(months_sce2)))
        scenario_text = "both"
    else:
        months_all = df1["snapshot"].dt.month.unique()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    key = f"shared_month_{config_plot['tab_name']}"
    label = f"{slider_id} Select Month:"
    months_names = {m: _convert_month_to_name(m) for m in months_all}

    # Month selection widget
    selected_month = st.segmented_control(
        label,
        options=months_all,
        format_func=lambda x: months_names[x],
        key=key,
        default=months_all[0],
        label_visibility="collapsed",
    )

    return selected_month


def setup_date_filter_complete(
    config_plot: Dict,
    df1_m: pd.DataFrame,
    df2_m: pd.DataFrame | None = None,
    shared_year: str = None,
    selected_month: str = None,
    is_dual_scenario: bool = False,
) -> tuple[dt.datetime, dt.datetime]:
    """Setup date range slider.

    This is called when data is complete (no discontinuous hours).

    Parameters
    ----------
    config_plot : Dict
        Configuration dictionary for the plot.
    df1_m : pd.DataFrame
        Data for sce1.
    df2_m : pd.DataFrame, optional
        Data for sce2 (only if sce2 is selected). Defaults to None.
    shared_year : str, optional
        The year selected in the year selection filter. Defaults to None.
    selected_month : str, optional
        The month selected in the month selection filter. Defaults to None.
    is_dual_scenario : bool, optional
        True if sce2 is selected, False otherwise. Defaults to False.

    Returns
    -------
    tuple[dt.datetime, dt.datetime]
        A tuple containing the start and end timestamps selected in the slider.
    """
    # Set widget configuration params based on one or two scenarios
    if is_dual_scenario:
        min_date = min(df1_m["snapshot"].min(), df2_m["snapshot"].min())
        max_date = max(df1_m["snapshot"].max(), df2_m["snapshot"].max())
        scenario_text = "both"
    else:
        min_date = df1_m["snapshot"].min()
        max_date = df1_m["snapshot"].max()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    label = f"{slider_id} Select Date Range:"

    # Date range slider widget
    selected_dates = st.slider(
        label=label,
        min_value=min_date,
        max_value=max_date,
        value=(
            dt.datetime(int(shared_year), int(selected_month), 1, 0, 0),
            dt.datetime(int(shared_year), int(selected_month), 14, 0, 0),
        ),
        format="DD/MM/YY HH:mm",
        step=dt.timedelta(hours=1),
        label_visibility="collapsed",
    )

    return selected_dates


def setup_date_filter_incomplete(
    config_plot: Dict,
    df1_m: pd.DataFrame,
    df2_m: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Setup integer range slider.

    This is called when there are missing hours in the data.

    Parameters
    ----------
    config_plot : Dict
        Configuration dictionary for the plot.
    df1_m : pd.DataFrame
        Data for sce1.
    df2_m : pd.DataFrame, optional
        Data for sce2 (only if sce2 is selected). Defaults to None.
    is_dual_scenario : bool, optional
        True if sce2 is selected, False otherwise.

    Returns
    -------
    tuple[pd.Timestamp, pd.Timestamp]
        A tuple containing the start and end timestamps selected in the slider.
    """
    # Set widget configuration params, including getting the number of unique
    # timestamps present to use in the slider's range
    if is_dual_scenario:
        timestamps_1 = set(df1_m["snapshot"].unique())
        timestamps_2 = set(df2_m["snapshot"].unique())
        all_timestamps = sorted(list(timestamps_1.union(timestamps_2)))
        scenario_text = "both"
    else:
        all_timestamps = df1_m["snapshot"].unique()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    label = f"{slider_id} Select Range:"
    num_timestamps = len(all_timestamps)

    # Integer slider widget representing hours present in the dataset
    row_range = st.slider(
        label=label,
        min_value=1,
        max_value=num_timestamps,
        value=(1, min(20, num_timestamps)),
        step=1,
        label_visibility="collapsed",
    )

    # Convert selected indices to corresponding timestamps
    selected_dates = (
        all_timestamps[row_range[0] - 1],
        all_timestamps[row_range[1] - 1],
    )

    return selected_dates


def setup_radio_button_filter(
    config_plot: Dict, is_dual_scenario: bool
) -> Optional[str]:
    """Setup the radio button filter that appears in bar_with_filter graphs.

    Parameters
    ----------
    config_plot : Dict
        Configuration dictionary for the plot.
    is_dual_scenario : bool
        True if sce2 is selected, False otherwise.

    Returns
    -------
    Optional[str]
        The selected filter option, or None if only sce1 is selected.
    """
    if config_plot.get("graph_type") != "bar_with_filter":
        return None

    df1 = read_result_csv(st.session_state.sce1, config_plot["tab_name"])

    if is_dual_scenario:
        df2 = read_result_csv(st.session_state.sce2, config_plot["tab_name"])
        fil_col = config_plot["fil_col"]
        filter_options = sorted(
            list(set(df1[fil_col].unique().tolist() + df2[fil_col].unique().tolist()))
        )
        slider_id = config_plot["slider_id"].format("both")

        spacer1, filter_col, spacer2 = st.columns([1, 2, 1])
        with filter_col:
            shared_filter = st.radio(
                f"{slider_id} Select {fil_col} (both):",
                options=[str(x) for x in filter_options],
                format_func=_prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
    else:
        shared_filter = None

    return shared_filter


def setup_hourly_data_filters(
    df1: pd.DataFrame,
    df2: pd.DataFrame | None,
    config_plot: Dict,
    is_dual_scenario: bool,
) -> Dict:
    """Setup relevant filters for graphs with hourly data.

    This will setup month and date selection filters if data is complete, and an
    integer range filter if data is incomplete. Also adds the region selection filter if
    relevant.

    Parameters
    ----------
    df1: pd.DataFrame
        Data for sce1.
    df2: pd.DataFrame, optional
        Data for sce2 (only if sce2 is selected). Defaults to None.
    config_plot : Dict
        Configuration dictionary for the graph.
    is_dual_scenario : bool
        True if sce2 is selected, False otherwise.

    Returns
    -------
    Dict
        A dictionary containing the selected parameters (year, month, date, and region
        if applicable) from the filters.

    """
    shared_year = config_plot.get("shared_year", None)
    shared_region = config_plot.get("shared_region", None)
    if shared_region is not None:
        df1 = df1[df1[config_plot["fil_col"]] == shared_region]
        if df2 is not None:
            df2 = df2[df2[config_plot["fil_col"]] == shared_region]

    # Check if data is complete (full year)
    is_complete = all(len(df) % 8760 == 0 for df in [df1, df2] if df is not None)
    is_empty = all(df.empty for df in [df1, df2] if df is not None)

    if not is_empty:
        if is_complete:
            # Month selection filter
            selected_month = setup_month_filter(config_plot, df1, df2, is_dual_scenario)

            # Filter by month
            df1_m = _filter_by_month(df1, selected_month)
            df2_m = _filter_by_month(df2, selected_month) if df2 is not None else None

            # Date range filter
            selected_dates = setup_date_filter_complete(
                config_plot, df1_m, df2_m, shared_year, selected_month, is_dual_scenario
            )
        else:
            # In the case of incomplete data, do not show the month selector
            selected_month = None
            df1_m, df2_m = df1, df2

            # Integer range slider
            selected_dates = setup_date_filter_incomplete(
                config_plot, df1_m, df2_m, is_dual_scenario
            )
    else:
        shared_year = selected_month = selected_dates = None

    return {
        "shared_years": shared_year,
        "shared_months": selected_month,
        "shared_dates": selected_dates,
        "shared_region": shared_region if shared_region else None,
    }


def plot_indicator(graph_type, config_plot: dict):
    # graph_type: function (plot_function)
    st.markdown(
        f"<div id='{config_plot['name'].replace(' ','-')}'></div>",
        unsafe_allow_html=True,
    )
    st.markdown(f"#### {config_plot['name']}")

    yaxis_scales = calculate_yaxis_scales(config_plot.get("graph_type"), config_plot)
    config_plot["yaxis_scales"] = yaxis_scales

    # These are the graph types that have the year+month+date filter
    graphs_with_date_filters = [
        "simple_bar_hourly",
        "simple_line_hourly",
        "line_with_secondary_y_hourly",
        "filtered_bar_hourly",
    ]

    # Track whether sce2 has been selected by the user or not
    is_dual_scenario = st.session_state.sce2 and st.session_state.sce2 != ""

    if is_dual_scenario:
        # Optional, but we can use this CSS to centre the radio buttons within the
        # column that will contain the filter menu
        apply_radio_menu_styles()

    # Setup country filter
    config_plot["shared_country"] = setup_country_filter(
        config_plot, is_dual_scenario, scenario_tag=st.session_state.sce1
    )

    # Setup filters based on graph type
    if config_plot.get("graph_type") == "bar_with_filter":
        # Two scenarios: plots with radio button filters
        shared_filter = setup_radio_button_filter(config_plot, is_dual_scenario)
        if shared_filter:
            config_plot["shared_filter"] = shared_filter

    elif config_plot.get("graph_type") in graphs_with_date_filters:
        # Setup year filter
        shared_year = setup_year_filter(config_plot, is_dual_scenario)
        config_plot["shared_year"] = str(shared_year)

        df1 = read_result_csv(
            st.session_state.sce1, config_plot["tab_name"], year=str(shared_year)
        )
        df1["snapshot"] = pd.to_datetime(df1["snapshot"])

        df2 = [
            (
                read_result_csv(
                    st.session_state.sce2,
                    config_plot["tab_name"],
                    year=str(shared_year),
                )
                if is_dual_scenario
                else None
            )
        ][0]

        # Setup region filter if needed
        config_plot["shared_region"] = setup_region_filter(
            config_plot, df1, df2, is_dual_scenario
        )

        if is_dual_scenario:
            # Two scenarios: plots with hourly data/year+month+date filters
            spacer1, filter_col, spacer2 = st.columns([1, 2, 1])
            with filter_col:
                df2["snapshot"] = pd.to_datetime(df2["snapshot"])

        # Setup all hourly data filters
        filter_results = setup_hourly_data_filters(
            df1, df2, config_plot, is_dual_scenario
        )
        config_plot.update(filter_results)

    # Render the plots
    if not is_dual_scenario:
        # Display the graph for single scenario
        config_plot["years"] = st.session_state.sce1_years
        st.markdown(f"#### {st.session_state.sce1} ")
        graph_type(sc_name=st.session_state.sce1, config_g=config_plot)

        # Display the data download part
        if config_plot.get("graph_type") in graphs_with_date_filters:
            display_data_download_button(st.session_state.sce1, config_plot)
        else:
            display_data_table(st.session_state.sce1, config_plot)

    else:
        # Display the graphs for each of the two scenarios
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            config_plot["years"] = st.session_state.sce1_years
            st.markdown(f"#### {st.session_state.sce1} ")
            graph_type(sc_name=st.session_state.sce1, config_g=config_plot)
        with col3:
            config_plot["years"] = st.session_state.sce2_years
            st.markdown(f"#### {st.session_state.sce2} ")
            graph_type(sc_name=st.session_state.sce2, config_g=config_plot)

        # Display the data download part
        col1, col2, col3 = st.columns([6, 1, 6])

        if config_plot.get("graph_type") in graphs_with_date_filters:
            with col1:
                display_data_download_button(st.session_state.sce1, config_plot)
            with col3:
                display_data_download_button(st.session_state.sce2, config_plot)
        else:
            with col1:
                display_data_table(st.session_state.sce1, config_plot)
            with col3:
                display_data_table(st.session_state.sce2, config_plot)

    st.divider()


def create_download_csv_button(csv_data, download_id):
    st.download_button(
        label="Download CSV",
        key=f"download_button_{download_id}",
        data=csv_data,
        file_name=f"{download_id}.csv",
        mime="text/csv",
    )


def display_data_download_button(sc_name: str, config_g: Dict[str, str]):
    """Generate and display the CSV download button for a graph.

    Parameters
    ----------
    sc_name : str
        The scenario name for the data that will be downloaded.
    config_g : Dict[str, str]
        The configuration dictionary for the graph.
    """
    leg_col = config_g["leg_col"]
    download_id = config_g["download_id"].format(sc_name)
    graph_type = config_g["graph_type"]

    df = read_result_csv(
        sc_name,
        config_g["tab_name"],
        year=str(config_g["shared_years"]),
        country=config_g["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, _ = _get_filtered_data(df, config_g)

        df_sel = _filter_between_dates(df_m, start_date=start_date, end_date=end_date)

        csv_data = (
            df_sel.to_csv().encode("utf-8")
            if graph_type == "filtered_bar_hourly"
            else df_sel[["snapshot", leg_col, "value"]].to_csv().encode("utf-8")
        )

        create_download_csv_button(csv_data, download_id)


def display_data_table(sc_name: str, config_g: Dict[str, str]):
    """Generate and display the data table and CSV download button for a graph.

    Parameters
    ----------
    sc_name : str
        The scenario name for the data being displayed in the table.
    config_g : Dict[str, str]
        The configuration dictionary for the graph.
    """
    leg_col = config_g["leg_col"]
    download_id = config_g["download_id"].format(sc_name)
    df = read_result_csv(
        sc_name, config_g["tab_name"], country=config_g["shared_country"]
    )
    base_group_cols = {"year", leg_col}

    # Add additional columns that may be present (for dfs that have additional columns
    # for filtering)
    additional_group_cols = [
        col
        for col in df.columns
        if col not in base_group_cols.union({"value"}) and df[col].dtype == "object"
    ]

    group_cols = ["year", leg_col] + additional_group_cols

    df_grouped = df.groupby(group_cols, as_index=False)["value"].sum(min_count=1)

    with st.expander(f":material/database: Data ({sc_name}):", expanded=False):
        # For bar_with_filter graphs with from/to columns, make sure 'from' is the first
        # column in the pivoted table
        pivot_index = (
            [leg_col] + additional_group_cols
            if leg_col == "from"
            else additional_group_cols + [leg_col]
        )
        df_pivot = pd.pivot_table(
            df_grouped,
            values="value",
            columns="year",
            index=pivot_index,
        )
        df_pivot = df_pivot.loc[~(df_pivot == 0).all(axis=1)].fillna(
            0
        )  # Drop all-0 rows
        # Add column names for index columns
        # Note that this only adds them to the downloaded csv, not the displayed table
        # in the UI, as st.table does not support displaying index level names
        df_pivot.index.names = pivot_index
        styled_df = df_pivot.style.apply(_highlight_diff, axis=None).format("{:.2f}")
        st.table(styled_df)

        csv_data = df_pivot.to_csv().encode("utf-8")

        create_download_csv_button(csv_data, download_id)


def _clean_df_for_plotting(leg_col: str, df: pd.DataFrame):
    """Clean the data used to plot the graph.

    1. Filter out rows from the raw data df where all values of that legend series are
    zero or NaN, in order to make them not appear in the graph.
    2. Convert all values <e-06 to 0.0

    Parameters
    ----------
    leg_col : str
        The legend column name as per the graph's configuration dictionary.
    df : pd.DataFrame
        The raw data returned by read_result_csv.

    Returns
    -------
    pd.DataFrame
        The input dataframe with zero and NaN rows removed and small values converted to
        0.0
    """
    df_pivoted = df.pivot_table(
        values="value",
        columns="year" if "year" in df.columns else "snapshot",
        index=leg_col,
        aggfunc="mean",
        dropna=False,
    )

    # Identify legends where all values are 0 across the row
    all_legends_to_remove = df_pivoted[
        (df_pivoted.isna() | (df_pivoted == 0)).all(axis=1)
    ].index

    # Exclude legends where all values are 0 from the original df
    df_filtered = df[~df[leg_col].isin(all_legends_to_remove)]

    df_filtered = _handle_small_values(df_filtered)

    return df_filtered


def _handle_small_values(df: pd.DataFrame) -> pd.DataFrame:
    """Converts <1e-6 values in the "value" column to 0.0.

    Parameters
    ----------
    df : pd.DataFrame
        The input dataframe containing a 'value' column.

    Returns
    -------
    pd.DataFrame
        The output dataframe.
    """
    if "value" in df.columns:
        df.loc[df["value"].abs() < 1e-6, "value"] = 0.0
    return df


def _highlight_diff(data):
    """Highlight differences between adjacent cells in a dataframe.

    Green/red arrows are shown for an increase/decrease in values. Colours are scaled
    relative to the first column's value in each row.

    Parameters
    ----------
    data : pd.DataFrame
        The dataframe containing values to be styled.
    """
    styled = pd.DataFrame("", index=data.index, columns=data.columns)
    numeric_data = data.apply(pd.to_numeric, errors="coerce")

    # Base padding in all cells to push numbers slightly to the left in order to create
    # space to accommodate arrows
    for idx in range(len(numeric_data)):
        for col in range(len(numeric_data.columns)):
            styled.iloc[
                idx, col
            ] = """
                padding-right: 20px;
                text-align: right;
            """

    for idx in range(len(numeric_data)):
        row_data = numeric_data.iloc[idx]
        base_value = row_data.iloc[0]  # Get the first column's value as reference

        # Use a small value if base value is 0 to avoid division by zero error
        safe_base_value = 0.01 if base_value == 0 else base_value

        for col in range(1, len(numeric_data.columns)):

            current_value = row_data.iloc[col]

            # If the current value is equal to the base value, skip calculating the
            # percentage change since there is no difference to visualise
            if (
                pd.isna(base_value)
                or pd.isna(current_value)
                or base_value == current_value
            ):
                continue

            # Calculate percentage change of the current value from the base value
            pct_change = (current_value - base_value) / abs(safe_base_value)

            # Calculate colour intensity based on change, ranges from 20% to 100%
            intensity = min(abs(pct_change), 1.0) * 0.8 + 0.2

            if pct_change > 0:
                # Up arrow from two triangles
                styled.iloc[
                    idx, col
                ] = f"""
                    background-color: rgba(100, 100, 100, 0.15);
                    position: relative;
                    padding-right: 20px;
                    background-image: 
                        linear-gradient(135deg, transparent 50%, rgb(40 167 69 / {intensity}) 30%),
                        linear-gradient(-135deg, transparent 50%, rgb(40 167 69 / {intensity}) 30%);
                    background-size: 6px 6px, 6px 6px;
                    background-position: 
                        right 8px top 2px,
                        right 2px top 2px;
                    background-repeat: no-repeat;
                """
            elif pct_change < 0:
                # Down arrow from two triangles
                styled.iloc[
                    idx, col
                ] = f"""
                    background-color: rgba(100, 100, 100, 0.15);
                    position: relative;
                    padding-right: 20px;
                    background-image: 
                        linear-gradient(45deg, transparent 50%, rgb(240 30 55 / {intensity}) 30%),
                        linear-gradient(-45deg, transparent 50%, rgb(240 30 55 / {intensity}) 30%);
                    background-size: 6px 6px, 6px 6px;
                    background-position: 
                        right 8px top 2px,
                        right 2px top 2px;
                    background-repeat: no-repeat;
                """

    return styled


def simple_bar_yearly(sc_name: str, config_g: Dict[str, str]) -> Optional[None]:
    """Generate a yearly stacked bar chart with a downloadable data table."""
    leg_col = config_g["leg_col"]
    download_id = config_g["download_id"].format(sc_name)
    tab_name = config_g["tab_name"]

    # Construct file path
    df = read_result_csv(
        sc_name, config_g["tab_name"], country=config_g["shared_country"]
    )
    df = _clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(tab_name)

    df_grouped = df.groupby(["year", leg_col], as_index=False)["value"].sum()
    df_grouped["nice_names"] = df_grouped[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else _prettify_label(x)
        )
    )

    unique_legends = df_grouped["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(tab_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    df1_grouped = df_grouped
    df2_grouped = None
    if config_g.get("shared_country") is not None and st.session_state.sce2:
        df1_grouped, df2_grouped = _get_yearly_dfs(
            config_g,
            lambda df: df.groupby(["year", leg_col], as_index=False)["value"].sum(),
        )

    y_range = _calculate_min_max_y_scale(df1_grouped, df2_grouped, "year")
    min_y = y_range["min"]
    max_y = y_range["max"]

    # Generate and display the chart
    try:
        fig = px.bar(
            df_grouped,
            x="year",
            y="value",
            color="nice_names",
            barmode="group" if tab_name == "pow_bats_ep_ratio" else "stack",
            color_discrete_map=colour_mapping,
        )

        # Add stacked bar total and update layout
        fig = _add_stackedbar_total(fig, df_grouped)
        fig = _update_layout(
            fig, df_grouped, {"max_scale": max_y, "min_scale": min_y}, config_g
        )

        # Display the chart with a unique key
        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{download_id}"
        )

    except ValueError as e:
        st.error(f"ValueError encountered: {e}")
        st.dataframe(df)


@st.fragment
def simple_line_yearly(sc_name: str, config_g):
    """_summary_

    Args:
        sc_name (str): _description_
        tab_name (str): _description_
    """
    # Read data from file
    tab_name = config_g["tab_name"]
    leg_col = config_g["leg_col"]
    download_id = config_g["download_id"].format(sc_name)

    df = read_result_csv(sc_name, tab_name, country=config_g["shared_country"])

    df = _clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(tab_name)

    df["nice_names"] = df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else _prettify_label(x)
        )
    )

    unique_legends = df["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(tab_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    df1 = df
    df2 = None
    if config_g.get("shared_country") is not None and st.session_state.sce2:
        df1, df2 = _get_yearly_dfs(config_g, None)

    y_range = _calculate_min_max_y_scale(df1, df2, None)
    min_y = y_range["min"]
    max_y = y_range["max"]

    fig = px.line(
        df, x="year", y="value", color="nice_names", color_discrete_map=colour_mapping
    )

    fig = _update_layout(fig, df, {"max_scale": max_y, "min_scale": min_y}, config_g)
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{sc_name}_{tab_name}"
    )


@st.fragment
def bar_with_filter(sc_name: str, config_g: dict):
    leg_col = config_g["leg_col"]
    fil_col = config_g["fil_col"]
    slider_id = config_g["slider_id"].format(sc_name)
    download_id = config_g["download_id"].format(sc_name)
    tab_name = config_g["tab_name"]
    df = read_result_csv(sc_name, tab_name, country=config_g["shared_country"])

    mapping_df = create_nice_names_and_color_mapping(tab_name)

    # This creates a shared filter that applies to both graphs if the entry is found in
    # the config dict, otherwise a local filter is generated for the single graph.
    if "shared_filter" in config_g:
        filter = config_g["shared_filter"]
    else:
        filter = st.radio(
            "{} Select {} ({})".format(slider_id, fil_col, sc_name) + ":",
            options=df[fil_col].unique(),
            format_func=_prettify_label,
            horizontal=True,
            label_visibility="collapsed",
        )

    df_reg = df.loc[df[fil_col] == filter]
    df_reg["nice_names"] = df_reg[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else _prettify_label(x)
        )
    )

    unique_legends = df_reg["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(tab_name, unique_legends)
        if mapping_df is not None
        else generate_default_colour_mapping(df, leg_col)
    )

    # If scenario2 exists, recalculate the maximum y tick value considering both
    # datasets, else calculate the maximum y for just the single dataset
    df1_reg = df_reg
    df2_reg = None
    if config_g.get("shared_filter") is not None and st.session_state.sce2:
        df1_reg, df2_reg = _get_yearly_dfs(
            config_g, lambda df: df.loc[df[fil_col] == filter]
        )

    y_range = _calculate_min_max_y_scale(df1_reg, df2_reg, "year")
    min_y = y_range["min"]
    max_y = y_range["max"]

    df_reg = _clean_df_for_plotting(leg_col, df_reg)

    fig = px.bar(
        df_reg,
        x="year",
        y="value",
        color="nice_names",
        color_discrete_map=colour_mapping,
    )

    fig = _add_stackedbar_total(fig, df_reg)
    fig = _update_layout(
        fig, df_reg, {"max_scale": max_y, "min_scale": min_y}, config_g
    )
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{sc_name}_{tab_name}"
    )


@st.fragment
def area_share_yearly(sc_name: str, config_g):
    """_summary_

    Args:
        sc_name (str): _description_
        tab_name (str): _description_
    """
    download_id = config_g["download_id"].format(sc_name)
    leg_col = config_g["leg_col"]
    tab_name = config_g["tab_name"]

    df = read_result_csv(sc_name, tab_name, country=config_g["shared_country"])
    df = _clean_df_for_plotting(leg_col, df)

    mapping_df = create_nice_names_and_color_mapping(tab_name)

    df["nice_names"] = df[leg_col].map(
        lambda x: (
            mapping_df.loc[x, "nice_names"]
            if (mapping_df is not None and x in mapping_df.index)
            else _prettify_label(x)
        )
    )
    unique_legends = df["nice_names"].unique().tolist()
    colour_mapping = (
        handle_color_mapping_for_chart(tab_name, unique_legends)
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

    fig = _update_layout(fig, df, None, config_g)
    st.plotly_chart(
        fig, use_container_width=True, key=f"plotly_chart_{sc_name}_{tab_name}"
    )


@st.fragment
def simple_bar_hourly(sc_name: str, config_g: Dict[str, str]) -> Optional[None]:
    """Generate hourly stacked bar chart with filters for datetime."""
    tab_name = config_g["tab_name"]
    leg_col = config_g["leg_col"]
    download_id = config_g["download_id"].format(sc_name)

    df = read_result_csv(
        sc_name,
        config_g["tab_name"],
        year=str(config_g["shared_years"]),
        country=config_g["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = _get_filtered_data(df, config_g)

        # Validate date range
        if start_date > end_date:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        mapping_df = create_nice_names_and_color_mapping(tab_name)

        # Filter data by date range
        df_sel = _filter_between_dates(df_m, start_date=start_date, end_date=end_date)
        df_sel["nice_names"] = df_sel[leg_col].map(
            lambda x: (
                mapping_df.loc[x, "nice_names"]
                if (mapping_df is not None and x in mapping_df.index)
                else _prettify_label(x)
            )
        )
        unique_legends = df_sel["nice_names"].unique().tolist()
        colour_mapping = (
            handle_color_mapping_for_chart(tab_name, unique_legends)
            if mapping_df is not None
            else generate_default_colour_mapping(df_sel, leg_col)
        )

        df_sel = _handle_small_values(df_sel)

        fig = px.bar(
            df_sel,
            x="snapshot",
            y="value",
            color="nice_names",
            color_discrete_map=colour_mapping,
        )

        fig = _update_hourly_plot_x_axis(fig, df_sel, start_date, end_date, is_complete)

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = _get_hourly_filtered_dfs(config_g)
            y_range = _calculate_min_max_y_scale(df1_sel, df2_sel, "snapshot")
        else:
            y_range = _calculate_min_max_y_scale(df_sel, None, "snapshot")

        fig = _update_layout(
            fig,
            df_sel,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            config_g,
        )
        # Display the chart with a unique key
        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{download_id}"
        )


@st.fragment
def simple_line_hourly(sc_name: str, config_g: dict):
    """_summary_

    Args:
        sc_name (str): _description_
        tab_name (str): _description_
        slider_id (str): _description_
    """
    tab_name = config_g["tab_name"]
    leg_col = config_g["leg_col"]
    download_id = config_g["download_id"].format(sc_name)

    df = read_result_csv(
        sc_name,
        config_g["tab_name"],
        year=str(config_g["shared_years"]),
        country=config_g["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = _get_filtered_data(df, config_g)

        mapping_df = create_nice_names_and_color_mapping(tab_name)

        if start_date <= end_date:
            pass
        else:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        df_sel = _filter_between_dates(df_m, start_date=start_date, end_date=end_date)
        df_sel["nice_names"] = df_sel[leg_col].map(
            lambda x: (
                mapping_df.loc[x, "nice_names"]
                if (mapping_df is not None and x in mapping_df.index)
                else _prettify_label(x)
            )
        )
        unique_legends = df_sel["nice_names"].unique().tolist()
        colour_mapping = (
            handle_color_mapping_for_chart(tab_name, unique_legends)
            if mapping_df is not None
            else generate_default_colour_mapping(df_sel, leg_col)
        )

        df_sel = _handle_small_values(df_sel)

        fig = px.line(
            df_sel,
            x="snapshot",
            y="value",
            color="nice_names",
            color_discrete_map=colour_mapping,
        )

        fig = _update_hourly_plot_x_axis(fig, df_sel, start_date, end_date, is_complete)

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = _get_hourly_filtered_dfs(config_g)
            y_range = _calculate_min_max_y_scale(df1_sel, df2_sel, None)
        else:
            y_range = _calculate_min_max_y_scale(df_sel, None, None)

        fig = _update_layout(
            fig,
            df_sel,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            config_g,
        )
        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{sc_name}_{tab_name}"
        )


@st.fragment
def filtered_bar_hourly(sc_name: str, config_g: dict):
    """_summary_

    Args:
        sc_name (str): _description_
        tab_name (str): _description_
        slider_id (str): _description_
    """
    tab_name = config_g["tab_name"]
    leg_col = config_g["leg_col"]
    fil_col = config_g["fil_col"]

    df = read_result_csv(
        sc_name,
        config_g["tab_name"],
        year=str(config_g["shared_years"]),
        country=config_g["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = _get_filtered_data(df, config_g)

        mapping_df = create_nice_names_and_color_mapping(tab_name)

        if start_date <= end_date:
            pass
        else:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        df_sel = _filter_between_dates(df_m, start_date=start_date, end_date=end_date)
        df_sel["nice_names"] = df_sel[leg_col].map(
            lambda x: (
                mapping_df.loc[x, "nice_names"]
                if (mapping_df is not None and x in mapping_df.index)
                else _prettify_label(x)
            )
        )

        df_sel = _handle_small_values(df_sel)

        unique_legends = df_sel["nice_names"].unique().tolist()
        colour_mapping = (
            handle_color_mapping_for_chart(tab_name, unique_legends)
            if mapping_df is not None
            else generate_default_colour_mapping(df_sel, leg_col)
        )

        fig = px.bar(
            df_sel[df_sel["value"] != 0],
            x="snapshot",
            y="value",
            color="nice_names",
            color_discrete_map=colour_mapping,
        )

        fig = _update_hourly_plot_x_axis(fig, df_sel, start_date, end_date, is_complete)

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = _get_hourly_filtered_dfs(config_g)
            y_range = _calculate_min_max_y_scale(df1_sel, df2_sel, "snapshot")
        else:
            y_range = _calculate_min_max_y_scale(df_sel, None, "snapshot")

        fig = _update_layout(
            fig,
            df_sel,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            config_g,
        )

        if fil_col in df_sel.columns:
            df_line = df_sel.groupby(["snapshot", fil_col])["value"].sum().reset_index()
            line_chart_trace = px.line(df_line, x="snapshot", y="value")
            line_chart_trace.update_traces(line=dict(color="blue", width=3))

            for trace in line_chart_trace.data:
                fig.add_trace(trace)

        df_line = df_sel.groupby(["snapshot", fil_col])["value"].sum().reset_index()
        line_chart_trace = px.line(df_line, x="snapshot", y="value")
        line_chart_trace.update_traces(line=dict(color="blue", width=3))

        for trace in line_chart_trace.data:
            fig.add_trace(trace)

        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{sc_name}_{tab_name}"
        )


@st.fragment
def line_with_secondary_y_hourly(sc_name: str, config_g: dict):
    """_summary_

    Args:
        sc_name (str): _description_
        tab_name (str): _description_
        slider_id (str): _description_
    """
    tab_name = config_g["tab_name"]
    leg_col = config_g["leg_col"]
    primary_y_lab = config_g["primary_y_lab"]
    secondary_y_lab = config_g["secondary_y_lab"]

    df = read_result_csv(
        sc_name,
        config_g["tab_name"],
        year=str(config_g["shared_years"]),
        country=config_g["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, is_complete = _get_filtered_data(df, config_g)

        mapping_df = create_nice_names_and_color_mapping(tab_name)
        colour_mapping = (
            handle_color_mapping_for_chart(tab_name)
            if mapping_df is not None
            else generate_default_colour_mapping(df_m, leg_col)
        )

        if start_date <= end_date:
            pass
        else:
            st.error("Error: End date must be greater than or equal to start date.")
            return

        df_sel = _filter_between_dates(df_m, start_date=start_date, end_date=end_date)

        # Create figure with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        for prim_y in primary_y_lab:
            nice_name = (
                mapping_df.loc[prim_y, "nice_names"]
                if (mapping_df is not None and prim_y in mapping_df.index)
                else _prettify_label(prim_y)
            )
            fig.add_trace(
                go.Line(
                    y=df_sel[df_sel[leg_col] == prim_y]["value"],
                    x=df_sel["snapshot"],
                    name=nice_name,
                    line=dict(color=colour_mapping.get(nice_name, "#a0a0a0")),
                ),
                secondary_y=False,
            )

        for secd_y in secondary_y_lab:
            nice_name = (
                mapping_df.loc[secd_y, "nice_names"]
                if (mapping_df is not None and secd_y in mapping_df.index)
                else _prettify_label(secd_y)
            )
            fig.add_trace(
                go.Line(
                    y=df_sel[df_sel[leg_col] == secd_y]["value"],
                    x=df_sel["snapshot"],
                    name=nice_name,
                    line=dict(color=colour_mapping.get(nice_name, "#a0a0a0")),
                ),
                secondary_y=True,
            )

        fig = _update_hourly_plot_x_axis(fig, df_sel, start_date, end_date, is_complete)

        if st.session_state.sce2 and st.session_state.sce2 != "":
            df1_sel, df2_sel = _get_hourly_filtered_dfs(config_g)
            y_range = _calculate_min_max_y_scale(df1_sel, df2_sel, None)
        else:
            y_range = _calculate_min_max_y_scale(df_sel, None, None)

        fig = _update_layout(
            fig,
            df_sel,
            {"max_scale": y_range["max"], "min_scale": y_range["min"]},
            config_g,
        )
        fig.update_yaxes(
            title_text=_handle_y_axis_list(primary_y_lab), secondary_y=False
        )
        fig.update_yaxes(
            title_text=_handle_y_axis_list(secondary_y_lab), secondary_y=True
        )
        st.plotly_chart(
            fig, use_container_width=True, key=f"plotly_chart_{sc_name}_{tab_name}"
        )


def _prettify_label(label: str) -> str:
    """
    Converts snake_case or camelCase string into readable text version for legends and
    filters.

    Parameters
    ----------
    label : str
        The input string to format.

    Returns
    -------
    str
        The formatted, human-readable string.
    """

    # Handle snake_case labels
    if "_" in label:
        if "_to_" in label:
            # Specific case of 'flow' strings
            parts = label.split("_")
            if len(parts) == 5 and parts[2] == "to":
                # Format XX_YY_to_AA_BB -> XX (YY) to AA (BB)
                return f"{parts[0]} ({parts[1]}) to {parts[3]} ({parts[4]})"
            # Fallback in case the string is malformed but still has _to_
            return " ".join(parts)
        else:
            return " ".join(label.split("_"))

    # Handle camelCase labels
    elif re.search(r"[a-z][A-Z]", label):
        spaced = re.sub(r"([a-z])([A-Z])", r"\1 \2", label)
        return spaced.capitalize()

    return label


def _handle_y_axis_list(title_list: list[str]) -> str:
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


def _filter_by_month(df: pd.DataFrame, month: int):
    """Filter input dataframe by month."""
    # Convert the 'Date' column to datetime format
    df["Date"] = pd.to_datetime(df["snapshot"])

    # Extract month from the 'Date' column
    df["Month"] = df["Date"].dt.month

    # Filter data by the selected month
    filtered_df = df[df["Month"] == month]

    return filtered_df


def _filter_by_year(df: pd.DataFrame, year: int):
    """Filter input dataframe by year."""
    # Convert the 'Date' column to datetime format
    df["Date"] = pd.to_datetime(df["snapshot"])

    # Extract month from the 'Date' column
    df["Year"] = df["Date"].dt.year

    # Filter data by the selected month
    filtered_df = df[df["Year"] == year]

    return filtered_df


def _filter_between_dates(
    df: pd.DataFrame, start_date: dt.datetime, end_date: dt.datetime
):
    """Filter input dataframe by specific date range."""
    # Convert the 'Date' column to datetime format
    df["Date"] = pd.to_datetime(df["snapshot"])

    # Filter data between the selected dates
    filtered_df = df[(df["Date"] >= start_date) & (df["Date"] <= end_date)]

    return filtered_df


def _get_filtered_data(df: pd.DataFrame, config_g: Dict[str, str]):
    """Get the filtered data, start date, and end date for graphs with date filters.
    Relevant graphs are:
    - simple_bar_hourly
    - simple_line_hourly
    - line_with_secondary_y_hourly
    - filtered_bar_hourly

    Args:
        df (DataFrame): DataFrame of the scenario
        config_g (Dict[str, str]): Configuration dictionary that may contain:
            - 'shared_years': int, shared year for both scenarios
            - 'shared_months': int, shared month for both scenarios
            - 'shared_dates': tuple of datetime objects, shared date range for both
            scenarios

    Returns:
        df_m: The data filtered by the selected month
        start_date: The start date selected by the user
        end_date: The end date selected by the user
        is_complete (bool): Whether all hours of the year are present

    """
    month = config_g["shared_months"]
    start_date, end_date = config_g["shared_dates"]

    # Check if data is complete (no discontinuous hours)
    is_complete = len(df) % 8760 == 0

    df["snapshot"] = pd.to_datetime(df["snapshot"])

    # Add region filtering to the config dict for filtered_bar_hourly
    if "shared_region" in config_g and "fil_col" in config_g:
        fil_col = config_g["fil_col"]
        df = df[df[fil_col] == config_g["shared_region"]]

    # Only filter by month if month selector exists
    df_m = _filter_by_month(df=df, month=month) if month is not None else df

    return df_m, start_date, end_date, is_complete


def _convert_month_to_name(month_num: int) -> str:
    """Helper function to convert a month number to the abbreviated month name."""
    return dt.datetime.strptime(str(month_num), "%m").strftime("%b")


def _add_stackedbar_total(fig: Figure, df: pd.DataFrame) -> Figure:
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


def _update_layout(
    fig: Figure, df: pd.DataFrame, yaxis_scales: dict = None, config_g: dict = None
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
    config_g : dict, optional
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
    if config_g and "units" in config_g:
        units = f" [{config_g['units']}]"

    # Handle y-axis scaling
    if yaxis_scales is None and config_g and "yaxis_scales" in config_g:
        yaxis_scales = config_g["yaxis_scales"]

    layout_dict = {
        "showlegend": True,
        "font": dict(family="Flexo, sans-serif", size=15),
        "legend": dict(
            orientation=legend_orientation,
            y=legend_y_pos,
            x=legend_x_pos,
            xanchor=legend_x_anchor,
            yanchor=legend_y_anchor,
            title_text=config_g["leg_col"].capitalize(),
            title_font_size=16,
        ),
        "margin": dict(t=50, b=margin_b),
        "xaxis": dict(tickfont=dict(size=x_tick_font_size)),
        "yaxis": dict(tickfont=dict(size=15)),
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
            )
        ],
    }

    # For the yearly bar charts, adjust the space between bars
    if config_g["graph_type"] in [
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
            dict(tickvals=list(years_with_data), tickmode="array")
        )

    fig.update_layout(**layout_dict)

    return fig


def _update_hourly_plot_x_axis(
    fig: Figure,
    df_sel: pd.DataFrame,
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
    df_sel : pd.DataFrame
        Dataframe containing the filtered data
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
            tickfont=dict(size=11),
        )
    else:
        unique_snapshots = df_sel["snapshot"].unique()
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
            tickfont=dict(size=13),
        )

    return fig


def plot_function(func_name: str = None):
    mapping = {
        "simple_bar_yearly": simple_bar_yearly,
        "simple_bar_yearly_2": simple_bar_yearly,
        "simple_line_yearly": simple_line_yearly,
        "bar_with_filter": bar_with_filter,
        "area_share_yearly": area_share_yearly,
        "simple_bar_hourly": simple_bar_hourly,
        "simple_line_hourly": simple_line_hourly,
        "filtered_bar_hourly": filtered_bar_hourly,
        "line_with_secondary_y_hourly": line_with_secondary_y_hourly,
    }
    func = mapping.get(func_name)
    if func is None:
        st.write(f"Function not mapped: {func_name}")
    return func


def _grouping_data(
    sc_name: str, tab_name: str, years: list = None, hourly: bool = False
) -> pd.DataFrame:
    """Helper function to read and concatenate data for hourly scenarios

    Parameters
    ----------
    sc_name : str
        scenario name
    tab_name : str
        tabe name
    years : list
        a list of the years

    Returns
    -------
    DataFrame
        DataFrame containing the concatenated data
    """
    if hourly and years:
        df_combined = pd.concat(
            [
                read_result_csv(
                    scenario_name=sc_name, table_name=tab_name, year=str(year)
                )
                for year in years
            ],
            ignore_index=True,
        )
    else:
        df_combined = read_result_csv(scenario_name=sc_name, table_name=tab_name)
    return df_combined


def _get_yearly_dfs(
    config_g: Dict, fn: Optional[Callable[[pd.DataFrame], pd.DataFrame]] = None
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the yearly dfs for two scenarios, and optionally process (before min/max y
    calculation takes place).

    Parameters
    ----------
    config_g : Dict
        Configuration dictionary for the current graph
    fn: Callable
        Optional function to further process each dataframe (e.g., groupby or filter)

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        The dataframes (optionally processed) for the two scenarios
    """
    dfs = []
    tab_name = config_g["tab_name"]
    shared_country = config_g["shared_country"]

    for scenario in [st.session_state.sce1, st.session_state.sce2]:
        df = read_result_csv(scenario, tab_name, shared_country)
        if fn:
            df = fn(df)
        dfs.append(df)

    return dfs


def _get_hourly_filtered_dfs(config_g: Dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the filtered hourly dataframes for two scenarios

    Parameters
    ----------
    config_g : Dict
        Configuration dictionary for the current graph

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        The filtered hourly dataframes for the two scenarios
    """
    df_sels = []
    for i, scenario in enumerate([st.session_state.sce1, st.session_state.sce2]):
        df = read_result_csv(
            scenario,
            config_g["tab_name"],
            year=str(config_g["shared_years"]),
            country=config_g["shared_country"],
        )

        if df is not None and not df.empty:
            df_m, s_date, e_date, _ = _get_filtered_data(df, config_g)
            if i == 0:
                start_date, end_date = s_date, e_date
            df_sel = _filter_between_dates(
                df_m, start_date=start_date, end_date=end_date
            )
            df_sel = _handle_small_values(df_sel)
            df_sels.append(df_sel)

    return df_sels


def _calculate_min_max_y_scale(
    df: pd.DataFrame, df2: pd.DataFrame, group_col: str = None
) -> dict:
    """Calculate manimum and maximum values to use on the y axis.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame. Assumes it contains a "value" column.
    df2 : pd.DataFrame
        Optional second input DataFrame to compare with to calculate a joint
    group_col : str, optional
        Column in the df to group by first, by default None

    Returns
    -------
    dict
        Dictionary with the following keys:
        - "min" : float
            Minimum value for the y-axis, scaled by 1.2 if negative
        - "max" : float
            Maximum value for the y-axis, scaled by 1.2
    """

    def compute_min_max(data: pd.DataFrame) -> tuple[float, float]:
        if group_col and group_col in data.columns:

            grouped = data.groupby(group_col)

            positive_sum = grouped["value"].apply(lambda x: x[x > 0].sum())
            negative_sum = grouped["value"].apply(lambda x: x[x < 0].sum())

            max_val = positive_sum.max() if not positive_sum.empty else 0
            min_val = negative_sum.min() if not negative_sum.empty else 0
        else:
            max_val = data["value"].max()
            min_val = data["value"].min()

        return min_val, max_val

    if df is None or (isinstance(df, pd.DataFrame) and df.empty):
        return {"min": 0, "max": 0}

    scaling_factor = 1.2  # Scaling factor to add headroom to the y-axis

    min_val, max_val = compute_min_max(df)

    # Compare and get overall min and max if scenario2 exists
    if df2 is not None and not df2.empty:
        min_val2, max_val2 = compute_min_max(df2)
        min_val = min(min_val, min_val2)
        max_val = max(max_val, max_val2)

    min_val_scaled = min_val * scaling_factor if min_val < 0 else 0
    max_val_scaled = max_val * scaling_factor

    return {"min": min_val_scaled, "max": max_val_scaled}


def _generate_min_max_scales(df: pd.DataFrame, config_g: dict) -> tuple:
    """Generate minimum and maximum values.

    Parameters
    ----------
    df : pd.DataFrame
        input Dataframe
    config_g : dict
        Configuration dictionary containing tab name and years.

    Returns
    -------
    tuple
        tuple of Maximum and minimum values
    """
    if df is None or (isinstance(df, pd.DataFrame) and df.empty):
        return 0, 0

    groupby_list = []
    if config_g.get("fil_col"):
        groupby_list.append(config_g["fil_col"])

    # Separate hourly and yearly charts to define groupby list
    if "hourly" in config_g["graph_type"]:
        groupby_list.append("snapshot")
        max_scale = df.groupby(groupby_list)["value"].sum().max()
        min_sum = df.groupby(groupby_list)["value"].sum().min()
        min_scale = [0 if min_sum >= 0 else min_sum][0]
        return max_scale * 1.2, min_scale * 1.2

    if "yearly" in config_g["graph_type"]:
        groupby_list.append("year")

    # For yearly charts, groupby activates when the input is not a number
    if "max" in str(config_g["yaxis_scale"]["max"]):
        if config_g.get("tab_name") == "pow_bats_ep_ratio":
            # For EP Ratio chart, do not sum all values to calculate the max y value
            # since the chart is unstacked
            max_scale = df.groupby(groupby_list)["value"].max().max() * 1.2
        else:
            max_scale = df.groupby(groupby_list)["value"].sum().max() * 1.2
    else:
        max_scale = config_g["yaxis_scale"]["max"]

    if "min" in str(config_g["yaxis_scale"]["min"]):
        min_scale = df.groupby(groupby_list)["value"].sum().min() * 1.2
    else:
        min_scale = config_g["yaxis_scale"]["min"]

    return max_scale, min_scale


def calculate_yaxis_scales(func_name: str, config_g: dict) -> dict:
    """Calculate the maximum and minimum scales of the y-axis.

    Parameters
    ----------
    func_name : str
        Name of the function determining the calculation type.
    config_g : dict
        Configuration dictionary containing tab name and years.

    Returns
    -------
    dict
        Maximum and minimum scales of y-axis.
    """
    if func_name is None:
        st.write(f"Function not mapped: {func_name}")
        return 0

    df_sc1 = _grouping_data(
        sc_name=st.session_state.sce1,
        tab_name=config_g.get("tab_name"),
        years=st.session_state.get("sce1_years"),
        hourly="hourly" in func_name,
    )

    max_scale, min_scale = _generate_min_max_scales(df_sc1, config_g)

    if st.session_state.sce2:
        df_sc2 = _grouping_data(
            sc_name=st.session_state.sce2,
            tab_name=config_g.get("tab_name"),
            years=st.session_state.get("sce2_years"),
            hourly="hourly" in func_name,
        )
        max_scale_sc2, min_scale_sc2 = _generate_min_max_scales(df_sc2, config_g)

        # Update global max and min scales
        max_scale = max(max_scale, max_scale_sc2)
        min_scale = min(min_scale, min_scale_sc2)

    return {"max_scale": max_scale, "min_scale": min_scale}
