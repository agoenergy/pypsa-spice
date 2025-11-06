# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Helper functions for handling Results section in visual app."""

import datetime as dt

import pandas as pd
import streamlit as st
from styles import apply_radio_menu_styles, apply_sidebar_chart_nav_styles, use_flexo

from scripts.data_utils import (
    convert_month_to_name,
    filter_dataframe_by_date_range,
    filter_dataframe_by_month,
    get_filtered_df_and_date_range,
    prettify_label,
    read_result_csv,
    slugify_text,
)
from scripts.plot_settings import (
    area_share_yearly,
    bar_with_filter,
    filtered_bar_hourly,
    line_with_secondary_y_hourly,
    simple_bar_hourly,
    simple_bar_yearly,
    simple_line_hourly,
    simple_line_yearly,
)

use_flexo()


# These are the graph types that have the year+month+date filter
GRAPHS_WITH_TIME_FILTERS = [
    "simple_bar_hourly",
    "simple_line_hourly",
    "line_with_secondary_y_hourly",
    "filtered_bar_hourly",
]


def generate_sidebar(table_of_content):
    """Generate sidebar with navigation links."""
    with st.sidebar:
        st.divider()

        apply_sidebar_chart_nav_styles()

        for section in table_of_content:
            anchor_id = slugify_text(section)
            st.markdown(
                f'<a href="#{anchor_id}" class="nav-link">{section}</a>',
                unsafe_allow_html=True,
            )


def setup_year_filter(config_plot: dict, is_dual_scenario: bool) -> str:
    """Set up the year filter that appears in graphs with hourly data.

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
        years = sorted(set(st.session_state.sce1_years + st.session_state.sce2_years))
        scenario_text = "both"
        key_prefix = "shared"
    else:
        years = st.session_state.sce1_years
        scenario_text = st.session_state.sce1
        key_prefix = "single"

    slider_id = config_plot["slider_id"].format(scenario_text)
    key = f"{key_prefix}_year_{config_plot['table_name']}"
    label = f"{slider_id} Select Year:"

    # Pills widget for the year filter
    pills_widget = lambda: st.pills(  # noqa:E731
        label,
        options=years,
        key=key,
        default=years[0],
        label_visibility="collapsed",
    )  # noqa:E731

    if is_dual_scenario:
        spacer1, filter_col, spacer2 = st.columns([1, 3, 1])
        with filter_col:
            shared_year = pills_widget()
    else:
        shared_year = pills_widget()

    return shared_year


def setup_country_filter(config_plot, is_dual_scenario=False, scenario_tag=None) -> str:
    """Set up the country filter that appears in filtered_bar_hourly graphs.

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
            country_options = sorted(set(df["country"].unique().tolist()))
            scenario_text = "both"  # noqa: F841
        else:
            country_options = df["country"].unique()
            scenario_text = st.session_state.sce1  # noqa: F841

        slider_id = config_plot["table_name"]
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
    config_plot: dict,
    df1: pd.DataFrame,
    df2: pd.DataFrame = None,
    is_dual_scenario: bool = False,
) -> str | None:
    """Set up the region filter that appears in filtered_bar_hourly graphs.

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
            set(df1[fil_col].unique().tolist() + df2[fil_col].unique().tolist())
        )
        scenario_text = "both"
    else:
        region_options = df1[fil_col].unique()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    key = f"shared_region_{config_plot['table_name']}"
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
    config_plot: dict,
    df1: pd.DataFrame,
    df2: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> int:
    """Set up month selection filter.

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
        months_all = sorted(months_sce1.union(months_sce2))
        scenario_text = "both"
    else:
        months_all = df1["snapshot"].dt.month.unique()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    key = f"shared_month_{config_plot['table_name']}"
    label = f"{slider_id} Select Month:"
    months_names = {m: convert_month_to_name(m) for m in months_all}

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
    config_plot: dict,
    df1_m: pd.DataFrame,
    df2_m: pd.DataFrame | None = None,
    shared_year: str = None,
    selected_month: str = None,
    is_dual_scenario: bool = False,
) -> tuple[dt.datetime, dt.datetime]:
    """Set up date range slider.

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
    config_plot: dict,
    df1_m: pd.DataFrame,
    df2_m: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Set up integer range slider.

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
        all_timestamps = sorted(timestamps_1.union(timestamps_2))
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


def setup_radio_button_filter(config_plot: dict, is_dual_scenario: bool) -> str | None:
    """Set up the radio button filter that appears in bar_with_filter graphs.

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

    df1 = read_result_csv(st.session_state.sce1, config_plot["table_name"])

    if is_dual_scenario:
        df2 = read_result_csv(st.session_state.sce2, config_plot["table_name"])
        fil_col = config_plot["fil_col"]
        filter_options = sorted(
            set(df1[fil_col].unique().tolist() + df2[fil_col].unique().tolist())
        )
        slider_id = config_plot["slider_id"].format("both")

        spacer1, filter_col, spacer2 = st.columns([1, 2, 1])  # noqa: F841
        with filter_col:
            shared_filter = st.radio(
                f"{slider_id} Select {fil_col} (both):",
                options=[str(x) for x in filter_options],
                format_func=prettify_label,
                horizontal=True,
                label_visibility="collapsed",
            )
    else:
        shared_filter = None

    return shared_filter


def setup_hourly_data_filters(
    df1: pd.DataFrame,
    df2: pd.DataFrame | None,
    config_plot: dict,
    is_dual_scenario: bool,
) -> dict:
    """Set up relevant filters for graphs with hourly data.

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
            df1_m = filter_dataframe_by_month(df1, selected_month)
            df2_m = (
                filter_dataframe_by_month(df2, selected_month)
                if df2 is not None
                else None
            )

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


def render_st_page_and_plot(graph_type, config_plot: dict):
    """Render and plot all graphs based on the provided graph type and configuration."""
    st.markdown(
        f"<div id='{config_plot['name'].replace(' ', '-')}'></div>",
        unsafe_allow_html=True,
    )
    st.markdown(f"#### {config_plot['name']}")

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

    elif config_plot.get("graph_type") in GRAPHS_WITH_TIME_FILTERS:
        # Setup year filter
        shared_year = setup_year_filter(config_plot, is_dual_scenario)
        config_plot["shared_year"] = str(shared_year)

        df1 = read_result_csv(
            st.session_state.sce1, config_plot["table_name"], year=str(shared_year)
        )
        df1["snapshot"] = pd.to_datetime(df1["snapshot"])

        df2 = [
            (
                read_result_csv(
                    st.session_state.sce2,
                    config_plot["table_name"],
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
        graph_type(scenario_name=st.session_state.sce1, graph_config=config_plot)

        # Display the data download part
        if config_plot.get("graph_type") in GRAPHS_WITH_TIME_FILTERS:
            display_download_button_without_data(st.session_state.sce1, config_plot)
        else:
            display_download_button_with_data(st.session_state.sce1, config_plot)

    else:
        # Display the graphs for each of the two scenarios
        col1, col2, col3 = st.columns([6, 1, 6])
        with col1:
            config_plot["years"] = st.session_state.sce1_years
            st.markdown(f"#### {st.session_state.sce1} ")
            graph_type(scenario_name=st.session_state.sce1, graph_config=config_plot)
        with col3:
            config_plot["years"] = st.session_state.sce2_years
            st.markdown(f"#### {st.session_state.sce2} ")
            graph_type(scenario_name=st.session_state.sce2, graph_config=config_plot)

        # Display the data download part
        col1, col2, col3 = st.columns([6, 1, 6])

        if config_plot.get("graph_type") in GRAPHS_WITH_TIME_FILTERS:
            with col1:
                display_download_button_without_data(st.session_state.sce1, config_plot)
            with col3:
                display_download_button_without_data(st.session_state.sce2, config_plot)
        else:
            with col1:
                display_download_button_with_data(st.session_state.sce1, config_plot)
            with col3:
                display_download_button_with_data(st.session_state.sce2, config_plot)

    st.divider()


def create_download_csv_button(csv_data, download_id):
    """Create a CSV download button for each data table."""
    st.download_button(
        label="Download CSV",
        key=f"download_button_{download_id}",
        data=csv_data,
        file_name=f"{download_id}.csv",
        mime="text/csv",
    )


def display_download_button_without_data(
    scenario_name: str, graph_config: dict[str, str]
):
    """Generate and display the CSV download button for a graph without showing data.

    Parameters
    ----------
    scenario_name : str
        The scenario name for the data that will be downloaded.
    graph_config : Dict[str, str]
        The configuration dictionary for the graph.
    """
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    graph_type = graph_config["graph_type"]

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )

    if df is not None and not df.empty:
        df_m, start_date, end_date, _ = get_filtered_df_and_date_range(df, graph_config)

        filtered_df = filter_dataframe_by_date_range(
            df_m, start_date=start_date, end_date=end_date
        )

        csv_data = (
            filtered_df.to_csv().encode("utf-8")
            if graph_type == "filtered_bar_hourly"
            else filtered_df[["snapshot", leg_col, "value"]].to_csv().encode("utf-8")
        )

        create_download_csv_button(csv_data, download_id)


def display_download_button_with_data(scenario_name: str, graph_config: dict[str, str]):
    """Generate and display the CSV download button for a with data table included.

    Parameters
    ----------
    scenario_name : str
        The scenario name for the data being displayed in the table.
    graph_config : Dict[str, str]
        The configuration dictionary for the graph.
    """
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        country=graph_config["shared_country"],
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

    with st.expander(f":material/database: Data ({scenario_name}):", expanded=False):
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
        styled_df = df_pivot.style.apply(generate_diff_arrows, axis=None).format(
            "{:.2f}"
        )
        st.table(styled_df)

        csv_data = df_pivot.to_csv().encode("utf-8")

        create_download_csv_button(csv_data, download_id)


def generate_diff_arrows(data):
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
                """  # noqa:E501
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
                            """  # noqa:E501

    return styled


def map_chart_to_plot_function(func_name: str = None):
    """Map plot function to corresponding chart type  from config yaml."""
    mapping = {
        "simple_bar_yearly": simple_bar_yearly,
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


def read_and_concatenate_hourly_data(
    scenario_name: str,
    table_name: str,
    list_of_years: list = None,
    is_hourly: bool = False,
) -> pd.DataFrame:
    """Read and concatenate data for hourly scenarios.

    Parameters
    ----------
    scenario_name : str
        scenario name
    table_name : str
        table name
    list_of_years : list
        a list of the years
    is_hourly: bool
        whether the table is hourly or not

    Returns
    -------
    DataFrame
        DataFrame containing the concatenated data
    """
    if is_hourly and list_of_years:
        df_combined = pd.concat(
            [
                read_result_csv(
                    scenario_name=scenario_name, table_name=table_name, year=str(year)
                )
                for year in list_of_years
            ],
            ignore_index=True,
        )
    else:
        df_combined = read_result_csv(
            scenario_name=scenario_name, table_name=table_name
        )
    return df_combined
