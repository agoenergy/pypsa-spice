# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Helper functions for handling Results section in visual app."""

import datetime as dt
from collections.abc import Callable
from typing import Any

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
    generate_default_colour_mapping_dict_for_chart,
    generate_color_mapping_dict_for_chart,
    line_with_secondary_y_hourly,
    simple_bar_hourly,
    simple_bar_yearly,
    simple_line_hourly,
    simple_line_yearly,
)

# pylint: disable=too-many-locals, too-many-branches

use_flexo()

GRAPHS_WITH_TIME_FILTERS = [
    "simple_bar_hourly",
    "simple_line_hourly",
    "line_with_secondary_y_hourly",
    "filtered_bar_hourly",
]

# =============================================================================
# Sidebar navigation + section headers
# =============================================================================


def generate_sidebar(table_of_content: list[str]) -> None:
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


def render_section_header(section_name: str) -> None:
    """Render a section header with an anchor for sidebar navigation."""
    anchor_id = slugify_text(section_name)
    st.markdown(f"<div id='{anchor_id}'></div>", unsafe_allow_html=True)
    st.markdown(f"#### {section_name}")


# =============================================================================
# Download helper functions
# =============================================================================


def create_download_csv_button(csv_data: bytes, download_id: str) -> None:
    """Create a CSV download button for each data table."""
    st.download_button(
        label="Download CSV",
        key=f"download_button_{download_id}",
        data=csv_data,
        file_name=f"{download_id}.csv",
        mime="text/csv",
    )


def render_download_with_data_table(
    df: pd.DataFrame,
    graph_config: dict[str, Any],
    scenario_name: str,
) -> None:
    """Render data table and CSV download for yearly charts."""
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    base_group_cols = {"year", leg_col}

    additional_group_cols = [
        col
        for col in df.columns
        if col not in base_group_cols.union({"value"}) and df[col].dtype == "object"
    ]
    group_cols = ["year", leg_col] + additional_group_cols
    df_grouped = df.groupby(group_cols, as_index=False)["value"].sum(min_count=1)

    with st.expander(f":material/database: Data ({scenario_name}):", expanded=False):
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
        df_pivot = df_pivot.loc[~(df_pivot == 0).all(axis=1)].fillna(0)
        df_pivot.index.names = pivot_index

        styled_df = df_pivot.style.apply(generate_diff_arrows, axis=None).format(
            "{:.1f}"
        )
        st.table(styled_df)

        csv_data = df_pivot.to_csv().encode("utf-8")
        create_download_csv_button(csv_data, download_id)


def render_download_without_data_table(
    filtered_df: pd.DataFrame,
    graph_config: dict[str, Any],
    scenario_name: str,
) -> None:
    """Render CSV download for hourly charts using pre-filtered data."""
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    graph_type = graph_config["graph_type"]

    if filtered_df is None or filtered_df.empty:
        return

    csv_data = (
        filtered_df.to_csv().encode("utf-8")
        if graph_type == "filtered_bar_hourly"
        else filtered_df[["snapshot", leg_col, "value"]].to_csv().encode("utf-8")
    )
    create_download_csv_button(csv_data, download_id)


def display_download_button_without_data(
    scenario_name: str,
    graph_config: dict[str, Any],
) -> None:
    """Display CSV download button for a graph without showing a data table."""
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)
    graph_type = graph_config["graph_type"]

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        year=str(graph_config["shared_years"]),
        country=graph_config["shared_country"],
    )
    if df is None or df.empty:
        return

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


def display_download_button_with_data(
    scenario_name: str,
    graph_config: dict[str, Any],
) -> None:
    """Display data table + CSV download button for yearly graphs."""
    leg_col = graph_config["leg_col"]
    download_id = graph_config["download_id"].format(scenario_name)

    df = read_result_csv(
        scenario_name,
        graph_config["table_name"],
        country=graph_config["shared_country"],
    )
    if df is None or df.empty:
        return

    base_group_cols = {"year", leg_col}
    additional_group_cols = [
        col
        for col in df.columns
        if col not in base_group_cols.union({"value"}) and df[col].dtype == "object"
    ]
    group_cols = ["year", leg_col] + additional_group_cols
    df_grouped = df.groupby(group_cols, as_index=False)["value"].sum(min_count=1)

    with st.expander(f":material/database: Data ({scenario_name}):", expanded=False):
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
        df_pivot = df_pivot.loc[~(df_pivot == 0).all(axis=1)].fillna(0)
        df_pivot.index.names = pivot_index

        styled_df = df_pivot.style.apply(generate_diff_arrows, axis=None).format(
            "{:.2f}"
        )
        st.table(styled_df)

        csv_data = df_pivot.to_csv().encode("utf-8")
        create_download_csv_button(csv_data, download_id)


# =============================================================================
# Colour mapping
# =============================================================================


def generate_colour_mapping_dict(
    table_name: str,
    mapping_df: pd.DataFrame | None,
    df: pd.DataFrame,
    leg_col: str,
) -> dict[str, str]:
    """Return colour mapping for a chart based on mapping files or defaults."""
    unique_legends = df["legend"].unique().tolist()
    if mapping_df is not None:
        return generate_color_mapping_dict_for_chart(table_name, unique_legends)
    return generate_default_colour_mapping_dict_for_chart(df, leg_col)


# =============================================================================
# Filter widgets
# =============================================================================


def setup_year_filter(config_plot: dict[str, Any], is_dual_scenario: bool) -> str:
    """Set up the year filter (pills) for graphs with hourly data."""
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

    def pills_widget() -> str:
        return st.pills(
            label,
            options=years,
            key=key,
            default=years[0],
            label_visibility="collapsed",
        )

    return pills_widget()


def setup_country_filter(
    config_plot: dict[str, Any],
    is_dual_scenario: bool = False,
    scenario_tag: str | None = None,
) -> str | None:
    """Set up the country filter (pills) for charts that support country selection."""
    df = read_result_csv(scenario_tag, "pow_cap_by_type_yearly")
    if df is None or df.empty or "country" not in df.columns:
        return None

    if is_dual_scenario:
        country_options = sorted(set(df["country"].unique().tolist()))
    else:
        country_options = df["country"].unique().tolist()

    units = config_plot.get("units")
    table_name = config_plot.get("table_name", "")
    is_regional_hourly = (
        "flow" in table_name
        or "charging" in table_name
        or ("region" in table_name and "price" not in table_name)
    )
    if units and not is_regional_hourly:
        country_options += ["ALL"]

    slider_id = config_plot["table_name"]
    country_id = config_plot.get("shared_country", "ALL")
    key = f"shared_country_{country_id}_{scenario_tag}_{slider_id}"
    label = f"{slider_id} Select country:"

    return st.pills(
        label,
        options=country_options,
        key=key,
        default=country_options[0],
        label_visibility="collapsed",
    )


def setup_region_filter(
    config_plot: dict[str, Any],
    df1: pd.DataFrame,
    df2: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> str | None:
    """Set up the region filter (pills) for filtered_bar_hourly graphs."""
    if "fil_col" not in config_plot:
        return None

    fil_col = config_plot["fil_col"]
    if is_dual_scenario and df2 is not None:
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

    return st.pills(
        label,
        options=region_options,
        key=key,
        default=region_options[0],
        label_visibility="collapsed",
    )


def setup_month_filter(
    config_plot: dict[str, Any],
    df1: pd.DataFrame,
    df2: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> int:
    """Set up month selection filter (for complete hourly datasets)."""
    if is_dual_scenario and df2 is not None:
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

    return st.segmented_control(
        label,
        options=months_all,
        format_func=lambda x: months_names[x],
        key=key,
        default=months_all[0],
        label_visibility="collapsed",
    )


def setup_date_filter_complete(
    config_plot: dict[str, Any],
    df1_m: pd.DataFrame,
    df2_m: pd.DataFrame | None = None,
    shared_year: str | None = None,
    selected_month: int | None = None,
    is_dual_scenario: bool = False,
) -> tuple[dt.datetime, dt.datetime]:
    """Set up date range slider for complete hourly datasets."""
    if is_dual_scenario and df2_m is not None:
        min_date = min(df1_m["snapshot"].min(), df2_m["snapshot"].min())
        max_date = max(df1_m["snapshot"].max(), df2_m["snapshot"].max())
        scenario_text = "both"
    else:
        min_date = df1_m["snapshot"].min()
        max_date = df1_m["snapshot"].max()
        scenario_text = st.session_state.sce1

    slider_id = config_plot["slider_id"].format(scenario_text)
    label = f"{slider_id} Select Date Range:"

    # Defaults: first half of the month
    assert shared_year is not None
    assert selected_month is not None

    return st.slider(
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


def setup_date_filter_incomplete(
    config_plot: dict[str, Any],
    df1_m: pd.DataFrame,
    df2_m: pd.DataFrame | None = None,
    is_dual_scenario: bool = False,
) -> tuple[pd.Timestamp, pd.Timestamp]:
    """Set up integer range slider for incomplete hourly datasets."""
    if is_dual_scenario and df2_m is not None:
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

    row_range = st.slider(
        label=label,
        min_value=1,
        max_value=num_timestamps,
        value=(1, min(20, num_timestamps)),
        step=1,
        label_visibility="collapsed",
    )
    return all_timestamps[row_range[0] - 1], all_timestamps[row_range[1] - 1]


def setup_radio_button_filter(
    config_plot: dict[str, Any], is_dual_scenario: bool
) -> str | None:
    """Set up the radio button filter that appears in bar_with_filter graphs."""
    if config_plot.get("graph_type") != "bar_with_filter":
        return None

    df1 = read_result_csv(st.session_state.sce1, config_plot["table_name"])
    if df1 is None or df1.empty:
        return None

    if not is_dual_scenario:
        return None

    df2 = read_result_csv(st.session_state.sce2, config_plot["table_name"])
    if df2 is None or df2.empty:
        return None

    fil_col = config_plot["fil_col"]
    filter_options = sorted(
        set(df1[fil_col].unique().tolist() + df2[fil_col].unique().tolist())
    )
    slider_id = config_plot["slider_id"].format("both")

    return st.radio(
        f"{slider_id} Select {fil_col} (both):",
        options=[str(x) for x in filter_options],
        format_func=prettify_label,
        horizontal=True,
        label_visibility="collapsed",
    )


def setup_hourly_filters(
    config_dict: dict[str, Any],
    scenario_1_raw: pd.DataFrame,
    scenario_2_raw: pd.DataFrame | None,
    has_dual: bool,
) -> None:
    """Set up region and date filters for hourly data and update config_dict."""
    config_dict["shared_region"] = setup_region_filter(
        config_dict,
        scenario_1_raw,
        scenario_2_raw,
        has_dual,
    )
    filter_results = setup_hourly_data_filters(
        scenario_1_raw,
        scenario_2_raw,
        config_dict,
        has_dual,
    )
    config_dict.update(filter_results)


def setup_hourly_data_filters(
    df1: pd.DataFrame,
    df2: pd.DataFrame | None,
    config_plot: dict[str, Any],
    is_dual_scenario: bool,
) -> dict[str, Any]:
    """Set up relevant filters for graphs with hourly data."""
    shared_year = config_plot.get("shared_year", None)
    shared_region = config_plot.get("shared_region", None)

    if shared_region is not None:
        df1 = df1[df1[config_plot["fil_col"]] == shared_region]
        if df2 is not None:
            df2 = df2[df2[config_plot["fil_col"]] == shared_region]

    is_complete = all(len(df) % 8760 == 0 for df in [df1, df2] if df is not None)
    is_empty = all(df.empty for df in [df1, df2] if df is not None)

    if is_empty:
        return {
            "shared_years": None,
            "shared_months": None,
            "shared_dates": None,
            "shared_region": shared_region if shared_region else None,
        }

    if is_complete:
        selected_month = setup_month_filter(config_plot, df1, df2, is_dual_scenario)

        df1_m = filter_dataframe_by_month(df1, selected_month)
        df2_m = (
            filter_dataframe_by_month(df2, selected_month) if df2 is not None else None
        )

        selected_dates = setup_date_filter_complete(
            config_plot,
            df1_m,
            df2_m,
            shared_year,
            selected_month,
            is_dual_scenario,
        )
    else:
        selected_month = None
        df1_m, df2_m = df1, df2
        selected_dates = setup_date_filter_incomplete(
            config_plot,
            df1_m,
            df2_m,
            is_dual_scenario,
        )

    return {
        "shared_years": shared_year,
        "shared_months": selected_month,
        "shared_dates": selected_dates,
        "shared_region": shared_region if shared_region else None,
    }


# =============================================================================
# Plot function mapping + misc helpers
# =============================================================================


def map_chart_to_plot_function(
    func_name: str | None = None,
) -> Callable[..., Any] | None:
    """Map plot function name to the corresponding chart function."""
    mapping: dict[str, Callable[..., Any]] = {
        "simple_bar_yearly": simple_bar_yearly,
        "simple_line_yearly": simple_line_yearly,
        "bar_with_filter": bar_with_filter,
        "area_share_yearly": area_share_yearly,
        "simple_bar_hourly": simple_bar_hourly,
        "simple_line_hourly": simple_line_hourly,
        "filtered_bar_hourly": filtered_bar_hourly,
        "line_with_secondary_y_hourly": line_with_secondary_y_hourly,
    }

    func = mapping.get(func_name or "")
    if func is None:
        st.write(f"Function not mapped: {func_name}")
    return func


# =============================================================================
# Layout renderers (single/dual chart layout)
# =============================================================================


def render_single_chart_layout(
    vis_display_data: pd.DataFrame,
    table_display_data: pd.DataFrame,
    config_dict: dict[str, Any],
    mapping_df: pd.DataFrame | None,
    y_range: dict[str, Any],
    plot_function: Callable[..., Any],
    render_download_function: Callable[..., Any],
    **plot_kwargs: Any,
) -> None:
    """Render single scenario chart with download button."""
    scenario_name = st.session_state.sce1
    table_name = config_dict["table_name"]
    legend_col = config_dict["leg_col"]

    st.markdown(f"{scenario_name}")

    colour_mapping = generate_colour_mapping_dict(
        table_name,
        mapping_df,
        vis_display_data,
        legend_col,
    )

    plot_key = plot_kwargs.pop("key", f"plotly_chart_{scenario_name}_{table_name}")
    plot_function(
        vis_display_data,
        config_dict,
        colour_mapping,
        y_range,
        key=plot_key,
        **plot_kwargs,
    )

    render_download_function(table_display_data, config_dict, scenario_name)


def render_dual_chart_layout(
    scenario_1_vis_display_data: pd.DataFrame,
    scenario_2_vis_display_data: pd.DataFrame | None,
    table_1_display_data: pd.DataFrame,
    table_2_display_data: pd.DataFrame | None,
    config_dict: dict[str, Any],
    mapping_df: pd.DataFrame | None,
    y_range: dict[str, Any],
    plot_function: Callable[..., Any],
    render_download_function: Callable[..., Any],
    **plot_kwargs: Any,
) -> None:
    """Render dual scenario charts side-by-side with download buttons."""
    if scenario_2_vis_display_data is None or table_2_display_data is None:
        raise ValueError("render_dual_chart_layout requires non-None scenario 2 data")

    scenario_1_name = st.session_state.sce1
    scenario_2_name = st.session_state.sce2
    table_name = config_dict["table_name"]
    legend_col = config_dict["leg_col"]

    plot_kwargs.pop("key", None)

    col1, col2 = st.columns([6, 6])
    with col1:
        st.caption(f"{scenario_1_name}")
        colour_mapping_1 = generate_colour_mapping_dict(
            table_name,
            mapping_df,
            scenario_1_vis_display_data,
            legend_col,
        )
        plot_function(
            scenario_1_vis_display_data,
            config_dict,
            colour_mapping_1,
            y_range,
            key=f"plotly_chart_{scenario_1_name}_{table_name}",
            **plot_kwargs,
        )

    with col2:
        st.caption(f"{scenario_2_name}")
        colour_mapping_2 = generate_colour_mapping_dict(
            table_name,
            mapping_df,
            scenario_2_vis_display_data,
            legend_col,
        )
        plot_function(
            scenario_2_vis_display_data,
            config_dict,
            colour_mapping_2,
            y_range,
            key=f"plotly_chart_{scenario_2_name}_{table_name}",
            **plot_kwargs,
        )

    col1, col2 = st.columns([6, 6])
    with col1:
        render_download_function(table_1_display_data, config_dict, scenario_1_name)
    with col2:
        render_download_function(table_2_display_data, config_dict, scenario_2_name)


# =============================================================================
# Styling helper (table arrows)
# =============================================================================


def generate_diff_arrows(data: pd.DataFrame) -> pd.DataFrame:
    """Highlight differences between adjacent cells in a dataframe."""
    styled = pd.DataFrame("", index=data.index, columns=data.columns)
    numeric_data = data.apply(pd.to_numeric, errors="coerce")

    for idx in range(len(numeric_data)):
        for col in range(len(numeric_data.columns)):
            styled.iloc[idx, col] = "padding-right: 20px; text-align: right;"

    for idx in range(len(numeric_data)):
        row_data = numeric_data.iloc[idx]
        base_value = row_data.iloc[0]
        safe_base_value = 0.01 if base_value == 0 else base_value

        for col in range(1, len(numeric_data.columns)):
            current_value = row_data.iloc[col]
            if (
                pd.isna(base_value)
                or pd.isna(current_value)
                or base_value == current_value
            ):
                continue

            pct_change = (current_value - base_value) / abs(safe_base_value)
            intensity = min(abs(pct_change), 1.0) * 0.8 + 0.2

            if pct_change > 0:
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
                """  # noqa: E501
            elif pct_change < 0:
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
                """  # noqa: E501

    return styled


# =============================================================================
# Main page renderer
# =============================================================================


def render_st_page_and_plot(
    graph_type_func: Callable[..., Any],
    config_plot: dict[str, Any],
) -> None:
    """Render and plot all graphs based on the provided graph type and configuration."""
    is_dual_scenario = bool(st.session_state.sce2)

    if is_dual_scenario:
        apply_radio_menu_styles()

    config_plot["shared_country"] = setup_country_filter(
        config_plot,
        is_dual_scenario,
        scenario_tag=st.session_state.sce1,
    )

    if config_plot.get("graph_type") == "bar_with_filter":
        shared_filter = setup_radio_button_filter(config_plot, is_dual_scenario)
        if shared_filter:
            config_plot["shared_filter"] = shared_filter

    elif config_plot.get("graph_type") in GRAPHS_WITH_TIME_FILTERS:
        shared_year = setup_year_filter(config_plot, is_dual_scenario)
        config_plot["shared_year"] = str(shared_year)

        df1 = read_result_csv(
            st.session_state.sce1,
            config_plot["table_name"],
            year=str(shared_year),
        )
        if df1 is None or df1.empty:
            return
        df1 = df1.copy()
        df1["snapshot"] = pd.to_datetime(df1["snapshot"])

        df2 = None
        if is_dual_scenario:
            df2 = read_result_csv(
                st.session_state.sce2,
                config_plot["table_name"],
                year=str(shared_year),
            )
            if df2 is not None and not df2.empty:
                df2 = df2.copy()
                df2["snapshot"] = pd.to_datetime(df2["snapshot"])

        config_plot["shared_region"] = setup_region_filter(
            config_plot,
            df1,
            df2,
            is_dual_scenario,
        )

        filter_results = setup_hourly_data_filters(
            df1, df2, config_plot, is_dual_scenario
        )
        config_plot.update(filter_results)

    if not is_dual_scenario:
        config_plot["years"] = st.session_state.sce1_years
        st.markdown(f"#### {st.session_state.sce1} ")
        graph_type_func(scenario_name=st.session_state.sce1, graph_config=config_plot)

        if config_plot.get("graph_type") in GRAPHS_WITH_TIME_FILTERS:
            display_download_button_without_data(st.session_state.sce1, config_plot)
        else:
            display_download_button_with_data(st.session_state.sce1, config_plot)
    else:
        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            config_plot["years"] = st.session_state.sce1_years
            st.markdown(f"#### {st.session_state.sce1} ")
            graph_type_func(
                scenario_name=st.session_state.sce1, graph_config=config_plot
            )

        with col3:
            config_plot["years"] = st.session_state.sce2_years
            st.markdown(f"#### {st.session_state.sce2} ")
            graph_type_func(
                scenario_name=st.session_state.sce2, graph_config=config_plot
            )

        col1, _, col3 = st.columns([6, 1, 6])
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
