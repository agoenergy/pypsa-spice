# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Store data utility functions used across handler modules."""

import datetime as dt
import os
import re
from typing import Any

import pandas as pd
import streamlit as st

# Scaling factor to add some headroom to the y-axis when comparing two scenarios
SCALING_FACTOR = 1.2
# Threshold below which values are set to 0.0 for cleaner plots
SMALL_VALUE_THRESHOLD = 1e-6

# =============================================================================
# General text helpers
# =============================================================================


def prettify_label(label: str) -> str:
    """Convert snake_case or camelCase strings into readable text.

    Parameters
    ----------
    label : str
        The input string to format.

    Returns
    -------
    str
        The formatted, human-readable string.
    """
    camel_re = re.compile(r"[a-z][A-Z]")
    split_camel_re = re.compile(r"([a-z])([A-Z])")

    if "_" in label:
        if "_to_" in label:
            parts = label.split("_")
            if len(parts) == 5 and parts[2] == "to":
                return f"{parts[0]} ({parts[1]}) to {parts[3]} ({parts[4]})"
            return " ".join(parts)
        return " ".join(label.split("_"))

    if camel_re.search(label):
        spaced = split_camel_re.sub(r"\1 \2", label)
        return spaced.capitalize()

    return label


def slugify_text(text: str) -> str:
    """Slugify text string for safe anchor IDs (URL fragments)."""
    text = text.lower()
    text = re.sub(r"[^a-z0-9]+", "-", text)
    return text.strip("-")


def convert_month_to_name(month_num: int) -> str:
    """Convert a month number to the abbreviated month name."""
    return dt.datetime.strptime(str(month_num), "%m").strftime("%b")


# =============================================================================
# I/O helpers (reading model results)
# =============================================================================


def read_result_csv(
    scenario_name: str,
    table_name: str,
    country: str | None = None,
    year: str | None = None,
) -> pd.DataFrame | None:
    """Read model output CSV for a given scenario and table name.

    Parameters
    ----------
    scenario_name : str
        Selected scenario in Streamlit UI.
    table_name : str
        Output table name.
    country : str | None, optional
        If not None, filter the CSV by inputted country.
    year : str | None, optional
        If not None, read CSV from year-specific folder else all_years folder.

    Returns
    -------
    pd.DataFrame | None
        DataFrame containing the CSV data, or None if file not found.
    """
    year_dir = year if year else "all_years"
    file_path = os.path.join(
        st.session_state.result_path,
        scenario_name,
        "csvs",
        st.session_state.sector,
        year_dir,
        f"{table_name}.csv",
    )

    try:
        df = pd.read_csv(os.path.abspath(file_path))
    except FileNotFoundError:
        with st.container(height=450, border=True):
            st.write(f":material/warning: File does not exist or is empty: {file_path}")
        return None

    if "country" in df.columns and country is not None and country != "ALL":
        df = df[df["country"] == country]

    return df.fillna(0)


def load_and_validate_hourly_data(
    scenario_name: str, table_name: str, year: str, country: str
) -> pd.DataFrame | None:
    """Load hourly data CSV and convert `snapshot` column to datetimes.

    Returns None if the file is missing or empty.
    """
    raw_data = read_result_csv(scenario_name, table_name, year=year, country=country)
    if raw_data is None or raw_data.empty:
        return None
    raw_data = raw_data.copy()
    raw_data["snapshot"] = pd.to_datetime(raw_data["snapshot"])
    return raw_data


# =============================================================================
# Data cleanup + normalization
# =============================================================================


def handle_small_values(df: pd.DataFrame) -> pd.DataFrame:
    """Convert <1e-6 values in the "value" column to 0.0.

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
        df.loc[df["value"].abs() < SMALL_VALUE_THRESHOLD, "value"] = 0.0
    return df


def normalize_dataframe(df: pd.DataFrame | pd.Series) -> pd.DataFrame:
    """Ensure groupby results are returned as a DataFrame, not a Series.

    If a Pandas Series is passed (common after .groupby(...).sum()), the
    function resets the index and returns a DataFrame.
    """
    return df.reset_index() if isinstance(df, pd.Series) else df


def add_nice_names(
    df: pd.DataFrame, leg_col: str, mapping_df: pd.DataFrame | None
) -> pd.DataFrame:
    """Add a 'legend' column using mapping_df or prettified labels."""
    df = df.copy()

    def _label(value: Any) -> str:
        if mapping_df is not None and value in mapping_df.index:
            return str(mapping_df.loc[value, "nice_names"])
        return prettify_label(str(value))

    df["legend"] = df[leg_col].map(_label)
    return df


def clean_df_for_plotting(leg_col: str, df: pd.DataFrame) -> pd.DataFrame:
    """Clean the data used to plot the graph.

    1. Filter out legend series where all values are zero or NaN (hide in plots).
    2. Convert all values < 1e-6 to 0.0.

    Parameters
    ----------
    leg_col : str
        The legend column name as per the graph's configuration dictionary.
    df : pd.DataFrame
        The raw data returned by read_result_csv.

    Returns
    -------
    pd.DataFrame
        Dataframe with zero/NaN legend series removed and small values converted.
    """
    df_pivoted = df.pivot_table(
        values="value",
        columns="year" if "year" in df.columns else "snapshot",
        index=leg_col,
        aggfunc="mean",
        dropna=False,
    )

    all_legends_to_remove = df_pivoted[
        (df_pivoted.isna() | (df_pivoted == 0)).all(axis=1)
    ].index

    df_filtered = df[~df[leg_col].isin(all_legends_to_remove)]
    return handle_small_values(df_filtered)


# =============================================================================
# Date filtering helpers (hourly data)
# =============================================================================


def filter_dataframe_by_date_range(
    df: pd.DataFrame,
    start_date: dt.datetime | None,
    end_date: dt.datetime | None,
) -> pd.DataFrame:
    """Filter input dataframe by a specific date range.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with a 'snapshot' column.
    start_date : dt.datetime | None
        Start date for filtering. If None, no lower bound is applied.
    end_date : dt.datetime | None
        End date for filtering. If None, no upper bound is applied.

    Returns
    -------
    pd.DataFrame
        Filtered dataframe.
    """
    df = df.copy()
    df["Date"] = pd.to_datetime(df["snapshot"])

    if start_date is not None and end_date is not None:
        return df[(df["Date"] >= start_date) & (df["Date"] <= end_date)]
    if start_date is not None:
        return df[df["Date"] >= start_date]
    if end_date is not None:
        return df[df["Date"] <= end_date]
    return df


def filter_dataframe_by_month(df: pd.DataFrame, month: int) -> pd.DataFrame:
    """Filter input dataframe by month.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with a 'snapshot' column.
    month : int
        Month number (1-12) to filter by.

    Returns
    -------
    pd.DataFrame
        Filtered dataframe containing only rows from the specified month.
    """
    df = df.copy()
    df["Date"] = pd.to_datetime(df["snapshot"])
    df["Month"] = df["Date"].dt.month
    return df[df["Month"] == month]


def get_filtered_df_and_date_range(
    df: pd.DataFrame, graph_config: dict[str, Any]
) -> tuple[pd.DataFrame, dt.datetime | None, dt.datetime | None, bool]:
    """Get filtered data and date range for graphs with date filters.

    Relevant graphs are:
    - simple_bar_hourly
    - simple_line_hourly
    - line_with_secondary_y_hourly
    - filtered_bar_hourly

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of the scenario.
    graph_config : dict[str, Any]
        Configuration dictionary that may contain:
        - 'shared_months': int | None
        - 'shared_dates': tuple[datetime | None, datetime | None]
        - 'shared_region': str (optional)
        - 'fil_col': str (optional)

    Returns
    -------
    tuple[pd.DataFrame, dt.datetime | None, dt.datetime | None, bool]
        (filtered_df, start_date, end_date, is_complete)
    """
    month = graph_config.get("shared_months")
    start_date, end_date = graph_config["shared_dates"]

    is_complete = len(df) % 8760 == 0
    df = df.copy()
    df["snapshot"] = pd.to_datetime(df["snapshot"])

    if "shared_region" in graph_config and "fil_col" in graph_config:
        fil_col = graph_config["fil_col"]
        df = df[df[fil_col] == graph_config["shared_region"]]

    df_m = filter_dataframe_by_month(df=df, month=month) if month is not None else df
    return df_m, start_date, end_date, is_complete


# pylint: disable=too-many-arguments
def filter_and_prepare_hourly_data(
    raw_data: pd.DataFrame | None,
    config_dict: dict[str, Any],
    legend_col: str,
    mapping_df: pd.DataFrame | None,
) -> tuple[pd.DataFrame | None, dt.datetime | None, dt.datetime | None, bool]:
    """Filter hourly raw data by date range and prepare for plotting.

    Returns (filtered_data, start_date, end_date, is_complete). If `raw_data`
    is None the function returns (None, None, None, False).
    """
    if raw_data is None:
        return None, None, None, False

    monthly_data, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        raw_data, config_dict
    )
    filtered_data = filter_dataframe_by_date_range(
        monthly_data, start_date=start_date, end_date=end_date
    )
    filtered_data = handle_small_values(filtered_data)
    filtered_data = add_nice_names(filtered_data, legend_col, mapping_df)
    return filtered_data, start_date, end_date, is_complete


# pylint: enable=too-many-arguments

# =============================================================================
# Y-axis scaling helpers
# =============================================================================


def _compute_min_max(data: pd.DataFrame, group_col: str | None) -> tuple[float, float]:
    """Compute min/max values for y-scale, optionally grouped by group_col."""
    if group_col and group_col in data.columns:
        grouped = data.groupby(group_col)

        positive_sum = grouped["value"].apply(lambda x: x[x > 0].sum())
        negative_sum = grouped["value"].apply(lambda x: x[x < 0].sum())

        max_val = float(positive_sum.max()) if not positive_sum.empty else 0.0
        min_val = float(negative_sum.min()) if not negative_sum.empty else 0.0
        return min_val, max_val

    return float(data["value"].min()), float(data["value"].max())


def calculate_min_max_y_scale(
    df: pd.DataFrame, df2: pd.DataFrame | None, group_col: str | None
) -> dict[str, float]:
    """Calculate minimum and maximum values to use on the y axis.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame. Assumes it contains a "value" column.
    df2 : pd.DataFrame | None
        Optional second input DataFrame to compare with for a joint scale.
    group_col : str | None
        Column in the df to group by first. If None, calculates overall min/max.

    Returns
    -------
    dict[str, float]
        Dictionary with keys:
        - "min": float
        - "max": float
    """
    if df is None or df.empty:
        return {"min": 0.0, "max": 0.0}

    min_val, max_val = _compute_min_max(df, group_col)

    if df2 is not None and not df2.empty:
        min_val2, max_val2 = _compute_min_max(df2, group_col)
        min_val = min(min_val, min_val2)
        max_val = max(max_val, max_val2)

    min_val_scaled = min_val * SCALING_FACTOR if min_val < 0 else 0.0
    max_val_scaled = max_val * SCALING_FACTOR
    return {"min": min_val_scaled, "max": max_val_scaled}


def prepare_y_range(
    scenario_1_df: pd.DataFrame,
    scenario_2_df: pd.DataFrame | None,
    x_col: str | None,
) -> dict[str, float]:
    """Calculate y-axis range for consistent scaling across scenarios.

    Parameters
    ----------
    scenario_1_df : pd.DataFrame
        First scenario dataframe.
    scenario_2_df : pd.DataFrame | None
        Optional second scenario dataframe.
    x_col : str | None
        Grouping column used when computing min/max (e.g., 'year' or 'snapshot').

    Returns
    -------
    dict[str, float]
        Dictionary with keys 'max_scale' and 'min_scale'.
    """
    y_range = calculate_min_max_y_scale(scenario_1_df, scenario_2_df, x_col)
    return {"max_scale": y_range["max"], "min_scale": y_range["min"]}
