# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Utility functions that are DataFrame or calculation related and used across handler 
modules.
"""

import os
import pandas as pd
from typing import Dict
import re
import streamlit as st
import datetime as dt


def handle_small_values(df: pd.DataFrame) -> pd.DataFrame:
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


def calculate_min_max_y_scale(
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


def clean_df_for_plotting(leg_col: str, df: pd.DataFrame):
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

    df_filtered = handle_small_values(df_filtered)

    return df_filtered


def convert_month_to_name(month_num: int) -> str:
    """Convert a month number to the abbreviated month name."""
    return dt.datetime.strptime(str(month_num), "%m").strftime("%b")


def filter_dataframe_by_date_range(
    df: pd.DataFrame, start_date: dt.datetime, end_date: dt.datetime
):
    """Filter input dataframe by specific date range."""
    df["Date"] = pd.to_datetime(df["snapshot"])
    filtered_df = df[(df["Date"] >= start_date) & (df["Date"] <= end_date)]

    return filtered_df


def filter_dataframe_by_month(df: pd.DataFrame, month: int):
    """Filter input dataframe by month."""
    df["Date"] = pd.to_datetime(df["snapshot"])
    df["Month"] = df["Date"].dt.month

    return df[df["Month"] == month]


def get_filtered_df_and_date_range(df: pd.DataFrame, graph_config: dict):
    """Get the filtered data, start date, and end date for graphs with date filters.

    Relevant graphs are:
    - simple_bar_hourly
    - simple_line_hourly
    - line_with_secondary_y_hourly
    - filtered_bar_hourly

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame of the scenario
    graph_config : dict
        Configuration dictionary that may contain:
        - 'shared_years': int, shared year for both scenarios
        - 'shared_months': int, shared month for both scenarios
        - 'shared_dates': tuple of datetime objects, shared date range for both
        scenarios

    Returns
    -------
    Pd.DataFrame
        The data filtered by the selected month
    """
    month = graph_config["shared_months"]
    start_date, end_date = graph_config["shared_dates"]

    # Check if data is complete (no discontinuous hours)
    is_complete = len(df) % 8760 == 0

    df["snapshot"] = pd.to_datetime(df["snapshot"])

    # Add region filtering to the config dict for filtered_bar_hourly
    if "shared_region" in graph_config and "fil_col" in graph_config:
        fil_col = graph_config["fil_col"]
        df = df[df[fil_col] == graph_config["shared_region"]]

    # Only filter by month if month selector exists
    df_m = filter_dataframe_by_month(df=df, month=month) if month is not None else df

    return df_m, start_date, end_date, is_complete


def get_hourly_dfs_for_both_scenarios(
    graph_config: Dict,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Get the filtered hourly dataframes for two scenarios

    Parameters
    ----------
    graph_config : Dict
        Configuration dictionary for the current graph

    Returns
    -------
    tuple[pd.DataFrame, pd.DataFrame]
        The filtered hourly dataframes for the two scenarios
    """
    filtered_dfs = []
    for i, scenario in enumerate([st.session_state.sce1, st.session_state.sce2]):
        df = read_result_csv(
            scenario,
            graph_config["table_name"],
            year=str(graph_config["shared_years"]),
            country=graph_config["shared_country"],
        )

        if df is not None and not df.empty:
            df_m, s_date, e_date, _ = get_filtered_df_and_date_range(df, graph_config)
            if i == 0:
                start_date, end_date = s_date, e_date
            filtered_df = filter_dataframe_by_date_range(
                df_m, start_date=start_date, end_date=end_date
            )
            filtered_df = handle_small_values(filtered_df)
            filtered_dfs.append(filtered_df)

    return filtered_dfs


def prettify_label(label: str) -> str:
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


def read_result_csv(
    scenario_name: str,
    table_name: str,
    country: str = None,
    year: str = None,
) -> pd.DataFrame:
    """Read model ouput csv files for a given scenario and table name.

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


def slugify_text(text: str):
    """Slugify text string.

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
