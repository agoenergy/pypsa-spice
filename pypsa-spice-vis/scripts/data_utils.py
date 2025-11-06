# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Utility functions that are independent of Streamlit and used across helper modules.
"""

import os
import pandas as pd
import re
import streamlit as st

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
