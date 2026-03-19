# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create an info page under the General section.

This page explains the technological nomenclature used in PyPSA-SPICE.
"""

import streamlit as st

from scripts.data_utils import load_tech_info_mapping_df

st.title("PyPSA-SPICE technological nomenclature")
st.divider()

tech_info = load_tech_info_mapping_df().reset_index()

tech_info = tech_info[["original_names", "nice_names", "sector"]]
tech_info.columns = ["Abbreviations", "Full names", "Sector (for technology only)"]
tech_info = tech_info.set_index("Abbreviations")

st.dataframe(tech_info, width="stretch")
