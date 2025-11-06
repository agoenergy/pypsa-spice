# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create an info page under the General section.

This page explains the technological nomenclature used in PyPSA-SPICE.
"""

import os

import pandas as pd
import streamlit as st

st.title("PyPSA-SPICE technological nomenclature")
st.divider()

with open(
    os.path.join(st.session_state.current_dir, "setting/tech_mapping.csv"),
    encoding="utf-8",
) as file:
    tech_info = pd.read_csv(file)[["original_names", "nice_names"]]
    tech_info.columns = ["Abbreviations", "Full names"]
    tech_info = tech_info.set_index("Abbreviations")

st.dataframe(tech_info, use_container_width=True)
