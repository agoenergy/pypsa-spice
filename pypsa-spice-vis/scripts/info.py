# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create info page under General section showing technological nomenclature used in
PyPSA-SPICE.
"""

import streamlit as st
import os
import pandas as pd

st.title("PyPSA-SPICE technological nomenclature")
st.divider()

with open(
    os.path.join(st.session_state.current_dir, "setting/tech_mapping.csv"), "r"
) as file:
    tech_info = pd.read_csv(file)[["original_names", "nice_names"]]
    tech_info.columns = ["Abbreviations", "Full names"]
    tech_info = tech_info.set_index("Abbreviations")

st.dataframe(tech_info, use_container_width=True)
