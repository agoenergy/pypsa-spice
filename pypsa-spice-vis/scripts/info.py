# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

import os

import pandas as pd
import streamlit as st

st.title("PyPSA-SPICE technological nomenclature")
st.divider()

# Tech table
with open(
    os.path.join(st.session_state.current_dir, "setting/tech_mapping.csv"), "r"
) as file:
    tech_info = pd.read_csv(file)[["original_names", "nice_names"]]
    tech_info.columns = ["Abbreviations", "Full names"]
    tech_info = tech_info.set_index("Abbreviations")

# Display DataFrame
st.dataframe(tech_info, use_container_width=True)
