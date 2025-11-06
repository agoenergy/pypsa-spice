# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Power page under Results section showing editable power related 
dataframes and visualisations from the modelling results.
"""

import streamlit as st
import os
import yaml
from scripts.output_st_handler import (
    plot_indicator,
    map_chart_to_plot_function,
    generate_sidebar,
)

st.title(":material/bolt: Power")


with open(
    os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"), "r"
) as file:
    config = yaml.safe_load(file)["power"]

table_of_content = [config[item]["name"] for item in config]

for item, values in config.items():
    plot_indicator(
        graph_type=map_chart_to_plot_function(values["graph_type"]),
        config_plot=values,
    )

generate_sidebar(table_of_content)
