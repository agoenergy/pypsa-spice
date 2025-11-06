# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Transport page under Results section showing editable transport related 
dataframes and visualisations from the modelling results.
"""

import yaml
import os
import streamlit as st
from scripts.output_st_handler import (
    plot_indicator,
    map_chart_to_plot_function,
    generate_sidebar,
)

st.title(":material/directions_car: Transport")

with open(
    os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"), "r"
) as file:
    config = yaml.safe_load(file)["transport"]

table_of_content = [config[item]["name"] for item in config]

for item, values in config.items():
    plot_indicator(
        graph_type=map_chart_to_plot_function(values["graph_type"]),
        config_plot=values,
    )


generate_sidebar(table_of_content)
