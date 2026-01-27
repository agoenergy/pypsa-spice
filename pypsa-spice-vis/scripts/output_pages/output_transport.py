# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Transport page under Results section.

Page shows editable transport related
dataframes and visualisations from the modelling results.
"""

import os

import streamlit as st
import yaml

from scripts.output_st_handler import (
    generate_sidebar,
    map_chart_to_plot_function,
    render_st_page_and_plot_settings,
)

st.title(":material/directions_car: Transport")

with open(
    os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"),
    encoding="utf-8",
) as file:
    config = yaml.safe_load(file)["transport"]

table_of_content = [config[item]["name"] for item in config]

for _item, values in config.items():
    render_st_page_and_plot_settings(
        graph_type=map_chart_to_plot_function(values["graph_type"]),
        config_plot=values,
    )


generate_sidebar(table_of_content)
