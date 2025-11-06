# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Costs page under Results section showing editable costs related
dataframes and visualisations from the modelling results.
"""

import os

import streamlit as st
import yaml

from scripts.output_st_handler import (
    generate_sidebar,
    map_chart_to_plot_function,
    plot_indicator,
)

st.title(":material/attach_money: Costs")

with open(
    os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"),
    encoding="utf-8",
) as file:
    config = yaml.safe_load(file)["costs"]

# table_of_content = [config[item]["name"] for item in config]
table_of_content = []


for item, values in config.items():
    if (
        values["incl_sector"] == "all"
        or values["incl_sector"] in st.session_state.sector
    ):
        plot_indicator(
            graph_type=map_chart_to_plot_function(values["graph_type"]),
            config_plot=values,
        )
        table_of_content.append(values["name"])
    else:
        pass

generate_sidebar(table_of_content)
