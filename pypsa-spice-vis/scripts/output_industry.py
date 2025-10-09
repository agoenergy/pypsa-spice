# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import streamlit as st
import os
import yaml
from scripts.output_helpers import (
    calculate_yaxis_scales, plot_indicator, plot_function, generate_sidebar
)

st.title(":material/construction: Industry")

with open(os.path.join(st.session_state.current_dir,
                       "setting/graph_settings.yaml"), "r") as file:
    config = yaml.safe_load(file)["industry"]

toc = [config[item]["name"] for item in config]
yaxis_scales_dict = {}

for item, values in config.items():
    if "with_filter" in values["graph_type"] and "hourly" in values["graph_type"]:
        yaxis_scales_dict=None
    else:
        yaxis_scales_dict[values["tab_name"]] = calculate_yaxis_scales(
            values["graph_type"], config_g=values
        )
    plot_indicator(
        graph_type=plot_function(values["graph_type"]), 
        config_plot=values,
        yaxis_scales_dict=yaxis_scales_dict
    )    

generate_sidebar(toc)