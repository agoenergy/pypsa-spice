# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import yaml
import os
import streamlit as st
from scripts.output_helpers import (
    calculate_yaxis_scales, plot_indicator, plot_function, generate_sidebar
)

st.title(":material/thermostat: Emissions")

with open(os.path.join(st.session_state.current_dir,
                       "setting/graph_settings.yaml"), "r") as file:
    config = yaml.safe_load(file)["emissions"]

toc = []
yaxis_scales_dict = {}

for item, values in config.items():
    if values["incl_sector"]=="all" or values["incl_sector"] in st.session_state.sector:
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
        toc.append(values["name"])
    else:
        pass


generate_sidebar(toc)