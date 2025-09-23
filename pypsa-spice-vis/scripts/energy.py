# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

import os

import streamlit as st
import yaml

from scripts.helpers import *

st.title(':material/bolt: Energy')

with open(
    os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"), "r"
) as file:
    config = yaml.safe_load(file)["energy"]

toc = [config[item]["name"] for item in config]
yaxis_scales_dict = {}

for item, values in config.items():
    yaxis_scales_dict = {}
    if "with_filter" in values["graph_type"] or "hourly" in values["graph_type"]:
        yaxis_scales_dict = None
    elif "sankey_diagram" in values["graph_type"]:     
        yaxis_scales_dict = None
    else:
        yaxis_scales_dict[values["tab_name"]] = calculate_yaxis_scales(
            values["graph_type"], config_g=values
        )
    plot_indicator(
        graph_type=plot_function(values["graph_type"]),
        config_plot=values,
        yaxis_scales_dict=yaxis_scales_dict,
    ) 

generate_sidebar(toc)