# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import yaml
import os
import streamlit as st
from scripts.output_helpers import (
    plot_indicator,
    plot_function,
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
        graph_type=plot_function(values["graph_type"]),
        config_plot=values,
    )


generate_sidebar(table_of_content)
