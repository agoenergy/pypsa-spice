# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Initialize the Streamlit application for PyPSA-SPICE visualization."""

import os
import sys

import streamlit as st
from styles import apply_sidebar_styles, apply_title_styles, use_flexo

from scripts.get_params import GetParams
from scripts.input_pages import input_static, input_timeseries
from scripts.output_pages import (output_costs, output_emissions,
                                  output_industry, output_power,
                                  output_transport)

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)


st.set_page_config(initial_sidebar_state="expanded", layout="wide")

use_flexo()

DEPLOY = False

get_params = GetParams()

current_dir = get_params.streamlit_base_dir
st.session_state.current_dir = current_dir

init_conf = get_params.init_config

st.session_state.input_data_folder_path = init_conf["input_folder_path"]

apply_sidebar_styles()

st.session_state.window_width = get_params.get_window_width(
    st.session_state.get("window_width")
)

with st.sidebar:
    apply_title_styles("Parameters for settings")

    st.session_state.project = st.sidebar.selectbox(
        ":material/globe: Project :",
        options=get_params.get_project_folder_list(init_conf["data_folder_path"]),
        index=0,
    )

    st.session_state.result_path = os.path.join(init_conf["results_folder_path"])

    st.session_state.sce1 = st.sidebar.selectbox(
        ":material/looks_one: Scenario 1:",
        options=get_params.get_output_scenario_list(
            selected_project_name=st.session_state.project
        ),
        index=0,
    )

    scenario_list = get_params.get_output_scenario_list(
        selected_project_name=st.session_state.project
    )
    if len(scenario_list) == 1:
        st.session_state.sce2 = ""
    else:
        scenario_list.append("None")
        st.session_state.sce2 = st.sidebar.selectbox(
            ":material/looks_two: Scenario 2:", options=scenario_list, index=1
        )
        if st.session_state.sce2 == "None":
            st.session_state.sce2 = ""

    st.session_state.sector = st.sidebar.selectbox(
        ":material/crossword: Sector:",
        options=get_params.get_sector_list(st.session_state.sce1),
        format_func=get_params.get_sector_display_name,
        index=0,
    )

    if st.session_state.sce1 == st.session_state.sce2:
        st.sidebar.error("⚠️ The two scenarios should not be the same!")
        st.stop()

try:
    st.session_state.sce1_years = get_params.get_year_list(
        st.session_state.sce1, st.session_state.sector
    )
    if st.session_state.sce2:
        st.session_state.sce2_years = get_params.get_year_list(
            st.session_state.sce2, st.session_state.sector
        )
except FileNotFoundError as e:
    st.write(e)


def input_main():
    """Render all Input views in one subpage."""
    st.title(":material/input: Input")

    input_sections = ["Static", "Timeseries"]
    selected_input_section = st.segmented_control(
        "Input section",
        options=input_sections,
        default=input_sections[0],
        selection_mode="single",
        label_visibility="collapsed",
        key="selected_input_section",
    )

    if selected_input_section is None:
        st.error("Please select a tab to view and edit the input data.")
        st.stop()

    if selected_input_section == "Static":
        input_static.main(get_params)
    else:
        input_timeseries.main(get_params)


def output_main():
    """Render all Output views in one subpage."""
    st.title(":material/monitoring: Output")

    output_tabs = ["Power", "Emissions", "Costs"]
    include_industry = "i" in st.session_state.sector
    include_transport = "t" in st.session_state.sector

    if include_industry:
        output_tabs.insert(1, "Industry")
    if include_transport:
        output_tabs.insert(2 if include_industry else 1, "Transport")

    selected_output_tab = st.segmented_control(
        "Output section",
        options=output_tabs,
        default=output_tabs[0],
        selection_mode="single",
        label_visibility="collapsed",
        key="selected_output_tab",
    )

    if selected_output_tab is None:
        st.error("Please select a tab to view the output charts.")
        st.stop()

    output_page_mapping = {
        "Power": output_power.main,
        "Industry": output_industry.main,
        "Transport": output_transport.main,
        "Emissions": output_emissions.main,
        "Costs": output_costs.main,
    }
    output_page_mapping[selected_output_tab]()


in_page = st.Page(input_main, title="Input", icon=":material/input:")
out_page = st.Page(output_main, title="Output", icon=":material/monitoring:")

pg = st.navigation([in_page, out_page], position="sidebar")
pg.run()
