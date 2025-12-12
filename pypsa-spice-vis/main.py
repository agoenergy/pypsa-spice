# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Initialize the Streamlit application for PyPSA-SPICE visualization."""

import os
import sys

import streamlit as st
from styles import apply_sidebar_styles, use_flexo

from scripts.getters import Getters

ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
if ROOT_DIR not in sys.path:
    sys.path.insert(0, ROOT_DIR)


st.set_page_config(initial_sidebar_state="expanded", layout="wide")

use_flexo()

DEPLOY = False

getters = Getters()

current_dir = getters.streamlit_base_dir
st.session_state.current_dir = current_dir

init_conf = getters.init_config

# Initialize input_data_folder_path in session state
st.session_state.input_data_folder_path = init_conf["input_folder_path"]

# st.logo( # noqa: E800
#     os.path.join( # noqa: E800
#       st.session_state.current_dir, "design/pypsa-spice-long.png" # noqa: E800
#     ),  # noqa: E800
#     icon_image=os.path.join( # noqa: E800
#       st.session_state.current_dir, "design/agora-icon.png" # noqa: E800
#     ),  # noqa: E800
# ) # noqa: E800

apply_sidebar_styles()

# Set window_width in session state
st.session_state.window_width = getters.get_window_width(
    st.session_state.get("window_width")
)

with st.sidebar:
    st.markdown(
        """<p style='font-size: 1.2em; font-weight: 600; margin-bottom: 12px;'>
                Parameters for settings
                </p>
                """,
        unsafe_allow_html=True,
    )

    # Set project in session state
    st.session_state.project = st.sidebar.selectbox(
        ":material/globe: Project :",
        options=getters.get_project_folder_list(init_conf["data_folder_path"]),
        index=0,
    )

    # Set result_path in session state (note this includes the project dir)
    st.session_state.result_path = os.path.join(init_conf["results_folder_path"])

    # Set sce1 name in session state
    st.session_state.sce1 = st.sidebar.selectbox(
        ":material/looks_one: Scenario 1:",
        options=getters.get_output_scenario_list(
            selected_project_name=st.session_state.project
        ),
        index=0,
    )

    # Set sce2 name in session state
    scenario_list = getters.get_output_scenario_list(
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

    # Set sector in session state
    st.session_state.sector = st.sidebar.selectbox(
        ":material/crossword: Sector:",  # init_conf["sector"]
        options=getters.get_sector_list(st.session_state.sce1),
        index=0,
    )

    if st.session_state.sce1 == st.session_state.sce2:
        st.sidebar.error("⚠️ The two scenarios should not be the same!")
        st.stop()

# Set year in session state
try:
    st.session_state.sce1_years = getters.get_year_list(
        st.session_state.sce1, st.session_state.sector
    )
    if st.session_state.sce2:
        st.session_state.sce2_years = getters.get_year_list(
            st.session_state.sce2, st.session_state.sector
        )
except FileNotFoundError as e:
    st.write(e)

p_page = st.Page(
    "scripts/output_pages/output_power.py", 
    title="Power", 
    icon=":material/bolt:"

)
ene_page = st.Page(
    "scripts/output_pages/output_energy.py", 
    title="Energy", 
    icon=':material/captive_portal:'
)
i_page = st.Page(
    "scripts/output_pages/output_industry.py",
    title="Industry",
    icon=":material/construction:",
)
t_page = st.Page(
    "scripts/output_pages/output_transport.py",
    title="Transport",
    icon=":material/directions_car:",
)
e_page = st.Page(
    "scripts/output_pages/output_emissions.py",
    title="Emissions",
    icon=":material/thermostat:",
)
c_page = st.Page(
    "scripts/output_pages/output_costs.py",
    title="Costs",
    icon=":material/attach_money:",
)
info_page = st.Page("scripts/info.py", title="Info", icon=":material/info:")
in_p_supply_page = st.Page(
    "scripts/input_pages/input_power_supply.py",
    title="Power - Supply",
    icon=":material/bolt:",
)

in_p_demand_page = st.Page(
    "scripts/input_pages/input_power_demand.py",
    title="Power - Demand",
    icon=":material/bolt:",
)

pages_list = {
    "General": [info_page],  # in_setting_page
    "Input": [in_p_supply_page, in_p_demand_page],
    "Results":[p_page, ene_page, i_page, t_page, e_page, c_page],
}

if "i" not in st.session_state.sector:
    pages_list["Results"].remove(i_page)
if "t" not in st.session_state.sector:
    pages_list["Results"].remove(t_page)

pg = st.navigation(pages_list, position="sidebar")
pg.run()
