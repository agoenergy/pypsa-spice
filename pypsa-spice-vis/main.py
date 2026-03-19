# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Initialize the Streamlit application for PyPSA-SPICE visualization."""

import os
import runpy
import sys

import streamlit as st
import yaml
from styles import apply_sidebar_styles, apply_title_styles, use_flexo

from scripts.get_params import GetParams

st.set_page_config(initial_sidebar_state="expanded", layout="wide")

use_flexo()

# Directory of app's entry point file (pypsa-spice/pypsa-spice-vis)
st.session_state.streamlit_base_dir = os.path.dirname(
    os.path.abspath(sys.modules["__main__"].__file__)  # pylint: disable=no-member
)

# Open all config files for the app session_state
try:
    with open(
        os.path.join(
            st.session_state.streamlit_base_dir, "setting/input_settings.yaml"
        ),
        encoding="utf-8",
    ) as file:
        st.session_state.input_config = yaml.safe_load(file)

    with open(
        os.path.join(
            st.session_state.streamlit_base_dir, "setting/graph_settings.yaml"
        ),
        encoding="utf-8",
    ) as file:
        st.session_state.output_config = yaml.safe_load(file)

    with open(
        os.path.join(st.session_state.streamlit_base_dir, "../base_config.yaml"),
        encoding="utf-8",
    ) as file:
        st.session_state.base_config = yaml.safe_load(file)

except FileNotFoundError as exc:
    raise FileNotFoundError(
        "Please ensure your working directory is at the pypsa-spice root level."
    ) from exc
except Exception as exc:
    raise Exception(
        f"Error loading configuration file: {str(exc)}. "
        + "Please ensure that the base_config.yaml file exists."
    ) from exc

st.session_state.input_data_folder_path = st.session_state.base_config["path_configs"][
    "input_scenario_name"
]

st.session_state.input_path = os.path.join(
    "data",
    st.session_state.base_config["path_configs"]["data_folder_name"],
    st.session_state.base_config["path_configs"]["project_name"],
    "input",
)

st.session_state.result_path = os.path.join(
    "data",
    st.session_state.base_config["path_configs"]["data_folder_name"],
    st.session_state.base_config["path_configs"]["project_name"],
    "results",
)

# st.logo( # noqa: E800
#     os.path.join( # noqa: E800
#       st.session_state.streamlit_base_dir, "design/pypsa-spice-long.png" # noqa: E800
#     ),  # noqa: E800
#     icon_image=os.path.join( # noqa: E800
#       st.session_state.streamlit_base_dir, "design/agora-icon.png" # noqa: E800
#     ),  # noqa: E800
# ) # noqa: E800

apply_sidebar_styles()


get_params = GetParams(st.session_state.base_config)

# Set window_width in session state
st.session_state.window_width = get_params.get_window_width(
    st.session_state.get("window_width")
)


def render_script_page(relative_script_path: str):
    """Execute a page script so it renders in the current Streamlit page."""
    script_path = os.path.join(
        st.session_state.streamlit_base_dir, relative_script_path
    )
    runpy.run_path(script_path, run_name="__main__")


def input_main():
    """Render all Input views in one page with tab-like controls."""
    st.title(":material/input: Input")

    input_page_mapping = {
        "Static": "scripts/input_pages/input_static.py",
        "Timeseries": "scripts/input_pages/input_timeseries.py",
    }
    selected_input_tab = st.segmented_control(
        "Input section",
        options=list(input_page_mapping.keys()),
        default=list(input_page_mapping.keys())[0],
        selection_mode="single",
        label_visibility="collapsed",
        key="selected_input_tab",
    )

    if selected_input_tab is None:
        st.error("Please select a tab to view and edit the input data.")
        st.stop()

    render_script_page(input_page_mapping[selected_input_tab])


def output_main():
    """Render all Output views in one page with tab-like controls."""
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
        "Power": "scripts/output_pages/output_power.py",
        "Industry": "scripts/output_pages/output_industry.py",
        "Transport": "scripts/output_pages/output_transport.py",
        "Emissions": "scripts/output_pages/output_emissions.py",
        "Costs": "scripts/output_pages/output_costs.py",
    }
    render_script_page(output_page_mapping[selected_output_tab])


with st.sidebar:
    apply_title_styles("Parameters for settings")

    # Set project in session state
    st.session_state.project = st.sidebar.selectbox(
        ":material/globe: Project :", options=get_params.get_project_folder_list()
    )

    # Set sce1 name in session state
    st.session_state.sce1 = st.sidebar.selectbox(
        ":material/looks_one: Scenario 1:",
        options=get_params.get_output_scenario_list(
            selected_project_name=st.session_state.project
        ),
        index=0,
    )

    # Set sce2 name in session state
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

    # Set sector in session state
    st.session_state.sector = st.sidebar.selectbox(
        ":material/crossword: Sector:",  # base_config["sector"]
        options=get_params.get_sector_list(st.session_state.sce1),
        format_func=get_params.get_sector_display_name,
        index=0,
    )

    if st.session_state.sce1 == st.session_state.sce2:
        st.sidebar.error("⚠️ The two scenarios should not be the same!")
        st.stop()

# Set year in session state
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


info_page = st.Page("scripts/info.py", title="Info", icon=":material/info:")
in_page = st.Page(input_main, title="Input", icon=":material/input:")
out_page = st.Page(output_main, title="Output", icon=":material/monitoring:")

pg = st.navigation([info_page, in_page, out_page], position="sidebar")
pg.run()
