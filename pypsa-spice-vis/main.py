# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

import os
import time

import pandas as pd
import streamlit as st
import yaml
from streamlit_js_eval import streamlit_js_eval
from styles import apply_sidebar_styles, use_flexo

st.set_page_config(initial_sidebar_state="expanded", layout="wide")

use_flexo()

DEPLOY = False

current_dir = os.path.dirname(os.path.abspath(__file__))
st.session_state.current_dir = current_dir

st.logo(
    os.path.join(st.session_state.current_dir, "design/pypsa-spice-long.png"),
    icon_image=os.path.join(st.session_state.current_dir, "design/agora-icon.png"),
)

if DEPLOY:
    with open(
        os.path.join(
            st.session_state.current_dir, "setting/initial_project_01_deploy.yaml"
        ),
        "r",
    ) as file:
        init_conf = yaml.safe_load(file)
else:
    config_path = "setting/initial.yaml"

with open(os.path.join(current_dir, config_path), "r") as file:
    init_conf = yaml.safe_load(file)


if "vis" in str(os.path.basename(os.getcwd())):
    init_conf["data_folder_path"] = "../results/"
else:
    init_conf["data_folder_path"] = "results/"


def get_window_width(max_attempts=5, delay=0.2):
    """
    Get the current window width, and set it in the session state. Since streamlit_js_eval
    is asynchronous and may return None on first pass if the js has not executed in the
    browser yet, we use a logic that tries for up to five attempts with a short delay
    in between, before falling back to a default width.

    The width is used to set the legend orientation later on (horizontal for narrow
    widths and two scenario cases).

    Args:
        max_attempts (int): Maximum no. of attempts to try and get the width (set to 5)
        delay (float): Delay between attempts (set to 0.2)

    """
    default_width = 1200  # Default fallback width

    if "window_width" not in st.session_state:
        st.session_state.window_width = None

    if st.session_state.window_width is not None and st.session_state.window_width > 0:
        return st.session_state.window_width

    for attempt in range(max_attempts):
        try:
            result = streamlit_js_eval(
                js_expressions="window.innerWidth",
                key=f"SCR_{attempt}",  # Use a different key for each attempt
                want_output=True,
            )
            # Note that this does not actually correspond to the true window.innerWidth
            # possibly because of streamlit's iframe context - it seems to be smaller

            # Check for a valid result
            if result is not None and isinstance(result, (int, float)) and result > 0:
                st.session_state.window_width = int(
                    result
                )  # Set width in session state

            # If previous attempt failed, wait a tiny bit before trying again
            if attempt < max_attempts - 1:
                time.sleep(delay)

        except Exception as e:
            st.write(f"Attempt {attempt + 1} failed: {e}")
            continue

    # Use the default fallback if all attempts fail
    st.session_state.window_width = default_width


def get_list_of_projects(result_path):
    """Getting a list of countries from the result folder."""
    if not os.path.exists(result_path):
        raise FileNotFoundError(f"folder not found: {result_path}")

    project_list = [
        project for project in os.listdir(result_path) if project not in [".DS_Store"]
    ]
    # make default project first option in the list if present
    if init_conf["project"] in project_list:
        project_list.remove(init_conf["project"])
        project_list.insert(0, init_conf["project"])
    return project_list


def get_list_of_scenarios(result_path, project):
    """Getting a list of scenarios from the result folder."""
    path = os.path.join(result_path, project)
    if not os.path.exists(path):
        raise FileNotFoundError(f"folder not found: {path}")
    scen_list = [
        scenario for scenario in os.listdir(path) if scenario not in [".DS_Store"]
    ]
    # make default scenario first option in the list if present
    if init_conf["sce1"] in scen_list:
        scen_list.remove(init_conf["sce1"])
        scen_list.insert(0, init_conf["sce1"])
    if init_conf["sce2"] in scen_list:
        scen_list.remove(init_conf["sce2"])
        scen_list.insert(1, init_conf["sce2"])
    # if len(scen_list) == 1:
    #     scen_list.insert(1, "")
    return scen_list


def get_list_of_sector(result_path, project, scenario):
    """Getting a list of sectors from the result/project/scenario folder."""
    path = os.path.join(result_path, project, scenario, "csvs")
    if not os.path.exists(path):
        raise FileNotFoundError(f"folder not found: {path}")
    return [sector for sector in os.listdir(path) if sector not in [".DS_Store"]]


def get_list_of_years(result_path, scenario, sector):
    """Getting a list of years from result/project/scenario folder."""
    path = os.path.join(result_path, scenario, "csvs", sector)
    if not os.path.exists(path):
        raise FileNotFoundError(f"folder not found: {path}")
    return sorted(int(year) for year in os.listdir(path) if year.isdigit())


apply_sidebar_styles()

get_window_width()

with st.sidebar:
    st.markdown(
        """<p style='font-size: 1.2em; font-weight: 600; margin-bottom: 12px;'>
                Parameters for settings
                </p>
                """,
        unsafe_allow_html=True,
    )

    st.session_state.project = st.sidebar.selectbox(
        ":material/globe: Project :",
        options=get_list_of_projects(init_conf["data_folder_path"]),
        index=0,
    )

    st.session_state.sce1 = st.sidebar.selectbox(
        ":material/looks_one: Scenario 1:",
        options=get_list_of_scenarios(
            init_conf["data_folder_path"], st.session_state.project
        ),
        index=0,
    )

    sc_list = get_list_of_scenarios(
        init_conf["data_folder_path"], st.session_state.project
    )
    if len(sc_list) == 1:
        st.session_state.sce2 = ""
    else:
        sc_list.append("None")
        st.session_state.sce2 = st.sidebar.selectbox(
            ":material/looks_two: Scenario 2:", options=sc_list, index=1
        )
        if st.session_state.sce2 == "None":
            st.session_state.sce2 = ""

    st.session_state.sector = st.sidebar.selectbox(
        ":material/crossword: Sector:",  # init_conf["sector"]
        options=get_list_of_sector(
            init_conf["data_folder_path"],
            st.session_state.project,
            st.session_state.sce1,
        ),
        index=0,
    )

    if st.session_state.sce1 == st.session_state.sce2:
        st.sidebar.error("⚠️ The two scenarios should not be the same!")
        st.stop()

st.session_state.result_path = os.path.join(
    init_conf["data_folder_path"], st.session_state.project
)

try:
    st.session_state.sce1_years = get_list_of_years(
        st.session_state.result_path, st.session_state.sce1, st.session_state.sector
    )
    if st.session_state.sce2:
        st.session_state.sce2_years = get_list_of_years(
            st.session_state.result_path, st.session_state.sce2, st.session_state.sector
        )
except FileNotFoundError as e:
    st.write(e)

p_page = st.Page("scripts/power.py", title="Power", icon=":material/bolt:")
i_page = st.Page(
    "scripts/industry.py", title="Industry", icon=":material/construction:"
)
t_page = st.Page(
    "scripts/transport.py", title="Transport", icon=":material/directions_car:"
)
e_page = st.Page(
    "scripts/emissions.py", title="Emissions", icon=":material/thermostat:"
)
c_page = st.Page("scripts/costs.py", title="Costs", icon=":material/attach_money:")
info_page = st.Page("scripts/info.py", title="Info", icon=":material/info:")

pages_list = [p_page, i_page, t_page, e_page, c_page, info_page]
if "i" not in st.session_state.sector:
    pages_list.remove(i_page)
if "t" not in st.session_state.sector:
    pages_list.remove(t_page)

pg = st.navigation(pages_list, position="sidebar")
pg.run()
