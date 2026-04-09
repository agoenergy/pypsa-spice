# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Render the landing page for the PyPSA-SPICE-Vis tool."""

import streamlit as st

st.title("Welcome to the PyPSA-SPICE-Vis tool")

st.markdown(
    """
    PyPSA-SPICE-Vis is the interactive visualisation interface for PyPSA-SPICE.
    It helps you inspect scenario inputs, review scenario configuration files,
    and compare model outputs across sectors and years from a single Streamlit app.

    The app is intended to run inside the PyPSA-SPICE repository so it can read
    project data, scenario configuration files, and generated result files directly.
    Use it when you want a faster way to explore a project than opening YAML, CSV,
    and result files one by one.
    """
)

st.markdown(
    "Detailed information can be found in: "
    f"[PyPSA-SPICE-Vis documentation]("
    "https://agoenergy.github.io/pypsa-spice/visualisation-tool/pypsa-spice-vis/)"
)

st.divider()
st.subheader("How To navigate the app")
st.markdown(
    """
    1. Select a project folder from the sidebar.
    2. Open the sub-page you want to work on from the top navigation.
    3. Choose the relevant scenario in the sidebar when those control buttons appear.
    4. Review/edit inputs and configuration, or compare model outputs across scenarios.
    """
)

st.divider()
st.subheader("Pages Overview")

info_col, config_col = st.columns(2)

with info_col:
    st.markdown(
        """
        #### Info
        Start here for a quick orientation to the tool.

        This page gives a short introduction to PyPSA-SPICE-Vis, explains how the
        interface is organised, and helps users understand what each section of the app is for.
        """
    )

with config_col:
    st.markdown(
        """
        #### Scenario config
        Edit the selected `scenario_config.yaml` file through a form-based interface.

        This page is used to review and update scenario settings such as snapshots,
        CO2 management, and custom constraints without editing YAML manually.
        """
    )

input_col, output_col = st.columns(2)

with input_col:
    st.markdown(
        """
        #### Input
        Inspect and edit model inputs for the selected sector.

        The Input page is split into two sections:

        - Global input: shared input data such as technology assumptions, profiles,
          and demand data used across scenarios.
        - Scenario specific input: scenario-level input tables such as loads and,
          for power, interconnections.

        The available sector tabs depend on the selected model scope and can include
        Power, Industry, and Transport.
        """
    )

with output_col:
    st.markdown(
        """
        #### Output
        Explore charts generated from PyPSA-SPICE results and compare scenarios.

        The Output page groups results into sector and summary views:

        - Power: electricity system capacities, generation, storage, and related metrics.
        - Industry: industrial energy system outputs for the selected scenario and year.
        - Transport: transport demand and technology results.
        - Emissions: emissions trends and comparisons across scenarios.
        - Costs: system cost breakdowns and scenario cost comparisons.
        """
    )

st.divider()
st.subheader("What Appears In The Sidebar")
st.markdown(
    """
    The sidebar changes with the page you open:

    - All pages show the project selector.
    - Scenario config and Input show an input scenario selector and sector selector.
    - Output shows one or two output scenario selectors, a sector selector, and then loads the available years for comparison.
    """
)

st.caption(
    "The exact tabs and charts shown in Input and Output depend on the sectors and result files available in the selected project."
)
