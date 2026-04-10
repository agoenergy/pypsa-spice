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
    """
)

st.markdown(
    "Detailed information can be found in: "
    "[PyPSA-SPICE-Vis documentation]("
    "https://agoenergy.github.io/pypsa-spice/visualisation-tool/pypsa-spice-vis/)"
)

st.divider()
st.subheader("How to use the app")
st.markdown(
    """
    1. Select a project folder from the sidebar.
    2. Open the sub-page you want to work on from the top navigation.
    3. Choose the relevant scenario in the sidebar when those control options appear.
    4. Review or edit inputs and configuration, or compare model outputs across
       scenarios.
    """
)

st.divider()
st.subheader("Pages overview (top navigation)")

info_col, config_col = st.columns(2)

with info_col:
    st.markdown(
        """
        #### :material/info: Info
        Start here for a quick orientation to the tool and a short introduction to \
        PyPSA-SPICE-Vis.
        """
    )

with config_col:
    st.markdown(
        """
        #### :material/settings: Scenario config
        Edit the selected `scenario_config.yaml` file through this interface to \
        review and update scenario settings such as snapshots, CO2 management, and \
        custom constraints.
        """
    )

input_col, output_col = st.columns(2)

with input_col:
    st.markdown(
        """
        #### :material/input: Input
        Inspect and edit model inputs for the selected sector in two sections:

        - Global input: shared input data such as technology assumptions, profiles,
          and demand data used across scenarios.
        - Scenario-specific input: scenario-level input tables such as loads and,
          for power, interconnections.

        The available sector tabs depend on the selected model scope and can include
        Power, Industry, and Transport.
        """
    )

with output_col:
    st.markdown(
        """
        #### :material/monitoring: Output
        Explore charts generated from PyPSA-SPICE results across the following sectors:

        - Power: electricity system capacities, generation, storage, etc.
        - Industry: industrial energy system outputs.
        - Transport: transport demand and technology results.
        - Emissions: emissions trends.
        - Costs: system cost breakdowns.
        """
    )

st.divider()
st.subheader("What appears in the sidebar")
st.markdown(
    """
    The sidebar changes with the page you open:

    - All pages show the project selector.
    - Scenario config and Input show an input scenario selector and the sector selector.
    - Output shows one or two output scenario selectors, a sector selector, and then \
    loads the available years for comparison. The number of selected scenarios
    can be defined by the user.
    """
)
