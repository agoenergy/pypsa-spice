# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Initialize the Streamlit application for PyPSA-SPICE visualization."""

import os
import runpy
import sys

import streamlit as st
import yaml
from styles import use_flexo

from scripts.get_params import GetParams


@st.cache_data
def load_yaml(path: str) -> dict:
    """Load a YAML file."""
    with open(path, encoding="utf-8") as file:
        return yaml.safe_load(file)


def render_script_page(relative_script_path: str) -> None:
    """Execute a page script in the current Streamlit page."""
    script_path = os.path.join(
        st.session_state.streamlit_base_dir, relative_script_path
    )
    runpy.run_path(script_path, run_name="__main__")


def get_sector_tabs(base_tabs: list[str], sector: str) -> list[str]:
    """Return tabs extended by enabled sectors."""
    tabs = base_tabs.copy()
    if "i" in sector:
        tabs.append("Industry")
    if "t" in sector:
        tabs.append("Transport")
    return tabs


def scenario_config_main() -> None:
    """Render the scenario configuration editor page."""
    st.session_state.input_sce1 = st.session_state.get("input_sce1", "")

    # Path to the selected scenario config file.
    scenario_config_path = os.path.join(
        st.session_state.input_path,
        st.session_state.input_sce1,
        "scenario_config.yaml",
    )
    st.session_state.scenario_config_path = scenario_config_path

    # Load the config file for the selected scenario.
    st.session_state.scenario_config = load_yaml(scenario_config_path)

    render_script_page("scripts/scenario_configuration_pages/scenario_config.py")


def input_main() -> None:
    """Render all input views in a single page."""
    st.title(":material/input: Input")

    st.session_state.input_sce1 = st.session_state.get("input_sce1", "")

    # Build the sector tabs dynamically based on enabled sectors.
    input_sector_tabs = get_sector_tabs(["Power"], st.session_state.sector)

    selected_input_sector = st.segmented_control(
        "Input sector",
        options=input_sector_tabs,
        default=input_sector_tabs[0],
        selection_mode="single",
        label_visibility="collapsed",
        key="selected_input_sector",
    )

    if selected_input_sector is None:
        st.error("Please select a sector to view and edit input data.")
        st.stop()

    input_page_mapping = {
        "Global input": {
            "Power": "scripts/input_pages/input_power_global.py",
            "Industry": "scripts/input_pages/input_industry_global.py",
            "Transport": "scripts/input_pages/input_transport_global.py",
        },
        "Scenario specific input": {
            "Power": "scripts/input_pages/input_power_scenario.py",
            "Industry": "scripts/input_pages/input_industry_scenario.py",
            "Transport": "scripts/input_pages/input_transport_scenario.py",
        },
    }

    selected_input_tab = st.segmented_control(
        "Input section",
        options=list(input_page_mapping),
        default=next(iter(input_page_mapping)),
        selection_mode="single",
        label_visibility="collapsed",
        key="selected_input_tab",
    )

    if selected_input_tab is None:
        st.error("Please select a tab to view and edit the input data.")
        st.stop()

    render_script_page(input_page_mapping[selected_input_tab][selected_input_sector])


def output_main() -> None:
    """Render all output views in a single page."""
    # Build output tabs dynamically from the active sector selection.
    output_tabs = ["Power"]
    if "i" in st.session_state.sector:
        output_tabs.append("Industry")
    if "t" in st.session_state.sector:
        output_tabs.append("Transport")
    output_tabs.extend(["Emissions", "Costs"])

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

    # Map each output tab to the script that renders it.
    output_page_mapping = {
        "Power": "scripts/output_pages/output_power.py",
        "Industry": "scripts/output_pages/output_industry.py",
        "Transport": "scripts/output_pages/output_transport.py",
        "Emissions": "scripts/output_pages/output_emissions.py",
        "Costs": "scripts/output_pages/output_costs.py",
    }
    render_script_page(output_page_mapping[selected_output_tab])


def initialize_paths() -> None:
    """Load app configs into session state and initialize input/output paths."""
    # Resolve the directory of the Streamlit app entry point.
    base_dir = os.path.dirname(
        os.path.abspath(sys.modules["__main__"].__file__)  # pylint: disable=no-member
    )
    st.session_state.streamlit_base_dir = base_dir

    try:
        # Load all YAML config files needed by the app.
        settings_dir = os.path.join(base_dir, "setting")
        st.session_state.input_config = load_yaml(
            os.path.join(settings_dir, "input_settings.yaml")
        )
        st.session_state.output_config = load_yaml(
            os.path.join(settings_dir, "graph_settings.yaml")
        )
        st.session_state.base_config = load_yaml(
            os.path.join(base_dir, "..", "base_config.yaml")
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            "Please ensure your working directory is at the pypsa-spice root level."
        ) from exc
    except yaml.YAMLError as exc:
        raise ValueError(f"Invalid YAML configuration: {exc}") from exc

    # Build path to project folder.
    path_configs = st.session_state.base_config["path_configs"]
    data_root = os.path.join(
        "data",
        path_configs["data_folder_name"],
        path_configs["project_name"],
    )

    # Store input/ output paths in session state.
    st.session_state.data_folder_name = path_configs["data_folder_name"]
    st.session_state.input_data_folder_path = path_configs["input_scenario_name"]
    st.session_state.result_path = os.path.join(data_root, "results")


if __name__ == "__main__":
    # Load configuration and initialize shared paths before rendering the UI.
    initialize_paths()

    # Set up Streamlit page layout and global styling.
    st.set_page_config(initial_sidebar_state="expanded", layout="wide")
    base_dir = os.path.dirname(os.path.abspath(__file__))
    logo_path = os.path.join(base_dir, "design", "pypsa-logo_rgb.png")

    st.logo(logo_path, size="large")

    use_flexo()

    # Define the main navigation pages shown in the sidebar.
    info_page = st.Page("scripts/info.py", title="Info", icon=":material/info:")
    config_page = st.Page(
        scenario_config_main, title="Scenario config", icon=":material/settings:"
    )
    in_page = st.Page(input_main, title="Input", icon=":material/input:")
    out_page = st.Page(output_main, title="Output", icon=":material/monitoring:")
    pages = st.navigation([info_page, config_page, in_page, out_page], position="top")
    current_page_title = getattr(pages, "title", "")

    # Helper object for reading project/scenario/sector metadata.
    get_params = GetParams(st.session_state.base_config)

    # Persist current browser window width in session state.
    st.session_state.window_width = get_params.get_window_width(
        st.session_state.get("window_width")
    )

    with st.sidebar:

        # Select project
        st.session_state.project = st.selectbox(
            ":material/globe: Project :",
            options=get_params.get_project_folder_list(),
        )

        input_scenario_list = get_params.get_input_scenario_list(
            selected_project_name=st.session_state.project
        )
        show_input_scenario: bool = current_page_title in {"Input", "Scenario config"}

        if show_input_scenario:
            # Input pages use an input scenario and sector selection.
            st.session_state.input_sce1 = st.selectbox(
                ":material/looks_one: Input scenario 1:",
                options=input_scenario_list,
                index=0,
            )
            st.session_state.sector = st.selectbox(
                ":material/crossword: Input Sector:",
                options=st.session_state.base_config["base_configs"]["sector"],
                format_func=get_params.get_sector_display_name,
                index=0,
            )
            data_root = os.path.join(
                "data",
                st.session_state.data_folder_name,
                st.session_state.project,
            )
            st.session_state.input_path = os.path.join(data_root, "input")
        else:
            st.session_state.input_sce1 = ""

        output_scenario_list = get_params.get_output_scenario_list(
            selected_project_name=st.session_state.project
        )
        show_output_scenario: bool = current_page_title == "Output"

        if show_output_scenario:
            # Output pages allow selecting one or two scenarios for comparison.
            st.session_state.output_sce1 = st.selectbox(
                ":material/looks_one: Output scenario 1:",
                options=output_scenario_list,
                index=0,
            )

            if len(output_scenario_list) == 1:
                st.session_state.output_sce2 = ""
            else:
                scenario_2_options = output_scenario_list + ["None"]
                st.session_state.output_sce2 = st.selectbox(
                    ":material/looks_two: Output scenario 2:",
                    options=scenario_2_options,
                    index=1,
                )
                if st.session_state.output_sce2 == "None":
                    st.session_state.output_sce2 = ""

            # Sector choices for output depend on the selected output scenario.
            st.session_state.sector = st.selectbox(
                ":material/crossword: Output Sector:",
                options=get_params.get_sector_list(st.session_state.output_sce1),
                format_func=get_params.get_sector_display_name,
                index=0,
            )

            # Prevent comparing a scenario with itself.
            if st.session_state.output_sce1 == st.session_state.output_sce2:
                st.error("⚠️ The two scenarios should not be the same!")
                st.stop()

            # Load available years for the selected output scenario(s).
            try:
                st.session_state.output_sce1_years = get_params.get_year_list(
                    st.session_state.output_sce1,
                    st.session_state.sector,
                )
                if st.session_state.output_sce2:
                    st.session_state.output_sce2_years = get_params.get_year_list(
                        st.session_state.output_sce2,
                        st.session_state.sector,
                    )
            except FileNotFoundError as exc:
                st.write(exc)
        else:
            st.session_state.output_sce1 = ""
            st.session_state.output_sce2 = ""

    # Render the selected page.
    pages.run()
