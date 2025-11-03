# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import sys
import os
import yaml
import time
import pandas as pd
import streamlit as st
from streamlit_js_eval import streamlit_js_eval


class Getters:
    """Functions that handle retrieval of app-related params"""

    def __init__(self):
        # Where app was launched from
        self.working_dir = os.path.basename(os.getcwd())

        # Directory of app's entry point file
        self.entry_dir = os.path.dirname(
            os.path.abspath(sys.modules["__main__"].__file__)
        )

        # Rel path to initial config file - modify as necessary if DEPLOY is true
        self.config_path = "../base_config.yaml"

        with open(os.path.join(self.entry_dir, self.config_path), "r") as file:
            self.init_conf = yaml.safe_load(file)

        if "vis" in str(self.working_dir):
            data_path = "../data/"
        else:
            data_path = "data/"

        self.init_conf["project_folder_path"] = (
            data_path + self.init_conf["path_configs"]["data_folder_name"] + "/"
        )

        self.init_conf["input_data_folder_path"] = (
            self.init_conf["project_folder_path"]
            + "/"
            + self.init_conf["path_configs"]["project_name"]
            + "/input/"
        )

        self.init_conf["data_folder_path"] = (
            self.init_conf["project_folder_path"]
            + "/"
            + self.init_conf["path_configs"]["project_name"]
            + "/results/"
        )

    def get_project_folder_list(self, folder_path: str) -> list[str]:
        """Get a list of project subfolders in a given folder.

        Parameters
        ----------
        folder_path : str
          The path to the input folder to look within.

        Returns
        -------
        list[str]
          The list of project subfolders.
        """
        if not os.path.exists(folder_path):
            raise FileNotFoundError(f"folder not found: {folder_path}")

        project_folders = [
            f
            for f in os.listdir(folder_path)
            if os.path.isdir(os.path.join(folder_path, f))
        ]

        # Remove hidden files and folders
        project_folders = [f for f in project_folders if not f.startswith(".")]

        # Make default project the first option in the list if present
        if self.init_conf["path_configs"]["project_name"] in project_folders:
            project_folders.remove(self.init_conf["path_configs"]["project_name"])
            project_folders.insert(0, self.init_conf["path_configs"]["project_name"])

        return project_folders

    def get_input_scenario_list(self) -> list[str]:
        """Get the list of input scenarios from a given project within the input/ folder.

        Parameters
        ----------
        project_dir : str
          Name of the project folder to look within for the scenarios.

        Returns
        -------
        list[str]
          The list of scenarios for this project.
        """
        data_folder_path = self.init_conf["input_data_folder_path"]

        if not os.path.exists(data_folder_path):
            raise FileNotFoundError(f"folder not found: {data_folder_path}")

        scenario_list = [
            scenario
            for scenario in os.listdir(data_folder_path)
            if scenario not in [".DS_Store"] and scenario != "global_input"
        ]

        # Make default scenario the first option in the list if present
        for sce in (self.init_conf["path_configs"]["input_scenario_name"], ""):
            if sce in scenario_list:
                scenario_list.insert(0, scenario_list.pop(scenario_list.index(sce)))

        return scenario_list

    def get_output_scenario_list(self) -> list[str]:
        """Get the list of output scenarios from a given project within the results/ folder.

        Parameters
        ----------
        project_dir : str
          Name of the project folder to look within for the scenarios.

        Returns
        -------
        list[str]
          The list of scenarios for this project.
        """
        data_folder_path = self.init_conf["data_folder_path"]

        if not os.path.exists(data_folder_path):
            raise FileNotFoundError(f"folder not found: {data_folder_path}")

        scenario_list = [
            scenario
            for scenario in os.listdir(data_folder_path)
            if scenario not in [".DS_Store"]
        ]

        # Make default scenario the first option in the list if present
        for sce in (self.init_conf["path_configs"]["output_scenario_name"], ""):
            if sce in scenario_list:
                scenario_list.insert(0, scenario_list.pop(scenario_list.index(sce)))

        return scenario_list

    def get_sector_list(self, scenario: str) -> list[str]:
        """Get the list of sectors from the scenario/ folder in a given project.

        Parameters
        ----------
        project_dir : str
          Name of the project folder in the results folder.
        scenario : str
          Scenario name.

        Returns
        -------
        list[str]
          The list of sectors for this scenario.
        """
        results_path = self.init_conf["data_folder_path"]
        csv_folder_path = os.path.join(results_path, scenario, "csvs")

        if not os.path.exists(csv_folder_path):
            raise FileNotFoundError(f"folder not found: {csv_folder_path}")

        sector_list = [
            sector
            for sector in os.listdir(csv_folder_path)
            if sector not in [".DS_Store"]
        ]

        return sector_list

    def get_year_list(self, scenario: str, sector: str) -> list[str]:
        """Get the list of years from the scenario/sector/ folder in a given project.

        Parameters
        ----------
        project_path: str
          Path to the project folder
        scenario : str
          Scenario name
        sector : str
          Sector (p/i/t)

        Returns
        -------
        list[str]
          The list of years for this sector.
        """
        results_path = self.init_conf["data_folder_path"]
        sector_folder_path = os.path.join(results_path, scenario, "csvs", sector)

        if not os.path.exists(sector_folder_path):
            raise FileNotFoundError(f"folder not found: {sector_folder_path}")

        years_list = sorted(
            int(year) for year in os.listdir(sector_folder_path) if year.isdigit()
        )

        return years_list

    def get_country_list(self, df: pd.DataFrame) -> list[str]:
        """Get the list of countries from a dataframe if it has a "country" column

        Parameters
        ----------
        df : pd.DataFrame
          The dataframe to check

        Returns
        -------
        list[str]
          The list of countries
        """
        if "country" in df.columns:
            return sorted(df["country"].unique().tolist())

        return []

    def get_mapping_list(self, *dfs: pd.DataFrame) -> list[str]:
        """Get the list of technologies to display in the input UI

        Parameters
        ----------
        dfs : pd.DataFrame
          Any number of dataframes to extract values from the "type" column from.

        Returns
        -------
        list[str]
          The sorted list of technology types.
        """
        type_set = set()

        for df in dfs:
            if "technology" in df.columns:
                type_set |= set(df["technology"].unique())
            if "profile_type" in df.columns:
                type_set |= set(df["profile_type"].unique())

        return sorted(type_set)

    def get_window_width(
        self, current_width: int, max_attempts: int = 5, delay: float = 0.2
    ) -> int:
        """Get the current window width.

        Since streamlit_js_eval is asynchronous and may return None on first pass if the js
        has not executed in the browser yet, we use a logic that tries for up to five
        attempts with a short delay in between, before falling back to a default width.

        This width is used to set window_width in the session state, in order to set the
        legend position and orientation later on (below the graph and horizontal for narrow
        widths + two scenario cases).

        Parameters
        ----------
        current_width : int
          The current window width.
        max_attempts : int, optional
          Maximum number of attempts to try and get the width, by default 5.
        delay : float, optional
          Delay between attempts, by default 0.2.

        Returns
        -------
        int
          The window width in px.
        """

        default_width = 1200  # Default fallback width

        # Return current width if it is valid
        if current_width is not None and current_width > 0:
            return current_width

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
                if (
                    result is not None
                    and isinstance(result, (int, float))
                    and result > 0
                ):
                    return int(result)  # Width to set in session_state later

                # If previous attempt failed, wait a tiny bit before trying again
                if attempt < max_attempts - 1:
                    time.sleep(delay)

            except Exception as e:
                st.write(f"Attempt {attempt + 1} failed: {e}")
                continue

        # Use the default fallback if all attempts fail
        return default_width
