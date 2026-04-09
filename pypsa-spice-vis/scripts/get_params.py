# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Get various parameters related to the visual app and data structure."""

import os
import time

import pandas as pd
import streamlit as st
from streamlit_js_eval import streamlit_js_eval


class GetParams:
    """Functions that handle retrieval of app-related params."""

    SECTOR_DISPLAY_NAMES = {
        "p-i-t": "Power+Industry+Transport",
        "p": "Power only",
        "p-i": "Power+Industry",
        "p-t": "Power+Transport",
    }

    def __init__(self, base_config: dict):
        self.base_config = base_config

        # data/data_folder_name
        data_folder_path = os.path.join(
            "data/", self.base_config["path_configs"]["data_folder_name"]
        )
        default_project_name = self.base_config["path_configs"]["project_name"]

        self.base_config["data_folder_path"] = data_folder_path
        # data/data_folder_name/project_name/results
        self.base_config["results_folder_path"] = os.path.join(
            data_folder_path, default_project_name, "results"
        )

    def get_input_scenario_list(self, selected_project_name: str) -> list[str]:
        """Get the list of input scenarios.

        From a given project within the input folder.

        Parameters
        ----------
        selected_project_name : str
          Name of the project folder selected by users in the app.

        Returns
        -------
        list[str]
          The list of scenarios for this project.
        """
        data_folder_path = os.path.join(
            self.base_config["data_folder_path"], selected_project_name, "input"
        )

        if not os.path.exists(data_folder_path):
            raise FileNotFoundError(f"folder not found: {data_folder_path}")

        scenario_list = [
            scenario
            for scenario in os.listdir(data_folder_path)
            if scenario not in [".DS_Store"] and scenario != "global_input"
        ]

        # Make default scenario the first option in the list if present
        default_scenarios = [
            self.base_config["path_configs"]["input_scenario_name"],
            "",
        ]
        for sce in default_scenarios:
            if sce in scenario_list:
                scenario_list.insert(0, scenario_list.pop(scenario_list.index(sce)))

        return scenario_list

    def get_output_scenario_list(self, selected_project_name: str) -> list[str]:
        """Get the list of output scenarios.

        From a given project within the results folder.

        Parameters
        ----------
        selected_project_name : str
          Name of the project folder selected by users in the app.

        Returns
        -------
        list[str]
          The list of scenarios for this project.
        """
        data_folder_path = os.path.join(
            self.base_config["data_folder_path"], selected_project_name, "results"
        )

        if not os.path.exists(data_folder_path):
            raise FileNotFoundError(f"folder not found: {data_folder_path}")

        scenario_list = [
            scenario
            for scenario in os.listdir(data_folder_path)
            if scenario not in [".DS_Store"]
        ]

        # Make default scenario the first option in the list if present
        default_scenarios = [
            self.base_config["path_configs"]["output_scenario_name"],
            "",
        ]
        for sce in default_scenarios:
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
        results_path = self.base_config["results_folder_path"]
        csv_folder_path = os.path.join(results_path, scenario, "csvs")

        if not os.path.exists(csv_folder_path):
            raise FileNotFoundError(f"folder not found: {csv_folder_path}")

        sector_list = [
            sector
            for sector in os.listdir(csv_folder_path)
            if sector not in [".DS_Store"]
        ]

        return sector_list

    def get_sector_display_name(self, sector_code: str) -> str:
        """Get display name for a sector code.

        Parameters
        ----------
        sector_code : str
            Sector code as used in folder names.

        Returns
        -------
        str
            Human-readable display name for the UI.
        """
        return self.SECTOR_DISPLAY_NAMES.get(sector_code, sector_code)

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
        results_path = self.base_config["results_folder_path"]
        sector_folder_path = os.path.join(results_path, scenario, "csvs", sector)

        if not os.path.exists(sector_folder_path):
            raise FileNotFoundError(f"folder not found: {sector_folder_path}")

        years_list = sorted(
            int(year) for year in os.listdir(sector_folder_path) if year.isdigit()
        )

        return years_list

    def get_country_list(self, df: pd.DataFrame) -> list[str]:
        """Get the list of countries from a dataframe if it has a "country" column.

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
        """Get the list of technologies to display in the input UI.

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

        Since streamlit_js_eval is asynchronous and may return None on first pass if
        the javascript has not executed in the browser yet, we use a logic that tries
        for up to five attempts with a short delay in between, before falling back to
        a default width.

        This width is used to set window_width in the session state, in order to set
        the legend position and orientation later on (below the graph and horizontal
        for narrow widths + two scenario cases).

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
        default_width = 1200

        if current_width is not None and current_width > 0:
            return current_width

        for attempt in range(max_attempts):
            try:
                result = streamlit_js_eval(
                    js_expressions="window.innerWidth",
                    key=f"SCR_{attempt}",
                    want_output=True,
                )

                if (
                    result is not None
                    and isinstance(result, (int, float))
                    and result > 0
                ):
                    return int(result)  # Width to set in session_state later

                # If previous attempt failed, wait a tiny bit before trying again
                if attempt < max_attempts - 1:
                    time.sleep(delay)

            except Exception as exc:
                st.write(f"Attempt {attempt + 1} failed: {exc}")
                continue

        return default_width

    def get_project_folder_list(self) -> list[str]:
        """Get a list of project subfolders in the data folder.

        Returns
        -------
        list[str]
          The list of project subfolders.
        """
        folder_path = self.base_config["data_folder_path"]

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
        if self.base_config["path_configs"]["project_name"] in project_folders:
            project_folders.remove(self.base_config["path_configs"]["project_name"])
            project_folders.insert(0, self.base_config["path_configs"]["project_name"])

        return project_folders
