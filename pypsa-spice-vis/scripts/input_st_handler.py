# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Helper functions for handling Input section in visual app."""

import csv
import datetime as dt
import hashlib
import os
import re
import shutil
import time
from collections.abc import Callable
from dataclasses import dataclass
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from scripts.get_params import GetParams


class DataFrameWidgetsHandler:
    """Handle the DataFrame widget."""

    def __init__(self):
        self.get_params = GetParams()
        self.input_ui_handler = InputUiHandler(
            year_list=self.get_params.init_config["base_configs"]["years"]
        )
        self.csvs_dict = self.input_ui_handler.csvs_dict
        self.data_folder_path = self.get_params.init_config["data_folder_path"]
        project_folders = self.get_params.get_project_folder_list(self.data_folder_path)
        if not project_folders:
            st.error(f"No valid project folders found in {self.data_folder_path}")

        self.input_data_project = project_folders[0]

        self.scenario_input_file_keys = [
            "decomission",
            "fuel_costs",
            "interconnector",
            "load",
            "generator",
            "links",
            "storageunit",
            "store",
        ]

        # Path parts to all the csvs
        self.csv_files_path_parts = {
            "technologies": ["technologies.csv"],
            "availability": ["availability.csv"],
            "demand": ["demand_profile.csv"],
            "pp_costs": ["power_plant_costs.csv"],
            "potentials": ["renewables_technical_potential.csv"],
            "storage_cost": ["storage_costs.csv"],
            "storage_inflows": ["storage_inflows.csv"],
            "decomission": ["power", "decomission_capacity.csv"],
            "fuel_costs": ["power", "fuel_supplies.csv"],
            "interconnector": ["power", "interconnector.csv"],
            "load": ["power", "loads.csv"],
            "generator": ["power", "power_generators.csv"],
            "links": ["power", "power_links.csv"],
            "storageunit": ["power", "storage_capacity.csv"],
            "store": ["power", "storage_energy.csv"],
        }

    def load_all_dfs(self) -> dict:
        """Load all the dataframes.

        Returns
        -------
        dict
            Contains all the loaded dataframes.
        """
        base_input_path = os.path.join(
            self.data_folder_path,
            self.input_data_project,
            "input",
        )

        global_input_path = os.path.join(
            base_input_path,
            "global_input",
        )

        sub_folder = (
            st.session_state["scenario"]
            if "scenario" in st.session_state
            else self.get_params.init_config["path_configs"]["input_scenario_name"]
        )

        scenario_input_path = os.path.join(base_input_path, sub_folder)

        # Update csvs_dict with paths to all csvs
        for key, parts in self.csv_files_path_parts.items():
            if key in self.csvs_dict:
                base_folder = (
                    scenario_input_path
                    if key in self.scenario_input_file_keys
                    else global_input_path
                )
                self.csvs_dict[key].path = os.path.join(base_folder, *parts)

        for key, csv_config in self.csvs_dict.items():
            if not os.path.exists(csv_config.path):
                st.error(f"{key} file not found at {csv_config.path}")
                return {}

        dfs = {
            "tech_df": pd.read_csv(self.csvs_dict["technologies"].path),
            "avail_df": pd.read_csv(self.csvs_dict["availability"].path),
            "demand_df": pd.read_csv(self.csvs_dict["demand"].path),
            "pp_costs_df": pd.read_csv(self.csvs_dict["pp_costs"].path),
            "potentials_df": pd.read_csv(self.csvs_dict["potentials"].path),
            "storage_cost_df": pd.read_csv(self.csvs_dict["storage_cost"].path),
            "storage_inflows_df": pd.read_csv(self.csvs_dict["storage_inflows"].path),
            "decomission_capacity_df": pd.read_csv(self.csvs_dict["decomission"].path),
            "fuel_costs_df": pd.read_csv(self.csvs_dict["fuel_costs"].path),
            "intercon_df": pd.read_csv(self.csvs_dict["interconnector"].path),
            "load_df": pd.read_csv(self.csvs_dict["load"].path),
            "generator_df": pd.read_csv(self.csvs_dict["generator"].path),
            "links_df": pd.read_csv(self.csvs_dict["links"].path),
            "storageunit_df": pd.read_csv(self.csvs_dict["storageunit"].path),
            "store_df": pd.read_csv(self.csvs_dict["store"].path),
        }

        return dfs

    def reload_scenario_dfs(self, dfs: dict, selected_scenario: str) -> dict:
        """Reload the scenario-related dataframes after user has selected a scenario.

        Parameters
        ----------
        dfs : dict
            Dictionary containing all the loaded dataframes.
        selected_scenario : str
            The scenario selected by the user.

        Returns
        -------
        dict
            Dictionary of all loaded dataframes with updated versions for the scenario-
            dependent dataframes.
        """
        # Update csvs_dict with new paths based on selected scenario
        for key in self.scenario_input_file_keys:
            parts = self.csv_files_path_parts[key]
            self.csvs_dict[key].path = os.path.join(
                self.data_folder_path,
                self.input_data_project,
                "input",
                selected_scenario,
                *parts,
            )

        for key in self.scenario_input_file_keys:
            dfs[f"{key}_df"] = pd.read_csv(self.csvs_dict[key].path)

        return dfs


@dataclass
class CsvDictConfig:
    """Configuration for different csvs."""

    identifier: str
    filter_col: str
    title: str
    filter_fn: Callable
    uniform_unit: str | None = None
    empty_df_fn: Callable | None = None
    empty_df_kwargs: dict | None = None
    path: str = ""  # Populate after user selections are made


class InputUiHandler:
    """Handler class for input UI."""

    def __init__(self, year_list: list):
        self.year_list = year_list
        self.csvs_dict = {
            "technologies": CsvDictConfig(
                identifier="tech",
                filter_col="technology",
                title="Technology Parameters",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "availability": CsvDictConfig(
                identifier="hourly_avail",
                filter_col="technology",
                title="Availability",
                uniform_unit="Dimensionless",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "demand": CsvDictConfig(
                identifier="hourly_demand",
                filter_col="profile_type",
                title="Demand Profiles",
                uniform_unit="MW",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "pp_costs": CsvDictConfig(
                identifier="powerplant_costs",
                filter_col="powerplant_type",
                title="Power Plant Costs",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "potentials": CsvDictConfig(
                identifier="potentials",
                filter_col="type",
                title="Renewable Technical Potentials",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "storage_cost": CsvDictConfig(
                identifier="storage_cost",
                filter_col="storage_type",
                title="Storage Costs",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "storage_inflows": CsvDictConfig(
                identifier="hourly_storage_inflows",
                filter_col="technology",
                title="Storage Inflows",
                uniform_unit="MW",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "decomission": CsvDictConfig(
                identifier="decomission",
                filter_col="name",
                title="Decomissioning",
                filter_fn=self.filter_df_decomission,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "fuel_costs": CsvDictConfig(
                identifier="fuel",
                filter_col="carrier",
                title="Fuel Costs",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "interconnector": CsvDictConfig(
                identifier="interconnector",
                filter_col="type",
                title="Interconnector",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "load": CsvDictConfig(
                identifier="load",
                filter_col="profile_type",
                title="Load",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "generator": CsvDictConfig(
                identifier="generator",
                filter_col="type",
                title="Assest Parameters - Generators",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "links": CsvDictConfig(
                identifier="links",
                filter_col="carrier",
                title="Assest Parameters - Links",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "storageunit": CsvDictConfig(
                identifier="storageunit",
                filter_col="carrier",
                title="Assest Parameters - StorageUnits",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
            "store": CsvDictConfig(
                identifier="store",
                filter_col="carrier",
                title="Assest Parameters - Stores",
                filter_fn=self.filter_df_generic,
                empty_df_fn=self.empty_df_message_generic,
            ),
        }

    def create_editable_df(
        self,
        filtered_df: pd.DataFrame,
        edited_df_key: str,
        has_changes_key: str,
    ) -> tuple:
        """Create widget with an editable dataframe.

        Parameters
        ----------
        filtered_df : pd.DataFrame
            The filtered dataframe based on user"s selection.
        edited_df_key : str
            Unique key that references the edited df.
        has_changes_key : str
            _description_

        Returns
        -------
            The edited dataframe and whether it is valid to save.
        """
        to_save = True
        editable_df = filtered_df.replace({np.inf: "inf"})

        editable_cols = filtered_df.select_dtypes(
            include=["number", float, int, "bool"]
        ).columns

        disabled_cols = [
            col
            for col in filtered_df.columns
            if col not in editable_cols and col != "max_supply [MWh/year]"
        ]
        edited_df = st.data_editor(
            editable_df,
            hide_index=True,
            key=edited_df_key,
            disabled=disabled_cols,
            on_change=lambda: st.session_state.update({has_changes_key: True}),
        )

        result_df = edited_df.replace({"inf": np.inf})
        # Validate float columns: only float or "inf" accepted
        for col in filtered_df.select_dtypes(include=[float]).columns:
            try:
                result_df[col] = result_df[col].astype(float)
            except Exception:
                invalid_mask = result_df[col].apply(
                    lambda x: not (isinstance(x, (float, int)) or x == np.inf)
                )

                if invalid_mask.any():
                    st.error(
                        f"Column '{col}' contains invalid entries. "
                        "Only numbers or 'inf' allowed."
                    )
                    to_save = False

                result_df[col] = result_df[col].astype(float, errors="ignore")

        return result_df, to_save

    def create_save_button(
        self,
        filtered_df: pd.DataFrame,
        edited_df: pd.DataFrame,
        has_changes: bool,
        has_changes_key: str,
        save_button_key: str,
        output_file_path: str,
        message_delay: float = 1,
    ):
        """Set up the save button that handles saving of the editable dataframe.

        Parameters
        ----------
        filtered_df : pd.DataFrame
            The filtered dataframe based on the user's selection.
        edited_df: pd.DataFrame
            The edited dataframe with changes made by the user.
        has_changes: bool
            True if df has been edited, False otherwise.
        has_changes_key : str
            Key to session_state var for whether dataframe has changed.
        save_button_key: str
            Button key.
        output_file_path : str
            Path to save the csv to.
        message_delay: float
            How long (s)to display the 'success message' after saving for.
            How long (s)to display the 'success message' after saving for.
        """
        if st.button(
            "Save Changes",
            key=save_button_key,
            type="primary" if has_changes else "secondary",
            disabled=not has_changes,
        ):
            success = True
            # Iterate through each row in the edited dataframe
            for idx in range(len(edited_df)):
                current_index = filtered_df.index[idx]
                for col in filtered_df.columns:
                    if filtered_df[col].iloc[idx] != edited_df[col].iloc[idx]:
                        success &= self.update_csv_file(
                            file_path=output_file_path,
                            row_identifier=str(current_index),
                            column_name=col,
                            new_value=str(edited_df[col].iloc[idx]),
                        )
            if success:
                st.success("Changes saved successfully!")
                st.session_state[has_changes_key] = False
                time.sleep(message_delay)
                st.rerun()
            else:
                st.error("Error saving some changes")

    def set_up_single_tab_widget(
        self,
        csv_dict_key: str,
        input_df: pd.DataFrame,
        selected_types: list,
        input_csv_path: str,
        selected_countries: list = None,
        selected_classes: list | None = None,
        secondary_df: pd.DataFrame | None = None,
    ):
        """Set up the widget with a single tab.

        The tab contains the editable df and save
        button that handles saving of the edited df.

        Parameters
        ----------
        csv_dict_key : str
            Key to the csv in csvs_dict.
        input_df : pd.DataFrame
            The input dataframe.
        selected_type : str
            Technology type selected by the user.
        input_csv_path : str
            Path to the input csv.
        selected_countries : list,
            Country(s) selected by the user in the global country select widget.
        selected_class : list, optional
            Class(es) associated with technologies selected by the user.
        secondary_df : pd.DataFrame, optional
            Additional df if needed (e.g., Availability may need Technologies.csv).
        """
        csv_config = self.csvs_dict[csv_dict_key]
        csv_identifier = csv_config.identifier

        # Create session state var keys based on the current csv
        list_key = self.list_to_key(selected_types)
        edited_df_key = f"{csv_identifier}_editor_{list_key}"
        has_changes_key = f"has_changes_{csv_identifier}_{list_key}"
        save_button_key = f"save_{csv_identifier}_{list_key}"

        filter_col = csv_config.filter_col
        filtered_df = pd.DataFrame()
        no_data_msg = None

        if csv_dict_key == "fuel_costs":
            if "Link" in selected_classes and secondary_df is not None:
                fuels = self.get_fuel_mapping(secondary_df, selected_types)
                filtered_df = csv_config.filter_fn(
                    input_df, filter_col, list(fuels.values())
                )
            else:
                no_data_msg = (
                    f"No fuel costs for technologies: {', '.join(selected_types)}"
                )

        else:
            filtered_df = csv_config.filter_fn(input_df, filter_col, selected_types)

        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        with st.container(border=True):
            st.write(f"### {csv_config.title}")

            path_to_display = os.path.normpath(input_csv_path)
            st.markdown(
                f"<small><i>{path_to_display}</i></small>", unsafe_allow_html=True
            )

            if not filtered_df.empty:
                edited_df, to_save = self.create_editable_df(
                    filtered_df, edited_df_key, has_changes_key
                )

                has_changes = st.session_state.get(has_changes_key, False)

                if to_save:
                    self.create_save_button(
                        filtered_df,
                        edited_df,
                        has_changes,
                        has_changes_key,
                        save_button_key,
                        input_csv_path,
                        message_delay=1,
                    )

            elif no_data_msg:
                st.info(no_data_msg)
            elif csv_config.empty_df_fn:
                # Fallback if defined in csvs_config
                kwargs = dict(csv_config.empty_df_kwargs or {})
                csv_config.empty_df_fn(**kwargs)

    def set_up_double_tab_widget(
        self,
        csv_dict_key: str,
        input_df: pd.DataFrame,
        selected_types: list,
        input_csv_path: str,
        selected_countries: list = None,
        secondary_df: pd.DataFrame | None = None,
    ):
        """Set up the widget with a double tab.

        First tab contains the editable df and save button. Second tab contains the
        visualisation of the selected data.

        Parameters
        ----------
        csv_dict_key : str
            Key to the csv in the csvs_dict.
        input_df : pd.DataFrame
            The input dataframe.
        selected_type : str
            Technology type selected by the user.
        input_csv_path : str
            Path to the input csv.
        selected_countries : list, optional
            Country(s) selected by the user in the global country select widget.
        secondary_df : pd.DataFrame, optional
            Additional df if needed (e.g., Availability may need technologies.csv).
        """
        csv_config = self.csvs_dict[csv_dict_key]

        csv_identifier = csv_config.identifier

        list_key = self.list_to_key(selected_types)
        edited_df_key = f"{csv_identifier}_editor_{list_key}"
        has_changes_key = f"has_changes_{csv_identifier}_{list_key}"
        save_button_key = f"save_{csv_identifier}_{list_key}"

        filter_col = csv_config.filter_col
        filtered_df = csv_config.filter_fn(
            input_df,
            filter_col,
            selected_types,
        )

        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        with st.container(border=True):
            st.write(f"### {csv_config.title}")

            path_to_display = os.path.normpath(input_csv_path)
            if csv_dict_key == "availability" and filtered_df.empty:
                path_to_display = path_to_display.replace(
                    "Availability", "Technologies"
                )

            st.markdown(
                f"<small><i>{path_to_display}</i></small>", unsafe_allow_html=True
            )

            if not filtered_df.empty:
                tab1, tab2 = st.tabs(["Table", "Visualisation"])

                with tab1:
                    edited_df, to_save = self.create_editable_df(
                        filtered_df,
                        edited_df_key,
                        has_changes_key,
                    )

                    has_changes = st.session_state.get(has_changes_key, False)

                with tab2:
                    self.visualise_data(filtered_df, csv_config)

                if to_save:
                    self.create_save_button(
                        filtered_df,
                        edited_df,
                        has_changes,
                        has_changes_key,
                        save_button_key,
                        input_csv_path,
                        message_delay=1,
                    )
            else:
                kwargs = dict(csv_config.empty_df_kwargs or {})
                if csv_dict_key == "availability":
                    kwargs.update(
                        {
                            "tech_df": secondary_df,
                            "selected_types": selected_types,
                        }
                    )
                csv_config.empty_df_fn()

    def get_tech_mapping(self):
        """Load the technology mapping csv."""
        current_dir = os.path.dirname(__file__)
        file_path = os.path.join(current_dir, "..", "setting", "tech_mapping.csv")
        return pd.read_csv(file_path)

    def list_to_key(self, selected_types: list) -> str:
        """Hash a list of selected types.

        This is used to generate short, unique widget keys as an alternative to joining
        a potentially long list of selected types.

        Parameters
        ----------
        selected_types : list
            The list of technology types selected by the user.

        Returns
        -------
        str
            An eight character hash
        """
        joined = ",".join(sorted(selected_types))
        return hashlib.md5(joined.encode(), usedforsecurity=False).hexdigest()[:8]

    def filter_df_generic(
        self,
        df: pd.DataFrame,
        filter_col: str,
        selected_types: list,
    ) -> pd.DataFrame:
        """Filter the input df based on the user's selection.

        Parameters
        ----------
        df : pd.DataFrame
            The input dataframe to filter.
        filter_col : str
            Name of column containing the filter values.
        selected_type : str
            Values selected by user to filter by.

        Returns
        -------
        pd.DataFrame
            The filtered dataframe.
        """
        return df[df[filter_col].isin(selected_types)]

    def filter_df_decomission(
        self,
        df: pd.DataFrame,
        filter_col: str,
        selected_types: list,
    ) -> pd.DataFrame:
        """Filter the decommission capacity df based on the user's selection.

        Parameters
        ----------
        df : pd.DataFrame
            The input dataframe to filter.
        filter_col : str
            Name of column containing the filter value.
        selected_type : str
            Value selected by user to filter by.

        Returns
        -------
        pd.DataFrame
            The filtered dataframe.
        """
        return df[df[filter_col].str.split("_").str[-1].isin(selected_types)]

    def empty_df_message_generic(self):
        """Display a generic info message when the dataframe is empty."""
        info_message = (
            "No data required in this table for the selected technology "
            "type(s) and country(ies)."
        )
        st.info(info_message)

    def visualise_data(self, df: pd.DataFrame, csv_config: CsvDictConfig):
        """Visualise cost data with a line graph.

        Parameters
        ----------
        csv_config : CsvDictConfig
            The csv configuration containing the csv identifier and other settings.
        df : pd.DataFrame
            The dataframe already filtered by technology type.
        """
        tech_mapping = self.get_tech_mapping()

        colour_map = dict(
            zip(tech_mapping["original_names"], tech_mapping["hex_codes"])
        )
        name_map = dict(zip(tech_mapping["original_names"], tech_mapping["nice_names"]))

        countries = df["country"].unique()
        country_select_key = f"country_select_key{csv_config.identifier}"

        selected_country = st.pills(
            "Select a country",
            options=countries,
            default=countries[0],
            selection_mode="single",
            key=country_select_key,
        )

        filtered_df = df[df["country"] == selected_country]

        if "hourly" in csv_config.identifier:
            filtered_df = self.resample_to_monthly(filtered_df, csv_config.filter_col)
            start_date = dt.date(self.year_list[0], 1, 1)
            end_date = dt.date(self.year_list[0], 12, 31)
            selected_range = st.date_input(
                "Select month range",
                value=(start_date, end_date),
                min_value=start_date,
                max_value=end_date,
                key=f"{csv_config.identifier}_date_range_{selected_country}",
            )

            selected_tech = st.selectbox(
                "Select specific technology",
                options=filtered_df[csv_config.filter_col].unique(),
                index=0,
                key=(
                    f"{csv_config.identifier}_tech_select_{csv_config.filter_col}"
                    f"_{selected_country}"
                ),
            )

            start_month, end_month = selected_range[0].month, selected_range[1].month
            filtered_df = filtered_df[
                (filtered_df["month"] >= start_month)
                & (filtered_df["month"] <= end_month)
                & (filtered_df[csv_config.filter_col] == selected_tech)
            ]

            x = "month"
            leg_col = "node"
        else:
            x = "year"
            leg_col = csv_config.filter_col

        y_options = [
            col
            for col in filtered_df.columns
            if col not in {"country", "node", x}
            and (
                pd.api.types.is_float_dtype(filtered_df[col])
                or pd.api.types.is_integer_dtype(filtered_df[col])
            )
        ]

        y = st.selectbox(
            "Legend column",
            options=y_options,
            index=0,
            key=f"{csv_config.identifier}_legend_col_{selected_country}",
        )

        labels = [csv_config.uniform_unit if csv_config.uniform_unit else y]

        fig = px.line(
            filtered_df,
            x=x,
            y=y,
            labels={y: labels[0]},
            color=leg_col,
            color_discrete_map=colour_map,
        )

        # Update legends to use nice names
        fig.for_each_trace(lambda t: t.update(name=name_map.get(t.name, t.name)))
        fig.update_layout(
            height=300, legend_title_text=re.sub(r"_+", " ", leg_col).capitalize()
        )
        fig.update_yaxes(range=[0, 1.2 * filtered_df[y].max()])

        chart_key = f"plot_{csv_config.identifier}_{selected_country}_{x}_{y}_{leg_col}"
        st.plotly_chart(fig, width="stretch", key=chart_key)

    def resample_to_monthly(self, df: pd.DataFrame, leg_col: str) -> pd.DataFrame:
        """Resamples dataframe from hourly to monthly.

        Parameters
        ----------
        df : pd.DataFrame
            Input df to resample. Assumes country, node, and techology columns in
            addition to the hourly data.

        Returns
        -------
        pd.DataFrame
            The resampled dataframe.
        """
        # Melt the df - wide to long for hours
        value_cols = [c for c in df.columns if c.isdigit()]
        df_melted = df.melt(
            id_vars=["country", "node", leg_col],
            value_vars=value_cols,
            var_name="hour",
            value_name="value",
        )
        df_melted["hour"] = df_melted["hour"].astype(int)

        start = pd.Timestamp("2000-01-01 00:00:00")  # Assume an arbitrary year
        df_melted["datetime"] = start + pd.to_timedelta(df_melted["hour"], unit="h")

        # Resample to monthly
        df_monthly = (
            df_melted.set_index("datetime")
            .groupby(["country", "node", leg_col])
            .resample("ME")["value"]
            .mean()
            .reset_index()
        )

        df_monthly["month"] = df_monthly["datetime"].dt.month
        df_monthly = df_monthly.drop(columns=["datetime"])

        return df_monthly

    def get_fuel_mapping(self, tech_df: pd.DataFrame, selected_types: str) -> dict:
        """Get mapping of technology types to their carriers.

        Parameters
        ----------
        tech_df : pd.DataFrame
            Dataframe containing technology data.
        selected_types : str
            Selected technology hash string from list_to_key.

        Returns
        -------
        dict
            Dictionary of technology type to carrier.
        """
        type_to_carrier = (
            tech_df.loc[
                tech_df["technology"].isin(selected_types), ["technology", "carrier"]
            ]
            .drop_duplicates(subset="technology")
            .set_index("technology")["carrier"]
            .to_dict()
        )

        return type_to_carrier

    def update_csv_file(
        self,
        file_path: str,
        row_identifier: str,
        column_name: str,
        new_value: str,
    ):
        """Make targeted changes to CSV file without loading entire file into memory."""
        temp_file = NamedTemporaryFile(mode="w", delete=False, newline="")

        try:
            with (
                open(file_path, encoding="utf-8") as csvfile,
                temp_file,
            ):
                reader = csv.DictReader(csvfile)
                fieldnames = reader.fieldnames
                writer = csv.DictWriter(temp_file, fieldnames=fieldnames)
                writer.writeheader()

                for i, row in enumerate(reader):
                    if str(i) == row_identifier:  # Compare with row index
                        row[column_name] = new_value

                    # Convert NaNs (which pandas reads in empty cells as) to back to
                    # empty strings before writing back to the file
                    for key, val in row.items():
                        if val is None or val == "" or str(val).lower() == "nan":
                            row[key] = ""

                    writer.writerow(row)

            shutil.move(temp_file.name, file_path)
            return True
        except Exception as exc:
            if os.path.exists(temp_file.name):
                os.unlink(temp_file.name)
            raise exc
