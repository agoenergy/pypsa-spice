# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

import os
import streamlit as st
import pandas as pd
import plotly.express as px
import time
from dataclasses import dataclass
from typing import Callable, Optional
import numpy as np
import re
import hashlib
from scripts.getters import Getters

pd.set_option("future.no_silent_downcasting", True)


class dfWidgetsHandler:
    def __init__(self):
        self.getters = Getters()
        self.input_ui_handler = InputUiHandler()
        self.csvs_dict = self.input_ui_handler.csvs_dict
        self.input_folder_path = st.session_state["input_folder_path"]
        project_folders = self.getters.get_project_folder_list(self.input_folder_path)
        if not project_folders:
            st.error(f"No valid project folders found in {self.input_folder_path}")

        self.input_data_project = project_folders[0]

        self.base_input_path = os.path.join(
            self.input_folder_path,
            self.input_data_project,
            "input",
        )

        self.global_input_path = os.path.join(
            self.base_input_path,
            "global_input",
        )
        sub_folder = (
            st.session_state["scenario"]
            if "scenario" in st.session_state
            else self.getters.init_config["path_configs"]["input_scenario_name"]
        )

        self.scenario_input_path = os.path.join(self.base_input_path, sub_folder)

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
        # Update csvs_dict with paths to all csvs
        for key, parts in self.csv_files_path_parts.items():
            if key in self.csvs_dict:
                base_folder = (
                    self.scenario_input_path
                    if key in self.scenario_input_file_keys
                    else self.global_input_path
                )
                self.csvs_dict[key].path = os.path.join(base_folder, *parts)

        # Check if files exist
        for key, csv_config in self.csvs_dict.items():
            if not os.path.exists(csv_config.path):
                st.error(f"{key} file not found at {csv_config.path}")
                return {}

        # Load all dataframes
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
        """Reload the scenario-related dataframes after the user has selected a scenario.

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
                self.input_folder_path,
                self.input_data_project,
                "input",
                selected_scenario,
                *parts,
            )

        # Reload the dataframes
        for key in self.scenario_input_file_keys:
            dfs[f"{key}_df"] = pd.read_csv(self.csvs_dict[key].path)

        return dfs


@dataclass
class CsvDictConfig:
    """Configuration for different csvs"""

    identifier: str
    filter_col: str
    title: str
    filter_fn: Callable
    empty_df_fn: Optional[Callable] = None
    empty_df_kwargs: Optional[dict] = None
    path: str = ""  # Populate after user selections are made


class InputUiHandler:
    def __init__(self):
        self.csvs_dict = {
            "technologies": CsvDictConfig(
                identifier="tech",
                filter_col="technology",
                title="Techonology Parameters",
                filter_fn=self._filter_df_generic,
            ),
            "links": CsvDictConfig(
                identifier="links",
                filter_col="type",
                title="Asset Parameters",
                filter_fn=self._filter_df_generic,
            ),
            "availability": CsvDictConfig(
                identifier="avail",
                filter_col="technology",
                title="Availability",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
                empty_df_kwargs={"msg": "Taking availability from availability.csv"},
            ),
            "demand": CsvDictConfig(
                identifier="demand",
                filter_col="profile_type",
                title="Demand Profiles",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
                empty_df_kwargs={"msg": "Taking demand from demand_profile.csv"},
            ),
            "pp_costs": CsvDictConfig(
                identifier="costs",
                filter_col="powerplant_type",
                title="Power Plant Costs",
                filter_fn=self._filter_df_generic,
            ),
            "potentials": CsvDictConfig(
                identifier="potentials",
                filter_col="type",
                title="Renewable Technical Potentials",
                filter_fn=self._filter_df_generic,
            ),
            "storage_cost": CsvDictConfig(
                identifier="storage_cost",
                filter_col="storage_type",
                title="Storage Costs",
                filter_fn=self._filter_df_generic,
            ),
            "storage_inflows": CsvDictConfig(
                identifier="storage_inflows",
                filter_col="technology",
                title="Storage Inflows",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
                empty_df_kwargs={"msg": "Taking inflows from storage_inflows.csv"},
            ),
            "decomission": CsvDictConfig(
                identifier="decomission",
                filter_col="name",
                title="Decomissioning",
                filter_fn=self._filter_df_decomission,
                empty_df_fn=self._empty_df_message_generic,
                empty_df_kwargs={"msg": "No decomissioning"},
            ),
            "fuel_costs": CsvDictConfig(
                identifier="fuel",
                filter_col="carrier",
                title="Fuel Costs",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
            ),
            "interconnector": CsvDictConfig(
                identifier="interconnector",
                filter_col="type",
                title="Interconnector",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
            ),
            "load": CsvDictConfig(
                identifier="load",
                filter_col="profile_type",
                title="Load",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
            ),
            "generator": CsvDictConfig(
                identifier="generator",
                filter_col="type",
                title="Assest Parameters - Generators",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
            ),
            "links": CsvDictConfig(
                identifier="links",
                filter_col="carrier",
                title="Assest Parameters - Links",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
            ),
            "storageunit": CsvDictConfig(
                identifier="storageunit",
                filter_col="carrier",
                title="Assest Parameters - StorageUnits",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
            ),
            "store": CsvDictConfig(
                identifier="store",
                filter_col="carrier",
                title="Assest Parameters - Stores",
                filter_fn=self._filter_df_generic,
                empty_df_fn=self._empty_df_message_generic,
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
            except:
                invalid_mask = result_df[col].apply(
                    lambda x: not (
                        isinstance(x, float) or isinstance(x, int) or x == np.inf
                    )
                )
                if invalid_mask.any():
                    st.error(
                        f"Column '{col}' contains invalid entries. Only numbers or 'inf' allowed."
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
        edited_df:
            The edited dataframe with changes made by the user.
        has_changes:
            True if df has been edited, False otherwise.
        has_changes_key : str
            Key to session_state var for whether dataframe has changed.
        save_button_key: str
            Button key.
        output_file_path : str
            Path to save the csv to.
        message_delay:
            How long (s)to display the 'success message' after saving for.
        """
        if st.button(
            "Save Changes",
            key=save_button_key,
            type="primary" if has_changes else "secondary",
            disabled=not has_changes,
        ):

            try:
                success = True
                # Iterate through each row in the edited dataframe
                for idx in range(len(edited_df)):
                    # Get the original dataframe index
                    current_index = filtered_df.index[idx]
                    for col in filtered_df.columns:
                        if filtered_df[col].iloc[idx] != edited_df[col].iloc[idx]:
                            success &= self.update_csv_file(
                                file_path=output_file_path,
                                row_identifier=str(current_index),
                                column_name=col,
                                new_value=str(edited_df[col].iloc[idx]),
                                # Use "Index" as indentifier column
                                identifier_column="Index",
                            )
                if success:
                    st.success("Changes saved successfully!")
                    st.session_state[has_changes_key] = False
                    time.sleep(message_delay)
                    # Force a complete reset cycle - resets widget internal states and
                    # makes sure our save button behaves normally (wait a bit first to
                    # allow the success message to stay on screen)
                    st.rerun()
                else:
                    st.error("Error saving some changes")
            except Exception as e:
                st.error(f"Error saving changes: {e}")

    def set_up_single_tab_widget(
        self,
        csv_dict_key: str,
        input_df: pd.DataFrame,
        selected_types: list,
        input_csv_path: str,
        selected_countries: list = None,
        selected_classes: Optional[list] = None,
        secondary_df: Optional[pd.DataFrame] = None,
    ):
        """Set up the widget with a single tab containing the editable df and save
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
        list_key = self._list_to_key(selected_types)
        edited_df_key = f"{csv_identifier}_editor_{list_key}"
        has_changes_key = f"has_changes_{csv_identifier}_{list_key}"
        save_button_key = f"save_{csv_identifier}_{list_key}"

        filter_col = csv_config.filter_col
        filtered_df = pd.DataFrame()
        no_data_msg = None

        # Filter the input df based on user's selection
        if csv_dict_key == "fuel_costs":
            if "Link" in selected_classes and secondary_df is not None:
                # Set up filtered_df and no_data_msg for the fuel_costs df
                try:
                    fuels = self._get_fuel_mapping(secondary_df, selected_types)
                    filtered_df = csv_config.filter_fn(
                        input_df, filter_col, list(fuels.values())
                    )
                except IndexError:
                    no_data_msg = f"No carrier found for {selected_types} in links data"
            else:
                no_data_msg = (
                    f"No fuel costs for technologies: {', '.join(selected_types)}"
                )

        else:
            # Apply the default filtering for all other csvs
            filtered_df = csv_config.filter_fn(input_df, filter_col, selected_types)

        # Apply country filter if available
        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        with st.container(border=True):
            st.write(f"### {csv_config.title}")

            path_to_display = os.path.normpath(input_csv_path)
            # link = f"file://{path_to_display}"
            # st.markdown(f'<small><i><a href="{link}">{path_to_display}</a></i></small>', unsafe_allow_html=True)
            st.markdown(
                f"<small><i>{path_to_display}</i></small>", unsafe_allow_html=True
            )

            if not filtered_df.empty:
                edited_df, to_save = self.create_editable_df(
                    filtered_df, edited_df_key, has_changes_key
                )

                # Check for changes to the dataframe
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
        secondary_df: Optional[pd.DataFrame] = None,
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

        # Create session state var keys based on the current csv
        list_key = self._list_to_key(selected_types)
        edited_df_key = f"{csv_identifier}_editor_{list_key}"
        has_changes_key = f"has_changes_{csv_identifier}_{list_key}"
        save_button_key = f"save_{csv_identifier}_{list_key}"

        filter_col = csv_config.filter_col
        filtered_df = csv_config.filter_fn(
            input_df,
            filter_col,
            selected_types,
        )

        # Apply country filter if available
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
                    self._visualise_data(filtered_df, csv_identifier)

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
                csv_config.empty_df_fn(**kwargs)

    def _get_tech_mapping(self):
        current_dir = os.path.dirname(__file__)
        file_path = os.path.join(current_dir, "..", "setting", "tech_mapping.csv")
        return pd.read_csv(file_path)

    def _list_to_key(self, selected_types: list) -> str:
        """Helper function to hash a list of selected types.

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
        return hashlib.md5(joined.encode()).hexdigest()[:8]

    def _filter_df_generic(
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

    def _filter_df_decomission(
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

    def _empty_df_message_generic(self, **kwargs):
        info_message = kwargs.get("msg")
        st.info(info_message)

    def _visualise_data(self, df: pd.DataFrame, csv_identifier: str):
        """Visualise cost data with a line graph.

        Parameters
        ----------
        csv_identifier : str
            The csv identifier from the csv_config dict (avail or costs).
        df : pd.DataFrame
            The dataframe already filtered by technology type.
        """
        tech_mapping = self._get_tech_mapping()

        colour_map = dict(
            zip(tech_mapping["original_names"], tech_mapping["hex_codes"])
        )
        name_map = dict(zip(tech_mapping["original_names"], tech_mapping["nice_names"]))

        countries = df["country"].unique()
        country_select_key = f"country_select_key{csv_identifier}"

        # Country filter within the vis tab
        selected_country = st.pills(
            "Select a country",
            options=countries,
            default=countries[0],  # Use the first country in the list as default
            selection_mode="single",
            key=country_select_key,
        )

        if csv_identifier == "avail":
            filtered_df = df[df["country"] == selected_country]
            resampled = self._resample_to_monthly(df, "technology")

            node_avail_toggle = st.toggle(
                "Include node filter in the Availability Profiles", value=False
            )

            if node_avail_toggle:
                nodes = df.loc[df["country"] == selected_country, "node"].unique()
                # # Node filter within the vis tab, only for Availability
                selected_avail_node = st.pills(
                    "Select a node",
                    options=nodes,
                    default=nodes[0],  # First node in the list as default
                    selection_mode="single",
                    key="node_select_key_avail",
                )
                filtered_df = resampled[
                    (resampled["country"] == selected_country)
                    & (resampled["node"] == selected_avail_node)
                ]
                leg_col = "technology"
            else:
                filtered_df = resampled[(resampled["country"] == selected_country)]
                leg_col = "node"

            x, y = "month", "value"
            labels = {"month": "Month", "value": "Value"}

        elif csv_identifier == "demand":
            filtered_df = df[df["country"] == selected_country]
            resampled = self._resample_to_monthly(df, "profile_type")

            node_demand_toggle = st.toggle(
                "Include node filter in the Demand Profile", value=False
            )

            if node_demand_toggle:
                nodes = df.loc[df["country"] == selected_country, "node"].unique()
                # # Node filter within the vis tab
                selected_demand_node = st.pills(
                    "Select a node",
                    options=nodes,
                    default=nodes[0],  # First node in the list as default
                    selection_mode="single",
                    key="node_select_key_avail",
                )
                filtered_df = resampled[
                    (resampled["country"] == selected_country)
                    & (resampled["node"] == selected_demand_node)
                ]
                leg_col = "profile_type"
            else:
                filtered_df = resampled[(resampled["country"] == selected_country)]
                leg_col = "node"

            x, y = "month", "value"
            labels = {"month": "Month", "value": "Value"}

        elif csv_identifier == "load":
            filtered_df = df[df["country"] == selected_country]

            node_load_toggle = st.toggle("Include node filter in the load", value=False)

            if node_load_toggle:
                nodes = df.loc[df["country"] == selected_country, "node"].unique()
                # # Node filter within the vis tab
                selected_load_node = st.pills(
                    "Select a node",
                    options=nodes,
                    default=nodes[0],  # First node in the list as default
                    selection_mode="single",
                    key="node_select_key_avail",
                )
                filtered_df = filtered_df[
                    (filtered_df["country"] == selected_country)
                    & (filtered_df["node"] == selected_load_node)
                ]
                leg_col = "profile_type"
            else:
                filtered_df = filtered_df[(filtered_df["country"] == selected_country)]
                leg_col = "node"

            x, y = "year", "total_load__mwh"
            labels = {"total_load__mwh": "Total Load (MW)", "year": "Year"}

        elif csv_identifier == "costs":
            filtered_df = df[df["country"] == selected_country]

            x, y = "year", "cap__usd_mw"
            leg_col = "powerplant_type"
            labels = {"cap__usd_mw": "Capital Cost (USD/MW)", "year": "Year"}
        else:
            raise ValueError(f"Unknown csv_identifier: {csv_identifier}")

        fig = px.line(
            filtered_df,
            x=x,
            y=y,
            labels=labels,
            color=leg_col,
            color_discrete_map=colour_map,
        )

        # Update legends to use nice names
        fig.for_each_trace(lambda t: t.update(name=name_map.get(t.name, t.name)))
        fig.update_layout(
            height=300, legend_title_text=re.sub(r"_+", " ", leg_col).capitalize()
        )
        fig.update_yaxes(range=[0, 1.2 * filtered_df[y].max()])

        st.plotly_chart(fig, use_container_width=True)

    def _resample_to_monthly(self, df: pd.DataFrame, leg_col: str) -> pd.DataFrame:
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

    def _get_fuel_mapping(self, tech_df: pd.DataFrame, selected_types):
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
        identifier_column: str,
    ):
        """Make targeted changes to CSV file without loading entire file into memory."""
        import csv
        from tempfile import NamedTemporaryFile
        import shutil

        temp_file = NamedTemporaryFile(mode="w", delete=False, newline="")

        try:
            with open(file_path, "r") as csvfile, temp_file:
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
        except Exception as e:
            if os.path.exists(temp_file.name):
                os.unlink(temp_file.name)
            raise e
