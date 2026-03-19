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
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

pd.set_option("future.no_silent_downcasting", True)


class DataFrameWidgetsHandler:
    """Handle the DataFrame widget."""

    def __init__(self, input_config: dict):
        self.input_config = input_config
        self.base_config = st.session_state.base_config
        self.base_input_path = st.session_state.input_path
        self.year_list = self.base_config["base_configs"]["years"]

        sectors = self.base_config["base_configs"]["sector"]
        self.has_industry_sector = (
            "i" in sectors
            if isinstance(sectors, str)
            else any("i" in sector for sector in sectors)
        )
        self.has_transport_sector = (
            "t" in sectors
            if isinstance(sectors, str)
            else any("t" in sector for sector in sectors)
        )

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
            Key to session_state var for whether dataframe has changed.

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

    def set_up_df_with_charts(
        self,
        title: str,
        selected_types: list,
        input_df: pd.DataFrame | None = None,
        sector: str | None = None,
        selected_classes: list | None = None,
        selected_countries: list | None = None,
        selected_scenario: str | None = None,
    ):
        """Set up the widget with a df and a chart (if it's enabled).

        First tab contains the editable df and save button. Second tab contains the
        visualisation of the selected data if enabled.

        Parameters
        ----------
        sector : str
            Global_input, Power, Industry, or Transport.
        title : str
            Title to the csv in the csvs_dict.
        input_df : pd.DataFrame
            The input dataframe.
        selected_types : list
            Technology types selected by the user.
        selected_classes : list
            Class(es) selected by the user in the global class select widget.
        selected_countries : list, optional
            Country(s) selected by the user in the global country select widget.
        """
        selected_classes = selected_classes or []

        table_config, input_csv_path = self.find_table_config_and_path(
            title=title,
            sector=sector,
            selected_scenario=selected_scenario,
        )

        csv_identifier = table_config["identifier"]

        widget_scope = self.build_widget_scope_key(
            sector or "Global_input",
            title,
            selected_types,
            selected_classes,
            selected_countries,
        )
        edited_df_key = f"{title}_{csv_identifier}_editor_{widget_scope}"
        has_changes_key = f"has_changes_{title}_{csv_identifier}_{widget_scope}"
        save_button_key = f"save_{title}_{csv_identifier}_{widget_scope}"

        with st.expander(f"{title}"):
            st.write(f"### {title}")

            path_to_display = os.path.normpath(input_csv_path)

            st.markdown(
                f"<small><i>{path_to_display}</i></small>", unsafe_allow_html=True
            )

            if input_df is None:
                if not os.path.exists(input_csv_path):
                    st.error(f"File not found: {input_csv_path}")
                    return
                input_df = pd.read_csv(input_csv_path)

            if table_config["tag_name"] == "decommission":
                filtered_df = self.filter_df_decommission(
                    df=input_df,
                    filter_col=table_config["filter_col"],
                    selected_types=selected_types,
                )
            elif table_config["tag_name"] == "fuel_costs":
                fuels = (
                    self.get_fuel_mapping(selected_types)
                    if "Link" in selected_classes
                    else {}
                )
                filtered_df = self.filter_df_generic(
                    df=input_df,
                    filter_col=table_config["filter_col"],
                    selected_types=list(fuels.values()),
                )
            elif table_config["tag_name"] == "direct_air_capture":
                filtered_df = input_df
            else:
                filtered_df = self.filter_df_generic(
                    df=input_df,
                    filter_col=table_config["filter_col"],
                    selected_types=selected_types,
                )

            if selected_countries and "country" in filtered_df.columns:
                filtered_df = filtered_df[
                    filtered_df["country"].isin(selected_countries)
                ]

            edited_df = pd.DataFrame()
            if not filtered_df.empty:
                if table_config["with_charts"] and not table_config.get("timeseries"):
                    tab1, tab2 = st.tabs(["Table", "Visualisation"])

                    with tab1:
                        edited_df, to_save = self.create_editable_df(
                            filtered_df,
                            edited_df_key,
                            has_changes_key,
                        )

                    with tab2:
                        self.visualise_data(filtered_df, table_config, widget_scope)
                elif table_config.get("timeseries"):
                    self.visualise_data(filtered_df, table_config, widget_scope)
                    to_save = False
                else:
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
            else:
                self.empty_df_message_generic()

    def find_table_config_and_path(
        self,
        title: str,
        sector: str | None,
        selected_scenario: str | None = None,
    ) -> tuple[dict, str]:
        """Find table config and CSV path for a section."""
        if not sector or sector == "Global_input":
            table_config = self.input_config["Global_input"][title]
            input_csv_path = os.path.join(
                self.base_input_path,
                "global_input",
                table_config["csv_name"],
            )
            return table_config, input_csv_path

        table_config = self.input_config[sector][title]

        if sector in {"Power", "Industry", "Transport"}:
            scenario_name = (
                selected_scenario
                or self.base_config["path_configs"]["input_scenario_name"]
            )
            input_csv_path = os.path.join(
                self.base_input_path,
                scenario_name,
                sector.lower(),
                table_config["csv_name"],
            )
        else:
            input_csv_path = os.path.join(
                self.base_input_path,
                sector.lower(),
                table_config["csv_name"],
            )

        return table_config, input_csv_path

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

    def build_widget_scope_key(self, *parts) -> str:
        """Hash widget context into a short stable key."""
        flattened_parts = []
        for part in parts:
            if part is None:
                flattened_parts.append("")
            elif isinstance(part, list):
                flattened_parts.append("|".join(map(str, sorted(part))))
            else:
                flattened_parts.append(str(part))

        joined = "::".join(flattened_parts)
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

    def filter_df_decommission(
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

    def visualise_data(
        self,
        df: pd.DataFrame,
        table_config: dict,
        widget_scope: str,
    ):
        """Visualise cost data with a line graph.

        Parameters
        ----------
        table_config : dict
            The table configuration containing the table identifier and other settings.
        df : pd.DataFrame
            The dataframe already filtered by technology type.
        """
        tech_mapping = self.get_tech_mapping()

        colour_map = dict(
            zip(tech_mapping["original_names"], tech_mapping["hex_codes"])
        )
        name_map = dict(zip(tech_mapping["original_names"], tech_mapping["nice_names"]))
        identifier = table_config["identifier"]
        countries = df["country"].unique()
        country_select_key = f"country_select_key_{identifier}_{widget_scope}"

        selected_country = st.pills(
            "Select a country",
            options=countries,
            default=countries[0],
            selection_mode="single",
            key=country_select_key,
        )

        filtered_df = df[df["country"] == selected_country]

        if "hourly" in identifier:
            selected_tech = st.selectbox(
                "Select specific technology",
                options=filtered_df[table_config["filter_col"]].unique(),
                index=0,
                key=(
                    f"{identifier}_tech_select_{table_config['filter_col']}"
                    f"_{selected_country}_{widget_scope}"
                ),
            )

            averaging_period = st.segmented_control(
                "Averaging period",
                options=["hourly", "daily", "monthly"],
                default="daily",
                selection_mode="single",
                key=f"{identifier}_avg_period_{selected_country}_{widget_scope}",
            )
            averaging_period = averaging_period or "daily"

            start_date = dt.date(self.year_list[0], 1, 1)
            end_date = dt.date(self.year_list[0], 12, 31)
            selected_range = st.date_input(
                "Select date range",
                value=(start_date, end_date),
                min_value=start_date,
                max_value=end_date,
                key=f"{identifier}_date_range_{selected_country}",
            )

            if isinstance(selected_range, tuple) and len(selected_range) == 2:
                start_datetime = pd.Timestamp(selected_range[0])
                end_datetime = pd.Timestamp(selected_range[1]) + pd.Timedelta(days=1)
            else:
                st.error(
                    "Please select a valid date range with both start and end dates."
                )
                return

            filtered_df = filtered_df[
                filtered_df[table_config["filter_col"]] == selected_tech
            ]
            filtered_df = self.resample_data_time_resolution(
                filtered_df,
                table_config["filter_col"],
                averaging_period,
            )
            filtered_df = filtered_df[
                (filtered_df["datetime"] >= start_datetime)
                & (filtered_df["datetime"] < end_datetime)
            ]

            if filtered_df.empty:
                st.info("No data available for the selected technology/date range.")
                return

            if averaging_period == "hourly":
                year_start = pd.Timestamp(f"{self.year_list[0]}-01-01 00:00:00")
                filtered_df["hour"] = (
                    (filtered_df["datetime"] - year_start).dt.total_seconds() // 3600
                ).astype(int)
                x = "hour"
            elif averaging_period == "daily":
                filtered_df["day"] = filtered_df["datetime"].dt.dayofyear
                x = "day"
            else:
                filtered_df["month"] = filtered_df["datetime"].dt.month
                x = "month"

            leg_col = "node"
        else:
            x = "year"
            leg_col = table_config["filter_col"]

        y_options = [
            col
            for col in filtered_df.columns
            if col not in {"country", "node", x, "datetime"}
            and (
                pd.api.types.is_float_dtype(filtered_df[col])
                or pd.api.types.is_integer_dtype(filtered_df[col])
            )
        ]

        y = st.selectbox(
            "Legend column",
            options=y_options,
            index=0,
            key=f"{identifier}_legend_col_{selected_country}_{widget_scope}",
        )
        label = table_config.get("uniform_unit") or y
        labels = ({y: label},)

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

        chart_key = (
            f"plot_{identifier}_{selected_country}_{x}_{y}_{leg_col}_{widget_scope}"
        )
        st.plotly_chart(fig, width="stretch", key=chart_key)

    def resample_data_time_resolution(
        self,
        df: pd.DataFrame,
        leg_col: str,
        averaging_period: str,
    ) -> pd.DataFrame:
        """Resample dataframe to selected averaging period.

        Parameters
        ----------
        df : pd.DataFrame
            Input df to resample. Assumes country, node, and technology columns in
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

        start = pd.Timestamp(f"{self.year_list[0]}-01-01 00:00:00")
        df_melted["datetime"] = start + pd.to_timedelta(df_melted["hour"], unit="h")

        resample_rule = {
            "hourly": "h",
            "daily": "D",
            "monthly": "ME",
        }.get(averaging_period, "D")

        df_resampled = (
            df_melted.set_index("datetime")
            .groupby(["country", "node", leg_col])
            .resample(resample_rule)["value"]
            .mean()
            .reset_index()
        )

        return df_resampled

    def get_fuel_mapping(self, selected_types: list) -> dict:
        """Get mapping of technology types to their carriers.

        Parameters
        ----------
        selected_types : list
            Selected technology hash string from list_to_key.

        Returns
        -------
        dict
            Dictionary of technology type to carrier.
        """
        tech_path = os.path.join(
            self.base_input_path,
            "global_input",
            self.input_config["Global_input"]["Technologies"]["csv_name"],
        )
        tech_df = pd.read_csv(tech_path)
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
            with open(file_path, encoding="utf-8") as csvfile, temp_file:
                reader = csv.DictReader(csvfile)
                fieldnames = reader.fieldnames or []
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
