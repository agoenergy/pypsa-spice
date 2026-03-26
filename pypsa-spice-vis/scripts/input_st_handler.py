# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Helper functions for handling streamlit input UI and CSV editing."""

import csv
import datetime as dt
import os
import re
import shutil
import time
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from scripts.data_utils import load_tech_info_mapping_df

SECTOR_TITLES = {
    "Power": ":material/bolt: Power",
    "Industry": ":material/construction: Industry",
    "Transport": ":material/directions_car: Transport",
}


# =============================================================================
# Sidebar navigation + section headers
# =============================================================================


def generate_sector_title(selected_sector: str) -> str:
    """Return the formatted title for a sector page."""
    return SECTOR_TITLES.get(selected_sector, selected_sector)


def generate_global_markdown_message() -> None:
    """Generate the warning message for global input pages."""
    st.markdown(
        """
        <p style="color:orange; font-weight:bold;">
        ⚠️ Changes made to the global input files will be automatically applied
        across all scenarios.
        </p>
        """,
        unsafe_allow_html=True,
    )


# =============================================================================
# Paths, variables, and mapping helper functions
# =============================================================================


def get_all_countries() -> list[str]:
    """Filter all countries from the base_configs. Used for country pills rendering."""
    return list(st.session_state.base_config["base_configs"]["regions"].keys())


def get_class_type_for_timeseries_data(
    tech_df: pd.DataFrame,
    selected_sector: str,
    selected_types: list[str],
) -> list[str]:
    """Define PyPSA component type for timeseries technologies."""
    if "class" not in tech_df.columns:
        return selected_types

    class_mapping = dict(zip(tech_df["technology"], tech_df["class"]))
    if selected_sector == "Power":
        allowed_classes = {"Generator", "StorageUnit", "Store", "Link"}
    else:  # Other sectors
        allowed_classes = {"Link", "StorageUnit", "Store"}

    filtered_types = [
        technology
        for technology in selected_types
        if class_mapping.get(technology) in allowed_classes
    ]
    return filtered_types or selected_types


def get_fuel_mapping(selected_types: list[str], input_config: dict) -> dict:
    """Return a mapping from technology type to fuel carrier."""
    tech_path = os.path.join(
        st.session_state.input_path,
        "global_input",
        input_config["Global_input"]["Technologies"]["csv_name"],
    )
    tech_df = pd.read_csv(tech_path)
    return (
        tech_df.loc[
            tech_df["technology"].isin(selected_types), ["technology", "carrier"]
        ]
        .drop_duplicates(subset="technology")
        .set_index("technology")["carrier"]
        .to_dict()
    )


def get_table_config_and_path(
    title: str,
    sector: str | None,
    input_config: dict,
    selected_scenario: str | None = None,
) -> tuple[dict, str]:
    """Get the table configuration and CSV path for a given input scenario."""
    base_input_path = st.session_state.input_path
    base_config = st.session_state.base_config

    if not sector or sector == "Global_input":
        table_config = input_config["Global_input"][title]
        input_csv_path = os.path.join(
            base_input_path,
            "global_input",
            table_config["csv_name"],
        )
        return table_config, input_csv_path

    table_config = input_config[sector][title]
    if sector in {"Power", "Industry", "Transport"}:
        scenario_name = (
            selected_scenario or base_config["path_configs"]["input_scenario_name"]
        )
        input_csv_path = os.path.join(
            base_input_path,
            scenario_name,
            sector.lower(),
            table_config["csv_name"],
        )
        return table_config, input_csv_path

    input_csv_path = os.path.join(
        base_input_path,
        sector.lower(),
        table_config["csv_name"],
    )
    return table_config, input_csv_path


def get_tech_mapping() -> pd.DataFrame:
    """Load the tech_mapping CSV used for visualising data."""
    current_dir = os.path.dirname(__file__)
    file_path = os.path.join(current_dir, "..", "setting", "tech_mapping.csv")
    return pd.read_csv(file_path)


def get_empty_df_notice_message() -> None:
    """Display a generic message when the filtered dataframe is empty."""
    st.info(
        "No data required in this table for the selected technology type(s) and "
        "country(ies)."
    )


def get_unique_type_key(*parts: object) -> str:
    """Create a unique streamlit key by combining all inputs from parts."""
    combine_all_parts_str = (
        (
            ""
            if part is None
            else (
                "|".join(map(str, sorted(part)))
                if isinstance(part, list)
                else str(part)
            )
        )
        for part in parts
    )
    return "::".join(combine_all_parts_str)


def resample_data_time_resolution(
    df: pd.DataFrame,
    leg_col: str,
    averaging_period: str,
) -> pd.DataFrame:
    """Resample an hourly table to the selected time resolution."""
    value_cols = [column for column in df.columns if column.isdigit()]
    df_melted = df.melt(
        id_vars=["country", "node", leg_col],
        value_vars=value_cols,
        var_name="hour",
        value_name="value",
    )
    df_melted["hour"] = df_melted["hour"].astype(int)

    year = st.session_state.base_config["base_configs"]["years"][0]
    start = pd.Timestamp(f"{year}-01-01 00:00:00")
    df_melted["datetime"] = start + pd.to_timedelta(df_melted["hour"], unit="h")

    resample_rule = {
        "hourly": "h",
        "daily": "D",
        "monthly": "ME",
    }.get(averaging_period, "D")

    return (
        df_melted.set_index("datetime")
        .groupby(["country", "node", leg_col])
        .resample(resample_rule)["value"]
        .mean()
        .reset_index()
    )


# =============================================================================
# Filter widgets or tables
# =============================================================================


def set_available_technology_df(
    selected_sector: str, input_config: dict
) -> pd.DataFrame:
    """Load, filter and create the technology table under selected sector."""
    tech_info_df = load_tech_info_mapping_df()
    # Select default technology tables from the selected sector
    tech_info_df = tech_info_df[tech_info_df["sector"] == selected_sector]

    tech_path = os.path.join(
        st.session_state.input_path,
        "global_input",
        input_config["Global_input"]["Technologies"]["csv_name"],
    )
    # technology tables from global input in the project folder
    tech_df = pd.read_csv(tech_path)

    # List all available technology names from tech_info_df
    allowed_technologies = set()
    if not tech_info_df.index.empty:
        allowed_technologies.update(tech_info_df.index.astype(str).str.strip().tolist())

    for column in ["original_names", "nice_names"]:
        if column in tech_info_df.columns:
            allowed_technologies.update(
                tech_info_df[column].dropna().astype(str).str.strip().tolist()
            )

    # map the technology table to only include technologies relevant to the selected
    # sector based on the tech info mapping
    if "technology" in tech_df.columns and allowed_technologies:
        tech_df["technology"] = tech_df["technology"].astype(str).str.strip()
        tech_df = tech_df[tech_df["technology"].isin(allowed_technologies)]

    return tech_df


def set_general_filter_df(
    df: pd.DataFrame,
    filter_col: str,
    selected_types: list[str],
) -> pd.DataFrame:
    """Set up a dataframe by the selected technology or profile values."""
    return df[df[filter_col].isin(selected_types)]


def set_decommission_filter_df(
    df: pd.DataFrame,
    filter_col: str,
    selected_types: list[str],
) -> pd.DataFrame:
    """Set up the decommission capacity table based on the selected technologies."""
    return df[df[filter_col].str.split("_").str[-1].isin(selected_types)]


# =============================================================================
# Table editing and update helper functions
# =============================================================================


def convert_inf_to_strings_in_df(df: pd.DataFrame) -> pd.DataFrame:
    """Convert infinite numeric values to strings for Streamlit editing display."""
    editable_df = df.copy()

    # find any cell with inf values in the df and convert it into "inf" string
    for column in editable_df.select_dtypes(include=[np.number]).columns:
        inf_mask = np.isinf(editable_df[column])
        if inf_mask.any():
            editable_df[column] = editable_df[column].astype(object)
            editable_df.loc[inf_mask, column] = "inf"

    return editable_df


def convert_strings_to_inf_in_df(
    edited_df: pd.DataFrame,
    reference_df: pd.DataFrame,
) -> pd.DataFrame:
    """Convert inf strings back to numeric inf after completion of Streamlit editing."""
    result_df = edited_df.copy()

    # find any cell with inf values in the df and convert it into "inf" string
    for column in reference_df.select_dtypes(include=[float]).columns:
        inf_mask = result_df[column] == "inf"
        if inf_mask.any():
            result_df.loc[inf_mask, column] = np.inf

    return result_df


def create_editable_df(
    filtered_df: pd.DataFrame,
    edited_df_key: str,
    has_changes_key: str,
) -> tuple[pd.DataFrame, bool]:
    """Render an editable dataframe and validate numeric columns."""
    to_save = True
    editable_df = convert_inf_to_strings_in_df(filtered_df)

    editable_cols = filtered_df.select_dtypes(
        include=["number", float, int, "bool"]
    ).columns
    disabled_cols = [
        column
        for column in filtered_df.columns
        if column not in editable_cols and column != "max_supply [MWh/year]"
    ]

    edited_df = st.data_editor(
        editable_df,
        hide_index=True,
        key=edited_df_key,
        disabled=disabled_cols,
        on_change=lambda: st.session_state.update({has_changes_key: True}),
    )
    result_df = convert_strings_to_inf_in_df(edited_df, filtered_df)

    for column in filtered_df.select_dtypes(include=[float]).columns:
        try:
            result_df[column] = result_df[column].astype(float)
        except Exception:
            invalid_mask = result_df[column].apply(
                lambda value: not (isinstance(value, (float, int)) or value == np.inf)
            )
            if invalid_mask.any():
                st.error(
                    f"Column '{column}' contains invalid entries. Only numbers or "
                    "'inf' allowed."
                )
                to_save = False

            result_df[column] = result_df[column].astype(float, errors="ignore")

    return result_df, to_save


# =============================================================================
# Layout renderers (single/dual chart layout)
# =============================================================================


def visualise_data(df: pd.DataFrame, table_config: dict, widget_scope: str) -> None:
    """Render the existing line-chart visualisation for a filtered input table."""
    tech_mapping = get_tech_mapping()
    colour_map = dict(zip(tech_mapping["original_names"], tech_mapping["hex_codes"]))
    name_map = dict(zip(tech_mapping["original_names"], tech_mapping["nice_names"]))
    identifier = table_config["identifier"]

    countries = df["country"].unique()
    selected_country = st.pills(
        "Select a country",
        options=countries,
        default=countries[0],
        selection_mode="single",
        key=f"country_select_key_{identifier}_{widget_scope}",
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
        averaging_period = (
            st.segmented_control(
                "Averaging period",
                options=["hourly", "daily", "monthly"],
                default="daily",
                selection_mode="single",
                key=f"{identifier}_avg_period_{selected_country}_{widget_scope}",
            )
            or "daily"
        )

        year = st.session_state.base_config["base_configs"]["years"][0]
        start_date = dt.date(year, 1, 1)
        end_date = dt.date(year, 12, 31)
        selected_range = st.date_input(
            "Select date range",
            value=(start_date, end_date),
            min_value=start_date,
            max_value=end_date,
            key=f"{identifier}_date_range_{selected_country}_{widget_scope}",
        )

        if not (isinstance(selected_range, tuple) and len(selected_range) == 2):
            st.error("Please select a valid date range with both start and end dates.")
            return

        start_datetime = pd.Timestamp(selected_range[0])
        end_datetime = pd.Timestamp(selected_range[1]) + pd.Timedelta(days=1)
        filtered_df = filtered_df[
            filtered_df[table_config["filter_col"]] == selected_tech
        ]
        filtered_df = resample_data_time_resolution(
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
            year_start = pd.Timestamp(f"{year}-01-01 00:00:00")
            filtered_df["hour"] = (
                (filtered_df["datetime"] - year_start).dt.total_seconds() // 3600
            ).astype(int)
            x_axis = "hour"
        elif averaging_period == "daily":
            filtered_df["day"] = filtered_df["datetime"].dt.dayofyear
            x_axis = "day"
        else:
            filtered_df["month"] = filtered_df["datetime"].dt.month
            x_axis = "month"

        legend_column = "node"
    else:
        x_axis = "year"
        legend_column = table_config["filter_col"]

    y_options = [
        column
        for column in filtered_df.columns
        if column not in {"country", "node", x_axis, "datetime"}
        and (
            pd.api.types.is_float_dtype(filtered_df[column])
            or pd.api.types.is_integer_dtype(filtered_df[column])
        )
    ]

    y_axis = st.selectbox(
        "Legend column",
        options=y_options,
        index=0,
        key=f"{identifier}_legend_col_{selected_country}_{widget_scope}",
    )
    label = table_config.get("uniform_unit") or y_axis

    fig = px.line(
        filtered_df,
        x=x_axis,
        y=y_axis,
        labels={y_axis: label},
        color=legend_column,
        color_discrete_map=colour_map,
    )
    fig.for_each_trace(
        lambda trace: trace.update(name=name_map.get(trace.name, trace.name))
    )
    fig.update_layout(
        height=300,
        legend_title_text=re.sub(r"_+", " ", legend_column).capitalize(),
    )
    fig.update_yaxes(range=[0, 1.2 * filtered_df[y_axis].max()])

    chart_key = (
        f"plot_{identifier}_{selected_country}_{x_axis}_{y_axis}_"
        f"{legend_column}_{widget_scope}"
    )
    st.plotly_chart(fig, width="stretch", key=chart_key)


def render_input_table_section(
    title: str,
    selected_types: list[str],
    input_config: dict,
    input_df: pd.DataFrame | None = None,
    sector: str | None = None,
    selected_classes: list[str] | None = None,
    selected_countries: list[str] | None = None,
    selected_scenario: str | None = None,
) -> None:
    """Render one input table section with editing and charts."""
    selected_classes = selected_classes or []
    table_config, input_csv_path = get_table_config_and_path(
        title=title,
        sector=sector,
        input_config=input_config,
        selected_scenario=selected_scenario,
    )
    csv_identifier = table_config["identifier"]
    unique_type_key = get_unique_type_key(
        sector or "Global_input",
        title,
        selected_types,
        selected_classes,
        selected_countries,
    )

    edited_df_key = f"{title}_{csv_identifier}_editor_{unique_type_key}"
    has_changes_key = f"has_changes_{title}_{csv_identifier}_{unique_type_key}"
    save_button_key = f"save_{title}_{csv_identifier}_{unique_type_key}"

    with st.expander(title):
        st.write(f"### {title}")
        st.markdown(
            f"<small><i>{os.path.normpath(input_csv_path)}</i></small>",
            unsafe_allow_html=True,
        )

        if input_df is None:
            if not os.path.exists(input_csv_path):
                st.error(f"File not found: {input_csv_path}")
                return
            input_df = pd.read_csv(input_csv_path)

        if table_config["tag_name"] == "decommission":
            filtered_df = set_decommission_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=selected_types,
            )
        elif table_config["tag_name"] == "fuel_costs":
            fuels = (
                get_fuel_mapping(selected_types, input_config)
                if "Link" in selected_classes
                else {}
            )
            filtered_df = set_general_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=list(fuels.values()),
            )
        elif table_config["tag_name"] == "direct_air_capture":
            filtered_df = input_df
        else:
            filtered_df = set_general_filter_df(
                df=input_df,
                filter_col=table_config["filter_col"],
                selected_types=selected_types,
            )

        if selected_countries and "country" in filtered_df.columns:
            filtered_df = filtered_df[filtered_df["country"].isin(selected_countries)]

        edited_df = pd.DataFrame()
        if filtered_df.empty:
            get_empty_df_notice_message()
            return

        if table_config["with_charts"] and not table_config.get("timeseries"):
            table_tab, visualisation_tab = st.tabs(["Table", "Visualisation"])
            with table_tab:
                edited_df, to_save = create_editable_df(
                    filtered_df,
                    edited_df_key,
                    has_changes_key,
                )
            with visualisation_tab:
                visualise_data(filtered_df, table_config, unique_type_key)
        elif table_config.get("timeseries"):
            visualise_data(filtered_df, table_config, unique_type_key)
            to_save = False
        else:
            edited_df, to_save = create_editable_df(
                filtered_df,
                edited_df_key,
                has_changes_key,
            )

        has_changes = st.session_state.get(has_changes_key, False)
        if to_save:
            render_save_button(
                filtered_df,
                edited_df,
                has_changes,
                has_changes_key,
                save_button_key,
                input_csv_path,
                message_delay=1,
            )


def render_demand_profiles_widget(
    selected_sector: str,
    sector_selected_countries: list[str] | None,
    input_config: dict,
) -> None:
    """Render the global demand profile widget for the active sector."""
    load_mapping = {
        "HV_LOAD": "Transmission/Wholesale market load (High voltage level)",
        "LV_LOAD": "Distribution/Building load (low/medium voltage level)",
        "IND_LOAD": "Industrial load (both high- and low-temperature heat)",
        "HPV_LOAD": "Transport load (high voltage level)",
        "LPV_LOAD": "Transport load (low/medium voltage level)",
    }
    reverse_load_mapping = {label: key for key, label in load_mapping.items()}
    defaults_by_sector = {
        "Power": ["HV_LOAD", "LV_LOAD"],
        "Industry": ["IND_LOAD"],
        "Transport": ["HPV_LOAD", "LPV_LOAD"],
    }

    default_profile_keys = [
        profile
        for profile in defaults_by_sector.get(selected_sector, [])
        if profile in load_mapping
    ]
    default_profile_selection = [
        load_mapping[profile] for profile in default_profile_keys
    ]

    selected_profile_label = st.selectbox(
        "Select type of demand/load profiles:",
        options=default_profile_selection,
        index=0,
        key=f"timeseries_demand_type_{selected_sector.lower()}",
        help="Choose the load profile variant to inspect for the active sector.",
    )
    st.markdown("Class: **Load**")

    selected_profile_key = reverse_load_mapping.get(selected_profile_label)

    if selected_profile_key is None:
        st.info("Select a demand profile type to load data.")
        return

    render_input_table_section(
        sector="Global_input",
        title="Demand_Profiles",
        selected_types=[selected_profile_key],
        selected_countries=sector_selected_countries,
        input_config=input_config,
    )


# =============================================================================
# Download helper functions
# =============================================================================


def update_csv_file(
    file_path: str,
    row_identifier: str,
    column_name: str,
    new_value: str,
) -> bool:
    """Update changes from Streamlit input editor back to the original CSV file."""
    temp_file = NamedTemporaryFile(mode="w", delete=False, newline="")

    try:
        with open(file_path, encoding="utf-8") as csv_file, temp_file:
            reader = csv.DictReader(csv_file)
            fieldnames = reader.fieldnames or []
            writer = csv.DictWriter(temp_file, fieldnames=fieldnames)
            writer.writeheader()

            for index, row in enumerate(reader):
                if str(index) == row_identifier:
                    row[column_name] = new_value

                for key, value in row.items():
                    if value is None or value == "" or str(value).lower() == "nan":
                        row[key] = ""

                writer.writerow(row)

        shutil.move(temp_file.name, file_path)
        return True
    except Exception:
        if os.path.exists(temp_file.name):
            os.unlink(temp_file.name)
        raise


def render_save_button(
    filtered_df: pd.DataFrame,
    edited_df: pd.DataFrame,
    has_changes: bool,
    has_changes_key: str,
    save_button_key: str,
    output_file_path: str,
    message_delay: float = 1,
) -> None:
    """Render the save button and save changes in the input data."""
    if st.button(
        "Save Changes",
        key=save_button_key,
        type="primary" if has_changes else "secondary",
        disabled=not has_changes,
    ):
        success = True
        for index in range(len(edited_df)):
            current_index = filtered_df.index[index]
            for column in filtered_df.columns:
                if filtered_df[column].iloc[index] != edited_df[column].iloc[index]:
                    success &= update_csv_file(
                        file_path=output_file_path,
                        row_identifier=str(current_index),
                        column_name=column,
                        new_value=str(edited_df[column].iloc[index]),
                    )

        if success:
            st.success("Changes saved successfully!")
            st.session_state[has_changes_key] = False
            time.sleep(message_delay)
            st.rerun()
        else:
            st.error("Error saving some changes")
