# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Helper functions for handling streamlit input UI and CSV editing."""

# pylint: disable=too-many-arguments,too-many-locals, too-many-positional-arguments
import csv
import datetime as dt
import os
import re
import shutil
import time
from datetime import date
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarfloat import ScalarFloat
from ruamel.yaml.scalarstring import DoubleQuotedScalarString

from scripts.data_utils import load_tech_info_mapping_df

pd.set_option("future.no_silent_downcasting", True)

SECTOR_TITLES = {
    "Power": ":material/bolt: Power",
    "Industry": ":material/construction: Industry",
    "Transport": ":material/directions_car: Transport",
}

SCENARIO_CONFIG_HEADER = """
# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

# coding: utf-8

# Scenario configuration includes model constraints, scenario and solver settings.

"""


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


@st.cache_data
def get_default_scenario_config() -> dict:
    """Load the default scenario configuration template used by the editor."""
    file_path = (
        Path(__file__).resolve().parents[2]
        / "data"
        / "scenario_config_template"
        / "scenario_config.default.yaml"
    )

    yaml_loader = YAML(typ="safe", pure=True)

    with file_path.open(encoding="utf-8") as file_handle:
        return yaml_loader.load(file_handle) or {}


def get_table_config_and_path(
    title: str,
    sector: str | None,
    input_config: dict,
    selected_scenario: str | None = None,
) -> tuple[dict, str]:
    """Get the table configuration and CSV path for a given input scenario."""
    base_input_path = st.session_state.input_path

    if not sector or sector == "Global_input":
        table_config = input_config["Global_input"][title]
        path_parts = ("global_input", table_config["csv_name"])
    else:
        table_config = input_config[sector][title]

        if sector in SECTOR_TITLES:
            if not selected_scenario:
                raise ValueError(f"selected_scenario is required for sector '{sector}'")
            path_parts = (
                selected_scenario,
                sector.lower(),
                table_config["csv_name"],
            )
        else:
            path_parts = (
                sector.lower(),
                table_config["csv_name"],
            )

    return table_config, os.path.join(base_input_path, *path_parts)


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
    return df[df[filter_col].str.contains("|".join(selected_types))]


def set_decommission_filter_df(
    df: pd.DataFrame,
    filter_col: str,
    selected_types: list[str],
) -> pd.DataFrame:
    """Set up the decommission capacity table based on the selected technologies."""
    return df[df[filter_col].str.split("_").str[-1].isin(selected_types)]


def set_filtered_timeseries_df(
    df: pd.DataFrame, table_config: dict, widget_scope: str
) -> tuple[pd.DataFrame, str]:
    """Set up a timeseries dataframe by the selected technology and date range."""
    identifier = table_config["identifier"]
    countries = df["country"].unique()
    filtered_df = df.copy()

    # Country filter
    selected_country = st.pills(
        "Select a country",
        options=countries,
        default=countries[0],
        selection_mode="single",
        key=f"country_select_key_{identifier}_{widget_scope}",
    )
    if selected_country is None:
        selected_country = countries[0]
    else:
        filtered_df = df[df["country"] == selected_country]

    return filtered_df, selected_country


# =============================================================================
# Table editing and update helper functions
# =============================================================================


def create_editable_df(
    filtered_df: pd.DataFrame,
    edited_df_key: str,
    has_changes_key: str,
) -> tuple[pd.DataFrame, bool]:
    """Render an editable dataframe and validate numeric columns."""
    to_save = True

    # convert inf values to string for streamlit data editor display
    editable_df = filtered_df.replace({np.inf: "inf"})

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

    # convert "inf" strings back to numeric inf
    result_df = edited_df.replace({"inf": np.inf})

    # validate numeric columns
    for column in filtered_df.select_dtypes(include=[float]).columns:
        try:
            result_df[column] = result_df[column].astype(float)
        except (ValueError, TypeError):
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
# Layout renderers (chart, table, config layout)
# =============================================================================


def render_line_chart(df: pd.DataFrame, table_config: dict, widget_scope: str) -> None:
    """Render a line chart for timeseries data or other yearly data."""
    tech_mapping = get_tech_mapping()
    colour_map = dict(zip(tech_mapping["original_names"], tech_mapping["hex_codes"]))
    name_map = dict(zip(tech_mapping["original_names"], tech_mapping["nice_names"]))
    identifier = table_config["identifier"]

    filtered_df, selected_country = set_filtered_timeseries_df(
        df, table_config, widget_scope
    )

    if "hourly" in identifier:
        # Hourly data
        # Technology filter
        selected_tech = st.selectbox(
            "Select specific technology",
            options=filtered_df[table_config["filter_col"]].unique(),
            index=0,
            key=(
                f"{identifier}_tech_select_{table_config['filter_col']}"
                f"_{selected_country}_{widget_scope}"
            ),
        )
        # Temporal resolution filter
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

        # Timeseries data only applies with the base year
        year = st.session_state.base_config["base_configs"]["years"][0]
        start_date = dt.date(year, 1, 1)
        end_date = dt.date(year, 12, 31)

        # Date range filter
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

        # Resample data to the selected time resolution for the selected date range
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

        # Set x and legend columns based on the selected averaging period
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
        # yearly data
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

    # visualisation
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


def render_demand_profiles_selectbox(selected_sector: str) -> list[str]:
    """Render the global demand profile selectbox."""
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

    # Filter the default profiles based on the selected sector and available load types
    default_profile_keys = [
        profile
        for profile in defaults_by_sector.get(selected_sector, [])
        if profile in load_mapping
    ]
    default_profile_selection = [
        load_mapping[profile] for profile in default_profile_keys
    ]

    # Render the selectbox with the filtered default options
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
        return [""]

    return [selected_profile_key]


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
        with open(file_path, encoding="utf-8") as csv_file:
            with temp_file:
                reader = csv.DictReader(csv_file)
                fieldnames = reader.fieldnames or []
                writer = csv.DictWriter(temp_file, fieldnames=fieldnames)
                writer.writeheader()

                # Iterate through CSV rows and update the target cell with the new value
                for index, row in enumerate(reader):
                    if str(index) == row_identifier:
                        row[column_name] = new_value

                    for key, value in row.items():
                        if value is None or value == "" or str(value).lower() == "nan":
                            row[key] = ""

                    writer.writerow(row)

        # Remove the original file and replace it with the updated temp file
        shutil.move(temp_file.name, file_path)
        return True
    except Exception:
        if os.path.exists(temp_file.name):
            os.unlink(temp_file.name)
        raise


# =============================================================================
# For input scenario config editing and saving
# =============================================================================


def format_keys_into_readable_titles(value: str) -> str:
    """Format snake_case keys into user-facing titles."""
    return value.replace("res_", "Renewable ").replace("_", " ").strip().title()


def convert_date_string_into_date_obj(raw_value: object, fallback: date) -> date:
    """Convert an ISO date string to a date object."""
    parsed_date = fallback

    if isinstance(raw_value, str):
        try:
            parsed_date = date.fromisoformat(raw_value)
        except ValueError:
            pass

    return parsed_date


def render_save_button_for_input_df(
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

        # Iterate through the edited df and update CSV file for any changes
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


def convert_scalar_value_to_yaml_value(value: str) -> object:
    """Convert a scalar string value to the appropriate YAML scalar type."""
    stripped_value = value.strip()

    # Check for scientific notation (e.g., "1.23e-4" or "-5E+6")
    if re.fullmatch(r"[+-]?\d+[eE][+-]?\d+", stripped_value):
        notation, exponent = stripped_value.lower().split("e", maxsplit=1)

        # Handle optional signs in both the notation and exponent parts
        notation_has_sign = notation.startswith(("+", "-"))
        notation_digits = notation[1:] if notation_has_sign else notation
        exponent_has_sign = exponent.startswith(("+", "-"))
        exponent_digits = exponent[1:] if exponent_has_sign else exponent

        # Validate that the notation and exponent parts contain only digits
        if not notation_digits.isdigit() or not exponent_digits.isdigit():
            raise ValueError(f"Invalid scientific notation: {stripped_value}")

        # Use ScalarFloat to dump to YAML
        return ScalarFloat(
            float(stripped_value),
            width=len(notation_digits),
            prec=-1,
            m_sign=notation_has_sign,
            exp="e",
            e_width=len(exponent_digits) + (1 if exponent_has_sign else 0),
            e_sign=exponent_has_sign,
        )

    # Check for integer values (e.g., "42" or "-7")
    if re.fullmatch(r"[+-]?\d+", stripped_value):
        return int(stripped_value)

    # Check for float values (e.g., "3.14", "-0.001", or ".5")
    if re.fullmatch(
        r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?",
        stripped_value,
    ):
        return float(stripped_value)

    return DoubleQuotedScalarString(value)


def convert_to_commented_yaml_value(value: object) -> object:
    """Convert a Python value to the corresponding ruamel YAML value."""
    if isinstance(value, str):
        return convert_scalar_value_to_yaml_value(value)

    # Initialize the root node
    root = CommentedMap() if isinstance(value, dict) else CommentedSeq()

    # Use a stack to traverse the data structure iteratively and convert scalar values
    pending_nodes: list[tuple[object, CommentedMap | CommentedSeq]] = [(value, root)]

    while pending_nodes:
        # source_node is the original Python data structure node (dict, list, or scalar)
        # target_node is the corresponding YAML node (CommentedMap or CommentedSeq)
        source_node, target_node = pending_nodes.pop()

        if isinstance(source_node, dict) and isinstance(target_node, CommentedMap):
            for key, item in source_node.items():
                if isinstance(item, dict):
                    # For dictionaries, create a CommentedMap and add it to the stack
                    child_node = CommentedMap()
                    target_node[key] = child_node
                    pending_nodes.append((item, child_node))
                elif isinstance(item, list):
                    # For lists, create a CommentedSeq and add it to the stack
                    child_node = CommentedSeq()
                    target_node[key] = child_node
                    pending_nodes.append((item, child_node))
                elif isinstance(item, str):
                    # For scalar string values, convert them to the YAML scalar type
                    target_node[key] = convert_scalar_value_to_yaml_value(item)
                else:
                    # For other types (e.g., numbers, booleans), assign them directly
                    target_node[key] = item
            continue

        if isinstance(source_node, list) and isinstance(target_node, CommentedSeq):
            for item in source_node:
                if isinstance(item, dict):
                    # For dictionaries, create a CommentedMap and add it to the stack
                    child_node = CommentedMap()
                    target_node.append(child_node)
                    pending_nodes.append((item, child_node))
                elif isinstance(item, list):
                    # For lists, create a CommentedSeq and add it to the stack
                    child_node = CommentedSeq()
                    target_node.append(child_node)
                    pending_nodes.append((item, child_node))
                elif isinstance(item, str):
                    # For scalar string values, convert them to the YAML scalar type
                    target_node.append(convert_scalar_value_to_yaml_value(item))
                else:
                    # For other types (e.g., numbers, booleans), assign them directly
                    target_node.append(item)

    return root


def save_input_config_into_yaml(
    scenario_section: dict,
    section_name: str,
    has_changes_key: str,
    message_delay: float = 1,
    change_type: str = "revert",
) -> None:
    """Save changes in the scenario config."""
    scenario_config_path = st.session_state.scenario_config_path

    # Initialize YAML instance
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.default_flow_style = False
    yaml.width = 4096  # Prevent line wrapping

    scenario_config_data = dict(st.session_state.scenario_config)
    scenario_config_data[section_name] = scenario_section
    st.session_state.scenario_config = scenario_config_data

    commented_map = CommentedMap()
    for key, item in scenario_config_data.items():
        # Convert each value to the appropriate YAML type while preserving comments
        commented_map[key] = convert_to_commented_yaml_value(item)

    scenario_config = commented_map

    # Save the updated scenario config back to the YAML file
    with open(scenario_config_path, "w", encoding="utf-8") as file_handle:
        file_handle.write(SCENARIO_CONFIG_HEADER)
        yaml.dump(scenario_config, file_handle)

    if change_type == "save":
        st.success("Changes saved successfully!")
    else:
        st.success("Changes reverted successfully!")

    st.session_state[has_changes_key] = False
    time.sleep(message_delay)
    st.rerun()


def add_key_for_save_or_revert_changes(section_name: str, widget_name: str) -> str:
    """Build a session key for saving or reverting a section's widgets."""
    reset_counter = st.session_state.get(
        f"scenario_config::{section_name}::reset_counter",
        0,
    )
    return f"scenario_config::{section_name}::{reset_counter}::{widget_name}"


def render_section_action_buttons(
    section_name: str,
    scenario_section: dict,
    edited_section: dict,
    has_changes_key: str,
) -> None:
    """Render save and revert buttons for changes in the scenario config sections."""
    has_changes = edited_section != (scenario_section or {})

    action_columns = st.columns(2)
    with action_columns[0]:
        if st.button(
            "Save Changes",
            key=add_key_for_save_or_revert_changes(section_name, "save"),
            type="primary" if has_changes else "secondary",
            disabled=not has_changes,
        ):
            save_input_config_into_yaml(
                scenario_section=edited_section,
                section_name=section_name,
                has_changes_key=has_changes_key,
                change_type="save",
            )

    with action_columns[1]:
        if st.button(
            "Revert Changes",
            key=add_key_for_save_or_revert_changes(section_name, "revert"),
            disabled=not has_changes,
        ):
            # reset every changes in the sections back to default
            reset_counter_key = f"scenario_config::{section_name}::reset_counter"
            st.session_state[reset_counter_key] = (
                st.session_state.get(reset_counter_key, 0) + 1
            )
            st.session_state[has_changes_key] = False

            save_input_config_into_yaml(
                scenario_section=scenario_section,
                section_name=section_name,
                has_changes_key=has_changes_key,
                change_type="revert",
            )
