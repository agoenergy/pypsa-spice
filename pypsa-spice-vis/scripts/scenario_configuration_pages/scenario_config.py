# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Helpers for rendering and editing scenario configuration YAML."""

import calendar
import json
from datetime import date

import streamlit as st

from scripts.input_st_handler import (
    convert_date_string_into_date_obj,
    create_inputbox_and_keep_nulls_for_empty_input_values,
    format_keys_into_readable_titles,
    render_decimal_text_input,
    render_save_button_for_input_config,
)

CO2_OPTIONS = ["co2_cap", "co2_price"]
INCLUSIVE_OPTIONS = ["both", "neither", "left", "right"]
MATH_SYMBOL_OPTIONS = ["<=", ">="]
RESERVE_MARGIN_METHODS = ["static", "dynamic"]
TEMPORAL_CLUSTERING_METHODS = ["nth_hour", "clustered"]
EXCLUDED_SECTIONS = {"version", "logging", "solving"}

PYPSA_SPICE_DOCS_URL = "https://agoenergy.github.io/pypsa-spice/"
SCENARIO_DOCS_PATH = "getting-started/input-data/model-builder-configuration/"

# =============================================================================
# Minor helper functions for rendering scenario config sections
# =============================================================================


def extract_default_snapshot(snapshots: dict) -> tuple[int, date, date, str]:
    """Extract default snapshot values from the default scenario config file.

    Parameters
    ----------
    snapshots : dict
        The snapshots section from the scenario config file.

    Returns
    -------
    tuple[int, date, date, str]
        A tuple containing the default year, start date, end date, and inclusive option.
    """
    fallback_start = date.today().replace(month=1, day=1)
    start_date = convert_date_string_into_date_obj(
        snapshots.get("start"), fallback_start
    )
    end_date = convert_date_string_into_date_obj(
        snapshots.get("end"),
        date(start_date.year + 1, 1, 1),
    )
    inclusive = snapshots.get("inclusive", "left")
    if inclusive not in INCLUSIVE_OPTIONS:
        inclusive = "left"
    return start_date.year, start_date, end_date, inclusive


def setup_json_editor(
    label: str,
    value: object,
    constraint_key: str,
    help_text: str | None = None,
) -> object:
    """Render a text area for editing JSON content.

    Parameters
    ----------
    label : str
        The label for the JSON editor.
    value : object
        The initial value to display in the editor.
    constraint_key : str
        The key to use for the Streamlit widget.
    help_text : str | None, optional
        The help text to display for the JSON editor, by default None

    Returns
    -------
    object
        The parsed JSON object or the previous value if parsing fails.
    """
    default_value = {} if value is None else value

    # Display the JSON editor with the current value as a formatted JSON string
    json_text = st.text_area(
        label,
        value=json.dumps(default_value, indent=2),
        key=constraint_key,
        help=help_text,
        height=180,
    )

    try:
        parsed_json = json.loads(json_text)
    except json.JSONDecodeError:
        st.warning(f"Invalid JSON in {label}. Keeping the previous value.")
        return value

    if not isinstance(parsed_json, (dict, list)):
        return value

    if isinstance(parsed_json, dict):
        normalized_root = {}

    if isinstance(parsed_json, list):
        normalized_root = []

    # Use a stack to traverse the parsed JSON & build the normalized structure
    stack: list[tuple[dict | list, dict | list]] = [(parsed_json, normalized_root)]
    while stack:
        source, target = stack.pop()

        if isinstance(source, dict):
            # If the source is a dict, process key-value pairs
            if not isinstance(target, dict):
                continue
            target_dict = target
            for key, child_value in source.items():
                normalized_key = (
                    int(key) if isinstance(key, str) and key.isdigit() else key
                )
                if isinstance(child_value, dict):
                    child_target: dict | list = {}
                    target_dict[normalized_key] = child_target
                    stack.append((child_value, child_target))
                elif isinstance(child_value, list):
                    child_target = []
                    target_dict[normalized_key] = child_target
                    stack.append((child_value, child_target))
                else:
                    target_dict[normalized_key] = child_value
        else:
            # If the source is a list, process each element
            if not isinstance(target, list):
                continue
            target_list = target
            for child_value in source:
                if isinstance(child_value, dict):
                    child_target = {}
                    target_list.append(child_target)
                    stack.append((child_value, child_target))
                elif isinstance(child_value, list):
                    child_target = []
                    target_list.append(child_target)
                    stack.append((child_value, child_target))
                else:
                    target_list.append(child_value)

    return normalized_root


def create_custom_constraint_field_inputbox(
    constraint_name: str,
    field_name: str,
    field_value: object,
    widget_prefix: str,
) -> object:
    """Render an input box for a custom constraint field.

    Parameters
    ----------
    constraint_name : str
        The name of the custom constraint.
    field_name : str
        The name of the field within the custom constraint.
    field_value : object
        The current value of the field.
    widget_prefix : str
        The prefix to use for the Streamlit widget key.

    Returns
    -------
    object
        The Streamlit widget for the custom constraint field.
    """
    constraint_key = f"{widget_prefix}_{field_name}"

    # Generate a user-friendly label for the field names
    label = (
        "Active"
        if field_name == "activate"
        else format_keys_into_readable_titles(field_name)
    )

    # Render a checkbox for "activate" fields, which are expected to be boolean
    if field_name == "activate":
        return st.checkbox(label, value=bool(field_value), key=constraint_key)

    options = None
    selected_value = str(field_value)
    if constraint_name == "reserve_margin" and field_name == "method":
        options = RESERVE_MARGIN_METHODS
        if selected_value not in options:
            selected_value = "static"
    elif field_name == "math_symbol":
        options = MATH_SYMBOL_OPTIONS
        if selected_value not in options:
            selected_value = options[0]

    if options is not None:
        return st.selectbox(
            label,
            options=options,
            index=options.index(selected_value),
            key=constraint_key,
        )

    if field_name in {"value", "values"} or isinstance(field_value, (dict, list)):
        help_text = "Use JSON for technology- or year-specific mappings."
        if isinstance(field_value, dict):
            help_text = "Use JSON to edit nested mappings."
        elif isinstance(field_value, list):
            help_text = "Use a JSON array for list values."

        return setup_json_editor(
            f"{label} (JSON)",
            field_value,
            constraint_key,
            help_text=help_text,
        )

    return create_inputbox_and_keep_nulls_for_empty_input_values(
        label, field_value, constraint_key
    )


def render_country_custom_constraints(country: str, country_constraints: dict) -> dict:
    """Render all custom constraints for one country and return edited values."""
    edited_country_constraints = {}

    with st.expander(country, expanded=False):
        for constraint_name, constraint_config in country_constraints.items():
            constraint_label = format_keys_into_readable_titles(constraint_name)
            constraint_key_prefix = f"custom_{country}_{constraint_name}"
            constraint_expander_key = f"{constraint_key_prefix}_expanded"

            with st.expander(
                constraint_label,
                expanded=False,
                key=constraint_expander_key,
                on_change="rerun",
            ):
                if not isinstance(constraint_config, dict):
                    edited_country_constraints[constraint_name] = (
                        create_inputbox_and_keep_nulls_for_empty_input_values(
                            constraint_label,
                            constraint_config,
                            constraint_key_prefix,
                        )
                    )
                    continue

                # Constraint activation: collapsed=False, expanded=True.
                edited_constraint: dict[str, object] = {
                    "activate": bool(
                        st.session_state.get(constraint_expander_key, False)
                    )
                }
                for field_name, field_value in constraint_config.items():
                    if field_name == "activate":
                        continue
                    edited_constraint[field_name] = (
                        create_custom_constraint_field_inputbox(
                            constraint_name,
                            field_name,
                            field_value,
                            constraint_key_prefix,
                        )
                    )

            edited_country_constraints[constraint_name] = edited_constraint

    return edited_country_constraints


# =============================================================================
# Render each scenario config section
# ==============================================================================


def render_scenario_settings_section(scenario_section: dict) -> None:
    """Extract default scenario setting values from the default config file.

    After extraction, render the editor and update the values.

    Parameters
    ----------
    scenario_section : dict
        Default value for scenario settings section from the default config file.
    """
    # Extract default scenario setting values from default config file
    scenario_settings = scenario_section or {}
    snapshots = scenario_settings.get("snapshots", {}) or {}
    resolution = scenario_settings.get("resolution", {}) or {}
    interest = scenario_settings.get("interest", {}) or {}

    # Extract default snapshot values based on the loaded configuration
    default_year, default_start, default_end, default_inclusive = (
        extract_default_snapshot(snapshots)
    )

    # Render the scenario settings editor with inputs
    with st.expander("Scenario settings", expanded=True):
        SETTING_PATH = (
            PYPSA_SPICE_DOCS_URL
            + SCENARIO_DOCS_PATH
            + "#scenario_configyaml-scenario-settings"
        )
        st.markdown(
            "Detailed explanation can be found in: "
            f"[config guides: global input]({SETTING_PATH})"
        )
        year_col, threshold_col = st.columns([1, 1])
        # Render model year and remove threshold inputs side by side
        with year_col:
            model_year = int(
                st.number_input(
                    "Model year",
                    min_value=2010,
                    max_value=3000,
                    value=default_year,
                    step=1,
                    key="scenario_configs_model_year",
                )
            )
            # Automatically adjust model year if it's a leap year
            # to avoid issues with snapshot ranges
            if calendar.isleap(model_year):
                model_year -= 1
                st.warning(
                    f"⚠️ The selected year is a leap year. "
                    f"Model year has been automatically adjusted to **{model_year}**.",
                )
        with threshold_col:
            remove_threshold = render_decimal_text_input(
                "Remove asset with capacity lower than (MW)",
                value=float(scenario_settings.get("remove_threshold", 0.0) or 0.0),
                constraint_key="scenario_configs_remove_threshold",
            )
            remove_threshold = max(remove_threshold or 0.0, 0.0)

        st.caption(
            "By default, the snapshot range is derived automatically from the model "
            "year and saved as a full-year hourly range."
        )

        # Allow users to edit snapshots manually
        edit_snapshots_manually = st.toggle(
            "Edit snapshot range manually",
            value=False,
            key="scenario_configs_manual_snapshots",
        )

        if edit_snapshots_manually:
            start_col, end_col, inclusive_col = st.columns(3)
            with start_col:
                start_date = st.date_input(
                    "Snapshot start",
                    value=default_start,
                    key="scenario_configs_snapshot_start",
                )
            with end_col:
                end_date = st.date_input(
                    "Snapshot end",
                    value=default_end,
                    key="scenario_configs_snapshot_end",
                )
            with inclusive_col:
                inclusive = st.selectbox(
                    "Inclusive",
                    options=INCLUSIVE_OPTIONS,
                    index=INCLUSIVE_OPTIONS.index(default_inclusive),
                    key="scenario_configs_snapshot_inclusive",
                )
        else:
            start_date = date(model_year, 1, 1)
            end_date = date(model_year + 1, 1, 1)
            inclusive = "left"
            st.caption(
                f"Saved snapshot range: {start_date.isoformat()} to "
                f"{end_date.isoformat()} with inclusive='{inclusive}'."
            )

        resolution_method = resolution.get("method", TEMPORAL_CLUSTERING_METHODS[0])
        if resolution_method not in TEMPORAL_CLUSTERING_METHODS:
            resolution_method = TEMPORAL_CLUSTERING_METHODS[0]

        method_col, detail_col = st.columns(2)
        # Render temporal clustering method selection
        with method_col:
            selected_method = st.selectbox(
                "Temporal clustering method",
                options=TEMPORAL_CLUSTERING_METHODS,
                index=TEMPORAL_CLUSTERING_METHODS.index(resolution_method),
                key="scenario_configs_resolution_method",
                help=(
                    "**nth_hour** - Select every Nth snapshot as representative\n\n"
                    "**clustered** - TSAM clustering for variable-duration segments"
                ),
            )

        number_of_days = int(resolution.get("number_of_days", 3) or 3)
        stepsize = int(resolution.get("stepsize", 25) or 25)
        with detail_col:
            if selected_method == "clustered":
                number_of_days = int(
                    st.number_input(
                        "Number of days",
                        min_value=1,
                        value=number_of_days,
                        step=1,
                        key="scenario_configs_number_of_days",
                    )
                )
            else:
                stepsize = int(
                    st.number_input(
                        "Step size",
                        min_value=1,
                        value=stepsize,
                        step=1,
                        key="scenario_configs_stepsize",
                    )
                )

        # Render interest rate inputs for each country
        st.markdown("#### Interest")
        edited_interest = {}
        interest_items = list(interest.items())
        interest_columns = st.columns(2) if interest_items else []

        # Render the country-specific interest values in two columns
        for index, (country, country_value) in enumerate(interest_items):
            with interest_columns[index % 2]:
                edited_interest[country] = (
                    create_inputbox_and_keep_nulls_for_empty_input_values(
                        country,
                        country_value,
                        f"scenario_configs_interest_{country}",
                        help_text="Interest rate in decimal form, e.g. 0.05 for 5%.",
                    )
                )
        if end_date < start_date:
            st.error("Snapshot end must be on or after snapshot start.")
            return

        # Compile the edited values into a new section dictionary
        edited_section = {
            "snapshots": {
                "start": start_date.isoformat(),
                "end": end_date.isoformat(),
                "inclusive": inclusive,
            },
            "resolution": {
                "method": selected_method,
                "number_of_days": number_of_days,
                "stepsize": stepsize,
            },
            "interest": edited_interest,
            "remove_threshold": remove_threshold,
        }

        # Check if there are changes to save and render the save button
        has_changes_key = f"has_changes_{st.title}_scenario_configs"
        has_changes = st.session_state.get(has_changes_key, False)

        render_save_button_for_input_config(
            scenario_section=edited_section,
            section_name="scenario_configs",
            save_button_key="save_scenario_configs",
            has_changes=has_changes,
            has_changes_key=has_changes_key,
        )


def render_co2_management_section(scenario_section: dict) -> None:
    """Render the CO2 management editor.

    Parameters
    ----------
    scenario_section : dict
        The current values for the scenario settings section.
    """
    # Set default CO2 management values based on the loaded configuration
    co2_management = scenario_section or {}

    # Render the CO2 management editor with inputs
    with st.expander("CO2 management", expanded=True):
        MANDATORY_SETTING_PATH = (
            PYPSA_SPICE_DOCS_URL
            + SCENARIO_DOCS_PATH
            + "#scenario_configyaml-mandatory-constraints"
        )
        st.markdown(
            "Detailed explanation can be found in: "
            f"[config guides: mandatory constraints]({MANDATORY_SETTING_PATH})"
        )

        st.caption(
            "Select the CO2 instrument per country and edit the year-specific values."
        )
        edited_section = {}

        for country, country_config in co2_management.items():
            country_config = country_config or {}
            option = country_config.get("option", CO2_OPTIONS[0])
            if option not in CO2_OPTIONS:
                option = CO2_OPTIONS[0]

            # Render the CO2 management options and values for each country
            with st.expander(country, expanded=False):
                # Render the CO2 management option selection
                selected_option = st.selectbox(
                    "CO2 management option",
                    options=CO2_OPTIONS,
                    index=CO2_OPTIONS.index(option),
                    key=f"co2_option_{country}",
                    help=(
                        "**co2_cap** - Maximum allowable CO2 emissions "
                        "(carbon budget constraint)\n\n"
                        "**co2_price** - Cost applied per unit of CO2 emitted "
                        "(carbon pricing mechanism)"
                    ),
                )

                edited_values = {}
                yearly_values = country_config.get("value", {}) or {}
                yearly_items = list(yearly_values.items())
                year_columns = st.columns(3) if yearly_items else []

                # Render the year-specific CO2 values in three columns
                for index, (year, year_value) in enumerate(yearly_items):
                    with year_columns[index % 3]:
                        edited_values[year] = (
                            create_inputbox_and_keep_nulls_for_empty_input_values(
                                str(year),
                                year_value,
                                f"co2_value_{country}_{year}",
                            )
                        )

                # Compile the edited values into a new country config dictionary
                edited_section[country] = {
                    "option": selected_option,
                    "value": edited_values,
                }

        # Check if there are changes to save and render the save button
        has_changes_key = f"has_changes_{st.title}_co2_management"
        has_changes = st.session_state.get(has_changes_key, False)

        render_save_button_for_input_config(
            scenario_section=edited_section,
            section_name="co2_management",
            save_button_key="save_co2_management",
            has_changes_key=has_changes_key,
            has_changes=has_changes,
        )


def render_custom_constraints_section(scenario_section: dict) -> None:
    """Render the custom constraints editor.

    Parameters
    ----------
    scenario_section : dict
        The current values for the scenario settings section.
    """
    # Set default custom constraints values based on the loaded configuration
    custom_constraints = scenario_section or {}

    # Render the custom constraints editor with inputs
    with st.expander("Custom constraints", expanded=True):
        CUSTOM_SETTING_PATH = (
            PYPSA_SPICE_DOCS_URL
            + SCENARIO_DOCS_PATH
            + "#scenario_configyaml-custom-constraints"
        )
        st.markdown(
            "Detailed explanation can be found in: "
            f"[config guides: custom constraints]({CUSTOM_SETTING_PATH})"
        )

        st.caption(
            "Mappings are configured in JSON format to enable easy customization of "
            "technology- and year-specific constraint values."
        )
        edited_section = {}

        # Render the custom constraints for each country
        for country, country_constraints in custom_constraints.items():
            edited_section[country] = render_country_custom_constraints(
                country,
                country_constraints or {},
            )

        # Check if there are changes to save and render the save button
        has_changes_key = f"has_changes_{st.title}_custom_constraints"
        has_changes = st.session_state.get(has_changes_key, False)

        render_save_button_for_input_config(
            scenario_section=edited_section,
            section_name="custom_constraints",
            save_button_key="save_custom_constraints",
            has_changes_key=has_changes_key,
            has_changes=has_changes,
        )


if __name__ == "__main__":
    st.title(":material/settings: Scenario config editor")
    st.markdown(
        "Scenario name: "
        + st.session_state.base_config["path_configs"]["input_scenario_name"]
    )
    DOCS_PATH = "getting-started/input-data/model-builder-configuration"
    st.markdown(
        "Detailed explanation can be found in: "
        f"[config guides](https://agoenergy.github.io/pypsa-spice/{DOCS_PATH})"
    )

    # Display the config file for the selected scenario.
    st.caption(st.session_state.scenario_config_path)

    scenario_config = st.session_state.scenario_config

    for section_name, scenario_section in scenario_config.items():
        if section_name in EXCLUDED_SECTIONS:
            continue

        if section_name == "scenario_configs":
            render_scenario_settings_section(scenario_section)
        elif section_name == "co2_management":
            render_co2_management_section(scenario_section)
        else:  # section_name == "custom_constraints"
            render_custom_constraints_section(scenario_section)
