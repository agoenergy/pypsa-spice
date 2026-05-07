# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Helpers for rendering and editing scenario configuration YAML."""

import calendar
from copy import deepcopy
from datetime import date
from typing import Any

import streamlit as st
from streamlit_extras.json_editor import json_editor

from scripts.input_st_handler import (
    convert_date_string_into_date_obj,
    format_keys_into_readable_titles,
    get_default_scenario_config,
    render_save_button_for_input_config,
)

CO2_OPTIONS = ["co2_cap", "co2_price"]
INCLUSIVE_OPTIONS = ["both", "neither", "left", "right"]
MATH_SYMBOL_OPTIONS = ["<=", ">="]
RESERVE_MARGIN_METHODS = ["static", "dynamic"]
TEMPORAL_CLUSTERING_METHODS = ["nth_hour", "clustered"]
EXCLUDED_SECTIONS = {"version", "logging", "solving"}
MISSING = object()

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


def normalize_json_like_value(value: Any) -> Any:
    """Normalize JSON-like data while restoring integer-looking mapping keys."""
    if not isinstance(value, (dict, list)):
        return value

    root = {} if isinstance(value, dict) else []
    pending_nodes: list[tuple[Any, Any]] = [(value, root)]

    while pending_nodes:
        source_node, target_node = pending_nodes.pop()

        if isinstance(source_node, dict) and isinstance(target_node, dict):
            for key, child_value in source_node.items():
                normalized_key = (
                    int(key) if isinstance(key, str) and key.isdigit() else key
                )
                if isinstance(child_value, dict):
                    child_target = {}
                    target_node[normalized_key] = child_target
                    pending_nodes.append((child_value, child_target))
                elif isinstance(child_value, list):
                    child_target = []
                    target_node[normalized_key] = child_target
                    pending_nodes.append((child_value, child_target))
                else:
                    target_node[normalized_key] = child_value
            continue

        if isinstance(source_node, list) and isinstance(target_node, list):
            for child_value in source_node:
                if isinstance(child_value, dict):
                    child_target = {}
                    target_node.append(child_target)
                    pending_nodes.append((child_value, child_target))
                elif isinstance(child_value, list):
                    child_target = []
                    target_node.append(child_target)
                    pending_nodes.append((child_value, child_target))
                else:
                    target_node.append(child_value)

    return root


def merge_nested_mappings(base_value: Any, override_value: Any) -> Any:
    """Recursively merge two nested mappings, preferring override values."""
    if not (isinstance(base_value, dict) and isinstance(override_value, dict)):
        return deepcopy(override_value)

    merged_value = deepcopy(base_value)
    pending_nodes: list[tuple[dict[str, Any], dict[str, Any]]] = [
        (merged_value, override_value)
    ]

    while pending_nodes:
        current_target, current_override = pending_nodes.pop()

        for key, value in current_override.items():
            existing_value = current_target.get(key)
            if isinstance(existing_value, dict) and isinstance(value, dict):
                pending_nodes.append((existing_value, value))
            else:
                current_target[key] = deepcopy(value)

    return merged_value


def count_nested_items(value: Any) -> int:
    """Return a simple complexity score for preferring richer templates."""
    item_count = 0
    pending_nodes = [value]

    while pending_nodes:
        current_value = pending_nodes.pop()
        item_count += 1

        if isinstance(current_value, dict):
            pending_nodes.extend(current_value.values())
        elif isinstance(current_value, list):
            pending_nodes.extend(current_value)

    return item_count


def build_custom_constraint_templates(
    default_custom_constraints: dict[str, Any],
) -> dict[str, Any]:
    """Build a template per custom constraint from the default scenario config."""
    constraint_templates: dict[str, Any] = {}

    for country_constraints in default_custom_constraints.values():
        if not isinstance(country_constraints, dict):
            continue

        for constraint_name, constraint_config in country_constraints.items():
            normalized_config = normalize_json_like_value(constraint_config)
            if isinstance(normalized_config, dict) and "activate" in normalized_config:
                normalized_config["activate"] = False
            current_template = constraint_templates.get(constraint_name)
            if current_template is None or count_nested_items(
                normalized_config
            ) > count_nested_items(current_template):
                constraint_templates[constraint_name] = deepcopy(normalized_config)

    return constraint_templates


def get_available_custom_constraints(
    current_custom_constraints: dict[str, Any],
) -> tuple[list[str], dict[str, Any]]:
    """Return country names and merged custom constraint templates."""
    default_scenario_config = get_default_scenario_config()
    default_custom_constraints = (
        default_scenario_config.get("custom_constraints", {}) or {}
    )
    constraint_templates = build_custom_constraint_templates(default_custom_constraints)

    countries = list(current_custom_constraints)
    if not countries:
        countries = list(default_custom_constraints)

    return countries, constraint_templates


def build_custom_constraint_editor_value(
    template_config: dict[str, Any],
    current_config: dict[str, Any],
) -> dict[str, Any]:
    """Build the editor state without letting template values override scenario data."""
    if current_config:
        return merge_nested_mappings(template_config, current_config)

    return deepcopy(template_config)


def strip_template_defaults_from_edit(
    edited_value: Any,
    current_value: Any = MISSING,
    template_value: Any = MISSING,
) -> Any:
    """Keep scenario values while dropping untouched defaults."""
    if not isinstance(edited_value, dict):
        current_exists = current_value is not MISSING
        template_exists = template_value is not MISSING

        if isinstance(edited_value, list):
            if current_exists or not template_exists or edited_value != template_value:
                return edited_value
            return MISSING

        if current_exists or not template_exists or edited_value != template_value:
            return edited_value

        return MISSING

    root_result: dict[str, Any] = {}
    pending_nodes: list[tuple[dict[str, Any], Any, Any, dict[str, Any]]] = [
        (edited_value, current_value, template_value, root_result)
    ]

    while pending_nodes:
        edited_dict, current_node, template_node, result_node = pending_nodes.pop()
        current_dict = current_node if isinstance(current_node, dict) else {}
        template_dict = template_node if isinstance(template_node, dict) else {}

        for key, child_value in edited_dict.items():
            child_current = current_dict.get(key, MISSING)
            child_template = template_dict.get(key, MISSING)
            current_exists = child_current is not MISSING
            template_exists = child_template is not MISSING

            if isinstance(child_value, dict):
                child_result: dict[str, Any] = {}
                result_node[key] = child_result
                pending_nodes.append(
                    (child_value, child_current, child_template, child_result)
                )
                continue

            if isinstance(child_value, list):
                if (
                    current_exists
                    or not template_exists
                    or child_value != child_template
                ):
                    result_node[key] = child_value
                continue

            if current_exists or not template_exists or child_value != child_template:
                result_node[key] = child_value

    cleanup_stack: list[
        tuple[dict[str, Any], Any, Any, dict[str, Any] | None, str | None]
    ] = [(root_result, current_value, template_value, None, None)]
    postorder_nodes: list[
        tuple[dict[str, Any], Any, Any, dict[str, Any] | None, str | None]
    ] = []

    while cleanup_stack:
        result_node, current_node, template_node, parent_node, parent_key = (
            cleanup_stack.pop()
        )
        postorder_nodes.append(
            (result_node, current_node, template_node, parent_node, parent_key)
        )
        current_dict = current_node if isinstance(current_node, dict) else {}
        template_dict = template_node if isinstance(template_node, dict) else {}

        for key, child_result in result_node.items():
            if isinstance(child_result, dict):
                cleanup_stack.append(
                    (
                        child_result,
                        current_dict.get(key, MISSING),
                        template_dict.get(key, MISSING),
                        result_node,
                        key,
                    )
                )

    root_removed = False
    for result_node, current_node, template_node, parent_node, parent_key in reversed(
        postorder_nodes
    ):
        current_exists = current_node is not MISSING
        template_exists = template_node is not MISSING

        if result_node or current_exists:
            continue

        if parent_node is None:
            root_removed = True
            continue

        if parent_key is None:
            continue

        if (
            template_exists
            and isinstance(template_node, dict)
            and result_node == template_node
        ):
            parent_node.pop(parent_key, None)
            continue

        parent_node.pop(parent_key, None)

    if root_removed and not root_result:
        return MISSING

    return root_result


def setup_json_editor(
    label: str,
    value: object,
    constraint_key: str,
    help_text: str | None = None,
    collapsed: bool | int = False,
    root_name: str | None = None,
) -> object:
    """Render a JSON editor for editing nested configuration content.

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
    collapsed : bool | int, optional
        Initial collapse state for nested values, by default False
    root_name : str | None, optional
        Root label shown by the JSON editor, by default None

    Returns
    -------
    object
        The edited JSON-compatible value.
    """
    default_value: dict[str, Any] | list[Any] | str
    if isinstance(value, (dict, list, str)):
        default_value = value
    else:
        default_value = {}

    st.markdown(f"**{label}**")
    if help_text:
        st.caption(help_text)

    editor_state = json_editor(
        default_value,
        name=root_name,
        collapsed=collapsed,
        display_object_size=False,
        editable=True,
        key=constraint_key,
    )

    return normalize_json_like_value(editor_state.get("data", default_value))


def render_interest_section(interest: dict[str, Any]) -> dict[str, Any]:
    """Render interest-rate mappings in a JSON editor."""
    edited_interest = setup_json_editor(
        "Interest rates",
        interest or {},
        "scenario_configs_interest",
        help_text="Edit country-specific interest rates in decimals, e.g. 0.05 for 5%.",
        collapsed=1,
        root_name=None,
    )
    return edited_interest if isinstance(edited_interest, dict) else interest


def render_co2_country_editor(
    country: str, country_config: dict[str, Any]
) -> dict[str, Any]:
    """Render one country's CO2 management configuration."""
    option = country_config.get("option", CO2_OPTIONS[0])
    if option not in CO2_OPTIONS:
        option = CO2_OPTIONS[0]

    with st.expander(country, expanded=False):
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
        edited_values = setup_json_editor(
            "Year-specific values",
            country_config.get("value", {}) or {},
            f"co2_value_{country}",
            help_text="Edit year-specific CO2 cap or price values as a mapping.",
            collapsed=1,
            root_name=None,
        )

    return {
        "option": selected_option,
        "value": edited_values if isinstance(edited_values, dict) else {},
    }


def render_country_custom_constraints(
    country: str,
    country_constraints: dict[str, Any],
    constraint_templates: dict[str, Any],
) -> dict[str, Any]:
    """Render all available custom constraints for one country."""
    edited_country_constraints = {}
    ordered_constraint_names = list(constraint_templates)
    for constraint_name in country_constraints:
        if constraint_name not in constraint_templates:
            ordered_constraint_names.append(constraint_name)

    with st.expander(country, expanded=False):
        for constraint_name in ordered_constraint_names:
            constraint_label = format_keys_into_readable_titles(constraint_name)
            constraint_key = f"custom_{country}_{constraint_name}"
            template_config = deepcopy(constraint_templates.get(constraint_name, {}))
            current_config = country_constraints.get(constraint_name, {})
            editor_config = build_custom_constraint_editor_value(
                template_config,
                current_config,
            )

            with st.expander(constraint_label, expanded=False):
                edited_constraint = setup_json_editor(
                    "Configuration",
                    editor_config,
                    constraint_key,
                    help_text=(
                        "Edit the full constraint configuration as JSON. "
                        "All available constraints are loaded from the default "
                        "scenario configuration."
                    ),
                    collapsed=False,
                    root_name=None,
                )

                filtered_constraint = strip_template_defaults_from_edit(
                    edited_constraint,
                    current_config if current_config else MISSING,
                    template_config if template_config else MISSING,
                )
                edited_country_constraints[constraint_name] = (
                    {} if filtered_constraint is MISSING else filtered_constraint
                )

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
            remove_threshold = st.number_input(
                "Remove asset with capacity lower than (MW)",
                min_value=0.0,
                value=float(scenario_settings.get("remove_threshold", 0.0) or 0.0),
                format="%.2f",
                step=0.1,
                key="scenario_configs_remove_threshold",
            )

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
        edited_interest = render_interest_section(interest)
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
        st.caption(
            "Select the CO2 instrument per country and edit the year-specific values."
        )
        edited_section = {}

        for country, country_config in co2_management.items():
            edited_section[country] = render_co2_country_editor(
                country,
                country_config or {},
            )

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
    countries, constraint_templates = get_available_custom_constraints(
        custom_constraints
    )

    # Render the custom constraints editor with inputs
    with st.expander("Custom constraints", expanded=True):
        st.caption(
            "Constraint templates are loaded from the default scenario configuration. "
            "Edit each constraint as JSON inside its expander."
        )
        edited_section = {}

        # Render the custom constraints for each country
        for country in countries:
            country_result = render_country_custom_constraints(
                country,
                custom_constraints.get(country, {}) or {},
                constraint_templates,
            )
            edited_section[country] = {
                constraint_name: constraint_value
                for constraint_name, constraint_value in country_result.items()
                if constraint_value not in ({}, [], None)
            }

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
