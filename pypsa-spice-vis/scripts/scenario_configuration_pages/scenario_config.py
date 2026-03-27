# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Helpers for rendering and editing scenario configuration YAML."""

import json
import os
from datetime import date

import streamlit as st

from scripts.input_st_handler import render_save_button_for_input_config

CO2_OPTIONS = ["co2_cap", "co2_price"]
INCLUSIVE_OPTIONS = ["both", "neither", "left", "right"]
MATH_SYMBOL_OPTIONS = ["<=", ">="]
RESERVE_MARGIN_METHODS = ["static", "dynamic"]
TEMPORAL_CLUSTERING_METHODS = ["nth_hour", "clustered"]
EXCLUDED_SECTIONS = {"version", "logging", "solving"}
STATUS_MESSAGE_KEY = "scenario_config_status_message"


def _safe_widget_key(raw_key: str) -> str:
    """Create stable Streamlit widget keys from nested YAML paths."""
    return "".join(character if character.isalnum() else "_" for character in raw_key)


def _title_from_key(value: str) -> str:
    """Convert snake_case keys into user-facing titles."""
    return value.replace("_", " ").strip().title()


def _parse_iso_date(raw_value: object, fallback: date) -> date:
    """Convert an ISO date string to a date object."""
    if isinstance(raw_value, str):
        try:
            return date.fromisoformat(raw_value)
        except ValueError:
            return fallback
    return fallback


def _parse_scalar_text(raw_text: str) -> object:
    """Convert free-form text input into a YAML-friendly scalar."""
    stripped = raw_text.strip()
    if stripped == "":
        return None

    lowered = stripped.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False

    try:
        if stripped.startswith("-"):
            if stripped[1:].isdigit():
                return int(stripped)
        elif stripped.isdigit():
            return int(stripped)
        return float(stripped)
    except ValueError:
        return stripped


def _normalize_json_value(value: object) -> object:
    """Normalize parsed JSON to better match YAML types used by the model."""
    if isinstance(value, dict):
        normalized = {}
        for key, child_value in value.items():
            normalized_key = int(key) if isinstance(key, str) and key.isdigit() else key
            normalized[normalized_key] = _normalize_json_value(child_value)
        return normalized

    if isinstance(value, list):
        return [_normalize_json_value(item) for item in value]

    return value


def resolve_scenario_config_file(scenario_folder: str) -> str:
    """Resolve the scenario config file inside the selected scenario folder."""
    preferred = os.path.join(scenario_folder, "scenario_config.yaml")
    if os.path.exists(preferred):
        return preferred

    candidates = sorted(
        file_name
        for file_name in os.listdir(scenario_folder)
        if file_name.startswith("scenario_config")
        and file_name.endswith(".yaml")
        and not file_name.endswith(".default.yaml")
    )
    if not candidates:
        raise FileNotFoundError(
            f"No `scenario_config*.yaml` file found in: {scenario_folder}"
        )

    return os.path.join(scenario_folder, candidates[0])


def _render_nullable_scalar_input(
    label: str,
    value: object,
    widget_key: str,
    help_text: str | None = None,
) -> object:
    """Render a scalar input that can be left empty to preserve null."""
    if isinstance(value, bool):
        return st.checkbox(label, value=value, key=widget_key, help=help_text)

    if isinstance(value, int) and not isinstance(value, bool):
        return st.number_input(
            label, value=value, step=1, key=widget_key, help=help_text
        )

    if isinstance(value, float):
        return st.number_input(
            label,
            value=float(value),
            format="%.6f",
            key=widget_key,
            help=help_text,
        )

    raw_value = "" if value is None else str(value)
    text_value = st.text_input(label, value=raw_value, key=widget_key, help=help_text)
    return _parse_scalar_text(text_value)


def _render_json_editor(
    label: str,
    value: object,
    widget_key: str,
    help_text: str | None = None,
) -> object:
    """Render a JSON editor and preserve the previous value on parse errors."""
    default_value = {} if value is None else value
    json_text = st.text_area(
        label,
        value=json.dumps(default_value, indent=2),
        key=widget_key,
        help=help_text,
        height=180,
    )

    try:
        return _normalize_json_value(json.loads(json_text))
    except json.JSONDecodeError:
        st.warning(f"Invalid JSON in {label}. Keeping the previous value.")
        return value


def _snapshot_defaults(snapshots: dict) -> tuple[int, date, date, str]:
    """Build default snapshot values from the loaded YAML section."""
    fallback_start = date.today().replace(month=1, day=1)
    start_date = _parse_iso_date(snapshots.get("start"), fallback_start)
    end_date = _parse_iso_date(
        snapshots.get("end"),
        date(start_date.year + 1, 1, 1),
    )
    inclusive = snapshots.get("inclusive", "left")
    if inclusive not in INCLUSIVE_OPTIONS:
        inclusive = "left"
    return start_date.year, start_date, end_date, inclusive


def render_scenario_settings_section(section_value: dict) -> None:
    """Render the scenario settings editor."""
    scenario_settings = section_value or {}
    snapshots = scenario_settings.get("snapshots", {}) or {}
    resolution = scenario_settings.get("resolution", {}) or {}
    interest = scenario_settings.get("interest", {}) or {}

    default_year, default_start, default_end, default_inclusive = _snapshot_defaults(
        snapshots
    )

    with st.expander("Scenario settings", expanded=True):
        year_col, threshold_col = st.columns([1, 1])
        with year_col:
            model_year = int(
                st.number_input(
                    "Model year",
                    min_value=1900,
                    max_value=3000,
                    value=default_year,
                    step=1,
                    key="scenario_configs_model_year",
                )
            )
        with threshold_col:
            remove_threshold = st.number_input(
                "Remove asset with capacity lower than (MW)",
                min_value=0.0,
                value=float(scenario_settings.get("remove_threshold", 0.0) or 0.0),
                format="%.6f",
                key="scenario_configs_remove_threshold",
            )

        st.caption(
            "By default, the snapshot range is derived automatically from the model "
            "year and saved as a full-year hourly range."
        )
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
        with method_col:
            selected_method = st.selectbox(
                "Temporal clustering method",
                options=TEMPORAL_CLUSTERING_METHODS,
                index=TEMPORAL_CLUSTERING_METHODS.index(resolution_method),
                key="scenario_configs_resolution_method",
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

        st.markdown("#### Interest")
        edited_interest = {}
        interest_items = list(interest.items())
        interest_columns = st.columns(2) if interest_items else []
        for index, (country, country_value) in enumerate(interest_items):
            with interest_columns[index % 2]:
                edited_interest[country] = _render_nullable_scalar_input(
                    country,
                    country_value,
                    _safe_widget_key(f"scenario_configs_interest_{country}"),
                    help_text="Interest rate in decimal form, for example 0.05 for 5%.",
                )

        if end_date < start_date:
            st.error("Snapshot end must be on or after snapshot start.")
            return

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

        has_changes_key = f"has_changes_{st.title}_scenario_configs"
        has_changes = st.session_state.get(has_changes_key, False)

        render_save_button_for_input_config(
            section_value=edited_section,
            section_name="scenario_configs",
            save_button_key="save_scenario_configs",
            has_changes=has_changes,
            has_changes_key=has_changes_key,
        )


def render_co2_management_section(section_value: dict) -> None:
    """Render the CO2 management editor."""
    co2_management = section_value or {}

    with st.expander("CO2 management", expanded=True):
        st.caption(
            "Select the CO2 instrument per country and edit the year-specific values."
        )
        edited_section = {}

        for country, country_config in co2_management.items():
            country_config = country_config or {}
            option = country_config.get("option", CO2_OPTIONS[0])
            if option not in CO2_OPTIONS:
                option = CO2_OPTIONS[0]

            with st.expander(country, expanded=False):
                selected_option = st.selectbox(
                    "CO2 management option",
                    options=CO2_OPTIONS,
                    index=CO2_OPTIONS.index(option),
                    key=_safe_widget_key(f"co2_option_{country}"),
                )

                edited_values = {}
                yearly_values = country_config.get("value", {}) or {}
                yearly_items = list(yearly_values.items())
                year_columns = st.columns(3) if yearly_items else []
                for index, (year, year_value) in enumerate(yearly_items):
                    with year_columns[index % 3]:
                        edited_values[year] = _render_nullable_scalar_input(
                            str(year),
                            year_value,
                            _safe_widget_key(f"co2_value_{country}_{year}"),
                        )

                edited_section[country] = {
                    "option": selected_option,
                    "value": edited_values,
                }

        has_changes_key = f"has_changes_{st.title}_co2_management"
        has_changes = st.session_state.get(has_changes_key, False)

        render_save_button_for_input_config(
            section_value=edited_section,
            section_name="co2_management",
            save_button_key="save_co2_management",
            has_changes_key=has_changes_key,
            has_changes=has_changes,
        )


def _render_custom_constraint_field(
    constraint_name: str,
    field_name: str,
    field_value: object,
    widget_prefix: str,
) -> object:
    """Render one field inside a custom constraint block."""
    widget_key = _safe_widget_key(f"{widget_prefix}_{field_name}")
    label = "Active" if field_name == "activate" else _title_from_key(field_name)

    if field_name == "activate":
        return st.checkbox(label, value=bool(field_value), key=widget_key)

    if constraint_name == "reserve_margin" and field_name == "method":
        method = str(field_value)
        if method not in RESERVE_MARGIN_METHODS:
            method = "static"
        return st.selectbox(
            label,
            options=RESERVE_MARGIN_METHODS,
            index=RESERVE_MARGIN_METHODS.index(method),
            key=widget_key,
        )

    if field_name == "math_symbol":
        symbol = str(field_value)
        if symbol not in MATH_SYMBOL_OPTIONS:
            symbol = MATH_SYMBOL_OPTIONS[0]
        return st.selectbox(
            label,
            options=MATH_SYMBOL_OPTIONS,
            index=MATH_SYMBOL_OPTIONS.index(symbol),
            key=widget_key,
        )

    if field_name in {"value", "values"}:
        return _render_json_editor(
            f"{label} (JSON)",
            field_value,
            widget_key,
            help_text="Use JSON for technology- or year-specific mappings.",
        )

    if isinstance(field_value, dict):
        return _render_json_editor(
            f"{label} (JSON)",
            field_value,
            widget_key,
            help_text="Use JSON to edit nested mappings.",
        )

    if isinstance(field_value, list):
        return _render_json_editor(
            f"{label} (JSON)",
            field_value,
            widget_key,
            help_text="Use a JSON array for list values.",
        )

    return _render_nullable_scalar_input(label, field_value, widget_key)


def render_custom_constraints_section(section_value: dict) -> None:
    """Render the custom constraints editor."""
    custom_constraints = section_value or {}

    with st.expander("Custom constraints", expanded=True):
        st.caption(
            "Nested mappings are edited as JSON so technology- and year-specific "
            "constraint values remain flexible."
        )
        edited_section = {}

        for country, country_constraints in custom_constraints.items():
            edited_country_constraints = {}
            country_constraints = country_constraints or {}

            with st.expander(country, expanded=False):
                for constraint_name, constraint_config in country_constraints.items():
                    constraint_config = constraint_config or {}
                    with st.expander(_title_from_key(constraint_name), expanded=False):
                        if isinstance(constraint_config, dict):
                            edited_constraint = {}
                            for field_name, field_value in constraint_config.items():
                                edited_constraint[field_name] = (
                                    _render_custom_constraint_field(
                                        constraint_name,
                                        field_name,
                                        field_value,
                                        f"custom_{country}_{constraint_name}",
                                    )
                                )
                            edited_country_constraints[constraint_name] = (
                                edited_constraint
                            )
                        else:
                            edited_country_constraints[constraint_name] = (
                                _render_nullable_scalar_input(
                                    _title_from_key(constraint_name),
                                    constraint_config,
                                    _safe_widget_key(
                                        f"custom_{country}_{constraint_name}"
                                    ),
                                )
                            )

            edited_section[country] = edited_country_constraints

        has_changes_key = f"has_changes_{st.title}_custom_constraints"
        has_changes = st.session_state.get(has_changes_key, False)

        render_save_button_for_input_config(
            section_value=edited_section,
            section_name="custom_constraints",
            save_button_key="save_custom_constraints",
            has_changes_key=has_changes_key,
            has_changes=has_changes,
        )


if __name__ == "__main__":
    scenario_config = st.session_state.scenario_config

    for section_name, section_value in scenario_config.items():
        if section_name in EXCLUDED_SECTIONS:
            continue

        if section_name == "scenario_configs":
            render_scenario_settings_section(section_value)
            continue

        if section_name == "co2_management":
            render_co2_management_section(section_value)
            continue

        if section_name == "custom_constraints":
            render_custom_constraints_section(section_value)
            continue
