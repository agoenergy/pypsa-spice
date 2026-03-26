# SPDX-FileCopyrightText: PyPSA-SPICE Developers
#
# SPDX-License-Identifier: GPL-2.0-or-later

"""Helpers for rendering and editing scenario configuration YAML."""

from __future__ import annotations

import os

import streamlit as st
import yaml


def _safe_widget_key(raw_key: str) -> str:
    """Create stable Streamlit widget keys from nested YAML paths."""
    return "".join(character if character.isalnum() else "_" for character in raw_key)


def _is_primitive(value: object) -> bool:
    """Check if a value is a primitive YAML type."""
    return value is None or isinstance(value, (str, int, float, bool))


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


def load_config(file_path: str) -> dict:
    """Load a YAML configuration file."""
    with open(file_path, encoding="utf-8") as file_handle:
        return yaml.safe_load(file_handle)


def save_config(new_config: dict, file_path: str) -> bool:
    """Persist the edited YAML configuration to disk."""
    try:
        with open(file_path, "w", encoding="utf-8") as file_handle:
            yaml.safe_dump(
                new_config,
                file_handle,
                default_flow_style=False,
                sort_keys=False,
            )
        return True
    except Exception as exc:
        raise RuntimeError(f"Error saving configuration: {exc}") from exc


def edit_config_value(label: str, value: object, key_prefix: str) -> object:
    """Render editable widgets for one YAML value recursively."""
    widget_key = _safe_widget_key(f"{key_prefix}_{label}")

    if isinstance(value, bool):
        return st.checkbox(label, value=value, key=widget_key)

    if isinstance(value, int):
        return st.number_input(label, value=value, step=1, key=widget_key)

    if isinstance(value, float):
        return st.number_input(label, value=float(value), format="%.6f", key=widget_key)

    if isinstance(value, str):
        return st.text_input(label, value=value, key=widget_key)

    if value is None:
        null_text = st.text_input(
            label,
            value="",
            key=widget_key,
            help="Leave empty to keep null, or type a value to convert it to string.",
        )
        return None if null_text.strip() == "" else null_text

    if isinstance(value, list):
        if all(_is_primitive(item) for item in value):
            list_yaml = st.text_area(
                label,
                value=yaml.safe_dump(value, sort_keys=False),
                key=widget_key,
                help="Edit as YAML list.",
            )
            try:
                parsed = yaml.safe_load(list_yaml)
            except yaml.YAMLError:
                st.warning(f"Invalid YAML list in `{label}`. Keeping original value.")
                return value

            if isinstance(parsed, list):
                return parsed

            st.warning(f"`{label}` must stay a list. Keeping original value.")
            return value

        edited_items = []
        st.caption(label)
        for index, item in enumerate(value):
            edited_items.append(
                edit_config_value(f"[{index}]", item, f"{key_prefix}_{label}_{index}")
            )
        return edited_items

    if isinstance(value, dict):
        edited_dict = {}
        st.caption(label)
        for child_key, child_value in value.items():
            edited_dict[child_key] = edit_config_value(
                str(child_key),
                child_value,
                f"{key_prefix}_{label}_{child_key}",
            )
        return edited_dict

    return st.text_input(label, value=str(value), key=widget_key)


def render_section_editor(section_name: str, section_value: object) -> object:
    """Render one top-level config section in an expander."""
    with st.expander(section_name, expanded=False):
        if isinstance(section_value, dict):
            edited_section = {}
            for key, value in section_value.items():
                edited_section[key] = edit_config_value(
                    str(key),
                    value,
                    f"section_{section_name}_{key}",
                )
            return edited_section

        return edit_config_value(section_name, section_value, f"section_{section_name}")


def render_scenario_config_editor() -> None:
    """Render the scenario configuration editor page."""
    scenario_config = st.session_state.scenario_config
    scenario_config_path = st.session_state.scenario_config_path

    edited_config = {}
    excluded_sections = {"version", "logging", "solving"}
    for section_name, section_value in scenario_config.items():
        if section_name not in excluded_sections:
            edited_config[section_name] = render_section_editor(
                section_name,
                section_value,
            )

    with st.container(border=True):
        if st.button("Save Scenario Config"):
            try:
                save_config(edited_config, scenario_config_path)
                st.success("Scenario config saved successfully.")
            except Exception as exc:
                st.error(f"Error saving scenario config: {exc}")


if __name__ == "__main__":
    render_scenario_config_editor()
