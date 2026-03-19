# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Create a generic editor for scenario_config files."""

import os

import streamlit as st
import yaml


def _safe_widget_key(raw_key):
    """Create stable Streamlit widget keys from nested yaml paths."""
    return "".join(ch if ch.isalnum() else "_" for ch in raw_key)


def _is_primitive(value):
    """Check if value is a primitive yaml type."""
    return value is None or isinstance(value, (str, int, float, bool))


def resolve_scenario_config_file(scenario_folder):
    """Find scenario_config yaml file inside the selected scenario folder."""
    preferred = os.path.join(scenario_folder, "scenario_config.yaml")
    if os.path.exists(preferred):
        return preferred

    candidates = sorted(
        f
        for f in os.listdir(scenario_folder)
        if f.startswith("scenario_config")
        and f.endswith(".yaml")
        and not f.endswith(".default.yaml")
    )

    if not candidates:
        raise FileNotFoundError(
            f"No `scenario_config*.yaml` file found in: {scenario_folder}"
        )

    return os.path.join(scenario_folder, candidates[0])


def load_config(file_path):
    """Load yaml configuration."""
    with open(file_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_config(new_config, file_path):
    """Save yaml configuration."""
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(new_config, f, default_flow_style=False, sort_keys=False)
        return True
    except Exception as e:
        raise RuntimeError(f"Error saving configuration: {str(e)}") from e


def edit_config_value(label, value, key_prefix):
    """Render editable widget(s) for one yaml value recursively."""
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
            item_label = f"[{index}]"
            edited_items.append(
                edit_config_value(item_label, item, f"{key_prefix}_{label}_{index}")
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

    raw_value = st.text_input(label, value=str(value), key=widget_key)
    return raw_value


def render_section_editor(section_name, section_value):
    """Render one top-level section in an expander."""
    with st.expander(f"{section_name}", expanded=False):
        if isinstance(section_value, dict):
            edited_section = {}
            for key, value in section_value.items():
                edited_section[key] = edit_config_value(
                    str(key), value, f"section_{section_name}_{key}"
                )
            return edited_section

        return edit_config_value(section_name, section_value, f"section_{section_name}")


if __name__ == "__main__":
    scenario_config = st.session_state.scenario_config
    scenario_config_path = st.session_state.scenario_config_path

    edited_config = {}
    exluded_sections = ["version", "logging", "solving"]
    for section_name, section_value in scenario_config.items():
        if section_name not in exluded_sections:
            edited_config[section_name] = render_section_editor(
                section_name, section_value
            )

    with st.container(border=True):
        if st.button("Save Scenario Config"):
            try:
                save_config(edited_config, scenario_config_path)
                st.success("Scenario config saved successfully.")
            except Exception as exc:
                st.error(f"Error saving scenario config: {exc}")
