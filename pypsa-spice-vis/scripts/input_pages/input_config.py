# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Config editor page under Input section.

Page shows editable input data from configuration files.
"""

import os

import pandas as pd
import streamlit as st
import yaml

st.write(
    """
# Config Editor
"""
)


def get_config_files(base_path):
    """Find configuration files from the specified base path."""
    # get files with .yaml extension in the base_path or subfolder as a list only give
    # file names
    folders = [
        f for f in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, f))
    ]
    config_files = []
    for folder in [""] + folders:
        folder_path = os.path.join(base_path, folder)
        yaml_files = [f for f in os.listdir(folder_path) if f.endswith(".yaml")]
        config_files.extend(yaml_files)
    # remove any file that startes with .
    config_files = [f for f in config_files if not f.startswith(".")]
    # remove files that end with .default.yaml
    config_files = [f for f in config_files if not f.endswith(".default.yaml")]
    return config_files


def load_config(file_path):
    """Load configuration while preserving comments and structure."""
    with open(file_path, encoding="utf-8") as f:
        return yaml.safe_load(f)


def save_config(new_config, file_path):
    """Save configuration while preserving original structure, comments, and lists."""
    try:
        with open(file_path, "w", encoding="utf-8") as f:
            yaml.safe_dump(new_config, f, default_flow_style=False, sort_keys=False)
        return True
    except Exception as e:
        raise RuntimeError(f"Error saving configuration: {str(e)}") from e


def handle_base_settings(config, file_path):
    """Handle base settings configuration with years and regions."""
    if "base_configs" not in config:
        return config

    base_config = config["base_configs"]

    with st.container(border=True):
        st.write("### Base Settings")

        # Years setting as a sequence
        st.write("Years (one per line):")
        years_text = st.text_area(
            "Edit Years",
            value="\n".join(map(str, base_config["years"])),
            help="Enter years (2019-2050), one per line. Will be saved as a sequence.",
        )

        # Convert text area input to sorted list of integers
        try:
            years_list = sorted(
                [int(year.strip()) for year in years_text.split("\n") if year.strip()]
            )
            valid_years = [year for year in years_list if 2019 <= year <= 2050]
            if len(valid_years) != len(years_list):
                st.warning(
                    "Some years were removed as they were outside the valid range "
                    + " (2019-2050)"
                )
            base_config["years"] = valid_years
        except ValueError as e:
            st.error(f"Please enter valid years (integers between 2019-2050): {str(e)}")
            base_config["years"] = base_config["years"]

        # Regions setting as a sequence
        st.write("Regions (one per line):")
        regions_text = st.text_area(
            "Edit Regions",
            value="\n".join(base_config["regions"]),
            help="Enter regions, one per line.",
        )
        regions_list = [
            region.strip() for region in regions_text.split("\n") if region.strip()
        ]
        if regions_list:
            base_config["regions"] = regions_list

        if st.button("Save Scenario Settings"):
            try:
                save_config(config, file_path)
                st.success("Scenario settings saved successfully!")
            except OSError as e:
                st.error(f"Error saving scenario settings: {e}")
            except yaml.YAMLError as e:
                st.error(f"Error formatting yaml file: {e}")

    return config


def display_co2_limits(config, years):
    """Display and edit CO2 limits for given years in the config editor."""
    st.write("CO2 Limits (Mt):")
    co2_df = pd.DataFrame(
        {
            "Year": years,
            "CO2 Limit": [
                config["co2_management"]["co2_limit"].get(str(year), 0.0)
                for year in years
            ],
        }
    )
    edited_co2_df = st.data_editor(
        co2_df,
        column_config={
            "Year": st.column_config.NumberColumn(
                "Year", min_value=2019, max_value=2050, step=1
            ),
            "CO2 Limit": st.column_config.NumberColumn(
                "CO2 Limit (Mt)", min_value=0.0, format="%.1f"
            ),
        },
        hide_index=True,
    )
    return {str(row["Year"]): row["CO2 Limit"] for _, row in edited_co2_df.iterrows()}


def display_emission_costs(config, years):
    """Display and edit emission costsfor given years in the config editor."""
    st.write("Emission Costs (€/tCO2):")
    emission_costs_df = pd.DataFrame(
        {
            "Year": years,
            "Emission Cost": [
                config["co2_management"]["emission_cost"].get(str(year), 0.0)
                for year in years
            ],
        }
    )
    edited_emission_costs_df = st.data_editor(
        emission_costs_df,
        column_config={
            "Year": st.column_config.NumberColumn(
                "Year", min_value=2019, max_value=2050, step=1
            ),
            "Emission Cost": st.column_config.NumberColumn(
                "Emission Cost (€/tCO2)", min_value=0.0, format="%.1f"
            ),
        },
        hide_index=True,
    )
    return {
        str(row["Year"]): row["Emission Cost"]
        for _, row in edited_emission_costs_df.iterrows()
    }


def handle_co2_management(config, file_path):
    """Handle CO2 management configuration with country-specific settings."""
    if "co2_management" not in config:
        return config

    management_config = config["co2_management"]
    with st.container(border=True):
        st.write("### CO2 Management Settings")

        # Get list of countries from co2_management (excluding comments)
        countries = [
            key
            for key in management_config.keys()
            if isinstance(management_config[key], dict)
        ]

        if not countries:
            st.info("No country-specific CO2 management settings found.")
            return config

        # Select country to edit
        selected_country = st.selectbox(
            "Select Country/Region for CO2 management:", countries, index=0
        )

        country_config = management_config[selected_country]

        # CO2 management option for selected country
        current_option = country_config.get("option", "co2_cap")
        co2_option = st.radio(
            f"CO2 Management Option for {selected_country}",
            ["co2_price", "co2_cap"],
            index=0 if current_option == "co2_price" else 1,
        )
        management_config[selected_country]["option"] = co2_option

        # Display and edit values for the selected country
        st.write(f"### {selected_country} CO2 Values")
        if "value" in country_config:
            value_dict = country_config["value"]
            years = list(value_dict.keys())
            values = list(value_dict.values())

            # Create dataframe for editing
            df = pd.DataFrame({"Year": [int(y) for y in years], "Value": values})

            unit = "€/tCO2" if co2_option == "co2_price" else "Mt CO2"

            edited_df = st.data_editor(
                df,
                column_config={
                    "Year": st.column_config.NumberColumn(
                        "Year", min_value=2019, max_value=2050, step=1
                    ),
                    "Value": st.column_config.NumberColumn(
                        f"Value ({unit})", min_value=0.0, format="%.1f"
                    ),
                },
                hide_index=True,
                use_container_width=True,
            )

            # Update the config with edited values
            updated_values = {
                str(row["Year"]): row["Value"] for _, row in edited_df.iterrows()
            }
            management_config[selected_country]["value"] = updated_values

        if st.button("Save CO2 Settings"):
            try:
                save_config(config, file_path)
                st.success("CO2 Management settings saved successfully!")
            except Exception as e:
                st.error(f"Error saving CO2 settings: {e}")

    return config


def handle_resolution(config, file_path):
    """Handle resolution configuration with method and interest rate."""
    if "scenario_configs" not in config:
        return config

    scenario_config = config["scenario_configs"]
    resolution_config = scenario_config["resolution"]

    with st.container(border=True):
        st.write("### Resolution Settings")

        resolution_method = st.radio(
            "Resolution Method",
            ["nth_hour", "clustered"],
            index=0 if resolution_config["method"] == "nth_hour" else 1,
        )
        resolution_config["method"] = resolution_method

        if resolution_method == "clustered":
            number_of_days = st.number_input(
                "Number of Representative Days",
                min_value=1,
                max_value=365,
                value=resolution_config.get("number_of_days", 3),
            )
            resolution_config["number_of_days"] = number_of_days

        if resolution_method == "nth_hour":
            stepsize = st.number_input(
                "Step Size",
                min_value=1,
                max_value=168,
                value=resolution_config.get("stepsize", 25),
            )
            resolution_config["stepsize"] = stepsize

        if st.button("Save Resolution Settings"):
            try:
                save_config(config, file_path)
                st.success("Resolution settings saved successfully!")
            except Exception as e:
                st.error(f"Error saving resolution settings: {e}")

        # Interest rate setting
        interest_rate = st.number_input(
            "Interest Rate",
            min_value=0.0,
            max_value=1.0,
            value=float(scenario_config["interest"]),
            format="%.3f",
            help="Interest rate value between 0 and 1",
        )
        scenario_config["interest"] = interest_rate

    return config


def handle_custom_constraints(config, file_path):
    """Handle custom constraints configuration with country-specific settings."""
    if "custom_constraints" not in config:
        return config

    custom_constraint = config["custom_constraints"]
    with st.container(border=True):
        st.write("### Custom Constraints Settings")

        # Get list of countries from co2_management (excluding comments)
        countries = [
            key
            for key in custom_constraint.keys()
            if isinstance(custom_constraint[key], dict)
        ]

        if not countries:
            st.info("No country-specific custom constraints found.")
            return config

        # Select country to edit
        selected_country = st.selectbox(
            "Select Country/Region for custom constraints:", countries, index=0
        )

        country_constraints = custom_constraint[selected_country]

        # RE Generation Constraint
        if country_constraints["res_generation"].get("activate", False):
            re_gen_constraint = country_constraints["res_generation"]
            re_gen_sense = st.radio(
                "RE Generation Sense",
                ["<=", ">=", "=="],
                index=0 if re_gen_constraint["math_symbol"] == "<=" else 1,
            )
            re_gen_constraint["math_symbol"] = re_gen_sense

            # RE Generation Fractions
            st.write("RE Generation Fractions:")
            re_gen_df = pd.DataFrame(
                {
                    "Year": list(
                        map(int, re_gen_constraint["res_generation_share"].keys())
                    ),
                    "Fraction": list(
                        re_gen_constraint["res_generation_share"].values()
                    ),
                }
            )
            edited_re_gen_df = st.data_editor(
                re_gen_df,
                column_config={
                    "Year": st.column_config.NumberColumn(
                        "Year", min_value=2020, max_value=2050, step=5
                    ),
                    "Fraction": st.column_config.NumberColumn(
                        "Fraction", min_value=0.0, max_value=1.0, format="%.2f"
                    ),
                },
                num_rows="dynamic",
                hide_index=True,
            )
            re_gen_constraint["res_generation_share"] = {
                str(row["Year"]): row["Fraction"]
                for _, row in edited_re_gen_df.iterrows()
            }

        if country_constraints["thermal_must_run"].get("activate", False):
            must_run_constraint = country_constraints["thermal_must_run"]
            must_run_frac = st.number_input(
                "Must Run Fraction",
                min_value=0.0,
                max_value=1.0,
                value=float(must_run_constraint["min_must_run_ratio"]),
                format="%.2f",
                help="Fraction of thermal capacity that must run",
            )
            must_run_constraint["min_must_run_ratio"] = must_run_frac

        # Reserve Margin
        if country_constraints["reserve_margin"].get("activate", False):
            margin_constraint = country_constraints["reserve_margin"]

            st.write("Reserve Parameters:")
            epsilon_load = st.number_input(
                "Epsilon Load",
                min_value=0.0,
                max_value=1.0,
                value=float(margin_constraint["epsilon_load"]),
                format="%.2f",
                help="Fraction of load considered as reserve",
            )
            epsilon_vre = st.number_input(
                "Epsilon VRE",
                min_value=0.0,
                max_value=1.0,
                value=float(margin_constraint["epsilon_vre"]),
                format="%.2f",
                help="Contribution of VRE to reserve",
            )
            contingency = st.number_input(
                "Contingency (MW)",
                min_value=0,
                value=int(margin_constraint["contingency"]),
                help="Extra contingency in MW",
            )

            margin_constraint.update(
                {
                    "epsilon_load": epsilon_load,
                    "epsilon_vre": epsilon_vre,
                    "contingency": contingency,
                }
            )

        if st.button("Save Custom Constraints"):
            try:
                save_config(config, file_path)
                st.success("Custom constraints saved successfully!")
            except Exception as e:
                st.error(f"Error saving custom constraints: {e}")

    return config


def get_project_folder_from_path(path):
    """Get the project folder name from the given path."""
    # get all folders that are project folders. only one folder should be there
    folders = [
        f
        for f in os.listdir(path)
        if os.path.isdir(os.path.join(path, f)) and not f.startswith(".")
    ]
    # Filter out common non-project folders
    excluded_folders = {"__pycache__", "logs", "results", "benchmarks"}
    project_folders = [f for f in folders if f not in excluded_folders]
    if len(project_folders) == 1:
        return project_folders[0]
    else:
        raise ValueError("There should be only one project folder in the path")


def config_editor():
    """Render the config editor page."""
    st.session_state.config_file = st.selectbox(
        "Select Configuration File",
        options=get_config_files(st.session_state.input_data_folder_path),
        index=0,
    )
    file_path = os.path.join(
        st.session_state.input_data_folder_path,
        get_project_folder_from_path(st.session_state.input_data_folder_path),
        st.session_state.config_file,
    )

    # Load and display current config
    config = load_config(file_path)
    st.json(config, expanded=False)

    functions = [
        handle_base_settings,
        handle_co2_management,
        handle_resolution,
        handle_custom_constraints,
    ]

    edited_config = config
    for func in functions:
        edited_config = func(edited_config, file_path)

    # Save all button at the bottom
    with st.container(border=True):
        if st.button("Save All Settings"):
            try:
                save_config(edited_config, file_path)
                st.success("All settings saved successfully!")
            except Exception as e:
                st.error(f"Error saving settings: {e}")


if __name__ == "__main__":
    config_editor()
