# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Emissions page under Results section.

Page shows editable emission related dataframes and visualisations
from the modelling results.
"""

import streamlit as st

from scripts.data_utils import (
    add_nice_names,
    clean_df_for_plotting,
    normalize_dataframe,
    prepare_y_range,
    read_result_csv,
)
from scripts.output_st_handler import (
    generate_sidebar,
    render_chart_layout,
    render_download_with_data_table,
    render_section_header,
    setup_country_filter,
)
from scripts.plot_functions import (
    create_nice_names_and_color_mapping,
    plot_simple_bar_yearly,
)

# =========================== Render functions for each chart =========================


def render_e1_pow_emi_by_carrier(graph_config: dict) -> None:
    """Render power emission by carrier (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_e2_ind_emi_by_carrier(graph_config: dict) -> None:
    """Render industry emission by carrier (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


def render_e3_ene_emi_by_carrier_by_sector(graph_config: dict) -> None:
    """Render emission by carrier for all sector (yearly bar chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "simple_bar_yearly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and process scenario 1 data
    scenario_1_raw = read_result_csv(
        st.session_state.output_sce1, table_name, country=graph_config["shared_country"]
    )
    if scenario_1_raw is None or scenario_1_raw.empty:
        st.divider()
        return
    scenario_1_grouped = clean_df_for_plotting(legend_col, scenario_1_raw)
    scenario_1_grouped = scenario_1_grouped.groupby(
        ["year", legend_col], as_index=False
    )["value"].sum()
    scenario_1_grouped = normalize_dataframe(scenario_1_grouped)
    scenario_1_grouped = add_nice_names(scenario_1_grouped, legend_col, mapping_df)

    # Load and process scenario 2 data (if dual mode)
    scenario_2_grouped = None
    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = read_result_csv(
            st.session_state.output_sce2,
            table_name,
            country=graph_config["shared_country"],
        )
        if scenario_2_raw is not None:
            scenario_2_grouped = clean_df_for_plotting(legend_col, scenario_2_raw)
            scenario_2_grouped = scenario_2_grouped.groupby(
                ["year", legend_col], as_index=False
            )["value"].sum()
            scenario_2_grouped = normalize_dataframe(scenario_2_grouped)
            scenario_2_grouped = add_nice_names(
                scenario_2_grouped, legend_col, mapping_df
            )

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_grouped, scenario_2_grouped, "year")

    # Render charts (only dual if scenario 2 data is available)
    has_dual_data = (
        is_dual and scenario_2_grouped is not None and scenario_2_raw is not None
    )
    render_chart_layout(
        is_dual_scenario=has_dual_data,
        scenario_1_vis_display_data=scenario_1_grouped,
        scenario_2_vis_display_data=scenario_2_grouped,
        table_1_display_data=scenario_1_raw,
        table_2_display_data=scenario_2_raw,
        config_dict=graph_config,
        mapping_df=mapping_df,
        y_range=y_range,
        plot_function=plot_simple_bar_yearly,
        render_download_function=render_download_with_data_table,
        key=(
            "plotly_chart_"
            f"{graph_config['download_id'].format(st.session_state.output_sce1)}"
        ),
    )

    st.divider()


if __name__ == "__main__":
    st.title(":material/thermostat: Emissions")
    docs_path = "visualisation-tool/vis-sections-and-charts/#emissions"
    st.markdown(
        "Detailed explanation can be found in: "
        f"[emissions guides](https://agoenergy.github.io/pypsa-spice/{docs_path})"
    )
    output_config = st.session_state.output_config["emissions"]
    # Always render power emission chart
    render_e1_pow_emi_by_carrier(output_config["e1"])

    EMISSION_CHART_KEYS = ["e1"]

    # Only show e2 and e3 chart if industry in session_state.sector
    show_industry = "i" in str(st.session_state.get("sector", "")).lower()
    if show_industry:
        render_e2_ind_emi_by_carrier(output_config["e2"])
        render_e3_ene_emi_by_carrier_by_sector(output_config["e3"])
        EMISSION_CHART_KEYS.extend(["e2", "e3"])

    # Generate table of content side bar
    emi_charts = [
        output_config[key] for key in EMISSION_CHART_KEYS if key in output_config
    ]
    table_of_content = [c["name"] for c in emi_charts]
    generate_sidebar(table_of_content)
