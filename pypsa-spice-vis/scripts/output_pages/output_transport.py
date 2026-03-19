# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Transport page under Results section.

Page shows editable transport related
dataframes and visualisations from the modelling results.
"""

import streamlit as st

from scripts.data_utils import (
    add_nice_names,
    filter_dataframe_by_date_range,
    get_filtered_df_and_date_range,
    handle_small_values,
    load_and_validate_hourly_data,
    prepare_y_range,
    prettify_label,
)
from scripts.output_st_handler import (
    generate_colour_mapping_dict,
    generate_sidebar,
    render_download_without_data_table,
    render_section_header,
    setup_country_filter,
    setup_hourly_filters,
    setup_year_filter,
)
from scripts.plot_functions import (
    create_nice_names_and_color_mapping,
    plot_line_with_secondary_y_hourly,
)

# =========================== Render functions for each chart =========================


def render_t1_ev_load_profile(graph_config: dict) -> None:
    """Render EV load profile with dual y-axes (hourly line chart).

    Parameters
    ----------
    graph_config : dict
        Chart configuration for this plot (from graph_settings.yaml). Expected
        keys include `name`, `table_name`, `leg_col`, and `download_id`.
    """
    # Inject graph_type for legacy compatibility
    graph_config = {**graph_config, "graph_type": "line_with_secondary_y_hourly"}
    render_section_header(graph_config["name"])

    # Setup filters
    is_dual = bool(st.session_state.output_sce2 and st.session_state.output_sce2 != "")
    shared_year = setup_year_filter(graph_config, is_dual)
    graph_config["shared_year"] = str(shared_year)
    graph_config["shared_country"] = setup_country_filter(
        graph_config, is_dual, scenario_tag=st.session_state.output_sce1
    )

    # Extract config values
    table_name = graph_config["table_name"]
    legend_col = graph_config["leg_col"]
    mapping_df = create_nice_names_and_color_mapping(table_name)

    # Load and validate scenario data
    scenario_1_raw = load_and_validate_hourly_data(
        st.session_state.output_sce1,
        table_name,
        str(shared_year),
        graph_config["shared_country"],
    )
    if scenario_1_raw is None:
        st.divider()
        return

    scenario_2_raw = None
    if is_dual:
        scenario_2_raw = load_and_validate_hourly_data(
            st.session_state.output_sce2,
            table_name,
            str(shared_year),
            graph_config["shared_country"],
        )

    # Setup hourly filters
    has_dual_filters = is_dual and scenario_2_raw is not None
    setup_hourly_filters(graph_config, scenario_1_raw, scenario_2_raw, has_dual_filters)

    # Filter and prepare data (without nice names for p14 special handling)
    monthly_1, start_date, end_date, is_complete = get_filtered_df_and_date_range(
        scenario_1_raw, graph_config
    )
    if start_date > end_date:
        st.error("Error: End date must be greater than or equal to start date.")
        st.divider()
        return

    scenario_1_filtered = filter_dataframe_by_date_range(
        monthly_1, start_date=start_date, end_date=end_date
    )
    scenario_1_filtered = handle_small_values(scenario_1_filtered)

    scenario_2_filtered = None
    if has_dual_filters and scenario_2_raw is not None:
        monthly_2, _, _, _ = get_filtered_df_and_date_range(
            scenario_2_raw, graph_config
        )
        scenario_2_filtered = filter_dataframe_by_date_range(
            monthly_2, start_date=start_date, end_date=end_date
        )
        scenario_2_filtered = handle_small_values(scenario_2_filtered)

    # Calculate common y-axis range
    y_range = prepare_y_range(scenario_1_filtered, scenario_2_filtered, None)

    # Create label mapping for primary and secondary y-axes
    label_map = {
        label: (
            mapping_df.loc[label, "nice_names"]
            if (mapping_df is not None and label in mapping_df.index)
            else prettify_label(label)
        )
        for label in graph_config["primary_y_lab"] + graph_config["secondary_y_lab"]
    }
    graph_config["label_map"] = label_map

    plot_kwargs = {
        "start_date": start_date,
        "end_date": end_date,
        "is_complete": is_complete,
    }

    if not is_dual or scenario_2_filtered is None:
        st.caption(f"{st.session_state.output_sce1}")
        colour_mapping = generate_colour_mapping_dict(
            table_name,
            mapping_df,
            add_nice_names(scenario_1_filtered, legend_col, mapping_df),
            legend_col,
        )
        plot_line_with_secondary_y_hourly(
            scenario_1_filtered,
            graph_config,
            colour_mapping,
            y_range,
            key=f"plotly_chart_{st.session_state.output_sce1}_{table_name}",
            **plot_kwargs,
        )
        render_download_without_data_table(
            scenario_1_filtered, graph_config, st.session_state.output_sce1
        )
    else:
        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            st.caption(f"{st.session_state.output_sce1}")
            colour_mapping_1 = generate_colour_mapping_dict(
                table_name,
                mapping_df,
                add_nice_names(scenario_1_filtered, legend_col, mapping_df),
                legend_col,
            )
            plot_line_with_secondary_y_hourly(
                scenario_1_filtered,
                graph_config,
                colour_mapping_1,
                y_range,
                key=f"plotly_chart_{st.session_state.output_sce1}_{table_name}",
                **plot_kwargs,
            )

        with col3:
            st.caption(f"{st.session_state.output_sce2}")
            colour_mapping_2 = generate_colour_mapping_dict(
                table_name,
                mapping_df,
                add_nice_names(scenario_2_filtered, legend_col, mapping_df),
                legend_col,
            )
            plot_line_with_secondary_y_hourly(
                scenario_2_filtered,
                graph_config,
                colour_mapping_2,
                y_range,
                key=f"plotly_chart_{st.session_state.output_sce2}_{table_name}",
                **plot_kwargs,
            )

        col1, _, col3 = st.columns([6, 1, 6])
        with col1:
            render_download_without_data_table(
                scenario_1_filtered, graph_config, st.session_state.output_sce1
            )
        with col3:
            render_download_without_data_table(
                scenario_2_filtered, graph_config, st.session_state.output_sce2
            )

    st.divider()


if __name__ == "__main__":
    st.title(":material/directions_car: Transport")
    output_config = st.session_state.output_config["transport"]
    TRANSPORT_CHART_KEYS = [
        "t1",
    ]
    show_transport = "t" in str(st.session_state.get("sector", "")).lower()
    if show_transport:
        render_t1_ev_load_profile(output_config["t1"])
        # Only include charts in the sidebar if they are present in the config
        transport_charts = [
            output_config[key] for key in TRANSPORT_CHART_KEYS if key in output_config
        ]
        table_of_content = [chart["name"] for chart in transport_charts]
        generate_sidebar(table_of_content)
    else:
        st.warning("Transport sector is not selected. No chart is displayed.")
