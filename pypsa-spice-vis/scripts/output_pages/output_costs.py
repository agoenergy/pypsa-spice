# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Create Costs page under Results section.

Page shows costs related dataframes and visualisations from the modelling results.
"""

import os

import streamlit as st
import yaml

from scripts.output_st_handler import (generate_sidebar,
                                       map_chart_to_plot_function,
                                       render_st_page_and_plot)


def main():
    """Render the Costs output page."""
    st.title(":material/attach_money: Costs")

    with open(
        os.path.join(st.session_state.current_dir, "setting/graph_settings.yaml"),
        encoding="utf-8",
    ) as file:
        config = yaml.safe_load(file)["costs"]

    table_of_content = []

    for _item, values in config.items():
        if (
            values["incl_sector"] == "all"
            or values["incl_sector"] in st.session_state.sector
        ):
            render_st_page_and_plot(
                graph_type_func=map_chart_to_plot_function(values["graph_type"]),
                config_plot=values,
            )
            table_of_content.append(values["name"])
        else:
            pass

    generate_sidebar(table_of_content)


if __name__ == "__main__":
    main()
