# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Render the unified Output subpage."""

import streamlit as st

from scripts.output_pages import (output_costs, output_emissions,
                                  output_industry, output_power,
                                  output_transport)


def main():
    """Render all Output views in one subpage."""
    st.title(":material/monitoring: Output")

    output_tabs = ["Power", "Emissions", "Costs"]
    include_industry = "i" in st.session_state.sector
    include_transport = "t" in st.session_state.sector

    if include_industry:
        output_tabs.insert(1, "Industry")
    if include_transport:
        output_tabs.insert(2 if include_industry else 1, "Transport")

    tabs = st.tabs(output_tabs)

    tab_idx = 0
    with tabs[tab_idx]:
        output_power.main()
    tab_idx += 1

    if include_industry:
        with tabs[tab_idx]:
            output_industry.main()
        tab_idx += 1

    if include_transport:
        with tabs[tab_idx]:
            output_transport.main()
        tab_idx += 1

    with tabs[tab_idx]:
        output_emissions.main()
    tab_idx += 1

    with tabs[tab_idx]:
        output_costs.main()


if __name__ == "__main__":
    main()
