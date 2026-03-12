# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Render the unified Input subpage."""

import streamlit as st

from scripts.getters import Getters
from scripts.input_pages import input_power_demand, input_power_supply


def main():
    """Render all Input views in one subpage."""
    st.title(":material/input: Input")

    getters = Getters()
    supply_tab, demand_tab = st.tabs(["Power - Supply", "Power - Demand"])

    with supply_tab:
        input_power_supply.main(getters)

    with demand_tab:
        input_power_demand.main(getters)


if __name__ == "__main__":
    main()
