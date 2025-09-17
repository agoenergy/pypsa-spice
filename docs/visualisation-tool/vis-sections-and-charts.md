<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# List of available sections and charts

The list of sections and available charts can be adjusted by `pypsa-spice-vis/setting/graph_settings.yaml` including maximum and minimum scales of the y-axis, units of the y-axis, etc.

## Power

| Chart name         | Unit | Description |
| ------------------ | ---- | ----------- |
| Capacity by type   | GW   | Installed capacity in the power sector by technology by modelling year |
| Capacity by region | GW   | Installed capacity in the power sector by region by modelling year |
| Generation by type | TWh  | Power Generation by technology by modelling year |
| Share category     | %    | Power Generation share by renewables & fossil fuels by modelling year |
| Transmission capacity between regions | GW | Capacity of transmission lines between different regions by modelling year |
| Hourly generation | MW | Hourly electricity generation by technology by modelling year |
| Regional hourly generation | MW | Hourly electricity generation by region by modelling year |
| Energy demand by carrier | TWh | Energy demand by carrier by modelling year |
| Hourly demand | MW | Hourly electricity demand by sector by modelling year |
| Hourly elec price | currency/MWh | Hourly marginal electricity price by region by modelling year |
| Hourly nodal flow between regions | MW | Hourly nodal exchange flow between regions by modelling year |
| Battery's E/P ratio | N/A | Hourly Energy-to-Power ratio of batteries by modelling year |
| Battery's charging profile | GW | Hourly charing/discharing status of batteries by modelling year |

## Industry

| Chart name           | Unit | Description |
| -------------------- | ---- | ----------- |
| Capacity by carrier  | GW   | Installed capacity in the industry sector by technology by carrier by modelling year |
| Capacity by region   | GW   | Installed capacity in the industry sector by region by modelling year |
| Generation by region | TWh  | Electricity generation used for industrial heat applications by region by modelling year |
| Generation by type   | TWh  | Electricity generation used for industrial heat applications by technology by modelling year |

## Transport

| Chart name           | Unit | Description |
| -------------------- | ---- | ----------- |
| EV load profile      | MW   | Hourly charging load of electric vehicles in the tranport sector by modelling year |

## Emissions

| Chart name              | Unit | Description |
| ----------------------- | ---- | ----------- |
| Emission (power sector) | MtCO<sub>2</sub> | Emission in the power sector by carrier by modelling year |
| Emission (industry sector) | MtCO<sub>2</sub> | Emission in the industry sector by carrier by modelling year |

## Costs

| Chart name           | Unit | Description |
| -------------------- | ---- | ----------- |
| CAPEX by type (power sector) | million currency | Capital Expenditure (CAPEX) in the power sector by technology by modelling year |
| Overnight investment by type (power sector) | million currency | Overnight investment in the power sector by technology by modelling year |
| CAPEX by type (industry sector) | million currency | Capital Expenditure (CAPEX) in the industry sector by technology by modelling year |
| CAPEX by type (transport sector) | million currency | Capital Expenditure (CAPEX) in the transport sector by technology by modelling year |
| CAPEX by type (all energy sectors) | million currency | Capital Expenditure (CAPEX) in all sectors by technology by modelling year |
| CAPEX and OPEX (all energy sectors) | million currency | Capital Expenditure (CAPEX) and Operational Expenditure in all sectors by modelling year |
| Average fuel costs (all energy sectors) | million currency/MWh<sub>thermal</sub> | Average fuel costs (thermal units) in all sectors by modelling year |

## Info

This section provides a list of all abbreviations and corresponding full names of all technologies and carriers that are used in the visualisation tool.