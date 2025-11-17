<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Standard output from the Snakemake workflow

The following NetCDF files, CSVs, and charts are generated automatically after the whole Snakemake workflow is successfully executed.

## NetCDF files (`*.nc` files)

[NetCDF (Network Common Data Form)](https://en.wikipedia.org/wiki/NetCDF){:target="_blank"} is a data format for efficiently storing multi-dimensional arrays. By default, the PyPSA-SPICE model creates:

- **Pre-solve networks** in the `results/pre-solve` directory in the `project_name` folder.
- **Pre-solve-brownfield networks** in the `results/pre-solve-brownfield` directory in the `project_name` folder.
- **Post networks** in the `results/post-solve` directory in the `project_name` folder.

These `*.nc` files allow users to examine the model’s structure, inputs, and results before and after optimisation.

## Excel-ready CSVs

The following CSV files are created automatically. Each file represents a key indicator, broken down by sector and modelling year or hour, depending on granularity.

| File Name                                | Sector    | Unit                                               | Description                                                                                |
| ----------------------------------- | --------- | -------------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `ene_avg_fuel_costs_fuel_yearly`      | Energy    | currency/MWh~th~    | Average fuel costs by modelling year                                                       |
| `ene_costs_opex_capex_yearly`         | Energy    | million currency   | Total CAPEX and OPEX in energy sectors by modelling year                                         |
| `ene_dmd_by_carrier_yearly`           | Energy    | TWh                                                | Energy demand by carrier (e.g., `Electricity`) by modelling year                              |
| `ene_emi_by_carrier_by_sector_yearly` | Energy    | MtCO~2~                                              | Total emissions by carrier and sector by modelling year                                  |
| `ene_fom_by_type_yearly`              | Energy    | million currency                                   | Total fixed operation and maintenance cost by technology (e.g., `CCGT`) by modelling year    |
| `ene_gen_by_carrier_yearly`           | Energy    | TWh                                                | Total generation by carrier (e.g., `Electricity`) by modelling year                              |
| `ene_opex_by_type_yearly`             | Energy    | million currency   | OPEX in energy sectors by technology (e.g., `CCGT`) by modelling year                         |
| `ind_cap_by_carrier_by_region_yearly` | Industry  | GW                                                 | Installed capacity in the industry sector by carrier and region/node, and by modelling year |
| `ind_cap_by_type_by_carrier_yearly`   | Industry  | GW                                                 | Installed capacity in the industry sector by technology and carrier, and by modelling year  |
| `ind_emi_by_carrier_yearly`           | Industry  | MtCO~2~                                              | Emissions in the industry sector by carrier by modelling year                              |
| `ind_gen_by_carrier_by_region_yearly` | Industry  | TWh                                                | Generation in the industry sector by carrier and region/node, and by modelling year         |
| `ind_gen_by_type_by_carrier_by_heatgroup_yearly`   | Industry  | TWh                                                | Generation in the industry sector by technology and carrier, heating group, and by modelling year          |
| `ind_gen_by_type_hourly`              | Industry  | MW                                                 | Generation in the industry sector by technology by modelling year                 |
| `ind_hh_gen_by_type_hourly`           | Industry  | MW                                                 | Generation in the industry sector (high-heat) by technology by modelling year                 |
| `ind_lh_gen_by_type_hourly`           | Industry  | MW                                                 | Generation in the industry sector (low-heat)  by technology by modelling year                 |
| `pow_bats_charging_hourly`            | Power     | MW                                                 | Hourly battery charging profile by modelling year                                          |
| `pow_bats_ep_ratio`                   | Power     | -                                                 | Battery Energy-to-Power ratio                                                              |
| `pow_battery_flows_by_region_hourly`  | Power     | MW                                                 | Hourly battery flows between different regions/nodes by modelling year                     |
| `pow_cap_by_region_yearly`            | Power     | GW                                                 | Installed capacity in the power sector by region/node by modelling year                    |
| `pow_cap_by_type_yearly`              | Power     | GW                                                 | Installed capacity in the power sector by technology (e.g., `CCGT`) by modelling year         |
| `pow_cap_by_type_by_region_yearly`    | Power     | GW                                                 | Installed capacity in the power sector by technology (e.g., `CCGT`) by region/node by modelling year         |
| `pow_capex_by_type_yearly`            | Power     | million currency   | CAPEX in the power sector by technology (e.g., `CCGT`) by modelling year                      |
| `pow_elec_load_by_sector_hourly`      | Power     | MW                                                 | Hourly electricity load by sector                                                          |
| `pow_emi_by_carrier_yearly`           | Power     | MtCO~2~                                              | Emissions in the power sector by carrier (e.g., `Gas`) by modelling year                      |
| `pow_flh_by_type_yearly`              | Power     | Hours                                              | Full load hours by technology (e.g., `CCGT`) by modelling year                                |
| `pow_fom_by_type_yearly`              | Power     | million currency                                   | Fixed operation and maintenance cost in the power sector by technology (e.g., `CCGT`) by modelling year    |
| `pow_gen_by_category_share_yearly`    | Power     | %                                                  | Fossil-Renewables share by modelling years                                                 |
| `pow_gen_by_type_hourly`              | Power     | MW                                                 | Hourly generation in the power sector by technology (e.g., `CCGT`)                            |
| `pow_gen_by_type_region_hourly`       | Power     | TWh                                                | Hourly generation in the power sector by region/node by modelling year                            |
| `pow_gen_by_type_yearly`              | Power     | TWh                                                | Generation in the power sector by technology (e.g., `CCGT`) by modelling year                 |
| `pow_hdam_flows_by_region_hourly`     | Power     | MW                                                 | Hourly hydro dam flow profile by modelling year                                            |
| `pow_hphs_flows_by_region_hourly`     | Power     | MW                                                 | Hourly hydro pumped storage profile by modelling year                                      |
| `pow_inter_cf_by_region_yearly`       | Power     | -                                                 | Capacity factor of the interconnectors by region/node by modelling year                   |
| `pow_intercap_by_region_yearly`       | Power     | GW                                                 | Installed capacity of the interconnectors by region/node by modelling year                 |
| `pow_marginal_price_by_region_hourly` | Power     | currency/MWh       | Hourly marginal price by region/node by modelling year                                     |
| `pow_nodal_flow_hourly`               | Power     | MW                                                 | Hourly exchange flow between different regions/nodes by modelling year                       |
| `pow_opex_by_type_yearly`             | Power     | million currency   | OPEX in the power sector by technology (e.g., `CCGT`) by modelling year                      |
| `pow_overnight_inv_by_type_yearly`    | Power     | million currency   | Overnight investment cost in the power sector by technology (e.g., `CCGT`) by modelling year                      |
| `pow_reserve_by_type_hourly`          | Power     | TWh                                                | Reserve by technology (e.g., onshore wind) by modelling year                                |
| `tran_capex_by_type_yearly`                | Transport | million currency   | CAPEX in the transport sector by technology by modelling year                              |
| `tran_charger_capacity_by_region_yearly`   | Transport | GW                 | Capacity of the chargers in the transport sector by region by modelling year                              |
| `tran_charger_capacity_by_type_yearly`     | Transport | GW                 | Capacity of the chargers in the transport sector by technology by modelling year                              |
| `tran_load_charging_all_hourly_by_region`  | Transport | MW                 | Hourly charging/discharging/load status in the transport sector by region/node by modelling year                                     |
| `tran_load_charging_all_hourly_by_type`    | Transport | MW                 | Hourly charging/discharging/load status in the transport sector by technology by modelling year                                     |
| `tran_load_charging_all_hourly`            | Transport | MW                 | Hourly charging/discharging/load status in the transport sector by modelling year                                     |
| `tran_storage_capacity_by_region_yearly`   | Transport | GW                 | Capacity of the storage in the transport sector by region by modelling year                              |
| `tran_storage_capacity_by_type_yearly`     | Transport | GW                 | Capacity of the storage in the transport sector by technology by modelling year                              |

## Jupyter Notebook

We provide a `post-analysis/post_analysis.ipynb` for manual exploration of the pre-solve, pre-solve-brownfield, or post-solve network files.

## Dynamic visualisation

To make it easier to explore and compare model outputs, we provide an open-source library [PyPSA-SPICE-Vis](../visualisation-tool/pypsa-spice-vis.md), an interactive visualisation library that generates dynamic charts using the outputs listed above.
