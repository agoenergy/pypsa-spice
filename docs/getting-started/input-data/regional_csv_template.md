<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input Data: Regional CSV Template

```text title="Structure of the regional CSV template files"
📦 data
 ┗ 📂 example
    ┗ 📂 project_01
     ┗ 📂 scenario_01
       ┣ 📂 Power
       ┃ ┣ 📜 buses.csv
       ┃ ┣ 📜 decommission_capacity.csv
       ┃ ┣ 📜 fuel_suppliers.csv
       ┃ ┣ 📜 interconnector.csv
       ┃ ┣ 📜 loads.csv
       ┃ ┣ 📜 power_generators.csv
       ┃ ┣ 📜 power_links.csv
       ┃ ┣ 📜 storage_capacity.csv
       ┃ ┗ 📜 storage_energy.csv
       ┣ 📂 Industry
       ┃ ┣ 📜 buses.csv
       ┃ ┣ 📜 decommission_capacity.csv
       ┃ ┣ 📜 direct_air_capture.csv
       ┃ ┣ 📜 fuel_conversion.csv
       ┃ ┣ 📜 heat_generators.csv
       ┃ ┣ 📜 heat_links.csv
       ┃ ┣ 📜 loads.csv
       ┃ ┣ 📜 storage_capacity.csv
       ┃ ┗ 📜 storage_energy.csv
       ┗ 📂 Transport
         ┣ 📜 buses.csv
         ┣ 📜 loads.csv
         ┣ 📜 pev_chargers.csv
         ┗ 📜 pev_storages.csv
```

!!! Tip
    The currency of all example data is `USD` defined in the `base_configs` section of `config.yaml`. You can refer to [Model Builder Configuration](model-builder-configuration.md#init-settings) for more information.

!!! Tip
    If there's a cell with `inf` in the csv files, it represents infinite value in `float` datatype when it is read into the network.

## Buses

This file defines the buses to be used in the model. All components need to be connected to one or more buses.

| Parameter  | definition                                                                 |
| ---------- | -------------------------------------------------------------------------- |
| `country`    | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `node`       | Node or region name within the country (can be the same as country if the model is not region-specific)  |
| `bus`        | PyPSA `bus` component. Format: `{NODE}_{TECHNOLOGY}N` or `{COUNTRY}_{TECHNOLOGY}N` (suffix `N` is not applied for technologies with `HVELEC`, `LVELEC`, `IND-LH`, `IND-HH`, or `ATMP`) |
| `carrier`    | Fuel or resource name. Unlike other entries that are in uppercase, here only the first letter is capitalized.                  |

## Decommission Capacity

`Decommission_capacity.csv` contains the installed capacity of power plants scheduled for decommissioning.

- For **Generators**, the `name` column must match the `name` column in `input_dir/project_name/input_scenario_name/Power/power_generators.csv`.
- For **Links**, the `name` column must match the `link` column in `input_dir/project_name/input_scenario_name/Power/power_links.csv`


| Parameter   | definition                                                |
| ----------- | --------------------------------------------------------- |
| `country`     | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format       |
| `name`        | Asset name (to be decommissioned)  |
| `class`       | Component type of the asset in PyPSA network. Only the first letter is capitalized |
| `years`| Decommission plan in each year [MW]             |

## Fuel Supplies

These are fuel supply generators that provide fuel in the thermal energy unit [MWh_th]. It is possible to put maximum supply constraint over a year for these fuel supplies. See model [schematic diagram](../../user-guide/model-builder-methodology.md/#model-builder-methodology)

| Parameter    | definition                                                                     |
| ------------ | ------------------------------------------------------------------------------ |
| `country`      | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                            |
| `bus`          | Fuel supply `bus`. Format: `{COUNTRY}_{TECHNOLOGY}N` (suffix `N` is not applied for technologies with `HVELEC`, `LVELEC`, `IND-LH`, `IND-HH`, or `ATMP`) |
| `supply_plant` | Fuel supply hub. Format: `TGEN_{TECHNOLOGY}N` |
| `carrier`   | Fuel or resource name                           |
| `max_supply`   | Annual fuel supply limit [MWh/year]                          |
| `fuel_cost`    | Fuel cost [CURRENCY/MWh]                                        |
| `year`         | Year of the fuel supply                                                        |

## Interconnector

Interconnectors connect different regions by their maximum power transfer capacity.

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                  | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `link`                     | Name of the interconnection link between two regions/countries. Format: `{NODE}_HVELEC_to_{NODE}_HVELEC` |
| `bus0`                     | Region/country exporting electricity to `bus1`. Used as the `bus` component in the PyPSA network Format: `{NODE}_HVELEC`       |
| `bus1`                     | Region/country importing electricity from `bus0`. Used as the `bus` component in the PyPSA network Format: `{NODE}_HVELEC`       |
| `carrier`                  | Energy carrier or resource (e.g., electricity, gas). Only the first letter is capitalized   |
| `type`                     | Interconnector technology (e.g., `ITCN`). All uppercase |
| `efficiency`               | Efficiency of the interconnector                                          |
| `p_max_pu`                 | Maximum availability per snapshot (per unit of `p_nom`)                           |
| `p_min_pu`                 | Minimum availability per snapshot (per unit of `p_nom`)                           |
| `p_nom`                    | Nominal capacity in the default year [MW]                        |
| `p_nom_extendable`         | Indicates if capacity can be expanded. Possible values: `TRUE` or `FALSE`    |
| `CAP`                      | Capital expenditure in USD/MW (currency based on input data)     |
| `FOM`                      | Fixed annual operation and maintenance cost in USD/MWa (currency based on input data)     |
| `marginal_cost`            | Marginal cost of the link in USD/MWh (currency based on input data)      |
| `p_nom_max_{YEAR}` | Maximum additional capacity allowed for the given year in MW  |
| `p_nom_min_{YEAR}` | Minimum additional capacity allowed for the given year in MW  |

## Loads

This file contains the **total** load per load type which is matched to `profile_type`. You can connect multiple load types to the same bus and they would be added to create total load for the bus.  

| Parameter    | definition                                                                 |
| ------------ | -------------------------------------------------------------------------- |
| `country`      | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `node`         | Name of the nodes or regions                 |
| `bus`          | PyPSA `bus` component. The values can be `{NODE}_HVELEC` or `{NODE}_LVELEC` for power sector, `{NODE}_IND-LH` or `{NODE}_IND-HH` for industry sector, and `{NODE}_TRAN-PRV` or `{NODE}_TRAN-PUB` for transport sector |
| `profile_type` | Load profile type                         |
| `name`         | Load name. Format: `{BUS}_{PROFILE_TYPE}`   |
| `total_load`   | Total annual load in MW                                      |
| `carrier`      | Energy carrier or resource. First letter capitalized                  |
| `year`         | Year of the load data                                                          |

## Generators

See details of implementation in [power](../../user-guide/power_sector.md/#power-generators) and [industry](../../user-guide/industry_sector.md/#heat-generators) sectors.

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                  | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `node`                     | Name of the nodes or regions                |
| `type`                     | Generator technology.                                    |
| `carrier`                  | Energy carrier or resource. First letter capitalized.                  |
| `bus`                      | PyPSA `bus` component. Format: `{NODE}_{TECHNOLOGY}N` (suffix `N` is not applied for technologies with `HVELEC`, `LVELEC`, `IND-LH`, `IND-HH`, or `ATMP`) |
| `name`                     | Generator name. Format: `{BUS}_{TECHNOLOGY}`               |
| `p_nom`                    | Nominal capacity in the default year in MW                        |
| `p_nom_extendable`         | Indicates if capacity can be expanded. Possible values: `TRUE` or `FALSE`     |
| `p_nom_max_{YEAR}` | Maximum additional capacity allowed in the given year in MW  |
| `p_nom_min_{YEAR}` | Minimum additional capacity allowed in the given year in MW  |

## Links

See details of implementation logic in [power](../../user-guide/power_sector.md/#power-links), [industry](../../user-guide/industry_sector.md/#heat-links), [transport](../../user-guide/transport_sector.md/#electric-vehicle-chargers) sectors.

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                  | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `bus0...3`                 | PyPSA `bus` components. Format: `{NODE}_{TECHNOLOGY}N` (suffix `N` is not applied for technologies with `HVELEC`, `LVELEC`, `IND-LH`, `IND-HH`, or `ATMP`) |
| `type`                     | Link technology                                    |
| `link`                     | Link name. Format: `{BUS0}_to_{BUS1}_by_{TECHNOLOGY}`        |
| `carrier`                  | Energy carrier or resource. First letter capitalized                  |
| `p_nom`                    | Nominal capacity in the default year in MW                        |
| `p_nom_extendable`         | Indicates if capacity can be expanded. Possible values: `TRUE` or `FALSE`     |
| `p_nom_max_{YEAR}` | Maximum additional capacity allowed in the given year in MW  |
| `p_nom_min_{YEAR}` | Minimum additional capacity allowed in the given year in MW  |

## Storage Capacity

See details of implementation logic for [power](../../user-guide/power_sector.md/#storage-capacity) and [industry](../../user-guide/industry_sector.md/#storage-capacity) sectors.

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                  | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `node`                     | Name of the nodes or regions                |
| `type`                     | Storage technology                                   |
| `carrier`                  | Energy carrier or resource. First letter capitalized                  |
| `bus`                      | PyPSA `bus` component. The values can be `{NODE}_HVELEC` or `{NODE}_LVELEC` for power sector, and `{NODE}_IND-LH` for industry sector |
| `name`                     | Storage name. Format: `{BUS}_{TECHNOLOGY}`              |
| `p_nom`                    | Nominal capacity in the default year in MW                        |
| `p_nom_extendable`         | Indicates if capacity can be expanded. Possible values: `TRUE` or `FALSE`      |
| `p_nom_max_{YEAR}` | Maximum additional capacity allowed in the given year in MW  |
| `p_nom_min_{YEAR}` | Minimum additional capacity allowed in the given year in MW  |

## Storage Energy

See details of description for use of storage energy in [power sector](../../user-guide/power_sector.md/#storage-energy) and [industry sectors](../../user-guide/industry_sector.md/#storage-energy). 

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                  | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `bus`                      | PyPSA `bus` component. Format: `{NODE}_{TECHNOLOGY}N` (suffix `N` is not applied for technologies with `HVELEC`, `LVELEC`, `IND-LH`, `IND-HH`, or `ATMP`) |
| `store`                    | Energy storage name. Format: `{BUS}_STOR`   |
| `type`                     | Storage technology                                    |
| `carrier`                  | Energy carrier or resource. First letter capitalized                  |
| `standing_loss`            | Hourly energy loss rate during storage, in %/hour.                       |
| `e_nom`                    | Nominal energy in the default year in MWh                         |
| `e_nom_extendable`         | Indicates if capacity can be expanded. Possible values: `TRUE` or `FALSE`      |
| `max_store_{YEAR}` | Maximum additional capacity allowed in the given year in MW      |

## EV Chargers

See details of implementation [here](../../user-guide/transport_sector.md/#electric-vehicle-chargers).

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                     |
| `link`                   | Link name. Format: `{BUS0}_to_{BUS1}_by_{TECHNOLOGY}`        |
| `bus0`                     | Region/country exporting electricity to `bus1`. Used as the `bus` component in the PyPSA network Format: `{NODE}_LVELEC`       |
| `bus1`                     | Region/country importing electricity from `bus0`. Used as the `bus` component in the PyPSA network Format: `{NODE}_TRAN-PRV` or `{NODE}_TRAN-PUB`       |
| `type`                     | Storage technology: private (`EVCH-PRV`) or public (`EVCH-PUB`)                    |
| `carrier`                  | Fixed as `Electricity`     |
| `p_max_pu`                 | Profile reading from `availability.csv` based on private (`EVCH-PRV`) or public (`EVCH-PUB`)   |
| `num_ch_{YEAR}`    | Number of electric vehicles charged in the given year      |

## EV Storages

See details of implementation [here](../../user-guide/transport_sector.md/#electric-vehicle-storage).

| Parameter                | definition                                                                 |
| ------------------------ | -------------------------------------------------------------------------- |
| `country`                  | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                        |
| `node`                     | Name of the nodes, regions, or countries |
| `type`                     | Storage technology: private (`EVST-PRV`) or public (`EVST-PUB`)                    |
| `carrier`                  | Fixed as `Electricity`     |
| `bus`                      | PyPSA `bus` component. Format: `{NODE}_TRAN-PRV` or `{NODE}_TRAN-PUB`   |
| `name`                     | Storage name. Format: `{BUS}_{TYPE}`   |
| `num_ev_{YEAR}`    | Number of electric vehicles in the given year      |
