<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input Data: Global CSV Template

Global csvs contain parameters that are typically kept constant accross scenarios. This is to maintain comparability of the scenarios. `global_input_template` folder is used for creating skeletons purposes, and thus it can be considered as hidden folder inside the template.

```text title="Structure of the global CSV template files"
📦 data
 ┗ 📂 global_input_template
 ┗ 📂 pypsa-spice-data
    ┗ 📂 project_01
       ┗ 📂 input
         ┗ 📂 global_input
           ┣ 📜 availability.csv
           ┣ 📜 demand_profile.csv
           ┣ 📜 ev_parameters.csv
           ┣ 📜 power_plant_costs.csv
           ┣ 📜 renewables_technical_potential.csv
           ┣ 📜 storage_costs.csv
           ┣ 📜 storage_inflows.csv
           ┗ 📜 technologies.csv
```

!!! Tip
    The currency of all example data is `USD` defined in the `base_configs` section of `base_config.yaml`. You can refer to [Model Builder Configuration](model-builder-configuration.md#base_configyaml) for more information.

## Availability

`availability.csv` contains time-series availability data, mainly for renewable plants. By default, availability is matched using renewables type (e.g., solar photovoltaic (`PHOT`), onshore wind (`WTON`), etc.) and their locations (e.g., region `XY_NO` in country `XY`, region `YZ_SO` in country `YZ`).

If a technology shares the same profile across the country (e.g., electric vehicle charger public (`EVCH-PUB`)), then both region and country fields use the same name (e.g., region `XY` in country `XY`). If the technology is listed and it requires availability profile, but the profile is not in this csv, then it will be defined as **constant 1** for all hours.

## Demand Profile

`demand_profile.csv` stores normalized hourly load profiles, which are scaled using total annual load values of each year to create time-series demand data.

By default, Load profiles are matched based on:

- **Profile type and location** for power sector loads such as wholesale market load (`HV_LOAD`) and building load (`LV_LOAD`).
- **Profile type only** for all other loads.

To add new load profiles (e.g., for a new project or country), insert a new row.

## EV Parameters

`ev_parameters.csv` stores the technical parameters relevant to the electric vehicles.

## Power Plant Costs

`power_plant_costs.csv` defines cost data for all technologies in each country. It includes:

- Capital expenditure (CAPEX) in USD/MW (currency based on input data)
- Fixed operation and maintenance cost (FOM) in USD/MW (currency based on input data)
- Variable operation and maintenance cost (VOM) in USD/MWh (currency based on input data)

Note: Currencies may vary depending on the source data.

This data applies to generators, storage, converters, and storage capacity expansion (In case of lithium battery, it refers to inverter costs).

## Renewables Technical Potential

`renewables_technical_potential.csv` defines maximum expansion limits (technical potential or land-use limits) for renewable technologies. It is currently only applied to solar photovoltaic (`PHOT`), hydro run-of-river (`HROR`), onshore wind (`WTON`), offshore wind (`WTOF`), rooftop PV (`RTPV`), solar hot water heater (`SWHT`) but can be modified to apply for other technologies.

The model builder does not allow higher expansion than what are specified in this CSV file. You can expand the file to include other technologies if needed.

## Storage Costs

`storage_costs.csv` covers the cost structure for all storage tanks (storage volume or energy capacity) in each country. It includes:

- Capital expenditure (CAPEX) in USD/MW (currency based on input data)
- Fixed operation and maintenance cost (FOM) in USD/MW (currency based on input data)
- Variable operation and maintenance cost (VOM) in USD/MWh (currency based on input data)

Note: Currencies may vary depending on the source data.

## Storage Inflows

`storage_inflows.csv` provides time-series inflow data [MW] for `StorageUnit` components. Inflow profiles are only designed for reservoir-based systems like hydropower and hydro pumped storage.

Matching the inflows is based on the technology (hydro dam (`HDAM`) and/or hydro pumped storage (`HPHS`)) and their location (e.g., region `XY_NO` in country `XY`). If the technology is listed and it requires inflow profile, but the profile is not in this csv, then it will be defined as **constant 0** for all hours.

## Technologies

technologies.csv defines typical technical parameters for each technology used in the model builder. It provides values like efficiency, ramp limits, and availability of various technologies. To add a new technology, insert a new row and fill in all required parameters.

Description of all technical parameters:

| Parameter                              | definition                                                                                                                 |
| -------------------------------------- | -------------------------------------------------------------------------------------------------------------------------- |
| `country`                                | 2-letter country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"} format                                                                                                |
| `technology`                             | Abbreviations of the technology                                                                                          |
| `technology_nomenclature`                | Full names of the technology                                                                                             |
| `carrier`                                | Resources used by the technologies                                                                                        |
| `class`                                  | Component class as defined in PyPSA                                                                                   |
| `efficiency`                | Energy conversion efficiency from primary energy to electricity for `Generators`, and to another form of energy for `Links`. For `StorageUnits`, this is the discharge efficiency.                                                                |
| `efficiency2` | **Positive** values represent emission factor and **negative** values correspond to the efficiency of generating the second product in a plant                                                                                  |
| `efficiency3`                    | Carbon capture efficiency (for CCS technologies)                                                                      |
| `efficiency_store`        | Efficiency of charging energy into storage                                                               |
| `max_hours`               | Maximum charge duration in hours (total storage volume / capacity)                                         |
| `cyclic_state_of_charge`  | If **True**, the final state of charge equals the initial state of charge                                 |
| `state_of_charge_initial` | Initial state of charge in MWh before the snapshots in the optimal Power Flow (MWh)                                                       |
| `p_max_pu`                               | The maximum availability per snapshot per unit of `p_nom`                                                                            |
| `p_min_pu`                               | The minimum availability per snapshot per unit of `p_nom`                                                                           |
| `ramp_limit_down`                        | Maximum active power decrease from one snapshot to the next (per unit)                                                     |
| `ramp_limit_up`                          | Maximum active power increase from one snapshot to the next (per unit)                                                     |
| `standing_loss`           | Hourly energy loss from storage                                                                              |
| `r_rating`                               | Contribution of reserve rating (if used)                                                                            |
