<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Power Sector

## Key Features

- **Co-optimisation** of generation and capacity expansion including interconnections.
- **Myopic** (Year-by-year) optimisation. Each year is optimised independently, without assuming knowledge of future developments.
- **Brownfield** modelling approach. The model builds on existing infrastructure, meaning capacity from previous years is retained and carried forward.

PyPSA-SPICE follows the component definitions from [PyPSA Components](https://docs.pypsa.org/latest/user-guide/design/){:target="_blank"}. The diagram below illustrates all components involved in energy flows at a single node in the power sector.

[![PyPSA-SPICE power sector energy flow](../assets/images/pypsa-spice_schema_power_sector.svg){ .img-center width="100%" }](../assets/images/pypsa-spice_schema_power_sector.svg){: target="_blank" }

## Power Generators

All the listed components are defined as `Generator` in PyPSA.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `CSP`         | Concentrated solar power plant           |
| `GEOT`        | Geothermal power plant                   |
| `HROR`        | Hydro run-of-river                       |
| `PHOT`        | Solar PV                                 |
| `RTPV`        | Rooftop PV                               |
| `WTOF`        | Offshore wind                            |
| `WTON`        | Onshore wind                             |

## Power Links

All the listed components are defined as `Link` in PyPSA.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `BIOT`        | Biomass power plant                      |
| `CCGT`        | Combined-cycle gas turbine power plant   |
| `CHP`         | Combined heat and power plant            |
| `ELTZ`        | Electrolyser (for hydrogen production)   |
| `NUCL`        | Nuclear power plant                      |
| `OCGT`        | Open-cycle gas turbine power plant       |
| `OCHT`        | Open-cycle hydrogen turbine power plant  |
| `OILT`        | Oil turbine power plant power plant                  |
| `SubC`        | Subcritical coal-fired power plant       |
| `SupC`        | Supercritical coal-fired power plant     |
| `WSTT`        | Waste-to-Energy power plant              |

## Storage Capacity

The following component is defined as `StorageUnit` in PyPSA.

Storages can be modelled with two approaches.

1. **Fixed Energy/Power Ratio:** In this case, the energy to power ratio for storage is predefined. You can use multiple storage type with different energy to power ratio. For example, `BATS` with E/P ratio of 4 and `BATS` with E/P ratio of 8 representing different energy to power ratio and the model will optimise the capacity of each of these technology. In this case, PyPSA type `Storage_units` can be used for modelling and defining the energy to power ratio in `technologies.csv`.
2. **Variable Energy/Power Ratio:** If you want the model to optimise the energy/power ratio of storage your have to model it using a combination of `links` + `store` component. This requires separate inputs like costs for capacity and energy component of the storage inputs.

!!! Tip
    In PyPSA components, `StorageUnit` is modelled as a storage asset with a fixed energy-to-power ratio defined by `max_hours` of the nominal power (you can also refer to [PyPSA Components - Storage Unit](https://docs.pypsa.org/latest/user-guide/components/storage-units/){:target="_blank"} for more information). Thus, in PyPSA-SPICE model builder, hydro dam `HDAM` is defined as a `StorageUnit` and it is given in storage capcaity only to represent nominal power-related params. <br><br>
    To model the storage energy separately from the power capacity, `store` + 2 `links` is a better combination. You can refer to [Storage Energy](power_sector.md#storage-energy) for more information. Technologies defined in the storage energy require storage capacity if the carrier is related to electricity (power).

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `HDAM`        | Hydro dam                                |
| `BATS`        | Utility-scale battery storage            |
| `HHBS`        | Household battery storage                |
| `HPHS`        | Hydro pumped storage                     |

## Storage Energy

All the listed components are defined as `Store` in PyPSA.

!!! Tip
    In PyPSA components, `Store` is modelled as a storage asset with only energy storage. It can optimise energy capacity separately from the power capacity with a combination of `store` + 2 `links`. The links represent charging and discharging characteristics to control the power output. Marginal cost and efficiency of charging and discharging can be defined in each link.<br><br>
    In PyPSA-SPICE model builder, technologies that are defined as storage energy, **they should also be included in [Storage Capacity](power_sector.md#storage-capacity) to describe charging and discharging processes. The links are created automatically , and hence it's not required to add charging and discharging links inside [Power Links](power_sector.md#power-links).**<br><br>
    Detailed information and example can be found in [PyPSA Components - Store](https://docs.pypsa.org/latest/user-guide/components/stores/){:target="_blank"} and [Replace StorageUnits with fundamental Links and Stores](https://docs.pypsa.org/latest/examples/replace-generator-storage-units-with-store/){:target="_blank"}.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `CO2STOR`     | CO~2~ storage                            |
| `BATS`        | Utility-scale battery storage            |
| `HHBS`        | Household battery storage                |
| `HPHS`        | Hydro pumped storage                     |

## Carriers

| Abbreviation  | Full Name                                                     |
| ------------- | ------------------------------------------------------------- |
| `Bio`         | Biomass                                                       |
| `Bit`         | Bituminous or brown coal                                      |
| `CO2`         | Carbon dioxide (in the atmosphere)                            |
| `Co2stor`     | Captured carbon dioxide                                       |
| `Electricity` | Electricity                                                   |
| `Gas`         | Domestic natural gas                                          |
| `Gas-imp`     | Imported natural gas                                          |
| `High_Heat`   | High-temperature heat (> 350°C)                               |
| `Hrdc`        | Anthracite or hard coal                                       |
| `Hyd`         | Hydrogen                                                      |
| `Lig`         | Lignite                                                       |
| `Lng`         | Liquefied natural gas                                         |
| `Low_Heat`    | Low/Medium-temperature heat (< 350°C)                         |
| `Oil`         | Oil                                                           |
| `Uranium`     | Uranium                                                       |
| `Waste`       | Waste                                                         |

## Buses

| Abbreviation  | Full Name                     |
| ------------- | ----------------------------- |
| `ATMP`        | Atmosphere                    |
| `BATSN`       | Lithium battery storage       |
| `BION`        | Biomass                       |
| `BITN`        | Bituminous                    |
| `CO2STORN`    | CO~2~ storage                 |
| `GASN`        | Gas                           |
| `HHBSN`       | Household battery storage     |
| `HPHSN`       | Hydro pumped storage          |
| `HRDCN`       | Anthracite or hard coal       |
| `HVELEC`      | High-voltage electricity      |
| `HYDN`        | Hydrogen                      |
| `LIGN`        | Lignite                       |
| `LNGN`        | liquefied natural gas         |
| `LVELEC`      | Low-voltage electricity       |
| `NUCLN`       | Uranium                       |
| `OILN`        | Oil                           |
| `WSTN`        | Waste                         |

## Other Components

| Abbreviation | Full Name                                               |
| ------------ | ------------------------------------------------------- |
| `co2Price`   | Price of emitting one unit of CO~2~ into the atmosphere |
| `r`          | Interest rate                                           |
| `HV_LOAD`   | Wholesale market load (High voltage level)              |
| `LV_LOAD`   | Building load (low/medium voltage level)                |

## Custom Constraints (Defined in the `config.yaml` File)

- CO</sub>2</sub> management
- Energy independence
- Fuel production constraint
- Reserve margin
- Renewable generation share constraint
- Must run constraint of thermal generators
- Capacity factor constraint

You can refer to [Model Builder Constraints](../getting-started/input-data/model-builder-configuration.md#custom-constraints) for more information.
