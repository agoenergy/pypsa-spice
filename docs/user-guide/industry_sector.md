<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Industry sector

## Key features

Optimisation of industry heat supply at two different temperature levels:

- High-temperature (above 350°C)
- Low-/medium-temperature (350°C or below)

The structure and functionality of components follow the [PyPSA Components](https://docs.pypsa.org/latest/user-guide/design/){:target="_blank"}. The diagram below shows the full set of components involved in energy flows for a single industrial node.

[![PyPSA-SPICE industry sector energy flow](../assets/images/pypsa-spice_schema_industrial_sector.svg){ .img-center width="100%" }](../assets/images/pypsa-spice_schema_industrial_sector.svg){: target="_blank" }

## Heat generators

The following component is defined as `Generator` in PyPSA.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `SWHT`        | Solar hot water heater                   |

## Heat links

All the listed components are defined as `Link` in PyPSA.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `EERH`        | Electric resistance heater               |
| `EDLH`        | Dielectric heating technology            |
| `EHPP`        | Industry heat pump                       |
| `EIDT`        | Induction heat boiler                    |
| `FITR`        | Fischer-Tropsch process                  |
| `IND_BOILER`  | Industry heat boiler                     |
| `INLHSTOR`    | Low-temperature heat storage             |
| `METH`        | Methanation                              |

## Fuel conversion

All the listed components are defined as `Link` in PyPSA.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `ELTZ`        | Electrolyser (for hydrogen production)   |
| `FITR`        | Fischer-Tropsch process                  |

## Direct air capture

The following components is defined as `Link` in PyPSA.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `DAC`         | Direct air capture                       |

## Storage capacity

The following component is defined as `StorageUnit` in PyPSA.

!!! Tip
    In PyPSA components, `StorageUnit` is modelled as a storage asset with a fixed energy-to-power ratio defined by `max_hours` of the nominal power (you can also refer to [PyPSA Components - StorageUnit](https://docs.pypsa.org/latest/user-guide/components/storage-units/){:target="_blank"} for more information).<br><br>
    To model the storage energy separately from the power capacity, `Store` + 2 `Links` is a better combination. You can refer to [Storage energy](power_sector.md#storage-energy) for more information. Technologies defined in the storage energy require storage capacity if the carrier is related to electricity (power).

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `INLHSTOR`    | Low-temperature heat storage             |

## Storage energy

All the listed components are defined as `Store` in PyPSA.

!!! Tip
    In PyPSA components, `Store` is modelled as a storage asset with only energy storage. It can optimise energy capacity separately from the power capacity with a combination of `Store` + 2 `Links`. The links represent charging and discharging characteristics to control the power output. Marginal cost and efficiency of charging and discharging can be defined in each link.<br><br>
    IIn the industry sector of PyPSA-SPICE model builder, the media `electricity` is replaced by `low-temperature heat` in the storage process. Technologies that are defined as storage energy, **they should also be included in [Storage capacity](power_sector.md#storage-capacity) to describe charging and discharging processes. The links are created automatically, and hence it's not required to add charging and discharging links inside [Heat links](industry_sector.md#heat-links), [Fuel conversion](industry_sector.md#fuel-conversion), or [Direct air capture](industry_sector.md#direct-air-capture).**<br><br>
    Detailed information and example can be found in [PyPSA Components - Store](https://docs.pypsa.org/latest/user-guide/components/stores/){:target="_blank"} and [Replace StorageUnits with fundamental Links and Stores](https://docs.pypsa.org/latest/examples/replace-generator-storage-units-with-store/){:target="_blank"}.

| Abbreviation  | Full Name                                |
| ------------- | ---------------------------------------- |
| `INLHSTOR`    | Low-temperature heat storage             |

## Carriers

| Abbreviation  | Full Name                                                     |
| ------------- | ------------------------------------------------------------- |
| `Bio`         | Biomass                                                       |
| `CO2`         | Carbon dioxide (in the atmosphere)                            |
| `Co2stor`     | Captured carbon dioxide                                       |
| `Electricity` | Electricity                                                   |
| `Gas`         | Domestic natural gas                                          |
| `High_Heat`   | High-temperature heat (> 350°C)                               |
| `Hrdc`        | Anthracite or hard coal                                       |
| `Hyd`         | Hydrogen                                                      |
| `Low_Heat`    | Low-/Medium-temperature heat (< 350°C)                        |
| `Oil`         | Oil                                                           |

## Buses

| Abbreviation  | Full Name                               |
| ------------- | --------------------------------------- |
| `CO2STORN`    | CO~2~ storage                           |
| `HVELEC`      | High-voltage electricity                |
| `HYDN`        | Hydrogen                                |
| `LVELEC`      | Low-voltage electricity                 |
| `IND_LH`      | Industrial low-temperature heat         |
| `IND_HH`      | Industrial high-temperature heat        |
| `INLHSTORN`   | Industrial low-temperature heat storage |

## Other components

| Abbreviation | Full Name                                               |
| ------------ | ------------------------------------------------------- |
| `co2Price`   | Price of emitting one unit of CO~2~ into the atmosphere |
| `IND_LOAD`   | Industrial load (both high- and low-temperature heat)   |

## Custom constraints (defined in the `config.yaml` file)

- coming soon...
