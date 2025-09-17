<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Transport Sector

## Key Features

- Optimal charging of electric vehicles (EVs), taking into account availability and charging constraints.
- Supply of fuels for transport.

The structure and function of each component follow the definitions from the [PyPSA components](https://pypsa.readthedocs.io/en/latest/user-guide/components.html){:target="_blank"}. The diagram below illustrates the components and energy flows at a single node in the transport sector.

[![PyPSA-SPICE transport sector energy flow](../assets/images/pypsa-spice_schema_transport_sector.svg){ .img-center width="60%" }](../assets/images/pypsa-spice_schema_transport_sector.svg){: target="_blank" }

!!! Tip
    In the transport sector, PyPSA-SPICE separated the demand and energy flow into two categories: private and public sector. The public sector (`PUB`) encompasses public transportation services, while the private sector (`PRV`) refers to personally owned vehicles.

## Electric Vehicle Chargers

All the listed components are defined as `Link` in PyPSA.

| Abbreviation  | Full Name                          |
| ------------- | ---------------------------------- |
| `EVCH-PRV`    | Electric vehicle charger (Private) |
| `EVCH-PUB`    | Electric vehicle charger (Public)  |

## Electric Vehicle Storage

All the listed components are defined as `Store` in PyPSA.

!!! Tip
    In PyPSA components, `Store` is modelled as a storage asset with only energy storage. It can optimise energy capacity separately from the power capacity with a combination of `store` + 2 `links`. The links represent charging and discharging characteristics to control the power output. Marginal cost and efficiency of charging and discharging can be defined in each link.<br><br>
    In the transport sector of PyPSA-SPICE model builder, technologies that are defined as storage energy, their links of charging and discharging links are defined in [Electric Vehicle Chargers](transport_sector.md#electric-vehicle-chargers).<br><br>
    Detailed information and example can be found in [PyPSA Components - Store](https://pypsa.readthedocs.io/en/latest/user-guide/components.html#store){:target="_blank"} and [Replace StorageUnits with fundamental Links and Stores](https://pypsa.readthedocs.io/en/latest/examples/replace-generator-storage-units-with-store.html){:target="_blank"}.

| Abbreviation  | Full Name                          |
| ------------- | ---------------------------------- |
| `EVST-PRV`    | Electric vehicle storage (Private) |
| `EVST-PUB`    | Electric vehicle storage (Public)  |

## Carriers

| Abbreviation  | Full Name    |
| ------------- | ------------ |
| `Electricity` | Electricity  |

## Buses

| Abbreviation  | Full Name                     |
| ------------- | ----------------------------- |
| `HVELEC`      | High-voltage electricity      |
| `LVELEC`      | Low-voltage electricity       |
| `TRAN-PUB`    | Public electric vehicle       |
| `TRAN-PRV`    | Private electric vehicle      |

## Other Components

| Abbreviation | Full Name                                   |
| ------------ | ------------------------------------------- |
| `HPV_LOAD`   | Transport load (High voltage level)         |
| `LPV_LOAD`   | Transport load (low/medium voltage level)   |

## Custom Constraints (Defined in the `config.yaml` File)

- Coming soon...
