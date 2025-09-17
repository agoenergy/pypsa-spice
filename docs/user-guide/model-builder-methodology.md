<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Model Builder Methodology

PyPSA-SPICE is a least-cost optimization model builder designed to evaluate long-term national energy scenarios at a nodal network level. Built on the [PyPSA](https://pypsa.org/){:target="_blank"} framework, it adopts a multi-sectoral cost optimization approach with a primary focus on the power system.

The diagram below illustrates the energy flows within a single node:

[![PyPSA-SPICE overview schema](../assets/images/pypsa-spice_schema_overview.svg){ .img-center width="100%" }](../assets/images/pypsa-spice_schema_overview.svg){: target="_blank" }

Each component in this model follows the definitions from [PyPSA components](https://pypsa.readthedocs.io/en/latest/user-guide/components.html){:target="_blank"} and the system is structured into three main sectors:

- [Power sector](power_sector.md)
- [Industry sector](industry_sector.md)
- [Transport sector](transport_sector.md)

You can find more detailed information within each sector’s dedicated section.
