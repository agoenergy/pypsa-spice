<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# PyPSA-SPICE: PyPSA-based Scenario Planning and Integrated Capacity Expansion

<!-- badges-begin -->
[![License][license badge]][license]{:target="_blank"}
[![PyPSA version][PyPSA version badge]][PyPSA version]{:target="_blank"}
[![Snakemake][Snakemake badge]][Snakemake]{:target="_blank"}
[![Code style][Code style badge]][Code style]{:target="_blank"}

[license badge]: https://eddelbuettel.github.io/badges/GPL2+.svg
[license]: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html

[PyPSA version badge]: https://img.shields.io/pypi/v/pypsa?label=pypsa
[PyPSA version]: https://pypi.org/project/pypsa/

[Snakemake badge]: https://img.shields.io/badge/snakemake-minimal==8.10.8-brightgreen.svg?style=flat
[Snakemake]: https://snakemake.readthedocs.io

[Code style badge]: https://img.shields.io/badge/code%20style-black-000000.svg
[Code style]: https://github.com/ambv/black

<!-- badges-end -->

!!! info
    If you are considering using this model builder, please reach out to us at [modelling@agora-thinktanks.org](mailto:modelling@agora-thinktanks.org){:target="_blank"}. We would be happy to help you get started.
    If you encounter a bug, please create a [new issue](https://github.com/agoenergy/pypsa-spice/issues){:target="_blank"}. For new ideas or feature requests, you can start a conversation in the [discussions](https://github.com/agoenergy/pypsa-spice/discussions){:target="_blank"} section of the repository.

PyPSA-SPICE (https://agoenergy.github.io/pypsa-spice/) is an open-source model builder for assessing national mid/long-term energy scenarios using a least-cost, multi-sectoral optimization approach based on the [PyPSA](https://pypsa.org/){:target="_blank"} framework. It can be used to build models that represent one or more countries across multiple interconnected nodes linked by electricity transmission, and within each region, it models the integration of power, heat, and transport sectors. The model focuses on the power sector, which is represented with a high level of detail.

The model workflow has been designed to be more accessible compared to other PyPSA-based models, though basic Python coding knowledge is required.

![PyPSA-SPICE_single_node_energy_flow](assets/images/pypsa-spice_intro.jpg)

## Key Features of PyPSA-SPICE

- Assessment of national or regional long-term energy scenarios using a least-cost optimization approach based on the [PyPSA](https://pypsa.org/){:target="_blank"} framework.
- Co-optimization of generation, capacity, and interconnector expansion at hourly resolution.
- Power plants are represented at technology resolution, with user-defined clustering within technologies as needed.
- Straightforward model creation for new countries and/or regions defined by the user.
- Easy integration of custom data into the model.
- Several pre-defined custom constraints including energy independence, reserve margin, and must-run constraints on thermal generators.
- Flexible sectoral coverage: base power sector model can be complemented with industry and transport sectors for full energy system investigation.
- [PyPSA-SPICE-Vis](visualisation-tool/pypsa-spice-vis.md) as a visual tool for easy visualization of model outputs.
- Extensive documentation to facilitate working with the model.

## Types of models outside the scope of PyPSA-SPICE

PyPSA-SPICE is not designed for modeling power or energy systems with very high geographic resolution, such as load flow modeling. Instead, it prioritizes ease of building country- or region-level models using manually provided custom data. For this reason, it is best suited for models with a maximum of 10–15 regional nodes.

- For higher geographic detail other open-source PyPSA-based frameworks can be used:
  - [PyPSA-Eur](https://github.com/pypsa/pypsa-eur){:target="_blank"}
  - [PyPSA-meets-Earth](https://github.com/pypsa-meets-earth/pypsa-earth){:target="_blank"}

- PyPSA-SPICE requires total energy demand as an input and does not optimize total demand, modal shifts, or other demand-side dynamics.  

## Studies conducted using PyPSA-SPICE

Please refer to [Publication](tutorials-and-examples/publication.md) to see the publications or studies
using the PyPSA-SPICE model builder.

## Visualisation of PyPSA-SPICE Results

To make it easy to visualise and compare scenario outputs, we provide an open-source
library [PyPSA-SPICE-Vis](visualisation-tool/pypsa-spice-vis.md) within the model builder.

## Citing PyPSA-SPICE

Please use the citation below:

- Agora Think Tanks (2025): PyPSA-SPICE: PyPSA-based Scenario Planning and Integrated Capacity Expansion

## Contributing

We welcome contributions from anyone interested in improving this project. Please take a moment to review our [Contributing Guide](contributing/contributing.md) and [Code of Conduct](contributing/code_of_conduct.md).
If you have ideas, suggestions, or encounter any issues, feel free to open an issue or submit a pull request on GitHub.

## Maintained by

[![AgoraEW](assets/images/Agora_EW.png)](https://www.agora-energiewende.org/)

## Supported by

[![CASE](assets/images/CASE.png)](https://caseforsea.org/) [![INETTT](assets/images/INETTT.png)](https://www.inettt.org/)

## License

Copyright &copy; [PyPSA-SPICE Developers](references/developers.md)

PyPSA-SPICE is licensed under the open source [GNU General Public License v2.0 or later](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
with the following information:

The documentation is licensed under [CC-BY-4.0](https://interoperable-europe.ec.europa.eu/licence/creative-commons-attribution-40-international-cc-40).

The repository uses [REUSE](https://reuse.software/) to expose the licenses of its files.
