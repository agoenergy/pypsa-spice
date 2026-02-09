# Adding a new technology logic

This note explains how to add a new technology in PyPSA-SPICE.
The procedure depends on the type of component, so we provide several common examples for reference:

!!! info
    Adding a new technology requires a solid understanding of both PyPSA and PyPSA-SPICE to avoid errors during model execution.
    If you are uncertain about any changes, please reach out to us at [modelling@agora-thinktanks.org](mailto:modelling@agora-thinktanks.org){:target="_blank"}. We would be happy to help you get started.

## A new fossil generator

In the default setting in `technologies.csv`, fossil fuels are defined as `links` which attached with the corresponding `fuel buses` as fossil fuel hubs. Taking `CCGT_Drived` as an exmple, you can:

- In `technologies.csv`, copy the rows of CCGT into new row, rename `CCGT` into `CCGT_Drived`, and adjust the other techical parameters.
- In `power_plant_costs.csv`, add new rows and enter the cost assumptions of `CCGT_Drived`.
- In `power_links.csv`, check which country/region that you would like to set this technology in, copy the rows of CCGT in the same country/region into new row, rename `CCGT` into `CCGT_Drived`, and adjust `p_nom_min` and `p_nom_max` in each year.
- In `decomission_capacity.csv`, add new decomission plan of `CCGT_Drived` if there's any.

## A new renewable generator

In the default setting in `technologies.csv`, renewales are defined as `Generators` (except hydro dam and hydro pumped storage which are `StorageUnits`). Taking `Wind-PV-Hybrid` as an exmple, you can:

- In `technologies.csv`, copy the rows of Wind or Solar into new row, rename the technology into `Wind-PV-Hybrid`, and adjust the other techical parameters.
- In `power_plant_costs.csv`, add new rows and enter the cost assumptions of `Wind-PV-Hybrid`.
- In `renewables_technical_potential.csv`, add the technical potentials for `Wind-PV-Hybrid`. This defines the overall p_nom_max of this technology in different countries/regions.
- In `availability.csv`, add the timeseries profile of capacity factor for `Wind-PV-Hybrid`.
- In `power_generator.csv`, check which country/region that you would like to set this technology in, copy the rows of CCGT in the same country/region into new row, rename technology into `Wind-PV-Hybrid`, and adjust `p_nom_min` and `p_nom_max` in each year.

## A battery with a fixed E/P ratio
