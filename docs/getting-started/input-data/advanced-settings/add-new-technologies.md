# Adding a new technology logic

This note explains how to add new or custom technologies in PyPSA-SPICE that are not included in the default set. In PyPSA-SPICE, some energy technologies are modeled by combining multiple PyPSA components to capture the characteristics of a specific plant type. Below are a few common examples for reference.

!!! info
    Adding new or custom technologies often requires a good understanding of both the PyPSA framework and PyPSA-SPICE to avoid errors during model execution.

    If you are uncertain about making these changes, please reach out to us at [modelling@agora-thinktanks.org](mailto:modelling@agora-thinktanks.org). We would be happy to help you get started.

## Example 1: Adding a multi-output plant (e.g., a new fossil plant with emissions)

In PyPSA-SPICE, a **multi-output plant** is modeled using a PyPSA `Link`, where the fuel bus (`bus0`) provides the input and the output buses (`bus1`, `bus2`, etc.) represent the plant’s different outputs. For example, a gas-fired power plant (CCGT) can be implemented as a `Link` that connects the **gas bus** (`bus0`) to the **electricity/power sector bus** (`bus1`) as its main output, and to an **atmosphere bus** (`bus2`) as a second output to represent CO₂ emissions released as a by-product of combustion; additional outputs can be included by connecting to further buses such as `bus3`, `bus4`, and so on.

To add a new gas-fired power plant type (e.g. `CCGT_Mod`), follow these steps:

- In **`technologies.csv`**: Add a new row for `CCGT_Mod` and fill in the technical parameters. Note: `efficiency` applies to the conversion from **`bus0 → bus1`**, `efficiency2` to **`bus0 → bus2`**, and so on for additional outputs.
- In **`power_plant_costs.csv`**: Add a new row for `CCGT_Mod` for each modeled year and enter the cost assumptions.
- In **`power_links.csv`**: Add a new row for `CCGT_Mod`, fill in required input, and make sure to set the buses correctly: **`bus0` = fuel**, **`bus1` = power sector**, **`bus2` = atmosphere**.
- In **`decomission_capacity.csv`**: Add a decommissioning plan for `CCGT_Mod`, if applicable.

!!! tip
    In most cases, the new plant type is only slightly different from an existing one (e.g. `CCGT_Mod` vs `CCGT`).

    To save time and reduce errors, copy the `CCGT` rows in all four files, change the relevant text to `CCGT_Mod`, and then modify only the parameters that differ.


## Example 2: Adding a new renewable plant with time-dependent availability

In PyPSA-SPICE, renewable plants with time-dependent availability are modeled as a PyPSA `Generator` (hydro dams are the exception and use `StorageUnit`). To add a new renewable plant (e.g. `Wind-PV-Hybrid`), do the following:

- In **`technologies.csv`**: Add a new row for `Wind-PV-Hybrid` and fill in the technical parameters. Leave `p_max_pu` and `p_min_pu` empty as renewable plant's availability is defined in `availability.csv`.
- In **`power_plant_costs.csv`**: Add a new row for each modeled year and enter the cost assumptions for `Wind-PV-Hybrid`.
- In **`renewables_technical_potential.csv`**: Add new rows for `Wind-PV-Hybrid` that define the maximum buildable capacity by region (per country).
- In **`availability.csv`**: Add the time-dependent availability profile (per unit) for `Wind-PV-Hybrid` by region/country.
- In **`power_generator.csv`**: Add a new row for `Wind-PV-Hybrid`, fill in the required inputs, and make sure to set the bus information correctly.

## Example 3: Adding a storage plant with fixed energy-to-power ratio

In PyPSA-SPICE, all storage plants are built from a combination of PyPSA `Store` and `Link` components, which lets the model optimize **energy capacity** and **power capacity** independently. If you want a custom storage type with a **fixed energy-to-power ratio**, use a PyPSA `StorageUnit` instead. For example, to add a hydro-pumped storage plant (`HPHS_Mod`) with an **energy-to-power ratio of 8**, follow these steps:

- In **`technologies.csv`**: Add a new row for `HPHS_Mod` and fill in the technical parameters. Set `max_hours = 8`. Note that storage can have different charge/discharge efficiencies: use `efficiency` for **discharging** and `efficiency_store` for **charging**.
- In **`power_plant_costs.csv`**: Add a new row for each modeled year and enter cost assumptions for `HPHS_Mod`. CAPEX and FOM in the csv are in **$/MW**. If your costs are in **$/MWh**, multiply by **8** to convert to **$/MW**.
- In **`storage_inflows.csv`**: Add the time-dependent inflow profile (MW) for `HPHS_Mod` by region/country (if applicable).
- In **`storage_capacity.csv`**: Add a new row for `HPHS_Mod`, fill in the required inputs, and make sure to set the bus information correctly.

