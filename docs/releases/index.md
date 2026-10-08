!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming
- Add legend ordering control for stacked column charts in PyPSA-SPICE-Vis.
- Add minimum curtailment support with curtailment penalty in the optimization objective.
- Support piecewise linear efficiency and rate curves for thermal power units.

### Fixed
- Fix `solve_network` checking the model status of the global network instead of the solved one, and raise an error when the model is infeasible. ([:material-source-pull:124](https://github.com/agoenergy/pypsa-spice/pull/124) by @nhlong2701)

### Changed
- Add ramp-up and ramp-down costs for thermal power units modelled as generators or links, activated per country via the new `ramp_costs` custom constraint. Ramp costs are included in the OPEX output tables and reported separately in `pow_ramp_cost_by_type_yearly`. ([:material-source-pull:124](https://github.com/agoenergy/pypsa-spice/pull/124) by @nhlong2701)

### Notes

- **New ramp cost columns in `technologies.csv`:** We have added the `ramp_up_cost` and `ramp_down_cost` columns (currency/MW) to `technologies.csv`. Existing files without these columns still work, with ramp costs set to 0. To use ramp costs, add the columns to your `technologies.csv`, rebuild the networks, and activate `ramp_costs` in the `custom_constraints` section of `scenario_config.yaml`.

--8<-- "releases/v2.1.0.md"
--8<-- "releases/v2.0.0.md"
--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
