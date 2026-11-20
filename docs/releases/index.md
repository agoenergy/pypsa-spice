!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming
- Add legend ordering control for stacked column charts in PyPSA-SPICE-Vis.
- Add minimum curtailment support with curtailment penalty in the optimization objective.
- Support piecewise linear efficiency and rate curves for thermal power units.

### Fixed
- Fix `build_skeleton` creating the global cost files (`power_plant_costs.csv` and `storage_costs.csv`) with the template years and USD columns instead of the years and currency from `base_config.yaml`. ([:material-source-pull:129](https://github.com/agoenergy/pypsa-spice/pull/129) by @sasa821)
- Fix `solve_network` checking the model status of the global network instead of the solved one, and raise an error when the model is infeasible. ([:material-source-pull:124](https://github.com/agoenergy/pypsa-spice/pull/124) by @nhlong2701)

### Changed
- Add ramp-up and ramp-down costs for thermal power units modelled as generators or links, activated per country via the new `ramp_costs` custom constraint. Ramp costs are included in the OPEX output tables and reported separately in `pow_ramp_cost_by_type_yearly`. ([:material-source-pull:124](https://github.com/agoenergy/pypsa-spice/pull/124) by @nhlong2701).
- Add load dumping generators to all buses (except CO~2~ buses) to avoid infeasibility when must-run generation exceeds the load, and report the buses with load shedding or load dumping in the post-analysis. ([:material-source-pull:126](https://github.com/agoenergy/pypsa-spice/pull/126) by @nhlong2701).
- Add optional demand side management (DSM) per country and load type in the new `demand_side_management` section of `scenario_config.yaml`. Load can be shed at a marginal cost (`shed`) or shifted within a fixed window (`shift`), capped at a share of the peak of the load on each bus. DSM is reported in `dsm_dmd_by_bus_yearly` and `dsm_cost_by_bus_yearly`, separate from load shedding and load dumping. ([:material-source-pull:PR_NUMBER](https://github.com/agoenergy/pypsa-spice/pull/PR_NUMBER) by @samarthiith)
- Add the PyPSA-SPICE web workspace (beta) in `pypsa-spice-web/`, a React and FastAPI interface to edit model inputs, compare scenario inputs, edit scenario configurations and solver settings, run the Snakemake workflow, and explore results and custom dashboards for a local project. ([:material-source-pull:127](https://github.com/agoenergy/pypsa-spice/pull/127) by @samarthiith and @nhlong2701).

### Notes

- **New ramp cost columns in `technologies.csv`:** We have added the `ramp_up_cost` and `ramp_down_cost` columns (currency/MW) to `technologies.csv`. Existing files without these columns still work, with ramp costs set to 0. To use ramp costs, add the columns to your `technologies.csv`, rebuild the networks, and activate `ramp_costs` in the `custom_constraints` section of `scenario_config.yaml`.

- **Load dumping generators:** Load dumping generators (`DUMPLOAD - <bus>`, type `LDMP`, carrier `EXS`) are added to the base year network next to the load shedding generators. Please rebuild the base year network to include them. The `test_energy_not_served_warning` function in the post-analysis has been renamed to `report_load_shedding_and_dumping`.

- **Demand side management:** DSM is only applied with hourly resolution (`nth_hour` with `stepsize: 1`). With any other resolution, a warning is logged and DSM is skipped. Loads now carry their `profile_type` as an attribute, and the carriers `DSM_SHED` and `DSM_SHIFT` are added in the base year, so please rebuild the base year network before using DSM. See [demand side management](../getting-started/input-data/model-builder-configuration.md#scenario_configyaml-demand-side-management) for the options and the formulation.

- **Web workspace (beta):** The web workspace is under active development and is not a fully released product, so features and layouts may change. It saves edits directly to your input files, so keep your data under version control. To try it, install [Node.js](https://nodejs.org/), update your environment (using `conda env update -f envs/environment.yaml --prune`), and start it with `python pypsa-spice-web/app.py` from the activated environment. See the [web workspace documentation](../visualisation-tool/web-workspace.md) for details.

--8<-- "releases/v2.1.0.md"
--8<-- "releases/v2.0.0.md"
--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
