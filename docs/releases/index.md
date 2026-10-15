!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming
- Add legend ordering control for stacked column charts in PyPSA-SPICE-Vis.
- Add minimum curtailment support with curtailment penalty in the optimization objective.
- Add ramp rate costs for thermal power units.
- Support piecewise linear efficiency and rate curves for thermal power units.

### Fixed


### Changed
- Add load dumping generators to all buses (except CO~2~ buses) to avoid infeasibility when must-run generation exceeds the load, and report the buses with load shedding or load dumping in the post-analysis. ([:material-source-pull:126](https://github.com/agoenergy/pypsa-spice/pull/126) by @nhlong2701)

### Notes

- **Load dumping generators:** Load dumping generators (`DUMPLOAD - <bus>`, type `LDMP`, carrier `EXS`) are added to the base year network next to the load shedding generators. Please rebuild the base year network to include them. The `test_energy_not_served_warning` function in the post-analysis has been renamed to `report_load_shedding_and_dumping`.

--8<-- "releases/v2.1.0.md"
--8<-- "releases/v2.0.0.md"
--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
