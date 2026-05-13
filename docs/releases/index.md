!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming

- Improve the input GUI’s handling of custom constraints and session state.
- Add minimum curtailment support via a curtailment penalty in the optimisation objective.

### Fixed

- Change website link display. ([:material-source-pull:110](https://github.com/agoenergy/pypsa-spice/pull/110) by @RichChang963)

### Changed

- Add customised legend order control for hourly stacked column charts, and add back-to-top button in PyPSA-SPICE-Vis ([:material-source-pull:104](https://github.com/agoenergy/pypsa-spice/pull/104) by @RichChang963)

### Notes

- **Add streamlit-sortables package:** We have added the `streamlit-sortables==0.3.1` to handle legend orders in PyPSA-SPICE-Vis. We encourage users to upgrade their environments (using `conda env update -f envs/environment.yaml --prune`) to ensure full compatibility with the latest features and improvements.

--8<-- "releases/v2.0.0.md"
--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
