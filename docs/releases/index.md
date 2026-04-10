!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming
- Add a legend order control for stacked column charts in PyPSA-SPICE-Vis.
- Improve the input GUI’s handling of custom constraints and session state.
- Add minimum curtailment support via a curtailment penalty in the optimisation objective.

### Fixed

- Fix build skeleton error and only use build skeleton when the project folder is not initialised. ([:material-source-pull:88](https://github.com/agoenergy/pypsa-spice/pull/88) by @RichChang963)
- Fix wrong definition for coal in documentation and remove site name `Model Builder Documentation` in the documentation website. ([:material-source-pull:101](https://github.com/agoenergy/pypsa-spice/pull/101) by @RichChang963)

### Changed

- Refactor PyPSA-SPICE-VIS output modules to improve readability and maintainability. ([:material-source-pull:75](https://github.com/agoenergy/pypsa-spice/pull/75) by @nhlong2701)
- Suppress PyPSA logging warnings for network export and import before solve. ([:material-source-pull:97](https://github.com/agoenergy/pypsa-spice/pull/97) by @nhlong2701)
- Add a new comparison bar chart if there is deviation between the two selected scenarios in each indicator. ([:material-source-pull:91](https://github.com/agoenergy/pypsa-spice/pull/91) by @samarthiith & @nhlong2701)
- Add descriptions in the landing page of the PyPSA-SPICE-Vis tool. ([:material-source-pull:100](https://github.com/agoenergy/pypsa-spice/pull/100) by @RichChang963)
- Update Streamlit to `v1.55.0` and introduce a GUI for input data and scenario configs. ([:material-source-pull:92](https://github.com/agoenergy/pypsa-spice/pull/92) by @RichChang963)

### Notes

- **PyPSA-SPICE-VIS remake:** Significantly refactored the PyPSA-SPICE-VIS Streamlit app for readability and maintainability. Reorganised modules, and now manage chart types and their configurations separately. We have also introduced two new pages for working directly with input data and scenario configs in the app.
- **Streamlit version upgrade:** We have upgraded streamlit to `v1.55.0` and encourage users to upgrade their environments (using `conda env update -f envs/environment.yaml --prune`) to ensure full compatibility with the latest features and improvements.

--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
