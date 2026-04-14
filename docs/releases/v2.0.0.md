## v2.0.0 (2026-04-14)

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
- **Correction of file name for decommission_capacity:** The file name for decommissioning capacity has been corrected to `decommission_capacity.csv` in the documentation and codebase. Please make sure to update your files accordingly to avoid any issues with the model execution.