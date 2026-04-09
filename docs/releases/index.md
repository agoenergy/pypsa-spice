!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming

### Fixed

- fix build skeleton error and only use build skeleton when the project folder is not initialised.([:material-source-pull:88](https://github.com/agoenergy/pypsa-spice/pull/88) by @RichChang963)
- fix wrong definition for coal in documentation and remove site name `Model Builder Documentation` in the documentation website.([:material-source-pull:101](https://github.com/agoenergy/pypsa-spice/pull/101) by @RichChang963)

### Changed

- suppress PyPSA logging warnings for network export and import before solve ([:material-source-pull:97](https://github.com/agoenergy/pypsa-spice/pull/97) by @nhlong2701)
- Add a new comparison bar chart if there is deviation between the two selected scenarios in each inidcator ([:material-source-pull:91](https://github.com/agoenergy/pypsa-spice/pull/91) by @samarthiith & @nhlong2701)
- Update streamlit to `v1.55.0` and clean up GUI for input data and scenario configs.([:material-source-pull:92](https://github.com/agoenergy/pypsa-spice/pull/92) by @RichChang963)

### Notes

- **Streamlit version upgrade:** We have upgraded PyPSA to `v1.55.0` and encourage users to upgrade their environments (using `conda env update -f envs/environment.yaml --prune`) to ensure full compatibility with the latest features and improvements.

--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
