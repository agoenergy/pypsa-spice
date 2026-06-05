!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming
- Add minimum curtailment support via a curtailment penalty in the optimisation objective.

### Fixed
- Remove polyfill package to enhance website security. ([:material-source-pull:113](https://github.com/agoenergy/pypsa-spice/pull/113) by @RichChang963)
- Fix correct barmode in plotly to display right stacked column charts. ([:material-source-pull:111](https://github.com/agoenergy/pypsa-spice/pull/111) by @RichChang963)
- Change website link display. ([:material-source-pull:110](https://github.com/agoenergy/pypsa-spice/pull/110) by @RichChang963)
- Fix error handling of missing p_max_pu for generators and links. ([:material-source-pull:107](https://github.com/agoenergy/pypsa-spice/pull/107) by @RichChang963)

### Changed

- Add customised legend order control for hourly stacked column charts, and add back-to-top button in PyPSA-SPICE-Vis ([:material-source-pull:104](https://github.com/agoenergy/pypsa-spice/pull/104) by @RichChang963)
- Add streamlit-extra pacakge and enhance addtional features in PyPSA-SPICE-Vis. ([:material-source-pull:105](https://github.com/agoenergy/pypsa-spice/pull/105) by @RichChang963)

### Notes

- **Add streamlit-sortables package:** We have added the `streamlit-sortables==0.3.1` to handle legend orders in PyPSA-SPICE-Vis. We encourage users to upgrade their environments (using `conda env update -f envs/environment.yaml --prune`) to ensure full compatibility with the latest features and improvements.
- **New Streamlit-extras package:** We have added a streamlit-extra `v1.5.0` and encourage users to upgrade their environments (using `conda env update -f envs/environment.yaml --prune`) to ensure full compatibility with the latest features and improvements.

--8<-- "releases/v2.0.0.md"
--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
