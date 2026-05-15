!!! info
    The features and bugfixes listed under **Upcoming** aren't released yet but will be included in the next version. If you'd like to try them early, you can switch to the `develop` branch. Just keep in mind that it's not stable and may contain issues. All previous releases are already available on the `main` branch.

## Upcoming
- Add a legend order control for stacked column charts in PyPSA-SPICE-Vis.
- Add minimum curtailment support via a curtailment penalty in the optimisation objective.

### Fixed
- Change website link display. ([:material-source-pull:110](https://github.com/agoenergy/pypsa-spice/pull/110) by @RichChang963)
- Fix error handling of missing p_max_pu in generators and links. ([:material-source-pull:107](https://github.com/agoenergy/pypsa-spice/pull/107) by @RichChang963)

### Changed
- Add streamlit-extra pacakge and enhance addtional features in PyPSA-SPICE-Vis. ([:material-source-pull:105](https://github.com/agoenergy/pypsa-spice/pull/105) by @RichChang963)

### Notes

- **New Streamlit-extras package:** We have added a streamlit-extra `v1.5.0` and encourage users to upgrade their environments (using `conda env update -f envs/environment.yaml --prune`) to ensure full compatibility with the latest features and improvements.

--8<-- "releases/v2.0.0.md"
--8<-- "releases/v1.1.1.md"
--8<-- "releases/v1.1.0.md"
--8<-- "releases/v1.0.0.md"
