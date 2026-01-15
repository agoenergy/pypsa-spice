# Unreleased

## Fixed

- [Bug] Fix issue with regional filters when adding load for non-electricity ones and missing interconnection by @nhlong2701 in [#51](https://github.com/agoenergy/pypsa-spice/pull/51)

## Changed

- Update pypsa version to v1.0.6 and refactor code to adapt to new changes by @nhlong2701 in [#54](https://github.com/agoenergy/pypsa-spice/pull/54)

- Add scenario config template and update build skeleton workflow by @RichChang963 in [#58](https://github.com/agoenergy/pypsa-spice/pull/58)

- Add custom maximum power generation constraint by @RichChang963 in [#59](https://github.com/agoenergy/pypsa-spice/pull/59)

## Notes

- PyPSA version has been updated to v1.0.6. Users are encouraged to upgrade their environments to ensure full compatibility with the latest features and improvements.
- The `build_skeleton` rule now automatically generates a `scenario_config.yaml` file under scenario data folder, simplifying initial project setup.
- A new boolean control field, `active`, has been introduced when defining new custom constraints within `scenario_config.yaml`. Users should update their configuration format as described in the [documentation](https://agoenergy.github.io/pypsa-spice/getting-started/input-data/model-builder-configuration/#scenario_configyaml-custom-constraints) to maintain compatibility and avoid configuration errors.
