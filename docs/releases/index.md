# Unreleased

## Fixed

- [Bug] Filtering issue at add_load for non-electricity loads and missing interconnection by @nhlong2701 in [#51](https://github.com/agoenergy/pypsa-spice/pull/51)

## Changed

- Update pypsa version to 1.0.6 and refactor code to adapt to new changes by @nhlong2701 in [#54](https://github.com/agoenergy/pypsa-spice/pull/54)

- Update scenario_config workflow and template structure by @RichChang963 in [#58](https://github.com/agoenergy/pypsa-spice/pull/58)

- Feature/add generation constraint by @RichChang963 in [#59](https://github.com/agoenergy/pypsa-spice/pull/59)

## Notes

- Upgraded PyPSA version: PyPSA has been updated to v1.0.6. Users are encouraged to upgrade their environments to ensure full compatibility with the latest features and improvements.
- Automatic `scenario_config.yaml` creation: The `build_skeleton` rule now automatically generates a `scenario_config.yaml` file, simplifying initial project setup.
- A  can be created automatically while executing `build_skeleton` rule.
- New `active` parameter in custom constraints: A new boolean control field, `active`, has been introduced for defining custom constraints within `scenario_config.yaml`. Users should update their configuration format as described in the [documentation](https://agoenergy.github.io/pypsa-spice/getting-started/input-data/model-builder-configuration/#scenario_configyaml-custom-constraints) to maintain compatibility and avoid configuration errors.
