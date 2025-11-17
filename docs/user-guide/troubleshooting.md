<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Troubleshooting

## Setup debugger

To test the workflow using Snakemake rules, a configuration file `launch.json` with the following content is required:

```json title="Debugger configurations"
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "mock snakemake",
            "type": "python",
            "request": "launch",
            "program": "${file}",
            "console": "integratedTerminal",
            "justMyCode": false,
            "cwd": "${cwd}${pathSeparator}scripts"
        }
    ]
}
```

Once this file is set up, you can easily run any Python file in the `scripts` folder under debugger mode for testing.

## Infeasibility issues

If your model is **infeasible** or **unbounded**, it often means that your input settings are leading to a situation where PyPSA can’t solve one or more objective functions. Common causes are:

- Insufficient generator capacity at a bus to meet the load across all snapshots.
- Must-run generators (`p_min_pu`) produce more power than the load at certain snapshots, causing excess power (i.e., power dumping).
- Negative values assigned to capacity parameters.
- Maximum capacity (`p_nom_max`) is smaller than minimum capacity (`p_nom_min`).
- Storage units has inflow and cyclic charging (`cyclic_state_of_charge = True`) but no discharging capability.

## Tips for diagnosing the problem

Try the following steps to identify and fix the issue:

- Run `pypsa.network.consistency_check()` to check if any warnings appear.
- Temporarily disable custom features or user extensions to isolate the cause.
- Ensure load-shedding generators are added to _every_ bus.
- Check `pypsa.network.generators.p_min_pu` to identify any must-run generators. Then verify if their generation exceeds the corresponding load in `pypsa.network.loads_t.p_set`.
- Compare the `p_nom_max` and `p_nom_min` values for each component to ensure they make sense (i.e., max ≥ min).
- Try disabling cyclic charging for storage units: `pypsa.network.storage_units.cyclic_state_of_charge = False`.
- Double-check for any negative values in:
    - `pypsa.network.{component}.p_nom`
    - `pypsa.network.{component}.p_nom_max`
    - `pypsa.network.{component}.p_nom_min`

## More detailed information

[PyPSA](https://pypsa--1250.org.readthedocs.build/en/1250/user-guide/troubleshooting/#optimisation-convergence-infeasibility) also provides official guidance on solving infeasiblity issues. You can also explore the link to seek for a good solution.
