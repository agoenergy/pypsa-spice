<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# PyPSA-SPICE web workspace (beta)

!!! warning "Beta version"
    The web workspace is under active development and is **not a fully released product**.
    Features, page layouts, and the API may change between versions, and some workflows are still incomplete.
    The workspace saves edits **directly to your input CSV and YAML files**, so keep your data under version
    control (for example Git) and check model results before using them.

    Please report bugs and ideas in the [GitHub issues](https://github.com/agoenergy/pypsa-spice/issues).

The web workspace is a browser interface for a local PyPSA-SPICE project. It runs on your own machine and works
directly on the files in the `data/` folder. It follows the modelling workflow:

- **Inputs & configuration**
    - **Input data**: view and edit the scenario's CSV inputs, grouped by technology.
    - **Compare inputs**: see configuration and input-data differences between two scenarios.
    - **Scenario configuration**: edit general settings, CO₂ management, and custom constraints in
      `scenario_config.yaml`.
- **Review & run**: check the run, choose the solver and its settings, and start the Snakemake workflow.
- **Results & analysis**: explore Power, Industry, Transport, Emissions, and Costs results, and compare two runs.
- **Dashboard builder**: combine result charts into custom dashboards, stored in your browser.

The existing [PyPSA-SPICE-Vis](pypsa-spice-vis.md) Streamlit tool remains available while the web workspace is in beta.

## How to start it

The web workspace needs the `hotpot` Conda environment and [Node.js](https://nodejs.org/) with `npm`.
From the repository root, install the dependencies once:

```bash
conda env update -n hotpot -f envs/environment.yaml
npm install --prefix pypsa_spice_web/frontend
```

Then start the workspace:

```bash
./run-web-locally.sh
```

Open <http://127.0.0.1:5173/ui/> in your browser.

Developer documentation for the web workspace is in `pypsa_spice_web/docs/overview.md`.
