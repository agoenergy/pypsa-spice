<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input data: model execution

The execution of the model is orchestrated using [Snakemake](https://snakemake.readthedocs.io/en/stable/){:target="_blank"}, as defined in the project’s `Snakefile` located in the root directory. Snakemake acts as the workflow engine, coordinating individual rules and chaining them into a single, automated execution process. This allows for scalable, parallel computation across multiple cores or threads.

## Execution of Snakemake workflow

Below is a simplified representation, showing the workflow for the years 2025-2050 with the power sector only.

![PyPSA-SPICE_snakemake_workflow_DAG](../../assets/images/snakemake_workflow_dag.jpg)

## Rules and functions

The table below explains the main Snakemake rules and helper functions used in the PyPSA-SPICE workflow:

| Name                    | Type           | Description                         |
| ----------------------- | -------------- | ----------------------------------- |
| `add_baseyear`          | Snakemake rule | First step of building up the model after input data preparation is done. This rule imports all input data in the base year into the model. The base year is defined in `config.yaml` (see [new-model](new-model.md) or [model-builder-configuration](model-builder-configuration.md)).|
| `add_brownfield`        | Snakemake rule | Imports all output data from the previous year, using the `previous_year_outputs` and `solved_previous_year` functions. The definition of the different years is described in `config.yaml` (see [new-model](new-model.md) or [model-builder-configuration](model-builder-configuration.md)). |
| `make_summary`          | Snakemake rule | Generates output CSVs and summary of a specific year. |
| `combine_summaries`     | Snakemake rule | Last step of the whole workflow aside from `solve_all_networks` rule. It generates summary of all years after all the networks from different years are solved.  |
| `solve_network`         | Snakemake rule | Solves the model builder for a specific year. This rule is only triggered after either `add_baseyear` or `add_brownfield`. Optimisation settings are defined in `config.yaml` (see [model-builder-configuration](model-builder-configuration.md)). |
| `solve_all_networks`    | Snakemake rule | Executes the entire model builder workflow at once. |
| `previous_year_outputs` | Python Function | Reads the output CSVs of the previous year. |
| `solved_previous_year`  | Python Function | Reads the network from the output of the previous year. |

## Running the model builder

### Run the entire workflow at once

To execute the entire workflow, modify the `configfile` parameter in the snakemake to point to your scenario `config_scenario.yaml` file. After this, you can run the whole workflow using:

```bash title="Run the entire model builder workflow at once."
snakemake -j1 -c4 solve_all_networks #(1)!
```

1. `-j1` means running only 1 job in parallel at a time and `-c4` means allowing each job to use up to 4 CPU threads. See the [Snakemake CLI documentation](https://snakemake.readthedocs.io/en/stable/executing/cli.html){:target="_blank"} to customise cores, threads, and parallel execution.

or

```bash title="Running the entire model builder workflow using this call command"
snakemake -call
```

### Run a single year

If you want to run a specific year instead of running multiple years, You can do it using the following command:

```bash title="Option 2: Run a single output file."
snakemake -j1 -c4 post-solve/elec_{SECTOR}_{YEAR}.nc
```

Since a particular year run needs previous year outputs (except base year), snakemake will run the all the required worksflows to generate the output for specified year.
