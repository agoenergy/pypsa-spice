<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input data: define a new model

This section explains how to set up a new model for a particular country/region using the PyPSA-SPICE model builder. Once your model is established you can run the model as described in [Model execution](model-builder-execution.md).

Recommended Steps of setting up a new model:

1. [Adjust the information inside `base_config.yaml`](#step-1-set-up-the-base-configuration-file).
2. [Run `snakemake -c1 build_skeleton`](#step-2-build-the-skeleton) to create a folder structure and template CSV files for your input data. All input folders and files shall be created inside the `data` folder after this command is executed.
3. [Save your input data systematically](#step-3-save-your-input-data-systematically)
4. [Fill in the skeleton CSVs](#step-4-fill-in-the-skeleton-csvs) with the required data manually or using available resources.
5. [Fill the scenario assumptions and constraints in `scenario_config.yaml` in your scenario folder](#step-5-fill-in-the-scenario-assumptions-and-constraints-in-scenario_configyaml).

An example structure created by ``build_skeleton`` is displayed below. The following sections will use this example to explain the settings.

## Step 1: set up the base configuration file

Setting up the base config requires defining the scope and resolution of the model. Specifically, defining which countries/regions will be represented in the model and for which year the model will be run. While it is possible to change these after initial model is created, it would require significant effort to add new regions/years in the input CSVs.  

The `base_config.yaml` contains two parts which will be used both for folder structure building and model executions:

```yaml title="Init settings in the config.yaml file"
path_configs: #(1)!
  data_folder_name: example
  project_name: project_01
  input_scenario_name: scenario_01 # (2)!
  output_scenario_name: scenario_01_tag1 # (3)!

base_configs:
  regions: 
    XY: ["NR","CE", "SO"] # (4)!
    YZ: ["NR","CE", "SO"] 
  years: [2025, 2030, 2035, 2040, 2045, 2050] # (5)!
  sector: ["p-i-t"] # (6)!
  currency: USD # (7)!
```

1. This section is for configuring directory structure for storing model inputs and results.
2. A custom name you define for the input scenario folder in your model.
3. A custom name you define for the output scenario folder to save your model results.
4. List of regions or nodes within each country. This defines the network’s nodal structure. The country list contains **2-letter** country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"}.
5. Modelled years should be provided as a list.
6. Options: [`p`, `p-i`, `p-t`, `p-i-t`], representing power (`p`), industry (`i`), and transport (`t`) sectors.
7. Currency usd in the model. The default setting is USD (also used in example data). Format shall be in all uppercases, [ISO4217](https://www.iso.org/iso-4217-currency-codes.html){:target="_blank"} format.

The final skeleton folder path will follow this structure:: `data`/`data_folder_name`/`project_name`.

By setting different ``input_scenario_name`` and country or regional settings in the ``base_configs`` section (see details in [Model configuration](model-builder-configuration.md)), a new skeleton structure under the same `data_folder_name` folder will be created.

## Step 2: build the skeleton

After modifying the configuration file, run the following command in your terminal.

```bash title="Generating the skeleton folder"
snakemake -c1 build_skeleton
```

This step creates your skeleton folder and files which can be feed with your data.

```text title="Structure of Folder and files created by build skeleton script"
📦 data
 ┗ 📂 example
    ┗ 📂 project_01
       ┣ 📂 input
       ┃ ┣ 📂 global_input
       ┃ ┃ ┣ 📜 availability.csv
       ┃ ┃ ┣ 📜 demand_profile.csv
       ┃ ┃ ┣ 📜 ev_parameters.csv
       ┃ ┃ ┣ 📜 power_plant_costs.csv
       ┃ ┃ ┣ 📜 renewables_technical_potential.csv
       ┃ ┃ ┣ 📜 storage_costs.csv
       ┃ ┃ ┣ 📜 storage_inflows.csv
       ┃ ┃ ┗ 📜 technologies.csv
       ┃ ┗ 📂 scenario_01
       ┃   ┣ 📂 industry
       ┃   ┃ ┣ 📜 buses.csv
       ┃   ┃ ┣ 📜 decommission_capacity.csv
       ┃   ┃ ┣ 📜 direct_air_capture.csv
       ┃   ┃ ┣ 📜 fuel_conversion.csv
       ┃   ┃ ┣ 📜 heat_generators.csv
       ┃   ┃ ┣ 📜 heat_links.csv
       ┃   ┃ ┣ 📜 loads.csv
       ┃   ┃ ┣ 📜 storage_capacity.csv
       ┃   ┃ ┗ 📜 storage_energy.csv
       ┃   ┣ 📂 power
       ┃   ┃ ┣ 📜 buses.csv
       ┃   ┃ ┣ 📜 decommission_capacity.csv
       ┃   ┃ ┣ 📜 fuel_suppliers.csv
       ┃   ┃ ┣ 📜 interconnector.csv
       ┃   ┃ ┣ 📜 loads.csv
       ┃   ┃ ┣ 📜 power_generators.csv
       ┃   ┃ ┣ 📜 storage_capacity.csv
       ┃   ┃ ┣ 📜 power_links.csv
       ┃   ┃ ┗ 📜 storage_energy.csv
       ┃   ┣ 📂 transport
       ┃   ┃ ┣ 📜 buses.csv
       ┃   ┃ ┣ 📜 loads.csv
       ┃   ┃ ┣ 📜 pev_chargers.csv
       ┃   ┃ ┗ 📜 pev_storages.csv
       ┃   ┗ 📜 scenario_config.yaml
       ┗ 📂 results (will be created when during the model execution)
```

!!! Tip
    Once you’ve created a skeleton data folder for one scenario, you can simply duplicate the scenario folder and rename the folder name to set up additional scenarios. However, we recommend doing this only after you’ve completed filling in the data for the first one.  

## Step 3: Save your input data systematically

If you change the parameters in Step 1—especially `data_folder_name`—and then run the ` build_skeleton` command in Step 2, the workflow creates a new folder under `data/` with that name (for example, `data/<data_folder_name>/`).

By default, any newly created folders inside `data/` are ignored via `.gitignore`, so they won’t be committed to the main PyPSA-SPICE repository. This helps protect sensitive or proprietary input data from being accidentally pushed.

We also recommend creating a new repository under your private GitHub account (or your team/organization). You can create a new empty repository with the same name as in `<data_folder_name>`, and push your data there by following Git instructions.

```bash title="push your data folder into your private repository"
cd data/<data_folder_name>/
git init -b main
git add <data_folder_name>/ # or: git add .
git commit -m "initialize data folder"
git remote add origin https://github.com/YOUR_ACCOUNT/YOUR_REPO_NAME.git
git push -u origin main
```

Use these Git commands to create your own data folder for working with PyPSA-SPICE. This data repository lives under `data/` in the PyPSA-SPICE working tree while remaining a separate Git repository, so it won’t interfere with the main PyPSA-SPICE Git workflow.
​
You can reuse this pattern whenever you need a new, self-contained data folder, making long-term maintenance straightforward. In the next steps, you’ll edit the input CSVs and configuration files, then commit and push your changes to your private GitHub repository (or repositories) with clear versioning and access control. An example multi–data-folder structure looks like this:

```text title="Structure of multi-data-folder"
📦 data
 ┣ 📂 example
 ┣ 📂 DATA_FOLDER_NAME_REPO1
 ┣ 📂 DATA_FOLDER_NAME_REPO2
 ┗ 📂 DATA_FOLDER_NAME_REPO3
```

## Step 4: fill in the skeleton CSVs

Once a new skeleton folder is created, project-specific CSV templates will be setup. Each CSV will include placeholders marked with `Please fill here`. These need to be completed with relevant data so the model can perform more accurate optimisations.

To help you fill these files:

- check [Global CSV template](global_csv_template.md) for default  file descriptions of country-level data.
- see [Regional CSV template](regional_csv_template.md) for detailed file descriptions of region-level data.

!!! Tip
    By default, the input structure considers large number of technologies represented in the model. If the particular technologies are not needed in your model, it is good practice to remove the input data for these technologies. You can also define your own technologies and customise the model accordingly.  

## Step 5: fill in the scenario assumptions and constraints in `scenario_config.yaml`

Once all the necessary input data is provided, you can adjust model and solver settings in [Model configuration](model-builder-configuration.md) and follow [Model execution](model-builder-execution.md) to understand the model logic and how to run the model.

!!! Tip
    Most of the time, after initializing the first project folder and the first scenario folder (such as the base or reference scenario), you can simply copy this scenario folder and rename it to create new scenarios. In most modeling work, new scenarios are compared against the base or reference scenario, and this approach ensures you don’t need to re-enter all the input data each time.
    Alternatively, you can use `build_skeleton` to create a new scenario folder within an existing project folder if you prefer. In this case, all input data will be empty, providing a clean template structure that requires you to refill each input file manually.
