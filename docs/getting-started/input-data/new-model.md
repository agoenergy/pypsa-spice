<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input Data: Define a New Model

This section explains how to set up a new model for a particular country/region using the PyPSA-SPICE model builder. Once your model is step you can run the model as described in [model execution](model-builder-execution.md).

Steps of setting up a new model:

1. Create a folder inside `data` folder. This name of this folder can be the name of your project or scenario (e.g. `test`).
2. Add and set up the configuration file (`config.yaml`). Users are encouraged to regularly compare their own ``config.yaml`` with ``config.default.yaml`` when pulling updates from the remote repository, to ensure compatibility with any new or changed configuration options.
3. Chang the input of `configfile` variable inside `Snakefile` to be the path of your config file. For example:

  ```yaml
    configfile: "data/test/config.yaml"
  ```

4. Run `snakemake -c1 build_skeleton` to create a folder structure and template CSV files for your input data. All input folders and files shall be created inside the `test` folder after this command is executed.
5. Fill in the skeleton CSVs with the required data manually or using available resources.

An example structure created by ``build_skeleton`` is displayed below. The following sections will use this example to explain the settings.

## Step 1-3: Set up the Configuration File

Setting up the config requires defining the scope and resolution of the model. Specifically, defining which countries/regions will be represented in the model and for which year the model will be run. While it is possible to change these after initial model is created, it would require significant effort to add new regions/years in the input CSVs.  

The `config.yaml` contains two parts:

- **Skeleton settings:** for defining the overall structure of the model that is used to create the input data files.
- **Model settings:** for model execution (see [Model Configuration](model-builder-configuration.md)) and turning on/off different conditions.

In the first part (named **Init settings**), following parameters should be defined as shown below.

```yaml title="Init settings in the config.yaml file"
path_configs: #(1)!
  input_dir: data/pypsa-spice-data/
  results_dir: results/ 
  project_name: project_01
  input_scenario_name: scenario_01 # (2)!
  output_scenario_name: scenario_01_tag1 # (3)!

base_configs:
  regions: 
    XY: ["NR","CE", "SO"] # (4)!
    YZ: ["NR","CE", "SO"] 
  years: [2025, 2030, 2035, 2040, 2045, 2050] # (5)!
  sector: ["p-i-t"] # (6)!
```

1. This section is for configuring directory structure for storing model inputs and results.
2. A custom name you define for the input scenario folder in your model.
3. A custom name you define for the output scenario folder to save your model results.
4. List of regions or nodes within each country. This defines the network’s nodal structure. The country list contains **2-letter** country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"}.
5. Modelled years should be provided as a list.
6. Options: [`p`, `p-i`, `p-t`, `p-i-t`], representing power (`p`), industry (`i`), and transport (`t`) sectors.

The final skeleton folder path will follow this structure:: `data`/`input_dir`/`countries`/`input_scenario_name`.

By setting different ``input_scenario_name`` and country or regional settings in the ``base_configs`` section (see details in [Model Configuration](model-builder-configuration.md)), a new skeleton structure under the same `input_dir` folder will be created.

## Step 4: Build the Skeleton

After modifying the configuration file, run the following command in your terminal.

```bash title="Generating the skeleton folder"
snakemake -c1 build_skeleton
```

This step creates your skeleton folder and files which can be feed with your data.

```text title="Structure of Folder and files created by build skeleton script"
📦 data
 ┗ 📂 pypsa-spice-data
    ┣ 📂 project_01
    ┃  ┣ 📂 scenario_01
    ┃  ┃ ┣ 📂 industry
    ┃  ┃ ┃ ┣ 📜 buses.csv
    ┃  ┃ ┃ ┣ 📜 decommission_capacity.csv
    ┃  ┃ ┃ ┣ 📜 direct_air_capture.csv
    ┃  ┃ ┃ ┣ 📜 fuel_conversion.csv
    ┃  ┃ ┃ ┣ 📜 heat_generators.csv
    ┃  ┃ ┃ ┣ 📜 heat_links.csv
    ┃  ┃ ┃ ┣ 📜 loads.csv
    ┃  ┃ ┃ ┣ 📜 storage_capacity.csv
    ┃  ┃ ┃ ┗ 📜 storage_energy.csv
    ┃  ┃ ┣ 📂 power
    ┃  ┃ ┃ ┣ 📜 buses.csv
    ┃  ┃ ┃ ┣ 📜 decommission_capacity.csv
    ┃  ┃ ┃ ┣ 📜 fuel_suppliers.csv
    ┃  ┃ ┃ ┣ 📜 interconnector.csv
    ┃  ┃ ┃ ┣ 📜 loads.csv
    ┃  ┃ ┃ ┣ 📜 power_generators.csv
    ┃  ┃ ┃ ┣ 📜 storage_capacity.csv
    ┃  ┃ ┃ ┣ 📜 power_links.csv
    ┃  ┃ ┃ ┗ 📜 storage_energy.csv
    ┃  ┃ ┗ 📂 transport
    ┃  ┃    ┣ 📜 buses.csv
    ┃  ┃    ┣ 📜 loads.csv
    ┃  ┃    ┣ 📜 pev_chargers.csv
    ┃  ┃    ┗ 📜 pev_storages.csv
    ┃  ┣ 📜 technologies.csv
    ┃  ┣ 📜 availability.csv
    ┃  ┣ 📜 demand_profile.csv
    ┃  ┣ 📜 power_plant_costs.csv
    ┃  ┣ 📜 renewables_technical_potential.csv
    ┃  ┣ 📜 storage_inflows.csv
    ┃  ┗ 📜 storage_costs.csv
    ┗ 📜 config.yaml
```

!!! Tip
    Once you’ve created a skeleton data folder for one scenario, you can simply duplicate it to set up additional scenarios. However, we recommend doing this only after you’ve completed filling in the data for the first one.  

## Step 5: Fill in the Skeleton CSVs

Once a new skeleton folder is created, project-specific CSV templates will be setup. Each CSV will include placeholders marked with `Please fill here`. These need to be completed with relevant data so the model can perform more accurate optimizations.

To help you fill these files:

- See [Regional CSV Template](regional_csv_template.md) for detailed file descriptions.
- Check [Global CSV Template](global_csv_template.md) for default global values. If needed, you can copy values from global templates into the regional ones. However, these values are based on the average assumption of global average, and might not be as accurate as valued obtained from regional sources. Global data is useful as a fallback, but using regional data is strongly recommended for more accurate results.

Once all the necessary input data is provided, adjust model and solver settings in [Model Configuration](model-builder-configuration.md) and follow [Model Execution](model-builder-execution.md) to understand the model logic and how to run the model.

!!! Tip
    By detault, the input structure considers large number of technologies represented in the model. If the particular technologies are not needed in your model, it is good practice to remove the input data for these technologies. You can also define your own technologies and customise the model accordingly.  
