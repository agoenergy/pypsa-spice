<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: 2020-2025 PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input Data: Define a New Model

This section explains how to set up a new model for a particular country/region using the PyPSA-SPICE model builder. Once your model is step you can run the model as described in [model execution](model-builder-execution.md).

Steps of setting up a new model:

1. Set up the configuration file (`config.yaml`).
2. Run `python script/build_skeleton.py` to create a folder structure and template CSV files for your input data.
3. Fill in the skeleton CSVs with the required data manually or using available resources.

An example structure created by ``build_skeleton`` is displayed below. The following sections will use this example to explain the settings.

## Step 1: Set up the Configuration File

Setting up the config requires defining the scope and resolution of the model. Specifically, defining which countries/regions will be represented in the model and for which year the model will be run. While it is possible to change these after initial model is created, it would require significant effort to add new regions/years in the input CSVs.  

The `config.yaml` contains two parts:

- **Skeleton settings:** for defining the overall structure of the model that is used to create the input data files.
- **Model settings:** for model execution (see [Model Configuration](model-builder-configuration.md)) and turning on/off different conditions.

In the first part (named **Init settings**), following parameters should be defined as shown below.

```yaml title="Init settings in the config.yaml file"
path_configs: #(1)!
  input_dir: data/example/
  results_dir: results/ 
  project_name: project_01
  scenario_name: scenario_01 # (2)!

base_configs:
  countries: ["XY", "YZ"] # (3)!
  regions: 
    XY: ["NR","CE", "SO"] # (4)!
    YZ: ["NR","CE", "SO"] 
  years: [2025, 2030, 2035, 2040, 2045, 2050] # (5)!
  sector: ["p-i-t"] # (6)!
```

1. This section is for configuring directory structure for storing model inputs and results.
2. A custom name you define for the scenario in your model.
3. A list containing **2-letter** country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"}.
4. List of subregions in each country.
5. Modelled years should be provided as a list.
6. Options: [`p`, `p-i`, `p-t`, `p-i-t`], representing power (`p`), industry (`i`), and transport (`t`) sectors.

The final skeleton folder path will follow this structure:: `input_dir`/`countries`/`scenario_name`.

By setting different ``scenario_name`` and country or regional settings in the ``base_configs`` section (see details in [Model Configuration](model-builder-configuration.md)), a new skeleton structure under the same `input_dir` folder will be created.

## Step 2: Build the Skeleton

After modifying the configuration file, run the following command in your terminal.

```bash title="Generating the skeleton folder"
python script/build_skeleton.py
```

This step creates your skeleton folder and files which can be feed with your data.

```text title="Structure of Folder and files created by build skeleton script"
📦 data
 ┣ 📂 example
 ┃  ┗ 📂 project_01
 ┃   ┗ 📂 scenario_01
 ┃     ┣ 📂 industry
 ┃     ┃ ┣ 📜 buses.csv
 ┃     ┃ ┣ 📜 decommission_capacity.csv
 ┃     ┃ ┣ 📜 direct_air_capture.csv
 ┃     ┃ ┣ 📜 fuel_conversion.csv
 ┃     ┃ ┣ 📜 heat_generators.csv
 ┃     ┃ ┣ 📜 heat_links.csv
 ┃     ┃ ┣ 📜 loads.csv
 ┃     ┃ ┣ 📜 storage_capacity.csv
 ┃     ┃ ┗ 📜 storage_energy.csv
 ┃     ┣ 📂 Power
 ┃     ┃ ┣ 📜 buses.csv
 ┃     ┃ ┣ 📜 decommission_capacity.csv
 ┃     ┃ ┣ 📜 fuel_suppliers.csv
 ┃     ┃ ┣ 📜 interconnector.csv
 ┃     ┃ ┣ 📜 loads.csv
 ┃     ┃ ┣ 📜 power_generators.csv
 ┃     ┃ ┣ 📜 storage_capacity.csv
 ┃     ┃ ┣ 📜 power_links.csv
 ┃     ┃ ┗ 📜 storage_energy.csv
 ┃     ┗ 📂 Transport
 ┃       ┣ 📜 buses.csv
 ┃       ┣ 📜 loads.csv
 ┃       ┣ 📜 pev_chargers.csv
 ┃       ┗ 📜 pev_storages.csv
 ┣ 📂 global_csv_templates
 ┃ ┣ 📜 Technologies.csv
 ┃ ┣ 📜 Availability.csv
 ┃ ┣ 📜 Demand_Profile.csv
 ┃ ┣ 📜 PowerPlant_costs.csv
 ┃ ┣ 📜 Renewables_technical_potential.csv
 ┃ ┣ 📜 Storage-Inflows.csv
 ┃ ┗ 📜 Storage_costs.csv
 ┗ 📂 override_component_attrs
```

**Note:** The ``override_component_attrs`` folder supports multi-output for PyPSA components. These files are used internally and don’t need to be edited.

!!! Tip
    You will only need to run `build_skeleton` ones while setting up the model. After this initial run you *should* turn of `build_skeleton` in `config.yaml`. Otherwise, the `build_skeleton` will overwrite your data from in the input files and create new empty files.
    Once you have created a skeleton for one of the scenarios, you can simply copy this folder to create additional scenarios. We recommend to do this after you have filled the data for the first scenario.  

## Step 3: Fill in the Skeleton CSVs

Once a new skeleton folder is created, project-specific CSV templates will be setup. Each CSV will include placeholders marked with `Please fill here`. These need to be completed with relevant data so the model can perform more accurate optimizations.

To help you fill these files:

- See [Regional CSV Template](regional_csv_template.md) for detailed file descriptions.
- Check [Global CSV Template](global_csv_template.md) for default global values. If needed, you can copy values from global templates into the regional ones. However, these values are based on the average assumption of global average, and might not be as accurate as valued obtained from regional sources. Global data is useful as a fallback, but using regional data is strongly recommended for more accurate results.

Once all the necessary input data is provided, adjust model and solver settings in [Model Configuration](model-builder-configuration.md) and follow [Model Execution](model-builder-execution.md) to understand the model logic and how to run the model.

!!! Tip
    By detault, the input structure considers large number of technologies represented in the model. If the particular technologies are not needed in your model, it is good practice to remove the input data for these technologies. You can also define your own technologies and customise the model accordingly.  
