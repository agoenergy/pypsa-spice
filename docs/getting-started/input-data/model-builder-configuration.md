<!--
-*- coding: utf-8 -*-
SPDX-FileCopyrightText: PyPSA-SPICE Developers
SPDX-License-Identifier: GPL-2.0-or-later
-->

# Input data: model builder configuration

PyPSA-SPICE requires two configuration files:

1. **`base_config.yaml`**: Located in the project root directory, this file contains general configurations needed to set up the initial input data structure.

2. **`scenario_config.yaml`**: Located inside each scenario folder, this file contains scenario-specific configurations. It is only used after the input data structure has been created.

To get started, configure `base_config.yaml` first, then run the data setup process. Once complete, you can configure individual scenarios using their respective `scenario_config.yaml` files.

## base_config.yaml

```yaml title="Path configurations"
path_configs: 
  data_folder_name: pypsa-spice-data #(1)!
  project_name: project_01 #(2)!
  input_scenario_name: scenario_01 # (3)!
  output_scenario_name: scenario_01_tag1 # (4)!
```

1. Directory containing all scenario data. Inside this folder, subfolders for `project_name` and `input_scenario_name` will be created.
2. Directory for project-related data, including the `input_scenario_name` folder.
3. Directory for storing scenario input CSVs.
4. Directory for storing scenario output CSVs.

The path for the skeleton folder follows the pattern: `data_folder_name`/`project_name`/input/`input_scenario_name`.
The path for the output folder follows the pattern: `data_folder_name`/`project_name`/results/`output_scenario_name`.

This config file is used for both creating a new model via `snakemake -c1 build_skeleton` (see the section on [Defining a new model](new-model.md)) and used for running different instances of the model.

!!! Tip
    Make sure your snakemake file points to correct config file. To run different scenarios, you just need to change the snakemake file to the corresponding scenario config file.

```yaml title="Base configurations"
base_configs:
  regions: # (1)!
    XY: ["NR","CE", "SO"] 
    YZ: ["NR","CE", "SO"] 
  years: [2025, 2030, 2035, 2040, 2045, 2050] # (2)!
  sector: ["p-i-t"] # (3)!
  currency: USD # (4)!
```

1. List of regions or nodes within each country. This defines the network’s nodal structure. The country list contains **2-letter** country codes according to [ISO 3166](https://www.iso.org/iso-3166-country-codes.html){:target="_blank"}.
2. List of years to be executed in the model builder.
3. List of sectors to include in model run. The power sector (`p`) needs to be included. Other available options are `p-i`, `p-t`, `p-i-t`, representing industry (`i`), and transport (`t`) sectors coupled with the power sector.
4. Currency usd in the model. The default setting is USD (also used in example data). Format shall be in all uppercases, [ISO4217](https://www.iso.org/iso-4217-currency-codes.html){:target="_blank"} format.

## scenario_config.yaml - scenario settings

```yaml title="Scenario configurations"
scenario_configs:
  snapshots: # (1)!
    start: "2025-01-01"
    end: "2026-01-01" # (2)!
    inclusive: "left" # (3)!
  resolution:
    method: "nth_hour" # (4)!
    number_of_days: 3 # (5)!
    stepsize: 25 # (6)!
  interest: 0.05  # (7)!
  remove_threshold: 0.1 # (8)!
```

1. Defines the start and end dates for the model’s time period. Dates are in **"YYYY-MM-DD"** format.
2. The model performs optimisation on a yearly basis, with each modelling year defined as 8,760 hours for hourly resolution. If the selected base year is a leap year, it is recommended to set the end date to `-12-31` of that year.
3. Defines which side of the selected snapshot should be included in the model builder. In the given example, if this parameter is set to `left`, the zeroth hour of the start time snapshot i.e. `2025-01-01 00:00` will be included in the model while the `2026-01-01 00:00` will not be included. We recommend to leave this as is.
4. Determines the method used by the model to deal with the time steps. For testing this reduces the compute time. Available options are `nth_hour` (recommended) and `clustered`. Depending on the selected option, one other parameter in the `resolution` section should be set.
5. If `method = "clustered"`, the number of representative days should be provided to group time steps and form the model builder timestamps. For example, if `number_of_days: 3`, the model builder will only solve 72 hours in the entire year.
6. This is used when `method = "nth_hour"`. In this case, the model builder will run at every **n-th** hour. Typical value to use for this would be 25, so every 25th hour is included in the model. To run the model at hourly resolution (the highest temporal resolution in the model builder), then it needs to be set to 1.
7. Interest rate in decimal form (e.g., 0.05 represents 5%).
8. Removes non expandable assets with a capacity below this threshold (in MW) to avoid numerical issues during optimisation.

## scenario_config.yaml - mandatory constraints

The CO~2~ management is the most important and mandatory constraint in the model. The model allows for two different instruments to decrease CO~2~ emissions: CO~2~ price or CO~2~ constraint. You can choose between each of these but not both.
The variables listed below should be filled out for each country individually.

```yaml title="Constraints for CO<sub>2</sub> management"
co2_management: # (1)!
  XY:
    option: "co2_cap" # (2)!
    value:
      2025: 100
      2030: 90
      2035: 130
      2040: 110
      2045: 100
      2050: 90
  YZ:
    option: "co2_price" # (3)!
    value:
      2025: 1
      2030: 1
      2035: 1
      2040: 1
      2045: 1
      2050: 1
```

1. Please indicate country/region specific CO~2~ management mode: `"co2_cap"` or `"co2_price"`.
2. Goal is to minimise the CO~2~ emissions in each year and target values are given in **mtCO~2~**.
3. Goal is to minimise the emission price (`"co2_price"`) in the system and target values are in **USD/tCO~2~**.

!!! Tip
    If you don't want to use any CO~2~ constraint, default to using `co2_price` with very small value for the years.  

## scenario_config.yaml - custom constraints

The custom constraints section allows you to apply additional rules or limits to the model’s behavior, tailoring it to specific scenario requirements. All custom constraints are listed below in the two countries as an example. These constraints can control various aspects of the model, such as renewable generation share, thermal power plant operation, reserve margins, energy independence, and production limitations. By adjusting these settings, you can implement assumptions or policies. The settings listed below should be configured for each country individually.

!!! Note
    By default these are not included, so if you need a custom constraint, the corresponding part needs to be included in your scenario config file.

```yaml title="Custom constraints"
custom_constraints:
  XY:
    energy_independence: # (1)!
      pe_conv_fraction: # (2)!
        Solar: 1
        Wind: 1
        Geothermal: 1
        Water: 1
      ei_fraction: # (3)!
        2025: 0.3
        2030: 0.4
        2035: 0.5
        2040: 0.6
        2045: 0.7
        2050: 0.8
    production_constraint_fuels: ["Bio", "Bit", "Gas", "Oil"] # (4)!
    reserve_margin: # (5)!
      epsilon_load: 0.1 # (6)!
      epsilon_vre: 0.1 # (7)!
      contingency: 1000 # (8)!
      method: static # (9)!
    res_generation: # (10)!
      math_symbol: "<=" # (11)!
      res_generation_share: # (12)!
        2030: 0.25
        2035: 0.35
        2040: 0.4
        2045: 0.45
        2050: 0.5
    thermal_must_run: # (13)!
      must_run_frac: 0.2 # (14)!
  YZ:
    capacity_factor_constraint: # (15)!
      "SubC": 0.6
      "SupC": 0.6
      "HDAM": 0.4
    production_constraint_fuels: ["Bio", "Bit", "Gas", "Oil"]
    res_generation: 
      math_symbol: "<="
      res_generation_share:
        2030: 0.1
        2035: 0.2
        2040: 0.3
        2045: 0.3
        2050: 0.3
```

1. The model adds constraints to ensure energy independence. This indicates how much of energy needs are met without relying on imports (by producing enough energy domestically). You can refer to Constraint - Energy Independence for more information.<br>To deactivate it, you can exclude them in the `custom_constraints`, and the model will identify it as deactivated.
2. Primary energy conversion factor (dimensionless) is used to convert electricity generation to `primary energy` to make renewables comparable to fossil at primary energy level. Different definitions can be used to arrive at the value of these.
3. **Minimum** energy independence fraction defined as: $$ \frac{\textit{locally produced energy}}{\textit{locally produced energy + imported energy}} $$ For details see `Energy independence constraint` below.
4. Maximum production limit of certain fuels can be defined here. Maximum values for these fuels are defined in `Power/fuel_supplies.csv`.<br>To deactivate it, you can remove them in the `custom_constraints`, and the model will identify it as dectivated.
5. The model adds reserve margin constraints based on `reserve_parameters`. See Constraint - Reserve Margin for more information.<br>To deactivate it, you can exclude them in the `custom_constraints`, and the model will identify it as dectivated.
6. Fraction of load considered as reserve.
7. Contribution of Variable Renewable Energy (VRE) to the reserve.
8. Extra contingency in MW. It is under `reserve_margin`. This is usually taken as a the size of largest individual power plants or defined by country specific regulations. 
9. Options: `static` (no VRE) or `dynamic` (includes VRE). See reserve margin definition below.
10. The model adds a constraint on renewable generation as a fraction of total electricity demand.
11. Defines the type of renewable constraint. If set to `<=` it means the fraction of renewable generation from the total electricity demand should be **less than or equal to** the given values or **greater than or equal to** if value set to `=>`
12. Fraction of renewable generation to the total electricity demand for each year.
13. The model forces combined thermal power plants to have minimum generation level as a fraction of load.
14. Fraction of thermal generation to the total electricity demand per snapshot providing the baseload.
15. Maximum capacity factor of certain technologies can be defined here.<br>To deactivate it, you can exclude them in the `custom_constraints`, and the model will identify it as dectivated.

In the following two sub-sections, we provide more information about the definition of energy independence and reserve margin.

### Energy independence: mathematical formulation

This constraint forces the model to keep the ratio of locally produced power to the sum of locally produced power and imported power to be more than the minimum energy independence factor:

```math
\frac{loc}{imp + loc} \ge \phi \qquad \Longrightarrow \qquad (1 -\phi) \cdot loc >= \phi \cdot imp
```

Where

| parameter  | description                            | mathematical formulation |
| ---------- | -------------------------------------- | ----- |
| $`imp`$  | Imported power:  Generation from the theoretical import fuel-based generators | $`\sum_{f \in F} \sum_{t \in T}{G_t^{TGEN-f-import}}`$ |
| $`loc`$  | Local power generation: Fuel-based generations from local resources + renewable generations x primary energy conversion factor | $`\sum_{f \in F} \sum_{t \in T} {G_t^{TGEN-f-local}} + \alpha_{res}*\sum_{t}G_{res}`$ |
| $`\phi`$ | minimum energy independence factor  | |
| $`\alpha_{res}`$ | Primary energy conversion factor used for renewable sources for electricity generation. This value can be 0-1 for renewables, or larger than 1 for other generation sources depending on the energy policy in the country. | |

### Reserve margin: mathematical formulation

The reserve margin constraint in PyPSA-SPICE is modeled similarly to the [GenX](https://genxproject.github.io/GenX.jl/stable/Model_Reference/core/#Operational-Reserves){:target="_blank"} approach.

```math
\sum_g r_{g,t} \ge \epsilon^{load} \cdot \sum_n d_{n,t} + \epsilon^{vres} \cdot \sum_{g \in \mathcal{G}^{VRES}} \bar g_{g,t} \cdot G_g + {contingency} \qquad \forall t
```

Where

| parameter                   | description                                      |
| ------------------------ | ------------------------------------------------ |
| $`r_{g,t}`$            | Reserve margin of generator `g` at time `t` |
| $`\epsilon^{load}`$ | Fraction of load considered for reserve |
| $`d_{n,t}`$            | Demand at node `n` and time `t` |
| $`\epsilon^{vres}`$ | Fraction of renewable energy for reserve |
| $`\bar g_{g,t}`$       | Forecasted capacity factor for renewable energy of generator `g` at time `t` |
| $`G_g`$                | Capacity of generator `g` |
| $`\mathcal{G}^{VRES}`$ | Set of renewable generators in the system |
| $`contingency`$ | Fixed contingency |

See [Linopy example](https://github.com/PyPSA/pypsa-eur/blob/7ac983e5b31bcaf3ae667ceec4fc9d5d91c18046/scripts/solve_network.py#L387-L454){:target="_blank"} of the reserve constraint implementation for more details.

## scenario_config.yaml - solver settings

Solving the optimisation model builder requires a good solver to boost the performance. PyPSA-SPICE supports solvers such as `gurobi`, `cplex`, and `highs`. A comparison of solver performance is available in [solver benchmarking results](https://openenergybenchmark.org/dashboard/main-results){:target="_blank"}.

```yaml title="Solver configurations"
solving:
  solver:
    name: highs #(1)!
    options: highs-default #(2)!
  oetc: # (3)!
    activate: false
    name: test-agora-job
    authentication_server_url: http://34.34.8.15:5050
    orchestrator_server_url: http://34.34.8.15:5000
    cpu_cores: 4
    disk_space_gb: 20
    delete_worker_on_error: false
  solver_options: #(4)!
    default: {}
    cbc-default:
      threads: 8 #(5)!
      cuts: 0 #(6)!
      maxsol: 1 #(7)!
      ratio: 0.1 #(8)!
      presolve: 1 #(9)!
      time_limit: 3600 #(10)!
    gurobi-default:
      threads: 8 #(11)!
      method: 2 #(12)!
      crossover: 0 #(13)!
      BarConvTol: 1.e-5 #(14)!
      AggFill: 0 #(15)!
      PreDual: 0 #(16)!
      GURO_PAR_BARDENSETHRESH: 200 #(17)!
    gurobi-numeric-focus:
      name: gurobi
      NumericFocus: 3 #(18)!
      method: 2 #(19)!
      crossover: 0 #(20)!
      BarHomogeneous: 1 #(21)!
      BarConvTol: 1.e-5 #(22)!
      FeasibilityTol: 1.e-4 #(23)!
      OptimalityTol: 1.e-4 #(24)!
      ObjScale: -0.5 #(25)!
      threads: 8 #(26)!
      Seed: 123 #(27)!
    cplex-default:
      threads: 4 #(28)!
      lpmethod: 4 #(29)!
      solutiontype: 2 #(30)!
      barrier_convergetol: 1.e-5 #(31)!
      feasopt_tolerance: 1.e-6
    highs-default: #(32)!
      threads: 4 #(33)!
      solver: "ipm"
      run_crossover: "off"
      small_matrix_value: 1e-6 #(34)!
      large_matrix_value: 1e9 #(35)!
      primal_feasibility_tolerance: 1e-5 #(36)!
      dual_feasibility_tolerance: 1e-5 #(37)!
      ipm_optimality_tolerance: 1e-4 #(38)!
      parallel: "on"
      random_seed: 123 #(39)!
    highs-simplex:
      solver: "simplex"
      parallel: "on"
      primal_feasibility_tolerance: 1e-5 #(40)!
      dual_feasibility_tolerance: 1e-5 #(41)!
      random_seed: 123 #(42)!
```

1. Compatible solvers are `gurobi`, `cplex`, `cbc`, and `highs`.
2. Depending on the selected solver, specify one of the corresponding options: `cbc-default`, `gurobi-default`, `gurobi-numeric-focus`, `cplex-default`, `highs-default`, or `highs-simplex`.
3. `oetc` is a colud computing service provided by [Open Energy Trainsition](https://openenergytransition.org/){:target="_blank"} Organisation. To access and activate the service. Please reach out to us for more details.
4. Solver options can be adjusted in the following list. If no value is provided in the `solver/option` section, the default `solver_options`, which is empty, will be considered.
5. Number of CPU threads to be used by the solver for parallel computation to speed up solving time.
6. Cutting planes are typically used to tighten the problem and improve performance (usually beneficial in MIP problems). By disabling them solution will be obtained faster but potentially at the cost of optimality.
7. Limits the solver to finding only one solution. The solver will stop once it finds a feasible solution (instead of finding all solutions).
8. Specifies that the solver should stop if it finds a solution within 10% of the best possible bound. This is useful for faster solutions when absolute optimality is not required.
9. Enabling `presolve` simplifies the problem before starting the optimisation to make it faster and more stable.
10. Sets a time limit for the solver (here 3600 seconds = 1 hour). The solver will stop if it exceeds this limit.
11. Number of CPU threads to be used by the solver for parallel computation to speed up solving time.
12. Algorithm used to solve continuous models or the initial root relaxation of a MIP model. `-1` chooses the algorithm automatically and other options are explained in [Gurobi documentation](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#method){:target="_blank"} for more information.
13. Barrier crossover strategy. `-1` chooses strategy automatically and `0` disables crossover, which will speed up the solution process but might reduce solution quality. Other options are explained in [Gurobi documentation](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#crossover){:target="_blank"}.
14. The barrier solver terminates when the relative difference between the primal and dual objective values is less than the specified tolerance. This parameter is in `[0, 1]` range and the default value is `1e-8`.
15. A parameter that controls the amount of **fill** allowed during the aggregation phase of presolve. `AggFill` determines how aggressively Gurobi merges constraints during aggregation. Higher values can potentially lead to more simplification but may also introduce numerical instability. `-1` chooses aggregation fill automatically and `0` disables it.
16. Determines whether the solver should dualize the problem during the presolve phase. Depending on the structure of the model, solving the dual can reduce overall solution time. The default setting (`-1`) decides about it automatically. Setting `0` forbids presolve from forming the dual, while setting `1` forces it to take the dual.
17. Sets the threshold for determining when a column in the constraint matrix is considered **dense** during barrier optimisation. When the constraint matrix is dense, it means its non-zero elements are more than the `GURO_PAR_BARDENSETHRESH` value.
18. With higher values, the model will spend more time checking the numerical accuracy of intermediate results.
19. Algorithm used to solve continuous models or the initial root relaxation of a MIP model. `-1` chooses the algorithm automatically and other options are explained in [Gurobi documentation](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#method){:target="_blank"} for more information.
20. Barrier crossover strategy. `-1` chooses strategy automatically and `0` disables crossover, which will speed up the solution process but might reduce solution quality. Other options are explained in [Gurobi documentation](https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#crossover){:target="_blank"}.
21. Setting the parameter to `0` turns it off, and setting it to `1` forces it on. The homogeneous algorithm is useful for recognizing infeasibility or unboundedness and is a bit slower than the default algorithm.
22. The barrier solver terminates when the relative difference between the primal and dual objective values is less than the specified tolerance. This parameter is in `[0, 1]` range and the default value is `1e-8`.
23. All constraints must be satisfied to a tolerance of `FeasibilityTol`. This parameter is in `[1e-9, 1e-2]` range and the default value is `1e-6`.
24. This parameter defines how close the solution needs to be to the best possible answer before the solver stops. It is in `[1e-9, 1e-2]` range and the default value is `1e-6`.
25. **Positive values:** divides the objective by the specified value to avoid numerical issues that may result from very large or very small objective coefficients. **Negative values:** uses the maximum coefficient **to the specified power** as the scaling (so ObjScale=-0.5 would scale by the square root of the largest objective coefficient). **Default:** `0` which means the model decides on the scaling automatically.
26. Number of CPU threads to be used by the solver for parallel computation to speed up solving time.
27. Fixed `Seed` values must be set if you want to get the same results (i.e., reproducibility of the optimisation process). The value does not matter.
28. Number of CPU threads to be used by the solver for parallel computation to speed up solving time.
29. This parameter changes the algorithm and accepts an integer from `0` to `6`, where `0` denotes automatic choice of the algorithm, `1` is for primal simplex, `2` is for dual simplex, and `4` is for barrier.
30. Crossover can be turned off with `solutiontype=2` that instructs CPLEX not to seek a basic solution. This can be useful for a quick insight of the approx. optimal solution, if crossover takes long time.
31. Sets the tolerance on complementarity for convergence. Values can be qny positive number greater than or equal to `1e-12`; default: `1e-8`.
32. Please visit [HiGHS documentation](https://ergo-code.github.io/HiGHS/dev/options/definitions/){:target="_blank"} for complete list of options.
33. Number of CPU threads to be used by the solver for parallel computation to speed up solving time.
34. Values less than or equal to this will be treated as zero.
35. Values greater than or equal to this will be treated as infinite.
36. Range: `[1e-10, inf]`, default: `1e-07`.
37. Range: `[1e-10, inf]`, default: `1e-07`.
38. Range: `[1e-10, inf]`, default: `1e-07`.
39. Fixed `random_seed` values must be set if you want to get the same results (i.e., reproducibility of the optimisation process). The value does not matter.
40. Range: `[1e-10, inf]`, default: `1e-07`.
41. Range: `[1e-12, inf]`, default: `1e-08`.
42. Fixed `random_seed` values must be set if you want to get the same results (i.e., reproducibility of the optimisation process). The value does not matter.
