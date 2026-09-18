import os
import pathlib

import sh
import yaml

SCENARIOS = [
    # "pdp_case1_infeasible",
    "pdp_case1_free_emi",
    "pdp_case1_free_emi_res",
    "pdp_case1_free_emi_res_bat",
    # "pdp_case3_infeasible",
    "pdp_case3_free_emi",
    "pdp_case3_free_emi_res",
    "pdp_case3_free_emi_res_bat",
    # "pdp_case4_infeasible",
    "pdp_case4_free_emi",
    "pdp_case4_free_emi_res",
    "pdp_case4_free_emi_res_bat",
    "res_case3_free_emi",
    "res_case3_free_gas_cap",
    "res_case3_free_emi_gas_cap",
    "res_case3_hIPS_free_emi",
    "res_case3_hIPS_free_gas_cap",
    "res_case3_hIPS_free_emi_gas_cap",
    "res_case4_free_emi",
    "res_case4_free_gas_cap",
    "res_case4_free_emi_gas_cap",
    "res_case4_hIPS_free_emi",
    "res_case4_hIPS_free_gas_cap",
    "res_case4_hIPS_free_emi_gas_cap",
]


def print_output(line: str) -> None:
    """Print subprocess output immediately."""
    print(line, end="")


def base_config() -> dict:
    with open(
        os.path.join(pathlib.Path(__file__).parent, "base_config.yaml"),
        encoding="utf-8",
    ) as file:
        data = yaml.safe_load(file)
    return data


if __name__ == "__main__":
    base_data = base_config()

    data_folder = base_data["path_configs"]["data_folder_name"]
    project_folder = base_data["path_configs"]["project_name"]

    for s in SCENARIOS:
        base_data["path_configs"]["input_scenario_name"] = s
        base_data["path_configs"]["output_scenario_name"] = s

        input_scenario_name = base_data["path_configs"]["input_scenario_name"]

        # Save the modified YAML back to the file
        with open(
            os.path.join(pathlib.Path(__file__).parent, "base_config.yaml"),
            "w",
            encoding="utf-8",
        ) as file:
            yaml.dump(base_data, file)

        try:
            print(f"Run scenario {s}")
            result = sh.conda(
                "run",
                "-n",
                "pypsa-spice",
                "snakemake",
                "-j1",
                "-c4",
                "solve_all_networks",
                "-F",
                _out=print_output,
                _err=print_output,
            )
        except sh.ErrorReturnCode as e:
            print(f"Error: {e}")


# run nohup python automation.py &
