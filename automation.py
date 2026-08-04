import sh
import os
import pathlib
import yaml

SCENARIOS = [
    "pdp_high_LOLE10", 
    "pdp_high_LOLE10_free_emi", 
    "pdp_high_LOLE10_shock_gas",
    "res_pdp_high_LOLE10",
    "res_pdp_high_LOLE10_shock_gas",
    "res_pdp_high_LOLE10_high_IPS",
]


def print_output(line: str) -> None:
    """Print subprocess output immediately."""
    print(line, end="")

def base_config() -> dict:
    with open(
        os.path.join(pathlib.Path(__file__).parent, "base_config.yaml"),
        "r",
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
