# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Generate summary output tables for a specified sector and year from the network data.

The generated tables are saved separately as CSV files in a structured directory based
on the sector and year.
"""

import os
from typing import Any

import pandas as pd
from _helpers import add_unit_column
from _post_analysis import OutputTables


def create_folder(sector: str, year: int, results_dir: str) -> str:
    """Create directory for saving output tables based on sector and year.

    Parameters
    ----------
    sector : str
        Sector identifier (e.g., 'power', 'industry').
    year : int
        Year for which the summary is generated.
    results_dir : str
        Directory for storing the results.

    Returns
    -------
    str
        Path to the created directory for saving output tables.
    """
    csv_folder = os.path.join(results_dir, "csvs")
    year_folder = os.path.join(csv_folder, sector, str(year))
    # If csv/{year} folder doesn't exist, make folders
    if os.path.exists(year_folder) is False:
        os.makedirs(year_folder)
    return year_folder  # return filepath to year folder


if __name__ == "__main__":
    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake  # pylint: disable=ungrouped-imports

        snakemake = mock_snakemake("make_summary", sector="p-i-t", years=2025)

    selected_year = int(snakemake.wildcards.years)
    selected_sector = snakemake.wildcards.sector

    save_output_folder = create_folder(
        year=selected_year,
        sector=selected_sector,
        results_dir=snakemake.params.results_dir,
    )
    network_list = []
    network_list.append(snakemake.input.network)
    ot = OutputTables(
        network_list=network_list,
        config=snakemake.config,
        scenario_configs=snakemake.params.scenario_configs,
    )
    # Get all non-default callable methods of output table class
    method_list = [
        func
        for func in dir(ot)
        if callable(getattr(ot, func)) and not func.startswith("__")
    ]
    # Extract only summary methods
    summary_methods = [
        func
        for func in method_list
        if "plot" not in func and "network_dict" not in func
    ]
    # If no industry in wildcard, remove industry methods
    if "i" not in snakemake.wildcards.sector:
        summary_methods = [func for func in summary_methods if "ind_" not in func]
    # If no transport in wildcard, remove transport methods
    if "t" not in snakemake.wildcards.sector:
        summary_methods = [func for func in summary_methods if "tran_" not in func]

    resolution = snakemake.params.scenario_configs["scenario_configs"]["resolution"]
    if resolution["method"] == "nth_hour":
        NTH_HOUR = resolution["stepsize"]
    else:  # clustered method
        NTH_HOUR = 1

    # Extracting output table of each method
    for method in summary_methods:
        # If hourly in method's name, define year parameter
        if "hourly" in method:
            df = getattr(ot, method)(year=selected_year, nth_hour=NTH_HOUR)
        else:
            df = getattr(ot, method)()
        # save table as csv to year folder with name same as the method's name
        if isinstance(df, pd.DataFrame):
            df["unit"] = add_unit_column(method, currency=snakemake.params.currency)
            df.to_csv(f"{save_output_folder}/{method}.csv")
    # Saving summary as
    with open(snakemake.output.summary, "w", encoding="utf-8") as f:
        f.write("all output tables are generated")
