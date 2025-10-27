# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Combine all summary output tables from different years into a single CSV file.

Each output table is combined across all specified years and saved in the designated
output directory. For example, if the input years are 2020, 2025, and 2030, the script
will generate combined output tables for each table type in the following format:
<< by_carrier_by_region.csv >> that includes data from all three years.
"""

import os
import warnings
from typing import Any

import pandas as pd
from pandas.errors import InvalidIndexError

warnings.filterwarnings("ignore")

if __name__ == "__main__":

    snakemake: Any = globals().get("snakemake")
    if snakemake is None:
        from _helpers import mock_snakemake

        snakemake = mock_snakemake("combine_summaries", sector="p-i-t", years=2025)

    sector = snakemake.wildcards.sector
    # Extracting list of all summary folder paths
    file_paths = [path.replace("summary.txt", "") for path in snakemake.input.summary]
    # Checking all output csv folders and list out all output files
    files = []
    for path in file_paths:
        files += os.listdir(path)
    # Get a set of output files (no duplicate) and ignore summary.txt
    output_tables = [
        file for file in set(files) if file != "summary.txt" and "hourly" not in file
    ]
    # Extract output paths
    # Output table paths
    ot_path = snakemake.output.combined_summaries.replace("combined_summary.txt", "")
    # Extract year specific output and concat into year_combined outputs csv
    for table in output_tables:
        combined_df = pd.DataFrame()
        for path in file_paths:
            index_cols = (
                [0, 1]
                if any(
                    keyword in table
                    for keyword in [
                        "intercap",
                        "by_carrier_by_region",
                        "by_type_by_carrier",
                        "by_type_by_region",
                    ]
                )
                else 0
            )
            df = pd.read_csv(f"{path}{table}", index_col=index_cols)
            try:
                combined_df = pd.concat([combined_df, df])
                # Export csvs
                combined_df.to_csv(f"{ot_path}{table}")
            except InvalidIndexError:
                print("not combining: ", table)
        with open(snakemake.output.combined_summaries, "w", encoding="utf-8") as f:
            f.write("all year combined output tables are generated")
