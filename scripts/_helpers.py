# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""Helper functions for pypsa-spice."""

import glob
import logging
import os
import urllib
from itertools import cycle
from os import PathLike
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import snakemake as sm
import yaml
from progressbar import ProgressBar
from snakemake.api import Workflow
from snakemake.common import SNAKEFILE_CHOICES
from snakemake.script import Snakemake
from snakemake.settings.types import (
    ConfigSettings,
    DAGSettings,
    ResourceSettings,
    StorageSettings,
    WorkflowSettings,
)

FilePath = Union[str, PathLike[str]]

agora_style = {
    "line_width": 2,
    "line_style": "dash",
    "font_family": "DejaVu Sans",
    "font_size": 14,
    "axes_facecolor": "#E3E4EA",
    "axes_edgecolor": "white",
    "axes_labelcolor": "black",
    "axes_labelweight": "bold",
    "xtick_color": "#000000",
    "ytick_color": "#000000",
    "legend_facecolor": "#E3E4EA",
    "figure_facecolor": "#E3E4EA",
    "figure_edgecolor": "#E3E4EA",
    "figure_figsize": [800, 450],
}

ag_cp = {
    "GAST": "#FF3300",
    "CCGT": "#FF3300",
    "OCGT": "#FF3300",
    "GASB": "#FF3300",
    "GAS": "#FF3300",
    "GAS_SUPPLY": "#FF3300",
    "HRDT_SubC": "#000000",
    "PHOT": "#FFD744",
    "RTPV": "#F1C40F",
    "WTON": "#1F82C0",
    "BIOT": "#467850",
    "OILT": "#1D565C",
    "GEOT": "#D557A6",
    "CSP": "#EF7F03",
    "WTOF": "#006099",
    "BATS": "#00ced1",
    "BATS_DISCHARGE": "#00ced1",
    "BATS_STORE": "#00ced1",
    "HHBS": "#00ced1",
    "HHBS_DISCHARGE": "#00ced1",
    "HHBS_STORE": "#00ced1",
    "HDAM": "#97BAD5",
    "HROR": "#97BAD5",
    "EVST": "#00ced1",
    "HPHS": "#008080",
    "HPHS_DISCHARGE": "#008080",
    "HPHS_STORE": "#008080",
    "EDLB": "#97BAD5",
    "EERH": "#D557A6",
    "EHPP": "#5E419B",
    "EIDT": "#F1C40F",
    "DAC": "#733e88",
    "EDLH": "#EF7F03",
    "INLHSTOR_DISCHARGE": "#8493C2",
    "INLHSTOR_STORE": "#8493C2",
    "INLHSTOR": "#8493C2",
    "IND_BOILER": "#5F605C",
    "OCHT": "#733e88",
    "Import": "#be0032",
    "SupC": "#B69E8C",
    "CHP": "#467850",
    "SubC": "#423429",
    "WSTT": "#0BDA51",
    "WST": "#0BDA51",
    "LSLO": "#0BDA51",
    "NUCL": "#5F605C",
    "SWHT": "#FFD744",
    "fossils": "#9ABBCA",
    "renewables": "#98B721",
    "Bio": "#408B2E",
    "BIO_SUPPLY": "#467850",
    "Gas": "#9ABBCA",
    "Hrdc": "#0C0C0C",
    "HRDC_SUPPLY": "#0C0C0C",
    "Lig": "#5E4745",
    "Bit": "#5E4745",
    "Lng": "#8493C2",
    "Oil": "#006057",
    "Waste": "#0BDA51",
    "Hyd": "#d2dde4",
    "Uranium": "#5F605C",
    "PUB": "#000000",
    "PRV": "#8393BE",
    "ITCN": "#733e88",
}

ag_st = {
    "#000000",
    "#8393BE",
    "#733e88",
    "#d05094",
    "#1e84b3",
    "#48a7ae",
    "#ad86b0",
    "#617494",
    "#64b9e4",
}


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined by
    either (in priority order): "snakemake.log.python", "snakemake.log[0]" or
    "logs/{rulename}.log". Additional keywords from logging.basicConfig are accepted via
    the snakemake configuration file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and
        file.
    """
    kwargs = snakemake.config.get("logging", {})
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ]
            }
        )
    logging.basicConfig(**kwargs)


def retrieve_snakemake_keys(snakemake):
    """Retrieve the most commonly used snakemake keys."""
    return (
        snakemake.input,
        snakemake.config,
        snakemake.wildcards,
        snakemake.log,
        snakemake.output,
    )


def aggregate_p_nom(n):
    """Aggregate nominal power and average load per carrier."""
    return pd.concat(
        [
            n.generators.groupby("carrier").p_nom_opt.sum(),
            n.storage_units.groupby("carrier").p_nom_opt.sum(),
            n.links.groupby("carrier").p_nom_opt.sum(),
            n.loads_t.p.groupby(n.loads.carrier, axis=1).sum().mean(),
        ]
    )


def aggregate_p(n):
    """Aggregate total power per carrier."""
    return pd.concat(
        [
            n.generators_t.p.sum().groupby(n.generators.carrier).sum(),
            n.storage_units_t.p.sum().groupby(n.storage_units.carrier).sum(),
            n.stores_t.p.sum().groupby(n.stores.carrier).sum(),
            -n.loads_t.p.sum().groupby(n.loads.carrier).sum(),
        ]
    )


def aggregate_e_nom(n):
    """Aggregate nominal energy capacity per carrier."""
    return pd.concat(
        [
            (n.storage_units["p_nom_opt"] * n.storage_units["max_hours"])
            .groupby(n.storage_units["carrier"])
            .sum(),
            n.stores["e_nom_opt"].groupby(n.stores.carrier).sum(),
        ]
    )


def aggregate_p_curtailed(n):
    """Aggregate curtailed power per carrier."""
    return pd.concat(
        [
            (
                (
                    n.generators_t.p_max_pu.sum().multiply(n.generators.p_nom_opt)
                    - n.generators_t.p.sum()
                )
                .groupby(n.generators.carrier)
                .sum()
            ),
            (
                (n.storage_units_t.inflow.sum() - n.storage_units_t.p.sum())
                .groupby(n.storage_units.carrier)
                .sum()
            ),
        ]
    )


def aggregate_costs(n, flatten=False, opts=None, existing_only=False):
    """Aggregate capital and marginal costs per carrier."""
    components = {
        "Link": ("p_nom", "p0"),
        "Generator": ("p_nom", "p"),
        "StorageUnit": ("p_nom", "p"),
        "Store": ("e_nom", "p"),
        "Line": ("s_nom", None),
        "Transformer": ("s_nom", None),
    }

    costs = {}
    for c, (p_nom, p_attr) in zip(
        n.iterate_components(components.keys(), skip_empty=False), components.values()
    ):
        if c.df.empty:
            continue
        if not existing_only:
            p_nom += "_opt"
        costs[(c.list_name, "capital")] = (
            (c.df[p_nom] * c.df.capital_cost).groupby(c.df.carrier).sum()
        )
        if p_attr is not None:
            p = c.pnl[p_attr].sum()
            if c.name == "StorageUnit":
                p = p.loc[p > 0]
            costs[(c.list_name, "marginal")] = (
                (p * c.df.marginal_cost).groupby(c.df.carrier).sum()
            )
    costs = pd.concat(costs)

    if flatten:
        assert opts is not None
        conv_techs = opts["conv_techs"]

        costs = costs.reset_index(level=0, drop=True)
        costs = costs["capital"].add(
            costs["marginal"].rename({t: t + " marginal" for t in conv_techs}),
            fill_value=0.0,
        )

    return costs


def progress_retrieve(url, file):
    """Download a file from a URL with a progress bar."""
    pbar = ProgressBar(0, 100)

    def dl_progress(count, block_size, total_size):
        pbar.update(int(count * block_size * 100 / total_size))

    urllib.request.urlretrieve(url, file, reporthook=dl_progress)


def get_reduced_hours(step_size: int, load: pd.Series):
    """Aggregate 8760 hours into xth hours based on stepSize and include peak hour.

    Parameters
    ----------
    step_size : int
        stepsize to aggregate
    load : pd.Series
        total electricity load

    Returns
    -------
    reduced_hours[list]:
        List of reduced hours
    weights:
        weight per reduced snapshot
    """
    hour_list = np.arange(0, 8760, 1)
    xth_hour = hour_list[::step_size]
    peak_hour = load.argmax()

    if peak_hour not in xth_hour:
        print(f"Adding peak hour: {peak_hour}th")
        xth_hour = np.append(xth_hour, peak_hour)
    reduced_hours = np.sort(xth_hour)
    weightings = [step_size if x != peak_hour else 1 for x in reduced_hours]
    weightings[-1] = 8760 - (step_size * (len(reduced_hours) - 2) + 1)
    weights = pd.DataFrame(weightings, index=list(range(1, len(reduced_hours) + 1, 1)))
    weights.rename(columns={0: "objective"}, inplace=True)
    weights["objective"] = step_size  # having uniform weighting
    weights["generators"] = weights["objective"]
    weights["stores"] = weights["objective"]

    return reduced_hours, weights


def mock_snakemake(
    rulename,
    root_dir=None,
    configfiles=None,
    **wildcards,
):
    """Create a snakemake.script.Snakemake object for a given rule in the Snakefile.

    This function is expected to be executed from the 'scripts'-directory of the
    snakemake project.

    If a rule has wildcards, you have to specify them in **wildcards.

    Parameters
    ----------
    rulename : str
        name of the rule for which the snakemake object should be generated
    root_dir : str, optional
        path to the root directory of the snakemake project, by default None
    configfiles : list, optional
        list of configfiles to be used to update the config, by default None

    **wildcards:
        keyword arguments fixing the wildcards. Only necessary if wildcards are needed.
    """
    script_dir = Path(__file__).parent.resolve()
    if root_dir is None:
        root_dir = script_dir.parent
    else:
        root_dir = Path(root_dir).resolve()

    workdir: Path | None = None
    user_in_script_dir = Path.cwd().resolve() == script_dir
    if user_in_script_dir:
        os.chdir(root_dir)
    elif Path.cwd().resolve() != root_dir:
        print(
            "Not in scripts or root directory, will assume this is a separate workdir"
        )
        workdir = Path.cwd()

    try:
        for p in SNAKEFILE_CHOICES:
            p = root_dir / p
            if os.path.exists(p):
                snakefile = p
                break
        if configfiles is None:
            configfiles = []
        elif isinstance(configfiles, str):
            configfiles = [configfiles]

        resource_settings = ResourceSettings()
        config_settings = ConfigSettings(configfiles=map(Path, configfiles))
        workflow_settings = WorkflowSettings()
        storage_settings = StorageSettings()
        dag_settings = DAGSettings(rerun_triggers=[])
        workflow = Workflow(
            config_settings,
            resource_settings,
            workflow_settings,
            storage_settings,
            dag_settings,
            storage_provider_settings={},
            overwrite_workdir=workdir,
        )
        workflow.include(snakefile)

        if configfiles:
            for f in configfiles:
                if not os.path.exists(f):
                    raise FileNotFoundError(f"Config file {f} does not exist.")
                workflow.configfile(f)

        workflow.global_resources = {}
        rule = workflow.get_rule(rulename)
        dag = sm.dag.DAG(workflow, rules=[rule])
        wc = dict(wildcards)
        job = sm.jobs.Job(rule, dag, wc)

        def make_accessable(*ios):
            for io in ios:
                for i, _ in enumerate(io):
                    io[i] = os.path.abspath(io[i])

        make_accessable(job.input, job.output, job.log)
        snakemake = Snakemake(
            job.input,
            job.output,
            job.params,
            job.wildcards,
            job.threads,
            job.resources,
            job.log,
            job.dag.workflow.config,
            job.rule.name,
            None,
        )
        # create log and output dir if not existent
        for path in list(snakemake.log) + list(snakemake.output):
            Path(path).parent.mkdir(parents=True, exist_ok=True)

    finally:
        if user_in_script_dir:
            os.chdir(script_dir)
    return snakemake


def scaling_conversion(
    input_df: pd.DataFrame, scaling_number: int, decimals: int
) -> pd.DataFrame:
    """Change scale of value column in a DataFrame and return the DataFrame.

    Parameters
    ----------
    input_df : pd.DataFrame
        Input DataFrame.
    scaling_number : int
        Scaling number used for conversion.
    decimals : int
        Round up decimal numbers after conversion.

    Returns
    -------
    pd.DataFrame
        Output DataFrame after scaling conversion.
    """
    input_df = input_df.reset_index()
    input_df = input_df.loc[:, ~input_df.columns.duplicated(keep="last")]
    col_list = input_df.columns.to_list()
    col_list.remove("value")
    input_df = input_df.set_index(col_list)
    input_df = (input_df / scaling_number).round(decimals)
    input_df = input_df.reset_index()

    return input_df.set_index(input_df.columns[0])


def get_buses(bus_df: pd.DataFrame) -> pd.DataFrame:
    """Get buses with carriers.

    Parameters
    ----------
    bus_df : pd.DataFrame
        table for all buses

    Returns
    -------
    pd.DataFrame
        table of buses with assigned carriers
    """
    bus_df = bus_df[["bus", "carrier", "node", "country"]].set_index("bus")
    return bus_df


def get_time_series_demands(
    load_df: pd.DataFrame, profile_dir: FilePath, year: int
) -> pd.DataFrame:
    """Convert total load per buses into timeseries for 8760 h based on profile pattern.

    Parameters
    ----------
    load_df : FilePath
        table of load to convert
    profile_dir : FilePath
        directory to database of profile patterns per bus type
    year : int
        load year to convert

    Returns
    -------
    pd.DataFrame
        DataFrame of load per buses converted into time series for 8760 hour based on
        profile pattern

    Raises
    ------
    ValueError
        If load_id has duplicate
    """
    # filter load table to only corresponding load year
    profile_df = pd.read_csv(profile_dir)
    load_df = load_df[load_df.year == year]

    # Start matching profile
    dfs = []
    load_df = load_df.to_dict("records")
    for row in load_df:
        # filter country first
        profile_country_df = profile_df[profile_df.country == row["country"]].drop(
            "country", axis=1
        )
        # If the demand type repeats in profile_df more than 1 row
        # meaning the demand type should be matched by regional specific profile.
        if list(profile_country_df.profile_type).count(row["profile_type"]) > 1:
            result = profile_country_df[
                (profile_country_df["profile_type"].isin([row["profile_type"]]))
                & (profile_country_df["node"].isin([row["node"]]))
            ]
        # Else all region will have similar profile.
        else:
            result = profile_country_df[
                profile_country_df["profile_type"].isin([row["profile_type"]])
            ]
        result = result.iloc[:, 2:].astype(float) * float(row["total_load__mwh"])
        if result.empty:
            print(f'{row["profile_type"]} is not in profile database')
        result["bus"] = row["bus"]
        result["node"] = row["node"]
        result["carrier"] = row["carrier"]
        result["load_id"] = row["name"]
        result["country"] = row["country"]
        if len(result) > 1:
            raise ValueError(
                f'{result["load_id"]} has duplicate. Please check input data again'
            )
        dfs.append(result)
    demand_df = pd.concat(dfs)
    demand_df.set_index(["load_id", "country", "bus", "carrier", "node"], inplace=True)

    return demand_df


def get_plant_availabilities(
    pp_df: pd.DataFrame, avail_dir: FilePath, arch_dir: FilePath
) -> pd.DataFrame:
    """Match each generator/storage unit to corresponding availability.

    Note: Matching is done based on location and technology.

    Parameters
    ----------
    pp_df : pd.DataFrame
        DataFrame of all power plants
    avail_dir : FilePath
        directory to database of availability per location and technology
    arch_dir : FilePath
        directory to database of techno-economic params per technology

    Returns
    -------
    pd.DataFrame
        DataFrame of all generators and storage units with assigned timeseries
        availability for 8760 snapshots

    Raises
    ------
    ValueError
        If plant name has duplicate
    """
    # Looping through each row in pp_df to find technology info and match to
    # corresponding availability
    dfs = []
    avail_df = pd.read_csv(avail_dir)
    arch_df = pd.read_csv(arch_dir)
    pp_df = pp_df.to_dict("records")
    for row in pp_df:
        # filter country first
        avail_country_df = avail_df[avail_df.country == row["country"]].drop(
            "country", axis=1
        )
        arch_country_df = arch_df[arch_df.country == row["country"]].drop(
            "country", axis=1
        )
        # If the technology type repeats in avail_df in more than 1 row
        # meaning the technology should be matched by regional specific profile.
        # e.g. renewable technologies like solar, wind, ror hydro, etc.
        if list(avail_country_df.technology).count(row["type"]) > 1:
            result = avail_country_df[
                (avail_country_df["technology"].isin([row["type"]]))
                & (avail_country_df["node"].isin([row["node"]]))
            ]
        # Else if the technology type appears in avail_df in only 1 row
        # all region will have similar availability profile for given technology.
        elif list(avail_country_df.technology).count(row["type"]) == 1:
            result = avail_country_df[
                (avail_country_df["technology"].isin([row["type"]]))
            ]
        # If technology type doesn't appear in avail_df
        # take availability from technology database.
        else:
            p_max_pu = arch_country_df[
                (arch_country_df["technology"].isin([row["type"]]))
                & (arch_country_df["carrier"].isin([row["carrier"]]))
            ]["p_max_pu"]
            # creating result to match format of other availability dataframes
            result = pd.DataFrame(
                np.repeat(p_max_pu.values, 8762, axis=0),
            ).T
            result.columns = ["node", "technology"] + [str(i) for i in range(0, 8760)]
            if result.empty:
                result = avail_df[(avail_df["technology"].isin(["Con10"]))]

        if len(result) > 1:
            raise ValueError(
                f'{row["name"]} has duplicate. Please check input data again'
            )
        result = result.iloc[:, 2:]
        result["plant"] = row["name"]
        dfs.append(result.set_index("plant"))
    availability_df = pd.concat(dfs)
    availability_df = availability_df.astype(float)
    return availability_df


def get_link_availabilities(
    link_df: pd.DataFrame, avail_dir: FilePath, arch_dir: FilePath
) -> pd.DataFrame:
    """Match each converters to corresponding availability.

    Note: Matching is done based on location and technology.

    Parameters
    ----------
    link_df : pd.DataFrame
        DataFrame of all links
    avail_dir : FilePath
        directory to database of availability per location and technology
    arch_dir : FilePath
        directory to database of techno-economic params per technology

    Returns
    -------
    pd.DataFrame
        DataFrame of all links with assigned time-series availability for 8760 snapshots

    Raises
    ------
    ValueError
        If link name has duplicate
    """
    # Looping through each row in pp_df to find technology info and match to
    # corresponding availability
    dfs = []
    avail_df = pd.read_csv(avail_dir)
    arch_df = pd.read_csv(arch_dir)
    link_df = link_df.to_dict("records")
    for row in link_df:
        # filter country first
        avail_country_df = avail_df[avail_df.country == row["country"]].drop(
            "country", axis=1
        )
        arch_country_df = arch_df[arch_df.country == row["country"]].drop(
            "country", axis=1
        )
        # If the technology type repeats in avail_df in more than 1 row
        # meaning the technology should be matched by regional specific profile.
        # e.g. renewable technologies like solar, wind, ror hydro, etc.
        if list(avail_country_df.technology).count(row["type"]) > 1:
            avail_country_df = avail_country_df[
                (avail_country_df["technology"].isin([row["type"]]))
                & (avail_country_df["node"].isin([row["node"]]))
            ]
        # Else if the technology type appears in avail_df in only 1 row
        # all region will have similar availability profile for given technology.
        elif list(avail_country_df.technology).count(row["type"]) == 1:
            result = avail_country_df[
                (avail_country_df["technology"].isin([row["type"]]))
            ]
        # If technology type doesn't appear in avail_df
        # get availability from technology database.
        else:
            p_max_pu = arch_country_df[
                (arch_country_df["technology"].isin([row["type"]]))
                & (arch_country_df["carrier"].isin([row["carrier"]]))
            ]["p_max_pu"]
            # creating result to match format of other availability dataframes
            result = pd.DataFrame(
                np.repeat(p_max_pu.values, 8762, axis=0),
            ).T
            # give column names to match with other availability dataframes
            result.columns = ["node", "technology"] + [str(i) for i in range(0, 8760)]
            if result.empty:
                result = avail_df[(avail_df["technology"].isin(["Con10"]))]

        if len(result) > 1:
            raise ValueError(
                f'{row["name"]} has duplicate. Please check input data again'
            )
        result = result.iloc[:, 2:]
        result["plant"] = row["link"]
        dfs.append(result.set_index("plant"))
    availability_df = pd.concat(dfs)
    availability_df = availability_df.astype(float)
    return availability_df


def get_store_min_availabilities(
    stores_df: pd.DataFrame, avail_dir: FilePath
) -> pd.DataFrame:
    """Match each storage unit to corresponding minimum availability.

    Parameters
    ----------
    stores_df : pd.DataFrame
        Table of all storage units
    avail_dir : FilePath
        directory to database of availability per location and technology

    Returns
    -------
    pd.DataFrame
        DataFrame of all Store components with assigned timeseries minimum availability
        for 8760 snapshots

    Raises
    ------
    ValueError
        If plant name has duplicate
    """
    dfs = []
    avail_df = pd.read_csv(avail_dir)
    stores_df = stores_df.to_dict("records")
    for row in stores_df:
        # filter country first
        avail_country_df = avail_df[avail_df.country == row["country"]].drop(
            "country", axis=1
        )
        # If the technology type repeats in avail_df in more than 1 row
        # meaning the technology should be matched by regional specific profile.
        if list(avail_country_df.technology).count(row["type"]) > 1:
            result = avail_country_df[
                (avail_country_df["technology"].isin([row["type"]]))
                & (avail_country_df["node"].isin([row["node"]]))
            ]
        # Else if the technology type appears in avail_df in only 1 row
        # all region will have similar availability profile for given technology.
        else:
            result = avail_country_df[
                (avail_country_df["technology"].isin([row["type"]]))
            ]

        if len(result) > 1:
            raise ValueError(
                f'{row["name"]} has duplicate. Please check input data again'
            )
        result = result.iloc[:, 2:]
        result["plant"] = row["name"]
        dfs.append(result.set_index("plant"))
    min_avail_df = pd.concat(dfs)
    min_avail_df = min_avail_df.astype(float)
    return min_avail_df


def get_storage_units_inflows(
    storage_capacity: pd.DataFrame,
    inflows_path: FilePath,
) -> pd.DataFrame:
    """Match each storage units to corresponding inflows.

    Note: Matching is done based on technology and region.

    Parameters
    ----------
    storage_capacity : pd.DataFrame
        DataFrame of all storage units
    inflows_path : FilePath
        inflow database

    Returns
    -------
    pd.DataFrame
        DataFrame of all storage units with assigned time-series inflows
        for 8760 snapshots
    """
    dfs = []
    storage_inflows = pd.read_csv(inflows_path)
    storage_unit_dict = storage_capacity.to_dict("records")

    for row in storage_unit_dict:
        # filter country first
        storage_inflows_country = storage_inflows[
            storage_inflows.country == row["country"]
        ].drop("country", axis=1)
        # If the technology type repeats in inflow_df in more than 1 row
        # meaning the technology should be matched by regional specific profile.
        if list(storage_inflows_country["technology"]).count(row["type"]) > 1:
            result = storage_inflows_country[
                (storage_inflows_country["technology"].isin([row["type"]]))
                & (storage_inflows_country["node"].isin([row["node"]]))
            ]
        # Else if the technology type appears in inflow_df in only 1 row
        # all region will have similar inflow profile for given technology.
        elif list(storage_inflows_country["technology"]).count(row["type"]) == 1:
            result = storage_inflows_country[
                (storage_inflows_country["technology"].isin([row["type"]]))
            ]
        # If technology type doesn't appear in inflow_df
        # a flat inflow profile of 0 (pu) is provided.
        else:
            # Country in row containing "Con00" is "WORLD"; thus we take
            # storage_inflows df instead
            result = storage_inflows[storage_inflows.technology == "Con00"].drop(
                "country", axis=1
            )

        result = result.iloc[:, 2:]
        result["node"] = row["node"]
        result["plant"] = row["name"]
        dfs.append(result)
    storage_inflows = pd.concat(dfs)
    storage_inflows = storage_inflows.drop(columns=["node"])
    storage_inflows.set_index(["plant"], inplace=True)
    return storage_inflows


def get_capital_cost(
    plant_type: list[str], tech_costs: pd.DataFrame, interest: float, currency: str
) -> float:
    """Calculate pypsa's capital cost from capital expenditure and fix O&M cost.

    Parameters
    ----------
    plant_type : list[str]
        list of plant types to calculate cost
    tech_costs : pd.DataFrame
        database of plant's costs
    interest : float
        discounted rate
    currency : str
        currency of model run

    Returns
    -------
    float
        pypsa capital cost for given technology
    """
    crf = (
        interest
        * (1 + interest) ** tech_costs.life__years
        / ((1 + interest) ** tech_costs.life__years - 1)
    )

    capex_column_name = f"cap__{str(currency).lower()}_mw"
    fom_column_name = f"fom__{str(currency).lower()}_mwa"
    # Handling incase of storage technologies with different unit
    if capex_column_name not in tech_costs.columns:
        capex_column_name = f"cap__{str(currency).lower()}_mwh"
    if fom_column_name not in tech_costs.columns:
        fom_column_name = f"fom__{str(currency).lower()}_mwh"

    return (
        (
            tech_costs[capex_column_name].astype(float) * crf
            + tech_costs[fom_column_name].astype(float)
        )
        .fillna(0)
        .loc[plant_type]
    )


def update_tech_fact_table(
    tech_table: pd.DataFrame,
    technologies_dir: FilePath,
    tech_costs_dir: FilePath,
    year: int,
    interest_dict: dict,
    currency: str,
) -> pd.DataFrame:
    """Process and add techno-economic columns to plants or links table.

    Parameters
    ----------
    tech_table : pd.DataFrame
        plants or links table
    technologies_dir : FilePath
        directory to database of technological params per technology
    tech_costs_dir : FilePath
        directory to database of cost params per technology
    year : int
        year of model run
    interest_dict : dict
        dict of interest rate
    currency : str
        currency of model run

    Returns
    -------
    pd.DataFrame
        Cleaned plant/link DataFrame with fact columns added
    """
    technology_df = pd.read_csv(technologies_dir)
    tech_costs = pd.read_csv(tech_costs_dir, index_col=["powerplant_type", "country"])
    # filter tech cost table to only corresponding model year
    tech_costs = tech_costs[tech_costs.year == year]

    tech_table = tech_table.to_dict("records")
    dfs = []
    for row in tech_table:
        # Matching plant type to technology to get technical params
        result = technology_df[
            (technology_df["technology"].isin([row["type"]]))
            & (technology_df["carrier"].isin([row["carrier"]]))
            & (technology_df["country"].isin([row["country"]]))
        ].loc[:, "class":]
        if result.empty:
            print(f'{row["type"]} is not in technology database')
        for column in row.keys():
            result[column] = row[column]

        # Processing capacity fact columns
        result["p_nom_min"] = row[f"p_nom_min_{year}"]
        result["p_nom_max"] = row[f"p_nom_max_{year}"]
        result["p_nom_extendable"] = row["p_nom_extendable"]
        result = result[
            [
                x
                for x in result.columns
                if "p_nom_max_" not in x and "p_nom_min_" not in x
            ]
        ]
        # Retrieve interest rate for given country
        try:
            interest = interest_dict[row["country"]]
        except KeyError:
            raise KeyError(f'No specific interest rate for {row["country"]} in config.')
        # Processing cost fact columns
        result["capital_cost"] = get_capital_cost(
            plant_type=(row["type"], row["country"]),
            tech_costs=tech_costs,
            interest=interest,
            currency=currency,
        )
        result["marginal_cost"] = (
            tech_costs[f"vom__{str(currency).lower()}_mwh"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["fom_cost"] = (
            tech_costs[f"fom__{str(currency).lower()}_mwa"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["inv_cost"] = (
            tech_costs[f"cap__{str(currency).lower()}_mw"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        if result["class"].values == "Link":
            for column in [
                "capital_cost",
                "marginal_cost",
                "fom_cost",
                "inv_cost",
            ]:
                result[column] = result[column] * result["efficiency"]
            for column in ["p_nom", "p_nom_max", "p_nom_min"]:
                result[column] = result[column] / result["efficiency"]

        result["lifetime"] = (
            tech_costs["life__years"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        # add noisy cost to marginal_cost if 0
        if result["marginal_cost"].values <= 0:
            result["marginal_cost"] += 1e-2 + 2e-3 * (
                np.random.random(len(result["marginal_cost"])) - 0.5
            )

        dfs.append(result)

    # combining df and sort columns
    plant_df = pd.concat(dfs)
    plant_df = plant_df[sorted(plant_df.columns)]
    plant_df.reset_index(drop=True)
    return plant_df


def update_storage_costs(
    storage_table: pd.DataFrame,
    storage_costs_dir: FilePath,
    year: int,
    interest_dict: dict,
    currency: str,
) -> pd.DataFrame:
    """Update storage cost parameters.

    capital, marginal, investment, and fixed O&M costs, as well as lifetime for each
    storage unit in the storage_table are updated based on the provided storage_costs
    database, model year, and interest rate.

    Parameters
    ----------
    storage_table : pd.DataFrame
        DataFrame containing storage units with at least 'type' and 'country' columns.
    storage_costs : FilePath
        Path to the CSV file containing storage cost data indexed by ['storage_type',
        'country'].
    year : int
        Model year to filter the storage_costs table.
    interest_dict : dict
        dict of interest rate.
    currency : str
        Currency in all uppercase, , ISO4217 format

    Returns
    -------
    pd.DataFrame
        The input storage_table with additional columns:
        'capital_cost', 'marginal_cost', 'lifetime', 'inv_cost', 'fom_cost'.
    """
    # Load and filter storage cost data for the specified year
    storage_costs_df = pd.read_csv(
        storage_costs_dir, index_col=["storage_type", "country"]
    )
    storage_costs_df = storage_costs_df[storage_costs_df.year == year]

    storage_table = storage_table.reset_index().to_dict("records")
    dfs = []
    for row in storage_table:
        result = pd.DataFrame([row])
        # Retrieve interest rate for given country
        try:
            interest = interest_dict[row["country"]]
        except KeyError:
            raise KeyError(f'No specific interest rate for {row["country"]} in config.')
        # Processing cost fact columns
        result["capital_cost"] = get_capital_cost(
            plant_type=(row["type"], row["country"]),
            tech_costs=storage_costs_df,
            interest=interest,
            currency=currency,
        )
        result["marginal_cost"] = (
            storage_costs_df[f"vom__{str(currency).lower()}_mwh"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["fom_cost"] = (
            storage_costs_df[f"fom__{str(currency).lower()}_mwh"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["inv_cost"] = (
            storage_costs_df[f"cap__{str(currency).lower()}_mwh"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["lifetime"] = (
            storage_costs_df["life__years"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        dfs.append(result)
    storage_table = pd.concat(dfs)
    storage_table = storage_table[sorted(storage_table.columns)]
    storage_table = storage_table.set_index("store")
    return storage_table


def update_ev_char_parameters(
    tech_df: pd.DataFrame,
    cost_df: pd.DataFrame,
    ev_param_dir: FilePath,
    year: int,
    interest_dict: dict,
    currency: str,
) -> pd.DataFrame:
    """Update ev storage and link parameters parameters.

    Parameters
    ----------
    tech_table : pd.DataFrame
        table of EV storage to update
    year : int
        year of model run
    ev_param_dir : FilePath
        directory of EV parameter database
    cost_df : pd.DataFrame
        DataFrame of cost parameters
    interest_dict : dict
        dict of interest rate.
    currency : str
        Currency in all uppercase, , ISO4217 format

    Returns
    -------
    pd.DataFrame
        Cleaned EV charger DataFrame with parameters updated
    """
    ev_param_df = pd.read_csv(ev_param_dir)
    cost_df_single_year = cost_df[cost_df.year == year]
    dfs = []
    for index, row in tech_df.iterrows():
        # Matching plant type to technology to get technical params
        result = ev_param_df[
            (ev_param_df["technology"].isin([row["type"]]))
            & (ev_param_df["carrier"].isin([row["carrier"]]))
            & (ev_param_df["country"].isin([row["country"]]))
        ].loc[:, "class":]
        if result.empty:
            print(f'{row["type"]} is not in technology database')
        for column in row.keys():
            result[column] = row[column]
        # Retrieve interest rate for given country
        try:
            interest = interest_dict[row["country"]]
        except KeyError:
            raise KeyError(f'No specific interest rate for {row["country"]} in config.')
        # processing costs for chargers
        result["capital_cost"] = get_capital_cost(
            plant_type=(row["type"], row["country"]),
            tech_costs=cost_df_single_year,
            interest=interest,
            currency=currency,
        )
        result["marginal_cost"] = (
            cost_df_single_year[f"vom__{str(currency).lower()}_mwh"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["fom_cost"] = (
            cost_df_single_year[f"fom__{str(currency).lower()}_mwa"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["inv_cost"] = (
            cost_df_single_year[f"cap__{str(currency).lower()}_mw"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )
        result["lifetime"] = (
            cost_df_single_year["life__years"]
            .astype(float)
            .fillna(0)
            .loc[(row["type"], row["country"])]
        )

        result = result.dropna(axis=1)
        result["p_nom_char"] = row[f"num_ch_{year}"] * result["p_nom_char"]
        result["link"] = index
        result["build_year"] = year

        if result["class"].values == "Link":
            for column in [
                "capital_cost",
                "marginal_cost",
                "fom_cost",
                "inv_cost",
            ]:
                result[column] = result[column] * result["efficiency"]
            for column in ["p_nom_char"]:
                result[column] = result[column] / result["efficiency"]

        dfs.append(result)

    plant_df = pd.concat(dfs)
    plant_df = plant_df[sorted(plant_df.columns)]
    plant_df = plant_df.set_index("link")
    return plant_df


def update_ev_store_parameters(
    tech_table: pd.DataFrame,
    year: int,
    ev_param_dir: FilePath,
) -> pd.DataFrame:
    """Calculate parameters for electric vehicle storage.

    Parameters
    ----------
    tech_table : pd.DataFrame
        table of EV storage to update
    year : int
        year of model run
    ev_param_dir : FilePath
        directory of EV parameter database

    Returns
    -------
    pd.DataFrame
        Cleaned EV storage table with parameters updated
    """
    ev_param_df = pd.read_csv(ev_param_dir)
    dfs = []
    for index, row in tech_table.iterrows():
        # matching the type of ev with the storage capacity
        result = ev_param_df[
            (
                ev_param_df["technology"].isin([row["type"]])
                & ev_param_df["country"].isin([row["country"]])
            )
        ].loc[:, "class":]
        if result.empty:
            print(f'{row["type"]} is not in ev parameter database')
        for column in row.keys():
            result[column] = row[column]
        result = result.dropna(axis=1)
        result["ev_e_nom"] = row[f"num_ev_{year}"] * result["e_nom_per_ev"]
        result["name"] = index
        result["build_year"] = year
        dfs.append(result)
    ev_store_clean_df = pd.concat(dfs)
    ev_store_clean_df = ev_store_clean_df[sorted(ev_store_clean_df.columns)]
    ev_store_clean_df = ev_store_clean_df[
        [x for x in ev_store_clean_df.columns if "ev_num_" not in x]
    ]
    ev_store_clean_df = ev_store_clean_df.set_index("name")
    return ev_store_clean_df


def get_link_min_availabilities(
    link_dir: FilePath, avail_dir: FilePath
) -> pd.DataFrame:
    """Match each links to corresponding maximum reverse flow based on technology.

    Parameters
    ----------
    link_dir : FilePath
        directory to a table of link to match
    avail_df_dir : FilePath
        directory to database of availability per location and technology

    Returns
    -------
    pd.DataFrame
        DataFrame of links with assigned time-series maximum reverse flow availability
        for 8760 snapshots
    """
    link_df = pd.read_csv(link_dir)
    avail_df = pd.read_csv(avail_dir)
    dfs = []
    link_df = link_df.to_dict("records")
    for row in link_df:
        # filter country first
        avail_country_df = avail_df[avail_df.country == row["country"]].drop(
            "country", axis=1
        )
        result = avail_country_df[
            (avail_country_df["technology"].isin([row["p_min_pu"]]))
        ]
        result = result.iloc[:, 2:]
        result["link"] = row["link"]
        dfs.append(result.set_index("link"))
    link_avail_df = pd.concat(dfs)

    return link_avail_df


def update_ev_link_cap(
    pev_load: pd.DataFrame, ev_link: pd.DataFrame, charge_ratio: float = 2
) -> pd.DataFrame:
    """Calculate EV charge capacities per buses based on PEV peak load and charge ratio.

    Parameters
    ----------
    pev_load : pd.DataFrame
        timeseries (8760 hours) for private electric vehicle loads
    ev_link : pd.DataFrame
        pev link table to update
    charge_ratio : float, optional
        ratio between charge capacity and peak load, by default 2

    Returns
    -------
    pd.DataFrame
        pev link table with p_nom updated and unused columns removed
    """
    print("########    PEV charging cap are updated")
    # Getting peak PEV for each bus
    pev_load_max = pev_load.reset_index(
        ["load_id", "node", "carrier", "sector"], drop=True
    ).max(axis=1)
    # Multiply peak load by charge ratio to get charging capacity for each bus
    for bus in pev_load_max.index:
        ev_link.loc[ev_link.bus1 == bus, "p_nom"] = pev_load_max[bus] * charge_ratio
    ev_link["p_nom_extendable"] = False
    ev_link = ev_link[
        [x for x in ev_link.columns if "p_nom_max" not in x and "p_nom_min" not in x]
    ]

    return ev_link


def update_ev_plant_cap(
    pev_load: pd.DataFrame, ev_plants: pd.DataFrame, store_ratio: float = 1
) -> pd.DataFrame:
    """Calculate EV storage capacities per buses based on PEV peak load and store ratio.

    Parameters
    ----------
    pev_load : pd.DataFrame
        timeseries (8760 hours) for private electric vehicle loads
    ev_plants : pd.DataFrame
        pev plants table to update
    store_ratio : float, optional
        ratio between storage capacity and peak load, by default 1

    Returns
    -------
    pd.DataFrame
        PEV plants table with p_nom updated and unused columns removed
    """
    print("########    PEV storage cap are updated")
    # Getting peak PEV for each bus
    pev_load_max = pev_load.reset_index(
        ["load_id", "node", "carrier", "sector"], drop=True
    ).max(axis=1)
    # Multiply peak load by charge ratio to get storage capacity for each bus
    for bus in pev_load_max.index:
        ev_plants.loc[ev_plants.bus == bus, "p_nom"] = pev_load_max[bus] * store_ratio
    ev_plants["p_nom_extendable"] = False
    ev_plants = ev_plants[
        [x for x in ev_plants.columns if "p_nom_max" not in x and "p_nom_min" not in x]
    ]
    return ev_plants


def plot_table(
    file_dir: FilePath,
    save_dir: FilePath,
    x: list,
    y: list,
    table: pd.DataFrame = None,
    days_range: int = 0,
):
    """Plot the charts based on output CSVs.

    Yearly --> stacked-column charts
    Hourly --> line charts

    Parameters
    ----------
    file_dir : FilePath
        Path of the CSVs which will be used for plotting.
    save_dir : FilePath
        Path to save the charts.
    x : list
        List of column name(s) used for x-axis.
    y : list
        List of column name(s) used for y-axis.
    table : pd.DataFrame, optional
        , by default None
    days_range : int, optional
        selected days range for hourly plot, by default 0
    """
    df = pd.read_csv(file_dir).fillna(0) if table is None else table.reset_index()

    if df.empty or "value" not in df.columns:
        return
    sort_list = [["snapshot", "value"] if "hourly" in file_dir else ["year", "value"]]
    df = df.sort_values(sort_list[0], ascending=[True, False])
    index_cols = list(set(df.columns.to_list()) - set(sort_list[0]))

    if len(index_cols) >= 2:
        df["index_cols"] = (
            "(" + df[index_cols].apply(lambda x: ", ".join(x), axis=1) + ")"
        )
        df = df.set_index("index_cols")
    elif len(index_cols) == 1:
        df = df.set_index(index_cols)
    else:  # len(index_cols) == 0
        df = df.set_index(df.columns[0])

    # Create Plotly figure
    data = []
    fig = go.Figure()
    try:
        if "hourly" in file_dir:
            # hourly data will be line charts
            df = df.reset_index().set_index("snapshot")
            df.index = pd.to_datetime(df.index)
            if days_range == 0:
                days = 7
            else:
                days = days_range

            type_col = df.columns[-2]
            df = df.iloc[: 24 * days * len(df[type_col].unique()), :]
            df = df.reset_index().set_index(index_cols)

            if "gen_by_type_hourly" in file_dir:
                # Create traces for each group
                for legend_name in df.index.unique():
                    trace = go.Scatter(
                        x=df[x][legend_name],
                        y=df["value"][legend_name].unique().tolist(),
                        mode="lines",
                        line={"dash": "dash", "width": 2},
                        name=legend_name,
                        marker_color=ag_cp[legend_name],
                    )
                    data.append(trace)
            else:
                for legend_name, color in zip(df.index.unique(), cycle(ag_st)):
                    trace = go.Scatter(
                        x=df[x][legend_name],
                        y=df["value"][legend_name].unique().tolist(),
                        mode="lines",
                        line={"dash": "dash", "width": 2},
                        name=str(legend_name),
                        marker_color=color,
                    )
                    data.append(trace)
        else:
            df[x] = df[x].astype(int)
            # yearly data will be column charts
            ag_cp_list = ["by_type_yearly", "avg", "emi_by_carrier", "share"]
            if any(x in file_dir for x in ag_cp_list):
                # Create traces for each group
                for legend_name in df.index.unique():
                    trace = go.Bar(
                        x=list(df.loc[df.index == legend_name][x]),
                        y=list(df.loc[df.index == legend_name][y]),
                        name=legend_name,
                        marker_color=ag_cp[legend_name],
                    )
                    data.append(trace)
            else:
                for legend_name, color in zip(df.index.unique(), cycle(ag_st)):
                    trace = go.Bar(
                        x=list(df.loc[df.index == legend_name][x]),
                        y=list(df.loc[df.index == legend_name][y]),
                        name=legend_name,
                        marker_color=color,
                    )
                    data.append(trace)

        fig = go.Figure(data=data)

        # Update layout properties
        fig.update_layout(
            barmode="stack",
            autosize=False,
            plot_bgcolor=agora_style["figure_edgecolor"],
            paper_bgcolor=agora_style["figure_facecolor"],
            font={
                "family": agora_style["font_family"],
                "size": agora_style["font_size"],
                "color": agora_style["axes_labelcolor"],
            },
            xaxis={
                "showline": True,  # Show axis line
                "linecolor": agora_style["axes_edgecolor"],  # Set axis line color
                "linewidth": agora_style["line_width"],  # Set axis line width
                "showticklabels": True,  # Show tick labels
                "ticks": "outside",  # Place ticks inside the plot
                "tickfont": {"color": agora_style["xtick_color"]},
                "tickwidth": agora_style["line_width"],  # Set tick width
                "showgrid": False,
            },
            yaxis={
                "showline": True,  # Show axis line
                "linecolor": agora_style["axes_edgecolor"],  # Set axis line color
                "linewidth": agora_style["line_width"],  # Set axis line width
                "showticklabels": True,  # Show tick labels
                "ticks": "outside",  # Place ticks inside the plot
                "tickfont": {"color": agora_style["ytick_color"]},
                "tickwidth": agora_style["line_width"],  # Set tick width
                "showgrid": False,
            },
            legend={
                "bgcolor": agora_style["legend_facecolor"],
            },
            showlegend=True,
        )
        fig.add_trace(
            go.Scatter(
                line={
                    "width": agora_style["line_width"],
                    "dash": agora_style["line_style"],
                }
            )
        )

        fig.update_xaxes(type="category")

        if "hourly" in file_dir:
            fig.update_xaxes(tickformat="%H:%M:%S")
        else:
            fig.update_yaxes(tickformat=",")  # not showing 3K but 3000

        if save_dir is not None:
            fig.write_image(save_dir, scale=3)
        else:
            min_value = [0 if min(df["value"]) >= 0 else min(df["value"])]
            max_value = [max(df["value"]) if "hourly" in file_dir else sum(df["value"])]
            fig.update_yaxes(range=[min_value, 1.2 * max_value[0]])
            fig.show()

    except (KeyError, TypeError) as e:
        print(e, "not plotting:", save_dir)


def type_to_bus_mapping(x: str) -> str:
    """Match pypsa type to bus type the technology should be located at.

    Parameters
    ----------
    x : str
        pypsa type

    Returns
    -------
    str
        bus type the technology should be located at
    """
    if x in ["PHOT", "FLOT", "HROR", "WTON", "WTOF"]:
        return "HVELEC"
    if x in ["RTPV"]:
        return "LVELEC"
    if x in ["SWHT"]:
        return "IND-LH"


def generation_type_mapping(x: str) -> str:
    """Match pypsa type to either renewable or fossil category.

    Parameters
    ----------
    x : str
        pypsa type

    Returns
    -------
    str
        "renewables"/"fossils"/"not defined"
    """
    if x in [
        "BIOT",
        "GEOT",
        "HROR",
        "PHOT",
        "CSP",
        "RTPV",
        "FLOT",
        "WTON",
        "WTOF",
        "HDAM",
        "WSTT",
        "OCHT",
        "HPHS",
    ]:
        return "renewables"
    if x in [
        "CCGT",
        "OCGT",
        "OILT",
        "SubC",
        "SupC",
        "NUCL",
        "SupC_RETRO",
        "SubC_RETRO",
    ]:
        return "fossils"

    return "not defined"


def add_unit_column(table_name: str, currency: str) -> str:
    """Return the appropriate unit string for a given output table name.

    Parameters
    ----------
    table_name : str
        the name of the table
    currency : str
        Currency in all uppercase, ISO4217 format

    Returns
    -------
    str
        The unit associated with the table.
    """
    table_name = table_name.lower()

    direct_mappings = [
        ("share", "%"),
        ("fuel_costs", f"{str(currency).lower()}/MWh_th"),
        ("flh", "hours"),
        ("emi_by", "MtCO2"),
    ]

    for keyword, unit in direct_mappings:
        if keyword in table_name:
            return unit

    # Generation or demand tables (yearly): use TWh, unless it's a share table
    if all(k in table_name for k in ["gen", "yearly"]) or ("dmd" in table_name):
        return "%" if "share" in table_name else "TWh"

    if "hourly" in table_name:
        return f"{str(currency).lower()}/MWh_el" if "price" in table_name else "MW"

    if any(k in table_name for k in ["cap_by", "capacity"]):
        return "GW"

    if any(
        k in table_name
        for k in ["capex", "opex", "fom", "overnight", "costs", "cost_by"]
    ):
        return f"Million {str(currency).lower()}"

    # Default fallback if no match found
    return "Dimensionless"


def filter_selected_countries_and_regions(
    df: pd.DataFrame,
    column: str,
    country_region: dict,
    currency: str,
    buses_csv: bool = False,
) -> pd.DataFrame:
    """Filter selected regions defined in the config.yaml.

    Parameters
    ----------
    df : pd.DataFrame
        Targeted DataFrame for filtering
    column: str
        Targeted column for filtering
    country_region : dict{str, list[str]}
        A dictionary with countries regions within those countries to filter by.
    currency: str
        Currency from the config file
    buses_csv : bool
        If True, apply filtering algorithms especially for buses.csv.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame
    """
    final_df = pd.DataFrame()
    for country, region in country_region.items():
        region_pattern = [country + "_" + x for x in region]
        # for buses cases
        if buses_csv:
            country_node_df = df[
                (df["country"].str.contains(country))
                & (~(df[column].str.contains("_")))
            ]
            region_df = df[(df[column].str.contains("|".join(region_pattern)))]
            filter_df = pd.concat([country_node_df, region_df])
        else:
            # for interconnector cases
            columns_to_check = [column, "bus0", "bus1", f"cap__{currency.lower()}_mw"]
            if (
                column == "link"
                and len([col for col in columns_to_check if col in df.columns]) == 4
            ):
                from_df = df[(df["bus0"].str.contains("|".join(region_pattern)))]
                filter_df = from_df[
                    (from_df["bus1"].str.contains("|".join(region_pattern)))
                ]
            # for storage_energy cases
            elif column == "store":
                region_df = df[(df[column].str.contains("|".join(region_pattern)))]
                extra_df = df[
                    (df["country"].str.contains(country))
                    & df[column].str.contains("|".join(["CO2STORN", "HYDN"]))
                ]
                filter_df = pd.concat([region_df, extra_df])
            else:
                filter_df = df[(df[column].str.contains("|".join(region_pattern)))]
        final_df = pd.concat([final_df, filter_df])
    return final_df


def load_scenario_config(path: str) -> dict:
    """Open and read scenario configuration from a YAML file.

    Parameters
    ----------
    path : str
        path to scenario_config.yaml

    Returns
    -------
    dict
        dictionary of scenario configuration
    """
    yaml_files = glob.glob(f"{path}/*.yaml")
    with open(yaml_files[0]) as f:
        return yaml.safe_load(f)


if __name__ == "__main__":
    type_to_bus_mapping("PHOT")
