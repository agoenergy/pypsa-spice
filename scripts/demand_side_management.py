# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""
Demand side management (DSM) for PyPSA-SPICE.

DSM is set per country and per load type (``profile_type``) in the
``demand_side_management`` section of ``scenario_config.yaml``. Two types exist:

- ``shed``: one Generator per bus (carrier ``DSM_SHED``, type ``DSMD``). Load can be
  left unserved at a marginal cost. The power is capped at a share of the peak of the
  load on the bus and the weighted energy at a share of the annual demand.
- ``shift``: one lossless StorageUnit per bus (carrier ``DSM_SHIFT``, type ``DSMF``).
  Charging means consuming earlier, discharging means consuming less. The state of
  charge is fixed at half of the energy capacity at every window boundary, so load can
  move earlier or later and the energy balances to zero in each window.

DSM is only applied with hourly resolution (``nth_hour`` with ``stepsize: 1``). With
any other resolution, a warning is logged and DSM is skipped.

The components are added here. The constraints are in ``custom_constraints.py``.
"""

import logging
import numbers

import numpy as np
import pandas as pd
import pypsa

logger = logging.getLogger(__name__)

DSM_SECTION = "demand_side_management"
DSM_TYPES = {
    "shed": {"carrier": "DSM_SHED", "type": "DSMD"},
    "shift": {"carrier": "DSM_SHIFT", "type": "DSMF"},
}
DSM_CARRIERS = [v["carrier"] for v in DSM_TYPES.values()]
DSM_COMPONENT_TYPES = [v["type"] for v in DSM_TYPES.values()]
REQUIRED_KEYS = {
    "shed": {"type", "max_peak_share", "max_energy_share", "marginal_cost"},
    "shift": {"type", "max_peak_share", "window_hours", "marginal_cost"},
}
OPTIONAL_KEYS = {"shed": set(), "shift": {"max_hours"}}
YEAR_KEYS = ["max_peak_share", "max_energy_share", "marginal_cost", "max_hours"]


def _is_number(value) -> bool:
    return isinstance(value, numbers.Number) and not isinstance(value, bool)


def value_for_year(value, year: int) -> float:
    """Return the value of a DSM parameter for a year (scalar or dict keyed by year)."""
    if isinstance(value, dict):
        return float(value[year])
    return float(value)


def validate_dsm_config(scenario_configs: dict, years: list, countries: list):
    """Validate the ``demand_side_management`` section of the scenario config.

    Validation always runs, also when the resolution does not support DSM, so that a
    broken config is found early.

    Parameters
    ----------
    scenario_configs : dict
        Scenario configurations.
    years : list
        All model years.
    countries : list
        Countries in the model.

    Raises
    ------
    ValueError
        If the section is malformed.
    """
    section = scenario_configs.get(DSM_SECTION) or {}
    errors = []
    if not isinstance(section, dict):
        raise ValueError(f"'{DSM_SECTION}' must be a mapping of countries.")
    for country, country_cfg in section.items():
        prefix = f"{DSM_SECTION}.{country}"
        if country not in countries:
            errors.append(f"{prefix}: country is not in base_configs.regions.")
            continue
        if not isinstance(country_cfg, dict):
            errors.append(f"{prefix}: must be a mapping with 'activate' and 'loads'.")
            continue
        unknown = set(country_cfg) - {"activate", "loads"}
        if unknown:
            errors.append(f"{prefix}: unknown keys {sorted(unknown)}.")
        if not isinstance(country_cfg.get("activate", False), bool):
            errors.append(f"{prefix}.activate: must be true or false.")
        loads = country_cfg.get("loads")
        # Inactive countries may leave 'loads' empty, but any settings given are
        # still validated so that errors show up before DSM is switched on.
        if not country_cfg.get("activate", False) and not loads:
            continue
        if not isinstance(loads, dict) or not loads:
            errors.append(f"{prefix}.loads: must list at least one load type.")
            continue
        for profile_type, load_cfg in loads.items():
            errors.extend(
                _validate_load_cfg(f"{prefix}.loads.{profile_type}", load_cfg, years)
            )
    if errors:
        raise ValueError(
            "Invalid demand_side_management configuration:\n  - "
            + "\n  - ".join(errors)
        )


def _validate_load_cfg(prefix: str, load_cfg, years: list) -> list:
    """Validate the DSM settings of one load type and return a list of errors."""
    if not isinstance(load_cfg, dict):
        return [f"{prefix}: must be a mapping."]
    dsm_type = load_cfg.get("type")
    if dsm_type not in DSM_TYPES:
        return [f"{prefix}.type: must be one of {sorted(DSM_TYPES)}, got {dsm_type!r}."]
    errors = []
    missing = REQUIRED_KEYS[dsm_type] - set(load_cfg)
    if missing:
        errors.append(f"{prefix}: missing keys {sorted(missing)} for type {dsm_type}.")
    unknown = set(load_cfg) - REQUIRED_KEYS[dsm_type] - OPTIONAL_KEYS[dsm_type]
    if unknown:
        errors.append(f"{prefix}: unknown keys {sorted(unknown)} for type {dsm_type}.")
    for key in YEAR_KEYS:
        if key not in load_cfg:
            continue
        value = load_cfg[key]
        values = value if isinstance(value, dict) else {None: value}
        if isinstance(value, dict):
            missing_years = [y for y in years if y not in value]
            if missing_years:
                errors.append(f"{prefix}.{key}: missing years {missing_years}.")
        for year, v in values.items():
            label = f"{prefix}.{key}" + (f"[{year}]" if year is not None else "")
            if not _is_number(v):
                errors.append(f"{label}: must be a number, got {v!r}.")
            elif key in ("max_peak_share", "max_energy_share") and not 0 <= v <= 1:
                errors.append(f"{label}: must be between 0 and 1, got {v}.")
            elif key == "marginal_cost" and v < 0:
                errors.append(f"{label}: must not be negative, got {v}.")
            elif key == "max_hours" and v <= 0:
                errors.append(f"{label}: must be positive, got {v}.")
    if "window_hours" in load_cfg:
        window = load_cfg["window_hours"]
        if not isinstance(window, int) or isinstance(window, bool) or window < 1:
            errors.append(
                f"{prefix}.window_hours: must be a positive integer, got {window!r}."
            )
    return errors


def dsm_resolution_supported(scenario_configs: dict) -> bool:
    """Return True if the temporal resolution supports DSM (hourly ``nth_hour``)."""
    resolution = scenario_configs["scenario_configs"]["resolution"]
    return (
        str(resolution.get("method", "")).lower() == "nth_hour"
        and resolution.get("stepsize") == 1
    )


def get_active_dsm_countries(
    scenario_configs: dict, countries: list, year: int
) -> dict:
    """Return the DSM settings of every country where DSM is active and supported.

    Logs a warning for every active country if the resolution does not support DSM.

    Returns
    -------
    dict
        ``{country: {profile_type: settings}}``, empty if DSM is skipped.
    """
    section = scenario_configs.get(DSM_SECTION) or {}
    active = {
        country: cfg["loads"]
        for country, cfg in section.items()
        if country in countries and cfg.get("activate", False)
    }
    if active and not dsm_resolution_supported(scenario_configs):
        resolution = scenario_configs["scenario_configs"]["resolution"]
        for country in active:
            logger.warning(
                "Demand side management is not implemented for resolution "
                "method=%s, stepsize=%s. Only nth_hour with stepsize 1 (hourly) is "
                "supported. DSM is skipped for %s in %s.",
                resolution.get("method"),
                resolution.get("stepsize"),
                country,
                year,
            )
        return {}
    return active


def remove_dsm_components(n: pypsa.Network):
    """Remove all DSM components, e.g. those carried over from the previous year."""
    for comp in ["Generator", "StorageUnit"]:
        df = n.components[comp].static
        if df.empty or "type" not in df:
            continue
        n.remove(comp, df.index[df["type"].isin(DSM_COMPONENT_TYPES)])


def add_dsm_components(
    n: pypsa.Network,
    full_loads: pd.DataFrame,
    dsm_settings: dict,
    year: int,
):
    """Add DSM shed generators and shift storage units for one year.

    Parameters
    ----------
    n : pypsa.Network
        Network to add the components to.
    full_loads : pd.DataFrame
        Loads with the full 8760 h profile before temporal aggregation. Index is the
        load name. Columns are ``bus``, ``country``, ``profile_type`` and one column
        per hour.
    dsm_settings : dict
        Output of :func:`get_active_dsm_countries`.
    year : int
        Model year.
    """
    meta_cols = ["bus", "country", "profile_type", "carrier", "node"]
    hour_cols = [c for c in full_loads.columns if c not in meta_cols]
    for country, loads_cfg in dsm_settings.items():
        for profile_type, cfg in loads_cfg.items():
            dsm_type = cfg["type"]
            sel = full_loads[
                (full_loads["country"] == country)
                & (full_loads["profile_type"] == profile_type)
            ]
            if sel.empty:
                logger.warning(
                    "DSM: no %s loads found in %s for %s. Nothing added.",
                    profile_type,
                    country,
                    year,
                )
                continue
            # peak of the load on each bus, from the full hourly profile
            peak = sel[hour_cols].astype(float).groupby(sel["bus"]).sum().max(axis=1)
            share = value_for_year(cfg["max_peak_share"], year)
            p_nom = (peak * share).clip(lower=0)
            p_nom = p_nom[p_nom > 0]
            if p_nom.empty:
                logger.warning(
                    "DSM: %s potential of %s in %s is zero in %s. Nothing added.",
                    dsm_type,
                    profile_type,
                    country,
                    year,
                )
                continue
            buses = p_nom.index
            carrier = DSM_TYPES[dsm_type]["carrier"]
            names = [f"{carrier} - {profile_type} - {bus}" for bus in buses]
            common = dict(
                bus=buses,
                country=country,
                carrier=carrier,
                type=DSM_TYPES[dsm_type]["type"],
                p_nom=p_nom.to_numpy(),
                p_nom_extendable=False,
                marginal_cost=value_for_year(cfg["marginal_cost"], year),
                build_year=year,
                lifetime=np.inf,
                dsm_profile_type=profile_type,
            )
            if dsm_type == "shed":
                n.add(
                    "Generator",
                    names,
                    dsm_max_energy_share=value_for_year(cfg["max_energy_share"], year),
                    **common,
                )
            else:
                window = int(cfg["window_hours"])
                max_hours = value_for_year(cfg.get("max_hours", window), year)
                if max_hours > window:
                    logger.warning(
                        "DSM: max_hours=%s of %s in %s is larger than "
                        "window_hours=%s and is clipped to %s.",
                        max_hours,
                        profile_type,
                        country,
                        window,
                        window,
                    )
                    max_hours = float(window)
                n.add(
                    "StorageUnit",
                    names,
                    max_hours=max_hours,
                    efficiency_store=1.0,
                    efficiency_dispatch=1.0,
                    standing_loss=0.0,
                    cyclic_state_of_charge=False,
                    state_of_charge_initial=p_nom.to_numpy() * max_hours / 2,
                    dsm_window_hours=float(window),
                    **common,
                )
            logger.info(
                "DSM: added %s %s on %d bus(es) for %s in %s (%.1f MW in total).",
                dsm_type,
                profile_type,
                len(buses),
                country,
                year,
                p_nom.sum(),
            )
