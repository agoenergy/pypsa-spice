"""Input-file discovery and direct editing for the web application.

The web editor intentionally uses the same files as the model.  All resolved paths are
confined to ``data/<dataset>/<project>/input`` and writes use an atomic replacement so
an interrupted request cannot leave a partially-written CSV or YAML file behind.
"""

from __future__ import annotations

import csv
import math
import os
import shutil
import tempfile
import threading
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Any, Literal

from fastapi import HTTPException
from pydantic import BaseModel, Field
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

EDITABLE_CONFIG_SECTIONS = ("scenario_configs", "co2_management", "custom_constraints")
SECTOR_LABELS = {"Power": "power", "Industry": "industry", "Transport": "transport"}
TABLE_FILE_FALLBACKS = {
    "Power_decommission": ("decommission_capacity.csv", "decomission_capacity.csv"),
    "Industry_decommission": ("decommission_capacity.csv", "decomission_capacity.csv"),
    "Direct_air_capture": ("direct_air_capture.csv",),
}
TABLE_COMPARE_KEYS = {
    "Power_decommission": ("country", "name", "class"),
    "Industry_decommission": ("country", "name", "class"),
    "Fuel_costs": ("country", "supply_plant", "year"),
    "Power_loads": ("country", "name", "year"),
    "Heat_loads": ("country", "name", "year"),
    "Transport_loads": ("country", "name", "year"),
    "Power_generators": ("country", "name"),
    "Heat_generators": ("country", "name"),
    "Power_links": ("country", "link"),
    "Heat_links": ("country", "link"),
    "Fuel_conversion": ("country", "link"),
    "Direct_air_capture": ("country", "link"),
    "Storage_units": ("country", "name"),
    "Industry_storage_units": ("country", "name"),
    "Stores": ("country", "store"),
    "Industry_stores": ("country", "store"),
    "Interconnectors": ("country", "link"),
    "Transport_chargers": ("country", "link"),
    "Transport_storages": ("country", "name"),
}
_WRITE_LOCK = threading.Lock()


class CellChange(BaseModel):
    """One cell mutation in an input CSV."""

    row: int = Field(ge=0)
    column: str = Field(min_length=1)
    value: Any


class TableUpdate(BaseModel):
    """Optimistic-concurrency request for a CSV table."""

    revision: str
    changes: list[CellChange] = Field(min_length=1, max_length=10000)


class ConfigSectionUpdate(BaseModel):
    """Optimistic-concurrency request for one scenario YAML section."""

    revision: str
    value: dict[str, Any]


class ConfigSectionsUpdate(BaseModel):
    """Optimistic-concurrency request for multiple scenario YAML sections."""

    revision: str
    sections: dict[str, dict[str, Any]]


def _load_yaml(path: Path, *, round_trip: bool = False) -> Any:
    yaml = YAML() if round_trip else YAML(typ="safe", pure=True)
    with path.open(encoding="utf-8") as handle:
        return yaml.load(handle) or {}


def _confine(candidate: Path, parent: Path, *, must_exist: bool = True) -> Path:
    """Resolve a user-selected path and guarantee it stays below ``parent``."""

    resolved_parent = parent.resolve()
    resolved = candidate.resolve()
    try:
        resolved.relative_to(resolved_parent)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="Invalid input path") from exc
    if must_exist and not resolved.exists():
        raise HTTPException(status_code=404, detail="Input file not found")
    return resolved


def _input_root(data_dir: Path, dataset: str, project: str) -> Path:
    project_root = _confine(data_dir / dataset / project, data_dir)
    input_root = _confine(project_root / "input", project_root)
    if not input_root.is_dir():
        raise HTTPException(status_code=404, detail="Project input folder not found")
    return input_root


def _revision(path: Path) -> str:
    stat = path.stat()
    return f"{stat.st_mtime_ns}:{stat.st_size}"


def _pretty_label(value: str) -> str:
    words = value.replace("_", " ").strip().split()
    return " ".join(word if word.isupper() else word.capitalize() for word in words)


def _comparison_value(value: Any) -> Any:
    """Normalise scalar input values without hiding meaningful differences."""

    if value is None:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return Decimal(str(value))
    if isinstance(value, str):
        stripped = value.strip()
        lowered = stripped.lower()
        if lowered in {"true", "false"}:
            return lowered == "true"
        if lowered in {"inf", "+inf", "-inf"}:
            return lowered.replace("+", "")
        try:
            return Decimal(stripped)
        except InvalidOperation:
            return stripped
    if isinstance(value, list):
        normalised = [_comparison_value(item) for item in value]
        return sorted(normalised, key=str)
    return value


def _display_delta(reference: Any, comparison: Any) -> int | float | None:
    left = _comparison_value(reference)
    right = _comparison_value(comparison)
    if not isinstance(left, Decimal) or not isinstance(right, Decimal):
        return None
    delta = right - left
    return int(delta) if delta == delta.to_integral_value() else float(delta)


def _change(
    *,
    status: str,
    item: str,
    country: str,
    parameter: str,
    reference: Any,
    comparison: Any,
) -> dict[str, Any]:
    return {
        "status": status,
        "item": item or "General",
        "country": country,
        "parameter": parameter,
        "reference": reference,
        "comparison": comparison,
        "delta": _display_delta(reference, comparison),
    }


def _technology_sector_mapping(settings_path: Path) -> dict[str, str]:
    """Return the Streamlit technology-to-sector mapping in web-friendly form."""

    mapping_path = settings_path.with_name("tech_mapping.csv")
    if not mapping_path.is_file():
        return {}
    with mapping_path.open(encoding="utf-8-sig", newline="") as handle:
        return {
            (row.get("original_names") or "").strip(): (
                row.get("sector") or ""
            ).strip().lower()
            for row in csv.DictReader(handle)
            if (row.get("original_names") or "").strip()
        }


def _technology_catalog(
    path: Path, sector_mapping: dict[str, str]
) -> list[dict[str, Any]]:
    """Summarise selectable technologies from a project's source table."""

    if not path.is_file():
        return []
    technologies: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            technology = (row.get("technology") or "").strip()
            if not technology:
                continue
            item = technologies.setdefault(
                technology,
                {
                    "id": technology,
                    "label": (row.get("technology_nomenclature") or technology).strip(),
                    "sector": sector_mapping.get(technology, "other") or "other",
                    "classes": set(),
                    "carriers": set(),
                },
            )
            if row.get("class"):
                item["classes"].add(row["class"].strip())
            if row.get("carrier"):
                item["carriers"].add(row["carrier"].strip())
    return [
        {
            **item,
            "classes": sorted(item["classes"]),
            "carriers": sorted(item["carriers"]),
        }
        for item in sorted(
            technologies.values(),
            key=lambda value: (value["label"].lower(), value["id"]),
        )
    ]


def input_catalog(data_dir: Path, settings_path: Path) -> dict[str, Any]:
    """Discover projects with model inputs and expose configured editable tables."""

    settings = _load_yaml(settings_path)
    datasets: list[dict[str, Any]] = []
    sector_mapping = _technology_sector_mapping(settings_path)
    if data_dir.is_dir():
        for dataset_path in sorted(path for path in data_dir.iterdir() if path.is_dir()):
            projects: list[dict[str, Any]] = []
            for project_path in sorted(path for path in dataset_path.iterdir() if path.is_dir()):
                input_root = project_path / "input"
                if not input_root.is_dir():
                    continue
                scenarios = sorted(
                    path.name
                    for path in input_root.iterdir()
                    if path.is_dir() and path.name != "global_input"
                )
                if not scenarios:
                    continue
                countries: set[str] = set()
                tech_path = input_root / "global_input" / "technologies.csv"
                if tech_path.is_file():
                    with tech_path.open(encoding="utf-8-sig", newline="") as handle:
                        countries.update(
                            row.get("country", "")
                            for row in csv.DictReader(handle)
                            if row.get("country")
                        )
                projects.append(
                    {
                        "name": project_path.name,
                        "scenarios": scenarios,
                        "countries": sorted(countries),
                        "technologies": _technology_catalog(tech_path, sector_mapping),
                    }
                )
            if projects:
                datasets.append({"name": dataset_path.name, "projects": projects})

    global_tables = [
        {
            "id": name,
            "label": _pretty_label(name),
            "scope": "global",
            "timeseries": bool(config.get("timeseries")),
        }
        for name, config in (settings.get("Global_input", {}) or {}).items()
    ]
    sector_tables: dict[str, list[dict[str, Any]]] = {}
    for label, slug in SECTOR_LABELS.items():
        sector_tables[slug] = [
            {
                "id": name,
                "label": _pretty_label(name),
                "scope": "scenario",
                "sector": slug,
            }
            for name, config in (settings.get(label, {}) or {}).items()
        ]
    grid_config = (settings.get("Grids", {}) or {}).get("Interconnectors")
    if grid_config:
        sector_tables.setdefault("power", []).append(
            {
                "id": "Interconnectors",
                "label": "Interconnectors",
                "scope": "scenario",
                "sector": "power",
            }
        )

    return {
        "table_query_version": 2,
        "datasets": datasets,
        "global_tables": global_tables,
        "sector_tables": sector_tables,
    }


def _table_config(
    settings_path: Path,
    scope: Literal["global", "scenario"],
    sector: str,
    table: str,
) -> dict[str, Any]:
    settings = _load_yaml(settings_path)
    if scope == "global":
        config = (settings.get("Global_input", {}) or {}).get(table)
    elif table == "Interconnectors" and sector == "power":
        config = (settings.get("Grids", {}) or {}).get("Interconnectors")
    else:
        label = next((name for name, slug in SECTOR_LABELS.items() if slug == sector), None)
        config = (settings.get(label, {}) or {}).get(table) if label else None
    if not config:
        raise HTTPException(status_code=404, detail="Input table is not configured")
    return dict(config)


def table_path(
    data_dir: Path,
    settings_path: Path,
    dataset: str,
    project: str,
    scenario: str,
    scope: Literal["global", "scenario"],
    sector: str,
    table: str,
) -> tuple[Path, dict[str, Any]]:
    """Resolve a configured CSV table within the selected input project."""

    input_root = _input_root(data_dir, dataset, project)
    config = _table_config(settings_path, scope, sector, table)
    if scope == "global":
        candidate = input_root / "global_input" / config["csv_name"]
    else:
        if not scenario:
            raise HTTPException(status_code=400, detail="Scenario is required")
        scenario_root = _confine(input_root / scenario, input_root)
        table_sector = "power" if table == "Interconnectors" else sector
        sector_root = scenario_root / table_sector
        candidate = sector_root / config["csv_name"]
        if not candidate.is_file():
            candidate = next(
                (
                    sector_root / filename
                    for filename in TABLE_FILE_FALLBACKS.get(table, ())
                    if (sector_root / filename).is_file()
                ),
                candidate,
            )
    path = _confine(candidate, input_root)
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Input CSV not found")
    return path, config


def _parse_number(raw: str) -> int | float | str:
    stripped = raw.strip()
    if stripped.lower() in {"inf", "+inf", "-inf"}:
        return stripped.lower().replace("+", "")
    number = float(stripped)
    return int(number) if number.is_integer() and not any(char in stripped.lower() for char in ".e") else number


def _column_kind(values: list[str]) -> Literal["boolean", "number", "string"]:
    populated = [value.strip() for value in values if value.strip()]
    if populated and all(value.lower() in {"true", "false"} for value in populated):
        return "boolean"
    if populated:
        try:
            for value in populated:
                _parse_number(value)
            return "number"
        except ValueError:
            pass
    return "string"


def _typed_value(raw: str, kind: str) -> Any:
    if raw == "":
        return ""
    if kind == "boolean":
        return raw.strip().lower() == "true"
    if kind == "number":
        return _parse_number(raw)
    return raw


def _technology_matches(
    row: dict[str, str],
    config: dict[str, Any],
    table: str,
    technology: str,
    technology_carriers: tuple[str, ...],
) -> bool:
    """Match one source row to a technology using model-data semantics."""

    if table == "Direct_air_capture":
        return technology == "DAC"
    filter_column = str(config.get("filter_col") or "")
    raw = (row.get(filter_column) or "").strip()
    if not raw:
        return False
    if technology == "PEVCH" and (raw.startswith("EVCH") or raw.startswith("EVST")):
        return True
    if "decommission" in table.lower():
        return raw.rsplit("_", 1)[-1] == technology
    if filter_column == "carrier":
        return raw in technology_carriers
    return raw == technology or technology in {
        value.strip() for value in raw.replace(";", ",").replace("|", ",").split(",")
    }


def read_table(
    path: Path,
    config: dict[str, Any],
    *,
    table: str = "",
    technology: str = "",
    technology_carriers: tuple[str, ...] = (),
    country: str = "ALL",
    filter_value: str = "ALL",
    query: str = "",
    offset: int = 0,
    limit: int = 100,
) -> dict[str, Any]:
    """Read a CSV with typed cells and legacy-compatible editability metadata."""

    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        raw_rows = list(reader)
    kinds = {
        column: _column_kind([row.get(column, "") or "" for row in raw_rows])
        for column in fieldnames
    }
    editable_table = not bool(config.get("timeseries"))
    columns = [
        {
            "name": column,
            "label": _pretty_label(column),
            "kind": kinds[column],
            "editable": editable_table
            and (kinds[column] in {"number", "boolean"} or column == "max_supply [MWh/year]"),
        }
        for column in fieldnames
    ]
    indexed_rows = list(enumerate(raw_rows))
    if technology:
        indexed_rows = [
            (index, row)
            for index, row in indexed_rows
            if _technology_matches(
                row,
                config,
                table,
                technology,
                technology_carriers,
            )
        ]

    country_column = next(
        (column for column in fieldnames if column.lower() == "country"), None
    )
    country_options = sorted(
        {
            (row.get(country_column) or "").strip()
            for _, row in indexed_rows
            if country_column and (row.get(country_column) or "").strip()
        }
    )
    if country != "ALL" and country_column:
        indexed_rows = [
            (index, row)
            for index, row in indexed_rows
            if (row.get(country_column) or "").strip() == country
        ]

    filter_column = str(config.get("filter_col") or "")
    filter_options = sorted(
        {
            (row.get(filter_column) or "").strip()
            for _, row in indexed_rows
            if filter_column
            and filter_column.lower() != "country"
            and (row.get(filter_column) or "").strip()
        }
    )
    if filter_value != "ALL" and filter_column:
        indexed_rows = [
            (index, row)
            for index, row in indexed_rows
            if (row.get(filter_column) or "").strip() == filter_value
        ]

    needle = query.strip().lower()
    if needle:
        indexed_rows = [
            (index, row)
            for index, row in indexed_rows
            if any(needle in str(value or "").lower() for value in row.values())
        ]

    total_filtered_rows = len(indexed_rows)
    page_rows = indexed_rows[offset : offset + limit]
    rows = [
        {
            "__row_id": index,
            **{
                column: _typed_value(row.get(column, "") or "", kinds[column])
                for column in fieldnames
            },
        }
        for index, row in page_rows
    ]
    return {
        "path": str(path.relative_to(path.parents[4])) if len(path.parents) > 4 else str(path),
        "revision": _revision(path),
        "columns": columns,
        "rows": rows,
        "total_rows": len(raw_rows),
        "total_filtered_rows": total_filtered_rows,
        "filter_column": filter_column or None,
        "country_options": country_options,
        "filter_options": filter_options,
    }


def _comparison_table_path(
    input_root: Path,
    scenario: str,
    sector: str,
    table: str,
    config: dict[str, Any],
) -> Path | None:
    """Resolve a scenario CSV for comparison while allowing a missing side."""

    scenario_root = _confine(input_root / scenario, input_root)
    table_sector = "power" if table == "Interconnectors" else sector
    sector_root = _confine(
        scenario_root / table_sector, scenario_root, must_exist=False
    )
    if table == "Direct_air_capture":
        direct_path = _confine(
            sector_root / "direct_air_capture.csv", input_root, must_exist=False
        )
        if direct_path.is_file():
            return direct_path
    candidate = _confine(
        sector_root / config["csv_name"], input_root, must_exist=False
    )
    if candidate.is_file():
        return candidate
    return next(
        (
            fallback
            for filename in TABLE_FILE_FALLBACKS.get(table, ())
            if (fallback := _confine(sector_root / filename, input_root, must_exist=False)).is_file()
        ),
        None,
    )


def _raw_csv(path: Path | None) -> tuple[list[str], list[dict[str, str]]]:
    if path is None:
        return [], []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        return reader.fieldnames or [], list(reader)


def _row_key(
    row: dict[str, str], table: str, fieldnames: list[str]
) -> tuple[str, ...]:
    return tuple(
        (row.get(column) or "").strip()
        for column in _row_key_columns(table, fieldnames)
    )


def _row_key_columns(table: str, fieldnames: list[str]) -> tuple[str, ...]:
    """Choose the stable CSV columns that identify a comparison row."""

    requested = TABLE_COMPARE_KEYS.get(table, ())
    columns = tuple(column for column in requested if column in fieldnames)
    if columns and columns != ("country",):
        return columns
    preferred = ("country", "name", "link", "store", "supply_plant", "bus", "year")
    return tuple(column for column in preferred if column in fieldnames)


def _index_rows(
    rows: list[dict[str, str]], table: str, fieldnames: list[str]
) -> tuple[dict[tuple[str, ...], dict[str, str]], set[tuple[str, ...]]]:
    indexed: dict[tuple[str, ...], dict[str, str]] = {}
    duplicates: set[tuple[str, ...]] = set()
    for row in rows:
        key = _row_key(row, table, fieldnames)
        if key in indexed:
            duplicates.add(key)
        indexed[key] = row
    return indexed, duplicates


def _row_item(
    key: tuple[str, ...], table: str, fieldnames: list[str]
) -> tuple[str, str]:
    values = dict(zip(_row_key_columns(table, fieldnames), key))
    country = values.pop("country", "")
    item = " · ".join(value for value in values.values() if value)
    return item, country


def _compare_csv_table(
    *,
    table: str,
    label: str,
    sector: str,
    reference_path: Path | None,
    comparison_path: Path | None,
) -> dict[str, Any] | None:
    reference_fields, reference_rows = _raw_csv(reference_path)
    comparison_fields, comparison_rows = _raw_csv(comparison_path)
    fieldnames = list(dict.fromkeys([*reference_fields, *comparison_fields]))
    reference_index, reference_duplicates = _index_rows(
        reference_rows, table, fieldnames
    )
    comparison_index, comparison_duplicates = _index_rows(
        comparison_rows, table, fieldnames
    )
    duplicate_keys = reference_duplicates | comparison_duplicates
    key_columns = set(_row_key_columns(table, fieldnames))
    changes: list[dict[str, Any]] = []
    for key in sorted(set(reference_index) | set(comparison_index)):
        item, country = _row_item(key, table, fieldnames)
        left = reference_index.get(key)
        right = comparison_index.get(key)
        if key in duplicate_keys:
            changes.append(
                _change(
                    status="ambiguous",
                    item=item,
                    country=country,
                    parameter="Duplicate row identity",
                    reference="Review source rows",
                    comparison="Review source rows",
                )
            )
            continue
        status = "added" if left is None else "removed" if right is None else "changed"
        for column in fieldnames:
            if column in key_columns:
                continue
            reference_value = left.get(column, "") if left is not None else None
            comparison_value = right.get(column, "") if right is not None else None
            if left is not None and right is not None and _comparison_value(
                reference_value
            ) == _comparison_value(comparison_value):
                continue
            if left is None and comparison_value in (None, ""):
                continue
            if right is None and reference_value in (None, ""):
                continue
            changes.append(
                _change(
                    status=status,
                    item=item,
                    country=country,
                    parameter=_pretty_label(column),
                    reference=reference_value,
                    comparison=comparison_value,
                )
            )
    if not changes:
        return None
    return {
        "id": f"input:{sector}:{table}",
        "label": label,
        "category": sector,
        "kind": "input",
        "changes": changes,
    }


def _flatten_config(
    value: Any, prefix: tuple[str, ...] = ()
) -> dict[tuple[str, ...], Any]:
    if isinstance(value, dict):
        flattened: dict[tuple[str, ...], Any] = {}
        for key, child in value.items():
            flattened.update(_flatten_config(child, (*prefix, str(key))))
        return flattened
    return {prefix: value}


def _compare_flat_values(
    reference: Any,
    comparison: Any,
    *,
    item_prefix: str = "",
) -> list[dict[str, Any]]:
    left = _flatten_config(reference)
    right = _flatten_config(comparison)
    changes: list[dict[str, Any]] = []
    for path in sorted(set(left) | set(right)):
        reference_value = left.get(path)
        comparison_value = right.get(path)
        if _comparison_value(reference_value) == _comparison_value(comparison_value):
            continue
        status = "added" if path not in left else "removed" if path not in right else "changed"
        country = path[0] if item_prefix == "country" and path else ""
        parameter_path = path[1:] if country else path
        if item_prefix == "settings" and len(path) > 1 and path[0] == "interest":
            country = path[1]
            parameter_path = (path[0], *path[2:])
        parameter = " · ".join(_pretty_label(part) for part in parameter_path)
        changes.append(
            _change(
                status=status,
                item=country or "General",
                country=country,
                parameter=parameter or "Value",
                reference=reference_value,
                comparison=comparison_value,
            )
        )
    return changes


def _effective_constraint(value: Any) -> dict[str, Any]:
    constraint = value if isinstance(value, dict) else {}
    if constraint.get("activate") is not True:
        return {"activate": False}
    return constraint


def _config_comparison_sections(
    reference: dict[str, Any], comparison: dict[str, Any]
) -> list[dict[str, Any]]:
    sections: list[dict[str, Any]] = []
    settings_changes = _compare_flat_values(
        reference.get("scenario_configs", {}),
        comparison.get("scenario_configs", {}),
        item_prefix="settings",
    )
    if settings_changes:
        sections.append(
            {
                "id": "config:scenario-settings",
                "label": "Scenario settings",
                "category": "configuration",
                "kind": "config",
                "changes": settings_changes,
            }
        )

    co2_changes = _compare_flat_values(
        reference.get("co2_management", {}),
        comparison.get("co2_management", {}),
        item_prefix="country",
    )
    if co2_changes:
        sections.append(
            {
                "id": "config:co2-management",
                "label": "CO₂ management",
                "category": "configuration",
                "kind": "config",
                "changes": co2_changes,
            }
        )

    left_constraints = reference.get("custom_constraints", {}) or {}
    right_constraints = comparison.get("custom_constraints", {}) or {}
    countries = sorted(set(left_constraints) | set(right_constraints))
    constraint_names = sorted(
        {
            name
            for country in countries
            for name in set((left_constraints.get(country) or {}))
            | set((right_constraints.get(country) or {}))
        }
    )
    for constraint_name in constraint_names:
        left = {
            country: _effective_constraint(
                (left_constraints.get(country) or {}).get(constraint_name)
            )
            for country in countries
        }
        right = {
            country: _effective_constraint(
                (right_constraints.get(country) or {}).get(constraint_name)
            )
            for country in countries
        }
        changes = _compare_flat_values(left, right, item_prefix="country")
        if changes:
            sections.append(
                {
                    "id": f"config:constraint:{constraint_name}",
                    "label": _pretty_label(constraint_name),
                    "category": "configuration",
                    "kind": "constraint",
                    "changes": changes,
                }
            )
    return sections


def compare_scenarios(
    data_dir: Path,
    settings_path: Path,
    dataset: str,
    project: str,
    reference: str,
    comparison: str,
) -> dict[str, Any]:
    """Return only effective model-input differences between two scenarios."""

    if reference == comparison:
        raise HTTPException(status_code=400, detail="Choose two different scenarios")
    input_root = _input_root(data_dir, dataset, project)
    reference_root = _confine(input_root / reference, input_root)
    comparison_root = _confine(input_root / comparison, input_root)
    if not reference_root.is_dir() or not comparison_root.is_dir():
        raise HTTPException(status_code=404, detail="Scenario input folder not found")

    reference_config = read_scenario_config(
        scenario_config_path(data_dir, dataset, project, reference)
    )["sections"]
    comparison_config = read_scenario_config(
        scenario_config_path(data_dir, dataset, project, comparison)
    )["sections"]
    sections = _config_comparison_sections(reference_config, comparison_config)

    settings = _load_yaml(settings_path)
    for sector_label, sector in SECTOR_LABELS.items():
        table_configs = dict(settings.get(sector_label, {}) or {})
        if sector == "power":
            interconnectors = (settings.get("Grids", {}) or {}).get("Interconnectors")
            if interconnectors:
                table_configs["Interconnectors"] = interconnectors
        for table, raw_config in table_configs.items():
            config = dict(raw_config)
            table_diff = _compare_csv_table(
                table=table,
                label=_pretty_label(table),
                sector=sector,
                reference_path=_comparison_table_path(
                    input_root, reference, sector, table, config
                ),
                comparison_path=_comparison_table_path(
                    input_root, comparison, sector, table, config
                ),
            )
            if table_diff:
                sections.append(table_diff)

    changes = [change for section in sections for change in section["changes"]]
    return {
        "dataset": dataset,
        "project": project,
        "reference": reference,
        "comparison": comparison,
        "global_inputs_shared": True,
        "summary": {
            "changes": len(changes),
            "groups": len(sections),
            "input_tables": sum(section["kind"] == "input" for section in sections),
            "configuration_groups": sum(
                section["kind"] != "input" for section in sections
            ),
        },
        "countries": sorted(
            {change["country"] for change in changes if change["country"]}
        ),
        "sections": sections,
    }


def _serialise_cell(value: Any, kind: str) -> str:
    if value is None or value == "":
        return ""
    if kind == "boolean":
        if isinstance(value, bool):
            return "TRUE" if value else "FALSE"
        if isinstance(value, str) and value.lower() in {"true", "false"}:
            return value.upper()
        raise HTTPException(status_code=422, detail="Boolean cells accept only true or false")
    if kind == "number":
        if isinstance(value, bool):
            raise HTTPException(status_code=422, detail="Numeric cells do not accept booleans")
        if isinstance(value, str) and value.strip().lower() in {"inf", "+inf", "-inf"}:
            return value.strip().lower().replace("+", "")
        try:
            number = float(value)
        except (TypeError, ValueError) as exc:
            raise HTTPException(status_code=422, detail="Numeric cells accept only numbers or inf") from exc
        if not math.isfinite(number):
            raise HTTPException(status_code=422, detail="Use inf or -inf for infinite values")
        return str(value)
    return str(value)


def _atomic_csv_write(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    descriptor, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)
            handle.flush()
            os.fsync(handle.fileno())
        shutil.copymode(path, temp_name)
        os.replace(temp_name, path)
    except Exception:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
        raise


def update_table(
    path: Path,
    config: dict[str, Any],
    update: TableUpdate,
    *,
    table: str = "",
    technology: str = "",
    technology_carriers: tuple[str, ...] = (),
    country: str = "ALL",
    filter_value: str = "ALL",
    query: str = "",
    offset: int = 0,
    limit: int = 100,
) -> dict[str, Any]:
    """Validate and atomically apply cell changes to a CSV."""

    if config.get("timeseries"):
        raise HTTPException(status_code=403, detail="Timeseries inputs are read-only in the web editor")
    with _WRITE_LOCK:
        if _revision(path) != update.revision:
            raise HTTPException(status_code=409, detail="This file changed on disk. Reload it before saving.")
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            rows = list(reader)
        kinds = {column: _column_kind([row.get(column, "") or "" for row in rows]) for column in fieldnames}
        editable = {
            column
            for column in fieldnames
            if kinds[column] in {"number", "boolean"} or column == "max_supply [MWh/year]"
        }
        seen: set[tuple[int, str]] = set()
        for change in update.changes:
            marker = (change.row, change.column)
            if marker in seen:
                raise HTTPException(status_code=422, detail="A cell may only be changed once per save")
            seen.add(marker)
            if change.row >= len(rows):
                raise HTTPException(status_code=422, detail="CSV row is out of range")
            if change.column not in editable:
                raise HTTPException(status_code=422, detail=f"Column '{change.column}' is read-only")
            rows[change.row][change.column] = _serialise_cell(change.value, kinds[change.column])
        _atomic_csv_write(path, fieldnames, rows)
    return read_table(
        path,
        config,
        table=table,
        technology=technology,
        technology_carriers=technology_carriers,
        country=country,
        filter_value=filter_value,
        query=query,
        offset=offset,
        limit=limit,
    )


def scenario_config_path(data_dir: Path, dataset: str, project: str, scenario: str) -> Path:
    input_root = _input_root(data_dir, dataset, project)
    scenario_root = _confine(input_root / scenario, input_root)
    path = _confine(scenario_root / "scenario_config.yaml", scenario_root)
    if not path.is_file():
        raise HTTPException(status_code=404, detail="Scenario configuration not found")
    return path


def _to_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _to_json(child) for key, child in value.items()}
    if isinstance(value, list):
        return [_to_json(child) for child in value]
    return value


def read_scenario_config(path: Path) -> dict[str, Any]:
    config = _load_yaml(path)
    return {
        "path": str(path),
        "revision": _revision(path),
        "sections": {name: _to_json(config.get(name, {}) or {}) for name in EDITABLE_CONFIG_SECTIONS},
    }


def _normalise_yaml_key(key: Any) -> Any:
    return int(key) if isinstance(key, str) and key.isdigit() else key


def _merge_yaml(existing: Any, incoming: Any) -> Any:
    """Update a round-trip YAML node while retaining comments on surviving keys."""

    if isinstance(incoming, dict):
        target = existing if isinstance(existing, CommentedMap) else CommentedMap()
        incoming_keys = {_normalise_yaml_key(key) for key in incoming}
        for key in list(target):
            if key not in incoming_keys:
                del target[key]
        for raw_key, child in incoming.items():
            key = _normalise_yaml_key(raw_key)
            target[key] = _merge_yaml(target.get(key), child)
        return target
    if isinstance(incoming, list):
        target = existing if isinstance(existing, CommentedSeq) else CommentedSeq()
        target.clear()
        target.extend(_merge_yaml(None, child) for child in incoming)
        return target
    return incoming


def _atomic_yaml_write(path: Path, document: Any) -> None:
    descriptor, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        yaml = YAML()
        yaml.preserve_quotes = True
        yaml.default_flow_style = False
        yaml.width = 4096
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            yaml.dump(document, handle)
            handle.flush()
            os.fsync(handle.fileno())
        shutil.copymode(path, temp_name)
        os.replace(temp_name, path)
    except Exception:
        if os.path.exists(temp_name):
            os.unlink(temp_name)
        raise


def update_scenario_section(path: Path, section: str, update: ConfigSectionUpdate) -> dict[str, Any]:
    """Atomically replace one editable YAML section while preserving comments."""

    if section not in EDITABLE_CONFIG_SECTIONS:
        raise HTTPException(status_code=404, detail="Scenario section is not editable")
    return update_scenario_sections(
        path,
        ConfigSectionsUpdate(revision=update.revision, sections={section: update.value}),
    )


def update_scenario_sections(path: Path, update: ConfigSectionsUpdate) -> dict[str, Any]:
    """Atomically replace multiple editable YAML sections while preserving comments."""

    if not update.sections:
        raise HTTPException(status_code=422, detail="At least one scenario section is required")
    if any(section not in EDITABLE_CONFIG_SECTIONS for section in update.sections):
        raise HTTPException(status_code=404, detail="Scenario section is not editable")
    with _WRITE_LOCK:
        if _revision(path) != update.revision:
            raise HTTPException(status_code=409, detail="This config changed on disk. Reload it before saving.")
        document = _load_yaml(path, round_trip=True)
        if not isinstance(document, CommentedMap):
            raise HTTPException(status_code=422, detail="Scenario config must be a YAML mapping")
        for section, value in update.sections.items():
            document[section] = _merge_yaml(document.get(section), value)
        _atomic_yaml_write(path, document)
    return read_scenario_config(path)
