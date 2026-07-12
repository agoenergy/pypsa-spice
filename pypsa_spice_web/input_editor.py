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
    return value.replace("_", " ").strip().title()


def input_catalog(data_dir: Path, settings_path: Path) -> dict[str, Any]:
    """Discover projects with model inputs and expose configured editable tables."""

    settings = _load_yaml(settings_path)
    datasets: list[dict[str, Any]] = []
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
                    }
                )
            if projects:
                datasets.append({"name": dataset_path.name, "projects": projects})

    global_tables = [
        {
            "id": name,
            "label": _pretty_label(name),
            "scope": "global",
            **config,
            "editable": not bool(config.get("timeseries")),
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
                **config,
                "editable": not bool(config.get("timeseries")),
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
                **grid_config,
                "editable": True,
            }
        )

    return {
        "datasets": datasets,
        "global_tables": global_tables,
        "sector_tables": sector_tables,
        "config_sections": list(EDITABLE_CONFIG_SECTIONS),
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


def read_table(path: Path, config: dict[str, Any], *, limit: int = 2000) -> dict[str, Any]:
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
    rows = [
        {"__row_id": index, **{column: _typed_value(row.get(column, "") or "", kinds[column]) for column in fieldnames}}
        for index, row in enumerate(raw_rows[:limit])
    ]
    return {
        "path": str(path.relative_to(path.parents[4])) if len(path.parents) > 4 else str(path),
        "revision": _revision(path),
        "columns": columns,
        "rows": rows,
        "total_rows": len(raw_rows),
        "truncated": len(raw_rows) > limit,
        "filter_column": config.get("filter_col"),
        "with_charts": bool(config.get("with_charts")),
        "timeseries": bool(config.get("timeseries")),
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


def update_table(path: Path, config: dict[str, Any], update: TableUpdate) -> dict[str, Any]:
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
    return read_table(path, config)


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
    with _WRITE_LOCK:
        if _revision(path) != update.revision:
            raise HTTPException(status_code=409, detail="This config changed on disk. Reload it before saving.")
        document = _load_yaml(path, round_trip=True)
        if not isinstance(document, CommentedMap):
            raise HTTPException(status_code=422, detail="Scenario config must be a YAML mapping")
        document[section] = _merge_yaml(document.get(section), update.value)
        _atomic_yaml_write(path, document)
    return read_scenario_config(path)
