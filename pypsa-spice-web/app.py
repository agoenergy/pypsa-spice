# SPDX-FileCopyrightText: PyPSA-SPICE Developers

# SPDX-License-Identifier: GPL-2.0-or-later

"""FastAPI application for the PyPSA-SPICE web workspace."""

from __future__ import annotations

import argparse
import csv
import io
import math
import shutil
import subprocess
import sys
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

import uvicorn
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from input_editor import (
    ConfigSectionsUpdate,
    ConfigSectionUpdate,
    TableUpdate,
    compare_scenarios,
    input_catalog,
    load_yaml,
    read_scenario_config,
    read_table,
    scenario_config_path,
    table_path,
    update_scenario_section,
    update_scenario_sections,
    update_table,
)
from pydantic import BaseModel, Field
from scenario_runner import (
    RunConflictError,
    RunValidationError,
    ScenarioRunManager,
)
from scenario_workspace import (
    ScenarioWorkspace,
    ScenarioWorkspaceError,
)

HOST = "127.0.0.1"
PORT = 8000
VITE_PORT = 5173  # frontend dev server, used only with --dev

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = Path(__file__).resolve().parent
FRONTEND_DIR = PACKAGE_DIR / "frontend"
FRONTEND_DIST = FRONTEND_DIR / "dist"
DATA_DIR = ROOT / "data"
GRAPH_CONFIG = PACKAGE_DIR / "graph_settings.yaml"
MAPPING_DIR = ROOT / "pypsa-spice-vis" / "setting"
INPUT_CONFIG = MAPPING_DIR / "input_settings.yaml"
RUN_MANAGER = ScenarioRunManager(ROOT)
SCENARIO_WORKSPACE = ScenarioWorkspace(ROOT, active_run=RUN_MANAGER.active)


class ModelRunRequest(BaseModel):
    """Editable base-config path values for a Snakemake model run."""

    dataset: str = Field(min_length=1, max_length=255)
    project: str = Field(min_length=1, max_length=255)
    input_scenario: str = Field(min_length=1, max_length=255)
    output_scenario: str = Field(min_length=1, max_length=255)
    cores: int = Field(default=1, ge=1, le=128)
    force: bool = False


class ScenarioCloneRequest(BaseModel):
    """Create a local scenario by duplicating an existing complete scenario."""

    dataset: str = Field(min_length=1, max_length=255)
    project: str = Field(min_length=1, max_length=255)
    source_scenario: str = Field(min_length=1, max_length=255)
    new_scenario: str = Field(min_length=1, max_length=255)


SECTION_META = {
    "power": {
        "label": "Power",
        "title": "Power Sector",
    },
    "industry": {
        "label": "Industry",
        "title": "Industry Sector",
    },
    "transport": {
        "label": "Transport",
        "title": "Transport Sector",
    },
    "emissions": {
        "label": "Emissions",
        "title": "Emissions",
    },
    "costs": {
        "label": "Costs",
        "title": "System Costs",
    },
}


def _read_mapping(filename: str) -> dict[str, dict[str, str]]:
    """Read a name mapping CSV into display labels and colours by name."""
    path = MAPPING_DIR / filename
    if not path.exists():
        return {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return {
            row["original_names"]: {
                "label": row.get("nice_names") or row["original_names"],
                "color": row.get("hex_codes") or "",
            }
            for row in csv.DictReader(handle)
            if row.get("original_names")
        }


def _chart_definitions() -> dict[str, list[dict[str, Any]]]:
    """Load the configured charts grouped by results section."""
    raw = load_yaml(GRAPH_CONFIG)
    sections: dict[str, list[dict[str, Any]]] = {}
    for section, charts in raw.items():
        if section not in SECTION_META:
            continue
        sections[section] = []
        for chart_id, values in charts.items():
            chart = {"id": chart_id, **values}
            chart["hourly"] = "hourly" in chart["type"]
            sections[section].append(chart)
    return sections


CHARTS = _chart_definitions()
MAPPINGS = {
    **_read_mapping("tech_mapping.csv"),
    **_read_mapping("carrier_mapping.csv"),
}


def _scenario_root(dataset: str, project: str, scenario: str) -> Path:
    """Resolve the CSV results folder of a scenario inside the data dir."""
    candidate = (DATA_DIR / dataset / project / "results" / scenario / "csvs").resolve()
    try:
        candidate.relative_to(DATA_DIR.resolve())
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="Invalid result path") from exc
    if not candidate.is_dir():
        raise HTTPException(status_code=404, detail="Result scenario not found")
    return candidate


def _sector_root(dataset: str, project: str, scenario: str, sector: str) -> Path:
    """Resolve the results folder of one sector within a scenario."""
    scenario_root = _scenario_root(dataset, project, scenario)
    root = (scenario_root / sector).resolve()
    if root.parent != scenario_root or not root.is_dir():
        raise HTTPException(status_code=404, detail="Sector result not found")
    return root


def _year_sort(value: str) -> tuple[int, str]:
    """Sort key that orders numeric years before other folder names."""
    return (0, f"{int(value):08d}") if value.isdigit() else (1, value)


def _table_paths(
    dataset: str,
    project: str,
    scenario: str,
    sector: str,
    table: str,
    year: str | None,
    hourly: bool,
) -> list[Path]:
    """Find the result CSV files for a configured table and year."""
    sector_path = _sector_root(dataset, project, scenario, sector)
    chart = next(
        (
            chart
            for charts in CHARTS.values()
            for chart in charts
            if chart["table_name"] == table
        ),
        None,
    )
    if chart is None:
        raise HTTPException(status_code=404, detail="Result table is not configured")
    if bool(chart["hourly"]) != hourly:
        raise HTTPException(status_code=400, detail="Invalid result table mode")

    year_dirs = sorted(
        [
            path
            for path in sector_path.iterdir()
            if path.is_dir() and path.name.isdigit()
        ],
        key=lambda path: _year_sort(path.name),
    )
    if hourly:
        available_years = {path.name for path in year_dirs}
        if year and year not in available_years:
            raise HTTPException(status_code=400, detail="Invalid result year")
        selected = year or next((path.name for path in year_dirs), None)
        paths = [sector_path / selected / f"{table}.csv"] if selected else []
    else:
        all_years = sector_path / "all_years" / f"{table}.csv"
        paths = (
            [all_years]
            if all_years.exists()
            else [p / f"{table}.csv" for p in year_dirs]
        )
    safe_paths: list[Path] = []
    for path in paths:
        resolved = path.resolve()
        try:
            resolved.relative_to(sector_path)
        except ValueError as exc:
            raise HTTPException(status_code=400, detail="Invalid result path") from exc
        if resolved.is_file():
            safe_paths.append(resolved)
    return safe_paths


@lru_cache(maxsize=256)
def _read_csv_cached(path_string: str, modified_ns: int) -> tuple[dict[str, Any], ...]:
    """Read a result CSV with typed values, cached per file version."""
    del modified_ns
    path = Path(path_string)
    rows: list[dict[str, Any]] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, Any] = {}
            for key, value in row.items():
                if key == "value":
                    try:
                        parsed[key] = float(value) if value not in (None, "") else None
                    except ValueError:
                        parsed[key] = None
                elif key == "year":
                    try:
                        parsed[key] = int(float(value))
                    except (TypeError, ValueError):
                        parsed[key] = value
                else:
                    parsed[key] = value
            if parsed.get("value") is not None:
                rows.append(parsed)
    return tuple(rows)


@lru_cache(maxsize=32)
def _read_paths_cached(
    versions: tuple[tuple[str, int], ...],
) -> tuple[dict[str, Any], ...]:
    """Concatenate result CSVs and drop duplicated rows, cached per file versions."""
    seen: set[tuple[Any, ...]] = set()
    unique_rows: list[dict[str, Any]] = []
    for path_string, modified_ns in versions:
        for row in _read_csv_cached(path_string, modified_ns):
            # Some exporters repeat the full yearly table in each year folder.
            marker = tuple(row.values())
            if marker not in seen:
                seen.add(marker)
                unique_rows.append(row)
    return tuple(unique_rows)


def _read_paths(paths: list[Path]) -> tuple[dict[str, Any], ...]:
    """Read result CSVs; returned rows are shared and must not be mutated."""
    return _read_paths_cached(tuple((str(p), p.stat().st_mtime_ns) for p in paths))


def _combine_csv_files(paths: list[Path]) -> str:
    """Stack result CSVs into one CSV, keeping every row and value as written."""
    columns: list[str] = []
    rows: list[dict[str, str]] = []
    for path in paths:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle)
            # Year files can differ in columns; keep them all, in order of appearance.
            columns += [name for name in reader.fieldnames or [] if name not in columns]
            rows += reader
    buffer = io.StringIO()
    # A column missing from a year file is left empty for that year's rows.
    writer = csv.DictWriter(buffer, fieldnames=columns, restval="")
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()


def _sample_hourly(
    rows: list[dict[str, Any]], legend_column: str, limit: int
) -> tuple[list[dict[str, Any]], bool]:
    """Downsample hourly rows evenly per series to stay within a limit."""
    if len(rows) <= limit:
        return rows, False
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(str(row.get(legend_column, "Series")), []).append(row)
    per_group, remainder = divmod(limit, len(groups))
    sampled: list[dict[str, Any]] = []
    for group_index, group in enumerate(groups.values()):
        target = per_group + (1 if group_index < remainder else 0)
        if target == 0:
            continue
        step = max(1, math.ceil(len(group) / target))
        sampled.extend(group[::step][:target])
    return sampled, True


def _chart_axis_extents(
    rows: list[dict[str, Any]], chart: dict[str, Any]
) -> dict[str, list[float] | None]:
    """Calculate exact axis extents before hourly rows are sampled."""
    x_column = "snapshot" if chart["hourly"] else "year"
    legend_column = chart["leg_col"]
    secondary_labels = set(chart.get("secondary_y_lab", []))
    aggregates: dict[tuple[str, str], float] = {}
    for row in rows:
        x_value = row.get(x_column)
        value = row.get("value")
        if x_value in (None, "") or not isinstance(value, (int, float)):
            continue
        legend = str(row.get(legend_column, "Series"))
        key = (str(x_value), legend)
        aggregates[key] = aggregates.get(key, 0.0) + float(value)

    should_stack = (chart["type"] != "grouped_bar" and "bar" in chart["type"]) or chart[
        "type"
    ] == "area_share"
    extents: dict[bool, list[float] | None] = {False: None, True: None}

    if should_stack:
        totals: dict[tuple[bool, str], list[float]] = {}
        for (x_value, legend), value in aggregates.items():
            axis_key = (legend in secondary_labels, x_value)
            total = totals.setdefault(axis_key, [0.0, 0.0])
            total[0 if value < 0 else 1] += value
        for (secondary, _), (negative, positive) in totals.items():
            current = extents[secondary]
            extents[secondary] = (
                [negative, positive]
                if current is None
                else [min(current[0], negative), max(current[1], positive)]
            )
    else:
        for (_, legend), value in aggregates.items():
            secondary = legend in secondary_labels
            current = extents[secondary]
            extents[secondary] = (
                [value, value]
                if current is None
                else [min(current[0], value), max(current[1], value)]
            )

    return {"primary": extents[False], "secondary": extents[True]}


def _parse_timestamp(value: str, field_name: str) -> datetime:
    """Parse an ISO timestamp supplied by the web chart controls."""
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise HTTPException(
            status_code=400, detail=f"Invalid {field_name} timestamp"
        ) from exc


def _catalog() -> dict[str, Any]:
    """List datasets, projects and result scenarios available on disk."""
    datasets: list[dict[str, Any]] = []
    if DATA_DIR.exists():
        for dataset_path in sorted(
            path for path in DATA_DIR.iterdir() if path.is_dir()
        ):
            projects: list[dict[str, Any]] = []
            for project_path in sorted(
                path for path in dataset_path.iterdir() if path.is_dir()
            ):
                results = project_path / "results"
                if not results.is_dir():
                    continue
                scenarios: list[dict[str, Any]] = []
                for scenario_path in sorted(
                    path for path in results.iterdir() if path.is_dir()
                ):
                    csv_root = scenario_path / "csvs"
                    if not csv_root.is_dir():
                        continue
                    sectors: list[dict[str, Any]] = []
                    for sector_path in sorted(
                        path for path in csv_root.iterdir() if path.is_dir()
                    ):
                        years = sorted(
                            [
                                path.name
                                for path in sector_path.iterdir()
                                if path.is_dir()
                            ],
                            key=_year_sort,
                        )
                        sectors.append(
                            {
                                "name": sector_path.name,
                                "years": [year for year in years if year.isdigit()],
                            }
                        )
                    if sectors:
                        scenarios.append(
                            {"name": scenario_path.name, "sectors": sectors}
                        )
                if scenarios:
                    projects.append({"name": project_path.name, "scenarios": scenarios})
            if projects:
                datasets.append({"name": dataset_path.name, "projects": projects})

    return {
        "datasets": datasets,
        "sections": [
            {"id": key, **SECTION_META[key], "charts": charts}
            for key, charts in CHARTS.items()
        ],
        "mappings": MAPPINGS,
    }


app = FastAPI(
    title="PyPSA-SPICE Web Workspace (beta)",
    description=(
        "Explore model results and directly edit PyPSA-SPICE CSV/YAML inputs. "
        "Beta: under active development, so endpoints and behaviour may change."
    ),
    version="0.2.0-beta",
)
app.add_middleware(GZipMiddleware, minimum_size=1000)


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> Response:
    """Serve the built frontend, or a hint if it has not been built."""
    del request
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(index_file)
    return HTMLResponse(
        "<h1>Frontend not built</h1><p>Start the app with "
        "<code>python pypsa-spice-web/app.py</code>, which builds it.</p>",
        status_code=503,
    )


@app.get("/health")
def health() -> dict[str, str]:
    """Report that the service is running."""
    return {"status": "ok", "service": "pypsa-spice-web"}


@app.get("/api/runs/options")
def model_run_options() -> JSONResponse:
    """Describe the repository base config used for web-launched runs."""
    try:
        return JSONResponse(RUN_MANAGER.options())
    except RunValidationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.get("/api/runs/latest")
def latest_model_run() -> JSONResponse:
    """Return the latest persisted web run, including its current log tail."""
    return JSONResponse(RUN_MANAGER.latest())


@app.get("/api/runs/{run_id}")
def model_run(run_id: str) -> JSONResponse:
    """Return current state for one web-launched Snakemake run."""
    try:
        return JSONResponse(RUN_MANAGER.get(run_id))
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Model run not found") from exc


@app.post("/api/runs", status_code=202)
def start_model_run(request: ModelRunRequest) -> JSONResponse:
    """Create a run-local base config and launch the full model workflow."""
    try:
        run = RUN_MANAGER.start(
            dataset=request.dataset,
            project=request.project,
            input_scenario=request.input_scenario,
            output_scenario=request.output_scenario,
            cores=request.cores,
            force=request.force,
        )
    except RunConflictError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc
    except RunValidationError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    return JSONResponse(run, status_code=202)


@app.delete("/api/runs/{run_id}")
def cancel_model_run(run_id: str) -> JSONResponse:
    """Request termination of the active Snakemake process."""
    try:
        return JSONResponse(RUN_MANAGER.cancel(run_id))
    except KeyError as exc:
        raise HTTPException(status_code=404, detail="Model run not found") from exc


@app.get("/api/catalog")
def catalog() -> JSONResponse:
    """Return the catalog of available result scenarios."""
    return JSONResponse(_catalog())


@app.get("/api/input/catalog")
def input_data_catalog() -> JSONResponse:
    """Return projects, scenarios, and configured input-table metadata."""
    return JSONResponse(input_catalog(DATA_DIR, INPUT_CONFIG))


@app.get("/api/input/workspace")
def input_workspace_status(dataset: str) -> JSONResponse:
    """Return web-mutation lock and dataset Git status for review and run."""
    try:
        return JSONResponse(SCENARIO_WORKSPACE.status(dataset))
    except ScenarioWorkspaceError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.post("/api/input/scenarios", status_code=201)
def create_input_scenario(request: ScenarioCloneRequest) -> JSONResponse:
    """Atomically clone a complete scenario inside the local data worktree."""
    try:
        scenario = SCENARIO_WORKSPACE.clone(
            dataset=request.dataset,
            project=request.project,
            source_scenario=request.source_scenario,
            new_scenario=request.new_scenario,
        )
    except ScenarioWorkspaceError as exc:
        detail = str(exc)
        status_code = 409 if "locked" in detail or "already exists" in detail else 422
        raise HTTPException(status_code=status_code, detail=detail) from exc
    return JSONResponse(scenario, status_code=201)


@app.get("/api/input/table")
def input_table_data(
    dataset: str,
    project: str,
    scenario: str = "",
    scope: Literal["global", "scenario"] = Query(),
    sector: str = Query(pattern="^(power|industry|transport)$"),
    table: str = Query(min_length=1),
    technology: str = "",
    technology_carrier: list[str] = Query(default=[]),
    country: str = "ALL",
    filter_value: str = "ALL",
    query: str = "",
    offset: int = Query(default=0, ge=0),
    limit: int = Query(default=100, ge=1, le=500),
) -> JSONResponse:
    """Read a configured model-input CSV with typed editing metadata."""
    path, config = table_path(
        DATA_DIR,
        INPUT_CONFIG,
        dataset,
        project,
        scenario,
        scope,
        sector,
        table,
    )
    return JSONResponse(
        read_table(
            path,
            config,
            table=table,
            technology=technology,
            technology_carriers=tuple(technology_carrier),
            country=country,
            filter_value=filter_value,
            query=query,
            offset=offset,
            limit=limit,
        )
    )


@app.put("/api/input/table")
def save_input_table(
    update: TableUpdate,
    dataset: str,
    project: str,
    scenario: str = "",
    scope: Literal["global", "scenario"] = Query(),
    sector: str = Query(pattern="^(power|industry|transport)$"),
    table: str = Query(min_length=1),
    technology: str = "",
    technology_carrier: list[str] = Query(default=[]),
    country: str = "ALL",
    filter_value: str = "ALL",
    query: str = "",
    offset: int = Query(default=0, ge=0),
    limit: int = Query(default=100, ge=1, le=500),
) -> JSONResponse:
    """Directly and atomically apply validated changes to an input CSV."""
    try:
        SCENARIO_WORKSPACE.ensure_mutable()
    except ScenarioWorkspaceError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc

    path, config = table_path(
        DATA_DIR,
        INPUT_CONFIG,
        dataset,
        project,
        scenario,
        scope,
        sector,
        table,
    )
    return JSONResponse(
        update_table(
            path,
            config,
            update,
            table=table,
            technology=technology,
            technology_carriers=tuple(technology_carrier),
            country=country,
            filter_value=filter_value,
            query=query,
            offset=offset,
            limit=limit,
        )
    )


@app.get("/api/input/scenario-config")
def scenario_configuration(dataset: str, project: str, scenario: str) -> JSONResponse:
    """Read the editable sections of a scenario configuration."""
    path = scenario_config_path(DATA_DIR, dataset, project, scenario)
    return JSONResponse(read_scenario_config(path))


@app.get("/api/input/compare")
def compare_scenario_inputs(
    dataset: str,
    project: str,
    reference: str = Query(min_length=1),
    comparison: str = Query(min_length=1),
) -> JSONResponse:
    """Show only configuration and scenario-input differences."""
    return JSONResponse(
        compare_scenarios(
            DATA_DIR,
            INPUT_CONFIG,
            dataset,
            project,
            reference,
            comparison,
        )
    )


@app.put("/api/input/scenario-config")
def save_scenario_configuration_sections(
    update: ConfigSectionsUpdate,
    dataset: str,
    project: str,
    scenario: str,
) -> JSONResponse:
    """Directly save multiple scenario-config sections in one atomic write."""
    try:
        SCENARIO_WORKSPACE.ensure_mutable()
    except ScenarioWorkspaceError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc

    path = scenario_config_path(DATA_DIR, dataset, project, scenario)
    return JSONResponse(update_scenario_sections(path, update))


@app.put("/api/input/scenario-config/{section}")
def save_scenario_configuration_section(
    section: str,
    update: ConfigSectionUpdate,
    dataset: str,
    project: str,
    scenario: str,
) -> JSONResponse:
    """Directly save one scenario-config section with YAML comments retained."""
    try:
        SCENARIO_WORKSPACE.ensure_mutable()
    except ScenarioWorkspaceError as exc:
        raise HTTPException(status_code=409, detail=str(exc)) from exc

    path = scenario_config_path(DATA_DIR, dataset, project, scenario)
    return JSONResponse(update_scenario_section(path, section, update))


@app.get("/api/chart")
def chart_data(
    dataset: str,
    project: str,
    scenario: str,
    sector: str,
    table: str,
    legend: str,
    year: str | None = None,
    country: str | None = None,
    filter_column: str | None = None,
    filter_value: str | None = None,
    start_time: str | None = None,
    end_time: str | None = None,
    hourly: bool = False,
    limit: int = Query(default=24000, ge=500, le=100000),
) -> JSONResponse:
    """Return filtered chart rows and filter options for a result table."""
    paths = _table_paths(dataset, project, scenario, sector, table, year, hourly)
    if not paths:
        raise HTTPException(status_code=404, detail=f"No data found for {table}")
    all_rows = _read_paths(paths)
    dimensions: dict[str, list[str]] = {}
    dimension_names = {"country", legend}
    if filter_column:
        dimension_names.add(filter_column)
    for name in dimension_names:
        dimensions[name] = sorted(
            {str(row[name]) for row in all_rows if row.get(name) not in (None, "")}
        )

    rows = all_rows
    if country and country != "ALL":
        rows = [row for row in rows if str(row.get("country")) == country]
    if filter_column and filter_value and filter_value != "ALL":
        rows = [row for row in rows if str(row.get(filter_column)) == filter_value]

    available_start: str | None = None
    available_end: str | None = None
    if hourly:
        for row in rows:
            if row.get("snapshot") in (None, ""):
                continue
            timestamp = str(row["snapshot"])
            available_start = (
                timestamp
                if available_start is None
                else min(available_start, timestamp)
            )
            available_end = (
                timestamp if available_end is None else max(available_end, timestamp)
            )

    start = _parse_timestamp(start_time, "start") if start_time else None
    end = _parse_timestamp(end_time, "end") if end_time else None
    if start and end and start > end:
        raise HTTPException(
            status_code=400, detail="Start time must be before end time"
        )
    if hourly and (start or end):
        filtered_rows: list[dict[str, Any]] = []
        for row in rows:
            snapshot_value = row.get("snapshot")
            if snapshot_value in (None, ""):
                continue
            try:
                snapshot = datetime.fromisoformat(
                    str(snapshot_value).replace("Z", "+00:00")
                )
            except ValueError:
                continue
            if start and snapshot < start:
                continue
            if end and snapshot > end:
                continue
            filtered_rows.append(row)
        rows = filtered_rows
    configured_chart = next(
        chart
        for charts in CHARTS.values()
        for chart in charts
        if chart["table_name"] == table
    )
    axis_extents = _chart_axis_extents(rows, configured_chart)
    source_count = len(rows)
    rows, truncated = _sample_hourly(rows, legend, limit) if hourly else (rows, False)
    return JSONResponse(
        {
            "rows": rows,
            "dimensions": dimensions,
            "meta": {
                "source_rows": source_count,
                "returned_rows": len(rows),
                "sampled": truncated,
                "files": len(paths),
                "available_start": available_start,
                "available_end": available_end,
                "axis_extents": axis_extents,
            },
        }
    )


@app.get("/api/download")
def download_data(
    dataset: str,
    project: str,
    scenario: str,
    sector: str,
    table: str,
    year: str | None = None,
    hourly: bool = False,
) -> Response:
    """Download a result table as CSV."""
    paths = _table_paths(dataset, project, scenario, sector, table, year, hourly)
    if not paths:
        raise HTTPException(status_code=404, detail=f"No data found for {table}")
    if len(paths) == 1:
        return FileResponse(
            paths[0], media_type="text/csv", filename=f"{scenario}_{table}.csv"
        )
    return Response(
        _combine_csv_files(paths),
        media_type="text/csv",
        headers={
            "Content-Disposition": f'attachment; filename="{scenario}_{table}.csv"'
        },
    )


@app.get("/vendor/plotly.min.js", include_in_schema=False)
def plotly_javascript() -> Response:
    """Serve the bundled Plotly JavaScript library."""
    from plotly.offline import get_plotlyjs

    return Response(get_plotlyjs(), media_type="application/javascript")


@app.get("/brand/pypsa-logo.svg", include_in_schema=False)
def pypsa_logo() -> FileResponse:
    """Serve the PyPSA logo."""
    return FileResponse(ROOT / "pypsa-spice-vis" / "design" / "pypsa-logo_rgb.svg")


@app.get("/brand/flexo.woff2", include_in_schema=False)
def flexo_font() -> FileResponse:
    """Serve the Flexo web font."""
    return FileResponse(
        ROOT / "pypsa-spice-vis" / "design" / "fonts" / "Flexo-Medium.woff2"
    )


app.mount(
    "/ui",
    StaticFiles(directory=FRONTEND_DIST, html=True, check_dir=False),
    name="frontend",
)


def _find_program(program: str) -> str:
    """Return the path of node or npm, or stop with install instructions."""
    path = shutil.which(program)
    if path is None:
        sys.exit(
            f"'{program}' was not found. The web workspace needs Node.js for its "
            "frontend: install it from https://nodejs.org, then start the app again."
        )
    return path


def _run_npm(*npm_args: str) -> None:
    """Run an npm command in the frontend folder; stop the app if it fails."""
    npm_command = [_find_program("npm"), "--prefix", str(FRONTEND_DIR), *npm_args]
    if subprocess.run(npm_command, check=False).returncode != 0:  # noqa: S603
        sys.exit(f"'npm {' '.join(npm_args)}' failed, see the messages above.")


def _frontend_needs_build() -> bool:
    """Return whether the frontend was never built or changed since its last build."""
    build_file = FRONTEND_DIST / "index.html"
    if not build_file.exists():
        return True
    source_files = [
        *(FRONTEND_DIR / "src").rglob("*"),
        *(FRONTEND_DIR / "public").rglob("*"),
        *FRONTEND_DIR.glob("*.*"),
    ]
    last_edit_time = max(
        path.stat().st_mtime for path in source_files if path.is_file()
    )
    return last_edit_time > build_file.stat().st_mtime


def main() -> None:
    """Start the web workspace: ``python pypsa-spice-web/app.py [--dev]``."""
    parser = argparse.ArgumentParser(description="PyPSA-SPICE web workspace.")

    # The --dev option is for development only, not for production deployment.
    # It runs the frontend with `Vite` and the backend with `auto-reload`, so that
    # changes show up immediately without restarting the app.
    parser.add_argument(
        "--dev",
        action="store_true",
        help="for development: serve the frontend with Vite and restart the backend, "
        "so frontend and backend changes show up right away!",
    )
    args = parser.parse_args()

    if not (FRONTEND_DIR / "node_modules").is_dir():
        _run_npm("install")

    # By default, the frontend is built once and served as static files.
    # In development mode, the frontend is served by Vite, which supports hot reloading.
    vite_server = None
    if args.dev:
        vite_command = [_find_program("node"), "node_modules/vite/bin/vite.js"]
        vite_command += ["--host", HOST, "--port", str(VITE_PORT), "--strictPort"]
        # pylint: disable-next=consider-using-with
        vite_server = subprocess.Popen(vite_command, cwd=FRONTEND_DIR)  # noqa: S603
        app_url = f"http://{HOST}:{VITE_PORT}/ui/"
    else:
        if _frontend_needs_build():
            print("Frontend changed, building it (takes a few seconds)...", flush=True)
            _run_npm("run", "build")
        app_url = f"http://{HOST}:{PORT}/ui/"

    print(71 * "-")
    print("Starting PyPSA-SPICE web workspace...\n")
    print(f"Open {app_url} in your browser. Press 'Ctrl+C' to stop.", flush=True)
    print(71 * "-")
    try:
        uvicorn.run(
            "app:app",
            app_dir=str(PACKAGE_DIR),
            host=HOST,
            port=PORT,
            reload=args.dev,
            reload_dirs=[str(PACKAGE_DIR)] if args.dev else None,
        )
    finally:
        if vite_server is not None:
            vite_server.terminate()
            vite_server.wait()


if __name__ == "__main__":
    main()
