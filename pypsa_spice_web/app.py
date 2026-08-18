"""FastAPI application for the PyPSA-SPICE web workspace."""

from __future__ import annotations

import csv
import io
import math
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Literal

import yaml
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from pypsa_spice_web.input_editor import (
    ConfigSectionUpdate,
    ConfigSectionsUpdate,
    TableUpdate,
    compare_scenarios,
    input_catalog,
    read_scenario_config,
    read_table,
    scenario_config_path,
    table_path,
    update_scenario_section,
    update_scenario_sections,
    update_table,
)
from pypsa_spice_web.scenario_runner import (
    RunConflictError,
    RunValidationError,
    ScenarioRunManager,
)
from pypsa_spice_web.scenario_workspace import (
    ScenarioWorkspace,
    ScenarioWorkspaceError,
)

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = Path(__file__).resolve().parent
FRONTEND_DIST = PACKAGE_DIR / "frontend" / "dist"
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


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _read_mapping(filename: str) -> dict[str, dict[str, str]]:
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
    raw = _load_yaml(GRAPH_CONFIG)
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
    candidate = (DATA_DIR / dataset / project / "results" / scenario / "csvs").resolve()
    try:
        candidate.relative_to(DATA_DIR.resolve())
    except ValueError as exc:
        raise HTTPException(status_code=400, detail="Invalid result path") from exc
    if not candidate.is_dir():
        raise HTTPException(status_code=404, detail="Result scenario not found")
    return candidate


def _sector_root(dataset: str, project: str, scenario: str, sector: str) -> Path:
    scenario_root = _scenario_root(dataset, project, scenario)
    root = (scenario_root / sector).resolve()
    if root.parent != scenario_root or not root.is_dir():
        raise HTTPException(status_code=404, detail="Sector result not found")
    return root


def _year_sort(value: str) -> tuple[int, str]:
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


def _read_paths(paths: list[Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for path in paths:
        rows.extend(
            dict(row) for row in _read_csv_cached(str(path), path.stat().st_mtime_ns)
        )
    # Some exporters repeat the full yearly table in each year folder.
    seen: set[tuple[Any, ...]] = set()
    unique_rows: list[dict[str, Any]] = []
    for row in rows:
        marker = tuple(row.items())
        if marker not in seen:
            seen.add(marker)
            unique_rows.append(row)
    return unique_rows


def _sample_hourly(
    rows: list[dict[str, Any]], legend_column: str, limit: int
) -> tuple[list[dict[str, Any]], bool]:
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


def _parse_timestamp(value: str, field_name: str) -> datetime:
    """Parse an ISO timestamp supplied by the web chart controls."""
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise HTTPException(
            status_code=400, detail=f"Invalid {field_name} timestamp"
        ) from exc


def _catalog() -> dict[str, Any]:
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
    title="PyPSA-SPICE Web Workspace",
    description="Explore model results and directly edit PyPSA-SPICE CSV/YAML inputs.",
    version="0.2.0",
)
app.add_middleware(GZipMiddleware, minimum_size=1000)


@app.get("/", response_class=HTMLResponse)
def index(request: Request) -> Response:
    del request
    index_file = FRONTEND_DIST / "index.html"
    if index_file.exists():
        return FileResponse(index_file)
    return HTMLResponse(
        "<h1>Frontend not built</h1><p>Run <code>npm run web:build</code> "
        "and restart the server.</p>",
        status_code=503,
    )


@app.get("/health")
def health() -> dict[str, str]:
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

    available_timestamps = sorted(
        str(row["snapshot"])
        for row in rows
        if hourly and row.get("snapshot") not in (None, "")
    )
    available_start = available_timestamps[0] if available_timestamps else None
    available_end = available_timestamps[-1] if available_timestamps else None

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
    paths = _table_paths(dataset, project, scenario, sector, table, year, hourly)
    if not paths:
        raise HTTPException(status_code=404, detail=f"No data found for {table}")
    if len(paths) == 1:
        return FileResponse(
            paths[0], media_type="text/csv", filename=f"{scenario}_{table}.csv"
        )
    rows = _read_paths(paths)
    if not rows:
        raise HTTPException(status_code=404, detail="The result table is empty")
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return Response(
        buffer.getvalue(),
        media_type="text/csv",
        headers={
            "Content-Disposition": f'attachment; filename="{scenario}_{table}.csv"'
        },
    )


@app.get("/vendor/plotly.min.js", include_in_schema=False)
def plotly_javascript() -> Response:
    from plotly.offline import get_plotlyjs

    return Response(get_plotlyjs(), media_type="application/javascript")


@app.get("/brand/pypsa-logo.svg", include_in_schema=False)
def pypsa_logo() -> FileResponse:
    return FileResponse(ROOT / "pypsa-spice-vis" / "design" / "pypsa-logo_rgb.svg")


@app.get("/brand/flexo.woff2", include_in_schema=False)
def flexo_font() -> FileResponse:
    return FileResponse(
        ROOT / "pypsa-spice-vis" / "design" / "fonts" / "Flexo-Medium.woff2"
    )


app.mount(
    "/ui",
    StaticFiles(directory=FRONTEND_DIST, check_dir=False),
    name="frontend",
)
