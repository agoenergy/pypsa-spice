"""FastAPI application for the PyPSA-SPICE output explorer."""

from __future__ import annotations

import csv
import io
import math
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml
from fastapi import FastAPI, HTTPException, Query, Request
from fastapi.middleware.gzip import GZipMiddleware
from fastapi.responses import FileResponse, HTMLResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles

ROOT = Path(__file__).resolve().parents[1]
PACKAGE_DIR = Path(__file__).resolve().parent
FRONTEND_DIST = PACKAGE_DIR / "frontend" / "dist"
DATA_DIR = ROOT / "data"
GRAPH_CONFIG = ROOT / "pypsa-spice-vis" / "setting" / "graph_settings.yaml"
MAPPING_DIR = ROOT / "pypsa-spice-vis" / "setting"

CHART_TYPES = {
    "p1": "bar",
    "p2": "filtered_bar",
    "p3": "bar",
    "p4": "area_share",
    "p6": "filtered_bar",
    "p7": "hourly_bar",
    "p8": "filtered_hourly_bar",
    "p9": "bar",
    "p10": "hourly_bar",
    "p11": "hourly_line",
    "p12": "hourly_line",
    "p13": "bar",
    "p14": "hourly_dual",
    "i1": "bar",
    "i2": "filtered_bar",
    "i3": "filtered_bar",
    "i4": "bar",
    "t1": "hourly_dual",
    "e1": "bar",
    "e2": "bar",
    "e3": "bar",
    "c1": "bar",
    "c2": "bar",
    "c3": "bar",
    "c4": "bar",
    "c5": "bar",
    "c6": "bar",
    "c7": "bar",
    "c8": "bar",
}

SECTION_META = {
    "power": {
        "label": "Power",
        "icon": "ϟ",
        "title": "Power Sector",
        "eyebrow": "Power system",
    },
    "industry": {
        "label": "Industry",
        "icon": "▦",
        "title": "Industry Sector",
        "eyebrow": "Industrial energy system",
    },
    "transport": {
        "label": "Transport",
        "icon": "→",
        "title": "Transport Sector",
        "eyebrow": "Transport energy system",
    },
    "emissions": {
        "label": "Emissions",
        "icon": "◌",
        "title": "Emissions",
        "eyebrow": "Cross-sector emissions",
    },
    "costs": {
        "label": "Costs",
        "icon": "$",
        "title": "System Costs",
        "eyebrow": "Cross-sector costs",
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
            chart["type"] = CHART_TYPES.get(chart_id, "bar")
            chart["hourly"] = "hourly" in chart["table_name"]
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
    root = (_scenario_root(dataset, project, scenario) / sector).resolve()
    if root.parent != _scenario_root(dataset, project, scenario) or not root.is_dir():
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
    year_dirs = sorted(
        [path for path in sector_path.iterdir() if path.is_dir()],
        key=lambda path: _year_sort(path.name),
    )
    if hourly:
        selected = year or next((p.name for p in year_dirs if p.name.isdigit()), None)
        paths = [sector_path / selected / f"{table}.csv"] if selected else []
    else:
        all_years = sector_path / "all_years" / f"{table}.csv"
        paths = (
            [all_years]
            if all_years.exists()
            else [p / f"{table}.csv" for p in year_dirs]
        )
    return [path for path in paths if path.is_file()]


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
    per_group = max(30, limit // max(1, len(groups)))
    sampled: list[dict[str, Any]] = []
    for group in groups.values():
        step = max(1, math.ceil(len(group) / per_group))
        sampled.extend(group[::step][:per_group])
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
                        files = {path.name for path in sector_path.glob("*/*.csv")}
                        section_counts = {
                            section: sum(
                                1
                                for chart in charts
                                if f"{chart['table_name']}.csv" in files
                            )
                            for section, charts in CHARTS.items()
                        }
                        chart_count = sum(section_counts.values())
                        sectors.append(
                            {
                                "name": sector_path.name,
                                "years": [year for year in years if year.isdigit()],
                                "chart_count": chart_count,
                                "section_counts": section_counts,
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
    title="PyPSA-SPICE Output Explorer",
    description="Interactive browser for PyPSA-SPICE result tables.",
    version="0.1.0",
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


@app.get("/api/catalog")
def catalog() -> JSONResponse:
    return JSONResponse(_catalog())


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
