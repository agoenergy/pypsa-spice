"""Local scenario workspace operations for the React/FastAPI application.

The model still reads its normal ``data/<dataset>/<project>/input`` tree.  This
module only adds safe, atomic workspace conveniences around that existing
layout; it does not introduce a web-only scenario format.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Callable


class ScenarioWorkspaceError(ValueError):
    """Raised when a local scenario workspace operation is invalid."""


def _safe_segment(value: str, label: str) -> str:
    cleaned = value.strip()
    if (
        not cleaned
        or cleaned in {".", "..", "global_input"}
        or "/" in cleaned
        or "\\" in cleaned
        or "\x00" in cleaned
        or cleaned.startswith(".")
    ):
        raise ScenarioWorkspaceError(f"Invalid {label}")
    return cleaned


def _confined(candidate: Path, parent: Path) -> Path:
    resolved_parent = parent.resolve()
    resolved = candidate.resolve()
    try:
        resolved.relative_to(resolved_parent)
    except ValueError as exc:
        raise ScenarioWorkspaceError("Invalid scenario path") from exc
    return resolved


def _git_value(repo: Path, *arguments: str) -> str:
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), *arguments],
            check=False,
            capture_output=True,
            text=True,
            timeout=5,
        )
    except (OSError, subprocess.TimeoutExpired):
        return ""
    return result.stdout.strip() if result.returncode == 0 else ""


def describe_data_repository(data_dir: Path, dataset: str) -> dict[str, Any]:
    """Describe a dataset-level Git worktree without changing it."""

    dataset_name = _safe_segment(dataset, "data folder")
    dataset_root = _confined(data_dir / dataset_name, data_dir)
    git_marker = dataset_root / ".git"
    if not git_marker.exists():
        return {
            "available": False,
            "root": str(dataset_root.relative_to(data_dir.parent)),
            "remote": "",
            "branch": "",
            "commit": "",
            "dirty": False,
            "changes": [],
        }

    changes = [line for line in _git_value(dataset_root, "status", "--short").splitlines() if line]
    return {
        "available": True,
        "root": str(dataset_root.relative_to(data_dir.parent)),
        "remote": _git_value(dataset_root, "remote", "get-url", "origin"),
        "branch": _git_value(dataset_root, "branch", "--show-current"),
        "commit": _git_value(dataset_root, "rev-parse", "--short=12", "HEAD"),
        "dirty": bool(changes),
        "changes": changes,
    }


def input_file_hashes(
    data_dir: Path, dataset: str, project: str, scenario: str
) -> dict[str, str]:
    """Hash the selected scenario and shared global inputs for run provenance."""

    dataset_name = _safe_segment(dataset, "data folder")
    project_name = _safe_segment(project, "project")
    scenario_name = _safe_segment(scenario, "scenario")
    project_root = _confined(data_dir / dataset_name / project_name, data_dir)
    input_root = _confined(project_root / "input", project_root)
    hashes: dict[str, str] = {}
    for folder_name in ("global_input", scenario_name):
        folder = _confined(input_root / folder_name, input_root)
        if not folder.is_dir():
            continue
        for path in sorted(item for item in folder.rglob("*") if item.is_file()):
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                while chunk := handle.read(1024 * 1024):
                    digest.update(chunk)
            hashes[str(path.relative_to(project_root))] = digest.hexdigest()
    return hashes


class ScenarioWorkspace:
    """Create local scenarios and expose their data-repository state."""

    def __init__(
        self,
        root: Path,
        *,
        active_run: Callable[[], dict[str, Any] | None] | None = None,
    ) -> None:
        self.root = root.resolve()
        self.data_dir = self.root / "data"
        self._active_run = active_run or (lambda: None)

    def active_run(self) -> dict[str, Any] | None:
        return self._active_run()

    def ensure_mutable(self) -> None:
        active = self.active_run()
        if active:
            raise ScenarioWorkspaceError(
                f"Scenario files are locked while run {active['id']} is {active['status']}"
            )

    def status(self, dataset: str) -> dict[str, Any]:
        active = self.active_run()
        return {
            "mutation_locked": bool(active),
            "active_run": {
                "id": active["id"],
                "status": active["status"],
            }
            if active
            else None,
            "repository": describe_data_repository(self.data_dir, dataset),
        }

    def clone(
        self,
        *,
        dataset: str,
        project: str,
        source_scenario: str,
        new_scenario: str,
    ) -> dict[str, Any]:
        """Atomically duplicate a complete existing scenario folder."""

        self.ensure_mutable()
        dataset_name = _safe_segment(dataset, "data folder")
        project_name = _safe_segment(project, "project")
        source_name = _safe_segment(source_scenario, "source scenario")
        destination_name = _safe_segment(new_scenario, "new scenario")
        project_root = _confined(
            self.data_dir / dataset_name / project_name, self.data_dir
        )
        input_root = _confined(project_root / "input", project_root)
        if not input_root.is_dir():
            raise ScenarioWorkspaceError("The selected project input folder is missing")
        source = _confined(input_root / source_name, input_root)
        destination = _confined(input_root / destination_name, input_root)
        if not source.is_dir() or not (source / "scenario_config.yaml").is_file():
            raise ScenarioWorkspaceError("The source scenario is incomplete or missing")
        if destination.exists():
            raise ScenarioWorkspaceError("A scenario with this name already exists")

        temporary = Path(
            tempfile.mkdtemp(prefix=f".{destination_name}-", dir=input_root)
        )
        try:
            shutil.copytree(source, temporary, dirs_exist_ok=True)
            os.replace(temporary, destination)
        except Exception:
            shutil.rmtree(temporary, ignore_errors=True)
            raise

        return {
            "dataset": dataset_name,
            "project": project_name,
            "scenario": destination_name,
            "source_scenario": source_name,
            "path": str(destination.relative_to(self.root)),
            "repository": describe_data_repository(self.data_dir, dataset_name),
        }
