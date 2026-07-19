"""Run the repository Snakemake workflow from the web workspace."""

from __future__ import annotations

import importlib.util
import json
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import yaml

from pypsa_spice_web.scenario_workspace import (
    describe_data_repository,
    input_file_hashes,
)


class RunConflictError(RuntimeError):
    """Raised when a second model run is requested while one is active."""


class RunValidationError(ValueError):
    """Raised when a requested run does not match the local data tree."""


class ScenarioRunManager:
    """Own one local Snakemake process and expose serialisable run state."""

    TARGET = "solve_all_networks"
    ACTIVE_STATUSES = {"queued", "running", "canceling"}
    STEP_PATTERN = re.compile(
        r"(?P<done>\d+)\s+of\s+(?P<total>\d+)\s+steps"
        r"(?:\s+\((?P<percent>\d+)%\)\s+done)?",
        re.IGNORECASE,
    )
    RULE_PATTERN = re.compile(r"^(?:localrule|rule)\s+([^:]+):", re.MULTILINE)

    def __init__(self, root: Path) -> None:
        self.root = root.resolve()
        self.base_config_path = self.root / "base_config.yaml"
        self.snakefile_path = self.root / "Snakefile"
        self.runs_root = self.root / "logs" / "web_runs"
        self._lock = threading.RLock()
        self._runs: dict[str, dict[str, Any]] = {}
        self._processes: dict[str, subprocess.Popen[Any]] = {}
        self._active_run_id: str | None = None
        self._load_existing_runs()

    @staticmethod
    def _now() -> str:
        return datetime.now(timezone.utc).isoformat()

    @staticmethod
    def _safe_segment(value: str, label: str) -> str:
        cleaned = value.strip()
        if (
            not cleaned
            or cleaned in {".", ".."}
            or "/" in cleaned
            or "\\" in cleaned
            or "\x00" in cleaned
        ):
            raise RunValidationError(f"Invalid {label}")
        return cleaned

    def _load_base_config(self) -> dict[str, Any]:
        if not self.base_config_path.is_file():
            raise RunValidationError("base_config.yaml was not found at the repository root")
        with self.base_config_path.open(encoding="utf-8") as handle:
            config = yaml.safe_load(handle) or {}
        if not isinstance(config.get("path_configs"), dict) or not isinstance(
            config.get("base_configs"), dict
        ):
            raise RunValidationError(
                "base_config.yaml must define path_configs and base_configs"
            )
        return config

    def options(self) -> dict[str, Any]:
        """Return the fixed config reference and model dimensions shown in the UI."""
        config = self._load_base_config()
        paths = config["path_configs"]
        base = config["base_configs"]
        return {
            "config_file": self.base_config_path.name,
            "target": self.TARGET,
            "defaults": {
                "dataset": str(paths.get("data_folder_name", "")),
                "project": str(paths.get("project_name", "")),
                "input_scenario": str(paths.get("input_scenario_name", "")),
                "output_scenario": str(paths.get("output_scenario_name", "")),
            },
            "years": [str(value) for value in base.get("years", [])],
            "sectors": [str(value) for value in base.get("sector", [])],
            "regions": [str(value) for value in (base.get("regions") or {})],
            "currency": str(base.get("currency", "")),
            "environment": "hotpot",
        }

    @staticmethod
    def _hotpot_active() -> bool:
        conda_environment = os.environ.get("CONDA_DEFAULT_ENV", "").strip()
        conda_prefix = Path(os.environ.get("CONDA_PREFIX", "")).name
        return conda_environment == "hotpot" or conda_prefix == "hotpot"

    def _validate_request(
        self,
        dataset: str,
        project: str,
        input_scenario: str,
        output_scenario: str,
    ) -> tuple[str, str, str, str]:
        values = (
            self._safe_segment(dataset, "data folder"),
            self._safe_segment(project, "project"),
            self._safe_segment(input_scenario, "input scenario"),
            self._safe_segment(output_scenario, "output scenario"),
        )
        input_path = self.root / "data" / values[0] / values[1] / "input" / values[2]
        if not input_path.is_dir():
            raise RunValidationError("The selected input scenario does not exist")
        if not (input_path / "scenario_config.yaml").is_file():
            raise RunValidationError(
                "The selected input scenario has no scenario_config.yaml"
            )
        if not self.snakefile_path.is_file():
            raise RunValidationError("The repository Snakefile was not found")
        return values

    def _write_run_config(
        self,
        run_dir: Path,
        dataset: str,
        project: str,
        input_scenario: str,
        output_scenario: str,
    ) -> Path:
        config = self._load_base_config()
        config["path_configs"].update(
            {
                "data_folder_name": dataset,
                "project_name": project,
                "input_scenario_name": input_scenario,
                "output_scenario_name": output_scenario,
            }
        )
        path = run_dir / "base_config.yaml"
        with path.open("w", encoding="utf-8") as handle:
            yaml.safe_dump(config, handle, sort_keys=False)
        return path

    def _resolve_command(self) -> list[str]:
        override = os.environ.get("SNAKEMAKE_COMMAND", "").strip()
        if override:
            command = shlex.split(override)
            if command:
                return command
        if not self._hotpot_active():
            conda = shutil.which("conda")
            if conda:
                return [conda, "run", "--no-capture-output", "-n", "hotpot", "snakemake"]
            raise RunValidationError(
                "Activate the hotpot Conda environment before running the model"
            )
        executable = shutil.which("snakemake")
        if executable:
            return [executable]
        if importlib.util.find_spec("snakemake") is not None:
            return [sys.executable, "-m", "snakemake"]
        raise RunValidationError(
            "Snakemake is not available in the hotpot Conda environment"
        )

    def active(self) -> dict[str, Any] | None:
        """Return a copy of the active web run without affecting CLI processes."""

        with self._lock:
            record = self._active_record()
            return dict(record) if record else None

    def _write_manifest(
        self,
        run_dir: Path,
        *,
        run_id: str,
        dataset: str,
        project: str,
        input_scenario: str,
        output_scenario: str,
        cores: int,
        config_path: Path,
    ) -> Path:
        """Record reproducible input and Git provenance beside the run log."""

        manifest = {
            "schema_version": 1,
            "run_id": run_id,
            "created_at": self._now(),
            "selection": {
                "dataset": dataset,
                "project": project,
                "input_scenario": input_scenario,
                "output_scenario": output_scenario,
            },
            "workflow": {
                "snakefile": str(self.snakefile_path.relative_to(self.root)),
                "target": self.TARGET,
                "cores": cores,
                "environment": "hotpot",
                "source_config": str(self.base_config_path.relative_to(self.root)),
                "run_config": str(config_path.relative_to(self.root)),
            },
            "data_repository": describe_data_repository(
                self.root / "data", dataset
            ),
            "input_sha256": input_file_hashes(
                self.root / "data", dataset, project, input_scenario
            ),
        }
        path = run_dir / "run_manifest.json"
        temporary = path.with_suffix(".tmp")
        temporary.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
        temporary.replace(path)
        return path

    @staticmethod
    def _signal_process(process: subprocess.Popen[Any], sig: signal.Signals) -> None:
        """Signal Snakemake and its child jobs when process groups are available."""
        if os.name == "posix":
            try:
                os.killpg(process.pid, sig)
                return
            except ProcessLookupError:
                return
            except OSError:
                pass
        if sig == signal.SIGKILL:
            process.kill()
        else:
            process.terminate()

    def _force_kill_later(self, run_id: str, process: subprocess.Popen[Any]) -> None:
        def force_kill() -> None:
            with self._lock:
                if self._processes.get(run_id) is process and process.poll() is None:
                    self._signal_process(process, signal.SIGKILL)

        timer = threading.Timer(10, force_kill)
        timer.daemon = True
        timer.start()

    def _state_path(self, run_id: str) -> Path:
        return self.runs_root / run_id / "state.json"

    def _persist(self, record: dict[str, Any]) -> None:
        path = self._state_path(record["id"])
        path.parent.mkdir(parents=True, exist_ok=True)
        temporary = path.with_suffix(".tmp")
        temporary.write_text(
            json.dumps(record, indent=2, sort_keys=True), encoding="utf-8"
        )
        temporary.replace(path)

    def _load_existing_runs(self) -> None:
        if not self.runs_root.is_dir():
            return
        for path in self.runs_root.glob("*/state.json"):
            try:
                record = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if not isinstance(record, dict) or not record.get("id"):
                continue
            if record.get("status") in self.ACTIVE_STATUSES:
                record.update(
                    {
                        "status": "failed",
                        "ended_at": self._now(),
                        "message": "The web service stopped before this run finished.",
                    }
                )
                self._persist(record)
            self._runs[str(record["id"])] = record

    def _active_record(self) -> dict[str, Any] | None:
        if not self._active_run_id:
            return None
        record = self._runs.get(self._active_run_id)
        if not record or record.get("status") not in self.ACTIVE_STATUSES:
            self._active_run_id = None
            return None
        return record

    def start(
        self,
        *,
        dataset: str,
        project: str,
        input_scenario: str,
        output_scenario: str,
        cores: int = 1,
    ) -> dict[str, Any]:
        """Validate a model selection, create its config copy, and queue Snakemake."""
        dataset, project, input_scenario, output_scenario = self._validate_request(
            dataset, project, input_scenario, output_scenario
        )
        if cores < 1 or cores > 128:
            raise RunValidationError("Cores must be between 1 and 128")

        with self._lock:
            active = self._active_record()
            if active:
                raise RunConflictError(
                    f"Run {active['id']} is already {active['status']}"
                )
            run_id = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:8]
            run_dir = self.runs_root / run_id
            run_dir.mkdir(parents=True, exist_ok=False)
            try:
                config_path = self._write_run_config(
                    run_dir, dataset, project, input_scenario, output_scenario
                )
                manifest_path = self._write_manifest(
                    run_dir,
                    run_id=run_id,
                    dataset=dataset,
                    project=project,
                    input_scenario=input_scenario,
                    output_scenario=output_scenario,
                    cores=cores,
                    config_path=config_path,
                )
            except Exception:
                shutil.rmtree(run_dir, ignore_errors=True)
                raise
            record: dict[str, Any] = {
                "id": run_id,
                "status": "queued",
                "progress": 0,
                "message": "Preparing Snakemake workflow",
                "current_rule": None,
                "created_at": self._now(),
                "started_at": None,
                "ended_at": None,
                "exit_code": None,
                "pid": None,
                "dataset": dataset,
                "project": project,
                "input_scenario": input_scenario,
                "output_scenario": output_scenario,
                "cores": cores,
                "target": self.TARGET,
                "config_file": str(config_path.relative_to(self.root)),
                "log_file": str((run_dir / "snakemake.log").relative_to(self.root)),
                "manifest_file": str(manifest_path.relative_to(self.root)),
            }
            self._runs[run_id] = record
            self._active_run_id = run_id
            self._persist(record)
            threading.Thread(
                target=self._execute,
                args=(run_id, config_path),
                name=f"snakemake-{run_id}",
                daemon=True,
            ).start()
            return self.get(run_id)

    def _execute(self, run_id: str, config_path: Path) -> None:
        log_path = self.runs_root / run_id / "snakemake.log"
        try:
            with self._lock:
                if self._runs[run_id].get("status") == "canceling":
                    self._runs[run_id].update(
                        {
                            "status": "canceled",
                            "ended_at": self._now(),
                            "message": "Run canceled before Snakemake started",
                        }
                    )
                    self._persist(self._runs[run_id])
                    return
            command = [
                *self._resolve_command(),
                "--snakefile",
                str(self.snakefile_path),
                "--configfile",
                str(config_path),
                "--cores",
                str(self._runs[run_id]["cores"]),
                "--printshellcmds",
                "--rerun-incomplete",
                self.TARGET,
            ]
            environment = os.environ.copy()
            environment["TMPDIR"] = str(self.runs_root / run_id / "tmp")
            Path(environment["TMPDIR"]).mkdir(parents=True, exist_ok=True)
            with log_path.open("a", encoding="utf-8", buffering=1) as log:
                log.write(f"$ {shlex.join(command)}\n\n")
                process = subprocess.Popen(
                    command,
                    cwd=self.root,
                    env=environment,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    start_new_session=os.name == "posix",
                )
                with self._lock:
                    self._processes[run_id] = process
                    record = self._runs[run_id]
                    if record.get("status") == "canceling":
                        self._signal_process(process, signal.SIGTERM)
                        return_code = process.wait()
                        record.update(
                            {
                                "status": "canceled",
                                "ended_at": self._now(),
                                "exit_code": return_code,
                                "message": "Run canceled",
                            }
                        )
                        self._persist(record)
                        return
                    record.update(
                        {
                            "status": "running",
                            "started_at": self._now(),
                            "pid": process.pid,
                            "message": "Snakemake workflow is running",
                        }
                    )
                    self._persist(record)
                exit_code = process.wait()

            with self._lock:
                record = self._runs[run_id]
                canceled = record.get("status") == "canceling"
                record.update(
                    {
                        "status": "canceled"
                        if canceled
                        else "succeeded"
                        if exit_code == 0
                        else "failed",
                        "progress": 100
                        if exit_code == 0 and not canceled
                        else record.get("progress", 0),
                        "message": "Run canceled"
                        if canceled
                        else "Model run completed"
                        if exit_code == 0
                        else f"Snakemake exited with code {exit_code}",
                        "ended_at": self._now(),
                        "exit_code": exit_code,
                    }
                )
                self._refresh_progress(record)
                self._persist(record)
        except Exception as exc:  # keep background failures visible to the UI
            with self._lock:
                record = self._runs[run_id]
                record.update(
                    {
                        "status": "failed",
                        "message": str(exc),
                        "ended_at": self._now(),
                    }
                )
                self._persist(record)
            try:
                with log_path.open("a", encoding="utf-8") as log:
                    log.write(f"\nWeb runner error: {exc}\n")
            except OSError:
                pass
        finally:
            with self._lock:
                self._processes.pop(run_id, None)
                if self._active_run_id == run_id:
                    self._active_run_id = None

    @staticmethod
    def _tail_text(path: Path, max_bytes: int = 96_000) -> str:
        if not path.is_file():
            return ""
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - max_bytes))
            data = handle.read()
        return data.decode("utf-8", errors="replace")

    def _refresh_progress(self, record: dict[str, Any]) -> str:
        log = self._tail_text(self.root / record["log_file"])
        matches = list(self.STEP_PATTERN.finditer(log))
        if matches:
            match = matches[-1]
            done = int(match.group("done"))
            total = max(1, int(match.group("total")))
            percent = int(match.group("percent") or round(done / total * 100))
            record["progress"] = max(int(record.get("progress", 0)), percent)
            if record.get("status") in self.ACTIVE_STATUSES:
                record["message"] = f"{done} of {total} workflow steps completed"
        rules = self.RULE_PATTERN.findall(log)
        if rules:
            record["current_rule"] = rules[-1].strip()
        return log

    def get(self, run_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._runs.get(run_id)
            if record is None:
                state_path = self._state_path(run_id)
                if not state_path.is_file():
                    raise KeyError(run_id)
                record = json.loads(state_path.read_text(encoding="utf-8"))
                self._runs[run_id] = record
            log = self._refresh_progress(record)
            response = dict(record)
            response["log"] = log
            return response

    def latest(self) -> dict[str, Any] | None:
        with self._lock:
            if not self._runs:
                return None
            record = max(
                self._runs.values(), key=lambda item: str(item.get("created_at", ""))
            )
            return self.get(str(record["id"]))

    def cancel(self, run_id: str) -> dict[str, Any]:
        with self._lock:
            record = self._runs.get(run_id)
            if record is None:
                raise KeyError(run_id)
            if record.get("status") not in self.ACTIVE_STATUSES:
                return self.get(run_id)
            record.update(
                {"status": "canceling", "message": "Stopping Snakemake workflow"}
            )
            self._persist(record)
            process = self._processes.get(run_id)
            if process and process.poll() is None:
                self._signal_process(process, signal.SIGTERM)
                self._force_kill_later(run_id, process)
            return self.get(run_id)
