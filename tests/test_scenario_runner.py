"""Regression tests for web-launched Snakemake model runs."""

from __future__ import annotations

import sys
import json
import tempfile
import time
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml
from fastapi.testclient import TestClient

from pypsa_spice_web.app import app
from pypsa_spice_web.scenario_runner import ScenarioRunManager


class ScenarioRunManagerTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)
        (self.root / "Snakefile").write_text(
            'configfile: "base_config.yaml"\nrule solve_all_networks:\n    input: []\n',
            encoding="utf-8",
        )
        self.base_config = {
            "version": "0.1.0",
            "path_configs": {
                "data_folder_name": "old-data",
                "project_name": "old-project",
                "input_scenario_name": "old-input",
                "output_scenario_name": "old-output",
            },
            "base_configs": {
                "regions": {"AA": ["one"]},
                "years": [2025, 2030],
                "sector": ["p-i-t"],
                "currency": "EUR",
            },
        }
        (self.root / "base_config.yaml").write_text(
            yaml.safe_dump(self.base_config, sort_keys=False), encoding="utf-8"
        )
        scenario = self.root / "data" / "dataset" / "project" / "input" / "scenario"
        scenario.mkdir(parents=True)
        (scenario / "scenario_config.yaml").write_text("scenario_configs: {}\n", encoding="utf-8")
        self.manager = ScenarioRunManager(self.root)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_options_are_read_from_repository_base_config(self) -> None:
        options = self.manager.options()
        self.assertEqual(options["config_file"], "base_config.yaml")
        self.assertEqual(options["target"], "solve_all_networks")
        self.assertEqual(options["years"], ["2025", "2030"])
        self.assertEqual(options["regions"], ["AA"])

    def test_run_uses_an_overridden_copy_and_preserves_source_config(self) -> None:
        fake_snakemake = (
            "print('localrule build_network:', flush=True);"
            "print('1 of 2 steps (50%) done', flush=True);"
            "print('2 of 2 steps (100%) done', flush=True)"
        )
        with patch.object(
            self.manager,
            "_resolve_command",
            return_value=[sys.executable, "-c", fake_snakemake],
        ):
            run = self.manager.start(
                dataset="dataset",
                project="project",
                input_scenario="scenario",
                output_scenario="scenario_web_run",
                cores=2,
            )
            deadline = time.monotonic() + 5
            while run["status"] in self.manager.ACTIVE_STATUSES and time.monotonic() < deadline:
                time.sleep(0.05)
                run = self.manager.get(run["id"])

        self.assertEqual(run["status"], "succeeded")
        self.assertEqual(run["progress"], 100)
        self.assertEqual(run["current_rule"], "build_network")
        self.assertIn("2 of 2 steps", run["log"])

        run_config = yaml.safe_load((self.root / run["config_file"]).read_text(encoding="utf-8"))
        self.assertEqual(
            run_config["path_configs"],
            {
                "data_folder_name": "dataset",
                "project_name": "project",
                "input_scenario_name": "scenario",
                "output_scenario_name": "scenario_web_run",
            },
        )
        self.assertEqual(run_config["base_configs"], self.base_config["base_configs"])
        source = yaml.safe_load((self.root / "base_config.yaml").read_text(encoding="utf-8"))
        self.assertEqual(source, self.base_config)
        manifest = json.loads(
            (self.root / run["manifest_file"]).read_text(encoding="utf-8")
        )
        self.assertEqual(manifest["workflow"]["environment"], "hotpot")
        self.assertEqual(manifest["selection"]["input_scenario"], "scenario")
        self.assertIn(
            "input/scenario/scenario_config.yaml", manifest["input_sha256"]
        )

    def test_run_api_starts_and_reports_a_managed_workflow(self) -> None:
        fake_snakemake = "print('1 of 1 steps (100%) done', flush=True)"
        with (
            patch("pypsa_spice_web.app.RUN_MANAGER", self.manager),
            patch.object(
                self.manager,
                "_resolve_command",
                return_value=[sys.executable, "-c", fake_snakemake],
            ),
        ):
            client = TestClient(app)
            response = client.post(
                "/api/runs",
                json={
                    "dataset": "dataset",
                    "project": "project",
                    "input_scenario": "scenario",
                    "output_scenario": "api-run",
                    "cores": 1,
                },
            )
            self.assertEqual(response.status_code, 202)
            run = response.json()
            deadline = time.monotonic() + 5
            while run["status"] in self.manager.ACTIVE_STATUSES and time.monotonic() < deadline:
                time.sleep(0.05)
                run = client.get(f"/api/runs/{run['id']}").json()

        self.assertEqual(run["status"], "succeeded")
        self.assertEqual(run["output_scenario"], "api-run")
        self.assertEqual(run["progress"], 100)


if __name__ == "__main__":
    unittest.main()
