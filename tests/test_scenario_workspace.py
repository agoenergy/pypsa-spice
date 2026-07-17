"""Regression tests for local scenario workspace operations."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from pypsa_spice_web.scenario_workspace import (
    ScenarioWorkspace,
    ScenarioWorkspaceError,
)


class ScenarioWorkspaceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.source = self.root / "data" / "dataset" / "project" / "input" / "reference"
        (self.source / "power").mkdir(parents=True)
        (self.source / "scenario_config.yaml").write_text(
            "scenario_configs: {}\n", encoding="utf-8"
        )
        (self.source / "power" / "power_generators.csv").write_text(
            "country,p_nom\nAA,1\n", encoding="utf-8"
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_clone_copies_complete_scenario_without_web_metadata(self) -> None:
        workspace = ScenarioWorkspace(self.root)
        created = workspace.clone(
            dataset="dataset",
            project="project",
            source_scenario="reference",
            new_scenario="policy_high",
        )

        destination = self.source.parent / "policy_high"
        self.assertEqual(created["scenario"], "policy_high")
        self.assertEqual(
            (destination / "power" / "power_generators.csv").read_text(
                encoding="utf-8"
            ),
            "country,p_nom\nAA,1\n",
        )
        self.assertFalse((destination / ".scenario-meta.yaml").exists())

    def test_clone_refuses_overwrite_and_active_web_run(self) -> None:
        workspace = ScenarioWorkspace(self.root)
        workspace.clone(
            dataset="dataset",
            project="project",
            source_scenario="reference",
            new_scenario="copy",
        )
        with self.assertRaises(ScenarioWorkspaceError):
            workspace.clone(
                dataset="dataset",
                project="project",
                source_scenario="reference",
                new_scenario="copy",
            )

        locked = ScenarioWorkspace(
            self.root, active_run=lambda: {"id": "run-1", "status": "running"}
        )
        with self.assertRaisesRegex(ScenarioWorkspaceError, "locked"):
            locked.clone(
                dataset="dataset",
                project="project",
                source_scenario="reference",
                new_scenario="blocked",
            )


if __name__ == "__main__":
    unittest.main()
