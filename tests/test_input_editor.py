"""Regression tests for direct model-input editing."""

from __future__ import annotations

import csv
import shutil
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi import HTTPException
from fastapi.testclient import TestClient

from pypsa_spice_web.app import ROOT, app
from pypsa_spice_web.input_editor import (
    CellChange,
    ConfigSectionUpdate,
    TableUpdate,
    read_scenario_config,
    read_table,
    update_scenario_section,
    update_table,
)


class InputEditorTests(unittest.TestCase):
    """Exercise writes only against temporary file copies."""

    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_input_catalog_endpoint_discovers_projects(self) -> None:
        response = TestClient(app).get("/api/input/catalog")
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json()["table_query_version"], 2)
        self.assertTrue(response.json()["datasets"])
        project = response.json()["datasets"][0]["projects"][0]
        self.assertTrue(project["technologies"])
        technology = project["technologies"][0]
        self.assertTrue(
            {"id", "label", "sector", "classes", "carriers"}.issubset(technology)
        )

    def test_csv_update_is_typed_atomic_and_conflict_checked(self) -> None:
        source = ROOT / "data/example/project_01/input/scenario_01/power/power_generators.csv"
        target = self.temp_path / source.name
        shutil.copy2(source, target)
        config = {"timeseries": False, "filter_col": "type", "with_charts": False}
        before = read_table(target, config, limit=3)

        after = update_table(
            target,
            config,
            TableUpdate(
                revision=before["revision"],
                changes=[CellChange(row=0, column="p_nom", value=12.5)],
            ),
        )
        self.assertNotEqual(before["revision"], after["revision"])
        with target.open(encoding="utf-8", newline="") as handle:
            self.assertEqual(next(csv.DictReader(handle))["p_nom"], "12.5")

        with self.assertRaises(HTTPException) as context:
            update_table(
                target,
                config,
                TableUpdate(
                    revision=before["revision"],
                    changes=[CellChange(row=0, column="p_nom", value=9)],
                ),
            )
        self.assertEqual(context.exception.status_code, 409)

    def test_table_filtering_happens_before_server_pagination(self) -> None:
        target = self.temp_path / "power_generators.csv"
        target.write_text(
            "country,technology,p_nom\n"
            "DE,solar,1\n"
            "FR,wind,2\n"
            "FR,solar,3\n",
            encoding="utf-8",
        )
        config = {
            "timeseries": False,
            "filter_col": "technology",
            "with_charts": False,
        }

        page = read_table(
            target,
            config,
            table="Power_generators",
            technology="solar",
            country="FR",
            offset=0,
            limit=1,
        )

        self.assertEqual(page["total_rows"], 3)
        self.assertEqual(page["total_filtered_rows"], 1)
        self.assertEqual(page["country_options"], ["DE", "FR"])
        self.assertEqual(page["rows"][0]["__row_id"], 2)
        self.assertEqual(page["rows"][0]["p_nom"], 3)

    def test_yaml_section_update_preserves_inline_comments(self) -> None:
        source = ROOT / "data/example/project_01/input/scenario_01/scenario_config.yaml"
        target = self.temp_path / source.name
        shutil.copy2(source, target)
        before = read_scenario_config(target)
        section = dict(before["sections"]["scenario_configs"])
        section["remove_threshold"] = 0.42

        after = update_scenario_section(
            target,
            "scenario_configs",
            ConfigSectionUpdate(revision=before["revision"], value=section),
        )
        self.assertEqual(after["sections"]["scenario_configs"]["remove_threshold"], 0.42)
        contents = target.read_text(encoding="utf-8")
        self.assertIn("# unit: MW", contents)
        self.assertIn("remove_threshold: 0.42", contents)

    def test_direct_write_endpoints_target_the_selected_project(self) -> None:
        data_root = self.temp_path / "data"
        input_root = data_root / "example" / "project_01" / "input"
        global_root = input_root / "global_input"
        scenario_root = input_root / "scenario_01"
        global_root.mkdir(parents=True)
        scenario_root.mkdir()
        csv_target = global_root / "technologies.csv"
        yaml_target = scenario_root / "scenario_config.yaml"
        shutil.copy2(
            ROOT / "data/example/project_01/input/global_input/technologies.csv",
            csv_target,
        )
        shutil.copy2(
            ROOT / "data/example/project_01/input/scenario_01/scenario_config.yaml",
            yaml_target,
        )
        table_params = {
            "dataset": "example",
            "project": "project_01",
            "scenario": "scenario_01",
            "scope": "global",
            "sector": "power",
            "table": "Technologies",
        }
        config_params = {
            "dataset": "example",
            "project": "project_01",
            "scenario": "scenario_01",
        }

        with patch("pypsa_spice_web.app.DATA_DIR", data_root):
            client = TestClient(app)
            table = client.get("/api/input/table", params=table_params)
            self.assertEqual(table.status_code, 200)
            saved_table = client.put(
                "/api/input/table",
                params=table_params,
                json={
                    "revision": table.json()["revision"],
                    "changes": [{"row": 0, "column": "efficiency", "value": 0.987}],
                },
            )
            self.assertEqual(saved_table.status_code, 200)

            config = client.get("/api/input/scenario-config", params=config_params)
            self.assertEqual(config.status_code, 200)
            section = config.json()["sections"]["scenario_configs"]
            section["remove_threshold"] = 0.321
            saved_config = client.put(
                "/api/input/scenario-config/scenario_configs",
                params=config_params,
                json={"revision": config.json()["revision"], "value": section},
            )
            self.assertEqual(saved_config.status_code, 200)

        with csv_target.open(encoding="utf-8", newline="") as handle:
            self.assertEqual(next(csv.DictReader(handle))["efficiency"], "0.987")
        self.assertIn("remove_threshold: 0.321", yaml_target.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
