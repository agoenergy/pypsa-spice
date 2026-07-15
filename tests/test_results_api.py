"""Regression tests for result-table path handling."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from pypsa_spice_web.app import app


class ResultsApiTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.data_dir = Path(self.temp_dir.name) / "data"
        self.result_dir = (
            self.data_dir
            / "dataset"
            / "project"
            / "results"
            / "run"
            / "csvs"
            / "p-i-t"
            / "2025"
        )
        self.result_dir.mkdir(parents=True)
        (self.result_dir / "pow_gen_by_type_hourly.csv").write_text(
            "snapshot,technology,value\n2025-01-01 00:00,solar,1\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_download_accepts_only_a_real_year_in_the_selected_sector(self) -> None:
        params = {
            "dataset": "dataset",
            "project": "project",
            "scenario": "run",
            "sector": "p-i-t",
            "table": "pow_gen_by_type_hourly",
            "hourly": "true",
        }
        with patch("pypsa_spice_web.app.DATA_DIR", self.data_dir):
            client = TestClient(app)
            valid = client.get("/api/download", params={**params, "year": "2025"})
            escaped = client.get(
                "/api/download",
                params={**params, "year": "../../../../input/scenario/power"},
            )

        self.assertEqual(valid.status_code, 200)
        self.assertEqual(escaped.status_code, 400)
        self.assertEqual(escaped.json()["detail"], "Invalid result year")

    def test_download_rejects_unconfigured_table_names(self) -> None:
        with patch("pypsa_spice_web.app.DATA_DIR", self.data_dir):
            response = TestClient(app).get(
                "/api/download",
                params={
                    "dataset": "dataset",
                    "project": "project",
                    "scenario": "run",
                    "sector": "p-i-t",
                    "table": "../../input/private",
                    "year": "2025",
                    "hourly": "true",
                },
            )

        self.assertEqual(response.status_code, 404)
        self.assertEqual(response.json()["detail"], "Result table is not configured")


if __name__ == "__main__":
    unittest.main()
