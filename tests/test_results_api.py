"""Regression tests for result-table path handling."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from fastapi.testclient import TestClient

from pypsa_spice_web.app import GRAPH_CONFIG, PACKAGE_DIR, app


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
        (self.result_dir / "pow_flh_by_type_yearly.csv").write_text(
            "country,technology,year,value,unit\n"
            "DE,solar,2025,1450,hours\n"
            "DE,wind,2025,2800,hours\n",
            encoding="utf-8",
        )
        (self.result_dir / "pow_cap_by_region_yearly.csv").write_text(
            "country,region,year,value,unit\n"
            "DE,North,2025,30.0,GW\n"
            "DE,South,2025,20.0,GW\n",
            encoding="utf-8",
        )
        (self.result_dir / "pow_inter_cf_by_region_yearly.csv").write_text(
            "country,from,to,year,value,unit\n"
            "DE,North,South,2025,0.75,Dimensionless\n"
            "DE,West,South,2025,0.5,Dimensionless\n",
            encoding="utf-8",
        )

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def test_web_uses_package_local_graph_settings(self) -> None:
        self.assertEqual(GRAPH_CONFIG, PACKAGE_DIR / "graph_settings.yaml")
        self.assertTrue(GRAPH_CONFIG.is_file())

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

    def test_full_load_hours_chart_reads_yearly_result_table(self) -> None:
        with patch("pypsa_spice_web.app.DATA_DIR", self.data_dir):
            response = TestClient(app).get(
                "/api/chart",
                params={
                    "dataset": "dataset",
                    "project": "project",
                    "scenario": "run",
                    "sector": "p-i-t",
                    "table": "pow_flh_by_type_yearly",
                    "legend": "technology",
                    "country": "DE",
                    "hourly": "false",
                },
            )

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["dimensions"]["technology"], ["solar", "wind"])
        self.assertEqual([row["value"] for row in payload["rows"]], [1450.0, 2800.0])

    def test_additional_yearly_charts_read_result_tables(self) -> None:
        cases = (
            ("pow_cap_by_region_yearly", "region", None, [30.0, 20.0]),
            ("pow_inter_cf_by_region_yearly", "from", "to", [0.75, 0.5]),
        )
        with patch("pypsa_spice_web.app.DATA_DIR", self.data_dir):
            client = TestClient(app)
            for table, legend, filter_column, expected_values in cases:
                with self.subTest(table=table):
                    params = {
                        "dataset": "dataset",
                        "project": "project",
                        "scenario": "run",
                        "sector": "p-i-t",
                        "table": table,
                        "legend": legend,
                        "country": "DE",
                        "hourly": "false",
                    }
                    if filter_column:
                        params["filter_column"] = filter_column
                        params["filter_value"] = "South"
                    response = client.get("/api/chart", params=params)

                    self.assertEqual(response.status_code, 200)
                    self.assertEqual(
                        [row["value"] for row in response.json()["rows"]],
                        expected_values,
                    )


if __name__ == "__main__":
    unittest.main()
