#!/usr/bin/env python3
"""Repository-level structural checks that do not execute scientific analyses."""

from __future__ import annotations

import csv
import json
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


class StaticPackageTests(unittest.TestCase):
    def test_active_code_has_no_project_absolute_paths(self) -> None:
        pattern = re.compile(r"/vol2/|/storage/|/public/|(?<![A-Za-z0-9_])[A-Za-z]:[\\/]")
        offenders: list[str] = []
        for path in (ROOT / "src").rglob("*"):
            if path.suffix not in {".py", ".sh"}:
                continue
            if pattern.search(path.read_text(encoding="utf-8")):
                offenders.append(str(path.relative_to(ROOT)))
        self.assertEqual(offenders, [])

    def test_notebook_is_valid_and_has_no_saved_outputs(self) -> None:
        path = ROOT / "notebooks" / "Fig8_final_plotting.ipynb"
        notebook = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(notebook["nbformat"], 4)
        for cell in notebook["cells"]:
            if cell["cell_type"] == "code":
                self.assertIsNone(cell["execution_count"])
                self.assertEqual(cell["outputs"], [])

    def test_manifest_contracts(self) -> None:
        contracts = {
            "input_manifest.tsv": ["config_key", "path_role", "required_for", "status", "note"],
            "output_contract.tsv": ["stage", "output_relative_to_run", "type", "required", "note"],
        }
        for name, expected in contracts.items():
            with (ROOT / "manifests" / name).open(encoding="utf-8", newline="") as handle:
                reader = csv.reader(handle, delimiter="\t")
                self.assertEqual(next(reader), expected, name)
                for index, row in enumerate(reader, start=2):
                    self.assertEqual(len(row), len(expected), f"{name}:{index}")

    def test_cts_header_fix_is_present(self) -> None:
        prepare = (ROOT / "src" / "gwas_ldsc" / "prepare_E1E9_reuse_existing_recalc.sh").read_text(encoding="utf-8")
        submit = (ROOT / "src" / "gwas_ldsc" / "submit_E1E9_reuse_existing_152.sh").read_text(encoding="utf-8")
        worker = (ROOT / "src" / "gwas_ldsc" / "run_A2_B2_cts_152.sh").read_text(encoding="utf-8")
        self.assertIn("ldcts_name\\ttrait_id\\tsumstats\\tldcts_file\\toutdir", prepare)
        self.assertIn("manifest header mismatch", submit)
        self.assertIn("GLOBAL_TASK+2", worker)


if __name__ == "__main__":
    unittest.main()
