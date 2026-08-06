"""Load the repository path table and expose values as environment variables."""

from __future__ import annotations

import csv
import os
from pathlib import Path


def load_paths(path: str | Path) -> dict[str, str]:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Path configuration not found: {path}")

    values: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"key", "value"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain tab-separated columns: key, value")
        for row in reader:
            key = (row.get("key") or "").strip()
            value = (row.get("value") or "").strip()
            if not key or key.startswith("#") or not value:
                continue
            if value.startswith("<") and value.endswith(">"):
                continue
            values[key] = value
    return values


def export_paths(path: str | Path) -> dict[str, str]:
    values = load_paths(path)
    os.environ.update(values)
    return values
