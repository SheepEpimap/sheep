#!/usr/bin/env python3
"""Create a configured repository copy by replacing legacy path prefixes."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


TEXT_SUFFIXES = {
    ".csv", ".env", ".ipynb", ".md", ".py", ".r", ".sbatch", ".sh",
    ".smk", ".toml", ".tsv", ".txt", ".yaml", ".yml",
}


def read_mapping(path: Path) -> list[tuple[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    mapping = []
    for index, row in enumerate(rows, start=2):
        source = (row.get("original_prefix") or "").rstrip("/")
        target = (row.get("replacement_prefix") or "").rstrip("/")
        if not source or not target or target.startswith("/path/to/"):
            raise ValueError(f"Complete both mapping columns at {path}:{index}")
        mapping.append((source, target))
    if not mapping:
        raise ValueError(f"No path mappings found in {path}")
    return sorted(mapping, key=lambda item: len(item[0]), reverse=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    if source == output or source in output.parents:
        raise ValueError("Output must be outside the source repository")

    mapping = read_mapping(args.mapping.resolve())
    shutil.copytree(source, output, ignore=shutil.ignore_patterns(".git", "__pycache__", "*.pyc"))
    changed_files = 0
    replacements = 0
    for path in output.rglob("*"):
        if not path.is_file() or (path.suffix.lower() not in TEXT_SUFFIXES and path.name != "Snakefile"):
            continue
        try:
            text = path.read_text(encoding="utf-8-sig")
        except UnicodeError:
            continue
        updated = text
        for old, new in mapping:
            count = updated.count(old)
            replacements += count
            updated = updated.replace(old, new)
        if updated != text:
            path.write_text(updated, encoding="utf-8", newline="\n")
            changed_files += 1
    print(f"Configured copy: {output}")
    print(f"Changed files: {changed_files}")
    print(f"Path replacements: {replacements}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
