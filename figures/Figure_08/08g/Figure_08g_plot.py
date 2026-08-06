#!/usr/bin/env python3
"""Generate Figure 8g: sfCRE trait-by-tissue GWAS enrichment dot plot."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paths", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    workflow = (
        Path(__file__).resolve().parents[1]
        / "Shared_Workflows"
        / "Figure_08_CRE_GWAS"
        / "workflow.py"
    )
    command = [
        sys.executable,
        str(workflow),
        "--paths",
        str(args.paths),
        "--stage",
        "fig8",
        "--panels",
        "g",
    ]
    if args.dry_run:
        command.append("--dry-run")
    subprocess.run(command, check=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
