"""Static path and executable preflight checks; no analysis is run."""

from __future__ import annotations

import argparse
import os
import shutil
from pathlib import Path


SCOPE_KEYS = {
    "cre": ["SRC_WORKDIR", "E1E9_WORKDIR", "LIFTOVER_BIN", "CHAIN_HG38_HG19"],
    "gwas": ["GWAS_BASE", "OLD_GWAS_BASE", "LDSC_ROOT", "LDSC_PYTHON", "LDSC_REF"],
    "fig8": [
        "E1E9_WORKDIR",
        "GWAS_BASE",
        "CONSERVATION_TABLE",
        "HUMAN_EMISSION",
        "SHEEP_EMISSION_R1",
        "SHEEP_EMISSION_R2",
        "TRAIT_38_XLSX",
        "TRAIT_152_TSV",
        "TISSUE_COLOR_FILE",
        "GWAS_STATS",
    ],
    "figure_inputs": [
        "E1E9_WORKDIR",
        "CONSERVATION_TABLE",
        "PHASTCONS_BW",
        "PHYLOP_BW",
        "GERP_BW",
    ],
}

EXPECTED_OUTPUT_KEYS = {
    "E1E9_WORKDIR",
    "GWAS_BASE",
    "CONSERVATION_TABLE",
    "GWAS_STATS",
    "TRAIT_152_TSV",
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=[*SCOPE_KEYS, "all"], default="all")
    args = parser.parse_args()

    scopes = list(SCOPE_KEYS) if args.scope == "all" else [args.scope]
    failed = False
    for scope in scopes:
        print(f"[{scope}]")
        for key in SCOPE_KEYS[scope]:
            value = os.environ.get(key)
            if not value:
                print(f"  MISSING_CONFIG\t{key}")
                failed = True
                continue
            path = Path(value)
            status = "exists" if path.exists() else "not-found"
            print(f"  {status}\t{key}\t{path}")
            if key not in EXPECTED_OUTPUT_KEYS and not path.exists():
                failed = True

    if args.scope in {"cre", "all"}:
        for executable in ("bash", "awk", "sort", "bedtools"):
            found = shutil.which(executable)
            print(f"  {'found' if found else 'not-found'}\texecutable\t{executable}")
            failed = failed or found is None

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
