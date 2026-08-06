#!/usr/bin/env python3
"""Unified launcher for the SheepEpimap E1-E9 workflow."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent
SRC = REPO / "src"
sys.path.insert(0, str(SRC))

from sheep_epimap.config import load_paths  # noqa: E402

CRE = SRC / "cre_pipeline"
GWAS = SRC / "gwas_ldsc"
FIG8 = SRC / "sheep_epimap" / "figures" / "fig8.py"
FIGURE_INPUTS = SRC / "sheep_epimap" / "figure_inputs"

STAGES: dict[str, list[str]] = {
    "preflight": [sys.executable, str(SRC / "sheep_epimap" / "preflight.py"), "--scope", "all"],
    "cre_01_reclassify": ["bash", str(CRE / "01_reclassify_e1_e9.sh")],
    "cre_02_hg38_to_hg19": ["bash", str(CRE / "02_hg38_to_hg19.sh")],
    "cre_03_collect_e5": ["bash", str(CRE / "03_collect_e5_hg19.sh")],
    "cre_04_reuse_human_tsr": ["bash", str(CRE / "04_reuse_human_tsr.sh")],
    "cre_05_label_tsr": ["bash", str(CRE / "05_label_cre_with_tsr.sh")],
    "cre_06_reuse_projection": ["bash", str(CRE / "06_reuse_sheep_projection.sh")],
    "cre_07_split_ldsc": ["bash", str(CRE / "07_split_ldsc_beds.sh")],
    "cre_08_collect_gwas_beds": ["bash", str(CRE / "08_collect_gwas_beds.sh")],
    "cre_submit_slurm": ["bash", str(CRE / "submit_stepwise_slurm.sh")],
    "traits_152": [sys.executable, str(GWAS / "expand_to_152.py")],
    "gwas_prepare": ["bash", str(GWAS / "prepare_E1E9_reuse_existing_recalc.sh")],
    "gwas_submit_slurm": ["bash", str(GWAS / "submit_E1E9_reuse_existing_152.sh")],
    "gwas_summary": ["bash", str(GWAS / "run_summaries.sh")],
    "figure_input_merge": [sys.executable, str(FIGURE_INPUTS / "build_merged_classification.py")],
    "figure_input_conservation": [sys.executable, str(FIGURE_INPUTS / "score_conservation.py")],
    "fig8": [sys.executable, str(FIG8)],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paths", type=Path, required=True, help="Tab-separated path configuration.")
    parser.add_argument("--stage", choices=sorted(STAGES), required=True)
    parser.add_argument("--run-id", help="Immutable run identifier used with RUNS_ROOT.")
    parser.add_argument("--panels", nargs="+", choices=list("bdefg"), default=list("bdefg"))
    parser.add_argument("--dry-run", action="store_true", help="Print the resolved command without running it.")
    return parser.parse_args()


def resolve_environment(paths: dict[str, str], run_id: str | None) -> tuple[dict[str, str], Path | None]:
    env = os.environ.copy()
    env.update(paths)
    run_root: Path | None = None
    if run_id:
        if not paths.get("RUNS_ROOT"):
            raise ValueError("RUNS_ROOT is required when --run-id is used")
        run_root = Path(paths["RUNS_ROOT"]) / run_id
        # A run ID always selects new output roots; configured existing result
        # roots are used only when --run-id is omitted for read/plot workflows.
        env["E1E9_WORKDIR"] = str(run_root / "cre_e1e9")
        env["GWAS_BASE"] = str(run_root / "gwas_ldsc")
        env["FIG8_FINAL_DIR"] = str(run_root / "figures" / "fig8")
        env["MERGED_CLASSIFICATION"] = str(
            run_root / "cre_e1e9" / "visualization_tables_E1E9_reuse_existing" / "results_hg38" / "merged_classification_summary.tsv"
        )
        env["CONSERVATION_TABLE"] = str(run_root / "figure_inputs" / "cre_conservation_scores_hg38_E1E9.tsv.gz")
        env["GWAS_STATS"] = str(
            run_root / "gwas_ldsc" / "summary" / "A2_B2_all_tissues_all_classes_enrichment_stats_long.csv"
        )
        env["TRAIT_152_TSV"] = str(run_root / "gwas_ldsc" / "all_152_traits.tsv")
    if env.get("E1E9_WORKDIR"):
        env["WORKDIR"] = env["E1E9_WORKDIR"]
        env["SRC"] = env["E1E9_WORKDIR"]
    return env, run_root


def write_run_snapshot(run_root: Path, args: argparse.Namespace, env: dict[str, str], command: list[str]) -> None:
    run_root.mkdir(parents=True, exist_ok=True)
    snapshots = run_root / "config_snapshot"
    snapshots.mkdir(exist_ok=True)
    snapshot = snapshots / args.paths.name
    if snapshot.exists():
        if snapshot.read_bytes() != args.paths.read_bytes():
            raise RuntimeError(f"Run configuration changed after creation: {snapshot}")
    else:
        shutil.copy2(args.paths, snapshot)
    configured_keys = set(load_paths(args.paths))
    metadata = {
        "stage": args.stage,
        "run_id": args.run_id,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "host": socket.gethostname(),
        "command": command,
        "configured_keys": sorted(k for k in env if k in configured_keys),
    }
    (run_root / f"command.{args.stage}.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def main() -> int:
    args = parse_args()
    paths = load_paths(args.paths)
    env, run_root = resolve_environment(paths, args.run_id)
    command = list(STAGES[args.stage])
    if args.stage == "fig8":
        command.extend(["--panels", *args.panels])

    print("stage:", args.stage)
    print("command:", shlex.join(command))
    for key in ("SRC_WORKDIR", "E1E9_WORKDIR", "GWAS_BASE", "FIG8_FINAL_DIR"):
        if env.get(key):
            print(f"{key}={env[key]}")

    if args.dry_run:
        return 0

    if args.stage not in {"preflight", "fig8"} and run_root is None:
        raise SystemExit("Execution of data-generating stages requires --run-id")

    marker: Path | None = None
    if run_root is not None:
        marker = run_root / f".completed.{args.stage}"
        if marker.exists():
            raise SystemExit(f"Stage already completed; refusing overwrite: {marker}")
        write_run_snapshot(run_root, args, env, command)

    subprocess.run(command, cwd=REPO, env=env, check=True)
    if marker is not None:
        marker.write_text(datetime.now(timezone.utc).isoformat() + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
