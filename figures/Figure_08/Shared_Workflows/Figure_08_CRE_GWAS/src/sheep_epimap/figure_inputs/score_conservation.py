#!/usr/bin/env python3
"""Score merged E1-E9 CREs with phastCons, phyloP and GERP bigWigs."""

from __future__ import annotations

import argparse
import os
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

STATES = [f"E{i}" for i in range(1, 10)]
CRE_TYPES = ["sfCRE", "sdCRE", "soCRE", "ssCRE"]


def env_path(name: str) -> Path | None:
    return Path(os.environ[name]) if os.environ.get(name) else None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    workdir = env_path("E1E9_WORKDIR")
    default_input = workdir / "visualization_tables_E1E9_reuse_existing" / "results_hg38" if workdir else None
    parser.add_argument("--input-dir", type=Path, default=default_input)
    parser.add_argument("--output", type=Path, default=env_path("CONSERVATION_TABLE"))
    parser.add_argument("--phastcons", type=Path, default=env_path("PHASTCONS_BW"))
    parser.add_argument("--phylop", type=Path, default=env_path("PHYLOP_BW"))
    parser.add_argument("--gerp", type=Path, default=env_path("GERP_BW"))
    parser.add_argument("--executable", default=os.environ.get("BIGWIG_AVERAGE_OVER_BED", "bigWigAverageOverBed"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    for key in ("input_dir", "output", "phastcons", "phylop", "gerp"):
        if getattr(args, key) is None:
            parser.error(f"--{key.replace('_', '-')} or its configured environment variable is required")
    return args


def load_intervals(input_dir: Path) -> pd.DataFrame:
    records: list[dict[str, int | str]] = []
    for state in STATES:
        for cre_type in CRE_TYPES:
            path = input_dir / f"{state}.{cre_type}.bed"
            if not path.is_file() or path.stat().st_size == 0:
                continue
            with path.open(encoding="utf-8") as handle:
                for line in handle:
                    columns = line.rstrip("\n").split("\t")
                    if len(columns) >= 3:
                        records.append(
                            {
                                "chrom": columns[0],
                                "start": int(columns[1]),
                                "end": int(columns[2]),
                                "state": columns[3] if len(columns) > 3 else state,
                                "cre_type": cre_type,
                            }
                        )
    if not records:
        raise ValueError(f"No merged BED intervals found under {input_dir}")
    return pd.DataFrame(records).drop_duplicates().reset_index(drop=True)


def write_named_bed(frame: pd.DataFrame, path: Path, strip_chr: bool) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        for index, row in frame.iterrows():
            chrom = str(row["chrom"])
            if strip_chr and chrom.startswith("chr"):
                chrom = chrom[3:]
            length = int(row["end"] - row["start"])
            handle.write(f"{chrom}\t{row['start']}\t{row['end']}\t{index}\t{length}\n")


def score(frame: pd.DataFrame, executable: str, bigwig: Path, bed: Path, name: str, output: Path) -> None:
    subprocess.run([executable, str(bigwig), str(bed), str(output)], check=True)
    values: dict[int, tuple[int, float]] = {}
    with output.open(encoding="utf-8") as handle:
        for line in handle:
            columns = line.rstrip("\n").split("\t")
            values[int(columns[0])] = (int(columns[2]), float(columns[5]))
    frame[f"{name}_covered_bp"] = frame.index.map(lambda index: values.get(index, (0, np.nan))[0])
    frame[f"{name}_mean"] = frame.index.map(lambda index: values.get(index, (0, np.nan))[1])


def main() -> int:
    args = parse_args()
    if args.output.exists() and not args.force:
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output}")
    for path in (args.phastcons, args.phylop, args.gerp):
        if not path.is_file():
            raise FileNotFoundError(path)

    frame = load_intervals(args.input_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="sheep_epimap_conservation_") as directory:
        temp = Path(directory)
        bed_chr = temp / "cre.hg38.bed"
        bed_no_chr = temp / "cre.no_chr.bed"
        write_named_bed(frame, bed_chr, strip_chr=False)
        write_named_bed(frame, bed_no_chr, strip_chr=True)
        score(frame, args.executable, args.phastcons, bed_chr, "phastCons", temp / "phastCons.tab")
        score(frame, args.executable, args.phylop, bed_chr, "phyloP", temp / "phyloP.tab")
        score(frame, args.executable, args.gerp, bed_no_chr, "GERP", temp / "GERP.tab")
    frame.to_csv(args.output, sep="\t", index=False, compression="gzip")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
