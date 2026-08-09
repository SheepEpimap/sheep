#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
match_homer_vs_deeplearning.py

 (Tissue) + motif   all_tables(HOMER/FIMO vs DeepLearning/Fi-NeMo).
 :
  -   (tissue, motif_base)   key
  - motif_base = motif.split("#", 1)[0]  #   deep learning   #1/#2
output:
  1) outdir/ALL.matched_homer_vs_deeplearning.tsv            #
  2) outdir/by_tissue/<tissue>.matched_homer_vs_deeplearning.tsv  #  file
  3) outdir/ALL.matched_motifs.summary.tsv                   # statistics  motif_base
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import defaultdict
from typing import Iterator, List, Tuple, Optional


def motif_base(m: str) -> str:
    """KLF5#1 -> KLF5;  #  """
    return m.split("#", 1)[0]


def split_fields(line: str) -> List[str]:
    """
     input:
      -   TSV(  \\t)
      -  (  \\t)
    """
    line = line.rstrip("\n").rstrip("\r")
    if not line:
        return []
    if "\t" in line:
        return line.split("\t")
    return line.split()


def read_header(path: str, encoding: str, errors: str) -> List[str]:
    """read """
    with open(path, "r", encoding=encoding, errors=errors) as f:
        for line in f:
            if line.strip():
                return split_fields(line)
    raise RuntimeError(f"Empty file: {path}")


def iter_rows(path: str, encoding: str, errors: str) -> Iterator[List[str]]:
    """ ,  yield  """
    with open(path, "r", encoding=encoding, errors=errors) as f:
        header_skipped = False
        for line in f:
            if not header_skipped:
                if line.strip():
                    header_skipped = True
                continue
            if not line.strip():
                continue
            fields = split_fields(line)
            if fields:
                yield fields


def pad(fields: List[str], n: int) -> List[str]:
    """ , output """
    if len(fields) >= n:
        return fields
    return fields + [""] * (n - len(fields))


def collect_keys(path: str, encoding: str, errors: str) -> set[Tuple[str, str]]:
    """  (tissue, motif_base)  """
    keys = set()
    for fields in iter_rows(path, encoding, errors):
        if len(fields) < 2:
            continue
        tissue = fields[0].strip()
        motif = fields[1].strip()
        if not tissue or not motif:
            continue
        keys.add((tissue, motif_base(motif)))
    return keys


def main():
    ap = argparse.ArgumentParser(
        description="Match HOMER vs DeepLearning all_tables by Tissue + Motif (ignore # variants)."
    )
    ap.add_argument("--homer", required=True, help="HOMER/FIMO all_tables.with_header.tsv")
    ap.add_argument("--dl", required=True, help="DeepLearning/Fi-NeMo ALL.all_tables.with_header.fullChain.txt")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--encoding", default="utf-8", help="File encoding (default: utf-8)")
    ap.add_argument("--errors", default="replace", help="Encoding errors handling (default: replace)")
    ap.add_argument("--no_per_tissue", action="store_true", help="Do not output per-tissue files")

    args = ap.parse_args()

    homer_path = args.homer
    dl_path = args.dl
    outdir = args.outdir
    encoding = args.encoding
    errors = args.errors
    per_tissue = not args.no_per_tissue

    for p in (homer_path, dl_path):
        if not os.path.isfile(p) or os.path.getsize(p) == 0:
            raise FileNotFoundError(f"Missing/empty input: {p}")

    os.makedirs(outdir, exist_ok=True)
    by_tissue_dir = os.path.join(outdir, "by_tissue")
    if per_tissue:
        os.makedirs(by_tissue_dir, exist_ok=True)

    homer_header = read_header(homer_path, encoding, errors)
    dl_header = read_header(dl_path, encoding, errors)
    max_cols = max(len(homer_header), len(dl_header))

    out_header = pad(homer_header, max_cols)
    for i in range(len(homer_header), max_cols):
        out_header[i] = f"Extra_{i+1}"
    out_header.append("Method")

    print("[INFO] Collecting keys...")
    keys_homer = collect_keys(homer_path, encoding, errors)
    keys_dl = collect_keys(dl_path, encoding, errors)
    keys_inter = keys_homer & keys_dl

    if not keys_inter:
        raise RuntimeError("No matched (Tissue, MotifBase) between HOMER and DeepLearning tables.")

    print(f"[INFO] Matched keys: {len(keys_inter)} (tissue,motif_base)")

    out_all = os.path.join(outdir, "ALL.matched_homer_vs_deeplearning.tsv")
    out_sum = os.path.join(outdir, "ALL.matched_motifs.summary.tsv")

    tissue_fhs: dict[str, csv.writer] = {}
    tissue_files: dict[str, object] = {}

    def get_writer_for_tissue(tissue: str) -> Optional[csv.writer]:
        if not per_tissue:
            return None
        if tissue in tissue_fhs:
            return tissue_fhs[tissue]
        path = os.path.join(by_tissue_dir, f"{tissue}.matched_homer_vs_deeplearning.tsv")
        fh = open(path, "w", encoding="utf-8", newline="")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_header)
        tissue_files[tissue] = fh
        tissue_fhs[tissue] = w
        return w

    counts = defaultdict(lambda: {"homer": 0, "deeplearning": 0})

    with open(out_all, "w", encoding="utf-8", newline="") as fh_all:
        w_all = csv.writer(fh_all, delimiter="\t", lineterminator="\n")
        w_all.writerow(out_header)

        def write_file(path: str, method: str):
            for fields in iter_rows(path, encoding, errors):
                if len(fields) < 2:
                    continue
                tissue = fields[0].strip()
                motif = fields[1].strip()
                if not tissue or not motif:
                    continue
                key = (tissue, motif_base(motif))
                if key not in keys_inter:
                    continue

                row = pad(fields, max_cols) + [method]
                w_all.writerow(row)

                wt = get_writer_for_tissue(tissue)
                if wt is not None:
                    wt.writerow(row)

                counts[key][method] += 1

        print("[INFO] Writing matched rows (HOMER)...")
        write_file(homer_path, "homer")
        print("[INFO] Writing matched rows (DeepLearning)...")
        write_file(dl_path, "deeplearning")

    for fh in tissue_files.values():
        fh.close()

    # 5) output summary
    with open(out_sum, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["Tissue", "MotifBase", "Homer_rows", "DeepLearning_rows"])
        for (tissue, mb) in sorted(counts.keys()):
            w.writerow([tissue, mb, counts[(tissue, mb)]["homer"], counts[(tissue, mb)]["deeplearning"]])

    print("[DONE]")
    print("  All matched:", out_all)
    print("  Summary    :", out_sum)
    if per_tissue:
        print("  By tissue  :", by_tissue_dir)


if __name__ == "__main__":
    main()
