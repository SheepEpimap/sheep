#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "DejaVu Sans",
    "text.usetex": False,
})

import matplotlib.pyplot as plt


def safe_name(motif: str) -> str:
    out = []
    for ch in motif:
        if ch.isalnum() or ch in "._-":
            out.append(ch)
        else:
            out.append("_")
    return "".join(out)


def read_meta(meta_tsv: Path) -> dict:
    if not meta_tsv.exists():
        return {}
    meta = {}
    with meta_tsv.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if "key" in row and "value" in row:
                meta[row["key"]] = row["value"]
    return meta


def read_gerp(gerp_tsv: Path, xcol="off", ycol="mean_gerp"):
    xs, ys = [], []
    with gerp_tsv.open("r", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if ycol not in row or xcol not in row:
                raise ValueError(f"Missing column {xcol}/{ycol} in {gerp_tsv}")
            yraw = row[ycol]
            if yraw is None or yraw == "" or str(yraw).upper() == "NA":
                continue
            try:
                x = int(float(row[xcol]))
                y = float(yraw)
            except Exception:
                continue
            if math.isnan(y):
                continue
            xs.append(x)
            ys.append(y)
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    xs = [xs[i] for i in order]
    ys = [ys[i] for i in order]
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results_dir", required=True,
                    help="OUTDIR/results (contains <SAFE>/<SAFE>.gerp_by_offset.w*.tsv)")
    ap.add_argument("--motif", required=True, help='motif name, e.g. "CTCF#2"')
    ap.add_argument("--window", type=int, default=25, help="which window file to load: w{window}")
    ap.add_argument("--box_width", type=int, required=True, help="blue box width in bp (for later logo alignment)")
    ap.add_argument("--box_center", type=float, default=0.0, help="box center on x-axis (default 0)")
    ap.add_argument("--out", required=True, help="output prefix (no extension)")
    ap.add_argument("--dpi", type=int, default=300)
    args = ap.parse_args()

    safe = safe_name(args.motif)
    mdir = Path(args.results_dir) / safe

    gerp_tsv = mdir / f"{safe}.gerp_by_offset.w{args.window}.tsv"
    meta_tsv = mdir / f"{safe}.meta.tsv"

    if not gerp_tsv.exists():
        raise SystemExit(f"[ERROR] not found: {gerp_tsv}")

    meta = read_meta(meta_tsv)
    xs, ys = read_gerp(gerp_tsv)

    left = args.box_center - (args.box_width // 2)
    right = left + args.box_width

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.plot(xs, ys, color="red", linewidth=2.0)
    ax.axvline(0, color="black", linewidth=1.0, alpha=0.8)
    ax.axvspan(left, right, ymin=0.0, ymax=1.0,
               facecolor="dodgerblue", edgecolor="dodgerblue",
               alpha=0.18, linewidth=1.0)

    ax.set_xlim(-args.window, args.window)
    ax.set_xlabel("Distance from motif center (bp)")
    ax.set_ylabel("Mean GERP (covered bases only)")

    title = args.motif
    if "n_hits" in meta:
        title += f" (N_hits={meta['n_hits']})"
    ax.set_title(title)

    out_prefix = Path(args.out)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(str(out_prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".png", dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    with (out_prefix.parent / (out_prefix.name + ".box.tsv")).open("w", encoding="utf-8") as f:
        f.write("motif\tsafe\tbox_left\tbox_right\tbox_width\twindow\n")
        f.write(f"{args.motif}\t{safe}\t{left}\t{right}\t{args.box_width}\t{args.window}\n")

    print("[DONE]", str(out_prefix) + ".png")
    print("[DONE]", str(out_prefix) + ".pdf")


if __name__ == "__main__":
    main()
