#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# export_footprints_onefolder_i0.py

import argparse
import os
import re
from typing import Dict, Optional, Tuple, List

def _prepare_hdf5_runtime() -> None:
    os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if conda_prefix:
        plugin_dir = os.path.join(conda_prefix, "lib", "hdf5", "plugin")
        try:
            os.makedirs(plugin_dir, exist_ok=True)
            os.environ.setdefault("HDF5_PLUGIN_PATH", plugin_dir)
        except Exception:
            pass

_prepare_hdf5_runtime()

try:
    import hdf5plugin  # noqa: F401
except Exception:
    pass

import h5py  # noqa: E402
import numpy as np  # noqa: E402

import matplotlib  # noqa: E402
matplotlib.use("Agg")

# =========================
# =========================
matplotlib.rcParams["pdf.fonttype"] = 42   # TrueType
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "DejaVu Sans"
matplotlib.rcParams["axes.unicode_minus"] = False

import matplotlib.pyplot as plt  # noqa: E402


def safe_token(s: str) -> str:
    s = s.strip()
    s = s.replace("/", "_")
    return re.sub(r"[^\w\.\-]+", "_", s)


def parse_ylim(s: str) -> Optional[Tuple[float, float]]:
    s = (s or "").strip()
    if not s:
        return None
    a, b = s.split(",")
    return float(a), float(b)


def load_i0_curves(h5_path: str) -> Dict[str, np.ndarray]:
    curves: Dict[str, np.ndarray] = {}
    try:
        fctx = h5py.File(h5_path, "r", locking="best-effort")
    except TypeError:
        fctx = h5py.File(h5_path, "r")

    with fctx as f:
        for motif in f.keys():
            g = f.get(motif, None)
            if g is None or not isinstance(g, h5py.Group):
                continue
            if "i0" not in g:
                continue
            arr = np.asarray(g["i0"])
            if arr.ndim != 1:
                continue
            curves[motif] = arr

    if not curves:
        raise RuntimeError(f"No <motif>/i0 1D arrays found in {h5_path}")
    return curves


def center_slice(y: np.ndarray, plot_bp: int) -> Tuple[np.ndarray, np.ndarray]:
    n = int(len(y))
    if plot_bp <= 0 or plot_bp >= n:
        half = n // 2
        if n % 2 == 0:
            x = np.arange(-half, half, dtype=int)
        else:
            x = np.arange(-half, half + 1, dtype=int)
        return x, y

    half = plot_bp // 2
    center = n // 2
    start = max(center - half, 0)
    end = min(start + plot_bp, n)

    y2 = y[start:end]
    x = np.arange(-half, -half + len(y2), dtype=int)
    return x, y2


def plot_one(tissue: str,
             motif: str,
             y: np.ndarray,
             out_pdf: str,
             ylim: Optional[Tuple[float, float]],
             plot_bp: int,
             dpi: int) -> None:
    x, y2 = center_slice(y, plot_bp)

    fig = plt.figure(figsize=(6.2, 3.6))
    ax = fig.add_subplot(111)

    ax.plot(x, y2, linewidth=1.5)
    ax.axvline(0, linewidth=1.0)

    ax.set_ylabel("Probability")
    if plot_bp <= 0:
        ax.set_xlabel("Position relative to motif center (bp)")
    else:
        ax.set_xlabel(f"{plot_bp}bp around motif center (bp)")

    ax.set_title(f"{tissue} | {motif}")

    if ylim is not None:
        ax.set_ylim(*ylim)

    fig.tight_layout()
    fig.savefig(out_pdf, dpi=dpi)   #   TrueType   PDF
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_h5", required=True, help="input *_footprints.h5")
    ap.add_argument("--tissue", required=True, help="tissue name used in output file name")
    ap.add_argument("--outdir", required=True, help="all PDFs go into this folder (no subdirs)")
    ap.add_argument("--plot_bp", type=int, default=200,
                    help="center window bp (default 200). Use 0 for full length.")
    ap.add_argument("--ylim", default="", help="optional y-range like '0,0.003'")
    ap.add_argument("--motifs", default="", help="optional comma-separated motif names to export")
    ap.add_argument("--dpi", type=int, default=200, help="save dpi (pdf is vector; keep reasonable)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    ylim = parse_ylim(args.ylim)

    keep: Optional[List[str]] = None
    if args.motifs.strip():
        keep = [m.strip() for m in args.motifs.split(",") if m.strip()]

    curves = load_i0_curves(args.in_h5)

    for motif, arr in curves.items():
        if keep is not None and motif not in keep:
            continue
        motif_token = safe_token(motif)
        out_pdf = os.path.join(args.outdir, f"{args.tissue}_{motif_token}_footprint.pdf")
        plot_one(args.tissue, motif, arr, out_pdf, ylim, args.plot_bp, args.dpi)


if __name__ == "__main__":
    main()
