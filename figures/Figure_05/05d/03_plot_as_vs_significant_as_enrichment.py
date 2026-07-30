#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from scipy.stats import wilcoxon
except Exception as e:
    raise ImportError("  scipy   Wilcoxon  , :pip install scipy") from e


# -----------------------------
# -----------------------------
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]


def read_input_table(input_file: Path) -> pd.DataFrame:
    """
     read;  1  , read.
     .
    """
    df = pd.read_csv(input_file, sep="\t", dtype=str)

    if df.shape[1] == 1:
        df = pd.read_csv(input_file, sep=r"\s+", engine="python", dtype=str)

    df.columns = [str(c).strip() for c in df.columns]
    return df


def read_tissue_colors(color_file: Path) -> dict:
    """
    read file.
     :tissue, color
    """
    df = pd.read_csv(color_file, sep="\t", dtype=str)

    if df.shape[1] == 1:
        df = pd.read_csv(color_file, sep=r"\s+", engine="python", dtype=str)

    df.columns = [str(c).strip() for c in df.columns]

    if "tissue" not in df.columns or "color" not in df.columns:
        raise ValueError(
            f"{color_file}  : tissue, color; : {', '.join(df.columns)}"
        )

    df["tissue"] = df["tissue"].astype(str).str.strip()
    df["color"] = df["color"].astype(str).str.strip()
    return dict(zip(df["tissue"], df["color"]))


def to_numeric_safe(series):
    return pd.to_numeric(series, errors="coerce")


def pvalue_to_stars(p):
    """
      p  .
    """
    if pd.isna(p):
        return "NA"
    if p < 1e-4:
        return "****"
    if p < 1e-3:
        return "***"
    if p < 1e-2:
        return "**"
    if p < 5e-2:
        return "*"
    return "ns"


def paired_wilcoxon_test(a, b):
    """
      Wilcoxon  .
      NA  .
      p  .
    """
    x = to_numeric_safe(pd.Series(a))
    y = to_numeric_safe(pd.Series(b))

    mask = x.notna() & y.notna()
    x = x[mask].to_numpy(dtype=float)
    y = y[mask].to_numpy(dtype=float)

    if len(x) < 2:
        return np.nan, len(x)

    diff = x - y
    if np.allclose(diff, 0):
        return 1.0, len(x)

    try:
        p = wilcoxon(x, y, alternative="two-sided", zero_method="pratt").pvalue
    except ValueError:
        p = 1.0

    return p, len(x)


def add_sig_bracket(ax, x1, x2, y, h, text, fontsize=12):
    """
     , .
    """
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=1.3, c="black")
    ax.text((x1 + x2) / 2, y + h, text, ha="center", va="bottom", fontsize=fontsize)


def ensure_required_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """
     inputfile .
     file :
      - p2_over_q1_t
      - p3_over_q1_t
     :
      - p1_over_q1_t

     :
      p1_over_q1_t = p1_t_testHet_over_allVCF / q1_t_motifAmongTestHet
    """
    df = df.copy()
    df.columns = [str(c).strip() for c in df.columns]

    required_base = ["tissue", "p2_over_q1_t", "p3_over_q1_t"]
    missing_base = [c for c in required_base if c not in df.columns]
    if missing_base:
        raise ValueError(
            "inputfileMissing required columns: "
            + ", ".join(missing_base)
            + f"\n :\n{', '.join(df.columns)}"
        )

    if "p1_over_q1_t" not in df.columns:
        need_for_compute = ["p1_t_testHet_over_allVCF", "q1_t_motifAmongTestHet"]
        missing_for_compute = [c for c in need_for_compute if c not in df.columns]
        if missing_for_compute:
            raise ValueError(
                "inputfile  p1_over_q1_t, , : "
                + ", ".join(missing_for_compute)
                + f"\n :\n{', '.join(df.columns)}"
            )

        p1 = to_numeric_safe(df["p1_t_testHet_over_allVCF"])
        q1 = to_numeric_safe(df["q1_t_motifAmongTestHet"])

        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = p1 / q1
        ratio = ratio.replace([np.inf, -np.inf], np.nan)

        df["p1_over_q1_t"] = ratio
        print("[INFO] inputfile  p1_over_q1_t, :")
        print("       p1_over_q1_t = p1_t_testHet_over_allVCF / q1_t_motifAmongTestHet")

    return df


def plot_paired_boxplot(
    df: pd.DataFrame,
    metrics,
    xlabels,
    color_map,
    title,
    ylabel,
    out_png: Path,
    out_pdf: Path,
    point_size=46,
    jitter_width=0.028,
):
    """
     :
    -
    -   tissue color map
    -
    -   Wilcoxon
    """

    need_cols = ["tissue"] + list(metrics)
    miss = [c for c in need_cols if c not in df.columns]
    if miss:
        raise ValueError(f"inputfile : {', '.join(miss)}")

    sub = df[need_cols].copy()
    sub["tissue"] = sub["tissue"].astype(str).str.strip()

    for m in metrics:
        sub[m] = to_numeric_safe(sub[m])

    sub = sub.dropna(subset=metrics).copy()

    if sub.empty:
        raise ValueError(f"{title}  , check : {metrics}")

    fig, ax = plt.subplots(figsize=(5.3, 6.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    positions = [1.00, 1.38]
    box_data = [sub[m].to_numpy(dtype=float) for m in metrics]

    ax.boxplot(
        box_data,
        positions=positions,
        widths=0.18,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.4),
        boxprops=dict(facecolor="white", edgecolor="black", linewidth=1.2),
        whiskerprops=dict(color="black", linewidth=1.0),
        capprops=dict(color="black", linewidth=1.0),
    )

    rng = np.random.default_rng(20260415)

    for _, row in sub.iterrows():
        tissue = row["tissue"]
        point_color = color_map.get(tissue, "#BDBDBD")

        for i, m in enumerate(metrics):
            xj = positions[i] + rng.uniform(-jitter_width, jitter_width)
            yv = float(row[m])

            ax.scatter(
                xj,
                yv,
                s=point_size,
                color=point_color,
                edgecolor="black",
                linewidth=0.6,
                zorder=3,
                clip_on=False,
            )

    ax.set_xticks(positions)
    ax.set_xticklabels(xlabels, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_title(title, fontsize=13)
    ax.grid(False)

    ax.set_xlim(0.84, 1.54)

    p, n = paired_wilcoxon_test(sub[metrics[0]], sub[metrics[1]])
    stars = pvalue_to_stars(p)

    all_vals = np.concatenate(box_data)
    finite_vals = all_vals[np.isfinite(all_vals)]

    y_data_min = finite_vals.min()
    y_data_max = finite_vals.max()
    yr = max(y_data_max - y_data_min, 1e-6)

    if y_data_min >= 0:
        y_lower = min(-0.45, y_data_min - 0.18 * yr)
    else:
        y_lower = y_data_min - 0.12 * yr

    bracket_y = y_data_max + 0.08 * yr
    bracket_h = max(0.03 * yr, 0.12)
    sig_text = f"{stars}\nP = {p:.2e}" if pd.notna(p) else "NA"

    add_sig_bracket(
        ax=ax,
        x1=positions[0],
        x2=positions[1],
        y=bracket_y,
        h=bracket_h,
        text=sig_text,
        fontsize=11,
    )

    y_upper = bracket_y + bracket_h + 0.10 * yr
    ax.set_ylim(y_lower, y_upper)

    for spine in ax.spines.values():
        spine.set_linewidth(1.0)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    print(f"[{title}] paired Wilcoxon p = {p:.6e}, n = {n}")
    print(f"[OUT] {out_png}")
    print(f"[OUT] {out_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description=" : , file, ."
    )
    parser.add_argument(
        "--input",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/p1t_q1t_p2_p3_ratios_strict_peak/p1t_q1t_p2_p3_ratios_strictPeak.tsv",
        help="inputresults "
    )
    parser.add_argument(
        "--color-file",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/tissue_colors.tsv",
        help=" ,  tissue,color  "
    )
    parser.add_argument(
        "--outdir",
        default="/data/home/sczd644/run/zsw_chrombpnet/snpscore/non_as_snpscore/non_AS_matched/p1t_q1t_p2_p3_ratios_strict_peak/plots_boxplot",
        help="Output directory"
    )
    args = parser.parse_args()

    input_file = Path(args.input)
    color_file = Path(args.color_file)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not input_file.exists():
        raise FileNotFoundError(f"Not foundinputfile: {input_file}")
    if not color_file.exists():
        raise FileNotFoundError(f"Not found file: {color_file}")

    df = read_input_table(input_file)
    print("[INFO] inputfile :")
    print(", ".join(df.columns))

    df = ensure_required_metrics(df)
    color_map = read_tissue_colors(color_file)

    plot_paired_boxplot(
        df=df,
        metrics=["p2_over_q1_t", "p3_over_q1_t"],
        xlabels=["p2_over_q1_t", "p3_over_q1_t"],
        color_map=color_map,
        title="Fold enrichment comparison: p2_over_q1_t vs p3_over_q1_t",
        ylabel="Fold enrichment",
        out_png=outdir / "boxplot_p2_over_q1_t_vs_p3_over_q1_t.png",
        out_pdf=outdir / "boxplot_p2_over_q1_t_vs_p3_over_q1_t.pdf",
    )

    plot_paired_boxplot(
        df=df,
        metrics=["p1_over_q1_t", "p2_over_q1_t"],
        xlabels=["p1_over_q1_t", "p2_over_q1_t"],
        color_map=color_map,
        title="Fold enrichment comparison: p1_over_q1_t vs p2_over_q1_t",
        ylabel="Fold enrichment",
        out_png=outdir / "boxplot_p1_over_q1_t_vs_p2_over_q1_t.png",
        out_pdf=outdir / "boxplot_p1_over_q1_t_vs_p2_over_q1_t.pdf",
    )

    print("[DONE]  output ")


if __name__ == "__main__":
    main()
