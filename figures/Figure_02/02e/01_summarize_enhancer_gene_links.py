#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-
"""
plot_topN_links_to_pdfs.py

 :
1)   Expression_tpm_39.tsv + Expression_tpm_40.tsv   86  (TPM),  expected_order  ;
2) read H3K27ac_*_counts_cpms.csv(  bam file ),canon   tissue_39/40,  expected_order  ;
3)   confident links file  topN(  LOC),  --pairs file  enhancer+gene;
4) output:
   - enhancer degree  ( ) PDF
   - gene degree  ( ) PDF
   -   pair   PDF( ,rep  ; ,r,p;  topX/topY  )

  Python 3.7/3.8(  “str | None”  ).
"""

import os
import re
import csv
import math
import argparse
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.stats import pearsonr


# -----------------------------
# -----------------------------
EXPECTED_ORDER = [
    "abomasum_39", "abomasum_40",
    "adipose_39", "adipose_40",
    "bone-marrow_39", "bone-marrow_40",
    "brainstem_39", "brainstem_40",
    "cecum_39", "cecum_40",
    "cerebellum_39", "cerebellum_40",
    "cerebral-cortex_39", "cerebral-cortex_40",
    "cervix_39", "cervix_40",
    "colon_39", "colon_40",
    "cornua-uteri_39", "cornua-uteri_40",
    "corpus-uteri_39", "corpus-uteri_40",
    "duodenum_39", "duodenum_40",
    "epididymis_39", "epididymis_40",
    "heart_39", "heart_40",
    "hippocampus_39", "hippocampus_40",
    "hypothalamus_39", "hypothalamus_40",
    "ileum_39", "ileum_40",
    "jejunum_39", "jejunum_40",
    "kidney_39", "kidney_40",
    "liver_39", "liver_40",
    "lung_39", "lung_40",
    "lymph-node_39", "lymph-node_40",
    "mammary-gland_39", "mammary-gland_40",
    "medulla-oblongata_39", "medulla-oblongata_40",
    "midbrain_39", "midbrain_40",
    "muscle_39", "muscle_40",
    "omasum_39", "omasum_40",
    "optic-chiasm_39", "optic-chiasm_40",
    "ovary_39", "ovary_40",
    "oviduct_39", "oviduct_40",
    "pineal_39", "pineal_40",
    "pituitary_39", "pituitary_40",
    "pons_39", "pons_40",
    "rectum_39", "rectum_40",
    "reticulum_39", "reticulum_40",
    "rumen_39", "rumen_40",
    "skin_39", "skin_40",
    "soft-horn_39", "soft-horn_40",
    "spleen_39", "spleen_40",
    "splenium_39", "splenium_40",
    "testis_39", "testis_40",
    "thymus_39", "thymus_40",
    "thyroid_39", "thyroid_40",
]


# -----------------------------
# -----------------------------
REPLACEMENTS = {
    "cornua.uteri": "cornua-uteri",
    "medulla.oblongata": "medulla-oblongata",
    "optic.chiasm": "optic-chiasm",
    "cerebral.cortex": "cerebral-cortex",
    "corpus.uteri": "corpus-uteri",
    "mammary.gland": "mammary-gland",
    "bone.marrow": "bone-marrow",
    "soft.horn": "soft-horn",
    "lymph.node": "lymph-node",
}


def apply_replacements(s: str) -> str:
    """  '.'   '-'  , ."""
    out = s
    for old, new in REPLACEMENTS.items():
        out = out.replace(old, new)
    out = out.replace(".", "-")
    return out


def tissue_from_sample(sample: str) -> str:
    """abomasum_39 -> abomasum"""
    return sample.rsplit("_", 1)[0]


def rep_from_sample(sample: str) -> str:
    """abomasum_39 -> 39"""
    return sample.rsplit("_", 1)[1]


def canon_sample_name(col: str) -> Optional[str]:
    """
      enhancer CPM  (  bam file )  tissue_39/40.
     :
      H3K27ac_abomasum_39.bowtie2.mapped.filtered.sort.bam -> abomasum_39
    """
    c = apply_replacements(col)

    if c in EXPECTED_ORDER:
        return c

    m = re.search(r"H3K27ac_(.+?)_(39|40)\b", c)
    if m:
        tissue = m.group(1)
        rep = m.group(2)
        tissue = apply_replacements(tissue)
        return f"{tissue}_{rep}"

    m2 = re.search(r"([A-Za-z0-9\-]+)_(39|40)\b", c)
    if m2:
        tissue = apply_replacements(m2.group(1))
        rep = m2.group(2)
        candidate = f"{tissue}_{rep}"
        if candidate in EXPECTED_ORDER:
            return candidate

    return None


def load_tissue_colors(path: str) -> Dict[str, str]:
    """
    tissue_colors.tsv  :
    tissue  color
    cerebral-cortex #dcd71a
    ...
    """
    m: Dict[str, str] = {}
    if not path:
        return m
    df = pd.read_csv(path, sep=r"\s+", header=0)
    for _, r in df.iterrows():
        t = apply_replacements(str(r["tissue"]))
        m[t] = str(r["color"])
    return m


def load_expression(expr39: str, expr40: str) -> pd.DataFrame:
    """
    read (axis=1),  EXPECTED_ORDER  .
    file :
      ID  abomasum_39  adipose_39 ...
      A1BG 1.83 0.88 ...
    """
    g39 = pd.read_csv(expr39, sep="\t", header=0, index_col=0)
    g40 = pd.read_csv(expr40, sep="\t", header=0, index_col=0)

    g39.columns = [apply_replacements(c) for c in g39.columns]
    g40.columns = [apply_replacements(c) for c in g40.columns]

    geneall = pd.concat([g39, g40], axis=1, join="inner")

    cols = [c for c in EXPECTED_ORDER if c in geneall.columns]
    if len(cols) != len(EXPECTED_ORDER):
        missing = [c for c in EXPECTED_ORDER if c not in geneall.columns]
        if missing:
            print(f"[WARN] expression missing {len(missing)} samples, first few: {missing[:8]}")
    geneall = geneall[cols]
    return geneall


def read_enhancer_header(enh_cpm_csv: str) -> List[str]:
    """read enhancer CPM file (header), ."""
    with open(enh_cpm_csv, "r", newline="") as f:
        line = f.readline().strip("\n")
    return line.split(",")


def build_enh_col_index(enh_cpm_csv: str) -> Dict[str, int]:
    """
      enhancer CPM file header  :canon_sample ->  (  0  , ).
     :enh_cpm file :
      enhancer_id, v1, v2, ... v86
      v1   header[0],  i   values[i].
    """
    hdr = read_enhancer_header(enh_cpm_csv)
    canon_list: List[Optional[str]] = [canon_sample_name(h) for h in hdr]

    idx: Dict[str, int] = {}
    for i, c in enumerate(canon_list):
        if c is None:
            continue
        idx[c] = i

    missing = [s for s in EXPECTED_ORDER if s not in idx]
    if missing:
        raise ValueError(
            f"Enhancer CPM missing samples ({len(missing)}): {missing[:10]} ...\n"
            f"Hint: header canon failed. Check enh_cpm header first few: {hdr[:3]}"
        )
    return idx


def load_enhancer_rows(enh_cpm_csv: str, enhancers: Set[str]) -> Dict[str, np.ndarray]:
    """
      enhancer CPM file,  enhancers   86  (  EXPECTED_ORDER  ).
     : enhancer_id -> np.array([86])
    """
    col_idx = build_enh_col_index(enh_cpm_csv)
    need = set(enhancers)
    found: Dict[str, np.ndarray] = {}

    with open(enh_cpm_csv, "r", newline="") as f:
        reader = csv.reader(f)
        _ = next(reader)  #   header
        for row in reader:
            if not row:
                continue
            enh = row[0]
            if enh not in need:
                continue

            vals = row[1:]

            vec = np.empty(len(EXPECTED_ORDER), dtype=float)
            for j, s in enumerate(EXPECTED_ORDER):
                i = col_idx[s]
                try:
                    vec[j] = float(vals[i])
                except Exception:
                    vec[j] = np.nan

            found[enh] = vec
            if len(found) == len(need):
                break

    missing_enh = need - set(found.keys())
    if missing_enh:
        print(f"[WARN] {len(missing_enh)} enhancers not found in enh_cpm file, first few: {list(missing_enh)[:5]}")
    return found


def read_links_tsv(links_tsv: str) -> pd.DataFrame:
    """
    read confident links:
      awk output :
      enhancer  pearson_r  pval  gene  distance  qval
    """
    df = pd.read_csv(
        links_tsv,
        sep=r"\s+",
        header=None,
        usecols=[0, 1, 2, 3, 4, 5],
        names=["enhancer", "r_file", "p_file", "gene", "distance", "qval"],
    )
    return df


def choose_pairs(df_links: pd.DataFrame,
                 topN: int,
                 exclude_loc: bool,
                 pairs_file: str) -> pd.DataFrame:
    """
      enhancer-gene:
    -   pairs_file  :  pairs_file( :enhancer gene)
    -  :  df_links   r_file   topN(  LOC)
     :enhancer gene r_file p_file ...
    """
    if pairs_file:
        p = pd.read_csv(pairs_file, sep=r"\s+", header=None, usecols=[0, 1], names=["enhancer", "gene"])
        if exclude_loc:
            p = p[~p["gene"].astype(str).str.startswith("LOC")]
        out = p.merge(df_links, on=["enhancer", "gene"], how="left")
        return out

    d = df_links.copy()
    if exclude_loc:
        d = d[~d["gene"].astype(str).str.startswith("LOC")]

    d = d.sort_values("r_file", ascending=False).head(topN)
    return d


def sanitize_filename(s: str, max_len: int = 160) -> str:
    s2 = re.sub(r"[^A-Za-z0-9\-\._]+", "_", s)
    if len(s2) > max_len:
        s2 = s2[:max_len]
    return s2


def compute_degrees(df_links: pd.DataFrame, exclude_loc: bool) -> Tuple[pd.Series, pd.Series]:
    """
    enhancer_degree:   enhancer   gene( )
    gene_degree:   gene   enhancer( )
    """
    d = df_links.copy()
    if exclude_loc:
        d = d[~d["gene"].astype(str).str.startswith("LOC")]

    enh_deg = d.groupby("enhancer")["gene"].nunique().sort_values(ascending=False)
    gene_deg = d.groupby("gene")["enhancer"].nunique().sort_values(ascending=False)
    return enh_deg, gene_deg


def save_degree_plots(enh_deg: pd.Series,
                      gene_deg: pd.Series,
                      out_dir: str,
                      prefix: str,
                      dpi: int) -> None:
    """
    output  PDF:
      enhancer_degree.pdf
      gene_degree.pdf
     .
    """
    fig = plt.figure(figsize=(12, 2.8))
    ax = fig.add_subplot(111)
    ax.bar(np.arange(len(enh_deg)), enh_deg.values, color="black", width=1.0)
    ax.set_title("Enhancer degree (sorted)")
    ax.set_ylabel("#genes per enhancer")
    ax.set_xlabel("Enhancers (ranked)")
    ax.set_xlim(-1, len(enh_deg))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{prefix}.enhancer_degree.pdf"))
    plt.close(fig)

    fig = plt.figure(figsize=(4.5, 12))
    ax = fig.add_subplot(111)
    ax.barh(np.arange(len(gene_deg)), gene_deg.values, color="black", height=1.0)
    ax.set_title("Gene degree (sorted)")
    ax.set_xlabel("#enhancers per gene")
    ax.set_ylabel("Genes (ranked)")
    ax.set_ylim(-1, len(gene_deg))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"{prefix}.gene_degree.pdf"))
    plt.close(fig)


def label_points(ax,
                 x: np.ndarray,
                 y: np.ndarray,
                 labels: List[str],
                 top_x: int,
                 top_y: int,
                 fontsize: int = 7) -> None:
    """
     :
    - top_x:  x   N
    - top_y:  y   N
    """
    idxs: Set[int] = set()
    if top_x and top_x > 0:
        idxs |= set(np.argsort(x)[-top_x:].tolist())
    if top_y and top_y > 0:
        idxs |= set(np.argsort(y)[-top_y:].tolist())

    for i in sorted(idxs):
        ax.annotate(
            labels[i],
            xy=(x[i], y[i]),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=fontsize,
            ha="left",
            va="bottom",
        )


def plot_one_pair(enh: str,
                  gene: str,
                  x: np.ndarray,
                  y: np.ndarray,
                  sample_names: List[str],
                  tissue_color: Dict[str, str],
                  out_pdf: str,
                  point_size: float,
                  alpha: float,
                  label_top_x: int,
                  label_top_y: int,
                  add_reg_line: bool = True,
                  title_fontsize: int = 9) -> Tuple[float, float]:
    """
      enhancer-gene  :
    - 39 =  ,40 =
    -   tissue_colors.tsv
    -   + r,p
    -  (box aspect = 1)
    """
    m = np.isfinite(x) & np.isfinite(y)
    x2 = x[m]
    y2 = y[m]
    sn2 = [sample_names[i] for i in np.where(m)[0]]

    if len(x2) < 3:
        r, p = np.nan, np.nan
    else:
        r, p = pearsonr(x2, y2)

    fig = plt.figure(figsize=(3.6, 3.6))
    ax = fig.add_subplot(111)

    for xi, yi, sname in zip(x2, y2, sn2):
        tissue = tissue_from_sample(sname)
        rep = rep_from_sample(sname)

        col = tissue_color.get(tissue, "#666666")  #
        marker = "o" if rep == "39" else "^"       # 39  ,40

        ax.scatter(
            [xi], [yi],
            s=point_size,           #  :matplotlib   s  “ (pt^2)”
            c=col,
            marker=marker,
            alpha=alpha,
            edgecolors="none"
        )

    if add_reg_line and len(x2) >= 2:
        try:
            slope, intercept = np.polyfit(x2, y2, 1)
            xx = np.linspace(np.nanmin(x2), np.nanmax(x2), 100)
            yy = slope * xx + intercept
            ax.plot(xx, yy, linewidth=1.2)
        except Exception:
            pass

    txt = f"Pearson r = {r:.4f}\nP-value = {p:.3e}\nN = {len(x2)}"
    ax.text(
        0.02, 0.98, txt,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=8,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, linewidth=0.4)
    )

    labels = [tissue_from_sample(s) for s in sn2]
    label_points(ax, x2, y2, labels, label_top_x, label_top_y, fontsize=7)

    ax.set_xlabel("Enhancer H3K27ac signal (CPM)")
    ax.set_ylabel("Target gene expression (TPM)")
    ax.set_title(f"{enh}  →  {gene}", fontsize=title_fontsize)

    if hasattr(ax, "set_box_aspect"):
        ax.set_box_aspect(1)

    fig.tight_layout()
    fig.savefig(out_pdf)
    plt.close(fig)
    return r, p


def main():
    ap = argparse.ArgumentParser(
        description="Plot topN enhancer-gene pairs (exclude LOC) as separate PDFs, plus degree barplots."
    )

    ap.add_argument("--links_tsv", required=True,
                    help="confident links TSV (columns: enhancer r p gene distance qval)")
    ap.add_argument("--enh_cpm", required=True,
                    help="Enhancer CPM/CPMS CSV (H3K27ac_*_counts_cpms.csv)")
    ap.add_argument("--expr39", required=True,
                    help="Expression_tpm_39.tsv")
    ap.add_argument("--expr40", required=True,
                    help="Expression_tpm_40.tsv")
    ap.add_argument("--tissue_colors", default="",
                    help="tissue_colors.tsv (two columns: tissue color)")
    ap.add_argument("--out_dir", required=True,
                    help="output directory (will be created)")
    ap.add_argument("--prefix", default="E5_toplinks",
                    help="prefix for outputs (default: E5_toplinks)")

    ap.add_argument("--topN", type=int, default=8,
                    help="top N pairs by r (after excluding LOC). Ignored if --pairs provided.")
    ap.add_argument("--pairs", default="",
                    help="Optional pairs file (two columns: enhancer gene). If given, use these pairs to plot.")
    ap.add_argument("--exclude_loc", action="store_true",
                    help="Exclude genes starting with 'LOC'")

    ap.add_argument("--label_top_x", type=int, default=0,
                    help="Label top N points with largest X (enhancer CPM)")
    ap.add_argument("--label_top_y", type=int, default=0,
                    help="Label top N points with largest Y (gene TPM)")
    ap.add_argument("--label_top_y2", type=int, default=None,
                    help="Alias of --label_top_y (if set, overrides)")

    ap.add_argument("--point_size", type=float, default=20.0,
                    help="scatter point size (matplotlib 's' in pt^2)")
    ap.add_argument("--alpha", type=float, default=0.85,
                    help="scatter alpha")
    ap.add_argument("--dpi", type=int, default=200,
                    help="(reserved) dpi for raster; PDFs are vector and ignore dpi mostly.")

    args = ap.parse_args()

    if args.label_top_y2 is not None:
        args.label_top_y = args.label_top_y2

    os.makedirs(args.out_dir, exist_ok=True)

    links_df = read_links_tsv(args.links_tsv)

    enh_deg, gene_deg = compute_degrees(links_df, exclude_loc=args.exclude_loc)
    save_degree_plots(enh_deg, gene_deg, args.out_dir, args.prefix, dpi=args.dpi)

    pairs_df = choose_pairs(links_df, args.topN, args.exclude_loc, args.pairs)
    if pairs_df.empty:
        raise RuntimeError("No pairs selected. Check --exclude_loc / --topN / --pairs.")

    expr = load_expression(args.expr39, args.expr40)

    tcol = load_tissue_colors(args.tissue_colors)

    enh_set = set(pairs_df["enhancer"].astype(str).tolist())
    enh_rows = load_enhancer_rows(args.enh_cpm, enh_set)

    sample_names = [c for c in EXPECTED_ORDER if c in expr.columns]  # expr
    if len(sample_names) < 10:
        raise RuntimeError(f"Too few expression samples loaded: {len(sample_names)}")

    out_rows = []
    for i, row in pairs_df.reset_index(drop=True).iterrows():
        enh = str(row["enhancer"])
        gene = str(row["gene"])

        if gene not in expr.index:
            print(f"[WARN] gene not found in expression: {gene}, skip")
            continue
        if enh not in enh_rows:
            print(f"[WARN] enhancer not found in enhancer CPM: {enh}, skip")
            continue

        x = enh_rows[enh]
        y = expr.loc[gene, EXPECTED_ORDER].to_numpy(dtype=float)

        # output pdf
        fname = f"{args.prefix}.pair_{i+1:02d}.{sanitize_filename(enh)}__{sanitize_filename(gene)}.pdf"
        out_pdf = os.path.join(args.out_dir, fname)

        r_calc, p_calc = plot_one_pair(
            enh=enh,
            gene=gene,
            x=x,
            y=y,
            sample_names=EXPECTED_ORDER,
            tissue_color=tcol,
            out_pdf=out_pdf,
            point_size=args.point_size,
            alpha=args.alpha,
            label_top_x=args.label_top_x,
            label_top_y=args.label_top_y,
            add_reg_line=True
        )

        out_rows.append({
            "rank": i + 1,
            "enhancer": enh,
            "gene": gene,
            "r_file": row.get("r_file", np.nan),
            "p_file": row.get("p_file", np.nan),
            "r_calc_86": r_calc,
            "p_calc_86": p_calc,
            "out_pdf": out_pdf
        })

    out_df = pd.DataFrame(out_rows)
    out_df.to_csv(os.path.join(args.out_dir, f"{args.prefix}.pairs_summary.tsv"),
                  sep="\t", index=False)
    print(f"[OK] outputs saved in: {args.out_dir}")
    print(f"[OK] summary: {os.path.join(args.out_dir, f'{args.prefix}.pairs_summary.tsv')}")


if __name__ == "__main__":
    main()
