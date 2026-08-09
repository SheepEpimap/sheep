#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob

missing = []
try:
    import numpy as np
except Exception:
    missing.append("numpy")

try:
    import pandas as pd
except Exception:
    missing.append("pandas")

try:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
except Exception:
    missing.append("matplotlib")

try:
    from scipy.stats import pearsonr, spearmanr
except Exception:
    missing.append("scipy")

if missing:
    sys.stderr.write(
        "Missing Python packages: %s\n"
        "Install with:\n"
        "conda activate Graph_py\n"
        "conda install -y numpy pandas matplotlib scipy\n"
        % ", ".join(missing)
    )
    sys.exit(1)

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["svg.fonttype"] = "none"
mpl.rcParams["figure.dpi"] = 150
mpl.rcParams["savefig.dpi"] = 300

PSEUDOCOUNT = 0.5
ROOT_DIR = "/vol2/wulingyun/AS_cert/xiangguanxing/"
OUT_DIR = os.path.join(ROOT_DIR, "merged_39_40_vs_all_samples")
MIN_DEPTH = 5


def safe_read_tsv(path):
    return pd.read_csv(path, sep="\t", dtype=str, encoding="utf-8", engine="python")


def to_num(s):
    return pd.to_numeric(s, errors="coerce")


def pick_first_nonempty(series):
    vals = []
    for x in series.tolist():
        if pd.isna(x):
            continue
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            continue
        vals.append(s)
    return vals[0] if vals else ""


def pick_unique_join(series):
    vals = []
    for x in series.tolist():
        if pd.isna(x):
            continue
        s = str(x).strip()
        if s == "" or s.lower() == "nan":
            continue
        vals.append(s)
    vals = list(dict.fromkeys(vals))
    return ",".join(vals)


def allele_relation(q_ref, q_alt, s_ref, s_alt):
    q_ref = str(q_ref).strip().upper()
    q_alt = str(q_alt).strip().upper()
    s_ref = str(s_ref).strip().upper()
    s_alt = str(s_alt).strip().upper()

    if not q_ref or not q_alt or not s_ref or not s_alt:
        return "unknown"
    if q_ref == s_ref and q_alt == s_alt:
        return "same"
    if q_ref == s_alt and q_alt == s_ref:
        return "swapped"
    return "incompatible"


def calc_log2fc(ref_count, alt_count, pseudocount=PSEUDOCOUNT):
    return np.log2((alt_count + pseudocount) / (ref_count + pseudocount))


def major_allele_fraction(ref_count, alt_count):
    total = ref_count + alt_count
    with np.errstate(divide="ignore", invalid="ignore"):
        frac = np.maximum(ref_count, alt_count) / total
    return frac


def calc_mse(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() == 0:
        return np.nan
    return np.mean((x[mask] - y[mask]) ** 2)


def pretty_p(p):
    if pd.isna(p):
        return "NA"
    if p < 1e-300:
        return "<1e-300"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4g}"


def get_square_limits(x, y, pad_ratio=0.06):
    vals = np.concatenate([x[np.isfinite(x)], y[np.isfinite(y)]])
    if vals.size == 0:
        return -1.0, 1.0

    max_abs = np.nanmax(np.abs(vals))
    if not np.isfinite(max_abs) or max_abs == 0:
        max_abs = 1.0

    lim = max_abs * (1.0 + pad_ratio)
    return -lim, lim


def add_reference_shading(ax, lim_lo, lim_hi):
    ax.set_facecolor("white")
    shade_color = "#cb96dd"
    shade_alpha = 0.22

    ax.add_patch(
        Rectangle(
            (0, 0),
            lim_hi,
            lim_hi,
            facecolor=shade_color,
            edgecolor="none",
            alpha=shade_alpha,
            zorder=0
        )
    )

    ax.add_patch(
        Rectangle(
            (lim_lo, lim_lo),
            -lim_lo,
            -lim_lo,
            facecolor=shade_color,
            edgecolor="none",
            alpha=shade_alpha,
            zorder=0
        )
    )


def draw_corr_plot(
    df,
    xcol,
    ycol,
    out_png,
    out_pdf,
    title,
    method="pearson",
    xlabel=None,
    ylabel=None,
    jitter=0.04,
    pad_ratio=0.06
):
    x = df[xcol].values.astype(float)
    y = df[ycol].values.astype(float)

    if len(df) < 3:
        print("Skip %s: fewer than 3 points" % title)
        return None

    if method == "pearson":
        stat, pval = pearsonr(x, y)
        stat_label = "Pearson r"
    elif method == "spearman":
        stat, pval = spearmanr(x, y)
        stat_label = "Spearman rho"
    else:
        raise ValueError("method must be 'pearson' or 'spearman'")

    mse = calc_mse(x, y)

    x_plot = x.copy()
    y_plot = y.copy()
    if jitter and len(x_plot) > 0:
        rng = np.random.default_rng(12345)
        x_plot = x_plot + rng.normal(0, jitter, size=len(x_plot))
        y_plot = y_plot + rng.normal(0, jitter, size=len(y_plot))

    lo, hi = get_square_limits(x_plot, y_plot, pad_ratio=pad_ratio)

    fig, ax = plt.subplots(figsize=(6.2, 6.2))

    add_reference_shading(ax, lo, hi)

    ax.scatter(
        x_plot,
        y_plot,
        s=18,
        alpha=0.78,
        color="#6e4b3a",
        edgecolors="none",
        zorder=3
    )

    ax.plot(
        [lo, hi],
        [lo, hi],
        linestyle=(0, (8, 8)),
        linewidth=1.2,
        color="black",
        alpha=0.70,
        zorder=2
    )

    if method == "pearson" and len(np.unique(x)) >= 2:
        coef = np.polyfit(x, y, 1)
        fit_x = np.linspace(lo, hi, 200)
        fit_y = coef[0] * fit_x + coef[1]
        ax.plot(fit_x, fit_y, linewidth=1.2, color="#1f77b4", zorder=4)

    ax.axhline(0, color="#d0d0d0", linewidth=0.8, zorder=1)
    ax.axvline(0, color="#d0d0d0", linewidth=0.8, zorder=1)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    ax.set_xlabel(xlabel if xlabel else xcol, fontsize=14)
    ax.set_ylabel(ylabel if ylabel else ycol, fontsize=14)
    ax.set_title(title, fontsize=12)

    txt = (
        "N = %d\n"
        "%s = %.4f\n"
        "P = %s\n"
        "MSE = %.4f"
    ) % (len(df), stat_label, stat, pretty_p(pval), mse)

    ax.text(
        0.03, 0.97, txt,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=11,
        bbox=dict(
            boxstyle="round,pad=0.35",
            facecolor="white",
            edgecolor="#999999",
            alpha=0.82
        )
    )

    plt.tight_layout()
    plt.savefig(out_png, bbox_inches="tight")
    plt.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)

    return {
        "n": len(df),
        "stat": stat,
        "pval": pval,
        "mse": mse,
        "method": method,
        "xcol": xcol,
        "ycol": ycol,
        "title": title,
    }


def reduce_binomial_file(binomial_path, sample_name, source_group):
    df = safe_read_tsv(binomial_path)

    required = [
        "match_contig", "match_position",
        "query_refAllele", "query_altAllele",
        "query_refCount", "query_altCount"
    ]
    for c in required:
        if c not in df.columns:
            raise ValueError("Missing column in %s: %s" % (binomial_path, c))

    keep_cols = [
        "match_contig", "match_position",
        "query_refAllele", "query_altAllele",
        "query_refCount", "query_altCount", "query_totalCount",
        "query_rawDepth", "query_lowMAPQDepth", "query_lowBaseQDepth",
        "query_otherBases", "query_improperPairs",
        "query_table_line_numbers", "binomial_line_number"
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]

    df = df[keep_cols].copy()
    df["match_position"] = to_num(df["match_position"]).astype("Int64")
    df["sample_name"] = sample_name
    df["source_group"] = source_group

    for col in [
        "query_refCount", "query_altCount", "query_totalCount",
        "query_rawDepth", "query_lowMAPQDepth", "query_lowBaseQDepth",
        "query_otherBases", "query_improperPairs"
    ]:
        if col in df.columns:
            df[col] = to_num(df[col])

    agg_dict = {
        "query_refAllele": pick_first_nonempty,
        "query_altAllele": pick_first_nonempty,
        "query_refCount": "max",
        "query_altCount": "max",
        "sample_name": pick_first_nonempty,
        "source_group": pick_first_nonempty,
    }

    if "query_totalCount" in df.columns:
        agg_dict["query_totalCount"] = "max"
    if "query_rawDepth" in df.columns:
        agg_dict["query_rawDepth"] = "max"
    if "query_lowMAPQDepth" in df.columns:
        agg_dict["query_lowMAPQDepth"] = "max"
    if "query_lowBaseQDepth" in df.columns:
        agg_dict["query_lowBaseQDepth"] = "max"
    if "query_otherBases" in df.columns:
        agg_dict["query_otherBases"] = "max"
    if "query_improperPairs" in df.columns:
        agg_dict["query_improperPairs"] = "max"
    if "query_table_line_numbers" in df.columns:
        agg_dict["query_table_line_numbers"] = pick_unique_join
    if "binomial_line_number" in df.columns:
        agg_dict["binomial_line_number"] = pick_unique_join

    out = df.groupby(["match_contig", "match_position"], as_index=False).agg(agg_dict).copy()

    out["query_refCount"] = to_num(out["query_refCount"])
    out["query_altCount"] = to_num(out["query_altCount"])
    if "query_totalCount" in out.columns:
        out["query_totalCount"] = to_num(out["query_totalCount"])
    else:
        out["query_totalCount"] = out["query_refCount"] + out["query_altCount"]

    for col in [
        "query_rawDepth", "query_lowMAPQDepth", "query_lowBaseQDepth",
        "query_otherBases", "query_improperPairs"
    ]:
        if col in out.columns:
            out[col] = to_num(out[col])

    return out


def reduce_srr_file(srr_path, sample_name, source_group):
    df = safe_read_tsv(srr_path)

    required = [
        "match_contig", "match_position",
        "refAllele", "altAllele",
        "refCount", "altCount"
    ]
    for c in required:
        if c not in df.columns:
            raise ValueError("Missing column in %s: %s" % (srr_path, c))

    keep_cols = [
        "match_contig", "match_position",
        "refAllele", "altAllele",
        "refCount", "altCount", "totalCount",
        "rawDepth", "lowMAPQDepth", "lowBaseQDepth",
        "otherBases", "improperPairs",
        "source_line_number"
    ]
    keep_cols = [c for c in keep_cols if c in df.columns]

    df = df[keep_cols].copy()
    df["match_position"] = to_num(df["match_position"]).astype("Int64")
    df["sample_name"] = sample_name
    df["source_group"] = source_group

    for col in [
        "refCount", "altCount", "totalCount",
        "rawDepth", "lowMAPQDepth", "lowBaseQDepth",
        "otherBases", "improperPairs"
    ]:
        if col in df.columns:
            df[col] = to_num(df[col])

    agg_dict = {
        "refAllele": pick_first_nonempty,
        "altAllele": pick_first_nonempty,
        "refCount": "max",
        "altCount": "max",
        "sample_name": pick_first_nonempty,
        "source_group": pick_first_nonempty,
    }

    if "totalCount" in df.columns:
        agg_dict["totalCount"] = "max"
    if "rawDepth" in df.columns:
        agg_dict["rawDepth"] = "max"
    if "lowMAPQDepth" in df.columns:
        agg_dict["lowMAPQDepth"] = "max"
    if "lowBaseQDepth" in df.columns:
        agg_dict["lowBaseQDepth"] = "max"
    if "otherBases" in df.columns:
        agg_dict["otherBases"] = "max"
    if "improperPairs" in df.columns:
        agg_dict["improperPairs"] = "max"
    if "source_line_number" in df.columns:
        agg_dict["source_line_number"] = pick_unique_join

    out = df.groupby(["match_contig", "match_position"], as_index=False).agg(agg_dict).copy()

    out["refCount"] = to_num(out["refCount"])
    out["altCount"] = to_num(out["altCount"])
    if "totalCount" in out.columns:
        out["totalCount"] = to_num(out["totalCount"])
    else:
        out["totalCount"] = out["refCount"] + out["altCount"]

    for col in [
        "rawDepth", "lowMAPQDepth", "lowBaseQDepth",
        "otherBases", "improperPairs"
    ]:
        if col in out.columns:
            out[col] = to_num(out[col])

    return out


def find_group_files(group_dir):
    binomial_files = sorted(glob.glob(os.path.join(group_dir, "*.binomial_overlap_plus_query_reads.tsv")))
    srr_files = sorted(glob.glob(os.path.join(group_dir, "*.srr_overlap_rows.tsv")))
    if len(binomial_files) == 0 or len(srr_files) == 0:
        return None, None
    return binomial_files[0], srr_files[0]


def collect_all_inputs(root_dir):
    query_parts = []
    sample_parts = []

    sample_dirs = sorted(glob.glob(os.path.join(root_dir, "SRR*")))
    sample_dirs = [x for x in sample_dirs if os.path.isdir(x)]

    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir)

        for group_name in ["39", "40"]:
            group_dir = os.path.join(sample_dir, group_name)
            if not os.path.isdir(group_dir):
                continue

            binomial_path, srr_path = find_group_files(group_dir)
            if binomial_path is None or srr_path is None:
                continue

            q_df = reduce_binomial_file(binomial_path, sample_name=sample_name, source_group=group_name)
            s_df = reduce_srr_file(srr_path, sample_name=sample_name, source_group=group_name)

            query_parts.append(q_df)
            sample_parts.append(s_df)

    return query_parts, sample_parts


def build_global_merged_plot_table(root_dir=ROOT_DIR, min_depth=MIN_DEPTH):
    query_parts, sample_parts = collect_all_inputs(root_dir)

    if len(query_parts) == 0 or len(sample_parts) == 0:
        raise FileNotFoundError("No usable 39/40 overlap files found under sample directories.")

    q_all = pd.concat(query_parts, axis=0, ignore_index=True)

    q_step1 = q_all.groupby(
        ["source_group", "match_contig", "match_position"],
        as_index=False
    ).agg({
        "query_refAllele": pick_first_nonempty,
        "query_altAllele": pick_first_nonempty,
        "query_refCount": "max",
        "query_altCount": "max",
        "query_totalCount": "max",
        **({"query_rawDepth": "max"} if "query_rawDepth" in q_all.columns else {}),
        **({"query_lowMAPQDepth": "max"} if "query_lowMAPQDepth" in q_all.columns else {}),
        **({"query_lowBaseQDepth": "max"} if "query_lowBaseQDepth" in q_all.columns else {}),
        **({"query_otherBases": "max"} if "query_otherBases" in q_all.columns else {}),
        **({"query_improperPairs": "max"} if "query_improperPairs" in q_all.columns else {}),
    }).copy()

    q_merged = q_step1.groupby(
        ["match_contig", "match_position"],
        as_index=False
    ).agg({
        "query_refAllele": pick_first_nonempty,
        "query_altAllele": pick_first_nonempty,
        "query_refCount": "sum",
        "query_altCount": "sum",
        "query_totalCount": "sum",
        "source_group": pick_unique_join,
        **({"query_rawDepth": "sum"} if "query_rawDepth" in q_step1.columns else {}),
        **({"query_lowMAPQDepth": "sum"} if "query_lowMAPQDepth" in q_step1.columns else {}),
        **({"query_lowBaseQDepth": "sum"} if "query_lowBaseQDepth" in q_step1.columns else {}),
        **({"query_otherBases": "sum"} if "query_otherBases" in q_step1.columns else {}),
        **({"query_improperPairs": "sum"} if "query_improperPairs" in q_step1.columns else {}),
    }).copy()

    s_all = pd.concat(sample_parts, axis=0, ignore_index=True)

    s_step1 = s_all.groupby(
        ["sample_name", "match_contig", "match_position"],
        as_index=False
    ).agg({
        "refAllele": pick_first_nonempty,
        "altAllele": pick_first_nonempty,
        "refCount": "max",
        "altCount": "max",
        "totalCount": "max",
        **({"rawDepth": "max"} if "rawDepth" in s_all.columns else {}),
        **({"lowMAPQDepth": "max"} if "lowMAPQDepth" in s_all.columns else {}),
        **({"lowBaseQDepth": "max"} if "lowBaseQDepth" in s_all.columns else {}),
        **({"otherBases": "max"} if "otherBases" in s_all.columns else {}),
        **({"improperPairs": "max"} if "improperPairs" in s_all.columns else {}),
    }).copy()

    s_merged = s_step1.groupby(
        ["match_contig", "match_position"],
        as_index=False
    ).agg({
        "refAllele": pick_first_nonempty,
        "altAllele": pick_first_nonempty,
        "refCount": "sum",
        "altCount": "sum",
        "totalCount": "sum",
        "sample_name": "nunique",
        **({"rawDepth": "sum"} if "rawDepth" in s_step1.columns else {}),
        **({"lowMAPQDepth": "sum"} if "lowMAPQDepth" in s_step1.columns else {}),
        **({"lowBaseQDepth": "sum"} if "lowBaseQDepth" in s_step1.columns else {}),
        **({"otherBases": "sum"} if "otherBases" in s_step1.columns else {}),
        **({"improperPairs": "sum"} if "improperPairs" in s_step1.columns else {}),
    }).copy().rename(columns={"sample_name": "n_samples_contributed"})

    merged = pd.merge(
        q_merged,
        s_merged,
        on=["match_contig", "match_position"],
        how="inner"
    ).copy()

    merged["relation"] = merged.apply(
        lambda r: allele_relation(
            r.get("query_refAllele", ""),
            r.get("query_altAllele", ""),
            r.get("refAllele", ""),
            r.get("altAllele", "")
        ),
        axis=1
    )

    merged = merged[merged["relation"].isin(["same", "swapped"])].copy()

    merged = merged[
        (merged["query_totalCount"] >= min_depth) &
        (merged["totalCount"] >= min_depth)
    ].copy()

    merged["query_log2FC"] = calc_log2fc(
        merged["query_refCount"].astype(float),
        merged["query_altCount"].astype(float)
    )
    merged["sample_log2FC_raw"] = calc_log2fc(
        merged["refCount"].astype(float),
        merged["altCount"].astype(float)
    )
    merged["sample_log2FC_aligned"] = np.where(
        merged["relation"] == "same",
        merged["sample_log2FC_raw"],
        -merged["sample_log2FC_raw"]
    )

    merged["query_abs_log2FC"] = np.abs(merged["query_log2FC"])
    merged["sample_abs_log2FC"] = np.abs(merged["sample_log2FC_aligned"])

    merged["query_major_frac"] = major_allele_fraction(
        merged["query_refCount"].astype(float),
        merged["query_altCount"].astype(float)
    )
    merged["sample_major_frac"] = major_allele_fraction(
        merged["refCount"].astype(float),
        merged["altCount"].astype(float)
    )

    merged = merged.replace([np.inf, -np.inf], np.nan)
    merged = merged.dropna(
        subset=[
            "query_log2FC", "sample_log2FC_aligned",
            "query_abs_log2FC", "sample_abs_log2FC",
            "query_major_frac", "sample_major_frac"
        ]
    ).copy()

    return merged, q_merged, s_merged


def save_all_outputs(plot_df, query_df, sample_df, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    query_path = os.path.join(out_dir, "merged_query_39_40_sumdedup.tsv")
    sample_path = os.path.join(out_dir, "merged_all_samples_sumdedup.tsv")
    plot_path = os.path.join(out_dir, "merged_query_vs_all_samples_plotting_table.tsv")

    query_df.to_csv(query_path, sep="\t", index=False)
    sample_df.to_csv(sample_path, sep="\t", index=False)
    plot_df.to_csv(plot_path, sep="\t", index=False)

    title_base = "Merged 39+40 vs merged all samples"

    results = []

    r = draw_corr_plot(
        plot_df, "query_log2FC", "sample_log2FC_aligned",
        os.path.join(out_dir, "pearson_signed_log2FC.png"),
        os.path.join(out_dir, "pearson_signed_log2FC.pdf"),
        title="%s | Pearson | signed log2FC" % title_base,
        method="pearson",
        xlabel=r"$\log_2$FC (This study)",
        ylabel=r"$\log_2$FC (Yan, Z. et al.)"
    )
    if r is not None:
        results.append(r)

    r = draw_corr_plot(
        plot_df, "query_log2FC", "sample_log2FC_aligned",
        os.path.join(out_dir, "spearman_signed_log2FC.png"),
        os.path.join(out_dir, "spearman_signed_log2FC.pdf"),
        title="%s | Spearman | signed log2FC" % title_base,
        method="spearman",
        xlabel=r"$\log_2$FC (This study)",
        ylabel=r"$\log_2$FC (Yan, Z. et al.)"
    )
    if r is not None:
        results.append(r)

    r = draw_corr_plot(
        plot_df, "query_abs_log2FC", "sample_abs_log2FC",
        os.path.join(out_dir, "pearson_abs_log2FC.png"),
        os.path.join(out_dir, "pearson_abs_log2FC.pdf"),
        title="%s | Pearson | abs(log2FC)" % title_base,
        method="pearson",
        xlabel=r"abs($\log_2$FC) (This study)",
        ylabel=r"abs($\log_2$FC) (Yan, Z. et al.)"
    )
    if r is not None:
        results.append(r)

    r = draw_corr_plot(
        plot_df, "query_abs_log2FC", "sample_abs_log2FC",
        os.path.join(out_dir, "spearman_abs_log2FC.png"),
        os.path.join(out_dir, "spearman_abs_log2FC.pdf"),
        title="%s | Spearman | abs(log2FC)" % title_base,
        method="spearman",
        xlabel=r"abs($\log_2$FC) (This study)",
        ylabel=r"abs($\log_2$FC) (Yan, Z. et al.)"
    )
    if r is not None:
        results.append(r)

    r = draw_corr_plot(
        plot_df, "query_major_frac", "sample_major_frac",
        os.path.join(out_dir, "pearson_major_frac.png"),
        os.path.join(out_dir, "pearson_major_frac.pdf"),
        title="%s | Pearson | major-allele fraction" % title_base,
        method="pearson",
        xlabel="Major-allele fraction (This study)",
        ylabel="Major-allele fraction (Yan, Z. et al.)"
    )
    if r is not None:
        results.append(r)

    r = draw_corr_plot(
        plot_df, "query_major_frac", "sample_major_frac",
        os.path.join(out_dir, "spearman_major_frac.png"),
        os.path.join(out_dir, "spearman_major_frac.pdf"),
        title="%s | Spearman | major-allele fraction" % title_base,
        method="spearman",
        xlabel="Major-allele fraction (This study)",
        ylabel="Major-allele fraction (Yan, Z. et al.)"
    )
    if r is not None:
        results.append(r)

    summary_path = os.path.join(out_dir, "summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("root_dir\t%s\n" % ROOT_DIR)
        f.write("out_dir\t%s\n" % OUT_DIR)
        f.write("min_depth\t%s\n" % MIN_DEPTH)
        f.write("n_query_sites_after_merge\t%s\n" % query_df.shape[0])
        f.write("n_sample_sites_after_merge\t%s\n" % sample_df.shape[0])
        f.write("n_sites_used_for_plot\t%s\n" % plot_df.shape[0])

    metric_summary_path = os.path.join(out_dir, "plot_metrics.tsv")
    if len(results) > 0:
        pd.DataFrame(results).to_csv(metric_summary_path, sep="\t", index=False)

    print("Done")
    print("query:", query_path)
    print("sample:", sample_path)
    print("plot:", plot_path)
    print("summary:", summary_path)
    if len(results) > 0:
        print("metrics:", metric_summary_path)


def main():
    plot_df, query_df, sample_df = build_global_merged_plot_table(
        root_dir=ROOT_DIR,
        min_depth=MIN_DEPTH
    )

    if plot_df.shape[0] < 3:
        raise ValueError("Too few plottable points after global merge: %d" % plot_df.shape[0])

    save_all_outputs(plot_df, query_df, sample_df, OUT_DIR)


if __name__ == "__main__":
    main()