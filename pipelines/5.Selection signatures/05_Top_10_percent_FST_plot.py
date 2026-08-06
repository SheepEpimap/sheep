#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pysam

# =========================
# 1. Input files and parameters
# =========================
vcf_file = "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/vcf/chrAuto_ancient.vcf.gz"
fst_file = "/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AA_selection_signature/vcf/chrAuto_neweurope_oldeurope_point.windowed.weir.fst"
region = "chr8:90340000-90380000"
out_prefix = "/vol2/mengzhu/jupyter/figure/Selection_signature"

# Smoothing window: rolling mean by number of SNPs
smooth_window = 15

# Label the positions with the highest FST values
top_n = 5


# =========================
# 2. Utility functions
# =========================
def normalize_chrom(chrom: str) -> str:
    """
    Normalize chromosome names by removing the chr prefix
    chr8 -> 8
    8    -> 8
    """
    chrom = str(chrom).strip()
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    return chrom


def parse_region(region_str: str):
    """
    Parse a region string
    Example: chr8:90340000-90380000
    Return: chrom, start, end
    """
    region_str = region_str.replace(",", "").strip()
    chrom_part, pos_part = region_str.split(":")
    start_str, end_str = pos_part.split("-")
    chrom = normalize_chrom(chrom_part)
    start = int(start_str)
    end = int(end_str)
    return chrom, start, end


def safe_float(x):
    """
    Safely convert a value to a floating-point number
    """
    if x is None:
        return np.nan
    if isinstance(x, (tuple, list)):
        if len(x) == 0:
            return np.nan
        x = x[0]
    try:
        return float(x)
    except Exception:
        return np.nan


# =========================
# 3. Read regional FST data
# =========================
def read_fst_region(fst_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Read the FST file and filter the target region
    File format:
    CHROM   POS     WEIR_AND_COCKERHAM_FST
    """
    df = pd.read_csv(
        fst_path,
        sep=r"\s+",
        engine="python",
        dtype={"CHROM": str}
    )

    required_cols = ["CHROM", "POS", "WEIR_AND_COCKERHAM_FST"]
    for col in required_cols:
        if col not in df.columns:
            raise ValueError(f"The FST file is missing column: {col}")

    df["CHROM_norm"] = df["CHROM"].map(normalize_chrom)
    df["POS"] = pd.to_numeric(df["POS"], errors="coerce")
    df["FST"] = pd.to_numeric(df["WEIR_AND_COCKERHAM_FST"], errors="coerce")

    region_df = df[
        (df["CHROM_norm"] == chrom) &
        (df["POS"] >= start) &
        (df["POS"] <= end)
    ][["CHROM_norm", "POS", "FST"]].copy()

    region_df.columns = ["CHROM", "POS", "FST"]
    region_df = region_df.sort_values("POS").reset_index(drop=True)
    return region_df


# =========================
# 4. Read regional VCF data using the .csi index
# =========================
def read_vcf_region(vcf_path: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    """
    Read the VCF by region with pysam
    Support both contig naming conventions: 8 and chr8
    """
    rows = []
    vcf = pysam.VariantFile(vcf_path)

    # Read contig names from the VCF header to determine whether names use 8 or chr8
    contigs_in_vcf = set(vcf.header.contigs.keys())

    candidates = [chrom, f"chr{chrom}"]
    fetch_contig = None
    for c in candidates:
        if c in contigs_in_vcf:
            fetch_contig = c
            break

    if fetch_contig is None:
        raise ValueError(
            f"The VCF does not contain chromosome {chrom} or chr{chrom}."
            f" Example contigs from the header: {list(contigs_in_vcf)[:10]}"
        )

    # pysam fetch: start is 0-based and stop defines a half-open interval
    for rec in vcf.fetch(fetch_contig, start - 1, end):
        rows.append({
            "CHROM": normalize_chrom(rec.chrom),
            "POS": rec.pos,
            "REF": rec.ref,
            "ALT": ",".join(rec.alts) if rec.alts else ".",
            "AF": safe_float(rec.info.get("AF")),
            "RAF": safe_float(rec.info.get("RAF")),
            "INFO_SCORE": safe_float(rec.info.get("INFO")),
            "MAF": safe_float(rec.info.get("MAF"))
        })

    vcf.close()

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("POS").reset_index(drop=True)
    return df


# =========================
# 5. Plotting
# =========================
def plot_region(df: pd.DataFrame, region_str: str, out_prefix: str,
                smooth_window: int = 15, top_n: int = 5):
    """
    Three-panel plot:
    1) FST scatter plot and smoothed curve
    2) INFO scatter plot
    3) MAF scatter plot
    """
    if df.empty:
        raise ValueError("No data are available for plotting.")

    df = df.sort_values("POS").copy()

    win = min(smooth_window, max(3, len(df)))
    df["FST_SMOOTH"] = df["FST"].rolling(
        window=win, center=True, min_periods=1
    ).mean()

    label_df = df.dropna(subset=["FST"]).nlargest(min(top_n, len(df)), "FST").copy()

    fig, axes = plt.subplots(
        3, 1,
        figsize=(12, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [3.5, 1.2, 1.2]}
    )

    ax1, ax2, ax3 = axes

    # ---------- panel 1: FST ----------
    ax1.scatter(df["POS"], df["FST"], s=22, alpha=0.8, label="Per-SNP FST")
    ax1.plot(df["POS"], df["FST_SMOOTH"], linewidth=2, label=f"Rolling mean (n={win})")
    ax1.axhline(0, linestyle="--", linewidth=1)
    ax1.set_ylabel("Weir & Cockerham FST")
    ax1.set_title(f"Local fine-scale selection signal: {region_str}")

    for _, row in label_df.iterrows():
        ax1.scatter(row["POS"], row["FST"], s=40, zorder=3)
        ax1.annotate(
            f"{int(row['POS'])}\n{row['FST']:.3f}",
            xy=(row["POS"], row["FST"]),
            xytext=(0, 8),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=8
        )

    ax1.legend(frameon=False, loc="upper right")

    # ---------- panel 2: INFO ----------
    if "INFO_SCORE" in df.columns and df["INFO_SCORE"].notna().any():
        ax2.scatter(df["POS"], df["INFO_SCORE"], s=18, alpha=0.8)
        ax2.set_ylabel("INFO")
        ymin = max(0, np.nanmin(df["INFO_SCORE"]) - 0.05)
        ymax = min(1.05, np.nanmax(df["INFO_SCORE"]) + 0.05)
        if np.isfinite(ymin) and np.isfinite(ymax) and ymin < ymax:
            ax2.set_ylim(ymin, ymax)
    else:
        ax2.text(0.5, 0.5, "No INFO field", transform=ax2.transAxes,
                 ha="center", va="center")
        ax2.set_ylabel("INFO")

    # ---------- panel 3: MAF ----------
    if "MAF" in df.columns and df["MAF"].notna().any():
        ax3.scatter(df["POS"], df["MAF"], s=18, alpha=0.8)
        ax3.set_ylabel("MAF")
        ax3.set_ylim(0, 0.55)
    else:
        ax3.text(0.5, 0.5, "No MAF field", transform=ax3.transAxes,
                 ha="center", va="center")
        ax3.set_ylabel("MAF")

    chrom, start, end = parse_region(region_str)
    ax3.set_xlim(start, end)
    ax3.set_xlabel("Genomic position")

    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300, bbox_inches="tight")
    plt.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    plt.close()


# =========================
# 6. Main workflow
# =========================
def main():
    chrom, start, end = parse_region(region)

    print("Reading FST...")
    fst_df = read_fst_region(fst_file, chrom, start, end)
    if fst_df.empty:
        raise ValueError(f"Region  {region}  has no sites in the FST file.")

    print("Reading VCF with CSI index...")
    vcf_df = read_vcf_region(vcf_file, chrom, start, end)

    print("Merging...")
    merged = fst_df.merge(vcf_df, on=["CHROM", "POS"], how="left")
    merged = merged.sort_values(["CHROM", "POS"]).reset_index(drop=True)

    # Save the regional data table
    merged.to_csv(f"{out_prefix}.region_data.tsv", sep="\t", index=False)

    print("Plotting...")
    plot_region(
        df=merged,
        region_str=region,
        out_prefix=out_prefix,
        smooth_window=smooth_window,
        top_n=top_n
    )

    print("Done!")
    print(f"Output table : {out_prefix}.region_data.tsv")
    print(f"Output figure: {out_prefix}.png")
    print(f"Output figure: {out_prefix}.pdf")


if __name__ == "__main__":
    main()
