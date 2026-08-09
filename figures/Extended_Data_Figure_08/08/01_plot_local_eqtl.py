#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
Local GWAS/eQTL Manhattan plot with optional per-panel LD coloring.

Key logic:
1. --plot-mode gwas / eqtl / both.
2. noLD mode:
   - ordinary local Manhattan plot;
   - user-specified SNPs are highlighted as red circles.
3. LD mode:
   - GWAS panel uses the minimum-P GWAS SNP within the window as its LD lead;
   - eQTL panel uses the minimum-P eQTL SNP within the window as its LD lead;
   - each panel is colored by LD r2 to its own lead SNP;
   - user-specified SNPs keep their LD color but are overdrawn with a black outline.
4. If --plot-mode eqtl is used, LD is automatically disabled.
5. The minimum-P LD lead SNP is shown as a purple diamond.
6. rsID labels can be filled from an Ensembl/EVA GVF file using --rsid-gvf.
7. LD mode removes SNPs with NA r2 directly; NA points and NA legend are not shown.
8. LD cache is saved automatically in <outdir>/LD_tmp as *.ld_cache.tsv.gz.
9. Output is editable PDF only.
"""

import argparse
import gzip
import hashlib
import math
import re
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
mpl.rcParams["axes.linewidth"] = 1.0
mpl.rcParams["xtick.major.width"] = 0.9
mpl.rcParams["ytick.major.width"] = 0.9
mpl.rcParams["xtick.direction"] = "out"
mpl.rcParams["ytick.direction"] = "out"

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MaxNLocator
from matplotlib.patches import Rectangle


# ============================================================
# Basic parsers
# ============================================================

def normalize_chr(chrom):
    chrom = str(chrom).strip()
    chrom = re.sub(r"^chr", "", chrom, flags=re.IGNORECASE)
    return chrom


def parse_site(site):
    site = str(site).strip().replace(",", "")
    m = re.match(r"^(chr)?([^:]+):(\d+)(?:-(\d+))?$", site, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Cannot parse site: {site}")
    return normalize_chr(m.group(2)), int(m.group(3))


def parse_region(region):
    region = str(region).strip().replace(",", "")
    m = re.match(r"^(chr)?([^:]+):(\d+)-(\d+)$", region, flags=re.IGNORECASE)
    if not m:
        raise ValueError(f"Cannot parse region: {region}")

    chrom = normalize_chr(m.group(2))
    start = int(m.group(3))
    end = int(m.group(4))

    if start > end:
        start, end = end, start

    return chrom, start, end


def split_multi_value(x):
    if x is None:
        return []

    x = str(x).strip()

    if x == "":
        return []

    return [p for p in re.split(r"[,\s;]+", x) if p != ""]


def parse_mark_sites(sites, rsids=None):
    site_list = split_multi_value(sites)
    rsid_list = split_multi_value(rsids)

    if len(site_list) == 0:
        raise ValueError("No sites were provided.")

    records = []

    for i, site in enumerate(site_list):
        chrom, pos = parse_site(site)
        label = rsid_list[i] if i < len(rsid_list) else f"chr{chrom}:{pos}"

        records.append({
            "chrom": chrom,
            "pos": pos,
            "label": label
        })

    return pd.DataFrame(records).drop_duplicates(subset=["chrom", "pos", "label"])


def update_mark_labels_from_gvf(mark_df, gvf_file):
    if mark_df is None or mark_df.empty or gvf_file is None:
        return mark_df

    gvf_file = Path(gvf_file)

    if not gvf_file.exists():
        raise FileNotFoundError(f"GVF file not found: {gvf_file}")

    targets = set(
        (normalize_chr(r["chrom"]), int(r["pos"]))
        for _, r in mark_df.iterrows()
    )

    found = {}

    opener = gzip.open if str(gvf_file).endswith(".gz") else open

    with opener(gvf_file, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")

            if len(parts) < 9:
                continue

            chrom = normalize_chr(parts[0])

            try:
                pos = int(parts[3])
            except ValueError:
                continue

            key = (chrom, pos)

            if key not in targets:
                continue

            m = re.search(r"(rs\d+)", parts[8])

            if m:
                found[key] = m.group(1)

            if len(found) == len(targets):
                break

    out = mark_df.copy()
    out["label"] = [
        found.get((normalize_chr(r["chrom"]), int(r["pos"])), r["label"])
        for _, r in out.iterrows()
    ]

    return out


def parse_variant_id_series(variant_series):
    s = variant_series.astype(str)
    split_df = s.str.rsplit("_", n=1, expand=True)

    if split_df.shape[1] != 2:
        raise ValueError("variant_id cannot be parsed by pattern chrom_pos")

    chrom = split_df[0].map(normalize_chr)
    pos = pd.to_numeric(split_df[1], errors="coerce").astype("Int64")

    return chrom, pos


def safe_name(x):
    x = str(x)

    for old in [":", "-", "/", "\\", " ", "|"]:
        x = x.replace(old, "_")

    return x


def add_logp(df, p_col="_p"):
    df = df.copy()
    df[p_col] = pd.to_numeric(df[p_col], errors="coerce")
    df = df.dropna(subset=[p_col])
    df = df[(df[p_col] >= 0) & (df[p_col] <= 1)].copy()

    if df.empty:
        return df

    df.loc[df[p_col] == 0, p_col] = 1e-300
    df["_logp"] = -np.log10(df[p_col].astype(float))
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(subset=["_logp"])

    return df


def resolve_column(columns, explicit_col, candidates, label):
    columns = list(columns)

    if explicit_col is not None and explicit_col != "auto":
        if explicit_col in columns:
            return explicit_col

        lower_map = {str(c).lower(): c for c in columns}

        if str(explicit_col).lower() in lower_map:
            return lower_map[str(explicit_col).lower()]

        raise ValueError(
            f"Column '{explicit_col}' for {label} not found. "
            f"Available columns: {columns}"
        )

    lower_map = {str(c).lower(): c for c in columns}

    for cand in candidates:
        if cand in columns:
            return cand

        if cand.lower() in lower_map:
            return lower_map[cand.lower()]

    raise ValueError(
        f"Cannot auto-detect column for {label}. "
        f"Candidates: {candidates}. Available columns: {columns}"
    )


# ============================================================
# Window
# ============================================================

def determine_window(
    sites=None,
    enhancer=None,
    region=None,
    flank_mult=5.0,
    window_bp=500000
):
    enhancer_info = None

    if region is not None:
        chrom, start, end = parse_region(region)
        return chrom, start, end, enhancer_info

    if enhancer is not None:
        chrom, enh_start, enh_end = parse_region(enhancer)

        enhancer_len = max(1, enh_end - enh_start)
        flank = int(round(enhancer_len * flank_mult))

        start = max(0, enh_start - flank)
        end = enh_end + flank

        enhancer_info = {
            "chrom": chrom,
            "start": enh_start,
            "end": enh_end
        }

        if sites is not None:
            mark_df = parse_mark_sites(sites)
            mark_df["chrom"] = mark_df["chrom"].map(normalize_chr)
            mark_df = mark_df[mark_df["chrom"] == normalize_chr(chrom)].copy()

            if not mark_df.empty:
                start = min(start, max(0, int(mark_df["pos"].min()) - flank))
                end = max(end, int(mark_df["pos"].max()) + flank)

        return chrom, start, end, enhancer_info

    if sites is not None:
        mark_df = parse_mark_sites(sites)
        chrom_set = set(mark_df["chrom"].map(normalize_chr).tolist())

        if len(chrom_set) != 1:
            raise ValueError(f"Marked sites contain multiple chromosomes: {chrom_set}")

        chrom = list(chrom_set)[0]
        start = max(0, int(mark_df["pos"].min()) - int(window_bp))
        end = int(mark_df["pos"].max()) + int(window_bp)

        return chrom, start, end, enhancer_info

    raise ValueError("At least one of --region, --enhancer, or --site/--sites is required.")


# ============================================================
# Read GWAS
# ============================================================

def read_gwas_local(
    gwas_file,
    chrom,
    window_start,
    window_end,
    gwas_chr_col="auto",
    gwas_pos_col="auto",
    gwas_p_col="auto",
    gwas_snp_col="auto",
    chunksize=1000000
):
    gwas_file = Path(gwas_file)

    if not gwas_file.exists():
        raise FileNotFoundError(f"GWAS file not found: {gwas_file}")

    chrom = normalize_chr(chrom)
    local_chunks = []

    reader = pd.read_csv(
        gwas_file,
        sep=r"\s+",
        compression="infer",
        chunksize=chunksize,
        low_memory=False
    )

    for chunk in reader:
        chr_col = resolve_column(
            chunk.columns,
            gwas_chr_col,
            ["CHR", "Chr", "chr", "CHROM", "chrom"],
            "GWAS chromosome"
        )

        pos_col = resolve_column(
            chunk.columns,
            gwas_pos_col,
            ["BP", "bp", "POS", "pos", "Position", "position"],
            "GWAS position"
        )

        p_col = resolve_column(
            chunk.columns,
            gwas_p_col,
            ["P", "p", "PVAL", "pval", "p_value", "pvalue"],
            "GWAS P value"
        )

        try:
            snp_col = resolve_column(
                chunk.columns,
                gwas_snp_col,
                ["SNP", "snp", "ID", "id", "variant_id", "MarkerName"],
                "GWAS SNP ID"
            )
        except ValueError:
            snp_col = None

        tmp = chunk.copy()
        tmp["_chr"] = tmp[chr_col].map(normalize_chr)
        tmp["_pos"] = pd.to_numeric(tmp[pos_col], errors="coerce").astype("Int64")
        tmp["_p"] = pd.to_numeric(tmp[p_col], errors="coerce")

        if snp_col is not None:
            tmp["_id"] = tmp[snp_col].astype(str)
        else:
            tmp["_id"] = tmp["_chr"].astype(str) + "_" + tmp["_pos"].astype(str)

        mask = (
            (tmp["_chr"] == chrom) &
            (tmp["_pos"].notna()) &
            (tmp["_pos"] >= window_start) &
            (tmp["_pos"] <= window_end)
        )

        sub = tmp.loc[mask, ["_chr", "_pos", "_p", "_id"]].copy()

        if not sub.empty:
            local_chunks.append(sub)

    if len(local_chunks) == 0:
        return pd.DataFrame(columns=["_chr", "_pos", "_p", "_id", "_logp"])

    return add_logp(pd.concat(local_chunks, ignore_index=True), p_col="_p")


# ============================================================
# Read eQTL
# ============================================================

def get_eqtl_file(
    eqtl_dir,
    tissue,
    chrom,
    eqtl_pattern="{tissue}.cis_qtl_pairs.chr{chrom}.txt.gz"
):
    filename = eqtl_pattern.format(tissue=tissue, chrom=normalize_chr(chrom))
    return Path(eqtl_dir) / filename


def read_eqtl_local(
    eqtl_dir,
    tissue,
    chrom,
    window_start,
    window_end,
    gene=None,
    gene_filter="exact",
    eqtl_pattern="{tissue}.cis_qtl_pairs.chr{chrom}.txt.gz",
    eqtl_gene_col="phenotype_id",
    eqtl_variant_col="variant_id",
    eqtl_p_col="pval_nominal",
    chunksize=1000000
):
    if eqtl_dir is None:
        raise ValueError("--eqtl-dir is required when plot-mode includes eQTL.")

    eqtl_file = get_eqtl_file(eqtl_dir, tissue, chrom, eqtl_pattern)

    if not eqtl_file.exists():
        raise FileNotFoundError(f"eQTL file not found: {eqtl_file}")

    chrom = normalize_chr(chrom)
    all_local_chunks = []
    gene_local_chunks = []

    reader = pd.read_csv(
        eqtl_file,
        sep=r"\s+",
        compression="infer",
        chunksize=chunksize,
        low_memory=False
    )

    for chunk in reader:
        for col in [eqtl_gene_col, eqtl_variant_col, eqtl_p_col]:
            if col not in chunk.columns:
                raise ValueError(
                    f"Column '{col}' not found in eQTL file. "
                    f"Available columns: {list(chunk.columns)}"
                )

        tmp = chunk.copy()
        var_chr, var_pos = parse_variant_id_series(tmp[eqtl_variant_col])

        tmp["_chr"] = var_chr
        tmp["_pos"] = var_pos
        tmp["_p"] = pd.to_numeric(tmp[eqtl_p_col], errors="coerce")
        tmp["_id"] = tmp[eqtl_variant_col].astype(str)
        tmp["_gene"] = tmp[eqtl_gene_col].astype(str)

        mask = (
            (tmp["_chr"] == chrom) &
            (tmp["_pos"].notna()) &
            (tmp["_pos"] >= window_start) &
            (tmp["_pos"] <= window_end)
        )

        local = tmp.loc[mask, ["_chr", "_pos", "_p", "_id", "_gene"]].copy()

        if not local.empty:
            all_local_chunks.append(local)

            if gene is not None:
                gene_local = local[local["_gene"].astype(str) == str(gene)].copy()

                if not gene_local.empty:
                    gene_local_chunks.append(gene_local)

    if len(all_local_chunks) == 0:
        return pd.DataFrame(columns=["_chr", "_pos", "_p", "_id", "_gene", "_logp"]), "none_found"

    all_local_df = add_logp(pd.concat(all_local_chunks, ignore_index=True), p_col="_p")

    if gene_filter == "none" or gene is None:
        return all_local_df, "none"

    if gene_filter == "exact":
        if len(gene_local_chunks) == 0:
            return pd.DataFrame(columns=["_chr", "_pos", "_p", "_id", "_gene", "_logp"]), "exact_no_gene_found"

        gene_df = add_logp(pd.concat(gene_local_chunks, ignore_index=True), p_col="_p")
        return gene_df, "exact"

    if gene_filter == "auto":
        if len(gene_local_chunks) > 0:
            gene_df = add_logp(pd.concat(gene_local_chunks, ignore_index=True), p_col="_p")
            return gene_df, "exact"

        return all_local_df, "auto_fallback_all_local"

    raise ValueError("gene_filter should be one of: exact, auto, none")


# ============================================================
# LD helpers
# ============================================================

def choose_min_p_lead(df, source_name):
    if df is None or df.empty:
        return None

    idx = df["_p"].astype(float).idxmin()
    row = df.loc[idx]

    return {
        "chrom": normalize_chr(row["_chr"]),
        "pos": int(row["_pos"]),
        "rsid": str(row["_id"]),
        "source": source_name,
        "p": float(row["_p"]),
        "id": str(row["_id"])
    }


def read_bim_maps(ld_bfile, chrom=None):
    bim_file = Path(str(ld_bfile) + ".bim")

    if not bim_file.exists():
        raise FileNotFoundError(f"BIM file not found: {bim_file}")

    bim = pd.read_csv(
        bim_file,
        sep=r"\s+",
        header=None,
        names=["CHR", "SNP", "CM", "BP", "A1", "A2"],
        dtype={"CHR": str, "SNP": str}
    )

    bim["CHR"] = bim["CHR"].map(normalize_chr)
    bim["BP"] = pd.to_numeric(bim["BP"], errors="coerce").astype("Int64")
    bim = bim.dropna(subset=["BP"]).copy()
    bim["BP"] = bim["BP"].astype(int)

    if chrom is not None:
        bim = bim[bim["CHR"] == normalize_chr(chrom)].copy()

    id_to_pos = dict(zip(bim["SNP"].astype(str), bim["BP"].astype(int)))

    pos_to_ids = {}
    for _, row in bim.iterrows():
        pos_to_ids.setdefault(int(row["BP"]), []).append(str(row["SNP"]))

    return bim, id_to_pos, pos_to_ids


def find_lead_snp_id(ld_bfile, chrom, lead_pos, lead_rsid=None):
    _, id_to_pos, pos_to_ids = read_bim_maps(ld_bfile, chrom=chrom)

    candidates = []

    if lead_rsid is not None and str(lead_rsid).strip() != "":
        candidates.append(str(lead_rsid).strip())

    chrom = normalize_chr(chrom)

    candidates.append(f"{chrom}_{lead_pos}")
    candidates.append(f"chr{chrom}_{lead_pos}")

    for cand in candidates:
        if cand in id_to_pos:
            return cand

    if int(lead_pos) in pos_to_ids:
        return pos_to_ids[int(lead_pos)][0]

    raise ValueError(
        f"Lead SNP not found in BIM. Tried lead={lead_rsid}, "
        f"{chrom}_{lead_pos}, chr{chrom}_{lead_pos}, and BP={lead_pos}."
    )


def run_plink_ld(
    plink_bin,
    ld_bfile,
    chrom,
    window_start,
    window_end,
    lead_pos,
    lead_rsid,
    out_prefix,
    ld_window_kb=None,
    keep_ld_files=False
):
    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    lead_id = find_lead_snp_id(
        ld_bfile=ld_bfile,
        chrom=chrom,
        lead_pos=lead_pos,
        lead_rsid=lead_rsid
    )

    if ld_window_kb is None:
        ld_window_kb = max(
            1000,
            int(math.ceil((window_end - window_start) / 1000.0)) + 1
        )

    cmd = [
        plink_bin,
        "--bfile", str(ld_bfile),
        "--chr", str(normalize_chr(chrom)),
        "--from-bp", str(int(window_start)),
        "--to-bp", str(int(window_end)),
        "--ld-snp", str(lead_id),
        "--r2",
        "--ld-window", "999999",
        "--ld-window-kb", str(int(ld_window_kb)),
        "--ld-window-r2", "0",
        "--out", str(out_prefix)
    ]

    print("[INFO] Running PLINK LD:")
    print("       " + " ".join(cmd))

    subprocess.run(cmd, check=True)

    ld_file = Path(str(out_prefix) + ".ld")

    if not ld_file.exists():
        raise FileNotFoundError(f"PLINK LD output not found: {ld_file}")

    if not keep_ld_files:
        for suffix in [".log", ".nosex"]:
            f = Path(str(out_prefix) + suffix)

            if f.exists():
                try:
                    f.unlink()
                except Exception:
                    pass

    return ld_file, lead_id


def read_plink_ld_as_maps(ld_file, ld_bfile, chrom, lead_id, lead_pos):
    ld_file = Path(ld_file)

    if not ld_file.exists():
        raise FileNotFoundError(f"LD file not found: {ld_file}")

    ld = pd.read_csv(ld_file, sep=r"\s+", engine="python")

    if "R2" not in ld.columns:
        raise ValueError(f"Column R2 not found in LD file: {ld_file}")

    if "SNP_A" not in ld.columns or "SNP_B" not in ld.columns:
        raise ValueError(f"SNP_A/SNP_B columns not found in LD file: {ld_file}")

    _, id_to_pos, _ = read_bim_maps(ld_bfile, chrom=chrom)

    id_to_r2 = {str(lead_id): 1.0}
    pos_to_r2 = {int(lead_pos): 1.0}

    for _, row in ld.iterrows():
        snp_a = str(row["SNP_A"])
        snp_b = str(row["SNP_B"])
        r2 = pd.to_numeric(row["R2"], errors="coerce")

        if pd.isna(r2):
            continue

        r2 = float(r2)

        if snp_a == str(lead_id):
            other = snp_b
        elif snp_b == str(lead_id):
            other = snp_a
        else:
            continue

        id_to_r2[other] = max(id_to_r2.get(other, 0.0), r2)

        if other in id_to_pos:
            pos = int(id_to_pos[other])
            pos_to_r2[pos] = max(pos_to_r2.get(pos, 0.0), r2)

    return id_to_r2, pos_to_r2


def add_ld_to_df(df, id_to_r2, pos_to_r2):
    if df is None or df.empty:
        return df

    df = df.copy()
    values = []

    for _, row in df.iterrows():
        snp_id = str(row["_id"])
        chrom = normalize_chr(row["_chr"])
        pos = int(row["_pos"])

        r2 = np.nan

        if snp_id in id_to_r2:
            r2 = id_to_r2[snp_id]
        else:
            key1 = f"{chrom}_{pos}"
            key2 = f"chr{chrom}_{pos}"

            if key1 in id_to_r2:
                r2 = id_to_r2[key1]
            elif key2 in id_to_r2:
                r2 = id_to_r2[key2]
            elif pos in pos_to_r2:
                r2 = pos_to_r2[pos]

        values.append(r2)

    df["_ld_r2"] = values
    df["_ld_r2"] = (
        df["_ld_r2"]
        .replace(["NA", "NaN", "nan", "NAN", ".", "", "None", "none", "null"], np.nan)
    )
    df["_ld_r2"] = pd.to_numeric(df["_ld_r2"], errors="coerce")

    return df


# ============================================================
# LD cache
# ============================================================

def build_ld_cache_path(
    outdir,
    out_prefix_name,
    source_name,
    chrom,
    window_start,
    window_end,
    lead_info,
    lead_id,
    ld_bfile
):
    cache_dir = Path(outdir) / "LD_tmp"
    cache_dir.mkdir(parents=True, exist_ok=True)

    key_text = (
        f"out_prefix={out_prefix_name}|"
        f"source={source_name}|"
        f"chr={normalize_chr(chrom)}|"
        f"start={int(window_start)}|"
        f"end={int(window_end)}|"
        f"lead_pos={int(lead_info['pos'])}|"
        f"lead_id={lead_info.get('id', '')}|"
        f"plink_lead_id={lead_id}|"
        f"ld_bfile={str(ld_bfile)}"
    )

    digest = hashlib.md5(key_text.encode("utf-8")).hexdigest()[:12]

    filename = (
        f"{safe_name(source_name)}."
        f"chr{normalize_chr(chrom)}_{int(window_start)}_{int(window_end)}."
        f"lead_{int(lead_info['pos'])}."
        f"{digest}.ld_cache.tsv.gz"
    )

    return cache_dir / filename


def save_ld_cache(
    cache_file,
    source_name,
    chrom,
    window_start,
    window_end,
    lead_info,
    lead_id,
    id_to_r2,
    pos_to_r2,
    ld_bfile
):
    rows = []

    for snp_id, r2 in id_to_r2.items():
        rows.append({
            "map_type": "id",
            "key": str(snp_id),
            "r2": float(r2)
        })

    for pos, r2 in pos_to_r2.items():
        rows.append({
            "map_type": "pos",
            "key": str(int(pos)),
            "r2": float(r2)
        })

    cache_df = pd.DataFrame(rows)

    if cache_df.empty:
        cache_df = pd.DataFrame(columns=["map_type", "key", "r2"])

    cache_df = cache_df.drop_duplicates(
        subset=["map_type", "key"],
        keep="first"
    ).copy()

    cache_df["source_name"] = str(source_name)
    cache_df["chrom"] = normalize_chr(chrom)
    cache_df["window_start"] = int(window_start)
    cache_df["window_end"] = int(window_end)
    cache_df["lead_chrom"] = normalize_chr(lead_info["chrom"])
    cache_df["lead_pos"] = int(lead_info["pos"])
    cache_df["lead_variant_id"] = str(lead_info.get("id", ""))
    cache_df["lead_rsid"] = str(lead_info.get("rsid", ""))
    cache_df["lead_p"] = lead_info.get("p", np.nan)
    cache_df["lead_plink_id"] = str(lead_id)
    cache_df["ld_bfile"] = str(ld_bfile)

    cache_file = Path(cache_file)
    cache_file.parent.mkdir(parents=True, exist_ok=True)

    cache_df.to_csv(
        cache_file,
        sep="\t",
        index=False,
        compression="gzip"
    )

    print(f"[INFO] LD cache saved: {cache_file}")


def read_ld_cache(cache_file):
    cache_file = Path(cache_file)

    if not cache_file.exists():
        raise FileNotFoundError(f"LD cache not found: {cache_file}")

    cache_df = pd.read_csv(
        cache_file,
        sep="\t",
        compression="infer",
        dtype={
            "map_type": str,
            "key": str,
            "source_name": str,
            "chrom": str,
            "lead_chrom": str,
            "lead_variant_id": str,
            "lead_rsid": str,
            "lead_plink_id": str,
            "ld_bfile": str
        }
    )

    if cache_df.empty:
        raise ValueError(f"LD cache is empty: {cache_file}")

    cache_df["r2"] = (
        cache_df["r2"]
        .replace(["NA", "NaN", "nan", "NAN", ".", "", "None", "none", "null"], np.nan)
    )
    cache_df["r2"] = pd.to_numeric(cache_df["r2"], errors="coerce")
    cache_df = cache_df.dropna(subset=["r2"]).copy()

    meta = cache_df.iloc[0]

    id_to_r2 = {}
    pos_to_r2 = {}

    id_df = cache_df[cache_df["map_type"] == "id"].copy()
    pos_df = cache_df[cache_df["map_type"] == "pos"].copy()

    for _, row in id_df.iterrows():
        id_to_r2[str(row["key"])] = float(row["r2"])

    for _, row in pos_df.iterrows():
        try:
            pos_to_r2[int(row["key"])] = float(row["r2"])
        except Exception:
            continue

    lead_info = {
        "chrom": normalize_chr(meta["lead_chrom"]),
        "pos": int(meta["lead_pos"]),
        "rsid": str(meta["lead_rsid"]),
        "source": str(meta["source_name"]),
        "p": float(meta["lead_p"]) if str(meta["lead_p"]) != "nan" else np.nan,
        "id": str(meta["lead_variant_id"]),
        "lead_id": str(meta["lead_plink_id"])
    }

    print(f"[INFO] LD cache loaded: {cache_file}")
    print(
        f"[INFO] Cached LD lead: "
        f"source={lead_info['source']}, "
        f"chr{lead_info['chrom']}:{lead_info['pos']}, "
        f"id={lead_info['id']}, "
        f"PLINK_ID={lead_info['lead_id']}"
    )

    return id_to_r2, pos_to_r2, lead_info


def add_panel_ld(
    df,
    source_name,
    ld_bfile,
    plink_bin,
    chrom,
    window_start,
    window_end,
    outdir,
    out_prefix_name,
    ld_window_kb=None,
    keep_ld_files=False
):
    if df is None or df.empty:
        return df, None

    lead_info = choose_min_p_lead(df, source_name=source_name)

    if lead_info is None:
        return df, None

    if normalize_chr(lead_info["chrom"]) != normalize_chr(chrom):
        raise ValueError(f"{source_name} LD lead chromosome differs from plot chromosome.")

    lead_id = find_lead_snp_id(
        ld_bfile=ld_bfile,
        chrom=chrom,
        lead_pos=lead_info["pos"],
        lead_rsid=lead_info["rsid"]
    )

    lead_info["lead_id"] = lead_id

    cache_file = build_ld_cache_path(
        outdir=outdir,
        out_prefix_name=out_prefix_name,
        source_name=source_name,
        chrom=chrom,
        window_start=window_start,
        window_end=window_end,
        lead_info=lead_info,
        lead_id=lead_id,
        ld_bfile=ld_bfile
    )

    if cache_file.exists():
        id_to_r2, pos_to_r2, lead_info = read_ld_cache(cache_file)

    else:
        ld_outdir = Path(outdir) / "LD_tmp"
        ld_prefix = ld_outdir / safe_name(out_prefix_name + f".{source_name}.LD")

        ld_file, lead_id_from_plink = run_plink_ld(
            plink_bin=plink_bin,
            ld_bfile=ld_bfile,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            lead_pos=lead_info["pos"],
            lead_rsid=lead_info["rsid"],
            out_prefix=ld_prefix,
            ld_window_kb=ld_window_kb,
            keep_ld_files=keep_ld_files
        )

        id_to_r2, pos_to_r2 = read_plink_ld_as_maps(
            ld_file=ld_file,
            ld_bfile=ld_bfile,
            chrom=chrom,
            lead_id=lead_id_from_plink,
            lead_pos=lead_info["pos"]
        )

        lead_info["lead_id"] = lead_id_from_plink

        save_ld_cache(
            cache_file=cache_file,
            source_name=source_name,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            lead_info=lead_info,
            lead_id=lead_id_from_plink,
            id_to_r2=id_to_r2,
            pos_to_r2=pos_to_r2,
            ld_bfile=ld_bfile
        )

        if not keep_ld_files:
            try:
                Path(ld_file).unlink()
            except Exception:
                pass

    print(
        f"[INFO] {source_name} LD lead selected by minimum P: "
        f"chr{lead_info['chrom']}:{lead_info['pos']}, "
        f"id={lead_info['id']}, P={lead_info['p']}, "
        f"PLINK_ID={lead_info.get('lead_id', lead_id)}"
    )

    df = add_ld_to_df(df, id_to_r2=id_to_r2, pos_to_r2=pos_to_r2)

    return df, lead_info


# ============================================================
# Colors and legends
# ============================================================

LD_BOUNDS = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

LD_COLORS = [
    "#2166AC",
    "#67A9CF",
    "#1A9850",
    "#FDAE61",
    "#D73027",
]


def clean_ld_dataframe(df, panel_name="LD"):
    """
      LD r2   SNP.

     :
      NaN,  NA / nan / . /  .
    """
    if df is None or df.empty:
        return df

    if "_ld_r2" not in df.columns:
        print(f"[WARN] {panel_name}: _ld_r2 column not found; no LD-colored points will be drawn.")
        return df.iloc[0:0].copy()

    out = df.copy()

    out["_ld_r2"] = (
        out["_ld_r2"]
        .replace(["NA", "NaN", "nan", "NAN", ".", "", "None", "none", "null"], np.nan)
    )
    out["_ld_r2"] = pd.to_numeric(out["_ld_r2"], errors="coerce")

    n_before = len(out)
    out = out[out["_ld_r2"].notna()].copy()
    n_after = len(out)

    print(f"[INFO] {panel_name}: removed {n_before - n_after} SNPs with NA LD; kept {n_after} SNPs.")

    return out


def ld_color_one(r2):
    if r2 < 0.2:
        return LD_COLORS[0]

    if r2 < 0.4:
        return LD_COLORS[1]

    if r2 < 0.6:
        return LD_COLORS[2]

    if r2 < 0.8:
        return LD_COLORS[3]

    return LD_COLORS[4]


def assign_ld_colors(r2_series):
    return [ld_color_one(float(r2)) for r2 in r2_series]


def draw_ld_scatter(
    ax,
    plot_df,
    point_size=8,
    alpha=0.75,
    ld_na_color="#C8C8C8"
):
    """
      LD  .

     :
      NA  .
    """
    if plot_df is None or plot_df.empty:
        return

    plot_df = clean_ld_dataframe(plot_df, panel_name="LD scatter")

    if plot_df.empty:
        return

    colors = assign_ld_colors(plot_df["_ld_r2"])

    ax.scatter(
        plot_df["_x_mb"],
        plot_df["_logp"],
        s=point_size,
        c=colors,
        alpha=alpha,
        edgecolors="none",
        rasterized=False,
        zorder=3
    )


def add_ld_combined_legend(
    ax,
    ld_na_color="#C8C8C8",
    fontsize=8
):
    """
    LD  .

    NA  ,  NA  .
    """
    leg_ax = ax.inset_axes([1.018, 0.20, 0.075, 0.66])

    for i, color in enumerate(LD_COLORS):
        y0 = LD_BOUNDS[i]
        y1 = LD_BOUNDS[i + 1]

        leg_ax.add_patch(
            Rectangle(
                (0.00, y0),
                0.42,
                y1 - y0,
                facecolor=color,
                edgecolor="black",
                linewidth=0.35
            )
        )

    for y in LD_BOUNDS:
        leg_ax.text(
            0.52,
            y,
            f"{y:.1f}",
            ha="left",
            va="center",
            fontsize=fontsize
        )

    leg_ax.text(
        0.21,
        1.08,
        r"$r^2$",
        ha="center",
        va="bottom",
        fontsize=fontsize + 1
    )

    leg_ax.set_xlim(0, 1.2)
    leg_ax.set_ylim(-0.02, 1.14)
    leg_ax.axis("off")


# ============================================================
# Annotation
# ============================================================

def annotate_multiple_sites(
    ax,
    plot_df,
    mark_df,
    chrom,
    ymax,
    label_fontsize=8,
    label_offset_step=9,
    label_sites=True,
    plot_ld=False,
    point_size=8,
    mark_point_size=None,
    no_ld_mark_color="#D73027",
    ld_na_color="#C8C8C8",
    mark_edge_color="black",
    mark_edge_width=1.1
):
    chrom = normalize_chr(chrom)

    if mark_point_size is None:
        mark_point_size = point_size + 2

    mark_df = mark_df.copy()
    mark_df["chrom"] = mark_df["chrom"].map(normalize_chr)
    mark_df = mark_df[mark_df["chrom"] == chrom].copy()

    if mark_df.empty:
        return

    mark_df = mark_df.sort_values("pos").reset_index(drop=True)

    for i, row in mark_df.iterrows():
        pos = int(row["pos"])
        label = str(row["label"])
        x_mb = pos / 1e6

        if plot_df is None or plot_df.empty:
            hit_df = pd.DataFrame()
        else:
            hit_df = plot_df[plot_df["_pos"].astype(int) == pos].copy()

        offset_y = 10 + (i % 5) * label_offset_step
        offset_x = 0

        if not hit_df.empty:
            best = hit_df.sort_values("_logp", ascending=False).iloc[0]

            if plot_ld:
                r2 = best.get("_ld_r2", np.nan)

                if pd.isna(r2):
                    continue

                face_color = ld_color_one(float(r2))

                ax.scatter(
                    [best["_x_mb"]],
                    [best["_logp"]],
                    s=mark_point_size,
                    marker="o",
                    facecolors=face_color,
                    edgecolors=mark_edge_color,
                    linewidths=mark_edge_width,
                    alpha=1.0,
                    rasterized=False,
                    zorder=25,
                    clip_on=False
                )

            else:
                ax.scatter(
                    [best["_x_mb"]],
                    [best["_logp"]],
                    s=mark_point_size,
                    marker="o",
                    facecolors=no_ld_mark_color,
                    edgecolors=mark_edge_color,
                    linewidths=0.6,
                    alpha=1.0,
                    rasterized=False,
                    zorder=25,
                    clip_on=False
                )

            if label_sites:
                ax.annotate(
                    label,
                    xy=(best["_x_mb"], best["_logp"]),
                    xytext=(offset_x, offset_y),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=label_fontsize,
                    color="black",
                    arrowprops=dict(
                        arrowstyle="-",
                        lw=0.7,
                        color="black",
                        shrinkA=0,
                        shrinkB=3
                    ),
                    zorder=30
                )

        else:
            if label_sites and not plot_ld:
                y_for_label = ymax * 0.80

                ax.text(
                    x_mb,
                    y_for_label,
                    label,
                    ha="center",
                    va="bottom",
                    fontsize=label_fontsize,
                    color="black",
                    zorder=30
                )


def draw_ld_lead_marker(
    ax,
    plot_df,
    lead_info,
    chrom,
    lead_marker_size,
    lead_marker_color="#7B3294"
):
    if lead_info is None:
        return

    if normalize_chr(lead_info["chrom"]) != normalize_chr(chrom):
        return

    if plot_df is None or plot_df.empty:
        return

    lead_pos = int(lead_info["pos"])
    hit_df = plot_df[plot_df["_pos"].astype(int) == lead_pos].copy()

    if hit_df.empty:
        return

    best = hit_df.sort_values("_logp", ascending=False).iloc[0]

    ax.scatter(
        [best["_x_mb"]],
        [best["_logp"]],
        s=lead_marker_size,
        marker="D",
        c=lead_marker_color,
        edgecolors="black",
        linewidths=0.55,
        zorder=28,
        clip_on=False
    )


# ============================================================
# Plot panels
# ============================================================

def draw_one_panel(
    ax,
    df,
    chrom,
    window_start,
    window_end,
    enhancer_info=None,
    threshold_logp=5.0,
    ylabel=r"$-\log_{10}(P)$",
    point_size=10,
    sig_point_size=None,
    mark_point_size=None,
    mark_df=None,
    alpha=0.75,
    label_fontsize=8,
    nonsig_color="#2F2A2E",
    sig_color="#93BFC8",
    enhancer_shade_color="#D9D9D9",
    enhancer_shade_alpha=0.30,
    show_threshold_line=True,
    threshold_line_color="#808080",
    threshold_line_width=0.9,
    empty_text="No data",
    plot_ld=False,
    show_ld_legend=True,
    ld_na_color="#C8C8C8",
    lead_info=None,
    lead_marker_size=None,
    lead_marker_color="#7B3294",
    no_ld_mark_color="#D73027",
    mark_edge_width=1.1
):
    if sig_point_size is None:
        sig_point_size = point_size

    if lead_marker_size is None:
        lead_marker_size = point_size + 1

    if mark_point_size is None:
        mark_point_size = point_size + 2

    ax.set_xlim(window_start / 1e6, window_end / 1e6)

    if enhancer_info is not None:
        if normalize_chr(enhancer_info["chrom"]) == normalize_chr(chrom):
            ax.axvspan(
                enhancer_info["start"] / 1e6,
                enhancer_info["end"] / 1e6,
                color=enhancer_shade_color,
                alpha=enhancer_shade_alpha,
                linewidth=0,
                zorder=0
            )

    if df is None or df.empty:
        ymax = max(threshold_logp + 1.0, 1.0)
        ax.set_ylim(0, ymax)

        if show_threshold_line:
            ax.axhline(
                threshold_logp,
                color=threshold_line_color,
                linestyle="-",
                linewidth=threshold_line_width,
                alpha=0.85,
                zorder=1
            )

        ax.text(
            0.5,
            0.55,
            empty_text,
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=10
        )

        ax.set_ylabel(ylabel, fontsize=11)
        return

    plot_df = df.copy()
    plot_df["_x_mb"] = plot_df["_pos"].astype(float) / 1e6

    if plot_ld:
        plot_df = clean_ld_dataframe(
            plot_df,
            panel_name=ylabel.replace("\n", " ")
        )

        if plot_df.empty:
            ymax = max(threshold_logp + 1.0, 1.0)
            ax.set_ylim(0, ymax)

            if show_threshold_line:
                ax.axhline(
                    threshold_logp,
                    color=threshold_line_color,
                    linestyle="-",
                    linewidth=threshold_line_width,
                    alpha=0.85,
                    zorder=1
                )

            ax.text(
                0.5,
                0.55,
                "No SNPs with LD",
                transform=ax.transAxes,
                ha="center",
                va="center",
                fontsize=10
            )

            ax.set_ylabel(ylabel, fontsize=11)
            return

        draw_ld_scatter(
            ax=ax,
            plot_df=plot_df,
            point_size=point_size,
            alpha=alpha,
            ld_na_color=ld_na_color
        )

        if show_ld_legend:
            add_ld_combined_legend(
                ax=ax,
                ld_na_color=ld_na_color,
                fontsize=8
            )

    else:
        sig_df = plot_df[plot_df["_logp"] >= threshold_logp].copy()
        nonsig_df = plot_df[plot_df["_logp"] < threshold_logp].copy()

        if not nonsig_df.empty:
            ax.scatter(
                nonsig_df["_x_mb"],
                nonsig_df["_logp"],
                s=point_size,
                c=nonsig_color,
                alpha=alpha,
                edgecolors="none",
                rasterized=False,
                zorder=2
            )

        if not sig_df.empty:
            ax.scatter(
                sig_df["_x_mb"],
                sig_df["_logp"],
                s=sig_point_size,
                c=sig_color,
                alpha=alpha,
                edgecolors="none",
                rasterized=False,
                zorder=3
            )

    if show_threshold_line:
        ax.axhline(
            threshold_logp,
            color=threshold_line_color,
            linestyle="-",
            linewidth=threshold_line_width,
            alpha=0.85,
            zorder=1
        )

    ymax_data = float(plot_df["_logp"].max())
    ymax = math.ceil(max(threshold_logp + 1.5, ymax_data + 2.5))
    ax.set_ylim(0, ymax)

    if mark_df is not None:
        annotate_multiple_sites(
            ax=ax,
            plot_df=plot_df,
            mark_df=mark_df,
            chrom=chrom,
            ymax=ymax,
            label_fontsize=label_fontsize,
            label_sites=True,
            plot_ld=plot_ld,
            point_size=point_size,
            mark_point_size=mark_point_size,
            no_ld_mark_color=no_ld_mark_color,
            ld_na_color=ld_na_color,
            mark_edge_width=mark_edge_width
        )

    if plot_ld:
        draw_ld_lead_marker(
            ax=ax,
            plot_df=plot_df,
            lead_info=lead_info,
            chrom=chrom,
            lead_marker_size=lead_marker_size,
            lead_marker_color=lead_marker_color
        )

    ax.set_ylabel(ylabel, fontsize=11)


def plot_gwas_eqtl_local(
    gwas_df,
    eqtl_df,
    chrom,
    window_start,
    window_end,
    out_prefix,
    trait="GWAS",
    tissue="Tissue",
    gene=None,
    mark_df=None,
    enhancer_info=None,
    plot_mode="both",
    gwas_threshold_logp=5.0,
    eqtl_threshold_logp=5.0,
    width=6.0,
    height=4.8,
    point_size=8,
    sig_point_size=12,
    mark_point_size=None,
    alpha=0.75,
    label_fontsize=8,
    enhancer_shade_color="#D9D9D9",
    enhancer_shade_alpha=0.30,
    show_threshold_line=True,
    threshold_line_color="#808080",
    threshold_line_width=0.9,
    x_decimals=2,
    x_nbins=6,
    plot_ld=False,
    show_ld_legend=True,
    ld_na_color="#C8C8C8",
    gwas_lead_info=None,
    eqtl_lead_info=None,
    lead_marker_size=None,
    lead_marker_color="#7B3294",
    no_ld_mark_color="#D73027",
    mark_edge_width=1.1
):
    if lead_marker_size is None:
        lead_marker_size = point_size + 1

    if mark_point_size is None:
        mark_point_size = point_size + 2

    panels = []

    if plot_mode in ["both", "gwas"]:
        panels.append({
            "name": "GWAS",
            "df": gwas_df,
            "threshold": gwas_threshold_logp,
            "ylabel": f"{trait}\n" + r"$-\log_{10}(P)$",
            "empty_text": "No GWAS data",
            "lead_info": gwas_lead_info
        })

    if plot_mode in ["both", "eqtl"]:
        panels.append({
            "name": "eQTL",
            "df": eqtl_df,
            "threshold": eqtl_threshold_logp,
            "ylabel": f"{tissue} eQTL\n" + r"$-\log_{10}(P)$",
            "empty_text": "No eQTL data",
            "lead_info": eqtl_lead_info
        })

    if len(panels) == 0:
        raise ValueError("plot_mode should be one of: both, gwas, eqtl")

    fig, axes = plt.subplots(
        nrows=len(panels),
        ncols=1,
        figsize=(width, height),
        sharex=True,
        gridspec_kw={"hspace": 0.08}
    )

    if len(panels) == 1:
        axes = [axes]

    for i, (ax, panel) in enumerate(zip(axes, panels)):
        draw_one_panel(
            ax=ax,
            df=panel["df"],
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            enhancer_info=enhancer_info,
            threshold_logp=panel["threshold"],
            ylabel=panel["ylabel"],
            point_size=point_size,
            sig_point_size=sig_point_size,
            mark_point_size=mark_point_size,
            mark_df=mark_df,
            alpha=alpha,
            label_fontsize=label_fontsize,
            enhancer_shade_color=enhancer_shade_color,
            enhancer_shade_alpha=enhancer_shade_alpha,
            show_threshold_line=show_threshold_line,
            threshold_line_color=threshold_line_color,
            threshold_line_width=threshold_line_width,
            empty_text=panel["empty_text"],
            plot_ld=plot_ld,
            show_ld_legend=(show_ld_legend and i == len(panels) - 1),
            ld_na_color=ld_na_color,
            lead_info=panel["lead_info"],
            lead_marker_size=lead_marker_size,
            lead_marker_color=lead_marker_color,
            no_ld_mark_color=no_ld_mark_color,
            mark_edge_width=mark_edge_width
        )

        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(axis="both", labelsize=10)

    axes[-1].xaxis.set_major_locator(MaxNLocator(nbins=x_nbins))
    axes[-1].xaxis.set_major_formatter(
        FuncFormatter(lambda x, pos: f"{x:.{x_decimals}f}")
    )
    axes[-1].set_xlabel(f"chr{chrom} position (Mb)", fontsize=11)

    gene_text = "" if gene is None else f" | {gene}"
    n_mark = 0 if mark_df is None else len(mark_df)
    ld_text = " | LD" if plot_ld else ""

    title = f"{trait} | {tissue}{gene_text} | {n_mark} marked SNPs{ld_text}"
    fig.suptitle(title, fontsize=12, y=0.98)

    if plot_ld and show_ld_legend:
        fig.tight_layout(rect=[0, 0, 0.94, 0.96])
    else:
        fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_prefix = Path(out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    out_file = f"{out_prefix}.pdf"
    fig.savefig(out_file, bbox_inches="tight")
    plt.close(fig)

    print(f"[INFO] PDF saved: {out_file}")


# ============================================================
# Main
# ============================================================

def run_one(
    gwas_file,
    eqtl_dir,
    tissue,
    trait,
    gene,
    sites=None,
    rsids=None,
    rsid_gvf=None,
    enhancer=None,
    shade_region=None,
    region=None,
    outdir="gwas_eqtl_local_manhattan",
    flank_mult=5.0,
    window_bp=500000,
    plot_mode="both",
    gwas_chr_col="auto",
    gwas_pos_col="auto",
    gwas_p_col="auto",
    gwas_snp_col="auto",
    eqtl_pattern="{tissue}.cis_qtl_pairs.chr{chrom}.txt.gz",
    eqtl_gene_col="phenotype_id",
    eqtl_variant_col="variant_id",
    eqtl_p_col="pval_nominal",
    gene_filter="exact",
    gwas_threshold_logp=5.0,
    eqtl_threshold_logp=5.0,
    width=6.0,
    height=4.8,
    point_size=8,
    sig_point_size=12,
    mark_point_size=None,
    alpha=0.75,
    label_fontsize=8,
    enhancer_shade_color="#D9D9D9",
    enhancer_shade_alpha=0.30,
    show_threshold_line=True,
    threshold_line_color="#808080",
    threshold_line_width=0.9,
    x_decimals=2,
    x_nbins=6,
    chunksize=1000000,
    plot_ld=False,
    ld_bfile=None,
    plink_bin="plink",
    ld_window_kb=None,
    keep_ld_files=False,
    show_ld_legend=True,
    ld_na_color="#C8C8C8",
    lead_marker_size=None,
    lead_marker_color="#7B3294",
    no_ld_mark_color="#D73027",
    mark_edge_width=1.1
):
    chrom, window_start, window_end, enhancer_info = determine_window(
        sites=sites,
        enhancer=enhancer,
        region=region,
        flank_mult=flank_mult,
        window_bp=window_bp
    )

    if shade_region is not None:
        sh_chrom, sh_start, sh_end = parse_region(shade_region)

        if normalize_chr(sh_chrom) == normalize_chr(chrom):
            enhancer_info = {
                "chrom": sh_chrom,
                "start": sh_start,
                "end": sh_end
            }
        else:
            print(f"[WARN] shade-region chromosome {sh_chrom} differs from plot chromosome {chrom}; ignored.")

    if sites is not None:
        mark_df = parse_mark_sites(sites=sites, rsids=rsids)

        if rsids is None and rsid_gvf is not None:
            mark_df = update_mark_labels_from_gvf(mark_df, rsid_gvf)
    else:
        mark_df = None

    if mark_df is not None:
        mark_chr_set = set(mark_df["chrom"].map(normalize_chr).tolist())

        if len(mark_chr_set) != 1:
            raise ValueError(f"Marked sites contain multiple chromosomes: {mark_chr_set}")

        if normalize_chr(list(mark_chr_set)[0]) != normalize_chr(chrom):
            raise ValueError(
                f"Marked site chromosome and plotting chromosome are different: "
                f"{list(mark_chr_set)[0]} vs {chrom}"
            )

    if plot_mode == "eqtl" and plot_ld:
        print("[INFO] plot_mode is eqtl only; LD calculation and LD coloring are disabled automatically.")
        plot_ld = False

    if lead_marker_size is None:
        lead_marker_size = point_size + 1

    if mark_point_size is None:
        mark_point_size = point_size + 2

    print(f"[INFO] Plot region: chr{chrom}:{window_start}-{window_end}")
    print(
        f"[INFO] point_size={point_size}, "
        f"mark_point_size={mark_point_size}, "
        f"lead_marker_size={lead_marker_size}, "
        f"mark_edge_width={mark_edge_width}"
    )

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if shade_region is not None:
        region_name = safe_name(shade_region)
    elif enhancer is not None:
        region_name = safe_name(enhancer)
    elif region is not None:
        region_name = safe_name(region)
    else:
        region_name = f"chr{chrom}_{window_start}_{window_end}"

    if mark_df is not None:
        if len(mark_df) == 1:
            site_part = f"{mark_df.iloc[0]['label']}.chr{chrom}_{mark_df.iloc[0]['pos']}"
        else:
            site_part = f"{len(mark_df)}SNPs"
    else:
        site_part = "noMarkedSNP"

    gene_part = gene if gene is not None else "allGenes"
    ld_part = "LD" if plot_ld else "noLD"

    out_prefix_name = (
        f"{safe_name(trait)}."
        f"{safe_name(tissue)}."
        f"{safe_name(gene_part)}."
        f"{safe_name(region_name)}."
        f"{safe_name(site_part)}."
        f"{ld_part}."
        f"GWAS_eQTL_localManhattan"
    )

    out_prefix = outdir / out_prefix_name

    gwas_df = pd.DataFrame()
    eqtl_df = pd.DataFrame()

    if plot_mode in ["both", "gwas"]:
        if gwas_file is None:
            raise ValueError("--gwas-file is required when plot-mode includes GWAS.")

        gwas_df = read_gwas_local(
            gwas_file=gwas_file,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            gwas_chr_col=gwas_chr_col,
            gwas_pos_col=gwas_pos_col,
            gwas_p_col=gwas_p_col,
            gwas_snp_col=gwas_snp_col,
            chunksize=chunksize
        )

        print(f"[INFO] GWAS points: {len(gwas_df)}")

    if plot_mode in ["both", "eqtl"]:
        eqtl_df, used_gene_filter = read_eqtl_local(
            eqtl_dir=eqtl_dir,
            tissue=tissue,
            chrom=chrom,
            window_start=window_start,
            window_end=window_end,
            gene=gene,
            gene_filter=gene_filter,
            eqtl_pattern=eqtl_pattern,
            eqtl_gene_col=eqtl_gene_col,
            eqtl_variant_col=eqtl_variant_col,
            eqtl_p_col=eqtl_p_col,
            chunksize=chunksize
        )

        print(f"[INFO] eQTL points: {len(eqtl_df)}")
        print(f"[INFO] eQTL gene_filter: {used_gene_filter}")

    gwas_lead_info = None
    eqtl_lead_info = None

    if plot_ld:
        if ld_bfile is None:
            raise ValueError("--ld-bfile is required when --plot-ld is used.")

        if plot_mode in ["both", "gwas"] and not gwas_df.empty:
            gwas_df, gwas_lead_info = add_panel_ld(
                df=gwas_df,
                source_name="GWAS",
                ld_bfile=ld_bfile,
                plink_bin=plink_bin,
                chrom=chrom,
                window_start=window_start,
                window_end=window_end,
                outdir=outdir,
                out_prefix_name=out_prefix_name,
                ld_window_kb=ld_window_kb,
                keep_ld_files=keep_ld_files
            )

        if plot_mode == "both" and not eqtl_df.empty:
            eqtl_df, eqtl_lead_info = add_panel_ld(
                df=eqtl_df,
                source_name="eQTL",
                ld_bfile=ld_bfile,
                plink_bin=plink_bin,
                chrom=chrom,
                window_start=window_start,
                window_end=window_end,
                outdir=outdir,
                out_prefix_name=out_prefix_name,
                ld_window_kb=ld_window_kb,
                keep_ld_files=keep_ld_files
            )

    plot_gwas_eqtl_local(
        gwas_df=gwas_df,
        eqtl_df=eqtl_df,
        chrom=chrom,
        window_start=window_start,
        window_end=window_end,
        out_prefix=out_prefix,
        trait=trait,
        tissue=tissue,
        gene=gene,
        mark_df=mark_df,
        enhancer_info=enhancer_info,
        plot_mode=plot_mode,
        gwas_threshold_logp=gwas_threshold_logp,
        eqtl_threshold_logp=eqtl_threshold_logp,
        width=width,
        height=height,
        point_size=point_size,
        sig_point_size=sig_point_size,
        mark_point_size=mark_point_size,
        alpha=alpha,
        label_fontsize=label_fontsize,
        enhancer_shade_color=enhancer_shade_color,
        enhancer_shade_alpha=enhancer_shade_alpha,
        show_threshold_line=show_threshold_line,
        threshold_line_color=threshold_line_color,
        threshold_line_width=threshold_line_width,
        x_decimals=x_decimals,
        x_nbins=x_nbins,
        plot_ld=plot_ld,
        show_ld_legend=show_ld_legend,
        ld_na_color=ld_na_color,
        gwas_lead_info=gwas_lead_info,
        eqtl_lead_info=eqtl_lead_info,
        lead_marker_size=lead_marker_size,
        lead_marker_color=lead_marker_color,
        no_ld_mark_color=no_ld_mark_color,
        mark_edge_width=mark_edge_width
    )


def main():
    parser = argparse.ArgumentParser(
        description="Draw local GWAS/eQTL Manhattan plots with optional per-panel LD coloring."
    )

    parser.add_argument("--gwas-file", default=None)
    parser.add_argument("--eqtl-dir", default=None)

    parser.add_argument("--tissue", default="NA")
    parser.add_argument("--trait", default="GWAS")
    parser.add_argument("--gene", default=None)

    parser.add_argument("--site", default=None)
    parser.add_argument("--sites", default=None)
    parser.add_argument("--rsid", default=None)
    parser.add_argument("--rsids", default=None)
    parser.add_argument("--rsid-gvf", default=None)

    parser.add_argument("--enhancer", default=None)
    parser.add_argument("--shade-region", default=None)
    parser.add_argument("--region", default=None)
    parser.add_argument("--outdir", default="gwas_eqtl_local_manhattan")

    parser.add_argument("--flank-mult", type=float, default=5.0)
    parser.add_argument("--window-bp", type=int, default=500000)

    parser.add_argument(
        "--plot-mode",
        default="both",
        choices=["both", "gwas", "eqtl"]
    )

    parser.add_argument("--gwas-chr-col", default="auto")
    parser.add_argument("--gwas-pos-col", default="auto")
    parser.add_argument("--gwas-p-col", default="auto")
    parser.add_argument("--gwas-snp-col", default="auto")

    parser.add_argument(
        "--eqtl-pattern",
        default="{tissue}.cis_qtl_pairs.chr{chrom}.txt.gz"
    )
    parser.add_argument("--eqtl-gene-col", default="phenotype_id")
    parser.add_argument("--eqtl-variant-col", default="variant_id")
    parser.add_argument("--eqtl-p-col", default="pval_nominal")

    parser.add_argument(
        "--gene-filter",
        default="exact",
        choices=["exact", "auto", "none"]
    )

    parser.add_argument("--gwas-threshold-logp", type=float, default=5.0)
    parser.add_argument("--eqtl-threshold-logp", type=float, default=5.0)
    parser.add_argument("--hide-threshold-line", action="store_true")
    parser.add_argument("--threshold-line-color", default="#808080")
    parser.add_argument("--threshold-line-width", type=float, default=0.9)

    parser.add_argument("--width", type=float, default=6.0)
    parser.add_argument("--height", type=float, default=4.8)

    parser.add_argument("--point-size", type=float, default=8)
    parser.add_argument("--sig-point-size", type=float, default=12)
    parser.add_argument("--mark-point-size", type=float, default=None)
    parser.add_argument("--mark-edge-width", type=float, default=1.1)
    parser.add_argument("--no-ld-mark-color", default="#D73027")

    parser.add_argument("--lead-marker-size", type=float, default=None)
    parser.add_argument("--lead-marker-color", default="#7B3294")

    parser.add_argument("--alpha", type=float, default=0.75)
    parser.add_argument("--label-fontsize", type=float, default=8)

    parser.add_argument("--enhancer-shade-color", default="#D9D9D9")
    parser.add_argument("--enhancer-shade-alpha", type=float, default=0.30)

    parser.add_argument("--x-decimals", type=int, default=2)
    parser.add_argument("--x-nbins", type=int, default=6)
    parser.add_argument("--chunksize", type=int, default=1000000)

    parser.add_argument("--formats", default="pdf", help="Deprecated compatibility option; output is always an editable PDF.")

    parser.add_argument("--plot-ld", action="store_true")
    parser.add_argument("--ld-bfile", default=None)
    parser.add_argument("--plink", default="plink")
    parser.add_argument("--ld-window-kb", type=int, default=None)
    parser.add_argument("--keep-ld-files", action="store_true")
    parser.add_argument("--hide-ld-legend", action="store_true")
    parser.add_argument("--ld-na-color", default="#C8C8C8")

    # Compatibility with older commands; ignored in current per-panel LD mode.
    parser.add_argument("--ld-lead-source", default="auto")
    parser.add_argument("--ld-lead-site", default=None)
    parser.add_argument("--ld-lead-rsid", default=None)

    args = parser.parse_args()

    sites = args.sites if args.sites is not None else args.site
    rsids = args.rsids if args.rsids is not None else args.rsid

    run_one(
        gwas_file=args.gwas_file,
        eqtl_dir=args.eqtl_dir,
        tissue=args.tissue,
        trait=args.trait,
        gene=args.gene,
        sites=sites,
        rsids=rsids,
        rsid_gvf=args.rsid_gvf,
        enhancer=args.enhancer,
        shade_region=args.shade_region,
        region=args.region,
        outdir=args.outdir,
        flank_mult=args.flank_mult,
        window_bp=args.window_bp,
        plot_mode=args.plot_mode,
        gwas_chr_col=args.gwas_chr_col,
        gwas_pos_col=args.gwas_pos_col,
        gwas_p_col=args.gwas_p_col,
        gwas_snp_col=args.gwas_snp_col,
        eqtl_pattern=args.eqtl_pattern,
        eqtl_gene_col=args.eqtl_gene_col,
        eqtl_variant_col=args.eqtl_variant_col,
        eqtl_p_col=args.eqtl_p_col,
        gene_filter=args.gene_filter,
        gwas_threshold_logp=args.gwas_threshold_logp,
        eqtl_threshold_logp=args.eqtl_threshold_logp,
        width=args.width,
        height=args.height,
        point_size=args.point_size,
        sig_point_size=args.sig_point_size,
        mark_point_size=args.mark_point_size,
        alpha=args.alpha,
        label_fontsize=args.label_fontsize,
        enhancer_shade_color=args.enhancer_shade_color,
        enhancer_shade_alpha=args.enhancer_shade_alpha,
        show_threshold_line=not args.hide_threshold_line,
        threshold_line_color=args.threshold_line_color,
        threshold_line_width=args.threshold_line_width,
        x_decimals=args.x_decimals,
        x_nbins=args.x_nbins,
        chunksize=args.chunksize,
        plot_ld=args.plot_ld,
        ld_bfile=args.ld_bfile,
        plink_bin=args.plink,
        ld_window_kb=args.ld_window_kb,
        keep_ld_files=args.keep_ld_files,
        show_ld_legend=not args.hide_ld_legend,
        ld_na_color=args.ld_na_color,
        lead_marker_size=args.lead_marker_size,
        lead_marker_color=args.lead_marker_color,
        no_ld_mark_color=args.no_ld_mark_color,
        mark_edge_width=args.mark_edge_width
    )


if __name__ == "__main__":
    main()
