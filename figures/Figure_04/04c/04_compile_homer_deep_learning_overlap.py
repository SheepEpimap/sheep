#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
02_compile_homer_dl_without_overlap_tsv.py


-   overlap.tsv
-   homer   + deeplearning  ,  overlap(  promoter  )
-  output results(  8  file):
  1) keepHash:DL motif   #1/#2( )
  2) dropHash:DL motif   #1/#2(  compendium   category=unresolved   motif   #, )
-  output:
    - raw  ( / )
    - dedup  (  Tissue+Motif+Enhancer+Promoter+Gene+Method  )
    - edge file( ):id type tissue method
    - node file( ):node1 node2 weight tissue method
   : :edge= ,node= , output.

overlap  (  overlap.tsv)
- overlap_key = (Tissue, motif_key_for_overlap, Enhancer, Gene)
- motif_key_for_overlap:
    - unresolved motif: (  unresolved#1),  split('#')
    -   motif:  base(  #  ),  KLF5#1 -> KLF5
- overlap  :  overlap_key   homer   dl
-  :
    -   overlap_key:
          homer  ,Method   "overlap"
        dl  (  dl motif   unresolved)
    -   overlap:
        homer   Method="homer"
        dl     Method="deeplearning"

dedup  (  chain key  )
- dedup_key = (Tissue, Motif_out, Enhancer, Promoter, Gene, Method)
-  :
    Motif_Enhancer:sum
    Enhancer_Gene:max
    Promoter_Gene:max(  1)

 input( )
- homer:      /vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/AARegulatory_module/AA_TSR_E5/motif/fimo_vertebrates_bg0.001/natwork/node/table/all_tables.with_header.tsv
- dl:         /data/home/sczd644/run/zsw_chrombpnet/network/TSR/ALL_fullChain/ALL.all_tables.with_header.fullChain.txt
- compendium: /data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/output/TABLE_motif_compendium.csv
- outdir:     /data/home/sczd644/run/zsw_chrombpnet/network/TSR/complied
"""

import argparse
import csv
import os
import shutil
import subprocess
import shlex
from collections import OrderedDict
from typing import Dict, Tuple, List, Optional


# -----------------------------
# -----------------------------
def split_fields(line: str) -> List[str]:
    """
     :
    -   \t
    -
    """
    line = line.rstrip("\n")
    if "\t" in line:
        return line.split("\t")
    return line.split()


def norm_col(s: str) -> str:
    """ :  +  / / """
    return s.strip().lower().replace(" ", "").replace("-", "").replace("_", "")


def require_idx(idx_map: Dict[str, int], colname: str) -> int:
    k = norm_col(colname)
    if k not in idx_map:
        raise RuntimeError(f"[ERROR] Missing column: {colname}; header={list(idx_map.keys())}")
    return idx_map[k]


def motif_base(m: str) -> str:
    """  motif   base:  #  """
    m = (m or "").strip()
    if "#" in m:
        return m.split("#", 1)[0]
    return m


def which(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def sort_u(infile: str, outfile: str):
    """ :  sort -u;  python set fallback"""
    if which("sort"):
        cmd = f"LC_ALL=C sort -u {shlex.quote(infile)} > {shlex.quote(outfile)}"
        r = subprocess.run(["bash", "-lc", cmd], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"[ERROR] sort -u failed:\n{r.stderr}")
    else:
        seen = set()
        with open(outfile, "w", encoding="utf-8") as w:
            with open(infile, "r", encoding="utf-8") as r:
                for line in r:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line in seen:
                        continue
                    seen.add(line)
                    w.write(line + "\n")


def to_float(x: str) -> float:
    try:
        return float(x)
    except Exception:
        return float("nan")


def fmt_num(x: float) -> str:
    if x != x:  # NaN
        return "NA"
    return repr(x)


def load_unresolved_set(compendium_csv: str) -> set:
    """
      compendium  read category==unresolved   motif_name
      motif:
      - overlap   split('#'),  unresolved#1   unresolved
      - dropHash output  '#'
    """
    unresolved = set()
    with open(compendium_csv, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        need = {"motif_name", "category"}
        got = set(reader.fieldnames or [])
        if not need.issubset(got):
            raise RuntimeError(f"[ERROR] compendium  :{need - got}; got={reader.fieldnames}")

        for row in reader:
            cat = (row.get("category") or "").strip().lower()
            if cat != "unresolved":
                continue
            m = (row.get("motif_name") or "").strip()
            if m:
                unresolved.add(m)
    return unresolved


# -----------------------------
# main workflow
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--homer",
        default="/vol2/mengzhu/SheepFANNG/04_ChromHMM_noblacklist_modif/Merge_chromatin_state/state_variability/"
                "AARegulatory_module/AA_TSR_E5/motif/fimo_vertebrates_bg0.001/natwork/node/table/all_tables.with_header.tsv"
    )
    ap.add_argument(
        "--dl",
        default="/data/home/sczd644/run/zsw_chrombpnet/network/TSR/ALL_fullChain/ALL.all_tables.with_header.fullChain.txt"
    )
    ap.add_argument(
        "--compendium",
        default="/data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/output/TABLE_motif_compendium.csv"
    )
    ap.add_argument(
        "--outdir",
        default="/data/home/sczd644/run/zsw_chrombpnet/network/TSR/complied"
    )
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    unresolved_set = load_unresolved_set(args.compendium)
    print(f"[INFO] unresolved motifs from compendium: {len(unresolved_set)}")

    out_keep_raw = os.path.join(args.outdir, "complied_all_tables.with_header.keepHash.raw.tsv")
    out_keep_dedup = os.path.join(args.outdir, "complied_all_tables.with_header.keepHash.dedup.tsv")
    out_keep_nodes = os.path.join(args.outdir, "complied.edge.fullChain.keepHash.txt")  #
    out_keep_edges = os.path.join(args.outdir, "complied.node.fullChain.keepHash.txt")  #

    out_drop_raw = os.path.join(args.outdir, "complied_all_tables.with_header.dropHash.raw.tsv")
    out_drop_dedup = os.path.join(args.outdir, "complied_all_tables.with_header.dropHash.dedup.tsv")
    out_drop_nodes = os.path.join(args.outdir, "complied.edge.fullChain.dropHash.txt")
    out_drop_edges = os.path.join(args.outdir, "complied.node.fullChain.dropHash.txt")

    tmp_keep_nodes_raw = os.path.join(args.outdir, "._tmp.keep.nodes.raw")
    tmp_keep_edges_raw = os.path.join(args.outdir, "._tmp.keep.edges.raw")
    tmp_keep_nodes_sorted = os.path.join(args.outdir, "._tmp.keep.nodes.sorted")
    tmp_keep_edges_sorted = os.path.join(args.outdir, "._tmp.keep.edges.sorted")

    tmp_drop_nodes_raw = os.path.join(args.outdir, "._tmp.drop.nodes.raw")
    tmp_drop_edges_raw = os.path.join(args.outdir, "._tmp.drop.edges.raw")
    tmp_drop_nodes_sorted = os.path.join(args.outdir, "._tmp.drop.nodes.sorted")
    tmp_drop_edges_sorted = os.path.join(args.outdir, "._tmp.drop.edges.sorted")

    for p in [
        tmp_keep_nodes_raw, tmp_keep_edges_raw, tmp_keep_nodes_sorted, tmp_keep_edges_sorted,
        tmp_drop_nodes_raw, tmp_drop_edges_raw, tmp_drop_nodes_sorted, tmp_drop_edges_sorted
    ]:
        if os.path.exists(p):
            os.remove(p)

    OUT_COLS = [
        "Tissue", "Motif", "Enhancer", "Promoter", "Gene",
        "Motif_Enhancer", "Enhancer_Gene", "Promoter_Gene", "Method"
    ]

    def motif_key_for_overlap(m: str) -> str:
        m = (m or "").strip()
        if m in unresolved_set:
            return m
        return motif_base(m)

    homer_keys = set()
    with open(args.homer, "r", encoding="utf-8") as f:
        hdr = split_fields(f.readline())
        idx = {norm_col(c): i for i, c in enumerate(hdr)}
        it = require_idx(idx, "Tissue")
        im = require_idx(idx, "Motif")
        ie = require_idx(idx, "Enhancer")
        ig = require_idx(idx, "Gene")

        for line in f:
            if not line.strip():
                continue
            a = split_fields(line)
            if len(a) <= max(it, im, ie, ig):
                continue
            tissue = a[it]
            motif = a[im]
            enh = a[ie]
            gene = a[ig]
            key = (tissue, motif_key_for_overlap(motif), enh, gene)
            homer_keys.add(key)

    overlap_keys = set()
    with open(args.dl, "r", encoding="utf-8") as f:
        hdr = split_fields(f.readline())
        idx = {norm_col(c): i for i, c in enumerate(hdr)}
        it = require_idx(idx, "Tissue")
        im = require_idx(idx, "Motif")
        ie = require_idx(idx, "Enhancer")
        ig = require_idx(idx, "Gene")

        for line in f:
            if not line.strip():
                continue
            a = split_fields(line)
            if len(a) <= max(it, im, ie, ig):
                continue
            tissue = a[it]
            motif = a[im]
            enh = a[ie]
            gene = a[ig]
            key = (tissue, motif_key_for_overlap(motif), enh, gene)
            if key in homer_keys:
                overlap_keys.add(key)

    print(f"[INFO] homer unique overlap-keys = {len(homer_keys)}")
    print(f"[INFO] inferred overlap keys (homer ∩ dl, promoter ignored) = {len(overlap_keys)}")

    def dl_motif_out_drop_hash(m: str) -> str:
        m = (m or "").strip()
        if m in unresolved_set:
            return m
        return motif_base(m)

    #    key = (tissue, motif_out, enh, pro, gene, method)
    AggVal = Tuple[float, float, float]  # sum_me, max_eg, max_pg

    keep_agg: "OrderedDict[Tuple[str, str, str, str, str, str], AggVal]" = OrderedDict()
    drop_agg: "OrderedDict[Tuple[str, str, str, str, str, str], AggVal]" = OrderedDict()

    def agg_add(agg: OrderedDict, key: Tuple[str, str, str, str, str, str],
                w_me_s: str, w_eg_s: str, w_pg_s: str):
        w_me = to_float(w_me_s)
        w_eg = to_float(w_eg_s)
        w_pg = to_float(w_pg_s)

        if key not in agg:
            agg[key] = (
                0.0 if w_me != w_me else w_me,  # sum_me
                float("-inf") if w_eg != w_eg else w_eg,  # max_eg
                float("-inf") if w_pg != w_pg else w_pg,  # max_pg
            )
        else:
            s_me, m_eg, m_pg = agg[key]
            if w_me == w_me:
                s_me += w_me
            if w_eg == w_eg:
                m_eg = max(m_eg, w_eg)
            if w_pg == w_pg:
                m_pg = max(m_pg, w_pg)
            agg[key] = (s_me, m_eg, m_pg)

    kept_homer = kept_overlap = kept_dl = 0

    with open(out_keep_raw, "w", encoding="utf-8") as w_keep_raw, \
         open(out_drop_raw, "w", encoding="utf-8") as w_drop_raw:

        w_keep_raw.write("\t".join(OUT_COLS) + "\n")
        w_drop_raw.write("\t".join(OUT_COLS) + "\n")

        with open(args.homer, "r", encoding="utf-8") as f:
            hdr = split_fields(f.readline())
            idx = {norm_col(c): i for i, c in enumerate(hdr)}

            it = require_idx(idx, "Tissue")
            im = require_idx(idx, "Motif")
            ie = require_idx(idx, "Enhancer")
            ip = require_idx(idx, "Promoter")
            ig = require_idx(idx, "Gene")
            i_me = require_idx(idx, "Motif_Enhancer")
            i_eg = require_idx(idx, "Enhancer_Gene")
            i_pg = require_idx(idx, "Promoter_Gene")

            for line in f:
                if not line.strip():
                    continue
                a = split_fields(line)
                if len(a) <= max(i_pg, it, im, ie, ip, ig):
                    continue

                tissue = a[it]
                motif = a[im]
                enh = a[ie]
                pro = a[ip]
                gene = a[ig]
                w_me = a[i_me]
                w_eg = a[i_eg]
                w_pg = a[i_pg]

                okey = (tissue, motif_key_for_overlap(motif), enh, gene)
                method = "overlap" if okey in overlap_keys else "homer"

                motif_keep = motif
                w_keep_raw.write("\t".join([tissue, motif_keep, enh, pro, gene, w_me, w_eg, w_pg, method]) + "\n")
                agg_add(keep_agg, (tissue, motif_keep, enh, pro, gene, method), w_me, w_eg, w_pg)

                motif_drop = motif
                w_drop_raw.write("\t".join([tissue, motif_drop, enh, pro, gene, w_me, w_eg, w_pg, method]) + "\n")
                agg_add(drop_agg, (tissue, motif_drop, enh, pro, gene, method), w_me, w_eg, w_pg)

                if method == "overlap":
                    kept_overlap += 1
                else:
                    kept_homer += 1

        with open(args.dl, "r", encoding="utf-8") as f:
            hdr = split_fields(f.readline())
            idx = {norm_col(c): i for i, c in enumerate(hdr)}

            it = require_idx(idx, "Tissue")
            im = require_idx(idx, "Motif")
            ie = require_idx(idx, "Enhancer")
            ip = require_idx(idx, "Promoter")
            ig = require_idx(idx, "Gene")
            i_me = require_idx(idx, "Motif_Enhancer")
            i_eg = require_idx(idx, "Enhancer_Gene")
            i_pg = require_idx(idx, "Promoter_Gene")

            for line in f:
                if not line.strip():
                    continue
                a = split_fields(line)
                if len(a) <= max(i_pg, it, im, ie, ip, ig):
                    continue

                tissue = a[it]
                motif = a[im]
                enh = a[ie]
                pro = a[ip]
                gene = a[ig]
                w_me = a[i_me]
                w_eg = a[i_eg]
                w_pg = a[i_pg]

                okey = (tissue, motif_key_for_overlap(motif), enh, gene)

                if (okey in overlap_keys) and (motif not in unresolved_set):
                    continue

                method = "deeplearning"

                motif_keep = motif
                w_keep_raw.write("\t".join([tissue, motif_keep, enh, pro, gene, w_me, w_eg, w_pg, method]) + "\n")
                agg_add(keep_agg, (tissue, motif_keep, enh, pro, gene, method), w_me, w_eg, w_pg)

                motif_drop = dl_motif_out_drop_hash(motif)
                w_drop_raw.write("\t".join([tissue, motif_drop, enh, pro, gene, w_me, w_eg, w_pg, method]) + "\n")
                agg_add(drop_agg, (tissue, motif_drop, enh, pro, gene, method), w_me, w_eg, w_pg)

                kept_dl += 1

    print(f"[INFO] wrote raw tables")
    print(f"[INFO] homer rows kept as homer    = {kept_homer}")
    print(f"[INFO] homer rows kept as overlap  = {kept_overlap}")
    print(f"[INFO] dl rows kept as deeplearning= {kept_dl}")

    def write_dedup_table(path: str, agg: OrderedDict):
        with open(path, "w", encoding="utf-8") as w:
            w.write("\t".join(OUT_COLS) + "\n")
            for (tissue, motif, enh, pro, gene, method), (sum_me, max_eg, max_pg) in agg.items():
                w_me = fmt_num(sum_me)
                w_eg = fmt_num(max_eg)
                w_pg = fmt_num(max_pg)
                w.write("\t".join([tissue, motif, enh, pro, gene, w_me, w_eg, w_pg, method]) + "\n")

    write_dedup_table(out_keep_dedup, keep_agg)
    write_dedup_table(out_drop_dedup, drop_agg)
    print(f"[INFO] wrote dedup tables")

    def build_nodes_edges_raw(agg: OrderedDict, nodes_raw: str, edges_raw: str):
        with open(nodes_raw, "w", encoding="utf-8") as wn, open(edges_raw, "w", encoding="utf-8") as we:
            for (tissue, motif, enh, pro, gene, method), (sum_me, max_eg, max_pg) in agg.items():
                w_me = fmt_num(sum_me)
                w_eg = fmt_num(max_eg)
                w_pg = fmt_num(max_pg)

                wn.write(f"{motif}\tMotif\t{tissue}\t{method}\n")
                wn.write(f"{enh}\tEnhancer\t{tissue}\t{method}\n")
                wn.write(f"{pro}\tPromoter\t{tissue}\t{method}\n")
                wn.write(f"{gene}\tGene\t{tissue}\t{method}\n")

                we.write(f"{motif}\t{enh}\t{w_me}\t{tissue}\t{method}\n")
                we.write(f"{enh}\t{gene}\t{w_eg}\t{tissue}\t{method}\n")
                we.write(f"{pro}\t{gene}\t{w_pg}\t{tissue}\t{method}\n")

    build_nodes_edges_raw(keep_agg, tmp_keep_nodes_raw, tmp_keep_edges_raw)
    build_nodes_edges_raw(drop_agg, tmp_drop_nodes_raw, tmp_drop_edges_raw)

    sort_u(tmp_keep_nodes_raw, tmp_keep_nodes_sorted)
    sort_u(tmp_keep_edges_raw, tmp_keep_edges_sorted)
    with open(out_keep_nodes, "w", encoding="utf-8") as w:
        w.write("id\ttype\ttissue\tmethod\n")
        with open(tmp_keep_nodes_sorted, "r", encoding="utf-8") as r:
            shutil.copyfileobj(r, w)
    with open(out_keep_edges, "w", encoding="utf-8") as w:
        w.write("node1\tnode2\tweight\ttissue\tmethod\n")
        with open(tmp_keep_edges_sorted, "r", encoding="utf-8") as r:
            shutil.copyfileobj(r, w)

    sort_u(tmp_drop_nodes_raw, tmp_drop_nodes_sorted)
    sort_u(tmp_drop_edges_raw, tmp_drop_edges_sorted)
    with open(out_drop_nodes, "w", encoding="utf-8") as w:
        w.write("id\ttype\ttissue\tmethod\n")
        with open(tmp_drop_nodes_sorted, "r", encoding="utf-8") as r:
            shutil.copyfileobj(r, w)
    with open(out_drop_edges, "w", encoding="utf-8") as w:
        w.write("node1\tnode2\tweight\ttissue\tmethod\n")
        with open(tmp_drop_edges_sorted, "r", encoding="utf-8") as r:
            shutil.copyfileobj(r, w)

    for p in [
        tmp_keep_nodes_raw, tmp_keep_edges_raw, tmp_keep_nodes_sorted, tmp_keep_edges_sorted,
        tmp_drop_nodes_raw, tmp_drop_edges_raw, tmp_drop_nodes_sorted, tmp_drop_edges_sorted
    ]:
        if os.path.exists(p):
            os.remove(p)

    print("[DONE] outputs written to:", args.outdir)
    print(" keepHash:")
    print("  -", out_keep_raw)
    print("  -", out_keep_dedup)
    print("  -", out_keep_nodes, "(nodes: id type tissue method)")
    print("  -", out_keep_edges, "(edges: node1 node2 weight tissue method)")
    print(" dropHash (DL only; unresolved keeps #):")
    print("  -", out_drop_raw)
    print("  -", out_drop_dedup)
    print("  -", out_drop_nodes)
    print("  -", out_drop_edges)


if __name__ == "__main__":
    main()
