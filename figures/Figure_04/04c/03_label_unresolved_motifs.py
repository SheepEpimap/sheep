#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

"""
01_unresolved_to_like_and_filter_fullChain.py

input:
  1) /data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/output/TABLE_motif_compendium.csv
     -
     - motif_name: unresolved motif
     - category: resolved/unresolved
     - TF:   motif( “ ”)

  2) /data/home/sczd644/run/zsw_chrombpnet/network/TSR/ALL_fullChain/ALL.all_tables.with_header.fullChain.txt
     - tab  (  tab  )
     -  :Tissue Motif Enhancer Promoter Gene Motif_Enhancer Enhancer_Gene Promoter_Gene

output(  7  file,  outdir):
  0) unresolved -> TF-like  ( )
     - unresolved_to_TF_like.tsv   # unresolved_motif  TF  TF_like

  A) unresolved(  motif  )
     - unresolved.all_tables.with_header.fullChain.txt
     - unresolved.edge.fullChain.txt   #  :id type tissue
     - unresolved.node.fullChain.txt   #  :node1 node2 weight tissue

  B) TF-like(  unresolved motif   TF-like / TF#i-like)
     - TF_like.all_tables.with_header.fullChain.txt
     - TF_like.edge.fullChain.txt
     - TF_like.node.fullChain.txt

 :
  -   unresolved motif   TF:
       1  unresolved -> TF#1-like
       2  unresolved -> TF#2-like...
    (  compendium  )
"""

import argparse
import csv
import os
import shutil
import subprocess
from collections import OrderedDict


def which(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def norm_col(s: str) -> str:
    """ :  +  / """
    return s.strip().lower().replace(" ", "").replace("-", "").replace("_", "")


def split_line(line: str):
    """  tab  , """
    line = line.rstrip("\n")
    if "\t" in line:
        return line.split("\t")
    return line.split()


def sh_quote(s: str) -> str:
    return "'" + s.replace("'", "'\"'\"'") + "'"


def load_unresolved_mapping(compendium_csv: str):
    """
     :
      unresolved_order: list(unresolved_motif)  #   compendium  ,
      unresolved_set: set(unresolved_motif)
      unresolved_to_like: dict(unresolved_motif -> TF-like / TF#i-like)
      unresolved_to_tf: dict(unresolved_motif -> TF)

     :
      -   category == 'unresolved'
      - TF   TF
      -   unresolved   TF:TF#1-like, TF#2-like...
    """
    unresolved_rows = []
    with open(compendium_csv, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        need = {"motif_name", "category", "TF"}
        got = set(reader.fieldnames or [])
        missing = need - got
        if missing:
            raise RuntimeError(f"[ERROR] compendium  :{missing}; :{reader.fieldnames}")

        for row in reader:
            cat = (row.get("category") or "").strip().lower()
            if cat != "unresolved":
                continue
            m = (row.get("motif_name") or "").strip()
            tf = (row.get("TF") or "").strip()
            if not m or not tf:
                continue
            unresolved_rows.append((m, tf))

    unresolved_order = []
    unresolved_to_tf = {}
    seen_motif = set()
    for m, tf in unresolved_rows:
        if m in seen_motif:
            continue
        seen_motif.add(m)
        unresolved_order.append(m)
        unresolved_to_tf[m] = tf

    tf_to_motifs = OrderedDict()
    for m in unresolved_order:
        tf = unresolved_to_tf[m]
        tf_to_motifs.setdefault(tf, [])
        tf_to_motifs[tf].append(m)

    # unresolved -> TF-like / TF#i-like
    unresolved_to_like = {}
    for tf, motifs in tf_to_motifs.items():
        if len(motifs) == 1:
            unresolved_to_like[motifs[0]] = f"{tf}-like"
        else:
            for i, m in enumerate(motifs, start=1):
                unresolved_to_like[m] = f"{tf}#{i}-like"

    unresolved_set = set(unresolved_to_like.keys())
    return unresolved_order, unresolved_set, unresolved_to_like, unresolved_to_tf


def dedup_with_sort_u(infile: str, outfile: str, header: str):
    """
      sort -u  ( ).
      sort,  Python set  ( ).
    """
    tmp_sorted = outfile + ".tmp.sorted"
    os.makedirs(os.path.dirname(outfile), exist_ok=True)

    if which("sort"):
        cmd = ["bash", "-lc", f"LC_ALL=C sort -u {sh_quote(infile)} > {sh_quote(tmp_sorted)}"]
        r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if r.returncode != 0:
            raise RuntimeError(f"[ERROR] sort -u  :\n{r.stderr}")

        with open(outfile, "w", encoding="utf-8") as w:
            w.write(header.rstrip("\n") + "\n")
            with open(tmp_sorted, "r", encoding="utf-8") as r2:
                shutil.copyfileobj(r2, w)

        os.remove(tmp_sorted)
    else:
        seen = set()
        with open(outfile, "w", encoding="utf-8") as w:
            w.write(header.rstrip("\n") + "\n")
            with open(infile, "r", encoding="utf-8") as r:
                for line in r:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    if line in seen:
                        continue
                    seen.add(line)
                    w.write(line + "\n")


def write_mapping_table(out_path: str, unresolved_order, unresolved_to_tf, unresolved_to_like):
    """output unresolved -> TF -> TF_like  """
    with open(out_path, "w", encoding="utf-8") as w:
        w.write("unresolved_motif\tTF\tTF_like\n")
        for m in unresolved_order:
            tf = unresolved_to_tf.get(m, "")
            like = unresolved_to_like.get(m, "")
            w.write(f"{m}\t{tf}\t{like}\n")


def main():
    ap = argparse.ArgumentParser(
        description="Filter unresolved motifs from ALL_fullChain and build TF-like renamed version (7 outputs)."
    )
    ap.add_argument("--compendium",
                    default="/data/home/sczd644/run/zsw_chrombpnet/03-syntax/01/output/TABLE_motif_compendium.csv")
    ap.add_argument("--fullchain",
                    default="/data/home/sczd644/run/zsw_chrombpnet/network/TSR/ALL_fullChain/ALL.all_tables.with_header.fullChain.txt")
    ap.add_argument("--outdir",
                    default="/data/home/sczd644/run/zsw_chrombpnet/network/TSR/ALL_fullChain/UNRESOLVED_6files")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    unresolved_order, unresolved_set, unresolved_to_like, unresolved_to_tf = load_unresolved_mapping(args.compendium)
    print(f"[INFO] unresolved motifs in compendium: {len(unresolved_set)}")

    mapping_out = os.path.join(args.outdir, "unresolved_to_TF_like.tsv")
    write_mapping_table(mapping_out, unresolved_order, unresolved_to_tf, unresolved_to_like)
    print(f"[INFO] wrote mapping: {mapping_out}")

    if len(unresolved_set) == 0:
        print("[WARN] compendium   category==unresolved   motif; output file( ).")

    raw_table = os.path.join(args.outdir, "unresolved.all_tables.with_header.fullChain.txt")
    raw_edge = os.path.join(args.outdir, "unresolved.edge.fullChain.txt")  #
    raw_node = os.path.join(args.outdir, "unresolved.node.fullChain.txt")  #

    like_table = os.path.join(args.outdir, "TF_like.all_tables.with_header.fullChain.txt")
    like_edge = os.path.join(args.outdir, "TF_like.edge.fullChain.txt")
    like_node = os.path.join(args.outdir, "TF_like.node.fullChain.txt")

    raw_edge_tmp = raw_edge + ".tmp"
    raw_node_tmp = raw_node + ".tmp"
    like_edge_tmp = like_edge + ".tmp"
    like_node_tmp = like_node + ".tmp"

    for p in [raw_edge_tmp, raw_node_tmp, like_edge_tmp, like_node_tmp]:
        if os.path.exists(p):
            os.remove(p)

    with open(args.fullchain, "r", encoding="utf-8") as f:
        header_line = f.readline()
        if not header_line:
            raise RuntimeError("[ERROR] fullChain file ")

        header = split_line(header_line)
        col_index = {norm_col(name): i for i, name in enumerate(header)}

        def must_col(name):
            k = norm_col(name)
            if k not in col_index:
                raise RuntimeError(f"[ERROR] fullChain  :{name}; :{header}")
            return col_index[k]

        i_tissue = must_col("Tissue")
        i_motif = must_col("Motif")
        i_enh = must_col("Enhancer")
        i_pro = must_col("Promoter")
        i_gene = must_col("Gene")
        i_w_me = must_col("Motif_Enhancer")
        i_w_eg = must_col("Enhancer_Gene")
        i_w_pg = must_col("Promoter_Gene")

        with open(raw_table, "w", encoding="utf-8") as w_raw, \
             open(like_table, "w", encoding="utf-8") as w_like, \
             open(raw_edge_tmp, "w", encoding="utf-8") as raw_edge_w, \
             open(raw_node_tmp, "w", encoding="utf-8") as raw_node_w, \
             open(like_edge_tmp, "w", encoding="utf-8") as like_edge_w, \
             open(like_node_tmp, "w", encoding="utf-8") as like_node_w:

            w_raw.write("\t".join(header) + "\n")
            w_like.write("\t".join(header) + "\n")

            kept = 0
            bad = 0
            for line in f:
                if not line.strip():
                    continue
                fields = split_line(line)

                if len(fields) < len(header):
                    bad += 1
                    continue

                motif = fields[i_motif]
                if motif not in unresolved_set:
                    continue

                tissue = fields[i_tissue]
                enh = fields[i_enh]
                pro = fields[i_pro]
                gene = fields[i_gene]

                w_me = fields[i_w_me]
                w_eg = fields[i_w_eg]
                w_pg = fields[i_w_pg]

                # --- raw table ---
                w_raw.write("\t".join(fields) + "\n")

                motif_like = unresolved_to_like.get(motif, motif)
                fields_like = list(fields)
                fields_like[i_motif] = motif_like
                w_like.write("\t".join(fields_like) + "\n")

                # raw
                raw_edge_w.write(f"{motif}\tMotif\t{tissue}\n")
                raw_edge_w.write(f"{enh}\tEnhancer\t{tissue}\n")
                raw_edge_w.write(f"{pro}\tPromoter\t{tissue}\n")
                raw_edge_w.write(f"{gene}\tGene\t{tissue}\n")

                # like
                like_edge_w.write(f"{motif_like}\tMotif\t{tissue}\n")
                like_edge_w.write(f"{enh}\tEnhancer\t{tissue}\n")
                like_edge_w.write(f"{pro}\tPromoter\t{tissue}\n")
                like_edge_w.write(f"{gene}\tGene\t{tissue}\n")

                # raw
                raw_node_w.write(f"{motif}\t{enh}\t{w_me}\t{tissue}\n")
                raw_node_w.write(f"{enh}\t{gene}\t{w_eg}\t{tissue}\n")
                raw_node_w.write(f"{pro}\t{gene}\t{w_pg}\t{tissue}\n")

                # like
                like_node_w.write(f"{motif_like}\t{enh}\t{w_me}\t{tissue}\n")
                like_node_w.write(f"{enh}\t{gene}\t{w_eg}\t{tissue}\n")
                like_node_w.write(f"{pro}\t{gene}\t{w_pg}\t{tissue}\n")

                kept += 1

    print(f"[INFO] kept unresolved rows from fullChain: {kept}")

    dedup_with_sort_u(raw_edge_tmp, raw_edge, header="id\ttype\ttissue")
    dedup_with_sort_u(raw_node_tmp, raw_node, header="node1\tnode2\tweight\ttissue")
    dedup_with_sort_u(like_edge_tmp, like_edge, header="id\ttype\ttissue")
    dedup_with_sort_u(like_node_tmp, like_node, header="node1\tnode2\tweight\ttissue")

    for p in [raw_edge_tmp, raw_node_tmp, like_edge_tmp, like_node_tmp]:
        if os.path.exists(p):
            os.remove(p)

    print("[DONE] outputs written to:", args.outdir)
    print("  - unresolved_to_TF_like.tsv")
    print("  - unresolved.all_tables.with_header.fullChain.txt")
    print("  - unresolved.edge.fullChain.txt   (nodes: id type tissue)")
    print("  - unresolved.node.fullChain.txt   (edges: node1 node2 weight tissue)")
    print("  - TF_like.all_tables.with_header.fullChain.txt")
    print("  - TF_like.edge.fullChain.txt")
    print("  - TF_like.node.fullChain.txt")


if __name__ == "__main__":
    main()
