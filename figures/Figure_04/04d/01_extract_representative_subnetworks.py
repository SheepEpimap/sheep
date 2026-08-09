#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
# -*- coding: utf-8 -*-

from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

REQUIRED_COLS = [
    "Tissue", "Motif", "Enhancer", "Promoter", "Gene",
    "Motif_Enhancer", "Enhancer_Gene", "Promoter_Gene", "Method"
]

METHODS = ["homer", "deeplearning"]


def read_table_auto(path: Path) -> pd.DataFrame:
    """
     input:
    1)   tab
    2)  /Tab  (  sep=r'\\s+')
    """
    try:
        df = pd.read_csv(path, sep="\t", dtype=str, engine="python", keep_default_na=False)
        if all(c in df.columns for c in REQUIRED_COLS):
            return df
    except Exception:
        pass

    df = pd.read_csv(path, sep=r"\s+", dtype=str, engine="python", keep_default_na=False)
    return df


def motif_base(m: str) -> str:
    """KLF5#1 -> KLF5; # """
    m = "" if m is None else str(m).strip()
    return m.split("#", 1)[0].strip()


def build_edge_nodes_table(df_method: pd.DataFrame) -> pd.DataFrame:
    """
    edge file( ):id type tissue
    """
    edge = pd.concat([
        df_method[["Motif", "Tissue"]].rename(columns={"Motif": "id", "Tissue": "tissue"}).assign(type="Motif"),
        df_method[["Enhancer", "Tissue"]].rename(columns={"Enhancer": "id", "Tissue": "tissue"}).assign(type="Enhancer"),
        df_method[["Promoter", "Tissue"]].rename(columns={"Promoter": "id", "Tissue": "tissue"}).assign(type="Promoter"),
        df_method[["Gene", "Tissue"]].rename(columns={"Gene": "id", "Tissue": "tissue"}).assign(type="Gene"),
    ], ignore_index=True).drop_duplicates()

    return edge[["id", "type", "tissue"]]


def build_node_edges_table(df_method: pd.DataFrame) -> pd.DataFrame:
    """
    node file( ):node1 node2 weight tissue
     :
      Motif -> Enhancer   weight = Motif_Enhancer
      Enhancer -> Gene    weight = Enhancer_Gene
      Promoter -> Gene    weight = Promoter_Gene
    """
    node_me = df_method[["Motif", "Enhancer", "Motif_Enhancer", "Tissue"]].rename(
        columns={"Motif": "node1", "Enhancer": "node2", "Motif_Enhancer": "weight", "Tissue": "tissue"}
    )
    node_eg = df_method[["Enhancer", "Gene", "Enhancer_Gene", "Tissue"]].rename(
        columns={"Enhancer": "node1", "Gene": "node2", "Enhancer_Gene": "weight", "Tissue": "tissue"}
    )
    node_pg = df_method[["Promoter", "Gene", "Promoter_Gene", "Tissue"]].rename(
        columns={"Promoter": "node1", "Gene": "node2", "Promoter_Gene": "weight", "Tissue": "tissue"}
    )

    node = pd.concat([node_me, node_eg, node_pg], ignore_index=True).drop_duplicates()
    return node[["node1", "node2", "weight", "tissue"]]


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Export overlap all_tables + split node/edge for homer and deeplearning."
    )
    ap.add_argument("-i", "--infile", required=True,
                    help="Input: ALL.matched_homer_vs_deeplearning.tsv")
    ap.add_argument("-o", "--outdir", default="overlap_outputs",
                    help="Output directory (default: overlap_outputs)")
    args = ap.parse_args()

    infile = Path(args.infile)
    if (not infile.exists()) or infile.stat().st_size == 0:
        raise SystemExit(f"[ERROR] missing/empty input: {infile}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_table_auto(infile)

    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise SystemExit(f"[ERROR] missing columns: {missing}\nGot columns: {list(df.columns)}")

    for c in REQUIRED_COLS:
        df[c] = df[c].astype(str).str.strip()

    df["Motif_base"] = df["Motif"].map(motif_base)

    # key:tissue + motif_base + enhancer + gene
    key_cols = ["Tissue", "Motif_base", "Enhancer", "Gene"]

    method_sets = df.groupby(key_cols)["Method"].agg(lambda x: set(x)).reset_index(name="MethodSet")
    keys_keep = method_sets[
        method_sets["MethodSet"].apply(lambda s: ("homer" in s) and ("deeplearning" in s))
    ][key_cols]

    df_keep = df.merge(keys_keep, on=key_cols, how="inner")

    out_all = outdir / "overlap_all_tables.with_header.tsv"
    df_keep[REQUIRED_COLS].to_csv(out_all, sep="\t", index=False)

    for method in METHODS:
        df_m = df_keep[df_keep["Method"] == method].copy()

        edge = build_edge_nodes_table(df_m)
        out_edge = outdir / f"overlap_{method}_edge.txt"
        edge.to_csv(out_edge, sep="\t", index=False)

        node = build_node_edges_table(df_m)
        out_node = outdir / f"overlap_{method}_node.txt"
        node.to_csv(out_node, sep="\t", index=False)

        print(f"[OUT] {method}: rows={len(df_m)} nodes={len(edge)} edges={len(node)}")

    print(f"[DONE] overlap(all) rows={len(df_keep)}")
    print(f"[DONE] outdir={outdir.resolve()}")
    print(f"[FILE] {out_all}")


if __name__ == "__main__":
    main()
