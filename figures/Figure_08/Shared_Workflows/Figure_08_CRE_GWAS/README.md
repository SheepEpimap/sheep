# Sheep Epimap E1–E9 workflow

This repository is the curated code release for the sheep regulatory-element analysis used to classify E1–E9 CREs, prepare LDSC annotations, run the 152-trait GWAS enrichment workflow, summarize results, and reproduce the code-derived panels of Figure 8.

The code was consolidated from the executed server lineage and the final plotting notebooks. Large data, references, LDSC software, historical results, and the manually assembled final `fig8.tif` are intentionally not included.

## Canonical workflow

1. Reuse the verified historical lifted human intervals and sheep reference intervals.
2. Reclassify E1–E9 into sfCRE, sdCRE and soCRE while retaining unmapped human intervals as ssCRE.
3. Convert classified human intervals from hg38 to hg19.
4. collect E5 CREs and reuse the validated human TSR labels.
5. label human E5 CREs with TSR classes.
6. validate E5 sfCRE identity before reusing the sheep projection.
7. create LDSC-ready human and sheep-projection BED slices and collect the 68 final BED files.
8. prepare normalized, header-bearing LDSC manifests and control annotations.
9. submit annotation, LD-score, partitioned-heritability and CTS arrays for 152 traits.
10. build the merged-deduplicated hg38 table and, when requested, score conservation bigWigs.
11. summarize LDSC results and generate Figure 8 panels b, d, e, f and g.

## Quick start

```bash
cp config/paths.example.tsv config/paths.server.tsv
# Edit config/paths.server.tsv with paths available on your system.

bash run_workflow.sh config/paths.server.tsv YYYYMMDD_label --dry-run
```

The executable launcher above runs the following three workflow entry points in order:

```bash
python workflow.py --paths config/paths.server.tsv --stage preflight \
  --run-id YYYYMMDD_label
python workflow.py --paths config/paths.server.tsv --stage cre_submit_slurm \
  --run-id YYYYMMDD_label --dry-run
python workflow.py --paths config/paths.server.tsv --stage gwas_submit_slurm \
  --run-id YYYYMMDD_label --dry-run
```

Figure generation is run separately:

```bash
python workflow.py --paths config/paths.server.tsv --stage fig8 \
  --panels b d e f g --dry-run
```

Replace `--dry-run` with `--execute` in `run_workflow.sh` only after the path preflight has passed and the intended run ID has been reviewed. Both Slurm submission stages require `--run-id`; completed stages receive immutable marker files and are not overwritten.

## Historical corrections encoded here

- The active E1–E9 chain uses three reuse scripts and five shared scripts from the historical blank branch; the blank branch is not treated as wholly obsolete.
- E10 and E11 are excluded from the E1–E9 "other state" evidence.
- The final GWAS manifests have explicit headers. Submission stops if the header does not match the registered schema. This removes the historical A2/B2 first-row skip and empty-array-task failure.
- Figure 8 panel g colors points only when `Coefficient > 0` and `P < 0.05`.
- The conservation cache is an explicit input because the final notebook used an older in-kernel `WORKDIR`; the code does not infer a replacement path.

## Reproducibility boundary

The source code has been structurally checked and its Python and shell syntax validated. The biological analyses were **not rerun** during packaging. Panels a, c, h and i and the final composite `fig8.tif` do not have a confirmed code-generation lineage and are therefore not claimed as reproducible by this repository.

See [docs/WORKFLOW_AND_PATHS.md](docs/WORKFLOW_AND_PATHS.md), [docs/FIGURE_PROVENANCE.md](docs/FIGURE_PROVENANCE.md), and the tables under `manifests/` before execution or public release.
