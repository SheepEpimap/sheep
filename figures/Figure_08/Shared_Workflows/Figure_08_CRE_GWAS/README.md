# SheepEpimap E1-E9 CRE and GWAS workflow

This project classifies E1-E9 cis-regulatory elements, prepares LDSC
annotations, runs the 152-trait GWAS enrichment workflow, summarizes results,
and generates Figure 8 panels b, d, e, f and g.

## Workflow

1. Prepare the human and sheep regulatory-element intervals.
2. Classify E1-E9 CREs into sfCRE, sdCRE, soCRE and ssCRE groups.
3. Convert human intervals from hg38 to hg19 for LDSC.
4. Prepare E5 CRE and TSR annotations.
5. Build LDSC annotations, manifests and control sets.
6. Run annotation, LD-score, partitioned-heritability and CTS jobs.
7. Summarize the enrichment results and generate the Figure 8 panels.

## Quick start

```bash
cp config/paths.example.tsv config/paths.server.tsv
# Edit config/paths.server.tsv for the local software, reference and data paths.

bash run_workflow.sh config/paths.server.tsv YYYYMMDD_label --dry-run
```

The launcher runs preflight checks followed by the CRE and GWAS Slurm stages.
Figure generation can be run separately:

```bash
python workflow.py --paths config/paths.server.tsv --stage fig8 \
  --panels b d e f g --dry-run
```

Replace `--dry-run` with `--execute` only after path validation and review of
the run ID. See [docs/WORKFLOW_AND_PATHS.md](docs/WORKFLOW_AND_PATHS.md) for
configuration and execution details.
