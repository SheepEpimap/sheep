# Workflow and paths

## Execution order

| Order | Stage | Canonical implementation | Main output |
|---:|---|---|---|
| 1 | `cre_01_reclassify` | `src/cre_pipeline/01_reclassify_e1_e9.sh` | classified CRE BEDs and classification summary |
| 2 | `cre_02_hg38_to_hg19` | `src/cre_pipeline/02_hg38_to_hg19.sh` | hg19 classified and ssCRE BEDs |
| 3 | `cre_03_collect_e5` | `src/cre_pipeline/03_collect_e5_hg19.sh` | E5 collection and summary |
| 4 | `cre_04_reuse_human_tsr` | `src/cre_pipeline/04_reuse_human_tsr.sh` | validated TSR source link or copy |
| 5 | `cre_05_label_tsr` | `src/cre_pipeline/05_label_cre_with_tsr.sh` | TSR-labelled E5 CRE BEDs |
| 6 | `cre_06_reuse_projection` | `src/cre_pipeline/06_reuse_sheep_projection.sh` | validated sheep projection link or copy |
| 7 | `cre_07_split_ldsc` | `src/cre_pipeline/07_split_ldsc_beds.sh` | LDSC-ready BED slices |
| 8 | `cre_08_collect_gwas_beds` | `src/cre_pipeline/08_collect_gwas_beds.sh` | final 68 BED collection |
| 9 | `gwas_prepare` | `src/gwas_ldsc/prepare_E1E9_reuse_existing_recalc.sh` | controls, normalized manifests and LDCTS files |
| 10 | `gwas_submit_slurm` | `src/gwas_ldsc/submit_E1E9_reuse_existing_152.sh` | LDSC arrays and dependent summaries |
| 11 | `gwas_summary` | `src/gwas_ldsc/run_summaries.sh` | summary/QC tables and final GWAS plotting input |
| 12 | `figure_input_merge` | `src/sheep_epimap/figure_inputs/build_merged_classification.py` | merged hg38 BEDs and classification summary |
| 13 | `figure_input_conservation` | `src/sheep_epimap/figure_inputs/score_conservation.py` | per-CRE phastCons, phyloP and GERP table |
| 14 | `fig8` | `src/sheep_epimap/figures/fig8.py` | panels b, d, e, f and g |

`cre_submit_slurm` submits steps 1–8 with `afterok` dependencies and the resource requests recorded for the completed historical jobs. `gwas_submit_slurm` chunks arrays at no more than 1,000 tasks and 80 concurrent tasks by default.

## Path rules

- Executable source contains no project-specific absolute path.
- Copy `config/paths.example.tsv` to an ignored local/server config and set paths there.
- `SRC_WORKDIR` and `OLD_GWAS_BASE` are read-only historical inputs.
- New run outputs are derived from `RUNS_ROOT/<run_id>/` unless an existing result tree is supplied for figure-only use.
- The conservation table is not derived by string substitution. Supply the verified table explicitly.
- Server mount paths reflect the original compute environment and are not portable defaults.

## Manifest contract

Every active LDSC array reads a header-bearing TSV. The accepted first columns are `annotation_id` for annotation, LD-score and h² manifests, and `ldcts_name` for the CTS manifest. The preparation step adds the registered header to a historical headerless manifest and the submitter validates it before calculating array size.

## Not executed during release preparation

No BED transformation, liftOver, bedtools intersection, LDSC calculation, notebook, summary script, Slurm submission or biological-result comparison was run while constructing this repository.
