# Human-to-sheep CRE classification

Mapping status: **supporting**.

Run scripts in numeric order after reviewing the cluster-specific paths and input schemas.

| Order | Script | Role |
|---:|---|---|
| 01 | `01_prepare_human_chromatin_states.sh` | Convert human hg38 E1-E9 intervals to sheep ARS-UI_Ramb_v3.0 (`GCF_016772045.2`) with `minMatch=0.1`, then classify them in sheep V3 coordinates. |
| 02 | `02_convert_hg38_hg19_coordinates.sh` | Convert human assemblies for cross-species comparison. |

Set `HUMAN_TO_SHEEP_CHAIN`, `SHEEP_V3_CHROM_MAP`,
`HUMAN_HG38_DIR`, `SHEEP_V3_DIR`, and `WORKDIR` before running step 01.
The reverse `SHEEP_TO_HG38_CHAIN` is used by Figure 8a and is deliberately a
separate configuration item.

The final Illustrator assembly, typography and panel labels are not generated
by these scripts unless explicitly stated.
