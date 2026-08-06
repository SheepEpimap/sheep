# Analysis workflows and figure relationships

The repository provides two complementary views of the project:

- `pipelines/` contains stage-oriented analysis code;
- `figures/` contains code that can be defensibly associated with a figure or
  panel.

The stage-oriented workflows provide upstream analyses, while the figure tree
provides panel-focused entry points.

| Workflow | Main purpose | Figure relationship |
|---|---|---|
| `1.Data preprocessing` | ATAC-seq and CUT&Tag processing, peak calling, QC, signal matrices, TSS profiles and assay PCA | Supports Fig. 1, Extended Data Fig. 1 and Supplementary Figs. 1-3; the PCA scripts directly support Supplementary Fig. 1 |
| `2.ChromImpute_ChromHMM` | Epigenome imputation and 11-state chromatin modelling | Supports Fig. 2 and Supplementary Figs. 4-8 |
| `3.TSR` | Tissue-specific regulatory regions, enhancer-gene links, enrichment, motifs and tissue-specificity | Supports Fig. 2 and regulatory-element characterization panels |
| `4.Network` | Regulatory-network construction and network statistics | Supports Fig. 4; panel-specific network tables are also available under `figures/Figure_04/` |
| `5.Selection signatures` | Selection scans, GAT enrichment and regional FST analysis | Supports Fig. 7; the exact tissue-specific enhancer enrichment panel is under `figures/Figure_07/07c/` |

## Recommended execution pattern

1. Configure project, reference and data paths.
2. Run the required upstream pipeline stages.
3. Read the target panel README and satisfy its documented input schema.
4. Run panel scripts in numeric order.
5. Run any validation script in the panel directory.
6. Compare the generated panel content with the final RGB PDF before final
   multipanel assembly.

Absolute paths and scheduler commands reflect the original HPC environment and must be
reviewed before use on another system.
