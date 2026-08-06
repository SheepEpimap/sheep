# Figure-code coverage

This table reports conservative coverage for all 34 supplied final RGB PDFs.
“Exact/partial” means at least one panel has defensible panel-specific code; it
does not mean the entire composite figure can be rebuilt with one command.

| Figure | Subject | Coverage |
|---|---|---|
| Figure 1 | Sheep multi-omics regulatory atlas | Upstream pipeline support; no exact final panel assembly |
| Figure 2 | Chromatin states and enhancer-associated regulatory features | Exact/partial: 2c, 2d, 2e and 2j-2k |
| Figure 3 | ChromBPNet sequence features, motifs and footprints | Exact/partial: 3a-3h |
| Figure 4 | ChromBPNet-derived regulatory networks | Inputs/partial panels: 4a-4c |
| Figure 5 | Allele-specific signals and trait relevance | Exact/partial: 5d and 5f-5l |
| Figure 6 | Variant and trait interpretation | Shared analysis support; no exact final panel command |
| Figure 7 | Selection signals | Exact panel code for 7c; upstream support for other panels |
| Figure 8 | Cross-species CRE conservation and human traits | Exact/input code for 8a-8g |
| Extended Data Fig. 1 | Multi-omics overview and quality metrics | Upstream preprocessing/QC support; no exact final panel assembly |
| Extended Data Fig. 2 | Super-enhancers and their functions | No exact supplied panel code |
| Extended Data Fig. 3 | Enhancer-target gene-link evaluation | Exact/partial: 3a-3n |
| Extended Data Fig. 4 | ChromBPNet performance and motif-instance features | Exact/partial: 4a-4h |
| Extended Data Fig. 5 | Additional regulatory networks and loci | Shared engines only; exact regions/layouts were not supplied |
| Extended Data Fig. 6 | Enhancer sharing, connectivity, conservation and variation | Exact/partial: 6a-6n |
| Extended Data Fig. 7 | PIK3CD regulatory locus | Shared engines only; exact invocation was not supplied |
| Extended Data Fig. 8 | GAA and CAPN5 selected-variant loci | Shared engines only; exact invocations were not supplied |
| Supplementary Fig. 1 | Multi-assay PCA across tissues | Direct upstream PCA workflow available |
| Supplementary Fig. 2 | Signals across gene classes | Upstream signal-processing support; no exact final panel script |
| Supplementary Fig. 3 | RXFP2 locus tracks | Shared track support; no exact final invocation |
| Supplementary Fig. 4 | Epigenome-imputation validation | Upstream analysis support; no exact final panel scripts |
| Supplementary Fig. 5 | Imputed-profile QC | Upstream analysis support; no exact final panel scripts |
| Supplementary Fig. 6 | ChromHMM model selection and characterization | Upstream analysis support; no exact final panel scripts |
| Supplementary Fig. 7 | Genomic enrichment, tissue sharing and state switching | Upstream analysis support; no exact final panel scripts |
| Supplementary Fig. 8 | Conservation, accessibility and methylation profiles | Exact code for panel 8d only |
| Supplementary Fig. 9 | Super-enhancer summary | No exact supplied panel code |
| Supplementary Fig. 10 | Tissue-specific regulatory-element distribution | Upstream analysis support; no exact final panel code |
| Supplementary Fig. 11 | Mouse phenotype and VISTA enhancer enrichment | Supporting enrichment/heatmap workflow; no complete exact composite code |
| Supplementary Fig. 12 | Additional motif-instance features and AGPAT3 | Exact panel code for 12a; no exact 12b invocation |
| Supplementary Fig. 13 | Footprint distributions and tissue-pair GO heatmaps | Exact panel code for 13a; no exact 13b-13c code |
| Supplementary Fig. 14 | Cross-assay allele-specific signal summaries | No exact supplied panel code |
| Supplementary Fig. 15 | GWAS enrichment in AS/non-AS enhancers | No exact supplied panel code |
| Supplementary Fig. 16 | CRE conservation and class-specific trait enrichment | Figure 8 workflows generate supporting inputs; final plots are incomplete |
| Supplementary Fig. 17 | Cross-species conserved loci | Shared track support; exact invocations were not supplied |
| Supplementary Fig. 18 | Regional human GWAS plots | Shared association plotting support; exact invocations were not supplied |

For script-level detail, use [`../FIGURE_CODE_MAP.tsv`](../FIGURE_CODE_MAP.tsv).
For explicit gaps, use [`UNMATCHED_PANELS.tsv`](UNMATCHED_PANELS.tsv).
