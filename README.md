# SheepEpimap figure code

This archive reorganizes only the code supplied during this figure-code review. It is arranged as:

`figures/<figure>/<panel>/<numbered script>`

The archive contains code mapped to 12 of the 34 submitted figure files. Figures without a defensible code match were intentionally omitted. See `FIGURE_CODE_MAP.tsv` for direct, supporting and partial relationships.

## Important conventions

- ChromBPNet compute paths use `/data/home/sczd644/run/zsw_chrombpnet`, matching the HPC workflow supplied by the author.
- The five uploaded `Pasted markdown*.md` files were byte-identical; only one source copy was used.
- The cross-tissue motif workflow is ordered from per-tissue TF-MoDISco outputs through clustering, merging, annotation and preparation of the final 118 unique motifs.
- Numbered scripts encode execution order, not panel order.
- Absolute paths are retained where they encode the original analysis location. Before public release, replace them with configuration variables or repository-relative paths.
- No code was invented for missing panels. `partial` means that a supplied script produces a verified component of the submitted panel but also contains outputs not used in that panel.
- Generic plotting engines without the submitted panel's exact coordinates or command line are kept only under `Shared_Workflows`; they are not assigned to a specific panel.
- Figure 2b and 2f-2g are intentionally absent because no supplied code reproduces those submitted panels. The target-gene connectivity GO scripts are Extended Data Figure 3 analyses and are not Figure 2f-2g code.
- Figure 3b contains only the JSON metric-extraction script because the submitted Pearson-r panel was drawn directly from its summary table. The observed-versus-predicted correlation script is assigned to Extended Data Figure 4b-4c.
- Figure 5a-5c are intentionally absent: no supplied script reproduced their Circos distribution, functional-annotation enrichment or independent-dataset concordance plots.
- Figure 5d is mapped only to the AS-versus-ChromBPNet-significant-AS enrichment workflow.
- All scripts, comments, messages, panel READMEs and audit tables in this archive use English only.

## Reproducibility notes

Nature Portfolio requires a Code availability statement for central custom code and expects sufficient access for evaluation and reuse. For publication, archive the final repository with a DOI, record exact software versions, and provide the data accessions required to regenerate each panel.

Official guidance:

- https://www.nature.com/nature-portfolio/editorial-policies/reporting-standards
- https://www.nature.com/documents/nr-editorial-policy-checklist-Apr-2023-flat.pdf

## Files

- `FIGURE_CODE_MAP.tsv`: figure/panel-to-script mapping and provenance.
- `UNMATCHED_PANELS.tsv`: panels intentionally omitted because no defensible supplied-code match was found.
- `MAPPING_CORRECTIONS.tsv`: material corrections made after checking the manuscript legends and submitted figures.
- `SOURCE_AUDIT.tsv`: disposition of each newly supplied source file.
- `NEW_SOURCE_BLOCK_AUDIT.tsv`: selected and excluded blocks from the enhancer-sharing notebook.
- `VALIDATION_REPORT.tsv`: static syntax-check results.
- `environment/`: unpinned dependency inventory inferred from the supplied code.
- `figures/`: only figures/panels with uploaded code support.
