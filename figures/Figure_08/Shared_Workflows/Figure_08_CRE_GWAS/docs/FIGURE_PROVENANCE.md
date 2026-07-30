# Figure 8 provenance

## Reproducible code-derived panels

The canonical module `src/sheep_epimap/figures/fig8.py` and companion notebook generate:

- b: human and mean-of-two-sheep ChromHMM emission probabilities for E1–E9;
- d: merged, deduplicated CRE-class composition;
- e: phyloP raincloud plot and pairwise Welch tests with BH correction;
- f: overlap of positive nominal trait × tissue hits among four CRE classes;
- g: sfCRE enrichment dot plot, with color requiring a positive coefficient and `P < 0.05`.

Panel b was extracted from the emission cells of the final cleaned cross-species notebook. Panels d–g were consolidated from the user-confirmed final visualization notebook and the later cleaned Figure 8 script/notebook.

## Inputs that must be configured

Three ChromHMM emission matrices, the merged classification summary, the verified conservation table, the final all-tissues enrichment table, the 38-trait workbook, the 152-trait metadata table and tissue colors are required. See `manifests/input_manifest.tsv` for the expected objects.

## Final composite boundary

The user-confirmed final composite was `fig8.tif` with SHA256 `5cb22295cfa48050dbcf29bcb499e3bd1148526c898c11319f605e3c034a3ce7`. Its panel-composition software, panels a/c/h, and export parameters were not established by code or an export log. The TIFF is therefore not included and this repository does not claim to recreate the final composite pixel-for-pixel.
