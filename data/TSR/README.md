# Tissue-specific regulatory elements (TSRs)

This directory contains tissue-specific regulatory element annotations derived
from the 11-state sheep chromatin-state model.

## Contents

- `AA_TSR_E1/` through `AA_TSR_E11/` contain BED annotations for chromatin
  states E1-E11. Each directory contains 53 tissue or tissue-system files.
- The 583 BED files use four columns: chromosome, start, end, and interval ID.
- `perm_enrich_tsr_gwas.py` contains the permutation-enrichment analysis code.
- `perm_enrich_tsr_gwas.tsv` contains the corresponding analysis results.
- `SHA256SUMS` records checksums for all 585 data and analysis files.

Coordinates should only be combined with other datasets after confirming that
the same sheep reference assembly is being used.

## Checksum verification

From this directory, run:

```bash
sha256sum -c SHA256SUMS
```
