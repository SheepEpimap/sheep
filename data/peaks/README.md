# Unblacklisted peak calls

This directory contains 430 study-generated peak files: 172 `*.broadPeak` and
258 `*narrowPeak`. The files were copied without modification from the analysis
server's `/vol2/mengzhu/snakemake_sheep/peak_unblacklist/` directory.

- `peak_unblacklist/` contains the original, uncompressed peak files.
- `SHA256SUMS` records a checksum for every peak file, relative to this directory.

To verify after cloning:

```bash
cd data/peaks
sha256sum -c SHA256SUMS
```

These are peak calls, not sequencing reads. The originating Snakemake workflow
and its reference-assembly configuration should be checked before combining
these coordinates with other datasets. BAM and BAM index files are not included
in this release.
