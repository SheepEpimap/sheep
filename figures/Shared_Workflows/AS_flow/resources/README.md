# Resources

This directory stores small metadata files and generated reference resources used by the workflow.

## Required user-maintained file

`sample_to_individual.txt` maps each sequencing sample to the individual represented in the phased genotype data. Use two whitespace-delimited columns without a header:

```text
SAMPLE_ID    INDIVIDUAL_ID
```

## Generated resources

The workflow may create the following items in this directory:

- chromosome-size and chromosome-list files;
- GenMap mappability BED files;
- Beagle genetic-map files.

Large raw sequencing data, reference genomes, external indices, and genotype VCF files should remain outside the repository and be referenced through `scripts/00_config.sh`.
