# Contributing

## Scope

Changes should preserve the numbered workflow, centralized configuration, assay-specific behavior, and reproducibility of the analysis.

## Development requirements

1. Create a dedicated branch for each change.
2. Keep external absolute paths in `scripts/00_config.sh` only.
3. Do not commit raw FASTQ, BAM, VCF, HDF5, index, log, temporary, or generated result files.
4. Document new inputs, outputs, parameters, and dependencies in `README.md`.
5. Update the relevant Conda environment file when adding software.
6. Preserve the following default peak rules unless a documented analysis change is intended:
   - ATAC-seq: `q = 0.01`
   - other peak-based assays: `q = 0.05`
   - RNA-seq: no peak calling
7. Run the static validation before committing:

```bash
bash scripts/22_self_test.sh
```

8. Run a representative end-to-end sample when changing scientific logic.
9. Record behavior-changing modifications in the release notes and create a new versioned release.

## Pull requests

A pull request should state:

- the scientific or engineering reason for the change;
- affected scripts and workflow stages;
- validation performed;
- expected changes to outputs;
- compatibility implications for previous releases.
