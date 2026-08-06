# Contributing

Contributions should preserve scientific provenance and avoid overstating
figure reproducibility.

1. Open an issue describing the affected workflow, input data and figure panel.
2. Keep upstream analysis code under `pipelines/` and exact panel workflows
   under `figures/`.
3. Add or update the panel README, `FIGURE_CODE_MAP.tsv` and figure coverage
   documentation whenever a panel mapping changes.
4. Do not commit controlled data, credentials, private URLs, reference genomes
   or large generated outputs.
5. Run `python tools/validate_repository.py .` before submitting changes.
6. State whether the change was tested end to end or only checked statically.

New panel assignments must include evidence from the final figure content and
the exact command or input schema required to reproduce it.
