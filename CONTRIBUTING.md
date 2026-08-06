# Contributing

Contributions should preserve scientific reproducibility and document the
inputs, software and parameters needed for each workflow.

1. Open an issue describing the affected workflow, input data and figure panel.
2. Keep upstream analysis code under `pipelines/` and exact panel workflows
   under `figures/`.
3. Update the panel README and `FIGURE_CODE_MAP.tsv` whenever a panel workflow
   changes.
4. Do not commit controlled data, credentials, private URLs, reference genomes
   or large generated outputs.
5. Run `python tools/validate_repository.py .` before submitting changes.
6. State whether the change was tested end to end or only checked statically.

New panel workflows should include the command and input schema required to run
them.
