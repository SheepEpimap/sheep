# ADR-001: Combine stage-oriented pipelines and panel-oriented figure code

## Status

Superseded by [ADR-002](0002-group-pipelines-by-author-headings.md)

## Date

2026-08-04

## Context

The supplied material consisted of a single 3,331-line Markdown notebook with
mixed Bash, R, Python and Snakemake blocks, plus a separately curated figure-code
directory organized by manuscript panel. A public repository needs both an
analysis-stage view for reproducibility and a panel view for manuscript review.

## Decision

Keep two complementary top-level trees:

- `pipelines/` contains executable units ordered by analytical stage;
- `figures/` preserves the existing figure/panel mapping.

Retain the original notebook under `archive/` and record every extracted unit in
`SOURCE_CODE_MAP.tsv` so cleanup does not erase provenance.

## Alternatives considered

### Keep only the original Markdown

Rejected because mixed code blocks do not expose reliable entry points and
contain several independent scripts in single fences.

### Reorganize everything only by manuscript panel

Rejected because preprocessing, imputation and regulatory-region construction
are shared upstream workflows that do not map cleanly to one panel.

### Merge all figure scripts into the stage-oriented pipelines

Rejected because that would discard the previously reviewed panel provenance
and make manuscript-code assessment harder.

## Consequences

- Some logic is intentionally represented in both stage and figure contexts.
- Users must choose the relevant entry point rather than run the repository as a
  single monolithic workflow.
- Source provenance remains auditable even after script-boundary cleanup.
- Future releases can progressively parameterize stages without changing the
  panel-oriented figure archive.
