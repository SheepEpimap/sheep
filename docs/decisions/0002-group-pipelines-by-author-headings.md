# ADR-002: Group pipeline code by the author's Markdown headings

## Status

Superseded by [ADR-004](0004-numbered-native-pipeline-code.md)

## Date

2026-08-04

## Context

ADR-001 extracted the original notebook into many small executable units. That
made individual languages and commands easier to inspect, but it separated code
that the author intended to read as one analysis section. The desired public
layout should mirror the visible outline of `SheepEpimap.md` and keep a section's
explanation, code blocks and lower-level headings together.

## Decision

- Each H1 heading becomes one directory under `pipelines/`.
- Each H2 heading becomes one Markdown code file containing the complete source
  slice through the next H2 or H1, including all H3/H4 content and code blocks.
- When an H1 has no H2, its content remains in one H1-named file. If that section
  contains H3 headings, the H1 preamble and each H3 section receive their own
  files because those are the next visible units in the author's outline.
- Original code fences are preserved exactly and in source order. A provenance
  comment and `SOURCE_CODE_MAP.tsv` retain the original line ranges.
- Characters that are unsafe in cross-platform filenames are replaced with an
  underscore. Therefore `2.ChromImpute/ChromHMM` is stored as
  `2.ChromImpute_ChromHMM` while the original heading remains inside the file.
- The previous split pipeline tree and its source map are retained in a local
  ZIP backup outside the public repository.

## Alternatives considered

### Keep the fine-grained executable split

Rejected because it fragments the author's analysis sections and makes the
repository harder to navigate using the original outline.

### Keep only the monolithic notebook

Rejected because a 3,331-line file does not provide direct folder-level access
to the main analytical stages.

### Organize only by programming language

Rejected because language is an implementation detail; it does not reflect the
scientific workflow or the author's titles.

## Consequences

- The pipeline tree now matches the author's conceptual outline and keeps each
  titled passage intact.
- Markdown code files are provenance-oriented records, not directly executable
  scripts; users must adapt and extract blocks before running them.
- Automated validation compares all 39 grouped code blocks against the archived
  notebook to detect omissions, edits or reordering.
- Future executable workflows should be added separately rather than fragmenting
  these source-preserving Markdown files again.
