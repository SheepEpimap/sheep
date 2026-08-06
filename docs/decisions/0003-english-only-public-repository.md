# ADR-003: Use English throughout the public repository

## Status

Accepted

## Date

2026-08-04

## Context

The source notebook and heading-grouped pipeline files contained Chinese
headings, explanatory prose, code comments and user-facing messages. Four
pipeline filenames were also in Chinese. An English-only public repository is
easier for an international audience to navigate, review and reuse.

## Decision

- Translate all Chinese repository text into technical English, including
  headings, prose, code comments, diagnostic messages and installation notes.
- Rename Chinese pipeline filenames and update all local Markdown links and
  `SOURCE_CODE_MAP.tsv` entries to match.
- Apply the same translations to `archive/SheepEpimap.md` and the grouped
  pipeline copies so provenance line ranges and block-level comparisons remain
  valid.
- Preserve computational commands, variable names, source order, data paths and
  scientific parameters unless a translated user-facing string requires a text
  change.
- Keep the original Chinese working source outside the public repository and
  retain the earlier split pipeline version in the local backup ZIP.

## Alternatives considered

### Translate documentation only

Rejected because reviewers would still encounter Chinese comments, headings and
filenames in the primary workflow record.

### Keep a bilingual repository

Rejected because duplicating every title and comment would make the large code
record harder to scan and maintain.

### Translate filenames but leave code comments unchanged

Rejected because the filenames would become accessible while the implementation
notes and scientific intent would remain inaccessible to many readers.

## Consequences

- The public repository contains no Chinese filenames or Han-script text.
- The archived notebook is an English translation rather than a byte-identical
  copy of the author's Chinese working source.
- Automated validation compares the grouped code blocks with the translated
  archive and continues to verify all 39 blocks in source order.
- Translation can improve accessibility but does not validate the scientific or
  computational correctness of the underlying workflows.
