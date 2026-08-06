# ADR-004: Export numbered native pipeline code files

## Status

Accepted

## Date

2026-08-05

## Context

ADR-002 grouped source code into one Markdown file per heading. This preserved
the author's outline but did not provide ordered, directly recognizable code
files: filenames lacked sequence numbers, contained spaces, and hid Shell,
Python, R and Snakemake code inside Markdown fences. Plotting code was therefore
not visible as native `.R` or `.py` files beside its analytical stage.

The source contains 39 fenced blocks. Thirty-seven contain code; the final two
contain only a single `#` placeholder. Several fence language labels are wrong,
including Snakemake blocks labeled `shell`, a Shell block labeled `R`, and an R
plot labeled `python`.

## Decision

- Retain the five H1-based pipeline directories established by ADR-002.
- Export every non-placeholder fenced block intact into one native code file.
- Number code files independently within each H1 directory using zero-padded
  source order: `01_`, `02_`, and so on.
- Replace spaces in code filenames with underscores and use concise descriptive
  English names.
- Infer the extension from the code content, producing 5 `.smk`, 29 `.sh`,
  2 `.R`, and 1 `.py` file instead of trusting incorrect fence labels.
- Store plotting files in their corresponding pipeline directories rather than
  in Markdown containers.
- Skip source blocks 38 and 39 because they contain no commands, only `#`.
- Record every exported file, source block, source and fence line range,
  heading, and inferred language in `SOURCE_CODE_MAP.tsv`.
- Keep the prior Markdown pipeline tree and source map in a local ZIP backup
  outside the public repository.

## Alternatives considered

### Prefix the existing Markdown files only

Rejected because numbering alone would not expose executable file types or make
plotting scripts directly accessible.

### Split each fenced block into multiple commands or scripts

Rejected because many blocks combine setup, execution and post-processing in a
single author-defined unit. Further splitting would repeat the fragmentation
problem that ADR-002 was intended to solve.

### Trust the Markdown fence language labels

Rejected because several labels conflict with the actual syntax and would
produce misleading or invalid extensions.

## Consequences

- The `pipelines/` tree contains 37 clearly ordered native code files with no
  spaces in code filenames.
- Source code within each fence stays together and remains traceable to the
  translated archive and original line ranges.
- Users can immediately distinguish Snakemake, Shell, R and Python files, but
  the historical paths and data dependencies still require configuration.
- R and Snakemake syntax cannot be verified in the current Windows workspace
  because `Rscript` and `snakemake` are unavailable.
- The English archive and exported HOMER script comment out one unmatched source
  `done` so static Bash parsing succeeds; the original Chinese source is not
  modified.
