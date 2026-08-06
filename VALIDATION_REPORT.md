# Validation report

Validation date: 2026-08-06

Command:

```bash
python tools/validate_repository.py .
```

## Results

| Check | Status | Result |
|---|---|---|
| Python AST | Pass | 73 Python files parsed |
| Notebook JSON | Pass | 1 notebook parsed |
| Bash syntax | Pass | 107 `.sh` and `.sbatch` files passed `bash -n` with GNU Bash 5.2 |
| R syntax | Pass | 20 R files parsed with R 4.4.1 |
| Snakemake syntax | Not run | Snakemake is unavailable in this workspace; 5 `.smk` files were found |
| Source map | Pass | 37 stage-oriented files exist, remain numbered and map to source blocks |
| Native code blocks | Pass | 37 non-placeholder notebook blocks match the archive; 2 comment-only placeholders are skipped |
| Figure map | Pass | 131 mapped scripts exist and use approved relationship labels |
| Figure coverage | Pass | 34 unique final RGB figure files are accounted for |
| Markdown links | Pass | 63 repository-relative links resolved |
| UTF-8 text | Pass | 322 text files decoded |
| GitHub file-size limit | Pass | 323 files checked; none is 100 MiB or larger |
| Credential-pattern scan | Pass | No high-confidence access-key, token or private-key patterns were found |

## Informational findings

The repository intentionally retains historical HPC paths for provenance. The
audit found `/vol2/`, `/public/home/`, `/storage/public/home/` and `/data/home/`
prefixes in 129 files. These paths are not credentials, but they prevent
portable execution and may expose local storage conventions. Follow
[`docs/PATH_CONFIGURATION.md`](docs/PATH_CONFIGURATION.md) before execution or a
privacy-sensitive public release.

## Scope limitation

Static parsing does not establish scientific reproducibility. Snakemake linting,
tool-version resolution, input-schema verification and end-to-end execution
still require the original datasets and a configured Linux bioinformatics/HPC
environment.
