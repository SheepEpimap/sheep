# Validation report

Validation date: 2026-08-06

```bash
python tools/validate_repository.py .
```

| Check | Status | Result |
|---|---|---|
| Python AST | Pass | 73 Python files parsed |
| Notebook JSON | Pass | 1 notebook parsed |
| Bash syntax | Pass | 107 `.sh` and `.sbatch` files passed `bash -n` |
| R syntax | Pass | 20 R files parsed with R 4.4.1 |
| Snakemake syntax | Not run | Snakemake is unavailable in the validation environment |
| Figure map | Pass | All mapped script paths exist and use approved relationship labels |
| Markdown links | Pass | All repository-relative links resolved |
| UTF-8 text | Pass | All checked text files decoded |
| GitHub file-size limit | Pass | No file is 100 MiB or larger |
| Credential-pattern scan | Pass | No high-confidence token, access-key or private-key pattern found |

The workflows retain environment-specific HPC paths that must be configured
before execution. Static checks do not establish end-to-end biological
reproducibility; full testing requires the original datasets, references and
bioinformatics software.
