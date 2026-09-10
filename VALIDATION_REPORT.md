# Validation report

Validation date: 2026-09-10

```bash
python tools/validate_repository.py .
```

| Check | Status | Result |
|---|---|---|
| Python AST | Pass | 93 Python files parsed |
| Notebook JSON | Pass | 1 notebook parsed |
| Bash syntax | Pass | 151 `.sh` and `.sbatch` files passed `bash -n` |
| R syntax | Not run | The detected local R installation is unusable; GitHub Actions performs the remaining repository checks |
| Snakemake syntax | Not run | Snakemake is unavailable in the local validation environment |
| Figure map | Pass | All 267 mapped script paths exist and use approved relationship labels |
| Markdown links | Pass | All 56 repository-relative links resolved |
| UTF-8 text | Pass | All 463 checked text files decoded |
| GitHub file-size limit | Pass | 533 files checked; none is 100 MiB or larger |
| Credential-pattern scan | Pass | No high-confidence token, access-key or private-key pattern found |
| Enhancer-target gene demo | Pass | Two expected H3K27ac-expression correlations reproduced in 1.3 seconds |

Legacy workflows retain environment-specific HPC path prefixes. Create a
configured working copy with `tools/relocate_paths.py` before running them.
Static checks and the synthetic demo do not establish full-data biological
reproducibility; that requires the original datasets, references and specialist
bioinformatics software.
