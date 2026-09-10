# Software environment

[`software_versions.tsv`](software_versions.tsv) records the versions reported
for the analyses in the manuscript. The two Trim Galore versions correspond to
different workflows and are intentionally listed separately.

[`environment.yml`](environment.yml) provides a convenient broad Conda
environment for inspecting and adapting the repository. Some large or
workflow-specific programs require separate environments or installations; use
the manuscript-recorded version in `software_versions.tsv` for reproduction.

The lightweight demo requires only Python 3.10 or newer and the Python standard
library. It was tested with Python 3.12.3 on Ubuntu 24.04.3 LTS under WSL2.

To capture a full environment on the analysis server, run:

```bash
conda env export --no-builds > environment/server_environment.yml
conda list --explicit > environment/server_conda_explicit.txt
```

Do not replace the manuscript-recorded versions with newer package versions
without validating the affected workflow.
