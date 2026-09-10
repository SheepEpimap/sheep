# Snakemake configuration

`config2.server.yaml` and `cluster.yaml` are the configurations recovered
from the production server at `/vol2/mengzhu/snakemake_sheep/`. They preserve
the 43-tissue RNA-seq/CUT&Tag/ATAC workflow settings and Slurm resource
profiles used for the SheepEpimap analysis.

Before running on another system, copy `config2.server.yaml` to a local config
file and replace the genome, annotation, and chromosome-size paths. Run, for
example:

```bash
snakemake --snakefile "pipelines/1.Data preprocessing/03_call_peak.smk" \
  --configfile pipelines/config/config2.server.yaml \
  --use-conda --cores 12 --cluster-config pipelines/config/cluster.yaml
```
