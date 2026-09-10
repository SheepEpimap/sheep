# Analysis pipelines

The numbered directories contain the main SheepEpimap analysis workflows:

- `1.Data preprocessing`: sequencing-data preprocessing and quality control
- `2.ChromImpute_ChromHMM`: chromatin-state imputation and annotation
- `3.TSR`: tissue-specific regulatory-region analyses
- `4.Network`: regulatory-network analyses
- `5.Selection signatures`: population-genetic and selection analyses

Scripts are organized in approximate execution order. The preprocessing
workflows use the versioned environments under `Envs/`. Start from the server
configuration templates in `config/`, replace the example reference paths and
sample locations, and then run Snakemake with `--use-conda`.

The archived production configuration was recovered from
`/vol2/mengzhu/snakemake_sheep/config2.yaml`; its portable counterpart is
`config/config2.example.yaml`. Slurm resources are recorded in
`config/cluster.yaml`.
