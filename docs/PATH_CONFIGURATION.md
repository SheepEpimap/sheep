# Path and cluster configuration

The code was developed across multiple Linux HPC storage roots. The most common
roots are:

| Original root | Typical role |
|---|---|
| `/vol2/mengzhu/` | Primary SheepEpimap analysis, references and software |
| `/vol2/zhangshiwen/` | ATAC-seq inputs and enhancer-gene correlation outputs |
| `/public/home/mengzhu/` | Reference index location |
| `/storage/public/home/.../` | Population-genetics source data |
| `/data/home/.../` | ChromBPNet and figure-workflow compute paths |

Before execution:

1. Copy `config/paths.example.env` to `config/paths.env` and fill in the local
   equivalents.
2. Search the target stage for absolute paths and map each one to the documented
   variable or a repository-relative path.
3. Update scheduler commands (`sbatch` or `jsub`) for the available cluster.
4. Confirm chromosome naming (`1` versus `chr1`), genome assembly, chromosome
   sizes and chain-file versions before using interval operations or liftOver.
5. Run syntax checks and a Snakemake dry run before submitting large jobs.

Useful audit command:

```bash
grep -RInE '/vol2/|/public/home/|/storage/public/home/|/data/home/' pipelines figures
```

`config/paths.example.env` is currently a migration aid, not a global runtime
configuration file. Individual legacy scripts must be refactored or edited to
consume those values.
