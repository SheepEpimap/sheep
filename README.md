# SheepEpimap

Analysis and figure-generation code for **“A multi-tissue epigenomic atlas
links the sheep non-coding genome to domestication and complex traits.”**

SheepEpimap integrates 516 RNA-seq, ATAC-seq and CUT&Tag datasets across 43
adult sheep tissues. The project includes workflows for raw-data processing,
epigenome imputation, chromatin-state annotation, enhancer-gene linking,
sequence-to-function modelling, regulatory-network construction,
allele-specific analyses, selection signatures, trait interpretation and
cross-species regulatory conservation.

## Repository structure

| Path | Contents |
|---|---|
| `pipelines/` | Stage-oriented upstream analysis workflows |
| `figures/` | Figure- and panel-oriented workflows with numbered execution order |
| `config/` | Example path configuration for adapting the original HPC layout |
| `environment/` | Python, R and command-line dependency inventories |
| `data/` | Data-access and local-data organization guidance |
| `docs/` | Workflow, path-configuration and release documentation |
| `tools/` | Repository validation utilities |
| `FIGURE_CODE_MAP.tsv` | Figure-panel to script index |

## Analysis workflows

1. `pipelines/1.Data preprocessing/` - ATAC-seq, CUT&Tag, peak calling, QC,
   signal matrices, TSS profiles and PCA.
2. `pipelines/2.ChromImpute_ChromHMM/` - ChromImpute and ChromHMM analyses.
3. `pipelines/3.TSR/` - tissue-specific regulatory regions, enhancer-gene
   links, enrichment, motif scans and tissue-specificity scores.
4. `pipelines/4.Network/` - regulatory-network construction and statistics.
5. `pipelines/5.Selection signatures/` - selection scans, GAT enrichment and
   regional FST analysis.

See [`docs/ANALYSIS_WORKFLOWS.md`](docs/ANALYSIS_WORKFLOWS.md) for workflow
relationships and recommended execution order.

## Getting started

```bash
git clone https://github.com/SheepEpimap/sheep.git
cd sheep

mamba env create -f environment/environment.yml
conda activate sheep-epimap

cp config/paths.example.env config/paths.env
python tools/validate_repository.py .
```

Read the README in the relevant pipeline or figure directory before execution.
Update all input paths, reference assemblies, scheduler directives and resource
requests for the target environment.

## Figure workflows

Figure code is organized as:

```text
figures/<figure>/<panel>/<numbered script>
```

Run numbered scripts in order. Panel READMEs describe required inputs, outputs
and execution details. [`FIGURE_CODE_MAP.tsv`](FIGURE_CODE_MAP.tsv) provides a
searchable index from figures and panels to scripts.

Final multipanel assembly, labels and typography are generally completed after
the computational workflows generate individual plots and tables.

## Data requirements

No FASTQ/BAM/VCF files, reference genomes, controlled-access cohorts or large
generated outputs are included. Before running a workflow, record the data
repository and accession, checksum, reference assembly, annotation release,
liftOver chain version and local path mapping. See [`data/README.md`](data/README.md).

## Software requirements

The environment files provide broad dependency inventories. ChromImpute,
ChromHMM, ChromBPNet, HOMER, MEME Suite, Juicer, Bismark, UCSC utilities and
cluster schedulers may require separate installation. The historical GAT
workflow requires a legacy Python 2.7 environment; see
[`environment/gat_legacy_install.txt`](environment/gat_legacy_install.txt).

## Validation

Run:

```bash
python tools/validate_repository.py .
```

The validator checks Python, R and Bash syntax when the corresponding tools are
available, notebook JSON, mapped paths, Markdown links, UTF-8 decoding, GitHub
file-size limits and common credential patterns. GitHub Actions runs the same
repository-level checks on pushes and pull requests.

Static validation does not replace end-to-end testing with the original data,
reference files and bioinformatics software.

## Citation and license

Add the final manuscript citation and DOI to `CITATION.cff` when available. The
code owners should select and add an appropriate software license before public
release.

See [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md) before creating a
versioned release or archival DOI.
