# SheepEpimap

Analysis and figure-reproduction code for **“A multi-tissue epigenomic atlas
links the sheep non-coding genome to domestication and complex traits.”**

SheepEpimap integrates 516 RNA-seq, ATAC-seq and CUT&Tag datasets across 43
adult sheep tissues. This repository organizes the supplied code for raw-data
processing, epigenome imputation, chromatin-state annotation, enhancer-gene
linking, sequence-to-function modelling, regulatory-network construction,
allele-specific analyses, selection signatures, trait interpretation and
cross-species regulatory conservation.

> **Scope.** This is a curated research-code archive, not a turnkey workflow.
> The original analyses were run on Linux/HPC systems and require datasets,
> reference files and specialist software that are not distributed here.
> Figure claims are conservative: code is linked to a panel only when the
> supplied scripts and final figure content support that relationship.

## Repository structure

| Path | Contents |
|---|---|
| `pipelines/` | Stage-oriented upstream analysis reconstructed from the supplied SheepEpimap notebook |
| `figures/` | Figure- and panel-oriented workflows with numbered execution order |
| `config/` | Example path configuration for adapting the original HPC layout |
| `environment/` | Python, R and command-line dependency inventories |
| `data/` | Data-access and local-data organization guidance; no research data are included |
| `docs/` | Figure coverage, caption reconciliation, path guidance and release notes |
| `archive/` | English-translated source notebook retained for provenance |
| `tools/` | Repository extraction and static-validation utilities |
| `SOURCE_CODE_MAP.tsv` | Mapping from stage-oriented pipeline files to source notebook blocks |
| `FIGURE_CODE_MAP.tsv` | Detailed mapping from figure panels to scripts and their relationship |

## Figure-code coverage

The supplied figure set contains 34 assembled PDFs: 8 main figures, 8 Extended
Data figures and 18 Supplementary figures. The PDFs are not committed to this
code repository. They were used only to verify figure identity and panel
content.

Use these files before running or citing a figure workflow:

- [`docs/FIGURE_COVERAGE.md`](docs/FIGURE_COVERAGE.md) summarizes all 34 figures.
- [`FIGURE_CODE_MAP.tsv`](FIGURE_CODE_MAP.tsv) lists each mapped script and
  labels it as panel code, panel input, partial panel code, validation or shared
  analysis support.
- [`docs/UNMATCHED_PANELS.tsv`](docs/UNMATCHED_PANELS.tsv) records panels for
  which no exact supplied code was found.
- [`docs/CAPTION_RECONCILIATION.md`](docs/CAPTION_RECONCILIATION.md) documents
  numbering conflicts among the manuscript, figure-caption document and final
  RGB PDFs.

Final multipanel assembly, typography and labels were generally completed
outside these scripts. A directory name therefore means that the contained
code supports the named panel; it does not imply that one command rebuilds the
publication-ready composite figure.

## Analysis workflows

The stage-oriented code is preserved in source order:

1. `pipelines/1.Data preprocessing/` - ATAC-seq, CUT&Tag, peak calling, QC,
   signal matrices, TSS profiles and PCA.
2. `pipelines/2.ChromImpute_ChromHMM/` - ChromImpute and ChromHMM analyses.
3. `pipelines/3.TSR/` - tissue-specific regulatory regions, enhancer-gene
   links, enrichment, motif scans and tissue-specificity scores.
4. `pipelines/4.Network/` - regulatory-network construction and statistics.
5. `pipelines/5.Selection signatures/` - selection scans, GAT enrichment and
   regional FST analysis.

See [`docs/ANALYSIS_WORKFLOWS.md`](docs/ANALYSIS_WORKFLOWS.md) for how these
upstream analyses relate to the figure-oriented code.

## Getting started

```bash
git clone https://github.com/SheepEpimap/sheep.git
cd sheep

mamba env create -f environment/environment.yml
conda activate sheep-epimap

cp config/paths.example.env config/paths.env
python tools/validate_repository.py .
```

Then read the README in the relevant pipeline or panel directory and adapt all
input paths, reference assemblies, scheduler directives and resource requests.
The example path file documents the required concepts, but legacy scripts do
not automatically load it.

## Data and software requirements

No FASTQ/BAM/VCF files, reference genomes, controlled-access cohorts or large
generated outputs are included. Before running a workflow, record the data
repository and accession, checksum, reference assembly, annotation release,
liftOver chain version and local path mapping. See [`data/README.md`](data/README.md).

The environment files are broad dependency inventories because exact package
versions were not present in every supplied source. ChromImpute, ChromHMM,
ChromBPNet, HOMER, MEME Suite, Juicer, Bismark, UCSC utilities and cluster
schedulers may require separate installation. The historical GAT workflow
requires a legacy Python 2.7 environment; see
[`environment/gat_legacy_install.txt`](environment/gat_legacy_install.txt).

## Validation and reproducibility status

Run the repository validator with:

```bash
python tools/validate_repository.py .
```

It checks Python syntax, notebook JSON, Bash syntax when Bash is available,
R/Snakemake syntax when those executables are installed, source and figure map
integrity, Markdown links, UTF-8 decoding, GitHub file-size limits and common
credential patterns. GitHub Actions runs the same repository-level checks on
push and pull requests.

Static checks do not demonstrate biological reproducibility. End-to-end
execution was not possible without the original data, references, HPC software
and cluster environment. Preserve this distinction when describing the code in
the manuscript or a public release.

## Provenance

The repository merges two supplied code collections:

- the stage-oriented `SheepEpimap.zip` project reconstructed from the author
  notebook; and
- the panel-oriented `SheepEpimap_figure_code` archive, including the updated
  Figure 2c methylation and Figure 2j-2k Hi-C workflows.

Source and mapping audits are retained under [`docs/provenance/`](docs/provenance/).
Code was not invented for missing panels. Generic plotting engines remain under
`figures/Shared_Workflows/` unless exact panel parameters were supplied.

## Citation and license

The manuscript citation, DOI and software license were not final in the supplied
materials. Before public release, the code owners should select a license and
add a complete `CITATION.cff` using the final published author list and DOI.
Until a license is added, reuse rights are not granted automatically.

See [`docs/RELEASE_CHECKLIST.md`](docs/RELEASE_CHECKLIST.md) before creating a
public release or archival DOI.
