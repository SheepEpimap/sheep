# Enhancer-target gene demo

This synthetic end-to-end example exercises the Figure 2d enhancer-gene
linking analysis. It combines two RNA-expression matrices, filters genes and
enhancers with sufficient cross-tissue dynamic range, selects enhancers within
500 kb of each transcription start site, correlates H3K27ac signal with gene
expression, and controls the false-discovery rate with Benjamini-Hochberg
correction.

The six columns represent three tissue pairs and are synthetic; they do not
contain study participant data. Two nearby enhancer-gene pairs are constructed
to have positive correlation, while an enhancer on another chromosome verifies
the genomic-window filter.

## Run

From the repository root, after activating `environment/environment.yml`:

```bash
python tests/demo/run_demo.py
```

The runner invokes the production script at
`figures/Figure_02/02d/08_correlate_enhancers_with_genes.py`, writes to a
temporary directory, compares the result with
`expected_output/enhancer_gene_correlations.tsv`, and prints `PASS` on success.

Typical runtime is under two seconds using one CPU core and less than 250 MB of
memory after dependencies are loaded.
