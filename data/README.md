# Data and reference resources

Unblacklisted peak calls are distributed in [`peaks/peak_unblacklist/`](peaks/peak_unblacklist/).
See [`peaks/README.md`](peaks/README.md) for provenance and checksums. Raw
sequencing files, reference genomes, population-genetics cohorts and other
generated analysis outputs are not distributed with this repository.

## Study-generated data

- Raw CUT&Tag (H3K4me3, H3K4me1, H3K27ac and H3K27me3), ATAC-seq and RNA-seq
  data across 43 sheep tissues: NCBI SRA BioProject `PRJNA1237432`.
- SheepGTEx RNA-seq data: NCBI SRA BioProject `PRJNA1304012`.
- Processed chromatin states, peaks, bigWig tracks, enhancer-gene links and
  motif-instance annotations: [SheepEpimap UCSC Genome Browser session](https://genome.ucsc.edu/s/mengzhu/SheepEpimap).

## Reference assemblies

- Primary sheep reference used for read alignment: ARS-UI_Ramb_v2.0,
  NCBI accession `GCF_016772045.1`.
- Cross-species CRE workflow intermediate sheep assembly: ARS-UI_Ramb_v3.0,
  NCBI accession `GCF_016772045.2`.
- Human coordinate systems: hg38 and hg19, as identified in the corresponding
  Figure 8 path configuration and liftOver scripts.

The Figure 8 comparison uses the `hg38ToGCF_016772045.2` chain with
`minMatch=0.1`; downstream human LDSC inputs are converted from hg38 to hg19.
Record the checksum and source URL for each downloaded FASTA, annotation and
chain file because locally generated chain files do not have a universal
accession.

## Additional published data

The analyses also use published Hi-C, WGBS, SheepGTEx molecular QTL, OMIA,
Animal QTLdb, ENCODE, Roadmap Epigenomics and human GWAS resources cited in the
manuscript Methods. Access restrictions and redistribution terms remain those
of the source repositories.

The root `.gitignore` excludes common large genomics formats. The redistributable
synthetic example under [`tests/demo/`](../tests/demo/) includes its expected
output and does not contain study participant or controlled-access data.
