# Public-release checklist

- [x] Add an MIT `LICENSE` file.
- [x] Add `CITATION.cff`; update the journal citation and DOI after publication.
- [x] Provide a repository-wide path-prefix relocation tool and example mapping.
- [x] Add study data accessions and access restrictions.
- [x] Record the analysis software versions reported in the manuscript.
- [x] Record the reference assemblies and liftOver chain requirements used in the manuscript.
- [ ] Test every claimed entry point on a clean Linux environment.
- [x] Add a synthetic Figure 2d enhancer-target gene correlation demo with expected output and runtime.
- [ ] Parse and dry-run the 5 `.smk` files with Snakemake on a configured Linux environment.
- [x] Check that no credentials, tokens, private URLs or controlled data are committed.
- [x] Confirm generated files and large genomics formats are excluded from Git.
- [ ] Review figure labels and captions against the final accepted manuscript.
- [ ] Tag the tested release and archive it in Zenodo or another DOI-granting repository.
