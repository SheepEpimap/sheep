# Public-release checklist

- [ ] Select and add a software license approved by all code owners.
- [ ] Add `CITATION.cff` with the final author list, manuscript citation and DOI.
- [ ] Replace or parameterize personal/HPC absolute paths.
- [ ] Add data repository accessions, checksums and access restrictions.
- [ ] Freeze exact Python, R, Java and command-line tool versions.
- [ ] Confirm reference genome, annotation and liftOver chain versions.
- [ ] Test every claimed entry point on a clean Linux environment.
- [ ] Add a small synthetic or redistributable test dataset where possible.
- [ ] Install Snakemake and lint the 5 `.smk` files. Python, R and Bash static syntax checks pass.
- [ ] Check that no credentials, tokens, private URLs or controlled data are committed.
- [ ] Confirm generated files and large genomics formats are excluded from Git.
- [ ] Review figure labels and captions against the final accepted manuscript.
- [ ] Tag the tested release and archive it in Zenodo or another DOI-granting repository.
