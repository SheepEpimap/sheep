# Public-release checklist

- [ ] Choose and add a software license approved by the authors/institution.
- [ ] Add manuscript citation, DOI/preprint and author list when final.
- [ ] Replace all placeholder values in a private `config/paths.server.tsv`.
- [ ] Confirm that no controlled GWAS summary statistics or large data are staged for Git.
- [ ] Record checksums and versions for liftOver, chain files, bedtools, LDSC and reference LD resources.
- [ ] Run `python tests/test_static.py`.
- [ ] Run `workflow.py --stage preflight` on the target compute system.
- [ ] Perform a small authorized test before a full Slurm submission.
- [ ] Compare regenerated tables and panels with the archived results.
