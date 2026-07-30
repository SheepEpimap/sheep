# Shared Figure 8 workflows

`Figure_08_CRE_GWAS/` is retained as one canonical executable project because
its launchers and common workers serve more than one analysis branch.

The `workflow.py --stage gwas_submit_slurm` command calls
`src/gwas_ldsc/submit_E1E9_reuse_existing_152.sh`, which submits:

- A1: marginal human CRE-class by TSR-group LDSC;
- B1: marginal sheep-projection by TSR-group LDSC;
- A2: human CRE-class by tissue CTS LDSC, used by Figure 8f, Figure 8g and
  Supplementary Figure 16c;
- B2: sheep-projection by tissue CTS LDSC for 10 paired tissues, not the
  approximately 43-tissue analysis shown in Figure 8a.

Keeping the controller here avoids duplicating one multi-branch workflow in
several panel folders. Panel-specific plotting and input-preparation scripts
remain in the numbered panel directories.
