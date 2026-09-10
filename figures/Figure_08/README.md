# Figure 8 code

This directory contains the cross-species CRE and GWAS/LDSC workflows used for
Figure 8. Panel-specific launchers are grouped in the `08a-08g` directories,
and shared components are under `Shared_Workflows`.

- `08a`: sheep ARS-UI_Ramb_v3.0 (`GCF_016772045.2`) to human hg38
  coordinate conversion and projected-EnhA preparation
- `08b-08d`: human hg38 to sheep ARS-UI_Ramb_v3.0
  (`GCF_016772045.2`) conversion (`minMatch=0.1`) followed by E1-E9
  cross-species CRE classification
- `08e`, `08f`, `08g`: panel analysis and plotting workflows
- `08f-08g`: shared input preparation for panels f and g

The Figure 8 analyses use the E1-E9 chromatin-state scope. Configure paths in
the shared workflow before execution.
