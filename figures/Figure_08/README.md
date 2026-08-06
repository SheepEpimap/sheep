# Figure 8 code

This directory contains the cross-species CRE and GWAS/LDSC workflows used for
Figure 8. Panel-specific launchers are grouped in the `08a-08g` directories,
and shared components are under `Shared_Workflows`.

- `08a`: sheep-to-human coordinate conversion and projected-EnhA preparation
- `08b`, `08d`, `08e`, `08f`, `08g`: panel analysis and plotting workflows
- `08c`: E1-E9 cross-species CRE-classification analysis
- `08f-08g`: shared input preparation for panels f and g

The Figure 8 analyses use the E1-E9 chromatin-state scope. Configure paths in
the shared workflow before execution.
