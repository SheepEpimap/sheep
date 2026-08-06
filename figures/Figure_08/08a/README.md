# Figure 8a

Figure 8a is the approximately 43-sheep-tissue by 37-human-trait projected-EnhA
LDSC heatmap.

The two numbered Bash scripts in this directory are Figure 8a analysis code:

1. `Figure_08a_01_liftover_sheep_states_to_hg38.sh` projects sheep
   chromatin-state intervals from Ramb v2 through Ramb v3 to hg38.
2. `Figure_08a_02_convert_hg38_to_hg19.sh` converts the projected intervals to
   hg19 for LDSC/GWAS analysis.

Run the scripts in numeric order after checking the configured paths. The
multi-branch `gwas_submit_slurm` controller remains under
`../Shared_Workflows/Figure_08_CRE_GWAS/` because it is shared by several
Figure 8 GWAS analyses.

The supplied coordinate scripts match `*_E*.bed` and use
`LIFTOVER_MINMATCH=0.8`; preserve these original parameters unless the final
analysis records specify otherwise.
