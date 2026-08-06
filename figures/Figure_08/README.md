# Figure 8 code

Every Python, Bash and notebook file in the supplied archive was inspected.
The complete translated, executable project is under `Shared_Workflows`; direct
panel launchers and confirmed panel-analysis scripts are named by figure in
08a-08g.

- 8a contains the sheep-to-hg38 and hg38-to-hg19 coordinate-conversion steps
  used before the projected-EnhA GWAS/LDSC analysis.
- Direct plotting code exists for 8b, 8d, 8e, 8f and 8g.
- 8c contains the E1-E9 cross-species CRE-classification analysis represented
  by its scheme.
- Shared input preparation for 8f and 8g is stored in the `08f-08g` panel
  directory.
- `gwas_submit_slurm` is retained once in `Shared_Workflows` because it submits
  A1/B1/A2/B2 together; its A2 branch serves 8f, 8g and Supplementary 16c.
- No exact plotting code was found for 8h or 8i.
- E1-E9 is intentionally retained for Figure 8. The E1-E11 scope used by other
  figures must not be substituted here.
