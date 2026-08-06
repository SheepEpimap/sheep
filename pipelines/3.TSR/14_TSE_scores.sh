#Remove rows with a maximum value below 0.1
####/vol2/mengzhu/snakemake_sheep/expressiondir/Average####
awk -F'\t' '
  NR==1 {print; next}
  {
    max = 0
    for (i=7; i<=NF; i++) {
      v = $i + 0
      if (v > max) max = v
    }
    if (max >= 0.1) print
  }' all_tisssues_expression_aveage_tpm.csv > filtered.tsv
cut -f1,7-49 filtered.tsv > expr.tsv

# SPM: one 0-1 score per tissue
vim 02_tspex_spm.sh
#!/bin/bash
tspex --log expr.tsv tspex_spm.tsv spm
#sbatch -p low -c 4 --mem=8G  02_tspex_spm.sh

# Tau: one 0-1 score per gene
vim 03_tspex_tau.sh
#!/bin/bash
tspex expr.tsv tspex_tau.tsv tau
#sbatch -p low -c 4 --mem=8G 03_tspex_tau.sh

# Tsi: one 0-1 score per gene
cat 04_tspex_tsi.sh
#!/bin/bash
tspex --log expr.tsv tspex_tsi.tsv tsi
#sbatch -p low -c 4 --mem=8G 04_tspex_tsi.sh

