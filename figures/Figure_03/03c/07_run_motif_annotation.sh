#!/bin/bash
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11
python 1.py \
    --compiled-h5 /data/home/sczd644/run/zsw_chrombpnet/compiled_results/modisco_compiled.h5 \
    --out-dir /data/home/sczd644/run/zsw_chrombpnet/final_report \
    --meme-db /data/home/sczd644/run/zsw_chrombpnet/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
    --n-matches 10 \
    --top-n-matches 3 \
    --verbose

#!/bin/bash
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11
    modisco report -i /data/home/sczd644/run/zsw_chrombpnet/compiled_results/modisco_compiled.h5 -o /data/home/sczd644/run/zsw_chrombpnet/final_report -m /data/home/sczd644/run/zsw_chrombpnet/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt

#!/bin/bash
module load cuda/11.8
module load cudnn/8.9.6.50_cuda11
python /data/home/sczd644/run/zsw_chrombpnet/HDMA/code/03-chrombpnet/02-compendium/04b-get_tomtom_matches_new.py --modisco-h5 /data/home/sczd644/run/zsw_chrombpnet/compiled_results/modisco_compiled.h5 \
    --out-dir /data/home/sczd644/run/zsw_chrombpnet/final_report \
    --meme-db /data/home/sczd644/run/zsw_chrombpnet/JASPAR2024_CORE_vertebrates_non-redundant_pfms_meme.txt \
    --verbose True
