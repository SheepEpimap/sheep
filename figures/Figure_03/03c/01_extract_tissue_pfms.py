#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
"""
 MoDIScoresults PFM  -  outputPFM file
"""

import os
import h5py
import numpy as np
import pandas as pd

def trim_cwm(cwm, trim_threshold=0.3, trim_min_length=3):
    """ CWM"""
    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * trim_threshold
    pass_inds = np.where(score >= trim_thresh)[0]

    if len(pass_inds) == 0 or (pass_inds[-1] - pass_inds[0] + 1) < trim_min_length:
        return None

    imp_start = np.min(pass_inds)
    imp_end = np.max(pass_inds) + 1
    return imp_start, imp_end

def main():
    modisco_dir = "/data/home/sczd644/run/zsw_chrombpnet/modisco_h5"
    output_pfm = "all_tissues_motifs.pfm"
    output_metadata = "all_tissues_motifs_metadata.tsv"

    with open("tissue.txt", "r") as f:
        tissues = [line.strip() for line in f if line.strip()]

    print(f"  {len(tissues)}  ...")

    all_pfms = []
    metadata_records = []

    for tissue in tissues:
        modisco_file = os.path.join(modisco_dir, f"{tissue}.modisco_results.h5")

        if not os.path.exists(modisco_file):
            print(f"  {tissue}: file ")
            continue

        print(f"  {tissue}...")

        try:
            with h5py.File(modisco_file, 'r') as f:
                if 'pos_patterns' not in f:
                    print(f"  {tissue}:  pos_patterns, ")
                    continue

                pattern_count = 0
                for pattern_name in f['pos_patterns'].keys():
                    pattern = f['pos_patterns'][pattern_name]

                    ppm = pattern['sequence'][()]
                    cwm = pattern['contrib_scores'][()]
                    num_seqlets = pattern['seqlets']['n_seqlets'][0]

                    trim_result = trim_cwm(cwm, 0.3, 3)
                    if trim_result:
                        imp_start, imp_end = trim_result
                        trimmed_ppm = ppm[imp_start:imp_end]
                        trimmed_length = trimmed_ppm.shape[0]

                        motif_id = f"{tissue}__pos_patterns.{pattern_name}"

                        pfm = trimmed_ppm * num_seqlets

                        pfm_str = f">{motif_id}\n"
                        for row in pfm:
                            pfm_str += f"{int(round(row[0]))}\t{int(round(row[1]))}\t{int(round(row[2]))}\t{int(round(row[3]))}\n"

                        all_pfms.append(pfm_str)

                        metadata_records.append({
                            'motif_id': motif_id,
                            'tissue': tissue,
                            'pattern_name': pattern_name,
                            'pattern_type': 'pos_patterns',
                            'num_seqlets': int(num_seqlets),
                            'trimmed_length': trimmed_length,
                            'original_length': ppm.shape[0],
                            'trim_start': int(imp_start),
                            'trim_end': int(imp_end)
                        })

                        pattern_count += 1

                print(f"  {tissue}:   {pattern_count}  motif")

        except Exception as e:
            print(f"  {tissue}  : {e}")
            continue

    with open(output_pfm, 'w') as f:
        for pfm in all_pfms:
            f.write(pfm)

    if metadata_records:
        df_metadata = pd.DataFrame(metadata_records)
        df_metadata.to_csv(output_metadata, sep='\t', index=False)
        print(f"\n file save: {output_metadata}")

    print(f"\n !   {len(all_pfms)}  motif")
    print(f"PFMfile: {output_pfm}")
    print(f" file: {output_metadata}")

if __name__ == '__main__':
    main()
