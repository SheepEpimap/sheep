#!/usr/bin/env python3
# Extracted and reorganized from user-supplied analysis notes.
# Review PROJECT_ROOT and all input paths before execution.
#egcorr_downsample_1.py
import sys
import pandas as pd
import numpy as np
import scipy.stats
from statsmodels.stats.multitest import multipletests
import multiprocessing
import traceback

print("Arguments received:", sys.argv)

print("Loading gene expression data...")
genes1 = pd.read_csv(sys.argv[1], sep='\t', index_col=0)
genes2 = pd.read_csv(sys.argv[2], sep='\t', index_col=0)

print(f"Genes1 shape: {genes1.shape}, columns: {genes1.columns[:3]}")
print(f"Genes2 shape: {genes2.shape}, columns: {genes2.columns[:3]}")

geneall = pd.concat([genes1, genes2], axis=1, join='inner')
print(f"Merged gene expression shape: {geneall.shape}")

expected_order = [
    "abomasum_39", "abomasum_40",
    "adipose_39", "adipose_40",
    "bone-marrow_39", "bone-marrow_40",
    "brainstem_39", "brainstem_40",
    "cecum_39", "cecum_40",
    "cerebellum_39", "cerebellum_40",
    "cerebral-cortex_39", "cerebral-cortex_40",
    "cervix_39", "cervix_40",
    "colon_39", "colon_40",
    "cornua-uteri_39", "cornua-uteri_40",
    "corpus-uteri_39", "corpus-uteri_40",
    "duodenum_39", "duodenum_40",
    "epididymis_39", "epididymis_40",
    "heart_39", "heart_40",
    "hippocampus_39", "hippocampus_40",
    "hypothalamus_39", "hypothalamus_40",
    "ileum_39", "ileum_40",
    "jejunum_39", "jejunum_40",
    "kidney_39", "kidney_40",
    "liver_39", "liver_40",
    "lung_39", "lung_40",
    "lymph-node_39", "lymph-node_40",
    "mammary-gland_39", "mammary-gland_40",
    "medulla-oblongata_39", "medulla-oblongata_40",
    "midbrain_39", "midbrain_40",
    "muscle_39", "muscle_40",
    "omasum_39", "omasum_40",
    "optic-chiasm_39", "optic-chiasm_40",
    "ovary_39", "ovary_40",
    "oviduct_39", "oviduct_40",
    "pineal_39", "pineal_40",
    "pituitary_39", "pituitary_40",
    "pons_39", "pons_40",
    "rectum_39", "rectum_40",
    "reticulum_39", "reticulum_40",
    "rumen_39", "rumen_40",
    "skin_39", "skin_40",
    "soft-horn_39", "soft-horn_40",
    "spleen_39", "spleen_40",
    "splenium_39", "splenium_40",
    "testis_39", "testis_40",
    "thymus_39", "thymus_40",
    "thyroid_39", "thyroid_40"
]

print("Reordering gene expression columns...")

common_columns = [col for col in expected_order if col in geneall.columns]

if len(common_columns) != len(geneall.columns):
    missing = set(geneall.columns) - set(common_columns)
    print(f"Warning: {len(missing)} columns not in expected order: {list(missing)[:3]}...")
    geneall = geneall[common_columns]
else:
    geneall = geneall[common_columns]

common_samples = geneall.columns.tolist()
print(f"Reordered gene expression columns: {len(common_samples)} samples")
print(f"Sample order (first 4): {common_samples[:4]}")

common_samples = geneall.columns.tolist()
print(f"Gene expression samples ({len(common_samples)}): {common_samples[:5]}...")

print("Loading regulatory elements data...")

signal_columns = ['enhancer'] + common_samples

try:
    regs = pd.read_csv(
        sys.argv[3],
        sep=',',
        skiprows=1,  #  ( )
        header=None,
        names=signal_columns,
        index_col=0,  #  enhancer
        dtype={col: float for col in signal_columns[1:]},  #
        low_memory=False
    )
    print(f"Initial regulatory elements shape: {regs.shape}")

    for col in regs.columns:
        regs[col] = pd.to_numeric(regs[col], errors='coerce')

    regs = regs.dropna(how='all')
    print(f"After dropping all-NaN rows: {regs.shape}")
except Exception as e:
    print(f"Error loading regulatory data: {e}")
    traceback.print_exc()
    sys.exit(1)

print("Filtering genes by dynamic range...")
genes = geneall.copy()
genes = genes[genes.apply(lambda row: row.max() / (row.min() + 0.0001) > 6, axis=1)]
print(f"Filtered gene expression shape: {genes.shape}")

print("Filtering regulatory elements by dynamic range...")
regs = regs.apply(pd.to_numeric, errors='coerce')
regs = regs[regs.apply(lambda row: row.max() / (row.min() + 0.0001) > 6, axis=1)]
print(f"Filtered regulatory elements shape: {regs.shape}")

print("Loading TSS positions...")
try:
    tss = pd.read_csv(
        sys.argv[4],
        sep='\t',
        header=None,
        names=['chr', 'start', 'end', 'gene_id1', 'gene_id2', 'strand'],
        index_col=3  #  gene_id1
    )
    print(f"TSS positions loaded. Shape: {tss.shape}")
    print(f"TSS index sample: {tss.index[:5].tolist()}")
    print(f"Gene expression index sample: {genes.index[:5].tolist()}")
except Exception as e:
    print(f"Error loading TSS data: {e}")
    traceback.print_exc()
    sys.exit(1)

def check_in_window(row, chrom, pos, window=500000):
    try:
        ch, rng = row.name.split(':')
        start, end = map(int, rng.split('-'))

        if ch != chrom:
            return None

        enhancer_midpoint = (start + end) // 2
        distance = abs(pos - enhancer_midpoint)

        if distance <= window:
            return distance
        return None
    except Exception as e:
        print(f"Error processing {row.name}: {e}")
        return None

def find_regulators(gene):
    try:
        gene_name = gene.name

        if gene_name not in tss.index:
            return None

        gene_tss = tss.loc[gene_name]
        chrom = gene_tss['chr']
        pos = int(gene_tss['start'])

        dists = regs.apply(check_in_window, axis=1, args=(chrom, pos))
        candidates = dists[dists.notnull()]
        if candidates.empty:
            return None

        results = []
        for reg, dist in candidates.items():
            try:
                gene_vector = gene[common_samples]
                enhancer_vector = regs.loc[reg, common_samples]

                df = pd.DataFrame({
                    'gene': gene_vector,
                    'enhancer': enhancer_vector
                }).dropna()

                if len(df) < 3:
                    continue

                r, pval = scipy.stats.pearsonr(df['gene'], df['enhancer'])

                results.append({
                    'gene': gene_name,
                    'regulator': reg,
                    'pearson_r': r,
                    'pval': pval,
                    'distance': dist,
                    'n_samples': len(df)
                })
            except Exception as e:
                continue

        if not results:
            return None
        return pd.DataFrame(results)

    except Exception as e:
        return None

def batch_process(subset):
    results = []
    for i in range(len(subset)):
        res = find_regulators(subset.iloc[i])
        if res is not None:
            results.append(res)
    return pd.concat(results, ignore_index=True) if results else None

print("Starting correlation analysis...")
print(f"Genes to process: {genes.shape[0]}, Regulatory elements: {regs.shape[0]}")
try:
    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        df_split = np.array_split(genes, multiprocessing.cpu_count())
        results = pool.map(batch_process, df_split)

    print("Merging results...")
    final = pd.concat([res for res in results if res is not None], ignore_index=True)

    if not final.empty:
        print(f"Found {len(final)} correlations before FDR correction")
        _, final['qval'], _, _ = multipletests(final['pval'], method='fdr_bh')
        print(f"After FDR correction: {len(final)} correlations")
    else:
        print("No significant correlations found.")
        final = pd.DataFrame(columns=['gene', 'regulator', 'pearson_r', 'pval', 'distance', 'n_samples', 'qval'])

    # saveresults
    output_file = sys.argv[6]
    final.to_csv(output_file, sep='\t', index=False)
    print(f"Analysis completed. Results saved to {output_file}")

except Exception as e:
    print(f"Error during correlation analysis: {e}")
    traceback.print_exc()
    sys.exit(1)
