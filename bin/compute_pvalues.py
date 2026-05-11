#!/usr/bin/env python3
"""Compute Wilcoxon rank-sum p-values for all perturbations in a filtered h5ad file.

For each perturbation, tests genes where |log2FC| > 0.25 (treated vs non-targeting).
Applies BH FDR correction per perturbation.

Usage: compute_pvalues.py <input.h5ad>
Output: pvalues_{cell_line}.parquet
"""
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gc

MIN_CELLS = 5
LOG2FC_THRESHOLD = 0.25
EPS = 1e-6


def main():
    input_file = sys.argv[1]
    print(f"Loading: {input_file}")

    adata = ad.read_h5ad(input_file)
    n_cells, n_genes = adata.shape

    cell_line = str(adata.obs['cell_line'].iloc[0])
    print(f"Cell line: {cell_line}")
    print(f"Shape: {n_cells:,} cells x {n_genes:,} genes")

    X = adata.X
    if issparse(X):
        X = X.toarray()
    X = np.asarray(X, dtype=np.float32)

    gene_names = adata.var_names.tolist()
    perts = adata.obs['gene'].values

    control_mask = perts == 'non-targeting'
    n_control = int(control_mask.sum())
    X_control = X[control_mask]
    mean_control = X_control.mean(axis=0)
    print(f"Control cells: {n_control:,}")

    unique_perts = [p for p in np.unique(perts) if p != 'non-targeting']
    print(f"Perturbations: {len(unique_perts)}")

    results = []

    for idx, pert in enumerate(unique_perts):
        pert_mask = perts == pert
        n_treated = int(pert_mask.sum())

        if n_treated < MIN_CELLS:
            continue

        X_treated = X[pert_mask]
        mean_treated = X_treated.mean(axis=0)
        log2fc = np.log2(mean_treated + EPS) - np.log2(mean_control + EPS)

        sig_indices = np.where(np.abs(log2fc) > LOG2FC_THRESHOLD)[0]

        if len(sig_indices) == 0:
            continue

        pvals = np.ones(len(sig_indices))
        for j, gene_idx in enumerate(sig_indices):
            treated_vals = X_treated[:, gene_idx]
            control_vals = X_control[:, gene_idx]

            if np.std(treated_vals) == 0 and np.std(control_vals) == 0:
                continue
            _, pval = stats.ranksums(treated_vals, control_vals)
            pvals[j] = pval

        _, fdr_vals, _, _ = multipletests(pvals, method='fdr_bh')

        for j, gene_idx in enumerate(sig_indices):
            results.append((
                cell_line,
                pert,
                gene_names[gene_idx],
                float(log2fc[gene_idx]),
                float(pvals[j]),
                float(fdr_vals[j])
            ))

        if (idx + 1) % 100 == 0:
            print(f"  {idx+1}/{len(unique_perts)} perturbations ({len(results):,} rows)")

    del X, X_control, adata
    gc.collect()

    df = pd.DataFrame(results, columns=[
        'cell_line', 'target_gene', 'readout_gene', 'log2fc_observed', 'pval', 'fdr'
    ])

    output_file = f"pvalues_{cell_line}.parquet"
    df.to_parquet(output_file, index=False)
    print(f"Written: {output_file} ({len(df):,} rows)")


if __name__ == '__main__':
    main()
