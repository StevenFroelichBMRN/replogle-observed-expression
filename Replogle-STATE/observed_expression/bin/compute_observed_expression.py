#!/usr/bin/env python3
"""Compute pseudobulk observed expression for all perturbations in a filtered h5ad file.

Usage: compute_observed_expression.py <input.h5ad>
Output: observed_expression_{cell_line}.parquet
"""
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import issparse
import gc

MIN_CELLS = 5
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
    mean_control = X[control_mask].mean(axis=0)
    print(f"Control cells: {n_control:,}")

    unique_perts = [p for p in np.unique(perts) if p != 'non-targeting']
    print(f"Perturbations: {len(unique_perts)}")

    results = []

    for idx, pert in enumerate(unique_perts):
        pert_mask = perts == pert
        n_treated = int(pert_mask.sum())

        if n_treated < MIN_CELLS:
            continue

        mean_treated = X[pert_mask].mean(axis=0)
        log2fc = np.log2(mean_treated + EPS) - np.log2(mean_control + EPS)

        for j, gene in enumerate(gene_names):
            results.append((
                cell_line,
                pert,
                gene,
                float(log2fc[j]),
                float(mean_treated[j]),
                float(mean_control[j]),
                n_treated,
                n_control
            ))

        if (idx + 1) % 100 == 0:
            print(f"  {idx+1}/{len(unique_perts)} perturbations ({len(results):,} rows)")

    del X, adata
    gc.collect()

    df = pd.DataFrame(results, columns=[
        'cell_line', 'target_gene', 'readout_gene', 'log2fc_observed',
        'mean_expr_treated', 'mean_expr_control', 'n_cells_treated', 'n_cells_control'
    ])

    output_file = f"observed_expression_{cell_line}.parquet"
    df.to_parquet(output_file, index=False)
    print(f"Written: {output_file} ({len(df):,} rows, {df.memory_usage(deep=True).sum()/1e6:.0f} MB)")

if __name__ == '__main__':
    main()
