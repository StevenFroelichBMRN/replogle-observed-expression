#!/usr/bin/env python3
"""
Stream-convert a large Matrix Market (MTX) to 10x HDF5 (.h5) for Seurat.

ASSUMPTION: MTX is CELLS x GENES (i.e., rows = cells, columns = genes).
We produce a CSC matrix where columns = cells (10x convention) and rows = genes.

Two passes:
  Pass 1: count per-cell nnz from MTX rows -> build CSC indptr
  Pass 2: stream MTX again in chunks; write data/indices into HDF5 at proper offsets

Outputs (10x v3-style):
  /matrix/barcodes       (S)
  /matrix/features/id    (S)  # gene IDs (we reuse names here)
  /matrix/features/name  (S)
  /matrix/features/feature_type (S)  # 'Gene Expression'
  /matrix/features/genome (S)
  /matrix/data    (float32)
  /matrix/indices (int32)   # gene row indices
  /matrix/indptr  (int64)   # length = n_cells + 1
  /matrix/shape   (int64)   # [n_genes, n_cells]
"""

import os, sys, gzip, argparse
import numpy as np
import pandas as pd
import h5py
from typing import Tuple

def opengz(path: str):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def read_mtx_header(fp) -> Tuple[int,int,int]:
    """
    Returns (n_rows, n_cols, nnz) directly from the size line.
    NOTE: We treat n_rows as #CELLS and n_cols as #GENES for this script.
    """
    line = fp.readline()
    while line.startswith('%'):
        line = fp.readline()
    n_rows, n_cols, nnz = map(int, line.strip().split())
    return n_rows, n_cols, nnz

def pass1_count_cells_as_rows(mtx_path: str, n_cells_expected: int) -> Tuple[np.ndarray, int, int, int]:
    """
    Count nnz per CELL (MTX row) and build CSC indptr for columns=cells.
    Returns (indptr, n_cells, n_genes, nnz_total).
    """
    with opengz(mtx_path) as f:
        n_cells, n_genes, nnz_declared = read_mtx_header(f)
        if n_cells_expected and n_cells_expected != n_cells:
            print(f"[warn] Barcodes={n_cells_expected}, MTX rows={n_cells}. Proceeding…", file=sys.stderr)

        counts = np.zeros(n_cells, dtype=np.int64)   # nnz per cell (per column of CSC)
        trip = 0
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            i_cell = int(parts[0]) - 1  # 0-based cell
            if 0 <= i_cell < n_cells:
                counts[i_cell] += 1
            trip += 1
            if (trip % 10_000_000) == 0:
                print(f"[pass1] scanned {trip:,} triples…", file=sys.stderr)

    indptr = np.zeros(n_cells + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    return indptr, n_cells, n_genes, int(counts.sum())

def create_10x_h5(out_h5: str, n_genes: int, n_cells: int, nnz: int,
                  genes: pd.Series, barcodes: pd.Series,
                  feature_type: str = "Gene Expression"):
    """
    Create 10x v3 HDF5 skeleton and write strings & shape. Return open file + datasets.
    """
    if os.path.exists(out_h5):
        os.remove(out_h5)
    h5 = h5py.File(out_h5, 'w')
    grp = h5.create_group('matrix')

    grp.create_dataset('barcodes', data=np.asarray(barcodes.astype(str).values, dtype='S'))
    grp.create_dataset('shape', data=np.array([n_genes, n_cells], dtype='int64'))

    data_ds    = grp.create_dataset('data',    shape=(nnz,),      dtype='float32', chunks=True)
    indices_ds = grp.create_dataset('indices', shape=(nnz,),      dtype='int32',  chunks=True)
    indptr_ds  = grp.create_dataset('indptr',  shape=(n_cells+1,),dtype='int64',  chunks=True)

    fgrp = grp.create_group('features')
    gnames = genes.astype(str).values
    fgrp.create_dataset('id',           data=np.asarray(gnames, dtype='S'))
    fgrp.create_dataset('name',         data=np.asarray(gnames, dtype='S'))
    fgrp.create_dataset('feature_type', data=np.asarray([feature_type]*n_genes, dtype='S'))
    fgrp.create_dataset('genome',       data=np.asarray(['']*n_genes, dtype='S'))

    return h5, data_ds, indices_ds, indptr_ds

def pass2_fill_cells_as_rows(mtx_path: str,
                             data_ds, indices_ds, indptr_ds,
                             indptr: np.ndarray,
                             buffer_nnz: int = 10_000_000,
                             sort_within_cell_by_gene: bool = False):
    """
    Fill CSC arrays when MTX rows=cells, cols=genes.
    - "columns" of CSC are cells, so write pointer is per cell.
    - indices are gene indices (rows of CSC).
    """
    n_cells = indptr.size - 1
    write_ptr = indptr.copy()

    buf_cell = np.empty(buffer_nnz, dtype=np.int64)  # which cell (column of CSC)
    buf_gene = np.empty(buffer_nnz, dtype=np.int32)  # which gene (row of CSC)
    buf_val  = np.empty(buffer_nnz, dtype=np.float32)

    def flush(n_used: int):
        if n_used == 0:
            return
        # Sort by cell to cluster writes
        order = np.argsort(buf_cell[:n_used], kind='mergesort')
        cells = buf_cell[:n_used][order]
        genes = buf_gene[:n_used][order]
        vals  = buf_val [:n_used][order]

        if sort_within_cell_by_gene:
            # stable boundaries per cell
            cuts = np.concatenate(([0], np.nonzero(cells[1:] != cells[:-1])[0] + 1, [n_used]))
            for k in range(cuts.size - 1):
                a, b = cuts[k], cuts[k+1]
                sub = np.argsort(genes[a:b], kind='mergesort')
                genes[a:b] = genes[a:b][sub]
                vals [a:b] = vals [a:b][sub]

        # write per cell block
        start = 0
        while start < n_used:
            c = cells[start]
            end = start + 1
            while end < n_used and cells[end] == c:
                end += 1
            L = end - start
            pos = write_ptr[c]
            indices_ds[pos:pos+L] = genes[start:end]
            data_ds[pos:pos+L]    = vals[start:end]
            write_ptr[c] += L
            start = end

    with opengz(mtx_path) as f:
        _ = read_mtx_header(f)
        used = 0
        seen = 0
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            i_cell = int(parts[0]) - 1
            j_gene = int(parts[1]) - 1
            x      = float(parts[2])

            buf_cell[used] = i_cell
            buf_gene[used] = j_gene
            buf_val [used] = x
            used += 1; seen += 1

            if used == buffer_nnz:
                flush(used)
                used = 0
                if (seen % 10_000_000) == 0:
                    print(f"[pass2] streamed {seen:,} triples…", file=sys.stderr)
        flush(used)

    indptr_ds[:] = indptr
    print("[done] pass2 wrote data/indices/indptr.", file=sys.stderr)

def make_unique_series(x: pd.Series) -> pd.Series:
    seen = {}
    out = []
    for v in x.astype(str):
        if v not in seen:
            seen[v] = 0; out.append(v)
        else:
            seen[v] += 1; out.append(f"{v}-dup{seen[v]}")
    return pd.Series(out, index=x.index)

def main():
    ap = argparse.ArgumentParser(description="Convert large CELLS×GENES MTX to 10x HDF5 (CSC) for Seurat.")
    ap.add_argument('--mtx', required=True, help='Path to MTX (.mtx or .mtx.gz), assumed CELLS×GENES')
    ap.add_argument('--genes', required=True, help='CSV with gene names; column = --gene-col')
    ap.add_argument('--barcodes', required=True, help='CSV with cell barcodes; column = --barcode-col')
    ap.add_argument('--gene-col', default='gene_name')
    ap.add_argument('--barcode-col', default='bc_wells')
    ap.add_argument('--out', required=True, help='Output 10x HDF5, e.g., matrix_10x.h5')
    ap.add_argument('--buffer-nnz', type=int, default=10_000_000, help='Triplets buffered per flush (default 10M)')
    ap.add_argument('--sort-within-cell', action='store_true', help='Sort by gene within each cell block (slower)')
    args = ap.parse_args()

    # Read metadata tables
    gdf = pd.read_csv(args.genes, header=0)
    if args.gene_col not in gdf.columns:
        raise SystemExit(f"Column '{args.gene_col}' not found in {args.genes}")
    genes = make_unique_series(gdf[args.gene_col])  # ensure unique IDs/names

    bdf = pd.read_csv(args.barcodes, header=0)
    if args.barcode_col not in bdf.columns:
        raise SystemExit(f"Column '{args.barcode_col}' not found in {args.barcodes}")
    barcodes = bdf[args.barcode_col].astype(str)
    n_cells_expected = len(barcodes)

    print(f"[info] expected cells (barcodes): {n_cells_expected:,}; genes: {len(genes):,}", file=sys.stderr)

    # Pass 1: build CSC indptr from MTX ROWS (cells)
    print("[info] pass1: counting per-cell nnz (rows=cells)…", file=sys.stderr)
    indptr, n_cells_mtx, n_genes_mtx, nnz_total = pass1_count_cells_as_rows(args.mtx, n_cells_expected)

    if len(genes) != n_genes_mtx:
        print(f"[warn] genes.csv ({len(genes)}) != MTX columns ({n_genes_mtx}). Proceeding…", file=sys.stderr)
    if n_cells_expected != n_cells_mtx:
        print(f"[warn] barcodes.csv ({n_cells_expected}) != MTX rows ({n_cells_mtx}). Proceeding…", file=sys.stderr)

    n_cells = n_cells_mtx
    n_genes = len(genes)  # prefer external list as feature catalog

    print(f"[pass1] cells={n_cells:,} genes={n_genes_mtx:,} nnz={nnz_total:,}", file=sys.stderr)
    print("[info] creating 10x HDF5 skeleton…", file=sys.stderr)
    h5, data_ds, indices_ds, indptr_ds = create_10x_h5(args.out, n_genes, n_cells, nnz_total, genes, barcodes)

    # Pass 2: stream-fill CSC arrays
    print("[info] pass2: streaming fill (cells-as-columns CSC)…", file=sys.stderr)
    pass2_fill_cells_as_rows(args.mtx, data_ds, indices_ds, indptr_ds, indptr,
                             buffer_nnz=args.buffer_nnz,
                             sort_within_cell_by_gene=args.sort_within_cell)

    h5.flush(); h5.close()
    print(f"[done] wrote 10x HDF5 -> {args.out}", file=sys.stderr)

if __name__ == '__main__':
    main()