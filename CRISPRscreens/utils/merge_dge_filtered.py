#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Out-of-core row-wise merge for Parse DGE_filtered matrices (cells = rows, genes = columns).

For each sublibrary, expects:
  <SUB>/MPRA/DGE_filtered/count_matrix.mtx(.gz)
  <SUB>/MPRA/DGE_filtered/all_genes.csv
  <SUB>/MPRA/DGE_filtered/cell_metadata.csv

This script:
  • Verifies that all_genes.csv are identical (line-by-line).
  • Reads only MTX headers to accumulate total rows/cols/nnz.
  • Validates: rows in cell_metadata.csv (excluding header) == MTX rows.
  • Streams each MTX and writes a merged Matrix Market by offsetting the ROW index (i).
  • Streams and concatenates cell metadata, adding a 'sublibrary' column.
  • Optionally gzips the merged MTX.

Usage:
  python3 merge_dge_filtered_rows_streamed.py \
    --sublibs output_sublibrary_01 ... output_sublibrary_16 \
    --outdir output_combined_filtered \
    --gzip_out
"""

import sys, os, gzip, argparse, csv
from pathlib import Path

# ---------- small I/O helpers ----------
def opengz(path, mode='rt'):
    path = str(path)
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)

def find_either(gz, plain):
    gz, plain = Path(gz), Path(plain)
    if gz.exists(): return gz
    if plain.exists(): return plain
    sys.exit(f"[ERROR] Not found: {gz} or {plain}")

# ---------- MatrixMarket helpers ----------
def mtx_sizes_from_header(path):
    """Return (nrows, ncols, nnz) from the MTX header only."""
    with opengz(path, 'rt') as fh:
        first = fh.readline()
        if not first:
            sys.exit(f"[ERROR] Empty MTX: {path}")
        if not first.startswith('%%MatrixMarket'):
            fh.seek(0)
        # sizes line = first non-comment
        for line in fh:
            if line.startswith('%'):
                continue
            parts = line.strip().split()
            if len(parts) != 3:
                sys.exit(f"[ERROR] Bad sizes line in MTX: {path}")
            try:
                r, c, nnz = map(int, parts)
            except Exception:
                sys.exit(f"[ERROR] Non-integer sizes in MTX: {path}")
            return r, c, nnz
    sys.exit(f"[ERROR] No sizes line found in MTX: {path}")

def stream_mtx_data_with_row_offset(in_path, out_fh, row_offset):
    """Skip header+sizes, then stream lines; offset ROW index (i)."""
    with opengz(in_path, 'rt') as fh:
        first = fh.readline()
        if not first:
            return
        if not first.startswith('%%MatrixMarket'):
            fh.seek(0)
        # skip sizes line
        for line in fh:
            if line.startswith('%'):
                continue
            break
        # now data lines: i j v
        for line in fh:
            if not line or line.startswith('%'):
                continue
            parts = line.split()
            if len(parts) != 3:
                continue
            i, j, v = parts
            out_fh.write(f"{int(i) + row_offset} {j} {v}\n")

# ---------- genes: verify equal & copy ----------
def copy_first_genes_and_verify(out_genes_path, first_genes_path, other_gene_paths):
    """Copy first all_genes.csv to output; verify others are identical line-by-line."""
    with open(first_genes_path, 'rt', newline='') as inf, open(out_genes_path, 'wt', newline='') as outf:
        header = inf.readline()
        if not header:
            sys.exit(f"[ERROR] Empty genes file: {first_genes_path}")
        outf.write(header)
        n_genes = 0
        for line in inf:
            outf.write(line)
            n_genes += 1
    for p in other_gene_paths:
        with open(out_genes_path, 'rt', newline='') as ref, open(p, 'rt', newline='') as test:
            if ref.readline() != test.readline():
                sys.exit(f"[ERROR] all_genes header differs in {p}")
            for i, (a, b) in enumerate(zip(ref, test), start=1):
                if a != b:
                    sys.exit(f"[ERROR] all_genes row {i} differs in {p}")
            extra = test.readline()
            if extra != "":
                sys.exit(f"[ERROR] all_genes row count differs (extra lines) in {p}")
    return n_genes

# ---------- metadata helpers ----------
def count_metadata_rows(path):
    """Count CSV rows excluding header; supports .csv and .csv.gz."""
    if str(path).endswith(".gz"):
        with opengz(path, 'rt') as fh:
            _ = fh.readline()
            return sum(1 for _ in fh)
    with open(path, 'rt', newline='') as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample) if sample else csv.get_dialect("excel")
        except Exception:
            dialect = csv.get_dialect("excel")
        reader = csv.reader(fh, dialect=dialect)
        try:
            next(reader)
        except StopIteration:
            return 0
        n = 0
        for _ in reader:
            n += 1
        return n

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Row-wise merge of DGE_filtered across sublibraries (cells=rows).")
    ap.add_argument("--sublibs", nargs='+', required=True, help="Sublibrary root folders")
    ap.add_argument("--outdir", required=True, help="Output directory for merged DGE")
    ap.add_argument("--rel_mtx",   default="MPRA/DGE_filtered/count_matrix.mtx",
                    help="Relative path to matrix.mtx (or .mtx.gz) under each sublibrary")
    ap.add_argument("--rel_genes", default="MPRA/DGE_filtered/all_genes.csv",
                    help="Relative path to all_genes.csv")
    ap.add_argument("--rel_cells", default="MPRA/DGE_filtered/cell_metadata.csv",
                    help="Relative path to cell_metadata.csv")
    ap.add_argument("--sublib_label", default=None,
                    help="Comma-separated labels; defaults to directory names")
    ap.add_argument("--gzip_out", action="store_true", help="Write merged MTX as count_matrix.mtx.gz")
    args = ap.parse_args()

    sublibs = [Path(s).resolve() for s in args.sublibs]
    labels = args.sublib_label.split(",") if args.sublib_label else [p.name for p in sublibs]
    if len(labels) != len(sublibs):
        sys.exit("[ERROR] --sublib_label count must match --sublibs")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Resolve inputs (support gz/plain MTX)
    mtx_files   = [find_either(s / (args.rel_mtx + ".gz"),  s / args.rel_mtx) for s in sublibs]
    genes_files = [s / args.rel_genes for s in sublibs]
    cells_files = [s / args.rel_cells for s in sublibs]

    # 1) Copy first genes and verify equality
    out_genes = outdir / "all_genes.csv"
    n_genes = copy_first_genes_and_verify(out_genes, genes_files[0], genes_files[1:])

    # 2) First pass: read MTX sizes, verify metadata rows == MTX rows
    total_rows, total_cols, total_nnz = 0, None, 0
    per_rows = []  # MTX rows per sublib (cells)
    for mtx_path, meta_path in zip(mtx_files, cells_files):
        r, c, nnz = mtx_sizes_from_header(mtx_path)
        meta_rows = count_metadata_rows(meta_path)
        if meta_rows != r:
            sys.exit(f"[ERROR] cell_metadata rows ({meta_rows}) != MTX rows ({r}) in {meta_path}")
        per_rows.append(r)
        total_rows += r
        total_nnz  += nnz
        if total_cols is None:
            total_cols = c
        elif total_cols != c:
            sys.exit(f"[ERROR] Gene (column) count differs across sublibraries: saw {total_cols} vs {c}.")

    # Additional sanity: total_cols should match n_genes from all_genes.csv
    if total_cols != n_genes:
        sys.exit(f"[ERROR] Gene count mismatch: MTX columns={total_cols} vs all_genes rows={n_genes}")

    # 3) Write merged metadata (streamed, prepend sublibrary)
    # Reuse header from first metadata; ensure all headers match
    with open(cells_files[0], 'rt', newline='') as fh0:
        hdr0 = fh0.readline().rstrip('\n')
    for p in cells_files[1:]:
        with open(p, 'rt', newline='') as fh:
            if fh.readline().rstrip('\n') != hdr0:
                sys.exit(f"[ERROR] cell_metadata header differs in {p}")

    out_meta = outdir / "cell_metadata.csv"
    with open(out_meta, 'wt', newline='') as outfh:
        outfh.write("sublibrary," + hdr0 + "\n")
        for lab, meta_path, r in zip(labels, cells_files, per_rows):
            written = 0
            with open(meta_path, 'rt', newline='') as fh:
                _ = fh.readline()  # skip header
                for line in fh:
                    outfh.write(f"{lab},{line.rstrip()}\n")
                    written += 1
            if written != r:
                sys.exit(f"[ERROR] Wrote {written} metadata rows but MTX rows={r} in {meta_path}")

    # 4) Stream and merge MTX row-wise (offset row index)
    mtx_out_path = outdir / ("count_matrix.mtx.gz" if args.gzip_out else "count_matrix.mtx")
    open_mtx_out = (lambda p: gzip.open(p, 'wt')) if args.gzip_out else (lambda p: open(p, 'wt'))

    with open_mtx_out(mtx_out_path) as outfh:
        outfh.write("%%MatrixMarket matrix coordinate real general\n")
        outfh.write(f"{total_rows} {total_cols} {total_nnz}\n")
        row_offset = 0
        for mtx_path, r in zip(mtx_files, per_rows):
            stream_mtx_data_with_row_offset(mtx_path, outfh, row_offset)
            row_offset += r

    print("[OK] Merged DGE_filtered (row-wise)")
    print(f"     outdir        : {outdir}")
    print(f"     cells (rows)  : {total_rows}")
    print(f"     genes (cols)  : {total_cols}")
    print(f"     nnz entries   : {total_nnz}")
    print(f"     matrix output : {mtx_out_path}")

if __name__ == "__main__":
    main()