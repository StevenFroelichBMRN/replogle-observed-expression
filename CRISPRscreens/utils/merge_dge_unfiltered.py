#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Out-of-core row-wise merge for Parse DGE_unfiltered matrices (cells = rows, genes = columns).

For each sublibrary, expects:
  <SUB>/MPRA/DGE_unfiltered/count_matrix.mtx(.gz)
  <SUB>/MPRA/DGE_unfiltered/all_genes.csv
  <SUB>/MPRA/DGE_unfiltered/cell_metadata.csv

This script:
  • Verifies that all_genes.csv are identical (line-by-line).
  • Reads only MTX headers to accumulate total rows/cols/nnz.
  • Aligns cell_metadata rows to MTX rows:
       - If rows == MTX rows: stream all rows.
       - Else: auto-detect a boolean/0-1 "keep" column whose sum == MTX rows, and stream only those rows.
       - Or force a specific keep column via --keep_col.
  • Streams each MTX and writes a merged Matrix Market by offsetting the ROW index (i).
  • Streams and concatenates cell metadata, adding a 'sublibrary' column.
  • Optionally gzips the merged MTX.

Usage:
  python3 merge_dge_unfiltered_rows_streamed.py \
    --sublibs output_sublibrary_01 ... output_sublibrary_16 \
    --outdir output_combined_unfiltered_rows \
    --gzip_out

If some cell_metadata files have extra rows, you can force the keep column, e.g.:
  --keep_col is_cell    (or pass_filter, classifier, etc.)
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

# ---------- metadata helpers (keep-column detection) ----------
COMMON_KEEP_COLS = [
    "is_cell","pass","pass_filter","keep","keep_cell","cell_selected",
    "final_cell_call","classifier","selected","kept","valid_cell"
]

def try_sum_keep_col_stream(meta_path, colname, dialect):
    """Return (sum_keep, col_index) if column exists and can be coerced to 0/1; else (None, idx_or_None)."""
    with open(meta_path, 'rt', newline='') as fh:
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)
        if colname not in header:
            return None, None
        idx = header.index(colname)
        total = 0
        for row in reader:
            try:
                v = row[idx].strip().lower()
                if v in ("1","true","t","yes","y"):
                    total += 1
                elif v in ("0","false","f","no","n",""):
                    pass
                else:
                    total += int(float(v))  # tolerate numeric strings
            except Exception:
                return None, idx
        return total, idx

def detect_keep_column(meta_path, nrows_target):
    """Auto-detect a keep column whose sum equals nrows_target; returns (name, idx, dialect) or (None, None, dialect)."""
    with open(meta_path, 'rt', newline='') as fh:
        sample = fh.read(4096)
        fh.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample) if sample else csv.get_dialect("excel")
        except Exception:
            dialect = csv.get_dialect("excel")
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)

    # Try common names first
    for name in COMMON_KEEP_COLS:
        if name in header:
            s, idx = try_sum_keep_col_stream(meta_path, name, dialect)
            if s == nrows_target:
                return name, idx, dialect

    # Brute-force every column (0/1/bool sums only)
    sums = [0] * len(header)
    valid = [True] * len(header)
    with open(meta_path, 'rt', newline='') as fh:
        reader = csv.reader(fh, dialect=dialect)
        _ = next(reader)
        for row in reader:
            for i, x in enumerate(row):
                if not valid[i]:
                    continue
                v = x.strip().lower()
                if v in ("1","true","t","yes","y"):
                    sums[i] += 1
                elif v in ("0","false","f","no","n",""):
                    pass
                else:
                    try:
                        sums[i] += int(float(v))
                    except Exception:
                        valid[i] = False
    for i, (ok, s) in enumerate(zip(valid, sums)):
        if ok and s == nrows_target:
            return header[i], i, dialect
    return None, None, dialect

def stream_write_metadata_subset(meta_path, out_fh, sublib_label, keep_idx, dialect):
    """Stream rows where keep==1/True; write with leading sublibrary column."""
    written = 0
    with open(meta_path, 'rt', newline='') as fh:
        reader = csv.reader(fh, dialect=dialect)
        header = next(reader)
        for row in reader:
            v = row[keep_idx].strip().lower()
            is_keep = (v in ("1","true","t","yes","y"))
            if not is_keep:
                try:
                    is_keep = int(float(v)) > 0
                except Exception:
                    is_keep = False
            if is_keep:
                out_fh.write(sublib_label + "," + ",".join(row) + "\n")
                written += 1
    return written

def stream_write_metadata_all(meta_path, out_fh, sublib_label, dialect):
    """Stream all rows; write with leading sublibrary column."""
    written = 0
    with open(meta_path, 'rt', newline='') as fh:
        reader = csv.reader(fh, dialect=dialect)
        _ = next(reader)
        for row in reader:
            out_fh.write(sublib_label + "," + ",".join(row) + "\n")
            written += 1
    return written

def count_metadata_rows_fast(meta_path):
    """Fast line-count excluding header (for sizing checks)."""
    if str(meta_path).endswith(".gz"):
        with opengz(meta_path, 'rt') as fh:
            _ = fh.readline()
            return sum(1 for _ in fh)
    with open(meta_path, 'rt', newline='') as fh:
        _ = fh.readline()
        return sum(1 for _ in fh)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="Row-wise merge of DGE_unfiltered across sublibraries (cells=rows).")
    ap.add_argument("--sublibs", nargs='+', required=True, help="Sublibrary root folders")
    ap.add_argument("--outdir", required=True, help="Output directory for merged DGE")
    ap.add_argument("--rel_mtx",   default="MPRA/DGE_unfiltered/count_matrix.mtx",
                    help="Relative path to matrix.mtx (or .mtx.gz) under each sublibrary")
    ap.add_argument("--rel_genes", default="MPRA/DGE_unfiltered/all_genes.csv",
                    help="Relative path to all_genes.csv")
    ap.add_argument("--rel_cells", default="MPRA/DGE_unfiltered/cell_metadata.csv",
                    help="Relative path to cell_metadata.csv")
    ap.add_argument("--sublib_label", default=None,
                    help="Comma-separated labels; defaults to directory names")
    ap.add_argument("--keep_col", default=None,
                    help="Force this metadata keep column (must sum to MTX rows)")
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

    # 2) First pass: read MTX sizes; plan metadata handling; validate columns (genes)
    total_rows, total_cols, total_nnz = 0, None, 0
    per_rows = []
    meta_plans = []  # dict per sublib: {"mode": "all"|"subset", "keep_idx": int or None, "dialect": dialect}
    hdr0 = None

    for mtx_path, meta_path in zip(mtx_files, cells_files):
        r, c, nnz = mtx_sizes_from_header(mtx_path)
        total_rows += r
        total_nnz  += nnz
        per_rows.append(r)
        if total_cols is None:
            total_cols = c
        elif total_cols != c:
            sys.exit(f"[ERROR] Gene (column) count differs across sublibraries: saw {total_cols} vs {c}.")

        # metadata header & dialect
        with open(meta_path, 'rt', newline='') as fh:
            sample = fh.read(4096); fh.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample) if sample else csv.get_dialect("excel")
            except Exception:
                dialect = csv.get_dialect("excel")
            reader = csv.reader(fh, dialect=dialect)
            header = next(reader)
        if hdr0 is None:
            hdr0 = header
        elif header != hdr0:
            sys.exit(f"[ERROR] cell_metadata header differs in {meta_path}")

        # Decide metadata mode
        m_rows = count_metadata_rows_fast(meta_path)
        if m_rows == r:
            meta_plans.append({"mode":"all", "keep_idx":None, "dialect":dialect})
        else:
            if args.keep_col:
                s, idx = try_sum_keep_col_stream(meta_path, args.keep_col, dialect)
                if s != r:
                    sys.exit(f"[ERROR] --keep_col '{args.keep_col}' sums to {s}, but MTX rows={r} in {meta_path}")
                meta_plans.append({"mode":"subset", "keep_idx":idx, "dialect":dialect})
            else:
                name, idx, dialect2 = detect_keep_column(meta_path, r)
                if name is None:
                    sys.exit(f"[ERROR] Could not auto-detect a 0/1 keep column in {meta_path} "
                             f"that sums to MTX rows ({r}). Use --keep_col <column_name>.")
                meta_plans.append({"mode":"subset", "keep_idx":idx, "dialect":dialect2})

    # sanity: columns==n_genes
    if total_cols != n_genes:
        sys.exit(f"[ERROR] Gene count mismatch: MTX columns={total_cols} vs all_genes rows={n_genes}")

    # 3) Write merged metadata (streamed; prepend 'sublibrary')
    out_meta = outdir / "cell_metadata.csv"
    with open(out_meta, 'wt', newline='') as outfh:
        outfh.write("sublibrary," + ",".join(hdr0) + "\n")
        for lab, meta_path, plan, r in zip(labels, cells_files, meta_plans, per_rows):
            if plan["mode"] == "all":
                written = stream_write_metadata_all(meta_path, outfh, lab, plan["dialect"])
            else:
                written = stream_write_metadata_subset(meta_path, outfh, lab, plan["keep_idx"], plan["dialect"])
            if written != r:
                sys.exit(f"[ERROR] Wrote {written} metadata rows but MTX rows={r} in {meta_path}")

    # 4) Stream & merge MTX by offsetting ROW index
    mtx_out_path = outdir / ("count_matrix.mtx.gz" if args.gzip_out else "count_matrix.mtx")
    open_mtx_out = (lambda p: gzip.open(p, 'wt')) if args.gzip_out else (lambda p: open(p, 'wt'))

    with open_mtx_out(mtx_out_path) as outfh:
        outfh.write("%%MatrixMarket matrix coordinate real general\n")
        outfh.write(f"{total_rows} {total_cols} {total_nnz}\n")
        row_offset = 0
        for mtx_path, r in zip(mtx_files, per_rows):
            stream_mtx_data_with_row_offset(mtx_path, outfh, row_offset)
            row_offset += r

    print("[OK] Merged DGE_unfiltered (row-wise)")
    print(f"     outdir        : {outdir}")
    print(f"     cells (rows)  : {total_rows}")
    print(f"     genes (cols)  : {total_cols}")
    print(f"     nnz entries   : {total_nnz}")
    print(f"     matrix output : {mtx_out_path}")

if __name__ == "__main__":
    main()