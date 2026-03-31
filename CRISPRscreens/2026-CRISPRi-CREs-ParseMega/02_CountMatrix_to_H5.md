# Process count matrices and converts them into H5 format.

Regular commands cannot handle such large matrices due to memory restrictions. Therefore, will convert to an H5 file format which can be read via BPCells.

See: https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette

```bash
python utils/mtx_to_h5.py \
  --mtx  output_combined_filtered/count_matrix.mtx \
  --genes output_combined_filtered/all_genes.csv \
  --barcodes output_combined_filtered/cell_metadata_sub.csv \
  --gene-col gene_name --barcode-col bc_wells_s \
  --out output_combined_filtered/count_matrix.h5 \
  --buffer-nnz 10000000
```