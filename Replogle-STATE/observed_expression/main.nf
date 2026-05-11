#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input_bucket = 's3://r6333-pep-nppc-oi-bmn333-dev/gwas_testing_nwk/Gene_Editing/replogle_nadig/filtered'
params.outdir = 's3://r6333-pep-nppc-oi-bmn333-dev/gwas_testing_nwk/Gene_Editing/replogle_nadig/observed_expression'

workflow {
    ch_files = Channel.of(
        file("${params.input_bucket}/GSE264667_hepg2_raw_singlecell_01_filtered.h5ad"),
        file("${params.input_bucket}/GSE264667_jurkat_raw_singlecell_01_filtered.h5ad"),
        file("${params.input_bucket}/K562_essential_normalized_singlecell_01_filtered.h5ad"),
        file("${params.input_bucket}/rpe1_normalized_singlecell_01_filtered.h5ad")
    )

    COMPUTE_EXPRESSION(ch_files)
}

process COMPUTE_EXPRESSION {
    publishDir params.outdir, mode: 'copy'

    input:
    path h5ad_file

    output:
    path "*.parquet"

    script:
    """
    compute_observed_expression.py ${h5ad_file}
    """
}
