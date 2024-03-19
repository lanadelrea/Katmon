#!/usr/bin/env nextflow

process report {
    tag "Generating Summary Report"
    container 'jmeigs1/rscript:latest'

    publishDir(
    path: "${params.out_dir}/06-report",
    mode: 'copy',
    overwrite: 'true'
    )

    input
    path (bammix_plot)
    path (freyja_plot)
    path (sample), path(aaf_tsv), path(aafplot_mut)
    path (sample), path(aafplot_amp)
    path 

    output:

    script:
    """
    Rscript summary-report.R ${bammix_plot} ${freyja_plot} ${aafplot_mut} ${aafplot_amp}
    """
}