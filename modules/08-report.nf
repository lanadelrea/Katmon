#!/usr/bin/env nextflow

process report {
    tag "Generating Summary Report"
    container 'stefanfritsch/rmd_statup:latest'

    publishDir(
    path: "${params.out_dir}/08-Report",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path (bammix_plot)
    path (freyja_plot)
    tuple val(sample), path(tsv), path(aafplot_mut)
    tuple val(sample), path(aafplot_amp)

    output:
    path (".pdf")

    script:
    """
    Rscript ${baseDir}/bin/summary-report.R ${bammix_plot} ${freyja_plot} ${aafplot_mut} ${aafplot_amp}
    """
}