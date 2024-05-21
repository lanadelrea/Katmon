#!/usr/bin/env nextflow

process report {
    tag "Generating Summary Report"
    container 'ufuomababatunde/rmarkdown:1.1.0'

    publishDir(
    path: "${params.out_dir}/08-Report",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path (lineage_assignment)
    path (bammix_plot)
    path (freyja_plot)
    tuple val(sample), path(tsv), path(aafplot_mut)
    tuple val(sample), path(aafplot_amp)
    path (report_rmd)
    path (virstrain_tsv)

    output:
    path ("*.html")

    script:
    """
    Rscript ${baseDir}/bin/summary-report.R ${lineage_assignment} ${bammix_plot} ${freyja_plot} ${aafplot_mut} ${aafplot_amp} ${report_rmd} ${sample} ${virstrain_tsv}
    """
}
