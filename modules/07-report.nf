#!/usr/bin/env nextflow

process report {
//    errorStrategy = 'ignore'
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
    tuple val (sample), path(tsv), path(aafplot_mut)
    path (aafplot_amp)
    path (virstrain_tsv)
    path (report_rmd)


    output:
    path ("*.html")

    script: 
    """
    Rscript ${baseDir}/bin/summary-report.R \
    ${lineage_assignment} \
    ${bammix_plot} \
    ${freyja_plot} \
    ${aafplot_mut} \
    ${aafplot_amp} \
    ${virstrain_tsv} \
    ${report_rmd}
    """
}

process report_no_flag {
//    errorStrategy = 'ignore'
    tag "Generating Summary Report"
    container 'ufuomababatunde/rmarkdown:1.1.0'

    publishDir(
    path: "${params.out_dir}/08-Report",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path (lineage_assignment)
    path (freyja_plot)
    path (virstrain_tsv)
    path (report_rmd_no_bammix)


    script: 
    """
    Rscript ${baseDir}/bin/summary-report-no-bammix.R \
    ${lineage_assignment} \
    ${freyja_plot} \
    ${virstrain_tsv} \
    ${report_rmd_no_bammix}
    """
}