#!/usr/bin/env nextflow

process report {
//    errorStrategy = 'ignore'
    tag "Generating Summary Report"
    container 'ufuomababatunde/rmarkdown:1.1.0'

    publishDir(
    path: "${params.out_dir}/07-Report",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path (report_r)
    path (lineage_assignment)
    path (bammix_plot)
    path (freyja_plot)
    path (aafplot_mut)
    path (aafplot_amp)
    path (virstrain_tsv)
    path (report_rmd)


    output:
    path ("*.html")

    script: 
    """
    Rscript ${report_r} \
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
    path: "${params.out_dir}/087-Report",
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