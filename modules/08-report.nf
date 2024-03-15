#!/usr/bin

process report {
    tag "Generating Summary Report"
    container ''

    publishDir(
    path: "${params.out_dir}/06-report",
    mode: 'copy',
    overwrite: 'true'
    )

    input:

    output:

    script:
    """
    summary-report.py  
    """
}