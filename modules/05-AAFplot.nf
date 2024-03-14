#!/usr/bin

process aafplot {
    tag "Plotting alternative allele fraction of ${sample}"
    container 'bitnami/python:latest'

    publishDir(
    path: "${params.out_dir}/05-AAFplot",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(vcf)
    path(bam)
    
    output:
    tuple val(sample), path("*.png")

    script:
    """
    aafplot.py ${vcf} ${baseDir}/assets/mutations.csv ${bam} ${sample} ${sample}_AAFplot.png
    """
}