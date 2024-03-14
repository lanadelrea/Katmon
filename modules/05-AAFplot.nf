#!/usr/bin

process aafplot {
    tag "Plotting alternative allele fraction of ${sample}"
//    container 'biocontainers/pandas:1.5.1_cv1'

    publishDir(
    path: "${params.out_dir}/05-AAFplot",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (filtered_bam), path (filtered_bam_bai), path(vcf)

    output:
    tuple val(sample), path("*.png"), emit: aafplot_1

    script:
    """
    aafplot.py ${vcf} ${baseDir}/assets/mutations.csv ${filtered_bam} ${sample} ${sample}_AAFplot.png
    """
}