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
    tuple val(sample), path (sorted_bam), path (sorted_bam_bai), path(vcf)
    path(bam)
    
    output:
    tuple val(sample), path("*.png")

    script:
    """
    aafplot.py ${vcf} ${baseDir}/assets/mutations.csv ${bam} ${sample} ${sample}_AAFplot.png
    """
}