#!/usr/bin

process aafplot {
    tag "Plotting alternative allele fraction per mutation of ${sample}"
//    container 'biocontainers/pandas:1.5.1_cv1'

    publishDir(
    path: "${params.out_dir}/05-AAFplot",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (filtered_bam), path (filtered_bam_bai), path(vcf)

    output:
    tuple val(sample), path ("*_read_depth.tsv"), path("*_AAFplot_mutations.png"), path("*_AAFplot_amplicons.png"), emit: aafplots

    script:
    """
    mkdir -p ${params.out_dir}/05-AAFplot/
    touch ${sample}_read_depth.tsv
    aafplot.py ${vcf} \
    ${baseDir}/assets/mutations.tsv \
    ${filtered_bam} \
    ${sample} \
    ${params.out_dir}/05-AAFplot/${sample}_read_depth.tsv \
    ${sample}_AAFplot_mutations.png \
    ${baseDir}/assets/primer-schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed \
    ${sample}_AAFplot_amplicons.png
    """
}

