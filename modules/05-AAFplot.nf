#!/usr/bin

process aafplot_mutations {
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
    tuple val(sample), path ("*.tsv"), path("*.png"), emit: aafplot_mut

    script:
    """
    touch ${sample}_read_depth.tsv
    
    aafplot.py ${vcf} \
    ${baseDir}/assets/mutations.tsv \
    ${filtered_bam} \
    ${sample} \
    ${params.out_dir}/05-AAFplot/${sample}_read_depth.tsv \
    ${sample}_AAFplot_mutations.png
    """
}

process aafplot_amplicons {
    tag "Plotting alternative allele fraction per mutation of ${sample}"

    publishDir(
    path: "${params.out_dir}/05-AAFplot",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (read_depth_tsv), path("*.png")

    output:
    tuple val(sample), path("*.png"), emit: aafplot_amp

    script:
    """
    aafplot_amplicon.py ${baseDir}/assets/primer-schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed ${read_depth_tsv} ${sample}_AAFplot_amplicons.png
    """
}