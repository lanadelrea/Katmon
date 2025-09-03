#!/usr/bin/env nextflow

process bammixplot {
//    errorStrategy = 'ignore'
    tag "Plotting bammix plot for ${sample}"
    container 'ufuomababatunde/bammix:v1.1.0'

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(csv)

    output:
    path ("*.png"), emit: bammix_plot

    script:
    """
    bammix_plot.py ${csv} \
    ${sample}_bammix_plot.png \
    ${sample}
    """
}


process aafplots {
//    errorStrategy = 'ignore'
    tag "Plotting alternative allele fraction per mutation of ${sample}"
    container 'ufuomababatunde/bammix:v1.1.0' // to fix

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path (primer_scheme)
    tuple val(sample), path (processed_mut_tsv), path(vcf)

    output:
    tuple val (sample), path ("*.tsv"), emit: aafplot_mut_tsv
    path ("*mutations.png"), emit: aafplot_mut
    path ("*amplicons.png"), emit: aafplot_amp

    script:
    """
    aafplot_mutation.py ${vcf} \
    ${processed_mut_tsv} \
    ${params.in_dir}/${sample}.bam \
    ${sample} \
    ${sample}_read_depth.tsv \
    ${sample}_AAFplot_mutations.png

    aafplot_amplicon.py ${primer_scheme} \
    ${sample}_read_depth.tsv \
    ${sample} \
    ${processed_mut_tsv} \
    ${sample}_AAFplot_amplicons.png
    """
}