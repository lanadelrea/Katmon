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
    tuple val(sample), path(vcf)

    output:
    path ("*.png"), emit: bammix_plot

    script:
    """
    bammix_plot.py ${params.out_dir}/02-Bammix/${sample}_position_base_counts.csv \
    ${sample}_bammix_plot.png \
    ${sample}
    """
}


process aafplot_mutations {
//    errorStrategy = 'ignore'
    tag "Plotting alternative allele fraction per mutation of ${sample}"
    container 'ufuomababatunde/bammix:v1.1.0' // to fix

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(vcf)

    output:
    tuple val (sample), path ("*.tsv"), emit: aafplot_mut_tsv
    path ("*.png"), emit: aafplot_mut

    script:
    """
    aafplot_mutation.py ${vcf} \
    ${params.out_dir}/04-Freyja/Mutations/${sample}_mutations.tsv \
    ${params.in_dir}/${sample}.bam \
    ${sample} \
    ${sample}_read_depth.tsv \
    ${sample}_AAFplot_mutations.png
    """
}

process aafplot_amplicons {
    tag "Plotting alternative allele fraction per mutation of ${sample}"
    container 'ufuomababatunde/bammix:v1.1.0'

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path primer_scheme
    tuple val (sample), path (read_depth_tsv)

    output:
    path ("*.tsv")
    path ("*.png"), emit: aafplot_amp

    script:
    """
    aafplot_amplicon.py ${primer_scheme} \
    ${read_depth_tsv} \
    ${sample} \
    ${params.out_dir}/04-Freyja/Mutations/${sample}_mutations.tsv \
    ${sample}_AAFplot_amplicons.png
    """
}