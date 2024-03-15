#!/usr/bin

process bammixplot {
    tag "Plotting bammix plot for ${sample}"

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (filtered_bam), path (filtered_bam_bai), path(vcf)

    output:
    path ("*.png")

    script:
    """
    bammixplot.py ${params.out_dir}/02-Bammix/${sample}_position_base_counts.csv \
    ${sample}_bammix_plot.png \
    ${sample}
    """
}

process aafplot_mutations {
    tag "Plotting alternative allele fraction per mutation of ${sample}"

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (filtered_bam), path (filtered_bam_bai), path(vcf)

    output:
    tuple val(sample), path ("*.tsv"), path("*.png"), emit: aafplot_mut

    script:
    """
    aafplot.py ${vcf} \
    ${baseDir}/assets/mutations.tsv \
    ${filtered_bam} \
    ${sample} \
    ${sample}_read_depth.tsv \
    ${sample}_AAFplot_mutations.png
    """
}

process aafplot_amplicons {
    tag "Plotting alternative allele fraction per mutation of ${sample}"

    publishDir(
    path: "${params.out_dir}/06-Plots",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (read_depth_tsv), path("*.png")

    output:
    tuple val(sample), path("*.png"), emit: aafplot_amp

    script:
    """
    aafplot_amplicon.py ${baseDir}/assets/primer-schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed \
    ${read_depth_tsv} \
    ${sample} \
    ${sample}_AAFplot_amplicons.png
    """
}