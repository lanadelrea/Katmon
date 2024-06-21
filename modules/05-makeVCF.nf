#!/usr/bin/env nextflow
process makevcf {
    tag "Making vcf file of high quality reads from bam file of ${filtered_bam.name}"
    container 'pegi3s/samtools_bcftools:latest'
    
    publishDir(
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
    )

    input:
    path filtered_bam
    path reference

    output:
    tuple val (filtered_bam.SimpleName), path ("*_filtered.sorted.bam"), path ("*_filtered.sorted.bam.bai"), path ("*.vcf"), emit: filtered_vcf

    script:
    """
    samtools sort ${filtered_bam} -o ${filtered_bam.SimpleName}_filtered.sorted.bam
    samtools index ${filtered_bam.SimpleName}_filtered.sorted.bam
    samtools mpileup -uf ${reference} ${filtered_bam.SimpleName}_filtered.sorted.bam > ${filtered_bam.SimpleName}.mpileup
    bcftools call -mv -O b -o ${filtered_bam.SimpleName}.bcf ${filtered_bam.SimpleName}.mpileup
    bcftools view -O v -o ${filtered_bam.SimpleName}.vcf ${filtered_bam.SimpleName}.bcf
    """
}