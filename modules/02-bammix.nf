#!/usr/bin nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures using bammix"
        container 'ufuomababatunde/bammix:v1.1.0'

        cpus 1

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path(nextclade_tsv)
        path(bam)
        path(bam_index)

        output:
        path '*.csv'
        path '*.pdf'
        path 'bammixFlags.csv', emit: bammixflagged_csv

        script:
        """
        bammix.py ${nextclade_tsv} ${params.in_dir}
        """
}

process bam_filter {
        errorStrategy = 'ignore'
        tag "Determining samples flagged by bammix then filtering the high quality reads"
        container 'pegi3s/samtools_bcftools:latest'

        publishDir (
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path bammixflagged_csv

        output:
        path ('*.bam'), emit: filtered_bam

        script:
        """
        bammix-flagged-sample-name.py ${bammixflagged_csv} ${params.in_dir}
        """
}

process makevcf {
    errorStrategy = 'ignore'
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