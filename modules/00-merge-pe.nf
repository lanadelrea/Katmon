#!/usr/bin/env nextflow

process merge{
    tag "Merging paired-end reads"
    container 'pipecraft/pandaseq:2.11'

    input:
    tuple val(sample), path(forward), path(reverse)

    output:
    tuple val(sample), path("${sample}.merged.fastq.gz"), emit: fastq

    script:
    """
    pandaseq -f ${forward} -r ${reverse} | gzip -c > ${sample}.merged.fastq.gz
    """
}