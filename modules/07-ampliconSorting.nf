#!/usr/bin/env nextflow

process ampliconsorting_DeltaReads {
        tag "Sorting Delta reads of ${sample} suspect for coinfection"
        container 'lindenb/jvarkit:1b2aedf24'

        publishDir (
        path: "${params.out_dir}/07-AmpliconSorting",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path(bam), path(bai), path(vcf)
        path jvarkit_jar
        path sort_delta_reads

        output:
        tuple val(sample), path ("*.bam"), emit: delta_bam 

        script:
        """
        java -jar ${jvarkit_jar} samjdk -f ${sort_delta_reads} ${bam} \
        -o ${sample}_sorted_delta_reads.bam
        """
}

process ampliconsorting_OmicronReads {
    tag "Sorting Omicron reads of ${sample} suspect for coinfection"
    container 'lindenb/jvarkit:1b2aedf24'

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(bam), path(bai), path(vcf)
    path jvarkit_jar
    path sort_delta_reads

    output:
    tuple val(sample), path ("*.bam"), emit: omicron_bam

    script:
    """
    java -jar ${jvarkit_jar} samjdk -f ${sort_omicron_reads} ${bam_file} \
    -o ${sample}_sorted_omicron_reads.bam
    """
}

process ampliconsorting_consensusFasta {
    tag "Creating consensus from sorted reads"

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path (delta_bam)
    path (omicron_bam)

    output: 
    tuple val(sample), path ("*delta.fasta"), path ("omicron.fasta"), emit: ampliconsorting_consensus_fasta

    script:
    """
    """
}

process ampliconsorting_lineageAssignment_Pangolin {
    tag "Lineage assignment of sorted reads using Pangolin tool"

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (delta_fasta), path (omicron_fasta)

    output:
    
    script:
    """
    """
}