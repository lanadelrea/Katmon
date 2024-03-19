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
    path sort_omicron_reads

    output:
    tuple val(sample), path ("*.bam"), emit: omicron_bam

    script:
    """
    java -jar ${jvarkit_jar} samjdk -f ${sort_omicron_reads} ${bam} \
    -o ${sample}_sorted_omicron_reads.bam
    """
}

process ampliconsorting_samtools {
    tag "Creating vcf from Delta and Omicron sorted reads"
    container 'pegi3s/samtools_bcftools:latest'

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (delta_bam)
    tuple val(sample), path (omicron_bam)
    path (reference)

    output:
    path("*.bai")
    tuple val(sample), path("*delta.vcf"), path("*omicron.vcf"), emit: vcf

    script:
    """
    samtools index ${delta_bam}
    bcftools mpileup -f ${reference} ${delta_bam} > ${sample}_delta.mpileup
    bcftools call -mv -O b -o ${sample}_delta.bcf ${sample}_delta.mpileup
    bcftools view -O v -o ${sample}_delta.vcf ${sample}_delta.bcf

    samtools index ${omicron_bam}
    bcftools mpileup -f ${reference} ${omicron_bam} > ${sample}_omicron.mpileup
    bcftools call -mv -O b -o ${sample}_omicron.bcf ${sample}_omicron.mpileup
    bcftools view -O v -o ${sample}_omicron.vcf ${sample}_omicron.bcf
    """
}

process ampliconsorting_bgzip {
    container 'vandhanak/bcftools:1.3.1'
    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(delta_vcf), path(omicron_vcf)

    output:
    tuple val(sample), path("*delta.vcf.gz"), path("*omicron.vcf.gz"), path("*.vcf.gz.tbi"), emit: vcfgz

    script:
    """
    bgzip ${delta_vcf}
    bgzip ${omicron_vcf}

    tabix -p vcf ${delta_vcf}.gz
    tabix -p vcf ${omicron_vcf}.gz
    """
}

process ampliconsorting_fasta {
    tag "Creating consensus from sorted reads"
    container 'pegi3s/samtools_bcftools:latest'

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(delta_vcfgz), path(omicron_vcfgz), path (vcfgz_tbi)
    path (reference)

    output:
    tuple val(sample), path("*delta_consensus.fasta"), path("*omicron_consensus.fasta"), emit: fasta

    script:
    """
    bcftools consensus -f ${reference} ${delta_vcfgz} > ${sample}_delta_consensus.fasta
    bcftools consensus -f ${reference} ${omicron_vcfgz} > ${sample}_omicron_consensus.fasta
    """
}

process ampliconsorting_fasta_rename {
    tag "Creating consensus from sorted reads"

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path(delta_fasta), path(omicron_fasta)

    output:
    tuple val(sample), path("*ampliconsorted.fasta"), emit: fasta

    script:
    """
    awk '/^>/{print ">${sample}_Delta_sorted_reads"; next}{print}' ${delta_fasta} > ${delta_fasta}
    awk '/^>/{print ">${sample}_Omicron_sorted_reads"; next}{print}' ${omicron_fasta} > ${omicron_fasta}
    
    cat ${delta_fasta} ${omicron_fasta} > ${sample}_ampliconsorted.fasta
    """
}

process ampliconsorting_lineageAssignment_Pangolin {
    tag "Lineage assignment of sorted reads using Pangolin tool"
    container 'staphb/pangolin:latest'

    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (fasta)

    output:
    path ("*.csv"), emit: pangolineageAssign
    
    script:
    """
    pangolin ${fasta} > ${sample}_ampliconsorted_lineage_assignment.csv
    """
}

process ampliconsorting_lineageAssignment_Nextclade {
    tag "Lineage assignment of sorted reads using Nextclade"
    container 'nextstrain/nextclade:latest'
    
    publishDir (
    path: "${params.out_dir}/07-AmpliconSorting",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    tuple val(sample), path (fasta)
    path (SC2_dataset)

    output:
    path ("*.tsv"), emit: nextcladeAssign

    script:
    """
    nextclade run --input-dataset ${SC2_dataset} --output-tsv=${sample}_ampliconsorted_nextclade.tsv ${fasta}
    """
}