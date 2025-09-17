#!/usr/bin/env nextflow

process get_pos_mut {
        tag "Getting list of positions and mutations for specific lineages in ${sample}"
        container 'ufuomababatunde/bammix:v1.1.0' // to fix

        publishDir (
        path: "${params.outdir}/07-AmpliconSorting/PosMut_txt",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val (sample), path (vcf), path (lin_mut_tsv)

        output:
        tuple val (sample), path ("*lineage_A.txt"), emit: pos_mut_lineage_A // ${sample}_lineage_A.txt
        tuple val (sample), path ("*lineage_B.txt"), emit: pos_mut_lineage_B // ${sample}_lineage_B.txt

        script:
        """
        get_pos_mut.py ${sample} ${lin_mut_tsv}
        """
}

process ampliconsorting {
        tag "Sorting reads of ${sample} suspect for co-infection"
        container 'maelanade/jvarkit:1.0'

        publishDir (
        path: "${params.outdir}/07-AmpliconSorting",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val (sample), path (pos_mut_lineage_A), path (pos_mut_lineage_B)
        tuple val (sample), path (bam)
        path (jvarkit_jar)
        path (sort_reads)

        output:
        tuple val(sample), path ("*lineage_A.bam"), emit: sorted_bam_lineage_A // ${sample}_ampsort_lineage_A.bam
        tuple val(sample), path ("*lineage_B.bam"), emit: sorted_bam_lineage_B // ${sample}_ampsort_lineage_B.bam

        script:
        """
        samjdk -DmutationFilePath=${pos_mut_lineage_A} \
        -f ${sort_reads} \
        ${bam} \
        -o ${sample}_lineage_A.bam

        samjdk -DmutationFilePath=${pos_mut_lineage_B} \
        -f ${sort_reads} \
        ${bam} \
        -o ${sample}_lineage_B.bam
        """
}

process consensus {
        tag "Creating consensus fasta from the sorted reads"
        container 'maelanade/consensus:1.0'

        publishDir (
        path: "${params.outdir}/07-AmpliconSorting",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path (lineage_A_bam), path (lineage_B_bam)
        path (reference)

        output:
        path ("*.bai")
        tuple val(sample), path ("*lineage_A.consensus.fasta"), emit: consensus_lineage_A
        tuple val(sample), path ("*lineage_B.consensus.fasta"), emit: consensus_lineage_B

        script:
        """
        samtools index ${lineage_A_bam}
        bcftools mpileup -f ${reference} ${lineage_A_bam} > ${sample}_lineage_A.mpileup
        bcftools call -mv -O b -o ${sample}_lineage_A.bcf ${sample}_lineage_A.mpileup
        bcftools view -O v -o ${sample}_lineage_A.vcf ${sample}_lineage_A.bcf

        samtools index ${lineage_B_bam}
        bcftools mpileup -f ${reference} ${lineage_B_bam} > ${sample}_lineage_B.mpileup
        bcftools call -mv -O b -o ${sample}_lineage_B.bcf ${sample}_lineage_B.mpileup
        bcftools view -O v -o ${sample}_lineage_B.vcf ${sample}_lineage_B.bcf

        bgzip ${sample}_lineage_A.vcf
        bgzip ${sample}_lineage_B.vcf

        tabix -p vcf ${sample}_lineage_A.vcf.gz
        tabix -p vcf ${sample}_lineage_B.vcf.gz

        bcftools consensus -f ${reference} ${sample}_lineage_A.vcf.gz > ${sample}_lineage_A.consensus.fasta
        bcftools consensus -f ${reference} ${sample}_lineage_B.vcf.gz > ${sample}_lineage_B.consensus.fasta
        """
}

process renamefasta {
	container 'nanozoo/seqkit:latest'
	tag "Renaming consensus fasta sequence from amplicon sorting"

	publishDir (
	path: "${params.outdir}/07-AmpliconSorting",
	mode: 'copy',
	overwrite: 'true',
	)

	input:
	tuple val(sample), path(fasta_lineage_A), path(fasta_lineage_B)

	output:
	path ("*ampsort.consensus.fasta"), emit: fasta

	script:
	"""
	seqkit replace -p 'MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome' -r ${sample}_lineage_A ${fasta_lineage_A} > ${sample}_lineage_A_consensus_renamed.fasta
        seqkit replace -p 'MN908947.3 Severe acute respiratory syndrome coronavirus 2 isolate Wuhan-Hu-1, complete genome' -r ${sample}_lineage_B ${fasta_lineage_B} > ${sample}_lineage_B_consensus_renamed.fasta
        cat ${sample}_lineage_A_consensus_renamed.fasta ${sample}_lineage_B_consensus_renamed.fasta > ${sample}.ampsort.consensus.fasta
	"""
}

process pangolin {
        tag "Lineage assignment of sorted reads using Pangolin tool"
        container 'staphb/pangolin:latest'

        publishDir (
        path: "${params.outdir}/07-AmpliconSorting/Lineage_Assignment",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (fasta) // ampsort_consensus_final

        output:
        path ('*.csv'), emit: csv

        script:
        """
        pangolin ${fasta} --outfile ampsort.pangolin.csv
        """
}

process nextclade {
        tag "Lineage assignment of sorted reads using Nextclade"
        container 'nextstrain/nextclade:latest'
    
        publishDir (
        path: "${params.outdir}/07-AmpliconSorting/Lineage_assignment",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (fasta) // ampsort_consensus_final
        path SC2_dataset

        output:
        path ('*.tsv'), emit: tsv

        script:
        """
        nextclade run --input-dataset ${SC2_dataset} --output-tsv=ampsort.nextclade.tsv ${fasta}
        """
}

process summary {
        tag "Creating summary table for lineage assignment of sorted reads for ${sample}"
        container 'ufuomababatunde/bammix:v1.1.0' // to fix

        publishDir (
        path: "${params.outdir}/07-AmpliconSorting/Lineage_assignment",
        mode: 'copy',
        overwrite: 'true',
        )

        input:
        path (ampsort_pangolin_csv)
        path (ampsort_nextclade_tsv)

        output:
        path ('*.tsv'), emit: table

        script:
        """
        lineageAssign_table.py ${ampsort_pangolin_csv} ${ampsort_nextclade_tsv} 'ampsort'
        """
}