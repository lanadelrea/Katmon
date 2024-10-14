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
//        errorStrategy = 'ignore'
        tag "Determining samples flagged by bammix then filtering the high quality reads"
        container 'staphb/samtools:latest'

        publishDir (
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path bammixflagged_csv

        output:
        stdout
//        path ('*.bam'), emit: filtered_bam

        script:
        """
        bammix-flagged-2.py ${bammixflagged_csv}
        """
}

process makevcf {
//    errorStrategy = 'ignore'
        tag "Making vcf file of high quality reads from bam file of ${sample}"
        container 'staphb/samtools:latest'

        publishDir(
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        val sample
        path reference

        output:
        tuple path ("*_filtered.sorted.bam"), path ("*_filtered.sorted.bam.bai"), emit: filtered_bam_bai
        path ("*.mpileup"), emit: mpileup

        script:
        """
        samtools view -q 20 -b ${params.in_dir}/${sample}.bam > ${sample}.filtered.bam 
        samtools sort ${sample}.filtered.bam -o ${sample}_filtered.sorted.bam
        samtools index ${sample}_filtered.sorted.bam
        samtools mpileup -uf ${reference} ${sample}_filtered.sorted.bam > ${sample}.mpileup
        """
}

process makevcf_2 {
//    errorStrategy = 'ignore'
        tag "Making vcf file of high quality reads from bam file of ${sample}"
        container 'staphb/samtools:latest'

        publishDir(
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        val sample
        path reference

        output:
        tuple path ("*_filtered.sorted.bam"), path ("*_filtered.sorted.bam.bai"), emit: filtered_bam_bai
        path ("*.mpileup"), emit: mpileup

        script:
        """
        samtools mpileup -f ${reference} ${params.in_dir}/${sample}.bam > ${sample}.mpileup
        """
}

process bcftools {
        tag "Making vcf file of high quality reads from bam file of ${sample}"
        container 'staphb/bcftools:latest'

        publishDir(
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path(mpileup)

        output:
        tuple val(mpileup.SimpleName), path("*.vcf"), emit: filtered_vcf

        script:
        """
        bcftools call -mv -O b -o ${mpileup.SimpleName}.bcf ${mpileup}
        bcftools view -O v -o ${mpileup.SimpleName}.vcf ${mpileup.SimpleName}.bcf
        """
}