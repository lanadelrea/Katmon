#!/usr/bin nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures"
        container 'ufuomababatunde/bammix:v1.1.0'

        cpus 1

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (nextclade_tsv)
        path (bam)
        path (bam_index)

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
        tag "Determining samples flagged for having nucleotide mixtures"

        publishDir (
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (bammixflagged_csv)

        output:
        stdout

        script:
        """
        bammix-flagged-2.py ${bammixflagged_csv}
        """
}

process makevcf {
        tag "Making vcf file from bam file of ${sample}"
        container 'staphb/bcftools:latest'

        publishDir(
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        val (sample)
        path GISAID_reference

        output:
        path ("*.mpileup"), emit: mpileup

        script:
        """
        bcftools mpileup -f ${GISAID_reference} ${params.in_dir}/${sample}.bam -o ${sample}.mpileup
        """
}

process bcftools {
        tag "Making vcf file from bam file of ${mpileup.BaseName}"
        container 'staphb/bcftools:latest'

        publishDir(
        path: "${params.out_dir}/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (mpileup)

        output:
        tuple val(mpileup.BaseName), path("*.vcf"), emit: filtered_vcf

        script:
        """
        bcftools call -mv -O b -o ${mpileup.BaseName}.bcf ${mpileup}
        bcftools view -O v -o ${mpileup.BaseName}.vcf ${mpileup.BaseName}.bcf
        """
}