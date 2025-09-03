#!/usr/bin nextflow

process bammix {
        tag "Looking for positions with nucleotide mixtures"
        container 'ufuomababatunde/bammix:v1.1.0'
//        conda '/home/bdmu/miniforge3/envs/bammix'

        cpus 1

        containerOptions = "-v ${params.in_dir}:/data"

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (nextclade_tsv)

        output:
        path '*.csv'
        path '*.pdf'
        path 'bammixFlags.csv', emit: bammixflagged_csv

        script:
        """
        bammix.py ${nextclade_tsv} ${params.in_dir} ${params.bammix_thresh}
        """
}

process bammix_01_edited {
        tag "Looking for positions with nucleotide mixtures"
        container 'ufuomababatunde/bammix:v1.1.0'
//        conda '/home/bdmu/miniforge3/envs/bammix'
        cpus 1

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (nextclade_tsv)


        output:
        stdout

        script:
        """
        bammix-01-edited.py ${nextclade_tsv}
        """
}

process bammix_02_edited {
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (bam)
        path (bai)
        val (snps)

        output:
        tuple val (bai.simpleName), path ('*.csv'), emit: bammix_csv
        path ('*.pdf')

        script:
        """
        bammix -b ${bam} \
        -o ${bai.simpleName} \
        -p ${snps.join(' ')}
        """
}

process bammix_03_edited {
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val (sample), path (bammix_csv)

        output:
        path ('*_bammix_flags.csv'), emit: bammix_flags_csv

        script:
        """
        bammix-02-edited.py ${sample} ${bammix_csv} ${params.bammix_thresh}
        """
}

process bam_filter {
        tag "Determining samples flagged for having nucleotide mixtures"

        publishDir (
        path: "${params.out_dir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (bammixflagged_csv)

        output:
        path ('*.txt'), emit: samples_txt

        script:
        """
        get-bammix-flagged-samples.py ${bammixflagged_csv}
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
        tuple val(sample), path (bam)
        path (reference)

        output:
        path ("*.mpileup"), emit: mpileup

        script:
        """
        bcftools mpileup \
        -f ${reference} \
        ${bam} \
        -o ${sample}.mpileup
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