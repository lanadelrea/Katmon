#!/usr/bin nextflow

process bammix {
        tag "This is the original bammix process, but has been split into 3 processes for better flow"
        container 'ufuomababatunde/bammix:v1.1.0'
        cpus 1

        containerOptions = "-v ${params.indir}:/data"

        publishDir (
        path: "${params.outdir}/02-Bammix",
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
        bammix.py ${nextclade_tsv} ${params.indir} ${params.bammix_thresh}
        """
}

process positions {
        tag "Getting relevant positions"
        container 'ufuomababatunde/bammix:v1.1.0'
        cpus 1

        publishDir (
        path: "${params.outdir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (nextclade_tsv)


        output:
        stdout

        script:
        """
        bammix-01.py ${nextclade_tsv}
        """
}

process bammix_process {
        tag "Looking for positions with nucleotide mixtures"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.outdir}/02-Bammix",
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

process flagged_positions {
        tag "Getting flagged positions with nucleotide mixtures per sample"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.outdir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val (sample), path (bammix_csv)

        output:
        path ('*.csv'), emit: bammix_flags_csv

        script:
        """
        bammix-02.py ${sample} ${bammix_csv} ${params.bammix_thresh}
        """
}

process flagged_samples {
        tag "Determining samples flagged for having nucleotide mixtures"
        container 'ufuomababatunde/bammix:v1.1.0'

        publishDir (
        path: "${params.outdir}/02-Bammix",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (bammixflagged_csv)

        output:
        path ('*.txt'), emit: samples_txt

        script:
        """
        bammix-03.py ${bammixflagged_csv}
        """
}

process makevcf {
        tag "Making vcf file from bam file of ${sample}"
        container 'staphb/bcftools:latest'

        publishDir(
        path: "${params.outdir}/05-makeVCF",
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
        path: "${params.outdir}/05-makeVCF",
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