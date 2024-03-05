#!/usr/bin

process freyja {
        tag "Identifying relative lineage abundances of sample ${bamfilePath.baseName} from potential mixed SARS-CoV-2 samples using Freyja"
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path bamfilePath

        output:
        tuple val(bamfilePath.baseName), path ("*variants.tsv"), path ("*depth.tsv"), emit: freyja_variants

        script:
        """
        freyja variants ${bamfilePath} \
        --variants ${bamfilePath.baseName}_variants.tsv \
        --depths ${bamfilePath.baseName}_depth.tsv \
        --ref ${baseDir}/assets/sars-cov-2/reference-sequence.fasta
        """
}

process freyja_demix {
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja/Demix_results",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val(sample), path ("*variants.tsv"), path ("*depth.tsv")

        script:
        """
        freyja demix ${sample}_variants.tsv ${sample}_depth.tsv --output ${params.out_dir}/04-Freyja/Demix_results
        """
}

process freyja_aggregate {
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )
        
        output:
        path ("*.tsv")

        script:
        """
        freyja aggregate ${params.out_dir}/04-Freyja/Demix_results --output aggregated-file.tsv
        """
}

process freyja_plot {
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path freyja_aggregate

        output:
        path ("*.png")

        script:
        """
        freyja plot aggregated-file.tsv --output freyja-lineage-abundance-plot.png
        """
}