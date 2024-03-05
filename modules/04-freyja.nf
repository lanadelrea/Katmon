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

        output:
        val true

        script:
        """
        mkdir -p ${params.out_dir}/04-Freyja/Demix_results
        freyja demix ${params.out_dir}/04-Freyja/${sample}_variants.tsv ${params.out_dir}/04-Freyja/${sample}_depth.tsv \
        --output ${params.out_dir}/04-Freyja/Demix_results \
        --depthcutoff 1
        """
}

process freyja_aggregate {
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )
        
        input:
        val ready

        output:
        path ("*.tsv"), emit: freyja_aggregated_file

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
        path aggregated_file

        output:
        path ("*.png")

        script:
        """
        freyja plot ${aggregated_file} --output freyja-lineage-abundance-plot.png
        """
}