#!/usr/bin

process freyja {
        tag "Identifying relative lineage abundances of sample ${bamfilePath.baseName} from potential mixed SARS-CoV-2 samples"
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja/VariantsDepth_results",
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
        tuple val(sample), path (variants_tsv), path (depth_tsv)

        output:
        path ("*.tsv"), emit: tsv_demix


        script:
        """
        freyja demix ${variants_tsv} ${depth_tsv} --output ${sample}_demix.tsv
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
        path tsv_demix

        output:
        path ("*.tsv"), emit: freyja_aggregated_file

        script:
        """
        mkdir Demix_results
        cp ${tsv_demix} Demix_results/
        freyja aggregate Demix_results/ --output aggregated-file.tsv
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
        path ("*.png"), emit: freyja_plot

        script:
        """
        freyja plot ${aggregated_file} --output freyja-lineage-abundance-plot.png
        """
}