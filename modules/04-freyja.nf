#!/usr/bin/env nextflow

process freyja {
        tag "Identifying relative lineage abundances of samples to see potential coinfection"
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja/VariantsDepth_results",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path reference
        path bamfilePath

        output:
        tuple val(bamfilePath.baseName), path ("*variants.tsv"), path ("*depth.tsv"), emit: freyja_variants

        script:
        """
        freyja variants ${bamfilePath} \
        --variants ${bamfilePath.baseName}_variants.tsv \
        --depths ${bamfilePath.baseName}_depth.tsv \
        --ref ${reference}
        """
}

process freyja_demix {
        tag "Identifying relative lineage abundances of sample ${sample} from potential mixed SARS-CoV-2 samples"
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
        freyja demix ${variants_tsv} ${depth_tsv} --output ${sample}_demix.tsv --depthcutoff 20
        """
}

process freyja_aggregate {
        tag "Aggregating Freyja demix results"
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
        tag "Plotting relative lineage abundances of samples"
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path aggregated_file

        output:
        path ("*.tsv")
        path ("*.png"), emit: freyja_plot

        script:
        """
        sed 's/_variants\\.tsv\s*//' '${aggregated_file}' > aggregated.tsv

        freyja plot aggregated.tsv --output freyja-lineage-abundance-plot.png
        """

}