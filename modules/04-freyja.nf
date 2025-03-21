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
        freyja demix ${variants_tsv} ${depth_tsv} --output ${sample}.tsv --depthcutoff 20
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
        path (tsv_demix)

        output:
        path ("*.tsv"), emit: freyja_aggregated_file

        script:
        """
        mkdir Demix_results
        cp ${tsv_demix} Demix_results/
        freyja aggregate Demix_results/ --output aggregated-file.tsv
        """
}

process freyja_plot_summarized {
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
        path ("*.tsv"), emit: aggregated_tsv
        path ("*.png"), emit: freyja_plot_sum

        script:
        """
        sed 's/_variants\\.tsv\s*//' '${aggregated_file}' > aggregated.tsv

        freyja plot aggregated.tsv --output freyja-summarized-plot.png
        """

}

process freyja_plot_lineage {
        tag "Plotting specific lineage abundances of samples"
        container 'staphb/freyja:latest'

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path aggregated_tsv

        output:
        path ("*.png"), emit: freyja_plot_lin

        script:
        """
        freyja plot ${aggregated_tsv} --lineages --output freyja-lineage-specific-plot.png
        """

}

process freyja_list_lineages {
        tag "Getting lineage defining mutations of variants detected by Freyja"
        container 'ufuomababatunde/bammix:v1.1.0' // to fix

        publishDir (
        path: "${params.out_dir}/04-Freyja/Lineages",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        path (tsv_demix)

        output:
        tuple val(tsv_demix.baseName), path ("*.tsv"), emit: freyja_list_lin

        script:
        """
        freyja-list-lineages.py ${tsv_demix.baseName}_lineages.tsv ${tsv_demix}
        """
}

process freyja_get_lineage_def {
        tag "Getting lineage defining mutations of variants detected by Freyja"
        container 'staphb/freyja:1.5.2-03_02_2025-02-03-2025-03-03'

        publishDir (
        path: "${params.out_dir}/04-Freyja/Mutations/${sample}",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val (sample), path (lineage_list)
        path annot
        path ref

        output:
        tuple val (sample), path ("*.tsv"), emit: lin_mut_tsv

        script:
        """
        freyja-get-mutations.sh ${lineage_list} ${sample} ${annot} ${ref}
        """
}

process mutations {
        tag "Creating mutations tsv file"
        container 'ufuomababatunde/bammix:v1.1.0' // to fix

        publishDir (
        path: "${params.out_dir}/04-Freyja/Mutations",
        mode: 'copy',
        overwrite: 'true'
        )

        input:
        tuple val (sample), path (lin_mut_tsv)

        output:
        tuple val (sample), path ("*.tsv"), emit: processed_mut_tsv

        script:
        """
        process-mutations.py ${params.out_dir}/04-Freyja/Mutations/${sample} ${sample}
        """
}