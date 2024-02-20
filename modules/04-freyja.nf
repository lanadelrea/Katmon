#!/usr/bin

process freyja {
        tag "Identifying relative lineage abundances from potential mixed SARS-CoV-2 samples using Freyja"

        publishDir (
        path: "${params.out_dir}/04-Freyja",
        mode: 'copy',
        overwrite: 'true'
        )    

        input:
        tuple val(sample), path(fasta)

        output:
        path '*.csv'
        
        script:
        """
        freyja update
        freyja variants ${bamfile} --variants ${sample}_variants.tsv --depths ${sample}_depth.tsv --ref ${reference}
        freyja demix ${sample}_variants.tsv ${sample}_depth.tsv --output ${params.out_dir}/04-Freyja/demix_results
        freyja aggregate ${params.out_dir}/04-Freyja/demix_results --output aggregated-file.tsv
        freyja plot aggregated-file.tsv --output freyja-lineage-abundance-plot.png
        """
}