#!/usr/bin

process aafplot {
    tag "Plotting alternative allele fraction of ${sample}"

    publishDir(
    path: "${params.out_dir}/05-AAFplot",
    mode: 'copy',
    overwrite: 'true'
    )

    input:
    path vcfPath
    
    output:
    tuple val(sample), path("*.png")

    script:
    """
    python aafplot.py [1] [2] [3] [4] [5]
    """
}