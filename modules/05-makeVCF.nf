#!/usr/bin 
process makevcf {
    tag "Making vcf file of high quality reads from bam file of ${sample}"
    container 'pegi3s/samtools_bcftools:latest'
    
    publishDir(
        path: "$params.out_dir/05-makeVCF",
        mode: 'copy',
        overwrite: 'true'
    )

    input:
    stdin sample

    output:
    path ("*.vcf.gz")

    script:
    """
    samtools view -b -q 30 ${params.in_dir}/${sample}.bam > ${sample}_filtered.bam
    samtools sort ${sample}_filtered.bam -o ${sample}_filtered.sorted.bam
    samtools index ${sample}_filtered.sorted.bam
    bcftools mpileup -Ou -f ./assets/sars-cov-2/reference-sequence.fasta ${sample}_filtered.sorted.bam | bcftools call -mv -Ob -o ${sample}.vcf.gz
    """
}