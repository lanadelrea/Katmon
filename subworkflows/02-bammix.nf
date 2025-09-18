// enable dsl2
nextflow.enable.dsl=2

// import modules
include { positions } from '../modules/02-bammix.nf'
include { bammix_process } from '../modules/02-bammix.nf'
include { flagged_positions } from '../modules/02-bammix.nf'
include { flagged_samples } from '../modules/02-bammix.nf'
include { makevcf } from '../modules/02-bammix.nf'
include { bcftools } from '../modules/02-bammix.nf'

workflow bammix {
    take:
        nextclade_tsv
        ch_bam_file
        ch_bam_index

    main:
        // Looking for positions with nucleotide mixtures
        positions( nextclade_tsv )
            snps = (positions.out).toList()

        ch_bam_bai = ch_bam_file.join(ch_bam_index) // Join bam and bai channel by sample name

        bammix_process( ch_bam_bai, snps)
        flagged_positions( bammix_process.out.bammix_csv )
        flagged_samples( flagged_positions.out.bammix_flags_csv.collect() )

        // Filter high quality reads from samples with nucleotide mixture and make VCF
        ch_samples_bammix_flagged = flagged_samples.out.samples_txt
            .splitText()
            .map { it.trim() }

        ch_samples_bammix_flagged_keyed = ch_samples_bammix_flagged.map { s -> tuple(s) }

        ch_samples_bammix_flagged_keyed.view()

        ch_flagged_bams = ch_samples_bammix_flagged_keyed.join(ch_bam_file)

//        ch_flagged_bams = ch_samples_bammix_flagged_keyed
//             .join(ch_bam_file)
//             .map { sample, bam -> tuple(sample, bam) }

        makevcf( ch_flagged_bams, params.reference )
        bcftools(makevcf.out.mpileup)

    emit:
        flagged = ch_samples_bammix_flagged_keyed
        bammix_csv = bammix_process.out.bammix_csv
        filtered_vcf = bcftools.out.filtered_vcf
        flagged_bams = ch_flagged_bams
}