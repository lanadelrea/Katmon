// enable dsl2
nextflow.enable.dsl=2

//import modules
include { nextclade } from './modules/01-lineageAssignment.nf'
include { bammix } from './modules/02-bammix.nf'
include { bam_filter } from './modules/02-bammix.nf'
include { pangolin } from './modules/01-lineageAssignment.nf'
include { lineage_assignment } from './modules/01-lineageAssignment.nf'

// import subworkflows
include { bammix_flagged } from './workflows/bammix_flagged.nf'
include { no_bammix_flag } from './workflows/no_bammix_flag.nf'
include { conditional } from './workflows/conditional.nf'


workflow {

    ch_bam_file = Channel
                    .fromPath("${params.in_dir}/**.bam", type: 'file')
                    .ifEmpty { error "Cannot find any BAM files on ${params.in_dir}"}
    ch_bam_index = Channel
                    .fromPath("${params.in_dir}/**.bai", type: 'file')
                    .ifEmpty { error "Cannot find any BAM index files on ${params.in_dir}"}
    ch_fastq = Channel
                    .fromPath("${params.in_dir}/**.fastq{.gz,}", type: 'file')
                    .ifEmpty { error "Cannot find any fastq files on ${params.in_dir}"}
    ch_fasta = Channel
                    .fromPath("${params.in_dir}/**.fasta", type: 'file')
                    .ifEmpty { error "Cannot find any fasta files on ${params.in_dir}"}

    main:
        ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
        ch_cat_fasta = ch_fasta.collectFile(name: 'all_sequences.fasta', newLine: true )

        nextclade( ch_cat_fasta, params.SC2_dataset )

        bammix( nextclade.out.nextclade_tsv, ch_bam_file, ch_bam_index )
        bam_filter( bammix.out.bammixflagged_csv )
        ch_count = bam_filter.out.flatMap { it.split("\n") }.count()
        result = ch_count.branch{ n ->
          present: n != 0 || null
          absent: n == 0 || null
        }
        conditional(result.present) | bammix_flagged()
        conditional(result.absent) | not_flagged()
}