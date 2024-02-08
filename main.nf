// enable dsl2
nextflow.enable.dsl=2

//import modules

include { concat } from './modules/01-lineageAssignment.nf'
include { pangolin } from './modules/01-lineageAssignment.nf'
include { nextclade } from './modules/01-lineageAssignment.nf'
include { bammix } from './modules/02-bammix.nf'

workflow {
        ch_bam_file = Channel
                    .fromPath("${params.in_dir}/*.bam", type: 'file')
                    .ifEmpty { error "Cannot find any BAM files on ${params.in_dir}"}
                    .view()
        ch_bam_index = Channel
                    .fromPath("${params.in_dir}/*.bai", type: 'file')
                    .ifEmpty { error "Cannot find any BAM index files on ${params.in_dir}"}
                    .view()

        main:
               concat()
               pangolin( concat.out.fasta )
               nextclade( concat.out.fasta )
               bammix ( nextclade.out.nextclade_tsv, ch_bam_file, ch_bam_index )
               virstrain ( concat.out.fasta, Virstrain_DB )
}