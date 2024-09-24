// enable dsl2
nextflow.enable.dsl=2

//import modules
include { nextclade } from './modules/01-lineageAssignment.nf'
include { bammix } from './modules/02-bammix.nf'
include { bam_filter } from './modules/02-bammix.nf'

// import subworkflows
include { bammix_flagged } from './workflows/bammix_flagged.nf'
include { no_bammix_flag } from './workflows/no_bammix_flag.nf'

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

// HOW IM CURRENTLY RUNNING THE PIPELINE: 
// nextflow run Katmon --in_dir <absolute path no "/" in the end> --out_dir <absolute path no "/" in the end > --run full -profile docker

// added the --run params
        if (params.run == 'full' ) {
            bammix_flagged()
        }
        else if (params.run == 'notfull' ){
            no_bammix_flag()
        }

// SOLUTIONS I'VE TRIED BUT DID NOT WORK:
////////////////////////////////////////////////////////////////
//        ch_count = bam_filter.out.flatMap { it.split("\n") }
//        ch_count_value = ch_count.count()
//       ch_count_value.view()
//        ch_count_value.branch {
//                bammix_flagged: it > 0
//                no_bammix_flag: it <= 0
//                }
//                .set { result }
//        result.bammix_flagged | bammix_flagged()
//        result.no_bammix_flag | no_bammix_flag()
//
//               if ( result.bammix_flagged ) {
//                    bammix_flagged()
//                }
//               else if ( result.no_bammix_flag ) {
//                    no_bammix_flag()
//                }
//////////////////////////////////////////////////////////////
//        ch_count_value = ch_count.count()
//        result = ch_count_value
//                 | branch { n ->
//                   TRUE: n > 0
//                   FALSE: n <= 0
//                   }
//        result.TRUE | bammix_flagged()
//        result.FALSE | no_bammix_flag()
//////////////////////////////////////////////////////////////
// SOLUTIONS TO TRY:
   // it, length of values 
   // convert ch_count to string / array and print out
}