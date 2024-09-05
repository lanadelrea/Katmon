// enable dsl2
nextflow.enable.dsl=2

//import modules
include { pangolin } from './modules/01-lineageAssignment.nf'
include { nextclade } from './modules/01-lineageAssignment.nf'
include { lineage_assignment } from './modules/01-lineageAssignment.nf'
include { bammix } from './modules/02-bammix.nf'
include { bam_filter } from './modules/02-bammix.nf'
include { makevcf } from './modules/02-bammix.nf'
include { virstrain } from './modules/03-virstrain.nf'
include { virstrain_summary } from './modules/03-virstrain.nf'
include { freyja } from './modules/04-freyja.nf'
include { freyja_demix } from './modules/04-freyja.nf'
include { freyja_aggregate } from './modules/04-freyja.nf'
include { freyja_plot } from './modules/04-freyja.nf'
include { bammixplot } from './modules/05-plots.nf'
include { aafplot_mutations } from './modules/05-plots.nf'
include { aafplot_amplicons } from './modules/05-plots.nf'
include { ampliconsorting_DeltaReads } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_OmicronReads } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_samtools } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_bgzip } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_fasta } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_lineageAssignment_Pangolin } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_lineageAssignment_Nextclade } from './modules/06-ampliconSorting.nf'
include { report } from './modules/07-report.nf'

workflow {
        ch_bam_file = Channel
                    .fromPath("${params.in_dir}/**.bam", type: 'file')
                    .ifEmpty { error "Cannot find any BAM files on ${params.in_dir}"}
        ch_bam_index = Channel
                    .fromPath("${params.in_dir}/**.bai", type: 'file')
                    .ifEmpty { error "Cannot find any BAM index files on ${params.in_dir}"}
        ch_fastq = Channel
                    .fromPath("${params.in_dir}/**.fastq{.gz,}", type: 'file')
                    .ifEmpty { error "Cannot find any fastq files" on ${params.in_dir}}
        ch_fasta = Channel
                    .fromPath("${params.in_dir}/**.fasta", type: 'file')
                    .ifEmpty { error "Cannot find any fasta files" on ${params.in_dir}}

        main:
               ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
               ch_fastq.map { fastqPath -> tuple(fastqPath) }
               ch_cat_fasta = ch_fasta.collectFile(name: 'all_sequences.fasta', newLine: true )

               pangolin(ch_cat_fasta)
               nextclade( ch_cat_fasta, params.SC2_dataset )
               lineage_assignment( pangolin.out.pangolin_csv, nextclade.out.nextclade_tsv )

               bammix ( nextclade.out.nextclade_tsv, ch_bam_file, ch_bam_index ) // TO-DO: will put operator here if there are no flags
               bam_filter ( bammix.out.bammixflagged_csv )

//               bam_filter_result = Channel.of(bam_filter.out.filtered_bam).map() // TO-DO: fixing to accomodate multiple flagged samples

               makevcf( bam_filter.out.filtered_bam.flatten(), params.reference ) 

               virstrain ( ch_fastq )
               virstrain_summary(virstrain.out.txt.collect())
//               ch_virstrain_txt = Channel.fromPath("${params.out_dir}/03-VirStrain", type: 'dir')

//               virstrain_summary( ch_virstrain_txt.collect() )
               
               freyja( ch_bam_file )
               freyja_demix( freyja.out.freyja_variants )
               freyja_aggregate( freyja_demix.out.tsv_demix.collect())
               freyja_plot( freyja_aggregate.out.freyja_aggregated_file )

               bammixplot ( makevcf.out.filtered_vcf)
               aafplot_mutations( makevcf.out.filtered_vcf)
               aafplot_amplicons( aafplot_mutations.out.aafplot_mut)

               ampliconsorting_DeltaReads( makevcf.out.filtered_vcf.collect(), params.jvarkit_jar, params.sort_delta_reads)
               ampliconsorting_OmicronReads( makevcf.out.filtered_vcf.collect(), params.jvarkit_jar, params.sort_omicron_reads)
               ampliconsorting_samtools( ampliconsorting_DeltaReads.out.delta_bam, ampliconsorting_OmicronReads.out.omicron_bam, params.reference)
               ampliconsorting_bgzip( ampliconsorting_samtools.out.vcf.collect())
               ampliconsorting_fasta( ampliconsorting_bgzip.out.vcfgz.collect(), params.reference )
               ampliconsorting_lineageAssignment_Pangolin( ampliconsorting_fasta.out.fasta.collect())
               ampliconsorting_lineageAssignment_Nextclade( ampliconsorting_fasta.out.fasta.collect(), params.SC2_dataset) // TO-DO: finding a way to show amplicon sorting in the final report

               report( lineage_assignment.out.lineageAssign_tsv, bammixplot.out.bammix_plot, freyja_plot.out.freyja_plot, aafplot_mutations.out.aafplot_mut, aafplot_amplicons.out.aafplot_amp, virstrain_summary.out.tsv, params.report_rmd )
               // Surely, there is a much better way to do this but for now you have to deal with this lenghty line of code. I'm just a girl lmao
}

if (params.help) {
        help = """The Katmon pipeline is designed to look for potential Omicron and Delta Co-infection from an NGS run.
                |
                | To run the pipeline, do:
                |       nextflow run Katmon --in_dir <input dir> --out_dir <out dir> -profile <docker or conda>
                |  
                | Required arguments:
                |
                |      --in_dir       Location of the input files containing fasta, fastq, bam and bam index files
                |      --out_dir      Location of the results directory
                | 
                | Optional arguments:
                |      -profile       Can be docker or conda
                |      -resume        To resume the pipeline
                |      -w             The NextFlow work directory. Delete this directory once the process is finished
                |                     Default: ${workDir} 
                |      --help         To view this help page
                |""".stripMargin()
        println(help)
        exit(0)
}