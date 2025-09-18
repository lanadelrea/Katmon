// enable dsl2
nextflow.enable.dsl=2

//import subworkflows
include { lineage_assignment } from '../subworkflows/01-lineage-assignment.nf'
include { bammix             } from '../subworkflows/02-bammix.nf'
include { virstrain          } from '../subworkflows/03-virstrain.nf'
include { freyja             } from '../subworkflows/04-freyja.nf'
include { amplicon_sorting   } from '../subworkflows/05-amplicon-sorting.nf'
include { plots              } from '../subworkflows/06-plots.nf'
include { report             } from '../subworkflows/07-report.nf'

workflow ont {

    ch_bam = Channel
                    .fromPath("${params.indir}/**.bam", type: 'file')
                    .ifEmpty { error "Cannot find any BAM files on ${params.indir}"}
    ch_bai = Channel
                    .fromPath("${params.indir}/**.bai", type: 'file')
                    .ifEmpty { error "Cannot find any BAM index files on ${params.indir}"}
    ch_fastq = Channel
                    .fromPath("${params.indir}/**.{fastq,fq}{,.gz}", type: 'file')
                    .ifEmpty { error "Cannot find any fastq files on ${params.indir}"}
    ch_fasta = Channel
                    .fromPath("${params.indir}/**.fasta", type: 'file')
                    .ifEmpty { error "Cannot find any fasta files on ${params.indir}"}

    main:

    // Input bam, fasta, and fastq files
        ch_bam_file = ch_bam
                 .map { bam ->
                   def sample_name = bam.baseName
                   tuple(sample_name, bam)
                }

        ch_bam_index = ch_bai
                 .map { bai ->
                   def sample_name = bai.baseName
                   tuple(sample_name, bai)
                }

        ch_fastq_files = ch_fastq
                 .map { fastq -> 
                   def sample_name = fastq.baseName 
                   tuple(sample_name, fastq) // Create tuple of sample name and path to the input fastq files
                }

        ch_cat_fasta = ch_fasta.collectFile(name: 'all_sequences.fasta', newLine: true )

    // Analysis pipeline
        lineage_assignment ( ch_cat_fasta )
        bammix ( 
            lineage_assignment.out.nextclade_tsv, 
            ch_bam_file,
            ch_bam_index )
        virstrain ( ch_fastq_files )
        freyja ( 
            ch_bam_file, 
            bammix.out.flagged )
        amplicon_sorting ( 
            bammix.out.filtered_vcf,
            bammix.out.flagged_bams,
            freyja.out.freyja_result )
        plots ( 
            bammix.out.bammix_csv, 
            bammix.out.filtered_vcf,
            bammix.out.flagged_bams, 
            freyja.out.mutations_tsv )
        report (
            plots.out.plot_bammix,
            plots.out.plot_mutation,
            plots.out.plot_amplicon,
            lineage_assignment.out.lineage_summary,
            freyja.out.freyja_summarized,
            freyja.out.freyja_lineage,
            virstrain.out.virstrain_summary,
            amplicon_sorting.out.table)
}