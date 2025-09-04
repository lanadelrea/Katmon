// enable dsl2
nextflow.enable.dsl=2

//import modules
include { pangolin } from './modules/01-lineageAssignment.nf'
include { nextclade } from './modules/01-lineageAssignment.nf'
include { lineage_assignment } from './modules/01-lineageAssignment.nf'

include { virstrain } from './modules/03-virstrain.nf'
include { virstrain_summary } from './modules/03-virstrain.nf'

//include { bammix } from './modules/02-bammix.nf'
include { bammix_positions } from './modules/02-bammix.nf'
include { bammix_process } from './modules/02-bammix.nf'
include { bammix_flagged_positions } from './modules/02-bammix.nf'
include { bammix_flagged_samples } from './modules/02-bammix.nf'
include { makevcf } from './modules/02-bammix.nf'
include { bcftools } from './modules/02-bammix.nf'

include { freyja } from './modules/04-freyja.nf'
include { freyja_demix } from './modules/04-freyja.nf'
include { freyja_aggregate } from './modules/04-freyja.nf'
include { freyja_plot_summarized } from './modules/04-freyja.nf'
include { freyja_plot_lineage } from './modules/04-freyja.nf'
include { freyja_list_lineages } from './modules/04-freyja.nf'
include { freyja_get_lineage_def } from './modules/04-freyja.nf'

include { mutations } from './modules/04-freyja.nf'

include { bammixplot } from './modules/05-plots.nf'
include { aafplots } from './modules/05-plots.nf'

include { get_pos_mut } from './modules/06-ampliconSorting.nf'
include { ampliconsorting } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_consensus } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_renamefasta } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_pangolin } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_nextclade } from './modules/06-ampliconSorting.nf'
include { ampliconsorting_table } from './modules/06-ampliconSorting.nf'

include { report } from './modules/07-report.nf'
//include { report_no_flag } from './modules/07-report.nf'


// import subworkflows
//include { bammix_flagged } from './workflows/bammix_flagged.nf'
//include { no_bammix_flag } from './workflows/no_bammix_flag.nf'
//include { conditional } from './workflows/conditional.nf'


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

        // Input bam, fastq, and fasta files
        ch_bam_file.map { bamfilePath -> tuple(bamfilePath) }
        ch_bam_index.map { baifilePath -> tuple(baifilePath)}
        ch_cat_fasta = ch_fasta.collectFile(name: 'all_sequences.fasta', newLine: true )
        ch_fastq.map { fastqPath -> tuple(fastqPath) }

        // Lineage assignment
           // through Pangolin and Nextclade using FASTA files
           pangolin( ch_cat_fasta )
           nextclade( ch_cat_fasta, params.SC2_dataset )
           lineage_assignment( pangolin.out.pangolin_csv, nextclade.out.nextclade_tsv )
           // through VirStrain using FASTQ files
           virstrain(params.virstrain_database, ch_fastq)
           virstrain_summary(params.virstrain_txt_dir, virstrain.out.txt)

        // Detecting of nucleotide mixtures from all samples
//           bammix( nextclade.out.nextclade_tsv)
           bammix_positions( nextclade.out.nextclade_tsv )
           snps = (bammix_positions.out).toList()
           bammix_process( ch_bam_file, ch_bam_index, snps)
           bammix_flagged_positions( bammix_process.out.bammix_csv )
           bammix_flagged_samples( bammix_flagged_positions.out.bammix_flags_csv.collect() )

        // Filter high quality reads from samples with nucleotide mixture and make VCF
           ch_samples_bammix_flagged = bammix_flagged_samples.out.samples_txt
               .splitText()
               .map { it.trim() }

           ch_bam_bammix_flagged = ch_bam_file
               .map { bam -> 
                   def sample_name = bam.baseName 
                   tuple(sample_name, bam)
                }

           ch_samples_bammix_flagged_keyed = ch_samples_bammix_flagged.map { s -> tuple(s) }

           ch_flagged_bams = ch_samples_bammix_flagged_keyed
               .join(ch_bam_bammix_flagged)
               .map { sample, bam -> tuple(sample, bam) }

           makevcf( ch_flagged_bams, params.reference )
           bcftools(makevcf.out.mpileup)

        // Lineage abundance estimation per sample through Freyja using BAM files
           freyja(params.reference, ch_bam_file)
           freyja_demix(freyja.out.freyja_variants)
           freyja_aggregate(freyja_demix.out.tsv_demix.collect())
           freyja_plot_summarized(freyja_aggregate.out.freyja_aggregated_file)
           freyja_plot_lineage(freyja_plot_summarized.out.aggregated_tsv)

           // Only for samples flagged by bammix
           ch_freyja_tsv = freyja_demix.out.tsv_demix.collect()
           ch_freyja_samples = ch_freyja_tsv.flatten()
               .map { tsv ->
                   def sample_name = tsv.baseName
                   tuple(sample_name, tsv)
               }


           ch_freyja_flagged_tsv = ch_samples_bammix_flagged_keyed
               .join (ch_freyja_samples)
               .map { sample, tsv -> tuple (sample, tsv) }
           
//           freyja_list_lineages(freyja_demix.out.tsv_demix)
           freyja_list_lineages(ch_freyja_flagged_tsv)
           freyja_get_lineage_def(freyja_list_lineages.out.freyja_list_lin, params.annot, params.ref)

           freyja_result = freyja_get_lineage_def.out.lin_mut_tsv

           mutations(freyja_result)
        
           // Alternative Allele Fraction (AAF) plotting
           bammixplot(bammix_process.out.bammix_csv)
           //Combine the input for aafplots
           mutations.out.processed_mut_tsv
               .join(bcftools.out.filtered_vcf)
               .set { aafplot_inputs }
           aafplots(params.primer_scheme, aafplot_inputs)

           // Amplicon sorting
           bcftools.out.filtered_vcf
               .join(freyja_get_lineage_def.out.lin_mut_tsv) // only get mutations list for bammix flagged samples
               .set { get_posmut_flagged }
           get_pos_mut( get_posmut_flagged ) // List mutations for the two most abundant lineages

           get_pos_mut.out.pos_mut_lineage_A
               .join(get_pos_mut.out.pos_mut_lineage_B)
               .set { pos_mut }
           ampliconsorting(pos_mut, params.jvarkit_jar, params.sort_reads)

           ampliconsorting.out.sorted_bam_lineage_A
               .join(ampliconsorting.out.sorted_bam_lineage_B)
               .set { ampsort_bam }
           ampliconsorting_consensus(ampsort_bam, params.reference)

           ampliconsorting_consensus.out.consensus_lineage_A
               .join(ampliconsorting_consensus.out.consensus_lineage_B)
               .set { ampsort_consensus }
           ampliconsorting_renamefasta(ampsort_consensus)
           ampliconsorting_renamefasta.out.ampsort_consensus_final
               .collectFile(name: 'all_sequences_ampsort.fasta', newLine: true )
               .set { ch_ampsort_cat }

           ampliconsorting_pangolin(ch_ampsort_cat)
           ampliconsorting_nextclade(ch_ampsort_cat, params.SC2_dataset)

           ampliconsorting_table(ampliconsorting_pangolin.out.ampsort_pangolin_csv, ampliconsorting_nextclade.out.ampsort_nextclade_tsv)


        // IF ch_count = 0 ; No samples flagged by bammix 

        // Report generation
           bammix_plot_tsv = bammixplot.out.bammix_plot
                              .collect() // Collect all paths to bam plots
                              .map{ it.join("\t") + "\n" } // Join paths with tabs and add a new line
                              .collectFile ( name: "bammix_plots.tsv")
           aafplot_mutations_tsv = aafplots.out.aafplot_mut
                              .collect() 
                              .map{ it.join("\t") + "\n" }
                              .collectFile ( name: "aafplot_mutations.tsv")
           aafplot_amplicons_tsv = aafplots.out.aafplot_amp
                              .collect() 
                              .map{ it.join("\t") + "\n" }
                              .collectFile( name: "aafplot_amplicons.tsv")

           report(
                 params.report_r,
                 lineage_assignment.out.lineageAssign_tsv,
                 bammix_plot_tsv,
                 freyja_plot_summarized.out.freyja_plot_sum,
                 freyja_plot_lineage.out.freyja_plot_lin,
                 aafplot_mutations_tsv,
                 aafplot_amplicons_tsv,
                 virstrain_summary.out.tsv,
                 ampliconsorting_table.out.ampsort_table,
                 params.report_rmd )
}

if ( params.help ) {
         help = """The Katmon pipeline is designed to look for potential SARS-CoV-2 Co-infection from an NGS run.
                 
                  To run the pipeline, do:
                        nextflow run Katmon --in_dir <input dir> --out_dir <out dir>
                   
                  Required arguments:
                 
                       --in_dir             Input directory containing FASTA, FASTQ, BAM, and BAM index files
                       --out_dir            Output directory for results
                  
                  Optional arguments:
                       --bammix_thresh      Set the bammix threshold for the proportion of the major allele
                                                    Default: 0.8
                       -profile             Can be docker or conda
                       -resume              To resume the pipeline
                       -w                   The NextFlow work directory. Delete this directory once the process is finished
                                                    Default: ${workDir} 
                       --help               To view this help message
                 """.stripMargin()
         println(help)
         exit(0)
}

// TO-DO:
   // Add branch for when there are no samples flagged by bammix for nucleotide mixtures/clean batch of samples.
   // Create one env/docker image for the main processes

// I know this is a very loooong main.nf
// I'm doing my best hahahu (╥ᆺ╥；)
// Will write this in a better way, I promise!