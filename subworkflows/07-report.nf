// enable dsl2
nextflow.enable.dsl=2

// import modules 
include { generation } from '../modules/07-report.nf'

workflow report {
    take:
        plot_bammix
        plot_mutation
        plot_amplicon
        lineage_summary
        freyja_summarized
        freyja_lineage
        virstrain_summary
        amplicon_sorting


        
    main:

// Report generation
   bammix_plot_tsv = plot_bammix
                      .collect() // Collect all paths to bam plots
                      .map{ it.join("\t") + "\n" } // Join paths with tabs and add a new line
                      .collectFile ( name: "bammix_plots.tsv")
   aafplot_mutations_tsv = plot_mutation
                      .collect() 
                      .map{ it.join("\t") + "\n" }
                      .collectFile ( name: "aafplot_mutations.tsv")
   aafplot_amplicons_tsv = plot_amplicon
                      .collect() 
                      .map{ it.join("\t") + "\n" }
                      .collectFile( name: "aafplot_amplicons.tsv")

   generation(
         params.report_r,
         lineage_summary,
         bammix_plot_tsv,
         freyja_summarized,
         freyja_lineage,
         aafplot_mutations_tsv,
         aafplot_amplicons_tsv,
         virstrain_summary,
         amplicon_sorting,
         params.report_rmd )
}