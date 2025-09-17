// enable dsl2
nextflow.enable.dsl=2

// import modules
include { bammixplot } from '../modules/06-plots.nf'
include { aafplots } from '../modules/06-plots.nf'

workflow plots {
    take:
        bammix_csv
        filtered_vcf
        flagged_bams
        mutations_tsv


    main:
        // Alternative Allele Fraction (AAF) plotting
        bammixplot(bammix_csv)
    
        mutations_tsv
           .join(filtered_vcf) //Combine the input for aafplots
           .set { aafplot_inputs }

        all_files_flagged = aafplot_inputs.join.flagged_bams.view() // Join the mutations tsv and bam file of flagged in one input

        aafplots(params.primer_scheme,
                 all_files_flagged)

    emit:
        plot_bammix   = bammixplot.out.bammix_plot
        plot_mutation = aafplots.out.aafplot_mut
        plot_amplicon = aafplots.out.aafplot_amp
}