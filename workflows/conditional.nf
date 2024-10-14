// enable dsl2
nextflow.enable.dsl=2

// import subworkflows
include { bammix_flagged } from '../workflows/bammix_flagged.nf'
include { no_bammix_flag } from '../workflows/no_bammix_flag.nf'

workflow conditional {
    take:
        result
    
    main:
    if ( result == true ){
        bammix_flagged()
    }
    else if ( result == false ){
        no_bammix_flag()
    }
}