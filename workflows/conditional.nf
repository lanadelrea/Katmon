// enable dsl2
nextflow.enable.dsl=2

// import subworkflows
include { bammix_flagged } from '../workflows/bammix_flagged.nf'
include { no_bammix_flag } from '../workflows/no_bammix_flag.nf'

workflow conditional {
    take:
        result
    
    main:
    result.view()
    if ( result != null ){
        bammix_flagged()
    }
    else if ( result == null || 0 ){
        no_bammix_flag()
    }
}