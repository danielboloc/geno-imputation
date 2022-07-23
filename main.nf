/*
 * Genotype imputation pipeline using Michigan Imputation Server locally
 * Using Nextflow DSL2
 *
 * Author: Daniel Boloc
 * E-mail: danielboloc@gmail.com
 */
nextflow.enable.dsl=2

// Intro
log.info """\
         G E N O T Y P E - I M P U T A T I O N    P I P E L I N E
         ========================================================
         ""                                                                   ""
         "---------------------------------------------------------------------"
         " Geno-Imputation v1.0.0: Imputation of a dataset against a reference "
         "---------------------------------------------------------------------"
         ""                                                                   ""
         vcf          : ${params.vcf}
         ref_assembly : ${params.ref_assembly}
         outdir       : ${params.outdir}
         """
         .stripIndent()

// LOG HELP
if(params.help){
    log.info ""
    log.info "-------------------------------------------------------------"
    log.info " USAGE                                                       "
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--eagle-ref-panel  FILE       chromosome reference file   "
    log.info "--ref-assembly     FILE       FASTA reference file   "
    log.info "-------------------------------------------------------------"
}

include { info } from './modules/utils.nf'

workflow {

    info()
}