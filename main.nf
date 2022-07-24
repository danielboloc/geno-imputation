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
         vcf          : ${params.vcf}
         ref_assembly : ${params.ref_assembly}
         outdir       : ${params.outdir}
         """
         .stripIndent()

// LOG HELP
if(params.help){
    log.info ""
    log.info "--------------------------------------------------------"
    log.info "USAGE                                                       "
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--eagle_ref_panel  FILE       chromosome reference file   "
    log.info "--ref_assembly     FILE       FASTA reference file   "
    log.info "--------------------------------------------------------"
    log.info ""
}

// Include
include { info; create_dir_structure } from './modules/utils.nf'
include { process_vcf; concatenate } from './workflows/QC_workflow.nf'
include { concatImputed } from './modules/imputation.nf'

// Chanels
// Get sample name
getPatient = """bcftools query -l ${params.vcf}"""
params.patient = getPatient.execute().text.replaceAll(/\s+$/, "")
channel
    .fromPath(params.vcf, checkIfExists:true)
    .map{ it -> tuple( params.patient, file(params.vcf), file(params.vcf + '.csi') ) }
    .set{ samples_ch }

// Testing only chr1 and chr2
chromosomes = channel.of(1..2)
chroms_persons_ch = chromosomes
    .combine(samples_ch)

workflow {

    // Print pipeline info
    info()

    // Create output folder structure
    create_dir_structure(samples_ch)

    // Filter, normalize and phase vcf then impute by chromosomes
    process_vcf(chroms_persons_ch)

    // Concatenate chromosomes from normalized, phased chromosomes
    concatenate( process_vcf.out[0] )

    // Concatenate imputed chromosomes
    concatImputed( process_vcf.out[0] )

}