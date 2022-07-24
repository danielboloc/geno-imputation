/*
 * Genotype imputation pipeline using Michigan Imputation Server locally
 * Using Nextflow DSL2
 *
 * Author: Daniel Boloc
 * E-mail: danielboloc@gmail.com
 */
nextflow.enable.dsl=2

// LOG HELP
if(params.help){
    log.info ""
    log.info "----------------------------------------------------------------------------------------"
    log.info "USAGE"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "-------------------"
    log.info "--vcf                 FILE    Input .vcf file (it need to have .csi index from bcftools)"
    log.info "--ref_assembly        FILE    FASTA reference file                                      "
    log.info "--outdir              PATH    where output should be stored                             "
    log.info "--eagle_ref_panel     PATH    where chromosome reference files are located              "
    log.info "--token               STR     Michigan Imputation Server (MIS) local token              "
    log.info "--imputation_job_dir  PATH    where MIS 'data' folder is located on disk                "
    log.info ""
    log.info "Optional arguments:"
    log.info "------------------"
    log.info "--max_memory          STR     specifies memory to allocate                              "
    log.info "--max_cpus            INT     specifies number of logical cpus required                 "
    log.info "--max_time            FILE    how long a process is allowed to run                      "
    log.info "--monochrome_logs     FILE    whether to display logs using just one color (white)      "
    log.info "----------------------------------------------------------------------------------------"
    log.info ""
    exit(1)
}

// Include
include { info; header; create_dir_structure } from './modules/utils.nf'
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

    // Print header and pipeline info
    log.info header()
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