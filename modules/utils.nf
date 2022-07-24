def info (){
    def summary = [:]
    summary['Pipeline Name']    = 'geno-imputation'
    summary['Max Memory']       = params.max_memory
    summary['Max CPUs']         = params.max_cpus
    summary['Max Time']         = params.max_time
    summary['Output dir']       = params.outdir
    summary['Working dir']      = workflow.workDir
    summary['Script dir']       = workflow.projectDir
    summary['Current path']     = "$PWD"

    log.info "========================================================"
    log.info "SUMMARY OF THE RUN                                      "
    log.info "------------------                                      "
    log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
    log.info "========================================================"
    log.info ""
}

process create_dir_structure {
    input:
    tuple val(patient), file(vcf), file(csi)

    """
    mkdir -p $baseDir/results/{concat,helper}/
    mkdir -p $baseDir/results/$patient/{concat,helper}/
    mkdir -p $baseDir/results/$patient/imputation/{raw,proc,concat,helper}/
    """
}

