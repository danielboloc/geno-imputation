def info (){
    def summary = [:]
    summary['Pipeline Name']    = 'geno-imputation'
    summary['Max Memory']       = params.max_memory
    summary['Max CPUs']         = params.max_cpus
    summary['Max Time']         = params.max_time
    summary['Output dir']       = params.outDir
    summary['Working dir']      = workflow.workDir
    summary['Script dir']       = workflow.projectDir
    summary['Current path']     = "$PWD"

    log.info summary.collect { k,v -> "${k.padRight(20)}: $v" }.join("\n")
    log.info "========================================================"
    log.info ""
}