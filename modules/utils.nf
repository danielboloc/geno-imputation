def info (){

    // Intro
    log.info """\
        A R G U M E N T S   U S E D
        ===========================
        vcf                : ${params.vcf}
        ref_assembly       : ${params.ref_assembly}
        outdir             : ${params.outdir}
        eagle_ref_panel    : ${params.eagle_ref_panel}
        token              : ${params.token}
        imputation_job_dir : ${params.imputation_job_dir}
        max_memory         : ${params.max_memory}
        max_cpus           : ${params.max_cpus}
        max_time           : ${params.max_time}
        monochrome_logs    : ${params.monochrome_logs}
        """.stripIndent()


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

def header(){
    // ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";

    return """
    -${c_dim}---------------------------------------------------------------------------------------${c_reset}-
    ${c_purple}   ______                       ____                      __        __  _            ${c_reset}
    ${c_purple}  / ____/__  ____  ____        /  _/___ ___  ____  __  __/ /_____ _/ /_(_)___  ____  ${c_reset}
    ${c_purple} / / __/ _ \\/ __ \\/ __ \\______ / // __ `__ \\/ __ \\/ / / / __/ __ `/ __/ / __ \\/ __ \\ ${c_reset}
    ${c_purple}/ /_/ /  __/ / / / /_/ /_____// // / / / / / /_/ / /_/ / /_/ /_/ / /_/ / /_/ / / / / ${c_reset}
    ${c_purple}\\____/\\___/_/ /_/\\____/     /___/_/ /_/ /_/ .___/\\__,_/\\__/\\__,_/\\__/_/\\____/_/ /_/  ${c_reset}
    ${c_purple}                                         /_/                                         ${c_reset}
    -${c_dim}---------------------------------------------------------------------------------------${c_reset}-
    """.stripIndent()
}