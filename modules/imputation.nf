process impute {

    errorStrategy 'retry'
    maxRetries 3

    input:
    val(patient)
    val(chromosome)
    file(phased_chr)

    output:
    val(patient)
    val(chromosome)
    stdout emit: job_id

    script:
    """
    # If error with permission denied appear can fix with 'chmod +x impute_on_MIS.py'
    impute_on_MIS.py \
        --file ${phased_chr} \
        --token ${params.token} \
        --patient ${patient}
    """
}

process extractZip {

    input:
    val(patient)
    val(chromosome)
    val(job_id)

    output:
    val(patient), emit: patient_name
    val(chromosome), emit: chromosome

    script:
    """
    7z x ${params.imputation_job_dir}/${job_id}/local/chr_${chromosome}.zip -aoa -p'password' -o${params.outdir}/results/${patient}/imputation/raw/
    """
}

process indexImputed {

    input:
    val(patient)
    val(chromosome)

    output:
    val(patient), emit: patient_name

    """
    bcftools index -f ${params.outdir}/results/${patient}/imputation/raw/chr${chromosome}.dose.vcf.gz
    """
}

process concatImputed {

    input:
    val(patient)

    output:
    val(patient)

    """
    bcftools concat -a -d none ${params.outdir}/results/${patient}/imputation/raw/chr*.dose.vcf.gz -Ob -o ${params.outdir}/results/${patient}/imputation/concat/${patient}_imputed.bcf.gz;
    bcftools index -f ${params.outdir}/results/${patient}/imputation/concat/${patient}_imputed.bcf.gz;
    """
}
