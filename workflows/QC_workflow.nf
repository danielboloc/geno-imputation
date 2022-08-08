/*
 * include requires tasks
 */
include { split_norm_phase; chr_concat } from '../modules/vcf.nf'
include { impute; extractZip; indexImputed } from '../modules/imputation.nf'

/*
 * define the data analysis workflow
 */
workflow process_vcf {

    // required inputs
    take:
        chroms_samples_ch

    // workflow implementation
    main:
        split_norm_phase(chroms_samples_ch)
        impute( split_norm_phase.out[0],
                split_norm_phase.out[1],
                split_norm_phase.out[2])
        extractZip( impute.out[0],
                    impute.out[1],
                    impute.out[2].flatMap { n -> n.replaceAll(/\s+$/, "") }) // replace whitespace
        indexImputed( extractZip.out.patient_name,
                        extractZip.out.chromosome )

    // results
    emit:
        indexImputed.out.patient_name
}

workflow concatenate {

    // normalized vcfs
    take:
        patient_dir_w_norm_vcfs

    // perform task
    main:
        chr_concat( patient_dir_w_norm_vcfs )
}
