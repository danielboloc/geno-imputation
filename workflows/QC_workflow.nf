/*
 * include requires tasks
 */
include { split_norm_phase } from '../modules/vcf.nf'
//include { split_norm; phase } from '../tasks/vcf_tasks.nf'
//include { impute; extractZip; indexImputed } from '../tasks/imputation.nf'

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
    //     impute( split_norm_phase.out[0],
    //             split_norm_phase.out[1],
    //             split_norm_phase.out[2] )
    //     extractZip( impute.out[0],
    //                 impute.out[1],
    //                 impute.out[2].flatMap { n -> n.replaceAll(/\s+$/, "") } )
    //     indexImputed( extractZip.out.patient_name,
    //                   extractZip.out.chromosome )

	// emit:
    //     indexImputed.out.patient_name
}

// workflow concatenate {
//     take:
//         patient_dir_w_norm_vcfs
//     main:
//         chr_concat( patient_dir_w_norm_vcfs )
// }
