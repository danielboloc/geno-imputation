process split_norm_phase {

	label 'big_mem'
	publishDir "${baseDir}/results/${patient}"

	input:
	tuple val(chromosome), val(patient), file(vcf), file(csi), file(gvcf)

    output:
    val(patient)
    val(chromosome)
    path("chr${chromosome}_phased.vcf.gz")
    path("chr${chromosome}_phased.vcf.gz.csi")
    path("chr${chromosome}_norm.vcf.gz")
    path("chr${chromosome}_norm.vcf.gz.csi")

	"""

	bcftools view -r chr${chromosome} ${vcf} -Ob -o chr${chromosome}.bcf.gz
	bcftools index -f chr${chromosome}.bcf.gz

    # Filter PASS, DP>10, GQ>20, GT not missing
    bcftools view -i'FORMAT/GQ>20 & FORMAT/DP>10 & GT!~"mis" & FILTER~"PASS"' \
        chr${chromosome}.bcf.gz \
        --threads 7 \
        -Ob \
        -o chr${chromosome}_only_PASS_DPmore10_GQmore20.bcf.gz
    # Normalize
	bcftools norm \
        -m \
        -any \
        --check-ref w \
        -f ${params.ref_assembly} chr${chromosome}_only_PASS_DPmore10_GQmore20.bcf.gz \
	| bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' \
	| bcftools norm --rm-dup all --threads 6 -Oz -o chr${chromosome}_norm.vcf.gz
	bcftools index -f chr${chromosome}_norm.vcf.gz

    # Phasing for imputation
    eagle \
	    --vcfRef ${params.coverage_30_dir}/bcf/chr${chromosome}.bcf.gz \
        --vcfTarget chr${chromosome}_norm.vcf.gz \
        --geneticMapFile ${params.coverage_30_dir}/genetic_map_hg38_withX.txt \
        --outPrefix chr${chromosome}_phased \
        --allowRefAltSwap \
        --vcfOutFormat z \
        --keepMissingPloidyX \
        --numThreads 6
    bcftools index -f chr${chromosome}_phased.vcf.gz
	"""
}