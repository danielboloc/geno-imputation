process split_norm_phase {

    label 'big_mem'
    publishDir "${params.outdir}/results/${patient}", overwrite: true, mode:'copy'

    input:
    tuple val(chromosome), val(patient), file(vcf), file(csi)

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
        --threads ${task.cpus} \
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
        --vcfRef ${params.eagle_ref_panel}/chr${chromosome}.bcf.gz \
        --vcfTarget chr${chromosome}_norm.vcf.gz \
        --geneticMapFile ${params.eagle_ref_panel}/genetic_map_hg38_withX.txt \
        --outPrefix chr${chromosome}_phased \
        --allowRefAltSwap \
        --vcfOutFormat z \
        --keepMissingPloidyX \
        --numThreads ${task.cpus}
    bcftools index -f chr${chromosome}_phased.vcf.gz
    """
}

process chr_concat {
    publishDir "${params.outdir}/results/${patient}", overwrite: true, mode:'copy'

    input:
    val(patient)

    output:
    path("${patient}_norm.bcf.gz")
    path("${patient}_norm.bcf.gz.csi")
    path("${patient}_phased.bcf.gz")
    path("${patient}_phased.bcf.gz.csi")

    """
    # concatenate normalized
    bcftools concat -a -d none ${params.outdir}/results/${patient}/chr*_norm.vcf.gz --threads 6 -Ob -o ${patient}_norm.bcf.gz;
    bcftools index -f ${patient}_norm.bcf.gz;
    # concatenate phased
    bcftools concat -a -d none ${params.outdir}/results/${patient}/chr*_phased.vcf.gz --threads 6 -Ob -o ${patient}_phased.bcf.gz;
    bcftools index -f ${patient}_phased.bcf.gz;
    """
}
