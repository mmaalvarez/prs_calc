process CALCULATE_PRS {

    tag "${meta.id}"

    label 'process_medium'
    label 'prs_calc_env'
    
    input:
    tuple val(meta), path(vcf), path(scorefile)
    val target_build

    output:
    tuple val(meta),
          path("${meta.id}.prs.tsv"),
          path("${meta.id}.variants.tsv"),
          emit: results

    script:
    def strictAlleles = params.strict_alleles
        .toString()
        .toLowerCase()

    def missingMode = params.missing_genotype
        .toString()
        .toLowerCase()

    def defaultPloidy = params.default_ploidy as Integer

    """
    bcftools query -l "${vcf}" > vcf_samples.txt

    n_samples=\$(awk 'NF { n += 1 } END { print n + 0 }' vcf_samples.txt)

    if [[ "\${n_samples}" -ne 1 ]]; then
        echo "ERROR: ${vcf} contains \${n_samples} samples." >&2
        echo "Each input VCF/BCF must contain exactly one sample." >&2
        exit 1
    fi

    vcf_sample=\$(head -n 1 vcf_samples.txt)

    bcftools view -h "${vcf}" > vcf_header.txt

    make_vcf_targets.py \
        --scorefile "${scorefile}" \
        --vcf-header "vcf_header.txt" \
        --output "targets.tsv" \
        --report "target_mapping.tsv"

    printf 'chrom\\tposition\\tref\\talt\\tgenotype\\n' \
        > queried_genotypes.tsv

    if [[ -s targets.tsv ]]; then
        bcftools view \
            --targets-file targets.tsv \
            --targets-overlap 0 \
            -Ou \
            "${vcf}" \
        | bcftools query \
            --samples "\${vcf_sample}" \
            --format '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT]\\n' \
        >> queried_genotypes.tsv
    fi

    calculate_prs.R \
        --scorefile "${scorefile}" \
        --genotypes "queried_genotypes.tsv" \
        --sample-id "${meta.id}" \
        --vcf-sample-file "vcf_samples.txt" \
        --source-vcf "${vcf.name}" \
        --target-build "${target_build}" \
        --missing-genotype "${missingMode}" \
        --default-ploidy "${defaultPloidy}" \
        --strict-alleles "${strictAlleles}" \
        --output-summary "${meta.id}.prs.tsv" \
        --output-details "${meta.id}.variants.tsv"
    """
}
