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

    def missingMode = params.missing_genotype == null
        ? ''
        : params.missing_genotype.toString().trim().toLowerCase()

    def missingArgument = missingMode
        ? "--missing-genotype '${missingMode}'"
        : ''

    def defaultPloidy = params.default_ploidy as Integer

    def gvcfMode = (params.gvcf_mode ?: 'auto')
        .toString()
        .toLowerCase()

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

    #
    # Decide whether this input carries hom-ref (reference) blocks.
    # GATK, DRAGEN and bcftools declare a symbolic non-ref ALT in the
    # header. Starling/Strelka-style gVCFs do not, and instead emit
    # ALT="." records carrying INFO/END, so fall back to inspecting the
    # first records of the file.
    #
    case "${gvcfMode}" in
        gvcf)
            is_gvcf="true"
            ;;
        plain)
            is_gvcf="false"

            if grep -qE '^##ALT=<ID=(NON_REF|\\*)[,>]' vcf_header.txt; then
                echo "WARNING: --gvcf_mode plain requested, but ${vcf} declares a symbolic non-ref ALT." >&2
            fi
            ;;
        *)
            is_gvcf="false"

            if grep -qE '^##ALT=<ID=(NON_REF|\\*)[,>]' vcf_header.txt; then
                is_gvcf="true"
            elif grep -q '^##INFO=<ID=END,' vcf_header.txt; then
                n_blocks=\$(
                    { bcftools view -H "${vcf}" 2>/dev/null || true; } \
                    | head -n 100000 \
                    | awk -F '\\t' '
                        (\$5 == "." || \$5 == "<NON_REF>" || \$5 == "<*>") &&
                        \$8 ~ /(^|;)END=/ { n += 1 }
                        END { print n + 0 }
                    '
                )

                if [[ "\${n_blocks:-0}" -gt 0 ]]; then
                    is_gvcf="true"
                fi
            fi
            ;;
    esac

    echo "INFO: ${vcf} treated as gVCF: \${is_gvcf}" >&2

    make_vcf_targets.py \
        --scorefile "${scorefile}" \
        --vcf-header "vcf_header.txt" \
        --output "targets.tsv" \
        --report "target_mapping.tsv"

    #
    # INFO/END is only queryable when declared in the header.
    #
    if grep -q '^##INFO=<ID=END,' vcf_header.txt; then
        end_field='%INFO/END'
    else
        end_field='.'
    fi

    printf 'chrom\\tposition\\tref\\talt\\tend\\tgenotype\\n' \
        > queried_genotypes.tsv

    #
    # --targets-overlap 1 (record) is required so that a reference block
    # starting before a score position is still retrieved.
    #
    if [[ -s targets.tsv ]]; then
        bcftools view \
            --targets-file targets.tsv \
            --targets-overlap 1 \
            -Ou \
            "${vcf}" \
        | bcftools query \
            --samples "\${vcf_sample}" \
            --format "%CHROM\\t%POS\\t%REF\\t%ALT\\t\${end_field}[\\t%GT]\\n" \
        >> queried_genotypes.tsv
    fi

    calculate_prs.R \
        --scorefile "${scorefile}" \
        --genotypes "queried_genotypes.tsv" \
        --sample-id "${meta.id}" \
        --vcf-sample-file "vcf_samples.txt" \
        --source-vcf "${vcf.name}" \
        --input-is-gvcf "\${is_gvcf}" \
        --gvcf-mode "${gvcfMode}" ${missingArgument} \
        --target-build "${target_build}" \
        --default-ploidy "${defaultPloidy}" \
        --strict-alleles "${strictAlleles}" \
        --output-summary "${meta.id}.prs.tsv" \
        --output-details "${meta.id}.variants.tsv"
    """
}
