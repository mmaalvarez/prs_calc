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
        .trim()
        .toLowerCase()

    def noCallsMode = params.no_calls
        .toString()
        .trim()
        .toLowerCase()

    def minCoveredArgument = params.min_covered_fraction == null
        ? ''
        : "--min_covered_fraction '${params.min_covered_fraction}'"

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


    # Detect using the header and, when needed, ALL records. Checking
    # only the first 100,000 records can miss later gVCF blocks.
    #
    # Operational definition:
    #   gvcf    : symbolic non-ref marker, or ALT=. with INFO/END
    #   plain   : no gVCF marker anywhere, but at least one concrete/
    #             other non-placeholder ALT record
    #   unknown : neither kind of evidence
    #
    if grep -qE '^##ALT=<ID=(NON_REF|\\*)[,>]' vcf_header.txt; then
        detected_format="gvcf"
    else
        detected_format=\$(
            bcftools view -H "${vcf}" |
                awk -F '\\t' '
                    {
                        n = split(\$5, alleles, ",")
                        for (i = 1; i <= n; i++) {
                            if (alleles[i] == "<NON_REF>" ||
                                alleles[i] == "<*>") {
                                has_gvcf = 1
                            } else if (alleles[i] != ".") {
                                has_plain = 1
                            }
                        }

                        if (\$4 ~ /^[ACGTNacgtn]\$/ &&
                            \$5 == "." &&
                            \$8 ~ /(^|;)END=[0-9]+(;|\$)/) {
                            has_gvcf = 1
                        }
                    }
                    END {
                        if (has_gvcf) print "gvcf"
                        else if (has_plain) print "plain"
                        else print "unknown"
                    }
                '
        )
    fi

    case "${gvcfMode}" in
        auto)
            if [[ "\${detected_format}" == "unknown" ]]; then
                echo "ERROR: Cannot determine whether ${vcf} is a plain VCF or gVCF. Specify the input format explicitly with --gvcf_mode plain or --gvcf_mode gvcf." >&2
                exit 1
            fi
            effective_format="\${detected_format}"
            ;;
        gvcf|plain)
            if [[ "\${detected_format}" != "unknown" &&
                  "\${detected_format}" != "${gvcfMode}" ]]; then
                echo "ERROR: Specified input format (--gvcf_mode ${gvcfMode}) does not match the detected actual input format (\${detected_format}) for ${vcf}." >&2
                exit 1
            fi

            effective_format="${gvcfMode}"

            if [[ "\${detected_format}" == "unknown" ]]; then
                echo "WARNING: The format of ${vcf} cannot be independently determined; trusting the explicit --gvcf_mode ${gvcfMode} declaration." >&2
            fi
            ;;
    esac

    is_gvcf="false"
    if [[ "\${effective_format}" == "gvcf" ]]; then
        is_gvcf="true"
    fi

    if [[ "\${is_gvcf}" == "true" &&
          "${missingMode}" == "reference" ]]; then
        echo "ERROR: --missing_genotype reference is allowed only for plain VCF input; ${vcf} is being treated as a gVCF." >&2
        exit 1
    fi

    echo "INFO: ${vcf}: detected=\${detected_format}; effective=\${effective_format}" >&2

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
    # --targets_overlap 1 (record) is required so that a reference block
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
        --sample_id "${meta.id}" \
        --vcf_sample_file "vcf_samples.txt" \
        --source_vcf "${vcf.name}" \
        --input_is_gvcf "\${is_gvcf}" \
        --gvcf_mode "${gvcfMode}" \
        --missing_genotype "${missingMode}" \
        --no_calls "${noCallsMode}" ${minCoveredArgument} \
        --target_build "${target_build}" \
        --default_ploidy "${defaultPloidy}" \
        --strict_alleles "${strictAlleles}" \
        --output_summary "${meta.id}.prs.tsv" \
        --output_details "${meta.id}.variants.tsv"
    """
}
