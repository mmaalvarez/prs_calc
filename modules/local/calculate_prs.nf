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

    tuple val(meta),
          path("${meta.id}.target_mapping.tsv"),
          emit: mapping

    script:
    def q = { value ->
        "'" + value.toString().replace("'", "'\"'\"'") + "'"
    }

    def vcfArg = q.call(vcf)
    def scoreArg = q.call(scorefile)
    def sampleIdArg = q.call(meta.id)
    def sourceVcfArg = q.call(vcf.name)
    def buildArg = q.call(target_build)
    def mappingFile = "${meta.id}.target_mapping.tsv"

    def strictAlleles = params.strict_alleles
        .toString().trim().toLowerCase()
    def missingMode = params.missing_genotype
        .toString().trim().toLowerCase()
    def noCallsMode = params.no_calls
        .toString().trim().toLowerCase()
    def gvcfMode = (params.gvcf_mode ?: 'auto')
        .toString().trim().toLowerCase()
    def defaultPloidy = params.default_ploidy as Integer
    def genotypeCalls = params.genotype_calls
        .toString().trim().toLowerCase()
    def nonAdditive = params.non_additive
        .toString().trim().toLowerCase()
    def minCoveredArgument = params.min_covered_fraction == null
        ? ''
        : "--min_covered_fraction ${q.call(params.min_covered_fraction)}"

    """
    bcftools query -l ${vcfArg} > vcf_samples.txt

    n_samples=\$(awk 'NF { n += 1 } END { print n + 0 }' vcf_samples.txt)
    if [[ "\${n_samples}" -ne 1 ]]; then
        printf 'ERROR: %s contains %s samples; expected exactly one.\\n' \
            ${vcfArg} "\${n_samples}" >&2
        exit 1
    fi

    vcf_sample=\$(head -n 1 vcf_samples.txt)
    if [[ "\${vcf_sample}" == "NA" ]]; then
        echo "ERROR: literal VCF sample ID 'NA' is reserved." >&2
        exit 1
    fi

    bcftools view -h ${vcfArg} > vcf_header.txt

    make_vcf_targets.py \
        --scorefile ${scoreArg} \
        --vcf-header vcf_header.txt \
        --output targets.tsv \
        --report ${q.call(mappingFile)}

    # One complete record scan performs overlap/END checks and
    # determines whether positive gVCF evidence exists.
    bcftools view -H ${vcfArg} |
        check_vcf_score_overlaps.py \
            --scorefile ${scoreArg} \
            --vcf-header vcf_header.txt \
            --format-report preflight_format.tsv \
            --expected-records expected_records.tsv

    detected_format=\$(cut -f1 preflight_format.tsv)
    end_declared=\$(cut -f2 preflight_format.tsv)

    case "${gvcfMode}" in
        auto)
            if [[ "\${detected_format}" != "gvcf" ]]; then
                echo "ERROR: No record-level gVCF evidence was found. Ordinary variant records cannot establish that this is a plain VCF rather than a gVCF subset. Specify --gvcf_mode plain or --gvcf_mode gvcf." >&2
                exit 1
            fi
            effective_format="gvcf"
            ;;
        gvcf|plain)
            if [[ "\${detected_format}" == "gvcf" &&
                  "${gvcfMode}" == "plain" ]]; then
                echo "ERROR: --gvcf_mode plain contradicts gVCF markers found in VCF records." >&2
                exit 1
            fi
            effective_format="${gvcfMode}"
            if [[ "\${detected_format}" == "ambiguous" ]]; then
                echo "WARNING: No record-level gVCF evidence was found; trusting explicit --gvcf_mode ${gvcfMode}." >&2
            fi
            ;;
    esac

    is_gvcf="false"
    if [[ "\${effective_format}" == "gvcf" ]]; then
        is_gvcf="true"
    fi

    if [[ "\${is_gvcf}" == "true" &&
          "${missingMode}" == "reference" ]]; then
        echo "ERROR: --missing_genotype reference is allowed only for declared plain VCF input." >&2
        exit 1
    fi

    printf 'INFO: %s: evidence=%s; effective=%s\\n' \
        ${vcfArg} "\${detected_format}" "\${effective_format}" >&2

    if [[ "\${end_declared}" == "true" ]]; then
        end_field='%INFO/END'
    else
        end_field='.'
    fi

    printf 'chrom\\tposition\\tref\\talt\\tend\\tgenotype\\tgp\\n' \
        > queried_genotypes.tsv

    # determine GP availability without allowing every undefined tag
    
    gp_field='.'
    if grep -q '^##FORMAT=<ID=GP[,>]' vcf_header.txt; then
        gp_field='%GP'
    fi

    gt_field='%GT'

    # Optional support for GP-only soft input:
    # represent an undeclared GT as unavailable, not as an explicit ".".
    # The R reader converts this NA marker to NA.
    if [[ "${genotypeCalls}" == "soft" ]] &&
       ! grep -q '^##FORMAT=<ID=GT[,>]' vcf_header.txt; then
        gt_field='NA'
    fi

    if [[ -s targets.tsv ]]; then
        bcftools view \
            --targets-file targets.tsv \
            --targets-overlap 1 \
            -Ou \
            ${vcfArg} \
        | bcftools query \
            --allow-undef-tag \
            --format "%CHROM\\t%POS\\t%REF\\t%ALT\\t\${end_field}[\\t%GT\\t%GP]\\n" \
            >> queried_genotypes.tsv
    fi

    verify_vcf_query.py \
        --expected expected_records.tsv \
        --queried queried_genotypes.tsv

    calculate_prs.R \
        --scorefile ${scoreArg} \
        --genotypes queried_genotypes.tsv \
        --sample_id ${sampleIdArg} \
        --vcf_sample_file vcf_samples.txt \
        --source_vcf ${sourceVcfArg} \
        --input_is_gvcf "\${is_gvcf}" \
        --gvcf_mode "${gvcfMode}" \
        --genotype_calls "${genotypeCalls}" \
        --non_additive "${nonAdditive}" \
        --missing_genotype "${missingMode}" \
        --no_calls "${noCallsMode}" ${minCoveredArgument} \
        --target_build ${buildArg} \
        --default_ploidy "${defaultPloidy}" \
        --strict_alleles "${strictAlleles}" \
        --output_summary "${meta.id}.prs.tsv" \
        --output_details "${meta.id}.variants.tsv"
    """
}
