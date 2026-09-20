#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PRSCALC } from './workflows/prs_calc'


def normaliseBuild(def rawBuild) {
    if (rawBuild == null) {
        throw new IllegalArgumentException(
            'A target genome build must be provided with --target_build'
        )
    }

    def build = rawBuild
        .toString()
        .trim()
        .toLowerCase()
        .replaceFirst(/^grch/, '')
        .replaceFirst(/^hg/, '')

    if (build == '38') {
        return 'hg38'
    }

    if (build == '37') {
        return 'hg37'
    }

    if (build == '19') {
        return 'hg19'
    }

    throw new IllegalArgumentException(
        "Unsupported target build '${rawBuild}'. " +
        'Supported values are hg38, hg37, hg19, 38, 37 and 19.'
    )
}


def stripOuterQuotes(def rawValue) {
    def value = rawValue.toString().trim()

    if (value.size() >= 2) {
        if (
            (value.startsWith('"') && value.endsWith('"')) ||
            (value.startsWith("'") && value.endsWith("'"))
        ) {
            return value.substring(1, value.size() - 1)
        }
    }

    return value
}


workflow {

    /*
     * Validate global parameters.
     */
    if (!params.input) {
        error 'Missing required parameter: --input'
    }

    if (!params.scorefile) {
        error 'Missing required parameter: --scorefile'
    }

    def requestedBuild = params.hg != null ? params.hg : params.target_build
    def targetBuild = normaliseBuild(requestedBuild)

    def missingMode = params.missing_genotype
        .toString()
        .trim()
        .toLowerCase()

    if (!(missingMode in ['reference', 'zero', 'error'])) {
        error(
            "--missing_genotype must be one of: reference, zero, error. " +
            "Received '${params.missing_genotype}'."
        )
    }

    def strictAlleles = params.strict_alleles
        .toString()
        .trim()
        .toLowerCase()

    if (!(strictAlleles in ['true', 'false'])) {
        error '--strict_alleles must be true or false'
    }

    def defaultPloidy
    def plotTopVariants

    try {
        defaultPloidy = params.default_ploidy as Integer
        plotTopVariants = params.plot_top_variants as Integer
    } catch (Exception _ignored) {
        error '--default_ploidy and --plot_top_variants must be integers'
    }

    if (defaultPloidy < 1) {
        error '--default_ploidy must be at least 1'
    }

    if (plotTopVariants < 1) {
        error '--plot_top_variants must be at least 1'
    }

    def inputSheet = file(params.input, checkIfExists: true)
    def scoreFile = file(params.scorefile, checkIfExists: true)
    def phenotypeParameter = params.phenotypes == null
        ? ''
        : params.phenotypes.toString().trim()

    def phenotypeFile = phenotypeParameter
        ? file(phenotypeParameter, checkIfExists: true)
        : null
    def inputSheetParent = inputSheet.parent ?: file(launchDir)

    log.info "Input file:      ${inputSheet}"
    log.info "Score file:      ${scoreFile}"
    log.info "Phenotypes:      ${phenotypeFile ?: 'not provided; ROC/OR plots disabled'}"
    log.info "Target build:    ${targetBuild}"
    log.info "Output directory:${params.outdir}"

    /*
     * Supported samplesheet layouts:
     *
     *   /path/to/sample.vcf.gz
     *
     * or:
     *
     *   sample_id<TAB>/path/to/sample.vcf.gz
     *
     * Blank lines, comments beginning with #, and a simple header are ignored.
     */
    ch_parsed_samples = channel
        .fromPath(inputSheet.toString(), checkIfExists: true)
        .splitText()
        .map { rawLine ->
            def line = rawLine.trim()

            if (!line || line.startsWith('#')) {
                return null
            }

            def fields = line
                .split('\t', -1)
                .collect { value -> stripOuterQuotes(value) }

            def firstField = fields[0].toLowerCase()

            if (
                fields.size() == 1 &&
                firstField in ['vcf', 'path', 'vcf_path']
            ) {
                return null
            }

            if (
                fields.size() >= 2 &&
                firstField in ['sample', 'sample_id', 'id'] &&
                fields[1].toLowerCase() in ['vcf', 'path', 'vcf_path']
            ) {
                return null
            }

            if (!(fields.size() in [1, 2])) {
                error(
                    "Invalid samplesheet line:\n${line}\n" +
                    'Expected either one VCF path or sample_id<TAB>VCF_path.'
                )
            }

            def explicitId = fields.size() == 2 ? fields[0] : null
            def vcfText = fields.size() == 2 ? fields[1] : fields[0]

            if (!vcfText) {
                error "Empty VCF path in samplesheet line: ${line}"
            }

            if (vcfText.startsWith('~/')) {
                vcfText = "${System.getProperty('user.home')}/${vcfText.substring(2)}"
            }

            def candidatePath = java.nio.file.Paths.get(vcfText)
            def resolvedPath = candidatePath.isAbsolute()
                ? candidatePath.normalize()
                : inputSheetParent.resolve(candidatePath).normalize()

            def vcf = file(resolvedPath.toString(), checkIfExists: true)
            def vcfName = vcf.name

            if (!(vcfName ==~ /(?i).+(\.vcf|gvcf|\.bcf)(\.gz|\.bgz)?$/)) {
                error(
                    "Unsupported input file '${vcf}'. " +
                    'Expected .vcf, .vcf.gz, .vcf.bgz, gvcf.gz, .bcf or .bcf.gz.'
                )
            }

            def inferredId = vcfName.replaceFirst(
                /(?i)\.(vcf|bcf)(\.gz|\.bgz)?$/,
                ''
            )

            def rawId = explicitId ?: inferredId

            if (!rawId) {
                error "Could not determine a sample ID for '${vcf}'"
            }

            def safeId = rawId
                .replaceAll(/[^A-Za-z0-9_.-]+/, '_')
                .replaceAll(/^[._-]+/, '')
                .replaceAll(/[._-]+$/, '')

            if (!safeId) {
                error "Sample ID '${rawId}' does not contain usable characters"
            }

            if (safeId != rawId) {
                log.warn "Sample ID '${rawId}' was normalised to '${safeId}'"
            }

            def meta = [
                id          : safeId,
                source_name : vcfName
            ]

            tuple(meta, vcf)
        }
        .filter { sample -> sample != null }

    /*
     * Collect once before scattering so that duplicate sample IDs can be
     * detected before tasks start publishing files.
     */
    ch_samples = ch_parsed_samples
        .toList()
        .flatMap { sampleRows ->
            if (!sampleRows) {
                error "No VCF/BCF files were found in '${inputSheet}'"
            }

            def duplicates = sampleRows
                .groupBy { row -> row[0].id }
                .findAll { _id, rows -> rows.size() > 1 }
                .keySet()

            if (duplicates) {
                error(
                    'Duplicate sample IDs were found: ' +
                    duplicates.sort().join(', ') +
                    '. Use the optional sample_id<TAB>VCF format to provide ' +
                    'unique IDs.'
                )
            }

            log.info "Found ${sampleRows.size()} sample VCF/BCF file(s)"
            sampleRows
        }

    ch_scorefile = channel.value(scoreFile)
    ch_target_build = channel.value(targetBuild)

    ch_phenotypes = phenotypeFile == null
        ? channel.empty()
        : channel.value(phenotypeFile)

    PRSCALC(
        ch_samples,
        ch_scorefile,
        ch_target_build,
        ch_phenotypes
    )

    workflow.onError = {
        log.error "Pipeline failed: ${workflow.errorMessage ?: 'No error message available'}"
    }

    workflow.onComplete = {
        log.info """
                 Pipeline completed
                 ------------------
                 Status:   ${workflow.success ? 'SUCCESS' : 'FAILED'}
                 Duration: ${workflow.duration}
                 Output:   ${params.outdir}
                 """.stripIndent()
    }
}
