// ============================================================
// modules/validate_association.nf
// Validate association CSV against detected FASTQ pairs
// ============================================================

process VALIDATE_ASSOCIATION {
    tag "validate_association"
    label 'process_single'

    publishDir "${params.outdir}/00_fastq_pairs", mode: 'copy'

    input:
    path association_csv
    path pairs_tsv

    output:
    path "association_validated.csv", emit: validated_csv
    path "validation_report.txt",     emit: report

    script:
    def allow_no_ctrl = params.allow_no_control ? "--allow_no_control" : ""
    """
    validate_association.py \\
        --association_csv "${association_csv}" \\
        --pairs_tsv       "${pairs_tsv}" \\
        --output          association_validated.csv \\
        --report          validation_report.txt \\
        ${allow_no_ctrl}
    """

    stub:
    """
    cp ${association_csv} association_validated.csv
    touch validation_report.txt
    """
}
