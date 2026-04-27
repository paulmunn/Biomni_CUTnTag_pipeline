// ============================================================
// modules/version_log.nf
// Log all tool versions for reproducibility
// ============================================================

process VERSION_LOG {
    tag "version_log"
    label 'process_single'

    publishDir "${params.outdir}/pipeline_info",
        mode: 'copy'

    output:
    path "software_versions.yml", emit: versions
    path "pipeline_info.txt",     emit: info

    script:
    """
    # Collect versions from all available tools
    {
    echo "pipeline:"
    echo "  name: cuttag-pipeline"
    echo "  version: 1.0.0"
    echo "  nextflow: \$(nextflow -version 2>&1 | grep version | head -1 | sed 's/.*version //')"
    echo "tools:"
    echo "  fastqc: \$(fastqc --version 2>&1 | sed 's/FastQC v//' || echo 'not found')"
    echo "  multiqc: \$(multiqc --version 2>&1 | sed 's/multiqc, version //' || echo 'not found')"
    echo "  cutadapt: \$(cutadapt --version 2>&1 || echo 'not found')"
    echo "  bowtie2: \$(bowtie2 --version 2>&1 | head -1 | sed 's/.*bowtie2-align-s version //' || echo 'not found')"
    echo "  samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //' || echo 'not found')"
    echo "  picard: \$(picard MarkDuplicates --version 2>&1 | grep -i version | head -1 || echo 'not found')"
    echo "  bedtools: \$(bedtools --version 2>&1 | sed 's/bedtools v//' || echo 'not found')"
    echo "  deeptools: \$(bamCoverage --version 2>&1 | sed 's/bamCoverage //' || echo 'not found')"
    echo "  macs3: \$(macs3 --version 2>&1 | sed 's/macs3 //' || echo 'not found')"
    echo "  python: \$(python --version 2>&1 | sed 's/Python //' || echo 'not found')"
    echo "  r: \$(R --version 2>&1 | head -1 | sed 's/R version //' | sed 's/ .*//' || echo 'not found')"
    } > software_versions.yml

    # Pipeline run info
    {
    echo "Pipeline: CUT&Tag Nextflow DSL2 Pipeline v1.0.0"
    echo "Run date: \$(date)"
    echo "Hostname: \$(hostname)"
    echo "Work dir: ${workflow.workDir}"
    echo "Launch dir: ${workflow.launchDir}"
    echo "Profile: ${workflow.profile}"
    echo "Command: ${workflow.commandLine}"
    } > pipeline_info.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -1 | sed 's/.*version //' | sed 's/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch software_versions.yml pipeline_info.txt versions.yml
    """
}
