// ============================================================
// modules/bowtie2.nf
// Align paired-end reads with Bowtie2
// ============================================================
// Recommendation: --local --very-sensitive-local is preferred for
// CUT&Tag because reads may be slightly shorter after trimming and
// local alignment tolerates soft-clipping at read ends. The fragment
// size range -I 10 -X 700 captures nucleosome-free and mono-
// nucleosomal fragments typical of CUT&Tag.
//
// NOTE: bowtie2_index is declared as val (not path) so that Nextflow
// does not attempt to stage the index prefix as a file.  Bowtie2
// index files share a common prefix (e.g. /data/hg38) and the actual
// files are hg38.1.bt2, hg38.2.bt2, etc.  Passing the prefix as a
// path causes Nextflow to look for a file/directory literally named
// "hg38", which does not exist, resulting in a broken symlink and
// bowtie2 exit code 2.  Passing it as val lets bowtie2 receive the
// exact absolute path string the user supplied via --bowtie2_index.
//
// NOTE: --rg-id and --rg tags are added so that every read in the
// output BAM carries an @RG header and RG:Z aux tag.  Picard 3.x
// MarkDuplicates requires read group information to be present;
// without it, getReadGroup() returns null and Picard throws a
// NullPointerException at buildSortedReadEndLists.

process BOWTIE2_ALIGN {
    tag "${sample_id}"
    label 'process_medium'

    publishDir "${params.outdir}/04_alignment/logs",
        mode: 'copy',
        pattern: "*.log"

    input:
    tuple val(sample_id), path(r1), path(r2)
    val  bowtie2_index   // absolute path prefix, e.g. /data/indexes/hg38

    output:
    tuple val(sample_id), path("${sample_id}.bam"), emit: bam
    path "${sample_id}_bowtie2.log",                emit: log
    path "versions.yml",                            emit: versions

    script:
    def bt2_args = params.bowtie2_args ?: "--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700"
    """
    bowtie2 \\
        -p ${task.cpus} \\
        ${bt2_args} \\
        --rg-id "${sample_id}" \\
        --rg "SM:${sample_id}" \\
        --rg "PL:ILLUMINA" \\
        --rg "LB:${sample_id}" \\
        -x "${bowtie2_index}" \\
        -1 "${r1}" \\
        -2 "${r2}" \\
        2> "${sample_id}_bowtie2.log" \\
    | samtools view -bS -F 4 - \\
    > "${sample_id}.bam"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie2: \$(bowtie2 --version | head -1 | sed 's/.*bowtie2-align-s version //')
        samtools: \$(samtools --version | head -1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch ${sample_id}.bam ${sample_id}_bowtie2.log versions.yml
    """
}
