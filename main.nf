#!/usr/bin/env nextflow
// ============================================================
// CUT&Tag Nextflow DSL2 Pipeline
// Author: Genomics Innovation Hub
// Version: 1.0.0
// ============================================================

nextflow.enable.dsl = 2

// ── Parameter defaults ──────────────────────────────────────
params.input_dir            = null
params.outdir               = "results"
params.genome               = null
params.bowtie2_index        = null
params.chrom_sizes          = null
params.blacklist            = null
params.annotation_gtf       = null
params.tss_bed              = null
params.association_csv      = null

// Adapter sequences (Illumina TruSeq / Nextera defaults)
params.adapter_fwd          = "CTGTCTCTTATACACATCT"   // Nextera Transposase Read 1
params.adapter_rev          = "CTGTCTCTTATACACATCT"   // Nextera Transposase Read 2
params.min_length           = 20

// Alignment
params.bowtie2_args         = "--local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700"
params.mapq                 = 20

// Mitochondrial chromosome names (pipe-separated for grep)
params.mito_name            = "chrM|MT|M|mitochondria|Mito"

// Deduplication
params.dedup_mode           = "picard"   // picard | umi_tools
params.umi_separator        = ":"        // UMI separator in read name

// Signal tracks
params.bigwig_binsize       = 10
params.normalization        = "CPM"      // CPM | RPGC | BPM | RPKM | None
params.effective_genome_size = 2913022398  // hg38 default; override for other genomes

// Peak calling
params.macs3_qvalue         = 0.05
params.macs3_genome         = "hs"      // hs | mm | ce | dm
params.allow_no_control     = false
params.allow_extra_samples  = false   // Allow association CSV rows with no matching FASTQ pair

// TSS / deepTools
params.tss_window           = 2000
params.peak_center_window   = 1000
params.threads              = 4

// FASTQ pairing
params.paired_pattern       = "auto"    // auto | _R1/_R2 | _1/_2 | .R1/.R2

// Annotation
params.annotation_tool      = "chipseeker"  // chipseeker | homer

// Misc
params.skip_umi             = true
params.skip_blacklist       = false
params.help                 = false

// ── Help message ────────────────────────────────────────────
def helpMessage() {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║          CUT&Tag Nextflow DSL2 Pipeline v1.0.0          ║
    ╚══════════════════════════════════════════════════════════╝

    Usage:
        nextflow run main.nf \\
            --input_dir     /path/to/fastqs \\
            --outdir        results \\
            --genome        hg38 \\
            --bowtie2_index /path/to/bowtie2/hg38 \\
            --chrom_sizes   /path/to/hg38.chrom.sizes \\
            --association_csv samples.csv \\
            --annotation_gtf /path/to/hg38.gtf \\
            --tss_bed       /path/to/hg38_tss.bed \\
            -profile        docker

    Required parameters:
        --input_dir         Directory containing paired-end FASTQ files
        --outdir            Output directory
        --genome            Genome identifier (hg38, mm10, rn6, custom)
        --bowtie2_index     Path to Bowtie2 index basename
        --chrom_sizes       Chromosome sizes file
        --association_csv   Sample association/metadata CSV

    Optional parameters:
        --blacklist         BED file of blacklisted regions
        --annotation_gtf    GTF annotation file (required for peak annotation)
        --tss_bed           TSS BED file (required for TSS enrichment)
        --adapter_fwd       Forward adapter sequence [${params.adapter_fwd}]
        --adapter_rev       Reverse adapter sequence [${params.adapter_rev}]
        --min_length        Minimum read length after trimming [${params.min_length}]
        --bowtie2_args      Bowtie2 alignment arguments
        --mapq              Minimum MAPQ score [${params.mapq}]
        --mito_name         Mitochondrial chromosome name pattern [${params.mito_name}]
        --dedup_mode        Deduplication mode: picard|umi_tools [${params.dedup_mode}]
        --bigwig_binsize    BigWig bin size [${params.bigwig_binsize}]
        --normalization     BigWig normalization: CPM|RPGC|BPM|RPKM|None [${params.normalization}]
        --effective_genome_size  Effective genome size for RPGC [${params.effective_genome_size}]
        --macs3_qvalue      MACS3 q-value threshold [${params.macs3_qvalue}]
        --macs3_genome      MACS3 genome size: hs|mm|ce|dm [${params.macs3_genome}]
        --allow_no_control  Allow peak calling without control [${params.allow_no_control}]
        --tss_window        TSS window size for enrichment [${params.tss_window}]
        --peak_center_window Window around peak centers [${params.peak_center_window}]
        --threads           Default number of threads [${params.threads}]
        --annotation_tool   Peak annotation tool: chipseeker|homer [${params.annotation_tool}]
        --skip_umi          Skip UMI-aware deduplication [${params.skip_umi}]
        --paired_pattern    FASTQ pairing pattern: auto|_R1/_R2|_1/_2 [${params.paired_pattern}]
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// ── Validate required parameters ────────────────────────────
def validateParams() {
    def errors = []
    if (!params.input_dir)       errors << "ERROR: --input_dir is required"
    if (!params.genome)          errors << "ERROR: --genome is required"
    if (!params.bowtie2_index)   errors << "ERROR: --bowtie2_index is required"
    if (!params.chrom_sizes)     errors << "ERROR: --chrom_sizes is required"
    if (!params.association_csv) errors << "ERROR: --association_csv is required"

    if (params.input_dir && !file(params.input_dir).exists())
        errors << "ERROR: --input_dir does not exist: ${params.input_dir}"
    if (params.association_csv && !file(params.association_csv).exists())
        errors << "ERROR: --association_csv does not exist: ${params.association_csv}"
    if (params.chrom_sizes && !file(params.chrom_sizes).exists())
        errors << "ERROR: --chrom_sizes does not exist: ${params.chrom_sizes}"
    if (params.blacklist && !file(params.blacklist).exists())
        errors << "ERROR: --blacklist does not exist: ${params.blacklist}"
    if (params.annotation_gtf && !file(params.annotation_gtf).exists())
        errors << "ERROR: --annotation_gtf does not exist: ${params.annotation_gtf}"
    if (params.tss_bed && !file(params.tss_bed).exists())
        errors << "ERROR: --tss_bed does not exist: ${params.tss_bed}"

    if (errors) {
        errors.each { log.error it }
        exit 1, "Parameter validation failed. See errors above."
    }
}

validateParams()

// ── Import modules ───────────────────────────────────────────
include { PAIR_FASTQS }              from './modules/pair_fastqs'
include { VALIDATE_ASSOCIATION }     from './modules/validate_association'
include { FASTQC as FASTQC_RAW }     from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { MULTIQC as MULTIQC_RAW }   from './modules/multiqc'
include { MULTIQC as MULTIQC_TRIM }  from './modules/multiqc'
include { MULTIQC as MULTIQC_ALIGN } from './modules/multiqc'
include { MULTIQC as MULTIQC_FINAL } from './modules/multiqc'
include { CUTADAPT }                 from './modules/cutadapt'
include { BOWTIE2_ALIGN }            from './modules/bowtie2'
include { SAMTOOLS_SORT_INDEX }      from './modules/samtools'
include { SAMTOOLS_FLAGSTAT }        from './modules/samtools'
include { SAMTOOLS_IDXSTATS }        from './modules/samtools'
include { SAMTOOLS_STATS }           from './modules/samtools'
include { FILTER_MITO }              from './modules/filter'
include { FILTER_BLACKLIST }         from './modules/filter'
include { PICARD_MARKDUPLICATES }    from './modules/filter'
include { UMI_TOOLS_DEDUP }          from './modules/filter'
include { BAMCOVERAGE_BIGWIG }       from './modules/deeptools'
include { BAMCOVERAGE_BEDGRAPH }     from './modules/deeptools'
include { COMPUTE_MATRIX_TSS }                    from './modules/deeptools'
include { COMPUTE_MATRIX_PEAKS }                  from './modules/deeptools'
include { PLOT_PROFILE as PLOT_PROFILE_TSS }      from './modules/deeptools'
include { PLOT_PROFILE as PLOT_PROFILE_PEAKS }    from './modules/deeptools'
include { PLOT_HEATMAP as PLOT_HEATMAP_TSS }      from './modules/deeptools'
include { PLOT_HEATMAP as PLOT_HEATMAP_PEAKS }    from './modules/deeptools'
include { MACS3_CALLPEAK_SAMPLE }    from './modules/macs3'
include { MACS3_CALLPEAK_GROUP }     from './modules/macs3'
include { MERGE_BAMS as MERGE_TREATMENT_BAMS } from './modules/merge_bams'
include { MERGE_BAMS as MERGE_CONTROL_BAMS }   from './modules/merge_bams'
include { ANNOTATE_PEAKS }           from './modules/annotate_peaks'
include { FRIP_SCORE }               from './modules/frip'
include { READ_RETENTION_SUMMARY }   from './modules/read_retention'
include { MAKE_MQC_CUSTOM }          from './modules/make_mqc_custom'
include { VERSION_LOG }              from './modules/version_log'

// ── Main workflow ────────────────────────────────────────────
workflow {

    // ── Step 0: Log pipeline info ──────────────────────────
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║          CUT&Tag Nextflow DSL2 Pipeline v1.0.0          ║
    ╚══════════════════════════════════════════════════════════╝
    input_dir       : ${params.input_dir}
    outdir          : ${params.outdir}
    genome          : ${params.genome}
    bowtie2_index   : ${params.bowtie2_index}
    association_csv : ${params.association_csv}
    dedup_mode      : ${params.dedup_mode}
    normalization   : ${params.normalization}
    annotation_tool : ${params.annotation_tool}
    """.stripIndent()

    VERSION_LOG()

    // ── Step 1: Detect and validate FASTQ pairs ────────────
    // Pass input_dir as a plain string (val), not file(), so that
    // pair_fastqs.py receives the real filesystem path and can
    // recursively search all subdirectories with Path.rglob().
    // Absolute path resolution is done here to guard against
    // relative paths being misinterpreted inside the work directory.
    PAIR_FASTQS(
        file(params.input_dir).toAbsolutePath().toString(),
        params.paired_pattern
    )

    // ── Step 2: Validate association CSV ──────────────────
    VALIDATE_ASSOCIATION(
        file(params.association_csv),
        PAIR_FASTQS.out.pairs_tsv
    )

    // ── Build sample channel from validated pairs TSV ─────
    // Channel: [sample_id, r1, r2]
    ch_samples = PAIR_FASTQS.out.pairs_tsv
        .splitCsv(header: true, sep: '\t')
        .filter { row -> row.status == 'OK' }
        .map { row ->
            def meta = [
                id:       row.sample_id,
                r1:       file(row.r1_path),
                r2:       file(row.r2_path)
            ]
            return meta
        }

    // ── Load association metadata ─────────────────────────
    // Filter to only the samples whose genome column matches --genome.
    // This allows a single association CSV to cover a multi-species
    // experiment; each pipeline invocation processes one genome at a time.
    ch_assoc = Channel.fromPath(params.association_csv)
        .splitCsv(header: true)
        .filter { row ->
            def keep = row.genome?.trim() == params.genome?.trim()
            if (!keep) {
                log.info "Skipping sample '${row.sample_id}' " +
                         "(genome='${row.genome}' != --genome '${params.genome}')"
            }
            keep
        }
        .map { row ->
            [
                row.sample_id,
                [
                    species:          row.species,
                    genome:           row.genome,
                    antibody:         row.antibody,
                    condition:        row.condition,
                    replicate:        row.replicate,
                    group_id:         row.group_id,
                    is_control:       row.is_control.toLowerCase() == 'true',
                    control_group_id: row.control_group_id,
                    merge_group_id:   row.merge_group_id,
                    peak_mode:        row.peak_calling_mode ?: 'narrow',
                    notes:            row.notes ?: ''
                ]
            ]
        }

    // Join sample channel with metadata
    ch_samples_meta = ch_samples
        .map { meta -> [meta.id, meta] }
        .join(ch_assoc, by: 0)
        .map { sample_id, meta, assoc ->
            [sample_id, meta.r1, meta.r2, assoc]
        }

    // ── Step 3: Raw FastQC ─────────────────────────────────
    ch_raw_reads = ch_samples_meta.map { id, r1, r2, assoc -> [id, r1, r2] }

    FASTQC_RAW(ch_raw_reads, "raw")

    // Extract bare zip paths from the FASTQC tuple output.
    // FASTQC now emits tuple(val(id), path(r1_zip), path(r2_zip)) —
    // two explicit named paths, NOT a glob list — so we extract both
    // path values and flatten into a single stream of bare path values.
    ch_fastqc_raw_paths = FASTQC_RAW.out.zip
        .flatMap { id, r1_zip, r2_zip -> [r1_zip, r2_zip] }

    MULTIQC_RAW(
        ch_fastqc_raw_paths.collect(),
        "raw",
        "${params.outdir}/13_multiqc/01_raw"
    )

    // ── Step 4: Adapter trimming ───────────────────────────
    CUTADAPT(ch_raw_reads)

    FASTQC_TRIMMED(CUTADAPT.out.trimmed_reads, "trimmed")

    // CUTADAPT.out.log emits a bare path; FASTQC_TRIMMED.out.zip emits tuples.
    // FASTQC now emits tuple(val(id), path(r1_zip), path(r2_zip)) — two
    // explicit named paths — so flatMap extracts both into a bare path stream.
    ch_fastqc_trim_paths = FASTQC_TRIMMED.out.zip
        .flatMap { id, r1_zip, r2_zip -> [r1_zip, r2_zip] }

    MULTIQC_TRIM(
        ch_fastqc_trim_paths.mix(CUTADAPT.out.log).collect(),
        "trimmed",
        "${params.outdir}/13_multiqc/02_trimmed"
    )

    // ── Step 5: Alignment ──────────────────────────────────
    ch_trimmed_meta = CUTADAPT.out.trimmed_reads
        .map { id, r1, r2 -> [id, r1, r2] }
        .join(ch_samples_meta.map { id, r1, r2, assoc -> [id, assoc] }, by: 0)

    // Pass bowtie2_index as a plain string (val), NOT file(), so that
    // Nextflow does not try to stage the prefix as a path.  The index
    // files (hg38.1.bt2, hg38.2.bt2, …) share a common prefix and
    // there is no file/directory literally named by that prefix.
    BOWTIE2_ALIGN(
        CUTADAPT.out.trimmed_reads,
        params.bowtie2_index
    )

    SAMTOOLS_SORT_INDEX(BOWTIE2_ALIGN.out.bam)
    SAMTOOLS_FLAGSTAT(SAMTOOLS_SORT_INDEX.out.bam_bai)
    SAMTOOLS_IDXSTATS(SAMTOOLS_SORT_INDEX.out.bam_bai)
    SAMTOOLS_STATS(SAMTOOLS_SORT_INDEX.out.bam_bai)

    // ── Step 6: Filtering ──────────────────────────────────
    // 6a. Remove mitochondrial reads
    FILTER_MITO(SAMTOOLS_SORT_INDEX.out.bam_bai)

    // 6b. Remove blacklisted regions (if provided)
    ch_after_mito = FILTER_MITO.out.bam_bai
    if (params.blacklist && !params.skip_blacklist) {
        FILTER_BLACKLIST(ch_after_mito, file(params.blacklist))
        ch_for_dedup = FILTER_BLACKLIST.out.bam_bai
    } else {
        ch_for_dedup = ch_after_mito
    }

    // 6c. Deduplication
    // ch_dedup_metrics is assigned INSIDE the if/else so that only the
    // invoked process's .out channel is ever referenced.  Using a ternary
    // operator outside the if/else would wire BOTH process outputs into the
    // dataflow graph, causing the un-invoked branch to emit a DataflowStream
    // containing only a PoisonPill — which then poisons any downstream .collect().
    if (params.dedup_mode == 'umi_tools' && !params.skip_umi) {
        UMI_TOOLS_DEDUP(ch_for_dedup)
        ch_final_bam    = UMI_TOOLS_DEDUP.out.bam_bai
        ch_dedup_metrics = UMI_TOOLS_DEDUP.out.metrics
    } else {
        PICARD_MARKDUPLICATES(ch_for_dedup)
        ch_final_bam    = PICARD_MARKDUPLICATES.out.bam_bai
        ch_dedup_metrics = PICARD_MARKDUPLICATES.out.metrics
    }

    // 6d. Read retention summary

    ch_read_counts = SAMTOOLS_FLAGSTAT.out.stats
        .join(FILTER_MITO.out.stats, by: 0)
        .join(ch_dedup_metrics, by: 0, remainder: true)

    READ_RETENTION_SUMMARY(ch_read_counts)

    // ── Step 7: Signal tracks ──────────────────────────────
    BAMCOVERAGE_BIGWIG(
        ch_final_bam,
        file(params.chrom_sizes)
    )

    BAMCOVERAGE_BEDGRAPH(
        ch_final_bam,
        file(params.chrom_sizes)
    )

    // ── Step 8: Per-sample peak calling ───────────────────
    // Join final BAMs with metadata to get control relationships
    ch_final_bam_meta = ch_final_bam
        .map { id, bam, bai -> [id, bam, bai] }
        .join(ch_samples_meta.map { id, r1, r2, assoc -> [id, assoc] }, by: 0)

    // Separate treatment and control samples
    ch_treatment_samples = ch_final_bam_meta
        .filter { id, bam, bai, assoc -> !assoc.is_control }

    // ── Step 9: Group-level BAM merging ───────────────────
    // Group treatment BAMs by merge_group_id
    // Emit: [merge_group_id, [bam_list], [bai_list], ctrl_group_id, peak_mode]
    ch_merge_groups = ch_final_bam_meta
        .filter { id, bam, bai, assoc -> !assoc.is_control }
        .map { id, bam, bai, assoc ->
            tuple(assoc.merge_group_id, bam, bai, assoc.control_group_id, assoc.peak_mode)
        }
        .groupTuple(by: 0)
        // After groupTuple: [merge_id, [bams], [bais], [ctrl_ids], [modes]]
        // ctrl_id and mode are the same for all members of a merge group (validated upstream)
        .map { merge_id, bams, bais, ctrl_ids, modes ->
            tuple(merge_id, bams, bais, ctrl_ids[0], modes[0])
        }

    // Group control BAMs by group_id
    // Emit: [ctrl_group_id, [bam_list], [bai_list]]
    ch_ctrl_groups = ch_final_bam_meta
        .filter { id, bam, bai, assoc -> assoc.is_control }
        .map { id, bam, bai, assoc -> tuple(assoc.group_id, bam, bai) }
        .groupTuple(by: 0)

    // Merge treatment BAMs per merge_group_id
    MERGE_TREATMENT_BAMS(
        ch_merge_groups.map { merge_id, bams, bais, ctrl_id, mode ->
            tuple(merge_id, bams, bais, 'treatment')
        }
    )

    // Merge control BAMs per control group_id
    MERGE_CONTROL_BAMS(
        ch_ctrl_groups.map { ctrl_id, bams, bais ->
            tuple(ctrl_id, bams, bais, 'control')
        }
    )

    // Build lookup: ctrl_group_id -> [merged_ctrl_bam, merged_ctrl_bai]
    ch_merged_ctrl_lookup = MERGE_CONTROL_BAMS.out.merged_bam
        .map { ctrl_id, bam, bai, type -> tuple(ctrl_id, bam, bai) }

    // ── Step 8: Per-sample peak calling ───────────────────
    // Use the MERGED control BAM (not individual control replicates) so
    // that each treatment sample is matched to exactly one control entry.
    // Using individual control BAMs with .combine() would produce one row
    // per control replicate, causing MACS3 to run multiple times per
    // treatment sample and generating duplicate *_peak_stats.txt filenames
    // that collide when staged into MAKE_MQC_CUSTOM.
    ch_sample_with_control = ch_treatment_samples
        .map { id, bam, bai, assoc ->
            [assoc.control_group_id, id, bam, bai, assoc]
        }
        .join(ch_merged_ctrl_lookup, by: 0)
        .map { ctrl_grp_id, id, bam, bai, assoc, ctrl_bam, ctrl_bai ->
            [id, bam, bai, ctrl_bam, ctrl_bai, assoc]
        }

    MACS3_CALLPEAK_SAMPLE(ch_sample_with_control)

    // Join merged treatment BAMs with their ctrl_group_id metadata,
    // then join with the merged control BAM on ctrl_group_id
    ch_group_peak_input = MERGE_TREATMENT_BAMS.out.merged_bam
        .map { merge_id, bam, bai, type -> tuple(merge_id, bam, bai) }
        // Re-attach ctrl_id and peak_mode from ch_merge_groups
        .join(
            ch_merge_groups.map { merge_id, bams, bais, ctrl_id, mode ->
                tuple(merge_id, ctrl_id, mode)
            },
            by: 0
        )
        // Now: [merge_id, merged_bam, merged_bai, ctrl_id, mode]
        // Join with merged control on ctrl_id
        .map { merge_id, bam, bai, ctrl_id, mode ->
            tuple(ctrl_id, merge_id, bam, bai, mode)
        }
        .join(ch_merged_ctrl_lookup, by: 0)
        // Now: [ctrl_id, merge_id, bam, bai, mode, ctrl_bam, ctrl_bai]
        .map { ctrl_id, merge_id, bam, bai, mode, ctrl_bam, ctrl_bai ->
            tuple(merge_id, bam, bai, ctrl_bam, ctrl_bai, mode)
        }

    MACS3_CALLPEAK_GROUP(ch_group_peak_input)

    // ── Step 10: Peak annotation ───────────────────────────
    if (params.annotation_gtf) {
        ch_peaks_for_annotation = MACS3_CALLPEAK_SAMPLE.out.peaks
            .mix(MACS3_CALLPEAK_GROUP.out.peaks)

        ANNOTATE_PEAKS(
            ch_peaks_for_annotation,
            file(params.annotation_gtf)
        )
    }

    // ── Step 11: FRiP scores ───────────────────────────────
    // Per-sample FRiP against own peaks
    ch_frip_input = ch_final_bam
        .map { id, bam, bai -> [id, bam, bai] }
        .join(MACS3_CALLPEAK_SAMPLE.out.peaks, by: 0)

    FRIP_SCORE(ch_frip_input, "per_sample")

    // ── Step 12: deepTools QC ─────────────────────────────
    if (params.tss_bed) {
        // Extract bare BigWig paths from tuple(val(id), path(bw)) channel
        ch_bigwigs_for_matrix = BAMCOVERAGE_BIGWIG.out.bigwig
            .map { id, bw -> bw }
            .collect()

        COMPUTE_MATRIX_TSS(
            ch_bigwigs_for_matrix,
            file(params.tss_bed),
            params.tss_window
        )

        PLOT_PROFILE_TSS(COMPUTE_MATRIX_TSS.out.matrix, "tss")
        PLOT_HEATMAP_TSS(COMPUTE_MATRIX_TSS.out.matrix, "tss")
    }

    // Peak-center matrix using group-level peaks
    // Extract bare paths from tuple channels before passing to computeMatrix
    ch_group_peaks_for_matrix = MACS3_CALLPEAK_GROUP.out.peaks
        .map { id, peak -> peak }
        .collect()

    COMPUTE_MATRIX_PEAKS(
        BAMCOVERAGE_BIGWIG.out.bigwig.map { id, bw -> bw }.collect(),
        ch_group_peaks_for_matrix,
        params.peak_center_window
    )

    PLOT_PROFILE_PEAKS(COMPUTE_MATRIX_PEAKS.out.matrix, "peaks")
    PLOT_HEATMAP_PEAKS(COMPUTE_MATRIX_PEAKS.out.matrix, "peaks")

    // ── Step 13: MultiQC custom content ───────────────────
    // MAKE_MQC_CUSTOM receives flat lists of bare path files.
    // All four source channels emit tuple(val(id), path(file)) —
    // extract the path element from each before collecting.
    // Deduplicate peak stats by sample_id before collecting.
    // ch_sample_with_control is built with .combine(), which can produce
    // multiple rows for the same treatment sample if the control-matching
    // logic yields more than one hit.  That causes MACS3_CALLPEAK_SAMPLE
    // to run (or cache) multiple times for the same sample_id, emitting
    // duplicate *_peak_stats.txt paths with identical basenames.
    // Nextflow then raises a "file name collision" error when staging them
    // all into MAKE_MQC_CUSTOM's flat input directory.
    // .unique { id, f -> id } keeps only the first emission per sample_id.
    ch_peak_stats_dedup = MACS3_CALLPEAK_SAMPLE.out.stats
        .unique { id, f -> id }
        .map    { id, f -> f }
        .collect()

    MAKE_MQC_CUSTOM(
        READ_RETENTION_SUMMARY.out.tsv.map    { id, f -> f }.collect(),
        FRIP_SCORE.out.tsv.map                { id, f -> f }.collect(),
        ch_peak_stats_dedup,
        SAMTOOLS_FLAGSTAT.out.stats.map       { id, f -> f }.collect()
    )

    // ── Step 14: Final MultiQC reports ────────────────────
    // ── MultiQC alignment report ───────────────────────────
    // All three samtools channels emit tuple(val(id), path(file)).
    // Extract the path element from each before mixing and collecting.
    ch_align_qc = SAMTOOLS_FLAGSTAT.out.stats.map { id, f -> f }
        .mix(SAMTOOLS_IDXSTATS.out.stats.map { id, f -> f })
        .mix(SAMTOOLS_STATS.out.stats.map    { id, f -> f })
        .collect()

    MULTIQC_ALIGN(
        ch_align_qc,
        "alignment",
        "${params.outdir}/13_multiqc/03_alignment"
    )

    // ── Final integrated MultiQC report ───────────────────
    // Build a flat channel of bare path values from every QC source.
    // Rule: every tuple channel needs .map { id, f -> f } (or
    // .map { id, files -> files }.flatten() for multi-file tuples)
    // before being mixed. Bare-path channels (CUTADAPT.out.log,
    // MAKE_MQC_CUSTOM.out.all_mqc) can be mixed in directly.

    // Dedup metrics for MultiQC: use ch_dedup_metrics which was already
    // set to whichever dedup process was actually invoked (Picard or UMI-tools).
    // NEVER reference PICARD_MARKDUPLICATES.out or UMI_TOOLS_DEDUP.out directly
    // here — referencing .out on an un-invoked process returns a DataflowStream
    // containing only a PoisonPill, which poisons the entire .collect() call.
    ch_dedup_mqc_paths = ch_dedup_metrics.map { id, f -> f }

    // Use individually-named MAKE_MQC_CUSTOM outputs instead of the glob
    // all_mqc emit, which produced a DataflowStream that leaked PoisonPill
    // into the collected list.  Each named output is a single path value.
    ch_final_mqc = ch_fastqc_raw_paths
        .mix( ch_fastqc_trim_paths )
        .mix( CUTADAPT.out.log )
        .mix( SAMTOOLS_FLAGSTAT.out.stats.map { id, f -> f } )
        .mix( SAMTOOLS_STATS.out.stats.map    { id, f -> f } )
        .mix( ch_dedup_mqc_paths )
        .mix( MAKE_MQC_CUSTOM.out.retention_mqc )
        .mix( MAKE_MQC_CUSTOM.out.frip_mqc )
        .mix( MAKE_MQC_CUSTOM.out.peaks_mqc )
        .mix( MAKE_MQC_CUSTOM.out.mito_mqc )
        .collect()

    MULTIQC_FINAL(
        ch_final_mqc,
        "final",
        "${params.outdir}/13_multiqc/04_final"
    )
}

