#!/usr/bin/env Rscript
# ============================================================
# annotate_peaks_chipseeker.R
# Annotate peaks using ChIPseeker + TxDb built from GTF
# ============================================================
# Usage:
#   Rscript annotate_peaks_chipseeker.R \
#     --peaks   sample_peaks.narrowPeak \
#     --gtf     annotation.gtf \
#     --prefix  sample_id \
#     --outdir  .
#
# Outputs:
#   {prefix}_annotated.tsv          - Full annotation table
#   {prefix}_annotation_summary.tsv - Category counts
#   {prefix}_annotation_pie.png     - Pie chart
#   {prefix}_annotation_bar.png     - Bar chart
# ============================================================

suppressPackageStartupMessages({
    library(optparse)
    library(ChIPseeker)
    library(GenomicFeatures)
    library(ggplot2)
    library(dplyr)
})

# ── Argument parsing ──────────────────────────────────────────
option_list <- list(
    make_option("--peaks",  type="character", help="Peak file (narrowPeak/broadPeak/BED)"),
    make_option("--gtf",    type="character", help="GTF annotation file"),
    make_option("--prefix", type="character", default="sample", help="Output prefix"),
    make_option("--outdir", type="character", default=".", help="Output directory"),
    make_option("--tss_upstream",   type="integer", default=3000,
                help="TSS upstream window [default: 3000]"),
    make_option("--tss_downstream", type="integer", default=3000,
                help="TSS downstream window [default: 3000]")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$peaks) || is.null(opt$gtf)) {
    stop("--peaks and --gtf are required.")
}

dir.create(opt$outdir, showWarnings=FALSE, recursive=TRUE)

# ── Build TxDb from GTF ───────────────────────────────────────
message("Building TxDb from GTF: ", opt$gtf)
txdb <- tryCatch(
    makeTxDbFromGFF(opt$gtf, format="GTF"),
    error = function(e) {
        message("ERROR building TxDb: ", e$message)
        quit(status=1)
    }
)

# ── Load peaks ────────────────────────────────────────────────
message("Loading peaks: ", opt$peaks)
peaks <- tryCatch(
    readPeakFile(opt$peaks, as="GRanges"),
    error = function(e) {
        # Fallback: read as BED
        message("Trying BED fallback for peak file...")
        df <- read.table(opt$peaks, header=FALSE, sep="\t",
                         stringsAsFactors=FALSE)
        GRanges(seqnames=df[,1],
                ranges=IRanges(start=df[,2]+1, end=df[,3]),
                score=if(ncol(df)>=5) df[,5] else rep(0, nrow(df)))
    }
)

message("Loaded ", length(peaks), " peaks.")

if (length(peaks) == 0) {
    message("WARNING: No peaks to annotate. Writing empty output files.")
    write.table(data.frame(), file.path(opt$outdir, paste0(opt$prefix, "_annotated.tsv")),
                sep="\t", quote=FALSE, row.names=FALSE)
    write.table(data.frame(annotation=character(), count=integer()),
                file.path(opt$outdir, paste0(opt$prefix, "_annotation_summary.tsv")),
                sep="\t", quote=FALSE, row.names=FALSE)
    quit(status=0)
}

# ── Annotate peaks ────────────────────────────────────────────
message("Annotating peaks...")
anno <- annotatePeak(
    peaks,
    tssRegion    = c(-opt$tss_upstream, opt$tss_downstream),
    TxDb         = txdb,
    annoDb       = NULL,   # No org.db required; avoids dependency issues
    verbose      = FALSE
)

# ── Write annotated TSV ───────────────────────────────────────
anno_df <- as.data.frame(anno)

# Clean up column names
colnames(anno_df) <- gsub("\\.", "_", colnames(anno_df))

out_annotated <- file.path(opt$outdir, paste0(opt$prefix, "_annotated.tsv"))
write.table(anno_df, out_annotated, sep="\t", quote=FALSE, row.names=FALSE)
message("Annotated peaks written to: ", out_annotated)

# ── Annotation summary ────────────────────────────────────────
# Simplify annotation categories
simplify_annotation <- function(ann) {
    ann <- as.character(ann)
    dplyr::case_when(
        grepl("Promoter",    ann, ignore.case=TRUE) ~ "Promoter",
        grepl("5' UTR",      ann, ignore.case=TRUE) ~ "5' UTR",
        grepl("3' UTR",      ann, ignore.case=TRUE) ~ "3' UTR",
        grepl("Exon",        ann, ignore.case=TRUE) ~ "Exon",
        grepl("Intron",      ann, ignore.case=TRUE) ~ "Intron",
        grepl("Downstream",  ann, ignore.case=TRUE) ~ "Downstream",
        grepl("Intergenic",  ann, ignore.case=TRUE) ~ "Intergenic",
        TRUE                                         ~ "Other"
    )
}

anno_df$simplified_annotation <- simplify_annotation(anno_df$annotation)

summary_df <- anno_df %>%
    dplyr::count(simplified_annotation, name="count") %>%
    dplyr::mutate(
        percentage = round(100 * count / sum(count), 2)
    ) %>%
    dplyr::arrange(dplyr::desc(count)) %>%
    dplyr::rename(annotation = simplified_annotation)

out_summary <- file.path(opt$outdir, paste0(opt$prefix, "_annotation_summary.tsv"))
write.table(summary_df, out_summary, sep="\t", quote=FALSE, row.names=FALSE)
message("Annotation summary written to: ", out_summary)

# ── Pie chart ─────────────────────────────────────────────────
annotation_colors <- c(
    "Promoter"    = "#E41A1C",
    "Exon"        = "#377EB8",
    "Intron"      = "#4DAF4A",
    "Intergenic"  = "#984EA3",
    "Downstream"  = "#FF7F00",
    "5' UTR"      = "#A65628",
    "3' UTR"      = "#F781BF",
    "Other"       = "#999999"
)

pie_colors <- annotation_colors[summary_df$annotation]
pie_colors[is.na(pie_colors)] <- "#CCCCCC"

out_pie <- file.path(opt$outdir, paste0(opt$prefix, "_annotation_pie.png"))
png(out_pie, width=800, height=700, res=150)
p_pie <- ggplot(summary_df, aes(x="", y=count, fill=annotation)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    scale_fill_manual(values=pie_colors) +
    labs(
        title = paste0("Peak Annotation Distribution\n", opt$prefix),
        fill  = "Genomic Feature"
    ) +
    theme_void() +
    theme(
        plot.title   = element_text(hjust=0.5, size=12, face="bold"),
        legend.title = element_text(size=10),
        legend.text  = element_text(size=9)
    )
print(p_pie)
dev.off()
message("Pie chart written to: ", out_pie)

# ── Bar chart ─────────────────────────────────────────────────
out_bar <- file.path(opt$outdir, paste0(opt$prefix, "_annotation_bar.png"))
png(out_bar, width=900, height=600, res=150)
p_bar <- ggplot(summary_df,
                aes(x=reorder(annotation, -count), y=percentage, fill=annotation)) +
    geom_bar(stat="identity", color="black", linewidth=0.3) +
    geom_text(aes(label=paste0(percentage, "%")),
              vjust=-0.3, size=3) +
    scale_fill_manual(values=pie_colors) +
    labs(
        title = paste0("Peak Annotation Distribution - ", opt$prefix),
        x     = "Genomic Feature",
        y     = "% of Peaks"
    ) +
    theme_bw(base_size=12) +
    theme(
        axis.text.x  = element_text(angle=45, hjust=1),
        legend.position = "none",
        plot.title   = element_text(hjust=0.5, face="bold")
    )
print(p_bar)
dev.off()
message("Bar chart written to: ", out_bar)

message("Peak annotation complete for: ", opt$prefix)
