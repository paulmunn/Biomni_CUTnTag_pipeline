> [!CAUTION]
> **Research prototype — not validated for scientific use**
>
> This pipeline is an experimental, AI-generated research artifact with known
> implementation and reporting flaws. **It should not be used to process,
> analyze, or interpret experimental data.**
>
> This repository was created as part of a study evaluating AI coding agents
> for bioinformatics pipeline construction. The pipeline required debugging
> during testing and did not fully reproduce the behavior or outputs of the
> hand-coded reference pipeline.
>
> The code is provided only for research, evaluation, and reproducibility of
> the study. Anyone adapting this code should independently validate every
> processing step and output against an appropriate validated workflow.

# CUT&Tag Nextflow DSL2 Pipeline

**Version:** 1.0.0  
**Author:** Genomics Innovation Hub  
**Nextflow DSL:** 2  
**Minimum Nextflow version:** 23.04.0

---

## Table of Contents

1. [Overview](#overview)
2. [Pipeline Architecture](#pipeline-architecture)
3. [Quick Start](#quick-start)
4. [Installation](#installation)
5. [Required Input Files](#required-input-files)
6. [Association CSV Specification](#association-csv-specification)
7. [Example Association CSV](#example-association-csv)
8. [Parameters Reference](#parameters-reference)
9. [Output Directory Structure](#output-directory-structure)
10. [Example Run Commands](#example-run-commands)
11. [Execution Profiles](#execution-profiles)
12. [Testing Strategy](#testing-strategy)
13. [Software Requirements](#software-requirements)
14. [Scientific Recommendations](#scientific-recommendations)
15. [Troubleshooting](#troubleshooting)
16. [Version History](#version-history)

---

## Overview

This pipeline processes paired-end CUT&Tag FASTQ data from raw reads through
quality control, adapter trimming, alignment, multi-stage filtering, signal
track generation, peak calling (per-sample and group-level), peak annotation,
FRiP calculation, deepTools QC visualizations, and integrated MultiQC reports.

**Key features:**
- Handles multi-species, multi-antibody, multi-condition experiments
- Flexible sample grouping via a user-defined association CSV
- Per-sample and group-level (merged replicate) peak calling
- Control-aware MACS3 peak calling with automatic treatment-control matching
- Picard or UMI-tools deduplication
- CPM/RPGC-normalized BigWig and BedGraph signal tracks
- ChIPseeker or HOMER peak annotation with distribution plots
- FRiP scores at per-sample and group levels
- deepTools TSS enrichment and peak-center heatmaps
- Four MultiQC reports (raw, trimmed, alignment, final integrated)
- Resume-friendly Nextflow design with full resource labeling
- Docker, Singularity, Conda, local, and SLURM execution profiles

---

## Pipeline Architecture

```
FASTQ files
    │
    ▼
[PAIR_FASTQS]          Detect R1/R2 pairs, validate, write manifest TSV
    │
    ▼
[VALIDATE_ASSOCIATION] Cross-validate association CSV vs FASTQ manifest
    │
    ▼
[FASTQC_RAW]           FastQC on raw reads
    │
    ▼
[CUTADAPT]             Trim Nextera/Illumina adapters (paired-end)
    │
    ▼
[FASTQC_TRIMMED]       FastQC on trimmed reads
    │
    ▼
[BOWTIE2_ALIGN]        Align to reference genome (--local --very-sensitive-local)
    │
    ▼
[SAMTOOLS_SORT_INDEX]  Sort and index BAM
[SAMTOOLS_FLAGSTAT]    Alignment statistics
[SAMTOOLS_IDXSTATS]    Per-chromosome read counts
[SAMTOOLS_STATS]       Comprehensive alignment stats
    │
    ▼
[FILTER_MITO]          Remove mitochondrial reads + MAPQ filter
    │
    ▼
[FILTER_BLACKLIST]     Remove blacklisted regions (optional)
    │
    ▼
[PICARD_MARKDUPLICATES] Remove PCR/optical duplicates
  or [UMI_TOOLS_DEDUP]  UMI-aware deduplication (if --dedup_mode umi_tools)
    │
    ├──▶ [BAMCOVERAGE_BIGWIG]    CPM/RPGC-normalized BigWig tracks
    ├──▶ [BAMCOVERAGE_BEDGRAPH]  BedGraph tracks
    │
    ├──▶ [MACS3_CALLPEAK_SAMPLE] Per-sample peak calling (with matched control)
    │
    ├──▶ [MERGE_BAMS]            Merge BAMs by merge_group_id
    │       │
    │       ▼
    │   [MACS3_CALLPEAK_GROUP]   Group-level peak calling (merged treatment vs merged control)
    │
    ├──▶ [ANNOTATE_PEAKS]        ChIPseeker or HOMER annotation
    │
    ├──▶ [FRIP_SCORE]            FRiP per sample and per group
    │
    ├──▶ [COMPUTE_MATRIX_TSS]    deepTools matrix around TSS
    │   [COMPUTE_MATRIX_PEAKS]   deepTools matrix around peak centers
    │   [PLOT_PROFILE]           Signal profile plots
    │   [PLOT_HEATMAP]           Signal heatmaps
    │
    └──▶ [MAKE_MQC_CUSTOM]       Generate MultiQC custom content TSVs
         [MULTIQC_RAW]           MultiQC: raw QC
         [MULTIQC_TRIM]          MultiQC: trimming + post-trim QC
         [MULTIQC_ALIGN]         MultiQC: alignment + filtering
         [MULTIQC_FINAL]         MultiQC: integrated final report
```

---

## Quick Start

```bash
# 1. Install Nextflow
curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

# 2. Clone the pipeline
git clone https://github.com/genomics-innovation-hub/cuttag-pipeline.git
cd cuttag-pipeline

# 3. Run with Docker (recommended)
nextflow run main.nf \
    --input_dir     /path/to/fastqs \
    --outdir        results \
    --genome        hg38 \
    --bowtie2_index /path/to/bowtie2/hg38 \
    --chrom_sizes   /path/to/hg38.chrom.sizes \
    --association_csv samples.csv \
    --annotation_gtf /path/to/hg38.gtf \
    --tss_bed       /path/to/hg38_tss.bed \
    --blacklist     /path/to/hg38_blacklist.bed \
    -profile docker
```

---

## Installation

### Option 1: Docker (Recommended)

```bash
# Install Docker: https://docs.docker.com/get-docker/
# Install Nextflow:
curl -s https://get.nextflow.io | bash

# Run pipeline - Docker images are pulled automatically
nextflow run main.nf -profile docker [parameters]
```

### Option 2: Singularity (HPC)

```bash
# Install Singularity: https://sylabs.io/guides/3.0/user-guide/installation.html
# Install Nextflow:
curl -s https://get.nextflow.io | bash

nextflow run main.nf -profile singularity,slurm [parameters]
```

### Option 3: Conda

```bash
# Install Miniconda: https://docs.conda.io/en/latest/miniconda.html
conda env create -f environment.yml
conda activate cuttag-pipeline

nextflow run main.nf -profile conda [parameters]
```

### Option 4: Manual installation

Install all tools listed in [Software Requirements](#software-requirements)
and ensure they are on your `$PATH`.

```bash
nextflow run main.nf -profile local [parameters]
```

### Reference genome files

You will need:
1. **Bowtie2 index** - Build with `bowtie2-build genome.fa index_prefix`
2. **Chromosome sizes** - Generate with `samtools faidx genome.fa && cut -f1,2 genome.fa.fai > genome.chrom.sizes`
3. **Blacklist BED** - Download from [ENCODE blacklists](https://github.com/Boyle-Lab/Blacklist/tree/master/lists)
4. **GTF annotation** - Download from Ensembl or GENCODE
5. **TSS BED** - Extract from GTF: `awk '$3=="transcript"' annotation.gtf | awk '{OFS="\t"; if($7=="+") print $1,$4-1,$4,$9,".",$7; else print $1,$5-1,$5,$9,".",$7}' | sort -k1,1 -k2,2n | uniq > tss.bed`

---

## Required Input Files

| Parameter | Required | Description |
|-----------|----------|-------------|
| `--input_dir` | Yes | Directory containing paired-end FASTQ files |
| `--outdir` | Yes | Output directory |
| `--genome` | Yes | Genome identifier (hg38, mm10, rn6, custom) |
| `--bowtie2_index` | Yes | Path to Bowtie2 index basename or directory |
| `--chrom_sizes` | Yes | Chromosome sizes file (2-column TSV) |
| `--association_csv` | Yes | Sample metadata and grouping table |
| `--annotation_gtf` | For annotation | GTF/GFF annotation file |
| `--tss_bed` | For TSS enrichment | BED file of transcription start sites |
| `--blacklist` | Recommended | BED file of blacklisted genomic regions |

---

## Association CSV Specification

The `--association_csv` file is a comma-separated table that defines sample
metadata, grouping relationships, and control assignments. It is the central
configuration file for the pipeline's group-level analysis logic.

### Column definitions

| Column | Required | Description |
|--------|----------|-------------|
| `sample_id` | Yes | Must exactly match the FASTQ-derived sample ID (filename prefix before `_R1`). |
| `species` | Yes | Species name: `human`, `mouse`, `rat`, `custom`, etc. |
| `genome` | Yes | Genome build: `hg38`, `mm10`, `rn6`, etc. |
| `antibody` | Yes | Target antibody: `H3K27ac`, `H3K4me3`, `IgG`, `Input`, etc. |
| `condition` | Yes | Biological condition, treatment, or cell type. |
| `replicate` | Yes | Biological replicate number or label. Must be unique within a `group_id`. |
| `group_id` | Yes | Unique identifier for this biological sample group (species + antibody + condition + replicate). |
| `is_control` | Yes | `true` if this is a control sample (IgG, Input, etc.), `false` otherwise. |
| `control_group_id` | Yes* | For non-control samples: the `group_id` of the control sample to use for MACS3. Leave empty for control samples. *Required unless `--allow_no_control` is set. |
| `merge_group_id` | Yes | Defines which samples are merged for group-level peak calling. All samples with the same `merge_group_id` must share `species`, `genome`, `antibody`, and `condition`. |
| `peak_calling_mode` | Yes | `narrow` (TFs, H3K4me3, H3K27ac), `broad` (H3K27me3, H3K9me3, H3K36me3), or `auto`. |
| `notes` | No | Free-text notes. |

### Validation rules

The pipeline enforces these rules at startup and will **fail with a clear error** if any are violated:

1. All required columns must be present.
2. No duplicate `sample_id` values.
3. Every `sample_id` in the FASTQ pairs manifest must exist in the association CSV.
4. No extra `sample_id` rows (unless `--allow_extra` is passed to the validator).
5. `is_control` must be exactly `true` or `false`.
6. `peak_calling_mode` must be `narrow`, `broad`, or `auto`.
7. All rows sharing a `merge_group_id` must have identical `species`, `genome`, `antibody`, and `condition`.
8. Every non-control sample must have a `control_group_id` that resolves to at least one control sample's `group_id` (unless `--allow_no_control`).
9. `replicate` values must be unique within each `group_id`.

---

## Example Association CSV

```csv
sample_id,species,genome,antibody,condition,replicate,group_id,is_control,control_group_id,merge_group_id,peak_calling_mode,notes
H3K27ac_MCF7_rep1,human,hg38,H3K27ac,MCF7,1,H3K27ac_MCF7_rep1,false,IgG_MCF7_group,human_H3K27ac_MCF7,narrow,H3K27ac replicate 1
H3K27ac_MCF7_rep2,human,hg38,H3K27ac,MCF7,2,H3K27ac_MCF7_rep2,false,IgG_MCF7_group,human_H3K27ac_MCF7,narrow,H3K27ac replicate 2
H3K27ac_MCF7_rep3,human,hg38,H3K27ac,MCF7,3,H3K27ac_MCF7_rep3,false,IgG_MCF7_group,human_H3K27ac_MCF7,narrow,H3K27ac replicate 3
H3K4me3_MCF7_rep1,human,hg38,H3K4me3,MCF7,1,H3K4me3_MCF7_rep1,false,IgG_MCF7_group,human_H3K4me3_MCF7,narrow,H3K4me3 replicate 1
H3K4me3_MCF7_rep2,human,hg38,H3K4me3,MCF7,2,H3K4me3_MCF7_rep2,false,IgG_MCF7_group,human_H3K4me3_MCF7,narrow,H3K4me3 replicate 2
H3K27me3_MCF7_rep1,human,hg38,H3K27me3,MCF7,1,H3K27me3_MCF7_rep1,false,IgG_MCF7_group,human_H3K27me3_MCF7,broad,H3K27me3 broad mark
H3K27me3_MCF7_rep2,human,hg38,H3K27me3,MCF7,2,H3K27me3_MCF7_rep2,false,IgG_MCF7_group,human_H3K27me3_MCF7,broad,H3K27me3 broad mark
IgG_MCF7_rep1,human,hg38,IgG,MCF7,1,IgG_MCF7_group,true,,human_IgG_MCF7,narrow,IgG negative control
H3K27ac_K562_rep1,human,hg38,H3K27ac,K562,1,H3K27ac_K562_rep1,false,IgG_K562_group,human_H3K27ac_K562,narrow,K562 H3K27ac rep 1
H3K27ac_K562_rep2,human,hg38,H3K27ac,K562,2,H3K27ac_K562_rep2,false,IgG_K562_group,human_H3K27ac_K562,narrow,K562 H3K27ac rep 2
IgG_K562_rep1,human,hg38,IgG,K562,1,IgG_K562_group,true,,human_IgG_K562,narrow,K562 IgG control
H3K27ac_MEF_rep1,mouse,mm10,H3K27ac,MEF,1,H3K27ac_MEF_rep1,false,IgG_MEF_group,mouse_H3K27ac_MEF,narrow,Mouse MEF H3K27ac
H3K27ac_MEF_rep2,mouse,mm10,H3K27ac,MEF,2,H3K27ac_MEF_rep2,false,IgG_MEF_group,mouse_H3K27ac_MEF,narrow,Mouse MEF H3K27ac
IgG_MEF_rep1,mouse,mm10,IgG,MEF,1,IgG_MEF_group,true,,mouse_IgG_MEF,narrow,Mouse MEF IgG control
```

**Key points illustrated:**
- Multiple antibodies per cell line (H3K27ac, H3K4me3, H3K27me3)
- Multiple cell lines (MCF7, K562) with separate IgG controls
- Multi-species experiment (human + mouse)
- Broad peak mode for H3K27me3
- All replicates of the same antibody+condition share a `merge_group_id`

---

## Parameters Reference

### Required

| Parameter | Description |
|-----------|-------------|
| `--input_dir` | Directory containing paired-end FASTQ files |
| `--outdir` | Output directory |
| `--genome` | Genome identifier (hg38, mm10, rn6, custom) |
| `--bowtie2_index` | Path to Bowtie2 index basename |
| `--chrom_sizes` | Chromosome sizes file |
| `--association_csv` | Sample association CSV |

### Optional

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--blacklist` | null | BED file of blacklisted regions |
| `--annotation_gtf` | null | GTF annotation file |
| `--tss_bed` | null | TSS BED file |
| `--adapter_fwd` | `CTGTCTCTTATACACATCT` | Forward adapter (Nextera) |
| `--adapter_rev` | `CTGTCTCTTATACACATCT` | Reverse adapter (Nextera) |
| `--min_length` | 20 | Minimum read length after trimming |
| `--bowtie2_args` | `--local --very-sensitive-local ...` | Bowtie2 alignment arguments |
| `--mapq` | 20 | Minimum MAPQ score |
| `--mito_name` | `chrM\|MT\|M\|mitochondria\|Mito` | Mitochondrial chromosome name pattern |
| `--dedup_mode` | `picard` | Deduplication: `picard` or `umi_tools` |
| `--skip_umi` | true | Skip UMI-aware deduplication |
| `--skip_blacklist` | false | Skip blacklist filtering |
| `--bigwig_binsize` | 10 | BigWig bin size (bp) |
| `--normalization` | `CPM` | BigWig normalization: CPM, RPGC, BPM, RPKM, None |
| `--effective_genome_size` | 2913022398 | Effective genome size (for RPGC; hg38 default) |
| `--macs3_qvalue` | 0.05 | MACS3 q-value threshold |
| `--macs3_genome` | `hs` | MACS3 genome size: hs, mm, ce, dm |
| `--allow_no_control` | false | Allow peak calling without control |
| `--tss_window` | 2000 | TSS window size (bp each side) |
| `--peak_center_window` | 1000 | Peak center window (bp each side) |
| `--threads` | 4 | Default number of threads |
| `--annotation_tool` | `chipseeker` | Peak annotation: `chipseeker` or `homer` |
| `--paired_pattern` | `auto` | FASTQ pairing pattern |

---

## Output Directory Structure

```
results/
├── 00_fastq_pairs/
│   ├── fastq_pairs.tsv              # Detected R1/R2 pairs manifest
│   ├── pairing_report.txt           # Human-readable pairing report
│   ├── association_validated.csv    # Validated association CSV
│   └── validation_report.txt        # Validation report
│
├── 01_fastqc_raw/                   # FastQC HTML/ZIP for raw reads
├── 02_trimmed/                      # Trimmed FASTQ files + Cutadapt logs
├── 03_fastqc_trimmed/               # FastQC HTML/ZIP for trimmed reads
│
├── 04_alignment/
│   ├── *.bam, *.bam.bai             # Sorted, indexed BAM files
│   ├── logs/                        # Bowtie2 alignment logs
│   └── stats/                       # flagstat, idxstats, stats files
│
├── 05_filtering/
│   ├── mito_removed/                # BAMs after mito removal
│   ├── blacklist_removed/           # BAMs after blacklist filtering
│   ├── dedup/                       # Deduplicated BAMs
│   ├── stats/                       # Per-stage filtering stats
│   └── read_retention/              # Per-sample read retention TSVs
│
├── 06_bigwig/                       # Normalized BigWig tracks
├── 07_bedgraph/                     # BedGraph tracks
│
├── 08_peaks_per_sample/
│   └── {sample_id}/
│       ├── *_peaks.narrowPeak       # Per-sample narrow peaks
│       ├── *_peaks.broadPeak        # Per-sample broad peaks (if applicable)
│       ├── *_summits.bed            # Peak summits
│       ├── *_peaks.xls              # MACS3 peak table
│       ├── *_macs3.log              # MACS3 log
│       └── *_peak_stats.txt         # Peak count summary
│
├── 09_peaks_merged_groups/
│   └── {merge_group_id}/
│       ├── *_merged.bam             # Merged treatment BAM
│       ├── *_peaks.narrowPeak       # Group-level peaks
│       ├── *_summits.bed
│       └── *_peak_stats.txt
│
├── 10_peak_annotation/
│   └── {peak_id}/
│       ├── *_annotated.tsv          # Full annotation table
│       ├── *_annotation_summary.tsv # Category distribution
│       ├── *_annotation_pie.png     # Pie chart
│       └── *_annotation_bar.png     # Bar chart
│
├── 11_frip/
│   └── *.frip.txt                   # FRiP scores per sample
│
├── 12_deeptools/
│   ├── matrices/                    # computeMatrix output (.gz)
│   ├── profiles/                    # plotProfile PNG/PDF
│   ├── heatmaps/                    # plotHeatmap PNG/PDF
│   └── tss_enrichment/              # TSS enrichment matrices
│
├── 13_multiqc/
│   ├── 01_raw/                      # Raw reads MultiQC report
│   ├── 02_trimmed/                  # Post-trimming MultiQC report
│   ├── 03_alignment/                # Alignment MultiQC report
│   ├── 04_final/                    # Final integrated MultiQC report
│   └── custom_content/              # MultiQC custom TSV files
│
└── pipeline_info/
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    ├── execution_trace_*.txt
    ├── pipeline_dag_*.html
    ├── software_versions.yml
    └── pipeline_info.txt
```

---

## Example Run Commands

### Basic run with Docker

```bash
nextflow run main.nf \
    --input_dir     /data/fastq \
    --outdir        results \
    --genome        hg38 \
    --bowtie2_index /ref/hg38/bowtie2/hg38 \
    --chrom_sizes   /ref/hg38/hg38.chrom.sizes \
    --blacklist     /ref/hg38/hg38_blacklist.bed \
    --annotation_gtf /ref/hg38/gencode.v44.annotation.gtf \
    --tss_bed       /ref/hg38/hg38_tss.bed \
    --association_csv samples.csv \
    -profile docker \
    -resume
```

### SLURM HPC with Singularity

```bash
nextflow run main.nf \
    --input_dir     /scratch/fastq \
    --outdir        /scratch/results \
    --genome        mm10 \
    --bowtie2_index /ref/mm10/bowtie2/mm10 \
    --chrom_sizes   /ref/mm10/mm10.chrom.sizes \
    --blacklist     /ref/mm10/mm10_blacklist.bed \
    --annotation_gtf /ref/mm10/gencode.vM33.annotation.gtf \
    --tss_bed       /ref/mm10/mm10_tss.bed \
    --association_csv samples.csv \
    --macs3_genome  mm \
    --effective_genome_size 2652783500 \
    -profile singularity,slurm \
    --slurm_queue   normal \
    -resume
```

### With UMI-aware deduplication

```bash
nextflow run main.nf \
    [standard parameters] \
    --dedup_mode    umi_tools \
    --skip_umi      false \
    --umi_separator ":" \
    -profile docker
```

### Broad peak mode (H3K27me3)

```bash
# Set peak_calling_mode=broad in association_csv, or override globally:
nextflow run main.nf \
    [standard parameters] \
    --macs3_qvalue  0.05 \
    -profile docker
# (peak_calling_mode is set per-sample in association_csv)
```

### Resume a failed run

```bash
nextflow run main.nf [parameters] -resume
```

### Dry run (stub mode)

```bash
nextflow run main.nf [parameters] -stub-run
```

---

## Execution Profiles

| Profile | Description |
|---------|-------------|
| `local` | Run on local machine, no containers |
| `docker` | Run with Docker containers (recommended for reproducibility) |
| `singularity` | Run with Singularity containers (recommended for HPC) |
| `conda` | Run with Conda environment |
| `slurm` | Submit jobs to SLURM scheduler (combine with docker/singularity) |
| `test` | Minimal test run with synthetic data |

Combine profiles with commas: `-profile singularity,slurm`

---

## Testing Strategy

### Minimal synthetic test

```bash
# Generate synthetic test data
python bin/make_test_data.py --outdir test_data --n_reads 1000

# Run pipeline in stub mode (no actual tool execution)
nextflow run main.nf -profile test -stub-run

# Run pipeline with real tools on synthetic data
nextflow run main.nf -profile test,docker
```

**Note:** Synthetic random reads will not produce meaningful alignments or peaks.
The test profile validates pipeline structure, channel logic, and file I/O.

### Biological validation test

For biological validation, use published CUT&Tag data:

```bash
# Download test data from GEO (Kaya-Okur et al. 2019, GSE145187)
# H3K27ac CUT&Tag in K562 cells (small subset recommended)
# SRR11074290 (H3K27ac), SRR11074291 (IgG control)

# Download with SRA tools:
fastq-dump --split-files --gzip SRR11074290
fastq-dump --split-files --gzip SRR11074291

# Rename to match expected pattern:
mv SRR11074290_1.fastq.gz H3K27ac_K562_rep1_R1.fastq.gz
mv SRR11074290_2.fastq.gz H3K27ac_K562_rep1_R2.fastq.gz
mv SRR11074291_1.fastq.gz IgG_K562_rep1_R1.fastq.gz
mv SRR11074291_2.fastq.gz IgG_K562_rep1_R2.fastq.gz
```

Expected results for H3K27ac K562 CUT&Tag:
- Alignment rate: >90%
- Mitochondrial fraction: <5%
- Duplicate rate: <30%
- FRiP score: >0.2
- Peak count: >10,000 narrow peaks

---

## Software Requirements

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | ≥23.04.0 | Workflow manager |
| FastQC | 0.12.1 | Raw and trimmed read QC |
| MultiQC | 1.21 | QC report aggregation |
| Cutadapt | 4.7 | Adapter trimming |
| Bowtie2 | 2.5.3 | Read alignment |
| SAMtools | 1.19.2 | BAM manipulation |
| Picard | 3.1.1 | Duplicate marking |
| UMI-tools | 1.1.5 | UMI-aware deduplication (optional) |
| BEDTools | 2.31.1 | Genomic interval operations |
| deepTools | 3.5.4 | Signal tracks and QC visualizations |
| MACS3 | 3.0.1 | Peak calling |
| R | 4.3.3 | Statistical computing |
| ChIPseeker | 1.38.0 | Peak annotation (Bioconductor) |
| HOMER | 4.11+ | Alternative peak annotation (optional) |
| Python | 3.11 | Helper scripts |

---

## Scientific Recommendations

### Adapter sequences
The default adapter sequence (`CTGTCTCTTATACACATCT`) is the Nextera transposase
adapter used in most CUT&Tag protocols. If you used a different library prep
(e.g., TruSeq), override with `--adapter_fwd` and `--adapter_rev`.

### Alignment parameters
`--local --very-sensitive-local` is recommended for CUT&Tag because:
- Local alignment tolerates soft-clipping at read ends (common after tagmentation)
- Fragment size range `-I 10 -X 700` captures nucleosome-free (NFR) and
  mono-nucleosomal fragments typical of CUT&Tag

### Deduplication
Picard MarkDuplicates is the default. For low-input CUT&Tag experiments with
very high duplication rates, consider whether deduplication is appropriate —
some protocols recommend keeping duplicates for low-input samples.

### Peak calling mode
- **Narrow peaks** (`--call-summits`): H3K4me3, H3K27ac, H3K4me1, TF binding
- **Broad peaks**: H3K27me3, H3K9me3, H3K36me3, H3K9ac (broad domains)
- Set `peak_calling_mode` per sample in the association CSV

### Normalization
- **CPM** (default): Good for comparing samples within an experiment
- **RPGC**: Better for cross-experiment comparisons; requires accurate
  `--effective_genome_size`
- Effective genome sizes: hg38=2913022398, mm10=2652783500, rn6=2729860805

### FRiP thresholds
ENCODE recommends FRiP ≥ 0.01 as a minimum. High-quality CUT&Tag experiments
typically achieve FRiP > 0.2 for sharp marks (H3K4me3, H3K27ac).

### Control samples
IgG is the standard negative control for CUT&Tag. Input DNA is less common
but also supported. Each treatment group should have a matched control from
the same cell line/condition.

---

## Troubleshooting

### "No FASTQ files found"
- Check that `--input_dir` contains `.fastq.gz` or `.fq.gz` files
- Verify file permissions
- Check that files follow a supported naming pattern (see `bin/pair_fastqs.py`)

### "Sample X found in FASTQ pairs but missing from association CSV"
- The sample ID derived from the FASTQ filename must exactly match `sample_id`
  in the association CSV
- The sample ID is derived by stripping the R1 suffix and Illumina lane tags
- Check `results/00_fastq_pairs/fastq_pairs.tsv` to see the derived sample IDs

### "control_group_id does not match any control sample"
- Verify that the `control_group_id` in treatment rows matches the `group_id`
  of a control sample (not the `sample_id`)
- Check spelling and case sensitivity

### Bowtie2 index not found
- Provide the full path to the index basename (without `.1.bt2` extension)
- Or provide the directory containing the index files

### Out of memory errors
- Increase memory in `conf/base.config` or use `--max_memory`
- For SLURM, adjust `conf/slurm.config`

### Resume not working
- Ensure the `-resume` flag is used
- Check that the work directory has not been deleted
- Verify that input files have not been modified

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0.0 | 2026-04-27 | Initial release |

---

## Citation

If you use this pipeline, please cite the tools it uses:

- **Bowtie2**: Langmead B, Salzberg SL. *Nat Methods* 2012
- **MACS3**: Zhang Y et al. *Genome Biol* 2008; updated MACS3 2023
- **deepTools**: Ramírez F et al. *Nucleic Acids Res* 2016
- **Picard**: Broad Institute. https://broadinstitute.github.io/picard/
- **ChIPseeker**: Yu G et al. *Bioinformatics* 2015
- **Cutadapt**: Martin M. *EMBnet.journal* 2011
- **FastQC**: Andrews S. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
- **MultiQC**: Ewels P et al. *Bioinformatics* 2016
- **SAMtools**: Li H et al. *Bioinformatics* 2009
- **BEDTools**: Quinlan AR, Hall IM. *Bioinformatics* 2010

For the CUT&Tag protocol itself:
- **CUT&Tag**: Kaya-Okur HS et al. *Nat Commun* 2019. doi:10.1038/s41467-019-09982-5
