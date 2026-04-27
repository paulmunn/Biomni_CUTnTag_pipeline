#!/usr/bin/env python3
"""
make_mqc_custom.py
==================
Aggregate per-sample QC files and generate MultiQC-compatible
custom content TSV files for:
  - Read retention through pipeline stages
  - FRiP scores
  - Peak counts
  - Mitochondrial read fractions
  - Annotation distributions (if available)

MultiQC custom content format:
  Files named *_mqc.tsv with a YAML header block are auto-detected
  by MultiQC as custom content sections.
"""

import argparse
import csv
import glob
import os
import sys
from pathlib import Path


# ── MultiQC custom content YAML headers ──────────────────────────────────────

READ_RETENTION_HEADER = """# id: 'read_retention'
# section_name: 'Read Retention Through Pipeline'
# description: 'Number of reads retained at each filtering stage.'
# format: 'tsv'
# plot_type: 'bargraph'
# pconfig:
#     id: 'read_retention_plot'
#     title: 'Read Retention'
#     ylab: 'Read Count'
#     stacking: 'normal'
"""

FRIP_HEADER = """# id: 'frip_scores'
# section_name: 'FRiP Scores'
# description: 'Fraction of Reads in Peaks (FRiP) per sample. ENCODE recommends >= 0.01.'
# format: 'tsv'
# plot_type: 'bargraph'
# pconfig:
#     id: 'frip_plot'
#     title: 'FRiP Scores'
#     ylab: 'FRiP Score'
#     ymax: 1.0
#     ymin: 0.0
"""

PEAK_COUNTS_HEADER = """# id: 'peak_counts'
# section_name: 'Peak Counts'
# description: 'Number of peaks called per sample.'
# format: 'tsv'
# plot_type: 'bargraph'
# pconfig:
#     id: 'peak_count_plot'
#     title: 'Peak Counts'
#     ylab: 'Number of Peaks'
"""

MITO_HEADER = """# id: 'mito_fraction'
# section_name: 'Mitochondrial Read Fraction'
# description: 'Percentage of reads mapping to mitochondrial chromosome.'
# format: 'tsv'
# plot_type: 'bargraph'
# pconfig:
#     id: 'mito_plot'
#     title: 'Mitochondrial Read Fraction'
#     ylab: '% Reads'
#     ymax: 100
"""

FLAGSTAT_SUMMARY_HEADER = """# id: 'alignment_summary'
# section_name: 'Alignment Summary'
# description: 'Summary of alignment statistics per sample.'
# format: 'tsv'
# plot_type: 'table'
"""


def load_tsv_files(directory: str, pattern: str = '*.tsv') -> list:
    """Load all TSV files matching pattern from directory."""
    rows = []
    for fpath in sorted(glob.glob(os.path.join(directory, pattern))):
        try:
            with open(fpath, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    rows.append(row)
        except Exception as e:
            print(f"WARNING: Could not read {fpath}: {e}", file=sys.stderr)
    return rows


def load_frip_files(directory: str) -> list:
    """Load all FRiP TSV files."""
    return load_tsv_files(directory, '*.frip.txt')


def load_peak_stats(directory: str) -> list:
    """Load all peak stats TSV files."""
    rows = []
    for fpath in sorted(glob.glob(os.path.join(directory, '**/*_peak_stats.txt'),
                                  recursive=True)):
        try:
            with open(fpath, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    rows.append(row)
        except Exception as e:
            print(f"WARNING: Could not read {fpath}: {e}", file=sys.stderr)
    return rows


def parse_flagstat_file(fpath: str) -> dict:
    """Parse a single flagstat file."""
    import re
    stats = {'total': 0, 'mapped': 0, 'properly_paired': 0}
    try:
        with open(fpath) as fh:
            for line in fh:
                m = re.match(r'^(\d+) \+ \d+ (.+)$', line.strip())
                if not m:
                    continue
                count = int(m.group(1))
                desc  = m.group(2).lower()
                if 'in total' in desc:
                    stats['total'] = count
                elif 'mapped' in desc and '%' not in desc and 'primary' not in desc:
                    stats['mapped'] = count
                elif 'properly paired' in desc:
                    stats['properly_paired'] = count
    except Exception:
        pass
    return stats


def write_mqc_tsv(header_yaml: str, column_headers: list, rows: list,
                  output_path: str):
    """Write a MultiQC-compatible TSV with YAML header."""
    with open(output_path, 'w') as fh:
        fh.write(header_yaml)
        fh.write('\t'.join(column_headers) + '\n')
        for row in rows:
            fh.write('\t'.join(str(row.get(col, '')) for col in column_headers) + '\n')
    print(f"Written: {output_path}", file=sys.stderr)


def make_read_retention_mqc(retention_dir: str, outdir: str):
    """Generate read retention MultiQC TSV."""
    rows = load_tsv_files(retention_dir, '*.read_retention.tsv')
    if not rows:
        print("WARNING: No read retention files found.", file=sys.stderr)
        return

    # MultiQC bargraph expects: Sample | Category1 | Category2 ...
    # We'll output stacked bar: properly_paired, mito_removed, dup_removed, final
    mqc_rows = []
    for row in rows:
        sid = row.get('sample_id', 'unknown')
        try:
            total    = int(row.get('total_reads', 0))
            pp       = int(row.get('properly_paired', 0))
            mito     = int(row.get('mito_reads', 0))
            dups     = int(row.get('duplicate_reads', 0))
            final    = int(row.get('final_reads', 0))
            # Reads removed at each stage (for stacked bar)
            unpaired = total - pp
            mqc_rows.append({
                'Sample':           sid,
                'Final_Reads':      final,
                'Duplicates_Removed': dups,
                'Mito_Removed':     mito,
                'Unpaired_Removed': unpaired,
            })
        except (ValueError, TypeError) as e:
            print(f"WARNING: Could not parse retention row for {sid}: {e}",
                  file=sys.stderr)

    if mqc_rows:
        write_mqc_tsv(
            READ_RETENTION_HEADER,
            ['Sample', 'Final_Reads', 'Duplicates_Removed', 'Mito_Removed', 'Unpaired_Removed'],
            mqc_rows,
            os.path.join(outdir, 'read_retention_mqc.tsv')
        )

    # Also write mito fraction separately
    mito_rows = []
    for row in rows:
        sid = row.get('sample_id', 'unknown')
        try:
            mito_pct = float(row.get('pct_mito', 0))
            mito_rows.append({'Sample': sid, 'Mito_Fraction_Pct': round(mito_pct, 2)})
        except (ValueError, TypeError):
            pass

    if mito_rows:
        write_mqc_tsv(
            MITO_HEADER,
            ['Sample', 'Mito_Fraction_Pct'],
            mito_rows,
            os.path.join(outdir, 'mito_fraction_mqc.tsv')
        )


def make_frip_mqc(frip_dir: str, outdir: str):
    """Generate FRiP scores MultiQC TSV."""
    rows = load_frip_files(frip_dir)
    if not rows:
        print("WARNING: No FRiP files found.", file=sys.stderr)
        return

    mqc_rows = []
    for row in rows:
        sid = row.get('sample_id', 'unknown')
        try:
            frip = float(row.get('frip_score', 0))
            mqc_rows.append({'Sample': sid, 'FRiP_Score': round(frip, 4)})
        except (ValueError, TypeError):
            pass

    if mqc_rows:
        write_mqc_tsv(
            FRIP_HEADER,
            ['Sample', 'FRiP_Score'],
            mqc_rows,
            os.path.join(outdir, 'frip_scores_mqc.tsv')
        )


def make_peak_counts_mqc(peaks_dir: str, outdir: str):
    """Generate peak counts MultiQC TSV."""
    rows = load_peak_stats(peaks_dir)
    if not rows:
        print("WARNING: No peak stats files found.", file=sys.stderr)
        return

    mqc_rows = []
    for row in rows:
        # Handle both per-sample (sample_id) and group-level (merge_group_id)
        sid = row.get('sample_id') or row.get('merge_group_id', 'unknown')
        try:
            n_peaks = int(row.get('n_peaks', 0))
            mqc_rows.append({'Sample': sid, 'N_Peaks': n_peaks})
        except (ValueError, TypeError):
            pass

    if mqc_rows:
        write_mqc_tsv(
            PEAK_COUNTS_HEADER,
            ['Sample', 'N_Peaks'],
            mqc_rows,
            os.path.join(outdir, 'peak_counts_mqc.tsv')
        )


def make_flagstat_summary_mqc(flagstat_dir: str, outdir: str):
    """Generate alignment summary MultiQC TSV from flagstat files."""
    import re
    mqc_rows = []
    for fpath in sorted(glob.glob(os.path.join(flagstat_dir, '*.flagstat'))):
        sample_id = Path(fpath).stem.replace('.flagstat', '')
        stats = parse_flagstat_file(fpath)
        total = stats['total']
        if total > 0:
            mqc_rows.append({
                'Sample':              sample_id,
                'Total_Reads':         total,
                'Mapped_Reads':        stats['mapped'],
                'Pct_Mapped':          round(100.0 * stats['mapped'] / total, 2),
                'Properly_Paired':     stats['properly_paired'],
                'Pct_Properly_Paired': round(100.0 * stats['properly_paired'] / total, 2),
            })

    if mqc_rows:
        write_mqc_tsv(
            FLAGSTAT_SUMMARY_HEADER,
            ['Sample', 'Total_Reads', 'Mapped_Reads', 'Pct_Mapped',
             'Properly_Paired', 'Pct_Properly_Paired'],
            mqc_rows,
            os.path.join(outdir, 'alignment_summary_mqc.tsv')
        )


def main():
    parser = argparse.ArgumentParser(
        description='Generate MultiQC custom content files.'
    )
    parser.add_argument('--retention_dir', required=True,
                        help='Directory with *.read_retention.tsv files')
    parser.add_argument('--frip_dir',      required=True,
                        help='Directory with *.frip.txt files')
    parser.add_argument('--peaks_dir',     required=True,
                        help='Directory with *_peak_stats.txt files')
    parser.add_argument('--flagstat_dir',  required=True,
                        help='Directory with *.flagstat files')
    parser.add_argument('--outdir',        required=True,
                        help='Output directory for MultiQC TSV files')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    make_read_retention_mqc(args.retention_dir, args.outdir)
    make_frip_mqc(args.frip_dir, args.outdir)
    make_peak_counts_mqc(args.peaks_dir, args.outdir)
    make_flagstat_summary_mqc(args.flagstat_dir, args.outdir)

    print("MultiQC custom content generation complete.", file=sys.stderr)


if __name__ == '__main__':
    main()
