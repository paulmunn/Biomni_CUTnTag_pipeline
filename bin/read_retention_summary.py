#!/usr/bin/env python3
"""
read_retention_summary.py
=========================
Parse per-sample QC files and produce a read-retention TSV showing
how many reads were retained at each filtering stage.

Stages tracked:
  raw_reads          - total reads from flagstat (pre-filter)
  properly_paired    - properly paired reads from flagstat
  mito_reads         - mitochondrial reads removed
  nomito_reads       - reads after mito removal
  blacklist_removed  - reads removed by blacklist (if available)
  after_blacklist    - reads after blacklist filtering
  duplicates         - duplicate reads removed (Picard or UMI-tools)
  final_reads        - final usable reads after all filtering
"""

import argparse
import re
import sys
from pathlib import Path


def parse_flagstat(flagstat_path: str) -> dict:
    """Parse samtools flagstat output."""
    stats = {
        'total_reads':      0,
        'mapped_reads':     0,
        'properly_paired':  0,
        'secondary':        0,
        'supplementary':    0,
    }
    try:
        with open(flagstat_path) as fh:
            for line in fh:
                line = line.strip()
                # Extract the primary count (before the '+')
                m = re.match(r'^(\d+) \+ \d+ (.+)$', line)
                if not m:
                    continue
                count = int(m.group(1))
                desc  = m.group(2).lower()

                if 'in total' in desc:
                    stats['total_reads'] = count
                elif 'mapped' in desc and 'primary mapped' not in desc and '%' not in desc:
                    stats['mapped_reads'] = count
                elif 'properly paired' in desc:
                    stats['properly_paired'] = count
                elif 'secondary' in desc:
                    stats['secondary'] = count
                elif 'supplementary' in desc:
                    stats['supplementary'] = count
    except FileNotFoundError:
        pass
    return stats


def parse_mito_stats(mito_stats_path: str) -> dict:
    """Parse mito_stats.txt produced by FILTER_MITO."""
    stats = {
        'total_reads':   0,
        'mito_reads':    0,
        'nomito_reads':  0,
        'mito_fraction': 0.0,
    }
    try:
        with open(mito_stats_path) as fh:
            lines = fh.readlines()
            if len(lines) >= 2:
                parts = lines[1].strip().split('\t')
                if len(parts) >= 5:
                    stats['total_reads']   = int(parts[1])
                    stats['mito_reads']    = int(parts[2])
                    stats['nomito_reads']  = int(parts[3])
                    try:
                        stats['mito_fraction'] = float(parts[4])
                    except ValueError:
                        stats['mito_fraction'] = 0.0
    except FileNotFoundError:
        pass
    return stats


def parse_picard_metrics(metrics_path: str) -> dict:
    """Parse Picard MarkDuplicates metrics file."""
    stats = {
        'estimated_library_size': 0,
        'percent_duplication':    0.0,
        'read_pairs_examined':    0,
        'read_pair_duplicates':   0,
        'reads_after_dedup':      0,
    }
    if not metrics_path or metrics_path == 'None':
        return stats
    try:
        with open(metrics_path) as fh:
            in_metrics = False
            header     = None
            for line in fh:
                line = line.strip()
                if line.startswith('## METRICS CLASS'):
                    in_metrics = True
                    continue
                if in_metrics and header is None and line:
                    header = line.split('\t')
                    continue
                if in_metrics and header and line:
                    values = line.split('\t')
                    row    = dict(zip(header, values))
                    try:
                        stats['read_pairs_examined']  = int(row.get('READ_PAIRS_EXAMINED', 0))
                        stats['read_pair_duplicates'] = int(row.get('READ_PAIR_DUPLICATES', 0))
                        stats['percent_duplication']  = float(row.get('PERCENT_DUPLICATION', 0))
                        stats['estimated_library_size'] = int(row.get('ESTIMATED_LIBRARY_SIZE', 0))
                    except (ValueError, KeyError):
                        pass
                    break
    except FileNotFoundError:
        pass
    return stats


def parse_umi_dedup_log(log_path: str) -> dict:
    """Parse UMI-tools dedup log for duplicate counts."""
    stats = {
        'read_pairs_examined':  0,
        'read_pair_duplicates': 0,
        'percent_duplication':  0.0,
    }
    if not log_path or log_path == 'None':
        return stats
    try:
        with open(log_path) as fh:
            for line in fh:
                m = re.search(r'Input Reads: (\d+)', line)
                if m:
                    stats['read_pairs_examined'] = int(m.group(1))
                m = re.search(r'Number of reads out: (\d+)', line)
                if m:
                    out = int(m.group(1))
                    stats['read_pair_duplicates'] = (
                        stats['read_pairs_examined'] - out
                    )
                    if stats['read_pairs_examined'] > 0:
                        stats['percent_duplication'] = (
                            stats['read_pair_duplicates'] /
                            stats['read_pairs_examined']
                        )
    except FileNotFoundError:
        pass
    return stats


def safe_pct(numerator: int, denominator: int) -> float:
    if denominator == 0:
        return 0.0
    return round(100.0 * numerator / denominator, 2)


def main():
    parser = argparse.ArgumentParser(
        description='Summarize read retention across filtering stages.'
    )
    parser.add_argument('--sample_id',   required=True)
    parser.add_argument('--flagstat',    required=True)
    parser.add_argument('--mito_stats',  required=True)
    parser.add_argument('--dup_metrics', default=None,
                        help='Picard MarkDuplicates metrics or UMI-tools log')
    parser.add_argument('--output',      required=True)
    args = parser.parse_args()

    flagstat  = parse_flagstat(args.flagstat)
    mito      = parse_mito_stats(args.mito_stats)
    dup_path  = args.dup_metrics

    # Detect Picard vs UMI-tools by file extension
    if dup_path and dup_path != 'None':
        if dup_path.endswith('.umi_dedup.log'):
            dup = parse_umi_dedup_log(dup_path)
        else:
            dup = parse_picard_metrics(dup_path)
    else:
        dup = {
            'read_pairs_examined':  0,
            'read_pair_duplicates': 0,
            'percent_duplication':  0.0,
            'estimated_library_size': 0,
        }

    # Compute final reads
    total_reads     = flagstat['total_reads']
    properly_paired = flagstat['properly_paired']
    mito_reads      = mito['mito_reads']
    nomito_reads    = mito['nomito_reads']
    dup_reads       = dup['read_pair_duplicates'] * 2  # pairs -> reads
    final_reads     = nomito_reads - dup_reads if nomito_reads > 0 else 0

    # Write TSV
    header = [
        'sample_id',
        'total_reads', 'pct_total',
        'properly_paired', 'pct_properly_paired',
        'mito_reads', 'pct_mito',
        'nomito_reads', 'pct_nomito',
        'duplicate_reads', 'pct_duplicates',
        'final_reads', 'pct_final',
        'estimated_library_size',
    ]

    row = [
        args.sample_id,
        total_reads,     100.0,
        properly_paired, safe_pct(properly_paired, total_reads),
        mito_reads,      safe_pct(mito_reads, total_reads),
        nomito_reads,    safe_pct(nomito_reads, total_reads),
        dup_reads,       round(dup['percent_duplication'] * 100, 2),
        final_reads,     safe_pct(final_reads, total_reads),
        dup.get('estimated_library_size', 0),
    ]

    with open(args.output, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        fh.write('\t'.join(str(v) for v in row) + '\n')

    print(
        f"[{args.sample_id}] total={total_reads:,} | "
        f"properly_paired={properly_paired:,} | "
        f"mito={mito_reads:,} ({safe_pct(mito_reads, total_reads):.1f}%) | "
        f"dups={dup_reads:,} ({dup['percent_duplication']*100:.1f}%) | "
        f"final={final_reads:,} ({safe_pct(final_reads, total_reads):.1f}%)",
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
