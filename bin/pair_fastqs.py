#!/usr/bin/env python3
"""
pair_fastqs.py
==============
Recursively scan an input directory (and all subdirectories) for
paired-end FASTQ files, detect R1/R2 pairs, validate that every R1
has a matching R2, and write a TSV manifest.

Directory traversal:
  Uses Path.rglob() — files at ANY nesting depth are discovered.
  Common layouts supported:
    input_dir/sample1_R1.fastq.gz          (flat)
    input_dir/sample1/sample1_R1.fastq.gz  (one level per sample)
    input_dir/run1/lane1/sample_R1.fq.gz   (deep nesting)

Supported R1/R2 naming patterns (auto-detected):
  _R1_001.fastq.gz / _R2_001.fastq.gz   (Illumina BCL2FASTQ default)
  _R1.fastq.gz     / _R2.fastq.gz
  _1.fastq.gz      / _2.fastq.gz
  .R1.fastq.gz     / .R2.fastq.gz
  _R1.fq.gz        / _R2.fq.gz
  _1.fq.gz         / _2.fq.gz

Output TSV columns:
  sample_id, r1_path, r2_path, input_dir, pattern, status, notes

Note: r1_path and r2_path are always written as absolute paths so that
downstream Nextflow channels can resolve them correctly regardless of
which work directory the process runs in.

IMPORTANT — sample_id derivation:
  The sample_id is derived by stripping ONLY the R1 read suffix from the
  filename (e.g. _R1.fastq.gz).  Lane suffixes (_L001, _L2, etc.) and
  sample index suffixes (_S1) are intentionally preserved so that the
  derived sample_id matches the sample_id column in the association CSV
  exactly.  Do NOT strip lane suffixes here — the association CSV
  validator handles normalisation if needed.
"""

import argparse
import os
import re
import sys
from pathlib import Path
from collections import defaultdict


# ── Supported R1/R2 pattern pairs ────────────────────────────────────────────
PATTERNS = [
    # (r1_suffix_regex, r2_suffix_regex, pattern_name, r1_replace, r2_replace)
    (r'_R1_001\.(fastq|fq)(\.gz)?$', r'_R2_001\.(fastq|fq)(\.gz)?$',
     '_R1_001/_R2_001', '_R1_001', '_R2_001'),
    (r'_R1\.(fastq|fq)(\.gz)?$',     r'_R2\.(fastq|fq)(\.gz)?$',
     '_R1/_R2',         '_R1',     '_R2'),
    (r'_1\.(fastq|fq)(\.gz)?$',      r'_2\.(fastq|fq)(\.gz)?$',
     '_1/_2',           '_1',      '_2'),
    (r'\.R1\.(fastq|fq)(\.gz)?$',    r'\.R2\.(fastq|fq)(\.gz)?$',
     '.R1/.R2',         '.R1',     '.R2'),
]


def find_fastq_files(input_dir: Path) -> list:
    """
    Recursively find all FASTQ files in input_dir and all subdirectories.
    Uses Path.rglob() so files nested at any depth are discovered.
    """
    fastq_files = []
    for ext in ['*.fastq.gz', '*.fq.gz', '*.fastq', '*.fq']:
        fastq_files.extend(input_dir.rglob(ext))
    fastq_files = sorted(set(fastq_files))

    # Report which subdirectories were found (helps users debug missing samples)
    subdirs = sorted(set(f.parent for f in fastq_files))
    if len(subdirs) == 1:
        print(f"  Searching: {subdirs[0]} (top-level only)", file=sys.stderr)
    else:
        print(f"  Searching {len(subdirs)} subdirectories:", file=sys.stderr)
        for d in subdirs:
            n = sum(1 for f in fastq_files if f.parent == d)
            rel = d.relative_to(input_dir) if d != input_dir else Path('.')
            print(f"    {rel}/  ({n} file(s))", file=sys.stderr)

    return fastq_files


def detect_pattern(filename: str) -> tuple:
    """
    Detect which R1 pattern a filename matches.
    Returns (pattern_name, r1_regex, r2_regex, r1_suf, r2_suf)
    or (None, None, None, None, None) if no R1 pattern matches.
    """
    for r1_regex, r2_regex, pattern_name, r1_suf, r2_suf in PATTERNS:
        if re.search(r1_regex, filename, re.IGNORECASE):
            return pattern_name, r1_regex, r2_regex, r1_suf, r2_suf
    return None, None, None, None, None


def derive_sample_id(filename: str, r1_suffix_regex: str) -> str:
    """
    Derive a sample ID from a filename by stripping only the R1 read suffix.

    Design decision: lane suffixes (_L001, _L2, etc.) and sample index
    suffixes (_S1) are intentionally PRESERVED so that the derived
    sample_id matches the sample_id column in the association CSV exactly.
    For example:
      A1_HEK293A_H3K27me3_10486710_23GVWTLT3_L2_R1.fastq.gz
      → A1_HEK293A_H3K27me3_10486710_23GVWTLT3_L2

    Do NOT strip lane suffixes here — the association CSV validator
    handles normalisation as a fallback if needed.
    """
    sample = re.sub(r1_suffix_regex, '', filename, flags=re.IGNORECASE)
    # Remove any trailing punctuation left by the suffix removal
    sample = sample.rstrip('._-')
    return sample


def pair_fastqs(input_dir: Path, pattern: str = 'auto') -> list:
    """
    Main pairing logic. Returns list of dicts with pairing results.
    """
    all_files = find_fastq_files(input_dir)
    if not all_files:
        print(f"ERROR: No FASTQ files found in {input_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(all_files)} FASTQ files in {input_dir}", file=sys.stderr)

    # Separate R1 files from R2/unclassified files
    r1_files  = []
    unmatched = []

    for f in all_files:
        fname = f.name
        pat_name, r1_regex, r2_regex, r1_suf, r2_suf = detect_pattern(fname)
        if pat_name:
            r1_files.append((f, pat_name, r1_regex, r2_regex, r1_suf, r2_suf))
        else:
            # Check if it's an R2 file (expected — don't warn about these)
            is_r2 = any(
                re.search(r2_pat, fname, re.IGNORECASE)
                for r2_pat in [
                    r'_R2_001\.(fastq|fq)(\.gz)?$',
                    r'_R2\.(fastq|fq)(\.gz)?$',
                    r'_2\.(fastq|fq)(\.gz)?$',
                    r'\.R2\.(fastq|fq)(\.gz)?$',
                ]
            )
            if not is_r2:
                unmatched.append(f)

    results = []

    for r1_path, pat_name, r1_regex, r2_regex, r1_suf, r2_suf in r1_files:
        sample_id = derive_sample_id(r1_path.name, r1_regex)

        # Construct expected R2 filename by simple suffix replacement
        r2_name = r1_path.name
        for old, new in [('_R1_001', '_R2_001'), ('_R1', '_R2'),
                         ('_1.', '_2.'), ('.R1.', '.R2.')]:
            if old in r2_name:
                r2_name = r2_name.replace(old, new, 1)
                break

        r2_path = r1_path.parent / r2_name

        if r2_path.exists():
            status = 'OK'
            notes  = ''
        else:
            # Fallback: search for any R2 in the same directory whose
            # derived sample_id matches
            r2_candidates = [
                f for f in r1_path.parent.iterdir()
                if re.search(r2_regex, f.name, re.IGNORECASE)
                and derive_sample_id(f.name, r2_regex) == sample_id
            ]
            if len(r2_candidates) == 1:
                r2_path = r2_candidates[0]
                status  = 'OK'
                notes   = 'R2 found by sample ID matching'
            elif len(r2_candidates) > 1:
                r2_path = r2_candidates[0]
                status  = 'WARNING'
                notes   = (f'Multiple R2 candidates: '
                           f'{[str(c) for c in r2_candidates]}')
            else:
                r2_path = None
                status  = 'ERROR'
                notes   = f'No matching R2 found for R1: {r1_path.name}'

        results.append({
            'sample_id': sample_id,
            'r1_path':   str(r1_path.resolve()),
            'r2_path':   str(r2_path.resolve()) if r2_path else 'MISSING',
            'input_dir': str(input_dir.resolve()),
            'pattern':   pat_name,
            'status':    status,
            'notes':     notes,
        })

    # Warn about files that couldn't be classified
    for f in unmatched:
        print(f"WARNING: Could not classify file (not R1 or R2): {f}",
              file=sys.stderr)

    return results


def write_tsv(results: list, output_path: str):
    """Write pairing results to TSV."""
    header = ['sample_id', 'r1_path', 'r2_path', 'input_dir',
              'pattern', 'status', 'notes']
    with open(output_path, 'w') as fh:
        fh.write('\t'.join(header) + '\n')
        for row in results:
            fh.write('\t'.join(str(row.get(col, '')) for col in header) + '\n')


def write_report(results: list, report_path: str):
    """Write a human-readable pairing report."""
    n_ok   = sum(1 for r in results if r['status'] == 'OK')
    n_warn = sum(1 for r in results if r['status'] == 'WARNING')
    n_err  = sum(1 for r in results if r['status'] == 'ERROR')

    subdirs = sorted(set(
        str(Path(r['r1_path']).parent)
        for r in results if r['r1_path'] != 'MISSING'
    ))

    with open(report_path, 'w') as fh:
        fh.write("=" * 60 + "\n")
        fh.write("FASTQ Pairing Report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Total samples detected : {len(results)}\n")
        fh.write(f"  OK                   : {n_ok}\n")
        fh.write(f"  Warnings             : {n_warn}\n")
        fh.write(f"  Errors               : {n_err}\n\n")

        fh.write(f"Subdirectories searched ({len(subdirs)}):\n")
        for d in subdirs:
            n = sum(1 for r in results
                    if str(Path(r['r1_path']).parent) == d)
            fh.write(f"  {d}/  ({n} sample(s))\n")
        fh.write("\n")

        if n_err > 0:
            fh.write("ERRORS (pipeline will fail for these samples):\n")
            for r in results:
                if r['status'] == 'ERROR':
                    fh.write(f"  {r['sample_id']}: {r['notes']}\n")
            fh.write("\n")

        if n_warn > 0:
            fh.write("WARNINGS:\n")
            for r in results:
                if r['status'] == 'WARNING':
                    fh.write(f"  {r['sample_id']}: {r['notes']}\n")
            fh.write("\n")

        fh.write("Sample pairs:\n")
        for r in results:
            fh.write(f"  [{r['status']:7s}] {r['sample_id']}\n")
            fh.write(f"           R1: {r['r1_path']}\n")
            fh.write(f"           R2: {r['r2_path']}\n")
            if r['notes']:
                fh.write(f"           Note: {r['notes']}\n")
            fh.write("\n")

    if n_err > 0:
        print(
            f"\nERROR: {n_err} sample(s) could not be paired. "
            "Check pairing_report.txt for details.",
            file=sys.stderr
        )
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Detect and validate paired-end FASTQ files.'
    )
    parser.add_argument('--input_dir', required=True,
                        help='Directory containing FASTQ files')
    parser.add_argument('--pattern',   default='auto',
                        help='Pairing pattern: auto|_R1/_R2|_1/_2|.R1/.R2')
    parser.add_argument('--output',    required=True,
                        help='Output TSV file path')
    parser.add_argument('--report',    required=True,
                        help='Output report text file path')
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"ERROR: Input directory does not exist: {input_dir}",
              file=sys.stderr)
        sys.exit(1)

    results = pair_fastqs(input_dir, args.pattern)
    write_tsv(results, args.output)
    write_report(results, args.report)

    print(f"Pairing complete. {len(results)} sample(s) detected.",
          file=sys.stderr)
    print(f"Output TSV: {args.output}", file=sys.stderr)
    print(f"Report:     {args.report}", file=sys.stderr)


if __name__ == '__main__':
    main()
