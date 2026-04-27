#!/usr/bin/env python3
"""
make_test_data.py
=================
Generate minimal synthetic test data for pipeline CI testing.

Creates:
  test_data/fastq/          - Tiny synthetic paired FASTQ files
  test_data/test_association.csv - Minimal association CSV

NOTE: This generates synthetic reads for pipeline structure testing only.
For biological validation, use real CUT&Tag data from a published dataset
(e.g., GEO GSE145187 - Kaya-Okur et al. 2019 CUT&Tag data).

Usage:
  python bin/make_test_data.py --outdir test_data --n_reads 1000
"""

import argparse
import gzip
import os
import random
import string
from pathlib import Path


BASES = 'ACGT'
QUAL  = 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'


def random_seq(length: int) -> str:
    return ''.join(random.choices(BASES, k=length))


def write_fastq_pair(r1_path: str, r2_path: str, n_reads: int = 1000,
                     read_length: int = 75):
    """Write a minimal paired FASTQ file pair."""
    with gzip.open(r1_path, 'wt') as r1, gzip.open(r2_path, 'wt') as r2:
        for i in range(n_reads):
            read_id = f"@SIM_READ_{i:06d}"
            seq1    = random_seq(read_length)
            seq2    = random_seq(read_length)
            qual    = QUAL[:read_length]

            r1.write(f"{read_id}/1\n{seq1}\n+\n{qual}\n")
            r2.write(f"{read_id}/2\n{seq2}\n+\n{qual}\n")


def write_association_csv(output_path: str, samples: list):
    """Write a minimal association CSV for testing."""
    header = (
        'sample_id,species,genome,antibody,condition,replicate,'
        'group_id,is_control,control_group_id,merge_group_id,'
        'peak_calling_mode,notes'
    )
    rows = []
    for s in samples:
        rows.append(','.join([
            s['sample_id'],
            s.get('species', 'human'),
            s.get('genome', 'hg38_chr21'),
            s.get('antibody', 'H3K27ac'),
            s.get('condition', 'treated'),
            s.get('replicate', '1'),
            s.get('group_id', s['sample_id']),
            s.get('is_control', 'false'),
            s.get('control_group_id', 'IgG_group'),
            s.get('merge_group_id', 'human_H3K27ac_treated'),
            s.get('peak_calling_mode', 'narrow'),
            s.get('notes', 'test sample'),
        ]))

    with open(output_path, 'w') as fh:
        fh.write(header + '\n')
        for row in rows:
            fh.write(row + '\n')


def main():
    parser = argparse.ArgumentParser(
        description='Generate minimal synthetic test data for pipeline testing.'
    )
    parser.add_argument('--outdir',   default='test_data',
                        help='Output directory for test data')
    parser.add_argument('--n_reads',  type=int, default=1000,
                        help='Number of synthetic reads per sample')
    parser.add_argument('--n_samples', type=int, default=4,
                        help='Number of test samples (including 1 control)')
    args = parser.parse_args()

    outdir = Path(args.outdir)
    fastq_dir = outdir / 'fastq'
    fastq_dir.mkdir(parents=True, exist_ok=True)

    # Define test samples: 3 treatment + 1 IgG control
    samples = []

    # Treatment samples
    for i in range(1, args.n_samples):
        sid = f"H3K27ac_rep{i}"
        samples.append({
            'sample_id':       sid,
            'species':         'human',
            'genome':          'hg38_chr21',
            'antibody':        'H3K27ac',
            'condition':       'treated',
            'replicate':       str(i),
            'group_id':        f'H3K27ac_treated_rep{i}',
            'is_control':      'false',
            'control_group_id': 'IgG_group',
            'merge_group_id':  'human_H3K27ac_treated',
            'peak_calling_mode': 'narrow',
            'notes':           f'Test H3K27ac replicate {i}',
        })
        r1 = str(fastq_dir / f"{sid}_R1.fastq.gz")
        r2 = str(fastq_dir / f"{sid}_R2.fastq.gz")
        print(f"Writing {sid}: {r1}, {r2}")
        write_fastq_pair(r1, r2, n_reads=args.n_reads)

    # IgG control
    ctrl_sid = 'IgG_ctrl'
    samples.append({
        'sample_id':       ctrl_sid,
        'species':         'human',
        'genome':          'hg38_chr21',
        'antibody':        'IgG',
        'condition':       'treated',
        'replicate':       '1',
        'group_id':        'IgG_group',
        'is_control':      'true',
        'control_group_id': '',
        'merge_group_id':  'human_IgG_ctrl',
        'peak_calling_mode': 'narrow',
        'notes':           'IgG negative control',
    })
    r1 = str(fastq_dir / f"{ctrl_sid}_R1.fastq.gz")
    r2 = str(fastq_dir / f"{ctrl_sid}_R2.fastq.gz")
    print(f"Writing {ctrl_sid}: {r1}, {r2}")
    write_fastq_pair(r1, r2, n_reads=args.n_reads)

    # Write association CSV
    assoc_path = str(outdir / 'test_association.csv')
    write_association_csv(assoc_path, samples)
    print(f"Association CSV written to: {assoc_path}")

    print(f"\nTest data generated in: {outdir}")
    print(f"  Samples: {len(samples)} ({len(samples)-1} treatment + 1 control)")
    print(f"  Reads per sample: {args.n_reads}")
    print("\nNOTE: These are random synthetic reads and will NOT produce")
    print("meaningful alignments or peaks. Use real data for biological testing.")
    print("\nFor real test data, download from:")
    print("  GEO: GSE145187 (Kaya-Okur et al. 2019 CUT&Tag)")
    print("  GEO: GSE124557 (Henikoff et al. 2019 CUT&RUN)")


if __name__ == '__main__':
    main()
