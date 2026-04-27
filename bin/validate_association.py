#!/usr/bin/env python3
"""
validate_association.py
=======================
Validate the association CSV against the detected FASTQ pairs TSV.

Validation rules enforced:
  1. Required columns are present.
  2. No duplicate sample_id rows.
  3. Every sample_id in the FASTQ pairs TSV exists in the association CSV.
  4. No extra sample_id rows in association CSV (unless --allow_extra).
  5. is_control must be 'true' or 'false'.
  6. peak_calling_mode must be 'narrow', 'broad', or 'auto'.
  7. All rows in the same merge_group_id share the same species, genome,
     antibody, and condition.
  8. Every non-control sample must have a valid control_group_id that
     resolves to at least one control sample (unless --allow_no_control).
  9. Every control_group_id referenced by a treatment sample must exist
     as a group_id in a control sample row.
 10. replicate values within a group_id must be unique.
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path


REQUIRED_COLUMNS = [
    'sample_id', 'species', 'genome', 'antibody', 'condition',
    'replicate', 'group_id', 'is_control', 'control_group_id',
    'merge_group_id', 'peak_calling_mode',
]

OPTIONAL_COLUMNS = ['notes']

VALID_PEAK_MODES = {'narrow', 'broad', 'auto', ''}
VALID_IS_CONTROL = {'true', 'false'}


def load_csv(path: str) -> list:
    rows = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})
    return rows


def load_tsv(path: str) -> list:
    rows = []
    with open(path, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            rows.append({k.strip(): v.strip() for k, v in row.items()})
    return rows


def validate(assoc_rows: list, pairs_rows: list, allow_no_control: bool,
             allow_extra: bool) -> tuple:
    """
    Returns (errors, warnings) as lists of strings.
    """
    errors   = []
    warnings = []

    # ── 1. Required columns ──────────────────────────────────────────────────
    if not assoc_rows:
        errors.append("Association CSV is empty.")
        return errors, warnings

    actual_cols = set(assoc_rows[0].keys())
    missing_cols = [c for c in REQUIRED_COLUMNS if c not in actual_cols]
    if missing_cols:
        errors.append(
            f"Association CSV is missing required columns: {missing_cols}\n"
            f"  Required: {REQUIRED_COLUMNS}\n"
            f"  Found:    {sorted(actual_cols)}"
        )
        return errors, warnings  # Cannot continue without required columns

    # ── 2. Duplicate sample_id ───────────────────────────────────────────────
    sample_ids = [r['sample_id'] for r in assoc_rows]
    seen = set()
    for sid in sample_ids:
        if sid in seen:
            errors.append(f"Duplicate sample_id in association CSV: '{sid}'")
        seen.add(sid)

    assoc_by_id = {r['sample_id']: r for r in assoc_rows}

    # ── 3. FASTQ pairs cross-check ───────────────────────────────────────────
    fastq_ids = {r['sample_id'] for r in pairs_rows if r.get('status') == 'OK'}
    assoc_ids = set(assoc_by_id.keys())

    missing_in_assoc = fastq_ids - assoc_ids
    for sid in sorted(missing_in_assoc):
        errors.append(
            f"Sample '{sid}' found in FASTQ pairs but missing from association CSV."
        )

    extra_in_assoc = assoc_ids - fastq_ids
    for sid in sorted(extra_in_assoc):
        msg = f"Sample '{sid}' in association CSV has no matching FASTQ pair."
        if allow_extra:
            warnings.append(msg + " (allowed by --allow_extra)")
        else:
            errors.append(msg + " Add --allow_extra to suppress this error.")

    # ── 4. Field value validation ────────────────────────────────────────────
    for i, row in enumerate(assoc_rows, start=2):
        sid = row['sample_id']

        # is_control
        if row['is_control'].lower() not in VALID_IS_CONTROL:
            errors.append(
                f"Row {i} (sample '{sid}'): is_control must be 'true' or 'false', "
                f"got '{row['is_control']}'"
            )

        # peak_calling_mode
        if row['peak_calling_mode'].lower() not in VALID_PEAK_MODES:
            errors.append(
                f"Row {i} (sample '{sid}'): peak_calling_mode must be "
                f"'narrow', 'broad', or 'auto', got '{row['peak_calling_mode']}'"
            )

        # Required non-empty fields
        for col in ['species', 'genome', 'antibody', 'condition', 'replicate',
                    'group_id', 'merge_group_id']:
            if not row.get(col, '').strip():
                errors.append(
                    f"Row {i} (sample '{sid}'): column '{col}' is empty."
                )

    # ── 5. Control group resolution ──────────────────────────────────────────
    control_group_ids = {
        r['group_id']
        for r in assoc_rows
        if r['is_control'].lower() == 'true'
    }

    for row in assoc_rows:
        if row['is_control'].lower() == 'false':
            ctrl_id = row.get('control_group_id', '').strip()
            if not ctrl_id:
                if allow_no_control:
                    warnings.append(
                        f"Sample '{row['sample_id']}' has no control_group_id. "
                        "Peak calling will proceed without control (--allow_no_control)."
                    )
                else:
                    errors.append(
                        f"Sample '{row['sample_id']}' has no control_group_id. "
                        "Set --allow_no_control to allow control-free peak calling."
                    )
            elif ctrl_id not in control_group_ids:
                errors.append(
                    f"Sample '{row['sample_id']}' references control_group_id "
                    f"'{ctrl_id}' which does not match any control sample's group_id. "
                    f"Available control group_ids: {sorted(control_group_ids)}"
                )

    # ── 6. merge_group_id consistency ────────────────────────────────────────
    merge_groups = defaultdict(list)
    for row in assoc_rows:
        merge_groups[row['merge_group_id']].append(row)

    for merge_id, rows in merge_groups.items():
        species_vals   = set(r['species']   for r in rows)
        genome_vals    = set(r['genome']    for r in rows)
        antibody_vals  = set(r['antibody']  for r in rows)
        condition_vals = set(r['condition'] for r in rows)

        if len(species_vals) > 1:
            errors.append(
                f"merge_group_id '{merge_id}' has inconsistent species: {species_vals}"
            )
        if len(genome_vals) > 1:
            errors.append(
                f"merge_group_id '{merge_id}' has inconsistent genome: {genome_vals}"
            )
        if len(antibody_vals) > 1:
            errors.append(
                f"merge_group_id '{merge_id}' has inconsistent antibody: {antibody_vals}"
            )
        if len(condition_vals) > 1:
            warnings.append(
                f"merge_group_id '{merge_id}' spans multiple conditions: {condition_vals}. "
                "Verify this is intentional."
            )

    # ── 7. Replicate uniqueness within group_id ───────────────────────────────
    group_replicates = defaultdict(list)
    for row in assoc_rows:
        group_replicates[row['group_id']].append(row['replicate'])

    for grp_id, reps in group_replicates.items():
        if len(reps) != len(set(reps)):
            dup_reps = [r for r in reps if reps.count(r) > 1]
            errors.append(
                f"group_id '{grp_id}' has duplicate replicate values: {dup_reps}"
            )

    return errors, warnings


def write_report(errors: list, warnings: list, report_path: str,
                 assoc_rows: list, pairs_rows: list):
    with open(report_path, 'w') as fh:
        fh.write("=" * 60 + "\n")
        fh.write("Association CSV Validation Report\n")
        fh.write("=" * 60 + "\n\n")
        fh.write(f"Samples in association CSV : {len(assoc_rows)}\n")
        fh.write(f"Samples in FASTQ pairs TSV : {len(pairs_rows)}\n")
        fh.write(f"Errors                     : {len(errors)}\n")
        fh.write(f"Warnings                   : {len(warnings)}\n\n")

        if errors:
            fh.write("ERRORS:\n")
            for i, e in enumerate(errors, 1):
                fh.write(f"  [{i}] {e}\n")
            fh.write("\n")

        if warnings:
            fh.write("WARNINGS:\n")
            for i, w in enumerate(warnings, 1):
                fh.write(f"  [{i}] {w}\n")
            fh.write("\n")

        if not errors and not warnings:
            fh.write("All validation checks passed.\n")

        # Summary table
        fh.write("\nSample Summary:\n")
        fh.write(f"{'sample_id':<30} {'group_id':<25} {'is_control':<12} "
                 f"{'merge_group_id':<30} {'control_group_id':<25}\n")
        fh.write("-" * 122 + "\n")
        for row in assoc_rows:
            fh.write(
                f"{row['sample_id']:<30} {row['group_id']:<25} "
                f"{row['is_control']:<12} {row['merge_group_id']:<30} "
                f"{row.get('control_group_id',''):<25}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description='Validate association CSV against FASTQ pairs TSV.'
    )
    parser.add_argument('--association_csv', required=True)
    parser.add_argument('--pairs_tsv',       required=True)
    parser.add_argument('--output',          required=True,
                        help='Validated association CSV output path')
    parser.add_argument('--report',          required=True,
                        help='Validation report text file path')
    parser.add_argument('--allow_no_control', action='store_true',
                        help='Allow peak calling without control samples')
    parser.add_argument('--allow_extra',      action='store_true',
                        help='Allow extra sample_ids in association CSV')
    args = parser.parse_args()

    assoc_rows = load_csv(args.association_csv)
    pairs_rows = load_tsv(args.pairs_tsv)

    errors, warnings = validate(
        assoc_rows, pairs_rows,
        allow_no_control=args.allow_no_control,
        allow_extra=args.allow_extra
    )

    write_report(errors, warnings, args.report, assoc_rows, pairs_rows)

    # Print summary to stderr
    for w in warnings:
        print(f"WARNING: {w}", file=sys.stderr)
    for e in errors:
        print(f"ERROR: {e}", file=sys.stderr)

    if errors:
        print(
            f"\nValidation FAILED with {len(errors)} error(s). "
            "See validation_report.txt for details.",
            file=sys.stderr
        )
        sys.exit(1)

    # Write validated CSV (pass-through with normalized fields)
    with open(args.output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=REQUIRED_COLUMNS + OPTIONAL_COLUMNS,
                                extrasaction='ignore')
        writer.writeheader()
        for row in assoc_rows:
            # Normalize is_control and peak_calling_mode
            row['is_control']       = row['is_control'].lower()
            row['peak_calling_mode'] = row['peak_calling_mode'].lower() or 'narrow'
            writer.writerow(row)

    print(
        f"Validation passed ({len(warnings)} warning(s)). "
        f"Validated CSV written to {args.output}",
        file=sys.stderr
    )


if __name__ == '__main__':
    main()
