#!/usr/bin/env python3
"""
RAPID — scripts/track_threshold_scan.py
Explore how many reference-specific genes you WOULD keep for different filter
thresholds, using featureCounts tables already on disk (no need to re-run).

For a grid of --other-max-count values (and optionally a few --ref-min-count
values), reports the number of genes that pass:
    ref count   >= ref_min      AND
    every other <= other_max

Usage:
    python3 track_threshold_scan.py \
        --ref counts_ref.txt \
        --others counts_other_1.txt counts_other_2.txt ... \
        [--ref-min 10] [--other-max-grid 0 1 2 5 10 20 50 100]
        [--ref-min-grid 5 10 20]     # optional: scan ref-min too

Prints a table to stdout.
"""

import argparse
import sys
from pathlib import Path


def read_counts(path):
    """Return {geneid: count} from a featureCounts table (count = last column)."""
    counts = {}
    header_seen = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True          # skip the column header row
                continue
            row = line.rstrip("\n")
            if not row:
                continue
            cols = row.split("\t")
            try:
                counts[cols[0]] = int(round(float(cols[-1])))
            except ValueError:
                counts[cols[0]] = 0
    return counts


def count_kept(ref, others, ref_min, other_max):
    """Number of genes passing the filter, plus how many are multi-exonic."""
    kept = 0
    for gene, rc in ref.items():
        if rc < ref_min:
            continue
        if all(o.get(gene, 0) <= other_max for o in others):
            kept += 1
    return kept


def main():
    p = argparse.ArgumentParser(description="Scan track filter thresholds.")
    p.add_argument("--ref", required=True)
    p.add_argument("--others", nargs="+", required=True)
    p.add_argument("--ref-min", type=int, default=10, dest="ref_min")
    p.add_argument("--other-max-grid", nargs="+", type=int,
                   default=[0, 1, 2, 5, 10, 20, 50, 100, 200, 500],
                   dest="other_grid")
    p.add_argument("--ref-min-grid", nargs="+", type=int, default=None,
                   dest="ref_grid",
                   help="Optional: also scan several ref-min values.")
    args = p.parse_args()

    ref = read_counts(args.ref)
    others = [read_counts(o) for o in args.others]

    n_ref_expr = sum(1 for v in ref.values() if v >= args.ref_min)
    print(f"# reference genes with count >= {args.ref_min}: {n_ref_expr}", file=sys.stderr)
    print(f"# other samples: {len(others)}\n", file=sys.stderr)

    ref_mins = args.ref_grid if args.ref_grid else [args.ref_min]

    if len(ref_mins) == 1:
        rm = ref_mins[0]
        print(f"{'other_max':>10} | {'genes_kept':>10}   (ref_min = {rm})")
        print("-" * 34)
        for om in args.other_grid:
            k = count_kept(ref, others, rm, om)
            print(f"{om:>10} | {k:>10}")
    else:
        # matrix: rows = other_max, cols = ref_min
        header = f"{'other_max':>10} | " + " ".join(f"rm={rm:<6}" for rm in ref_mins)
        print(header)
        print("-" * len(header))
        for om in args.other_grid:
            cells = " ".join(f"{count_kept(ref, others, rm, om):<9}" for rm in ref_mins)
            print(f"{om:>10} | {cells}")

    print("\nInterpretation: pick the smallest other_max that yields a workable "
          "number of candidate genes (e.g. a few dozen to a few hundred). "
          "Higher other_max = more permissive = keeps genes with some leaky "
          "expression in other stages.", file=sys.stderr)


if __name__ == "__main__":
    main()
