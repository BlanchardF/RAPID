#!/usr/bin/env python3
"""
RAPID — scripts/track_cpm_scan.py
Explore a CPM + fold-enrichment filter for `track` mode, using featureCounts
tables already on disk (no need to re-run the pipeline).

For a grid of fold-change values, reports how many genes would be kept with:
    CPM_ref >= --min-cpm-ref                      (expressed in reference), AND
    CPM_ref >= fold * CPM_other for EVERY other   (enriched vs each other stage)

CPM = count / (library size / 1e6), where library size = sum of assigned counts
in that sample. This normalises for the different sequencing depths of your
samples, so the comparison is fair (unlike a raw-count threshold).

Optionally traces specific genes (--gene) to show their per-sample CPM and the
exact fold-change threshold below which they are kept.

Usage:
    python3 track_cpm_scan.py \
        --ref counts_ref.txt \
        --others counts_other_1.txt counts_other_2.txt ... \
        [--min-cpm-ref 1] [--fold-grid 2 3 5 10 20 50] \
        [--gene Smp_169190 --gene Smp_032670]
"""

import argparse
import sys


def read_counts(path):
    """Return {geneid: count} from a featureCounts table (last column = count)."""
    counts = {}
    header_seen = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if not header_seen:
                header_seen = True          # skip column header row
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


def to_cpm(counts):
    """Return ({geneid: cpm}, library_size)."""
    lib = sum(counts.values())
    if lib == 0:
        return {g: 0.0 for g in counts}, 0
    factor = 1e6 / lib
    return {g: v * factor for g, v in counts.items()}, lib


def main():
    p = argparse.ArgumentParser(description="Scan CPM + fold-enrichment thresholds (track mode).")
    p.add_argument("--ref", required=True)
    p.add_argument("--others", nargs="+", required=True)
    p.add_argument("--min-cpm-ref", type=float, default=1.0, dest="min_cpm",
                   help="Minimum CPM in the reference to consider a gene expressed. Default: 1.")
    p.add_argument("--fold-grid", nargs="+", type=float,
                   default=[2, 3, 5, 10, 20, 50], dest="fold_grid",
                   help="Fold-change values to scan. Default: 2 3 5 10 20 50.")
    p.add_argument("--gene", action="append", default=[],
                   help="Trace a specific gene (repeatable). Matches with or without the "
                        "'gene:' prefix.")
    args = p.parse_args()

    ref_raw = read_counts(args.ref)
    others_raw = [read_counts(o) for o in args.others]
    ref_cpm, ref_lib = to_cpm(ref_raw)
    others_cpm, libs = [], []
    for o in others_raw:
        cpm, lib = to_cpm(o)
        others_cpm.append(cpm)
        libs.append(lib)

    print("# library sizes (assigned reads):", file=sys.stderr)
    print(f"#   ref     : {ref_lib:,}", file=sys.stderr)
    for i, lib in enumerate(libs, 1):
        print(f"#   other_{i} : {lib:,}", file=sys.stderr)
    print(file=sys.stderr)

    def kept_count(fold):
        n = 0
        for g, rc in ref_cpm.items():
            if rc < args.min_cpm:
                continue
            if all(rc >= fold * oc.get(g, 0.0) for oc in others_cpm):
                n += 1
        return n

    print(f"{'fold_change':>12} | {'genes_kept':>10}   (min_cpm_ref = {args.min_cpm})")
    print("-" * 42)
    for fold in args.fold_grid:
        print(f"{fold:>12g} | {kept_count(fold):>10}")

    # ── Gene tracing ──────────────────────────────────────────────────────────
    def find_gene(q):
        cands = [q, f"gene:{q}"]
        for g in ref_cpm:
            if g in cands or g.split(":")[-1] == q:
                return g
        for g in ref_cpm:        # fallback: substring
            if q in g:
                return g
        return None

    for q in args.gene:
        g = find_gene(q)
        print(f"\n=== {q} ===")
        if g is None:
            print("  not found in the reference table.")
            continue
        print(f"  full id             : {g}")
        print(f"  CPM ref             : {ref_cpm[g]:10.2f}   (raw {ref_raw.get(g, 0)})")
        max_other, max_i = 0.0, None
        for i, oc in enumerate(others_cpm, 1):
            v = oc.get(g, 0.0)
            flag = "  <== max" if v > max_other else ""
            if v > max_other:
                max_other, max_i = v, i
            print(f"  CPM other_{i}        : {v:10.2f}   (raw {others_raw[i-1].get(g, 0)}){flag}")
        passes_min = ref_cpm[g] >= args.min_cpm
        if max_other > 0:
            enr = ref_cpm[g] / max_other
            print(f"  min enrichment      : x{enr:.1f}  (vs other_{max_i}, the most expressed of the others)")
            print(f"  -> kept if --fold-change <= {enr:.1f}  AND  CPM ref >= {args.min_cpm}"
                  f"  ({'OK' if passes_min else 'CPM ref TOO LOW'})")
        else:
            print(f"  min enrichment      : infinite (absent from all others)")
            print(f"  -> kept for any --fold-change  ({'OK' if passes_min else 'CPM ref TOO LOW'})")


if __name__ == "__main__":
    main()
