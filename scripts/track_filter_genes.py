#!/usr/bin/env python3
"""
RAPID — scripts/track_filter_genes.py
Keep genes specific to the reference RNA-seq sample (track mode).

Two filtering modes:

  absolute  (default, backward-compatible):
      keep gene iff  ref count >= --ref-min  AND  every other count <= --other-max
      Works on raw counts. Simple, but depth-dependent and blind to enrichment.

  cpm:
      keep gene iff  CPM_ref >= --min-cpm-ref
                     AND  CPM_ref >= --fold-change * CPM_other  for EVERY other
      CPM = count / (library size / 1e6), library size = sum of assigned counts
      in that sample. This normalises for sequencing depth and keeps genes that
      are ENRICHED in the reference (not merely present/absent), which is the
      right notion for stage-specific markers.

Output (both modes): the reference featureCounts table subset to the kept
genes, in the exact same format, so the downstream R step reads it unchanged.

Usage:
  # absolute (raw counts)
  python3 track_filter_genes.py --mode absolute \
      --ref counts_ref.txt --others counts_a.txt counts_b.txt \
      --ref-min 10 --other-max 0 --out counts_ref_specific.txt

  # cpm (depth-normalised enrichment)
  python3 track_filter_genes.py --mode cpm \
      --ref counts_ref.txt --others counts_a.txt counts_b.txt \
      --min-cpm-ref 1 --fold-change 10 --out counts_ref_specific.txt
"""

import argparse
import sys


def read_featurecounts(path):
    """
    Parse a featureCounts file.
    Returns (comment_lines, header_line, data_lines, counts_dict, lib_size)
    where counts_dict maps Geneid -> int count (last column) and lib_size is the
    sum of all counts (used for CPM).
    """
    comments = []
    header   = None
    data     = []
    counts   = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
                continue
            if header is None:
                header = line.rstrip("\n")
                continue
            row = line.rstrip("\n")
            if not row:
                continue
            cols = row.split("\t")
            gene = cols[0]
            try:
                c = int(round(float(cols[-1])))
            except ValueError:
                c = 0
            counts[gene] = c
            data.append(row)
    lib_size = sum(counts.values())
    return comments, header, data, counts, lib_size


def cpm_dict(counts, lib_size):
    if lib_size == 0:
        return {g: 0.0 for g in counts}
    factor = 1e6 / lib_size
    return {g: v * factor for g, v in counts.items()}


def main():
    p = argparse.ArgumentParser(description="Keep reference-specific genes (RAPID track mode).")
    p.add_argument("--ref", required=True, help="featureCounts table of the reference sample.")
    p.add_argument("--others", nargs="+", required=True,
                   help="featureCounts tables of the other samples (>= 1).")
    p.add_argument("--out", required=True,
                   help="Output featureCounts table restricted to reference-specific genes.")
    p.add_argument("--mode", choices=["absolute", "cpm"], default="absolute",
                   help="Filtering mode. Default: absolute.")
    # absolute-mode parameters
    p.add_argument("--ref-min", type=int, default=10, dest="ref_min",
                   help="[absolute] Min raw count in the reference. Default: 10.")
    p.add_argument("--other-max", type=int, default=0, dest="other_max",
                   help="[absolute] Max raw count allowed in an 'other' sample. Default: 0.")
    # cpm-mode parameters
    p.add_argument("--min-cpm-ref", type=float, default=1.0, dest="min_cpm_ref",
                   help="[cpm] Min CPM in the reference. Default: 1.0.")
    p.add_argument("--fold-change", type=float, default=10.0, dest="fold_change",
                   help="[cpm] Required fold-enrichment of ref CPM over EACH other. Default: 10.")
    args = p.parse_args()

    comments, header, data, ref_counts, ref_lib = read_featurecounts(args.ref)
    if header is None:
        print(f"[track_filter] ERROR: no header found in {args.ref}", file=sys.stderr)
        sys.exit(1)

    others = [read_featurecounts(o) for o in args.others]
    other_counts = [o[3] for o in others]
    other_libs   = [o[4] for o in others]

    kept = set()

    if args.mode == "absolute":
        n_ref_expressed = sum(1 for v in ref_counts.values() if v >= args.ref_min)
        for gene, rc in ref_counts.items():
            if rc < args.ref_min:
                continue
            if all(oc.get(gene, 0) <= args.other_max for oc in other_counts):
                kept.add(gene)
        criterion = (f"ref count >= {args.ref_min} AND every other <= {args.other_max} "
                     f"(raw counts)")
    else:  # cpm
        ref_cpm = cpm_dict(ref_counts, ref_lib)
        other_cpm = [cpm_dict(c, l) for c, l in zip(other_counts, other_libs)]
        n_ref_expressed = sum(1 for v in ref_cpm.values() if v >= args.min_cpm_ref)
        for gene, rcpm in ref_cpm.items():
            if rcpm < args.min_cpm_ref:
                continue
            if all(rcpm >= args.fold_change * oc.get(gene, 0.0) for oc in other_cpm):
                kept.add(gene)
        criterion = (f"CPM_ref >= {args.min_cpm_ref} AND CPM_ref >= "
                     f"{args.fold_change}x every other CPM (depth-normalised)")

    with open(args.out, "w") as out:
        for c in comments:
            out.write(c + "\n")
        out.write(header + "\n")
        for row in data:
            gene = row.split("\t", 1)[0]
            if gene in kept:
                out.write(row + "\n")

    print(f"[track_filter] mode                 : {args.mode}", file=sys.stderr)
    print(f"[track_filter] criterion            : {criterion}", file=sys.stderr)
    print(f"[track_filter] reference sample     : {args.ref}", file=sys.stderr)
    print(f"[track_filter] other samples        : {len(other_counts)}", file=sys.stderr)
    if args.mode == "cpm":
        print(f"[track_filter] library sizes (ref/others): {ref_lib:,} / "
              f"{', '.join(f'{l:,}' for l in other_libs)}", file=sys.stderr)
    print(f"[track_filter] genes expressed in ref: {n_ref_expressed}", file=sys.stderr)
    print(f"[track_filter] reference-specific kept: {len(kept)} -> {args.out}", file=sys.stderr)

    if not kept:
        hint = ("lower --fold-change or --min-cpm-ref" if args.mode == "cpm"
                else "raise --other-max or lower --ref-min")
        print(f"[track_filter] WARNING: no gene passed the filter. Consider: {hint}.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
