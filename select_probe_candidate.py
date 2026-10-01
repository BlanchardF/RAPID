#!/usr/bin/env python3
"""
RAPID — scripts/select_probe_candidate.py
Only used when `rapid auto --probe` is active.

Primer3 was asked (via PRIMER_NUM_RETURN in common.smk) to return several
candidate primer PAIRS per gene/junction, each pair carrying its own
internal oligo (probe), picked by Primer3 to fit the size/Tm window given
in common.smk (18-30 nt, Tm ~ primers +8/+10°C).

For every SEQUENCE_ID block, this script walks the candidate pairs in
Primer3's own rank order (0 = best penalty) and keeps the first one whose
probe does NOT start with a 5' G. If none qualifies, it falls back to the
best pair (index 0) even if its probe starts with a G.

The chosen candidate is rewritten as index "0" so that the rest of the
pipeline (primer3_top_n's grep/awk, primer3_summary.py) keeps working
unchanged, exactly as if Primer3 had picked it directly.

Usage:
    python3 select_probe_candidate.py <primer3_raw_results.txt> <output.txt>
"""

import re
import sys

PAIR_INDEX_RE = re.compile(r"^PRIMER_PAIR_(\d+)_PENALTY$")

# Fields copied verbatim from the chosen candidate index, renamed to "_0".
LEFT_RIGHT_FIELDS = [
    "{side}_{i}",
    "{side}_{i}_SEQUENCE",
    "{side}_{i}_TM",
    "{side}_{i}_GC_PERCENT",
    "{side}_{i}_SELF_ANY_TH",
    "{side}_{i}_SELF_END_TH",
    "{side}_{i}_HAIRPIN_TH",
    "{side}_{i}_END_STABILITY",
]
PAIR_FIELDS = [
    "PRIMER_PAIR_{i}_PENALTY",
    "PRIMER_PAIR_{i}_COMPL_ANY_TH",
    "PRIMER_PAIR_{i}_COMPL_END_TH",
    "PRIMER_PAIR_{i}_PRODUCT_SIZE",
    "PRIMER_PAIR_{i}_PRODUCT_TM",
]
INTERNAL_FIELDS = [
    "PRIMER_INTERNAL_{i}",
    "PRIMER_INTERNAL_{i}_SEQUENCE",
    "PRIMER_INTERNAL_{i}_TM",
    "PRIMER_INTERNAL_{i}_GC_PERCENT",
    "PRIMER_INTERNAL_{i}_SELF_ANY_TH",
    "PRIMER_INTERNAL_{i}_SELF_END_TH",
    "PRIMER_INTERNAL_{i}_HAIRPIN_TH",
]
# Fields kept as-is, independent of candidate index.
PASSTHROUGH_KEYS = ["SEQUENCE_ID", "SEQUENCE_TEMPLATE", "SEQUENCE_OVERLAP_JUNCTION_LIST"]


def parse_blocks(path):
    """Yield (dict, raw_lines) for every '='-terminated block, in order."""
    current = {}
    raw_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "=":
                yield current, raw_lines
                current, raw_lines = {}, []
            else:
                raw_lines.append(line)
                if "=" in line:
                    key, _, val = line.partition("=")
                    current[key] = val


def candidate_indices(rec):
    """All candidate pair indices present in this block, in Primer3 rank order."""
    idxs = set()
    for key in rec:
        m = PAIR_INDEX_RE.match(key)
        if m:
            idxs.add(int(m.group(1)))
    return sorted(idxs)


def choose_index(rec, idxs):
    """First candidate whose probe doesn't start with a 5' G; else idx 0."""
    for i in idxs:
        seq = rec.get(f"PRIMER_INTERNAL_{i}_SEQUENCE")
        if seq and not seq.upper().startswith("G"):
            return i
    return idxs[0] if idxs else None


def write_block(out_fh, rec, chosen):
    for key in PASSTHROUGH_KEYS:
        if key in rec:
            out_fh.write(f"{key}={rec[key]}\n")

    if chosen is None:
        out_fh.write("=\n")
        return

    for side in ("PRIMER_LEFT", "PRIMER_RIGHT"):
        for tmpl in LEFT_RIGHT_FIELDS:
            src_key = tmpl.format(side=side, i=chosen)
            if src_key in rec:
                dst_key = tmpl.format(side=side, i=0)
                out_fh.write(f"{dst_key}={rec[src_key]}\n")

    for tmpl in PAIR_FIELDS:
        src_key = tmpl.format(i=chosen)
        if src_key in rec:
            dst_key = tmpl.format(i=0)
            out_fh.write(f"{dst_key}={rec[src_key]}\n")

    for tmpl in INTERNAL_FIELDS:
        src_key = tmpl.format(i=chosen)
        if src_key in rec:
            dst_key = tmpl.format(i=0)
            out_fh.write(f"{dst_key}={rec[src_key]}\n")

    out_fh.write("=\n")


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <primer3_raw_results.txt> <output.txt>", file=sys.stderr)
        sys.exit(1)

    in_path, out_path = sys.argv[1], sys.argv[2]

    n_blocks = 0
    n_probe_ok = 0
    n_probe_fallback_g = 0
    n_no_pair = 0

    with open(out_path, "w") as out_fh:
        for rec, raw_lines in parse_blocks(in_path):
            if not rec:
                continue
            n_blocks += 1
            idxs = candidate_indices(rec)

            if not idxs:
                # No candidate pair at all for this gene (Primer3 found nothing
                # or failed) — keep the original lines untouched for debugging.
                n_no_pair += 1
                for line in raw_lines:
                    out_fh.write(line + "\n")
                out_fh.write("=\n")
                continue

            chosen = choose_index(rec, idxs)
            probe_seq = rec.get(f"PRIMER_INTERNAL_{chosen}_SEQUENCE")
            if probe_seq:
                if probe_seq.upper().startswith("G"):
                    n_probe_fallback_g += 1
                else:
                    n_probe_ok += 1

            write_block(out_fh, rec, chosen)

    print(
        f"[select_probe_candidate] {n_blocks} gene/junction block(s) processed — "
        f"{n_probe_ok} with a G-free 5' probe, "
        f"{n_probe_fallback_g} fell back to a probe starting with G, "
        f"{n_no_pair} had no candidate pair.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
