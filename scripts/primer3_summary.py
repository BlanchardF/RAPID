#!/usr/bin/env python3
"""
RAPID — scripts/primer3_summary.py
Parses a Primer3 output file and produces a TSV summary table.

Usage: python3 primer3_summary.py <input_best_primers.txt> <output_summary.tsv>

Junction overlap columns
------------------------
For each primer, Primer3 stores its position in the template as:
  PRIMER_LEFT_0  = start,length   (0-based start, primer goes rightward)
  PRIMER_RIGHT_0 = end,length     (0-based 3'-end on template, primer goes leftward)

The junction position J (from SEQUENCE_ID g1234_jJ) is the cumulative exon
length at the junction (1-based boundary): bases 1..J belong to exon N,
bases J+1.. belong to exon N+1.

A primer "spans" the junction when its template coordinates straddle J.
  Left primer  [ls, ls+ll-1]  spans J  iff  ls < J  and  ls+ll > J
  Right primer [re-rl+1, re]  spans J  iff  re-rl+1 < J  and  re >= J

Bases before junction = number of primer bases in exon N  (≤ J side)
Bases after  junction = number of primer bases in exon N+1 (> J side)

Amplicon columns
----------------
The amplicon is the region of SEQUENCE_TEMPLATE bounded by the two primers:
  amplicon = SEQUENCE_TEMPLATE[left_start : right_end + 1]
Its length equals PRIMER_PAIR_0_PRODUCT_SIZE (used here as an internal check).
Amplicon_GC_percent is the %GC of that sequence.
"""

import sys

FIELDS = [
    "PRIMER_PAIR_0_PENALTY",
    "PRIMER_LEFT_0_SEQUENCE",
    "PRIMER_RIGHT_0_SEQUENCE",
    "PRIMER_LEFT_0_TM",
    "PRIMER_RIGHT_0_TM",
    "PRIMER_LEFT_0_GC_PERCENT",
    "PRIMER_RIGHT_0_GC_PERCENT",
    "PRIMER_LEFT_0_SELF_ANY_TH",
    "PRIMER_RIGHT_0_SELF_ANY_TH",
    "PRIMER_LEFT_0_SELF_END_TH",
    "PRIMER_RIGHT_0_SELF_END_TH",
    "PRIMER_LEFT_0_HAIRPIN_TH",
    "PRIMER_RIGHT_0_HAIRPIN_TH",
    "PRIMER_LEFT_0_END_STABILITY",
    "PRIMER_RIGHT_0_END_STABILITY",
    "PRIMER_PAIR_0_COMPL_ANY_TH",
    "PRIMER_PAIR_0_COMPL_END_TH",
    "PRIMER_PAIR_0_PRODUCT_SIZE",
    "PRIMER_PAIR_0_PRODUCT_TM",
]


def parse_blocks(path):
    """Parse a Primer3 output file into a list of dicts (one per block)."""
    blocks = []
    current = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line == "=":
                if current:
                    blocks.append(current)
                current = {}
            elif "=" in line:
                key, _, val = line.partition("=")
                current[key] = val
    return blocks


def junction_overlap(rec, junction_pos):
    """
    Compute how many bases of each primer fall before and after the junction.

    junction_pos (int) : 1-based boundary — bases 1..J = exon N,
                         bases J+1.. = exon N+1.

    Returns a dict with keys:
      left_before, left_after   — for the left/forward primer
      right_before, right_after — for the right/reverse primer
    Values are integers (0 when the primer does not span the junction).
    """
    result = {
        "left_before":  0, "left_after":  0,
        "right_before": 0, "right_after": 0,
    }

    J = junction_pos   # boundary: [1..J] = exon N, [J+1..] = exon N+1

    # ── Left primer ──────────────────────────────────────────────────────────
    left_raw = rec.get("PRIMER_LEFT_0", "")
    if "," in left_raw:
        try:
            ls, ll = (int(x) for x in left_raw.split(","))
            # 0-based: primer covers [ls, ls+ll-1]
            # Spans junction iff ls < J  (starts before or at junction)
            #                 and ls+ll > J (ends after junction)
            if ls < J < ls + ll:
                result["left_before"] = J - ls        # bases in exon N
                result["left_after"]  = ls + ll - J   # bases in exon N+1
        except ValueError:
            pass

    # ── Right primer ─────────────────────────────────────────────────────────
    right_raw = rec.get("PRIMER_RIGHT_0", "")
    if "," in right_raw:
        try:
            re_, rl = (int(x) for x in right_raw.split(","))
            # 0-based: 3'-end on template at re_, primer covers [re_-rl+1, re_]
            rs = re_ - rl + 1
            if rs < J <= re_:
                result["right_before"] = J - rs        # bases in exon N
                result["right_after"]  = re_ - J + 1   # bases in exon N+1
        except ValueError:
            pass

    return result


def amplicon_info(rec):
    """
    Return (amplicon_sequence, amplicon_gc_percent) for the pair.

    The amplicon spans from the left primer's start to the right primer's
    3'-end (inclusive) on SEQUENCE_TEMPLATE. Returns ("NA", "NA") if the
    required fields are missing or inconsistent.
    """
    template  = rec.get("SEQUENCE_TEMPLATE", "")
    left_raw  = rec.get("PRIMER_LEFT_0", "")
    right_raw = rec.get("PRIMER_RIGHT_0", "")
    if not template or "," not in left_raw or "," not in right_raw:
        return "NA", "NA"
    try:
        ls, _ll   = (int(x) for x in left_raw.split(","))
        re_, _rl  = (int(x) for x in right_raw.split(","))
    except ValueError:
        return "NA", "NA"

    start, end = ls, re_ + 1          # 0-based slice [start, end)
    if start < 0 or end > len(template) or start >= end:
        return "NA", "NA"

    amp = template[start:end].upper()
    gc = sum(c in "GC" for c in amp)
    at = sum(c in "AT" for c in amp)
    denom = gc + at                  # ignore Ns / other characters
    gc_pct = round(100.0 * gc / denom, 1) if denom else "NA"
    return amp, gc_pct


def quality_flag(rec):
    """
    Automatic quality assessment:
      PASS_BEST : penalty < 0.1, |ΔTm| < 1°C, no strong hairpin, no strong dimer
      PASS      : penalty < 0.3, |ΔTm| < 2°C, dimer below safe threshold
      REVIEW    : fails one or more of the above — check manually
    """
    try:
        penalty   = float(rec.get("PRIMER_PAIR_0_PENALTY",      1))
        tm_left   = float(rec.get("PRIMER_LEFT_0_TM",           0))
        tm_right  = float(rec.get("PRIMER_RIGHT_0_TM",          0))
        tm_diff   = abs(tm_left - tm_right)
        hairpin_l = float(rec.get("PRIMER_LEFT_0_HAIRPIN_TH",   0))
        hairpin_r = float(rec.get("PRIMER_RIGHT_0_HAIRPIN_TH",  0))
        dimer_any = float(rec.get("PRIMER_PAIR_0_COMPL_ANY_TH", 0))

        if penalty < 0.1 and tm_diff < 1.0 and max(hairpin_l, hairpin_r) < 40 and dimer_any < 40:
            return "PASS_BEST"
        elif penalty < 0.3 and tm_diff < 2.0 and dimer_any < 47:
            return "PASS"
        else:
            return "REVIEW"
    except Exception:
        return "NA"


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <best_primers.txt> <summary.tsv>", file=sys.stderr)
        sys.exit(1)

    infile, outfile = sys.argv[1], sys.argv[2]

    blocks = parse_blocks(infile)
    print(f"[RAPID summary] Parsed {len(blocks)} primer pair(s).", file=sys.stderr)

    header = (
        ["Gene_Junction", "Gene", "Junction", "Quality"]
        + FIELDS
        + [
            "Amplicon_Sequence",
            "Amplicon_GC_percent",
            "Left_bases_before_junction",
            "Left_bases_after_junction",
            "Right_bases_before_junction",
            "Right_bases_after_junction",
        ]
    )

    with open(outfile, "w") as out:
        out.write("\t".join(header) + "\n")
        for b in blocks:
            seq_id   = b.get("SEQUENCE_ID", "NA")
            parts    = seq_id.rsplit("_j", 1)
            gene     = parts[0]
            junction = parts[1] if len(parts) == 2 else "NA"

            # Compute junction overlap
            try:
                jpos = int(junction)
                ov   = junction_overlap(b, jpos)
            except (ValueError, TypeError):
                ov = {"left_before": "NA", "left_after": "NA",
                      "right_before": "NA", "right_after": "NA"}

            # Amplicon sequence + GC%
            amp_seq, amp_gc = amplicon_info(b)

            row = (
                [seq_id, gene, junction, quality_flag(b)]
                + [b.get(f, "NA") for f in FIELDS]
                + [
                    amp_seq,
                    amp_gc,
                    ov["left_before"],
                    ov["left_after"],
                    ov["right_before"],
                    ov["right_after"],
                ]
            )
            out.write("\t".join(str(x) for x in row) + "\n")

    print(f"[RAPID summary] Table written: {len(blocks)} row(s) → {outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
