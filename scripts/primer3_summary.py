#!/usr/bin/env python3
"""
RAPID — scripts/primer3_summary.py
Parses a Primer3 output file and produces a TSV summary table.
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
    # Probe / internal oligo — only present when `rapid auto --probe` was used.
    "PRIMER_INTERNAL_0_SEQUENCE",
    "PRIMER_INTERNAL_0_TM",
    "PRIMER_INTERNAL_0_GC_PERCENT",
    "PRIMER_INTERNAL_0_SELF_ANY_TH",
    "PRIMER_INTERNAL_0_SELF_END_TH",
    "PRIMER_INTERNAL_0_HAIRPIN_TH",
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
    Return amplicon_sequence, amplicon_gc_percent for the pair.
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


def probe_info(rec):
    """
    Return (Probe_Size, Probe_Tm_Offset_vs_Primers, Probe_Starts_With_G) for
    the internal oligo (probe) of this record, or ("NA", "NA", "NA") when no
    probe was designed (e.g. `rapid auto` was run without --probe, or Primer3
    could not find a suitable internal oligo for this gene/junction).

    Probe_Tm_Offset_vs_Primers = probe Tm - mean(left Tm, right Tm), i.e. how
    many °C above the primers the probe's Tm actually landed — this is what
    should sit in the +8 to +10°C window targeted by --probe-tm-offset-*.
    """
    seq = rec.get("PRIMER_INTERNAL_0_SEQUENCE")
    if not seq:
        return "NA", "NA", "NA"

    size = len(seq)
    starts_with_g = "yes" if seq.upper().startswith("G") else "no"

    try:
        tm_probe = float(rec.get("PRIMER_INTERNAL_0_TM", ""))
        tm_left  = float(rec.get("PRIMER_LEFT_0_TM", ""))
        tm_right = float(rec.get("PRIMER_RIGHT_0_TM", ""))
        offset = round(tm_probe - (tm_left + tm_right) / 2.0, 1)
    except ValueError:
        offset = "NA"

    return size, offset, starts_with_g


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
            "Probe_Size",
            "Probe_Tm_Offset_vs_Primers",
            "Probe_Starts_With_G",
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

            # Probe (internal oligo), if any — see probe_info() docstring.
            probe_size, probe_tm_offset, probe_starts_g = probe_info(b)

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
                    probe_size,
                    probe_tm_offset,
                    probe_starts_g,
                ]
            )
            out.write("\t".join(str(x) for x in row) + "\n")

    print(f"[RAPID summary] Table written: {len(blocks)} row(s) → {outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
