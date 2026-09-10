#!/usr/bin/env python3
"""
RAPID — scripts/prepare_ipcress_input.py
Parse a primers TSV (rapid auto output OR generic 2-column TSV) and write
an ipcress input file.
"""

import sys


def load_rows(tsv_path):
    with open(tsv_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]
    return header, rows


def load_primers(tsv_path, top_c):
    header, rows = load_rows(tsv_path)
    is_rapid = "PRIMER_PAIR_0_PENALTY" in header

    if is_rapid:
        try:
            rows.sort(key=lambda r: float(r.get("PRIMER_PAIR_0_PENALTY", 9999)))
        except Exception:
            pass
        rows = rows[:top_c]
        pairs = []
        for r in rows:
            pid     = r.get("Gene_Junction", f"pair_{len(pairs) + 1}")
            forward = r.get("PRIMER_LEFT_0_SEQUENCE", "").upper()
            reverse = r.get("PRIMER_RIGHT_0_SEQUENCE", "").upper()
            if forward and reverse:
                pairs.append((pid, forward, reverse))
    else:
        pairs = []
        for i, r in enumerate(rows[:top_c]):
            vals = list(r.values())
            if len(vals) < 2:
                continue
            pairs.append((f"pair_{i + 1}", vals[0].upper(), vals[1].upper()))

    return pairs, is_rapid


def main():
    if len(sys.argv) != 6:
        print(f"Usage: {sys.argv[0]} <tsv> <out_ipcress> <top_c> <min_size> <max_size>",
              file=sys.stderr)
        sys.exit(1)

    tsv_path, out_path, top_c, min_size, max_size = sys.argv[1:6]
    top_c, min_size, max_size = int(top_c), int(min_size), int(max_size)

    pairs, is_rapid = load_primers(tsv_path, top_c)
    if not pairs:
        print("ERROR: no valid primer pairs found in input TSV.", file=sys.stderr)
        sys.exit(1)

    with open(out_path, "w") as fh:
        for pid, fwd, rev in pairs:
            fh.write(f"{pid}\t{fwd}\t{rev}\t{min_size}\t{max_size}\n")

    mode = "rapid auto output" if is_rapid else "generic TSV"
    print(f"[prepare_ipcress_input] Detected: {mode}")
    print(f"[prepare_ipcress_input] {len(pairs)} primer pair(s) written to {out_path}")


if __name__ == "__main__":
    main()
