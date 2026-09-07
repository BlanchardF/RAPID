#!/usr/bin/env python3
"""
find_primer_amplicon.py  (standalone utility — NOT part of the RAPID pipeline)
Locate a primer pair on a sequence and report the amplicon (product) size.

For junction-spanning primers (RAPID's output), search the SPLICED transcript
(top_genes/merge_exons.fasta), not the genome: the primers cross an intron and
will not match the genomic sequence contiguously.

Primer strings may contain annotation brackets like GTAC{G}GAAT — the brackets
are stripped before searching, and their position within the primer is reported
(useful e.g. to mark an LNA base or a junction base).

Usage:
    python3 find_primer_amplicon.py <fasta> "<forward>" "<reverse>" [--id gene:Smp_169190]
"""

import argparse
import re
import sys

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


def clean_primer(p):
    """Keep only nucleotides (removes { } spaces etc.), uppercase, U->T."""
    return re.sub(r"[^ACGTUacgtu]", "", p).upper().replace("U", "T")


def bracket_offset(primer_raw):
    """Number of nucleotides before the first '{' (position of the marked base)."""
    before = primer_raw.split("{", 1)[0]
    return len(clean_primer(before)) if "{" in primer_raw else None


def read_fasta(path):
    recs, name, seq = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            recs[name] = "".join(seq)
    return recs


def find_all(hay, needle):
    idxs, start = [], 0
    while needle:
        i = hay.find(needle, start)
        if i < 0:
            break
        idxs.append(i)
        start = i + 1
    return idxs


def main():
    ap = argparse.ArgumentParser(description="Locate a primer pair and compute amplicon size.")
    ap.add_argument("fasta")
    ap.add_argument("forward")
    ap.add_argument("reverse")
    ap.add_argument("--id", default=None, help="Restrict to this FASTA record id "
                                               "(matches with or without the 'gene:' prefix).")
    a = ap.parse_args()

    F_raw, R_raw = a.forward, a.reverse
    F, R = clean_primer(F_raw), clean_primer(R_raw)
    rc = revcomp(R)
    jF, jR = bracket_offset(F_raw), bracket_offset(R_raw)

    recs = read_fasta(a.fasta)
    if a.id:
        want = a.id.split(":")[-1]
        keys = [k for k in recs if k == a.id or k.split(":")[-1] == want]
        if not keys:
            print(f"Record '{a.id}' not found. Present: {list(recs)[:8]}...", file=sys.stderr)
            sys.exit(1)
        recs = {k: recs[k] for k in keys}

    print(f"Forward : {F}  ({len(F)} nt)" + (f"  [marked base after position {jF}]" if jF else ""))
    print(f"Reverse : {R}  ({len(R)} nt)" + (f"  [marked base after position {jR}]" if jR else ""))
    print(f"Reverse (rev-comp, as expected on the + strand): {rc}\n")

    any_found = False
    for name, seq in recs.items():
        s = seq.upper()
        f_hits = find_all(s, F)
        r_hits = find_all(s, rc)
        if not f_hits and not r_hits:
            continue
        print(f">>> {name}  (spliced length: {len(seq)} nt)")
        print(f"    forward found at position : {[i+1 for i in f_hits] or 'NO'}")
        print(f"    reverse (rev-comp) at pos.: {[j+1 for j in r_hits] or 'NO'}")
        for i in f_hits:
            for j in r_hits:
                if j + len(rc) > i:            # reverse downstream of forward
                    any_found = True
                    product = (j + len(rc)) - i
                    amp = seq[i:j + len(rc)]
                    print(f"\n    > AMPLICON")
                    print(f"      forward : {i+1}..{i+len(F)}")
                    print(f"      reverse : {j+1}..{j+len(rc)}")
                    print(f"      SIZE    : {product} bp")
                    if jF is not None:
                        print(f"      marked base (forward) at position {i+jF+1} of the spliced sequence")
                    print(f"      amplified sequence:\n      {amp}")
        if not any_found and f_hits and r_hits:
            print("    (both primers found, but reverse upstream of forward "
                  "— check the F/R orientation)")
        print()

    if not any_found:
        print("No amplicon found on this/these sequence(s).", file=sys.stderr)
        print("Make sure you are searching the SPLICED sequence "
              "(top_genes/merge_exons.fasta), not the genome — junction-spanning "
              "primers do not match genomic DNA contiguously.", file=sys.stderr)


if __name__ == "__main__":
    main()
