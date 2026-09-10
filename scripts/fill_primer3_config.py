#!/usr/bin/env python3
"""
RAPID — scripts/fill_primer3_config.py
Fill PLACEHOLDER sequences in a Primer3 pre-config file AND split the result
into ONE config file per gene.
"""

import re
import shutil
import sys
from pathlib import Path


def load_fasta(fasta_path):
    """
    Load a FASTA file into a dict {header_id: uppercase_sequence}.
    header_id is the first word after '>'.
    Sequences are uppercased to remove soft-masking artefacts.
    """
    seqs = {}
    current_id = None
    current_parts = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_parts).upper()
                current_id    = line[1:].split()[0]   # first word = ID
                current_parts = []
            else:
                current_parts.append(line)

    if current_id is not None:
        seqs[current_id] = "".join(current_parts).upper()

    return seqs


def gene_of(seq_id):
    """Gene id = SEQUENCE_ID with the trailing '_j<digits>' stripped."""
    parts = seq_id.rsplit("_j", 1)
    return parts[0] if len(parts) == 2 else seq_id


def sanitize_gene_name(gene, used_names):
    """Filesystem-safe, unique file stem for a gene id."""
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", gene).strip("_") or "gene"
    name = safe
    i = 1
    while name in used_names:
        i += 1
        name = f"{safe}_{i}"
    return name


def fill_and_split(preconfig_path, fasta_seqs, extra_params, out_dir):
    """
    Read the preconfig, replace PLACEHOLDER with actual sequences, inject
    extra_params before each '=' separator, and write one file per gene into
    out_dir. Returns (stats, gene_paths) where gene_paths is a dict
    {gene_id: Path}.
    """
    out_dir = Path(out_dir)
    if out_dir.exists():
        shutil.rmtree(out_dir)          # fresh start; avoid stale files
    out_dir.mkdir(parents=True)

    stats = {"total": 0, "filled": 0, "missing": 0}
    used_names = set()
    gene_paths = {}                     # gene -> Path (first-seen order)

    def path_for(gene):
        if gene not in gene_paths:
            name = sanitize_gene_name(gene, used_names)
            used_names.add(name)
            gene_paths[gene] = out_dir / f"{name}.txt"
        return gene_paths[gene]

    current_gene = None
    cur_handle   = None

    with open(preconfig_path) as fin:
        for line in fin:
            line = line.rstrip("\n")

            if line.startswith("SEQUENCE_ID="):
                seq_id = line[len("SEQUENCE_ID="):]
                gene   = gene_of(seq_id)
                stats["total"] += 1
                # switch output file when the gene changes
                if gene != current_gene:
                    if cur_handle is not None:
                        cur_handle.close()
                    cur_handle   = open(path_for(gene), "a")
                    current_gene = gene
                cur_handle.write(line + "\n")

            elif line == "SEQUENCE_TEMPLATE=PLACEHOLDER":
                seq = fasta_seqs.get(current_gene, "")
                if seq:
                    cur_handle.write(f"SEQUENCE_TEMPLATE={seq}\n")
                    stats["filled"] += 1
                else:
                    # Empty template — Primer3 will report no primers found
                    cur_handle.write("SEQUENCE_TEMPLATE=\n")
                    stats["missing"] += 1
                    print(f"  [WARNING] No sequence found for gene: {current_gene}",
                          file=sys.stderr)

            elif line == "=":
                # Inject extra params before the block separator
                for param in extra_params:
                    cur_handle.write(param + "\n")
                cur_handle.write("=\n")

            else:
                if cur_handle is not None:
                    cur_handle.write(line + "\n")

    if cur_handle is not None:
        cur_handle.close()

    stats["genes"] = len(gene_paths)
    return stats, gene_paths


def main():
    if len(sys.argv) < 4:
        print(
            f"Usage: {sys.argv[0]} <preconfig> <merged_fasta> <output_dir> "
            "[EXTRA_PARAM=VALUE ...]",
            file=sys.stderr,
        )
        sys.exit(1)

    preconfig_path = sys.argv[1]
    fasta_path     = sys.argv[2]
    out_dir        = sys.argv[3]
    extra_params   = sys.argv[4:]   # list of "KEY=VALUE" strings

    print(f"[fill_primer3_config] Loading {fasta_path} …", file=sys.stderr)
    fasta_seqs = load_fasta(fasta_path)
    print(f"[fill_primer3_config] {len(fasta_seqs)} sequences loaded.", file=sys.stderr)

    if extra_params:
        print(f"[fill_primer3_config] Extra Primer3 params: {extra_params}", file=sys.stderr)

    stats, gene_paths = fill_and_split(preconfig_path, fasta_seqs, extra_params, out_dir)

    print(
        f"[fill_primer3_config] Done — {stats['filled']}/{stats['total']} blocks filled"
        + (f", {stats['missing']} missing" if stats["missing"] else "")
        + f" → {stats['genes']} per-gene config file(s) in {out_dir}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
