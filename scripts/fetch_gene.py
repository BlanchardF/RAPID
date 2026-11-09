#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fetch_gene.py - backend of `rapid fetch_gene`.

Extracts exon sequences, writes an analysis report and a FASTA file ready to be
aligned in SeaView (spliced mRNA, exons, amplicon, primers and probe), for one
or all gene/junction pairs of primers_summary.tsv.
"""

import argparse
import os
import re
import sys
import subprocess
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("Error: PyYAML is required. Install it using 'pip install pyyaml'")


# Config file names (auto mode, then track mode)
CONFIG_NAMES = ("rapid_config.yaml", "rapid_track_config.yaml")


class AnalysisError(Exception):
    """Recoverable error: lets batch mode continue with the next junction."""


def parse_args(argv=None, prog=None):
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Extracts exon sequences, generates an analysis report and a FASTA "
            "alignment file for SeaView (exons, amplicon, primers and probe). "
            "Without -g/-j, processes every gene/junction found in "
            "primers_summary.tsv. Works with RAPID auto and track outputs."
        )
    )
    parser.add_argument("-d", "--dir", required=True, help="Results directory (e.g., Results)")
    parser.add_argument(
        "-g", "--gene", default=None,
        help="Gene ID (e.g., g9316, 9316 or Smp_123456). If omitted, all genes are processed."
    )
    parser.add_argument(
        "-j", "--junction", default=None,
        help="Junction ID/position (e.g., 183). If omitted, all junctions "
             "of the selected gene(s) are processed."
    )
    args = parser.parse_args(argv)
    if args.junction is not None and args.gene is None:
        parser.error("-j/--junction requires -g/--gene")
    return args


def clean_gene(gene_id):
    """Removes the 'gene:' prefix (safe for file names and FASTA headers)."""
    gene_id = str(gene_id)
    return gene_id[5:] if gene_id.startswith("gene:") else gene_id


def canon_gene(gene_id):
    """'gene:g9316', 'g9316' and '9316' are the same gene; 'Smp_...' IDs are kept as is."""
    gene_id = clean_gene(gene_id)
    if gene_id.startswith("g") and gene_id[1:].isdigit():
        return gene_id[1:]
    return gene_id


def reverse_complement(seq):
    complement = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
    return seq.translate(complement)[::-1]


def load_config(results_dir):
    for name in CONFIG_NAMES:
        candidate = results_dir / name
        if candidate.exists():
            with open(candidate, "r") as f:
                config = yaml.safe_load(f) or {}
            print(f"[*] Using config: {candidate.name}")
            return config
    sys.exit(f"Error: Unable to find any of {', '.join(CONFIG_NAMES)} in {results_dir}")


def get_genome_path(config):
    genome_path = config.get("genome")
    if not genome_path or not os.path.exists(genome_path):
        sys.exit(f"Error: Genome file not found in config: {genome_path}")
    return genome_path


def load_primers(results_dir):
    """Returns all primer rows (dicts) from primers_summary.tsv."""
    tsv_file = results_dir / "primer3" / "results" / "primers_summary.tsv"
    if not tsv_file.exists():
        sys.exit(f"Error: {tsv_file} does not exist.")

    rows = []
    with open(tsv_file, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            rows.append(dict(zip(header, parts)))
    return rows, tsv_file


def select_primers(rows, tsv_file, gene_id, junc_id):
    """Selects the rows to process according to -g / -j."""
    selected = []
    for row in rows:
        values = list(row.values())
        row_key, row_gene, row_junc = values[0], values[1], values[2]
        if gene_id is None:
            selected.append((row_gene, row_junc, row))
        elif canon_gene(row_gene) == canon_gene(gene_id):
            if junc_id is None:
                selected.append((row_gene, row_junc, row))
            elif row_junc == str(junc_id) or row_key.endswith(f"_j{junc_id}"):
                selected.append((row_gene, row_junc, row))
                break

    if not selected:
        target = "any gene" if gene_id is None else gene_id
        if junc_id is not None:
            target += f" and junction {junc_id}"
        sys.exit(f"Error: No primer found for {target} in {tsv_file}")
    return selected


def find_counts_files(results_dir, config):
    """
    Returns the featureCounts files usable for gene/exon coordinates, by priority.
    - auto mode : featurecounts/counts.txt
    - track mode: counts_<reference sample>.txt, counts_ref.txt, counts_ref_specific.txt
    """
    fc_dir = results_dir / "featurecounts"
    candidates = [fc_dir / "counts.txt"]

    for sample in config.get("samples", []) or []:
        if sample.get("role") == "reference":
            candidates.append(fc_dir / f"counts_{sample.get('name')}.txt")
    candidates += [fc_dir / "counts_ref.txt", fc_dir / "counts_ref_specific.txt"]

    found, seen = [], set()
    for c in candidates:
        if c.exists() and c not in seen:
            found.append(c)
            seen.add(c)
    if not found:
        sys.exit(f"Error: no featureCounts file found in {fc_dir}")
    return found


def load_counts(counts_files):
    """Parses the featureCounts files and returns {canonical_gene_id: gene_data}."""
    genes = {}
    for counts_file in counts_files:
        print(f"[*] Reading gene coordinates from: {counts_file.name}")
        with open(counts_file, "r") as f:
            for line in f:
                if line.startswith("#") or line.startswith("Geneid"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                key = canon_gene(parts[0])
                if key in genes:  # priority to the first file
                    continue
                chrs = parts[1].split(";")
                starts = [int(x) for x in parts[2].split(";")]
                ends = [int(x) for x in parts[3].split(";")]
                strands = parts[4].split(";")
                genes[key] = {
                    "gene_id": parts[0],
                    "chr": chrs[0],
                    "strand": strands[0],
                    "exons": sorted(zip(starts, ends), key=lambda x: x[0]),
                    "gene_start": min(starts),
                    "gene_end": max(ends),
                }
    return genes


def extract_sequence_samtools(genome_fasta, chrom, start, end):
    region = f"{chrom}:{start}-{end}"
    cmd = ["samtools", "faidx", genome_fasta, region]
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        lines = res.stdout.strip().split("\n")
        return "".join(lines[1:])
    except (subprocess.CalledProcessError, FileNotFoundError):
        return extract_sequence_pyfaidx(genome_fasta, chrom, start, end)


def extract_sequence_pyfaidx(genome_fasta, chrom, start, end):
    try:
        from pyfaidx import Fasta
        genes_seq = Fasta(genome_fasta)
        return str(genes_seq[chrom][start - 1:end])
    except ImportError:
        sys.exit("Error: Please install 'samtools' in your PATH or run 'pip install pyfaidx'.")


def process_junction(out_dir, genome_path, gene_id, junc_id, primer_info, gene_data, use_probe):
    name = clean_gene(gene_id)
    report_file = out_dir / f"{name}_j{junc_id}_report.txt"
    align_file = out_dir / f"{name}_j{junc_id}_for_align.txt"

    chrom = gene_data["chr"]
    strand = gene_data["strand"]

    ordered_exons = (
        gene_data["exons"]
        if strand == "+"
        else sorted(gene_data["exons"], key=lambda x: x[0], reverse=True)
    )

    exon_seqs = []
    for i, (s, e) in enumerate(ordered_exons, 1):
        seq = extract_sequence_samtools(genome_path, chrom, s, e)
        if strand == "-":
            seq = reverse_complement(seq)
        exon_seqs.append((i, s, e, seq))

    full_mrna_seq = "".join(seq for _, _, _, seq in exon_seqs)
    amplicon_seq = primer_info.get("Amplicon_Sequence", "")

    left_primer = primer_info.get("PRIMER_LEFT_0_SEQUENCE", "")
    right_primer = primer_info.get("PRIMER_RIGHT_0_SEQUENCE", "")
    probe_seq = primer_info.get("PRIMER_INTERNAL_0_SEQUENCE", "") if use_probe else ""
    probe_tm = primer_info.get("PRIMER_INTERNAL_0_TM")

    # ---------------- Report ----------------
    with open(report_file, "w") as f:
        f.write("=" * 80 + "\n")
        f.write(f"ANALYSIS REPORT: {name} | JUNCTION {junc_id}\n")
        f.write("=" * 80 + "\n\n")

        f.write("--- GENERAL INFORMATION ---\n")
        f.write(f"Chromosome: {chrom}\n")
        f.write(f"Strand: {strand}\n")
        f.write(f"Genomic coordinates: {gene_data['gene_start']} - {gene_data['gene_end']}\n")
        f.write(f"Total exons identified: {len(gene_data['exons'])}\n\n")

        f.write("--- PRIMER INFORMATION ---\n")
        f.write(f"Left Primer  : {left_primer} (Tm: {primer_info.get('PRIMER_LEFT_0_TM')}°C)\n")
        f.write(f"Right Primer : {right_primer} (Tm: {primer_info.get('PRIMER_RIGHT_0_TM')}°C)\n")
        if probe_seq:
            f.write(f"Probe        : {probe_seq} (Tm: {probe_tm}°C, "
                    f"size: {primer_info.get('Probe_Size')} bp, "
                    f"starts with G: {primer_info.get('Probe_Starts_With_G')})\n")
        f.write(f"Product Size : {primer_info.get('PRIMER_PAIR_0_PRODUCT_SIZE')} bp\n")
        f.write(f"Amplicon Seq : {amplicon_seq}\n\n")

        f.write("--- INDIVIDUAL EXON SEQUENCES ---\n")
        for num, s, e, seq in exon_seqs:
            f.write(f">Exon_{num} | {chrom}:{s}-{e} (Strand {strand}) | Length: {len(seq)} bp\n")
            f.write(seq + "\n")
        f.write("\n")

        f.write("--- RECONSTRUCTED mRNA SEQUENCE ---\n")
        f.write(f">mRNA_{name} | Length: {len(full_mrna_seq)} bp\n")
        f.write(full_mrna_seq + "\n")

    # ---------------- FASTA for alignment (SeaView) ----------------
    with open(align_file, "w") as f:
        f.write(f">{name}_spliced_mRNA | {chrom}:{gene_data['gene_start']}-{gene_data['gene_end']} (Strand {strand})\n")
        f.write(full_mrna_seq + "\n")

        for num, s, e, seq in exon_seqs:
            f.write(f">Exon_{num}_{name} | {chrom}:{s}-{e}\n")
            f.write(seq + "\n")

        if amplicon_seq:
            f.write(f">Amplicon_{name}_j{junc_id}\n")
            f.write(amplicon_seq + "\n")

        # Primers / probe, written in the same orientation as the mRNA/amplicon
        if left_primer:
            f.write(f">Left_primer_{name}_j{junc_id} | Tm {primer_info.get('PRIMER_LEFT_0_TM')}\n")
            f.write(left_primer + "\n")
        if right_primer:
            f.write(f">Right_primer_revcomp_{name}_j{junc_id} | Tm {primer_info.get('PRIMER_RIGHT_0_TM')}\n")
            f.write(reverse_complement(right_primer) + "\n")
        if probe_seq:
            f.write(f">Probe_{name}_j{junc_id} | Tm {probe_tm}\n")
            f.write(probe_seq + "\n")

    return report_file, align_file


def main(argv=None, prog=None):
    args = parse_args(argv, prog)
    results_dir = Path(args.dir)
    gene_id = args.gene
    junc_id = args.junction

    out_dir = results_dir / "Find_genes"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load inputs once
    config = load_config(results_dir)
    genome_path = get_genome_path(config)
    use_probe = bool(config.get("design_probe", False))
    rows, tsv_file = load_primers(results_dir)
    targets = select_primers(rows, tsv_file, gene_id, junc_id)
    all_genes = load_counts(find_counts_files(results_dir, config))

    batch = len(targets) > 1
    if gene_id is None:
        print(f"[*] No gene/junction specified: processing all {len(targets)} entries.")
    print(f"[*] Probe included: {'yes' if use_probe else 'no (design_probe is false)'}")

    n_ok, failed = 0, []
    for g, j, primer_info in targets:
        print(f"[*] Analyzing gene {clean_gene(g)} (Junction {j})...")
        try:
            gene_data = all_genes.get(canon_gene(g))
            if gene_data is None:
                raise AnalysisError(f"Gene {g} not found in featureCounts files")
            report_file, align_file = process_junction(
                out_dir, genome_path, g, j, primer_info, gene_data, use_probe
            )
            print(f"[✔] Report successfully generated: {report_file}")
            print(f"[✔] SeaView alignment file generated: {align_file}")
            n_ok += 1
        except AnalysisError as err:
            if not batch:
                sys.exit(f"Error: {err}")
            print(f"[!] Skipped {clean_gene(g)}_j{j}: {err}", file=sys.stderr)
            failed.append(f"{clean_gene(g)}_j{j}")

    if batch:
        print(f"\n[*] Done: {n_ok}/{len(targets)} succeeded.")
        if failed:
            print(f"[!] Failed: {', '.join(failed)}", file=sys.stderr)


if __name__ == "__main__":
    main()
