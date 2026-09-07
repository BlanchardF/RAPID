# =============================================================================
# RAPID — track.smk   (mode: track)
# Like `auto`, but with one REFERENCE RNA-seq sample and one or more OTHER
# samples (different stages / treatments of the same species). Each sample is
# aligned and quantified independently, then a reference-specificity filter
# keeps only genes specific to the reference. From there on, the pipeline is
# identical to `auto` (common.smk), plus an enriched primer report at the end.
#
# Filter modes (config["filter_mode"]):
#   absolute : ref count >= ref_min_count AND every other <= other_max_count
#   cpm      : CPM_ref >= min_cpm_ref AND CPM_ref >= fold_change * CPM_other
#              for every other sample (depth-normalised enrichment)
# =============================================================================

import os as _os
import re
from pathlib import Path

# ── Config ────────────────────────────────────────────────────────────────────
OUT         = config["output_dir"]
GENOME      = config["genome"]
ANNOT_IN    = config.get("annotation")
RUN_BRAKER  = config["run_braker"]
THREADS     = config["threads"]
EX_MODE     = config.get("ex_mode", "max")
PRODUCT_MIN = config.get("primer_product_min", 100)
PRODUCT_MAX = config.get("primer_product_max", 200)

# Filter parameters
FILTER_MODE = config.get("filter_mode", "absolute")
REF_MIN     = config.get("ref_min_count", 10)
OTHER_MAX   = config.get("other_max_count", 0)
MIN_CPM_REF = config.get("min_cpm_ref", 1.0)
FOLD_CHANGE = config.get("fold_change", 10.0)

SCRIPTS_DIR = Path(workflow.basedir)
R_SCRIPT    = str(SCRIPTS_DIR / "01_expression_filter_and_primer3_prep.R")

# Annotation format, detected once (see Snakefile for the rationale).
_annot_source = (ANNOT_IN if ANNOT_IN else "").lower()
IS_GTF        = (not ANNOT_IN) or _annot_source.endswith(".gtf")

# Reproducible HISAT2 (see Snakefile note): supply known splice sites from the
# annotation and disable the dynamic table so every sample aligns
# deterministically regardless of thread count. GTF only.
if IS_GTF:
    _SS_FILE  = f"{OUT}/hisat2_index/splicesites.txt"
    _SS_INPUT = _SS_FILE
    _SS_FLAGS = f"--known-splicesite-infile {_SS_FILE} --no-temp-splicesite"
else:
    _SS_INPUT = []
    _SS_FLAGS = ""

# ── Samples ───────────────────────────────────────────────────────────────────
# config["samples"] = list of dicts: {name, role ("reference"|"other"), r1, r2}
SAMPLES = {}
for _s in config["samples"]:
    SAMPLES[_s["name"]] = _s

REF_NAME    = next(n for n, s in SAMPLES.items() if s["role"] == "reference")
OTHER_NAMES = [n for n, s in SAMPLES.items() if s["role"] == "other"]
ALL_NAMES   = list(SAMPLES.keys())

# Constrain the {sample} wildcard to the known sample names.
wildcard_constraints:
    sample = "|".join(re.escape(n) for n in ALL_NAMES)


def sample_reads(name):
    """List of read files for a sample: [r1] (SE) or [r1, r2] (PE)."""
    s = SAMPLES[name]
    reads = [s["r1"]]
    if s.get("r2"):
        reads.append(s["r2"])
    return reads


# ── Final targets ─────────────────────────────────────────────────────────────
rule all:
    input:
        f"{OUT}/primer3/results/best_primers.txt",
        f"{OUT}/primer3/results/primers_summary.tsv",
        f"{OUT}/primer3/results/primer_report.tsv",
        f"{OUT}/primer3/results/primer_report.html",


# =============================================================================
# STEP 1 — Annotation: Braker3 (all samples as evidence) OR link provided
# =============================================================================
if RUN_BRAKER:
    _BRAKER_WORKDIR = f"{OUT}/braker3"
    _link_lines = [f"mkdir -p {_BRAKER_WORKDIR}/rnaseq"]
    _set_ids    = []
    for _n in ALL_NAMES:
        _s = SAMPLES[_n]
        _link_lines.append(
            f"ln -sf $(realpath {_s['r1']}) {_BRAKER_WORKDIR}/rnaseq/{_n}_R1.fastq.gz 2>/dev/null || true"
        )
        if _s.get("r2"):
            _link_lines.append(
                f"ln -sf $(realpath {_s['r2']}) {_BRAKER_WORKDIR}/rnaseq/{_n}_R2.fastq.gz 2>/dev/null || true"
            )
        _set_ids.append(_n)
    _link_block  = "\n        ".join(_link_lines)
    _set_ids_str = ",".join(_set_ids)

    rule braker3:
        input:
            genome = GENOME,
            reads  = [r for n in ALL_NAMES for r in sample_reads(n)],
        output: annot = f"{OUT}/braker3/braker.gtf"
        log:    f"{OUT}/logs/braker3.log"
        threads: THREADS
        shell:
            _link_block + "\n" + f"""
        braker.pl \\
            --genome={{input.genome}} \\
            --softmasking \\
            --cores={{threads}} \\
            --workingdir={_BRAKER_WORKDIR} \\
            --rnaseq_sets_ids={_set_ids_str} \\
            --rnaseq_sets_dir={_BRAKER_WORKDIR}/rnaseq \\
            &> {{log}}
        """
    ANNOTATION = f"{OUT}/braker3/braker.gtf"

else:
    _annot_ext  = _os.path.splitext(ANNOT_IN)[1]
    _annot_link = f"{OUT}/annotation/annotation{_annot_ext}"

    rule link_annotation:
        input:  annot = ANNOT_IN
        output: annot = _annot_link
        shell:  "mkdir -p {OUT}/annotation && ln -sf $(realpath {input.annot}) {output.annot}"

    ANNOTATION = _annot_link


# =============================================================================
# STEP 2 — Index genome with HISAT2 (once)
# =============================================================================
rule hisat2_build:
    input:  genome = GENOME
    output: touch(f"{OUT}/hisat2_index/index.done")
    log:    f"{OUT}/logs/hisat2_build.log"
    threads: THREADS
    params: prefix = f"{OUT}/hisat2_index/genome"
    shell:
        """
        mkdir -p {OUT}/hisat2_index
        hisat2-build -p {threads} {input.genome} {params.prefix} &> {log}
        """


# =============================================================================
# STEP 2b — Extract known splice sites from the annotation (GTF only, once)
# =============================================================================
if IS_GTF:
    rule hisat2_extract_splicesites:
        input:  annot = ANNOTATION
        output: ss    = _SS_FILE
        log:    f"{OUT}/logs/hisat2_splicesites.log"
        shell:
            """
            mkdir -p {OUT}/hisat2_index
            hisat2_extract_splice_sites.py {input.annot} > {output.ss} 2> {log}
            echo "[RAPID] Known splice sites: $(wc -l < {output.ss})" >> {log}
            """


# =============================================================================
# STEP 3 — Align each sample with HISAT2 (deterministic)
# =============================================================================
rule hisat2_align:
    input:
        index_done  = f"{OUT}/hisat2_index/index.done",
        splicesites = _SS_INPUT,                       # [] when annotation is GFF3
        reads       = lambda wc: sample_reads(wc.sample),
    output: sam = f"{OUT}/hisat2/{{sample}}.sam"
    log:    f"{OUT}/logs/hisat2_align_{{sample}}.log"
    threads: THREADS
    params:
        prefix     = f"{OUT}/hisat2_index/genome",
        ss_flags   = _SS_FLAGS,
        reads_flag = lambda wc: (
            f"-1 {SAMPLES[wc.sample]['r1']} -2 {SAMPLES[wc.sample]['r2']}"
            if SAMPLES[wc.sample].get("r2")
            else f"-U {SAMPLES[wc.sample]['r1']}"
        ),
    shell:
        """
        mkdir -p {OUT}/hisat2
        hisat2 -p {threads} -x {params.prefix} {params.ss_flags} \\
            {params.reads_flag} -S {output.sam} &> {log}
        """


# =============================================================================
# STEP 4 — SAM → sorted BAM (per sample)
# =============================================================================
rule samtools_sort:
    input:  sam = f"{OUT}/hisat2/{{sample}}.sam"
    output:
        bam = f"{OUT}/samtools/{{sample}}_sorted.bam",
        bai = f"{OUT}/samtools/{{sample}}_sorted.bam.bai",
    log:    f"{OUT}/logs/samtools_sort_{{sample}}.log"
    threads: THREADS
    shell:
        """
        mkdir -p {OUT}/samtools
        samtools sort -@ {threads} -o {output.bam} {input.sam} &> {log}
        samtools index {output.bam} &>> {log}
        """


# =============================================================================
# STEP 5 — featureCounts (per sample)
# =============================================================================
if IS_GTF:
    _fc_fmt  = ""
    _fc_attr = "-g gene_id"
else:
    _fc_fmt  = "-F GFF"
    _fc_attr = "-g Parent"

rule featurecounts:
    input:
        bam   = f"{OUT}/samtools/{{sample}}_sorted.bam",
        annot = ANNOTATION,
    output: counts = f"{OUT}/featurecounts/counts_{{sample}}.txt"
    log:    f"{OUT}/logs/featurecounts_{{sample}}.log"
    threads: THREADS
    params:
        paired_flag = lambda wc: "-p -B" if SAMPLES[wc.sample].get("r2") else "",
        fmt         = _fc_fmt,
        attr        = _fc_attr,
    shell:
        """
        mkdir -p {OUT}/featurecounts
        featureCounts -T {threads} {params.paired_flag} \\
            {params.fmt} -t exon {params.attr} \\
            --primary \\
            -a {input.annot} -o {output.counts} {input.bam} &> {log}
        """


# =============================================================================
# STEP 6 — Reference-specificity filter (absolute or cpm)
# =============================================================================
rule track_filter_genes:
    input:
        ref    = f"{OUT}/featurecounts/counts_{REF_NAME}.txt",
        others = expand(f"{OUT}/featurecounts/counts_{{s}}.txt", s=OTHER_NAMES),
    output: counts = f"{OUT}/featurecounts/counts_ref_specific.txt"
    log:    f"{OUT}/logs/track_filter_genes.log"
    params:
        script      = str(SCRIPTS_DIR / "track_filter_genes.py"),
        mode        = FILTER_MODE,
        ref_min     = REF_MIN,
        other_max   = OTHER_MAX,
        min_cpm_ref = MIN_CPM_REF,
        fold_change = FOLD_CHANGE,
    shell:
        """
        mkdir -p {OUT}/logs
        python3 {params.script} \\
            --mode {params.mode} \\
            --ref {input.ref} --others {input.others} \\
            --ref-min {params.ref_min} --other-max {params.other_max} \\
            --min-cpm-ref {params.min_cpm_ref} --fold-change {params.fold_change} \\
            --out {output.counts} &> {log}
        cat {log}
        """


# =============================================================================
# Downstream (shared with auto): gene selection → primer design → summary
# =============================================================================
COUNTS = f"{OUT}/featurecounts/counts_ref_specific.txt"
include: "common.smk"


# =============================================================================
# STEP 7 — Enriched primer report (track-only: needs per-stage CPM)
# =============================================================================
# Merges the kept primer pairs (best_primers.txt) with per-stage expression
# (CPM + enrichment) into a sortable TSV and an interactive HTML explorer where
# primers can be hidden/kept while reviewing. Uses all per-sample counts and the
# track config (for stage labels). Reference stage first.
rule primer_report:
    input:
        best   = f"{OUT}/primer3/results/best_primers.txt",
        counts = expand(f"{OUT}/featurecounts/counts_{{s}}.txt", s=ALL_NAMES),
        config = f"{OUT}/rapid_track_config.yaml",
    output:
        tsv  = f"{OUT}/primer3/results/primer_report.tsv",
        html = f"{OUT}/primer3/results/primer_report.html",
    log:    f"{OUT}/logs/primer_report.log"
    params:
        script = str(SCRIPTS_DIR / "primer_report.py"),
        outdir = OUT,
    shell:
        """
        python3 {params.script} {params.outdir} \\
            --out-tsv {output.tsv} --out-html {output.html} &> {log}
        cat {log}
        """
