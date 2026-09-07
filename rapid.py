#!/usr/bin/env python3
"""
RAPID - RNA Automated Primer Identification and Design
=======================================================
Main dispatcher — usage: rapid <mode> [options]
"""

import argparse
import sys
import os
import re
import shutil
import yaml
import subprocess
import textwrap
from pathlib import Path


# =============================================================================
# ASCII art + banner
# =============================================================================

ASCII_ART = r"""
░                                                                 ░
  ░                                                             ▒  
   █▒                                                         ░░   
    ░░▒                                                     ░█     
      ▒░▒                                                 ░░░      
       ▒▓░▒                      ▒                      ░░▒        
   ░▒▒▒▒▓▓▒▒▒░░░░░▒▒▒▒░░░░░▒▒░░░▒▒▒░░▒▒▒▒▒▒▒▒▒▒░░░░░░░░▒▓▓▓▒▒▒▒░   
   ▓░▓██▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓████████████████████▓▓▓▓▓▓▓▓▓▓▓▓▓████▓░   
    ░▒▓▓▓▓▓███░░░░░░█████░░░█████░░░░░░███░░██░░░░░░░██▓▓▓▓█░▒     
     ▒░█▓▓▓▓▓░░░███░░███░░█░░███░░░███░░██░░██░░███▓░░▓▓▓▓██░▒     
     ░░█▓▓▓▓▓░░░░░░░███░░░░░░░██░░░░░░░░██░░██░░███▓░░▓▓▓▓█░░      
      ░░█▓▓▓▓▓░░███░░█░░█████░░█░░████████░░██░░░░░░░▓▓▓▓██░░      
      ░▒▒█▓▓█████████████████████████████████████████████▒█▒       
         ▓▒▓▓▒▒░▒▓███████████████████████████████▒░░▒▓▓▓▓░         
                 ░▒▓▓▒░░░░▓█████████████▒░▒▒▒▓▓█▓░                 
                      ▒░░░████▓▒▒ ▒▓▓██▓▒▒▒░░                      
                       ░░░░▒▓██  █  ▓▓▓░▒░░                        
                         ▒░░▓▒▓▓█ █▓▓▒▒░░▒                         
                          ░░░▒░▒▓▒▒▒▓▒░░                           
                            ▓░░▒░▒░▒░░▒                            
                             ░░░▒▒░░░                              
                               ▒░░░▓                               
                                ░░                                 

"""

BANNER = """
╔══════════════════════════════════════════════════════╗
║   RAPID — RNA Automated Primer Identification        ║
║             and Design                               ║
╚══════════════════════════════════════════════════════╝"""

MODES = {
    "auto":    "Full automated pipeline: alignment → expression → primer design",
    "track":   "Like auto, but keep only genes specific to a reference RNA-seq sample",
    "check":   "Validate primer candidates with ipcress against a sequence database",
    "install": "Create (or update) the 'rapid' conda environment",
    "clean":   "Clean Snakemake metadata to allow a fresh re-run",
}


def print_banner():
    print(ASCII_ART)
    print(BANNER)


def print_mode_help():
    print_banner()
    print("\n  Usage:  rapid <mode> [options]\n")
    print("  Available modes:\n")
    for mode, desc in MODES.items():
        print(f"    {mode:<12}  {desc}")
    print()
    print("  Run 'rapid <mode> --help' for mode-specific options.\n")


# =============================================================================
# Shared helpers (config, Snakemake launch, tool checks)
# =============================================================================

def write_yaml_config(config, output_dir, filename):
    config_path = Path(output_dir) / filename
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    with open(config_path, "w") as fh:
        yaml.dump(config, fh, default_flow_style=False, sort_keys=False)
    return config_path


def find_snakemake():
    if shutil.which("snakemake"):
        return shutil.which("snakemake"), None
    python_bin = Path(sys.executable).parent
    candidate = python_bin / "snakemake"
    if candidate.is_file():
        return str(candidate), None
    conda_exe = shutil.which("conda") or shutil.which("mamba")
    if conda_exe:
        try:
            r = subprocess.run([conda_exe, "run", "-n", "rapid", "which", "snakemake"],
                               capture_output=True, text=True, timeout=10)
            if r.returncode == 0 and r.stdout.strip():
                return r.stdout.strip(), "rapid"
        except Exception:
            pass
    return None, None


def check_tools(required):
    missing = [t for t in required if not shutil.which(t)]
    if missing:
        print("\n[RAPID] Tools not found in PATH:", file=sys.stderr)
        for t in missing:
            print(f"  ✗ {t}", file=sys.stderr)
        print("\n  Activate the rapid conda environment first:\n"
              "    conda activate rapid\n", file=sys.stderr)
        sys.exit(1)


def run_snakemake(config_path, extra_args, threads, snakefile_name="Snakefile"):
    """
    Launch Snakemake against <snakefile_name>, looked up in scripts/ first
    then at the project root (legacy flat layout).
    """
    base = Path(__file__).parent
    if (base / "scripts" / snakefile_name).is_file():
        snakefile_path = base / "scripts" / snakefile_name
    elif (base / snakefile_name).is_file():
        snakefile_path = base / snakefile_name
    else:
        print(f"[RAPID] ERROR: {snakefile_name} not found.", file=sys.stderr)
        sys.exit(1)

    sm_exe, conda_env = find_snakemake()
    if sm_exe is None:
        print("\n[RAPID] ERROR: snakemake not found.\n"
              "  conda activate rapid && rapid ...\n", file=sys.stderr)
        sys.exit(1)

    cmd = ([shutil.which("conda") or "conda", "run", "-n", conda_env,
             "--no-capture-output"] if conda_env else [])
    cmd += [sm_exe, "--snakefile", str(snakefile_path),
            "--configfile", str(config_path), "--cores", str(threads)]
    if extra_args:
        cmd += extra_args.split()

    print("\n[RAPID] Launching Snakemake pipeline …")
    print(f"  Command: {' '.join(cmd)}\n")
    sys.exit(subprocess.run(cmd).returncode)


# Tools required by the alignment → primer design pipeline (auto and track).
PIPELINE_TOOLS = ["hisat2", "samtools", "featureCounts",
                  "bedtools", "primer3_core", "parallel"]


def validate_expression_mode(ex):
    """Return an error string if ex is not a valid -ex mode, else None."""
    ex_num = None
    try:
        ex_num = int(ex)
    except ValueError:
        pass
    if ex_num is None and ex not in ("max", "+", "mid", "min"):
        return (f"Invalid -ex mode: '{ex}'. Use: max, +, mid, min, an integer "
                "(top N), or a negative integer (bottom N).")
    return None


# =============================================================================
# MODE: auto
# =============================================================================

def parse_auto_args():
    parser = argparse.ArgumentParser(
        prog="rapid auto",
        description="RAPID auto — Full automated primer design pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              rapid auto -g genome.fa -r R1.fq.gz R2.fq.gz -o results/
              rapid auto -g genome.fa -r reads.fq.gz -o results/ \\
                         -a braker.gtf --top-primers 20
              rapid auto ... -p "PRIMER_MIN_SIZE=18" -p "PRIMER_MAX_SIZE=27"
              rapid auto ... --primer3-timeout 3600   # 1h per-gene limit
        """),
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument("-g", "--genome",     required=True, metavar="FA",
                     help="Reference genome (FASTA).")
    req.add_argument("-r", "--rna",        required=True, nargs="+", metavar="READS",
                     help="RNA-seq reads: 1 file (single-end) or 2 files (paired-end R1 R2).")
    req.add_argument("-o", "--output",     required=True, metavar="DIR",
                     help="Output directory.")

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument("-a", "--annot",      metavar="GFF/GTF", default=None,
                     help="Existing annotation. If omitted, Braker3 runs automatically.")
    opt.add_argument("-p", "--primer3-param", metavar="KEY=VALUE", action="append",
                     default=[], dest="primer3_params",
                     help=("Extra Primer3 global parameter (repeatable). "
                           "E.g. -p \"PRIMER_MIN_SIZE=18\" -p \"PRIMER_MAX_SIZE=27\""))
    opt.add_argument("-ex", "--expression", default="max", metavar="MODE",
                     help=(
                         "Gene expression selection mode:\n"
                         "  max / +   : top 33%% most expressed (default)\n"
                         "  mid       : middle 33%% (33–66%%)\n"
                         "  min       : bottom 33%% least expressed\n"
                         "  N (int)   : top N most expressed genes  e.g. -ex 50\n"
                         "  -N (int)  : bottom N least expressed  e.g. -ex -20"
                     ))
    opt.add_argument("--primer-product-min", type=int, default=100, metavar="INT",
                     help="Minimum amplicon size for Primer3. Default: 100.")
    opt.add_argument("--primer-product-max", type=int, default=200, metavar="INT",
                     help="Maximum amplicon size for Primer3. Default: 200.")
    opt.add_argument("--top-primers",     type=int, default=10,  metavar="INT",
                     help="Number of best primer pairs to keep (ranked by penalty). Default: 10.")
    opt.add_argument("--primer3-jobs",    type=int, default=None, metavar="INT",
                     help="Parallel Primer3 jobs (default: same as --threads). "
                          "Primer3 is single-threaded; each job processes one gene.")
    opt.add_argument("--primer3-timeout", type=int, default=7200, metavar="SECONDS",
                     help="Per-gene wall-clock time limit for Primer3, in seconds. "
                          "A gene exceeding it is excluded and reported as failed, so a "
                          "pathological gene (e.g. a very long fused transcript with "
                          "thousands of junctions) can never stall the whole run. "
                          "Default: 7200 (2 hours). Use 0 to disable.")
    opt.add_argument("--threads",         type=int, default=8,   metavar="INT",
                     help="CPU threads per tool. Default: 8.")
    opt.add_argument("--snakemake-args",  metavar="'...'", default="",
                     help="Extra Snakemake flags (e.g. '--dryrun --forceall'). Quote the string.")

    return parser.parse_args(sys.argv[2:])


def validate_auto(args):
    errors = []
    if not Path(args.genome).is_file():
        errors.append(f"Genome not found: {args.genome}")
    for f in args.rna:
        if not Path(f).is_file():
            errors.append(f"RNA-seq file not found: {f}")
    if len(args.rna) > 2:
        errors.append("Provide 1 (single-end) or 2 (paired-end) RNA-seq files.")
    if args.annot and not Path(args.annot).is_file():
        errors.append(f"Annotation not found: {args.annot}")
    for param in args.primer3_params:
        if "=" not in param:
            errors.append(f"Invalid Primer3 param (expected KEY=VALUE): '{param}'")
    if args.primer_product_min >= args.primer_product_max:
        errors.append(
            f"--primer-product-min ({args.primer_product_min}) must be smaller than "
            f"--primer-product-max ({args.primer_product_max}). "
            "Primer3 requires PRIMER_PRODUCT_SIZE_RANGE=min-max with min < max."
        )
    if args.primer3_timeout < 0:
        errors.append(
            f"--primer3-timeout ({args.primer3_timeout}) must be >= 0 "
            "(0 disables the per-gene timeout)."
        )
    ex_err = validate_expression_mode(args.expression)
    if ex_err:
        errors.append(ex_err)
    if errors:
        print("\n[RAPID auto] Input errors:", file=sys.stderr)
        for e in errors:
            print(f"  ✗ {e}", file=sys.stderr)
        sys.exit(1)


def build_auto_config(args):
    paired = len(args.rna) == 2
    return {
        "genome":               str(Path(args.genome).resolve()),
        "rna_r1":               str(Path(args.rna[0]).resolve()),
        "rna_r2":               str(Path(args.rna[1]).resolve()) if paired else None,
        "paired":               paired,
        "output_dir":           str(Path(args.output).resolve()),
        "annotation":           str(Path(args.annot).resolve()) if args.annot else None,
        "run_braker":           args.annot is None,
        "ex_mode":              args.expression,
        "primer_product_min":   args.primer_product_min,
        "primer_product_max":   args.primer_product_max,
        "top_primers":          args.top_primers,
        "primer3_extra_params": args.primer3_params,
        "primer3_jobs":         args.primer3_jobs if args.primer3_jobs else args.threads,
        "primer3_timeout":      args.primer3_timeout,
        "threads":              args.threads,
    }


def run_auto():
    args = parse_auto_args()
    print_banner()

    paired_label = "paired-end" if len(args.rna) == 2 else "single-end"
    print(f"\n  Mode         : auto")
    print(f"  Genome       : {args.genome}")
    print(f"  RNA-seq      : {paired_label} — {', '.join(args.rna)}")
    print(f"  Output dir   : {args.output}")
    print(f"  Annotation   : {args.annot or 'not provided → Braker3 will run'}")
    print(f"  Threads      : {args.threads}")
    print(f"  Expression   : {args.expression}")
    timeout_label = f"{args.primer3_timeout}s" if args.primer3_timeout > 0 else "disabled"
    print(f"  Primer3 t/o  : {timeout_label} (per gene)")
    if args.primer3_params:
        print(f"  Primer3 extra: {', '.join(args.primer3_params)}")
    print()

    validate_auto(args)
    config = build_auto_config(args)
    config_path = write_yaml_config(config, args.output, "rapid_config.yaml")
    print(f"[RAPID] Config written to: {config_path}")
    check_tools(PIPELINE_TOOLS)
    run_snakemake(config_path, args.snakemake_args, args.threads, snakefile_name="Snakefile")


# =============================================================================
# MODE: track
# =============================================================================

def parse_track_args():
    parser = argparse.ArgumentParser(
        prog="rapid track",
        description=(
            "RAPID track — Like auto, but keep only genes specific to a reference\n"
            "RNA-seq sample.\n\n"
            "You provide one REFERENCE sample (--ref) and one or more OTHER samples\n"
            "(--other, repeatable), corresponding to different stages / treatments of\n"
            "the same species. Two filters are available (--filter):\n"
            "  absolute : keep genes expressed in ref (>= --ref-min-count) and absent\n"
            "             from every other (<= --other-max-count). Raw counts.\n"
            "  cpm      : keep genes with CPM_ref >= --min-cpm-ref and enriched at\n"
            "             least --fold-change x over EVERY other stage (depth-\n"
            "             normalised). Recommended for stage-specific markers.\n"
            "After filtering, the pipeline is identical to `rapid auto`."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            Examples:
              # depth-normalised enrichment (recommended): >=10x over every other stage
              rapid track -g genome.fa -o results/ -a annot.gtf \\
                  --filter cpm --min-cpm-ref 1 --fold-change 10 \\
                  --ref cerc_R1.fq.gz cerc_R2.fq.gz \\
                  --other stageB_R1.fq.gz stageB_R2.fq.gz \\
                  --other stageC_R1.fq.gz stageC_R2.fq.gz

              # simple presence/absence on raw counts
              rapid track -g genome.fa -o results/ -a annot.gtf \\
                  --filter absolute --other-max-count 5 \\
                  --ref ref.fq.gz --other other1.fq.gz --other other2.fq.gz
        """),
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument("-g", "--genome", required=True, metavar="FA",
                     help="Reference genome (FASTA).")
    req.add_argument("-o", "--output", required=True, metavar="DIR",
                     help="Output directory.")
    req.add_argument("--ref", required=True, nargs="+", metavar="READS",
                     help="Reference RNA-seq sample: 1 file (SE) or 2 files (R1 R2). "
                          "We keep genes expressed / enriched HERE.")
    req.add_argument("--other", required=True, nargs="+", action="append", metavar="READS",
                     help="An 'other' RNA-seq sample: 1 file (SE) or 2 files (R1 R2). "
                          "Repeat --other for each additional sample.")

    filt = parser.add_argument_group("filter options")
    filt.add_argument("--filter", choices=["absolute", "cpm"], default="absolute",
                      help="Filtering strategy. Default: absolute. "
                           "Use 'cpm' for depth-normalised enrichment.")
    filt.add_argument("--ref-min-count", type=int, default=10, metavar="INT",
                      help="[absolute] Min raw count in the reference. Default: 10.")
    filt.add_argument("--other-max-count", type=int, default=0, metavar="INT",
                      help="[absolute] Max raw count allowed in an 'other' sample. Default: 0.")
    filt.add_argument("--min-cpm-ref", type=float, default=1.0, metavar="FLOAT",
                      help="[cpm] Min CPM in the reference. Default: 1.0.")
    filt.add_argument("--fold-change", type=float, default=10.0, metavar="FLOAT",
                      help="[cpm] Required fold-enrichment of ref CPM over EACH other "
                           "stage. Default: 10.")

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument("-a", "--annot", metavar="GFF/GTF", default=None,
                     help="Existing annotation. If omitted, Braker3 runs automatically "
                          "using all provided samples as RNA-seq evidence.")
    opt.add_argument("-p", "--primer3-param", metavar="KEY=VALUE", action="append",
                     default=[], dest="primer3_params",
                     help="Extra Primer3 global parameter (repeatable).")
    opt.add_argument("-ex", "--expression", default="max", metavar="MODE",
                     help="Gene expression selection mode (max/+, mid, min, N, -N). "
                          "Ranking is on the reference sample. Default: max.")
    opt.add_argument("--primer-product-min", type=int, default=100, metavar="INT",
                     help="Minimum amplicon size for Primer3. Default: 100.")
    opt.add_argument("--primer-product-max", type=int, default=200, metavar="INT",
                     help="Maximum amplicon size for Primer3. Default: 200.")
    opt.add_argument("--top-primers", type=int, default=10, metavar="INT",
                     help="Number of best primer pairs to keep. Default: 10.")
    opt.add_argument("--primer3-jobs", type=int, default=None, metavar="INT",
                     help="Parallel Primer3 jobs (default: same as --threads).")
    opt.add_argument("--primer3-timeout", type=int, default=7200, metavar="SECONDS",
                     help="Per-gene Primer3 wall-clock limit in seconds. Default: 7200 (2h). "
                          "Use 0 to disable.")
    opt.add_argument("--threads", type=int, default=8, metavar="INT",
                     help="CPU threads per tool. Default: 8.")
    opt.add_argument("--snakemake-args", metavar="'...'", default="",
                     help="Extra Snakemake flags. Quote the string.")

    return parser.parse_args(sys.argv[2:])


def validate_track(args):
    errors = []
    if not Path(args.genome).is_file():
        errors.append(f"Genome not found: {args.genome}")

    # Reference sample
    if len(args.ref) > 2:
        errors.append("--ref takes 1 (single-end) or 2 (paired-end R1 R2) files.")
    for f in args.ref:
        if not Path(f).is_file():
            errors.append(f"Reference RNA-seq file not found: {f}")

    # Other samples (list of lists)
    for i, other in enumerate(args.other, start=1):
        if len(other) > 2:
            errors.append(f"--other #{i} takes 1 or 2 files (got {len(other)}).")
        for f in other:
            if not Path(f).is_file():
                errors.append(f"'other' RNA-seq file not found: {f}")

    if args.annot and not Path(args.annot).is_file():
        errors.append(f"Annotation not found: {args.annot}")
    for param in args.primer3_params:
        if "=" not in param:
            errors.append(f"Invalid Primer3 param (expected KEY=VALUE): '{param}'")
    if args.primer_product_min >= args.primer_product_max:
        errors.append(
            f"--primer-product-min ({args.primer_product_min}) must be smaller than "
            f"--primer-product-max ({args.primer_product_max})."
        )
    if args.primer3_timeout < 0:
        errors.append(f"--primer3-timeout ({args.primer3_timeout}) must be >= 0.")

    # Filter parameters
    if args.filter == "absolute":
        if args.ref_min_count < 0:
            errors.append(f"--ref-min-count ({args.ref_min_count}) must be >= 0.")
        if args.other_max_count < 0:
            errors.append(f"--other-max-count ({args.other_max_count}) must be >= 0.")
    else:  # cpm
        if args.min_cpm_ref < 0:
            errors.append(f"--min-cpm-ref ({args.min_cpm_ref}) must be >= 0.")
        if args.fold_change <= 0:
            errors.append(f"--fold-change ({args.fold_change}) must be > 0.")

    ex_err = validate_expression_mode(args.expression)
    if ex_err:
        errors.append(ex_err)

    if errors:
        print("\n[RAPID track] Input errors:", file=sys.stderr)
        for e in errors:
            print(f"  ✗ {e}", file=sys.stderr)
        sys.exit(1)


def build_track_config(args):
    def resolve(p):
        return str(Path(p).resolve())

    samples = []
    ref_paired = len(args.ref) == 2
    samples.append({
        "name": "ref",
        "role": "reference",
        "r1":   resolve(args.ref[0]),
        "r2":   resolve(args.ref[1]) if ref_paired else None,
    })
    for i, other in enumerate(args.other, start=1):
        paired = len(other) == 2
        samples.append({
            "name": f"other_{i}",
            "role": "other",
            "r1":   resolve(other[0]),
            "r2":   resolve(other[1]) if paired else None,
        })

    return {
        "genome":               resolve(args.genome),
        "output_dir":           resolve(args.output),
        "annotation":           resolve(args.annot) if args.annot else None,
        "run_braker":           args.annot is None,
        "samples":              samples,
        # filter
        "filter_mode":          args.filter,
        "ref_min_count":        args.ref_min_count,
        "other_max_count":      args.other_max_count,
        "min_cpm_ref":          args.min_cpm_ref,
        "fold_change":          args.fold_change,
        # downstream
        "ex_mode":              args.expression,
        "primer_product_min":   args.primer_product_min,
        "primer_product_max":   args.primer_product_max,
        "top_primers":          args.top_primers,
        "primer3_extra_params": args.primer3_params,
        "primer3_jobs":         args.primer3_jobs if args.primer3_jobs else args.threads,
        "primer3_timeout":      args.primer3_timeout,
        "threads":              args.threads,
    }


def run_track():
    args = parse_track_args()
    print_banner()

    ref_label = "paired-end" if len(args.ref) == 2 else "single-end"
    print(f"\n  Mode         : track")
    print(f"  Genome       : {args.genome}")
    print(f"  Reference    : {ref_label} — {', '.join(args.ref)}")
    print(f"  Other samples: {len(args.other)}")
    for i, other in enumerate(args.other, start=1):
        lbl = "paired-end" if len(other) == 2 else "single-end"
        print(f"     other_{i}  : {lbl} — {', '.join(other)}")
    print(f"  Output dir   : {args.output}")
    print(f"  Annotation   : {args.annot or 'not provided → Braker3 will run'}")
    if args.filter == "absolute":
        print(f"  Filter       : absolute — keep if ref >= {args.ref_min_count} "
              f"AND every other <= {args.other_max_count} (raw counts)")
    else:
        print(f"  Filter       : cpm — keep if CPM_ref >= {args.min_cpm_ref} "
              f"AND ref enriched >= {args.fold_change}x over every other")
    print(f"  Threads      : {args.threads}")
    print(f"  Expression   : {args.expression}")
    timeout_label = f"{args.primer3_timeout}s" if args.primer3_timeout > 0 else "disabled"
    print(f"  Primer3 t/o  : {timeout_label} (per gene)")
    if args.primer3_params:
        print(f"  Primer3 extra: {', '.join(args.primer3_params)}")
    print()

    validate_track(args)
    config = build_track_config(args)
    config_path = write_yaml_config(config, args.output, "rapid_track_config.yaml")
    print(f"[RAPID] Config written to: {config_path}")
    check_tools(PIPELINE_TOOLS)
    run_snakemake(config_path, args.snakemake_args, args.threads, snakefile_name="track.smk")


# =============================================================================
# MODE: check
# =============================================================================

PRIMERS_SUMMARY_RELPATH = "primer3/results/primers_summary.tsv"


def find_primers_summary(path_arg):
    """
    Accept either:
      - a directory (rapid auto/track output) → look for primers_summary.tsv inside
      - a TSV file directly
    Returns a Path to the TSV file.
    """
    p = Path(path_arg).resolve()
    if p.is_dir():
        candidate = p / PRIMERS_SUMMARY_RELPATH
        if candidate.is_file():
            return candidate
        raise FileNotFoundError(
            f"Could not find {PRIMERS_SUMMARY_RELPATH} inside '{p}'.\n"
            "  Run rapid auto/track first, or provide the TSV file directly with -i."
        )
    if p.is_file():
        return p
    raise FileNotFoundError(f"Input not found: {p}")


def sanitize_db_name(path, used_names):
    """Make a filesystem/wildcard-safe unique name from a database file path."""
    stem = Path(path).name
    stem = re.sub(r"\.gz$", "", stem)
    stem = re.sub(r"\.(fsa_nt|fasta|fa|fna)$", "", stem)
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", stem).strip("_") or "db"
    name = safe
    i = 1
    while name in used_names:
        i += 1
        name = f"{safe}_{i}"
    return name


def resolve_db_files(args):
    """
    Resolve the list of ipcress database files from -tsa or -d.
    Returns (db_files_dict, use_tsa_bool) where db_files_dict maps a
    sanitized wildcard-safe name to the resolved absolute path.
    """
    if args.tsa:
        tsa_dir = Path(args.tsa)
        if not tsa_dir.is_dir():
            print(f"[RAPID check] ERROR: TSA directory not found: {tsa_dir}", file=sys.stderr)
            sys.exit(1)
        files = sorted(tsa_dir.glob("tsa.*.fsa_nt.gz"))
        if not files:
            files = sorted(list(tsa_dir.glob("*.gz")) + list(tsa_dir.glob("*.fasta")))
        if not files:
            print(f"[RAPID check] ERROR: no database files found in {tsa_dir}", file=sys.stderr)
            sys.exit(1)
        use_tsa = True
    else:
        db = Path(args.database)
        if not db.is_file():
            print(f"[RAPID check] ERROR: database file not found: {db}", file=sys.stderr)
            sys.exit(1)
        files = [db]
        use_tsa = False

    used_names = set()
    db_files = {}
    for f in files:
        name = sanitize_db_name(f, used_names)
        used_names.add(name)
        db_files[name] = str(f.resolve())

    return db_files, use_tsa


def build_check_config(args, tsv_path, out_dir, db_files, use_tsa):
    return {
        "output_dir":      str(out_dir),
        "input_tsv":       str(tsv_path),
        "top_c":           args.c,
        "min_size":        args.min_size,
        "max_size":        args.max_size,
        "mismatch":        args.mismatch,
        "db_files":        db_files,
        "extract_species": use_tsa,
        "generate_maps":   args.maps,
    }


def run_check():
    parser = argparse.ArgumentParser(
        prog="rapid check",
        description=(
            "RAPID check — Validate primer candidates with ipcress (Snakemake-backed).\n\n"
            "Input: primers_summary.tsv from rapid auto/track (auto-detected if you\n"
            "pass the output directory), or any TSV with forward primer in\n"
            "column 1 and reverse primer in column 2.\n\n"
            "This mode is backed by Snakemake, like rapid auto: if ipcress is\n"
            "slow against a large database (e.g. TSA) and gets interrupted,\n"
            "simply re-run the same command — completed database files are\n"
            "cached and skipped. Re-running later with --maps added will only\n"
            "generate the maps, without re-running ipcress."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  rapid check -i ../Results/ -tsa /ceph/banques/blast/db/tsa/nucl/\n"
            "  rapid check -i ../Results/ -tsa /path/to/tsa/ --maps\n"
            "  rapid check -i my_primers.tsv -d sequences.fasta -c 5 -min 80 -max 300\n"
        ),
    )

    req = parser.add_argument_group("required arguments")
    req.add_argument("-i", "--input", required=True, metavar="DIR_OR_TSV",
                     help="rapid auto/track output directory OR a TSV file with primers.")

    db_grp = parser.add_mutually_exclusive_group(required=True)
    db_grp.add_argument("-tsa", metavar="TSA_DIR",
                        help="Directory containing TSA database files (tsa.*.fsa_nt.gz).")
    db_grp.add_argument("-d", "--database", metavar="FILE",
                        help="Single FASTA / FASTA.gz database file.")

    opt = parser.add_argument_group("optional arguments")
    opt.add_argument("-o", "--output", metavar="DIR", default=None,
                     help="Output directory. Default: <input_dir>/check/ or ./rapid_check/")
    opt.add_argument("-c", type=int, default=10, metavar="INT",
                     help="Number of best primer pairs to test. Default: 10.")
    opt.add_argument("-min", "--min-size", type=int, default=70, metavar="INT",
                     help="Minimum expected amplicon size (bp). Default: 70.")
    opt.add_argument("-max", "--max-size", type=int, default=250, metavar="INT",
                     help="Maximum expected amplicon size (bp). Default: 250.")
    opt.add_argument("--mismatch", type=int, default=3, metavar="INT",
                     help="Max mismatches allowed by ipcress. Default: 3.")
    opt.add_argument("--maps", action="store_true",
                     help="Generate GBIF biogeography maps after ipcress (requires internet).")
    opt.add_argument("--threads", type=int, default=4, metavar="INT",
                     help="Parallel ipcress jobs (one per database file). Default: 4.")
    opt.add_argument("--snakemake-args", metavar="'...'", default="",
                     help="Extra Snakemake flags (e.g. '--dryrun'). Quote the string.")

    args = parser.parse_args(sys.argv[2:])

    print_banner()
    print(f"\n  Mode       : check")
    print(f"  Input      : {args.input}")
    print(f"  Database   : {args.tsa or args.database}")
    print(f"  Pairs      : top {args.c}")
    print(f"  Amplicon   : {args.min_size}–{args.max_size} bp")
    print(f"  Mismatches : {args.mismatch}")
    print(f"  Maps       : {'yes' if args.maps else 'no'}")
    print()

    check_tools(["ipcress"])

    try:
        tsv_path = find_primers_summary(args.input)
    except FileNotFoundError as e:
        print(f"[RAPID check] ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    print(f"[RAPID check] Primers file : {tsv_path}")

    if args.output:
        out_dir = Path(args.output).resolve()
    elif Path(args.input).is_dir():
        out_dir = Path(args.input).resolve() / "check"
    else:
        out_dir = Path("rapid_check").resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    db_files, use_tsa = resolve_db_files(args)
    print(f"[RAPID check] {len(db_files)} database file(s) to process.")

    config = build_check_config(args, tsv_path, out_dir, db_files, use_tsa)
    config_path = write_yaml_config(config, out_dir, "rapid_check_config.yaml")
    print(f"[RAPID check] Config written to: {config_path}")

    run_snakemake(config_path, args.snakemake_args, args.threads, snakefile_name="check.smk")


# =============================================================================
# MODE: install
# =============================================================================

def run_install():
    parser = argparse.ArgumentParser(
        prog="rapid install",
        description=(
            "Create or update the 'rapid' conda environment.\n\n"
            "Looks for rapid.yaml in the following locations (in order):\n"
            "  1. envs/rapid.yaml  (inside the Rapid project directory)\n"
            "  2. rapid.yaml       (at the project root — flat layout)\n\n"
            "If the 'rapid' environment already exists it will be updated.\n"
            "After installation, activate with:  conda activate rapid"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--conda", default=None, metavar="CONDA_EXE",
                        help="Path to the conda/mamba executable. "
                             "Auto-detected if not specified.")
    parser.add_argument("--update-only", action="store_true",
                        help="Force update mode even if the environment does not exist.")
    args = parser.parse_args(sys.argv[2:])

    print_banner()
    print("\n  Mode : install\n")

    conda_exe = args.conda
    if conda_exe is None:
        for candidate in ("mamba", "conda"):
            found = shutil.which(candidate)
            if found:
                conda_exe = found
                break
    if conda_exe is None:
        print("  ERROR: conda/mamba not found in PATH.\n"
              "  Install Miniconda or Anaconda first, then re-run 'rapid install'.",
              file=sys.stderr)
        sys.exit(1)
    print(f"  Using conda : {conda_exe}")

    base = Path(__file__).parent
    env_file = None
    for candidate in [base / "envs" / "rapid.yaml", base / "rapid.yaml"]:
        if candidate.is_file():
            env_file = candidate
            break
    if env_file is None:
        print("  ERROR: rapid.yaml not found. Expected at:\n"
              f"    {base / 'envs' / 'rapid.yaml'}\n"
              f"    {base / 'rapid.yaml'}",
              file=sys.stderr)
        sys.exit(1)
    print(f"  Environment file : {env_file}\n")

    result = subprocess.run(
        [conda_exe, "env", "list"],
        capture_output=True, text=True
    )
    env_exists = any(
        line.split()[0] == "rapid"
        for line in result.stdout.splitlines()
        if line and not line.startswith("#")
    )

    if env_exists or args.update_only:
        action = "Updating"
        cmd = [conda_exe, "env", "update", "-n", "rapid",
               "-f", str(env_file), "--prune"]
    else:
        action = "Creating"
        cmd = [conda_exe, "env", "create", "-f", str(env_file)]

    print(f"  {action} conda environment 'rapid' …\n")
    ret = subprocess.run(cmd).returncode

    if ret == 0:
        print("\n  ✓ Done!\n")
        print("  Activate with:\n    conda activate rapid\n")
        print("  Then run the pipeline:\n    rapid auto --help\n")
    else:
        print("\n  ✗ Installation failed (see error above).", file=sys.stderr)
        sys.exit(ret)


# =============================================================================
# MODE: clean
# =============================================================================

def run_clean():
    parser = argparse.ArgumentParser(
        prog="rapid clean",
        description=(
            "Clean Snakemake metadata so the pipeline can be re-run cleanly.\n\n"
            "Removes .snakemake/ directories (metadata, locks, incomplete file markers)\n"
            "from the current working directory AND from the Rapid project directory.\n"
            "Works for rapid auto, track and check.\n"
            "Use --outputs to also wipe a results directory."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-o", "--output", metavar="DIR", default=None,
                        help="Output directory to optionally wipe (requires --outputs).")
    parser.add_argument("--outputs", action="store_true",
                        help="Also delete the results directory entirely (irreversible).")
    args = parser.parse_args(sys.argv[2:])

    print_banner()
    print(f"\n  Mode : clean\n")

    cleaned = False

    candidates = list({
        Path.cwd() / ".snakemake",
        Path(__file__).parent / ".snakemake",
    })

    for meta in candidates:
        if meta.exists():
            shutil.rmtree(meta)
            print(f"  ✓ Removed Snakemake metadata : {meta}")
            cleaned = True
        else:
            print(f"  ℹ Not found : {meta}")

    if args.outputs:
        if not args.output:
            print("\n  ERROR: --outputs requires -o <dir>.", file=sys.stderr)
            sys.exit(1)
        out = Path(args.output).resolve()
        if out.exists():
            confirm = input(f"\n  WARNING: delete '{out}' and ALL its contents? [yes/N] ").strip()
            if confirm.lower() == "yes":
                shutil.rmtree(out)
                print(f"  ✓ Removed results directory : {out}")
                cleaned = True
            else:
                print("  Aborted — results directory kept.")
        else:
            print(f"  ℹ Results directory not found: {out}")

    if cleaned:
        print("\n  Done. Re-run your pipeline normally — Snakemake will resume from scratch.\n")
    else:
        print("\n  Nothing to clean.\n")


# =============================================================================
# Entry point
# =============================================================================

def main():
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        print_mode_help()
        sys.exit(0)

    mode = sys.argv[1]

    if mode == "auto":
        run_auto()
    elif mode == "track":
        run_track()
    elif mode == "check":
        run_check()
    elif mode == "install":
        run_install()
    elif mode == "clean":
        run_clean()
    elif mode in ("-v", "--version"):
        print("RAPID v0.1.0 — RNA Automated Primer Identification and Design")
    else:
        print(f"\n[RAPID] Unknown mode: '{mode}'\n", file=sys.stderr)
        print_mode_help()
        sys.exit(1)


if __name__ == "__main__":
    main()
