# =============================================================================
# RAPID — check.smk
# Validate primer candidates with ipcress, split results per pair,
# extract per-pair taxonomy, build a cross-pair species report, and
# (optionally) generate biogeography maps.
#
# Snakemake-backed for the same reason as the main pipeline: ipcress against
# a large database (e.g. TSA) can take a long time. If it is interrupted or
# the connection drops, simply re-run the same command — completed database
# files are cached and skipped. Re-running later with --maps only triggers
# map generation, without re-running ipcress.
# =============================================================================

from pathlib import Path

OUT             = config["output_dir"]
INPUT_TSV       = config["input_tsv"]
TOP_C           = config.get("top_c", 10)
MIN_SIZE        = config.get("min_size", 70)
MAX_SIZE        = config.get("max_size", 250)
MISMATCH        = config.get("mismatch", 3)
DB_FILES        = config["db_files"]            # dict {sanitized_name: path}
EXTRACT_SPECIES = config.get("extract_species", False)
GENERATE_MAPS   = config.get("generate_maps", False)

SCRIPTS_DIR = Path(workflow.basedir)
DB_NAMES    = list(DB_FILES.keys())


def _final_targets():
    # split is always produced; maps only when requested.
    targets = [f"{OUT}/per_pair/.split_done"]
    # The cross-pair species report is only meaningful for TSA-style databases,
    # where species are parsed from the ipcress "TSA: Genus species" lines.
    if EXTRACT_SPECIES:
        targets += [
            f"{OUT}/species_by_pair.tsv",
            f"{OUT}/species_explorer.html",
        ]
    if GENERATE_MAPS:
        targets.append(f"{OUT}/maps/.maps_done")
    return targets


rule all:
    input:
        _final_targets()


# =============================================================================
# STEP 1 — Parse primers TSV and write the ipcress input file
# =============================================================================
rule prepare_ipcress_input:
    input:  tsv = INPUT_TSV
    output: ipcress_in = f"{OUT}/primers.ipcress"
    log:    f"{OUT}/logs/prepare_ipcress_input.log"
    params:
        script   = str(SCRIPTS_DIR / "prepare_ipcress_input.py"),
        top_c    = TOP_C,
        min_size = MIN_SIZE,
        max_size = MAX_SIZE,
    shell:
        """
        mkdir -p {OUT}/logs
        python3 {params.script} {input.tsv} {output.ipcress_in} \\
            {params.top_c} {params.min_size} {params.max_size} &> {log}
        cat {log}
        """


# =============================================================================
# STEP 2 — Run ipcress against each database file independently
# =============================================================================
# One rule instance per database file (wildcard {db}). This means:
#   - completed database files are cached and skipped on re-run
#   - a crash partway through a large TSA only loses the in-progress file
rule ipcress_per_db:
    input:
        ipcress_in = f"{OUT}/primers.ipcress",
        db         = lambda wc: DB_FILES[wc.db],
    output: raw = f"{OUT}/ipcress_raw/{{db}}.txt"
    log:    f"{OUT}/logs/ipcress_{{db}}.log"
    params: mismatch = MISMATCH
    shell:
        """
        mkdir -p {OUT}/ipcress_raw {OUT}/tmp
        DBFILE="{input.db}"
        case "$DBFILE" in
            *.gz)
                TMP="{OUT}/tmp/{wildcards.db}.fasta"
                zcat "$DBFILE" > "$TMP" 2> {log}
                ipcress --mismatch {params.mismatch} {input.ipcress_in} "$TMP" \\
                    > {output.raw} 2>> {log}
                rm -f "$TMP"
                ;;
            *)
                ipcress --mismatch {params.mismatch} {input.ipcress_in} "$DBFILE" \\
                    > {output.raw} 2> {log}
                ;;
        esac
        """


# =============================================================================
# STEP 3 — Combine all per-database results
# =============================================================================
rule combine_ipcress:
    input:  expand(f"{OUT}/ipcress_raw/{{db}}.txt", db=DB_NAMES)
    output: combined = f"{OUT}/ipcress_raw_combined.txt"
    shell:  "cat {input} > {output.combined}"


# =============================================================================
# STEP 4 — Split combined results by primer pair + extract taxonomy
# =============================================================================
rule split_by_pair:
    input:  combined = f"{OUT}/ipcress_raw_combined.txt"
    output: marker    = touch(f"{OUT}/per_pair/.split_done")
    log:    f"{OUT}/logs/split_by_pair.log"
    params:
        per_pair_dir    = f"{OUT}/per_pair",
        extract_species = "true" if EXTRACT_SPECIES else "false",
        script          = str(SCRIPTS_DIR / "split_ipcress_by_pair.py"),
    shell:
        """
        python3 {params.script} {input.combined} {params.per_pair_dir} \\
            {params.extract_species} &> {log}
        cat {log}
        """


# =============================================================================
# STEP 5 — Cross-pair species report (matrix TSV + interactive HTML explorer)
# =============================================================================
# Reads the per-pair directory produced above and builds:
#   - species_by_pair.tsv    : pair × species matrix (amplicon counts)
#   - species_explorer.html  : self-contained page to filter species out and
#                              see which primer pairs remain specific.
# Only meaningful for TSA-style databases (species parsed from TSA lines), so
# it is wired into the targets only when extract_species is set. The rule is
# always defined, so it can also be requested manually if needed.
rule check_species_report:
    input:  marker = f"{OUT}/per_pair/.split_done"
    output:
        tsv  = f"{OUT}/species_by_pair.tsv",
        html = f"{OUT}/species_explorer.html",
    log:    f"{OUT}/logs/check_species_report.log"
    params:
        per_pair_dir = f"{OUT}/per_pair",
        script       = str(SCRIPTS_DIR / "check_species_report.py"),
    shell:
        """
        python3 {params.script} {params.per_pair_dir} {output.tsv} {output.html} &> {log}
        cat {log}
        """


# =============================================================================
# STEP 6 — Biogeography maps (only when --maps is requested)
# =============================================================================
if GENERATE_MAPS:
    rule generate_maps:
        input:  marker = f"{OUT}/per_pair/.split_done"
        output: marker = touch(f"{OUT}/maps/.maps_done")
        log:    f"{OUT}/logs/generate_maps.log"
        params:
            per_pair_dir = f"{OUT}/per_pair",
            maps_dir     = f"{OUT}/maps",
            script       = str(SCRIPTS_DIR / "biogeography.py"),
        shell:
            """
            python3 {params.script} {params.per_pair_dir} {params.maps_dir} &> {log}
            cat {log}
            """
