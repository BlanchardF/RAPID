# =============================================================================
# RAPID — common.smk
# Shared downstream rules, from gene selection through primer design and the
# final summary table. Included by BOTH:
#   - Snakefile   
#   - track.smk  
# =============================================================================

PRIMER3_EXTRA   = config.get("primer3_extra_params", [])   # list of "KEY=VALUE"
PRIMER3_JOBS    = config.get("primer3_jobs", THREADS)       # parallel primer3 processes
PRIMER3_TIMEOUT = config.get("primer3_timeout", 7200)       # per-gene timeout (s), 0 = off


# =============================================================================
# DESeq2 normalization, top-N gene selection, BED + Primer3 preconfig
# =============================================================================
rule select_genes_and_prepare_primer3:
    input:
        counts = COUNTS,
    output:
        gene_list  = f"{OUT}/top_genes/liste_ID_genes_cibles.txt",
        gene_table = f"{OUT}/top_genes/table_multi_exon_genes_cibles.txt",
        bed        = f"{OUT}/top_genes/coords/exons_coords.bed",
        preconfig  = f"{OUT}/primer3/primer3_preconfig.txt",
    log:   f"{OUT}/logs/select_genes.log"
    params:
        ex_mode     = EX_MODE,
        product_min = PRODUCT_MIN,
        product_max = PRODUCT_MAX,
    script: R_SCRIPT


# =============================================================================
# Extract individual exon sequences with BEDTools
# =============================================================================
rule bedtools_getfasta:
    input:
        bed    = f"{OUT}/top_genes/coords/exons_coords.bed",
        genome = GENOME,
    output: fasta = f"{OUT}/top_genes/every_single_exons.fasta"
    log:    f"{OUT}/logs/bedtools_getfasta.log"
    shell:
        """
        bedtools getfasta \\
            -fi {input.genome} -bed {input.bed} \\
            -name -s -fo {output.fasta} &> {log}
        """


# =============================================================================
# Merge exons per gene into full spliced transcript sequences
# =============================================================================
rule merge_exons:
    input:  single_exons = f"{OUT}/top_genes/every_single_exons.fasta"
    output:
        merged   = f"{OUT}/top_genes/merge_exons.fasta",
        combined = f"{OUT}/top_genes/brut_sequences_for_primer3.fasta",
    log:   f"{OUT}/logs/merge_exons.log"
    shell:
        r"""
        awk '
        /^>/ {{
            split(substr($0,2), a, "|");
            gene = a[1];
            if (gene != current) {{
                if (current != "") printf "\n";
                current = gene;
                printf ">%s\n", gene;
            }}
            next;
        }}
        {{ printf "%s", toupper($0) }}
        END {{ if (current != "") printf "\n" }}
        ' {input.single_exons} > {output.merged} 2> {log}

        # Also uppercase the individual exon file before combining
        awk '/^>/ {{print}} !/^>/ {{print toupper($0)}}' \
            {input.single_exons} > {input.single_exons}.upper 2>> {log}
        cat {input.single_exons}.upper {output.merged} > {output.combined} 2>> {log}
        rm -f {input.single_exons}.upper
        """


# =============================================================================
# Fill PLACEHOLDER sequences AND split into one config file per gene
# =============================================================================
rule fill_primer3_config:
    input:
        preconfig = f"{OUT}/primer3/primer3_preconfig.txt",
        merged    = f"{OUT}/top_genes/merge_exons.fasta",
    output: config_dir = directory(f"{OUT}/primer3/gene_configs")
    log:    f"{OUT}/logs/fill_primer3_config.log"
    params:
        script      = str(SCRIPTS_DIR / "fill_primer3_config.py"),
        extra_flags = " ".join(
            f'"{p}"' for p in PRIMER3_EXTRA
        ) if PRIMER3_EXTRA else "",
    shell:
        """
        python3 {params.script} \\
            {input.preconfig} {input.merged} {output.config_dir} \\
            {params.extra_flags} &> {log}
        """


# =============================================================================
# Run Primer3 (one process per gene, in parallel via GNU Parallel)
# =============================================================================
# A per-gene wall-clock timeout (PRIMER3_TIMEOUT seconds, default 7200 = 2h)
# excludes any pathological gene so it can never stall the whole run.
rule primer3:
    input:  config_dir = f"{OUT}/primer3/gene_configs"
    output: raw        = f"{OUT}/primer3/results/primer3_raw_results.txt"
    log:    f"{OUT}/logs/primer3.log"
    threads: PRIMER3_JOBS
    params:
        script  = str(SCRIPTS_DIR / "run_primer3_parallel.py"),
        n_jobs  = PRIMER3_JOBS,
        timeout = PRIMER3_TIMEOUT,
    shell:
        """
        mkdir -p {OUT}/primer3/results
        python3 {params.script} \\
            {input.config_dir} {output.raw} {params.n_jobs} {params.timeout} &> {log}
        """


# =============================================================================
# Select top N primer pairs ranked by penalty
# =============================================================================
rule primer3_top_n:
    input:  raw  = f"{OUT}/primer3/results/primer3_raw_results.txt"
    output: best = f"{OUT}/primer3/results/best_primers.txt"
    log:    f"{OUT}/logs/primer3_top_n.log"
    params: top_n = config.get("top_primers", 10)
    shell:
        """
        echo "Extracting group-0 blocks and ranking by penalty..." > {log}

        grep -E "(^SEQUENCE_ID=|^SEQUENCE_TEMPLATE=|_0_|PRIMER_LEFT_0=|PRIMER_RIGHT_0=|^=$)" {input.raw} \\
            > {output.best}.tmp 2>> {log}

        awk -v top_n={params.top_n} '
        {{
            buffer = buffer $0 "\\n"
            if ($0 ~ /PRIMER_PAIR_0_PENALTY=/) {{
                split($0, part, "=")
                penalty = part[2] + 0
            }}
            if ($0 == "=") {{
                if (penalty != "") {{
                    count++
                    blocks[count]    = buffer
                    penalties[count] = penalty
                }}
                buffer  = ""
                penalty = ""
            }}
        }}
        END {{
            for (i = 1; i <= count; i++) {{
                for (j = i + 1; j <= count; j++) {{
                    if (penalties[i] > penalties[j]) {{
                        tmp_p = penalties[i]; penalties[i] = penalties[j]; penalties[j] = tmp_p
                        tmp_b = blocks[i];    blocks[i]    = blocks[j];    blocks[j]    = tmp_b
                    }}
                }}
            }}
            limit = (top_n < count) ? top_n : count
            for (i = 1; i <= limit; i++) {{
                printf "%s", blocks[i]
            }}
        }}
        ' {output.best}.tmp > {output.best} 2>> {log}

        rm -f {output.best}.tmp
        total=$(grep -c "^SEQUENCE_ID=" {output.best} 2>/dev/null || echo 0)
        echo "Done — ${{total}} pair(s) written (top {params.top_n} requested)." >> {log}
        """


# =============================================================================
# Summary table of best primer pairs
# =============================================================================
rule primer3_summary:
    input:  best   = f"{OUT}/primer3/results/best_primers.txt"
    output: tsv    = f"{OUT}/primer3/results/primers_summary.tsv"
    log:    f"{OUT}/logs/primer3_summary.log"
    params:
        script = str(SCRIPTS_DIR / "primer3_summary.py"),
    shell:
        """
        python3 {params.script} {input.best} {output.tsv} &> {log}
        """
