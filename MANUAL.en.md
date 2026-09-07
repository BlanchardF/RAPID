# RAPID — Manual & Tutorial

**R**NA **A**utomated **P**rimer **I**dentification and **D**esign

This manual covers concepts, installation, every mode and option, choosing filter thresholds, the interactive reports, reproducibility, parallel runs, and troubleshooting. For a short overview see `README.en.md`.

---

## 1. Concept

RAPID designs primers whose 3′ end straddles an **exon–exon junction** of a spliced transcript. Such a primer can only anneal on mature mRNA (where the two exons are contiguous), not on genomic DNA (where an intron separates them). This gives RNA-specific assays without relying solely on a DNase step.

The design is driven by RNA-seq: genes are selected by their expression, and — in `track` mode — by their **specificity to one reference stage/condition**.

Pipeline (mode `auto`):

1. **Annotation** — provided (`-a`) or produced de novo by BRAKER3 from the reads.
2. **Alignment** — HISAT2 (splice-aware), then SAMtools sort/index.
3. **Quantification** — featureCounts, per gene, at the exon level.
4. **Gene selection & junctions** — DESeq2 size-factor normalisation, selection by expression stratum, computation of internal exon–exon junctions, export of exon BED + a Primer3 preconfig with the junction directives.
5. **Spliced templates** — BEDTools extracts exons; they are concatenated per gene.
6. **Primer3** — one config per gene, run in parallel (GNU Parallel), with a per-gene timeout.
7. **Ranking & summary** — best pairs by penalty → `best_primers.txt` + `primers_summary.tsv`.

`track` inserts a reference-specificity filter between steps 3 and 4; everything downstream is identical. `check` is a separate validation workflow.

---

## 2. Installation

```bash
rapid install                 # create the 'rapid' env from rapid.yaml
rapid install --update-only   # update an existing env
conda activate rapid
```

The `rapid` environment contains: Snakemake, BRAKER3 + AUGUSTUS, HISAT2, SAMtools, Subread (featureCounts), BEDTools, Primer3, GNU Parallel, exonerate (provides `ipcress`), R with DESeq2/dplyr/tidyr, and folium (biogeography maps).

**AGAT (for GFF3 → GTF conversion) is intentionally NOT in this environment** — it pulls a large Perl stack that conflicts with the pinned dependencies. Install it in its own environment when needed (see §9).

---

## 3. Inputs

**Genome** — FASTA, **decompressed** (HISAT2 and BEDTools do not read `.gz`) and ideally **soft-masked** (repeats in lowercase). Soft-masking is required by BRAKER3 and is transparent to alignment; RAPID upper-cases sequences before Primer3, so lowercase has no downstream effect. Do **not** hard-mask (Ns break alignment and primer design).

```bash
gunzip -k genome.fa.gz        # keep the .gz, produce genome.fa
```

**Annotation** *(optional but recommended)* — **GTF preferred**. RAPID detects the format from the file extension:
- `.gtf` → featureCounts uses `-g gene_id`, and reproducible splice sites are extracted for HISAT2 (see §7).
- `.gff/.gff3` → featureCounts uses `-g Parent` (aggregates per **transcript**, not per gene) and the HISAT2 reproducibility safeguard is skipped.

If your annotation is GFF3, convert it to GTF first (see §9). If you omit `-a`, BRAKER3 runs automatically — but its predictions are stochastic, so for reproducible results run it once and reuse the resulting `braker.gtf` via `-a`.

**RNA-seq reads** — FASTQ(.gz). One file = single-end; two files (R1 R2) = paired-end.

---

## 4. Mode `auto`

Full pipeline on a single RNA-seq sample.

```bash
rapid auto -g genome.fa -r R1.fq.gz R2.fq.gz -o results/ -a annotation.gtf
```

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `-g, --genome` | — | Reference genome (FASTA). **Required.** |
| `-r, --rna` | — | 1 file (SE) or 2 files (PE R1 R2). **Required.** |
| `-o, --output` | — | Output directory. **Required.** |
| `-a, --annot` | — | Annotation (GTF/GFF3). If omitted, BRAKER3 runs. |
| `-ex, --expression` | `max` | Expression stratum: `max`/`+` (top 33%), `mid` (33–66%), `min` (bottom 33%), `N` (top N), `-N` (bottom N). |
| `--primer-product-min` | `100` | Minimum amplicon size (bp). |
| `--primer-product-max` | `200` | Maximum amplicon size (bp). Must be > min. |
| `--top-primers` | `10` | Number of best pairs to keep (by penalty). |
| `-p, --primer3-param` | — | Extra Primer3 global parameter `KEY=VALUE`, repeatable (e.g. `-p PRIMER_MIN_SIZE=18`). |
| `--primer3-jobs` | = `--threads` | Parallel Primer3 processes (one gene each). |
| `--primer3-timeout` | `7200` | Per-gene wall-clock limit (s); a gene exceeding it is excluded. `0` disables. |
| `--threads` | `8` | CPU threads per tool. |
| `--snakemake-args` | — | Extra flags passed to Snakemake, e.g. `'--dryrun'`. Quote the string. |

### Expression modes explained

Genes must be **multi-exonic** (a junction requires ≥ 2 exons). Among those, RAPID ranks by normalised expression and keeps the chosen stratum. Use `-ex N` (e.g. `-ex 50`) to design primers for a fixed number of the most-expressed genes.

---

## 5. Mode `track` — stage/condition-specific primers

`track` takes one **reference** sample (`--ref`) and any number of **other** samples (`--other`, repeatable), typically different stages or treatments of the same species. It keeps only genes specific to the reference, then designs primers exactly as in `auto`.

```bash
rapid track -g genome.fa -o results_track/ -a annotation.gtf \
    --filter cpm --min-cpm-ref 1 --fold-change 10 \
    --ref   cercariae_R1.fq.gz cercariae_R2.fq.gz \
    --other stageA_R1.fq.gz stageA_R2.fq.gz \
    --other stageB_R1.fq.gz stageB_R2.fq.gz
```

Each `--other` accepts 1 (SE) or 2 (PE) files and can be repeated for as many samples as you like. Samples are named `ref`, `other_1`, `other_2`, … internally.

### The two filters

**`--filter absolute`** *(default, raw counts)* — keep a gene if `ref count ≥ --ref-min-count` **and** `every other count ≤ --other-max-count`.
- `--ref-min-count` (default 10), `--other-max-count` (default 0).
- Simple presence/absence, but depth-dependent and blind to enrichment. With real data, `--other-max-count 0` is very strict; you will usually need to raise it.

**`--filter cpm`** *(recommended)* — keep a gene if `CPM_ref ≥ --min-cpm-ref` **and** `CPM_ref ≥ --fold-change × CPM_other` for **every** other sample.
- `--min-cpm-ref` (default 1.0), `--fold-change` (default 10.0).
- CPM = count / (library size / 1e6). This normalises for sequencing depth and captures **enrichment**, which is the right notion for a stage-specific marker.

All other `auto` options (`-ex`, product size, `--top-primers`, Primer3 params, timeout, threads) apply. Expression ranking is computed on the reference sample.

### Choosing thresholds without prior knowledge

Use the helper `track_cpm_scan.py` on the per-sample count tables produced by a run (no need to re-run the pipeline). It reports how many genes you would keep at each fold-change:

```bash
python3 Rapid/scripts/track_cpm_scan.py \
    --ref    results_track/featurecounts/counts_ref.txt \
    --others results_track/featurecounts/counts_other_*.txt \
    --min-cpm-ref 1 --fold-grid 2 3 5 10 20 50
```

Pick the smallest fold-change that yields a workable number of candidates (a few dozen to a few hundred). If you have a known good gene, add `--gene Smp_169190` to see its per-stage CPM and the exact fold-change below which it is kept — a useful sanity check.

### Output specific to `track`

In addition to the standard outputs, `track` writes an **enriched primer report** (`primer3/results/primer_report.tsv` and `.html`). See §6.

---

## 6. Choosing the best primers — the primer report

`track` runs generate `primer_report.tsv` (a table) and `primer_report.html` (a self-contained interactive page — open it in any browser). For each kept primer pair the report gives:

- amplicon (product) size, Primer3 penalty, primer Tm and GC%
- the gene's CPM in **every stage** (reference first)
- the **enrichment** fold-change (reference vs. the most-expressed other stage) and the limiting stage

In the HTML page you can:

- **Sort** by enrichment ↓, amplicon size ↑, penalty ↑, or gene.
- **Filter** by minimum enrichment and amplicon size range.
- **Expand** a row to see a per-stage CPM bar chart (reference highlighted).
- **Discard** unwanted pairs with the ✕ button; the list narrows as you review. Your selection persists across reloads. Use *View hidden* to review/restore, *Restore all* to reset.
- **Download** your kept selection as a TSV.

Recommended workflow: sort by enrichment, then among the top pairs pick the lowest penalty and an amplicon size in your target window; expand to confirm no other stage is expressed.

The stage labels (e.g. `cercariae`, `stageA`) are read from `rapid_track_config.yaml` — they are the folder names of your reads. If reads are not organised in named folders, generic labels (`other_1`, …) are used.

You can regenerate the report at any time without re-running the pipeline:

```bash
python3 Rapid/scripts/primer_report.py results_track/
```

---

## 7. Mode `check` — in-silico validation

`check` runs electronic PCR (ipcress, from exonerate) with your candidate primers against a sequence database, to assess specificity and possible cross-amplification.

```bash
# against an NCBI TSA nucleotide database (per-species transcriptomes)
rapid check -i results_track/ -tsa /path/to/tsa/nucl/ --maps

# against a single custom FASTA database
rapid check -i my_primers.tsv -d sequences.fasta -c 5 -min 80 -max 300
```

### Options

| Option | Default | Meaning |
|--------|---------|---------|
| `-i, --input` | — | A `rapid` output directory (auto-finds `primers_summary.tsv`) or a TSV with forward in col 1, reverse in col 2. **Required.** |
| `-tsa` / `-d, --database` | — | TSA directory **or** a single FASTA(.gz). One is **required**. |
| `-o, --output` | `<input>/check/` | Output directory. |
| `-c` | `10` | Number of top pairs to test. |
| `-min, --min-size` | `70` | Minimum expected amplicon (bp). |
| `-max, --max-size` | `250` | Maximum expected amplicon (bp). |
| `--mismatch` | `3` | Maximum mismatches allowed by ipcress. |
| `--maps` | off | Also generate GBIF biogeography maps (needs internet). |
| `--threads` | `4` | Parallel ipcress jobs (one per database file). |

`check` is resumable: each database file is processed independently and cached, so an interrupted run resumes where it stopped. Re-running later with `--maps` only generates the maps.

### Species report (TSA)

With a TSA database, `check` extracts the species behind each hit and writes:

- `check/species_by_pair.tsv` — a pair × species matrix (amplicon counts).
- `check/species_explorer.html` — an interactive page. Uncheck species you don't care about; a pair turns green once it no longer amplifies any checked species. This makes it easy to see which pairs stay specific once you disregard acceptable off-targets.

---

## 8. Modes `install` and `clean`

```bash
rapid install [--update-only] [--conda /path/to/conda]
rapid clean                       # remove Snakemake metadata (.snakemake/) to re-run cleanly
rapid clean -o results/ --outputs # also delete a results directory (asks for confirmation)
```

Run `rapid clean` **only when no run is active** — it removes `.snakemake/` from the current directory and from the project directory.

---

## 9. GFF3 → GTF conversion (AGAT)

RAPID works best with GTF (reproducible splice sites + per-gene aggregation). If you have a GFF3 (e.g. from WormBase ParaSite), convert it in a **dedicated** environment:

```bash
conda create -n agat -c bioconda -c conda-forge agat
conda activate agat
gunzip -k annotation.gff3.gz
agat_convert_sp_gff2gtf.pl --gff annotation.gff3 -o annotation.gtf
conda deactivate && conda activate rapid
# sanity check: exons must carry gene_id
grep -P "\texon\t" annotation.gtf | head -3
```

If `agat_*` errors with `Can't locate AGAT/AGAT.pm`, your shell is mixing environments (a different Perl is in `PATH`/`PERL5LIB`). Use `conda run -n agat agat_convert_sp_gff2gtf.pl …`, or clear `PERL5LIB` for the command.

---

## 10. Reproducibility

With a **fixed GTF annotation**, RAPID is deterministic — the same inputs give the same primers, even across machines and different `--threads`. Two points:

- **HISAT2** normally reuses splice sites discovered on the fly, which makes multi-threaded runs slightly non-deterministic. RAPID avoids this by extracting the known splice sites from the annotation and passing `--known-splicesite-infile … --no-temp-splicesite`. This is applied for **GTF** annotations only (the extraction script reads GTF). With GFF3, convert to GTF (§9) for reproducible alignment.
- **BRAKER3** (used when no `-a` is given) is stochastic. Run it once, keep `braker.gtf`, and always reuse it via `-a` for reproducible downstream results.

---

## 11. Running several analyses at once

`rapid` launches Snakemake in the **current working directory**, so Snakemake's lock lives in `./.snakemake/`. Two runs started from the **same** directory will clash on that lock even if their `-o` differ. To run in parallel:

```bash
mkdir job1 job2
( cd job1 && rapid auto -g genome.fa -r A_R1.fq.gz A_R2.fq.gz -o results/ -a annot.gtf ) &
( cd job2 && rapid auto -g genome.fa -r B_R1.fq.gz B_R2.fq.gz -o results/ -a annot.gtf ) &
```

On a cluster with a scheduler (e.g. SLURM), each job already runs in its own working directory, so there is no conflict. Give each run a distinct `-o`, size `--threads` to what the scheduler allocates, and do not run `rapid clean` while jobs are active.

---

## 12. Troubleshooting

**`hisat2-build`/`bedtools` fail on the genome** — the genome is still gzipped. Decompress it: `gunzip -k genome.fa.gz`.

**featureCounts: "failed to find the gene identifier attribute … 'gene_id'"** — your GTF has no `gene_id` on exons (common after some GFF3 conversions). Reconvert with AGAT (§9), then re-run after `rapid clean`.

**`track` keeps 0 genes** — the filter is too strict. In `absolute` mode raise `--other-max-count`; better, switch to `--filter cpm` and pick `--fold-change` with `track_cpm_scan.py` (§5).

**A known marker gene does not appear** — trace it: check its counts in each `counts_*.txt`, then whether it is in `counts_ref_specific.txt`, `top_genes/liste_ID_genes_cibles.txt`, the preconfig, and `best_primers.txt`. It is usually excluded because it is expressed in another stage (correct behaviour), because it is mono-exonic (no junction, cannot be targeted), or because its only junctions differ from those of an external reference (different annotation version).

**Empty `SEQUENCE_TEMPLATE` / Primer3 "SEQUENCE_OVERLAP_JUNCTION_LIST beyond end of sequence"** — this was a gene-name parsing issue with underscore/colon gene IDs; the current scripts handle it. Make sure you are running the up-to-date `01_expression_filter_and_primer3_prep.R` and `common.smk`.

**Two runs give slightly different counts** — see §10 (HISAT2 splice-site determinism; use a GTF).

**Results differ from a published primer set for the same gene** — RAPID is indexed on the annotation you provide, and it only designs across junctions. A published intra-exonic primer, or one based on a different annotation version, will not be reproduced. Use the same annotation version to compare, or accept RAPID's junction-based design.

---

## 13. Output reference

```
results/
├── rapid_config.yaml                (auto)  /  rapid_track_config.yaml (track)
├── annotation/annotation.gtf        symlink to the provided annotation
├── hisat2_index/                    HISAT2 index (+ splicesites.txt for GTF)
├── hisat2/                          alignments (SAM)
├── samtools/                        sorted, indexed BAM
├── featurecounts/
│   ├── counts.txt                   (auto)
│   ├── counts_<sample>.txt          (track, one per sample)
│   └── counts_ref_specific.txt      (track, after the filter)
├── top_genes/
│   ├── liste_ID_genes_cibles.txt    selected gene IDs
│   ├── table_multi_exon_genes_cibles.txt
│   ├── coords/exons_coords.bed
│   ├── every_single_exons.fasta
│   └── merge_exons.fasta            spliced transcripts (Primer3 templates)
├── primer3/
│   ├── primer3_preconfig.txt        per-junction Primer3 configs (unfilled)
│   ├── gene_configs/                one filled config per gene
│   └── results/
│       ├── primer3_raw_results.txt
│       ├── best_primers.txt         top pairs (Primer3 format)
│       ├── primers_summary.tsv      one row per pair
│       ├── primer_report.tsv        (track) enriched metrics + per-stage CPM
│       └── primer_report.html       (track) interactive selection page
├── logs/                            per-rule logs
└── check/                           (check mode)
    ├── ipcress_raw/<db>.txt
    ├── ipcress_raw_combined.txt
    ├── per_pair/<pair>.txt, <pair>_species.txt
    ├── species_by_pair.tsv
    ├── species_explorer.html
    └── maps/                        (with --maps)
```

---

## 14. Tools & citation

RAPID orchestrates: BRAKER3, HISAT2, SAMtools, Subread/featureCounts, DESeq2, BEDTools, Primer3, GNU Parallel, exonerate (ipcress), and folium, via Snakemake. Please cite these tools and RAPID itself.

*(Add the RAPID citation / DOI here.)*
