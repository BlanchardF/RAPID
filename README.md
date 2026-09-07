# RAPID

**R**NA **A**utomated **P**rimer **I**dentification and **D**esign


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


RAPID designs exon–exon **junction-spanning** PCR primers directly from a genome and RNA-seq data. Because each primer sits across a spliced junction, it amplifies mature mRNA but not the corresponding genomic (intron-containing) DNA — ideal for RNA-based assays (RT-PCR, RT-qPCR, ddPCR, eRNA).

The pipeline is built on [Snakemake](https://snakemake.readthedocs.io) and ships as a single command-line tool, `rapid`, with all dependencies pinned in one conda environment.

---

## What it does

Starting from a reference genome + RNA-seq reads, RAPID aligns the reads, quantifies expression, selects genes, computes their exon–exon junctions, reconstructs the spliced transcripts, and runs Primer3 to design junction-spanning primer pairs — then ranks and reports the best candidates.

## Modes

| Mode | Purpose |
|------|---------|
| `auto`    | Full pipeline from one RNA-seq sample → ranked junction-spanning primers |
| `track`   | Like `auto`, but keeps only genes **specific to a reference** stage/condition (one reference sample vs. several others) |
| `check`   | Validates candidate primers by in-silico PCR (ipcress) against a sequence database, with an optional species/biogeography report |
| `install` | Creates or updates the `rapid` conda environment |
| `clean`   | Clears Snakemake metadata for a fresh re-run |

## Requirements

- Linux, [conda](https://docs.conda.io) or [mamba](https://mamba.readthedocs.io)
- A reference genome in **FASTA** (decompressed, ideally **soft-masked**)
- RNA-seq reads (FASTQ, single- or paired-end)
- Optionally, a structural annotation. **GTF is strongly recommended** (see the manual): GFF3 works but disables the reproducibility safeguard and aggregates per transcript.

## Installation

```bash
git clone <your-repo-url> Rapid
cd Rapid
rapid install          # or: python rapid.py install
conda activate rapid
```

`rapid install` looks for `rapid.yaml` and creates the `rapid` environment (HISAT2, SAMtools, Subread/featureCounts, BEDTools, Primer3, GNU Parallel, exonerate/ipcress, BRAKER3, R + DESeq2, …).

## Quick start

```bash
# 1) Full pipeline, using a provided annotation (recommended)
rapid auto -g genome.fa -r reads_R1.fq.gz reads_R2.fq.gz -o results/ -a annotation.gtf

# 2) Stage-specific primers (reference vs. others), depth-normalised enrichment
rapid track -g genome.fa -o results_track/ -a annotation.gtf \
    --filter cpm --min-cpm-ref 1 --fold-change 10 \
    --ref  ref_R1.fq.gz ref_R2.fq.gz \
    --other otherA_R1.fq.gz otherA_R2.fq.gz \
    --other otherB_R1.fq.gz otherB_R2.fq.gz

# 3) In-silico validation of the designed primers against a TSA database
rapid check -i results_track/ -tsa /path/to/tsa/nucl/ --maps
```

Run `rapid <mode> --help` for the full option list of each mode.

## Key outputs

- `primer3/results/best_primers.txt` — top primer pairs (Primer3 format)
- `primer3/results/primers_summary.tsv` — one row per pair (sequences, Tm, GC%, product size, penalty)
- `primer3/results/primer_report.tsv` and `primer_report.html` *(track mode)* — primers enriched with per-stage expression (CPM), enrichment fold-change and amplicon size, plus an **interactive explorer** to sort, filter and shortlist candidates
- `check/species_by_pair.tsv` and `species_explorer.html` *(check mode, TSA)* — which species each pair amplifies, with an interactive filter

## Documentation

See **MANUAL.en.md** for the full tutorial: concepts, every option, choosing thresholds, the interactive reports, reproducibility notes, running several analyses in parallel, and troubleshooting.

## License / citation

*(Add your license and citation here.)*
