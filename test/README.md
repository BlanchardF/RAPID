# Synthetic dataset for RAPID (`auto`, `track` and `check` modes)

Fake genome, annotation, RNA-seq reads and TSA databases, built so that the
pipeline runs from start to finish (no real sequences).

## Contents

```
synthetic_data/
├── synthetic_genome.fasta        # 2 chromosomes (~134 kb), soft-masked
├── synthetic_annot.gtf           # 30 genes (3-4 exons), 1 transcript/gene, no isoforms
├── Juvenil_RNAseq/R1_test.fastq.gz, R2_test.fastq.gz   # paired-end 2x100
├── Adult_RNAseq/R1_test.fastq.gz,   R2_test.fastq.gz
├── Tsa_synthetic/                # 2 fake TSA files for `rapid check`
│   ├── tsa.GKQB.1.fsa_nt.gz      #   Vicia sativa
│   └── tsa.GLTR.1.fsa_nt.gz      #   Lupinus albus
├── expected_genes.tsv            # ground truth (category of each gene)
├── make_synthetic_dataset.py     # regenerates genome/GTF/reads (--seed 42)
└── make_synthetic_tsa.py         # regenerates the TSA files from primers_summary.tsv
```

| Category | Genes | Juvenil | Adult |
|---|---|---|---|
| juvenil_specific | 7 | high | zero |
| adult_specific | 7 | zero | high |
| shared | 8 | same | same |
| shared_diff | 4 | x3 | baseline |
| silent | 4 | zero | zero |

## Commands (run from `Rapid/test`, with this folder copied to `test/synthetic_data/`)

```bash
cd path/to/Rapid/test

rapid auto \
    -g synthetic_data/synthetic_genome.fasta \
    -r synthetic_data/Adult_RNAseq/R1_test.fastq.gz synthetic_data/Adult_RNAseq/R2_test.fastq.gz \
    -a synthetic_data/synthetic_annot.gtf \
    -o results_auto/ \
    --probe \
    --threads 4

rapid track \
    -g synthetic_data/synthetic_genome.fasta \
    -a synthetic_data/synthetic_annot.gtf \
    -o results_track/ \
    --filter cpm --min-cpm-ref 1 --fold-change 10 \
    --ref   synthetic_data/Juvenil_RNAseq/R1_test.fastq.gz synthetic_data/Juvenil_RNAseq/R2_test.fastq.gz \
    --other synthetic_data/Adult_RNAseq/R1_test.fastq.gz   synthetic_data/Adult_RNAseq/R2_test.fastq.gz \
    --threads 4
```

### In silico validation (`check`)

```bash
rapid check -i results_track/ -tsa synthetic_data/Tsa_synthetic/ -c 10 --maps --threads 2
rapid check -i results_auto/  -tsa synthetic_data/Tsa_synthetic/ -c 10 --maps --threads 2
```

(`--maps` needs an internet connection; it was tested without it.)

## Expected results (tested)

- `auto`: ~99.98% alignment rate, 10 primer pairs with probe in `primers_summary.tsv`.
- `track`: `counts_ref_specific.txt` contains exactly the 7 `juvenil_specific` genes
  (Gene006, Gene007, Gene017, Gene020, Gene023, Gene028, Gene029).
- `primer3` is the slow part (a few minutes depending on the machine).
- `check`: each TSA file = 150 random transcripts + "homologs" of the genes
  (6% divergence, 0 to 3 mismatches in the primer sites). With the 10 best pairs:
  - `results_auto`: Gene028 (x2 in *Vicia*) and Gene029 amplify in both species;
    Gene004, Gene006 and Gene011 only in *Vicia sativa*; Gene008 and Gene018 only in
    *Lupinus albus*; Gene002 amplifies nowhere.
  - `results_track`: Gene020, Gene028 and Gene029 in both species, Gene006 only in *Vicia*.
- The same primers are produced at every run (Primer3 is deterministic). If you change
  the parameters, or apply the minus-strand fix to RAPID (which changes the primers of
  '-' strand genes), re-run `make_synthetic_tsa.py` with your own `primers_summary.tsv`.
