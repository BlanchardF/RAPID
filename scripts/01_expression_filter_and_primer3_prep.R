## =============================================================================
## RAPID — scripts/01_expression_filter_and_primer3_prep.R
##
## Purpose:
##   1. Normalise with DESeq2 (estimateSizeFactors)
##   2. Select multi-exonic genes according to the expression mode:
##        "max" / "+"  -> top 33%   (most expressed)
##        "mid"        -> 33-66%    (intermediate expression)
##        "min"        -> bottom 33% (least expressed)
##        integer N    -> the top N most expressed genes
##   3. Export the exon BED coordinates
##   4. Compute the exon-exon junctions and write the Primer3 preconfig
##
## Note: the conservation (Diamond) filter has been removed from the pipeline.
##       All genes from the featureCounts table are now used.
##
## Called by Snakemake via the script: directive (access to the snakemake@* objects)
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(DESeq2)
})

log_msg <- function(...) message("[RAPID R] ", ...)

# Separator used to build the Exon_ID (BED name). It MUST be a character that
# never appears in gene identifiers OR chromosome names, because merge_exons
# recovers the gene name by splitting the header on this separator. WormBase/AGAT
# gene_id values contain "_" and ":" (e.g. "gene:Smp_155860"), so "_" cannot be
# used -- we use "|", which merge_exons splits on to recover the full gene name.
EXON_ID_SEP <- "|"

# ── Parameters from Snakemake ─────────────────────────────────────────────────
counts_file       <- snakemake@input$counts
out_gene_list     <- snakemake@output$gene_list
out_gene_table    <- snakemake@output$gene_table
out_bed           <- snakemake@output$bed
out_preconfig     <- snakemake@output$preconfig
ex_mode           <- snakemake@params$ex_mode    # "max", "mid", "min", or integer string
product_min       <- as.integer(snakemake@params$product_min)
product_max       <- as.integer(snakemake@params$product_max)

# Create output directories if needed
for (f in c(out_gene_list, out_gene_table, out_bed, out_preconfig)) {
  dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
}

## ============================================================================
## PART 1 — Normalisation, top-N selection
## ============================================================================

log_msg("Loading data...")

count_table <- read.delim(counts_file, comment.char = "#")

log_msg(nrow(count_table), " genes in the featureCounts table (no filtering).")

if (nrow(count_table) == 0) {
  stop("No gene in the featureCounts table. Check your input files.")
}

# ── DESeq2 normalisation ─────────────────────────────────────────────────────
# Use the 7th column: first sample in the standard featureCounts format
# (columns 1-6 = GeneID, Chr, Start, End, Strand, Length ; col 7 = first BAM)
count_matrix <- count_table[, 7, drop = FALSE]
rownames(count_matrix) <- count_table$Geneid

sample_info <- data.frame(condition = "sample", row.names = colnames(count_matrix))

dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData   = sample_info,
                               design    = ~ 1)
dds <- estimateSizeFactors(dds)
norm_counts <- counts(dds, normalized = TRUE)

df_expr <- data.frame(
  Geneid     = rownames(norm_counts),
  Expression = norm_counts[, 1]
)

# ── Select multi-exonic genes according to the expression mode ───────────────
# First: all multi-exonic genes sorted by decreasing expression
all_multiexon <- df_expr %>%
  arrange(desc(Expression)) %>%
  inner_join(count_table, by = "Geneid") %>%
  filter(grepl(";", Strand))

n_total <- nrow(all_multiexon)
log_msg(n_total, " multi-exonic genes available.")

if (n_total == 0) {
  stop("No multi-exonic gene. Check the annotation / featureCounts.")
}

# Tercile boundaries
low_cut  <- max(1L, floor(n_total * 0.33))   # at least 1 gene
high_cut <- max(2L, floor(n_total * 0.66))   # at least 2 genes

# Determine the mode
# as.integer("max") -> NA with a warning; suppressWarnings avoids noise
ex_num <- suppressWarnings(as.integer(as.character(ex_mode)))

# Safety: if ex_mode is exactly a known keyword, force ex_num to NA so that
# stray characters are never interpreted as an integer
if (is.character(ex_mode) && ex_mode %in% c("max", "+", "mid", "min")) {
  ex_num <- NA_integer_
}

if (!is.na(ex_num) && ex_num > 0) {
  # Positive integer mode: N most expressed genes
  top_candidates <- head(all_multiexon, ex_num)
  mode_label <- paste0("top ", ex_num, " genes (numeric mode)")

} else if (!is.na(ex_num) && ex_num < 0) {
  # Negative integer mode: |N| least expressed genes
  n_bottom <- abs(ex_num)
  top_candidates <- tail(all_multiexon, n_bottom)
  mode_label <- paste0("bottom ", n_bottom, " genes (negative numeric mode)")

} else if (is.na(ex_num) && ex_mode %in% c("max", "+")) {
  # Top 33%: most expressed
  top_candidates <- head(all_multiexon, low_cut)
  mode_label <- paste0("top 33% (", low_cut, " genes, max mode)")

} else if (is.na(ex_num) && ex_mode == "mid") {
  # Middle tercile: 33% to 66%
  top_candidates <- slice(all_multiexon, (low_cut + 1):high_cut)
  mode_label <- paste0("33-66% (", nrow(top_candidates), " genes, mid mode)")

} else if (is.na(ex_num) && ex_mode == "min") {
  # Bottom 33%: least expressed
  top_candidates <- tail(all_multiexon, n_total - high_cut)
  mode_label <- paste0("bottom 33% (", n_total - high_cut, " genes, min mode)")

} else {
  stop(paste0("Unknown expression mode: '", ex_mode,
              "'. Accepted values: max, +, mid, min, or an integer."))
}

log_msg(nrow(top_candidates), " genes selected — ", mode_label, ".")

if (nrow(top_candidates) == 0) {
  stop("No gene selected. Check the -ex mode and your data.")
}

# Retrieve the full featureCounts rows
count_table_clean_multiexon <- count_table %>%
  filter(Geneid %in% top_candidates$Geneid)

# ── Export ────────────────────────────────────────────────────────────────────
write.table(count_table_clean_multiexon$Geneid, out_gene_list,
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(count_table_clean_multiexon, out_gene_table,
            row.names = TRUE, col.names = TRUE, quote = FALSE, sep = "\t")

log_msg("Gene list -> ", out_gene_list)
log_msg("Multi-exon table -> ", out_gene_table)

# ── BED file (0-based coordinates) ───────────────────────────────────────────
# Exon_ID = Geneid | Chr | Start(1-based). The "|" separator lets merge_exons
# recover the full Geneid even when it contains "_" or ":".
exons_bed <- count_table_clean_multiexon %>%
  separate_rows(Chr, Start, End, Strand, sep = ";") %>%
  mutate(
    Start   = as.numeric(Start) - 1,     # BED is 0-based
    End     = as.numeric(End),
    Exon_ID = paste(Geneid, Chr, Start + 1, sep = EXON_ID_SEP)
  ) %>%
  select(Chr, Start, End, Exon_ID, Geneid, Strand)

write.table(exons_bed, out_bed,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

log_msg("BED file (", nrow(exons_bed), " exons) -> ", out_bed)

## ============================================================================
## PART 2 — Exon-exon junction computation + Primer3 preconfig
## ============================================================================

log_msg("Computing exon-exon junctions...")

primer3_data <- count_table_clean_multiexon %>%
  separate_rows(Chr, Start, End, Strand, sep = ";") %>%
  mutate(
    Start  = as.numeric(Start),
    End    = as.numeric(End),
    Length = End - Start + 1
  ) %>%
  group_by(Geneid) %>%
  # Sort by strand: ascending (+) or descending (-)
  arrange(if_else(Strand == "+", Start, -Start), .by_group = TRUE) %>%
  mutate(
    Junction_Point = cumsum(Length),
    Seq_ID         = paste0(Geneid, "_j", Junction_Point)
  ) %>%
  # Exclude the last exon (no junction after it)
  filter(row_number() < n()) %>%
  ungroup()

log_msg(nrow(primer3_data), " exon-exon junctions identified.")

# ── Write the Primer3 preconfig ──────────────────────────────────────────────
# PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4 forces the 3' end of the primer to be
# at least 4 bases INTO the next exon, guaranteeing a true junction overlap
# (the primer cannot anneal to genomic DNA).
file_conn <- file(out_preconfig, "w")
for (i in seq_len(nrow(primer3_data))) {
  cat("SEQUENCE_ID=",                          primer3_data$Seq_ID[i],        "\n", file = file_conn, sep = "")
  cat("SEQUENCE_TEMPLATE=PLACEHOLDER\n",                                             file = file_conn)
  cat("SEQUENCE_OVERLAP_JUNCTION_LIST=",       primer3_data$Junction_Point[i],"\n", file = file_conn, sep = "")
  cat("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4\n",                                 file = file_conn)
  cat("PRIMER_TASK=generic\n",                                                       file = file_conn)
  cat("PRIMER_PICK_LEFT_PRIMER=1\n",                                                 file = file_conn)
  cat("PRIMER_PICK_RIGHT_PRIMER=1\n",                                                file = file_conn)
  cat("PRIMER_PRODUCT_SIZE_RANGE=", product_min, "-", product_max, "\n",             file = file_conn, sep = "")
  cat("=\n",                                                                          file = file_conn)
}
close(file_conn)

log_msg("Primer3 preconfig -> ", out_preconfig)
log_msg("R step finished.")
