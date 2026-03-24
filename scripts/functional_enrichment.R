#!/usr/bin/env Rscript

# =============================================================================
# Functional Enrichment Module for LncPro (Snakemake version)
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(Hmisc)
  library(enrichR)
  library(AnnotationDbi)
  library(org.At.tair.db)
})

# -----------------------------------------------------------------------------
# Get inputs/outputs/params from Snakemake
# -----------------------------------------------------------------------------
gene_counts_file    <- snakemake@input[["gene_counts"]]
de_results_file     <- snakemake@input[["deseq2_results"]]
lncrna_gtf_file     <- snakemake@input[["lncrna_gtf"]]

corr_output         <- snakemake@output[["corr_table"]]
go_output           <- snakemake@output[["go_enrich"]]
reactome_output     <- snakemake@output[["reactome_enrich"]]

cor_threshold       <- snakemake@params[["cor_threshold"]]
p_threshold         <- snakemake@params[["p_threshold"]]
databases           <- snakemake@params[["databases"]]

# -----------------------------------------------------------------------------
# 1. Read data
# -----------------------------------------------------------------------------
counts <- read.csv(gene_counts_file, row.names = 1, check.names = FALSE) %>% as.matrix()
de_results <- read.csv(de_results_file, row.names = 1)

padj_cutoff <- 0.05
lfc_cutoff  <- 1

sig_genes <- de_results %>%
  filter(padj <= padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) %>%
  rownames()
cat("Number of significantly DE genes:", length(sig_genes), "\n")

# -----------------------------------------------------------------------------
# 2. Identify lncRNAs from GTF
# -----------------------------------------------------------------------------
if (!file.exists(lncrna_gtf_file)) stop("lncRNA GTF file not found.")
library(stringr)
gtf_lines <- read_lines(lncrna_gtf_file)
gene_ids <- str_extract(gtf_lines, '(?<=gene_id ")[^"]+')
lncrna_genes <- unique(na.omit(gene_ids))
cat("Number of high-confidence lncRNA genes:", length(lncrna_genes), "\n")

sig_lncrna <- intersect(sig_genes, lncrna_genes)
sig_mrna   <- setdiff(sig_genes, lncrna_genes)
cat("Significant DE lncRNAs:", length(sig_lncrna), "\n")
cat("Significant DE mRNAs:", length(sig_mrna), "\n")

if (length(sig_lncrna) < 2 || length(sig_mrna) < 2) {
  stop("Not enough DE lncRNAs or mRNAs for correlation analysis.")
}

# -----------------------------------------------------------------------------
# 3. Subset counts and compute correlations
# -----------------------------------------------------------------------------
common_genes <- intersect(sig_genes, rownames(counts))
counts_sub <- counts[common_genes, , drop = FALSE]

lncrna_counts <- counts_sub[intersect(sig_lncrna, common_genes), , drop = FALSE]
mrna_counts   <- counts_sub[intersect(sig_mrna,   common_genes), , drop = FALSE]

cat("Computing correlations...\n")
cor_results <- rcorr(t(lncrna_counts), t(mrna_counts), type = "pearson")
cor_matrix <- cor_results$r
p_matrix   <- cor_results$P

cor_df <- as.data.frame(as.table(cor_matrix)) %>%
  set_names("lncRNA", "mRNA", "correlation") %>%
  mutate(p_value = as.vector(p_matrix))

significant_pairs <- cor_df %>%
  filter(abs(correlation) >= cor_threshold, p_value <= p_threshold) %>%
  arrange(desc(abs(correlation)))

cat("Number of significant lncRNA-mRNA pairs:", nrow(significant_pairs), "\n")
write_csv(significant_pairs, corr_output)

# -----------------------------------------------------------------------------
# 4. ID conversion: Extract TAIR pattern and map to gene symbols
# -----------------------------------------------------------------------------
mrna_ids <- unique(significant_pairs$mRNA)
cat("Unique mRNAs in significant pairs:", length(mrna_ids), "\n")

if (length(mrna_ids) == 0) {
  cat("No mRNAs to test for enrichment.\n")
  quit(save = "no")
}

# Remove version suffixes (e.g., .1) first
mrna_ids_clean <- str_remove(mrna_ids, "\\..*$")

# Extract the TAIR locus ID using a regex pattern (matches AT1G01010, AT2G12345, etc.)
tair_pattern <- "AT[1-5C]G\\d+"
extracted_tair <- str_extract(mrna_ids_clean, tair_pattern)

# Keep only those that had a match
tair_ids <- extracted_tair[!is.na(extracted_tair)]
cat("TAIR IDs extracted after pattern matching:", length(tair_ids), "\n")

if (length(tair_ids) == 0) {
  stop("No TAIR IDs found among significant mRNAs. Enrichment cannot be performed.")
}

cat("First 10 extracted TAIR IDs:", paste(head(tair_ids, 10), collapse = ", "), "\n")

# Convert using org.At.tair.db
gene_symbols <- mapIds(org.At.tair.db,
                       keys = tair_ids,
                       column = "SYMBOL",
                       keytype = "TAIR",
                       multiVals = "first")
gene_symbols <- gene_symbols[!is.na(gene_symbols)]
cat("Final gene symbols for enrichment:", length(gene_symbols), "\n")

if (length(gene_symbols) == 0) {
  stop("No gene symbols could be converted for enrichment analysis.")
}

# Save conversion for reference
conversion_df <- data.frame(TAIR_id = names(gene_symbols),
                            symbol = gene_symbols,
                            row.names = NULL)
write.csv(conversion_df, file = "tair_to_symbol_conversion.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. Enrichment analysis using enrichR
# -----------------------------------------------------------------------------
cat("Running enrichment with", length(gene_symbols), "gene symbols...\n")
enriched <- enrichr(gene_symbols, databases)

# Write enrichment results
for (db in names(enriched)) {
  df <- enriched[[db]]
  if (nrow(df) > 0) {
    if (db == databases[1]) write.csv(df, go_output, row.names = FALSE)
    if (db == databases[2]) write.csv(df, reactome_output, row.names = FALSE)
  }
}

cat("Enrichment results saved.\n")
