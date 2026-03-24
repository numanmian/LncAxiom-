#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ---- 1. Get arguments from Snakemake ----
counts_file   <- snakemake@input[["counts"]]
metadata_file <- snakemake@input[["metadata"]]

csv_out      <- snakemake@output[["csv"]]
volcano_out  <- snakemake@output[["volcano"]]

design_formula <- as.formula(snakemake@params[["design"]])
factor_name    <- snakemake@params[["factor"]]
case_level     <- snakemake@params[["case"]]
control_level  <- snakemake@params[["control"]]
padj_cutoff    <- as.numeric(snakemake@params[["padj_cutoff"]])
lfc_cutoff     <- as.numeric(snakemake@params[["lfc_cutoff"]])
level          <- snakemake@params[["level"]]   # "gene" or "transcript"

# Decide feature column name just for nicer output
feature_col_name <- if (level == "gene") "gene_id" else "tx_id"

# ---- 2. Read counts ----
counts_df <- read_csv(counts_file, show_col_types = FALSE)
feature_ids <- counts_df[[1]]
counts_mat <- as.data.frame(counts_df[, -1])
rownames(counts_mat) <- feature_ids

# ---- 3. Read metadata ----
coldata <- read_tsv(metadata_file, show_col_types = FALSE) %>% as.data.frame()
rownames(coldata) <- coldata$sample
coldata$sample <- NULL

# ---- 4. Ensure matching samples ----
common_samples <- intersect(colnames(counts_mat), rownames(coldata))
if (length(common_samples) == 0) {
    stop("No common samples between counts matrix and metadata!")
}
counts_mat <- counts_mat[, common_samples, drop = FALSE]
coldata <- coldata[common_samples, , drop = FALSE]

# ---- 5. Create DESeq2 object ----
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_mat)),
  colData   = coldata,
  design    = design_formula
)

dds[[factor_name]] <- relevel(dds[[factor_name]], ref = control_level)

# ---- 6. Run DESeq2 ----
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c(factor_name, case_level, control_level)
)

res_df <- as.data.frame(res)
res_df[[feature_col_name]] <- rownames(res_df)

res_df <- res_df %>%
  dplyr::select(all_of(feature_col_name), everything()) %>%
  arrange(padj)

# ---- 7. Write CSV output ----
write_csv(res_df, csv_out)

# ---- 8. Volcano plot ----
plot_df <- res_df %>%
  filter(!is.na(padj), !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(padj),
    regulation = case_when(
      padj < padj_cutoff & log2FoldChange >= lfc_cutoff  ~ "Up",
      padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Down",
      TRUE                                               ~ "NS"
    )
  )

p <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(aes(color = regulation), alpha = 0.6, size = 1.2) +
  geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
  labs(
    title = paste("Volcano plot -", level, "level"),
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw(base_size = 12)

ggsave(
  filename = volcano_out,
  plot = p,
  width = 6,
  height = 5,
  dpi = 300
)
