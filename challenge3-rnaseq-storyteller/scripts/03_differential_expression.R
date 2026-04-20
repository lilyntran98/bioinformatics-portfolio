# =============================================================================
# 03_differential_expression.R
# DESeq2 differential expression: tamoxifen-treated vs. control MCF-7 cells.
#
# Outputs:
#   results/de_genes.csv          — full DE results table
#   results/de_significant.csv    — filtered (padj<0.05, |LFC|>1)
#   results/figures/volcano_plot.pdf
#   results/figures/ma_plot.pdf
#   results/figures/top50_heatmap.pdf
#   results/figures/top_gene_boxplots.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

set.seed(42)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat(">>> Loading DESeq2 object from QC step...\n")

if (file.exists("results/dds_qc.rds")) {
  dds <- readRDS("results/dds_qc.rds")
} else {
  # Rebuild from scratch if needed
  counts   <- read.table("data/counts.txt",   header = TRUE, row.names = 1, sep = "\t")
  metadata <- read.table("data/metadata.txt", header = TRUE, sep = "\t")
  rownames(metadata) <- metadata$sample
  metadata$condition <- factor(metadata$condition, levels = c("control", "treated"))
  dds <- DESeqDataSetFromMatrix(counts, metadata, ~ condition)
  keep <- rowSums(counts(dds)) >= 10
  dds  <- dds[keep, ]
}

# ---------------------------------------------------------------------------
# Add gene symbols via AnnotationDbi
# ---------------------------------------------------------------------------
cat(">>> Mapping Ensembl IDs to gene symbols...\n")

gene_ids <- rownames(dds)

# Handle both Ensembl IDs and gene symbols
if (any(grepl("^ENSG", gene_ids))) {
  symbols <- mapIds(org.Hs.eg.db,
                    keys    = gene_ids,
                    column  = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")
  entrez  <- mapIds(org.Hs.eg.db,
                    keys    = gene_ids,
                    column  = "ENTREZID",
                    keytype = "ENSEMBL",
                    multiVals = "first")
} else {
  symbols <- gene_ids
  entrez  <- mapIds(org.Hs.eg.db,
                    keys    = gene_ids,
                    column  = "ENTREZID",
                    keytype = "SYMBOL",
                    multiVals = "first")
}

# ---------------------------------------------------------------------------
# Run DESeq2
# ---------------------------------------------------------------------------
cat(">>> Running DESeq2...\n")
dds <- DESeq(dds)

# Extract results with LFC shrinkage (apeglm)
res_shrunk <- lfcShrink(dds,
                         coef    = "condition_treated_vs_control",
                         type    = "apeglm",
                         quiet   = TRUE)

res_df <- as.data.frame(res_shrunk) |>
  rownames_to_column("gene_id") |>
  mutate(
    symbol  = symbols[gene_id],
    entrez  = entrez[gene_id],
    symbol  = coalesce(symbol, gene_id),   # fallback to ID if no symbol
    sig     = case_when(
      padj < 0.05 & log2FoldChange >  1 ~ "Up",
      padj < 0.05 & log2FoldChange < -1 ~ "Down",
      TRUE                               ~ "NS"
    )
  ) |>
  arrange(padj)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
up_genes   <- filter(res_df, sig == "Up")
down_genes <- filter(res_df, sig == "Down")
ns_genes   <- filter(res_df, sig == "NS")

cat("\n")
cat(strrep("=", 55), "\n")
cat("DIFFERENTIAL EXPRESSION SUMMARY\n")
cat(strrep("=", 55), "\n")
cat(sprintf("  Total genes tested:    %d\n", nrow(res_df)))
cat(sprintf("  Significant (padj<0.05, |LFC|>1):\n"))
cat(sprintf("    Up-regulated:        %d\n", nrow(up_genes)))
cat(sprintf("    Down-regulated:      %d\n", nrow(down_genes)))
cat(sprintf("    Total significant:   %d\n", nrow(up_genes) + nrow(down_genes)))

cat("\n  Top 10 up-regulated genes:\n")
up_genes |> select(symbol, log2FoldChange, padj) |>
  mutate(log2FoldChange = round(log2FoldChange, 2),
         padj = formatC(padj, format = "e", digits = 2)) |>
  head(10) |> print()

cat("\n  Top 10 down-regulated genes:\n")
down_genes |> select(symbol, log2FoldChange, padj) |>
  mutate(log2FoldChange = round(log2FoldChange, 2),
         padj = formatC(padj, format = "e", digits = 2)) |>
  head(10) |> print()
cat(strrep("=", 55), "\n")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------
write.csv(res_df, "results/de_genes.csv", row.names = FALSE)
write.csv(
  filter(res_df, sig != "NS"),
  "results/de_significant.csv",
  row.names = FALSE
)
cat("\n  Saved: results/de_genes.csv\n")
cat("  Saved: results/de_significant.csv\n")

# ---------------------------------------------------------------------------
# Plot 1: Volcano plot
# ---------------------------------------------------------------------------
cat("\n>>> Plotting volcano...\n")

# Label top 15 up and top 15 down by padj
top_labels <- bind_rows(
  up_genes   |> head(15),
  down_genes |> head(15)
) |> pull(symbol)

sig_colors <- c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey70")

p_volcano <- ggplot(res_df,
                    aes(x = log2FoldChange, y = -log10(padj),
                        color = sig)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_point(data = filter(res_df, symbol %in% top_labels),
             size = 2, alpha = 0.9) +
  geom_text_repel(
    data          = filter(res_df, symbol %in% top_labels),
    aes(label     = symbol),
    size          = 2.8,
    max.overlaps  = 20,
    box.padding   = 0.3,
    show.legend   = FALSE
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray40", linewidth = 0.4) +
  scale_color_manual(values = sig_colors,
                     labels = c(Up   = paste0("Up (n=", nrow(up_genes), ")"),
                                Down = paste0("Down (n=", nrow(down_genes), ")"),
                                NS   = paste0("NS (n=", nrow(ns_genes), ")"))) +
  labs(
    title    = "Volcano plot — Tamoxifen vs. Control (MCF-7)",
    subtitle = "padj < 0.05 and |log2FC| > 1",
    x        = "log2 Fold Change (treated / control)",
    y        = "-log10(adjusted p-value)",
    color    = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("results/figures/volcano_plot.pdf", p_volcano, width = 8, height = 7)
cat("  Saved: results/figures/volcano_plot.pdf\n")

# ---------------------------------------------------------------------------
# Plot 2: MA plot
# ---------------------------------------------------------------------------
cat(">>> Plotting MA plot...\n")

p_ma <- ggplot(res_df, aes(x = log10(baseMean + 1),
                            y = log2FoldChange,
                            color = sig)) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_hline(yintercept = c(-1, 0, 1),
             linetype   = c("dashed", "solid", "dashed"),
             color      = c("gray40", "black", "gray40"),
             linewidth  = 0.4) +
  scale_color_manual(values = sig_colors) +
  labs(
    title = "MA plot — Tamoxifen vs. Control",
    x     = "log10(mean expression + 1)",
    y     = "log2 Fold Change",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(face = "bold"), legend.position = "top")

ggsave("results/figures/ma_plot.pdf", p_ma, width = 7, height = 5)
cat("  Saved: results/figures/ma_plot.pdf\n")

# ---------------------------------------------------------------------------
# Plot 3: Top 50 DE genes heatmap
# ---------------------------------------------------------------------------
cat(">>> Plotting top-50 heatmap...\n")

top50 <- res_df |>
  filter(sig != "NS") |>
  arrange(padj) |>
  head(50) |>
  pull(gene_id)

# Get VST matrix
vsd <- vst(dds, blind = FALSE)
vsd_mat <- assay(vsd)

# Z-score per gene
heat_mat <- vsd_mat[top50, ]
heat_mat <- t(scale(t(heat_mat)))
rownames(heat_mat) <- res_df$symbol[match(rownames(heat_mat), res_df$gene_id)]

# Sample annotation
ann_df <- data.frame(
  Condition = dds$condition,
  row.names = colnames(dds)
)
ann_colors <- list(
  Condition = c(control = "#4dac26", treated = "#d01c8b")
)

pdf("results/figures/top50_heatmap.pdf", width = 8, height = 12)
pheatmap(
  heat_mat,
  annotation_col  = ann_df,
  annotation_colors = ann_colors,
  color           = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  breaks          = seq(-2.5, 2.5, length.out = 101),
  cluster_rows    = TRUE,
  cluster_cols    = TRUE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  fontsize_row    = 8,
  main            = "Top 50 DE genes — z-scored VST counts",
  border_color    = NA
)
dev.off()
cat("  Saved: results/figures/top50_heatmap.pdf\n")

# ---------------------------------------------------------------------------
# Plot 4: Box plots for top 6 key genes
# ---------------------------------------------------------------------------
cat(">>> Plotting key gene box plots...\n")

key_genes_symbols <- c("GREB1", "TFF1", "MKI67", "CCNB1", "ESR1", "BCL2")

# Match to gene_ids
key_ids <- res_df |>
  filter(symbol %in% key_genes_symbols) |>
  pull(gene_id)

norm_counts <- counts(dds, normalized = TRUE)

box_data <- norm_counts[key_ids, , drop = FALSE] |>
  as.data.frame() |>
  rownames_to_column("gene_id") |>
  left_join(select(res_df, gene_id, symbol), by = "gene_id") |>
  pivot_longer(cols = -c(gene_id, symbol),
               names_to = "sample",
               values_to = "norm_count") |>
  left_join(as.data.frame(colData(dds)) |>
              rownames_to_column("sample") |>
              select(sample, condition),
            by = "sample")

p_box <- ggplot(box_data,
                aes(x = condition, y = log2(norm_count + 1),
                    fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.7) +
  scale_fill_manual(values = c(control = "#4dac26", treated = "#d01c8b")) +
  facet_wrap(~symbol, scales = "free_y", ncol = 3) +
  labs(
    title = "Key gene expression — tamoxifen vs. control",
    x     = NULL, y = "log2(normalized counts + 1)",
    fill  = "Condition"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )

ggsave("results/figures/key_gene_boxplots.pdf", p_box, width = 9, height = 6)
cat("  Saved: results/figures/key_gene_boxplots.pdf\n")

# Save dds and results for downstream scripts
saveRDS(dds,    "results/dds_final.rds")
saveRDS(res_df, "results/res_df.rds")
cat("\n  dds_final.rds and res_df.rds saved to results/\n")
