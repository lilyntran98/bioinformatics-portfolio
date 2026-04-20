# =============================================================================
# 02_qc_analysis.R
# Quality control and exploratory analysis of the count matrix.
#
# Outputs (saved to results/figures/):
#   library_sizes.pdf        — bar chart of per-sample read depth
#   pca_plot.pdf             — PCA colored by condition
#   correlation_heatmap.pdf  — sample-to-sample Pearson correlation
#   boxplot_normalized.pdf   — normalized count distributions
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(tidyverse)
})

set.seed(42)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
cat(">>> Loading count matrix and metadata...\n")

counts   <- read.table("data/counts.txt",   header = TRUE, row.names = 1,
                       sep = "\t", check.names = FALSE)
metadata <- read.table("data/metadata.txt", header = TRUE, sep = "\t")
rownames(metadata) <- metadata$sample

# Align sample order
counts <- counts[, metadata$sample]

cat("  Counts:", nrow(counts), "genes ×", ncol(counts), "samples\n")
cat("  Conditions:", table(metadata$condition), "\n")

# Condition as factor — control is reference
metadata$condition <- factor(metadata$condition,
                              levels = c("control", "treated"))

# ---------------------------------------------------------------------------
# Build DESeq2 object (for normalization only at this step)
# ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ condition
)

# Filter: keep genes with at least 10 counts across all samples
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat("  After low-count filter:", nrow(dds), "genes retained\n")

# Normalize
dds <- estimateSizeFactors(dds)

# Variance-stabilizing transformation for visualization
vsd <- vst(dds, blind = TRUE)

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------
condition_colors <- c(control = "#4dac26", treated = "#d01c8b")
sample_colors    <- condition_colors[metadata$condition]
names(sample_colors) <- metadata$sample

# ---------------------------------------------------------------------------
# Plot 1: Library sizes
# ---------------------------------------------------------------------------
cat(">>> Plotting library sizes...\n")

lib_sizes <- data.frame(
  sample    = colnames(counts),
  total_counts = colSums(counts),
  condition = metadata$condition
)

p_lib <- ggplot(lib_sizes, aes(x = reorder(sample, -total_counts),
                                y = total_counts / 1e6,
                                fill = condition)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_hline(yintercept = 10, linetype = "dashed",
             color = "gray40", linewidth = 0.5) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title    = "Library size per sample",
    subtitle = "Dashed line = 10M reads (minimum recommended)",
    x        = "Sample",
    y        = "Total mapped reads (millions)",
    fill     = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("results/figures/library_sizes.pdf", p_lib,
       width = 7, height = 5)
cat("  Saved: results/figures/library_sizes.pdf\n")

# ---------------------------------------------------------------------------
# Plot 2: PCA
# ---------------------------------------------------------------------------
cat(">>> Plotting PCA...\n")

pca_data  <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var   <- round(100 * attr(pca_data, "percentVar"), 1)

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2,
                               color = condition,
                               label = name)) +
  geom_point(size = 5, alpha = 0.9) +
  ggrepel::geom_text_repel(size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = condition_colors) +
  labs(
    title    = "PCA — VST-normalized counts",
    subtitle = paste0("PC1 explains ", pct_var[1],
                      "% of variance; PC2 explains ", pct_var[2], "%"),
    x        = paste0("PC1 (", pct_var[1], "%)"),
    y        = paste0("PC2 (", pct_var[2], "%)"),
    color    = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title  = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("results/figures/pca_plot.pdf", p_pca, width = 7, height = 6)
cat("  Saved: results/figures/pca_plot.pdf\n")

# ---------------------------------------------------------------------------
# Plot 3: Sample-to-sample correlation heatmap
# ---------------------------------------------------------------------------
cat(">>> Plotting sample correlation heatmap...\n")

vsd_mat <- assay(vsd)
cor_mat  <- cor(vsd_mat, method = "pearson")

# Annotation for heatmap columns/rows
ann_df <- data.frame(
  condition = metadata$condition,
  row.names = metadata$sample
)
ann_colors <- list(condition = condition_colors)

pdf("results/figures/correlation_heatmap.pdf", width = 7, height = 6)
pheatmap(
  cor_mat,
  annotation_col  = ann_df,
  annotation_row  = ann_df,
  annotation_colors = ann_colors,
  color           = colorRampPalette(c("#f7f7f7", "#4393c3", "#053061"))(100),
  breaks          = seq(0.85, 1.0, length.out = 101),
  display_numbers = TRUE,
  number_format   = "%.3f",
  fontsize_number = 9,
  main            = "Sample-to-sample Pearson correlation (VST counts)",
  border_color    = "white",
  treeheight_row  = 30,
  treeheight_col  = 30
)
dev.off()
cat("  Saved: results/figures/correlation_heatmap.pdf\n")

# ---------------------------------------------------------------------------
# Plot 4: Normalized count distributions
# ---------------------------------------------------------------------------
cat(">>> Plotting count distributions...\n")

norm_long <- as.data.frame(log2(counts(dds, normalized = TRUE) + 1)) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample", values_to = "log2_norm_count") |>
  left_join(metadata |> select(sample, condition), by = "sample")

p_box <- ggplot(norm_long, aes(x = sample, y = log2_norm_count,
                                fill = condition)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  scale_fill_manual(values = condition_colors) +
  labs(
    title = "Normalized count distributions",
    x     = "Sample",
    y     = "log2(normalized counts + 1)",
    fill  = "Condition"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("results/figures/boxplot_normalized.pdf", p_box,
       width = 8, height = 5)
cat("  Saved: results/figures/boxplot_normalized.pdf\n")

# ---------------------------------------------------------------------------
# QC summary
# ---------------------------------------------------------------------------
cat("\n")
cat(strrep("=", 55), "\n")
cat("QC SUMMARY\n")
cat(strrep("=", 55), "\n")
cat(sprintf("  Samples:               %d\n", ncol(counts)))
cat(sprintf("  Conditions:            control=%d, treated=%d\n",
            sum(metadata$condition == "control"),
            sum(metadata$condition == "treated")))
cat(sprintf("  Genes (raw):           %d\n", nrow(counts)))
cat(sprintf("  Genes (after filter):  %d\n", nrow(dds)))
cat(sprintf("  Median library size:   %.1fM reads\n",
            median(colSums(counts)) / 1e6))

# Check sample clustering — PC1 should separate conditions
pca_check <- pca_data |>
  group_by(condition) |>
  summarise(mean_PC1 = mean(PC1), .groups = "drop")
pc1_sep <- abs(diff(pca_check$mean_PC1))
cat(sprintf("  PC1 condition separation: %.1f units\n", pc1_sep))
if (pc1_sep > 20) {
  cat("  ✓ Samples cluster clearly by condition on PC1\n")
} else {
  cat("  ⚠ Weak PC1 separation — check for batch effects\n")
}

# Check correlations
diag(cor_mat) <- NA
within_ctrl  <- mean(cor_mat[metadata$condition == "control",
                              metadata$condition == "control"], na.rm = TRUE)
within_treat <- mean(cor_mat[metadata$condition == "treated",
                              metadata$condition == "treated"], na.rm = TRUE)
between      <- mean(cor_mat[metadata$condition == "control",
                              metadata$condition == "treated"], na.rm = TRUE)
cat(sprintf("  Within-control correlation:   %.4f\n", within_ctrl))
cat(sprintf("  Within-treated correlation:   %.4f\n", within_treat))
cat(sprintf("  Between-condition correlation: %.4f\n", between))
cat(strrep("=", 55), "\n")

# Save dds for next script
saveRDS(dds, "results/dds_qc.rds")
cat("\n  dds object saved to results/dds_qc.rds\n")
