# =============================================================================
# 04_pathway_enrichment.R
# GO, KEGG, and GSEA enrichment analysis on DE genes from tamoxifen experiment.
#
# Outputs:
#   results/enrichment_results/go_up.csv
#   results/enrichment_results/go_down.csv
#   results/enrichment_results/kegg_results.csv
#   results/enrichment_results/gsea_results.csv
#   results/figures/go_dotplot_up.pdf
#   results/figures/go_dotplot_down.pdf
#   results/figures/kegg_dotplot.pdf
#   results/figures/gsea_plot.pdf
#   results/figures/gsea_hallmarks_barplot.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(msigdbr)
  library(tidyverse)
  library(RColorBrewer)
})

set.seed(42)
dir.create("results/enrichment_results", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures",            recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load DE results
# ---------------------------------------------------------------------------
cat(">>> Loading DE results...\n")

if (file.exists("results/res_df.rds")) {
  res_df <- readRDS("results/res_df.rds")
} else {
  res_df <- read.csv("results/de_genes.csv")
}

up_genes   <- filter(res_df, sig == "Up")
down_genes <- filter(res_df, sig == "Down")

cat("  Up-regulated:   ", nrow(up_genes), "genes\n")
cat("  Down-regulated: ", nrow(down_genes), "genes\n")

# ---------------------------------------------------------------------------
# Gene ID conversion — need Entrez IDs for enrichment
# ---------------------------------------------------------------------------
cat(">>> Converting gene IDs to Entrez...\n")

# Use SYMBOL if available, otherwise ENSEMBL
id_type <- if (any(grepl("^ENSG", res_df$gene_id))) "ENSEMBL" else "SYMBOL"
id_col  <- if (id_type == "ENSEMBL") "gene_id" else "symbol"

convert_to_entrez <- function(gene_vec, from_type = id_type) {
  bitr(gene_vec,
       fromType = from_type,
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) |>
    pull(ENTREZID) |>
    unique()
}

up_entrez   <- convert_to_entrez(up_genes[[id_col]])
down_entrez <- convert_to_entrez(down_genes[[id_col]])
all_entrez  <- convert_to_entrez(res_df[[id_col]])

cat("  Up Entrez IDs:   ", length(up_entrez), "\n")
cat("  Down Entrez IDs: ", length(down_entrez), "\n")

# ---------------------------------------------------------------------------
# GO Enrichment — Biological Process
# ---------------------------------------------------------------------------
cat("\n>>> Running GO enrichment (Biological Process)...\n")

run_go <- function(entrez_ids, label) {
  enrichGO(
    gene          = entrez_ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    readable      = TRUE
  )
}

go_up   <- run_go(up_entrez,   "up")
go_down <- run_go(down_entrez, "down")

# Simplify to remove redundant terms
go_up_simp   <- simplify(go_up,   cutoff = 0.7, by = "p.adjust")
go_down_simp <- simplify(go_down, cutoff = 0.7, by = "p.adjust")

cat("  GO up   — significant terms:", nrow(as.data.frame(go_up_simp)), "\n")
cat("  GO down — significant terms:", nrow(as.data.frame(go_down_simp)), "\n")

# Save
write.csv(as.data.frame(go_up_simp),
          "results/enrichment_results/go_up.csv",   row.names = FALSE)
write.csv(as.data.frame(go_down_simp),
          "results/enrichment_results/go_down.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# GO dot plots
# ---------------------------------------------------------------------------
plot_go_dot <- function(go_obj, title, color_high) {
  dotplot(go_obj, showCategory = 20) +
    scale_color_gradient(low = "grey80", high = color_high) +
    labs(title = title) +
    theme_bw(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 12),
      axis.text.y  = element_text(size = 9)
    )
}

p_go_up <- plot_go_dot(go_up_simp,
                        "GO Biological Process — Up-regulated genes\n(tamoxifen treatment)",
                        "#d73027")
ggsave("results/figures/go_dotplot_up.pdf", p_go_up, width = 10, height = 9)
cat("  Saved: results/figures/go_dotplot_up.pdf\n")

p_go_down <- plot_go_dot(go_down_simp,
                          "GO Biological Process — Down-regulated genes\n(tamoxifen treatment)",
                          "#4575b4")
ggsave("results/figures/go_dotplot_down.pdf", p_go_down, width = 10, height = 9)
cat("  Saved: results/figures/go_dotplot_down.pdf\n")

# ---------------------------------------------------------------------------
# KEGG Pathway enrichment
# ---------------------------------------------------------------------------
cat("\n>>> Running KEGG enrichment...\n")

kegg_up <- enrichKEGG(
  gene          = up_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

kegg_down <- enrichKEGG(
  gene          = down_entrez,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

# Combine for saving
kegg_combined <- bind_rows(
  as.data.frame(kegg_up)   |> mutate(direction = "Up"),
  as.data.frame(kegg_down) |> mutate(direction = "Down")
)
write.csv(kegg_combined, "results/enrichment_results/kegg_results.csv",
          row.names = FALSE)

cat("  KEGG up   — pathways:", nrow(as.data.frame(kegg_up)), "\n")
cat("  KEGG down — pathways:", nrow(as.data.frame(kegg_down)), "\n")

# KEGG dot plot (combined up + down)
p_kegg <- bind_rows(
    as.data.frame(kegg_up)   |> mutate(direction = "Up-regulated"),
    as.data.frame(kegg_down) |> mutate(direction = "Down-regulated")
  ) |>
  mutate(
    GeneRatio_num = sapply(GeneRatio, function(x) {
      parts <- strsplit(x, "/")[[1]]
      as.numeric(parts[1]) / as.numeric(parts[2])
    }),
    Description = str_wrap(Description, 40)
  ) |>
  group_by(direction) |>
  slice_min(p.adjust, n = 10) |>
  ungroup() |>
  ggplot(aes(x = GeneRatio_num, y = reorder(Description, GeneRatio_num),
             color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "#d73027", high = "grey70",
                       name = "Adj. p-value") +
  scale_size_continuous(name = "Gene count", range = c(3, 8)) +
  facet_wrap(~direction, scales = "free_y") +
  labs(
    title = "KEGG pathway enrichment — Tamoxifen vs. Control",
    x     = "Gene ratio",
    y     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title  = element_text(face = "bold"),
    axis.text.y = element_text(size = 9),
    strip.text  = element_text(face = "bold")
  )

ggsave("results/figures/kegg_dotplot.pdf", p_kegg, width = 14, height = 8)
cat("  Saved: results/figures/kegg_dotplot.pdf\n")

# ---------------------------------------------------------------------------
# GSEA — MSigDB Hallmark gene sets
# ---------------------------------------------------------------------------
cat("\n>>> Running GSEA with MSigDB Hallmarks...\n")

# Build ranked gene list (all genes, ranked by sign(LFC) * -log10(pvalue))
ranked_list <- res_df |>
  filter(!is.na(padj), !is.na(log2FoldChange)) |>
  mutate(
    entrez = bitr(
      if (id_type == "ENSEMBL") gene_id else symbol,
      fromType = id_type,
      toType   = "ENTREZID",
      OrgDb    = org.Hs.eg.db
    )$ENTREZID[match(
      if (id_type == "ENSEMBL") gene_id else symbol,
      bitr(if (id_type == "ENSEMBL") gene_id else symbol,
           fromType = id_type, toType = "ENTREZID",
           OrgDb = org.Hs.eg.db)[[id_type]]
    )],
    rank_stat = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300))
  ) |>
  filter(!is.na(entrez), !is.na(rank_stat)) |>
  arrange(desc(rank_stat)) |>
  distinct(entrez, .keep_all = TRUE)

gene_ranks <- setNames(ranked_list$rank_stat, ranked_list$entrez)

# Get Hallmark gene sets
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") |>
  select(gs_name, entrez_gene) |>
  mutate(entrez_gene = as.character(entrez_gene))

gsea_res <- GSEA(
  geneList     = gene_ranks,
  TERM2GENE    = hallmarks,
  pvalueCutoff = 0.05,
  eps          = 0,
  seed         = TRUE
)

gsea_df <- as.data.frame(gsea_res) |>
  mutate(ID = str_replace(ID, "HALLMARK_", "") |> str_replace_all("_", " ") |>
           str_to_title())

write.csv(gsea_df, "results/enrichment_results/gsea_results.csv",
          row.names = FALSE)
cat("  GSEA significant hallmarks:", nrow(gsea_df), "\n")

# GSEA enrichment plot — top activated and suppressed hallmarks
top_activated   <- gsea_df |> filter(NES > 0) |> arrange(pvalue) |> head(3) |> pull(ID)
top_suppressed  <- gsea_df |> filter(NES < 0) |> arrange(pvalue) |> head(3) |> pull(ID)

if (length(c(top_activated, top_suppressed)) > 0) {
  plot_ids <- gsea_res@result$ID[
    gsea_res@result$ID %in% c(
      paste0("HALLMARK_", toupper(str_replace_all(top_activated, " ", "_"))),
      paste0("HALLMARK_", toupper(str_replace_all(top_suppressed, " ", "_")))
    )
  ]

  if (length(plot_ids) > 0) {
    pdf("results/figures/gsea_enrichment_plots.pdf", width = 10, height = 5)
    for (pid in plot_ids[seq_len(min(6, length(plot_ids)))]) {
      print(gseaplot2(gsea_res, geneSetID = pid,
                      title = str_replace(pid, "HALLMARK_", "") |>
                        str_replace_all("_", " ") |> str_to_title()))
    }
    dev.off()
    cat("  Saved: results/figures/gsea_enrichment_plots.pdf\n")
  }
}

# GSEA summary barplot
p_gsea_bar <- gsea_df |>
  arrange(NES) |>
  mutate(
    ID        = factor(ID, levels = ID),
    direction = ifelse(NES > 0, "Activated", "Suppressed")
  ) |>
  ggplot(aes(x = NES, y = ID, fill = direction)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c(Activated = "#d73027", Suppressed = "#4575b4")) +
  labs(
    title    = "GSEA — MSigDB Hallmark gene sets",
    subtitle = "Tamoxifen-treated vs. control MCF-7 cells (padj < 0.05)",
    x        = "Normalized Enrichment Score (NES)",
    y        = NULL,
    fill     = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold"),
    legend.position = "top",
    axis.text.y     = element_text(size = 9)
  )

ggsave("results/figures/gsea_hallmarks_barplot.pdf",
       p_gsea_bar, width = 9, height = max(5, nrow(gsea_df) * 0.35))
cat("  Saved: results/figures/gsea_hallmarks_barplot.pdf\n")

# ---------------------------------------------------------------------------
# Enrichment summary
# ---------------------------------------------------------------------------
cat("\n")
cat(strrep("=", 60), "\n")
cat("ENRICHMENT SUMMARY\n")
cat(strrep("=", 60), "\n")
cat("\n  Top GO BP terms — UP-regulated:\n")
as.data.frame(go_up_simp) |>
  select(Description, Count, p.adjust) |>
  mutate(p.adjust = formatC(p.adjust, format = "e", digits = 2)) |>
  head(5) |> print()

cat("\n  Top GO BP terms — DOWN-regulated:\n")
as.data.frame(go_down_simp) |>
  select(Description, Count, p.adjust) |>
  mutate(p.adjust = formatC(p.adjust, format = "e", digits = 2)) |>
  head(5) |> print()

cat("\n  Top GSEA Hallmarks (by |NES|):\n")
gsea_df |>
  arrange(desc(abs(NES))) |>
  select(ID, NES, p.adjust) |>
  mutate(NES = round(NES, 3),
         p.adjust = formatC(p.adjust, format = "e", digits = 2)) |>
  head(8) |> print()
cat(strrep("=", 60), "\n")
