# =============================================================================
# 05_network_analysis.R
# STRING protein interaction network for top DE genes.
# Identifies hub genes by degree centrality and visualizes subnetworks
# for the most significantly enriched biological processes.
#
# Outputs:
#   results/figures/string_network_up.pdf
#   results/figures/string_network_down.pdf
#   results/figures/hub_genes_barplot.pdf
#   results/enrichment_results/hub_genes.csv
# =============================================================================

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(tidygraph)
  library(tidyverse)
  library(RColorBrewer)
})

set.seed(42)
dir.create("results/figures",            recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment_results", recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Load DE results
# ---------------------------------------------------------------------------
cat(">>> Loading DE results...\n")

if (file.exists("results/res_df.rds")) {
  res_df <- readRDS("results/res_df.rds")
} else {
  res_df <- read.csv("results/de_genes.csv")
}

up_genes   <- filter(res_df, sig == "Up")   |> arrange(padj)
down_genes <- filter(res_df, sig == "Down") |> arrange(padj)

# Use top 100 per direction for network (manageable size)
top_up   <- head(up_genes,   100)
top_down <- head(down_genes, 100)

cat("  Top up genes for network:   ", nrow(top_up), "\n")
cat("  Top down genes for network: ", nrow(top_down), "\n")

# ---------------------------------------------------------------------------
# STRING database setup
# ---------------------------------------------------------------------------
cat("\n>>> Connecting to STRING database (human, v11.5)...\n")
string_db <- STRINGdb$new(
  version    = "11.5",
  species    = 9606,         # Homo sapiens
  score_threshold = 400,     # medium confidence
  network_type = "functional",
  input_directory = "."
)

# ---------------------------------------------------------------------------
# Map genes to STRING IDs
# ---------------------------------------------------------------------------
map_to_string <- function(gene_df, symbol_col = "symbol") {
  mapped <- string_db$map(
    as.data.frame(gene_df),
    symbol_col,
    removeUnmappedRows = TRUE
  )
  return(mapped)
}

up_mapped   <- map_to_string(top_up)
down_mapped <- map_to_string(top_down)

cat("  Mapped up:   ", nrow(up_mapped),   "of", nrow(top_up),   "genes\n")
cat("  Mapped down: ", nrow(down_mapped), "of", nrow(top_down), "genes\n")

# ---------------------------------------------------------------------------
# Build networks and compute hub genes
# ---------------------------------------------------------------------------
build_network <- function(mapped_df, label) {
  cat("\n>>> Building", label, "network...\n")

  # Get interactions from STRING
  interactions <- string_db$get_interactions(mapped_df$STRING_id)

  if (nrow(interactions) == 0) {
    cat("  No interactions found for", label, "\n")
    return(NULL)
  }

  # Build igraph object
  g <- graph_from_data_frame(
    d        = interactions |> select(from, to, combined_score),
    directed = FALSE,
    vertices = mapped_df |>
      select(STRING_id, symbol, log2FoldChange, padj) |>
      rename(name = STRING_id)
  )

  # Compute centrality metrics
  V(g)$degree      <- degree(g)
  V(g)$betweenness <- betweenness(g, normalized = TRUE)
  V(g)$hub_score   <- hub_score(g)$vector

  cat("  Nodes:", vcount(g), "| Edges:", ecount(g), "\n")

  # Hub genes — top 10 by degree
  hub_df <- data.frame(
    symbol          = V(g)$symbol,
    degree          = V(g)$degree,
    betweenness     = round(V(g)$betweenness, 4),
    hub_score       = round(V(g)$hub_score, 4),
    log2FoldChange  = round(V(g)$log2FoldChange, 3),
    direction       = label
  ) |>
    arrange(desc(degree)) |>
    head(20)

  cat("  Top hub genes:\n")
  print(head(hub_df, 10))

  return(list(graph = g, hubs = hub_df))
}

net_up   <- build_network(up_mapped,   "Up-regulated")
net_down <- build_network(down_mapped, "Down-regulated")

# ---------------------------------------------------------------------------
# Visualize networks with ggraph
# ---------------------------------------------------------------------------
plot_network <- function(net_obj, title, node_color_high, output_file) {
  if (is.null(net_obj)) return(invisible(NULL))

  g <- net_obj$graph

  # Keep only the largest connected component for clarity
  comp  <- components(g)
  giant <- induced_subgraph(g, which(comp$membership == which.max(comp$csize)))

  # Limit to top 60 nodes by degree for readability
  if (vcount(giant) > 60) {
    top_nodes <- order(degree(giant), decreasing = TRUE)[1:60]
    giant     <- induced_subgraph(giant, top_nodes)
  }

  tg <- as_tbl_graph(giant) |>
    mutate(
      degree     = centrality_degree(),
      hub_score  = node_attr(giant, "hub_score"),
      lfc        = node_attr(giant, "log2FoldChange")
    )

  p <- ggraph(tg, layout = "fr") +
    geom_edge_link(aes(alpha = combined_score / 1000),
                   color = "grey60", linewidth = 0.4,
                   show.legend = FALSE) +
    geom_node_point(aes(size = degree, color = lfc)) +
    geom_node_text(
      aes(label = symbol,
          filter = degree >= quantile(degree, 0.75)),
      repel    = TRUE,
      size     = 2.8,
      fontface = "bold"
    ) +
    scale_color_gradient2(
      low      = "#4575b4",
      mid      = "white",
      high     = "#d73027",
      midpoint = 0,
      name     = "log2FC"
    ) +
    scale_size_continuous(range = c(2, 10), name = "Degree") +
    scale_edge_alpha_continuous(range = c(0.1, 0.6)) +
    labs(
      title    = title,
      subtitle = paste0("STRING v11.5 · score ≥ 400 · ",
                        vcount(giant), " nodes · ", ecount(giant), " edges")
    ) +
    theme_graph(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "grey50"),
      legend.position = "right"
    )

  ggsave(output_file, p, width = 11, height = 9)
  cat("  Saved:", output_file, "\n")
}

plot_network(net_up,
             "STRING protein network — Up-regulated genes (tamoxifen)",
             "#d73027",
             "results/figures/string_network_up.pdf")

plot_network(net_down,
             "STRING protein network — Down-regulated genes (tamoxifen)",
             "#4575b4",
             "results/figures/string_network_down.pdf")

# ---------------------------------------------------------------------------
# Hub gene barplot (combined up + down)
# ---------------------------------------------------------------------------
cat("\n>>> Plotting hub gene summary...\n")

all_hubs <- bind_rows(
  if (!is.null(net_up))   net_up$hubs   |> head(10),
  if (!is.null(net_down)) net_down$hubs |> head(10)
)

if (nrow(all_hubs) > 0) {
  write.csv(all_hubs, "results/enrichment_results/hub_genes.csv",
            row.names = FALSE)

  p_hubs <- all_hubs |>
    mutate(symbol = reorder(symbol, degree)) |>
    ggplot(aes(x = degree, y = symbol, fill = direction)) +
    geom_bar(stat = "identity", width = 0.7) +
    scale_fill_manual(
      values = c("Up-regulated"   = "#d73027",
                 "Down-regulated" = "#4575b4")
    ) +
    labs(
      title    = "Hub genes by network degree — STRING network",
      subtitle = "Top 10 per direction (tamoxifen vs. control, MCF-7)",
      x        = "Degree (number of interactions)",
      y        = NULL,
      fill     = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      legend.position = "top"
    )

  ggsave("results/figures/hub_genes_barplot.pdf", p_hubs,
         width = 8, height = 7)
  cat("  Saved: results/figures/hub_genes_barplot.pdf\n")
  cat("  Saved: results/enrichment_results/hub_genes.csv\n")
}

# ---------------------------------------------------------------------------
# Final summary
# ---------------------------------------------------------------------------
cat("\n")
cat(strrep("=", 60), "\n")
cat("NETWORK ANALYSIS SUMMARY\n")
cat(strrep("=", 60), "\n")

if (!is.null(net_up)) {
  cat("\n  Up-regulated network top hubs:\n")
  net_up$hubs |> select(symbol, degree, log2FoldChange) |> head(5) |> print()
}
if (!is.null(net_down)) {
  cat("\n  Down-regulated network top hubs:\n")
  net_down$hubs |> select(symbol, degree, log2FoldChange) |> head(5) |> print()
}
cat(strrep("=", 60), "\n")
