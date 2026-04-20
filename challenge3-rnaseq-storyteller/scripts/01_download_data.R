# =============================================================================
# 01_download_data.R
# Download GSE183947 from GEO using GEOquery, extract count matrix,
# and prepare metadata for DESeq2.
#
# Output:
#   data/counts.txt     — raw count matrix (genes × samples)
#   data/metadata.txt   — sample condition table
# =============================================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
})

set.seed(42)
dir.create("data", showWarnings = FALSE)

GEO_ID <- "GSE183947"

# ---------------------------------------------------------------------------
# Download GEO series
# ---------------------------------------------------------------------------
cat(">>> Downloading", GEO_ID, "from GEO...\n")

gse <- getGEO(GEO_ID, GSEMatrix = TRUE, getGPL = FALSE)

# GEO sometimes returns a list — take first element
if (is.list(gse)) gse <- gse[[1]]

# ---------------------------------------------------------------------------
# Extract expression / count matrix
# ---------------------------------------------------------------------------
cat(">>> Extracting count matrix...\n")

expr <- exprs(gse)
sample_info <- pData(phenoData(gse))

cat("  Dimensions:", nrow(expr), "genes ×", ncol(expr), "samples\n")
cat("  Sample IDs:\n")
print(colnames(expr))

# ---------------------------------------------------------------------------
# Build metadata table
# ---------------------------------------------------------------------------
# GEO characteristic fields vary by study — adjust column names if needed
cat("\n>>> Building metadata...\n")
cat("  Available phenotype columns:\n")
print(colnames(sample_info))

# Extract condition from characteristics — common GEO fields
# Adjust `condition_col` to match the actual column in this dataset
condition_col <- grep("treatment|condition|group|characteristics",
                      colnames(sample_info),
                      ignore.case = TRUE, value = TRUE)[1]

if (!is.na(condition_col)) {
  conditions_raw <- sample_info[[condition_col]]
  cat("  Condition column found:", condition_col, "\n")
  cat("  Raw values:", unique(conditions_raw), "\n")

  # Standardize to "control" / "treated"
  conditions <- case_when(
    grepl("control|DMSO|vehicle|ctrl", conditions_raw, ignore.case = TRUE) ~ "control",
    grepl("tamoxifen|treated|treatment|TAM", conditions_raw, ignore.case = TRUE) ~ "treated",
    TRUE ~ conditions_raw
  )
} else {
  # Fallback: manually assign based on known sample order
  # GSM5571001-003 = control, GSM5571004-006 = treated
  cat("  No condition column detected — using positional assignment.\n")
  n_samples <- ncol(expr)
  conditions <- rep(c("control", "treated"), each = n_samples / 2)
}

metadata <- data.frame(
  sample    = colnames(expr),
  condition = conditions,
  row.names = colnames(expr)
)

cat("\n  Metadata summary:\n")
print(metadata)

# Verify balance
if (!all(c("control", "treated") %in% metadata$condition)) {
  warning("Could not detect both 'control' and 'treated' conditions. ",
          "Please manually inspect data/metadata.txt and correct the condition column.")
}

# ---------------------------------------------------------------------------
# If GEO provides normalized data (not raw counts), attempt to get
# supplementary count files instead
# ---------------------------------------------------------------------------
# Check if values look like counts (integers) or normalized (floats)
sample_vals <- expr[1:10, 1]
if (any(sample_vals != round(sample_vals))) {
  cat("\n  NOTE: Expression matrix appears to contain normalized values, not raw counts.\n")
  cat("  Attempting to download supplementary count files...\n")

  supp_files <- getGEOSuppFiles(GEO_ID, makeDirectory = FALSE)
  cat("  Supplementary files available:\n")
  print(rownames(supp_files))

  # Look for count matrix files
  count_file <- rownames(supp_files)[grep("count|raw|matrix",
                                           rownames(supp_files),
                                           ignore.case = TRUE)[1]]

  if (!is.na(count_file) && file.exists(count_file)) {
    cat("  Loading count file:", count_file, "\n")
    if (grepl("\\.gz$", count_file)) {
      expr <- read.table(gzfile(count_file), header = TRUE,
                         row.names = 1, sep = "\t", check.names = FALSE)
    } else {
      expr <- read.table(count_file, header = TRUE,
                         row.names = 1, sep = "\t", check.names = FALSE)
    }
    cat("  Count matrix dimensions:", nrow(expr), "×", ncol(expr), "\n")
  } else {
    cat("  Could not auto-detect count file.\n")
    cat("  Please manually download counts from:\n")
    cat("  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", GEO_ID, "\n")
  }
}

# ---------------------------------------------------------------------------
# Clean up gene IDs
# ---------------------------------------------------------------------------
# Remove version suffix from Ensembl IDs if present (ENSG00000000.12 → ENSG00000000)
if (any(grepl("^ENSG", rownames(expr)))) {
  rownames(expr) <- sub("\\.\\d+$", "", rownames(expr))
  cat("  Stripped Ensembl version suffixes from gene IDs.\n")
}

# Remove rows with all zeros
expr <- expr[rowSums(expr) > 0, ]
cat("  After removing all-zero rows:", nrow(expr), "genes remaining\n")

# Ensure integer counts for DESeq2
expr <- round(expr)
mode(expr) <- "integer"

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
write.table(expr, "data/counts.txt",
            sep = "\t", quote = FALSE, col.names = NA)

write.table(metadata, "data/metadata.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n>>> Done.\n")
cat("    data/counts.txt    —", nrow(expr), "genes ×", ncol(expr), "samples\n")
cat("    data/metadata.txt  —", nrow(metadata), "samples\n")
cat("\nCondition breakdown:\n")
print(table(metadata$condition))
