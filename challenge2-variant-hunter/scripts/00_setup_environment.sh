#!/usr/bin/env bash
# =============================================================================
# 00_setup_environment.sh
# Set up conda environment and install all required tools
# =============================================================================

set -euo pipefail

ENV_NAME="variant-hunter"
echo ">>> Creating conda environment: ${ENV_NAME}"

conda create -y -n "${ENV_NAME}" python=3.10 r-base=4.3

conda run -n "${ENV_NAME}" conda install -y -c bioconda -c conda-forge \
    bcftools \
    tabix \
    htslib \
    ensembl-vep \
    igv \
    samtools

conda run -n "${ENV_NAME}" conda install -y -c conda-forge \
    pandas \
    numpy \
    openpyxl \
    pysam \
    requests

echo ">>> Installing R packages"
conda run -n "${ENV_NAME}" Rscript - <<'EOF'
pkgs <- c(
  "tidyverse",
  "data.table",
  "readxl",
  "openxlsx",
  "shiny",
  "shinydashboard",
  "DT",
  "plotly",
  "ggplot2",
  "scales",
  "RColorBrewer"
)
install.packages(pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "GenomicRanges", "IRanges"), ask = FALSE)

cat("R packages installed successfully.\n")
EOF

echo ""
echo ">>> Setup complete. Activate environment with:"
echo "    conda activate ${ENV_NAME}"
