# Shiny Variant Browser — Setup & Usage

## Overview

An interactive R Shiny app to explore, filter, and present candidate variants from the Variant Hunter pipeline. Designed for the clinical review step after running the filtering and scoring scripts.

## Features

| Tab | What it shows |
|---|---|
| **Candidate Browser** | Ranked table with ACMG color-coding; click rows for full variant details |
| **ACMG Evidence** | Criteria heatmap, CADD vs AF scatter, phenotype match chart |
| **Filtering Funnel** | Log-scale bar chart showing variant reduction at each step |
| **Clinical Report** | Auto-generated clinical interpretation for the top candidate |

## Requirements

```r
install.packages(c(
  "shiny", "shinydashboard", "DT", "ggplot2",
  "plotly", "dplyr", "tidyr", "scales", "readr"
))
```

## Running the App

### Option A — From the repo root

```r
shiny::runApp("shiny_app")
```

### Option B — From inside shiny_app/

```bash
cd shiny_app
Rscript -e "shiny::runApp('.')"
```

### Option C — RStudio

Open `shiny_app/app.R` in RStudio, then click **Run App**.

## Loading Your Own Data

By default the app uses built-in demo data (5 candidate variants).

To load your real pipeline results, make sure `results/candidates_scored.tsv` exists (produced by `scripts/04_prioritize_candidates.R`). The app detects this file automatically:

```
challenge2-variant-hunter/
├── results/
│   └── candidates_scored.tsv   ← app reads this if present
└── shiny_app/
    └── app.R
```

## Screenshots

> _Add IGV screenshots and app screenshots to `results/igv_screenshots/` after running the full pipeline._

## Deploying to shinyapps.io (optional)

```r
install.packages("rsconnect")
library(rsconnect)
rsconnect::deployApp("shiny_app")
```
