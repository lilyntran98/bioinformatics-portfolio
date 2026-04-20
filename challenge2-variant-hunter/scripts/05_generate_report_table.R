#!/usr/bin/env Rscript
# =============================================================================
# 05_generate_report_table.R
# Generate a polished HTML + TSV summary table of the top candidates
# for inclusion in the clinical report and GitHub README.
#
# Output:
#   results/report_table.html  — embeddable HTML table
#   results/report_table.tsv   — plain-text version
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(knitr)
  library(kableExtra)
})

# ---------------------------------------------------------------------------
# Load scored candidates
# ---------------------------------------------------------------------------
scored_path <- "results/candidates_scored.tsv"

if (file.exists(scored_path)) {
  scored <- read_tsv(scored_path, show_col_types = FALSE)
} else {
  message("results/candidates_scored.tsv not found — using built-in example data.")

  scored <- tibble(
    chrom     = c("chr2",        "chr14",       "chrX"),
    pos       = c(166182386,     23884527,      153296777),
    ref       = c("G",           "C",           "G"),
    alt       = c("A",           "T",           "A"),
    gene      = c("SCN1A",       "FOXG1",       "MECP2"),
    hgvsp     = c("p.Arg1648His","p.Arg461Ter", "p.Arg306Cys"),
    consequence = c("missense_variant","stop_gained","missense_variant"),
    max_af    = c(0.000001, 0.0, 0.000003),
    cadd      = c(28.4, 44.0, 25.6),
    clinvar   = c("Pathogenic","Likely_pathogenic","Pathogenic"),
    acmg_score = c(30L, 26L, 24L),
    acmg_class = c("Pathogenic","Likely Pathogenic","Pathogenic"),
    acmg_criteria = c(
      "PS1, PS2, PM2, PP2, PP3, PP4",
      "PVS1, PS2, PM2, PP4",
      "PS1, PS2, PM2, PP2, PP3, PP4"
    ),
    disease   = c("Dravet syndrome","Congenital Rett syndrome","Rett syndrome"),
    omim_id   = c("607208","613454","312750"),
    phenotype_seizures = c(TRUE, FALSE, TRUE),
    phenotype_id       = c(TRUE, TRUE,  TRUE),
    phenotype_dd       = c(TRUE, TRUE,  TRUE),
  )
}

# ---------------------------------------------------------------------------
# Format for display
# ---------------------------------------------------------------------------
top3 <- scored |>
  slice_head(n = 3) |>
  mutate(
    Rank        = row_number(),
    Variant     = paste0(chrom, ":", formatC(pos, format="d", big.mark=","),
                         " ", ref, ">", alt),
    `Gene`      = gene,
    `Protein change` = hgvsp,
    `Consequence` = str_replace_all(consequence, "_", " "),
    `gnomAD AF` = case_when(
      max_af == 0 ~ "Novel",
      TRUE        ~ formatC(max_af, format = "e", digits = 1)
    ),
    `CADD`      = round(cadd, 1),
    `ClinVar`   = clinvar,
    `ACMG`      = acmg_class,
    `Criteria`  = acmg_criteria,
    `Disease (OMIM)` = paste0(disease, " (", omim_id, ")"),
    `Phenotype` = paste0(
      ifelse(phenotype_seizures, "Sz ✓", "Sz ✗"), " | ",
      ifelse(phenotype_id,       "ID ✓", "ID ✗"), " | ",
      ifelse(phenotype_dd,       "DD ✓", "DD ✗")
    )
  ) |>
  select(Rank, Variant, Gene, `Protein change`, `Consequence`,
         `gnomAD AF`, CADD, ClinVar, ACMG, Criteria, `Disease (OMIM)`, Phenotype)

# ---------------------------------------------------------------------------
# Write TSV
# ---------------------------------------------------------------------------
write_tsv(top3, "results/report_table.tsv")
cat("TSV saved to: results/report_table.tsv\n")

# ---------------------------------------------------------------------------
# Write formatted HTML table
# ---------------------------------------------------------------------------
acmg_bg <- c(
  "Pathogenic"        = "#fce8e8",
  "Likely Pathogenic" = "#fff3cd",
  "VUS"               = "#e8f4f8"
)

html_table <- top3 |>
  kbl(format = "html", escape = FALSE,
      caption = "Top Candidate Variants — Variant Hunter Pipeline") |>
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "bordered"),
    full_width = TRUE,
    font_size  = 13
  ) |>
  column_spec(1, bold = TRUE, width = "3em") |>
  column_spec(9, bold = TRUE,
              background = case_when(
                top3$ACMG == "Pathogenic"        ~ "#fce8e8",
                top3$ACMG == "Likely Pathogenic" ~ "#fff3cd",
                TRUE                             ~ "#e8f4f8"
              )) |>
  row_spec(1, background = "#fff8f8") |>
  add_header_above(c(
    " " = 1,
    "Variant" = 3,
    "Evidence" = 4,
    "Classification" = 2,
    "Clinical" = 2
  ))

writeLines(html_table, "results/report_table.html")
cat("HTML table saved to: results/report_table.html\n")

# ---------------------------------------------------------------------------
# Print to console
# ---------------------------------------------------------------------------
cat("\n")
top3 |>
  select(Rank, Variant, Gene, `Protein change`, `gnomAD AF`, CADD, ACMG) |>
  print(width = Inf)
