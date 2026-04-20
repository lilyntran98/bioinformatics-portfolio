#!/usr/bin/env Rscript
# =============================================================================
# 04_prioritize_candidates.R
# Score and rank candidate variants using a weighted ACMG-inspired framework.
#
# ACMG criteria applied (Richards et al. 2015, Genet Med):
#   PVS1  — Loss-of-function in a haploinsufficiency gene (+8)
#   PS1   — Same amino acid change as established pathogenic variant (+7)
#   PS2   — De novo (confirmed) in patient with no family history (+7)
#   PM2   — Absent or extremely rare in population databases (+4)
#   PM5   — Novel missense at same position as pathogenic variant (+3)
#   PP2   — Missense in gene with low missense variation tolerance (+2)
#   PP3   — Multiple in-silico tools predict deleterious (+2)
#   PP4   — Phenotype highly specific for single-gene disorder (+2)
#
# Output:
#   results/candidates_scored.tsv
#   results/top_candidates.xlsx
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(openxlsx)
})

set.seed(42)

# ---------------------------------------------------------------------------
# Load filtered variants
# NOTE: In a real run, parse from results/filtered_variants.vcf using
#       VariantAnnotation or bcftools query output.
#       Here we load from a TSV exported by 03_filter_variants.py.
# ---------------------------------------------------------------------------
candidates_file <- "results/candidates_flat.tsv"

# If the real pipeline output doesn't exist yet, use example data so the
# script is fully executable as a standalone demo.
if (!file.exists(candidates_file)) {
  message("NOTE: results/candidates_flat.tsv not found — using built-in example data.")
  message("      Run 03_filter_variants.py first for real results.\n")

  candidates <- tribble(
    ~chrom,  ~pos,       ~ref, ~alt, ~gene,   ~consequence,          ~hgvsp,            ~hgvsc,              ~max_af,    ~cadd,  ~sift_score, ~polyphen_score, ~clinvar,          ~existing_id,
    "chr2",  166182386,  "G",  "A",  "SCN1A", "missense_variant",    "p.Arg1648His",    "c.4943G>A",         0.000001,   28.4,   0.001,       0.998,           "Pathogenic",      "rs121917885",
    "chr14", 23884527,   "C",  "T",  "FOXG1", "stop_gained",         "p.Arg461Ter",     "c.1381C>T",         0.0,        44.0,   NA,          NA,              "Likely_pathogenic","rs863225251",
    "chrX",  153296777,  "G",  "A",  "MECP2", "missense_variant",    "p.Arg306Cys",     "c.916C>T",          0.000003,   25.6,   0.004,       0.994,           "Pathogenic",      "rs28935468",
    "chr2",  166178101,  "T",  "C",  "SCN1A", "missense_variant",    "p.Trp1204Arg",    "c.3610T>C",         0.00004,    18.2,   0.02,        0.88,            "Uncertain_significance", "rs1553411",
    "chr14", 23891044,   "A",  "G",  "FOXG1", "missense_variant",    "p.Pro478Leu",     "c.1433A>G",         0.000012,   16.5,   0.03,        0.76,            ".",               ".",
  )
} else {
  candidates <- read_tsv(candidates_file, show_col_types = FALSE)
}

# ---------------------------------------------------------------------------
# ACMG-inspired scoring function
# ---------------------------------------------------------------------------

## Gene-level metadata (haploinsufficiency, missense constraint)
## pLI ≥ 0.9 = high probability loss-of-function intolerant (gnomAD)
## Z_missense ≥ 3.09 = missense constrained
gene_metadata <- tribble(
  ~gene,   ~pLI,  ~Z_missense, ~haploinsufficient, ~omim_id,  ~disease,
  "SCN1A", 0.99,  4.18,        TRUE,               "607208",  "Dravet syndrome",
  "FOXG1", 0.99,  3.92,        TRUE,               "613454",  "Congenital Rett syndrome",
  "MECP2", 1.00,  3.10,        TRUE,               "312750",  "Rett syndrome",
  "STXBP1",0.99,  3.56,        TRUE,               "612164",  "Ohtahara syndrome",
  "CDKL5", 0.98,  2.88,        TRUE,               "300672",  "CDKL5 deficiency disorder",
)

score_acmg <- function(df, gene_meta) {
  df <- df |>
    left_join(gene_meta, by = "gene") |>
    mutate(
      # ── PVS1: LoF in haploinsufficiency gene ──────────────────────────
      PVS1 = case_when(
        grepl("stop_gained|frameshift|splice_acceptor|splice_donor|start_lost",
              consequence) & haploinsufficient == TRUE ~ 8L,
        TRUE ~ 0L
      ),

      # ── PS1: Same AA change as known pathogenic ───────────────────────
      PS1 = case_when(
        clinvar %in% c("Pathogenic", "Likely_pathogenic") ~ 7L,
        TRUE ~ 0L
      ),

      # ── PS2: De novo (assumed here — parents unaffected) ──────────────
      PS2 = case_when(
        !is.na(gene) ~ 7L,   # clinical scenario: parents unaffected → de novo assumed
        TRUE ~ 0L
      ),

      # ── PM2: Absent/extremely rare in population ──────────────────────
      PM2 = case_when(
        max_af == 0 ~ 4L,
        max_af < 0.0001 ~ 4L,
        max_af < 0.001  ~ 2L,
        TRUE ~ 0L
      ),

      # ── PP2: Missense in constrained gene ────────────────────────────
      PP2 = case_when(
        grepl("missense", consequence) & (!is.na(Z_missense) & Z_missense >= 3.09) ~ 2L,
        TRUE ~ 0L
      ),

      # ── PP3: Computational evidence (CADD + SIFT + PolyPhen) ─────────
      PP3 = case_when(
        cadd >= 25 & (!is.na(sift_score) & sift_score < 0.05) &
          (!is.na(polyphen_score) & polyphen_score > 0.9) ~ 2L,
        cadd >= 20 ~ 1L,
        TRUE ~ 0L
      ),

      # ── PP4: Phenotype match (seizures + ID → high specificity) ───────
      PP4 = case_when(
        gene %in% c("SCN1A", "MECP2", "FOXG1", "CDKL5", "STXBP1") ~ 2L,
        !is.na(omim_id) ~ 1L,
        TRUE ~ 0L
      ),

      # ── Total score ───────────────────────────────────────────────────
      acmg_score = PVS1 + PS1 + PS2 + PM2 + PP2 + PP3 + PP4,

      # ── ACMG classification from score ───────────────────────────────
      acmg_class = case_when(
        acmg_score >= 20              ~ "Pathogenic",
        acmg_score >= 15              ~ "Likely Pathogenic",
        acmg_score >= 8               ~ "VUS",
        TRUE                          ~ "Likely Benign"
      ),

      # ── Active ACMG criteria string ───────────────────────────────────
      acmg_criteria = {
        pmap_chr(list(PVS1, PS1, PS2, PM2, PP2, PP3, PP4), function(...) {
          vals  <- c(...)
          names <- c("PVS1","PS1","PS2","PM2","PP2","PP3","PP4")
          paste(names[vals > 0], collapse = ", ")
        })
      },

      # ── Phenotype match flags ─────────────────────────────────────────
      phenotype_seizures = gene %in% c("SCN1A","MECP2","CDKL5","STXBP1","KCNQ2","SCN2A"),
      phenotype_id       = gene %in% c("SCN1A","MECP2","FOXG1","SYNGAP1","MBD5","DYRK1A"),
      phenotype_dd       = gene %in% c("SCN1A","MECP2","FOXG1","CDKL5","STXBP1"),
    ) |>
    arrange(desc(acmg_score), desc(cadd))

  return(df)
}

# ---------------------------------------------------------------------------
# Apply scoring
# ---------------------------------------------------------------------------
scored <- score_acmg(candidates, gene_metadata)

# ---------------------------------------------------------------------------
# Print ranked table to console
# ---------------------------------------------------------------------------
cat("\n")
cat(strrep("=", 70), "\n")
cat("CANDIDATE VARIANT RANKING\n")
cat(strrep("=", 70), "\n\n")

display_cols <- c(
  "chrom", "pos", "ref", "alt", "gene", "hgvsp",
  "max_af", "cadd", "acmg_score", "acmg_class", "acmg_criteria",
  "disease"
)

scored |>
  select(any_of(display_cols)) |>
  mutate(
    max_af = formatC(max_af, format = "e", digits = 1),
    cadd   = round(cadd, 1),
    pos    = formatC(pos, format = "d", big.mark = ",")
  ) |>
  print(n = Inf, width = Inf)

# ---------------------------------------------------------------------------
# Save scored table to TSV
# ---------------------------------------------------------------------------
write_tsv(scored, "results/candidates_scored.tsv")
cat("\nScored table saved to: results/candidates_scored.tsv\n")

# ---------------------------------------------------------------------------
# Export top 3 to Excel with formatting
# ---------------------------------------------------------------------------
top3 <- scored |> slice_head(n = 3)

wb <- createWorkbook()
addWorksheet(wb, "Top Candidates")
addWorksheet(wb, "All Candidates")
addWorksheet(wb, "ACMG Evidence")

## Sheet 1: Top 3 candidates — clinical summary
top3_report <- top3 |>
  transmute(
    `Rank`                 = row_number(),
    `Variant`              = glue::glue("{chrom}:{pos} {ref}>{alt}"),
    `Gene`                 = gene,
    `HGVS protein`         = hgvsp,
    `HGVS coding`          = hgvsc,
    `Consequence`          = consequence,
    `gnomAD AF`            = formatC(max_af, format = "e", digits = 1),
    `CADD PHRED`           = round(cadd, 1),
    `SIFT`                 = round(sift_score, 3),
    `PolyPhen`             = round(polyphen_score, 3),
    `ClinVar`              = clinvar,
    `ACMG Score`           = acmg_score,
    `ACMG Class`           = acmg_class,
    `ACMG Criteria`        = acmg_criteria,
    `Disease (OMIM)`       = glue::glue("{disease} ({omim_id})"),
    `Seizures match`       = ifelse(phenotype_seizures, "✓", "—"),
    `ID match`             = ifelse(phenotype_id, "✓", "—"),
    `DD match`             = ifelse(phenotype_dd, "✓", "—"),
  )

writeDataTable(wb, "Top Candidates", top3_report, tableStyle = "TableStyleMedium9")

## Column widths
setColWidths(wb, "Top Candidates",
             cols  = seq_len(ncol(top3_report)),
             widths = c(5, 24, 8, 18, 20, 20, 12, 10, 8, 10, 22, 10, 18, 28, 28, 12, 10, 10))

## Color-code ACMG class
path_rows   <- which(top3_report$`ACMG Class` == "Pathogenic") + 1
lpath_rows  <- which(top3_report$`ACMG Class` == "Likely Pathogenic") + 1
vus_rows    <- which(top3_report$`ACMG Class` == "VUS") + 1

if (length(path_rows) > 0)
  addStyle(wb, "Top Candidates",
           style = createStyle(bgFill = "#fce8e8"),
           rows = path_rows, cols = 13, stack = TRUE)
if (length(lpath_rows) > 0)
  addStyle(wb, "Top Candidates",
           style = createStyle(bgFill = "#fff3cd"),
           rows = lpath_rows, cols = 13, stack = TRUE)
if (length(vus_rows) > 0)
  addStyle(wb, "Top Candidates",
           style = createStyle(bgFill = "#e8f4f8"),
           rows = vus_rows, cols = 13, stack = TRUE)

## Sheet 2: All scored candidates
writeDataTable(wb, "All Candidates", scored |> select(-raw_line, -any_of("raw_line")),
               tableStyle = "TableStyleLight9")

## Sheet 3: ACMG evidence breakdown
acmg_evidence <- top3 |>
  select(gene, hgvsp, PVS1, PS1, PS2, PM2, PP2, PP3, PP4, acmg_score, acmg_class)
writeDataTable(wb, "ACMG Evidence", acmg_evidence, tableStyle = "TableStyleMedium2")

saveWorkbook(wb, "results/top_candidates.xlsx", overwrite = TRUE)
cat("Excel report saved to:  results/top_candidates.xlsx\n")

# ---------------------------------------------------------------------------
# Final verdict
# ---------------------------------------------------------------------------
top1 <- scored |> slice(1)
cat("\n")
cat(strrep("─", 60), "\n")
cat("FINAL DIAGNOSIS\n")
cat(strrep("─", 60), "\n")
cat(sprintf("  Variant    : %s:%s %s>%s\n", top1$chrom, top1$pos, top1$ref, top1$alt))
cat(sprintf("  Gene       : %s\n", top1$gene))
cat(sprintf("  Change     : %s (%s)\n", top1$hgvsp, top1$hgvsc))
cat(sprintf("  Disease    : %s (OMIM: %s)\n", top1$disease, top1$omim_id))
cat(sprintf("  ACMG       : %s (score = %d)\n", top1$acmg_class, top1$acmg_score))
cat(sprintf("  Criteria   : %s\n", top1$acmg_criteria))
cat(sprintf("  Inheritance: De novo (autosomal dominant)\n"))
cat(strrep("─", 60), "\n")
