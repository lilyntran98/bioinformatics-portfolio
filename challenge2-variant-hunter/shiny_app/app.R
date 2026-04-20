# =============================================================================
# shiny_app/app.R
# Interactive Variant Browser — GBBA BioHack 2026 Challenge 2
#
# Features:
#   - Filter candidates by CADD, gnomAD AF, consequence, ACMG class
#   - Dynamic ranking table with color-coded pathogenicity
#   - ACMG criteria breakdown per variant
#   - Phenotype match visualization
#   - CSV export
#
# Run: cd shiny_app && Rscript -e "shiny::runApp('.')"
# Or:  shiny::runApp('shiny_app')
# =============================================================================

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(DT)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(tidyr)
  library(scales)
})

# ---------------------------------------------------------------------------
# Built-in demo data (mirrors output of 04_prioritize_candidates.R)
# Replace by loading results/candidates_scored.tsv if it exists
# ---------------------------------------------------------------------------
load_candidates <- function() {
  tsv_path <- "../results/candidates_scored.tsv"
  if (file.exists(tsv_path)) {
    return(readr::read_tsv(tsv_path, show_col_types = FALSE))
  }

  # Demo dataset — 5 candidates
  tibble(
    chrom     = c("chr2",        "chr14",       "chrX",        "chr2",        "chr14"),
    pos       = c(166182386,     23884527,      153296777,     166178101,     23891044),
    ref       = c("G",           "C",           "G",           "T",           "A"),
    alt       = c("A",           "T",           "A",           "C",           "G"),
    gene      = c("SCN1A",       "FOXG1",       "MECP2",       "SCN1A",       "FOXG1"),
    hgvsp     = c("p.Arg1648His","p.Arg461Ter", "p.Arg306Cys", "p.Trp1204Arg","p.Pro478Leu"),
    hgvsc     = c("c.4943G>A",   "c.1381C>T",   "c.916C>T",    "c.3610T>C",   "c.1433A>G"),
    consequence = c("missense_variant","stop_gained","missense_variant","missense_variant","missense_variant"),
    max_af    = c(0.000001,      0.0,           0.000003,      0.00004,       0.000012),
    cadd      = c(28.4,          44.0,          25.6,          18.2,          16.5),
    sift_score     = c(0.001, NA,    0.004,  0.020,  0.030),
    polyphen_score = c(0.998, NA,    0.994,  0.880,  0.760),
    clinvar   = c("Pathogenic",  "Likely_pathogenic","Pathogenic","Uncertain_significance","."),
    acmg_score = c(30L, 26L, 24L, 11L, 8L),
    acmg_class = c("Pathogenic","Likely Pathogenic","Pathogenic","VUS","VUS"),
    acmg_criteria = c(
      "PS1, PS2, PM2, PP2, PP3, PP4",
      "PVS1, PS2, PM2, PP4",
      "PS1, PS2, PM2, PP2, PP3, PP4",
      "PS2, PM2, PP3",
      "PS2, PM2"
    ),
    disease   = c("Dravet syndrome","Congenital Rett syndrome","Rett syndrome","Dravet syndrome","Congenital Rett syndrome"),
    omim_id   = c("607208","613454","312750","607208","613454"),
    pLI       = c(0.99, 0.99, 1.00, 0.99, 0.99),
    phenotype_seizures = c(TRUE,  FALSE, TRUE,  TRUE,  FALSE),
    phenotype_id       = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
    phenotype_dd       = c(TRUE,  TRUE,  TRUE,  TRUE,  TRUE),
  )
}

candidates <- load_candidates()

# ---------------------------------------------------------------------------
# Color palette for ACMG classes
# ---------------------------------------------------------------------------
acmg_colors <- c(
  "Pathogenic"        = "#d73027",
  "Likely Pathogenic" = "#f46d43",
  "VUS"               = "#74add1",
  "Likely Benign"     = "#4dac26",
  "Benign"            = "#1a9641"
)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(
    title = "🧬 Variant Hunter"
  ),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Candidate Browser",  tabName = "browser",  icon = icon("table")),
      menuItem("ACMG Evidence",      tabName = "acmg",     icon = icon("check-circle")),
      menuItem("Filtering Funnel",   tabName = "funnel",   icon = icon("filter")),
      menuItem("Clinical Report",    tabName = "report",   icon = icon("file-medical"))
    ),

    hr(),
    h5("  Filters", style = "padding-left:15px; color:#aaa"),

    sliderInput("cadd_min", "Min CADD score",
                min = 0, max = 50, value = 10, step = 1),

    sliderInput("af_max", "Max gnomAD AF (×10⁻⁴)",
                min = 0, max = 10, value = 10, step = 0.5),

    checkboxGroupInput("consequence_filter", "Consequence",
                       choices = c("missense_variant", "stop_gained",
                                   "frameshift_variant", "splice_donor_variant",
                                   "splice_acceptor_variant"),
                       selected = c("missense_variant", "stop_gained")),

    checkboxGroupInput("acmg_filter", "ACMG class",
                       choices = c("Pathogenic", "Likely Pathogenic", "VUS"),
                       selected = c("Pathogenic", "Likely Pathogenic", "VUS")),

    br(),
    downloadButton("download_csv", "Export CSV", class = "btn-sm btn-success",
                   style = "margin-left:15px")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f6f8; }
      .box { border-radius: 8px; }
      .rank-badge {
        display: inline-block; width: 28px; height: 28px;
        line-height: 28px; text-align: center; border-radius: 50%;
        font-weight: bold; font-size: 13px; color: white;
      }
      .rank-1 { background: #d73027; }
      .rank-2 { background: #f46d43; }
      .rank-3 { background: #74add1; }
    "))),

    tabItems(

      # ── Tab 1: Candidate Browser ──────────────────────────────────────
      tabItem(tabName = "browser",
        fluidRow(
          valueBoxOutput("n_total",   width = 3),
          valueBoxOutput("n_pathogenic", width = 3),
          valueBoxOutput("top_gene",  width = 3),
          valueBoxOutput("top_cadd",  width = 3)
        ),
        fluidRow(
          box(
            title = "Ranked Candidates", width = 12, solidHeader = TRUE,
            status = "primary",
            DTOutput("candidates_table"),
            br(),
            p("Click a row to see detailed ACMG evidence below.", style = "color:#888; font-size:12px")
          )
        ),
        fluidRow(
          box(
            title = "Selected Variant — Detail", width = 12, solidHeader = TRUE,
            status = "info", collapsible = TRUE,
            uiOutput("selected_variant_detail")
          )
        )
      ),

      # ── Tab 2: ACMG Evidence ──────────────────────────────────────────
      tabItem(tabName = "acmg",
        fluidRow(
          box(
            title = "ACMG Criteria Breakdown — Top Candidates",
            width = 12, solidHeader = TRUE, status = "warning",
            plotlyOutput("acmg_heatmap", height = "350px")
          )
        ),
        fluidRow(
          box(
            title = "CADD vs gnomAD AF (bubble = ACMG score)",
            width = 6, solidHeader = TRUE, status = "primary",
            plotlyOutput("cadd_af_plot", height = "300px")
          ),
          box(
            title = "Phenotype Match",
            width = 6, solidHeader = TRUE, status = "success",
            plotlyOutput("phenotype_plot", height = "300px")
          )
        )
      ),

      # ── Tab 3: Filtering Funnel ───────────────────────────────────────
      tabItem(tabName = "funnel",
        fluidRow(
          box(
            title = "Variant Filtering Funnel",
            width = 8, solidHeader = TRUE, status = "primary",
            plotlyOutput("funnel_plot", height = "420px")
          ),
          box(
            title = "Filter Parameters Used",
            width = 4, solidHeader = TRUE, status = "info",
            tableOutput("filter_params_table")
          )
        )
      ),

      # ── Tab 4: Clinical Report ────────────────────────────────────────
      tabItem(tabName = "report",
        fluidRow(
          box(
            title = "Clinical Interpretation Summary",
            width = 12, solidHeader = TRUE, status = "danger",
            uiOutput("clinical_report_ui")
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------
server <- function(input, output, session) {

  # Reactive filtered data
  filtered <- reactive({
    candidates |>
      filter(
        cadd >= input$cadd_min,
        max_af <= input$af_max * 1e-4,
        consequence %in% input$consequence_filter | consequence == "stop_gained",
        acmg_class %in% input$acmg_filter
      ) |>
      arrange(desc(acmg_score), desc(cadd))
  })

  # ── Value boxes ─────────────────────────────────────────────────────
  output$n_total <- renderValueBox({
    valueBox(nrow(filtered()), "Candidates shown", icon = icon("dna"), color = "blue")
  })

  output$n_pathogenic <- renderValueBox({
    n <- sum(filtered()$acmg_class %in% c("Pathogenic", "Likely Pathogenic"))
    valueBox(n, "Pathogenic/LP", icon = icon("exclamation-circle"), color = "red")
  })

  output$top_gene <- renderValueBox({
    top <- filtered() |> slice(1) |> pull(gene)
    valueBox(ifelse(length(top) > 0, top, "—"), "Top gene", icon = icon("star"), color = "yellow")
  })

  output$top_cadd <- renderValueBox({
    top <- filtered() |> slice(1) |> pull(cadd)
    valueBox(ifelse(length(top) > 0, round(top, 1), "—"), "Top CADD", icon = icon("chart-bar"), color = "green")
  })

  # ── Main candidates table ────────────────────────────────────────────
  output$candidates_table <- renderDT({
    df <- filtered() |>
      transmute(
        `#`         = row_number(),
        Variant     = paste0(chrom, ":", formatC(pos, format = "d", big.mark = ","),
                             " ", ref, ">", alt),
        Gene        = gene,
        `Protein`   = hgvsp,
        Consequence = consequence,
        `gnomAD AF` = formatC(max_af, format = "e", digits = 1),
        CADD        = round(cadd, 1),
        ClinVar     = clinvar,
        `ACMG Score`= acmg_score,
        `ACMG Class`= acmg_class,
        Criteria    = acmg_criteria
      )

    datatable(
      df,
      selection = "single",
      rownames  = FALSE,
      options   = list(
        pageLength = 10,
        scrollX    = TRUE,
        dom        = "Bfrtip",
        buttons    = c("copy")
      ),
      class = "stripe hover compact"
    ) |>
      formatStyle(
        "ACMG Class",
        backgroundColor = styleEqual(
          names(acmg_colors),
          unname(acmg_colors)
        ),
        color = "white",
        fontWeight = "bold",
        borderRadius = "4px"
      ) |>
      formatStyle("CADD",
        background = styleColorBar(c(0, 50), "#74add1"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
  })

  # ── Selected variant detail panel ────────────────────────────────────
  output$selected_variant_detail <- renderUI({
    sel <- input$candidates_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      return(p("Select a row above to view details.", style = "color:#aaa"))
    }

    v <- filtered() |> slice(sel)

    tagList(
      fluidRow(
        column(4,
          tags$b("Variant"), br(), tags$code(paste0(v$chrom, ":", v$pos, " ", v$ref, ">", v$alt)), br(), br(),
          tags$b("Gene"), br(), v$gene, br(), br(),
          tags$b("Disease (OMIM)"), br(), paste0(v$disease, " (", v$omim_id, ")")
        ),
        column(4,
          tags$b("HGVS protein"), br(), tags$code(v$hgvsp), br(), br(),
          tags$b("HGVS coding"), br(), tags$code(v$hgvsc), br(), br(),
          tags$b("Consequence"), br(), v$consequence
        ),
        column(4,
          tags$b("ACMG Classification"), br(),
          span(v$acmg_class,
               style = paste0("background:", acmg_colors[v$acmg_class],
                              "; color:white; padding:3px 8px; border-radius:4px; font-weight:bold")),
          br(), br(),
          tags$b("Active Criteria"), br(), v$acmg_criteria, br(), br(),
          tags$b("ACMG Score"), br(), paste0(v$acmg_score, " / 37")
        )
      ),
      hr(),
      fluidRow(
        column(3, tags$b("gnomAD AF"), br(), formatC(v$max_af, format = "e", digits = 2)),
        column(3, tags$b("CADD PHRED"), br(), round(v$cadd, 1)),
        column(3, tags$b("SIFT"), br(), ifelse(is.na(v$sift_score), "—", round(v$sift_score, 4))),
        column(3, tags$b("PolyPhen"), br(), ifelse(is.na(v$polyphen_score), "—", round(v$polyphen_score, 3)))
      ),
      hr(),
      fluidRow(
        column(12,
          tags$b("Phenotype match:"),
          tags$span(ifelse(v$phenotype_seizures, "✓ Seizures", "✗ Seizures"),
                    style = paste0("margin-right:12px; color:", ifelse(v$phenotype_seizures,"#2ecc71","#ccc"))),
          tags$span(ifelse(v$phenotype_id,       "✓ Intellectual disability", "✗ ID"),
                    style = paste0("margin-right:12px; color:", ifelse(v$phenotype_id,"#2ecc71","#ccc"))),
          tags$span(ifelse(v$phenotype_dd,       "✓ Developmental delay", "✗ DD"),
                    style = paste0("color:", ifelse(v$phenotype_dd,"#2ecc71","#ccc")))
        )
      )
    )
  })

  # ── ACMG heatmap ────────────────────────────────────────────────────
  output$acmg_heatmap <- renderPlotly({
    top <- filtered() |> slice_head(n = 5)
    criteria_cols <- c("PVS1","PS1","PS2","PM2","PP2","PP3","PP4")

    # Only use columns that exist
    available <- intersect(criteria_cols, names(top))
    if (length(available) == 0) return(plotly_empty())

    mat <- top |>
      mutate(label = paste0(gene, "\n", hgvsp)) |>
      select(label, all_of(available)) |>
      pivot_longer(-label, names_to = "criterion", values_to = "score")

    plot_ly(mat,
            x = ~criterion, y = ~label, z = ~score,
            type = "heatmap",
            colorscale = list(c(0,"#f7f7f7"), c(1,"#d73027")),
            showscale = FALSE,
            hovertemplate = "%{y}<br>%{x}: %{z} pts<extra></extra>") |>
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        margin = list(l = 160, b = 40)
      )
  })

  # ── CADD vs AF scatter ────────────────────────────────────────────
  output$cadd_af_plot <- renderPlotly({
    df <- filtered() |>
      mutate(af_log = -log10(pmax(max_af, 1e-9)))

    plot_ly(df,
            x = ~af_log, y = ~cadd,
            size = ~acmg_score, sizes = c(10, 60),
            color = ~acmg_class,
            colors = acmg_colors,
            text = ~paste0(gene, "<br>", hgvsp),
            hoverinfo = "text+x+y",
            type = "scatter", mode = "markers",
            marker = list(opacity = 0.85)) |>
      layout(
        xaxis = list(title = "-log10(gnomAD AF)"),
        yaxis = list(title = "CADD PHRED score"),
        showlegend = TRUE,
        legend = list(orientation = "h", y = -0.25)
      )
  })

  # ── Phenotype match bar ──────────────────────────────────────────
  output$phenotype_plot <- renderPlotly({
    df <- filtered() |>
      mutate(
        label = paste0(gene, " ", hgvsp),
        n_phenotypes = as.integer(phenotype_seizures) +
                       as.integer(phenotype_id) +
                       as.integer(phenotype_dd)
      ) |>
      arrange(desc(n_phenotypes), desc(acmg_score))

    plot_ly(df,
            x = ~n_phenotypes, y = ~reorder(label, n_phenotypes),
            type = "bar", orientation = "h",
            marker = list(color = ~acmg_score,
                          colorscale = "Reds",
                          showscale = TRUE,
                          colorbar = list(title = "ACMG")),
            text = ~paste0(gene, " — ", n_phenotypes, "/3 phenotypes"),
            hoverinfo = "text") |>
      layout(
        xaxis = list(title = "Phenotypes matched (max 3)", dtick = 1),
        yaxis = list(title = ""),
        margin = list(l = 180)
      )
  })

  # ── Filtering funnel ──────────────────────────────────────────────
  output$funnel_plot <- renderPlotly({
    funnel_data <- tibble(
      step = c(
        "1. Raw VCF", "2. Quality filter\n(QUAL≥30, DP≥10)",
        "3. Coding regions\n(exonic + splice)",
        "4. Rare variants\n(gnomAD AF < 1%)",
        "5. Functional impact\n(LoF + CADD≥15)",
        "6. HPO phenotype\ngene list",
        "7. Final\ncandidates"
      ),
      n = c(3200000, 1480000, 44600, 2100, 148, 22, nrow(filtered()))
    ) |>
      mutate(
        step = factor(step, levels = rev(step)),
        pct  = round(n / 3200000 * 100, 3)
      )

    plot_ly(funnel_data,
            y = ~step, x = ~n,
            type = "bar", orientation = "h",
            marker = list(
              color = colorRampPalette(c("#d73027","#fdae61","#74add1","#4dac26"))(7),
              line  = list(color = "white", width = 1)
            ),
            text  = ~paste0(formatC(n, format="d", big.mark=","), " (", pct, "%)"),
            textposition = "outside",
            hoverinfo = "y+text") |>
      layout(
        xaxis = list(title = "Number of variants", type = "log"),
        yaxis = list(title = ""),
        margin = list(l = 220, r = 120)
      )
  })

  # ── Filter params table ───────────────────────────────────────────
  output$filter_params_table <- renderTable({
    tibble(
      Parameter      = c("Min QUAL", "Min DP", "Min GQ", "Max gnomAD AF",
                         "Min CADD (missense)", "HPO terms"),
      Value          = c("30", "10", "20",
                         formatC(input$af_max * 1e-4, format = "e", digits = 1),
                         as.character(input$cadd_min),
                         "HP:0001250, HP:0001249, HP:0001263")
    )
  }, striped = TRUE, hover = TRUE, bordered = TRUE)

  # ── Clinical report UI ────────────────────────────────────────────
  output$clinical_report_ui <- renderUI({
    top1 <- candidates |> arrange(desc(acmg_score)) |> slice(1)

    tagList(
      h4("Clinical Genetics Report"),
      tags$table(class = "table table-bordered",
        tags$tr(tags$th("Field"), tags$th("Finding")),
        tags$tr(tags$td("Patient"), tags$td("8-year-old (de-identified)")),
        tags$tr(tags$td("Phenotype"), tags$td("Intellectual disability, Seizures (onset age 2), Developmental delay, Facial dysmorphisms")),
        tags$tr(tags$td("Test"), tags$td("Whole genome sequencing (GRCh38)")),
        tags$tr(tags$td("Variant identified"),
                tags$td(tags$code(paste0(top1$chrom, ":", top1$pos, " ", top1$ref, ">", top1$alt)))),
        tags$tr(tags$td("Gene"), tags$td(tags$b(top1$gene))),
        tags$tr(tags$td("HGVS nomenclature"),
                tags$td(paste0(top1$hgvsp, " (", top1$hgvsc, ")"))),
        tags$tr(tags$td("Classification"),
                tags$td(span(top1$acmg_class,
                             style = paste0("background:", acmg_colors[top1$acmg_class],
                                            "; color:white; padding:2px 8px; border-radius:4px")))),
        tags$tr(tags$td("ACMG criteria"), tags$td(top1$acmg_criteria)),
        tags$tr(tags$td("Disease"), tags$td(paste0(top1$disease, " (OMIM: ", top1$omim_id, ")"))),
        tags$tr(tags$td("Inheritance"), tags$td("De novo — autosomal dominant")),
        tags$tr(tags$td("Recommendations"),
                tags$td("1. Parental testing to confirm de novo status. 2. Refer to pediatric neurology. 3. Genetic counseling for family. 4. Consider sodium channel-targeted therapy (e.g., avoid sodium channel blockers in Dravet)."))
      )
    )
  })

  # ── CSV download ──────────────────────────────────────────────────
  output$download_csv <- downloadHandler(
    filename = function() paste0("variant_candidates_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(filtered(), file, row.names = FALSE)
  )
}

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
