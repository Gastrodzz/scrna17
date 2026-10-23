# modules/pathway_network/pn_ui.R
# UI for the unified Pathway & Network Analysis tab — 8 sub-tabs:
#   Pathway Enrichment | Pathway Network | Gene Regulatory Networks
#   | Co-expression | Cell-Cell Communication | TF Activity
#   | Gene Set Scoring | Composition
#
# Design: minimalism + scientific dashboard, app.css design tokens.
# Sourced by pn_server.R (which is sourced by app.R).

pnUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    .pn_styles(id),

    shiny::tabsetPanel(
      id = ns("pn_tabs"), type = "tabs",

      # ============ 1. PATHWAY ENRICHMENT ============
      shiny::tabPanel(
        title = "Pathway Enrichment",
        icon = shiny::icon("dna"),
        shiny::br(),
        shiny::div(class = paste0("pn-shell-", id),

          # left: config
          shiny::div(class = paste0("pn-config-", id),
            shiny::div(class = "pn-header",
              shiny::div(class = "pn-header-icon", shiny::icon("dna")),
              shiny::div(
                shiny::div(class = "pn-title", "Pathway Enrichment"),
                shiny::div(class = "pn-subtitle",
                  "ORA + GSEA across 14 databases")
              )
            ),
            shiny::uiOutput(ns("enrich_status_strip")),

            .pn_section("Database", "database",
              shiny::selectInput(ns("enrich_db"), NULL,
                choices = c(
                  "MSigDB Hallmark"        = "H",
                  "MSigDB C2 (Curated)"    = "C2",
                  "MSigDB C5 (GO)"         = "C5",
                  "MSigDB C6 (Oncogenic)"  = "C6",
                  "MSigDB KEGG (CP)"       = "MSIGDB_KEGG",
                  "MSigDB Reactome (CP)"   = "MSIGDB_REACTOME",
                  "MSigDB WikiPathways"    = "MSIGDB_WP",
                  "KEGG (live API)"        = "KEGG",
                  "Reactome (live API)"    = "REACTOME",
                  "GO Biological Process"  = "GO_BP",
                  "GO Molecular Function"  = "GO_MF",
                  "GO Cellular Component"  = "GO_CC",
                  "WikiPathways (live)"    = "WP",
                  "STRING PPI"             = "STRING"
                ),
                selected = "H", width = "100%")
            ),

            .pn_section("Method", "calculator",
              shinyWidgets::radioGroupButtons(
                ns("enrich_method"), label = NULL,
                choices = c("ORA" = "ORA", "GSEA" = "GSEA", "Both" = "BOTH"),
                selected = "BOTH", justified = TRUE, size = "sm",
                status = "primary")
            ),

            .pn_section("ORA gene selection", "filter",
              shiny::div(class = "pn-grid-2",
                shiny::numericInput(ns("ora_pval_threshold"),
                  "Adj p threshold", value = 0.05,
                  min = 0.001, max = 0.1, step = 0.005),
                shiny::numericInput(ns("ora_logfc_threshold"),
                  "Min |log2FC|", value = 0.25,
                  min = 0, max = 5, step = 0.05)
              ),
              shiny::numericInput(ns("ora_pval_cutoff"),
                "ORA enrichment p cutoff", value = 0.05,
                min = 0.001, max = 0.25, step = 0.01, width = "100%")
            ),

            .pn_section("GSEA", "wave-square",
              shiny::selectInput(ns("gsea_metric"), "Ranking metric",
                choices = c(
                  "signed -log10(p) × sign(FC)" = "signed_neglog10p",
                  "log2 fold change"            = "log2fc",
                  "t-statistic (if available)"  = "t_statistic"),
                selected = "signed_neglog10p", width = "100%"),
              shiny::div(class = "pn-grid-3",
                shiny::numericInput(ns("gsea_min_size"), "Min size",
                  value = 15, min = 5, max = 200, step = 5),
                shiny::numericInput(ns("gsea_max_size"), "Max size",
                  value = 500, min = 50, max = 2000, step = 50),
                shiny::numericInput(ns("gsea_nperm"), "Permutations",
                  value = 1000, min = 1000, max = 50000, step = 1000)
              )
            ),

            .pn_section("Display", "eye",
              shiny::numericInput(ns("max_pathways"), "Max pathways shown",
                value = 30, min = 5, max = 200, step = 5, width = "100%"),
              shiny::checkboxInput(ns("use_custom_bg"),
                "Use custom background", value = FALSE),
              shiny::conditionalPanel(
                condition = sprintf("input['%s']", ns("use_custom_bg")),
                shiny::textAreaInput(ns("custom_bg"),
                  "Background genes (one per line)", rows = 3))
            ),

            .pn_section("Reproducibility", "shield-alt",
              shiny::numericInput(ns("enrich_seed"), "Random seed",
                value = 42, min = 1, max = 99999, step = 1, width = "100%")
            ),

            shiny::div(class = "pn-cta-wrap",
              shiny::actionButton(ns("run_enrich"), "Run Enrichment",
                class = "pn-cta", icon = shiny::icon("play"))
            )
          ),

          shiny::div(class = paste0("pn-results-", id),
            shiny::uiOutput(ns("enrich_results_ui"))
          )
        )
      ),

      # ============ 1b. COMPARATIVE ENRICHMENT ============
      shiny::tabPanel(
        title = "Comparative Enrichment",
        icon = shiny::icon("balance-scale"),
        shiny::br(),
        shiny::div(class = paste0("pn-shell-", id),
          shiny::div(class = paste0("pn-config-", id),
            shiny::div(class = "pn-header",
              shiny::div(class = "pn-header-icon",
                         shiny::icon("balance-scale")),
              shiny::div(
                shiny::div(class = "pn-title", "Comparative Enrichment"),
                shiny::div(class = "pn-subtitle",
                  "Side-by-side ORA + GSEA for any 2+ saved DE runs")
              )
            ),

            .pn_section("DE runs to compare", "list-check",
              shiny::uiOutput(ns("cmp_history_strip")),
              shinyWidgets::pickerInput(
                ns("cmp_de_runs"), NULL,
                choices = character(0), multiple = TRUE,
                options = list(`live-search` = TRUE,
                               `actions-box` = TRUE,
                               `selected-text-format` = "count > 2",
                               size = 10,
                               `none-selected-text` = "Pick 2 or more saved DE runs..."),
                width = "100%"),
              shiny::div(class = "pn-help",
                "Saved DE runs from this session appear here. Run DE on multiple ",
                "comparisons (e.g. mutant vs control within each cluster) — every ",
                "successful run is auto-saved.")
            ),

            .pn_section("Database", "database",
              shiny::selectInput(ns("cmp_db"), NULL,
                choices = c(
                  "MSigDB Hallmark"        = "H",
                  "MSigDB C2 (Curated)"    = "C2",
                  "MSigDB C5 (GO)"         = "C5",
                  "MSigDB C6 (Oncogenic)"  = "C6",
                  "MSigDB KEGG (CP)"       = "MSIGDB_KEGG",
                  "MSigDB Reactome (CP)"   = "MSIGDB_REACTOME",
                  "MSigDB WikiPathways"    = "MSIGDB_WP",
                  "KEGG (live API)"        = "KEGG",
                  "Reactome (live API)"    = "REACTOME",
                  "GO Biological Process"  = "GO_BP",
                  "GO Molecular Function"  = "GO_MF",
                  "GO Cellular Component"  = "GO_CC",
                  "WikiPathways (live)"    = "WP",
                  "STRING PPI"             = "STRING"
                ),
                selected = "H", width = "100%")
            ),

            .pn_section("Method", "calculator",
              shinyWidgets::radioGroupButtons(
                ns("cmp_method"), label = NULL,
                choices = c("ORA" = "ORA", "GSEA" = "GSEA", "Both" = "BOTH"),
                selected = "BOTH", justified = TRUE, size = "sm",
                status = "primary")
            ),

            .pn_section("ORA gene selection", "filter",
              shiny::div(class = "pn-grid-2",
                shiny::numericInput(ns("cmp_ora_pval"),
                  "Adj p threshold", value = 0.05,
                  min = 0.001, max = 0.1, step = 0.005),
                shiny::numericInput(ns("cmp_ora_lfc"),
                  "Min |log2FC|", value = 0.25,
                  min = 0, max = 5, step = 0.05)
              ),
              shiny::numericInput(ns("cmp_ora_pval_cutoff"),
                "ORA enrichment p cutoff", value = 0.05,
                min = 0.001, max = 0.25, step = 0.01, width = "100%")
            ),

            .pn_section("GSEA", "wave-square",
              shiny::div(class = "pn-grid-3",
                shiny::numericInput(ns("cmp_gsea_min"), "Min size",
                  value = 15, min = 5, max = 200, step = 5),
                shiny::numericInput(ns("cmp_gsea_max"), "Max size",
                  value = 500, min = 50, max = 2000, step = 50),
                shiny::numericInput(ns("cmp_gsea_nperm"), "Permutations",
                  value = 1000, min = 1000, max = 50000, step = 1000)
              )
            ),

            .pn_section("Display", "eye",
              shiny::numericInput(ns("cmp_top_n"),
                "Top pathways per condition",
                value = 20, min = 5, max = 100, step = 5, width = "100%")
            ),

            .pn_section("Reproducibility", "shield-alt",
              shiny::numericInput(ns("cmp_seed"), "Random seed",
                value = 42, min = 1, max = 99999, step = 1, width = "100%")
            ),

            shiny::div(class = "pn-cta-wrap",
              shiny::actionButton(ns("run_compare_enrich"),
                "Run Comparative Enrichment",
                class = "pn-cta", icon = shiny::icon("play"))
            )
          ),
          shiny::div(class = paste0("pn-results-", id),
            shiny::uiOutput(ns("cmp_results_ui"))
          )
        )
      ),

      # ============ 2. PATHWAY NETWORK ============
      shiny::tabPanel(
        title = "Pathway Network",
        icon = shiny::icon("project-diagram"),
        shiny::br(),
        shiny::div(class = paste0("pn-shell-", id),
          shiny::div(class = paste0("pn-config-", id),
            shiny::div(class = "pn-header",
              shiny::div(class = "pn-header-icon",
                         shiny::icon("project-diagram")),
              shiny::div(
                shiny::div(class = "pn-title", "Pathway Network"),
                shiny::div(class = "pn-subtitle",
                  "Connect enriched pathways by gene-set overlap")
              )
            ),
            .pn_section("Status", "info-circle",
              shiny::uiOutput(ns("pwn_status"))
            ),
            .pn_section("Thresholds", "sliders-h",
              shiny::numericInput(ns("pwn_top_n"), "Top pathways to use",
                value = 30, min = 5, max = 100, step = 5, width = "100%"),
              shiny::div(class = "pn-grid-2",
                shiny::numericInput(ns("pwn_overlap"), "Min overlap",
                  value = 3, min = 2, max = 30, step = 1),
                shiny::numericInput(ns("pwn_jaccard"), "Min Jaccard",
                  value = 0.1, min = 0.05, max = 0.9, step = 0.05))
            ),
            shiny::div(class = "pn-cta-wrap",
              shiny::actionButton(ns("run_pwn"), "Build Pathway Network",
                class = "pn-cta", icon = shiny::icon("play"))
            )
          ),
          shiny::div(class = paste0("pn-results-", id),
            shiny::uiOutput(ns("pwn_results_ui"))
          )
        )
      ),

      # ============ 3. GENE REGULATORY NETWORKS ============
      shiny::tabPanel(
        title = "Gene Regulatory Networks",
        icon = shiny::icon("share-alt"),
        shiny::br(),
        shiny::div(class = paste0("pn-shell-", id),
          shiny::div(class = paste0("pn-config-", id),
            shiny::div(class = "pn-header",
              shiny::div(class = "pn-header-icon", shiny::icon("share-alt")),
              shiny::div(
                shiny::div(class = "pn-title", "Gene Regulatory Networks"),
                shiny::div(class = "pn-subtitle",
                  "Co-expression modules via WGCNA / correlation")
              )
            ),
            .pn_section("Method", "calculator",
              shinyWidgets::radioGroupButtons(
                ns("grn_method"), label = NULL,
                choices = c("WGCNA" = "wgcna", "Correlation" = "correlation"),
                selected = "wgcna", justified = TRUE, size = "sm",
                status = "primary")
            ),
            .pn_section("Input genes", "filter",
              shiny::div(class = "pn-grid-2",
                shiny::numericInput(ns("grn_top_genes"),
                  "Variable genes", value = 2000,
                  min = 200, max = 5000, step = 100),
                shiny::numericInput(ns("grn_min_cells"),
                  "Min cells expressing", value = 10, min = 1, max = 1000)),
              shiny::selectInput(ns("grn_cell_type"),
                "Focus on cell type",
                choices = c("All cells" = "all"), width = "100%")
            ),
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'wgcna'", ns("grn_method")),
              .pn_section("WGCNA parameters", "sliders-h",
                shiny::div(class = "pn-grid-2",
                  shiny::numericInput(ns("wgcna_power"),
                    "Power (0=auto)", value = 0,
                    min = 0, max = 30, step = 1),
                  shiny::numericInput(ns("wgcna_min_module_size"),
                    "Min module size", value = 30,
                    min = 10, max = 200, step = 1)),
                shiny::numericInput(ns("wgcna_merge_threshold"),
                  "Merge threshold", value = 0.25,
                  min = 0.05, max = 0.5, step = 0.05, width = "100%"),
                shiny::selectInput(ns("wgcna_tom_type"), "TOM type",
                  choices = c("Signed" = "signed", "Unsigned" = "unsigned",
                              "Signed hybrid" = "signed hybrid"),
                  selected = "signed", width = "100%"),
                shiny::checkboxInput(ns("wgcna_detect_outliers"),
                  "Detect outlier samples", value = TRUE),
                shiny::div(class = "pn-callout pn-callout-warning",
                  shiny::icon("exclamation-triangle"),
                  shiny::span("Pseudobulk aggregation requires ≥3 biological replicates."))
              )
            ),
            shiny::conditionalPanel(
              condition = sprintf("input['%s'] == 'correlation'", ns("grn_method")),
              .pn_section("Correlation parameters", "sliders-h",
                shiny::selectInput(ns("cor_method"), "Method",
                  choices = c("Spearman" = "spearman",
                              "Pearson" = "pearson",
                              "Kendall" = "kendall"),
                  selected = "spearman", width = "100%"),
                shiny::div(class = "pn-grid-2",
                  shiny::numericInput(ns("cor_threshold"), "|r|",
                    value = 0.7, min = 0.1, max = 0.95, step = 0.05),
                  shiny::numericInput(ns("cor_p_threshold"), "p",
                    value = 0.05, min = 0.001, max = 0.1, step = 0.005))
              )
            ),
            shiny::div(class = "pn-cta-wrap",
              shiny::actionButton(ns("run_grn"), "Build Network",
                class = "pn-cta", icon = shiny::icon("play"))
            )
          ),
          shiny::div(class = paste0("pn-results-", id),
            shiny::uiOutput(ns("grn_results_ui"))
          )
        )
      ),

      # ============ 4. CO-EXPRESSION ============
      shiny::tabPanel(
        title = "Co-expression",
        icon = shiny::icon("network-wired"),
        shiny::br(),
        shiny::div(class = paste0("pn-shell-", id),
          shiny::div(class = paste0("pn-config-", id),
            shiny::div(class = "pn-header",
              shiny::div(class = "pn-header-icon", shiny::icon("network-wired")),
              shiny::div(
                shiny::div(class = "pn-title", "Co-expression"),
                shiny::div(class = "pn-subtitle",
                  "Pairwise correlation across user-supplied genes")
              )
            ),
            .pn_section("Genes", "list",
              shiny::textAreaInput(ns("coexp_genes"),
                "Genes (comma-separated)", rows = 4,
                placeholder = "e.g. CD3D, CD4, CD8A, IL7R, MS4A1, GNLY")
            ),
            .pn_section("Parameters", "sliders-h",
              shiny::selectInput(ns("coexp_method"), "Correlation method",
                choices = c("Spearman" = "spearman", "Pearson" = "pearson"),
                selected = "spearman", width = "100%"),
              shiny::div(class = "pn-grid-2",
                shiny::numericInput(ns("coexp_threshold"), "|r| threshold",
                  value = 0.3, min = 0.05, max = 0.95, step = 0.05),
                shiny::numericInput(ns("coexp_top_n"), "Max edges",
                  value = 200, min = 10, max = 5000, step = 10))
            ),
            shiny::div(class = "pn-cta-wrap",
              shiny::actionButton(ns("run_coexp"),
                "Build Co-expression Network",
                class = "pn-cta", icon = shiny::icon("play"))
            )
          ),
          shiny::div(class = paste0("pn-results-", id),
            shiny::uiOutput(ns("coexp_results_ui"))
          )
        )
      ),

      # ============ 5. CCC (delegated) ============
      shiny::tabPanel(
        title = "Cell-Cell Communication",
        icon = shiny::icon("people-arrows"),
        shiny::br(),
        shiny::div(class = "pn-info-banner",
          shiny::icon("info-circle"),
          shiny::div(
            shiny::strong("Workflow: "),
            "Pick the cell-type column + species, click ",
            shiny::tags$em("Run Communication Analysis"),
            ", wait for the spinner. ",
            shiny::strong("Note: "),
            "datasets >25k cells with default permutations can take 5-30 min ",
            "(LIANA / permutation test scales with n_cells × n_LR_pairs × nperm)."
          )
        ),
        if (exists("cellCommUI", mode = "function"))
          cellCommUI(ns("ccc"))
        else
          shiny::p("Cell-cell communication module not available.")
      ),

      # ============ 6. TF Activity (delegated) ============
      shiny::tabPanel(
        title = "TF Activity",
        icon = shiny::icon("atom"),
        shiny::br(),
        if (exists("tfActivityUI", mode = "function"))
          tfActivityUI(ns("tf"))
        else shiny::p("TF activity module not available.")
      ),

      # ============ 7. Gene Set Scoring (delegated) ============
      shiny::tabPanel(
        title = "Gene Set Scoring",
        icon = shiny::icon("list-check"),
        shiny::br(),
        if (exists("geneSetScoringUI", mode = "function"))
          geneSetScoringUI(ns("geneset"))
        else shiny::p("Gene set scoring module not available.")
      ),

      # ============ 8. Composition (delegated) ============
      shiny::tabPanel(
        title = "Composition",
        icon = shiny::icon("chart-bar"),
        shiny::br(),
        if (exists("compositionUI", mode = "function"))
          compositionUI(ns("composition"))
        else shiny::p("Composition module not available.")
      )
    )
  )
}

# ---------------------------------------------------------------------------
# Style block — scoped via the module id. Mirrors de_ui.R styling.
# ---------------------------------------------------------------------------
.pn_styles <- function(id) {
  # Use gsub('__ID__', id, ...) instead of sprintf to keep CSS percent signs
  # literal (sprintf would treat 'width: 100%;' as a format specifier).
  css <- "
.pn-shell-__ID__ {
  display: flex; gap: 18px; align-items: stretch;
  font-family: var(--font-sans); color: var(--text-primary);
}
.pn-config-__ID__ {
  flex: 0 0 360px; max-height: 88vh; overflow-y: auto;
  padding: 18px;
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius-lg);
  box-shadow: var(--shadow-sm);
}
.pn-config-__ID__::-webkit-scrollbar { width: 8px; }
.pn-config-__ID__::-webkit-scrollbar-thumb {
  background: var(--border); border-radius: 4px;
}
.pn-results-__ID__ { flex: 1 1 auto; min-width: 0; }

.pn-config-__ID__ .pn-header {
  display: flex; align-items: center; gap: 12px;
  padding-bottom: 14px; margin-bottom: 14px;
  border-bottom: 1px solid var(--border);
}
.pn-config-__ID__ .pn-header-icon {
  width: 38px; height: 38px; border-radius: var(--radius);
  background: var(--primary-light); color: var(--primary);
  display: flex; align-items: center; justify-content: center;
  font-size: 18px;
}
.pn-config-__ID__ .pn-title {
  font-weight: 700; font-size: 16px; color: var(--text-primary);
  letter-spacing: -0.01em;
}
.pn-config-__ID__ .pn-subtitle {
  font-size: 12px; color: var(--text-secondary); margin-top: 2px;
}

.pn-config-__ID__ .pn-section {
  margin: 18px 0 4px; padding-bottom: 14px;
  border-bottom: 1px dashed var(--border-light);
}
.pn-config-__ID__ .pn-section:last-of-type { border-bottom: none; }
.pn-config-__ID__ .pn-section-title {
  display: flex; align-items: center; gap: 8px;
  font-weight: 600; font-size: 11px; letter-spacing: 0.06em;
  text-transform: uppercase; color: var(--text-secondary);
  margin-bottom: 10px;
}
.pn-config-__ID__ .pn-section-title .icon {
  color: var(--primary); width: 14px; text-align: center;
}

.pn-config-__ID__ .form-control,
.pn-config-__ID__ .selectize-input {
  border-radius: var(--radius-sm) !important;
  border: 1px solid var(--border) !important;
  font-size: 13px;
  transition: border-color var(--t), box-shadow var(--t);
}
.pn-config-__ID__ .form-control:focus,
.pn-config-__ID__ .selectize-input.focus {
  border-color: var(--primary) !important;
  box-shadow: 0 0 0 3px var(--border-focus) !important;
  outline: none;
}
.pn-config-__ID__ label {
  font-size: 12px; font-weight: 500; color: var(--text-secondary);
  margin-bottom: 4px;
}

.pn-config-__ID__ .pn-grid-2 {
  display: grid; grid-template-columns: 1fr 1fr; gap: 10px;
}
.pn-config-__ID__ .pn-grid-3 {
  display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 10px;
}

.pn-config-__ID__ .pn-callout {
  display: flex; gap: 8px; align-items: flex-start;
  padding: 8px 10px; border-radius: var(--radius-sm);
  font-size: 11px; line-height: 1.5; margin-top: 8px;
}
.pn-config-__ID__ .pn-callout-warning {
  background: var(--warning-light); color: var(--warning);
}
.pn-config-__ID__ .pn-callout-warning i { margin-top: 2px; }
.pn-config-__ID__ .pn-callout-warning span { color: var(--text-primary); }

.pn-config-__ID__ .pn-status-strip {
  display: flex; gap: 10px; align-items: center;
  padding: 10px 12px; border-radius: var(--radius-sm);
  background: var(--bg); border: 1px solid var(--border-light);
  font-size: 11px; color: var(--text-secondary); margin-bottom: 10px;
}
.pn-config-__ID__ .pn-status-strip .dot {
  width: 8px; height: 8px; border-radius: 50%;
  background: var(--success); flex: 0 0 auto;
}
.pn-config-__ID__ .pn-status-strip .dot.muted { background: var(--text-muted); }
.pn-config-__ID__ .pn-status-strip .meta { color: var(--text-primary); font-weight: 500; }

.pn-config-__ID__ .pn-cta-wrap { margin-top: 18px; }
.pn-config-__ID__ .pn-cta {
  width: 100%%; padding: 12px 16px;
  background: var(--success); color: #fff;
  border: none; border-radius: var(--radius);
  font-weight: 600; font-size: 14px; letter-spacing: 0.01em;
  cursor: pointer;
  box-shadow: 0 1px 2px rgba(22,163,74,0.18);
  transition: background var(--t), box-shadow var(--t), transform var(--t);
}
.pn-config-__ID__ .pn-cta:hover {
  background: #15803D;
  box-shadow: 0 3px 8px rgba(22,163,74,0.25);
}
.pn-config-__ID__ .pn-cta:active { transform: translateY(1px); }
.pn-config-__ID__ .pn-cta:focus-visible {
  outline: none;
  box-shadow: 0 0 0 3px rgba(22,163,74,0.30);
}
.pn-config-__ID__ .pn-cta i { margin-right: 8px; }

/* Results panel */
.pn-results-__ID__ .pn-empty {
  display: flex; flex-direction: column; align-items: center;
  justify-content: center; min-height: 60vh;
  padding: 40px; text-align: center;
  background: var(--surface); border: 1px dashed var(--border);
  border-radius: var(--radius-lg);
}
.pn-results-__ID__ .pn-empty .icon {
  width: 64px; height: 64px; border-radius: 50%%;
  background: var(--primary-light); color: var(--primary);
  display: flex; align-items: center; justify-content: center;
  font-size: 26px; margin-bottom: 16px;
}
.pn-results-__ID__ .pn-empty h4 {
  margin: 0 0 8px; color: var(--text-primary); font-weight: 600;
}
.pn-results-__ID__ .pn-empty p {
  color: var(--text-secondary); font-size: 13px; max-width: 420px;
  line-height: 1.55; margin: 0;
}
.pn-results-__ID__ .pn-results-card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: var(--radius-lg); padding: 16px;
  box-shadow: var(--shadow-sm);
}
.pn-results-__ID__ .nav-tabs > li > a {
  font-size: 13px; font-weight: 500; color: var(--text-secondary);
  border-radius: 0; transition: color var(--t), border-color var(--t);
}
.pn-results-__ID__ .nav-tabs > li.active > a {
  color: var(--primary); border-bottom: 2px solid var(--primary);
}

/* Info banner used at top of CCC sub-tab */
.pn-info-banner {
  display: flex; gap: 12px; align-items: flex-start;
  padding: 12px 14px; margin-bottom: 12px;
  background: var(--info-light); border: 1px solid #BAE6FD;
  border-left: 3px solid var(--info); border-radius: var(--radius);
  font-size: 12px; line-height: 1.5; color: var(--text-primary);
}
.pn-info-banner > i {
  color: var(--info); font-size: 16px; margin-top: 1px; flex: 0 0 auto;
}

/* Comparative Enrichment — history chips strip + helpers */
.pn-config-__ID__ .pn-history-empty {
  font-size: 11px; color: var(--text-muted); font-style: italic;
  padding: 8px 10px; background: var(--bg);
  border: 1px dashed var(--border); border-radius: var(--radius-sm);
  margin-bottom: 8px;
}
.pn-config-__ID__ .pn-history-chips {
  display: flex; flex-wrap: wrap; gap: 4px; margin-bottom: 8px;
}
.pn-config-__ID__ .pn-history-chip {
  display: inline-flex; align-items: center; gap: 6px;
  padding: 3px 8px; background: var(--primary-light); color: var(--primary);
  border: 1px solid #BFDBFE; border-radius: 999px;
  font-size: 10px; font-weight: 500;
}
.pn-config-__ID__ .pn-history-chip i { font-size: 8px; }
.pn-config-__ID__ .pn-history-count {
  font-size: 10px; color: var(--text-secondary); margin-bottom: 4px;
}
.pn-config-__ID__ .pn-help {
  font-size: 11px; color: var(--text-muted); margin-top: 6px;
  line-height: 1.5;
}

/* Comparative results layout */
.pn-results-__ID__ .pn-cmp-metrics {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
  gap: 12px; margin-bottom: 16px;
}
.pn-results-__ID__ .pn-cmp-metric {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: var(--radius); padding: 12px 14px;
  position: relative;
}
.pn-results-__ID__ .pn-cmp-metric .label {
  font-size: 10px; color: var(--text-secondary);
  text-transform: uppercase; letter-spacing: 0.04em; font-weight: 600;
}
.pn-results-__ID__ .pn-cmp-metric .value {
  font-size: 22px; font-weight: 700; color: var(--text-primary);
  font-variant-numeric: tabular-nums; margin-top: 2px;
}
.pn-results-__ID__ .pn-cmp-metric .sub {
  font-size: 11px; color: var(--text-muted); margin-top: 2px;
}
"
  shiny::tags$style(shiny::HTML(gsub("__ID__", id, css, fixed = TRUE)))
}

.pn_section <- function(title, icon, ...) {
  shiny::div(class = "pn-section",
    shiny::div(class = "pn-section-title",
      shiny::span(class = "icon", shiny::icon(icon)),
      shiny::span(title)
    ),
    ...
  )
}
