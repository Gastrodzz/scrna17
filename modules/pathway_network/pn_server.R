# modules/pathway_network/pn_server.R
# Server for the unified Pathway & Network Analysis tab.
# Sources core_robustness.R, pn_core.R, pn_enrichment.R, pn_network.R, pn_ui.R.
# Also sources the existing CCC / TF / geneset / composition sub-modules from
# modules/downstream/ so they can be re-mounted under this tab.

local({
  this_dir <- tryCatch({
    f <- sys.frame(1)$ofile
    if (!is.null(f) && nzchar(f)) normalizePath(dirname(f), mustWork = TRUE)
    else file.path(getwd(), "modules", "pathway_network")
  }, error = function(e) file.path(getwd(), "modules", "pathway_network"))

  src <- function(rel, base = this_dir) {
    p <- normalizePath(file.path(base, rel), mustWork = FALSE)
    if (file.exists(p)) source(p, local = FALSE)
  }
  modules_root <- dirname(this_dir)
  src("../core_robustness.R",                       base = modules_root)
  src("../downstream/downstream_core.R",            base = modules_root)
  src("../downstream/downstream_ccc.R",             base = modules_root)
  src("../downstream/downstream_tf.R",              base = modules_root)
  src("../downstream/downstream_geneset.R",         base = modules_root)
  src("../downstream/downstream_composition.R",     base = modules_root)
  src("pn_core.R")
  src("pn_enrichment.R")
  src("pn_network.R")
  src("pn_ui.R")
})

# Null-coalescing operator
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# ---------------------------------------------------------------------------
# Plotly palette — mirrors app.css design tokens. Plotly takes raw hex/rgba,
# not CSS variables, so we name them here in one place. Any plot helper that
# needs a color references .PN_PAL[[...]] instead of inlining a hex literal.
# ---------------------------------------------------------------------------
.PN_PAL <- list(
  primary       = "#2563EB",   # --primary
  primary_dark  = "#1E3A8A",   # darker shade for label/text
  danger        = "#DC2626",   # --danger
  danger_dark   = "#7F1D1D",
  success       = "#16A34A",   # --success
  success_dark  = "#14532D",
  text          = "#0F172A",   # --text-primary
  font_family   = "Inter"
)
#' Build an rgba() string from a hex + alpha (helper for fillcolor).
.pn_rgba <- function(hex, alpha = 0.30) {
  rgb_vals <- grDevices::col2rgb(hex)
  sprintf("rgba(%d,%d,%d,%.2f)", rgb_vals[1, 1], rgb_vals[2, 1], rgb_vals[3, 1], alpha)
}

pnServer <- function(id, values) {
  shiny::moduleServer(id, function(input, output, session) {

    # ----- Local reactive state -----
    enrichment_results <- shiny::reactiveVal(NULL)
    enrichment_audit   <- shiny::reactiveVal(NULL)
    show_enrich        <- shiny::reactiveVal(FALSE)

    pwn_results        <- shiny::reactiveVal(NULL)
    show_pwn           <- shiny::reactiveVal(FALSE)

    grn_results        <- shiny::reactiveVal(NULL)
    show_grn           <- shiny::reactiveVal(FALSE)

    coexp_results      <- shiny::reactiveVal(NULL)
    show_coexp         <- shiny::reactiveVal(FALSE)

    cmp_results        <- shiny::reactiveVal(NULL)
    cmp_audit          <- shiny::reactiveVal(NULL)
    show_cmp           <- shiny::reactiveVal(FALSE)

    # ----- Mount delegated sub-modules (CCC / TF / Geneset / Composition) -----
    # Each mount is wrapped so a failure in one doesn't take down the others.
    # On error we surface a Shiny notification AND write to the console so the
    # user sees something concrete rather than a silent blank panel.
    .pn_safe_mount <- function(name, fn) {
      tryCatch(fn(), error = function(e) {
        msg <- sprintf("[pn] %s init failed: %s", name, conditionMessage(e))
        message(msg)
        shiny::showNotification(msg, type = "error", duration = 12)
      })
    }
    if (exists("cellCommServer", mode = "function"))
      .pn_safe_mount("cellCommServer",
                     function() cellCommServer("ccc", values))
    if (exists("tfActivityServer", mode = "function"))
      .pn_safe_mount("tfActivityServer",
                     function() tfActivityServer("tf", values))
    if (exists("geneSetScoringServer", mode = "function"))
      .pn_safe_mount("geneSetScoringServer",
                     function() geneSetScoringServer("geneset", values))
    if (exists("compositionServer", mode = "function"))
      .pn_safe_mount("compositionServer",
                     function() compositionServer("composition", values))

    # =====================================================================
    # 1. PATHWAY ENRICHMENT
    # =====================================================================
    shiny::observeEvent(input$run_enrich, {
      de_data <- values$last_de_results
      if (is.null(de_data) || !is.data.frame(de_data) || nrow(de_data) == 0) {
        shiny::showNotification(
          "No DE results available — run Differential Expression first.",
          type = "error", duration = 8)
        return()
      }
      .sp_resolved <- pn_resolve_species(values)
      species <- .sp_resolved$species
      if (isTRUE(.sp_resolved$fell_back))
        shiny::showNotification(
          sprintf("No species set in project settings — falling back to '%s'. Set it for correct results on non-human data.",
                  species),
          type = "warning", duration = 8)

      withProgress_msg <- "Loading gene sets..."
      shiny::withProgress(message = withProgress_msg, value = 0.05, {

        # 1. Load gene sets (or sentinel for live API)
        gs <- tryCatch(
          load_gene_sets(database = input$enrich_db, species = species),
          error = function(e) e)
        if (inherits(gs, "error")) {
          shiny::showNotification(
            paste("Gene-set load failed:", conditionMessage(gs)),
            type = "error", duration = 10)
          return()
        }

        shiny::incProgress(0.2, detail = "Preparing gene lists...")
        custom_bg <- if (isTRUE(input$use_custom_bg) &&
                         !is.null(input$custom_bg) &&
                         nzchar(trimws(input$custom_bg))) {
          trimws(unlist(strsplit(input$custom_bg, "[,\n\\s]+")))
        } else NULL
        bg <- pn_background_universe(values$active_dataset, custom = custom_bg)

        gl <- prepare_gene_lists(
          de_data,
          ora_pval = input$ora_pval_threshold %||% 0.05,
          ora_lfc  = input$ora_logfc_threshold %||% 0.25,
          gsea_metric = input$gsea_metric %||% "signed_neglog10p",
          universe = bg,
          seed = as.integer(input$enrich_seed %||% 42L))

        params <- list(
          ora_pval_cutoff = input$ora_pval_cutoff %||% 0.05,
          gsea_pval_cutoff = 0.05,
          gsea_min_size = input$gsea_min_size %||% 15L,
          gsea_max_size = input$gsea_max_size %||% 500L,
          gsea_nperm    = input$gsea_nperm %||% 1000L,
          max_pathways  = input$max_pathways %||% 30L,
          seed          = as.integer(input$enrich_seed %||% 42L)
        )

        out <- list(ora = NULL, gsea = NULL)
        if (input$enrich_method %in% c("ORA", "BOTH")) {
          shiny::incProgress(0.25, detail = "Running ORA...")
          out$ora <- tryCatch(run_ora_analysis(gl, gs, params),
                              error = function(e) {
                                list(up = NULL, down = NULL,
                                     error = conditionMessage(e))
                              })
        }
        if (input$enrich_method %in% c("GSEA", "BOTH")) {
          shiny::incProgress(0.5, detail = "Running GSEA...")
          out$gsea <- tryCatch(run_gsea_analysis(gl$ranked_genes, gs, params),
                               error = function(e) NULL)
        }

        # Stash results for cross-module + Pathway Network sub-tab
        enrichment_results(out)
        values$last_enrich_results <- out
        n_results <- (if (!is.null(out$ora$up))   nrow(out$ora$up)   else 0L) +
                     (if (!is.null(out$ora$down)) nrow(out$ora$down) else 0L) +
                     (if (!is.null(out$gsea))     nrow(out$gsea)     else 0L)
        enrichment_audit(build_enrichment_audit(
          database = input$enrich_db, method = input$enrich_method,
          species  = species,
          n_input  = length(gl$up_genes) + length(gl$down_genes),
          n_background = length(gl$universe),
          seed     = as.integer(input$enrich_seed %||% 42L),
          n_results = n_results))
        show_enrich(TRUE)
        # --- Auto-push figures + CSVs to gallery ---
        .pn_push_enrichment_to_gallery(values, out, input$enrich_db,
                                       input$enrich_method,
                                       input$max_pathways %||% 30L)
        shiny::incProgress(1, detail = "Done.")
        shiny::showNotification(
          sprintf("Enrichment complete -> %d total results across %s. Saved to Figure Gallery.",
                  n_results, input$enrich_db),
          type = "message", duration = 5)
      })
    })

    output$enrich_results_ui <- shiny::renderUI({
      shiny::req(isTRUE(show_enrich()))
      ns <- session$ns
      shiny::tabsetPanel(
        shiny::tabPanel("Dot plot",
          shiny::br(), plotly::plotlyOutput(ns("enrich_dotplot"), height = "560px")),
        shiny::tabPanel("Bar plot",
          shiny::br(), plotly::plotlyOutput(ns("enrich_barplot"), height = "560px")),
        shiny::tabPanel("GSEA",
          shiny::br(), plotly::plotlyOutput(ns("gsea_plot"), height = "560px")),
        shiny::tabPanel("Tables",
          shiny::br(),
          shiny::h5("ORA results"),
          DT::DTOutput(ns("ora_table")),
          shiny::br(),
          shiny::h5("GSEA results"),
          DT::DTOutput(ns("gsea_table")),
          shiny::br(),
          shiny::downloadButton(ns("download_ora_csv"),  "Download ORA (CSV)",
            class = "btn-success"),
          shiny::downloadButton(ns("download_gsea_csv"), "Download GSEA (CSV)",
            class = "btn-success")),
        shiny::tabPanel("Audit Log",
          shiny::br(), DT::DTOutput(ns("enrich_audit_table")))
      )
    })

    output$enrich_dotplot <- plotly::renderPlotly({
      shiny::req(enrichment_results())
      .pn_render_dotplot(enrichment_results(), input$max_pathways %||% 30L)
    })
    output$enrich_barplot <- plotly::renderPlotly({
      shiny::req(enrichment_results())
      .pn_render_barplot(enrichment_results(), input$max_pathways %||% 30L)
    })
    output$gsea_plot <- plotly::renderPlotly({
      shiny::req(enrichment_results())
      .pn_render_gsea(enrichment_results())
    })
    output$ora_table <- DT::renderDT({
      shiny::req(enrichment_results())
      .pn_combined_ora_df(enrichment_results())
    }, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    output$gsea_table <- DT::renderDT({
      shiny::req(enrichment_results())
      df <- enrichment_results()$gsea
      if (is.null(df) || nrow(df) == 0)
        return(data.frame(message = "No GSEA results"))
      df
    }, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    output$enrich_audit_table <- DT::renderDT({
      shiny::req(enrichment_audit())
      a <- enrichment_audit()
      data.frame(
        Field = names(a),
        Value = vapply(a, function(v) {
          if (length(v) == 0) "(none)"
          else if (length(v) > 1) paste(v, collapse = ", ")
          else as.character(v)
        }, character(1)),
        stringsAsFactors = FALSE)
    }, options = list(pageLength = 20, dom = "t"), rownames = FALSE)

    output$download_ora_csv <- shiny::downloadHandler(
      filename = function() sprintf("ora_%s.csv",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        .pn_combined_ora_df(enrichment_results()), file, row.names = FALSE))
    output$download_gsea_csv <- shiny::downloadHandler(
      filename = function() sprintf("gsea_%s.csv",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        df <- enrichment_results()$gsea
        utils::write.csv(if (is.null(df)) data.frame() else df,
                         file, row.names = FALSE)
      })

    # =====================================================================
    # 1b. COMPARATIVE ENRICHMENT
    # =====================================================================

    # Live-update the picker + history strip from values$de_history.
    shiny::observe({
      hist <- values$de_history %||% list()
      ids  <- names(hist)
      shinyWidgets::updatePickerInput(
        session = session, inputId = "cmp_de_runs",
        choices = ids,
        selected = intersect(input$cmp_de_runs, ids))
    })

    output$cmp_history_strip <- shiny::renderUI({
      hist <- values$de_history %||% list()
      if (length(hist) == 0) {
        return(shiny::div(class = "pn-history-empty",
          shiny::icon("info-circle"),
          " No saved DE runs yet — run a few DE comparisons in the ",
          shiny::strong("Differential Expression"),
          " tab. Each successful run is auto-saved."))
      }
      # Show count + recent chips (last 8)
      recent <- utils::tail(hist, 8L)
      shiny::tagList(
        shiny::div(class = "pn-history-count",
          sprintf("%d saved run%s · showing most recent %d",
                  length(hist),
                  if (length(hist) == 1) "" else "s",
                  length(recent))),
        shiny::div(class = "pn-history-chips",
          lapply(rev(names(recent)), function(nm)
            shiny::span(class = "pn-history-chip",
              shiny::icon("bookmark"), shiny::span(nm))))
      )
    })

    shiny::observeEvent(input$run_compare_enrich, {
      picked <- input$cmp_de_runs
      hist   <- values$de_history %||% list()
      if (length(picked) < 2) {
        shiny::showNotification(
          "Pick at least 2 saved DE runs to compare.",
          type = "error", duration = 8); return()
      }
      de_runs <- hist[picked]
      .sp_resolved <- pn_resolve_species(values)
      species <- .sp_resolved$species
      if (isTRUE(.sp_resolved$fell_back))
        shiny::showNotification(
          sprintf("No species set in project settings — falling back to '%s'. Set it for correct results on non-human data.",
                  species),
          type = "warning", duration = 8)
      bg <- pn_background_universe(values$active_dataset, custom = NULL)

      params <- list(
        ora_pval         = input$cmp_ora_pval %||% 0.05,
        ora_lfc          = input$cmp_ora_lfc  %||% 0.25,
        ora_pval_cutoff  = input$cmp_ora_pval_cutoff %||% 0.05,
        gsea_pval_cutoff = 0.05,
        gsea_metric      = "signed_neglog10p",
        gsea_min_size    = input$cmp_gsea_min   %||% 15L,
        gsea_max_size    = input$cmp_gsea_max   %||% 500L,
        gsea_nperm       = input$cmp_gsea_nperm %||% 1000L,
        max_pathways     = input$cmp_top_n %||% 20L,
        universe         = bg,
        seed             = as.integer(input$cmp_seed %||% 42L)
      )

      shiny::withProgress(message = "Running comparative enrichment...",
                          value = 0.1, {
        res <- tryCatch(
          run_compare_enrichment(de_runs = de_runs,
                                 database = input$cmp_db,
                                 species  = species,
                                 method   = input$cmp_method %||% "BOTH",
                                 params   = params),
          error = function(e) {
            shiny::showNotification(
              paste("Comparative enrichment failed:", conditionMessage(e)),
              type = "error", duration = 12); NULL
          })
        shiny::incProgress(0.9)
        if (is.null(res)) { show_cmp(FALSE); return() }

        cmp_results(res); cmp_audit(res$audit); show_cmp(TRUE)
        values$last_compare_enrich <- res
        n_ora  <- if (is.null(res$combined_ora))  0L else nrow(res$combined_ora)
        n_gsea <- if (is.null(res$combined_gsea)) 0L else nrow(res$combined_gsea)
        # --- Auto-push to gallery ---
        .pn_push_compare_to_gallery(values, res,
                                    input$cmp_db, input$cmp_method %||% "BOTH",
                                    input$cmp_top_n %||% 20L)
        shiny::showNotification(
          sprintf("Comparative enrichment complete -> %d ORA + %d GSEA rows across %d conditions. Saved to Figure Gallery.",
                  n_ora, n_gsea, length(picked)),
          type = "message", duration = 6)
      })
    })

    output$cmp_results_ui <- shiny::renderUI({
      if (!isTRUE(show_cmp())) {
        return(shiny::div(class = "pn-empty",
          shiny::div(class = "icon", shiny::icon("balance-scale")),
          shiny::h4("Comparative results will appear here"),
          shiny::p("Pick 2+ saved DE runs from the history on the left, ",
                   "choose a database and method, and click ",
                   shiny::strong("Run Comparative Enrichment"),
                   ". You'll see a side-by-side dotplot, a Venn / UpSet of ",
                   "overlapping significant pathways, a combined long-form ",
                   "table you can sort + download, and an audit log.")))
      }
      ns <- session$ns
      shiny::div(class = "pn-results-card",
        shiny::uiOutput(ns("cmp_metrics")),
        shiny::tabsetPanel(
          shiny::tabPanel("Dotplot",
            shiny::br(),
            plotly::plotlyOutput(ns("cmp_dotplot"), height = "640px")),
          shiny::tabPanel("Overlap (Venn / UpSet)",
            shiny::br(),
            plotly::plotlyOutput(ns("cmp_overlap"), height = "440px"),
            shiny::br(),
            DT::DTOutput(ns("cmp_overlap_table"))),
          shiny::tabPanel("Combined ORA",
            shiny::br(),
            DT::DTOutput(ns("cmp_ora_table")),
            shiny::br(),
            shiny::downloadButton(ns("cmp_dl_ora"),
              "Download (CSV)", class = "btn-success")),
          shiny::tabPanel("Combined GSEA",
            shiny::br(),
            DT::DTOutput(ns("cmp_gsea_table")),
            shiny::br(),
            shiny::downloadButton(ns("cmp_dl_gsea"),
              "Download (CSV)", class = "btn-success")),
          shiny::tabPanel("Audit Log",
            shiny::br(),
            DT::DTOutput(ns("cmp_audit_table")))
        )
      )
    })

    output$cmp_metrics <- shiny::renderUI({
      shiny::req(cmp_results(), cmp_audit())
      a <- cmp_audit()
      res <- cmp_results()
      n_sig_per <- vapply(res$sig_pathways,
        function(x) length(x$union), integer(1))
      ora_rows <- if (is.null(res$combined_ora)) 0 else nrow(res$combined_ora)
      gsea_rows <- if (is.null(res$combined_gsea)) 0 else nrow(res$combined_gsea)
      .cmp_metric <- function(label, value, sub = "") {
        shiny::div(class = "pn-cmp-metric",
          shiny::span(class = "label", label),
          shiny::div(class = "value", value),
          if (nzchar(sub)) shiny::div(class = "sub", sub))
      }
      shiny::div(class = "pn-cmp-metrics",
        .cmp_metric("Conditions",  format(a$n_conditions, big.mark = ","),
                    "compared side-by-side"),
        .cmp_metric("Database",    a$database,
                    sprintf("%s · %s", a$method, a$species)),
        .cmp_metric("ORA rows",    format(ora_rows, big.mark = ",")),
        .cmp_metric("GSEA rows",   format(gsea_rows, big.mark = ",")),
        .cmp_metric("Sig pathways (min)",
                    if (length(n_sig_per) > 0) min(n_sig_per) else 0,
                    "lowest across conditions"),
        .cmp_metric("Sig pathways (max)",
                    if (length(n_sig_per) > 0) max(n_sig_per) else 0,
                    "highest across conditions")
      )
    })

    output$cmp_dotplot <- plotly::renderPlotly({
      shiny::req(cmp_results())
      .render_compare_dotplot(cmp_results(),
                              top_n = input$cmp_top_n %||% 20L,
                              method = input$cmp_method %||% "BOTH")
    })
    output$cmp_overlap <- plotly::renderPlotly({
      shiny::req(cmp_results())
      .render_compare_overlap(cmp_results())
    })
    output$cmp_overlap_table <- DT::renderDT({
      shiny::req(cmp_results())
      .build_overlap_table(cmp_results())
    }, options = list(pageLength = 12, scrollX = TRUE), rownames = FALSE)
    output$cmp_ora_table <- DT::renderDT({
      shiny::req(cmp_results())
      df <- cmp_results()$combined_ora
      if (is.null(df) || nrow(df) == 0)
        return(data.frame(message = "No ORA results"))
      keep <- intersect(c("condition","Direction","ID","Description",
                          "GeneRatio","BgRatio","pvalue","p.adjust","Count","geneID"),
                        colnames(df))
      df[, keep, drop = FALSE]
    }, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    output$cmp_gsea_table <- DT::renderDT({
      shiny::req(cmp_results())
      df <- cmp_results()$combined_gsea
      if (is.null(df) || nrow(df) == 0)
        return(data.frame(message = "No GSEA results"))
      df
    }, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    output$cmp_audit_table <- DT::renderDT({
      shiny::req(cmp_audit())
      a <- cmp_audit()
      data.frame(
        Field = names(a),
        Value = vapply(a, function(v) {
          if (length(v) == 0) "(none)"
          else if (length(v) > 1) paste(v, collapse = ", ")
          else as.character(v)
        }, character(1)),
        stringsAsFactors = FALSE)
    }, options = list(pageLength = 20, dom = "t"), rownames = FALSE)

    output$cmp_dl_ora <- shiny::downloadHandler(
      filename = function() sprintf("compare_ora_%s.csv",
        format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        df <- cmp_results()$combined_ora
        utils::write.csv(if (is.null(df)) data.frame() else df,
                         file, row.names = FALSE)
      })
    output$cmp_dl_gsea <- shiny::downloadHandler(
      filename = function() sprintf("compare_gsea_%s.csv",
        format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) {
        df <- cmp_results()$combined_gsea
        utils::write.csv(if (is.null(df)) data.frame() else df,
                         file, row.names = FALSE)
      })

    # =====================================================================
    # 2. PATHWAY NETWORK
    # =====================================================================
    output$pwn_status <- shiny::renderText({
      if (!is.null(values$last_enrich_results) &&
          !is.null(values$last_enrich_results$ora) &&
          (!is.null(values$last_enrich_results$ora$up) ||
           !is.null(values$last_enrich_results$ora$down))) {
        "Enrichment results loaded. Set thresholds and build."
      } else {
        "No enrichment results yet — run Pathway Enrichment (ORA) first."
      }
    })

    shiny::observeEvent(input$run_pwn, {
      enrich <- values$last_enrich_results %||% enrichment_results()
      if (is.null(enrich) || is.null(enrich$ora) ||
          (is.null(enrich$ora$up) && is.null(enrich$ora$down))) {
        shiny::showNotification(
          "Run Pathway Enrichment (ORA) before building a pathway network.",
          type = "error", duration = 8)
        return()
      }
      shiny::withProgress(message = "Building pathway network...", value = 0.3, {
        res <- build_pathway_network(enrich, list(
          pathway_top_n = input$pwn_top_n %||% 30L,
          pathway_overlap_threshold = input$pwn_overlap %||% 3L,
          pathway_jaccard_threshold = input$pwn_jaccard %||% 0.1
        ))
        shiny::incProgress(0.7)
        if (isTRUE(res$success)) {
          pwn_results(res); show_pwn(TRUE)
          # --- Auto-push to gallery ---
          .pn_push_pwn_to_gallery(values, res)
          shiny::showNotification(
            "Pathway network built. Saved to Figure Gallery.",
            type = "message")
        } else {
          show_pwn(FALSE)
          shiny::showNotification(paste("Pathway network failed:", res$error),
            type = "error", duration = 10)
        }
      })
    })

    output$pwn_results_ui <- shiny::renderUI({
      shiny::req(isTRUE(show_pwn()))
      ns <- session$ns
      shiny::tagList(
        shiny::fluidRow(
          shiny::column(9, visNetwork::visNetworkOutput(ns("pwn_network"),
            height = "600px")),
          shiny::column(3,
            shiny::h5("Pathway clusters"),
            DT::DTOutput(ns("pwn_clusters")),
            shiny::br(),
            shiny::downloadButton(ns("download_pwn_csv"),
              "Download Edges (CSV)", class = "btn-success"))
        )
      )
    })

    output$pwn_network <- visNetwork::renderVisNetwork({
      shiny::req(pwn_results())
      render_pathway_network(pwn_results())
    })
    output$pwn_clusters <- DT::renderDT({
      shiny::req(pwn_results())
      pwn_results()$clusters
    }, options = list(pageLength = 15, dom = "t"), rownames = FALSE)
    output$download_pwn_csv <- shiny::downloadHandler(
      filename = function() sprintf("pathway_network_%s.csv",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        pwn_results()$network_data, file, row.names = FALSE))

    # =====================================================================
    # 3. GENE REGULATORY NETWORKS
    # =====================================================================
    shiny::observe({
      shiny::req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[vapply(meta,
        function(x) is.factor(x) || is.character(x), logical(1))]
      ct_col <- intersect(c("cell_type","celltype","seurat_clusters","cluster"),
                          cat_cols)
      ct_col <- if (length(ct_col) > 0) ct_col[1]
                else if (length(cat_cols) > 0) cat_cols[1] else NULL
      if (!is.null(ct_col)) {
        cts <- sort(unique(as.character(meta[[ct_col]])))
        cts <- cts[!is.na(cts)]
        shiny::updateSelectInput(session, "grn_cell_type",
          choices = c("All cells" = "all", cts))
      }
    })

    shiny::observeEvent(input$run_grn, {
      shiny::req(values$active_dataset)
      shiny::withProgress(message = "Building GRN...", value = 0.2, {
        params <- list(
          grn_method = input$grn_method,
          grn_top_genes = input$grn_top_genes %||% 2000L,
          grn_min_cells = input$grn_min_cells %||% 10L,
          grn_cell_type = input$grn_cell_type %||% "all",
          wgcna_power = input$wgcna_power %||% 0L,
          wgcna_min_module_size = input$wgcna_min_module_size %||% 30L,
          wgcna_merge_threshold = input$wgcna_merge_threshold %||% 0.25,
          wgcna_tom_type = input$wgcna_tom_type %||% "signed",
          wgcna_detect_outliers = isTRUE(input$wgcna_detect_outliers),
          cor_method = input$cor_method %||% "spearman",
          cor_threshold = input$cor_threshold %||% 0.7,
          cor_p_threshold = input$cor_p_threshold %||% 0.05,
          seed = 42L
        )
        res <- build_gene_regulatory_network(values$active_dataset, params)
        shiny::incProgress(0.8)
        if (isTRUE(res$success)) {
          grn_results(res); show_grn(TRUE)
          # --- Auto-push to gallery ---
          .pn_push_grn_to_gallery(values, res)
          shiny::showNotification(
            sprintf("GRN built: %d edges across %d genes (method: %s). Saved to Figure Gallery.",
                    igraph::ecount(res$network), igraph::vcount(res$network),
                    res$method), type = "message", duration = 5)
        } else {
          show_grn(FALSE)
          shiny::showNotification(paste("GRN failed:", res$error),
            type = "error", duration = 10)
        }
      })
    })

    output$grn_results_ui <- shiny::renderUI({
      shiny::req(isTRUE(show_grn()))
      ns <- session$ns
      shiny::tagList(
        shiny::tabsetPanel(
          shiny::tabPanel("Network",
            shiny::br(),
            visNetwork::visNetworkOutput(ns("grn_network"), height = "640px")),
          shiny::tabPanel("Modules",
            shiny::br(),
            DT::DTOutput(ns("grn_modules_table"))),
          shiny::tabPanel("Edges",
            shiny::br(),
            DT::DTOutput(ns("grn_edges_table")))
        ),
        shiny::br(),
        shiny::downloadButton(ns("download_grn_edges_csv"),
          "Download Edges (CSV)", class = "btn-success"),
        shiny::downloadButton(ns("download_grn_modules_csv"),
          "Download Modules (CSV)", class = "btn-info")
      )
    })

    output$grn_network <- visNetwork::renderVisNetwork({
      shiny::req(grn_results())
      render_grn_network(grn_results())
    })
    output$grn_modules_table <- DT::renderDT({
      shiny::req(grn_results())
      df <- grn_results()$modules
      if (is.null(df)) data.frame(message = "(method has no module info)") else df
    }, options = list(pageLength = 15, scrollX = TRUE), rownames = FALSE)
    output$grn_edges_table <- DT::renderDT({
      shiny::req(grn_results())
      g <- grn_results()$network
      igraph::as_data_frame(g, what = "edges")
    }, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    output$download_grn_edges_csv <- shiny::downloadHandler(
      filename = function() sprintf("grn_edges_%s.csv",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        igraph::as_data_frame(grn_results()$network, what = "edges"),
        file, row.names = FALSE))
    output$download_grn_modules_csv <- shiny::downloadHandler(
      filename = function() sprintf("grn_modules_%s.csv",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        grn_results()$modules %||% data.frame(), file, row.names = FALSE))

    # =====================================================================
    # 4. CO-EXPRESSION
    # =====================================================================
    shiny::observeEvent(input$run_coexp, {
      shiny::req(values$active_dataset)
      shiny::withProgress(message = "Building co-expression network...",
                          value = 0.4, {
        res <- build_coexpression_network(values$active_dataset, list(
          coexp_genes = input$coexp_genes %||% "",
          coexp_method = input$coexp_method %||% "spearman",
          coexp_threshold = input$coexp_threshold %||% 0.3,
          coexp_top_n = input$coexp_top_n %||% 200L
        ))
        shiny::incProgress(0.6)
        if (isTRUE(res$success)) {
          coexp_results(res); show_coexp(TRUE)
          # --- Auto-push to gallery ---
          .pn_push_coexp_to_gallery(values, res)
          shiny::showNotification(
            "Co-expression network built. Saved to Figure Gallery.",
            type = "message", duration = 4)
        } else {
          show_coexp(FALSE)
          shiny::showNotification(paste("Co-expression failed:", res$error),
            type = "error", duration = 10)
        }
      })
    })

    output$coexp_results_ui <- shiny::renderUI({
      shiny::req(isTRUE(show_coexp()))
      ns <- session$ns
      shiny::tagList(
        shiny::tabsetPanel(
          shiny::tabPanel("Network",
            shiny::br(),
            visNetwork::visNetworkOutput(ns("coexp_network"), height = "560px")),
          shiny::tabPanel("Correlation matrix",
            shiny::br(),
            plotly::plotlyOutput(ns("coexp_matrix"), height = "560px")),
          shiny::tabPanel("Edges",
            shiny::br(), DT::DTOutput(ns("coexp_edges_table")))
        ),
        shiny::br(),
        shiny::downloadButton(ns("download_coexp_csv"),
          "Download Edges (CSV)", class = "btn-success")
      )
    })

    output$coexp_network <- visNetwork::renderVisNetwork({
      shiny::req(coexp_results())
      render_coexp_network(coexp_results())
    })
    output$coexp_matrix <- plotly::renderPlotly({
      shiny::req(coexp_results())
      render_correlation_matrix(coexp_results())
    })
    output$coexp_edges_table <- DT::renderDT({
      shiny::req(coexp_results())
      coexp_results()$edges
    }, options = list(pageLength = 25, scrollX = TRUE), rownames = FALSE)
    output$download_coexp_csv <- shiny::downloadHandler(
      filename = function() sprintf("coexp_%s.csv",
                                    format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        coexp_results()$edges, file, row.names = FALSE))
  })
}

# ---- Plot helpers (file-local) ----------------------------------------------
.pn_combined_ora_df <- function(results) {
  if (is.null(results) || is.null(results$ora)) return(data.frame())
  do.call(rbind, c(
    list(if (!is.null(results$ora$up))   results$ora$up),
    list(if (!is.null(results$ora$down)) results$ora$down)
  ))
}

.pn_render_dotplot <- function(results, top_n = 20L) {
  df <- .pn_combined_ora_df(results)
  if (is.null(df) || nrow(df) == 0)
    return(plotly::plotly_empty() |>
           plotly::layout(title = "No ORA results"))
  df <- utils::head(df[order(df$p.adjust), , drop = FALSE], top_n)
  df$log_p <- -log10(pmax(as.numeric(df$p.adjust), 1e-300))
  if ("GeneRatio" %in% colnames(df)) {
    gr <- vapply(as.character(df$GeneRatio), function(x) {
      parts <- strsplit(x, "/")[[1]]
      if (length(parts) == 2) {
        n <- suppressWarnings(as.numeric(parts[1]))
        d <- suppressWarnings(as.numeric(parts[2]))
        if (!is.na(n) && !is.na(d) && d != 0) return(n / d)
      }
      val <- suppressWarnings(as.numeric(x))
      if (is.na(val)) 0 else val
    }, numeric(1))
    df$gene_ratio <- gr
  } else df$gene_ratio <- df$Count / max(df$Count, 1)
  df$desc <- ifelse(nchar(as.character(df$Description)) > 55,
                    paste0(substr(df$Description, 1, 52), "..."),
                    as.character(df$Description))
  df <- df[order(df$log_p), ]
  plotly::plot_ly(
    df, x = ~gene_ratio, y = ~desc, color = ~log_p,
    size = ~Count, type = "scatter", mode = "markers",
    text = ~paste0(Description,
                   "<br>Count: ", Count,
                   "<br>adj.p: ", signif(p.adjust, 3),
                   "<br>Direction: ", Direction),
    hoverinfo = "text",
    marker = list(sizemode = "diameter", sizeref = max(1, max(df$Count) / 30))
  ) |>
    plotly::layout(
      title = "Top enriched pathways (ORA)",
      xaxis = list(title = "Gene ratio"),
      yaxis = list(title = "", categoryorder = "array", categoryarray = df$desc))
}

.pn_render_barplot <- function(results, top_n = 20L) {
  df <- .pn_combined_ora_df(results)
  if (is.null(df) || nrow(df) == 0)
    return(plotly::plotly_empty() |>
           plotly::layout(title = "No ORA results"))
  df <- utils::head(df[order(df$p.adjust), , drop = FALSE], top_n)
  df$log_p <- -log10(pmax(as.numeric(df$p.adjust), 1e-300))
  df$desc <- ifelse(nchar(as.character(df$Description)) > 55,
                    paste0(substr(df$Description, 1, 52), "..."),
                    as.character(df$Description))
  df <- df[order(df$log_p), ]
  plotly::plot_ly(df, x = ~log_p, y = ~desc, color = ~Direction,
                  type = "bar", orientation = "h") |>
    plotly::layout(
      title = "Top enriched pathways (ORA)",
      xaxis = list(title = "-log10 adj p"),
      yaxis = list(title = "", categoryorder = "array", categoryarray = df$desc))
}

.pn_render_gsea <- function(results) {
  df <- results$gsea
  if (is.null(df) || nrow(df) == 0)
    return(plotly::plotly_empty() |>
           plotly::layout(title = "No GSEA results"))
  pad_col <- if ("padj" %in% colnames(df)) "padj"
             else if ("p.adjust" %in% colnames(df)) "p.adjust" else NA_character_
  pw_col  <- if ("pathway" %in% colnames(df)) "pathway"
             else if ("Description" %in% colnames(df)) "Description" else NA_character_
  if (is.na(pad_col) || is.na(pw_col) || !"NES" %in% colnames(df))
    return(plotly::plotly_empty() |>
           plotly::layout(title = "GSEA result schema unrecognised"))
  d <- df[order(as.numeric(df[[pad_col]])), , drop = FALSE]
  d <- utils::head(d, 30)
  d$desc <- ifelse(nchar(as.character(d[[pw_col]])) > 55,
                   paste0(substr(d[[pw_col]], 1, 52), "..."),
                   as.character(d[[pw_col]]))
  d$log_p <- -log10(pmax(as.numeric(d[[pad_col]]), 1e-300))
  d <- d[order(d$NES), ]
  plotly::plot_ly(d, x = ~NES, y = ~desc, color = ~NES,
                  size = ~log_p, type = "scatter", mode = "markers",
                  hoverinfo = "text",
                  text = ~paste0(.data[[pw_col]],
                                 "<br>NES: ", round(NES, 3),
                                 "<br>adj.p: ", signif(.data[[pad_col]], 3)),
                  marker = list(sizemode = "diameter",
                                sizeref = max(1, max(d$log_p) / 30))) |>
    plotly::layout(
      title = "GSEA — NES vs pathway",
      xaxis = list(title = "Normalised Enrichment Score (NES)"),
      yaxis = list(title = "", categoryorder = "array", categoryarray = d$desc))
}

# ---------------------------------------------------------------------------
# Comparative Enrichment plot helpers (rendered in pn_server's renderers).
# All take the run_compare_enrichment() return value.
# ---------------------------------------------------------------------------

#' Side-by-side dotplot: rows = top pathways, cols = conditions.
#' Dot color = adj p (darker = more significant), dot size = gene count (ORA)
#' or |NES| (GSEA). Top pathways picked by smallest adj p across all conditions.
.render_compare_dotplot <- function(res, top_n = 20L, method = "BOTH") {
  use_gsea <- method == "GSEA" ||
              (method == "BOTH" && !is.null(res$combined_gsea) &&
               nrow(res$combined_gsea) > 0)
  if (use_gsea && !is.null(res$combined_gsea) && nrow(res$combined_gsea) > 0) {
    df <- res$combined_gsea
    pad_col <- if ("padj" %in% colnames(df)) "padj"
               else if ("p.adjust" %in% colnames(df)) "p.adjust" else NA_character_
    pw_col  <- if ("pathway" %in% colnames(df)) "pathway"
               else if ("Description" %in% colnames(df)) "Description"
               else if ("ID" %in% colnames(df)) "ID" else NA_character_
    if (is.na(pad_col) || is.na(pw_col) || !"condition" %in% colnames(df))
      return(plotly::plotly_empty() |>
             plotly::layout(title = "GSEA result has unexpected schema"))
    df$padj_n <- as.numeric(df[[pad_col]])
    df$pw     <- as.character(df[[pw_col]])
    df$cond   <- as.character(df$condition)
    df$abs_nes <- if ("NES" %in% colnames(df)) abs(as.numeric(df$NES)) else 1
    # Pick top-N pathways by minimum padj across conditions
    pw_min <- tapply(df$padj_n, df$pw, min, na.rm = TRUE)
    top_pw <- names(sort(pw_min))[seq_len(min(top_n, length(pw_min)))]
    df <- df[df$pw %in% top_pw, , drop = FALSE]
    df$pw_short <- ifelse(nchar(df$pw) > 50,
                          paste0(substr(df$pw, 1, 47), "..."), df$pw)
    df$log_p <- -log10(pmax(df$padj_n, 1e-300))
    return(
      plotly::plot_ly(df,
        x = ~cond, y = ~pw_short,
        size = ~abs_nes, color = ~log_p,
        type = "scatter", mode = "markers",
        marker = list(sizemode = "diameter",
                      sizeref = max(1, max(df$abs_nes, na.rm = TRUE) / 18),
                      line = list(color = "white", width = 0.5)),
        text = ~sprintf("Pathway: %s<br>Condition: %s<br>NES: %s<br>adj.p: %s",
                        pw, cond,
                        if ("NES" %in% colnames(df))
                          signif(as.numeric(df$NES), 3) else NA,
                        signif(padj_n, 3)),
        hoverinfo = "text"
      ) |>
        plotly::colorbar(title = "-log10\nadj p") |>
        plotly::layout(
          title = sprintf("Top %d GSEA pathways across %d conditions",
                          length(top_pw), length(unique(df$cond))),
          xaxis = list(title = "", tickangle = -30),
          yaxis = list(title = "", automargin = TRUE),
          margin = list(l = 200, b = 100)
        )
    )
  }

  # ORA branch
  df <- res$combined_ora
  if (is.null(df) || nrow(df) == 0)
    return(plotly::plotly_empty() |>
           plotly::layout(title = "No ORA results to plot"))
  if (!all(c("Description","p.adjust","Count","condition") %in% colnames(df)))
    return(plotly::plotly_empty() |>
           plotly::layout(title = "ORA result has unexpected schema"))
  df$padj_n <- as.numeric(df$p.adjust)
  df$pw     <- as.character(df$Description)
  df$cond   <- as.character(df$condition)
  if ("Direction" %in% colnames(df)) {
    df$cond <- paste0(df$cond, " · ", df$Direction)
  }
  pw_min <- tapply(df$padj_n, df$pw, min, na.rm = TRUE)
  top_pw <- names(sort(pw_min))[seq_len(min(top_n, length(pw_min)))]
  df <- df[df$pw %in% top_pw, , drop = FALSE]
  df$pw_short <- ifelse(nchar(df$pw) > 50,
                        paste0(substr(df$pw, 1, 47), "..."), df$pw)
  df$log_p <- -log10(pmax(df$padj_n, 1e-300))
  df$count_n <- as.numeric(df$Count)
  plotly::plot_ly(df,
    x = ~cond, y = ~pw_short,
    size = ~count_n, color = ~log_p,
    type = "scatter", mode = "markers",
    marker = list(sizemode = "diameter",
                  sizeref = max(1, max(df$count_n, na.rm = TRUE) / 28),
                  line = list(color = "white", width = 0.5)),
    text = ~sprintf("Pathway: %s<br>Condition: %s<br>Count: %d<br>adj.p: %s",
                    pw, cond, as.integer(count_n), signif(padj_n, 3)),
    hoverinfo = "text"
  ) |>
    plotly::colorbar(title = "-log10\nadj p") |>
    plotly::layout(
      title = sprintf("Top %d ORA pathways across %d conditions",
                      length(top_pw), length(unique(df$cond))),
      xaxis = list(title = "", tickangle = -30),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 200, b = 100)
    )
}

#' Overlap visualization. For 2 conditions: 2-circle Venn. For 3 conditions:
#' 3-circle Venn. For 4+: a "presence" heatmap (rows = pathways, cols = cond).
.render_compare_overlap <- function(res) {
  sets <- lapply(res$sig_pathways, function(x) unique(x$union))
  sets <- sets[lengths(sets) > 0]
  k <- length(sets)
  if (k < 2) {
    return(plotly::plotly_empty() |>
           plotly::layout(title =
             "Need at least 2 conditions with significant pathways for overlap view"))
  }
  if (k == 2) return(.cmp_venn_2(sets))
  if (k == 3) return(.cmp_venn_3(sets))
  .cmp_heatmap_presence(sets)
}

.cmp_venn_2 <- function(sets) {
  a <- sets[[1]]; b <- sets[[2]]
  only_a <- length(setdiff(a, b))
  only_b <- length(setdiff(b, a))
  shared <- length(intersect(a, b))
  alpha <- 0.30
  shapes <- list(
    list(type = "circle", xref = "x", yref = "y",
         x0 = -1.4, x1 = 1.0, y0 = -1.2, y1 = 1.2,
         line = list(color = .PN_PAL$danger, width = 2),
         fillcolor = .pn_rgba(.PN_PAL$danger, alpha)),
    list(type = "circle", xref = "x", yref = "y",
         x0 = -1.0, x1 = 1.4, y0 = -1.2, y1 = 1.2,
         line = list(color = .PN_PAL$primary, width = 2),
         fillcolor = .pn_rgba(.PN_PAL$primary, alpha))
  )
  short <- function(s, n = 22) if (nchar(s) > n) paste0(substr(s, 1, n - 1), "…") else s
  fam <- .PN_PAL$font_family
  annotations <- list(
    list(x = -0.85, y = 0,    text = format(only_a, big.mark = ","),
         showarrow = FALSE, font = list(size = 20, color = .PN_PAL$danger_dark,  family = fam)),
    list(x =  0.85, y = 0,    text = format(only_b, big.mark = ","),
         showarrow = FALSE, font = list(size = 20, color = .PN_PAL$primary_dark, family = fam)),
    list(x =  0.00, y = 0,    text = format(shared, big.mark = ","),
         showarrow = FALSE, font = list(size = 20, color = .PN_PAL$text,         family = fam)),
    list(x = -1.05, y = 1.45, text = short(names(sets)[1]),
         showarrow = FALSE, font = list(size = 11, color = .PN_PAL$danger,       family = fam)),
    list(x =  1.05, y = 1.45, text = short(names(sets)[2]),
         showarrow = FALSE, font = list(size = 11, color = .PN_PAL$primary,      family = fam))
  )
  plotly::plot_ly() |>
    plotly::layout(
      shapes = shapes, annotations = annotations,
      title = "Significant pathways — overlap",
      xaxis = list(range = c(-2.3, 2.3), zeroline = FALSE, showgrid = FALSE,
                   showticklabels = FALSE, fixedrange = TRUE),
      yaxis = list(range = c(-2.0, 2.0), zeroline = FALSE, showgrid = FALSE,
                   showticklabels = FALSE, scaleanchor = "x", scaleratio = 1,
                   fixedrange = TRUE),
      margin = list(l = 10, r = 10, t = 50, b = 10),
      plot_bgcolor = "rgba(0,0,0,0)", paper_bgcolor = "rgba(0,0,0,0)",
      showlegend = FALSE
    ) |>
    plotly::config(displayModeBar = FALSE)
}

.cmp_venn_3 <- function(sets) {
  a <- sets[[1]]; b <- sets[[2]]; c_ <- sets[[3]]
  ab  <- intersect(a, b); ac <- intersect(a, c_); bc <- intersect(b, c_)
  abc <- intersect(ab, c_)
  cnt <- list(
    onlyA = length(setdiff(a, union(b, c_))),
    onlyB = length(setdiff(b, union(a, c_))),
    onlyC = length(setdiff(c_, union(a, b))),
    AB    = length(setdiff(ab, c_)),
    AC    = length(setdiff(ac, b)),
    BC    = length(setdiff(bc, a)),
    ABC   = length(abc)
  )
  alpha <- 0.26
  shapes <- list(
    list(type = "circle", x0 = -1.6, x1 =  0.4, y0 = -0.6, y1 = 1.4,
         line = list(color = .PN_PAL$danger, width = 2),
         fillcolor = .pn_rgba(.PN_PAL$danger, alpha)),
    list(type = "circle", x0 = -0.4, x1 =  1.6, y0 = -0.6, y1 = 1.4,
         line = list(color = .PN_PAL$primary, width = 2),
         fillcolor = .pn_rgba(.PN_PAL$primary, alpha)),
    list(type = "circle", x0 = -1.0, x1 =  1.0, y0 = -1.6, y1 = 0.4,
         line = list(color = .PN_PAL$success, width = 2),
         fillcolor = .pn_rgba(.PN_PAL$success, alpha))
  )
  short <- function(s, n = 16) if (nchar(s) > n) paste0(substr(s, 1, n - 1), "…") else s
  fam <- .PN_PAL$font_family
  ann <- function(x, y, txt, col = .PN_PAL$text, sz = 16) list(
    x = x, y = y, text = txt, showarrow = FALSE,
    font = list(size = sz, color = col, family = fam))
  annotations <- list(
    ann(-1.05,  0.85, format(cnt$onlyA, big.mark = ","), .PN_PAL$danger_dark),
    ann( 1.05,  0.85, format(cnt$onlyB, big.mark = ","), .PN_PAL$primary_dark),
    ann( 0.00, -1.10, format(cnt$onlyC, big.mark = ","), .PN_PAL$success_dark),
    ann( 0.00,  0.95, format(cnt$AB,    big.mark = ",")),
    ann(-0.65, -0.30, format(cnt$AC,    big.mark = ",")),
    ann( 0.65, -0.30, format(cnt$BC,    big.mark = ",")),
    ann( 0.00,  0.15, format(cnt$ABC,   big.mark = ","), .PN_PAL$text, 18),
    ann(-1.30,  1.65, short(names(sets)[1]), .PN_PAL$danger,  11),
    ann( 1.30,  1.65, short(names(sets)[2]), .PN_PAL$primary, 11),
    ann( 0.00, -1.85, short(names(sets)[3]), .PN_PAL$success, 11)
  )
  plotly::plot_ly() |>
    plotly::layout(
      shapes = shapes, annotations = annotations,
      title = "Significant pathways — 3-way overlap",
      xaxis = list(range = c(-2.6, 2.6), zeroline = FALSE, showgrid = FALSE,
                   showticklabels = FALSE, fixedrange = TRUE),
      yaxis = list(range = c(-2.4, 2.4), zeroline = FALSE, showgrid = FALSE,
                   showticklabels = FALSE, scaleanchor = "x", scaleratio = 1,
                   fixedrange = TRUE),
      margin = list(l = 10, r = 10, t = 50, b = 10),
      plot_bgcolor = "rgba(0,0,0,0)", paper_bgcolor = "rgba(0,0,0,0)",
      showlegend = FALSE
    ) |>
    plotly::config(displayModeBar = FALSE)
}

#' Presence/absence heatmap (≥4 conditions) — rows = pathways, cols = cond,
#' cell = 1 if pathway is significantly enriched in that condition.
.cmp_heatmap_presence <- function(sets) {
  all_pw <- sort(unique(unlist(sets)))
  if (length(all_pw) == 0)
    return(plotly::plotly_empty() |>
           plotly::layout(title = "No significant pathways across conditions"))
  m <- vapply(sets, function(s) as.integer(all_pw %in% s),
              integer(length(all_pw)))
  rownames(m) <- all_pw
  # Order rows by frequency (most-shared at top)
  m <- m[order(-rowSums(m)), , drop = FALSE]
  short <- function(s) ifelse(nchar(s) > 50, paste0(substr(s, 1, 47), "..."), s)
  plotly::plot_ly(z = m, x = colnames(m), y = short(rownames(m)),
                  type = "heatmap", colorscale = "Blues",
                  showscale = FALSE,
                  hovertemplate = "Pathway: %{y}<br>Condition: %{x}<br>Significant: %{z}<extra></extra>") |>
    plotly::layout(
      title = sprintf("Pathway presence across %d conditions", ncol(m)),
      xaxis = list(title = "", tickangle = -30),
      yaxis = list(title = "", automargin = TRUE),
      margin = list(l = 200, b = 100))
}

#' Long-form table: pathway × condition × {significant ✓/✗}
.build_overlap_table <- function(res) {
  sets <- lapply(res$sig_pathways, function(x) unique(x$union))
  all_pw <- sort(unique(unlist(sets)))
  if (length(all_pw) == 0)
    return(data.frame(message = "No significant pathways"))
  out <- data.frame(pathway = all_pw, stringsAsFactors = FALSE)
  for (cond in names(sets)) {
    out[[cond]] <- ifelse(all_pw %in% sets[[cond]], "✓", "")
  }
  out$n_conditions <- rowSums(out[, names(sets), drop = FALSE] == "✓")
  out[order(-out$n_conditions, out$pathway), , drop = FALSE]
}

# ---------------------------------------------------------------------------
# Gallery integration — auto-capture figures + register CSVs after each PN
# sub-tab run. All calls are guarded with exists() so the module degrades
# cleanly if the umap_gallery module isn't loaded. Wired into the run-
# observers above so the user never has to click a "save" button.
# ---------------------------------------------------------------------------

.pn_safe_capture <- function(values, plot_obj, stage, label, source_module) {
  if (!exists("capture_figure", mode = "function")) return(invisible(NULL))
  tryCatch(capture_figure(values, plot_obj, stage = stage,
                          label = label, source_module = source_module,
                          figure_type = "plot"),
           error = function(e) invisible(NULL))
}
.pn_safe_register <- function(values, df, description, source_module) {
  if (!exists("register_csv_export", mode = "function")) return(invisible(NULL))
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
    return(invisible(NULL))
  tryCatch(register_csv_export(values, data = df,
                               description = description,
                               source_module = source_module),
           error = function(e) invisible(NULL))
}

.pn_push_enrichment_to_gallery <- function(values, result, database,
                                           method, top_n) {
  src   <- "Pathway Enrichment"
  stage <- paste0("Enrichment: ", database, " · ", method)

  # Figures
  .pn_safe_capture(values, .pn_render_dotplot(result, top_n),
                   stage, paste0(database, " — Top pathways (dot plot)"), src)
  .pn_safe_capture(values, .pn_render_barplot(result, top_n),
                   stage, paste0(database, " — Top pathways (bar plot)"), src)
  if (!is.null(result$gsea) && nrow(result$gsea) > 0)
    .pn_safe_capture(values, .pn_render_gsea(result),
                     stage, paste0(database, " — GSEA NES"), src)

  # CSVs
  ora_df <- .pn_combined_ora_df(result)
  .pn_safe_register(values, ora_df,
                    paste0(database, " — ORA results"), src)
  .pn_safe_register(values, result$gsea,
                    paste0(database, " — GSEA results"), src)
}

.pn_push_compare_to_gallery <- function(values, result, database, method, top_n) {
  src   <- "Comparative Enrichment"
  cond_count <- length(result$audit$conditions %||% character(0))
  stage <- sprintf("Compare: %s · %s · %d conditions",
                   database, method, cond_count)

  # Figures
  .pn_safe_capture(values,
                   .render_compare_dotplot(result, top_n = top_n, method = method),
                   stage,
                   paste0(database, " — Comparative dotplot (", cond_count, " cond)"),
                   src)
  .pn_safe_capture(values,
                   .render_compare_overlap(result),
                   stage,
                   paste0(database, " — Pathway overlap (", cond_count, " cond)"),
                   src)

  # CSVs
  .pn_safe_register(values, result$combined_ora,
                    paste0(database, " — Comparative ORA"), src)
  .pn_safe_register(values, result$combined_gsea,
                    paste0(database, " — Comparative GSEA"), src)
  .pn_safe_register(values, .build_overlap_table(result),
                    paste0(database, " — Pathway-by-condition overlap"), src)
}

.pn_push_pwn_to_gallery <- function(values, result) {
  src   <- "Pathway Network"
  stage <- "Pathway-network graph"
  .pn_safe_capture(values, render_pathway_network(result),
                   stage, "Pathway network (Jaccard / Louvain communities)", src)
  .pn_safe_register(values, result$network_data,
                    "Pathway network — edges (CSV)", src)
  .pn_safe_register(values, result$clusters,
                    "Pathway network — communities (CSV)", src)
}

.pn_push_grn_to_gallery <- function(values, result) {
  src   <- "Gene Regulatory Networks"
  stage <- paste0("GRN: ", result$method %||% "?")
  .pn_safe_capture(values, render_grn_network(result),
                   stage, paste0("GRN (", result$method, ") — network"), src)
  if (!is.null(result$modules))
    .pn_safe_register(values, result$modules,
                      paste0("GRN (", result$method, ") — gene modules"), src)
  # Edges
  if (!is.null(result$network) && inherits(result$network, "igraph"))
    .pn_safe_register(values,
                      igraph::as_data_frame(result$network, what = "edges"),
                      paste0("GRN (", result$method, ") — edges"), src)
}

.pn_push_coexp_to_gallery <- function(values, result) {
  src   <- "Co-expression Network"
  stage <- "Co-expression"
  .pn_safe_capture(values, render_coexp_network(result),
                   stage, "Co-expression — network", src)
  .pn_safe_capture(values, render_correlation_matrix(result),
                   stage, "Co-expression — correlation matrix", src)
  .pn_safe_register(values, result$edges,
                    "Co-expression — edges (CSV)", src)
}
