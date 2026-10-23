# modules/de/de_server.R
# Sourced by app.R. Sources core_robustness.R, de_core.R, de_methods.R, de_ui.R.

local({
  this_dir <- tryCatch({
    f <- sys.frame(1)$ofile
    if (!is.null(f) && nzchar(f)) normalizePath(dirname(f), mustWork = TRUE)
    else file.path(getwd(), "modules", "de")
  }, error = function(e) file.path(getwd(), "modules", "de"))

  src <- function(rel, base = this_dir) {
    p <- normalizePath(file.path(base, rel), mustWork = FALSE)
    if (file.exists(p)) source(p, local = FALSE)
  }
  # core_robustness.R provides prepare_for_de + safe_get_assay_data
  src("../core_robustness.R", base = dirname(this_dir))
  src("de_core.R")
  src("de_methods.R")
  src("de_ui.R")
})

# Null-coalescing operator (used heavily below)
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# Plotly palette — mirrors app.css design tokens. Plotly takes raw hex/rgba,
# not CSS variables, so we name them here in one place.
# ---------------------------------------------------------------------------
.DE_PAL <- list(
  primary       = "#2563EB", primary_dark   = "#1E3A8A",
  danger        = "#DC2626", danger_dark    = "#7F1D1D",
  success       = "#16A34A",
  text          = "#0F172A",
  text_secondary= "#475569",   # --text-secondary
  font_family   = "Inter"
)
.de_rgba <- function(hex, alpha = 0.32) {
  rgb_vals <- grDevices::col2rgb(hex)
  sprintf("rgba(%d,%d,%d,%.2f)",
          rgb_vals[1, 1], rgb_vals[2, 1], rgb_vals[3, 1], alpha)
}

deServer <- function(id, values) {
  shiny::moduleServer(id, function(input, output, session) {

    de_results       <- shiny::reactiveVal(NULL)
    de_audit         <- shiny::reactiveVal(NULL)
    show_results     <- shiny::reactiveVal(FALSE)

    # ---- Populate dropdowns when active_dataset changes ----
    shiny::observe({
      shiny::req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[vapply(meta,
                    function(x) is.factor(x) || is.character(x), logical(1))]
      num_cols <- colnames(meta)[vapply(meta, is.numeric, logical(1))]

      shiny::updateSelectInput(session, "groupby",   choices = sort(cat_cols))
      shiny::updateSelectInput(session, "subset_col",
                               choices = c("(none)" = "", sort(cat_cols)))
      shiny::updateSelectInput(session, "sample_col",
                               choices = sort(cat_cols),
                               selected = detect_sample_col(values$active_dataset))
      shiny::updateSelectInput(session, "latent_vars", choices = num_cols)
    })

    shiny::observeEvent(input$groupby, {
      shiny::req(values$active_dataset, input$groupby)
      meta <- values$active_dataset@meta.data
      if (input$groupby %in% colnames(meta)) {
        groups <- sort(unique(as.character(meta[[input$groupby]])))
        groups <- groups[!is.na(groups)]
        shiny::updateSelectInput(session, "ident1", choices = groups)
        shiny::updateSelectInput(session, "ident2",
          choices = c("rest", "all", groups))
      }
    })

    shiny::observeEvent(input$subset_col, {
      shiny::req(values$active_dataset)
      col <- input$subset_col
      if (!is.null(col) && nzchar(col) &&
          col %in% colnames(values$active_dataset@meta.data)) {
        vals <- sort(unique(as.character(
          values$active_dataset@meta.data[[col]])))
        vals <- vals[!is.na(vals)]
        shiny::updateSelectInput(session, "subset_val", choices = vals)
      } else {
        shiny::updateSelectInput(session, "subset_val", choices = character(0))
      }
    })

    # ---- Run DE ----
    shiny::observeEvent(input$run, {
      shiny::req(values$active_dataset, input$groupby, input$ident1)

      method <- if (input$analysis_type == "pseudobulk") input$pb_method
                else input$sc_method

      shiny::withProgress(
        message = paste("Running", method, "..."), value = 0.1, {
        result <- run_de_analysis(
          obj = values$active_dataset,
          group_col = input$groupby,
          ident1 = input$ident1,
          ident2 = input$ident2,
          test_method = method,
          analysis_type = input$analysis_type,
          sample_col = if (input$analysis_type == "pseudobulk")
                         input$sample_col else NULL,
          subset_col = if (nzchar(input$subset_col)) input$subset_col else NULL,
          subset_val = if (!is.null(input$subset_val) &&
                           nzchar(input$subset_val)) input$subset_val else NULL,
          latent_vars = if (length(input$latent_vars) > 0)
                          input$latent_vars else NULL,
          seed = as.integer(input$seed %||% 42L),
          min_pct_display = input$min_pct_display %||% 0.10,
          logfc_threshold = input$logfc_threshold %||% 0,
          only_pos = isTRUE(input$only_pos)
        )
        shiny::incProgress(0.9)

        if (!is.null(result$data)) {
          de_results(result$data)
          de_audit(attr(result$data, "de_audit"))
          values$last_de_results <- result$data
          # --- Save to history (named list, FIFO eviction at .HISTORY_CAP) ---
          label_used <- .save_de_to_history(values,
            label = .auto_de_label(
              ident1 = input$ident1, ident2 = input$ident2,
              method = method,
              subset_col = if (nzchar(input$subset_col)) input$subset_col else NULL,
              subset_val = if (!is.null(input$subset_val) &&
                               nzchar(input$subset_val)) input$subset_val else NULL,
              custom = if (!is.null(input$run_label) &&
                           nzchar(trimws(input$run_label)))
                         trimws(input$run_label) else NULL),
            data = result$data,
            audit = attr(result$data, "de_audit"))
          show_results(TRUE)
          # --- Push figures + CSVs to gallery (auto, no clicks needed) ---
          .de_push_to_gallery(values, result$data,
                              audit = attr(result$data, "de_audit"),
                              input = input, label = label_used)
          n_sig <- sum(result$data$p_val_adj < (input$pval_cutoff %||% 0.05),
                       na.rm = TRUE)
          shiny::showNotification(
            sprintf("DE complete -> %d significant genes (adj p < %.3f). Saved to history + Figure Gallery.",
                    n_sig, input$pval_cutoff %||% 0.05),
            type = "message", duration = 5)
        } else {
          show_results(FALSE)
          shiny::showNotification(
            paste("DE failed:", result$error),
            type = "error", duration = 10)
        }
      })
    })

    # ---- Run DE in every cell type (loop) ----
    # Loops the comparison currently configured across every unique value of
    # the user's "Subset by cell-type column" — appends one DE result per
    # cell type to the history.
    shiny::observeEvent(input$run_per_cluster, {
      shiny::req(values$active_dataset, input$groupby, input$ident1)
      col <- input$subset_col
      if (is.null(col) || !nzchar(col)) {
        shiny::showNotification(
          "Pick a 'Subset by cell-type column' first — the loop iterates over its values.",
          type = "error", duration = 10); return()
      }
      vals <- sort(unique(as.character(
        values$active_dataset@meta.data[[col]])))
      vals <- vals[!is.na(vals)]
      if (length(vals) < 2) {
        shiny::showNotification(
          sprintf("Column '%s' has fewer than 2 values — nothing to loop over.",
                  col),
          type = "error", duration = 10); return()
      }
      # Safety cap
      if (length(vals) > 50) {
        shiny::showNotification(
          sprintf("Column '%s' has %d values — capping to 50.",
                  col, length(vals)),
          type = "warning", duration = 8)
        vals <- vals[seq_len(50)]
      }

      method <- if (input$analysis_type == "pseudobulk") input$pb_method
                else input$sc_method

      shiny::withProgress(
        message = sprintf("Looping DE across %d cell types...",
                          length(vals)),
        value = 0, {
        n_ok <- 0L; n_err <- 0L
        for (i in seq_along(vals)) {
          v <- vals[i]
          shiny::incProgress(1 / length(vals),
            detail = sprintf("%d / %d : %s", i, length(vals), v))
          result <- run_de_analysis(
            obj = values$active_dataset,
            group_col = input$groupby,
            ident1 = input$ident1, ident2 = input$ident2,
            test_method = method,
            analysis_type = input$analysis_type,
            sample_col = if (input$analysis_type == "pseudobulk")
                           input$sample_col else NULL,
            subset_col = col, subset_val = v,
            latent_vars = if (length(input$latent_vars) > 0)
                            input$latent_vars else NULL,
            seed = as.integer(input$seed %||% 42L),
            min_pct_display = input$min_pct_display %||% 0.10,
            logfc_threshold = input$logfc_threshold %||% 0,
            only_pos = isTRUE(input$only_pos))
          if (!is.null(result$data)) {
            label <- .auto_de_label(
              ident1 = input$ident1, ident2 = input$ident2,
              method = method, subset_col = col, subset_val = v)
            saved_label <- .save_de_to_history(values, label = label,
              data = result$data,
              audit = attr(result$data, "de_audit"))
            n_ok <- n_ok + 1L
            # Keep the last successful one visible in the result panel
            de_results(result$data)
            de_audit(attr(result$data, "de_audit"))
            values$last_de_results <- result$data
            show_results(TRUE)
            # Each iteration also pushes its plots + CSVs to the gallery
            .de_push_to_gallery(values, result$data,
                                audit = attr(result$data, "de_audit"),
                                input = input, label = saved_label)
          } else {
            n_err <- n_err + 1L
            message(sprintf("[de loop] %s -> %s", v,
                            result$error %||% "(no message)"))
          }
        }
        shiny::showNotification(
          sprintf("Per-cell-type loop done — %d saved, %d skipped (errors). Open Comparative Enrichment to compare.",
                  n_ok, n_err),
          type = if (n_ok > 0) "message" else "error", duration = 10)
      })
    })

    # ---- Dataset status strip ----
    output$dataset_status <- shiny::renderUI({
      obj <- values$active_dataset
      if (is.null(obj)) {
        return(shiny::div(class = "de-dataset-strip",
          shiny::div(class = "dot muted"),
          shiny::span("No dataset loaded yet — pick one in Data Upload.")))
      }
      sp <- if (!is.null(values$species_global) &&
                nzchar(values$species_global))
              values$species_global else "—"
      shiny::div(class = "de-dataset-strip",
        shiny::div(class = "dot"),
        shiny::span(class = "meta", format(ncol(obj), big.mark = ",")),
        shiny::span("cells"),
        shiny::span("·"),
        shiny::span(class = "meta", format(nrow(obj), big.mark = ",")),
        shiny::span("genes"),
        shiny::span("·"),
        shiny::span(class = "meta", sp))
    })

    # ---- Results UI scaffold ----
    output$results_ui <- shiny::renderUI({
      ns <- session$ns
      if (!isTRUE(show_results())) {
        return(shiny::div(class = "de-empty",
          shiny::div(class = "icon", shiny::icon("chart-area")),
          shiny::h4("Results will appear here"),
          shiny::p("Configure the comparison on the left, pick a method, ",
                   "then click ", shiny::strong("Run DE Analysis"), ". ",
                   "You'll see a volcano plot, MA plot, sortable results table, ",
                   "expression heatmap, gene plots, and the audit log.")))
      }
      shiny::div(class = "de-results-card",
      shiny::tabsetPanel(
        shiny::tabPanel("Volcano",
          shiny::br(),
          plotly::plotlyOutput(ns("volcano"), height = "560px"),
          shiny::uiOutput(ns("volcano_summary"))),
        shiny::tabPanel("MA Plot",
          shiny::br(), plotly::plotlyOutput(ns("ma"), height = "560px")),
        shiny::tabPanel("Results Table",
          shiny::br(), DT::DTOutput(ns("table")),
          shiny::br(),
          shiny::downloadButton(ns("download_csv"),
                                "Download all (CSV)", class = "btn-success")),
        shiny::tabPanel("Heatmap",
          shiny::br(), plotly::plotlyOutput(ns("heatmap"), height = "640px")),
        shiny::tabPanel("Gene plots",
          shiny::br(),
          shiny::textAreaInput(ns("genes_to_plot"),
                               "Genes (comma-separated)",
                               rows = 2,
                               placeholder = "e.g. Gene001, Gene002"),
          shiny::actionButton(ns("plot_genes"), "Plot",
                              class = "btn-primary"),
          shiny::br(), shiny::br(),
          plotly::plotlyOutput(ns("gene_plots"), height = "560px")),
        shiny::tabPanel("Audit Log",
          shiny::br(), DT::DTOutput(ns("audit_table")))
      ))
    })

    # ---- Renderers ----
    output$volcano <- plotly::renderPlotly({
      shiny::req(de_results())
      .render_volcano(de_results(), input)
    })
    output$ma <- plotly::renderPlotly({
      shiny::req(de_results())
      .render_ma(de_results(), input)
    })
    output$table <- DT::renderDT({
      shiny::req(de_results())
      df <- de_results()
      keep_cols <- c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2",
                     "AUC","cohens_d","method")
      keep_cols <- intersect(keep_cols, colnames(df))
      DT::datatable(
        df[, keep_cols, drop = FALSE],
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE) |>
        DT::formatSignif(c("p_val","p_val_adj"), 3) |>
        DT::formatRound(intersect(c("avg_log2FC","pct.1","pct.2","AUC","cohens_d"),
                                  keep_cols), 3)
    })
    output$heatmap <- plotly::renderPlotly({
      shiny::req(de_results(), values$active_dataset)
      .render_heatmap(de_results(), values$active_dataset, input)
    })

    shiny::observeEvent(input$plot_genes, {
      output$gene_plots <- plotly::renderPlotly({
        shiny::req(values$active_dataset, input$genes_to_plot)
        genes <- trimws(strsplit(input$genes_to_plot, ",")[[1]])
        genes <- genes[nzchar(genes)]
        genes <- intersect(genes, rownames(values$active_dataset))
        shiny::validate(shiny::need(length(genes) > 0,
                                    "No valid genes found in dataset."))
        .render_gene_box(values$active_dataset, genes, input$groupby)
      })
    })

    output$audit_table <- DT::renderDT({
      shiny::req(de_audit())
      a <- de_audit()
      df <- data.frame(
        Field = names(a),
        Value = vapply(a, function(v) {
          if (length(v) == 0) "(none)"
          else if (length(v) > 1) paste(v, collapse = ", ")
          else as.character(v)
        }, character(1)),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, options = list(pageLength = 20, dom = "t"),
                    rownames = FALSE)
    })

    output$download_csv <- shiny::downloadHandler(
      filename = function()
        sprintf("de_%s_%s.csv",
                attr(de_results(), "de_audit")$test %||% "result",
                format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        de_results(), file, row.names = FALSE)
    )

    # ---- Volcano summary panel (cards + venn + stats) ----
    output$volcano_summary <- shiny::renderUI({
      shiny::req(de_results(), de_audit())
      .render_volcano_summary(de_results(), de_audit(), input, session$ns)
    })
    output$de_venn <- plotly::renderPlotly({
      shiny::req(de_results())
      .render_de_venn(de_results(), input)
    })
  })
}

# ---------------------------------------------------------------------------
# DE result history — persistent list of all DE runs in this session so that
# the Comparative Enrichment sub-tab can pick any 2+ to enrich side-by-side.
# Capped at .HISTORY_CAP entries (FIFO) to prevent memory blowup.
# ---------------------------------------------------------------------------
.HISTORY_CAP <- 50L

#' Build a human-readable label from the DE inputs.
.auto_de_label <- function(ident1, ident2 = NULL, method, subset_col = NULL,
                           subset_val = NULL, custom = NULL) {
  if (!is.null(custom) && nzchar(custom)) return(custom)
  ident2_lbl <- if (is.null(ident2) || ident2 %in% c("rest", "", "all"))
                  "rest" else as.character(ident2)
  base <- sprintf("%s vs %s · %s",
                  as.character(ident1), ident2_lbl, method)
  if (!is.null(subset_col) && nzchar(subset_col) &&
      !is.null(subset_val) && nzchar(subset_val)) {
    base <- paste0(base, " · ", subset_col, "=", subset_val)
  }
  base
}

#' Append a DE run to the history. Collisions resolved by suffixing "(2)" etc.
#' Caps total entries at .HISTORY_CAP (FIFO drop).
.save_de_to_history <- function(values, label, data, audit) {
  hist <- values$de_history %||% list()
  # Resolve label collisions
  final_label <- label
  i <- 2L
  while (final_label %in% names(hist)) {
    final_label <- sprintf("%s (%d)", label, i); i <- i + 1L
  }
  hist[[final_label]] <- list(
    label     = final_label,
    data      = data,
    audit     = audit,
    timestamp = Sys.time()
  )
  # FIFO eviction
  if (length(hist) > .HISTORY_CAP) {
    drop_n <- length(hist) - .HISTORY_CAP
    hist <- hist[-seq_len(drop_n)]
  }
  values$de_history <- hist
  invisible(final_label)
}

# ---- Plot helpers (file-local) ----------------------------------------------
.render_volcano <- function(df, input) {
  pcut <- input$pval_cutoff %||% 0.05
  fcut <- input$logfc_threshold %||% 0.10
  df <- df[!is.na(df$avg_log2FC) & !is.na(df$p_val_adj), ]
  df$neglog10p <- -log10(pmax(df$p_val_adj, 1e-300))
  df$status <- ifelse(df$p_val_adj < pcut & df$avg_log2FC >  fcut, "Up",
               ifelse(df$p_val_adj < pcut & df$avg_log2FC < -fcut, "Down",
                                                                   "NS"))
  pal <- c(Up = "#E64B35", Down = "#3C5488", NS = "#9aa3af")
  top <- head(df[order(df$p_val_adj), "gene"],
              input$label_top_n %||% 20)
  df$label <- ifelse(df$gene %in% top, df$gene, "")
  plotly::plot_ly(df, x = ~avg_log2FC, y = ~neglog10p,
                  type = "scatter", mode = "markers+text",
                  color = ~status, colors = pal,
                  text = ~label, textposition = "top center",
                  hoverinfo = "text",
                  hovertext = ~paste0(gene, "<br>log2FC=", round(avg_log2FC, 2),
                                      "<br>adj.p=", signif(p_val_adj, 3))) |>
    plotly::layout(
      title = "Volcano",
      xaxis = list(title = "avg log2 fold change"),
      yaxis = list(title = "-log10(adj p)"))
}

.render_ma <- function(df, input) {
  df <- df[!is.na(df$avg_log2FC) & !is.na(df$pct.1) & !is.na(df$pct.2), ]
  if (nrow(df) == 0) {
    return(plotly::plotly_empty() |>
           plotly::layout(title = "MA plot requires pct.1/pct.2 (single-cell methods)"))
  }
  df$mean_expr <- (df$pct.1 + df$pct.2) / 2
  pcut <- input$pval_cutoff %||% 0.05
  df$status <- ifelse(df$p_val_adj < pcut, "sig", "ns")
  plotly::plot_ly(df, x = ~mean_expr, y = ~avg_log2FC,
                  type = "scatter", mode = "markers",
                  color = ~status,
                  colors = c(sig = "#E64B35", ns = "#9aa3af"),
                  text = ~gene,
                  hoverinfo = "text") |>
    plotly::layout(
      title = "MA plot",
      xaxis = list(title = "mean detection rate"),
      yaxis = list(title = "avg log2 fold change"))
}

.render_heatmap <- function(df, obj, input) {
  top_n <- as.integer(input$label_top_n %||% 20)
  pcut  <- input$pval_cutoff %||% 0.05
  top <- head(df[df$p_val_adj < pcut, "gene"], top_n)
  if (length(top) < 2) {
    return(plotly::plotly_empty() |>
           plotly::layout(title = "Not enough significant genes for heatmap"))
  }
  expr <- safe_get_assay_data(obj, "data")[top, , drop = FALSE]
  group <- as.character(obj@meta.data[[input$groupby]])
  ord <- order(group)
  expr <- expr[, ord, drop = FALSE]; group <- group[ord]
  plotly::plot_ly(z = as.matrix(expr), x = colnames(expr), y = rownames(expr),
                  type = "heatmap", colors = viridis::viridis(50)) |>
    plotly::layout(title = sprintf("Top %d significant genes", length(top)))
}

.render_gene_box <- function(obj, genes, group_col) {
  expr <- safe_get_assay_data(obj, "data")[genes, , drop = FALSE]
  long <- do.call(rbind, lapply(genes, function(g) {
    data.frame(gene = g,
               expr = as.numeric(expr[g, ]),
               group = as.character(obj@meta.data[[group_col]]))
  }))
  plotly::plot_ly(long, x = ~group, y = ~expr,
                  color = ~gene, type = "box") |>
    plotly::layout(boxmode = "group",
                   yaxis = list(title = "log-normalised expression"),
                   xaxis = list(title = group_col))
}

# ---------------------------------------------------------------------------
# Volcano-summary panel: metric cards + conditions card + Venn diagram +
# detailed stats table. Rendered under the volcano plot in the same tab.
# All values are derived from the DE result + audit log — no hard-coding.
# ---------------------------------------------------------------------------

#' Build the summary UI block. Returns a single shiny::div().
.render_volcano_summary <- function(df, audit, input, ns) {
  pcut <- input$pval_cutoff      %||% 0.05
  fcut <- input$logfc_threshold  %||% 0.10

  # Significance partitioning
  sig <- !is.na(df$p_val_adj) & df$p_val_adj < pcut &
         !is.na(df$avg_log2FC) & abs(df$avg_log2FC) > fcut
  up   <- sig & df$avg_log2FC >  fcut
  down <- sig & df$avg_log2FC < -fcut
  n_up   <- sum(up,   na.rm = TRUE)
  n_down <- sum(down, na.rm = TRUE)
  n_sig  <- n_up + n_down

  # Top hits by p-value and by magnitude
  top_p   <- if (any(!is.na(df$p_val_adj)))
    df$gene[order(df$p_val_adj, na.last = NA)][1] else NA_character_
  top_fc  <- if (any(!is.na(df$avg_log2FC)))
    df$gene[order(-abs(df$avg_log2FC), na.last = NA)][1] else NA_character_
  min_padj <- if (any(sig)) min(df$p_val_adj[sig], na.rm = TRUE) else NA_real_
  med_fc   <- if (any(sig)) stats::median(abs(df$avg_log2FC[sig]),
                                          na.rm = TRUE) else NA_real_
  n_gt1 <- sum(!is.na(df$avg_log2FC) & abs(df$avg_log2FC) > 1, na.rm = TRUE)
  n_gt2 <- sum(!is.na(df$avg_log2FC) & abs(df$avg_log2FC) > 2, na.rm = TRUE)

  ident1 <- audit$ident1   %||% "Group 1"
  ident2 <- audit$ident2   %||% "rest"
  if (is.na(ident2) || !nzchar(as.character(ident2))) ident2 <- "rest"

  fmt <- function(x) format(x, big.mark = ",", scientific = FALSE)
  sci <- function(x) if (is.na(x) || !is.finite(x)) "—"
                     else signif(x, 3)

  # --- 1. Metric cards row ---
  cards <- shiny::div(class = "de-metrics",
    .de_metric_card("Cells analysed", fmt(audit$n_cells_total %||% NA),
                    sprintf("%s vs %s",
                            fmt(audit$n_group1 %||% NA),
                            fmt(audit$n_group2 %||% NA)),
                    accent = "neutral", icon = "users"),
    .de_metric_card("Genes tested", fmt(audit$n_genes_tested %||% NA),
                    sprintf("of %s in dataset",
                            fmt(audit$n_genes_tested %||% NA)),
                    accent = "success", icon = "dna"),
    .de_metric_card("Up-regulated",  fmt(n_up),
                    sprintf("adj p < %g  &  log2FC > %g", pcut, fcut),
                    accent = "up",   icon = "arrow-up"),
    .de_metric_card("Down-regulated", fmt(n_down),
                    sprintf("adj p < %g  &  log2FC < -%g", pcut, fcut),
                    accent = "down", icon = "arrow-down")
  )

  # --- 2. Conditions + Venn row ---
  cond_card <- shiny::div(class = "de-card",
    shiny::h5(class = "de-card-title", "Conditions compared"),
    shiny::div(class = "de-cond-block",
      shiny::div(class = "de-cond-row",
        shiny::div(class = "badge-grp g1", "1"),
        shiny::div(class = "meta",
          shiny::span(class = "name", ident1),
          shiny::span(class = "sub",
            sprintf("%s cells", fmt(audit$n_group1 %||% NA))))
      ),
      shiny::div(class = "de-cond-vs", "vs"),
      shiny::div(class = "de-cond-row",
        shiny::div(class = "badge-grp g2", "2"),
        shiny::div(class = "meta",
          shiny::span(class = "name", ident2),
          shiny::span(class = "sub",
            sprintf("%s cells", fmt(audit$n_group2 %||% NA))))
      )
    ),
    shiny::tags$dl(class = "de-cond-meta",
      shiny::tags$dt("Group column"),  shiny::tags$dd(audit$group_col),
      shiny::tags$dt("Method"),        shiny::tags$dd(audit$test),
      shiny::tags$dt("Latent vars"),
        shiny::tags$dd(if (length(audit$latent_vars) > 0)
                         paste(audit$latent_vars, collapse = ", ") else "(none)"),
      shiny::tags$dt("Seed"),          shiny::tags$dd(as.character(audit$seed)),
      shiny::tags$dt("Significant"),
        shiny::tags$dd(sprintf("%s / %s (%.1f%%)",
                               fmt(n_sig),
                               fmt(audit$n_genes_tested %||% nrow(df)),
                               100 * n_sig /
                                 max(1, audit$n_genes_tested %||% nrow(df))))
    )
  )

  venn_card <- shiny::div(class = "de-card",
    shiny::h5(class = "de-card-title",
              sprintf("Gene expression overlap (≥ %g%% cells)",
                      100 * (input$min_pct_display %||% 0.10))),
    plotly::plotlyOutput(ns("de_venn"), height = "300px"),
    shiny::div(style = "font-size: 11px; color: var(--text-secondary); margin-top: 6px;",
      shiny::icon("info-circle"),
      " Genes expressed in each group above the display threshold. The overlap shows ",
      "commonly-expressed genes (often housekeeping). Unique regions are candidate markers.")
  )

  split_row <- shiny::div(class = "de-split-row", cond_card, venn_card)

  # --- 3. Detailed stats table ---
  rows <- list(
    list("Top hit (smallest adj p)",
         if (!is.na(top_p)) sprintf("%s  (adj p = %s)",
                                    top_p, sci(min_padj)) else "—"),
    list("Top hit (largest |log2FC|)",
         if (!is.na(top_fc))
           sprintf("%s  (log2FC = %s)", top_fc,
                   sci(df$avg_log2FC[df$gene == top_fc][1])) else "—"),
    list("Median |log2FC| (significant)",          sci(med_fc)),
    list("Smallest adj p value",                    sci(min_padj)),
    list("Genes with |log2FC| > 1",                 fmt(n_gt1)),
    list("Genes with |log2FC| > 2",                 fmt(n_gt2)),
    list("Significance thresholds",
         sprintf("adj p < %g  &  |log2FC| > %g", pcut, fcut)),
    list("Total genes tested",                      fmt(audit$n_genes_tested %||% nrow(df))),
    list("Reproducibility (seed)",                  as.character(audit$seed)),
    list("Seurat version",                          audit$seurat_version %||% "—"),
    list("R version",                               audit$r_version %||% "—"),
    list("Analysis timestamp",
         if (!is.null(audit$timestamp))
           format(audit$timestamp, "%Y-%m-%d %H:%M:%S %Z") else "—")
  )
  stat_card <- shiny::div(class = "de-card",
    shiny::h5(class = "de-card-title", "Detailed statistics"),
    shiny::tags$table(class = "de-stat-table",
      shiny::tags$tbody(
        lapply(rows, function(r)
          shiny::tags$tr(
            shiny::tags$td(r[[1]]),
            shiny::tags$td(r[[2]])))
      )
    )
  )

  shiny::div(class = "de-summary-section",
    shiny::div(class = "de-section-h",
      shiny::span(class = "icon", shiny::icon("chart-pie")),
      shiny::span("Analysis Summary")),
    cards,
    split_row,
    stat_card
  )
}

#' Single metric card.
.de_metric_card <- function(label, value, sub, accent = "neutral", icon = NULL) {
  shiny::div(class = paste("de-metric-card", accent),
    shiny::div(class = "accent"),
    shiny::span(class = "label",
      if (!is.null(icon)) list(shiny::icon(icon), " ", label) else label),
    shiny::span(class = "value", value),
    shiny::span(class = "sub", sub)
  )
}

#' 2-circle Venn of genes expressed in ident1 vs ident2 (via pct.1 / pct.2).
.render_de_venn <- function(df, input) {
  thr <- input$min_pct_display %||% 0.10
  d <- df[!is.na(df$pct.1) & !is.na(df$pct.2), , drop = FALSE]
  if (nrow(d) == 0) {
    return(plotly::plotly_empty() |>
           plotly::layout(title = list(text =
             "Venn requires pct.1 / pct.2 (single-cell methods only)",
             font = list(size = 12))))
  }
  in_g1   <- d$pct.1 >= thr
  in_g2   <- d$pct.2 >= thr
  only1   <- sum(in_g1 & !in_g2, na.rm = TRUE)
  only2   <- sum(in_g2 & !in_g1, na.rm = TRUE)
  shared  <- sum(in_g1 &  in_g2, na.rm = TRUE)
  total_e <- only1 + only2 + shared

  # Use design-token palette (matches app.css --danger / --primary)
  alpha <- 0.32
  shapes <- list(
    list(type = "circle", xref = "x", yref = "y",
         x0 = -1.4, x1 =  1.0, y0 = -1.2, y1 = 1.2,
         line = list(color = .DE_PAL$danger, width = 2),
         fillcolor = .de_rgba(.DE_PAL$danger, alpha)),
    list(type = "circle", xref = "x", yref = "y",
         x0 = -1.0, x1 =  1.4, y0 = -1.2, y1 = 1.2,
         line = list(color = .DE_PAL$primary, width = 2),
         fillcolor = .de_rgba(.DE_PAL$primary, alpha))
  )

  # Pull ident names from the de_audit attribute attached to df
  audit <- attr(df, "de_audit") %||% list()
  ident1 <- audit$ident1 %||% "Group 1"
  ident2 <- audit$ident2
  if (is.null(ident2) || is.na(ident2) || !nzchar(as.character(ident2)))
    ident2 <- "rest"
  # Truncate long labels so they fit above the circles
  short <- function(s, n = 18) if (nchar(s) > n) paste0(substr(s, 1, n - 1), "…") else s

  # Counts: positioned in each region
  fam <- .DE_PAL$font_family
  annotations <- list(
    list(x = -0.85, y = 0, text = format(only1, big.mark = ","),
         showarrow = FALSE, font = list(size = 18, color = .DE_PAL$danger_dark,  family = fam)),
    list(x =  0.85, y = 0, text = format(only2, big.mark = ","),
         showarrow = FALSE, font = list(size = 18, color = .DE_PAL$primary_dark, family = fam)),
    list(x =  0.00, y = 0, text = format(shared, big.mark = ","),
         showarrow = FALSE, font = list(size = 18, color = .DE_PAL$text,         family = fam)),
    # Group labels above circles
    list(x = -1.05, y = 1.45, text = short(ident1),
         showarrow = FALSE, font = list(size = 11, color = .DE_PAL$danger,  family = fam)),
    list(x =  1.05, y = 1.45, text = short(as.character(ident2)),
         showarrow = FALSE, font = list(size = 11, color = .DE_PAL$primary, family = fam)),
    list(x =  0.00, y = -1.55, text = sprintf("Total expressed: %s",
                                              format(total_e, big.mark = ",")),
         showarrow = FALSE,
         font = list(size = 11, color = .DE_PAL$text_secondary,
                     family = .DE_PAL$font_family))
  )

  plotly::plot_ly() |>
    plotly::layout(
      shapes = shapes, annotations = annotations,
      xaxis = list(range = c(-2.3, 2.3), zeroline = FALSE,
                   showgrid = FALSE, showticklabels = FALSE,
                   fixedrange = TRUE),
      yaxis = list(range = c(-2.0, 2.0), zeroline = FALSE,
                   showgrid = FALSE, showticklabels = FALSE,
                   scaleanchor = "x", scaleratio = 1, fixedrange = TRUE),
      margin = list(l = 10, r = 10, t = 10, b = 10),
      plot_bgcolor = "rgba(0,0,0,0)",
      paper_bgcolor = "rgba(0,0,0,0)",
      showlegend = FALSE
    ) |>
    plotly::config(displayModeBar = FALSE)
}

# ---------------------------------------------------------------------------
# Gallery integration — auto-capture figures + register CSVs after a DE run.
# All calls are guarded with exists() so the module degrades cleanly if the
# umap_gallery module isn't loaded.
# ---------------------------------------------------------------------------
.de_push_to_gallery <- function(values, df, audit, input, label) {
  src <- "Differential Expression"
  stage <- paste0("DE: ", label)

  # --- Figures ---
  if (exists("capture_figure", mode = "function")) {
    tryCatch(capture_figure(values, .render_volcano(df, input),
      stage = stage, label = paste0(label, " — Volcano"),
      source_module = src, figure_type = "plot"),
      error = function(e) NULL)
    # MA needs pct.1/pct.2 — guard
    if (all(c("pct.1","pct.2") %in% colnames(df))) {
      tryCatch(capture_figure(values, .render_ma(df, input),
        stage = stage, label = paste0(label, " — MA plot"),
        source_module = src, figure_type = "plot"),
        error = function(e) NULL)
    }
    # Heatmap requires the Seurat object
    if (!is.null(values$active_dataset)) {
      tryCatch(capture_figure(values,
        .render_heatmap(df, values$active_dataset, input),
        stage = stage, label = paste0(label, " — Top-gene heatmap"),
        source_module = src, figure_type = "plot"),
        error = function(e) NULL)
    }
    # Venn (only meaningful for single-cell methods with pct.* columns)
    if (all(c("pct.1","pct.2") %in% colnames(df))) {
      tryCatch(capture_figure(values, .render_de_venn(df, input),
        stage = stage,
        label = paste0(label, " — Expression overlap (Venn)"),
        source_module = src, figure_type = "plot"),
        error = function(e) NULL)
    }
  }

  # --- CSVs ---
  if (exists("register_csv_export", mode = "function")) {
    tryCatch(register_csv_export(values, data = df,
      description = paste0(label, " — DE results"),
      source_module = src),
      error = function(e) NULL)
    # Audit as a 2-column table
    if (!is.null(audit)) {
      audit_df <- data.frame(
        Field = names(audit),
        Value = vapply(audit, function(v) {
          if (length(v) == 0) "(none)"
          else if (length(v) > 1) paste(v, collapse = ", ")
          else as.character(v)
        }, character(1)), stringsAsFactors = FALSE)
      tryCatch(register_csv_export(values, data = audit_df,
        description = paste0(label, " — DE audit log"),
        source_module = src),
        error = function(e) NULL)
    }
  }
  invisible(NULL)
}
