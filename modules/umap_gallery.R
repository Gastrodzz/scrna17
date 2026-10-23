# =============================================================================
# FIGURE & EXPORT GALLERY MODULE
# =============================================================================
# Combines UMAP snapshot gallery + universal figure gallery + CSV export library
# into a single browsable interface.

# -----------------------------------------------------------------------------
# Helper: capture a UMAP snapshot into values$umap_gallery
# -----------------------------------------------------------------------------
capture_umap_snapshot <- function(values, plot_obj, stage, coloring = "", method = "", label = "") {
  tryCatch({
    if (is.null(plot_obj)) return(invisible(NULL))

    n_cells <- tryCatch(ncol(isolate(values$active_dataset)), error = function(e) NA)
    entry <- list(
      plot = plot_obj,
      stage = stage,
      coloring = coloring,
      method = method,
      label = if (nzchar(label)) label else stage,
      timestamp = Sys.time(),
      n_cells = n_cells,
      id = paste0(stage, "_", as.numeric(Sys.time()))
    )

    gallery <- isolate(values$umap_gallery)
    if (is.null(gallery)) gallery <- list()

    # Dedup: skip if same stage added within last 5 seconds
    now <- Sys.time()
    for (g in gallery) {
      if (g$stage == stage && abs(difftime(now, g$timestamp, units = "secs")) < 5) {
        return(invisible(NULL))
      }
    }

    if (length(gallery) >= 50) gallery <- gallery[-1]
    gallery[[length(gallery) + 1]] <- entry
    values$umap_gallery <- gallery
    invisible(NULL)
  }, error = function(e) invisible(NULL))
}

# -----------------------------------------------------------------------------
# Helper: capture any figure into values$figure_gallery
# -----------------------------------------------------------------------------
# Works with plotly, ggplot, or grob objects. Auto-converts ggplot to plotly.
capture_figure <- function(values, plot_obj, stage, figure_type = "plot",
                           label = "", source_module = "") {
  tryCatch({
    if (is.null(plot_obj)) return(invisible(NULL))

    # Auto-convert ggplot to plotly
    if (inherits(plot_obj, "gg") || inherits(plot_obj, "ggplot")) {
      plot_obj <- tryCatch(plotly::ggplotly(plot_obj), error = function(e) plot_obj)
    }

    entry <- list(
      plot = plot_obj,
      stage = stage,
      figure_type = figure_type,
      label = if (nzchar(label)) label else stage,
      source_module = source_module,
      timestamp = Sys.time(),
      id = paste0("fig_", as.numeric(Sys.time()))
    )

    gallery <- isolate(values$figure_gallery)
    if (is.null(gallery)) gallery <- list()

    # Dedup: skip if same stage+label added within last 5 seconds
    now <- Sys.time()
    for (g in gallery) {
      if (g$stage == stage && g$label == entry$label &&
          abs(difftime(now, g$timestamp, units = "secs")) < 5) {
        return(invisible(NULL))
      }
    }

    if (length(gallery) >= 100) gallery <- gallery[-1]
    gallery[[length(gallery) + 1]] <- entry
    values$figure_gallery <- gallery
    invisible(NULL)
  }, error = function(e) invisible(NULL))
}

# -----------------------------------------------------------------------------
# Helper: register a CSV export into values$csv_library
# -----------------------------------------------------------------------------
register_csv_export <- function(values, data = NULL, filename = NULL, description = "",
                                source_module = "", n_rows = NA, n_cols = NA) {
  tryCatch({
    if (is.null(data) && is.null(filename)) return(invisible(NULL))

    # Extract data.frame from DT::datatable objects
    if (!is.null(data) && inherits(data, "datatables")) {
      data <- tryCatch(data$x$data, error = function(e) NULL)
      if (is.null(data)) return(invisible(NULL))
    }

    # Auto-detect dimensions from data
    if (!is.null(data) && is.data.frame(data)) {
      if (is.na(n_rows)) n_rows <- nrow(data)
      if (is.na(n_cols)) n_cols <- ncol(data)
      # Write to tempdir so gallery can offer downloads
      if (is.null(filename)) {
        fname <- paste0(gsub("[^a-zA-Z0-9_]", "_", description), "_",
                        format(Sys.time(), "%H%M%S"), ".csv")
        filename <- file.path(tempdir(), fname)
      }
      write.csv(data, filename, row.names = FALSE)
    }

    # If we still don't have a filename, skip registration
    if (is.null(filename)) return(invisible(NULL))

    entry <- list(
      filename = filename,
      description = description,
      source_module = source_module,
      n_rows = n_rows,
      n_cols = n_cols,
      timestamp = Sys.time(),
      id = paste0("csv_", as.numeric(Sys.time()))
    )

    lib <- isolate(values$csv_library)
    if (is.null(lib)) lib <- list()

    # Dedup: if same description+module already exists, update it
    for (i in seq_along(lib)) {
      if (lib[[i]]$description == description &&
          lib[[i]]$source_module == source_module) {
        lib[[i]] <- entry
        values$csv_library <- lib
        return(invisible(NULL))
      }
    }

    if (length(lib) >= 100) lib <- lib[-1]
    lib[[length(lib) + 1]] <- entry
    values$csv_library <- lib
    invisible(NULL)
  }, error = function(e) invisible(NULL))
}

# -----------------------------------------------------------------------------
# Helper: wrap plotlyOutput with "Add to Gallery" button
# -----------------------------------------------------------------------------
# Usage in UI: plotlyWithSave(ns, "my_plot", height = "500px", label = "QC Violins")
plotlyWithSave <- function(ns, outputId, height = "400px", label = outputId) {
  div(style = "position: relative;",
    plotlyOutput(ns(outputId), height = height),
    actionLink(ns(paste0("save_fig_", outputId)), NULL,
      icon = icon("plus-circle"),
      style = "position:absolute;top:4px;right:8px;font-size:16px;color:#4caf50;z-index:10;opacity:0.7;",
      title = paste0("Add '", label, "' to Figure Gallery"),
      class = "save-to-gallery-btn")
  )
}

# Helper: wrap DT output with "Add to Library" button
dtWithSave <- function(ns, outputId, label = outputId) {
  div(style = "position: relative;",
    DT::dataTableOutput(ns(outputId)),
    actionLink(ns(paste0("save_csv_", outputId)), NULL,
      icon = icon("plus-circle"),
      style = "position:absolute;top:4px;right:8px;font-size:16px;color:#2196f3;z-index:10;opacity:0.7;",
      title = paste0("Add '", label, "' to CSV Library"),
      class = "save-to-library-btn")
  )
}

# Server helper: set up save-to-gallery observers for a module's plots
# plot_registry: named list, e.g. list(my_plot = list(getter = function() { ... }, label = "My Plot", module = "QC"))
# Each getter should return the plotly object (or NULL)
setup_gallery_buttons <- function(input, output, session, values, plot_registry) {
  lapply(names(plot_registry), function(plot_id) {
    btn_id <- paste0("save_fig_", plot_id)
    info <- plot_registry[[plot_id]]
    observeEvent(input[[btn_id]], {
      tryCatch({
        p <- info$getter()
        if (!is.null(p)) {
          capture_figure(values, p, info$module,
                         figure_type = info$label,
                         label = info$label,
                         source_module = info$module)
          showNotification(paste0("'", info$label, "' added to Figure Gallery"),
                           type = "message", duration = 2)
        } else {
          showNotification("No figure to save — run the analysis first", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Could not save figure:", e$message), type = "error")
      })
    }, ignoreInit = TRUE)
  })
}

# Server helper: set up save-to-library observers for a module's tables
# table_registry: named list, e.g. list(my_table = list(getter = function() { ... }, label = "My Table", module = "QC"))
# Each getter should return a data.frame (or NULL)
setup_csv_buttons <- function(input, output, session, values, table_registry) {
  lapply(names(table_registry), function(tbl_id) {
    btn_id <- paste0("save_csv_", tbl_id)
    info <- table_registry[[tbl_id]]
    observeEvent(input[[btn_id]], {
      tryCatch({
        df <- info$getter()
        if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
          fname <- paste0(gsub("[^a-zA-Z0-9_]", "_", info$label), "_",
                          format(Sys.time(), "%H%M%S"), ".csv")
          fpath <- file.path(tempdir(), fname)
          write.csv(df, fpath, row.names = FALSE)
          register_csv_export(values, filename = fpath,
                              description = info$label,
                              source_module = info$module,
                              n_rows = nrow(df), n_cols = ncol(df))
          showNotification(paste0("'", info$label, "' added to CSV Library (",
                                  nrow(df), " rows)"),
                           type = "message", duration = 2)
        } else {
          showNotification("No data to save — run the analysis first", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Could not save table:", e$message), type = "error")
      })
    }, ignoreInit = TRUE)
  })
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------

umapGalleryUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        div(
          style = "display: flex; align-items: center; justify-content: space-between; margin-bottom: 15px;",
          h4(icon("images"), "Figure & Export Gallery", style = "margin: 0;")
        )
      )
    ),
    tabsetPanel(id = ns("gallery_tabs"),
      tabPanel("UMAP Gallery",
        fluidRow(
          column(12,
            div(style = "margin-top: 10px; margin-bottom: 5px;",
              actionButton(ns("clear_umap_gallery"), "Clear UMAPs",
                icon = icon("trash"), class = "btn-sm btn-warning")
            ),
            uiOutput(ns("umap_summary")),
            uiOutput(ns("umap_grid"))
          )
        )
      ),
      tabPanel("Figure Gallery",
        fluidRow(
          column(12,
            div(style = "margin-top: 10px; margin-bottom: 5px;",
              actionButton(ns("clear_figure_gallery"), "Clear Figures",
                icon = icon("trash"), class = "btn-sm btn-warning")
            ),
            uiOutput(ns("figure_summary")),
            uiOutput(ns("figure_grid"))
          )
        )
      ),
      tabPanel("CSV Library",
        fluidRow(
          column(12,
            div(style = "margin-top: 10px; margin-bottom: 5px;",
              actionButton(ns("clear_csv_library"), "Clear CSV Library",
                icon = icon("trash"), class = "btn-sm btn-warning")
            ),
            uiOutput(ns("csv_summary")),
            DT::DTOutput(ns("csv_table"))
          )
        )
      )
    )
  )
}

# -----------------------------------------------------------------------------
# Server
# -----------------------------------------------------------------------------

umapGalleryServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    umap_registered <- reactiveVal(integer(0))
    fig_registered <- reactiveVal(integer(0))

    # ======================== UMAP Gallery ========================

    output$umap_summary <- renderUI({
      gallery <- values$umap_gallery
      n <- if (is.null(gallery)) 0L else length(gallery)
      stages <- if (n > 0) length(unique(sapply(gallery, `[[`, "stage"))) else 0L

      div(class = "alert alert-info", style = "padding: 8px 15px; margin-bottom: 12px;",
        icon("chart-area"),
        sprintf(" %d snapshot%s from %d pipeline stage%s",
                n, if (n != 1) "s" else "",
                stages, if (stages != 1) "s" else ""),
        if (n == 0) span(" -- UMAPs will appear here as you run the pipeline.",
                         style = "color: #888;")
      )
    })

    output$umap_grid <- renderUI({
      gallery <- values$umap_gallery
      if (is.null(gallery) || length(gallery) == 0) {
        return(div(
          style = "text-align: center; padding: 40px; color: #999;",
          icon("images", class = "fa-3x"), br(), br(),
          p("No UMAP snapshots yet."),
          p("Run pipeline steps to automatically capture UMAP plots here.")
        ))
      }

      cards <- lapply(rev(seq_along(gallery)), function(i) {
        entry <- gallery[[i]]
        plot_output_id <- ns(paste0("uthumb_", i))
        remove_btn_id <- ns(paste0("uremove_", i))

        div(
          class = "col-md-4 col-sm-6",
          style = "margin-bottom: 15px;",
          div(
            style = "background: white; border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.1); overflow: hidden;",
            div(style = "height: 260px; overflow: hidden;",
              plotlyOutput(plot_output_id, height = "260px")
            ),
            div(style = "padding: 8px 12px; border-top: 1px solid #eee;",
              div(style = "display: flex; justify-content: space-between; align-items: center;",
                div(
                  strong(entry$label, style = "font-size: 13px;"), br(),
                  span(style = "font-size: 11px; color: #888;",
                    if (!is.na(entry$n_cells)) paste0(format(entry$n_cells, big.mark = ","), " cells | "),
                    format(entry$timestamp, "%H:%M:%S")
                  )
                ),
                actionButton(remove_btn_id, "", icon = icon("times"),
                  class = "btn-xs btn-default", style = "padding: 2px 6px; opacity: 0.5;")
              )
            )
          )
        )
      })
      do.call(tagList, cards)
    })

    observe({
      gallery <- values$umap_gallery
      if (is.null(gallery) || length(gallery) == 0) return()
      already <- umap_registered()
      new_indices <- setdiff(seq_along(gallery), already)
      if (length(new_indices) == 0) return()

      lapply(new_indices, function(i) {
        local({
          idx <- i
          output[[paste0("uthumb_", idx)]] <- renderPlotly({
            g <- isolate(values$umap_gallery)
            if (is.null(g) || idx > length(g)) return(plotly_empty())
            p <- g[[idx]]$plot
            if (is.null(p)) return(plotly_empty())
            p %>% layout(
              margin = list(l = 30, r = 10, t = 30, b = 30),
              showlegend = FALSE,
              xaxis = list(showticklabels = FALSE),
              yaxis = list(showticklabels = FALSE)
            ) %>% config(displayModeBar = FALSE)
          })
          observeEvent(input[[paste0("uremove_", idx)]], {
            g <- values$umap_gallery
            if (!is.null(g) && idx <= length(g)) {
              g[[idx]] <- NULL
              values$umap_gallery <- g
            }
          }, ignoreInit = TRUE)
        })
      })
      umap_registered(union(already, new_indices))
    })

    observeEvent(input$clear_umap_gallery, {
      values$umap_gallery <- list()
      umap_registered(integer(0))
      showNotification("UMAP gallery cleared.", type = "warning", duration = 3)
    })

    # ======================== Figure Gallery ========================

    output$figure_summary <- renderUI({
      gallery <- values$figure_gallery
      n <- if (is.null(gallery)) 0L else length(gallery)
      modules <- if (n > 0) length(unique(sapply(gallery, function(x) x$source_module))) else 0L

      div(class = "alert alert-info", style = "padding: 8px 15px; margin-bottom: 12px;",
        icon("chart-bar"),
        sprintf(" %d figure%s from %d module%s",
                n, if (n != 1) "s" else "",
                modules, if (modules != 1) "s" else ""),
        if (n == 0) span(" -- Figures from analyses will appear here.",
                         style = "color: #888;")
      )
    })

    output$figure_grid <- renderUI({
      gallery <- values$figure_gallery
      if (is.null(gallery) || length(gallery) == 0) {
        return(div(
          style = "text-align: center; padding: 40px; color: #999;",
          icon("chart-bar", class = "fa-3x"), br(), br(),
          p("No figures captured yet."),
          p("Run analyses to auto-capture figures here.")
        ))
      }

      cards <- lapply(rev(seq_along(gallery)), function(i) {
        entry <- gallery[[i]]
        plot_output_id <- ns(paste0("fthumb_", i))
        remove_btn_id <- ns(paste0("fremove_", i))

        div(
          class = "col-md-4 col-sm-6",
          style = "margin-bottom: 15px;",
          div(
            style = "background: white; border-radius: 8px; box-shadow: 0 2px 6px rgba(0,0,0,0.1); overflow: hidden;",
            div(style = "height: 260px; overflow: hidden;",
              plotlyOutput(plot_output_id, height = "260px")
            ),
            div(style = "padding: 8px 12px; border-top: 1px solid #eee;",
              div(style = "display: flex; justify-content: space-between; align-items: center;",
                div(
                  strong(entry$label, style = "font-size: 13px;"), br(),
                  span(style = "font-size: 11px; color: #888;",
                    entry$source_module, " | ",
                    format(entry$timestamp, "%H:%M:%S")
                  )
                ),
                actionButton(remove_btn_id, "", icon = icon("times"),
                  class = "btn-xs btn-default", style = "padding: 2px 6px; opacity: 0.5;")
              )
            )
          )
        )
      })
      do.call(tagList, cards)
    })

    observe({
      gallery <- values$figure_gallery
      if (is.null(gallery) || length(gallery) == 0) return()
      already <- fig_registered()
      new_indices <- setdiff(seq_along(gallery), already)
      if (length(new_indices) == 0) return()

      lapply(new_indices, function(i) {
        local({
          idx <- i
          output[[paste0("fthumb_", idx)]] <- renderPlotly({
            g <- isolate(values$figure_gallery)
            if (is.null(g) || idx > length(g)) return(plotly_empty())
            p <- g[[idx]]$plot
            if (is.null(p)) return(plotly_empty())
            if (inherits(p, "plotly")) {
              p %>% layout(
                margin = list(l = 30, r = 10, t = 30, b = 30),
                showlegend = FALSE
              ) %>% config(displayModeBar = FALSE)
            } else {
              plotly_empty()
            }
          })
          observeEvent(input[[paste0("fremove_", idx)]], {
            g <- values$figure_gallery
            if (!is.null(g) && idx <= length(g)) {
              g[[idx]] <- NULL
              values$figure_gallery <- g
            }
          }, ignoreInit = TRUE)
        })
      })
      fig_registered(union(already, new_indices))
    })

    observeEvent(input$clear_figure_gallery, {
      values$figure_gallery <- list()
      fig_registered(integer(0))
      showNotification("Figure gallery cleared.", type = "warning", duration = 3)
    })

    # ======================== CSV Library ========================

    output$csv_summary <- renderUI({
      lib <- values$csv_library
      n <- if (is.null(lib)) 0L else length(lib)

      div(class = "alert alert-info", style = "padding: 8px 15px; margin-bottom: 12px;",
        icon("file-csv"),
        sprintf(" %d CSV export%s recorded", n, if (n != 1) "s" else ""),
        if (n == 0) span(" -- CSV exports will be logged here as you run analyses.",
                         style = "color: #888;")
      )
    })

    output$csv_table <- DT::renderDT({
      lib <- values$csv_library
      if (is.null(lib) || length(lib) == 0) {
        return(DT::datatable(
          data.frame(Message = "No CSV exports recorded yet. Run analyses to auto-capture tables."),
          rownames = FALSE, options = list(dom = "t")
        ))
      }

      df <- data.frame(
        Description = vapply(lib, function(x) if (is.null(x$description)) "" else as.character(x$description), character(1)),
        Module = vapply(lib, function(x) if (is.null(x$source_module)) "" else as.character(x$source_module), character(1)),
        Rows = vapply(lib, function(x) if (is.null(x$n_rows) || is.na(x$n_rows)) "-" else format(x$n_rows, big.mark = ","), character(1)),
        Columns = vapply(lib, function(x) if (is.null(x$n_cols) || is.na(x$n_cols)) "-" else as.character(x$n_cols), character(1)),
        Time = vapply(lib, function(x) if (is.null(x$timestamp)) "" else format(x$timestamp, "%H:%M:%S"), character(1)),
        File = vapply(lib, function(x) if (is.null(x$filename)) "" else basename(x$filename), character(1)),
        stringsAsFactors = FALSE
      )

      DT::datatable(df,
        options = list(pageLength = 25, scrollX = TRUE, dom = "tip",
                       order = list(list(4, "desc")),
                       columnDefs = list(list(visible = FALSE, targets = 5))),
        rownames = FALSE
      )
    })

    observeEvent(input$clear_csv_library, {
      values$csv_library <- list()
      showNotification("CSV library cleared.", type = "warning", duration = 3)
    })

  })
}
