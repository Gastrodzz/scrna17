# Export and Reporting Module (enhanced and comprehensive)
# modules/export.R

suppressPackageStartupMessages({
  library(Seurat)
  library(openxlsx)
  library(rmarkdown)
  library(knitr)
  library(jsonlite)
  library(zip)
})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R

# ----------------------------- UI --------------------------------------------

dataExportUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Advanced Data Export"),
          p("Export your analysis results in multiple formats with comprehensive options."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("export"),

          # Export type selection
          fluidRow(
            column(6,
              selectInput(ns("export_type"), "Export Type",
                choices = c(
                  "Seurat Object" = "seurat",
                  "Analysis Results" = "results", 
                  "Complete Analysis Package" = "package",
                  "Custom Selection" = "custom"
                ),
                selected = "results")
            ),
            column(6,
              selectInput(ns("export_format"), "Format",
                choices = c(
                  "Excel Workbook (.xlsx)" = "excel",
                  "ZIP Archive" = "zip",
                  "Seurat RDS" = "rds",
                  "H5Seurat" = "h5seurat", 
                  "H5AD (AnnData)" = "h5ad",
                  "CSV Files" = "csv",
                  "JSON Package" = "json"
                ),
                selected = "excel")
            )
          ),
          
          # Seurat object export options
          conditionalPanel(
            condition = "input.export_type == 'seurat'",
            ns = ns,
            wellPanel(
              h5("Seurat Object Export Options"),
              fluidRow(
                column(12,
                  checkboxGroupInput(ns("seurat_components"), "Components to Include",
                    choices = c(
                      "Raw counts" = "counts",
                      "Normalized data" = "data", 
                      "Scaled data" = "scaled",
                      "Variable features" = "hvgs",
                      "Dimension reductions" = "reductions",
                      "Neighbor graphs" = "graphs",
                      "Metadata" = "metadata"
                    ),
                    selected = c("data", "metadata", "reductions", "hvgs"),
                    inline = TRUE
                  )
                )
              ),
              fluidRow(
                column(6,
                  textInput(ns("export_filename"), "Filename (without extension)", 
                          value = paste0("scRNA_analysis_", Sys.Date()))
                ),
                column(6,
                  checkboxInput(ns("compress_export"), "Compress output", TRUE)
                )
              )
            )
          ),
          
          # Results export options
          conditionalPanel(
            condition = "input.export_type == 'results'",
            ns = ns,
            wellPanel(
              h5("Analysis Results Export Options"),
              fluidRow(
                column(12,
                  checkboxGroupInput(ns("results_components"), "Results to Include",
                    choices = c(
                      "QC Metrics" = "qc",
                      "Clustering Results" = "clusters",
                      "Cell Type Annotations" = "annotations",
                      "Differential Expression" = "de",
                      "Pathway Enrichment" = "pathways",
                      "Dimension Reduction Coordinates" = "embeddings",
                      "Marker Genes" = "markers",
                      "Summary Statistics" = "stats"
                    ),
                    selected = c("qc", "clusters", "annotations", "de", "embeddings"),
                    inline = TRUE
                  )
                )
              )
            )
          ),
          
          # Custom export options
          conditionalPanel(
            condition = "input.export_type == 'custom'",
            ns = ns,
            wellPanel(
              h5("Custom Export Options"),
              fluidRow(
                column(6,
                  h6("Data Components"),
                  checkboxGroupInput(ns("custom_data"), "Data to Export",
                    choices = c(
                      "Expression Matrix" = "expression",
                      "Metadata" = "metadata",
                      "Variable Features" = "hvgs",
                      "PCA Embeddings" = "pca",
                      "UMAP Embeddings" = "umap",
                      "t-SNE Embeddings" = "tsne"
                    )
                  )
                ),
                column(6,
                  h6("Analysis Results"),
                  checkboxGroupInput(ns("custom_results"), "Results to Export",
                    choices = c(
                      "DE Results" = "de_results",
                      "Enrichment Results" = "enrichment",
                      "Cluster Markers" = "cluster_markers",
                      "Cell Type Scores" = "annotation_scores"
                    )
                  )
                )
              ),
              fluidRow(
                column(12,
                  h6("Visualization Exports"),
                  checkboxGroupInput(ns("custom_plots"), "Plots to Export",
                    choices = c(
                      "Dimension Reduction Plots" = "dimplots",
                      "Feature Plots" = "feature_plots", 
                      "QC Plots" = "qc_plots",
                      "Volcano Plots" = "volcano_plots",
                      "Heatmaps" = "heatmaps"
                    )
                  )
                )
              )
            )
          ),
          
          # Export execution
          fluidRow(
            column(12,
              div(style = "margin-top: 20px;",
                actionButton(ns("start_export"), "Start Export", 
                           class = "btn-success btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    # Export progress and results
    fluidRow(
      column(12,
        conditionalPanel(
          condition = sprintf("output['%s']", ns("show_export_progress")),
          ns = ns,
          wellPanel(
            h4("Export Progress"),
            verbatimTextOutput(ns("export_status")),
            br(),
            uiOutput(ns("download_links"))
          )
        )
      )
    ),
    
    # Report generation section
    fluidRow(
      column(12,
        wellPanel(
          h4("Analysis Report Generation"),
          p("Generate comprehensive HTML or PDF reports of your analysis."),
          
          fluidRow(
            column(4,
              selectInput(ns("report_type"), "Report Type",
                choices = c(
                  "Quick Summary" = "quick",
                  "Detailed Analysis" = "detailed", 
                  "Publication Ready" = "publication",
                  "Custom Report" = "custom_report"
                ),
                selected = "detailed")
            ),
            column(4,
              selectInput(ns("report_format"), "Output Format",
                choices = c("HTML" = "html", "PDF" = "pdf", "Word" = "docx"),
                selected = "html")
            ),
            column(4,
              checkboxInput(ns("include_code"), "Include R code", FALSE)
            )
          ),
          
          fluidRow(
            column(6,
              textInput(ns("report_title"), "Report Title", 
                      value = paste("scRNA-seq Analysis Report -", Sys.Date()))
            ),
            column(6,
              textInput(ns("report_author"), "Author", value = "")
            )
          ),
          
          # Report sections
          fluidRow(
            column(12,
              h5("Report Sections to Include"),
              checkboxGroupInput(ns("report_sections"), "",
                choices = c(
                  "Executive Summary" = "summary",
                  "Data Overview" = "overview",
                  "Quality Control" = "qc_section",
                  "Preprocessing Steps" = "preprocessing", 
                  "Dimensionality Reduction" = "dimred_section",
                  "Clustering Analysis" = "clustering_section",
                  "Cell Type Annotation" = "annotation_section",
                  "Differential Expression" = "de_section",
                  "Pathway Analysis" = "pathway_section",
                  "Conclusions" = "conclusions"
                ),
                selected = c("summary", "overview", "qc_section", "clustering_section", "de_section"),
                inline = TRUE
              )
            )
          ),
          
          fluidRow(
            column(12,
              div(style = "margin-top: 20px;",
                actionButton(ns("generate_report"), "Generate Report", 
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    # AI Report Generation
    fluidRow(
      column(12,
        wellPanel(
          h4("AI-Powered Report Generation"),
          p("Use AI to generate publication-ready Methods, Results, and Discussion sections based on your analysis."),
          fluidRow(
            column(8,
              checkboxGroupInput(ns("ai_report_sections"), "Report Sections:",
                choices = c("Methods" = "methods", "Results" = "results", "Discussion" = "discussion"),
                selected = c("methods", "results"),
                inline = TRUE)
            ),
            column(4,
              div(style = "margin-top: 20px;",
                actionButton(ns("ai_report_generate"), "AI Generate Report",
                             icon = icon("robot"), class = "ai-btn btn-sm", width = "100%")
              )
            )
          ),
          uiOutput(ns("ai_report_output"))
        )
      )
    ),

    # Provenance and reproducibility bundle
    fluidRow(
      column(12,
        wellPanel(
          h4("Reproducibility and Provenance"),
          p("Download the full provenance bundle, methods inventory, and audit trail for publication and peer review."),
          fluidRow(
            column(4,
              downloadButton(ns("download_provenance"), "Download Provenance JSON", class = "btn-success", width = "100%")
            ),
            column(4,
              downloadButton(ns("download_methods"), "Download Methods CSV", class = "btn-info", width = "100%")
            ),
            column(4,
              downloadButton(ns("download_audit"), "Download Audit Log CSV", class = "btn-warning", width = "100%")
            )
          )
        )
      )
    ),

    # Plot export section
    fluidRow(
      column(12,
        wellPanel(
          h4("Plot Export (PNG / PDF)"),
          p("Export the current plot as a static image file."),
          fluidRow(
            column(4,
              selectInput(ns("plot_export_format"), "Image Format",
                choices = c("PNG" = "png", "PDF" = "pdf"),
                selected = "png")
            ),
            column(4,
              numericInput(ns("plot_export_width"), "Width (inches)", value = 8, min = 2, max = 30, step = 0.5)
            ),
            column(4,
              numericInput(ns("plot_export_height"), "Height (inches)", value = 6, min = 2, max = 30, step = 0.5)
            )
          ),
          fluidRow(
            column(4,
              numericInput(ns("plot_export_dpi"), "DPI (for PNG)", value = 300, min = 72, max = 600, step = 50)
            ),
            column(8,
              div(style = "margin-top: 25px;",
                actionButton(ns("export_current_plot"), "Save Current Plot",
                             class = "btn-info", width = "100%")
              )
            )
          )
        )
      )
    ),

    # Session management
    fluidRow(
      column(12,
        wellPanel(
          h4("Session Management"),
          p("Save and restore your complete analysis session."),
          
          fluidRow(
            column(6,
              h5("Save Session"),
              textInput(ns("session_name"), "Session Name", 
                      value = paste0("scRNA_session_", Sys.Date())),
              br(),
              actionButton(ns("save_session"), "Save Complete Session", 
                         class = "btn-info", width = "100%")
            ),
            column(6,
              h5("Load Session"),
              fileInput(ns("session_file"), "Choose Session File (.rds)",
                       accept = ".rds"),
              br(),
              actionButton(ns("load_session"), "Load Session", 
                         class = "btn-warning", width = "100%")
            )
          )
        )
      )
    )
  )
}

# --------------------------- Server ------------------------------------------

dataExportServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Export progress tracking
    export_progress <- reactiveVal("")
    show_export_progress <- reactiveVal(FALSE)
    export_files <- reactiveVal(list())
    
    output$show_export_progress <- reactive({ show_export_progress() })
    outputOptions(output, "show_export_progress", suspendWhenHidden = FALSE)

    # ---- AI Report Generation Handler ----
    observeEvent(input$ai_report_generate, {
      if (is.null(values$ai_api_key) || !nzchar(values$ai_api_key)) {
        showNotification("Please set your Claude API key in the Data Management tab first.",
                         type = "warning", duration = 6)
        return()
      }

      sections <- input$ai_report_sections
      if (is.null(sections) || length(sections) == 0) {
        showNotification("Please select at least one report section.", type = "warning", duration = 4)
        return()
      }
      output$ai_report_output <- renderUI({ ai_thread_loading("AI generating report sections...") })

      # Build dataset context
      dataset_ctx <- if (exists("build_dataset_context", mode = "function")) {
        tryCatch(
          build_dataset_context(values,
            include = c("summary", "clusters", "annotations", "de", "enrichment")),
          error = function(e) "No dataset context available."
        )
      } else {
        "No dataset context available."
      }

      # Get method inventory if available
      method_info <- ""
      if (!is.null(values$method_inventory)) {
        inv <- values$method_inventory
        if (is.data.frame(inv) && nrow(inv) > 0) {
          method_info <- paste0("\n\nMETHOD INVENTORY:\n",
            paste(utils::capture.output(print(inv)), collapse = "\n"))
        }
      } else if (exists("build_method_inventory", mode = "function")) {
        inv <- tryCatch(build_method_inventory(), error = function(e) NULL)
        if (!is.null(inv) && is.data.frame(inv) && nrow(inv) > 0) {
          method_info <- paste0("\n\nMETHOD INVENTORY:\n",
            paste(utils::capture.output(print(inv)), collapse = "\n"))
        }
      }

      # Build section-specific prompts
      section_prompts <- list(
        methods = paste0(
          "Write a Methods section for a publication describing this scRNA-seq analysis. ",
          "Include: data processing, normalization, dimensionality reduction, clustering, ",
          "annotation methods. Cite tools with versions where possible."
        ),
        results = paste0(
          "Write a Results section describing the key findings from this scRNA-seq analysis. ",
          "Include: cell type composition, DE results, pathway enrichment, any trajectory findings."
        ),
        discussion = paste0(
          "Write a brief Discussion section interpreting the findings in biological context."
        )
      )

      selected_prompts <- vapply(sections, function(s) {
        paste0("## ", toupper(s), " SECTION\n", section_prompts[[s]])
      }, character(1))

      full_prompt <- paste0(
        "DATASET CONTEXT:\n", dataset_ctx, "\n",
        method_info, "\n\n",
        "Please write the following sections for a publication-quality scRNA-seq report:\n\n",
        paste(selected_prompts, collapse = "\n\n")
      )

      withProgress(message = "AI is generating report sections...", value = 0.5, {
        result <- call_ai("full_report", full_prompt, values)
      })

      if (result$ok && !is.null(result$text)) {
        report_html <- .mdtohtml(result$text)
        report_id <- paste0("ai_report_text_", as.integer(Sys.time()))
        copy_btn <- actionButton(session$ns("copy_ai_report"), "Copy to Clipboard",
                     icon = icon("clipboard"), class = "ai-thread-apply-btn",
                     onclick = paste0(
                       "var el = this.closest('.ai-thread-content');",
                       "navigator.clipboard.writeText(el.innerText).then(function(){",
                       "alert('Report copied to clipboard!');",
                       "}).catch(function(){ document.execCommand('copy'); });"
                     ))
        output$ai_report_output <- renderUI({
          ai_thread_message(
            content_html = report_html,
            title = "AI Generated Report",
            model = get_model_for_task("full_report", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out,
            apply_btn = copy_btn
          )
        })
      } else {
        output$ai_report_output <- renderUI({ ai_thread_error(result$error %||% "Unknown error") })
      }
    })

    # ---- AI Helper (generic Q&A modal for Export & Report) ----
    if (exists("ai_helper_server", mode = "function")) {
      ai_helper_server("exporter", values, "Export & Report", context_fn = function() {
        ctx <- ""
        if (!is.null(values$active_dataset)) {
          obj <- values$active_dataset
          ctx <- paste0(
            "Dataset has ", ncol(obj), " cells and ", nrow(obj), " genes.\n",
            "Available assays: ", paste(names(obj@assays), collapse = ", "), ".\n",
            "Available reductions: ", paste(names(obj@reductions), collapse = ", "), "."
          )
        }
        if (!is.null(values$last_de_results) && is.data.frame(values$last_de_results)) {
          ctx <- paste0(ctx, "\nDE results available with ", nrow(values$last_de_results), " rows.")
        }
        if (!is.null(values$last_enrich_results)) {
          ctx <- paste0(ctx, "\nEnrichment results available.")
        }
        ctx
      })
    }

    output$download_provenance <- downloadHandler(
      filename = function() paste0("provenance_bundle_", Sys.Date(), ".json"),
      content = function(file) {
        if (exists("write_provenance_bundle", mode = "function")) {
          write_provenance_bundle(values, file)
        } else if (requireNamespace("jsonlite", quietly = TRUE)) {
          jsonlite::write_json(
            list(
              generated_at = Sys.time(),
              dataset_version = values$dataset_version %||% NA,
              dataset_signature = values$dataset_signature %||% NA,
              environment = list(r_version = R.version.string)
            ),
            path = file,
            pretty = TRUE,
            auto_unbox = TRUE
          )
        } else {
          writeLines("{}", con = file)
        }
      }
    )

    output$download_methods <- downloadHandler(
      filename = function() paste0("methods_inventory_", Sys.Date(), ".csv"),
      content = function(file) {
        inv <- values$method_inventory %||% if (exists("build_method_inventory", mode = "function")) build_method_inventory() else data.frame()
        write.csv(inv, file, row.names = FALSE)
      }
    )

    output$download_audit <- downloadHandler(
      filename = function() paste0("analysis_audit_log_", Sys.Date(), ".csv"),
      content = function(file) {
        log_df <- values$analysis_audit %||% data.frame()
        write.csv(log_df, file, row.names = FALSE)
      }
    )
    
    # Main export handler
    observeEvent(input$start_export, {
      req(values$active_dataset)
      
      show_export_progress(TRUE)
      export_progress("Starting export process...")
      
      withProgress(message = "Exporting data...", value = 0, {
        
        export_result <- perform_export(values, input)
        
        if (export_result$success) {
          export_files(export_result$files)
          export_progress(paste("✅ Export completed successfully!", 
                               length(export_result$files), "files generated"))
        } else {
          export_progress(paste("❌ Export failed:", export_result$error))
        }
        
        incProgress(1)
      })
    })
    
    # Export status display
    output$export_status <- renderText({
      export_progress()
    })
    
    # Download links
    output$download_links <- renderUI({
      files <- export_files()
      if (length(files) > 0) {
        download_buttons <- lapply(seq_along(files), function(i) {
          file_info <- files[[i]]
          downloadButton(
            session$ns(paste0("download_", i)), 
            paste("📥 Download", file_info$name),
            class = "btn-success",
            style = "margin: 5px; width: 300px;"
          )
        })
        
        do.call(tagList, download_buttons)
      }
    })
    
    # Dynamic download handlers
    observe({
      files <- export_files()
      if (length(files) > 0) {
        lapply(seq_along(files), function(i) {
          file_info <- files[[i]]
          output_name <- paste0("download_", i)
          
          output[[output_name]] <- downloadHandler(
            filename = function() file_info$name,
            content = function(file) {
              file.copy(file_info$path, file, overwrite = TRUE)
            }
          )
        })
      }
    })
    
    # Report generation
    observeEvent(input$generate_report, {
      req(values$active_dataset)
      
      withProgress(message = "Generating report...", value = 0, {
        
        incProgress(0.2, detail = "Preparing report template...")
        
        report_result <- generate_analysis_report(values, input)
        
        if (report_result$success) {
          # Add report to download files
          current_files <- export_files()
          current_files[[length(current_files) + 1]] <- list(
            name = report_result$filename,
            path = report_result$filepath
          )
          export_files(current_files)
          
          export_progress(paste("✅ Report generated successfully:", report_result$filename))
          show_export_progress(TRUE)
          
          showNotification("Report generated successfully!", type = "message")
        } else {
          showNotification(paste("Report generation failed:", report_result$error), type = "error")
        }
        
        incProgress(1)
      })
    })
    
    # ---- Plot Export Handler ----
    observeEvent(input$export_current_plot, {
      # Attempt to export the last rendered plot
      fmt    <- input$plot_export_format %||% "png"
      width  <- input$plot_export_width  %||% 8
      height <- input$plot_export_height %||% 6
      dpi    <- input$plot_export_dpi    %||% 300

      out_file <- file.path(tempdir(), paste0("exported_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".", fmt))

      tryCatch({
        plot_obj <- values$last_plot  # assume modules store their last plot here

        if (!is.null(plot_obj)) {
          saved <- FALSE

          # 1. If it is a plotly object, try plotly::save_image (requires orca / kaleido)
          if (inherits(plot_obj, "plotly")) {
            if (requireNamespace("plotly", quietly = TRUE)) {
              tryCatch({
                plotly::save_image(plot_obj, file = out_file, width = width * dpi, height = height * dpi)
                saved <- TRUE
              }, error = function(e) {
                message("plotly::save_image failed: ", e$message)
              })
            }
            # Fallback: try webshot2 for plotly
            if (!saved && requireNamespace("webshot2", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
              tmp_html <- tempfile(fileext = ".html")
              htmlwidgets::saveWidget(plotly::as_widget(plot_obj), tmp_html, selfcontained = TRUE)
              webshot2::webshot(tmp_html, file = out_file, vwidth = width * dpi, vheight = height * dpi)
              unlink(tmp_html)
              saved <- TRUE
            }
          }

          # 2. If it is a ggplot, use ggsave
          if (!saved && inherits(plot_obj, "gg")) {
            ggplot2::ggsave(out_file, plot = plot_obj, width = width, height = height, dpi = dpi, device = fmt)
            saved <- TRUE
          }

          # 3. Recordedplot or other base-R plot -- use the appropriate device
          if (!saved && inherits(plot_obj, "recordedplot")) {
            if (fmt == "png") {
              grDevices::png(out_file, width = width, height = height, units = "in", res = dpi)
            } else {
              grDevices::pdf(out_file, width = width, height = height)
            }
            grDevices::replayPlot(plot_obj)
            grDevices::dev.off()
            saved <- TRUE
          }

          if (saved && file.exists(out_file)) {
            current_files <- export_files()
            current_files[[length(current_files) + 1]] <- list(name = basename(out_file), path = out_file)
            export_files(current_files)
            show_export_progress(TRUE)
            export_progress(paste("Plot saved:", basename(out_file)))
            showNotification(paste("Plot exported as", toupper(fmt)), type = "message")
          } else {
            showNotification("Could not export the plot. Ensure plotly (with kaleido), ggplot2, or webshot2 is installed.", type = "warning", duration = 8)
          }
        } else {
          showNotification("No plot available to export. Generate a plot first.", type = "warning")
        }
      }, error = function(e) {
        showNotification(paste("Plot export failed:", e$message), type = "error")
      })
    })

    # Session management
    observeEvent(input$save_session, {
      req(values$active_dataset)
      
      withProgress(message = "Saving session...", {
        session_result <- save_analysis_session(values, input$session_name)
        
        if (session_result$success) {
          # Add session file to downloads
          current_files <- export_files()
          current_files[[length(current_files) + 1]] <- list(
            name = session_result$filename,
            path = session_result$filepath
          )
          export_files(current_files)
          
          export_progress(paste("💾 Session saved:", session_result$filename))
          show_export_progress(TRUE)
          showNotification("Session saved successfully!", type = "message")
        } else {
          showNotification(paste("Session save failed:", session_result$error), type = "error")
        }
      })
    })
    
    observeEvent(input$load_session, {
      req(input$session_file)
      
      withProgress(message = "Loading session...", {
        session_result <- load_analysis_session(input$session_file$datapath)
        
        if (session_result$success) {
          # Restore session data to values
          names_to_restore <- names(session_result$data)
          for (name in names_to_restore) {
            values[[name]] <- session_result$data[[name]]
          }
          
          showNotification("Session loaded successfully.", type = "message")
        } else {
          showNotification(paste("Session load failed:", session_result$error), type = "error")
        }
      })
    })
  })
}

# ============================================================================
# EXPORT HELPER FUNCTIONS
# ============================================================================

perform_export <- function(values, input) {
  tryCatch({
    export_type <- input$export_type
    export_format <- input$export_format

    # Create temporary directory for export files
    temp_dir <- tempfile()
    dir.create(temp_dir, recursive = TRUE)
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
    
    files_created <- list()
    
    if (export_type == "seurat") {
      # Export Seurat object
      result <- export_seurat_object(values$active_dataset, input, temp_dir)
      files_created <- c(files_created, result$files)
      
    } else if (export_type == "results") {
      # Export analysis results
      result <- export_analysis_results(values, input, temp_dir)
      files_created <- c(files_created, result$files)
      
    } else if (export_type == "package") {
      # Export complete analysis package
      result <- export_complete_package(values, input, temp_dir)
      files_created <- c(files_created, result$files)
      
    } else if (export_type == "custom") {
      # Export custom selection
      result <- export_custom_selection(values, input, temp_dir)
      files_created <- c(files_created, result$files)
    }
    
    list(success = TRUE, files = files_created, error = NULL)
    
  }, error = function(e) {
    list(success = FALSE, files = list(), error = e$message)
  })
}

export_seurat_object <- function(seurat_obj, params, output_dir) {
  files_created <- list()
  
  # Prepare object based on selected components
  export_obj <- seurat_obj
  
  # Remove unwanted components (Seurat v4/v5 compatible)
  if (!"counts" %in% params$seurat_components) {
    for (assay_name in names(export_obj@assays)) {
      assay_obj <- export_obj[[assay_name]]
      if (inherits(assay_obj, "Assay5")) {
        lyrs <- SeuratObject::Layers(assay_obj)
        count_lyrs <- grep("^counts", lyrs, value = TRUE)
        for (cl in count_lyrs) {
          tryCatch(SeuratObject::LayerData(export_obj[[assay_name]], layer = cl) <- NULL, error = function(e) NULL)
        }
      } else {
        export_obj@assays[[assay_name]]@counts <- new("dgCMatrix")
      }
    }
  }

  if (!"scaled" %in% params$seurat_components) {
    for (assay_name in names(export_obj@assays)) {
      assay_obj <- export_obj[[assay_name]]
      if (inherits(assay_obj, "Assay5")) {
        lyrs <- SeuratObject::Layers(assay_obj)
        scale_lyrs <- grep("^scale\\.data", lyrs, value = TRUE)
        for (sl in scale_lyrs) {
          tryCatch(SeuratObject::LayerData(export_obj[[assay_name]], layer = sl) <- NULL, error = function(e) NULL)
        }
      } else {
        export_obj@assays[[assay_name]]@scale.data <- matrix()
      }
    }
  }
  
  if (!"reductions" %in% params$seurat_components) {
    export_obj@reductions <- list()
  }
  
  if (!"graphs" %in% params$seurat_components) {
    export_obj@graphs <- list()
  }
  
  if (!"hvgs" %in% params$seurat_components) {
    for (assay_name in names(export_obj@assays)) {
      VariableFeatures(export_obj[[assay_name]]) <- character(0)
    }
  }
  
  # Export in requested format
  filename <- paste0(params$export_filename, 
                    switch(params$export_format,
                           "rds" = ".rds",
                           "h5seurat" = ".h5seurat", 
                           "h5ad" = ".h5ad",
                           "csv" = "_data.csv",
                           "json" = ".json",
                           ".rds"))
  
  filepath <- file.path(output_dir, filename)
  
  if (params$export_format == "rds") {
    saveRDS(export_obj, filepath)
    
  } else if (params$export_format == "h5seurat") {
    if (requireNamespace("SeuratDisk", quietly = TRUE)) {
      SeuratDisk::SaveH5Seurat(export_obj, filename = filepath, overwrite = TRUE)
    } else {
      stop("SeuratDisk package required for H5Seurat export")
    }
    
  } else if (params$export_format == "h5ad") {
    if (requireNamespace("SeuratDisk", quietly = TRUE)) {
      temp_h5seurat <- gsub("\\.h5ad$", ".h5seurat", filepath)
      SeuratDisk::SaveH5Seurat(export_obj, filename = temp_h5seurat, overwrite = TRUE)
      SeuratDisk::Convert(temp_h5seurat, dest = "h5ad", overwrite = TRUE)
      unlink(temp_h5seurat)
    } else {
      stop("SeuratDisk package required for H5AD export")
    }
    
  } else if (params$export_format == "csv") {
    # Export expression data -- prefer sparse Matrix Market format
    expr_data <- safe_get_assay_data(export_obj, "data")

    # Default: write as Matrix Market (sparse) + gzip
    mm_file <- sub("\\.csv$", ".mtx", filepath)
    Matrix::writeMM(expr_data, mm_file)

    # Also write gene and barcode name files for completeness
    writeLines(rownames(expr_data), sub("\\.csv$", "_genes.tsv", filepath))
    writeLines(colnames(expr_data), sub("\\.csv$", "_barcodes.tsv", filepath))

    # Optionally write HDF5 if HDF5Array is available
    h5_file <- sub("\\.csv$", ".h5", filepath)
    if (requireNamespace("HDF5Array", quietly = TRUE)) {
      tryCatch({
        HDF5Array::writeHDF5Array(methods::as(expr_data, "HDF5Array"), filepath = h5_file, name = "expression")
        message("HDF5 export written to ", h5_file)
      }, error = function(e) {
        message("HDF5 export skipped: ", e$message)
      })
    }

    # If user explicitly requested CSV, limit to top N variable features with a warning
    n_cells <- ncol(expr_data)
    n_genes <- nrow(expr_data)
    est_mem_gb <- (n_cells * n_genes * 8) / (1024^3)
    top_n <- min(2000, n_genes)
    var_feats <- VariableFeatures(export_obj)
    if (length(var_feats) > 0) {
      top_n <- min(top_n, length(var_feats))
      selected_genes <- var_feats[seq_len(top_n)]
    } else {
      selected_genes <- rownames(expr_data)[seq_len(top_n)]
    }
    warning(sprintf(
      "CSV export limited to top %d variable features (full matrix would be ~%.1f GB in memory). Sparse .mtx file contains the complete matrix.",
      top_n, est_mem_gb
    ))
    write.csv(as.matrix(expr_data[selected_genes, ]), filepath)

    # Update filepath to return the mtx as the primary file
    filepath <- mm_file

  } else if (params$export_format == "json") {
    # Export as JSON
    obj_list <- list(
      metadata = export_obj@meta.data,
      reductions = lapply(export_obj@reductions, function(x) Embeddings(x)),
      assays = sapply(names(export_obj@assays), function(a) {
        list(
          counts = if("counts" %in% params$seurat_components) as.matrix(safe_get_assay_data(export_obj, "counts", assay = a)) else NULL,
          data = if("data" %in% params$seurat_components) as.matrix(safe_get_assay_data(export_obj, "data", assay = a)) else NULL
        )
      }, simplify = FALSE)
    )
    write_json(obj_list, filepath, auto_unbox = TRUE)
  }
  
  files_created[[1]] <- list(name = filename, path = filepath)
  
  list(files = files_created)
}

export_analysis_results <- function(values, params, output_dir) {
  files_created <- list()
  
  if (params$export_format == "excel") {
    # Create Excel workbook
    wb <- createWorkbook()
    
    # Add worksheets based on selected components
    if ("qc" %in% params$results_components && !is.null(values$active_dataset)) {
      add_qc_worksheet(wb, values$active_dataset)
    }
    
    if ("clusters" %in% params$results_components && !is.null(values$active_dataset)) {
      add_clustering_worksheet(wb, values$active_dataset)
    }
    
    if ("annotations" %in% params$results_components && !is.null(values$active_dataset)) {
      add_annotations_worksheet(wb, values$active_dataset)
    }
    
    if ("de" %in% params$results_components && !is.null(values$last_de_results)) {
      add_de_worksheet(wb, values$last_de_results)
    }
    
    if ("pathways" %in% params$results_components && !is.null(values$last_enrich_results)) {
      add_enrichment_worksheet(wb, values$last_enrich_results)
    }
    
    if ("embeddings" %in% params$results_components && !is.null(values$active_dataset)) {
      add_embeddings_worksheet(wb, values$active_dataset)
    }
    
    if ("stats" %in% params$results_components && !is.null(values$active_dataset)) {
      add_summary_worksheet(wb, values$active_dataset, values)
    }
    
    # Save workbook
    filename <- "analysis_results.xlsx"
    filepath <- file.path(output_dir, filename)
    saveWorkbook(wb, filepath, overwrite = TRUE)
    
    files_created[[1]] <- list(name = filename, path = filepath)
    
  } else if (params$export_format == "csv") {
    # Export as separate CSV files
    csv_files <- export_results_as_csv(values, params, output_dir)
    files_created <- c(files_created, csv_files)
    
  } else if (params$export_format == "zip") {
    # Create ZIP archive
    zip_files <- export_results_as_zip(values, params, output_dir)
    files_created <- c(files_created, zip_files)
  }
  
  list(files = files_created)
}

export_complete_package <- function(values, params, output_dir) {
  # Export everything as a comprehensive package
  files_created <- list()
  
  # Create subdirectories
  data_dir <- file.path(output_dir, "data")
  results_dir <- file.path(output_dir, "results")
  plots_dir <- file.path(output_dir, "plots")
  
  dir.create(data_dir, recursive = TRUE)
  dir.create(results_dir, recursive = TRUE) 
  dir.create(plots_dir, recursive = TRUE)
  
  # Export Seurat object
  seurat_file <- file.path(data_dir, "seurat_object.rds")
  saveRDS(values$active_dataset, seurat_file)
  files_created[[length(files_created) + 1]] <- list(name = "data/seurat_object.rds", path = seurat_file)

  # Export provenance and methods for reproducibility
  prov_file <- file.path(results_dir, "provenance_bundle.json")
  if (exists("write_provenance_bundle", mode = "function")) {
    write_provenance_bundle(values, prov_file)
  } else if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(
      list(
        generated_at = Sys.time(),
        dataset_version = values$dataset_version %||% NA,
        dataset_signature = values$dataset_signature %||% NA
      ),
      path = prov_file,
      pretty = TRUE,
      auto_unbox = TRUE
    )
  } else {
    writeLines("{}", prov_file)
  }
  files_created[[length(files_created) + 1]] <- list(name = "results/provenance_bundle.json", path = prov_file)

  methods_file <- file.path(results_dir, "methods_inventory.csv")
  methods_tbl <- values$method_inventory %||% if (exists("build_method_inventory", mode = "function")) build_method_inventory() else data.frame()
  write.csv(methods_tbl, methods_file, row.names = FALSE)
  files_created[[length(files_created) + 1]] <- list(name = "results/methods_inventory.csv", path = methods_file)

  audit_file <- file.path(results_dir, "analysis_audit_log.csv")
  write.csv(values$analysis_audit %||% data.frame(), audit_file, row.names = FALSE)
  files_created[[length(files_created) + 1]] <- list(name = "results/analysis_audit_log.csv", path = audit_file)
  
  # Export all analysis results
  if (!is.null(values$last_de_results)) {
    de_file <- file.path(results_dir, "differential_expression.csv")
    write.csv(values$last_de_results, de_file, row.names = FALSE)
    files_created[[length(files_created) + 1]] <- list(name = "results/differential_expression.csv", path = de_file)
  }
  
  if (!is.null(values$last_enrich_results)) {
    if (!is.null(values$last_enrich_results$ora)) {
      ora_file <- file.path(results_dir, "pathway_enrichment_ora.csv")
      combined_ora <- rbind(
        if (!is.null(values$last_enrich_results$ora$up)) values$last_enrich_results$ora$up,
        if (!is.null(values$last_enrich_results$ora$down)) values$last_enrich_results$ora$down
      )
      write.csv(combined_ora, ora_file, row.names = FALSE)
      files_created[[length(files_created) + 1]] <- list(name = "results/pathway_enrichment_ora.csv", path = ora_file)
    }
    
    if (!is.null(values$last_enrich_results$gsea)) {
      gsea_file <- file.path(results_dir, "pathway_enrichment_gsea.csv")
      write.csv(values$last_enrich_results$gsea, gsea_file, row.names = FALSE)
      files_created[[length(files_created) + 1]] <- list(name = "results/pathway_enrichment_gsea.csv", path = gsea_file)
    }
  }
  
  # Create final ZIP archive
  zip_filename <- "complete_analysis_package.zip"
  zip_filepath <- file.path(dirname(output_dir), zip_filename)
  
  zip::zip(zip_filepath, files = output_dir, mode = "cherry-pick")
  
  files_created <- list(list(name = zip_filename, path = zip_filepath))
  
  list(files = files_created)
}

export_custom_selection <- function(values, params, output_dir) {
  files_created <- list()
  
  obj <- values$active_dataset
  
  # Export selected data components
  if ("expression" %in% params$custom_data) {
    expr_data <- safe_get_assay_data(obj, "data")

    # Primary export: sparse Matrix Market format
    mtx_file <- file.path(output_dir, "expression_matrix.mtx")
    Matrix::writeMM(expr_data, mtx_file)
    writeLines(rownames(expr_data), file.path(output_dir, "expression_matrix_genes.tsv"))
    writeLines(colnames(expr_data), file.path(output_dir, "expression_matrix_barcodes.tsv"))
    files_created[[length(files_created) + 1]] <- list(name = "expression_matrix.mtx", path = mtx_file)

    # Optional HDF5 export
    if (requireNamespace("HDF5Array", quietly = TRUE)) {
      h5_file <- file.path(output_dir, "expression_matrix.h5")
      tryCatch({
        HDF5Array::writeHDF5Array(methods::as(expr_data, "HDF5Array"), filepath = h5_file, name = "expression")
        files_created[[length(files_created) + 1]] <- list(name = "expression_matrix.h5", path = h5_file)
      }, error = function(e) message("HDF5 export skipped: ", e$message))
    }

    # Also write a CSV limited to top variable features with memory warning
    csv_file <- file.path(output_dir, "expression_matrix.csv")
    n_cells <- ncol(expr_data)
    n_genes <- nrow(expr_data)
    est_mem_gb <- (n_cells * n_genes * 8) / (1024^3)
    top_n <- min(2000, n_genes)
    var_feats <- VariableFeatures(obj)
    if (length(var_feats) > 0) {
      top_n <- min(top_n, length(var_feats))
      selected_genes <- var_feats[seq_len(top_n)]
    } else {
      selected_genes <- rownames(expr_data)[seq_len(top_n)]
    }
    warning(sprintf(
      "CSV export limited to top %d variable features (full matrix ~%.1f GB). Use .mtx or .h5 for complete data.",
      top_n, est_mem_gb
    ))
    write.csv(as.matrix(expr_data[selected_genes, ]), csv_file)
    files_created[[length(files_created) + 1]] <- list(name = "expression_matrix.csv", path = csv_file)
  }
  
  if ("metadata" %in% params$custom_data) {
    meta_file <- file.path(output_dir, "cell_metadata.csv")
    write.csv(obj@meta.data, meta_file, row.names = TRUE)
    files_created[[length(files_created) + 1]] <- list(name = "cell_metadata.csv", path = meta_file)
  }
  
  # Export embeddings
  for (reduction in c("pca", "umap", "tsne")) {
    if (reduction %in% params$custom_data && reduction %in% names(obj@reductions)) {
      emb_file <- file.path(output_dir, paste0(reduction, "_coordinates.csv"))
      emb_data <- Embeddings(obj, reduction)
      write.csv(emb_data, emb_file, row.names = TRUE)
      files_created[[length(files_created) + 1]] <- list(name = paste0(reduction, "_coordinates.csv"), path = emb_file)
    }
  }
  
  # Export analysis results
  if ("de_results" %in% params$custom_results && !is.null(values$last_de_results)) {
    de_file <- file.path(output_dir, "de_results.csv")
    write.csv(values$last_de_results, de_file, row.names = FALSE)
    files_created[[length(files_created) + 1]] <- list(name = "de_results.csv", path = de_file)
  }
  
  if ("enrichment" %in% params$custom_results && !is.null(values$last_enrich_results)) {
    enrich_file <- file.path(output_dir, "enrichment_results.json")
    write_json(values$last_enrich_results, enrich_file, auto_unbox = TRUE)
    files_created[[length(files_created) + 1]] <- list(name = "enrichment_results.json", path = enrich_file)
  }
  
  list(files = files_created)
}

# ============================================================================
# REPORT GENERATION FUNCTIONS
# ============================================================================

generate_analysis_report <- function(values, params) {
  tryCatch({
    # Create temporary directory for report generation
    temp_dir <- tempfile()
    dir.create(temp_dir, recursive = TRUE)
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
    
    # Generate report based on type
    if (params$report_type == "quick") {
      report_result <- generate_quick_report(values, params, temp_dir)
    } else if (params$report_type == "detailed") {
      report_result <- generate_detailed_report(values, params, temp_dir)
    } else if (params$report_type == "publication") {
      report_result <- generate_publication_report(values, params, temp_dir)
    } else {
      report_result <- generate_custom_report(values, params, temp_dir)
    }
    
    report_result
    
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

generate_quick_report <- function(values, params, output_dir) {
  # Generate a quick summary report
  obj <- values$active_dataset
  
  # Create HTML content
  html_content <- paste(
    "<html><head><title>", params$report_title, "</title>",
    "<style>body{font-family:Arial,sans-serif;margin:40px;} table{border-collapse:collapse;width:100%;} th,td{border:1px solid #ddd;padding:8px;} th{background-color:#f2f2f2;}</style>",
    "</head><body>",
    "<h1>", params$report_title, "</h1>",
    if(nzchar(params$report_author)) paste("<p><strong>Author:</strong>", params$report_author, "</p>"),
    "<p><strong>Generated:</strong>", Sys.time(), "</p>",
    "<hr>",
    
    # Dataset overview
    "<h2>Dataset Overview</h2>",
    "<table>",
    "<tr><th>Metric</th><th>Value</th></tr>",
    "<tr><td>Total Cells</td><td>", ncol(obj), "</td></tr>",
    "<tr><td>Total Genes</td><td>", nrow(obj), "</td></tr>",
    "<tr><td>Median UMI count</td><td>", round(median(obj$nCount_RNA, na.rm = TRUE)), "</td></tr>",
    "<tr><td>Median Gene count</td><td>", round(median(obj$nFeature_RNA, na.rm = TRUE)), "</td></tr>",
    if ("percent.mt" %in% colnames(obj@meta.data)) paste("<tr><td>Median MT %</td><td>", round(median(obj$percent.mt, na.rm = TRUE), 1), "</td></tr>"),
    "</table>",
    
    # Analysis summary
    "<h2>Analysis Summary</h2>",
    "<ul>",
    "<li>Assays available: ", paste(names(obj@assays), collapse = ", "), "</li>",
    "<li>Dimensionality reductions: ", paste(names(obj@reductions), collapse = ", "), "</li>",
    if (!is.null(values$last_de_results)) paste("<li>Differential expression: ", nrow(values$last_de_results), " genes tested</li>"),
    "</ul>",
    
    "</body></html>",
    sep = ""
  )
  
  # Write HTML file
  filename <- paste0("quick_report_", Sys.Date(), ".html")
  filepath <- file.path(output_dir, filename)
  writeLines(html_content, filepath)
  
  list(success = TRUE, filename = filename, filepath = filepath)
}

generate_detailed_report <- function(values, params, output_dir) {

  # Guard: rmarkdown must be available

  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    message("rmarkdown package not available -- falling back to quick HTML report.")
    return(generate_quick_report(values, params, output_dir))
  }

  # Generate detailed analysis report using parameterized R Markdown

  # Build a parameterized Rmd template
  rmd_content <- create_detailed_rmd_template(values, params)

  # Write Rmd file
  rmd_filename <- "detailed_report.Rmd"
  rmd_filepath <- file.path(output_dir, rmd_filename)
  writeLines(rmd_content, rmd_filepath)

  # Serialize values needed by the template so the render environment can load them
  params_rds <- file.path(output_dir, ".report_params.rds")
  saveRDS(list(
    n_cells        = ncol(values$active_dataset),
    n_genes        = nrow(values$active_dataset),
    median_umi     = round(median(values$active_dataset$nCount_RNA, na.rm = TRUE)),
    median_features = round(median(values$active_dataset$nFeature_RNA, na.rm = TRUE)),
    assays         = paste(names(values$active_dataset@assays), collapse = ", "),
    reductions     = paste(names(values$active_dataset@reductions), collapse = ", "),
    n_clusters     = if ("seurat_clusters" %in% colnames(values$active_dataset@meta.data))
                       length(unique(values$active_dataset$seurat_clusters)) else NA,
    de_n_genes     = if (!is.null(values$last_de_results)) nrow(values$last_de_results) else NA,
    qc_summary     = summary(values$active_dataset@meta.data[, intersect(
                       c("nCount_RNA", "nFeature_RNA", "percent.mt"), colnames(values$active_dataset@meta.data)),
                       drop = FALSE]),
    report_title   = params$report_title,
    report_author  = params$report_author,
    report_sections = params$report_sections,
    dataset_version = values$dataset_version %||% "NA",
    dataset_signature = values$dataset_signature %||% "NA"
  ), params_rds)

  # Render report
  output_format <- switch(params$report_format,
    "html" = "html_document",
    "pdf" = "pdf_document",
    "docx" = "word_document",
    "html_document"
  )

  output_filename <- paste0("detailed_report_", Sys.Date(), ".",
                           switch(params$report_format, "html" = "html", "pdf" = "pdf", "docx" = "docx", "html"))

  tryCatch({
    rmarkdown::render(
      rmd_filepath,
      output_format = output_format,
      output_file = output_filename,
      output_dir = output_dir,
      quiet = TRUE,
      envir = new.env()
    )

    list(success = TRUE, filename = output_filename, filepath = file.path(output_dir, output_filename))

  }, error = function(e) {
    # Fallback to HTML generation if RMarkdown fails
    html_result <- generate_quick_report(values, params, output_dir)
    html_result$filename <- paste0("report_fallback_", Sys.Date(), ".html")
    html_result
  })
}

create_detailed_rmd_template <- function(values, params) {
  obj <- values$active_dataset
  
  rmd_template <- paste0(
    "---\n",
    "title: '", params$report_title, "'\n",
    "author: '", params$report_author, "'\n", 
    "date: '", Sys.Date(), "'\n",
    "output: html_document\n",
    "---\n\n",
    
    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(echo = ", if(params$include_code) "TRUE" else "FALSE", ", warning = FALSE, message = FALSE)\n",
    "library(Seurat)\n",
    "library(ggplot2)\n",
    "library(plotly)\n",
    "library(DT)\n",
    "```\n\n",
    
    if ("summary" %in% params$report_sections) {
      paste(
        "# Executive Summary\n\n",
        "This report presents the analysis of single-cell RNA sequencing data containing **", ncol(obj), " cells** and **", nrow(obj), " genes**.\n\n",
        "Key findings:\n",
        "- Dataset contains ", ncol(obj), " cells across ", nrow(obj), " genes\n",
        "- Median UMI count: ", round(median(obj$nCount_RNA, na.rm = TRUE)), "\n",
        "- Available assays: ", paste(names(obj@assays), collapse = ", "), "\n",
        if (length(obj@reductions) > 0) paste("- Dimensionality reductions computed: ", paste(names(obj@reductions), collapse = ", "), "\n"),
        "\n"
      )
    } else "",
    
    if ("overview" %in% params$report_sections) {
      paste(
        "# Data Overview\n\n",
        "```{r data-overview}\n",
        "# Dataset summary table\n",
        "summary_data <- data.frame(\n",
        "  Metric = c('Total Cells', 'Total Genes', 'Median UMI', 'Median Genes'),\n",
        "  Value = c(", ncol(obj), ", ", nrow(obj), ", ", 
        round(median(obj$nCount_RNA, na.rm = TRUE)), ", ",
        round(median(obj$nFeature_RNA, na.rm = TRUE)), ")\n",
        ")\n",
        "knitr::kable(summary_data, caption = 'Dataset Summary')\n",
        "```\n\n"
      )
    } else "",
    
    if ("qc_section" %in% params$report_sections) {
      "# Quality Control\n\nQuality control metrics were computed for all cells in the dataset.\n\n"
    } else "",
    
    if ("clustering_section" %in% params$report_sections) {
      paste(
        "# Clustering Analysis\n\n",
        if ("seurat_clusters" %in% colnames(obj@meta.data)) {
          paste("The dataset was clustered into **", length(unique(obj$seurat_clusters)), "** clusters using graph-based clustering.\n\n")
        } else {
          "Clustering analysis was not performed or results are not available.\n\n"
        }
      )
    } else "",
    
    if ("de_section" %in% params$report_sections && !is.null(values$last_de_results)) {
      paste(
        "# Differential Expression Analysis\n\n",
        "Differential expression analysis identified **", nrow(values$last_de_results), "** genes tested.\n\n"
      )
    } else "",

    if (!is.null(values$method_inventory) && nrow(values$method_inventory) > 0) {
      method_lines <- unique(paste0(values$method_inventory$method, " - ", values$method_inventory$citation))
      method_lines <- utils::head(method_lines, 25)
      paste(
        "# Methods and Citations\n\n",
        paste0("- ", method_lines, collapse = "\n"),
        "\n\n"
      )
    } else "",

    "# Reproducibility Metadata\n\n",
    "- Dataset version: ", values$dataset_version %||% "NA", "\n",
    "- Dataset signature: ", values$dataset_signature %||% "NA", "\n",
    "- Report generated at: ", Sys.time(), "\n\n",
    
    "# Session Information\n\n",
    "```{r session-info}\n",
    "sessionInfo()\n",
    "```\n"
  )
  
  return(rmd_template)
}

generate_publication_report <- function(values, params, output_dir) {
  # Generate publication-ready report
  # This would be similar to detailed but with more formal formatting
  generate_detailed_report(values, params, output_dir)
}

generate_custom_report <- function(values, params, output_dir) {
  # Generate custom report based on selected sections
  generate_detailed_report(values, params, output_dir)
}

# ============================================================================
# SESSION MANAGEMENT FUNCTIONS
# ============================================================================

save_analysis_session <- function(values, session_name) {
  tryCatch({
    # Helper: join Seurat v5 split layers so saveRDS works
    .prep_for_save <- function(obj) {
      if (is.null(obj)) return(NULL)
      tryCatch({
        for (assay_name in names(obj@assays)) {
          a <- obj[[assay_name]]
          if (inherits(a, "Assay5") && length(SeuratObject::Layers(a)) > 1) {
            obj[[assay_name]] <- JoinLayers(a)
          }
        }
        obj
      }, error = function(e) obj)
    }

    # Create session data object
    session_data <- list(
      active_dataset = .prep_for_save(values$active_dataset),
      active_name = values$active_name,
      seurat_list = lapply(values$seurat_list, .prep_for_save),
      last_de_results = values$last_de_results,
      last_enrich_results = values$last_enrich_results,
      normalization_done = values$normalization_done,
      features_selected = values$features_selected,
      clustering_done = values$clustering_done,
      cluster_colors = values$cluster_colors,
      dataset_version = values$dataset_version,
      dataset_signature = values$dataset_signature,
      dataset_history = values$dataset_history,
      analysis_audit = values$analysis_audit,
      method_inventory = values$method_inventory,
      session_info = list(
        saved_at = Sys.time(),
        session_name = session_name,
        r_version = R.version.string
      )
    )

    # Save to temporary file
    filename <- paste0(session_name, ".rds")
    filepath <- file.path(tempdir(), filename)
    saveRDS(session_data, filepath, compress = "gzip")
    gc()

    list(success = TRUE, filename = filename, filepath = filepath)

  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

load_analysis_session <- function(session_filepath) {
  tryCatch({
    # Load session data
    session_data <- readRDS(session_filepath)
    
    # Validate session data
    required_fields <- c("active_dataset", "seurat_list", "session_info")
    if (!all(required_fields %in% names(session_data))) {
      stop("Invalid session file format")
    }
    
    list(success = TRUE, data = session_data)
    
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

# ============================================================================
# WORKSHEET HELPER FUNCTIONS
# ============================================================================

add_qc_worksheet <- function(wb, seurat_obj) {
  addWorksheet(wb, "QC_Metrics")
  
  qc_data <- seurat_obj@meta.data %>%
    select(any_of(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")))
  
  if (ncol(qc_data) > 0) {
    writeData(wb, "QC_Metrics", qc_data, rowNames = TRUE)
  }
}

add_clustering_worksheet <- function(wb, seurat_obj) {
  addWorksheet(wb, "Clustering")
  
  cluster_cols <- grep("cluster|leiden|louvain", colnames(seurat_obj@meta.data), 
                      value = TRUE, ignore.case = TRUE)
  
  if (length(cluster_cols) > 0) {
    cluster_data <- seurat_obj@meta.data[, cluster_cols, drop = FALSE]
    writeData(wb, "Clustering", cluster_data, rowNames = TRUE)
  }
}

add_annotations_worksheet <- function(wb, seurat_obj) {
  addWorksheet(wb, "Cell_Annotations")
  
  ann_cols <- grep("cell_type|annotation", colnames(seurat_obj@meta.data), 
                  value = TRUE, ignore.case = TRUE)
  
  if (length(ann_cols) > 0) {
    ann_data <- seurat_obj@meta.data[, ann_cols, drop = FALSE]
    writeData(wb, "Cell_Annotations", ann_data, rowNames = TRUE)
  }
}

add_de_worksheet <- function(wb, de_results) {
  addWorksheet(wb, "Differential_Expression")
  writeData(wb, "Differential_Expression", de_results)
}

add_enrichment_worksheet <- function(wb, enrich_results) {
  if (!is.null(enrich_results$ora)) {
    addWorksheet(wb, "ORA_Results")
    combined_ora <- rbind(
      if (!is.null(enrich_results$ora$up)) enrich_results$ora$up,
      if (!is.null(enrich_results$ora$down)) enrich_results$ora$down
    )
    if (!is.null(combined_ora)) {
      writeData(wb, "ORA_Results", combined_ora)
    }
  }
  
  if (!is.null(enrich_results$gsea)) {
    addWorksheet(wb, "GSEA_Results")
    writeData(wb, "GSEA_Results", enrich_results$gsea)
  }
}

add_embeddings_worksheet <- function(wb, seurat_obj) {
  for (reduction in names(seurat_obj@reductions)) {
    sheet_name <- paste0(toupper(reduction), "_Coordinates")
    addWorksheet(wb, sheet_name)
    
    emb_data <- Embeddings(seurat_obj, reduction)
    writeData(wb, sheet_name, emb_data, rowNames = TRUE)
  }
}

add_summary_worksheet <- function(wb, seurat_obj, values) {
  addWorksheet(wb, "Analysis_Summary")
  
  summary_data <- data.frame(
    Metric = c("Total Cells", "Total Genes", "Assays", "Reductions", 
               "DE Genes Tested", "Analysis Date"),
    Value = c(
      ncol(seurat_obj),
      nrow(seurat_obj), 
      paste(names(seurat_obj@assays), collapse = ", "),
      paste(names(seurat_obj@reductions), collapse = ", "),
      if (!is.null(values$last_de_results)) nrow(values$last_de_results) else "N/A",
      as.character(Sys.Date())
    )
  )
  
  writeData(wb, "Analysis_Summary", summary_data)
}

export_results_as_csv <- function(values, params, output_dir) {
  files_created <- list()
  
  # Export each component as separate CSV
  if ("qc" %in% params$results_components && !is.null(values$active_dataset)) {
    qc_file <- file.path(output_dir, "qc_metrics.csv")
    qc_data <- values$active_dataset@meta.data %>%
      select(any_of(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")))
    write.csv(qc_data, qc_file, row.names = TRUE)
    files_created[[length(files_created) + 1]] <- list(name = "qc_metrics.csv", path = qc_file)
  }
  
  if ("de" %in% params$results_components && !is.null(values$last_de_results)) {
    de_file <- file.path(output_dir, "differential_expression.csv")
    write.csv(values$last_de_results, de_file, row.names = FALSE)
    files_created[[length(files_created) + 1]] <- list(name = "differential_expression.csv", path = de_file)
  }
  
  return(files_created)
}

export_results_as_zip <- function(values, params, output_dir) {
  # Create temporary directory for ZIP contents
  zip_contents_dir <- file.path(output_dir, "zip_contents")
  dir.create(zip_contents_dir, recursive = TRUE)
  
  # Export all components to temporary directory
  csv_files <- export_results_as_csv(values, params, zip_contents_dir)
  
  # Create ZIP file
  zip_filename <- "analysis_results.zip"
  zip_filepath <- file.path(output_dir, zip_filename)
  
  zip::zip(zip_filepath, files = zip_contents_dir, mode = "cherry-pick")
  
  list(list(name = zip_filename, path = zip_filepath))
}
