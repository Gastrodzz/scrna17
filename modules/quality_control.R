# Quality Control Module (improved and robust)
# modules/quality_control.R

suppressPackageStartupMessages({
  library(shiny)
  library(Seurat)
  library(ggplot2)
  library(plotly)
  library(dplyr)
  library(viridis)
})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R

map_ensembl_to_symbol <- function(ids, species) {
  # Try to map Ensembl IDs to gene symbols using installed org packages, if available.
  # Returns a named vector mapping input ids -> SYMBOL (or NA if missing).
  suppressWarnings({
    if (species == "human" && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
      keytype <- "ENSEMBL"
      orgpkg <- "org.Hs.eg.db"
    } else if (species == "mouse" && requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
      keytype <- "ENSEMBL"
      orgpkg <- "org.Mm.eg.db"
    } else {
      return(setNames(rep(NA_character_, length(ids)), ids))
    }
    suppressMessages({
      res <- try(AnnotationDbi::select(get(orgpkg), keys = ids, keytype = keytype, columns = "SYMBOL"), silent = TRUE)
    })
    if (inherits(res, "try-error") || is.null(res) || nrow(res) == 0) {
      return(setNames(rep(NA_character_, length(ids)), ids))
    }
    out <- setNames(res$SYMBOL, res[[keytype]])
    out[ids]
  })
}

# ============================================================================
# QC METRICS UI & SERVER (Enhanced)
# ============================================================================

qcMetricsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(tagList(icon("dna"), " Quality Control Metrics")),
          p("Compute comprehensive QC metrics including mitochondrial gene percentage, ribosomal gene percentage, and complexity scores."),
          
          fluidRow(
            column(4,
              selectInput(ns("species_qc"), "Species",
                choices = c("Human" = "human", "Mouse" = "mouse", "Custom" = "custom"),
                selected = "human")
            ),
            column(4,
              conditionalPanel(
                condition = "input.species_qc == 'custom'",
                ns = ns,
                textInput(ns("custom_mt_pattern"), "MT gene pattern", value = "^MT-")
              )
            ),
            column(4,
              div(style = "margin-top: 25px;",
                actionButton(ns("compute_qc"), tagList(icon("play"), " Compute QC Metrics"), 
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          ),
          
          # Additional QC options
          fluidRow(
            column(6,
              checkboxInput(ns("compute_ribo"), "Compute ribosomal gene %", TRUE)
            ),
            column(6,
              checkboxInput(ns("compute_complexity"), "Compute complexity scores", TRUE)
            )
          )
        )
      )
    ),
    
    # Results display - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("qc_results_ui"))
      )
    )
  )
}

qcMetricsServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_qc_results <- reactiveVal(FALSE)
    qc_data <- reactiveVal(NULL)

    # Dynamic UI rendering for results (replaces unreliable conditionalPanel)
    output$qc_results_ui <- renderUI({
      req(isTRUE(show_qc_results()))
      ns <- session$ns
      wellPanel(
        h4("QC Results"),
        # Batch grouping selector
        fluidRow(
          column(4,
            selectInput(ns("qc_group_by"), "Group by (batch/sample)",
              choices = c("None" = "none"), selected = "none")
          ),
          column(8,
            p(style = "margin-top: 30px; color: #888; font-size: 12px;",
              "Select a metadata column to split QC plots by batch, sample, or condition.")
          )
        ),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(6,
              h5("Summary Statistics"),
              downloadButton(ns("dl_qc_summary_table"), "CSV", class = "btn-xs btn-default"),
              div(style = "max-height: 400px; overflow-y: auto;",
                DT::dataTableOutput(ns("qc_summary_table"))
              )
            ),
            column(6,
              h5("Distribution Plots"),
              plotlyOutput(ns("qc_violins"), height = "400px")
            )
          )
        ),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(12,
              h5("QC Correlations"),
              plotlyOutput(ns("qc_correlations"), height = "400px")
            )
          )
        ),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(12,
              h5("Interactive QC Scatter Plots"),
              tabsetPanel(
                tabPanel("UMI vs Genes", plotlyOutput(ns("scatter_umi_genes"), height = "450px")),
                tabPanel("UMI vs MT%", plotlyOutput(ns("scatter_umi_mt"), height = "450px")),
                tabPanel("Genes vs MT%", plotlyOutput(ns("scatter_genes_mt"), height = "450px"))
              )
            )
          )
        ),
        downloadButton(ns("download_qc"), "Download QC Metrics", class = "btn-success")
      )
    })
    
    observeEvent(input$compute_qc, {
      req(values$active_dataset)
      
      withProgress(message = "Computing QC metrics...", value = 0, {
        
        tryCatch({
          obj <- values$active_dataset
          # Ensure assay + layer defaults (Seurat v5 friendly)
          if ("RNA" %in% names(obj@assays)) {
            DefaultAssay(obj) <- "RNA"
            rn <- tryCatch(obj[["RNA"]], error = function(e) NULL)
            if (!is.null(rn) && inherits(rn, "Assay5")) {
              # Join split layers first (avoids "No layers found" warnings)
              lyrs <- tryCatch(Layers(rn), error = function(e) NULL)
              if (!is.null(lyrs) && any(grepl("\\.", lyrs))) {
                obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
                rn <- obj[["RNA"]]
                lyrs <- Layers(rn)
              }
              if (!is.null(lyrs)) {
                if ("counts" %in% lyrs) {
                  DefaultLayer(rn) <- "counts"
                } else if ("data" %in% lyrs) {
                  DefaultLayer(rn) <- "data"
                }
              }
            }
          }
# Ensure we have the right assay
          if ("RNA" %in% names(obj@assays)) {
            DefaultAssay(obj) <- "RNA"
          }
          
          incProgress(0.2, detail = "Computing mitochondrial gene %...")
          
          # Resolve MT pattern from the dropdown, then verify against the data.
          mt_pattern <- switch(input$species_qc,
            "human"  = "^MT-",
            "mouse"  = "^mt-",
            "custom" = input$custom_mt_pattern,
            "^MT-"
          )
          mt_genes <- grep(mt_pattern, rownames(obj), value = TRUE)

          # Self-correct: if the chosen pattern finds nothing, auto-detect from data
          if (length(mt_genes) == 0) {
            det <- detect_species(rownames(obj))
            if (det$species != "unknown" && det$mt_pattern != mt_pattern) {
              mt_pattern <- det$mt_pattern
              mt_genes   <- grep(mt_pattern, rownames(obj), value = TRUE)
              if (length(mt_genes) > 0) {
                showNotification(sprintf(
                  "Species dropdown found no MT genes; auto-detected %s (pattern '%s', %d genes).",
                  det$species, mt_pattern, length(mt_genes)),
                  type = "warning", duration = 8)
              }
            }
          }

          if (length(mt_genes) > 0) {
            obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
            message("Found ", length(mt_genes), " mitochondrial genes")
          } else {
            # assign a named zero vector to metadata
            tmp <- rep(0, ncol(obj)); names(tmp) <- colnames(obj)
            obj[["percent.mt"]] <- tmp
            showNotification(
              "No mitochondrial genes found. If this is expected (pre-filtered matrix), ignore; otherwise check gene naming (symbols vs Ensembl IDs).",
              type = "warning", duration = 10)
          }
          
          # FIXED: Better ribosomal gene detection
          if (input$compute_ribo) {
            incProgress(0.3, detail = "Computing ribosomal gene %...")
            # Ribosomal genes (human: RPS/RPL, mouse: Rps/Rpl)
            ribo_pattern <- "^RP[SL]|^Rp[sl]"
            ribo_genes <- rownames(obj)[grepl(ribo_pattern, rownames(obj), ignore.case = FALSE)]
            if (length(ribo_genes) > 0) {
              obj[["percent.ribo"]] <- PercentageFeatureSet(obj, features = ribo_genes)
              message("Found ", length(ribo_genes), " ribosomal genes")
            } else {
              obj[["percent.ribo"]] <- 0
              showNotification("No ribosomal genes (RPS/RPL) detected in dataset.", type = "warning", duration = 6)
            }
          }
          
          # Cell cycle scoring lives in the Preprocessing tab (needs normalized data)

          # FIXED: Always compute basic complexity scores
          if (input$compute_complexity) {
            incProgress(0.7, detail = "Computing complexity scores...")
            obj <- compute_complexity_scores(obj)
          }
          
          # FIXED: Ensure basic QC metrics exist
          incProgress(0.8, detail = "Computing additional metrics...")
          obj <- compute_additional_qc_metrics(obj)
          
          # FIXED: Update values and set flags
          values$active_dataset <- obj
          values$qc_done <- TRUE  # Set QC flag
          
          # FIXED: Store QC data and show results
          qc_metadata <- obj@meta.data
          qc_data(qc_metadata)
          show_qc_results(TRUE)

          # Populate the batch/group selector with factor-like metadata columns
          md <- obj@meta.data
          factor_cols <- names(md)[sapply(md, function(x) {
            is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 30)
          })]
          # Prioritize common batch columns
          batch_cands <- grep("batch|orig.ident|sample|donor|condition|group|dataset",
                              factor_cols, value = TRUE, ignore.case = TRUE)
          all_choices <- c("None" = "none", setNames(unique(c(batch_cands, factor_cols)),
                                                      unique(c(batch_cands, factor_cols))))
          default_sel <- if ("orig.ident" %in% factor_cols) "orig.ident" else "none"
          updateSelectInput(session, "qc_group_by", choices = all_choices, selected = default_sel)

          incProgress(1, detail = "Complete!")
          
          # Use valid shiny types: "message" (not "success")
          showNotification("QC metrics computed successfully", type = "message")
          
        }, error = function(e) {
          showNotification(paste("QC computation failed:", e$message), type = "error")
          message("QC Error: ", e$message)
        })
      })
    })
    
    # QC Summary Table - FIXED
    output$qc_summary_table <- DT::renderDataTable({
      req(qc_data())

      tryCatch({
        group_col <- input$qc_group_by
        tbl <- if (!is.null(group_col) && group_col != "none" && group_col %in% colnames(qc_data())) {
          create_qc_summary_table_grouped(qc_data(), group_col)
        } else {
          create_qc_summary_table(qc_data())
        }
        if (exists("register_csv_export", mode = "function")) {
          register_csv_export(values, data = tbl, description = "QC Summary", source_module = "Quality Control")
        }
        tbl
      }, error = function(e) {
        message("Summary table error: ", e$message)
        return(DT::datatable(data.frame(Error = "Failed to create summary table")))
      })
    })

    # QC Violin plots - with optional batch grouping
    output$qc_violins <- renderPlotly({
      req(qc_data())
      group_col <- input$qc_group_by

      tryCatch({
        p <- if (!is.null(group_col) && group_col != "none" && group_col %in% colnames(qc_data())) {
          create_qc_violin_plot_grouped(qc_data(), group_col)
        } else {
          create_qc_violin_plot(qc_data())
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "QC Violin Plots",
                         label = "QC Violin Plots", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("Violin plot error: ", e$message)
        plotly_empty() %>% layout(title = "Error creating violin plot")
      })
    })
    
    # Cell cycle pie chart - FIXED
    output$phase_pie <- renderPlotly({
      req(qc_data())
      
      tryCatch({
        p <- create_cell_cycle_pie(qc_data())
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "Cell Phase Distribution",
                         label = "Cell Phase Distribution", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("Phase pie error: ", e$message)
        plotly_empty() %>% layout(title = "Cell cycle data not available")
      })
    })
    
    # QC correlations - FIXED
    output$qc_correlations <- renderPlotly({
      req(qc_data())
      
      tryCatch({
        p <- create_qc_correlation_plot(qc_data())
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "QC Correlations",
                         label = "QC Correlations", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("Correlation plot error: ", e$message)
        plotly_empty() %>% layout(title = "Error creating correlation plot")
      })
    })
    
    # Scatter plots - with optional batch coloring
    output$scatter_umi_genes <- renderPlotly({
      req(qc_data())
      group_col <- input$qc_group_by
      color_by <- if (!is.null(group_col) && group_col != "none" && group_col %in% colnames(qc_data())) group_col else NULL

      tryCatch({
        p <- create_qc_scatter(qc_data(), "nCount_RNA", "nFeature_RNA", "UMI Count", "Gene Count", color_by = color_by)
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "UMI vs Genes Scatter",
                         label = "UMI vs Genes Scatter", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("UMI-Genes scatter error: ", e$message)
        plotly_empty() %>% layout(title = "Error creating scatter plot")
      })
    })

    output$scatter_umi_mt <- renderPlotly({
      req(qc_data())
      group_col <- input$qc_group_by
      color_by <- if (!is.null(group_col) && group_col != "none" && group_col %in% colnames(qc_data())) group_col else NULL

      tryCatch({
        p <- create_qc_scatter(qc_data(), "nCount_RNA", "percent.mt", "UMI Count", "Mitochondrial %", color_by = color_by)
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "UMI vs MT% Scatter",
                         label = "UMI vs MT% Scatter", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("UMI-MT scatter error: ", e$message)
        plotly_empty() %>% layout(title = "Error creating scatter plot")
      })
    })
    
    output$scatter_genes_mt <- renderPlotly({
      req(qc_data())
      group_col <- input$qc_group_by
      color_by <- if (!is.null(group_col) && group_col != "none" && group_col %in% colnames(qc_data())) group_col else NULL

      tryCatch({
        p <- create_qc_scatter(qc_data(), "nFeature_RNA", "percent.mt", "Gene Count", "Mitochondrial %", color_by = color_by)
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "Genes vs MT% Scatter",
                         label = "Genes vs MT% Scatter", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("Genes-MT scatter error: ", e$message)
        plotly_empty() %>% layout(title = "Error creating scatter plot")
      })
    })
    
    # Download handler
    output$download_qc <- downloadHandler(
      filename = function() paste0("qc_metrics_", Sys.Date(), ".csv"),
      content = function(file) {
        if (!is.null(qc_data())) {
          write.csv(qc_data(), file, row.names = TRUE)
        }
      }
    )

    output$dl_qc_summary_table <- downloadHandler(
      filename = function() paste0("qc_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(qc_data())
        tryCatch({
          tbl <- create_qc_summary_table(qc_data())
          write.csv(tbl$x$data, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )

    # AI helper for QC metrics module
    ai_helper_server("qcMetrics", values, "Quality Control", context_fn = function() {
      obj <- values$active_dataset
      if (is.null(obj)) return("")
      md <- obj@meta.data
      lines <- character(0)
      lines <- c(lines, paste0("Cells: ", ncol(obj), ", Genes: ", nrow(obj)))
      for (metric in c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo")) {
        if (metric %in% colnames(md)) {
          vals <- md[[metric]]
          vals <- vals[!is.na(vals)]
          lines <- c(lines, paste0(metric, ": median=", round(median(vals), 2),
                                    ", mean=", round(mean(vals), 2),
                                    ", sd=", round(sd(vals), 2)))
        }
      }
      if ("Phase" %in% colnames(md)) {
        phase_tbl <- table(md$Phase)
        lines <- c(lines, paste0("Cell cycle phases: ",
                                  paste(names(phase_tbl), "=", phase_tbl, collapse = ", ")))
      }
      paste(lines, collapse = "\n")
    })
  })
}

# ============================================================================
# QC FILTERING UI & SERVER (Enhanced)
# ============================================================================

qcFilteringUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(tagList(icon("bullseye"), " Advanced Quality Control Filtering")),
          p("Apply sophisticated filtering based on QC metrics to remove low-quality cells and outliers."),
          
          # Filtering method selection
          fluidRow(
            column(12,
              selectInput(ns("filter_method"), "Filtering Method",
                choices = c(
                  "Manual Thresholds" = "manual",
                  "Median Absolute Deviation (MAD)" = "mad",
                  "Adaptive Outlier Detection" = "adaptive"
                ),
                selected = "manual"),
              helpText(style = "font-size: 11px; color: #888;",
                tags$b("MAD (parametric, symmetric):"), " Assumes approximately symmetric distributions. ",
                "Detects outliers using median \u00b1 N\u00d7MAD, where MAD = median absolute deviation. ",
                "Best for homogeneous samples with unimodal QC distributions. ",
                tags$b("Adaptive (non-parametric, skew-aware):"), " Uses percentile-based outlier detection ",
                "that adapts to the actual shape of each QC metric distribution without assuming symmetry. ",
                "Better for heterogeneous, multi-modal, or skewed datasets where MAD may be too aggressive.")
            )
          ),
          
          # Manual threshold controls
          conditionalPanel(
            condition = "input.filter_method == 'manual'",
            ns = ns,
            fluidRow(
              column(3,
                sliderInput(ns("nfeature_range"), "Gene Count Range",
                  min = 0, max = 10000, value = c(200, 5000), step = 50)
              ),
              column(3,
                sliderInput(ns("ncount_range"), "UMI Count Range",
                  min = 0, max = 100000, value = c(500, 50000), step = 500)
              ),
              column(3,
                sliderInput(ns("mt_max"), "Max Mitochondrial %",
                  min = 0, max = 100, value = 20, step = 1),
                helpText(style = "font-size: 11px; color: #888;",
                  "Tissue-specific guidance: PBMC ~5-10%, Brain ~5%, Kidney ~30%, Heart ~30-40%. ",
                  "Use MAD-based auto-detect for data-driven threshold."),
                actionButton(ns("auto_mt_threshold"), "Auto-detect MT threshold (MAD)",
                  class = "btn-info btn-xs", style = "margin-top: 2px; font-size: 11px;")
              ),
              column(3,
                sliderInput(ns("ribo_max"), "Max Ribosomal %",
                  min = 0, max = 100, value = 50, step = 1)
              )
            )
          ),
          
          # MAD-based filtering
          conditionalPanel(
            condition = "input.filter_method == 'mad'", 
            ns = ns,
            fluidRow(
              column(4,
                numericInput(ns("mad_threshold"), "MAD Threshold", value = 3, min = 1, max = 5, step = 0.5)
              ),
              column(4,
                checkboxInput(ns("log_transform_mad"), "Log transform for MAD", TRUE)
              ),
              column(4,
                checkboxInput(ns("bidirectional_mad"), "Bidirectional outliers", TRUE)
              )
            )
          ),
          
          # Adaptive filtering
          conditionalPanel(
            condition = "input.filter_method == 'adaptive'",
            ns = ns,
            fluidRow(
              column(6,
                sliderInput(ns("outlier_percentile"), "Outlier Percentile",
                  min = 0.01, max = 0.1, value = 0.05, step = 0.01)
              ),
              column(6,
                checkboxInput(ns("use_robust_stats"), "Use robust statistics", TRUE)
              )
            )
          ),
          
          # Additional filtering options
          fluidRow(
            column(4,
              checkboxInput(ns("filter_doublets"), "Filter likely doublets", FALSE)
            ),
            column(4,
              conditionalPanel(
                condition = "input.filter_doublets",
                ns = ns,
                numericInput(ns("doublet_threshold"), "Doublet score threshold",
                           value = 0.25, min = 0, max = 1, step = 0.05),
                actionButton(ns("auto_doublet_threshold"), "Auto-detect from distribution",
                  class = "btn-info btn-xs", style = "margin-top: 2px; font-size: 11px;")
              )
            ),
            column(4,
              checkboxInput(ns("filter_low_complexity"), "Filter low complexity cells", FALSE)
            )
          ),
          
          # AI-assisted threshold suggestion
          fluidRow(
            column(6,
              textInput(ns("tissue_context"), "Tissue Context for AI",
                        placeholder = "e.g., PBMC, brain cortex")
            ),
            column(6,
              div(style = "margin-top: 25px;",
                actionButton(ns("ai_qc_thresholds"), "AI Suggest Thresholds",
                             icon = icon("robot"), class = "ai-btn btn-sm")
              )
            )
          ),
          # AI thread output
          fluidRow(
            column(12, uiOutput(ns("ai_thread_qc_thresholds")))
          ),

          # Preview and apply
          fluidRow(
            column(6,
              actionButton(ns("preview_filter"), tagList(icon("search"), " Preview Filtering"), 
                         class = "btn-info", width = "100%")
            ),
            column(6,
              actionButton(ns("apply_filter"), tagList(icon("cut"), " Apply Filters"), 
                         class = "btn-primary", width = "100%")
            )
          )
        )
      )
    ),
    
    # Filtering results - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("filter_results_ui"))
      )
    )
  )
}

qcFilteringServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_filter_results <- reactiveVal(FALSE)
    filter_preview <- reactiveVal(NULL)
    original_data <- reactiveVal(NULL)

    # --- AI Suggest Thresholds ---
    ai_qc_params <- reactiveVal(NULL)

    observeEvent(input$ai_qc_thresholds, {
      if (is.null(values$ai_api_key) || !nzchar(trimws(values$ai_api_key))) {
        showNotification("Please set your Claude API key in the Data Management tab first.",
                         type = "warning", duration = 8)
        return()
      }
      if (is.null(values$active_dataset)) {
        showNotification("No dataset loaded.", type = "warning")
        return()
      }

      # Show loading state
      output$ai_thread_qc_thresholds <- renderUI({
        ai_thread_loading("AI is analyzing QC metrics...")
      })

      tryCatch({
        obj <- values$active_dataset
        md  <- obj@meta.data

        # Compute QC summary statistics
        qc_stats <- list()
        for (metric in c("nCount_RNA", "nFeature_RNA", "percent.mt")) {
          if (metric %in% colnames(md)) {
            vals <- md[[metric]]
            vals <- vals[!is.na(vals)]
            qc_stats[[metric]] <- list(
              median = round(median(vals), 2),
              mean   = round(mean(vals), 2),
              sd     = round(sd(vals), 2),
              min    = round(min(vals), 2),
              max    = round(max(vals), 2)
            )
          }
        }
        qc_summary_text <- paste(sapply(names(qc_stats), function(m) {
          s <- qc_stats[[m]]
          paste0(m, ": median=", s$median, ", mean=", s$mean, ", sd=", s$sd,
                 ", min=", s$min, ", max=", s$max)
        }), collapse = "\n")

        tissue <- input$tissue_context
        if (is.null(tissue) || !nzchar(trimws(tissue))) {
          tissue <- values$tissue_global
        }
        tissue_str <- if (!is.null(tissue) && nzchar(trimws(tissue))) tissue else "unknown"

        prompt <- build_expert_prompt(
          task_prompt = paste0(
            "TASK: Suggest optimal QC filtering thresholds for this scRNA-seq dataset.\n\n",
            "TISSUE: ", tissue_str, "\n",
            "NUMBER OF CELLS: ", ncol(obj), "\n",
            "QC METRICS SUMMARY:\n", qc_summary_text, "\n\n",
            "Based on these distributions and the tissue type, suggest appropriate filtering thresholds.\n",
            "Return your answer as a JSON object with these keys:\n",
            "  min_genes (integer), max_genes (integer), max_mt_percent (number), min_counts (integer)\n",
            "Wrap the JSON in ```json ... ``` code fences.\n",
            "After the JSON, provide a brief explanation of your reasoning."
          ),
          values = values,
          include_project = TRUE,
          include_dataset = FALSE
        )

        withProgress(message = "AI suggesting thresholds...", value = 0.5, {
          result <- call_ai("qc_thresholds", prompt, values)
        })

        if (!result$ok) {
          output$ai_thread_qc_thresholds <- renderUI({ ai_thread_error(result$error) })
          return()
        }

        ai_text <- result$text

        # Try to parse JSON from response
        json_match <- regmatches(ai_text, regexpr("```json\\s*\\n?(\\{.*?\\})\\s*```",
                                                    ai_text, perl = TRUE))
        if (length(json_match) == 0) {
          json_match <- regmatches(ai_text, regexpr("\\{[^{}]*\"min_genes\"[^{}]*\\}",
                                                      ai_text, perl = TRUE))
        }

        if (length(json_match) > 0) {
          clean_json <- gsub("```json|```", "", json_match[1])
          parsed <- tryCatch(jsonlite::fromJSON(clean_json), error = function(e) NULL)
          if (!is.null(parsed)) {
            ai_qc_params(parsed)
          }
        }

        # Show AI reasoning as thread message with Apply button
        apply_btn <- if (!is.null(ai_qc_params())) {
          actionButton(session$ns("apply_ai_qc"), "Apply Thresholds",
                       class = "ai-thread-apply-btn", icon = icon("check"))
        } else NULL

        output$ai_thread_qc_thresholds <- renderUI({
          ai_thread_message(
            content_html = .mdtohtml(ai_text),
            title = "AI QC Threshold Suggestions",
            model = get_model_for_task("qc_thresholds", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out,
            apply_btn = apply_btn
          )
        })

      }, error = function(e) {
        output$ai_thread_qc_thresholds <- renderUI({ ai_thread_error(e$message) })
      })
    })

    # Apply AI-suggested QC thresholds
    observeEvent(input$apply_ai_qc, {
      params <- ai_qc_params()
      req(!is.null(params))
      if (!is.null(params$min_genes) && !is.null(params$max_genes)) {
        tryCatch(updateSliderInput(session, "nfeature_range",
          value = c(as.integer(params$min_genes), as.integer(params$max_genes))),
          error = function(e) NULL)
      }
      if (!is.null(params$max_mt_percent)) {
        tryCatch(updateSliderInput(session, "mt_max",
          value = as.numeric(params$max_mt_percent)),
          error = function(e) NULL)
      }
      if (!is.null(params$min_counts)) {
        tryCatch(updateSliderInput(session, "ncount_range",
          value = c(as.integer(params$min_counts), input$ncount_range[2])),
          error = function(e) NULL)
      }
      showNotification("AI thresholds applied to sliders.", type = "message", duration = 4)
    })

    # Auto-detect doublet threshold from score distribution inflection point
    observeEvent(input$auto_doublet_threshold, {
      req(values$active_dataset)
      obj <- values$active_dataset
      if (!"doublet_score" %in% colnames(obj@meta.data)) {
        showNotification("Run doublet detection first (no doublet_score found).", type = "warning")
        return()
      }
      scores <- obj@meta.data$doublet_score
      scores <- scores[!is.na(scores)]
      if (length(scores) < 50) {
        showNotification("Too few cells for reliable doublet threshold estimation.", type = "warning")
        return()
      }
      # Compute inflection point from density of doublet scores
      # Use kernel density estimation and find the local minimum between the two modes
      tryCatch({
        d <- density(scores, n = 512, from = 0, to = 1)
        # Find local minima in density
        dy <- diff(d$y)
        sign_changes <- which(dy[-length(dy)] < 0 & dy[-1] >= 0)
        if (length(sign_changes) > 0) {
          # Pick the first local minimum after the main singlet peak
          threshold <- d$x[sign_changes[1] + 1]
        } else {
          # Fallback: use median + 2*MAD of scores
          threshold <- median(scores) + 2 * mad(scores)
        }
        threshold <- round(min(max(threshold, 0.05), 0.95), 3)
        updateNumericInput(session, "doublet_threshold", value = threshold)
        showNotification(
          paste0("Auto-detected doublet threshold: ", threshold, " (from score distribution inflection)"),
          type = "message", duration = 6
        )
      }, error = function(e) {
        showNotification(paste("Doublet threshold auto-detection failed:", e$message), type = "warning")
      })
    })

    # Auto-detect MT threshold using MAD
    observeEvent(input$auto_mt_threshold, {
      req(values$active_dataset)
      obj <- values$active_dataset
      if (!"percent.mt" %in% colnames(obj@meta.data)) {
        showNotification("Compute QC metrics first (percent.mt not found).", type = "warning")
        return()
      }
      mt_vals <- obj@meta.data$percent.mt
      mt_vals <- mt_vals[!is.na(mt_vals)]
      med <- median(mt_vals)
      mad_val <- mad(mt_vals, constant = 1.4826)
      # Use median + 3 * MAD as adaptive threshold, clamped to [1, 100]
      adaptive_threshold <- min(100, max(1, round(med + 3 * mad_val, 1)))
      updateSliderInput(session, "mt_max", value = adaptive_threshold)
      showNotification(
        paste0("MAD-based MT threshold: ", adaptive_threshold,
               "% (median=", round(med, 1), "%, MAD=", round(mad_val, 1), "%)"),
        type = "message", duration = 6
      )
    })

    # Dynamic UI rendering for filter results
    output$filter_results_ui <- renderUI({
      req(isTRUE(show_filter_results()))
      ns <- session$ns
      wellPanel(
        h4("Filtering Results"),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(6,
              h5("Before vs After"),
              downloadButton(ns("dl_filter_summary"), "CSV", class = "btn-xs btn-default"),
              DT::dataTableOutput(ns("filter_summary"))
            ),
            column(6,
              h5("Filter Effects"),
              plotlyOutput(ns("filter_effects_plot"), height = "350px")
            )
          )
        ),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(12,
              h5("QC Distributions: Before vs After Filtering"),
              plotlyOutput(ns("before_after_plot"), height = "450px")
            )
          )
        )
      )
    })
    
    # Preview filtering
    observeEvent(input$preview_filter, {
      req(values$active_dataset)
      
      withProgress(message = "Previewing filters...", value = 0, {
        obj <- values$active_dataset
        original_data(obj@meta.data)
        
        # Apply selected filtering method
        filter_results <- apply_qc_filters(obj, input, preview_only = TRUE)
        filter_preview(filter_results)
        show_filter_results(TRUE)
        
        incProgress(1)
      })
      
      showNotification("Filter preview generated", type = "message")
    })
    
    # Apply filtering
    observeEvent(input$apply_filter, {
      req(values$active_dataset, filter_preview())
      tryCatch({
        withProgress(message = "Applying filters...", value = 0, {
          obj <- values$active_dataset

          cells_to_keep <- filter_preview()$cells_to_keep
          filtered_obj <- subset(obj, cells = cells_to_keep)

          values$active_dataset <- filtered_obj
          incProgress(1)
        })

        n_original <- ncol(values$active_dataset) + length(filter_preview()$cells_to_remove)
        n_filtered <- ncol(values$active_dataset)
        n_removed <- n_original - n_filtered

        showNotification(
          paste("Filtering complete!", n_removed, "cells removed,", n_filtered, "cells retained"),
          type = "message"
        )
      }, error = function(e) {
        showNotification(paste("Filtering error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Filter summary table
    output$filter_summary <- DT::renderDataTable({
      req(filter_preview(), original_data())
      tbl <- create_filter_summary_table(filter_preview(), original_data())
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Filter Summary", source_module = "Quality Control")
      }
      tbl
    }, options = list(dom = 't'))
    
    # Filter effects plot
    output$filter_effects_plot <- renderPlotly({
      req(filter_preview())
      p <- create_filter_effects_plot(filter_preview())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Filter Effects",
                       label = "Filter Effects", source_module = "Quality Control")
      }
      p
    })
    
    # Before/after comparison plot
    output$before_after_plot <- renderPlotly({
      req(filter_preview(), original_data())
      p <- create_before_after_plot(filter_preview(), original_data())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Before/After Filtering",
                       label = "Before/After Filtering", source_module = "Quality Control")
      }
      p
    })

    output$dl_filter_summary <- downloadHandler(
      filename = function() paste0("filter_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(filter_preview(), original_data())
        tryCatch({
          tbl <- create_filter_summary_table(filter_preview(), original_data())
          write.csv(tbl$x$data, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )
  })
}

# ============================================================================
# DOUBLET DETECTION UI & SERVER (Enhanced)
# ============================================================================

doubletDetectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(tagList(icon("users"), " Doublet Detection")),
          p("Identify potential doublets using multiple complementary methods."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("quality_control"),
          
          # Method selection
          fluidRow(
            column(4,
              selectInput(ns("doublet_method"), "Detection Method",
                choices = c(
                  "DoubletFinder" = "doubletfinder",
                  "scDblFinder" = "scdblfinder",
                  "Consensus (DoubletFinder + scDblFinder)" = "consensus"
                ),
                selected = "doubletfinder")
            ),
            column(4,
              numericInput(ns("expected_doublet_rate"), "Expected Doublet Rate (%)",
                value = 8, min = 1, max = 20, step = 0.5)
            ),
            column(4,
              checkboxInput(ns("auto_estimate_rate"), "Auto-estimate rate from cell count", TRUE)
            )
          ),
          
          # Method-specific parameters
          conditionalPanel(
            condition = "input.doublet_method == 'doubletfinder'",
            ns = ns,
            fluidRow(
              column(4,
                numericInput(ns("pk_value"), "pK value", value = 0.09, min = 0.001, max = 1, step = 0.001)
              ),
              column(4,
                checkboxInput(ns("find_optimal_pk"), "Find optimal pK", TRUE)
              ),
              column(4,
                numericInput(ns("pn_value"), "pN value", value = 0.25, min = 0.05, max = 0.5, step = 0.05)
              )
            )
          ),
          
          conditionalPanel(
            condition = "input.doublet_method == 'scdblfinder'",
            ns = ns,
            fluidRow(
              column(6,
                checkboxInput(ns("use_clusters"), "Use existing clusters", TRUE)
              ),
              column(6,
                numericInput(ns("dbr_sd"), "Doublet rate SD", value = 0.015, min = 0.005, max = 0.05, step = 0.005)
              )
            )
          ),
          
          # Execution
          fluidRow(
            column(6,
              actionButton(ns("run_doublet_detection"), tagList(icon("play"), " Run Doublet Detection"),
                         class = "btn-primary btn-lg", width = "100%")
            ),
            column(6,
              actionButton(ns("filter_doublets"), tagList(icon("trash"), " Remove Detected Doublets"),
                         class = "btn-danger", width = "100%")
            )
          )
        )
      )
    ),
    
    # Results - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("doublet_results_ui"))
      )
    )
  )
}

doubletDetectionServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_doublet_results <- reactiveVal(FALSE)
    doublet_results <- reactiveVal(NULL)

    # Dynamic UI rendering for doublet results
    output$doublet_results_ui <- renderUI({
      req(isTRUE(show_doublet_results()))
      ns <- session$ns
      wellPanel(
        h4("Doublet Detection Results"),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(6,
              h5("Detection Summary"),
              downloadButton(ns("dl_doublet_summary"), "CSV", class = "btn-xs btn-default"),
              DT::dataTableOutput(ns("doublet_summary"))
            ),
            column(6,
              h5("Score Distribution"),
              plotlyOutput(ns("doublet_scores"), height = "350px")
            )
          )
        ),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(12,
              h5("Doublets on Embedding"),
              plotlyOutput(ns("doublet_embedding"), height = "500px")
            )
          )
        )
      )
    })
    
    observeEvent(input$run_doublet_detection, {
      req(values$active_dataset)
      job_id <- begin_qc_job(values, "run_doublet_detection", meta = list(
        method = input$doublet_method,
        expected_rate = input$expected_doublet_rate
      ))

      if (exists("validate_seurat_contract", mode = "function")) {
        contract <- validate_seurat_contract(values$active_dataset, stage = "generic", min_cells = 20, min_features = 50)
        if (!isTRUE(contract$ok)) {
          msg <- if (exists("format_contract_feedback", mode = "function")) {
            format_contract_feedback(contract)
          } else {
            paste(contract$errors, collapse = " | ")
          }
          showNotification(msg, type = "error", duration = 12)
          end_qc_job(values, job_id, status = "failed", message = "contract validation failed")
          return()
        }
      }

      tryCatch({
        withProgress(message = "Running doublet detection...", value = 0, {
          obj <- values$active_dataset

          if (input$auto_estimate_rate) {
            estimated_rate <- estimate_doublet_rate(ncol(obj))
            updateNumericInput(session, "expected_doublet_rate", value = estimated_rate)
          }

          result <- run_doublet_detection_method(obj, input)

          if (!is.null(result)) {
            values$active_dataset <- result$object
            doublet_results(result)
            show_doublet_results(TRUE)
            checkpoint_qc(
              values,
              action = "doublet_detection",
              obj = result$object,
              params = list(
                method = input$doublet_method,
                expected_doublet_rate = input$expected_doublet_rate
              ),
              method = result$method
            )
            end_qc_job(values, job_id, status = "done", message = "doublet detection completed")

            n_doublets <- sum(result$object$doublet_classification == "Doublet", na.rm = TRUE)
            showNotification(
              paste("Doublet detection complete!", n_doublets, "potential doublets identified"),
              type = "message"
            )
          } else {
            end_qc_job(values, job_id, status = "failed", message = "doublet detection returned NULL")
            showNotification("Doublet detection failed", type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        end_qc_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("Doublet detection error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Remove doublets
    observeEvent(input$filter_doublets, {
      req(values$active_dataset, doublet_results())
      job_id <- begin_qc_job(values, "filter_doublets", meta = list())
      tryCatch({
        obj <- values$active_dataset
        if ("doublet_classification" %in% colnames(obj@meta.data)) {
          singlets <- obj$doublet_classification == "Singlet"
          filtered_obj <- subset(obj, cells = colnames(obj)[singlets])

          n_removed <- ncol(obj) - ncol(filtered_obj)
          values$active_dataset <- filtered_obj
          checkpoint_qc(
            values,
            action = "filter_doublets",
            obj = filtered_obj,
            params = list(removed_cells = n_removed),
            method = "Metadata-based singlet filtering"
          )
          end_qc_job(values, job_id, status = "done", message = paste("removed", n_removed, "doublets"))

          showNotification(paste("Removed", n_removed, "doublets"), type = "warning")
        }
      }, error = function(e) {
        end_qc_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("Doublet filtering error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$doublet_summary <- DT::renderDataTable({
      req(doublet_results())
      tbl <- create_doublet_summary_table(doublet_results())
      if (exists("register_csv_export", mode = "function")) {
        export_df <- tryCatch({
          if (is.list(tbl) && !is.null(tbl$x$data)) as.data.frame(tbl$x$data) else data.frame()
        }, error = function(e) data.frame())
        register_csv_export(values, data = export_df, description = "Doublet Summary", source_module = "Quality Control")
      }
      tbl
    }, options = list(dom = 't'))
    
    output$doublet_scores <- renderPlotly({
      req(doublet_results())
      p <- create_doublet_score_plot(doublet_results())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Doublet Scores",
                       label = "Doublet Scores", source_module = "Quality Control")
      }
      p
    })
    
    output$doublet_embedding <- renderPlotly({
      req(doublet_results(), values$active_dataset)
      p <- create_doublet_embedding_plot(values$active_dataset)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Doublet Embedding",
                       label = "Doublet Embedding", source_module = "Quality Control")
      }
      p
    })

    output$dl_doublet_summary <- downloadHandler(
      filename = function() paste0("doublet_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(doublet_results())
        tryCatch({
          tbl <- create_doublet_summary_table(doublet_results())
          write.csv(tbl$x$data, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )
  })
}

# ============================================================================
# AMBIENT RNA & CELL CYCLE (Simplified but robust)
# ============================================================================

ambientRNAUI <- function(id) {
  ns <- NS(id)
  tagList(
    wellPanel(
      h4("Ambient RNA Correction"),
      p("Remove ambient RNA contamination using SoupX or decontX."),
      if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("quality_control"),

      fluidRow(
        column(4,
          selectInput(ns("ambient_method"), "Method",
            choices = c("SoupX (Young et al. 2020)" = "soupx",
                        "decontX (Yang et al. 2020)" = "decontx"),
            selected = "soupx")
        ),
        column(4,
          conditionalPanel(
            condition = "input.ambient_method == 'soupx'", ns = ns,
            div(class = "alert alert-warning",
              icon("exclamation-triangle"),
              strong(" SoupX requires BOTH the filtered AND raw (unfiltered) 10X output directories."),
              br(),
              "The raw matrix contains empty droplets used to estimate the ambient RNA profile. ",
              "Upload the raw barcode matrix below."
            ),
            fileInput(ns("raw_matrix_soupx"), "Raw (Unfiltered) Matrix",
                      multiple = TRUE, accept = c(".h5", ".gz")),
            fileInput(ns("raw_matrix_path"), "Raw unfiltered matrix (optional)",
              accept = c(".h5", ".gz", ".mtx"), multiple = FALSE,
              placeholder = "e.g., raw_feature_bc_matrix/"),
            helpText(style = "font-size: 11px; color: #999;",
              "Upload the raw_feature_bc_matrix directory or .h5 file. ",
              "If not provided, SoupX will be disabled and decontX used instead.")
          ),
          conditionalPanel(
            condition = "input.ambient_method == 'decontx'", ns = ns,
            p(style = "margin-top: 25px; font-size: 12px; color: #666;",
              "decontX uses a Bayesian model. Works with any filtered matrix.")
          )
        ),
        column(4,
          div(style = "margin-top: 25px;",
            actionButton(ns("run_soupx"), "Run Ambient RNA Correction", class = "btn-primary", width = "100%")
          )
        )
      ),

      verbatimTextOutput(ns("soupx_status")),
      br(),
      uiOutput(ns("ambient_results_ui"))
    )
  )
}

ambientRNAServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ambient_results <- reactiveVal(NULL)

    # Handle raw unfiltered matrix upload for SoupX
    observeEvent(input$raw_matrix_path, {
      req(input$raw_matrix_path)
      tryCatch({
        fpath <- input$raw_matrix_path$datapath
        if (grepl("\\.h5$", fpath, ignore.case = TRUE)) {
          values$raw_counts_matrix <- Seurat::Read10X_h5(fpath)
        } else {
          # Attempt to read as 10X directory or MatrixMarket
          values$raw_counts_matrix <- Seurat::Read10X(dirname(fpath))
        }
        showNotification("Raw unfiltered matrix loaded for SoupX.", type = "message")
      }, error = function(e) {
        values$raw_counts_matrix <- NULL
        showNotification(paste("Failed to load raw matrix:", e$message), type = "warning")
      })
    })

    output$ambient_results_ui <- renderUI({
      req(ambient_results())
      ns <- session$ns
      wellPanel(
        h5("Ambient RNA Diagnostics (Pre/Post)"),
        fluidRow(
          column(6, DT::dataTableOutput(ns("ambient_diagnostics_table"))),
          column(6, plotlyOutput(ns("ambient_diagnostics_plot"), height = "320px"))
        )
      )
    })

    output$ambient_diagnostics_table <- DT::renderDataTable({
      req(ambient_results())
      DT::datatable(ambient_results()$summary, options = list(dom = "t"), rownames = FALSE)
    })

    output$ambient_diagnostics_plot <- renderPlotly({
      req(ambient_results())
      df <- ambient_results()$per_cell
      n <- nrow(df)
      set.seed(42)
      idx <- if (n > 10000) sample.int(n, 10000) else seq_len(n)
      plot_ly(
        data = df[idx, , drop = FALSE],
        x = ~Pre_UMI,
        y = ~Post_UMI,
        type = "scattergl",
        mode = "markers",
        marker = list(size = 4, opacity = 0.4, color = "#2c7fb8"),
        hovertemplate = "Pre UMI: %{x}<br>Post UMI: %{y}<extra></extra>"
      ) %>%
        layout(
          title = "Per-cell UMI counts before vs after SoupX",
          xaxis = list(title = "Pre-correction UMI"),
          yaxis = list(title = "Post-correction UMI")
        )
    })

    output$soupx_status <- renderText({
      method <- input$ambient_method %||% "soupx"
      if (method == "decontx") {
        if (!requireNamespace("celda", quietly = TRUE)) {
          "celda package (decontX) not installed. Run: BiocManager::install('celda')"
        } else {
          "decontX available. Click to run ambient RNA correction."
        }
      } else if (!requireNamespace("SoupX", quietly = TRUE)) {
        "SoupX package not installed. Run: install.packages('SoupX')"
      } else {
        "SoupX available. Click to run ambient RNA correction."
      }
    })
    
    observeEvent(input$run_soupx, {
      req(values$active_dataset)
      ambient_method <- input$ambient_method %||% "soupx"

      # ---- decontX path (Yang et al., Genome Biology, 2020) ----
      if (ambient_method == "decontx") {
        if (!requireNamespace("celda", quietly = TRUE)) {
          showNotification("celda package (decontX) required. Install: BiocManager::install('celda')", type = "error")
          return()
        }
        job_id <- begin_qc_job(values, "ambient_rna_correction", meta = list(method = "decontX"))
        withProgress(message = "Running decontX...", {
          tryCatch({
            obj <- values$active_dataset
            counts_mat <- safe_get_assay_data(obj, "counts")
            if (is.null(counts_mat) || nrow(counts_mat) == 0) stop("No raw counts found.")
            pre_umi <- Matrix::colSums(counts_mat)
            pre_features <- Matrix::colSums(counts_mat > 0)

            incProgress(0.2, detail = "Converting to SCE...")
            sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_mat))
            cluster_col <- grep("seurat_clusters|cluster|Cluster", colnames(obj@meta.data), value = TRUE)[1]
            if (!is.na(cluster_col)) {
              SummarizedExperiment::colData(sce)$clusters <- as.character(obj@meta.data[[cluster_col]])
            }

            incProgress(0.5, detail = "Running decontX model...")
            sce <- celda::decontX(sce)

            incProgress(0.8, detail = "Extracting corrected counts...")
            corrected_counts <- SummarizedExperiment::assay(sce, "decontXcounts")
            contamination <- SummarizedExperiment::colData(sce)$decontX_contamination
            post_umi <- Matrix::colSums(corrected_counts)
            post_features <- Matrix::colSums(corrected_counts > 0)
            removed_frac <- ifelse(pre_umi > 0, (pre_umi - post_umi) / pre_umi, 0)

            corrected_obj <- CreateSeuratObject(counts = corrected_counts, meta.data = obj@meta.data)
            corrected_obj$decontX_contamination <- as.numeric(contamination)
            for (red_name in names(obj@reductions)) {
              corrected_obj@reductions[[red_name]] <- obj@reductions[[red_name]]
            }
            values$active_dataset <- corrected_obj

            summary_df <- data.frame(
              Metric = c("Total UMI (pre)", "Total UMI (post)", "Median UMI/cell (pre)",
                         "Median UMI/cell (post)", "Median genes/cell (pre)", "Median genes/cell (post)",
                         "Mean removed fraction", "Median contamination"),
              Value = c(sum(pre_umi), sum(post_umi), median(pre_umi), median(post_umi),
                        median(pre_features), median(post_features),
                        mean(removed_frac, na.rm = TRUE), median(contamination, na.rm = TRUE)),
              stringsAsFactors = FALSE
            )
            ambient_results(list(summary = summary_df, per_cell = data.frame(
              Pre_UMI = as.numeric(pre_umi), Post_UMI = as.numeric(post_umi),
              Pre_Features = as.numeric(pre_features), Post_Features = as.numeric(post_features),
              Removed_Fraction = as.numeric(removed_frac), stringsAsFactors = FALSE
            )))
            checkpoint_qc(values, action = "ambient_rna_correction", obj = corrected_obj,
                          params = list(method = "decontX"), method = "celda::decontX")
            end_qc_job(values, job_id, status = "done", message = "decontX correction completed")
            showNotification("decontX correction applied.", type = "message")
          }, error = function(e) {
            end_qc_job(values, job_id, status = "error", message = e$message)
            showNotification(paste("decontX failed:", e$message), type = "error", duration = 10)
          })
        })
        return()
      }

      # ---- SoupX path (Young et al., GigaScience, 2020) ----
      if (!requireNamespace("SoupX", quietly = TRUE)) {
        showNotification("SoupX package required", type = "error")
        return()
      }

      # Check if a raw unfiltered matrix is available
      raw_counts_mat <- NULL
      if (!is.null(values$raw_counts_matrix)) {
        raw_counts_mat <- values$raw_counts_matrix
      } else if (!is.null(values$raw_unfiltered_path) && file.exists(values$raw_unfiltered_path)) {
        tryCatch({
          raw_counts_mat <- Seurat::Read10X(values$raw_unfiltered_path)
        }, error = function(e) {
          raw_counts_mat <<- NULL
        })
      }

      if (is.null(raw_counts_mat)) {
        showNotification(
          "SoupX requires raw unfiltered droplet matrix. Use decontX instead.",
          type = "warning", duration = 10
        )
        # Auto-switch to decontX
        updateSelectInput(session, "ambient_method", selected = "decontx")
        showNotification(
          "Switched to decontX which works with filtered matrix only. Click 'Run' again to proceed.",
          type = "message", duration = 6
        )
        return()
      }

      job_id <- begin_qc_job(values, "ambient_rna_correction", meta = list(method = "SoupX"))

      withProgress(message = "Running SoupX...", {
        tryCatch({
          obj <- values$active_dataset

          # SoupX needs raw count matrix
          counts_mat <- safe_get_assay_data(obj, "counts")
          if (is.null(counts_mat) || nrow(counts_mat) == 0) {
            stop("No raw counts found. SoupX requires raw count data.")
          }
          pre_umi <- Matrix::colSums(counts_mat)
          pre_features <- Matrix::colSums(counts_mat > 0)

          incProgress(0.2, detail = "Creating SoupChannel...")

          # Build SoupChannel with proper tod (raw unfiltered) and toc (filtered)
          sc <- SoupX::SoupChannel(tod = raw_counts_mat, toc = counts_mat, metaData = obj@meta.data)

          # Add cluster info if available
          cluster_col <- grep("seurat_clusters|cluster|Cluster", colnames(obj@meta.data), value = TRUE)[1]
          if (!is.na(cluster_col)) {
            sc <- SoupX::setClusters(sc, setNames(as.character(obj@meta.data[[cluster_col]]), colnames(obj)))
          }

          incProgress(0.5, detail = "Estimating contamination...")
          sc <- SoupX::autoEstCont(sc, doPlot = FALSE)

          incProgress(0.8, detail = "Adjusting counts...")
          corrected_counts <- SoupX::adjustCounts(sc)
          post_umi <- Matrix::colSums(corrected_counts)
          post_features <- Matrix::colSums(corrected_counts > 0)
          removed_frac <- ifelse(pre_umi > 0, (pre_umi - post_umi) / pre_umi, 0)

          # Create corrected Seurat object
          corrected_obj <- CreateSeuratObject(counts = corrected_counts, meta.data = obj@meta.data)

          # Preserve reductions
          for (red_name in names(obj@reductions)) {
            corrected_obj@reductions[[red_name]] <- obj@reductions[[red_name]]
          }

          values$active_dataset <- corrected_obj
          rho <- tryCatch(sc$metaData$rho, error = function(e) NA_real_)
          if (all(is.na(rho))) {
            rho <- tryCatch(sc$metaData$contamination, error = function(e) NA_real_)
          }

          summary_df <- data.frame(
            Metric = c(
              "Total UMI (pre)",
              "Total UMI (post)",
              "Median UMI per cell (pre)",
              "Median UMI per cell (post)",
              "Median genes per cell (pre)",
              "Median genes per cell (post)",
              "Mean UMI removed fraction",
              "Median estimated contamination (rho)"
            ),
            Value = c(
              sum(pre_umi),
              sum(post_umi),
              median(pre_umi),
              median(post_umi),
              median(pre_features),
              median(post_features),
              mean(removed_frac, na.rm = TRUE),
              median(rho, na.rm = TRUE)
            ),
            stringsAsFactors = FALSE
          )

          ambient_results(list(
            summary = summary_df,
            per_cell = data.frame(
              Pre_UMI = as.numeric(pre_umi),
              Post_UMI = as.numeric(post_umi),
              Pre_Features = as.numeric(pre_features),
              Post_Features = as.numeric(post_features),
              Removed_Fraction = as.numeric(removed_frac),
              stringsAsFactors = FALSE
            )
          ))

          checkpoint_qc(
            values,
            action = "ambient_rna_correction",
            obj = corrected_obj,
            params = list(method = "SoupX"),
            method = "SoupX::autoEstCont + SoupX::adjustCounts"
          )
          end_qc_job(values, job_id, status = "done", message = "ambient RNA correction completed")

          showNotification("SoupX correction applied. Counts replaced with decontaminated counts.", type = "message")
        }, error = function(e) {
          end_qc_job(values, job_id, status = "error", message = e$message)
          showNotification(paste("SoupX failed:", e$message,
            "\nNote: SoupX works best with 10X data. Ensure clusters exist before running."),
            type = "error", duration = 10)
        })
      })
    })
  })
}

cellCycleUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(tagList(icon("sync"), " Cell Cycle Analysis")),
          p("Visualize cell cycle phase assignments (S / G2M / Phase)."),
          div(class = "alert alert-info", style = "margin-bottom: 10px;",
            icon("circle-info"),
            tags$b(" Scoring lives in the Preprocessing tab."),
            " Cell cycle scoring requires normalized data, so run",
            tags$b(" 'Cell Cycle Scoring'"), " under Preprocessing first; this tab then",
            " displays the results and includes them in the QC report."),
          fluidRow(
            column(4,
              actionButton(ns("run_cc"), tagList(icon("sync"), " Refresh Results"),
                class = "btn-primary", width = "100%")
            )
          )
        )
      )
    ),
    fluidRow(
      column(6,
        h5(tagList(icon("chart-bar"), " Cell Cycle Phase Distribution")),
        plotlyOutput(ns("cc_pie"), height = "400px")
      ),
      column(6,
        h5(tagList(icon("map"), " Cell Cycle on Embedding")),
        plotlyOutput(ns("cc_dimplot"), height = "400px")
      )
    ),
    fluidRow(
      column(12,
        h5(tagList(icon("clipboard"), " Phase Counts")),
        downloadButton(ns("dl_cc_table"), "CSV", class = "btn-xs btn-default"),
        DT::dataTableOutput(ns("cc_table"))
      )
    )
  )
}

cellCycleServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Cell-cycle SCORING now lives in the Preprocessing tab (it requires
    # normalized data). This tab is the visualization/inspection surface: it
    # renders the S/G2M/Phase results computed in Preprocessing. The button here
    # no longer scores raw data — it guides the user to the correct place.
    observeEvent(input$run_cc, {
      req(values$active_dataset)
      md <- values$active_dataset@meta.data
      if (all(c("S.Score", "G2M.Score", "Phase") %in% colnames(md))) {
        showNotification(
          "Cell cycle already scored (in Preprocessing). Phase plots are shown below.",
          type = "message", duration = 6)
      } else {
        showNotification(
          paste("Cell cycle scoring moved to the Preprocessing tab because it requires",
                "normalized data. Run 'Cell Cycle Scoring' there, then return here to",
                "explore the S/G2M scores and phase distribution."),
          type = "warning", duration = 12)
      }
    })

    output$cc_pie <- renderPlotly({
      req(values$active_dataset)
      req("Phase" %in% colnames(values$active_dataset@meta.data))

      phase_counts <- table(values$active_dataset$Phase)
      p <- plot_ly(
        labels = names(phase_counts),
        values = as.numeric(phase_counts),
        type = "pie",
        marker = list(colors = c("#2ecc71", "#e74c3c", "#3498db"))
      ) %>% layout(title = "Cell Cycle Phase Distribution", autosize = TRUE)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Cell Cycle Pie",
                       label = "Cell Cycle Pie", source_module = "Quality Control")
      }
      p
    })

    output$cc_dimplot <- renderPlotly({
      req(values$active_dataset)
      req("Phase" %in% colnames(values$active_dataset@meta.data))
      obj <- values$active_dataset

      reds <- names(obj@reductions)
      red <- if ("umap" %in% reds) "umap" else if ("tsne" %in% reds) "tsne" else if ("pca" %in% reds) "pca" else NULL

      if (is.null(red)) {
        return(plotly_empty() %>% layout(title = "No reduction available"))
      }

      emb <- Embeddings(obj, red)[, 1:2]
      plot_data <- data.frame(
        Dim1 = emb[, 1], Dim2 = emb[, 2],
        Phase = obj$Phase,
        Cell = colnames(obj)
      )

      p <- plot_ly(data = plot_data, x = ~Dim1, y = ~Dim2, color = ~Phase,
              colors = c("G1" = "#2ecc71", "S" = "#e74c3c", "G2M" = "#3498db"),
              type = "scatter", mode = "markers",
              marker = list(size = 3, opacity = 0.8),
              text = ~Cell) %>%
        layout(title = paste(toupper(red), "colored by Cell Cycle Phase"),
               xaxis = list(title = paste0(toupper(red), "_1")),
               yaxis = list(title = paste0(toupper(red), "_2")))
      if (exists("capture_umap_snapshot", mode = "function"))
        capture_umap_snapshot(values, p, "QC - Cell Cycle", coloring = "Phase", method = red)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Cell Cycle DimPlot",
                       label = "Cell Cycle DimPlot", source_module = "Quality Control")
      }
      p
    })

    output$cc_table <- DT::renderDataTable({
      req(values$active_dataset)
      req("Phase" %in% colnames(values$active_dataset@meta.data))

      phase_df <- data.frame(
        Phase = names(table(values$active_dataset$Phase)),
        Count = as.numeric(table(values$active_dataset$Phase))
      )
      phase_df$Percentage <- round(phase_df$Count / sum(phase_df$Count) * 100, 1)
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = phase_df, description = "Cell Cycle Results", source_module = "Quality Control")
      }
      DT::datatable(phase_df, options = list(dom = 't'))
    })

    output$dl_cc_table <- downloadHandler(
      filename = function() paste0("cell_cycle_phases_", Sys.Date(), ".csv"),
      content = function(file) {
        req(values$active_dataset, "Phase" %in% colnames(values$active_dataset@meta.data))
        phase_df <- data.frame(
          Phase = names(table(values$active_dataset$Phase)),
          Count = as.numeric(table(values$active_dataset$Phase))
        )
        phase_df$Percentage <- round(phase_df$Count / sum(phase_df$Count) * 100, 1)
        write.csv(phase_df, file, row.names = FALSE)
      }
    )
  })
}

batchEffectsUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(tagList(icon("chart-bar"), " Batch Effect Diagnostics")),
          p("Visualize and assess batch effects in your dataset."),
          fluidRow(
            column(4,
              selectInput(ns("batch_col"), "Batch Column", choices = NULL)
            ),
            column(4,
              selectInput(ns("batch_reduction"), "Reduction",
                choices = c("pca", "umap", "tsne"), selected = "pca")
            ),
            column(4,
              actionButton(ns("run_batch_diag"), tagList(icon("chart-bar"), " Run Batch Diagnostics"),
                class = "btn-primary btn-lg", width = "100%")
            )
          )
        )
      )
    ),
    fluidRow(
      column(6,
        h5(tagList(icon("map"), " Embedding Colored by Batch")),
        plotlyOutput(ns("batch_dimplot"), height = "450px")
      ),
      column(6,
        h5(tagList(icon("chart-bar"), " Cells per Batch")),
        plotlyOutput(ns("batch_bar"), height = "450px")
      )
    ),
    fluidRow(
      column(12,
        h5(tagList(icon("clipboard"), " Batch Summary Statistics")),
        downloadButton(ns("dl_batch_stats"), "CSV", class = "btn-xs btn-default"),
        DT::dataTableOutput(ns("batch_stats"))
      )
    ),
    fluidRow(
      column(12,
        h5(tagList(icon("music"), " Per-Batch QC Violin Plots")),
        plotlyOutput(ns("batch_violins"), height = "500px")
      )
    )
  )
}

batchEffectsServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    observe({
      req(values$active_dataset)
      md <- values$active_dataset@meta.data
      factor_cols <- colnames(md)[sapply(md, function(x) is.factor(x) || is.character(x))]
      batch_cands <- c(factor_cols, grep("batch|orig.ident|sample|donor", colnames(md), value = TRUE, ignore.case = TRUE))
      batch_cands <- sort(unique(batch_cands))
      updateSelectInput(session, "batch_col", choices = batch_cands,
                        selected = if ("orig.ident" %in% batch_cands) "orig.ident" else batch_cands[1])

      reds <- names(values$active_dataset@reductions)
      updateSelectInput(session, "batch_reduction", choices = reds,
                        selected = if ("umap" %in% reds) "umap" else reds[1])
    })

    observeEvent(input$run_batch_diag, {
      req(values$active_dataset, input$batch_col)
      showNotification("Batch diagnostics computed", type = "message")
    })

    output$batch_dimplot <- renderPlotly({
      req(values$active_dataset, input$batch_col, input$batch_reduction)
      obj <- values$active_dataset
      red <- input$batch_reduction

      if (!red %in% names(obj@reductions)) {
        return(plotly_empty() %>% layout(title = "Reduction not available"))
      }

      emb <- Embeddings(obj, red)[, 1:2]
      plot_data <- data.frame(
        Dim1 = emb[, 1], Dim2 = emb[, 2],
        Batch = as.character(obj@meta.data[[input$batch_col]]),
        Cell = colnames(obj)
      )

      n_batches <- length(unique(plot_data$Batch))
      pal <- if (n_batches <= 12) RColorBrewer::brewer.pal(max(3, n_batches), "Set3") else grDevices::hcl.colors(n_batches, "Dynamic")

      p <- plot_ly(data = plot_data, x = ~Dim1, y = ~Dim2, color = ~Batch,
              colors = pal, type = "scatter", mode = "markers",
              marker = list(size = 3, opacity = 0.7), text = ~Cell) %>%
        layout(title = paste(toupper(red), "colored by", input$batch_col),
               xaxis = list(title = paste0(toupper(red), "_1")),
               yaxis = list(title = paste0(toupper(red), "_2")))
      if (exists("capture_umap_snapshot", mode = "function"))
        capture_umap_snapshot(values, p, "QC - Batch", coloring = input$batch_col, method = red)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Batch DimPlot",
                       label = "Batch DimPlot", source_module = "Quality Control")
      }
      p
    })

    output$batch_bar <- renderPlotly({
      req(values$active_dataset, input$batch_col)
      batch_counts <- table(values$active_dataset@meta.data[[input$batch_col]])

      p <- plot_ly(x = names(batch_counts), y = as.numeric(batch_counts),
              type = "bar", marker = list(color = "#9b59b6")) %>%
        layout(title = "Cells per Batch",
               xaxis = list(title = input$batch_col),
               yaxis = list(title = "Number of Cells"))
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Quality Control", figure_type = "Batch Distribution",
                       label = "Batch Distribution", source_module = "Quality Control")
      }
      p
    })

    output$batch_stats <- DT::renderDataTable({
      req(values$active_dataset, input$batch_col)
      obj <- values$active_dataset
      batch_vec <- obj@meta.data[[input$batch_col]]

      stats_list <- lapply(unique(batch_vec), function(b) {
        cells <- which(batch_vec == b)
        data.frame(
          Batch = b,
          N_Cells = length(cells),
          Median_UMI = round(median(obj$nCount_RNA[cells], na.rm = TRUE)),
          Median_Genes = round(median(obj$nFeature_RNA[cells], na.rm = TRUE)),
          stringsAsFactors = FALSE
        )
      })

      tbl <- do.call(rbind, stats_list)
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Batch Statistics", source_module = "Quality Control")
      }
      DT::datatable(tbl, options = list(dom = 't', pageLength = 20))
    })

    output$batch_violins <- renderPlotly({
      req(values$active_dataset, input$batch_col)
      tryCatch({
        p <- create_per_batch_qc_violins(values$active_dataset, input$batch_col)
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Quality Control", figure_type = "Batch Violins",
                         label = "Batch Violins", source_module = "Quality Control")
        }
        p
      }, error = function(e) {
        message("Batch violin plot error: ", e$message)
        plotly_empty() %>% layout(title = "Error creating per-batch violin plots")
      })
    })

    output$dl_batch_stats <- downloadHandler(
      filename = function() paste0("batch_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        req(values$active_dataset, input$batch_col)
        obj <- values$active_dataset
        batch_vec <- obj@meta.data[[input$batch_col]]
        stats_list <- lapply(unique(batch_vec), function(b) {
          cells <- which(batch_vec == b)
          data.frame(Batch = b, N_Cells = length(cells),
            Median_UMI = round(median(obj$nCount_RNA[cells], na.rm = TRUE)),
            Median_Genes = round(median(obj$nFeature_RNA[cells], na.rm = TRUE)),
            stringsAsFactors = FALSE)
        })
        write.csv(do.call(rbind, stats_list), file, row.names = FALSE)
      }
    )
  })
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# DEPRECATED: CellCycleScoring should only be run post-normalization (preprocessing.R).
# This helper is retained as a shared utility but should NOT be called during QC.
# Cell cycle scoring in QC would use unnormalized data, producing unreliable results.
safe_cell_cycle_scoring <- function(seurat_obj, species) {
  tryCatch({
    # Seurat v5: join split layers before CellCycleScoring (avoids S4 [[<- error)
    assay <- DefaultAssay(seurat_obj)
    assay_obj <- seurat_obj[[assay]]
    if (inherits(assay_obj, "Assay5")) {
      lyrs <- Layers(assay_obj)
      if (any(grepl("\\.", lyrs))) {
        seurat_obj[[assay]] <- JoinLayers(seurat_obj[[assay]])
      }
      # Ensure normalized data layer exists
      lyrs <- Layers(seurat_obj[[assay]])
      if (!"data" %in% lyrs && "counts" %in% lyrs) {
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      }
    }

    if (species == "human") {
      if (exists("cc.genes.updated.2019", envir = asNamespace("Seurat"))) {
        cc_genes <- get("cc.genes.updated.2019", envir = asNamespace("Seurat"))
        CellCycleScoring(seurat_obj, s.features = cc_genes$s.genes,
                        g2m.features = cc_genes$g2m.genes, set.ident = FALSE)
      } else {
        NULL
      }
    } else if (species == "mouse") {
      if (exists("cc.genes.updated.2019", envir = asNamespace("Seurat"))) {
        cc_genes <- get("cc.genes.updated.2019", envir = asNamespace("Seurat"))
        s_genes_mouse <- paste0(toupper(substr(cc_genes$s.genes, 1, 1)),
                               tolower(substr(cc_genes$s.genes, 2, nchar(cc_genes$s.genes))))
        g2m_genes_mouse <- paste0(toupper(substr(cc_genes$g2m.genes, 1, 1)),
                                 tolower(substr(cc_genes$g2m.genes, 2, nchar(cc_genes$g2m.genes))))

        CellCycleScoring(seurat_obj, s.features = s_genes_mouse,
                        g2m.features = g2m_genes_mouse, set.ident = FALSE)
      } else {
        NULL
      }
    } else {
      NULL
    }
  }, error = function(e) {
    warning("Cell cycle scoring failed: ", e$message)
    NULL
  })
}

compute_complexity_scores <- function(seurat_obj) {
  # Guard against log10(0) = -Inf by clamping minimum values to 1
  safe_nCount <- pmax(seurat_obj$nCount_RNA, 1)
  safe_nFeature <- pmax(seurat_obj$nFeature_RNA, 1)

  # Log10 genes per UMI (complexity)
  seurat_obj$log10GenesPerUMI <- log10(safe_nFeature) / log10(safe_nCount)

  # Replace any remaining NaN/Inf with 0
  seurat_obj$log10GenesPerUMI[is.nan(seurat_obj$log10GenesPerUMI) | is.infinite(seurat_obj$log10GenesPerUMI)] <- 0

  # Novelty score (unique genes / total genes detected)
  seurat_obj$novelty_score <- seurat_obj$nFeature_RNA / safe_nCount

  return(seurat_obj)
}

compute_additional_qc_metrics <- function(seurat_obj) {
  # Ensure required metrics exist
  if (!"nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nCount_RNA <- Matrix::colSums(safe_get_assay_data(seurat_obj, "counts"))
  }
  if (!"nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$nFeature_RNA <- Matrix::colSums(safe_get_assay_data(seurat_obj, "counts") > 0)
  }
  
  return(seurat_obj)
}

create_qc_summary_table <- function(metadata) {
  numeric_cols <- sapply(metadata, is.numeric)
  qc_metrics <- metadata[, numeric_cols, drop = FALSE]
  
  summary_stats <- do.call(rbind, lapply(qc_metrics, function(x) {
    data.frame(
      Median = round(median(x, na.rm = TRUE), 2),
      Mean = round(mean(x, na.rm = TRUE), 2),
      SD = round(sd(x, na.rm = TRUE), 2),
      Min = round(min(x, na.rm = TRUE), 2),
      Max = round(max(x, na.rm = TRUE), 2)
    )
  }))
  
  summary_stats$Metric <- rownames(summary_stats)
  summary_stats <- summary_stats[, c("Metric", "Median", "Mean", "SD", "Min", "Max")]
  
  DT::datatable(summary_stats, options = list(pageLength = 10, dom = 't'))
}

create_qc_summary_table_grouped <- function(metadata, group_col) {
  qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                       colnames(metadata))
  if (length(qc_cols) == 0) return(DT::datatable(data.frame(Error = "No QC metrics found")))

  groups <- as.character(metadata[[group_col]])
  unique_groups <- sort(unique(groups))

  rows <- lapply(unique_groups, function(g) {
    sub <- metadata[groups == g, , drop = FALSE]
    row <- data.frame(Batch = g, N_Cells = nrow(sub), stringsAsFactors = FALSE)
    for (col in qc_cols) {
      vals <- sub[[col]]
      row[[paste0(col, "_median")]] <- round(median(vals, na.rm = TRUE), 2)
      row[[paste0(col, "_mean")]] <- round(mean(vals, na.rm = TRUE), 2)
    }
    row
  })
  summary_df <- do.call(rbind, rows)
  DT::datatable(summary_df, options = list(pageLength = 20, dom = 'tp', scrollX = TRUE),
                rownames = FALSE)
}

create_qc_violin_plot <- function(metadata) {
  # Select key QC metrics
  qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                      colnames(metadata))

  if (length(qc_cols) == 0) return(plotly_empty())

  # Use native plotly violin plots instead of ggplotly to avoid
  # serializing every data point (which crashes with large datasets)
  metric_colors <- c("nCount_RNA" = "#3498db", "nFeature_RNA" = "#2ecc71",
                     "percent.mt" = "#e74c3c", "percent.ribo" = "#9b59b6")
  metric_labels <- c("nCount_RNA" = "UMI Count", "nFeature_RNA" = "Gene Count",
                     "percent.mt" = "MT %", "percent.ribo" = "Ribo %")

  n_cells <- nrow(metadata)
  set.seed(42)
  p <- plot_ly()
  for (col in qc_cols) {
    vals <- metadata[[col]]
    # Downsample for rendering if too many cells
    if (length(vals) > 5000) {
      vals <- sample(vals, 5000)
    }
    p <- p %>% add_trace(
      y = vals,
      x = rep(metric_labels[col], length(vals)),
      type = "violin",
      name = metric_labels[col],
      box = list(visible = TRUE),
      meanline = list(visible = TRUE),
      fillcolor = metric_colors[col],
      line = list(color = metric_colors[col]),
      opacity = 0.7,
      # Add jitter points overlay for small datasets (n < 5000)
      points = if (n_cells < 5000) "all" else FALSE,
      jitter = if (n_cells < 5000) 0.3 else 0,
      pointpos = if (n_cells < 5000) -1.5 else 0,
      marker = if (n_cells < 5000) list(size = 2, opacity = 0.4, color = metric_colors[col]) else list()
    )
  }

  p %>% layout(
    title = "QC Metrics Distribution",
    yaxis = list(title = "Value"),
    xaxis = list(title = ""),
    showlegend = FALSE,
    autosize = TRUE
  )
}

create_qc_violin_plot_grouped <- function(metadata, group_col) {
  qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                      colnames(metadata))
  if (length(qc_cols) == 0 || !group_col %in% colnames(metadata)) return(plotly_empty())

  metric_labels <- c("nCount_RNA" = "UMI Count", "nFeature_RNA" = "Gene Count",
                     "percent.mt" = "MT %", "percent.ribo" = "Ribo %")

  groups <- as.character(metadata[[group_col]])
  unique_groups <- sort(unique(groups))
  n_groups <- length(unique_groups)
  pal <- if (n_groups <= 8) {
    RColorBrewer::brewer.pal(max(3, n_groups), "Set2")[seq_len(n_groups)]
  } else {
    grDevices::hcl.colors(n_groups, "Dynamic")
  }

  # Create subplots - one per QC metric, each with violins per batch
  set.seed(42)
  plots <- lapply(seq_along(qc_cols), function(i) {
    col <- qc_cols[i]
    p <- plot_ly()
    for (gi in seq_along(unique_groups)) {
      grp <- unique_groups[gi]
      vals <- metadata[[col]][groups == grp]
      if (length(vals) > 3000) vals <- sample(vals, 3000)
      p <- p %>% add_trace(
        y = vals, x = rep(grp, length(vals)),
        type = "violin", name = grp,
        legendgroup = grp, showlegend = (i == 1),
        box = list(visible = TRUE), meanline = list(visible = TRUE),
        fillcolor = pal[gi], line = list(color = pal[gi]), opacity = 0.7
      )
    }
    p %>% layout(
      annotations = list(list(
        text = metric_labels[col], xref = "paper", yref = "paper",
        x = 0.5, y = 1.05, showarrow = FALSE, font = list(size = 12)
      )),
      yaxis = list(title = metric_labels[col])
    )
  })

  subplot(plots, nrows = 1, shareY = FALSE, titleX = TRUE, titleY = TRUE) %>%
    layout(
      title = paste("QC Metrics by", group_col),
      showlegend = TRUE, autosize = TRUE
    )
}

create_cell_cycle_pie <- function(metadata) {
  if (!"Phase" %in% colnames(metadata)) {
    return(plotly_empty() %>% layout(title = "Cell cycle data not available"))
  }
  
  phase_counts <- table(metadata$Phase)
  
  plot_ly(
    labels = names(phase_counts),
    values = as.numeric(phase_counts),
    type = "pie",
    hole = 0.3,
    marker = list(colors = viridis::viridis(length(phase_counts)))
  ) %>%
    layout(title = "Cell Cycle Phase Distribution")
}

create_qc_correlation_plot <- function(metadata) {
  # Select only key QC numeric metrics to keep the heatmap readable
  qc_cols <- intersect(
    c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
      "log10GenesPerUMI", "novelty_score", "S.Score", "G2M.Score"),
    colnames(metadata)
  )

  if (length(qc_cols) < 2) return(plotly_empty())

  qc_data <- metadata[, qc_cols, drop = FALSE]

  # Compute correlation matrix
  cor_matrix <- cor(qc_data, use = "complete.obs")

  # Use native plotly heatmap instead of ggplotly
  plot_ly(
    x = colnames(cor_matrix),
    y = rownames(cor_matrix),
    z = cor_matrix,
    type = "heatmap",
    colorscale = list(c(0, "blue"), c(0.5, "white"), c(1, "red")),
    zmin = -1, zmax = 1,
    hovertemplate = "%{y} vs %{x}<br>Correlation: %{z:.3f}<extra></extra>"
  ) %>%
    layout(
      title = "QC Metrics Correlation",
      xaxis = list(title = "", tickangle = 45),
      yaxis = list(title = "")
    )
}

create_qc_scatter <- function(metadata, x_var, y_var, x_label, y_label, color_by = NULL) {
  if (!all(c(x_var, y_var) %in% colnames(metadata))) {
    return(plotly_empty() %>% layout(title = "Required metrics not available"))
  }

  # Downsample to prevent browser crashes with large datasets
  set.seed(42)
  n_cells <- nrow(metadata)
  if (n_cells > 10000) {
    idx <- sample(n_cells, 10000)
    plot_data <- metadata[idx, , drop = FALSE]
    subtitle <- paste0(" (", formatC(n_cells, format = "d", big.mark = ","),
                       " cells, showing 10,000)")
  } else {
    plot_data <- metadata
    subtitle <- ""
  }

  if (!is.null(color_by) && color_by %in% colnames(plot_data)) {
    grps <- as.character(plot_data[[color_by]])
    unique_grps <- sort(unique(grps))
    n_grps <- length(unique_grps)
    pal <- if (n_grps <= 8) RColorBrewer::brewer.pal(max(3, n_grps), "Set2")[seq_len(n_grps)]
           else grDevices::hcl.colors(n_grps, "Dynamic")

    p <- plot_ly()
    for (gi in seq_along(unique_grps)) {
      g <- unique_grps[gi]
      mask <- grps == g
      p <- p %>% add_trace(
        x = plot_data[[x_var]][mask], y = plot_data[[y_var]][mask],
        type = "scattergl", mode = "markers", name = g,
        marker = list(size = 3, opacity = 0.5, color = pal[gi]),
        hovertemplate = paste0(color_by, ": ", g, "<br>", x_label, ": %{x}<br>", y_label, ": %{y}<extra></extra>")
      )
    }
    p %>% layout(
      title = paste0(x_label, " vs ", y_label, " (by ", color_by, ")", subtitle),
      xaxis = list(title = x_label), yaxis = list(title = y_label),
      legend = list(title = list(text = color_by))
    )
  } else {
    plot_ly(
      x = plot_data[[x_var]],
      y = plot_data[[y_var]],
      type = "scattergl",
      mode = "markers",
      marker = list(size = 3, opacity = 0.5),
      hovertemplate = paste0(x_label, ": %{x}<br>", y_label, ": %{y}<extra></extra>")
    ) %>%
      layout(
        title = paste0(x_label, " vs ", y_label, subtitle),
        xaxis = list(title = x_label),
        yaxis = list(title = y_label)
      )
  }
}

apply_qc_filters <- function(seurat_obj, input_params, preview_only = FALSE) {
  metadata <- seurat_obj@meta.data
  
  if (input_params$filter_method == "manual") {
    # Manual threshold filtering
    keep_cells <- (
      metadata$nFeature_RNA >= input_params$nfeature_range[1] &
      metadata$nFeature_RNA <= input_params$nfeature_range[2] &
      metadata$nCount_RNA >= input_params$ncount_range[1] &
      metadata$nCount_RNA <= input_params$ncount_range[2] &
      metadata$percent.mt <= input_params$mt_max
    )
    # Also filter by ribosomal % if the column exists and threshold is provided
    if ("percent.ribo" %in% colnames(metadata) && !is.null(input_params$ribo_max)) {
      keep_cells <- keep_cells & (metadata$percent.ribo <= input_params$ribo_max)
    }
  } else if (input_params$filter_method == "mad") {
    # MAD-based filtering
    keep_cells <- apply_mad_filtering(metadata, input_params)
  } else if (input_params$filter_method == "adaptive") {
    # Adaptive filtering
    keep_cells <- apply_adaptive_filtering(metadata, input_params)
  } else {
    keep_cells <- rep(TRUE, nrow(metadata))
  }
  
  cells_to_keep <- rownames(metadata)[keep_cells]
  cells_to_remove <- rownames(metadata)[!keep_cells]
  
  list(
    cells_to_keep = cells_to_keep,
    cells_to_remove = cells_to_remove,
    n_kept = length(cells_to_keep),
    n_removed = length(cells_to_remove),
    filter_summary = create_filter_summary(metadata, keep_cells)
  )
}

apply_mad_filtering <- function(metadata, params) {
  # MAD-based outlier detection
  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  keep_all <- rep(TRUE, nrow(metadata))
  
  for (metric in metrics) {
    if (metric %in% colnames(metadata)) {
      values <- metadata[[metric]]
      if (params$log_transform_mad && metric %in% c("nFeature_RNA", "nCount_RNA")) {
        values <- log10(values + 1)
      }
      
      median_val <- median(values, na.rm = TRUE)
      mad_val <- mad(values, na.rm = TRUE)
      
      if (params$bidirectional_mad) {
        outliers <- abs(values - median_val) > params$mad_threshold * mad_val
      } else {
        # For MT%, only remove high outliers
        if (metric == "percent.mt") {
          outliers <- values > median_val + params$mad_threshold * mad_val
        } else {
          outliers <- values < median_val - params$mad_threshold * mad_val
        }
      }
      
      keep_all <- keep_all & !outliers
    }
  }
  
  keep_all
}

apply_adaptive_filtering <- function(metadata, params) {
  # Adaptive percentile-based filtering
  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  keep_all <- rep(TRUE, nrow(metadata))
  
  for (metric in metrics) {
    if (metric %in% colnames(metadata)) {
      values <- metadata[[metric]]
      
      if (metric == "percent.mt") {
        # Remove high MT% cells
        threshold <- quantile(values, 1 - params$outlier_percentile, na.rm = TRUE)
        keep_all <- keep_all & values <= threshold
      } else {
        # Remove low count/feature cells
        threshold <- quantile(values, params$outlier_percentile, na.rm = TRUE)
        keep_all <- keep_all & values >= threshold
      }
    }
  }
  
  keep_all
}


create_filter_summary <- function(metadata, keep_cells) {
  n_total <- nrow(metadata)
  n_kept <- sum(keep_cells, na.rm = TRUE)
  n_removed <- n_total - n_kept
  metrics <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), colnames(metadata))
  if (length(metrics) == 0) {
    return(data.frame(Metric=character(0), Kept_Median=numeric(0), Kept_Mean=numeric(0),
                      Removed_Median=numeric(0), Removed_Mean=numeric(0)))
  }
  rows <- lapply(metrics, function(m) {
    vals <- metadata[[m]]
    data.frame(
      Metric = m,
      Kept_Median = round(stats::median(vals[keep_cells], na.rm = TRUE), 2),
      Kept_Mean = round(mean(vals[keep_cells], na.rm = TRUE), 2),
      Removed_Median = round(stats::median(vals[!keep_cells], na.rm = TRUE), 2),
      Removed_Mean = round(mean(vals[!keep_cells], na.rm = TRUE), 2)
    )
  })
  do.call(rbind, rows)
}
create_filter_summary_table <- function(filter_results, original_metadata) {
  summary_df <- data.frame(
    Metric = c("Total Cells", "Cells Kept", "Cells Removed", "Retention Rate"),
    Value = c(
      nrow(original_metadata),
      filter_results$n_kept,
      filter_results$n_removed,
      paste0(round(filter_results$n_kept / nrow(original_metadata) * 100, 1), "%")
    )
  )
  
  DT::datatable(summary_df, options = list(dom = 't'))
}

create_filter_effects_plot <- function(filter_results) {
  effect_data <- data.frame(
    Category = c("Retained", "Removed"),
    Count = c(filter_results$n_kept, filter_results$n_removed),
    Percentage = c(
      filter_results$n_kept / (filter_results$n_kept + filter_results$n_removed) * 100,
      filter_results$n_removed / (filter_results$n_kept + filter_results$n_removed) * 100
    )
  )
  
  plot_ly(
    data = effect_data,
    x = ~Category,
    y = ~Count,
    type = "bar",
    marker = list(color = c("green", "red"))
  ) %>%
    layout(
      title = "Filter Effects",
      xaxis = list(title = ""),
      yaxis = list(title = "Number of Cells")
    )
}

create_before_after_plot <- function(filter_results, original_metadata) {
  # Create comparison of key metrics before and after filtering
  kept_metadata <- original_metadata[filter_results$cells_to_keep, ]

  # Downsample for rendering to prevent browser crashes
  set.seed(42)
  max_points <- 5000
  if (nrow(original_metadata) > max_points) {
    before_idx <- sample(nrow(original_metadata), max_points)
  } else {
    before_idx <- seq_len(nrow(original_metadata))
  }
  if (nrow(kept_metadata) > max_points) {
    after_idx <- sample(nrow(kept_metadata), max_points)
  } else {
    after_idx <- seq_len(nrow(kept_metadata))
  }

  # Use scattergl for WebGL-accelerated rendering
  p <- plot_ly() %>%
    add_trace(
      x = original_metadata$nCount_RNA[before_idx],
      y = original_metadata$nFeature_RNA[before_idx],
      type = "scattergl",
      mode = "markers",
      marker = list(size = 3, opacity = 0.4, color = "red"),
      name = "Before"
    ) %>%
    add_trace(
      x = kept_metadata$nCount_RNA[after_idx],
      y = kept_metadata$nFeature_RNA[after_idx],
      type = "scattergl",
      mode = "markers",
      marker = list(size = 3, opacity = 0.4, color = "blue"),
      name = "After"
    ) %>%
    layout(
      title = "QC Distribution: Before vs After Filtering",
      xaxis = list(title = "UMI Count"),
      yaxis = list(title = "Gene Count")
    )

  p
}

estimate_doublet_rate <- function(n_cells) {
  # Adaptive doublet rate estimation using linear formula

  # Formula: ~0.8% per 1000 cells loaded, capped at 20%
  # Reference: 10X Genomics multiplet rate guidelines
  doublet_rate <- min(0.20, n_cells * 0.0008)

  # Return as percentage (0-100 scale)
  doublet_rate * 100
}

begin_qc_job <- function(values, action, meta = list()) {
  if (!exists("start_job", mode = "function")) return(NULL)
  start_job(values, module = "quality_control", action = action, meta = meta)
}

end_qc_job <- function(values, job_id, status = "done", message = NULL, details = list()) {
  if (is.null(job_id) || !exists("finish_job", mode = "function")) return(invisible(NULL))
  finish_job(values, job_id = job_id, status = status, message = message, details = details)
}

checkpoint_qc <- function(values, action, obj, params = list(), method = NULL) {
  if (!exists("register_dataset_checkpoint", mode = "function")) return(invisible(NULL))
  register_dataset_checkpoint(
    values,
    module = "quality_control",
    action = action,
    obj = obj,
    params = params,
    methods = list(method = method)
  )
}

.normalize01 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) <= 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

run_doublet_detection_method <- function(seurat_obj, input_params) {
  method <- input_params$doublet_method
  
  if (method == "doubletfinder") {
    return(run_doubletfinder(seurat_obj, input_params))
  } else if (method == "scdblfinder") {
    return(run_scdblfinder(seurat_obj, input_params))
  } else if (method == "consensus") {
    return(run_consensus_doublet_detection(seurat_obj, input_params))
  } else {
    showNotification(paste("Unsupported doublet method:", method), type = "error")
    return(NULL)
  }
}

MAX_DOUBLETFINDER_CELLS <- 15000L

run_doubletfinder_safe <- function(seurat_obj, pK, nExp_poi, PCs = 1:20, pN = 0.25) {
  n_cells <- ncol(seurat_obj)
  subsample <- n_cells > MAX_DOUBLETFINDER_CELLS

  if (subsample) {
    set.seed(42)
    sample_idx <- sample.int(n_cells, MAX_DOUBLETFINDER_CELLS)
    sub_obj <- seurat_obj[, sample_idx]
    showNotification(
      paste0("DoubletFinder: subsampling to ", MAX_DOUBLETFINDER_CELLS,
             " cells (dataset has ", n_cells, " cells)."),
      type = "warning", duration = 6
    )
  } else {
    sub_obj <- seurat_obj
    sample_idx <- seq_len(n_cells)
  }

  sub_obj <- DoubletFinder::doubletFinder_v3(sub_obj, PCs = PCs, pN = pN, pK = pK, nExp = nExp_poi)
  df_col <- grep("DF.classifications", colnames(sub_obj@meta.data), value = TRUE)[1]
  doublet_labels <- as.character(sub_obj@meta.data[[df_col]])

  if (subsample) {
    full_labels <- rep("Singlet", n_cells)
    full_labels[sample_idx] <- doublet_labels
    # KNN propagation for unsampled cells
    if (requireNamespace("FNN", quietly = TRUE) && "pca" %in% names(seurat_obj@reductions)) {
      pca_all <- Embeddings(seurat_obj, "pca")[, 1:min(20, ncol(Embeddings(seurat_obj, "pca")))]
      pca_sampled <- pca_all[sample_idx, , drop = FALSE]
      pca_unsampled <- pca_all[-sample_idx, , drop = FALSE]
      if (nrow(pca_unsampled) > 0) {
        knn_res <- FNN::get.knnx(pca_sampled, pca_unsampled, k = 5)
        for (i in seq_len(nrow(pca_unsampled))) {
          nn_labels <- doublet_labels[knn_res$nn.index[i, ]]
          full_labels[which(!(seq_len(n_cells) %in% sample_idx))[i]] <-
            names(sort(table(nn_labels), decreasing = TRUE))[1]
        }
      }
    }
  } else {
    full_labels <- doublet_labels
  }

  seurat_obj$doublet_label <- full_labels
  seurat_obj
}

run_doubletfinder <- function(seurat_obj, params) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    showNotification("DoubletFinder package required", type = "error")
    return(NULL)
  }
  
  tryCatch({
    # Ensure normalized data exists
    assay_name <- DefaultAssay(seurat_obj)
    assay_obj <- seurat_obj[[assay_name]]
    layers <- if (inherits(assay_obj, "Assay5")) {
      Layers(assay_obj)
    } else {
      c("counts", "data", "scale.data")
    }
    if (!"data" %in% layers) {
      seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
      seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    }

    # Ensure PCA exists
    if (!"pca" %in% names(seurat_obj@reductions)) {
      if (length(VariableFeatures(seurat_obj)) == 0) {
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
      }
      seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
    }
    
    # Find optimal pK if requested
    # Handle DoubletFinder API change: v3 functions renamed in newer versions
    has_v3_api <- exists("paramSweep_v3", where = asNamespace("DoubletFinder"), mode = "function")

    if (params$find_optimal_pk) {
      sweep.res <- if (has_v3_api) {
        DoubletFinder::paramSweep_v3(seurat_obj, PCs = 1:20, sct = FALSE)
      } else {
        DoubletFinder::paramSweep(seurat_obj, PCs = 1:20, sct = FALSE)
      }
      sweep.stats <- DoubletFinder::summarizeSweep(sweep.res, GT = FALSE)
      bcmvn <- DoubletFinder::find.pK(sweep.stats)
      optimal_pk <- bcmvn$pK[which.max(bcmvn$BCmetric)]
      pk_to_use <- as.numeric(as.character(optimal_pk))
    } else {
      pk_to_use <- params$pk_value
    }

    # Calculate expected number of doublets
    nExp_poi <- round(params$expected_doublet_rate / 100 * ncol(seurat_obj))

    # Run DoubletFinder with safe subsampling for large datasets
    seurat_obj <- run_doubletfinder_safe(
      seurat_obj, pK = pk_to_use, nExp_poi = nExp_poi,
      PCs = 1:20, pN = params$pn_value
    )
    
    # Extract results
    pann_col <- grep("^pANN", colnames(seurat_obj@meta.data), value = TRUE)[1]
    class_col <- grep("^DF.classifications", colnames(seurat_obj@meta.data), value = TRUE)[1]
    
    seurat_obj$doublet_score <- seurat_obj@meta.data[[pann_col]]
    seurat_obj$doublet_classification <- seurat_obj@meta.data[[class_col]]
    
    list(
      object = seurat_obj,
      method = "DoubletFinder",
      scores = seurat_obj$doublet_score,
      classifications = seurat_obj$doublet_classification
    )
    
  }, error = function(e) {
    showNotification(paste("DoubletFinder failed:", e$message), type = "error")
    NULL
  })
}

run_scdblfinder <- function(seurat_obj, params) {
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    showNotification("scDblFinder package required", type = "error")
    return(NULL)
  }
  
  tryCatch({
    # Convert to SingleCellExperiment
    sce <- as.SingleCellExperiment(seurat_obj)
    
    # Run scDblFinder
    sce <- scDblFinder::scDblFinder(
      sce,
      dbr = params$expected_doublet_rate / 100,
      dbr.sd = params$dbr_sd
    )
    
    # Extract results
    doublet_scores <- sce$scDblFinder.score
    doublet_class <- ifelse(sce$scDblFinder.class == "doublet", "Doublet", "Singlet")
    
    seurat_obj$doublet_score <- doublet_scores
    seurat_obj$doublet_classification <- doublet_class
    
    list(
      object = seurat_obj,
      method = "scDblFinder",
      scores = doublet_scores,
      classifications = doublet_class
    )
    
  }, error = function(e) {
    showNotification(paste("scDblFinder failed:", e$message), type = "error")
    NULL
  })
}

run_consensus_doublet_detection <- function(seurat_obj, params) {
  df_res <- run_doubletfinder(seurat_obj, params)
  sc_res <- run_scdblfinder(seurat_obj, params)

  if (is.null(df_res) && is.null(sc_res)) {
    showNotification("Consensus mode failed: both DoubletFinder and scDblFinder are unavailable.", type = "error", duration = 10)
    return(NULL)
  }

  if (is.null(df_res)) {
    sc_res$method <- "Consensus fallback: scDblFinder"
    return(sc_res)
  }
  if (is.null(sc_res)) {
    df_res$method <- "Consensus fallback: DoubletFinder"
    return(df_res)
  }

  cells <- colnames(seurat_obj)
  df_scores <- setNames(as.numeric(df_res$scores), colnames(df_res$object))
  sc_scores <- setNames(as.numeric(sc_res$scores), colnames(sc_res$object))
  df_class <- setNames(as.character(df_res$classifications), colnames(df_res$object))
  sc_class <- setNames(as.character(sc_res$classifications), colnames(sc_res$object))

  df_scores <- df_scores[cells]
  sc_scores <- sc_scores[cells]
  df_class <- df_class[cells]
  sc_class <- sc_class[cells]

  df_norm <- .normalize01(df_scores)
  sc_norm <- .normalize01(sc_scores)

  # Weighted consensus: scDblFinder is more robust in benchmarks (Xi & Li, Cell Syst 2021)
  # so we give it slightly more weight (0.55 vs 0.45)
  w_df <- 0.45
  w_sc <- 0.55
  combined_score <- w_df * df_norm + w_sc * sc_norm
  combined_score[!is.finite(combined_score)] <- 0

  votes <- as.integer(df_class == "Doublet") + as.integer(sc_class == "Doublet")
  expected_n <- max(1L, round(as.numeric(params$expected_doublet_rate) / 100 * length(cells)))
  ranked_idx <- order(combined_score, decreasing = TRUE)
  top_expected <- rep(FALSE, length(cells))
  top_expected[utils::head(ranked_idx, expected_n)] <- TRUE

  final_class <- ifelse(votes == 2L | (votes == 1L & top_expected), "Doublet", "Singlet")

  # Compute estimated FPR based on score distribution of singlets
  singlet_scores <- combined_score[final_class == "Singlet"]
  doublet_scores <- combined_score[final_class == "Doublet"]
  estimated_fpr <- if (length(singlet_scores) > 0 && length(doublet_scores) > 0) {
    # Estimate FPR as fraction of singlets with score above min doublet score
    min_doublet <- min(doublet_scores, na.rm = TRUE)
    sum(singlet_scores >= min_doublet, na.rm = TRUE) / length(singlet_scores)
  } else {
    NA_real_
  }

  out_obj <- seurat_obj
  out_obj$doublet_score <- combined_score
  out_obj$doublet_classification <- final_class
  out_obj$doublet_score_doubletfinder <- df_scores
  out_obj$doublet_score_scdblfinder <- sc_scores
  out_obj$doublet_vote_count <- votes
  out_obj$doublet_class_doubletfinder <- df_class
  out_obj$doublet_class_scdblfinder <- sc_class

  agreement <- mean((df_class == sc_class), na.rm = TRUE)
  disagreement <- mean((df_class != sc_class), na.rm = TRUE)
  n_doublets <- sum(final_class == "Doublet", na.rm = TRUE)
  pct_doublets <- round(100 * n_doublets / length(cells), 2)

  list(
    object = out_obj,
    method = "Consensus (DoubletFinder + scDblFinder)",
    scores = combined_score,
    classifications = final_class,
    agreement_rate = agreement,
    disagreement_rate = disagreement,
    n_doublets = n_doublets,
    pct_doublets = pct_doublets,
    estimated_fpr = estimated_fpr,
    summary_table = data.frame(
      Method = c("DoubletFinder", "scDblFinder", "Consensus"),
      Doublets = c(sum(df_class == "Doublet", na.rm = TRUE),
                   sum(sc_class == "Doublet", na.rm = TRUE),
                   n_doublets),
      Percentage = c(round(100 * sum(df_class == "Doublet", na.rm = TRUE) / length(cells), 2),
                     round(100 * sum(sc_class == "Doublet", na.rm = TRUE) / length(cells), 2),
                     pct_doublets),
      Estimated_FPR = c(NA, NA, round(estimated_fpr, 4)),
      stringsAsFactors = FALSE
    )
  )
}

create_doublet_summary_table <- function(doublet_results) {
  classifications <- doublet_results$classifications
  summary_stats <- table(classifications)
  
  summary_df <- data.frame(
    Classification = names(summary_stats),
    Count = as.numeric(summary_stats),
    Percentage = round(as.numeric(summary_stats) / length(classifications) * 100, 2)
  )

  if (!is.null(doublet_results$method)) {
    summary_df$Method <- as.character(doublet_results$method)
  }
  if (!is.null(doublet_results$agreement_rate)) {
    summary_df$Agreement_Rate_Pct <- round(100 * as.numeric(doublet_results$agreement_rate), 2)
  }
  if (!is.null(doublet_results$disagreement_rate)) {
    summary_df$Disagreement_Rate_Pct <- round(100 * as.numeric(doublet_results$disagreement_rate), 2)
  }
  
  DT::datatable(summary_df, options = list(dom = 't'))
}

create_doublet_score_plot <- function(doublet_results) {
  scores <- doublet_results$scores
  classifications <- doublet_results$classifications
  
  plot_data <- data.frame(
    Score = scores,
    Classification = classifications
  )
  
  p <- ggplot(plot_data, aes(x = Score, fill = Classification)) +
    geom_histogram(alpha = 0.7, bins = 50) +
    scale_fill_manual(values = c("Singlet" = "blue", "Doublet" = "red")) +
    theme_minimal() +
    labs(title = "Doublet Score Distribution", x = "Doublet Score", y = "Count")
  
  ggplotly(p)
}

create_doublet_embedding_plot <- function(seurat_obj) {
  # Find best available reduction
  reductions <- names(seurat_obj@reductions)
  reduction <- if ("umap" %in% reductions) "umap" else if ("tsne" %in% reductions) "tsne" else "pca"
  
  if (!reduction %in% reductions) {
    return(plotly_empty() %>% layout(title = "No embedding available"))
  }
  
  embedding <- Embeddings(seurat_obj, reduction)[, 1:2]
  
  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Classification = seurat_obj$doublet_classification,
    Score = seurat_obj$doublet_score
  )
  
  plot_ly(
    data = plot_data,
    x = ~Dim1,
    y = ~Dim2,
    color = ~Classification,
    colors = c("Singlet" = "blue", "Doublet" = "red"),
    size = ~Score,
    type = "scatter",
    mode = "markers",
    marker = list(opacity = 0.7),
    hovertemplate = paste0(
      toupper(reduction), "1: %{x}<br>",
      toupper(reduction), "2: %{y}<br>",
      "Classification: %{color}<br>",
      "Score: %{marker.size:.3f}<extra></extra>"
    )
  ) %>%
    layout(
      title = paste("Doublets on", toupper(reduction), "embedding"),
      xaxis = list(title = paste0(toupper(reduction), "1")),
      yaxis = list(title = paste0(toupper(reduction), "2"))
    )
}

# ============================================================================
# PER-BATCH QC VIOLIN PLOTS
# ============================================================================

create_per_batch_qc_violins <- function(seurat_obj, batch_col) {
  # Creates faceted violin plots of QC metrics grouped by batch using plotly subplots
  metadata <- seurat_obj@meta.data

  if (!batch_col %in% colnames(metadata)) {
    return(plotly_empty() %>% layout(title = "Batch column not found"))
  }

  batch_vec <- as.character(metadata[[batch_col]])
  qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                        colnames(metadata))

  if (length(qc_cols) == 0) {
    return(plotly_empty() %>% layout(title = "No QC metrics available"))
  }

  metric_labels <- c(
    "nCount_RNA" = "UMI Count",
    "nFeature_RNA" = "Gene Count",
    "percent.mt" = "Mitochondrial %",
    "percent.ribo" = "Ribosomal %"
  )

  batches <- sort(unique(batch_vec))
  n_batches <- length(batches)
  pal <- if (n_batches <= 12) {
    RColorBrewer::brewer.pal(max(3, n_batches), "Set3")[1:n_batches]
  } else {
    grDevices::hcl.colors(n_batches, "Dynamic")
  }
  names(pal) <- batches

  subplot_list <- list()

  for (i in seq_along(qc_cols)) {
    metric <- qc_cols[i]
    label <- metric_labels[metric]

    p <- plot_ly()

    for (j in seq_along(batches)) {
      b <- batches[j]
      vals <- metadata[[metric]][batch_vec == b]

      p <- p %>% add_trace(
        y = vals,
        x = rep(b, length(vals)),
        type = "violin",
        name = b,
        legendgroup = b,
        showlegend = (i == 1),
        box = list(visible = TRUE),
        meanline = list(visible = TRUE),
        fillcolor = pal[j],
        line = list(color = pal[j]),
        opacity = 0.7
      )
    }

    p <- p %>% layout(
      yaxis = list(title = label),
      xaxis = list(title = "")
    )

    subplot_list[[i]] <- p
  }

  subplot(subplot_list, nrows = length(qc_cols), shareX = TRUE, titleY = TRUE) %>%
    layout(
      title = list(text = paste("Per-Batch QC Metrics by", batch_col)),
      showlegend = TRUE,
      autosize = TRUE
    )
}

# ============================================================================
# QC REPORT GENERATION
# ============================================================================

generate_qc_report <- function(seurat_obj, output_file, species = "human") {
  # Generate comprehensive HTML QC report
  # Includes: dataset summary, QC metric distributions (violin/histogram),
  # MT%/ribo% stats, cell cycle phase distribution, batch composition,
  # scatter plots (genes vs UMIs, MT% vs genes), filtering thresholds applied,
  # complexity scores

  tryCatch({
    metadata <- seurat_obj@meta.data
    n_cells <- ncol(seurat_obj)
    n_genes <- nrow(seurat_obj)
    assay_names <- paste(names(seurat_obj@assays), collapse = ", ")

    # ---- CSS Styling (matching enrichment report style) ----
    css <- '
    <style>
      body { font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #f9f9f9; color: #333; }
      .report-container { max-width: 1100px; margin: 0 auto; background: white; padding: 40px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
      h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
      h2 { color: #2980b9; border-bottom: 1px solid #ddd; padding-bottom: 8px; margin-top: 30px; }
      h3 { color: #34495e; margin-top: 20px; }
      table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 13px; }
      th { background-color: #3498db; color: white; padding: 10px 12px; text-align: left; }
      td { padding: 8px 12px; border-bottom: 1px solid #eee; }
      tr:nth-child(even) { background-color: #f2f6fc; }
      tr:hover { background-color: #e8f0fe; }
      .param-table td:first-child { font-weight: bold; width: 250px; }
      .summary-box { background: #eaf2f8; border-left: 4px solid #3498db; padding: 15px; margin: 15px 0; border-radius: 4px; }
      .metric-box { background: #f0f7ee; border-left: 4px solid #27ae60; padding: 15px; margin: 15px 0; border-radius: 4px; }
      .warning-box { background: #fdf2f2; border-left: 4px solid #e74c3c; padding: 15px; margin: 15px 0; border-radius: 4px; }
      .plot-container { text-align: center; margin: 20px 0; }
      .plot-container img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }
      .footer { margin-top: 40px; padding-top: 15px; border-top: 1px solid #ddd; font-size: 11px; color: #999; text-align: center; }
      .grid-2col { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
      .badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: bold; }
      .badge-info { background: #d6eaf8; color: #2471a3; }
      .badge-warn { background: #fce4e4; color: #c0392b; }
    </style>'

    # ---- Helper: build HTML table from data.frame ----
    df_to_html_table <- function(df, max_rows = 50, extra_class = "") {
      if (is.null(df) || nrow(df) == 0) return("<p><em>No data available.</em></p>")
      df <- head(df, max_rows)
      header <- paste0("<th>", colnames(df), "</th>", collapse = "")
      rows <- apply(df, 1, function(row) {
        paste0("<td>", htmltools::htmlEscape(as.character(row)), "</td>", collapse = "")
      })
      rows_html <- paste0("<tr>", rows, "</tr>", collapse = "\n")
      paste0('<table class="', extra_class, '"><thead><tr>', header, '</tr></thead><tbody>\n', rows_html, '\n</tbody></table>')
    }

    # ---- Helper: generate base64-encoded PNG from a plot expression ----
    plot_to_base64 <- function(plot_expr, width = 700, height = 450) {
      tmp <- tempfile(fileext = ".png")
      on.exit(unlink(tmp), add = TRUE)
      grDevices::png(tmp, width = width, height = height, res = 120)
      tryCatch({
        eval(plot_expr)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Plot error:", e$message), cex = 0.8)
      })
      grDevices::dev.off()

      if (file.exists(tmp) && file.size(tmp) > 0) {
        encoded <- base64enc::base64encode(tmp)
        return(paste0("data:image/png;base64,", encoded))
      }
      return(NULL)
    }

    # ---- Header ----
    header_html <- paste0(
      '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n<title>QC Report</title>\n',
      css, '\n</head>\n<body>\n',
      '<div class="report-container">\n',
      '<h1>Single-Cell RNA-seq QC Report</h1>\n',
      '<p><strong>Generated:</strong> ', format(Sys.time(), "%B %d, %Y at %H:%M:%S"), '</p>\n',
      '<p><strong>Species:</strong> ', htmltools::htmlEscape(species), '</p>\n'
    )

    # ---- Dataset Summary ----
    summary_html <- paste0(
      '<h2>Dataset Summary</h2>\n',
      '<div class="summary-box">\n',
      '<strong>Number of Cells:</strong> ', format(n_cells, big.mark = ","), '<br>\n',
      '<strong>Number of Genes:</strong> ', format(n_genes, big.mark = ","), '<br>\n',
      '<strong>Assays:</strong> ', htmltools::htmlEscape(assay_names), '<br>\n',
      '<strong>Metadata Columns:</strong> ', ncol(metadata), '\n',
      '</div>\n'
    )

    # ---- QC Metric Summary Table ----
    qc_metrics_cols <- intersect(
      c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
        "log10GenesPerUMI", "novelty_score", "S.Score", "G2M.Score"),
      colnames(metadata)
    )

    metrics_html <- '<h2>QC Metrics Summary</h2>\n'
    if (length(qc_metrics_cols) > 0) {
      qc_summary <- do.call(rbind, lapply(qc_metrics_cols, function(col) {
        vals <- metadata[[col]]
        data.frame(
          Metric = col,
          Mean = round(mean(vals, na.rm = TRUE), 3),
          Median = round(median(vals, na.rm = TRUE), 3),
          SD = round(sd(vals, na.rm = TRUE), 3),
          Min = round(min(vals, na.rm = TRUE), 3),
          Max = round(max(vals, na.rm = TRUE), 3),
          Q1 = round(quantile(vals, 0.25, na.rm = TRUE), 3),
          Q3 = round(quantile(vals, 0.75, na.rm = TRUE), 3),
          stringsAsFactors = FALSE
        )
      }))
      metrics_html <- paste0(metrics_html, df_to_html_table(qc_summary))
    } else {
      metrics_html <- paste0(metrics_html, '<p><em>No QC metrics computed yet.</em></p>\n')
    }

    # ---- MT% and Ribo% Statistics ----
    mt_ribo_html <- '<h2>Mitochondrial &amp; Ribosomal Statistics</h2>\n'
    if ("percent.mt" %in% colnames(metadata)) {
      mt_vals <- metadata$percent.mt
      mt_high <- sum(mt_vals > 20, na.rm = TRUE)
      mt_ribo_html <- paste0(mt_ribo_html,
        '<div class="metric-box">\n',
        '<h3>Mitochondrial %</h3>\n',
        '<strong>Mean:</strong> ', round(mean(mt_vals, na.rm = TRUE), 2), '% | ',
        '<strong>Median:</strong> ', round(median(mt_vals, na.rm = TRUE), 2), '% | ',
        '<strong>Max:</strong> ', round(max(mt_vals, na.rm = TRUE), 2), '%<br>\n',
        '<strong>Cells with MT% > 20:</strong> ', mt_high,
        ' (', round(mt_high / n_cells * 100, 1), '%)\n',
        '</div>\n'
      )
    }
    if ("percent.ribo" %in% colnames(metadata)) {
      ribo_vals <- metadata$percent.ribo
      mt_ribo_html <- paste0(mt_ribo_html,
        '<div class="metric-box">\n',
        '<h3>Ribosomal %</h3>\n',
        '<strong>Mean:</strong> ', round(mean(ribo_vals, na.rm = TRUE), 2), '% | ',
        '<strong>Median:</strong> ', round(median(ribo_vals, na.rm = TRUE), 2), '% | ',
        '<strong>Max:</strong> ', round(max(ribo_vals, na.rm = TRUE), 2), '%\n',
        '</div>\n'
      )
    }

    # ---- Distribution Histograms (base64 PNGs) ----
    plots_html <- '<h2>QC Metric Distributions</h2>\n'

    hist_metrics <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                              colnames(metadata))
    hist_labels <- c("nCount_RNA" = "UMI Count", "nFeature_RNA" = "Gene Count",
                     "percent.mt" = "Mitochondrial %", "percent.ribo" = "Ribosomal %")
    hist_colors <- c("nCount_RNA" = "#3498db", "nFeature_RNA" = "#2ecc71",
                     "percent.mt" = "#e74c3c", "percent.ribo" = "#9b59b6")

    if (length(hist_metrics) > 0) {
      # Combined histogram plot
      hist_b64 <- plot_to_base64(quote({
        n_plots <- length(hist_metrics)
        par(mfrow = c(1, n_plots), mar = c(4, 4, 3, 1))
        for (m in hist_metrics) {
          hist(metadata[[m]], breaks = 50, col = hist_colors[m], border = "white",
               main = hist_labels[m], xlab = hist_labels[m], ylab = "Frequency")
          abline(v = median(metadata[[m]], na.rm = TRUE), col = "black", lwd = 2, lty = 2)
        }
      }), width = 250 * length(hist_metrics), height = 350)

      if (!is.null(hist_b64)) {
        plots_html <- paste0(plots_html,
          '<div class="plot-container">\n',
          '<img src="', hist_b64, '" alt="QC Metric Histograms">\n',
          '<p><em>Dashed line indicates the median value.</em></p>\n',
          '</div>\n'
        )
      }
    }

    # ---- Scatter Plots (base64 PNGs) ----
    scatter_html <- '<h2>QC Scatter Plots</h2>\n'

    # Genes vs UMIs
    if (all(c("nCount_RNA", "nFeature_RNA") %in% colnames(metadata))) {
      scatter1_b64 <- plot_to_base64(quote({
        par(mar = c(4, 4, 3, 1))
        smoothScatter(metadata$nCount_RNA, metadata$nFeature_RNA,
                      xlab = "UMI Count", ylab = "Gene Count",
                      main = "Genes vs UMIs",
                      colramp = colorRampPalette(c("white", "#3498db", "#2c3e50")))
      }), width = 500, height = 420)

      if (!is.null(scatter1_b64)) {
        scatter_html <- paste0(scatter_html,
          '<div class="grid-2col">\n<div class="plot-container">\n',
          '<img src="', scatter1_b64, '" alt="Genes vs UMIs">\n</div>\n'
        )
      }
    }

    # MT% vs Genes
    if (all(c("nFeature_RNA", "percent.mt") %in% colnames(metadata))) {
      scatter2_b64 <- plot_to_base64(quote({
        par(mar = c(4, 4, 3, 1))
        plot(metadata$nFeature_RNA, metadata$percent.mt,
             pch = 16, cex = 0.3, col = adjustcolor("#e74c3c", alpha.f = 0.3),
             xlab = "Gene Count", ylab = "Mitochondrial %",
             main = "MT% vs Genes")
        abline(h = 20, col = "red", lty = 2, lwd = 2)
      }), width = 500, height = 420)

      if (!is.null(scatter2_b64)) {
        scatter_html <- paste0(scatter_html,
          '<div class="plot-container">\n',
          '<img src="', scatter2_b64, '" alt="MT% vs Genes">\n</div>\n</div>\n'
        )
      } else {
        scatter_html <- paste0(scatter_html, '</div>\n')
      }
    } else {
      scatter_html <- paste0(scatter_html, '</div>\n')
    }

    # ---- Cell Cycle Phase Distribution ----
    cc_html <- '<h2>Cell Cycle Phase Distribution</h2>\n'
    if ("Phase" %in% colnames(metadata)) {
      phase_counts <- table(metadata$Phase)
      phase_df <- data.frame(
        Phase = names(phase_counts),
        Count = as.numeric(phase_counts),
        Percentage = paste0(round(as.numeric(phase_counts) / n_cells * 100, 1), "%"),
        stringsAsFactors = FALSE
      )
      cc_html <- paste0(cc_html, df_to_html_table(phase_df))

      # Mean phase scores (when available from Preprocessing scoring)
      if (all(c("S.Score", "G2M.Score") %in% colnames(metadata))) {
        cc_html <- paste0(cc_html,
          sprintf('<p>Mean S.Score = %.3f &nbsp;|&nbsp; Mean G2M.Score = %.3f</p>\n',
                  mean(metadata$S.Score, na.rm = TRUE),
                  mean(metadata$G2M.Score, na.rm = TRUE)))
      }

      # Pie chart
      cc_b64 <- plot_to_base64(quote({
        par(mar = c(1, 1, 3, 1))
        pie(as.numeric(phase_counts),
            labels = paste0(names(phase_counts), "\n", round(as.numeric(phase_counts)/sum(phase_counts)*100, 1), "%"),
            col = c("#2ecc71", "#e74c3c", "#3498db"),
            main = "Cell Cycle Phase Distribution")
      }), width = 450, height = 400)

      if (!is.null(cc_b64)) {
        cc_html <- paste0(cc_html,
          '<div class="plot-container">\n',
          '<img src="', cc_b64, '" alt="Cell Cycle Phases">\n</div>\n'
        )
      }
    } else {
      cc_html <- paste0(cc_html,
        '<p><em>Cell cycle not scored. Run "Cell Cycle Scoring" in the Preprocessing tab ',
        '(requires normalized data) to include phase distribution here.</em></p>\n')
    }

    # ---- Batch Composition ----
    batch_html <- '<h2>Batch Composition</h2>\n'
    batch_col_candidates <- grep("batch|orig.ident|sample|donor", colnames(metadata),
                                  value = TRUE, ignore.case = TRUE)
    if ("orig.ident" %in% colnames(metadata)) {
      batch_col_candidates <- unique(c("orig.ident", batch_col_candidates))
    }

    if (length(batch_col_candidates) > 0) {
      bcol <- batch_col_candidates[1]
      batch_vec <- as.character(metadata[[bcol]])
      batch_counts <- table(batch_vec)

      batch_df_list <- lapply(names(batch_counts), function(b) {
        cells_idx <- which(batch_vec == b)
        df <- data.frame(
          Batch = b,
          N_Cells = length(cells_idx),
          Pct_Cells = paste0(round(length(cells_idx) / n_cells * 100, 1), "%"),
          stringsAsFactors = FALSE
        )
        if ("nCount_RNA" %in% colnames(metadata)) {
          df$Median_UMI <- round(median(metadata$nCount_RNA[cells_idx], na.rm = TRUE))
        }
        if ("nFeature_RNA" %in% colnames(metadata)) {
          df$Median_Genes <- round(median(metadata$nFeature_RNA[cells_idx], na.rm = TRUE))
        }
        if ("percent.mt" %in% colnames(metadata)) {
          df$Median_MT_Pct <- round(median(metadata$percent.mt[cells_idx], na.rm = TRUE), 2)
        }
        if ("percent.ribo" %in% colnames(metadata)) {
          df$Median_Ribo_Pct <- round(median(metadata$percent.ribo[cells_idx], na.rm = TRUE), 2)
        }
        df
      })
      batch_df <- do.call(rbind, batch_df_list)
      batch_html <- paste0(batch_html,
        '<p><strong>Batch column used:</strong> ', htmltools::htmlEscape(bcol), '</p>\n',
        df_to_html_table(batch_df)
      )
    } else {
      batch_html <- paste0(batch_html, '<p><em>No batch information found.</em></p>\n')
    }

    # ---- Complexity Scores ----
    complexity_html <- '<h2>Complexity Scores</h2>\n'
    if ("log10GenesPerUMI" %in% colnames(metadata)) {
      comp_vals <- metadata$log10GenesPerUMI
      low_complexity <- sum(comp_vals < 0.8, na.rm = TRUE)
      complexity_html <- paste0(complexity_html,
        '<div class="metric-box">\n',
        '<strong>log10(Genes/UMI) - Mean:</strong> ', round(mean(comp_vals, na.rm = TRUE), 4), ' | ',
        '<strong>Median:</strong> ', round(median(comp_vals, na.rm = TRUE), 4), '<br>\n',
        '<strong>Low complexity cells (< 0.8):</strong> ', low_complexity,
        ' (', round(low_complexity / n_cells * 100, 1), '%)\n',
        '</div>\n'
      )
    }
    if ("novelty_score" %in% colnames(metadata)) {
      nov_vals <- metadata$novelty_score
      complexity_html <- paste0(complexity_html,
        '<div class="metric-box">\n',
        '<strong>Novelty Score - Mean:</strong> ', round(mean(nov_vals, na.rm = TRUE), 4), ' | ',
        '<strong>Median:</strong> ', round(median(nov_vals, na.rm = TRUE), 4), '\n',
        '</div>\n'
      )
    }
    if (!"log10GenesPerUMI" %in% colnames(metadata) && !"novelty_score" %in% colnames(metadata)) {
      complexity_html <- paste0(complexity_html, '<p><em>Complexity scores not computed.</em></p>\n')
    }

    # ---- Footer ----
    footer_html <- paste0(
      '<div class="footer">\n',
      '<p>Generated by scRNA4 QC Module | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>\n',
      '</div>\n',
      '</div>\n</body>\n</html>'
    )

    # ---- Assemble and write report ----
    full_html <- paste0(
      header_html,
      summary_html,
      metrics_html,
      mt_ribo_html,
      plots_html,
      scatter_html,
      cc_html,
      batch_html,
      complexity_html,
      footer_html
    )

    writeLines(full_html, output_file)
    return(invisible(TRUE))

  }, error = function(e) {
    # Write a minimal error report
    writeLines(paste0(
      '<!DOCTYPE html><html><body><h1>QC Report Error</h1><p>',
      htmltools::htmlEscape(e$message), '</p></body></html>'
    ), output_file)
    return(invisible(FALSE))
  })
}

# ============================================================================
# QC REPORT UI & SERVER MODULES
# ============================================================================

qcReportUI <- function(id) {
  ns <- NS(id)
  tagList(
    wellPanel(
      h4("QC Report Generation"),
      p("Generate a comprehensive HTML report of QC metrics and visualizations."),

      fluidRow(
        column(12,
          h5("Include in Report:"),
          fluidRow(
            column(4,
              checkboxInput(ns("include_metrics_summary"), "Metrics Summary", TRUE)
            ),
            column(4,
              checkboxInput(ns("include_violin_plots"), "Violin Plots", TRUE)
            ),
            column(4,
              checkboxInput(ns("include_scatter_plots"), "Scatter Plots", TRUE)
            )
          ),
          fluidRow(
            column(4,
              checkboxInput(ns("include_cell_cycle"), "Cell Cycle Info", TRUE)
            ),
            column(4,
              checkboxInput(ns("include_batch_info"), "Batch Info", TRUE)
            )
          )
        )
      ),

      br(),
      downloadButton(ns("download_qc_report"), "Generate QC Report", class = "btn-primary btn-lg")
    )
  )
}

qcReportServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    output$download_qc_report <- downloadHandler(
      filename = function() {
        paste0("qc_report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
      },
      content = function(file) {
        req(values$active_dataset)

        withProgress(message = "Generating QC report...", value = 0, {

          seurat_obj <- values$active_dataset
          metadata <- seurat_obj@meta.data
          n_cells <- ncol(seurat_obj)
          n_genes <- nrow(seurat_obj)
          assay_names <- paste(names(seurat_obj@assays), collapse = ", ")

          species <- if (!is.null(values$species)) values$species else "human"

          incProgress(0.1, detail = "Building HTML structure...")

          # ---- CSS Styling ----
          css <- '
          <style>
            body { font-family: "Segoe UI", Tahoma, Geneva, Verdana, sans-serif; margin: 40px; background: #f9f9f9; color: #333; }
            .report-container { max-width: 1100px; margin: 0 auto; background: white; padding: 40px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
            h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
            h2 { color: #2980b9; border-bottom: 1px solid #ddd; padding-bottom: 8px; margin-top: 30px; }
            h3 { color: #34495e; margin-top: 20px; }
            table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 13px; }
            th { background-color: #3498db; color: white; padding: 10px 12px; text-align: left; }
            td { padding: 8px 12px; border-bottom: 1px solid #eee; }
            tr:nth-child(even) { background-color: #f2f6fc; }
            tr:hover { background-color: #e8f0fe; }
            .param-table td:first-child { font-weight: bold; width: 250px; }
            .summary-box { background: #eaf2f8; border-left: 4px solid #3498db; padding: 15px; margin: 15px 0; border-radius: 4px; }
            .metric-box { background: #f0f7ee; border-left: 4px solid #27ae60; padding: 15px; margin: 15px 0; border-radius: 4px; }
            .warning-box { background: #fdf2f2; border-left: 4px solid #e74c3c; padding: 15px; margin: 15px 0; border-radius: 4px; }
            .plot-container { text-align: center; margin: 20px 0; }
            .plot-container img { max-width: 100%; border: 1px solid #ddd; border-radius: 4px; }
            .footer { margin-top: 40px; padding-top: 15px; border-top: 1px solid #ddd; font-size: 11px; color: #999; text-align: center; }
            .grid-2col { display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }
            .badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 11px; font-weight: bold; }
            .badge-info { background: #d6eaf8; color: #2471a3; }
          </style>'

          # ---- Helper: build HTML table from data.frame ----
          df_to_html_table <- function(df, max_rows = 50, extra_class = "") {
            if (is.null(df) || nrow(df) == 0) return("<p><em>No data available.</em></p>")
            df <- head(df, max_rows)
            header <- paste0("<th>", colnames(df), "</th>", collapse = "")
            rows <- apply(df, 1, function(row) {
              paste0("<td>", htmltools::htmlEscape(as.character(row)), "</td>", collapse = "")
            })
            rows_html <- paste0("<tr>", rows, "</tr>", collapse = "\n")
            paste0('<table class="', extra_class, '"><thead><tr>', header, '</tr></thead><tbody>\n', rows_html, '\n</tbody></table>')
          }

          # ---- Helper: generate base64-encoded PNG ----
          plot_to_base64 <- function(plot_fn, width = 700, height = 450) {
            tmp <- tempfile(fileext = ".png")
            on.exit(unlink(tmp), add = TRUE)
            grDevices::png(tmp, width = width, height = height, res = 120)
            tryCatch({
              plot_fn()
            }, error = function(e) {
              plot.new()
              text(0.5, 0.5, paste("Plot error:", e$message), cex = 0.8)
            })
            grDevices::dev.off()

            if (file.exists(tmp) && file.size(tmp) > 0) {
              encoded <- base64enc::base64encode(tmp)
              return(paste0("data:image/png;base64,", encoded))
            }
            return(NULL)
          }

          # ---- Header ----
          header_html <- paste0(
            '<!DOCTYPE html>\n<html lang="en">\n<head>\n<meta charset="UTF-8">\n<title>QC Report</title>\n',
            css, '\n</head>\n<body>\n',
            '<div class="report-container">\n',
            '<h1>Single-Cell RNA-seq QC Report</h1>\n',
            '<p><strong>Generated:</strong> ', format(Sys.time(), "%B %d, %Y at %H:%M:%S"), '</p>\n',
            '<p><strong>Species:</strong> ', htmltools::htmlEscape(species), '</p>\n'
          )

          incProgress(0.2, detail = "Dataset summary...")

          # ---- Dataset Info ----
          info_html <- paste0(
            '<h2>Dataset Information</h2>\n',
            '<div class="summary-box">\n',
            '<strong>Number of Cells:</strong> ', format(n_cells, big.mark = ","), '<br>\n',
            '<strong>Number of Genes:</strong> ', format(n_genes, big.mark = ","), '<br>\n',
            '<strong>Assays:</strong> ', htmltools::htmlEscape(assay_names), '<br>\n',
            '<strong>Metadata Columns:</strong> ', ncol(metadata), '\n',
            '</div>\n'
          )

          # ---- QC Metric Summary Table (conditional) ----
          metrics_html <- ""
          if (input$include_metrics_summary) {
            qc_metrics_cols <- intersect(
              c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo",
                "log10GenesPerUMI", "novelty_score", "S.Score", "G2M.Score"),
              colnames(metadata)
            )

            metrics_html <- '<h2>QC Metrics Summary</h2>\n'
            if (length(qc_metrics_cols) > 0) {
              qc_summary <- do.call(rbind, lapply(qc_metrics_cols, function(col) {
                vals <- metadata[[col]]
                data.frame(
                  Metric = col,
                  Mean = round(mean(vals, na.rm = TRUE), 3),
                  Median = round(median(vals, na.rm = TRUE), 3),
                  SD = round(sd(vals, na.rm = TRUE), 3),
                  Min = round(min(vals, na.rm = TRUE), 3),
                  Max = round(max(vals, na.rm = TRUE), 3),
                  Q1 = round(quantile(vals, 0.25, na.rm = TRUE), 3),
                  Q3 = round(quantile(vals, 0.75, na.rm = TRUE), 3),
                  stringsAsFactors = FALSE
                )
              }))
              metrics_html <- paste0(metrics_html, df_to_html_table(qc_summary))

              # MT% and Ribo% detailed stats
              if ("percent.mt" %in% colnames(metadata)) {
                mt_vals <- metadata$percent.mt
                mt_high <- sum(mt_vals > 20, na.rm = TRUE)
                metrics_html <- paste0(metrics_html,
                  '<div class="metric-box">\n',
                  '<h3>Mitochondrial %</h3>\n',
                  '<strong>Mean:</strong> ', round(mean(mt_vals, na.rm = TRUE), 2), '% | ',
                  '<strong>Median:</strong> ', round(median(mt_vals, na.rm = TRUE), 2), '% | ',
                  '<strong>Max:</strong> ', round(max(mt_vals, na.rm = TRUE), 2), '%<br>\n',
                  '<strong>Cells with MT% > 20:</strong> ', mt_high,
                  ' (', round(mt_high / n_cells * 100, 1), '%)\n',
                  '</div>\n'
                )
              }
              if ("percent.ribo" %in% colnames(metadata)) {
                ribo_vals <- metadata$percent.ribo
                metrics_html <- paste0(metrics_html,
                  '<div class="metric-box">\n',
                  '<h3>Ribosomal %</h3>\n',
                  '<strong>Mean:</strong> ', round(mean(ribo_vals, na.rm = TRUE), 2), '% | ',
                  '<strong>Median:</strong> ', round(median(ribo_vals, na.rm = TRUE), 2), '% | ',
                  '<strong>Max:</strong> ', round(max(ribo_vals, na.rm = TRUE), 2), '%\n',
                  '</div>\n'
                )
              }
            } else {
              metrics_html <- paste0(metrics_html, '<p><em>No QC metrics computed yet.</em></p>\n')
            }
          }

          incProgress(0.4, detail = "Generating plots...")

          # ---- Distribution Histograms (conditional) ----
          violin_html <- ""
          if (input$include_violin_plots) {
            hist_metrics <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                                      colnames(metadata))
            hist_labels <- c("nCount_RNA" = "UMI Count", "nFeature_RNA" = "Gene Count",
                             "percent.mt" = "Mitochondrial %", "percent.ribo" = "Ribosomal %")
            hist_colors <- c("nCount_RNA" = "#3498db", "nFeature_RNA" = "#2ecc71",
                             "percent.mt" = "#e74c3c", "percent.ribo" = "#9b59b6")

            violin_html <- '<h2>QC Metric Distributions</h2>\n'

            if (length(hist_metrics) > 0) {
              # Histograms
              hist_b64 <- plot_to_base64(function() {
                n_plots <- length(hist_metrics)
                par(mfrow = c(1, n_plots), mar = c(4, 4, 3, 1))
                for (m in hist_metrics) {
                  hist(metadata[[m]], breaks = 50, col = hist_colors[m], border = "white",
                       main = hist_labels[m], xlab = hist_labels[m], ylab = "Frequency")
                  abline(v = median(metadata[[m]], na.rm = TRUE), col = "black", lwd = 2, lty = 2)
                }
              }, width = 250 * length(hist_metrics), height = 350)

              if (!is.null(hist_b64)) {
                violin_html <- paste0(violin_html,
                  '<div class="plot-container">\n',
                  '<img src="', hist_b64, '" alt="QC Metric Histograms">\n',
                  '<p><em>Dashed line indicates the median value.</em></p>\n',
                  '</div>\n'
                )
              }

              # Violin plots
              violin_b64 <- plot_to_base64(function() {
                n_plots <- length(hist_metrics)
                par(mfrow = c(1, n_plots), mar = c(4, 4, 3, 1))
                for (m in hist_metrics) {
                  vals <- metadata[[m]]
                  if (length(vals) > 0) {
                    vioplot_data <- list(vals)
                    boxplot(vals, col = hist_colors[m], main = hist_labels[m],
                            ylab = hist_labels[m], outline = FALSE, border = "#333")
                    stripchart(vals, vertical = TRUE, method = "jitter", add = TRUE,
                               pch = 16, cex = 0.2, col = adjustcolor(hist_colors[m], alpha.f = 0.3))
                  }
                }
              }, width = 250 * length(hist_metrics), height = 350)

              if (!is.null(violin_b64)) {
                violin_html <- paste0(violin_html,
                  '<div class="plot-container">\n',
                  '<img src="', violin_b64, '" alt="QC Metric Box Plots">\n',
                  '</div>\n'
                )
              }
            }
          }

          incProgress(0.5, detail = "Scatter plots...")

          # ---- Scatter Plots (conditional) ----
          scatter_html <- ""
          if (input$include_scatter_plots) {
            scatter_html <- '<h2>QC Scatter Plots</h2>\n<div class="grid-2col">\n'

            # Genes vs UMIs
            if (all(c("nCount_RNA", "nFeature_RNA") %in% colnames(metadata))) {
              s1_b64 <- plot_to_base64(function() {
                par(mar = c(4, 4, 3, 1))
                smoothScatter(metadata$nCount_RNA, metadata$nFeature_RNA,
                              xlab = "UMI Count", ylab = "Gene Count",
                              main = "Genes vs UMIs",
                              colramp = colorRampPalette(c("white", "#3498db", "#2c3e50")))
              }, width = 500, height = 420)

              if (!is.null(s1_b64)) {
                scatter_html <- paste0(scatter_html,
                  '<div class="plot-container">\n',
                  '<img src="', s1_b64, '" alt="Genes vs UMIs">\n</div>\n'
                )
              }
            }

            # MT% vs Genes
            if (all(c("nFeature_RNA", "percent.mt") %in% colnames(metadata))) {
              s2_b64 <- plot_to_base64(function() {
                par(mar = c(4, 4, 3, 1))
                plot(metadata$nFeature_RNA, metadata$percent.mt,
                     pch = 16, cex = 0.3, col = adjustcolor("#e74c3c", alpha.f = 0.3),
                     xlab = "Gene Count", ylab = "Mitochondrial %",
                     main = "MT% vs Genes")
                abline(h = 20, col = "red", lty = 2, lwd = 2)
              }, width = 500, height = 420)

              if (!is.null(s2_b64)) {
                scatter_html <- paste0(scatter_html,
                  '<div class="plot-container">\n',
                  '<img src="', s2_b64, '" alt="MT% vs Genes">\n</div>\n'
                )
              }
            }

            # UMI vs MT%
            if (all(c("nCount_RNA", "percent.mt") %in% colnames(metadata))) {
              s3_b64 <- plot_to_base64(function() {
                par(mar = c(4, 4, 3, 1))
                plot(metadata$nCount_RNA, metadata$percent.mt,
                     pch = 16, cex = 0.3, col = adjustcolor("#9b59b6", alpha.f = 0.3),
                     xlab = "UMI Count", ylab = "Mitochondrial %",
                     main = "MT% vs UMIs")
                abline(h = 20, col = "red", lty = 2, lwd = 2)
              }, width = 500, height = 420)

              if (!is.null(s3_b64)) {
                scatter_html <- paste0(scatter_html,
                  '<div class="plot-container">\n',
                  '<img src="', s3_b64, '" alt="MT% vs UMIs">\n</div>\n'
                )
              }
            }

            # Complexity plot
            if ("log10GenesPerUMI" %in% colnames(metadata) && "nCount_RNA" %in% colnames(metadata)) {
              s4_b64 <- plot_to_base64(function() {
                par(mar = c(4, 4, 3, 1))
                plot(metadata$nCount_RNA, metadata$log10GenesPerUMI,
                     pch = 16, cex = 0.3, col = adjustcolor("#27ae60", alpha.f = 0.3),
                     xlab = "UMI Count", ylab = "log10(Genes/UMI)",
                     main = "Complexity vs UMIs")
                abline(h = 0.8, col = "red", lty = 2, lwd = 2)
              }, width = 500, height = 420)

              if (!is.null(s4_b64)) {
                scatter_html <- paste0(scatter_html,
                  '<div class="plot-container">\n',
                  '<img src="', s4_b64, '" alt="Complexity vs UMIs">\n</div>\n'
                )
              }
            }

            scatter_html <- paste0(scatter_html, '</div>\n')
          }

          incProgress(0.7, detail = "Cell cycle and batch info...")

          # ---- Cell Cycle Phase (conditional) ----
          cc_html <- ""
          if (input$include_cell_cycle) {
            cc_html <- '<h2>Cell Cycle Phase Distribution</h2>\n'
            if ("Phase" %in% colnames(metadata)) {
              phase_counts <- table(metadata$Phase)
              phase_df <- data.frame(
                Phase = names(phase_counts),
                Count = as.numeric(phase_counts),
                Percentage = paste0(round(as.numeric(phase_counts) / n_cells * 100, 1), "%"),
                stringsAsFactors = FALSE
              )
              cc_html <- paste0(cc_html, df_to_html_table(phase_df))

              cc_b64 <- plot_to_base64(function() {
                par(mar = c(1, 1, 3, 1))
                pie(as.numeric(phase_counts),
                    labels = paste0(names(phase_counts), "\n",
                                    round(as.numeric(phase_counts)/sum(phase_counts)*100, 1), "%"),
                    col = c("#2ecc71", "#e74c3c", "#3498db"),
                    main = "Cell Cycle Phase Distribution")
              }, width = 450, height = 400)

              if (!is.null(cc_b64)) {
                cc_html <- paste0(cc_html,
                  '<div class="plot-container">\n',
                  '<img src="', cc_b64, '" alt="Cell Cycle Phases">\n</div>\n'
                )
              }

              # S.Score and G2M.Score details
              if (all(c("S.Score", "G2M.Score") %in% colnames(metadata))) {
                cc_html <- paste0(cc_html,
                  '<div class="metric-box">\n',
                  '<strong>S.Score</strong> - Mean: ', round(mean(metadata$S.Score, na.rm = TRUE), 4),
                  ' | Median: ', round(median(metadata$S.Score, na.rm = TRUE), 4), '<br>\n',
                  '<strong>G2M.Score</strong> - Mean: ', round(mean(metadata$G2M.Score, na.rm = TRUE), 4),
                  ' | Median: ', round(median(metadata$G2M.Score, na.rm = TRUE), 4), '\n',
                  '</div>\n'
                )
              }
            } else {
              cc_html <- paste0(cc_html, '<p><em>Cell cycle scoring not performed.</em></p>\n')
            }
          }

          # ---- Batch Info (conditional) ----
          batch_html <- ""
          if (input$include_batch_info) {
            batch_html <- '<h2>Batch Composition &amp; Per-Batch Statistics</h2>\n'
            batch_col_candidates <- grep("batch|orig.ident|sample|donor", colnames(metadata),
                                          value = TRUE, ignore.case = TRUE)
            if ("orig.ident" %in% colnames(metadata)) {
              batch_col_candidates <- unique(c("orig.ident", batch_col_candidates))
            }

            if (length(batch_col_candidates) > 0) {
              bcol <- batch_col_candidates[1]
              batch_vec <- as.character(metadata[[bcol]])
              batch_counts <- table(batch_vec)

              batch_df_list <- lapply(names(batch_counts), function(b) {
                cells_idx <- which(batch_vec == b)
                df <- data.frame(
                  Batch = b,
                  N_Cells = length(cells_idx),
                  Pct_Cells = paste0(round(length(cells_idx) / n_cells * 100, 1), "%"),
                  stringsAsFactors = FALSE
                )
                if ("nCount_RNA" %in% colnames(metadata)) {
                  df$Median_UMI <- round(median(metadata$nCount_RNA[cells_idx], na.rm = TRUE))
                }
                if ("nFeature_RNA" %in% colnames(metadata)) {
                  df$Median_Genes <- round(median(metadata$nFeature_RNA[cells_idx], na.rm = TRUE))
                }
                if ("percent.mt" %in% colnames(metadata)) {
                  df$Median_MT_Pct <- round(median(metadata$percent.mt[cells_idx], na.rm = TRUE), 2)
                }
                if ("percent.ribo" %in% colnames(metadata)) {
                  df$Median_Ribo_Pct <- round(median(metadata$percent.ribo[cells_idx], na.rm = TRUE), 2)
                }
                df
              })
              batch_df <- do.call(rbind, batch_df_list)
              batch_html <- paste0(batch_html,
                '<p><strong>Batch column:</strong> ', htmltools::htmlEscape(bcol), '</p>\n',
                df_to_html_table(batch_df)
              )

              # Per-batch QC violin plot
              batch_violin_b64 <- plot_to_base64(function() {
                qc_cols_batch <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"),
                                           colnames(metadata))
                n_qc <- length(qc_cols_batch)
                if (n_qc > 0) {
                  par(mfrow = c(n_qc, 1), mar = c(3, 4, 2, 1))
                  batch_labels <- c("nCount_RNA" = "UMI Count", "nFeature_RNA" = "Gene Count",
                                    "percent.mt" = "MT %", "percent.ribo" = "Ribo %")
                  for (m in qc_cols_batch) {
                    vals_by_batch <- split(metadata[[m]], batch_vec)
                    boxplot(vals_by_batch, col = grDevices::hcl.colors(length(vals_by_batch), "Dynamic"),
                            main = batch_labels[m], ylab = batch_labels[m], las = 2, outline = FALSE)
                  }
                }
              }, width = 800, height = 250 * max(1, length(intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), colnames(metadata)))))

              if (!is.null(batch_violin_b64)) {
                batch_html <- paste0(batch_html,
                  '<div class="plot-container">\n',
                  '<img src="', batch_violin_b64, '" alt="Per-Batch QC Distributions">\n</div>\n'
                )
              }
            } else {
              batch_html <- paste0(batch_html, '<p><em>No batch information found.</em></p>\n')
            }
          }

          incProgress(0.85, detail = "Complexity scores and filtering history...")

          # ---- Complexity Scores ----
          complexity_html <- '<h2>Complexity Scores</h2>\n'
          if ("log10GenesPerUMI" %in% colnames(metadata)) {
            comp_vals <- metadata$log10GenesPerUMI
            low_complexity <- sum(comp_vals < 0.8, na.rm = TRUE)
            complexity_html <- paste0(complexity_html,
              '<div class="metric-box">\n',
              '<strong>log10(Genes/UMI) - Mean:</strong> ', round(mean(comp_vals, na.rm = TRUE), 4), ' | ',
              '<strong>Median:</strong> ', round(median(comp_vals, na.rm = TRUE), 4), '<br>\n',
              '<strong>Low complexity cells (< 0.8):</strong> ', low_complexity,
              ' (', round(low_complexity / n_cells * 100, 1), '%)\n',
              '</div>\n'
            )
          }
          if ("novelty_score" %in% colnames(metadata)) {
            nov_vals <- metadata$novelty_score
            complexity_html <- paste0(complexity_html,
              '<div class="metric-box">\n',
              '<strong>Novelty Score - Mean:</strong> ', round(mean(nov_vals, na.rm = TRUE), 4), ' | ',
              '<strong>Median:</strong> ', round(median(nov_vals, na.rm = TRUE), 4), '\n',
              '</div>\n'
            )
          }
          if (!"log10GenesPerUMI" %in% colnames(metadata) && !"novelty_score" %in% colnames(metadata)) {
            complexity_html <- paste0(complexity_html, '<p><em>Complexity scores not computed.</em></p>\n')
          }

          # ---- Filtering History ----
          filter_html <- '<h2>Filtering History</h2>\n'
          if (!is.null(values$filtering_history) && length(values$filtering_history) > 0) {
            filter_html <- paste0(filter_html,
              '<div class="summary-box">\n',
              '<ul>\n'
            )
            for (entry in values$filtering_history) {
              filter_html <- paste0(filter_html, '<li>', htmltools::htmlEscape(entry), '</li>\n')
            }
            filter_html <- paste0(filter_html, '</ul>\n</div>\n')
          } else {
            filter_html <- paste0(filter_html,
              '<div class="summary-box">\n',
              '<p>No filtering steps recorded. Current dataset has <strong>',
              format(n_cells, big.mark = ","), '</strong> cells and <strong>',
              format(n_genes, big.mark = ","), '</strong> genes.</p>\n',
              '</div>\n'
            )
          }

          # ---- Footer ----
          footer_html <- paste0(
            '<div class="footer">\n',
            '<p>Generated by scRNA4 QC Module | ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"), '</p>\n',
            '</div>\n',
            '</div>\n</body>\n</html>'
          )

          incProgress(0.95, detail = "Writing report...")

          # ---- Assemble and write ----
          full_html <- paste0(
            header_html,
            info_html,
            metrics_html,
            violin_html,
            scatter_html,
            cc_html,
            batch_html,
            complexity_html,
            filter_html,
            footer_html
          )

          writeLines(full_html, file)

          incProgress(1, detail = "Complete!")
        })
      }
    )
  })
}
