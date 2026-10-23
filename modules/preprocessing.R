# Preprocessing & Integration Module (improved and robust)
# modules/preprocessing.R

suppressPackageStartupMessages({
  library(Seurat)
  library(plotly)
  library(ggplot2)
  library(dplyr)
  library(future)
  # Optional packages — loaded on demand, not required at startup
  # sctransform, harmony, batchelor, SingleCellExperiment
})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R

# -----------------------------------------------------------------------------
# NORMALIZATION (Enhanced with better error handling)
# -----------------------------------------------------------------------------

normalizationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Data Normalization & Transformation"),
          p("Choose appropriate normalization method based on your data characteristics and downstream analysis goals."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("preprocessing"),

          fluidRow(
            column(4,
              selectInput(ns("norm_method"), "Normalization Method",
                choices = c(
                  "Log Normalization (Standard)" = "LogNormalize",
                  "SCTransform (Recommended)" = "SCT",
                  "CLR (Centered Log Ratio)" = "CLR",
                  "RC (Relative Counts)" = "RC",
                  "scran (Sum Factors)" = "scran"
                ),
                selected = "LogNormalize"
              )
            ),
            column(4,
              conditionalPanel(
                condition = "input.norm_method == 'LogNormalize'",
                ns = ns,
                numericInput(ns("scale_factor"), "Scale Factor", 
                           value = 10000, min = 1000, max = 1000000, step = 1000)
              )
            ),
            column(4,
              conditionalPanel(
                condition = "input.norm_method == 'SCT'",
                ns = ns,
                checkboxInput(ns("use_glmgampoitype"), "Use glmGamPoi (faster)", TRUE)
              )
            )
          ),

          # CLR margin selector
          conditionalPanel(
            condition = "input.norm_method == 'CLR'",
            ns = ns,
            radioButtons(ns("clr_margin"), "CLR Margin",
              choices = c("Normalize across cells (margin=1)" = "1",
                          "Normalize across features (margin=2)" = "2"),
              selected = "2"),
            helpText("Margin=1: normalize each cell's feature counts. Margin=2 (default for ADT): normalize each feature across cells.")
          ),

          # Advanced options
          conditionalPanel(
            condition = "input.norm_method == 'SCT'",
            ns = ns,
            wellPanel(
              h5("SCTransform Options"),
              fluidRow(
                column(4,
                  selectInput(ns("vars_to_regress"), "Variables to Regress",
                    choices = NULL, multiple = TRUE)
                ),
                column(4,
                  numericInput(ns("n_hvg"), "Number of HVGs", 
                             value = 3000, min = 1000, max = 10000, step = 500)
                ),
                column(4,
                  checkboxInput(ns("return_only_var_genes"), "Return only variable genes", FALSE)
                )
              ),
              fluidRow(
                column(6,
                  numericInput(ns("min_cells"), "Min cells per gene", 
                             value = 5, min = 1, max = 50)
                ),
                column(6,
                  checkboxInput(ns("conserve_memory"), "Conserve memory", TRUE)
                )
              )
            )
          ),
          
          # Execution controls
          fluidRow(
            column(4,
              checkboxInput(ns("auto_scale"), "Auto-scale after normalization", TRUE)
            ),
            column(4,
              checkboxInput(ns("set_default_assay"), "Set as default assay", TRUE)
            ),
            column(4,
              div(style = "margin-top: 25px;",
                actionButton(ns("normalize"), "Normalize Data",
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          ),

          # AI normalization recommendation
          fluidRow(
            column(12,
              div(style = "margin-top: 8px;",
                actionButton(ns("ai_norm_advisor"), "AI Recommend Method",
                             icon = icon("robot"), class = "ai-btn btn-sm")
              ),
              uiOutput(ns("ai_thread_norm"))
            )
          )
        )
      )
    ),
    
    # Results display - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("norm_results_ui"))
      )
    )
  )
}

normalizationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_norm_results <- reactiveVal(FALSE)
    norm_results <- reactiveVal(NULL)

    # Dynamic UI rendering for normalization results
    output$norm_results_ui <- renderUI({
      req(isTRUE(show_norm_results()))
      ns <- session$ns
      wellPanel(
        h4("Normalization Results"),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(6,
              h5("Assay Information"),
              downloadButton(ns("dl_assay_info_table"), "CSV", class = "btn-xs btn-default"),
              div(style = "max-height: 350px; overflow-y: auto;",
                DT::DTOutput(ns("assay_info_table"))
              )
            ),
            column(6,
              h5("Normalization Quality"),
              plotlyOutput(ns("norm_quality_plot"), height = "350px")
            )
          )
        ),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(12,
              h5("HVG Selection Results"),
              plotlyOutput(ns("hvg_plot"), height = "400px")
            )
          )
        )
      )
    })
    
    # Populate regression variables
    observe({
      req(values$active_dataset)
      numeric_cols <- colnames(values$active_dataset@meta.data)[
        sapply(values$active_dataset@meta.data, is.numeric)
      ]
      # For SCT: exclude nCount_RNA/nFeature_RNA (SCT already models depth internally)
      suggested_vars <- intersect(c("percent.mt", "percent.ribo", "S.Score", "G2M.Score"), numeric_cols)
      updateSelectInput(session, "vars_to_regress", choices = numeric_cols, selected = suggested_vars)
    })
    
    observeEvent(input$normalize, {
      req(values$active_dataset)
      tryCatch({
        withProgress(message = "Normalizing data...", value = 0, {
          obj <- values$active_dataset
          method <- input$norm_method

          incProgress(0.2, detail = paste("Applying", method, "normalization..."))

          if (method == "RC") {
            showNotification("RC normalization produces relative counts without log-transformation. Many downstream methods assume log-normalized data. Consider using LogNormalize instead.", type = "warning", duration = 10)
          }

          result <- safely_normalize_data(obj, method, input)

          if (result$success) {
            obj <- result$object
            values$active_dataset <- obj
            values$normalization_done <- TRUE

            # SCTransform performs normalization + feature selection + scaling in one step
            if (method == "SCT") {
              values$features_selected <- TRUE
              values$scaling_done <- TRUE
            }
            # Auto-scale marks scaling as done
            if (isTRUE(input$auto_scale) && method != "SCT") {
              values$scaling_done <- TRUE
            }
            # Aggregate preprocessing flag: once normalisation has succeeded we
            # consider preprocessing underway; sidebar ticks it as soon as the
            # full trio is done.
            if (isTRUE(values$normalization_done) &&
                isTRUE(values$features_selected) &&
                isTRUE(values$scaling_done)) {
              values$preprocessing_done <- TRUE
            }

            norm_results(result)
            show_norm_results(TRUE)

            showNotification("Normalization completed successfully", type = "message")
          } else {
            showNotification(paste("Normalization failed:", result$error), type = "error")
            return()
          }

          incProgress(1)
        })
      }, error = function(e) {
        showNotification(paste("Normalization error:", e$message), type = "error", duration = 10)
      })
    })

    # Disable auto_scale when SCT is selected (ignoreInit prevents firing on app launch)
    observeEvent(input$norm_method, {
      if (input$norm_method == "SCT") {
        updateCheckboxInput(session, "auto_scale", value = FALSE)
        showNotification(
          "ScaleData will be skipped with SCTransform. SCT Pearson residuals are already variance-stabilized — no additional scaling needed.",
          type = "message", duration = 6
        )
      }
    }, ignoreInit = TRUE)

    # Display results
    output$assay_info_table <- DT::renderDT({
      req(norm_results())
      create_assay_info_table(values$active_dataset)
    })
    
    output$norm_quality_plot <- renderPlotly({
      req(norm_results())
      p <- create_normalization_quality_plot(values$active_dataset)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Preprocessing", figure_type = "Normalization Quality",
                       label = "Normalization Quality", source_module = "Preprocessing")
      }
      p
    })
    
    output$hvg_plot <- renderPlotly({
      req(norm_results())
      p <- create_hvg_plot(values$active_dataset)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Preprocessing", figure_type = "HVG Selection",
                       label = "HVG Selection", source_module = "Preprocessing")
      }
      p
    })

    output$dl_assay_info_table <- downloadHandler(
      filename = function() paste0("assay_info_", Sys.Date(), ".csv"),
      content = function(file) {
        req(values$active_dataset)
        tryCatch({
          tbl <- create_assay_info_table(values$active_dataset)
          write.csv(tbl$x$data, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )

    # --- AI Normalization Advisor ---
    observeEvent(input$ai_norm_advisor, {
      if (is.null(values$ai_api_key) || !nzchar(trimws(values$ai_api_key))) {
        showNotification("Please set your Claude API key in the Data Management tab first.",
                         type = "warning", duration = 8)
        return()
      }
      if (is.null(values$active_dataset)) {
        showNotification("No dataset loaded.", type = "warning")
        return()
      }
      output$ai_thread_norm <- renderUI({ ai_thread_loading("AI analyzing normalization options...") })

      tryCatch({
        obj <- values$active_dataset
        md  <- obj@meta.data

        cell_count <- ncol(obj)
        umi_vals   <- if ("nCount_RNA" %in% colnames(md)) md$nCount_RNA else NULL
        median_umi <- if (!is.null(umi_vals)) round(median(umi_vals, na.rm = TRUE), 0) else NA
        umi_var    <- if (!is.null(umi_vals)) round(var(umi_vals, na.rm = TRUE), 0) else NA

        # Detect batch-like columns
        batch_cols <- grep("batch|sample", colnames(md), ignore.case = TRUE, value = TRUE)
        n_batches  <- if (length(batch_cols) > 0) {
          max(sapply(batch_cols, function(col) length(unique(md[[col]]))))
        } else {
          1L
        }

        prompt <- build_expert_prompt(
          task_prompt = paste0(
            "TASK: Recommend the best normalization method for this scRNA-seq dataset.\n\n",
            "DATASET CHARACTERISTICS:\n",
            "- Cell count: ", cell_count, "\n",
            "- Median UMI per cell: ", median_umi, "\n",
            "- UMI variance: ", umi_var, "\n",
            "- Number of batches/samples: ", n_batches, "\n",
            "- Available assays: ", paste(names(obj@assays), collapse = ", "), "\n\n",
            "Should I use LogNormalize, SCTransform, or scran normalization?\n",
            "Explain the trade-offs for this specific dataset and provide a clear recommendation ",
            "with reasoning."
          ),
          values = values,
          include_project = TRUE,
          include_dataset = FALSE
        )

        withProgress(message = "AI analyzing normalization options...", value = 0.5, {
          result <- call_ai("integration_rec", prompt, values)
        })

        if (!result$ok) {
          output$ai_thread_norm <- renderUI({ ai_thread_error(result$error) })
          return()
        }

        output$ai_thread_norm <- renderUI({
          ai_thread_message(
            content_html = .mdtohtml(result$text),
            title = "AI Normalization Recommendation",
            model = get_model_for_task("integration_rec", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out
          )
        })

      }, error = function(e) {
        output$ai_thread_norm <- renderUI({ ai_thread_error(e$message) })
      })
    })

    # AI helper for Preprocessing module
    ai_helper_server("normalization", values, "Preprocessing", context_fn = function() {
      obj <- values$active_dataset
      if (is.null(obj)) return("")
      md <- obj@meta.data
      lines <- c(
        paste0("Cells: ", ncol(obj), ", Genes: ", nrow(obj)),
        paste0("Assays: ", paste(names(obj@assays), collapse = ", ")),
        paste0("Default assay: ", DefaultAssay(obj))
      )
      if ("nCount_RNA" %in% colnames(md)) {
        lines <- c(lines, paste0("Median UMI: ", round(median(md$nCount_RNA, na.rm = TRUE), 0)))
      }
      if (!is.null(values$normalization_done) && values$normalization_done) {
        lines <- c(lines, "Normalization: completed")
      }
      paste(lines, collapse = "\n")
    })
  })
}

# -----------------------------------------------------------------------------
# FEATURE SELECTION (Enhanced)
# -----------------------------------------------------------------------------

featureSelectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Highly Variable Feature Selection"),
          p("Identify genes with high biological variability for downstream analysis."),
          
          fluidRow(
            column(4,
              numericInput(ns("n_features"), "Number of Features", 
                         value = 2000, min = 500, max = 10000, step = 100)
            ),
            column(4,
              selectInput(ns("selection_method"), "Selection Method",
                choices = c(
                  "VST (Variance Stabilizing)" = "vst",
                  "Mean-Variance Plot" = "mean.var.plot", 
                  "Dispersion" = "dispersion"
                ),
                selected = "vst")
            ),
            column(4,
              numericInput(ns("loess_span"), "Loess Span (VST)", 
                         value = 0.3, min = 0.05, max = 1, step = 0.05)
            )
          ),
          
          # Advanced options
          fluidRow(
            column(4,
              numericInput(ns("mean_cutoff"), "Mean Expression Cutoff",
                         value = 0.1, min = 0, max = 2, step = 0.1)
            ),
            column(4,
              numericInput(ns("dispersion_cutoff"), "Dispersion Cutoff",
                         value = 1, min = 0, max = 5, step = 0.1)
            ),
            column(4,
              checkboxInput(ns("clip_max"), "Clip maximum variance", TRUE)
            )
          ),
          
          actionButton(ns("find_features"), "Find Variable Features", 
                     class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),
    
    # Results - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("hvg_results_ui"))
      )
    )
  )
}

featureSelectionServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_hvg_results <- reactiveVal(FALSE)
    hvg_info <- reactiveVal(NULL)

    # Dynamic UI rendering for feature selection results
    output$hvg_results_ui <- renderUI({
      req(isTRUE(show_hvg_results()))
      ns <- session$ns
      wellPanel(
        h4("Variable Features Results"),
        fluidRow(
          column(8,
            plotlyOutput(ns("var_features_plot"), height = "500px")
          ),
          column(4,
            h5("Selection Summary"),
            downloadButton(ns("dl_hvg_summary_table"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(ns("hvg_summary_table")),
            br(),
            downloadButton(ns("download_features"), "Download HVG List", class = "btn-success")
          )
        ),
        fluidRow(
          column(12,
            h5("Top Variable Features"),
            downloadButton(ns("dl_top_hvg_table"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(ns("top_hvg_table"))
          )
        )
      )
    })
    
    observeEvent(input$find_features, {
      req(values$active_dataset)

      # Warn if running FindVariableFeatures on SCT assay
      if (DefaultAssay(values$active_dataset) == "SCT") {
        showNotification(
          paste0("SCTransform already selected HVGs during normalization. ",
                 "Re-running FindVariableFeatures on SCT data replaces SCT-derived HVGs ",
                 "with VST-ranked ones. Switch DefaultAssay to RNA first if you need manual HVG selection."),
          type = "warning", duration = 12
        )
        return()
      }

      tryCatch({
        withProgress(message = "Finding variable features...", value = 0, {
          obj <- values$active_dataset

          incProgress(0.5, detail = "Computing feature variance...")

          result <- safely_find_variable_features(obj, input)

          if (result$success) {
            values$active_dataset <- result$object
            values$features_selected <- TRUE
            if (isTRUE(values$normalization_done) &&
                isTRUE(values$features_selected) &&
                isTRUE(values$scaling_done)) {
              values$preprocessing_done <- TRUE
            }
            hvg_info(result$hvg_info)
            show_hvg_results(TRUE)

            n_features <- length(VariableFeatures(result$object))
            showNotification(paste("Found", n_features, "variable features"), type = "message")
          } else {
            showNotification(paste("Feature selection failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        showNotification(paste("Feature selection error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$var_features_plot <- renderPlotly({
      req(hvg_info())
      p <- create_enhanced_hvg_plot(hvg_info(), values$active_dataset)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Preprocessing", figure_type = "Variable Features",
                       label = "Variable Features", source_module = "Preprocessing")
      }
      p
    })
    
    output$hvg_summary_table <- DT::renderDataTable({
      req(hvg_info())
      tbl <- create_hvg_summary_table(hvg_info(), values$active_dataset)
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "HVG Summary", source_module = "Preprocessing")
      }
      tbl
    }, options = list(dom = 't'))
    
    output$top_hvg_table <- DT::renderDataTable({
      req(values$active_dataset)
      tbl <- create_top_hvg_table(values$active_dataset)
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Top Variable Genes", source_module = "Preprocessing")
      }
      tbl
    }, options = list(pageLength = 15))
    
    output$download_features <- downloadHandler(
      filename = function() paste0("variable_features_", Sys.Date(), ".csv"),
      content = function(file) {
        hvg_list <- VariableFeatures(values$active_dataset)
        hvg_data <- hvg_info()
        if (!is.null(hvg_data)) {
          hvg_df <- hvg_data[hvg_list, , drop = FALSE]
          hvg_df$gene <- rownames(hvg_df)
          write.csv(hvg_df, file, row.names = FALSE)
        } else {
          write.csv(data.frame(Gene = hvg_list), file, row.names = FALSE)
        }
      }
    )

    output$dl_hvg_summary_table <- downloadHandler(
      filename = function() paste0("hvg_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(hvg_info())
        tryCatch({
          tbl <- create_hvg_summary_table(hvg_info(), values$active_dataset)
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )

    output$dl_top_hvg_table <- downloadHandler(
      filename = function() paste0("top_hvg_", Sys.Date(), ".csv"),
      content = function(file) {
        req(values$active_dataset)
        tryCatch({
          tbl <- create_top_hvg_table(values$active_dataset)
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )
  })
}

# -----------------------------------------------------------------------------
# SCALING (Enhanced)
# -----------------------------------------------------------------------------

scalingUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Data Scaling & Regression"),
          p("Scale data to unit variance and optionally regress out confounding factors."),
          
          fluidRow(
            column(4,
              selectInput(ns("scale_assay"), "Assay to Scale", choices = NULL)
            ),
            column(4,
              selectInput(ns("features_to_scale"), "Features to Scale",
                choices = c("All genes" = "all", "Variable features only" = "hvg"),
                selected = "hvg")
            ),
            column(4,
              selectInput(ns("vars_to_regress_scale"), "Variables to Regress",
                         choices = NULL, multiple = TRUE)
            )
          ),
          
          # Scaling parameters
          fluidRow(
            column(3,
              checkboxInput(ns("do_center"), "Center data", TRUE)
            ),
            column(3,
              checkboxInput(ns("do_scale"), "Scale to unit variance", TRUE)
            ),
            column(3,
              numericInput(ns("clip_max"), "Clip maximum value", 
                         value = 10, min = 1, max = 50)
            ),
            column(3,
              checkboxInput(ns("use_umi"), "Use UMI scaling", FALSE)
            )
          ),
          
          # Advanced options
          fluidRow(
            column(6,
              checkboxInput(ns("split_by_groups"), "Split scaling by groups", FALSE)
            ),
            column(6,
              conditionalPanel(
                condition = "input.split_by_groups",
                ns = ns,
                selectInput(ns("split_by_column"), "Split by column", choices = NULL)
              )
            )
          ),
          
          actionButton(ns("run_scaling"), "Scale Data", 
                     class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),
    
    # Results - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("scale_results_ui"))
      )
    )
  )
}

scalingServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_scale_results <- reactiveVal(FALSE)
    scale_results <- reactiveVal(NULL)

    # Dynamic UI rendering for scaling results
    output$scale_results_ui <- renderUI({
      req(isTRUE(show_scale_results()))
      ns <- session$ns
      wellPanel(
        h4("Scaling Results"),
        fluidRow(
          column(6,
            h5("Scaling Statistics"),
            downloadButton(ns("dl_scale_stats_table"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(ns("scale_stats_table"))
          ),
          column(6,
            h5("Before vs After"),
            plotlyOutput(ns("scaling_effect_plot"), height = "350px")
          )
        )
      )
    })
    
    # Update UI choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      
      # Assay choices
      assays <- names(obj@assays)
      updateSelectInput(session, "scale_assay", choices = assays, selected = DefaultAssay(obj))
      
      # Regression variables
      numeric_cols <- colnames(obj@meta.data)[sapply(obj@meta.data, is.numeric)]
      suggested_vars <- intersect(c("percent.mt", "S.Score", "G2M.Score", "CC.Difference"), numeric_cols)
      updateSelectInput(session, "vars_to_regress_scale", choices = numeric_cols, selected = suggested_vars)
      
      # Split by columns
      factor_cols <- colnames(obj@meta.data)[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "split_by_column", choices = factor_cols)
    })
    
    observeEvent(input$run_scaling, {
      req(values$active_dataset)

      # Warn against ScaleData on SCT assay — overwrites SCT-corrected scale.data
      if (DefaultAssay(values$active_dataset) == "SCT") {
        showNotification(
          paste0("Running ScaleData on SCT assay will overwrite SCTransform's regularised ",
                 "scale.data with naive z-scores, breaking downstream PCA/UMAP/clustering. ",
                 "Switch to RNA assay if you need manual scaling."),
          type = "error", duration = 12
        )
        return()
      }

      tryCatch({
        withProgress(message = "Scaling data...", value = 0, {
          obj <- values$active_dataset

          incProgress(0.3, detail = "Preparing scaling parameters...")

          result <- safely_scale_data(obj, input)

          if (result$success) {
            values$active_dataset <- result$object
            values$scaling_done <- TRUE
            if (isTRUE(values$normalization_done) &&
                isTRUE(values$features_selected) &&
                isTRUE(values$scaling_done)) {
              values$preprocessing_done <- TRUE
            }
            scale_results(result)
            show_scale_results(TRUE)

            showNotification("Data scaling completed", type = "message")
          } else {
            showNotification(paste("Scaling failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        showNotification(paste("Scaling error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$scale_stats_table <- DT::renderDataTable({
      req(scale_results())
      create_scaling_stats_table(scale_results())
    }, options = list(dom = 't'))
    
    output$scaling_effect_plot <- renderPlotly({
      req(scale_results())
      p <- create_scaling_effect_plot(scale_results())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Preprocessing", figure_type = "Scaling Effect",
                       label = "Scaling Effect", source_module = "Preprocessing")
      }
      p
    })

    output$dl_scale_stats_table <- downloadHandler(
      filename = function() paste0("scaling_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        req(scale_results())
        tryCatch({
          tbl <- create_scaling_stats_table(scale_results())
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )
  })
}

# -----------------------------------------------------------------------------
# CELL CYCLE SCORING
# -----------------------------------------------------------------------------

cellCycleUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Cell Cycle Scoring"),
          p("Score cells for S and G2M phase markers using Seurat's CellCycleScoring."),
          fluidRow(
            column(4,
              selectInput(ns("cc_species"), "Species",
                choices = c("Human" = "human", "Mouse" = "mouse"),
                selected = "human")
            ),
            column(4,
              checkboxInput(ns("cc_regress"), "Regress out cell cycle effects", FALSE)
            ),
            column(4,
              checkboxInput(ns("cc_diff_only"), "Regress only the difference (S - G2M)", TRUE)
            )
          ),
          actionButton(ns("run_cc_scoring"), "Run Cell Cycle Scoring",
            class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),
    fluidRow(
      column(12, uiOutput(ns("cc_results_ui")))
    )
  )
}

cellCycleServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$run_cc_scoring, {
      req(values$active_dataset)
      tryCatch({
        withProgress(message = "Scoring cell cycle phases...", value = 0, {
          obj <- values$active_dataset

          # Get cell cycle genes
          if (input$cc_species == "human") {
            s_genes <- Seurat::cc.genes.updated.2019$s.genes
            g2m_genes <- Seurat::cc.genes.updated.2019$g2m.genes
          } else {
            # Convert to mouse orthologs (capitalize first letter)
            s_genes <- sapply(Seurat::cc.genes.updated.2019$s.genes, function(g) {
              paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
            }, USE.NAMES = FALSE)
            g2m_genes <- sapply(Seurat::cc.genes.updated.2019$g2m.genes, function(g) {
              paste0(toupper(substr(g, 1, 1)), tolower(substr(g, 2, nchar(g))))
            }, USE.NAMES = FALSE)
          }

          # Filter to genes present in dataset
          all_genes <- rownames(obj)
          s_genes <- intersect(s_genes, all_genes)
          g2m_genes <- intersect(g2m_genes, all_genes)

          if (length(s_genes) < 3 || length(g2m_genes) < 3) {
            showNotification(
              paste0("Too few cell cycle genes found in dataset (S: ", length(s_genes),
                     ", G2M: ", length(g2m_genes), "). Need at least 3 each."),
              type = "error"
            )
            return()
          }

          incProgress(0.3, detail = "Computing phase scores...")

          obj <- CellCycleScoring(obj, s.features = s_genes, g2m.features = g2m_genes,
                                  set.ident = FALSE)

          # Always compute CC.Difference for downstream use
          obj$CC.Difference <- obj$S.Score - obj$G2M.Score

          if (input$cc_regress) {
            incProgress(0.6, detail = "Regressing cell cycle effects...")
            if (input$cc_diff_only) {
              obj <- ScaleData(obj, vars.to.regress = "CC.Difference", verbose = FALSE)
            } else {
              obj <- ScaleData(obj, vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
            }
          }

          values$active_dataset <- obj
          incProgress(1, detail = "Done!")

          phase_counts <- table(obj$Phase)
          showNotification(
            paste0("Cell cycle scoring complete: ",
                   "G1=", phase_counts["G1"] %||% 0, ", S=", phase_counts["S"] %||% 0,
                   ", G2M=", phase_counts["G2M"] %||% 0),
            type = "message"
          )

          output$cc_results_ui <- renderUI({
            ns <- session$ns
            wellPanel(
              h4("Cell Cycle Results"),
              fluidRow(
                column(6, plotlyOutput(ns("cc_phase_pie"), height = "300px")),
                column(6, plotlyOutput(ns("cc_score_scatter"), height = "300px"))
              )
            )
          })
        })
      }, error = function(e) {
        showNotification(paste("Cell cycle scoring failed:", e$message), type = "error")
      })
    })

    output$cc_phase_pie <- renderPlotly({
      req(values$active_dataset, "Phase" %in% colnames(values$active_dataset@meta.data))
      phase_counts <- as.data.frame(table(values$active_dataset$Phase))
      colnames(phase_counts) <- c("Phase", "Count")
      plot_ly(phase_counts, labels = ~Phase, values = ~Count, type = "pie",
              marker = list(colors = c("#66c2a5", "#fc8d62", "#8da0cb"))) %>%
        layout(title = "Cell Cycle Phase Distribution")
    })

    output$cc_score_scatter <- renderPlotly({
      req(values$active_dataset,
          "S.Score" %in% colnames(values$active_dataset@meta.data),
          "G2M.Score" %in% colnames(values$active_dataset@meta.data))
      md <- values$active_dataset@meta.data
      plot_ly(md, x = ~S.Score, y = ~G2M.Score, color = ~Phase,
              type = "scatter", mode = "markers",
              marker = list(size = 3, opacity = 0.5),
              colors = c("#66c2a5", "#fc8d62", "#8da0cb")) %>%
        layout(title = "S vs G2M Scores",
               xaxis = list(title = "S Score"),
               yaxis = list(title = "G2M Score"))
    })
  })
}

# -----------------------------------------------------------------------------
# INTEGRATION (Enhanced with multiple methods)
# -----------------------------------------------------------------------------

integrationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Dataset Integration & Batch Correction"),
          p("Integrate multiple datasets or correct for batch effects within a single dataset."),
          
          # Integration strategy
          fluidRow(
            column(6,
              selectInput(ns("integration_strategy"), "Integration Strategy",
                choices = c(
                  "Multi-dataset Integration" = "multi",
                  "Batch Correction (Single Dataset)" = "batch"
                ),
                selected = "multi")
            ),
            column(6,
              selectInput(ns("integration_method"), "Method",
                choices = c(
                  "Seurat SCT Integration" = "sct",
                  "Seurat CCA Integration" = "cca",
                  "Seurat RPCA Integration" = "rpca",
                  "Harmony" = "harmony",
                  "FastMNN" = "fastmnn",
                  "Combat" = "combat"
                ),
                selected = "sct")
            )
          ),
          
          # Dataset/batch selection
          conditionalPanel(
            condition = "input.integration_strategy == 'multi'",
            ns = ns,
            fluidRow(
              column(12,
                h5("Select Datasets to Integrate"),
                uiOutput(ns("dataset_selection_ui"))
              )
            )
          ),
          
          conditionalPanel(
            condition = "input.integration_strategy == 'batch'",
            ns = ns,
            fluidRow(
              column(6,
                selectInput(ns("batch_column"), "Batch Column", choices = NULL)
              ),
              column(6,
                checkboxInput(ns("add_batch_prefix"), "Add batch prefix to cell names", TRUE)
              )
            )
          ),
          
          # Method-specific parameters
          fluidRow(
            column(4,
              numericInput(ns("n_features_integrate"), "Integration Features",
                         value = 3000, min = 1000, max = 8000, step = 500)
            ),
            column(4,
              numericInput(ns("n_dims_integrate"), "Dimensions to Use",
                         value = 30, min = 10, max = 100, step = 5)
            ),
            column(4,
              numericInput(ns("k_anchor"), "k.anchor (neighbors)",
                         value = 5, min = 3, max = 20)
            )
          ),

          # Harmony-specific parameters
          conditionalPanel(
            condition = "input.integration_method == 'harmony'",
            ns = ns,
            fluidRow(
              column(6,
                sliderInput(ns("harmony_theta"), "Harmony diversity penalty (theta)",
                            min = 0, max = 4, value = 2, step = 0.5),
                helpText("Lower theta (0.5-1): preserves subtle biology. Higher theta (2-4): stronger batch removal. Default 2.")
              )
            )
          ),
          
          # Post-integration options
          fluidRow(
            column(4,
              checkboxInput(ns("run_pca_after"), "Run PCA after integration", TRUE)
            ),
            column(4,
              checkboxInput(ns("run_umap_after"), "Run UMAP after integration", TRUE)
            ),
            column(4,
              checkboxInput(ns("find_clusters_after"), "Find clusters after integration", TRUE)
            )
          ),
          
          actionButton(ns("run_integration"), "Run Integration",
                     class = "btn-success btn-lg", width = "100%"),

          # AI integration strategy recommendation
          div(style = "margin-top: 8px;",
            actionButton(ns("ai_integration_advisor"), "AI Recommend Strategy",
                         icon = icon("robot"), class = "ai-btn btn-sm")
          ),
          uiOutput(ns("ai_thread_integration"))
        )
      )
    ),

    # Results - using uiOutput for reliable rendering
    fluidRow(
      column(12,
        uiOutput(ns("integration_results_ui"))
      )
    )
  )
}

integrationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_integration_results <- reactiveVal(FALSE)
    integration_results <- reactiveVal(NULL)

    # Dynamic UI rendering for integration results
    output$integration_results_ui <- renderUI({
      req(isTRUE(show_integration_results()))
      ns <- session$ns
      wellPanel(
        h4("Integration Results"),
        fluidRow(
          column(6,
            h5("Integration Summary"),
            downloadButton(ns("dl_integration_summary"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(ns("integration_summary"))
          ),
          column(6,
            h5("Integration Quality Metrics"),
            downloadButton(ns("dl_integration_metrics"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(ns("integration_metrics"))
          )
        ),
        fluidRow(
          column(12,
            h5("Before vs After Integration"),
            plotlyOutput(ns("integration_comparison"), height = "500px")
          )
        )
      )
    })
    
    # Update UI based on available data
    observe({
      req(values$active_dataset)
      
      # Batch columns for single dataset batch correction
      meta_cols <- colnames(values$active_dataset@meta.data)
      factor_cols <- meta_cols[sapply(values$active_dataset@meta.data, function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "batch_column", choices = factor_cols)
    })
    
    # Dataset selection UI for multi-dataset integration
    output$dataset_selection_ui <- renderUI({
      if (length(values$seurat_list) < 2) {
        return(p("❌ At least 2 datasets required for multi-dataset integration."))
      }
      
      dataset_names <- names(values$seurat_list)
      checkboxGroupInput(session$ns("selected_datasets"), "Datasets to integrate:",
                        choices = dataset_names, selected = dataset_names)
    })
    
    observeEvent(input$run_integration, {
      req(values$active_dataset)
      tryCatch({
        withProgress(message = "Running integration...", value = 0, {

          if (input$integration_strategy == "multi" && length(values$seurat_list) < 2) {
            showNotification("At least 2 datasets required for integration", type = "error")
            return()
          }

          incProgress(0.2, detail = "Preparing integration...")

          result <- safely_run_integration(values, input)

          if (result$success) {
            integrated <- result$object

            # Reset DefaultAssay to RNA for downstream DE analysis
            if ("RNA" %in% names(integrated@assays)) {
              DefaultAssay(integrated) <- "RNA"
              integrated <- NormalizeData(integrated, verbose = FALSE)
              showNotification(
                "DefaultAssay reset to 'RNA' for downstream analysis. The integrated assay is used only for clustering/visualization.",
                type = "message", duration = 8
              )
            }

            values$active_dataset <- integrated
            result$object <- integrated
            integration_results(result)
            show_integration_results(TRUE)

            values$preprocessing_done <- TRUE
            showNotification("Integration completed successfully", type = "message")
          } else {
            showNotification(paste("Integration failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        showNotification(paste("Integration error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$integration_summary <- DT::renderDataTable({
      req(integration_results())
      create_integration_summary_table(integration_results())
    }, options = list(dom = 't'))
    
    output$integration_metrics <- DT::renderDataTable({
      req(integration_results())
      create_integration_metrics_table(integration_results())
    }, options = list(dom = 't'))
    
    output$integration_comparison <- renderPlotly({
      req(integration_results())
      p <- create_integration_comparison_plot(integration_results())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Preprocessing", figure_type = "Integration Comparison",
                       label = "Integration Comparison", source_module = "Preprocessing")
      }
      p
    })

    output$dl_integration_summary <- downloadHandler(
      filename = function() paste0("integration_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(integration_results())
        tryCatch({
          tbl <- create_integration_summary_table(integration_results())
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )

    output$dl_integration_metrics <- downloadHandler(
      filename = function() paste0("integration_metrics_", Sys.Date(), ".csv"),
      content = function(file) {
        req(integration_results())
        tryCatch({
          tbl <- create_integration_metrics_table(integration_results())
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )

    # --- AI Integration Strategy Advisor ---
    observeEvent(input$ai_integration_advisor, {
      if (is.null(values$ai_api_key) || !nzchar(trimws(values$ai_api_key))) {
        showNotification("Please set your Claude API key in the Data Management tab first.",
                         type = "warning", duration = 8)
        return()
      }
      if (is.null(values$active_dataset)) {
        showNotification("No dataset loaded.", type = "warning")
        return()
      }
      output$ai_thread_integration <- renderUI({ ai_thread_loading("AI analyzing integration strategies...") })

      tryCatch({
        obj <- values$active_dataset
        md  <- obj@meta.data

        cell_count <- ncol(obj)
        assays     <- names(obj@assays)

        # Detect batch columns
        batch_cols <- grep("batch|sample", colnames(md), ignore.case = TRUE, value = TRUE)
        n_batches  <- 1L
        batch_sizes_str <- "single batch"
        if (length(batch_cols) > 0) {
          # Use the column with the most groups
          best_col <- batch_cols[which.max(sapply(batch_cols, function(col) length(unique(md[[col]]))))]
          batch_tbl <- table(md[[best_col]])
          n_batches <- length(batch_tbl)
          batch_sizes_str <- paste(paste0(names(batch_tbl), ": ", batch_tbl), collapse = ", ")
        }

        prompt <- build_expert_prompt(
          task_prompt = paste0(
            "TASK: Recommend the best integration/batch correction strategy for this scRNA-seq dataset.\n\n",
            "DATASET CHARACTERISTICS:\n",
            "- Total cell count: ", cell_count, "\n",
            "- Number of batches: ", n_batches, "\n",
            "- Batch sizes: ", batch_sizes_str, "\n",
            "- Available assays: ", paste(assays, collapse = ", "), "\n",
            "- Batch columns detected: ", paste(batch_cols, collapse = ", "), "\n\n",
            "Which method should I use: Harmony, Seurat CCA, Seurat RPCA, or FastMNN?\n",
            "Consider the dataset size, number of batches, and computational requirements.\n",
            "Provide a clear recommendation with reasoning and any important caveats."
          ),
          values = values,
          include_project = TRUE,
          include_dataset = FALSE
        )

        withProgress(message = "AI analyzing integration strategies...", value = 0.5, {
          result <- call_ai("integration_rec", prompt, values)
        })

        if (!result$ok) {
          output$ai_thread_integration <- renderUI({ ai_thread_error(result$error) })
          return()
        }

        output$ai_thread_integration <- renderUI({
          ai_thread_message(
            content_html = .mdtohtml(result$text),
            title = "AI Integration Strategy Recommendation",
            model = get_model_for_task("integration_rec", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out
          )
        })

      }, error = function(e) {
        output$ai_thread_integration <- renderUI({ ai_thread_error(e$message) })
      })
    })
  })
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

safely_normalize_data <- function(seurat_obj, method, params) {
  tryCatch({
    # Force garbage collection before heavy operations
    gc(verbose = FALSE)

    # Seurat v5: join split layers BEFORE normalization to avoid "layer not found" errors
    for (a in names(seurat_obj@assays)) {
      if (inherits(seurat_obj[[a]], "Assay5")) {
        lyrs <- Layers(seurat_obj[[a]])
        if (any(grepl("\\.", lyrs))) {
          tryCatch(
            seurat_obj[[a]] <- JoinLayers(seurat_obj[[a]]),
            error = function(e) NULL
          )
        }
      }
    }

    if (method == "SCT") {
      # SCTransform with robust error handling
      vars_to_regress <- if (length(params$vars_to_regress) > 0) params$vars_to_regress else NULL

      # Remove nCount_RNA and nFeature_RNA from SCT regression variables —
      # SCTransform already models sequencing depth via its NB regression.
      # Regressing them out again causes catastrophic over-correction.
      if (!is.null(vars_to_regress)) {
        depth_vars <- c("nCount_RNA", "nFeature_RNA")
        removed <- intersect(vars_to_regress, depth_vars)
        if (length(removed) > 0) {
          vars_to_regress <- setdiff(vars_to_regress, depth_vars)
          if (length(vars_to_regress) == 0) vars_to_regress <- NULL
          warning("Removed ", paste(removed, collapse = ", "),
                  " from SCTransform regression — SCT already models sequencing depth.")
        }
      }

      # conserveMemory + autoScale guard
      use_conserve <- isTRUE(params$conserve_memory)
      if (use_conserve && isTRUE(params$auto_scale)) {
        warning("conserveMemory=TRUE: scale.data not retained by SCTransform. ",
                "Skipping ScaleData to avoid memory crash.")
      }

      # Check if glmGamPoi is available and requested
      use_glmgampoi <- isTRUE(params$use_glmgampoitype) && requireNamespace("glmGamPoi", quietly = TRUE)

      # Always run SCTransform sequentially to avoid memory duplication
      # Multisession copies the full Seurat object to each worker which can
      # exhaust system memory and crash the machine
      old_plan <- future::plan()
      future::plan("sequential")
      on.exit(future::plan(old_plan), add = TRUE)

      # Auto-enable memory conservation for large datasets
      sct_conserve_mem <- ncol(seurat_obj) > 50000L
      use_conserve <- use_conserve || sct_conserve_mem

      normalized_obj <- SCTransform(
        seurat_obj,
        vars.to.regress = vars_to_regress,
        variable.features.n = params$n_hvg %||% 3000L,
        return.only.var.genes = isTRUE(params$return_only_var_genes),
        min_cells = params$min_cells %||% 5L,
        method = if (use_glmgampoi) "glmGamPoi" else "poisson",
        conserve.memory = use_conserve,
        verbose = FALSE
      )

      if (isTRUE(params$set_default_assay)) {
        DefaultAssay(normalized_obj) <- "SCT"
      }

    } else if (method == "scran") {
      # scran sum-factor normalization via SingleCellExperiment
      if (!requireNamespace("scran", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("scran normalization requires the 'scran' and 'SingleCellExperiment' packages. ",
             "Install them with: BiocManager::install(c('scran', 'SingleCellExperiment'))")
      }
      sce <- as.SingleCellExperiment(seurat_obj)
      clusters <- scran::quickCluster(sce)
      sce <- scran::computeSumFactors(sce, clusters = clusters)
      sce <- scater::logNormCounts(sce)
      normalized_obj <- as.Seurat(sce, data = "logcounts")

      # Auto-scale if requested
      if (isTRUE(params$auto_scale)) {
        normalized_obj <- ScaleData(normalized_obj, verbose = FALSE)
      }

    } else {
      # Standard normalization methods (LogNormalize, CLR, RC)
      norm_args <- list(
        object = seurat_obj,
        normalization.method = method,
        scale.factor = if (method == "LogNormalize") params$scale_factor else 10000,
        verbose = FALSE
      )

      # Pass CLR margin if applicable
      if (method == "CLR" && !is.null(params$clr_margin)) {
        norm_args$margin <- as.integer(params$clr_margin)
      }

      normalized_obj <- do.call(NormalizeData, norm_args)

      # Auto-scale if requested and not SCT
      if (isTRUE(params$auto_scale) && method != "SCT") {
        normalized_obj <- ScaleData(normalized_obj, verbose = FALSE)
      } else if (isTRUE(params$auto_scale) && method == "SCT") {
        message("ScaleData skipped: SCT Pearson residuals are already variance-stabilized.")
      }
    }

    # Seurat v5: join any split layers so downstream steps work cleanly
    for (a in names(normalized_obj@assays)) {
      if (inherits(normalized_obj[[a]], "Assay5")) {
        lyrs <- Layers(normalized_obj[[a]])
        if (any(grepl("\\.", lyrs))) {
          tryCatch({
            normalized_obj[[a]] <- JoinLayers(normalized_obj[[a]])
          }, error = function(e) NULL)
        }
      }
    }

    # Free memory after normalization
    gc(verbose = FALSE)

    list(
      success = TRUE,
      object = normalized_obj,
      method = method,
      parameters = params,
      error = NULL
    )

  }, error = function(e) {
    gc(verbose = FALSE)
    list(
      success = FALSE,
      object = NULL,
      method = method,
      parameters = params,
      error = e$message
    )
  })
}

safely_find_variable_features <- function(seurat_obj, params) {
  tryCatch({
    hvg_obj <- FindVariableFeatures(
      seurat_obj,
      selection.method = params$selection_method,
      nfeatures = params$n_features,
      loess.span = params$loess_span,
      clip.max = if (params$clip_max) "auto" else Inf,
      mean.cutoff = c(params$mean_cutoff, Inf),
      dispersion.cutoff = c(params$dispersion_cutoff, Inf),
      verbose = FALSE
    )
    
    hvg_info <- HVFInfo(hvg_obj)
    
    list(
      success = TRUE,
      object = hvg_obj,
      hvg_info = hvg_info,
      n_features = length(VariableFeatures(hvg_obj)),
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      hvg_info = NULL,
      n_features = 0,
      error = e$message
    )
  })
}

safely_scale_data <- function(seurat_obj, params) {
  tryCatch({
    # Determine features to scale
    if (params$features_to_scale == "hvg") {
      features <- VariableFeatures(seurat_obj)
      if (length(features) == 0) {
        features <- rownames(seurat_obj)
      }
    } else {
      features <- rownames(seurat_obj)
    }
    
    # Set assay
    current_assay <- DefaultAssay(seurat_obj)
    if (params$scale_assay != current_assay) {
      DefaultAssay(seurat_obj) <- params$scale_assay
    }
    
    # Prepare regression variables
    vars_to_regress <- if (length(params$vars_to_regress_scale) > 0) params$vars_to_regress_scale else NULL
    
    # Memory budget guard before ScaleData (creates dense features x cells matrix)
    est_elements <- as.numeric(length(features)) * as.numeric(ncol(seurat_obj))
    if (est_elements > 1.2e9) {
      stop("ScaleData would create a ", round(est_elements / 1e9, 1),
           "B-element dense matrix (estimated ", round(est_elements * 8 / 1e9, 1),
           " GB RAM). Use 'Variable features only' instead of 'All genes', ",
           "or switch to SCTransform which avoids storing scale.data entirely.")
    }

    # Scale data
    scaled_obj <- ScaleData(
      seurat_obj,
      features = features,
      vars.to.regress = vars_to_regress,
      do.center = params$do_center,
      do.scale = params$do_scale,
      clip.max = params$clip_max,
      verbose = FALSE
    )
    
    # Calculate scaling statistics
    scaled_data <- safe_get_assay_data(scaled_obj, "scale.data")
    
    list(
      success = TRUE,
      object = scaled_obj,
      n_features_scaled = nrow(scaled_data),
      features_scaled = rownames(scaled_data),
      mean_variance = mean(apply(scaled_data, 1, var, na.rm = TRUE), na.rm = TRUE),
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      n_features_scaled = 0,
      features_scaled = NULL,
      mean_variance = NA,
      error = e$message
    )
  })
}

safely_run_integration <- function(values, params) {
  tryCatch({
    method <- params$integration_method
    strategy <- params$integration_strategy
    
    if (strategy == "multi") {
      # Multi-dataset integration
      selected_datasets <- params$selected_datasets
      if (is.null(selected_datasets) || length(selected_datasets) < 2) {
        stop("Select at least 2 datasets for integration")
      }
      
      object_list <- values$seurat_list[selected_datasets]
      result <- integrate_multiple_datasets(object_list, method, params)
      
    } else {
      # Single dataset batch correction
      batch_col <- params$batch_column
      if (is.null(batch_col) || !batch_col %in% colnames(values$active_dataset@meta.data)) {
        stop("Valid batch column required for batch correction")
      }
      
      result <- correct_batch_effects(values$active_dataset, method, batch_col, params)
    }
    
    result
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      method = params$integration_method,
      strategy = params$integration_strategy,
      error = e$message
    )
  })
}

integrate_multiple_datasets <- function(object_list, method, params) {
  # Seurat integration passes full Seurat objects to future workers — raise limit
  old_limit <- getOption("future.globals.maxSize")
  options(future.globals.maxSize = 16 * 1024^3)
  on.exit(options(future.globals.maxSize = old_limit), add = TRUE)

  # Detect Seurat version for workflow compatibility
  seurat_major <- as.integer(strsplit(as.character(packageVersion("Seurat")), "\\.")[[1]][1])

  # Seurat v5: use IntegrateLayers for CCA/RPCA/SCT instead of deprecated v4 workflow
  if (seurat_major >= 5 && method %in% c("cca", "rpca", "sct")) {
    # Merge all objects first (v5 works with layers)
    merged <- merge(object_list[[1]], y = object_list[-1])
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- FindVariableFeatures(merged, nfeatures = params$n_features_integrate, verbose = FALSE)
    merged <- ScaleData(merged, verbose = FALSE)
    merged <- RunPCA(merged, npcs = max(50, params$n_dims_integrate), verbose = FALSE)

    integration_method <- tryCatch({
      switch(method,
        "cca" = CCAIntegration,
        "rpca" = RPCAIntegration,
        "sct"  = CCAIntegration  # SCT uses CCA in v5 IntegrateLayers
      )
    }, error = function(e) NULL)

    if (!is.null(integration_method)) {
      integrated <- tryCatch({
        IntegrateLayers(
          merged,
          method = integration_method,
          orig.reduction = "pca",
          new.reduction = "integrated",
          verbose = FALSE
        )
      }, error = function(e) {
        warning("IntegrateLayers failed (", e$message, "), falling back to v4 workflow.")
        NULL
      })
      if (!is.null(integrated)) {
        integrated <- JoinLayers(integrated)
        # Post-integration processing
        if (params$run_pca_after) {
          integrated <- RunPCA(integrated, npcs = max(50, params$n_dims_integrate), verbose = FALSE)
        }
        if (params$run_umap_after) {
          red <- if ("integrated" %in% names(integrated@reductions)) "integrated" else "pca"
          integrated <- RunUMAP(integrated, reduction = red, dims = 1:params$n_dims_integrate, umap.method = "uwot", verbose = FALSE)
        }
        return(integrated)
      }
    }
    # If IntegrateLayers is not available (old Seurat 5 build), fall through to v4 workflow
  }

  if (method == "sct") {
    # SCT Integration
    object_list <- lapply(object_list, function(obj) {
      SCTransform(obj, variable.features.n = params$n_features_integrate, verbose = FALSE)
    })
    
    features <- SelectIntegrationFeatures(object_list, nfeatures = params$n_features_integrate)
    object_list <- PrepSCTIntegration(object_list, anchor.features = features, verbose = FALSE)
    
    anchors <- FindIntegrationAnchors(
      object.list = object_list,
      normalization.method = "SCT",
      anchor.features = features,
      dims = 1:params$n_dims_integrate,
      k.anchor = params$k_anchor,
      verbose = FALSE
    )
    
    integrated <- IntegrateData(
      anchorset = anchors,
      normalization.method = "SCT",
      dims = 1:params$n_dims_integrate,
      verbose = FALSE
    )
    
    DefaultAssay(integrated) <- "integrated"
    
  } else if (method == "cca") {
    # CCA Integration
    object_list <- lapply(object_list, function(obj) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, nfeatures = params$n_features_integrate, verbose = FALSE)
      obj
    })
    
    features <- SelectIntegrationFeatures(object_list, nfeatures = params$n_features_integrate)
    object_list <- lapply(object_list, function(obj) {
      ScaleData(obj, features = features, verbose = FALSE)
    })
    
    anchors <- FindIntegrationAnchors(
      object.list = object_list,
      anchor.features = features,
      dims = 1:params$n_dims_integrate,
      k.anchor = params$k_anchor,
      verbose = FALSE
    )
    
    integrated <- IntegrateData(
      anchorset = anchors,
      dims = 1:params$n_dims_integrate,
      verbose = FALSE
    )
    
    DefaultAssay(integrated) <- "integrated"
    
  } else if (method == "rpca") {
    # RPCA Integration
    object_list <- lapply(object_list, function(obj) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, nfeatures = params$n_features_integrate, verbose = FALSE)
      obj
    })
    
    features <- SelectIntegrationFeatures(object_list, nfeatures = params$n_features_integrate)
    object_list <- lapply(object_list, function(obj) {
      obj <- ScaleData(obj, features = features, verbose = FALSE)
      obj <- RunPCA(obj, features = features, npcs = max(50, params$n_dims_integrate), verbose = FALSE)
      obj
    })
    
    anchors <- FindIntegrationAnchors(
      object.list = object_list,
      anchor.features = features,
      reduction = "rpca",
      dims = 1:params$n_dims_integrate,
      k.anchor = params$k_anchor,
      verbose = FALSE
    )
    
    integrated <- IntegrateData(
      anchorset = anchors,
      dims = 1:params$n_dims_integrate,
      verbose = FALSE
    )
    
    DefaultAssay(integrated) <- "integrated"
    
  } else if (method == "fastmnn") {
    # FastMNN Integration (batchelor)
    if (!requireNamespace("batchelor", quietly = TRUE)) {
      stop("batchelor package required for FastMNN. Install: BiocManager::install('batchelor')")
    }
    # Normalize and find variable features
    object_list <- lapply(object_list, function(obj) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, nfeatures = params$n_features_integrate, verbose = FALSE)
      obj
    })
    features <- SelectIntegrationFeatures(object_list, nfeatures = params$n_features_integrate)
    # Convert to SCE for batchelor
    sce_list <- lapply(object_list, function(obj) {
      as.SingleCellExperiment(obj)
    })
    mnn_out <- do.call(batchelor::fastMNN, c(sce_list, list(
      subset.row = features,
      d = params$n_dims_integrate,
      BSPARAM = BiocSingular::RandomParam()
    )))
    # Convert back to Seurat
    integrated <- as.Seurat(mnn_out, counts = "counts", data = "logcounts")
    # Store corrected low-dim embedding with explicit cell alignment
    mnn_embed <- SingleCellExperiment::reducedDim(mnn_out, "corrected")
    # Ensure rownames match the integrated Seurat object's cell order
    if (!identical(rownames(mnn_embed), colnames(integrated))) {
      common_cells <- intersect(rownames(mnn_embed), colnames(integrated))
      mnn_embed <- mnn_embed[match(colnames(integrated), rownames(mnn_embed)), , drop = FALSE]
      rownames(mnn_embed) <- colnames(integrated)
    }
    integrated[["mnn"]] <- CreateDimReducObject(
      embeddings = mnn_embed,
      key = "MNN_", assay = DefaultAssay(integrated)
    )

  } else if (method == "harmony") {
    # Harmony for multi-dataset: merge first, then Harmony
    merged <- merge(object_list[[1]], y = object_list[-1], project = "harmony_integrated")
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- FindVariableFeatures(merged, nfeatures = params$n_features_integrate, verbose = FALSE)
    merged <- ScaleData(merged, verbose = FALSE)
    merged <- RunPCA(merged, npcs = max(50, params$n_dims_integrate), verbose = FALSE)
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("harmony package required. Install: install.packages('harmony')")
    }
    # Check for small batches
    batch_sizes <- table(merged@meta.data[["orig.ident"]])
    small_batches <- names(batch_sizes)[batch_sizes < 20]
    if (length(small_batches) > 0) {
      warning("Batches with < 20 cells: ", paste(small_batches, collapse = ", "),
              ". Harmony may produce NaN embeddings.")
    }
    harmony_theta <- params$harmony_theta %||% 2
    integrated <- harmony::RunHarmony(merged, group.by.vars = "orig.ident",
                                       theta = harmony_theta,
                                       dims.use = 1:params$n_dims_integrate, verbose = FALSE)

  } else {
    stop("Unsupported integration method for multi-dataset: ", method)
  }

  # Post-integration processing
  if (params$run_pca_after) {
    integrated <- RunPCA(integrated, npcs = max(50, params$n_dims_integrate), verbose = FALSE)
  }

  # Determine the correct corrected reduction for downstream steps
  corr_red <- if ("harmony" %in% names(integrated@reductions)) "harmony"
              else if ("integrated.dr" %in% names(integrated@reductions)) "integrated.dr"
              else if ("mnn" %in% names(integrated@reductions)) "mnn"
              else "pca"

  if (params$run_umap_after) {
    integrated <- RunUMAP(integrated, reduction = corr_red, dims = 1:params$n_dims_integrate, umap.method = "uwot", verbose = FALSE)
  }

  if (params$find_clusters_after) {
    cluster_dims <- min(
      max(2L, as.integer(params$n_dims_integrate %||% 30L)),
      ncol(Embeddings(integrated, corr_red))
    )
    cluster_k <- min(max(20L, floor(sqrt(ncol(integrated)) / 2)), ncol(integrated) - 1L)

    if (exists("run_standard_graph_clustering", mode = "function")) {
      integrated <- run_standard_graph_clustering(
        integrated,
        reduction = corr_red,
        dims = cluster_dims,
        k_param = cluster_k,
        resolution = 0.5,
        algorithm = "louvain",
        cluster_name = "seurat_clusters",
        graph_prefix = "integration_cluster",
        compute_quality = FALSE,
        random_seed = 42L
      )$object
    } else {
      integrated <- FindNeighbors(
        integrated,
        reduction = corr_red,
        dims = seq_len(cluster_dims),
        k.param = cluster_k,
        graph.name = c("integration_cluster_nn", "integration_cluster_snn"),
        verbose = FALSE
      )
      integrated <- FindClusters(
        integrated,
        graph.name = "integration_cluster_snn",
        resolution = 0.5,
        random.seed = 42L,
        verbose = FALSE
      )
    }
  }
  
  list(
    success = TRUE,
    object = integrated,
    method = method,
    strategy = "multi",
    n_datasets = length(object_list),
    error = NULL
  )
}

correct_batch_effects <- function(seurat_obj, method, batch_col, params) {
  if (method == "harmony") {
    # Harmony batch correction
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("harmony package required for Harmony integration")
    }

    # Check for small batches that can produce NaN embeddings
    batch_sizes <- table(seurat_obj@meta.data[[batch_col]])
    small_batches <- names(batch_sizes)[batch_sizes < 20]
    if (length(small_batches) > 0) {
      warning("Batches with < 20 cells detected: ", paste(small_batches, collapse = ", "),
              ". Harmony may produce NaN embeddings. Consider removing or merging these batches.")
    }

    # Ensure PCA exists
    if (!"pca" %in% names(seurat_obj@reductions)) {
      seurat_obj <- RunPCA(seurat_obj, npcs = max(50, params$n_dims_integrate), verbose = FALSE)
    }

    harmony_theta <- params$harmony_theta %||% 2
    corrected <- RunHarmony(
      seurat_obj,
      group.by.vars = batch_col,
      theta = harmony_theta,
      dims.use = 1:params$n_dims_integrate,
      verbose = FALSE
    )
    
    if (params$run_umap_after) {
      corrected <- RunUMAP(corrected, reduction = "harmony", dims = 1:params$n_dims_integrate, umap.method = "uwot", verbose = FALSE)
    }
    
  } else if (method == "combat") {
    # ComBat batch correction
    if (!requireNamespace("sva", quietly = TRUE)) {
      stop("sva package required for ComBat")
    }
    
    expr_data <- as.matrix(safe_get_assay_data(seurat_obj, "data"))
    batch_vector <- seurat_obj@meta.data[[batch_col]]
    
    corrected_data <- sva::ComBat(dat = expr_data, batch = batch_vector, verbose = FALSE)
    
    seurat_obj[["RNA_combat"]] <- CreateAssayObject(data = corrected_data)
    DefaultAssay(seurat_obj) <- "RNA_combat"
    corrected <- seurat_obj
    
  } else if (method == "fastmnn") {
    # FastMNN for single-dataset batch correction
    if (!requireNamespace("batchelor", quietly = TRUE)) {
      stop("batchelor package required for FastMNN. Install: BiocManager::install('batchelor')")
    }
    sce <- as.SingleCellExperiment(seurat_obj)
    mnn_out <- batchelor::fastMNN(sce, batch = seurat_obj@meta.data[[batch_col]],
                                   d = params$n_dims_integrate,
                                   BSPARAM = BiocSingular::RandomParam())
    corrected <- seurat_obj
    corrected[["mnn"]] <- CreateDimReducObject(
      embeddings = SingleCellExperiment::reducedDim(mnn_out, "corrected"),
      key = "MNN_", assay = DefaultAssay(corrected)
    )
    if (params$run_umap_after) {
      corrected <- RunUMAP(corrected, reduction = "mnn", dims = 1:params$n_dims_integrate, umap.method = "uwot", verbose = FALSE)
    }

  } else {
    stop("Unsupported batch correction method: ", method)
  }

  list(
    success = TRUE,
    object = corrected,
    method = method,
    strategy = "batch",
    batch_column = batch_col,
    error = NULL
  )
}

# Display helper functions
create_assay_info_table <- function(seurat_obj) {
  assay_info <- lapply(names(seurat_obj@assays), function(assay_name) {
    assay <- seurat_obj@assays[[assay_name]]
    data.frame(
      Assay = assay_name,
      Genes = nrow(assay),
      Cells = ncol(assay),
      Slots = paste(tryCatch(slotNames(assay), error = function(e) {
        if (isGeneric("Layers") && !is.null(tryCatch(Layers(assay), error = function(e2) NULL))) {
          Layers(assay)
        } else {
          "data"
        }
      }), collapse = ", "),
      Default = ifelse(assay_name == DefaultAssay(seurat_obj), "✓", "")
    )
  })
  
  DT::datatable(do.call(rbind, assay_info), options = list(dom = 't'))
}

create_normalization_quality_plot <- function(seurat_obj) {
  tryCatch({
    set.seed(42)
    expr_data <- safe_get_assay_data(seurat_obj, "data")
    sampled_genes <- sample(rownames(expr_data), min(100, nrow(expr_data)))
    sampled_cells <- sample(ncol(expr_data), min(500, ncol(expr_data)))
    sampled_data <- as.matrix(expr_data[sampled_genes, sampled_cells])

    # Calculate CV for sampled genes
    gene_means <- rowMeans(sampled_data)
    gene_vars <- apply(sampled_data, 1, var)
    gene_cv <- sqrt(gene_vars) / pmax(gene_means, 1e-10)
    gene_cv[is.infinite(gene_cv) | is.na(gene_cv)] <- 0

    plot_data <- data.frame(
      Mean = gene_means,
      CV = gene_cv,
      Gene = rownames(sampled_data)
    )

    plot_ly(
      data = plot_data,
      x = ~Mean,
      y = ~CV,
      type = "scattergl",
      mode = "markers",
      marker = list(size = 4, opacity = 0.7),
      text = ~Gene,
      hovertemplate = "Gene: %{text}<br>Mean: %{x:.3f}<br>CV: %{y:.3f}<extra></extra>"
    ) %>%
      layout(
        title = "Normalization Quality: Mean vs CV",
        xaxis = list(title = "Mean Expression", type = "log"),
        yaxis = list(title = "Coefficient of Variation")
      )
  }, error = function(e) {
    plotly_empty() %>% layout(title = paste("Quality plot unavailable:", e$message))
  })
}

create_hvg_plot <- function(seurat_obj) {
  if (length(VariableFeatures(seurat_obj)) == 0) {
    return(plotly_empty() %>% layout(title = "No variable features found"))
  }

  hvf_info <- tryCatch(HVFInfo(seurat_obj), error = function(e) NULL)
  if (is.null(hvf_info) || nrow(hvf_info) == 0) {
    return(plotly_empty() %>% layout(title = "HVG info not available"))
  }

  hvf_info$gene <- rownames(hvf_info)
  hvf_info$variable <- hvf_info$gene %in% VariableFeatures(seurat_obj)

  # Seurat v5 may use different column names; detect them
  mean_col <- if ("mean" %in% colnames(hvf_info)) "mean" else
              if ("gmean" %in% colnames(hvf_info)) "gmean" else NULL
  var_col <- if ("variance.standardized" %in% colnames(hvf_info)) "variance.standardized" else
             if ("residual_variance" %in% colnames(hvf_info)) "residual_variance" else
             if ("variance" %in% colnames(hvf_info)) "variance" else NULL

  if (is.null(mean_col) || is.null(var_col)) {
    # Fallback: just show which genes are variable using rank
    hvf_info$rank <- seq_len(nrow(hvf_info))
    colors <- ifelse(hvf_info$variable, "red", "gray70")
    return(plot_ly(
      x = hvf_info$rank, y = rep(1, nrow(hvf_info)),
      type = "scattergl", mode = "markers",
      marker = list(size = 3, opacity = 0.7, color = colors),
      text = hvf_info$gene,
      hovertemplate = "Gene: %{text}<br>Rank: %{x}<extra></extra>"
    ) %>% layout(title = paste(sum(hvf_info$variable), "Variable Features Selected"),
                 xaxis = list(title = "Gene Rank"), yaxis = list(title = ""), autosize = TRUE))
  }

  colors <- ifelse(hvf_info$variable, "red", "gray70")

  p <- plot_ly(
    x = hvf_info[[mean_col]],
    y = hvf_info[[var_col]],
    type = "scattergl",
    mode = "markers",
    marker = list(size = 3, opacity = 0.7, color = colors),
    text = hvf_info$gene,
    hovertemplate = "Gene: %{text}<br>Mean: %{x:.3f}<br>Variance: %{y:.3f}<extra></extra>"
  ) %>%
    layout(
      title = "Highly Variable Features",
      xaxis = list(title = "Average Expression (log)"),
      yaxis = list(title = "Standardized Variance")
    )

  # Add VST trend line overlay if data is available
  tryCatch({
    mean_vals <- hvf_info[[mean_col]]
    var_vals <- hvf_info[[var_col]]
    valid <- is.finite(mean_vals) & is.finite(var_vals) & mean_vals > 0
    if (sum(valid) > 20) {
      loess_fit <- loess(var_vals[valid] ~ mean_vals[valid], span = 0.3)
      ord <- order(mean_vals[valid])
      p <- p %>% add_trace(
        x = mean_vals[valid][ord],
        y = fitted(loess_fit)[ord],
        type = "scattergl",
        mode = "lines",
        line = list(color = "blue", width = 2, dash = "dash"),
        name = "VST trend",
        hoverinfo = "skip",
        inherit = FALSE
      )
    }
  }, error = function(e) NULL)

  p
}

create_enhanced_hvg_plot <- function(hvg_info, seurat_obj) {
  hvg_info$gene <- rownames(hvg_info)
  hvg_info$variable <- hvg_info$gene %in% VariableFeatures(seurat_obj)

  # Detect column names (Seurat v4 vs v5, LogNorm vs SCT)
  mean_col <- if ("mean" %in% colnames(hvg_info)) "mean" else
              if ("gmean" %in% colnames(hvg_info)) "gmean" else NULL
  var_col <- if ("variance.standardized" %in% colnames(hvg_info)) "variance.standardized" else
             if ("residual_variance" %in% colnames(hvg_info)) "residual_variance" else
             if ("variance" %in% colnames(hvg_info)) "variance" else NULL

  if (is.null(mean_col) || is.null(var_col)) {
    hvg_info$rank <- seq_len(nrow(hvg_info))
    colors <- ifelse(hvg_info$variable, "red", "gray70")
    return(plot_ly(x = hvg_info$rank, y = rep(1, nrow(hvg_info)),
      type = "scattergl", mode = "markers",
      marker = list(size = 3, opacity = 0.7, color = colors),
      text = hvg_info$gene,
      hovertemplate = "Gene: %{text}<br>Rank: %{x}<extra></extra>"
    ) %>% layout(title = paste(sum(hvg_info$variable), "Variable Features"),
                 xaxis = list(title = "Gene Rank"), yaxis = list(title = ""), autosize = TRUE))
  }

  # Add labels for top variable genes
  top_genes <- head(VariableFeatures(seurat_obj), 10)
  hvg_info$label <- ifelse(hvg_info$gene %in% top_genes, hvg_info$gene, "")

  p <- plot_ly(
    data = hvg_info,
    x = hvg_info[[mean_col]],
    y = hvg_info[[var_col]],
    color = ~variable,
    colors = c("FALSE" = "gray70", "TRUE" = "red"),
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.8),
    text = ~gene,
    hovertemplate = "Gene: %{text}<br>Mean: %{x:.3f}<br>Std Variance: %{y:.3f}<extra></extra>"
  ) %>%
    layout(
      title = "Highly Variable Features",
      xaxis = list(title = "Average Expression (log)"),
      yaxis = list(title = "Standardized Variance")
    )

  # Add VST trend line overlay if data is available
  tryCatch({
    mean_vals <- hvg_info[[mean_col]]
    var_vals <- hvg_info[[var_col]]
    valid <- is.finite(mean_vals) & is.finite(var_vals) & mean_vals > 0
    if (sum(valid) > 20) {
      loess_fit <- loess(var_vals[valid] ~ mean_vals[valid], span = 0.3)
      ord <- order(mean_vals[valid])
      p <- p %>% add_trace(
        x = mean_vals[valid][ord],
        y = fitted(loess_fit)[ord],
        type = "scatter",
        mode = "lines",
        line = list(color = "blue", width = 2, dash = "dash"),
        name = "VST trend",
        hoverinfo = "skip",
        inherit = FALSE
      )
    }
  }, error = function(e) NULL)

  p
}

create_hvg_summary_table <- function(hvg_info, seurat_obj) {
  var_col <- if ("variance.standardized" %in% colnames(hvg_info)) "variance.standardized" else
             if ("residual_variance" %in% colnames(hvg_info)) "residual_variance" else
             if ("variance" %in% colnames(hvg_info)) "variance" else NULL

  var_mean <- if (!is.null(var_col)) round(mean(hvg_info[[var_col]], na.rm = TRUE), 3) else "N/A"
  var_max  <- if (!is.null(var_col)) round(max(hvg_info[[var_col]], na.rm = TRUE), 3) else "N/A"

  summary_stats <- data.frame(
    Metric = c("Total Genes", "Variable Genes", "Mean Std Variance", "Max Std Variance"),
    Value = c(
      nrow(hvg_info),
      sum(rownames(hvg_info) %in% VariableFeatures(seurat_obj)),
      var_mean,
      var_max
    )
  )

  DT::datatable(summary_stats, options = list(dom = 't'))
}

create_top_hvg_table <- function(seurat_obj) {
  hvgs <- VariableFeatures(seurat_obj)
  if (length(hvgs) == 0) return(data.frame())

  hvg_info <- tryCatch(HVFInfo(seurat_obj), error = function(e) NULL)
  if (is.null(hvg_info)) return(DT::datatable(data.frame(Gene = hvgs), options = list(pageLength = 15)))

  top_hvg_info <- hvg_info[intersect(hvgs, rownames(hvg_info)), , drop = FALSE]
  top_hvg_info$Gene <- rownames(top_hvg_info)

  var_col <- if ("variance.standardized" %in% colnames(top_hvg_info)) "variance.standardized" else
             if ("residual_variance" %in% colnames(top_hvg_info)) "residual_variance" else
             if ("variance" %in% colnames(top_hvg_info)) "variance" else NULL

  if (!is.null(var_col)) {
    top_hvg_info <- top_hvg_info[order(top_hvg_info[[var_col]], decreasing = TRUE), ]
  }

  DT::datatable(head(top_hvg_info, 50), options = list(pageLength = 15, scrollX = TRUE))
}

create_scaling_stats_table <- function(scale_results) {
  stats_df <- data.frame(
    Metric = c("Features Scaled", "Mean Variance (Post-scaling)", "Scaling Method"),
    Value = c(
      scale_results$n_features_scaled,
      round(scale_results$mean_variance, 3),
      "ScaleData"
    )
  )
  
  DT::datatable(stats_df, options = list(dom = 't'))
}

create_scaling_effect_plot <- function(scale_results) {
  tryCatch({
    obj <- scale_results$object

    # Get scaled and unscaled data
    scaled_data   <- safe_get_assay_data(obj, "scale.data")
    unscaled_data <- safe_get_assay_data(obj, "data")

    # Sample ~20 genes from the scaled features for visualization
    set.seed(42)
    available_genes <- rownames(scaled_data)
    n_genes <- min(20, length(available_genes))
    sampled_genes <- sample(available_genes, n_genes)

    # Sample cells to keep plot manageable
    n_cells <- min(500, ncol(scaled_data))
    sampled_cells <- sample(colnames(scaled_data), n_cells)

    # Build long-format data for unscaled (before)
    before_vals <- as.matrix(unscaled_data[sampled_genes, sampled_cells, drop = FALSE])
    before_df <- data.frame(
      Gene = rep(rownames(before_vals), each = ncol(before_vals)),
      Expression = as.vector(t(before_vals))
    )

    # Build long-format data for scaled (after)
    after_vals <- as.matrix(scaled_data[sampled_genes, sampled_cells, drop = FALSE])
    after_df <- data.frame(
      Gene = rep(rownames(after_vals), each = ncol(after_vals)),
      Expression = as.vector(t(after_vals))
    )

    # Left panel: box plots before scaling
    p_before <- plot_ly(
      data = before_df,
      x = ~Gene,
      y = ~Expression,
      type = "box",
      marker = list(size = 2),
      line = list(width = 1),
      fillcolor = "rgba(70, 130, 180, 0.5)",
      name = "Before Scaling"
    ) %>%
      layout(
        xaxis = list(title = "Gene", tickangle = 45, tickfont = list(size = 8)),
        yaxis = list(title = "Expression (log-normalized)"),
        annotations = list(
          list(text = "Before Scaling", xref = "paper", yref = "paper",
               x = 0.5, y = 1.05, showarrow = FALSE,
               font = list(size = 14, weight = "bold"))
        )
      )

    # Right panel: box plots after scaling
    p_after <- plot_ly(
      data = after_df,
      x = ~Gene,
      y = ~Expression,
      type = "box",
      marker = list(size = 2),
      line = list(width = 1),
      fillcolor = "rgba(178, 34, 34, 0.5)",
      name = "After Scaling"
    ) %>%
      layout(
        xaxis = list(title = "Gene", tickangle = 45, tickfont = list(size = 8)),
        yaxis = list(title = "Expression (scaled)"),
        annotations = list(
          list(text = "After Scaling", xref = "paper", yref = "paper",
               x = 0.5, y = 1.05, showarrow = FALSE,
               font = list(size = 14, weight = "bold"))
        ),
        shapes = list(
          list(type = "line", x0 = 0, x1 = 1, xref = "paper",
               y0 = 0, y1 = 0,
               line = list(color = "red", width = 1, dash = "dash"))
        )
      )

    plotly::subplot(p_before, p_after, nrows = 1, shareY = FALSE,
                    titleX = TRUE, titleY = TRUE, margin = 0.08) %>%
      layout(
        title = list(
          text = paste("Scaling Effect on", n_genes, "Sampled Genes"),
          y = 0.98
        ),
        showlegend = FALSE
      )
  }, error = function(e) {
    plotly_empty() %>%
      layout(
        title = paste("Scaling completed for", scale_results$n_features_scaled, "features"),
        annotations = list(
          text = paste("Could not render before/after plot:", e$message),
          x = 0.5, y = 0.5,
          showarrow = FALSE,
          font = list(size = 14)
        )
      )
  })
}

create_integration_summary_table <- function(integration_results) {
  summary_df <- data.frame(
    Metric = c("Integration Method", "Strategy", "Datasets Integrated", "Total Cells", "Total Genes"),
    Value = c(
      integration_results$method,
      integration_results$strategy,
      integration_results$n_datasets %||% 1,
      ncol(integration_results$object),
      nrow(integration_results$object)
    )
  )
  
  DT::datatable(summary_df, options = list(dom = 't'))
}

create_integration_metrics_table <- function(integration_results) {
  obj <- integration_results$object
  
  metrics_df <- data.frame(
    Metric = c("Active Assay", "Available Reductions", "Default Assay Genes", "Integrated Assay Available"),
    Value = c(
      DefaultAssay(obj),
      paste(names(obj@reductions), collapse = ", "),
      nrow(obj@assays[[DefaultAssay(obj)]]),
      "integrated" %in% names(obj@assays)
    )
  )
  
  DT::datatable(metrics_df, options = list(dom = 't'))
}

create_integration_comparison_plot <- function(integration_results) {
  obj <- integration_results$object
  
  # Find best available reduction
  reductions <- names(obj@reductions)
  reduction <- if ("umap" %in% reductions) "umap" else if ("tsne" %in% reductions) "tsne" else "pca"
  
  if (!reduction %in% reductions) {
    return(plotly_empty() %>% layout(title = "No reduction available for visualization"))
  }
  
  embedding <- Embeddings(obj, reduction)[, 1:2]
  
  # Color by batch if available
  batch_cols <- grep("batch|orig.ident|dataset", colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
  color_by <- if (length(batch_cols) > 0) batch_cols[1] else NULL
  
  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2]
  )
  
  if (!is.null(color_by)) {
    plot_data$Batch <- obj@meta.data[[color_by]]
  }
  
  plot_args <- list(
    data = plot_data,
    x = ~Dim1,
    y = ~Dim2,
    type = "scatter",
    mode = "markers",
    marker = list(size = 3, opacity = 0.7)
  )
  if (!is.null(color_by)) {
    plot_args$color <- ~Batch
  }
  p <- do.call(plot_ly, plot_args)
  
  p %>% layout(
    title = paste("Integration Results on", toupper(reduction)),
    xaxis = list(title = paste0(toupper(reduction), "1")),
    yaxis = list(title = paste0(toupper(reduction), "2")),
    autosize = TRUE
  )
}
