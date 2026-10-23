# Dimensionality Reduction Module (enhanced and robust)
# modules/dimensionality_reduction.R

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(plotly)
  library(patchwork)
  # Optional packages — loaded on demand, not required at startup
  # destiny (diffusion maps), uwot (UMAP), Rtsne (t-SNE), irlba (fast PCA)
  for (.pkg in c("uwot", "Rtsne", "irlba")) {
    if (requireNamespace(.pkg, quietly = TRUE)) library(.pkg, character.only = TRUE)
  }
})

# get_color_palette(), best_reduction(), safe_default_assay(),
# detect_cell_type_col() are defined in core_robustness.R

# safe_get_assay_data() is defined in core_robustness.R

validate_dimred_contract <- function(values, obj, stage = "dimensionality_reduction", require_reductions = NULL) {
  if (!exists("validate_seurat_contract", mode = "function")) return(TRUE)

  contract <- validate_seurat_contract(
    obj,
    stage = stage,
    require_reductions = require_reductions,
    min_cells = 20,
    min_features = 50
  )
  if (!isTRUE(contract$ok)) {
    msg <- if (exists("format_contract_feedback", mode = "function")) {
      format_contract_feedback(contract)
    } else {
      paste(contract$errors, collapse = " | ")
    }
    showNotification(msg, type = "error", duration = 12)
    if (exists("append_audit_log", mode = "function")) {
      append_audit_log(
        values,
        module = "dimensionality_reduction",
        action = paste0("contract_failed:", stage),
        level = "error",
        details = list(errors = contract$errors, warnings = contract$warnings, suggestions = contract$suggestions)
      )
    }
    return(FALSE)
  }
  if (length(contract$warnings) > 0) {
    showNotification(paste("Contract warnings:", paste(contract$warnings, collapse = " | ")), type = "warning", duration = 8)
  }
  TRUE
}

begin_dimred_job <- function(values, action, meta = list()) {
  if (!exists("start_job", mode = "function")) return(NULL)
  start_job(values, module = "dimensionality_reduction", action = action, meta = meta)
}

end_dimred_job <- function(values, job_id, status = "done", message = NULL, details = list()) {
  if (is.null(job_id) || !exists("finish_job", mode = "function")) return(invisible(NULL))
  finish_job(values, job_id = job_id, status = status, message = message, details = details)
}

checkpoint_dimred <- function(values, action, obj, params = list(), method = NULL) {
  if (!exists("register_dataset_checkpoint", mode = "function")) return(invisible(NULL))
  register_dataset_checkpoint(
    values,
    module = "dimensionality_reduction",
    action = action,
    obj = obj,
    params = params,
    methods = list(method = method)
  )
}

# -----------------------------------------------------------------------------
# ENHANCED PCA UI & SERVER
# -----------------------------------------------------------------------------

pcaUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Principal Component Analysis (PCA)"),
          p("Perform linear dimensionality reduction to identify major sources of variation."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("dimensionality_reduction"),
          
          fluidRow(
            column(4,
              numericInput(ns("n_pcs"), "Number of PCs", 
                         value = 50, min = 10, max = 200, step = 5)
            ),
            column(4,
              selectInput(ns("pca_features"), "Features to Use",
                choices = c("Variable features" = "hvg", "All features" = "all"),
                selected = "hvg")
            ),
            column(4,
              selectInput(ns("pca_assay"), "Assay", choices = NULL)
            )
          ),
          
          fluidRow(
            column(3,
              checkboxInput(ns("use_irlba"), "Use IRLBA (faster)", TRUE)
            ),
            column(3,
              numericInput(ns("pca_seed"), "Random seed", value = 42, min = 1, max = 10000)
            ),
            column(3,
              checkboxInput(ns("run_jackstraw"), "Run JackStraw (slow)", FALSE),
              helpText("Optional: statistical PC selection (p < 0.01). Computationally expensive.")
            ),
            column(3,
              div(style = "margin-top: 25px;",
                actionButton(ns("run_pca"), "Run PCA",
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          ),

          fluidRow(
            column(12,
              div(style = "margin-top: 8px;",
                actionButton(ns("ai_pc_suggest"), "AI Suggest PCs",
                           icon = icon("robot"), class = "ai-btn btn-sm")
              ),
              uiOutput(ns("ai_thread_pc_suggest"))
            )
          )
        )
      )
    ),
    
    # Results display
    fluidRow(
      column(12,
        uiOutput(ns("pca_results_ui"))
      )
    )
  )
}

pcaServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_pca_results <- reactiveVal(FALSE)
    pca_results <- reactiveVal(NULL)
    
    output$pca_results_ui <- renderUI({
      req(isTRUE(show_pca_results()))
      ns <- session$ns
      wellPanel(
        h4("PCA Results"),

        tabsetPanel(
          tabPanel("Elbow Plot",
            fluidRow(
              column(8,
                plotlyOutput(ns("elbow_plot"), height = "450px")
              ),
              column(4,
                h5("PC Statistics"),
                DT::dataTableOutput(ns("pc_stats_table"))
              )
            )
          ),

          tabPanel("PC Heatmaps",
            fluidRow(
              column(3,
                numericInput(ns("heatmap_dims"), "PCs to show", value = 6, min = 1, max = 30)
              ),
              column(3,
                numericInput(ns("heatmap_cells"), "Cells to sample", value = 200, min = 50, max = 2000)
              ),
              column(3,
                numericInput(ns("heatmap_top_genes"), "Top genes per PC", value = 20, min = 5, max = 50)
              ),
              column(3,
                checkboxInput(ns("balanced_heatmap"), "Balanced sampling", TRUE)
              )
            ),
            uiOutput(ns("pc_heatmap_panels"))
          ),

          tabPanel("PC Plots",
            fluidRow(
              column(4,
                selectInput(ns("pc_x"), "PC X-axis", choices = NULL, selected = "PC_1")
              ),
              column(4,
                selectInput(ns("pc_y"), "PC Y-axis", choices = NULL, selected = "PC_2")
              ),
              column(4,
                selectInput(ns("color_by_pca"), "Color by", choices = NULL)
              )
            ),
            plotlyOutput(ns("pca_scatter"), height = "500px")
          ),

          tabPanel("Loadings",
            fluidRow(
              column(4,
                selectInput(ns("loading_pc"), "Select PC", choices = NULL)
              ),
              column(4,
                numericInput(ns("top_loadings"), "Top loadings to show", value = 20, min = 10, max = 50)
              ),
              column(4,
                radioButtons(ns("loading_type"), "Show",
                  choices = c("Positive" = "pos", "Negative" = "neg", "Both" = "both"),
                  selected = "both")
              )
            ),
            plotlyOutput(ns("loadings_plot"), height = "500px")
          )
        )
      )
    })
    
    # Update UI choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      
      # Update assay choices
      assays <- names(obj@assays)
      updateSelectInput(session, "pca_assay", choices = assays, selected = safe_default_assay(obj))
      
      # Update metadata for coloring
      meta_cols <- c("None" = "none", colnames(obj@meta.data))
      updateSelectInput(session, "color_by_pca", choices = meta_cols)
    })
    
    # Update PC choices after PCA is run
    observeEvent(pca_results(), {
      tryCatch({
        req(pca_results())
        pc_names <- paste0("PC_", 1:pca_results()$n_pcs)
        updateSelectInput(session, "pc_x", choices = pc_names, selected = "PC_1")
        updateSelectInput(session, "pc_y", choices = pc_names, selected = "PC_2")
        updateSelectInput(session, "loading_pc", choices = pc_names, selected = "PC_1")
      }, error = function(e) NULL)
    })
    
    # Run PCA
    observeEvent(input$run_pca, {
      req(values$active_dataset)
      job_id <- begin_dimred_job(values, "run_pca", meta = list(
        n_pcs = input$n_pcs,
        pca_features = input$pca_features,
        assay = input$pca_assay
      ))

      if (!validate_dimred_contract(values, values$active_dataset, stage = "dimensionality_reduction")) {
        end_dimred_job(values, job_id, status = "failed", message = "contract validation failed")
        return()
      }

      tryCatch({
        withProgress(message = "Running PCA...", value = 0, {
          obj <- values$active_dataset

          if (exists("runtime_budget_guard", mode = "function")) {
            budget <- runtime_budget_guard(obj, task = "pca")
            if (!isTRUE(budget$ok)) {
              showNotification(paste(budget$message, budget$suggestion), type = "error", duration = 10)
              end_dimred_job(values, job_id, status = "failed", message = budget$message)
              return()
            }
          }

          incProgress(0.3, detail = "Preparing features...")

          DefaultAssay(obj) <- input$pca_assay

          if (input$pca_features == "hvg") {
            features <- VariableFeatures(obj)
            if (length(features) == 0) {
              showNotification("No variable features found, using all features", type = "warning")
              features <- rownames(obj)
            }
          } else {
            features <- rownames(obj)
          }

          incProgress(0.6, detail = "Computing PCA...")

          result <- safely_run_pca(obj, input, features)

          if (result$success) {
            values$active_dataset <- result$object
            pca_results(result)
            show_pca_results(TRUE)
            checkpoint_dimred(
              values,
              action = "run_pca",
              obj = result$object,
              params = list(
                n_pcs = input$n_pcs,
                features_mode = input$pca_features,
                assay = input$pca_assay,
                seed = input$pca_seed
              ),
              method = "Seurat::RunPCA"
            )
            end_dimred_job(values, job_id, status = "done", message = "PCA completed")

            # Issue 4.1: Auto-detect elbow and notify user
            elbow_pc <- detect_elbow_point(result$stdev)
            values$.detected_elbow_pc <- elbow_pc
            values$dimred_done <- TRUE
            showNotification(
              paste0("PCA completed with ", result$n_pcs, " components. ",
                     "Suggested elbow point: PC ", elbow_pc,
                     " (auto-detected via second derivative of scree plot)."),
              type = "message", duration = 10)

            # Issue 4.2: Optional JackStraw statistical PC selection
            if (isTRUE(input$run_jackstraw)) {
              showNotification("Running JackStraw analysis (this may take a while)...", type = "message", duration = 6)
              tryCatch({
                js_dims <- min(input$n_pcs, 20)
                if (input$n_pcs > 20) {
                  showNotification("JackStraw evaluation capped at 20 PCs due to computational cost.", type = "warning", duration = 8)
                }
                js_obj <- JackStraw(values$active_dataset, num.replicate = 100,
                                    dims = js_dims)
                js_obj <- ScoreJackStraw(js_obj, dims = 1:js_dims)
                values$active_dataset <- js_obj

                # Identify significant PCs (p < 0.01)
                js_pvalues <- js_obj@reductions$pca@jackstraw$overall.p.values
                sig_pcs <- which(js_pvalues[, "Score"] < 0.01)
                if (length(sig_pcs) > 0) {
                  max_sig_pc <- max(sig_pcs)
                  showNotification(
                    paste0("JackStraw: ", length(sig_pcs), " significant PCs (p < 0.01). ",
                           "Highest significant PC: ", max_sig_pc, "."),
                    type = "message", duration = 12)
                } else {
                  showNotification("JackStraw: No PCs reached significance at p < 0.01.",
                                   type = "warning", duration = 8)
                }
              }, error = function(e) {
                showNotification(paste("JackStraw failed:", e$message), type = "warning", duration = 8)
              })
            }
          } else {
            end_dimred_job(values, job_id, status = "failed", message = result$error)
            showNotification(paste("PCA failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        end_dimred_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("PCA error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$elbow_plot <- renderPlotly({
      req(pca_results())
      p <- create_enhanced_elbow_plot(pca_results())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "Elbow Plot",
                       label = "Elbow Plot", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    output$pc_stats_table <- DT::renderDataTable({
      req(pca_results())
      tbl <- create_pc_stats_table(pca_results())
      if (exists("register_csv_export", mode = "function")) {
        pc_stats <- data.frame(
          PC = paste0("PC", 1:length(pca_results()$var_explained)),
          Variance_Explained = round(pca_results()$var_explained, 2),
          Cumulative_Variance = round(pca_results()$cumvar_explained, 2),
          Standard_Deviation = round(pca_results()$stdev, 3)
        )
        register_csv_export(values, data = head(pc_stats, 20), description = "PC Statistics", source_module = "Dimensionality Reduction")
      }
      tbl
    }, options = list(dom = 't', pageLength = 15))
    
    # Dynamic PC heatmap panels — one interactive plotly heatmap per PC
    output$pc_heatmap_panels <- renderUI({
      req(pca_results(), values$active_dataset)
      n_dims <- min(input$heatmap_dims, ncol(values$active_dataset@reductions$pca))
      ns <- session$ns

      panels <- lapply(seq_len(n_dims), function(i) {
        plot_id <- ns(paste0("pc_heatmap_", i))
        dl_id   <- ns(paste0("dl_pc_heatmap_", i))
        wellPanel(
          style = "margin-bottom: 10px; padding: 8px;",
          div(style = "display: flex; justify-content: space-between; align-items: center;",
            h5(paste0("PC ", i), style = "margin: 0;"),
            downloadButton(dl_id, "CSV", class = "btn-xs btn-default")
          ),
          plotlyOutput(plot_id, height = "350px")
        )
      })

      do.call(tagList, panels)
    })

    # Create individual PC heatmap plotly outputs and download handlers
    observe({
      req(pca_results(), values$active_dataset)
      n_dims <- min(input$heatmap_dims, ncol(values$active_dataset@reductions$pca))
      top_n  <- input$heatmap_top_genes
      n_cells <- input$heatmap_cells

      obj <- values$active_dataset
      loadings_mat <- Loadings(obj, reduction = "pca")
      embedding   <- Embeddings(obj, "pca")

      for (i in seq_len(n_dims)) {
        local({
          pc_i <- i

          output[[paste0("pc_heatmap_", pc_i)]] <- renderPlotly({
            tryCatch({
              pc_loads <- loadings_mat[, pc_i]
              top_pos <- names(sort(pc_loads, decreasing = TRUE))[seq_len(min(top_n, length(pc_loads)))]
              top_neg <- names(sort(pc_loads, decreasing = FALSE))[seq_len(min(top_n, length(pc_loads)))]
              sel_genes <- unique(c(top_pos, rev(top_neg)))

              # Sample cells ordered by PC score
              pc_scores <- embedding[, pc_i]
              cell_order <- order(pc_scores)
              n_sample <- min(n_cells, length(cell_order))
              set.seed(42)
              sample_idx <- sort(sample(cell_order, n_sample))

              expr <- tryCatch({
                  sd <- safe_get_assay_data(obj, "scale.data")
                  as.matrix(sd[sel_genes, sample_idx, drop = FALSE])
                },
                error = function(e) {
                  # Fallback to data layer
                  mat <- tryCatch(
                    as.matrix(safe_get_assay_data(obj, "data")[sel_genes, sample_idx, drop = FALSE]),
                    error = function(e2) {
                      as.matrix(GetAssayData(obj, layer = "data")[sel_genes, sample_idx, drop = FALSE])
                    }
                  )
                  t(scale(t(mat)))
                }
              )
              # Clamp for visualization
              expr[expr >  2.5] <-  2.5
              expr[expr < -2.5] <- -2.5

              plot_ly(
                z = expr,
                x = seq_len(ncol(expr)),
                y = rownames(expr),
                type = "heatmap",
                colorscale = list(c(0, "#313695"), c(0.5, "#FFFFBF"), c(1, "#A50026")),
                hovertemplate = "Gene: %{y}<br>Cell: %{x}<br>Expr: %{z:.2f}<extra></extra>"
              ) %>% layout(
                title = paste0("PC ", pc_i, " — Top ", top_n, " loading genes"),
                xaxis = list(title = "Cells (ordered by PC score)", showticklabels = FALSE),
                yaxis = list(title = "", tickfont = list(size = 9)),
                autosize = TRUE,
                margin = list(l = 100)
              ) %>% plotly::config(
                toImageButtonOptions = list(
                  format = "png", filename = paste0("PC", pc_i, "_heatmap"),
                  width = 1200, height = 500
                )
              )
            }, error = function(e) {
              plotly_empty() %>% layout(title = paste0("PC ", pc_i, " error: ", e$message))
            })
          })

          output[[paste0("dl_pc_heatmap_", pc_i)]] <- downloadHandler(
            filename = function() paste0("PC", pc_i, "_top_genes_", Sys.Date(), ".csv"),
            content = function(file) {
              tryCatch({
                pc_loads <- loadings_mat[, pc_i]
                top_pos <- names(sort(pc_loads, decreasing = TRUE))[seq_len(min(top_n, length(pc_loads)))]
                top_neg <- names(sort(pc_loads, decreasing = FALSE))[seq_len(min(top_n, length(pc_loads)))]
                sel_genes <- unique(c(top_pos, rev(top_neg)))
                df <- data.frame(Gene = sel_genes, Loading = pc_loads[sel_genes])
                df <- df[order(abs(df$Loading), decreasing = TRUE), ]
                write.csv(df, file, row.names = FALSE)
              }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
            }
          )
        })
      }
    })
    
    output$pca_scatter <- renderPlotly({
      req(pca_results(), values$active_dataset)
      p <- create_pca_scatter_plot(values$active_dataset, input$pc_x, input$pc_y, input$color_by_pca)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "PCA Scatter",
                       label = "PCA Scatter", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    output$loadings_plot <- renderPlotly({
      req(pca_results(), values$active_dataset)
      p <- create_loadings_plot(values$active_dataset, input$loading_pc, input$top_loadings, input$loading_type)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "PCA Loadings",
                       label = "PCA Loadings", source_module = "Dimensionality Reduction")
      }
      p
    })

    # --- AI: Suggest optimal number of PCs ---
    observeEvent(input$ai_pc_suggest, {
      if (is.null(values$ai_api_key) || !nzchar(trimws(values$ai_api_key))) {
        showNotification("Please set your Claude API key in the Data Management tab first.",
                         type = "warning", duration = 8)
        return()
      }
      if (is.null(values$active_dataset)) {
        showNotification("No active dataset loaded.", type = "warning", duration = 6)
        return()
      }
      obj <- values$active_dataset
      if (!("pca" %in% names(obj@reductions))) {
        showNotification("PCA has not been computed yet. Please run PCA first.",
                         type = "warning", duration = 6)
        return()
      }

      stdev <- obj@reductions$pca@stdev
      n_pcs <- min(length(stdev), 50)
      stdev_subset <- round(stdev[seq_len(n_pcs)], 4)

      prompt <- paste0(
        "Here are the PCA standard deviations for ", n_pcs, " PCs: ",
        paste(stdev_subset, collapse = ", "),
        ". Based on the elbow plot, suggest the optimal number of PCs to use ",
        "for downstream analysis (UMAP, clustering). Return a single number with brief reasoning."
      )

      output$ai_thread_pc_suggest <- renderUI({ ai_thread_loading("AI analyzing PCA variance...") })
      result <- call_ai("pc_selection", prompt, values)

      if (result$ok) {
        # Try to extract a number from the response for Apply
        pc_num <- as.integer(regmatches(result$text, regexpr("\\b([12345][0-9]|[1-9])\\b", result$text)))
        ai_pc_val <- if (length(pc_num) > 0 && !is.na(pc_num[1])) pc_num[1] else NULL

        apply_btn <- if (!is.null(ai_pc_val)) {
          actionButton(session$ns("apply_ai_pc"), paste0("Apply ", ai_pc_val, " PCs"),
                       class = "ai-thread-apply-btn", icon = icon("check"))
        } else NULL

        output$ai_thread_pc_suggest <- renderUI({
          ai_thread_message(
            content_html = .mdtohtml(result$text),
            title = "AI PC Suggestion",
            model = get_model_for_task("pc_selection", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out,
            apply_btn = apply_btn
          )
        })

        # Store value for apply handler
        values$.ai_pc_suggestion <- ai_pc_val
      } else {
        output$ai_thread_pc_suggest <- renderUI({ ai_thread_error(result$error) })
      }
    })

    # Apply AI-suggested PC count
    observeEvent(input$apply_ai_pc, {
      pc_val <- values$.ai_pc_suggestion
      req(!is.null(pc_val))
      updateNumericInput(session, "n_pcs", value = pc_val)
      showNotification(paste0("PC count set to ", pc_val, "."), type = "message", duration = 4)
    })

    # --- AI Helper: generic Q&A for Dimensionality Reduction ---
    ai_helper_server("pca", values, "Dimensionality Reduction", context_fn = function() {
      obj <- values$active_dataset
      if (is.null(obj) || !inherits(obj, "Seurat")) return("")
      parts <- character(0)
      reds <- names(obj@reductions)
      if (length(reds) > 0) {
        parts <- c(parts, paste0("Available reductions: ", paste(reds, collapse = ", ")))
      }
      if ("pca" %in% reds) {
        stdev <- obj@reductions$pca@stdev
        n_pcs <- length(stdev)
        total_var <- sum(stdev^2)
        cum_var <- cumsum(stdev^2) / total_var * 100
        parts <- c(parts, paste0("PCA: ", n_pcs, " components computed"))
        parts <- c(parts, paste0("Variance explained by PC1: ", round(stdev[1]^2 / total_var * 100, 1), "%"))
        top_n <- min(10, n_pcs)
        parts <- c(parts, paste0("Cumulative variance (first ", top_n, " PCs): ",
                                  round(cum_var[top_n], 1), "%"))
      }
      paste(parts, collapse = "\n")
    })
  })
}

# -----------------------------------------------------------------------------
# ENHANCED UMAP/t-SNE UI & SERVER
# -----------------------------------------------------------------------------

umapTsneUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Non-linear Dimensionality Reduction"),
          p("Generate 2D embeddings using UMAP or t-SNE for visualization and clustering."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("dimensionality_reduction"),
          
          fluidRow(
            column(4,
              selectInput(ns("reduction_method"), "Method",
                choices = c("UMAP" = "umap", "t-SNE" = "tsne"),
                selected = "umap")
            ),
            column(4,
              selectInput(ns("reduction_input"), "Input Reduction",
                choices = c("PCA" = "pca", "Harmony" = "harmony", "Integrated" = "integrated"),
                selected = "pca")
            ),
            column(4,
              numericInput(ns("dims_to_use"), "Dimensions to Use", 
                         value = 30, min = 5, max = 100, step = 5)
            )
          ),
          
          # UMAP-specific parameters
          conditionalPanel(
            condition = "input.reduction_method == 'umap'",
            ns = ns,
            wellPanel(
              h5("UMAP Parameters"),
              fluidRow(
                column(3,
                  numericInput(ns("n_neighbors"), "Neighbors", value = 30, min = 5, max = 200)
                ),
                column(3,
                  numericInput(ns("min_dist"), "Min Distance", value = 0.3, min = 0.001, max = 2, step = 0.1)
                ),
                column(3,
                  numericInput(ns("spread"), "Spread", value = 1.0, min = 0.1, max = 5, step = 0.1)
                ),
                column(3,
                  selectInput(ns("umap_metric"), "Distance Metric",
                    choices = c("euclidean", "cosine", "manhattan", "correlation"),
                    selected = "euclidean")
                )
              ),
              fluidRow(
                column(4,
                  numericInput(ns("a"), "a parameter", value = NULL, min = 0.1, max = 5)
                ),
                column(4,
                  numericInput(ns("b"), "b parameter", value = NULL, min = 0.1, max = 5)
                ),
                column(4,
                  checkboxInput(ns("densmap"), "Use DensMAP", FALSE)
                )
              )
            )
          ),
          
          # t-SNE-specific parameters
          conditionalPanel(
            condition = "input.reduction_method == 'tsne'",
            ns = ns,
            wellPanel(
              h5("t-SNE Parameters"),
              fluidRow(
                column(3,
                  numericInput(ns("perplexity"), "Perplexity", value = 30, min = 5, max = 100)
                ),
                column(3,
                  numericInput(ns("theta"), "Theta (speed)", value = 0.5, min = 0, max = 1, step = 0.1)
                ),
                column(3,
                  numericInput(ns("max_iter"), "Max Iterations", value = 1000, min = 250, max = 5000, step = 250)
                ),
                column(3,
                  numericInput(ns("learning_rate"), "Learning Rate", value = 200, min = 10, max = 1000, step = 50)
                )
              )
            )
          ),
          
          # General options
          fluidRow(
            column(4,
              numericInput(ns("reduction_seed"), "Random Seed", value = 42, min = 1, max = 10000)
            ),
            column(4,
              checkboxInput(ns("return_model"), "Return model (UMAP)", TRUE)
            ),
            column(4,
              div(style = "margin-top: 25px;",
                actionButton(ns("run_reduction"), "Run Reduction", 
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    # Results
    fluidRow(
      column(12,
        uiOutput(ns("reduction_results_ui"))
      )
    )
  )
}

umapTsneServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_reduction_results <- reactiveVal(FALSE)
    reduction_results <- reactiveVal(NULL)
    
    output$reduction_results_ui <- renderUI({
      req(isTRUE(show_reduction_results()))
      ns <- session$ns
      wellPanel(
        h4("Reduction Results"),

        tabsetPanel(
          tabPanel("Embedding Plot",
            fluidRow(
              column(9,
                plotlyOutput(ns("reduction_plot"), height = "600px")
              ),
              column(3,
                h5("Plot Options"),
                selectInput(ns("color_by_reduction"), "Color by", choices = NULL),
                selectInput(ns("shape_by_reduction"), "Shape by", choices = c("None" = "none")),
                sliderInput(ns("point_size"), "Point size", min = 0.1, max = 5, value = 1, step = 0.1),
                sliderInput(ns("point_alpha"), "Point opacity", min = 0.1, max = 1, value = 0.8, step = 0.1),
                checkboxInput(ns("show_legend"), "Show legend", TRUE)
              )
            )
          ),

          tabPanel("Quality Metrics",
            fluidRow(
              column(6,
                h5("Embedding Quality"),
                DT::dataTableOutput(ns("reduction_quality_table"))
              ),
              column(6,
                h5("Parameter Summary"),
                DT::dataTableOutput(ns("reduction_params_table"))
              )
            )
          ),

          tabPanel("Density Plot",
            plotlyOutput(ns("density_plot"), height = "500px")
          )
        )
      )
    })
    
    # Update UI choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset

      # Available input reductions
      available_reductions <- names(obj@reductions)
      updateSelectInput(session, "reduction_input",
                       choices = available_reductions,
                       selected = if ("pca" %in% available_reductions) "pca" else available_reductions[1])

      # Issue 4.5: Auto-compute t-SNE perplexity based on dataset size
      auto_perp <- max(5, min(50, floor(sqrt(ncol(obj)) / 3)))
      updateNumericInput(session, "perplexity", value = auto_perp)
      showNotification(
        paste0("t-SNE perplexity auto-set to ", auto_perp,
               " (based on ", ncol(obj), " cells: floor(sqrt(n)/3))."),
        type = "message", duration = 6)

      meta_all <- colnames(obj@meta.data)
factor_cols <- meta_all[sapply(obj@meta.data[meta_all], function(x) is.factor(x) || is.character(x))]
updateSelectInput(session, "color_by_reduction", choices = c("None" = "none", stats::setNames(meta_all, meta_all)))
updateSelectInput(session, "shape_by_reduction", choices = c("None" = "none", stats::setNames(factor_cols, factor_cols)))
    })
    
    # Run dimensionality reduction
    observeEvent(input$run_reduction, {
      req(values$active_dataset, input$reduction_input)
      job_id <- begin_dimred_job(values, "run_nonlinear_reduction", meta = list(
        method = input$reduction_method,
        input_reduction = input$reduction_input,
        dims = input$dims_to_use
      ))
      tryCatch({
        obj <- values$active_dataset
        selected_input <- input$reduction_input
        available_reductions <- names(obj@reductions)

        if (!selected_input %in% available_reductions) {
          fallback <- if (exists("find_best_reduction", mode = "function")) {
            find_best_reduction(obj, preferred = c("harmony", "integrated", "mnn", "pca", "umap", "tsne"))
          } else {
            NULL
          }
          if (!is.null(fallback) && fallback %in% available_reductions) {
            selected_input <- fallback
            updateSelectInput(session, "reduction_input", selected = fallback)
            showNotification(
              paste("Selected input reduction missing. Falling back to", fallback, "automatically."),
              type = "warning",
              duration = 8
            )
          } else {
            showNotification("No valid input reduction available. Run PCA first.", type = "error")
            end_dimred_job(values, job_id, status = "failed", message = "missing input reduction")
            return()
          }
        }

        if (!validate_dimred_contract(values, obj, stage = "dimensionality_reduction", require_reductions = selected_input)) {
          end_dimred_job(values, job_id, status = "failed", message = "contract validation failed")
          return()
        }

        # t-SNE perplexity validation: must be < (n_cells - 1) / 3
        if (input$reduction_method == "tsne") {
          max_perplexity <- floor((ncol(obj) - 1) / 3)
          if (input$perplexity >= max_perplexity) {
            safe_perp <- max(5, max_perplexity - 1)
            showNotification(
              paste0("Perplexity (", input$perplexity, ") must be < ", max_perplexity,
                     " for ", ncol(obj), " cells. Auto-adjusting to ", safe_perp, "."),
              type = "warning", duration = 8
            )
            updateNumericInput(session, "perplexity", value = safe_perp)
            end_dimred_job(values, job_id, status = "failed", message = "perplexity too high")
            return()
          }
        }

        withProgress(message = paste("Running", toupper(input$reduction_method), "..."), value = 0, {
          incProgress(0.3, detail = "Preparing input data...")

          params <- reactiveValuesToList(input)
          params$reduction_input <- selected_input
          result <- safely_run_nonlinear_reduction(obj, params)

          if (result$success) {
            values$active_dataset <- result$object
            reduction_results(result)
            show_reduction_results(TRUE)
            checkpoint_dimred(
              values,
              action = paste0("run_", input$reduction_method),
              obj = result$object,
              params = list(
                method = input$reduction_method,
                input_reduction = selected_input,
                dims = input$dims_to_use,
                seed = input$reduction_seed
              ),
              method = paste0("Seurat::Run", toupper(input$reduction_method))
            )
            end_dimred_job(values, job_id, status = "done", message = paste(input$reduction_method, "completed"))

            showNotification(paste(toupper(input$reduction_method), "completed successfully"), type = "message")
          } else {
            end_dimred_job(values, job_id, status = "failed", message = result$error)
            showNotification(paste(toupper(input$reduction_method), "failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        end_dimred_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("Reduction error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$reduction_plot <- renderPlotly({
      req(reduction_results(), values$active_dataset)
      p <- create_enhanced_reduction_plot(values$active_dataset, reduction_results(), input)
      if (exists("capture_umap_snapshot", mode = "function"))
        capture_umap_snapshot(values, p, "Dimensionality Reduction",
          coloring = tryCatch(input$color_by, error = function(e) ""),
          method = tryCatch(reduction_results()$method, error = function(e) ""))
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "UMAP/tSNE Plot",
                       label = "UMAP/tSNE Plot", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    output$reduction_quality_table <- DT::renderDataTable({
      req(reduction_results())
      create_reduction_quality_table(reduction_results())
    }, options = list(dom = 't'))
    
    output$reduction_params_table <- DT::renderDataTable({
      req(reduction_results())
      create_reduction_params_table(reduction_results(), input)
    }, options = list(dom = 't'))
    
    output$density_plot <- renderPlotly({
      req(reduction_results(), values$active_dataset)
      p <- create_density_plot(values$active_dataset, reduction_results()$method)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "Density Plot",
                       label = "Density Plot", source_module = "Dimensionality Reduction")
      }
      p
    })
  })
}

# -----------------------------------------------------------------------------
# 3D EMBEDDINGS UI & SERVER (Enhanced)
# -----------------------------------------------------------------------------

embedding3dUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("3D Embeddings"),
          p("Create interactive 3D visualizations for exploring cell relationships."),
          
          fluidRow(
            column(4,
              selectInput(ns("method3d"), "Method",
                choices = c("UMAP" = "umap", "t-SNE" = "tsne", "PCA" = "pca"),
                selected = "umap")
            ),
            column(4,
              selectInput(ns("input_reduction_3d"), "Input Reduction", choices = NULL)
            ),
            column(4,
              numericInput(ns("dims3d"), "Dimensions to Use", value = 30, min = 5, max = 100)
            )
          ),
          
          # 3D-specific parameters
          conditionalPanel(
            condition = "input.method3d == 'umap'",
            ns = ns,
            fluidRow(
              column(6,
                numericInput(ns("n_neighbors_3d"), "Neighbors", value = 30, min = 5, max = 200)
              ),
              column(6,
                numericInput(ns("min_dist_3d"), "Min Distance", value = 0.3, min = 0.001, max = 2, step = 0.1)
              )
            )
          ),
          
          conditionalPanel(
            condition = "input.method3d == 'tsne'",
            ns = ns,
            fluidRow(
              column(6,
                numericInput(ns("perplexity_3d"), "Perplexity", value = 30, min = 5, max = 100)
              ),
              column(6,
                numericInput(ns("theta_3d"), "Theta", value = 0.5, min = 0, max = 1, step = 0.1)
              )
            )
          ),
          
          fluidRow(
            column(4,
              selectInput(ns("color_by_3d"), "Color by", choices = NULL)
            ),
            column(4,
              numericInput(ns("seed_3d"), "Random Seed", value = 42, min = 1, max = 10000)
            ),
            column(4,
              div(style = "margin-top: 25px;",
                actionButton(ns("run_3d"), "Generate 3D Embedding", 
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    # Results
    fluidRow(
      column(12,
        uiOutput(ns("embed3d_results_ui"))
      )
    )
  )
}

embedding3dServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_3d_results <- reactiveVal(FALSE)
    embedding_3d_results <- reactiveVal(NULL)
    
    output$embed3d_results_ui <- renderUI({
      req(isTRUE(show_3d_results()))
      ns <- session$ns
      wellPanel(
        h4("3D Visualization"),

        fluidRow(
          column(9,
            plotlyOutput(ns("plot3d"), height = "700px")
          ),
          column(3,
            h5("3D Plot Controls"),
            sliderInput(ns("point_size_3d"), "Point Size", min = 1, max = 10, value = 3),
            sliderInput(ns("point_opacity_3d"), "Opacity", min = 0.1, max = 1, value = 0.8, step = 0.1),
            selectInput(ns("camera_projection"), "Projection",
              choices = c("perspective", "orthographic"),
              selected = "perspective"),
            checkboxInput(ns("show_axes"), "Show axes", TRUE),
            br(),
            downloadButton(ns("download_3d"), "Save 3D Plot", class = "btn-success")
          )
        )
      )
    })
    
    # Update UI choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      
      # Available reductions
      reductions <- names(obj@reductions)
      updateSelectInput(session, "input_reduction_3d", 
                       choices = reductions,
                       selected = if ("pca" %in% reductions) "pca" else reductions[1])
      
      # Color options
      meta_cols <- c("None" = "none", colnames(obj@meta.data))
      updateSelectInput(session, "color_by_3d", choices = meta_cols)
    })
    
    # Generate 3D embedding
    observeEvent(input$run_3d, {
      req(values$active_dataset, input$input_reduction_3d)
      job_id <- begin_dimred_job(values, "run_3d_embedding", meta = list(
        method = input$method3d,
        input_reduction = input$input_reduction_3d,
        dims = input$dims3d
      ))
      tryCatch({
        obj <- values$active_dataset
        if (!input$input_reduction_3d %in% names(obj@reductions)) {
          showNotification("Selected input reduction not available", type = "error")
          end_dimred_job(values, job_id, status = "failed", message = "missing input reduction")
          return()
        }
        if (!validate_dimred_contract(values, obj, stage = "dimensionality_reduction", require_reductions = input$input_reduction_3d)) {
          end_dimred_job(values, job_id, status = "failed", message = "contract validation failed")
          return()
        }

        withProgress(message = paste("Generating 3D", toupper(input$method3d), "..."), value = 0, {
          result <- safely_run_3d_embedding(obj, input)

          if (result$success) {
            values$active_dataset <- result$object
            embedding_3d_results(result)
            show_3d_results(TRUE)
            checkpoint_dimred(
              values,
              action = paste0("run_3d_", input$method3d),
              obj = result$object,
              params = list(
                method = input$method3d,
                input_reduction = input$input_reduction_3d,
                dims = input$dims3d,
                seed = input$seed_3d
              ),
              method = paste0("Seurat::Run", toupper(input$method3d), " (3D)")
            )
            end_dimred_job(values, job_id, status = "done", message = "3D embedding completed")

            showNotification(paste("3D", toupper(input$method3d), "completed"), type = "message")
          } else {
            end_dimred_job(values, job_id, status = "failed", message = result$error)
            showNotification(paste("3D embedding failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        end_dimred_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("3D embedding error:", e$message), type = "error", duration = 10)
      })
    })
    
    # 3D plot
    output$plot3d <- renderPlotly({
      req(embedding_3d_results(), values$active_dataset)
      p <- create_3d_plot(values$active_dataset, embedding_3d_results(), input)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "3D Reduction",
                       label = "3D Reduction", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    # Download handler
    output$download_3d <- downloadHandler(
      filename = function() paste0("3d_embedding_", input$method3d, "_", Sys.Date(), ".html"),
      content = function(file) {
        plot <- create_3d_plot(values$active_dataset, embedding_3d_results(), input)
        htmlwidgets::saveWidget(plot, file)
      }
    )
  })
}

# -----------------------------------------------------------------------------
# DIFFUSION MAPS UI & SERVER (Enhanced)
# -----------------------------------------------------------------------------

diffusionUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Diffusion Maps"),
          p("Compute diffusion maps for trajectory analysis and manifold learning."),
          
          fluidRow(
            column(4,
              numericInput(ns("dm_k"), "k (neighbors)", value = 30, min = 5, max = 200)
            ),
            column(4,
              numericInput(ns("dm_ndc"), "Diffusion Components", value = 10, min = 2, max = 50)
            ),
            column(4,
              numericInput(ns("dm_sigma"), "Sigma (bandwidth)", value = NA, min = 0.1, max = 10, step = 0.1)
            )
          ),
          
          fluidRow(
            column(4,
              selectInput(ns("dm_distance"), "Distance Metric",
                choices = c("euclidean", "cosine", "correlation", "manhattan"),
                selected = "euclidean")
            ),
            column(4,
              checkboxInput(ns("dm_density_norm"), "Density Normalization", TRUE)
            ),
            column(4,
              selectInput(ns("dm_features"), "Features to Use",
                choices = c("Variable features" = "hvg", "All features" = "all"),
                selected = "hvg")
            )
          ),
          
          fluidRow(
            column(4,
              numericInput(ns("dm_n_local"), "Local transitions", value = NULL, min = 1, max = 10)
            ),
            column(4,
              checkboxInput(ns("dm_censor_val"), "Censor transitions", FALSE)
            ),
            column(4,
              div(style = "margin-top: 25px;",
                actionButton(ns("run_dm"), "Run Diffusion Maps", 
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    # Results
    fluidRow(
      column(12,
        uiOutput(ns("dm_results_ui"))
      )
    )
  )
}

diffusionServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_dm_results <- reactiveVal(FALSE)
    dm_results <- reactiveVal(NULL)
    dpt_results <- reactiveVal(NULL)
    
    output$dm_results_ui <- renderUI({
      req(isTRUE(show_dm_results()))
      ns <- session$ns
      wellPanel(
        h4("Diffusion Maps Results"),

        tabsetPanel(
          tabPanel("Diffusion Components",
            fluidRow(
              column(9,
                plotlyOutput(ns("dm_plot"), height = "600px")
              ),
              column(3,
                h5("Plot Options"),
                selectInput(ns("dm_x"), "X-axis", choices = NULL),
                selectInput(ns("dm_y"), "Y-axis", choices = NULL),
                selectInput(ns("color_by_dm"), "Color by", choices = NULL),
                br(),
                downloadButton(ns("download_dm"), "Download DM", class = "btn-success")
              )
            )
          ),

          tabPanel("Eigenvalues",
            fluidRow(
              column(8,
                plotlyOutput(ns("eigenvalue_plot"), height = "450px")
              ),
              column(4,
                h5("Component Statistics"),
                DT::dataTableOutput(ns("eigenvalue_table"))
              )
            )
          ),

          tabPanel("Pseudotime",
            fluidRow(
              column(4,
                h5("Diffusion Pseudotime"),
                selectInput(ns("dpt_root"), "Root cells", choices = NULL),
                actionButton(ns("compute_dpt"), "Compute DPT", class = "btn-info")
              ),
              column(8,
                plotlyOutput(ns("dpt_plot"), height = "500px")
              )
            )
          )
        )
      )
    })
    
    # Update UI choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      
      # Color options
      meta_all <- colnames(obj@meta.data)
updateSelectInput(session, "color_by_dm", choices = c("None" = "none", stats::setNames(meta_all, meta_all)))
# Root cell options (clusters/cell types)
factor_cols <- meta_all[sapply(obj@meta.data[meta_all], function(x) is.factor(x) || is.character(x))]
updateSelectInput(session, "dpt_root", choices = stats::setNames(factor_cols, factor_cols))
    })
    
    # Update DC choices after diffusion maps
    observeEvent(dm_results(), {
      tryCatch({
        req(dm_results())
        dc_names <- paste0("DC_", 1:dm_results()$n_components)
        updateSelectInput(session, "dm_x", choices = dc_names, selected = "DC_1")
        updateSelectInput(session, "dm_y", choices = dc_names, selected = "DC_2")
      }, error = function(e) NULL)
    })
    
    # Run diffusion maps
    observeEvent(input$run_dm, {
      req(values$active_dataset)
      job_id <- begin_dimred_job(values, "run_diffusion_map", meta = list(
        k = input$dm_k,
        n_components = input$dm_ndc,
        features = input$dm_features
      ))
      if (!validate_dimred_contract(values, values$active_dataset, stage = "dimensionality_reduction")) {
        end_dimred_job(values, job_id, status = "failed", message = "contract validation failed")
        return()
      }
      tryCatch({
        withProgress(message = "Computing diffusion maps...", value = 0, {
          obj <- values$active_dataset

          if (exists("runtime_budget_guard", mode = "function")) {
            budget <- runtime_budget_guard(obj, task = "diffusion_maps")
            if (!isTRUE(budget$ok)) {
              showNotification(paste(budget$message, budget$suggestion), type = "error", duration = 10)
              end_dimred_job(values, job_id, status = "failed", message = budget$message)
              return()
            }
          }

          incProgress(0.2, detail = "Preparing gene expression matrix...")

          result <- safely_run_diffusion_maps(obj, input)

          if (result$success) {
            values$active_dataset <- result$object
            dm_results(result)
            show_dm_results(TRUE)
            checkpoint_dimred(
              values,
              action = "run_diffusion_map",
              obj = result$object,
              params = list(k = input$dm_k, n_components = input$dm_ndc, features = input$dm_features),
              method = "destiny::DiffusionMap"
            )
            end_dimred_job(values, job_id, status = "done", message = "Diffusion maps completed")

            showNotification("Diffusion maps completed", type = "message")
          } else {
            end_dimred_job(values, job_id, status = "failed", message = result$error)
            showNotification(paste("Diffusion maps failed:", result$error), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        end_dimred_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("Diffusion maps error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Compute diffusion pseudotime
    observeEvent(input$compute_dpt, {
      req(dm_results(), input$dpt_root)
      job_id <- begin_dimred_job(values, "compute_dpt", meta = list(root = input$dpt_root))
      tryCatch({
        withProgress(message = "Computing diffusion pseudotime...", value = 0, {
          result <- safely_compute_dpt(values$active_dataset, dm_results(), input$dpt_root)

          if (result$success) {
            values$active_dataset <- result$object
            dpt_results(result)
            checkpoint_dimred(
              values,
              action = "compute_dpt",
              obj = result$object,
              params = list(root = input$dpt_root),
              method = "Diffusion pseudotime"
            )
            end_dimred_job(values, job_id, status = "done", message = "DPT computed")
            showNotification("Diffusion pseudotime computed", type = "message")
          } else {
            end_dimred_job(values, job_id, status = "failed", message = result$error)
            showNotification(paste("DPT failed:", result$error), type = "error")
          }
        })
      }, error = function(e) {
        end_dimred_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("DPT error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Results displays
    output$dm_plot <- renderPlotly({
      req(dm_results(), values$active_dataset)
      p <- create_diffusion_plot(values$active_dataset, input$dm_x, input$dm_y, input$color_by_dm)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "Diffusion Map",
                       label = "Diffusion Map", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    output$eigenvalue_plot <- renderPlotly({
      req(dm_results())
      p <- create_eigenvalue_plot(dm_results())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "Eigenvalue Plot",
                       label = "Eigenvalue Plot", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    output$eigenvalue_table <- DT::renderDataTable({
      req(dm_results())
      create_eigenvalue_table(dm_results())
    }, options = list(dom = 't', pageLength = 10))
    
    output$dpt_plot <- renderPlotly({
      req(dpt_results(), values$active_dataset)
      p <- create_dpt_plot(values$active_dataset, dpt_results())
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Dimensionality Reduction", figure_type = "DPT Plot",
                       label = "DPT Plot", source_module = "Dimensionality Reduction")
      }
      p
    })
    
    # Download handler
    output$download_dm <- downloadHandler(
      filename = function() paste0("diffusion_maps_", Sys.Date(), ".csv"),
      content = function(file) {
        req(dm_results(), values$active_dataset)
        dm_embedding <- Embeddings(values$active_dataset, "diffmap")
        write.csv(dm_embedding, file, row.names = TRUE)
      }
    )
  })
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Issue 4.1: Auto-detect PCA elbow via second derivative of scree plot
detect_elbow_point <- function(stdev) {
  # Compute second derivative of the standard deviations
  # The elbow is where the second derivative is maximized (greatest curvature change)
  if (length(stdev) < 4) return(2)

  d1 <- diff(stdev)          # first derivative

  d2 <- diff(d1)             # second derivative

  # The elbow point is where the second derivative is largest (most positive),

  # meaning the rate of decrease slows down the most
  elbow_idx <- which.max(d2) + 1  # +1 to offset diff() index shift

  # Clamp to reasonable range
  elbow_idx <- max(2, min(elbow_idx, length(stdev) - 1))

  return(elbow_idx)
}

safely_run_pca <- function(seurat_obj, params, features) {
  tryCatch({
    set.seed(params$pca_seed)
    
    # Ensure data is scaled if needed (v4/v5 compatible check)
    has_scaled <- tryCatch({
      assay_obj <- seurat_obj@assays[[params$pca_assay]]
      if (inherits(assay_obj, "Assay5")) {
        "scale.data" %in% Layers(assay_obj)
      } else {
        "scale.data" %in% slotNames(assay_obj)
      }
    }, error = function(e) FALSE)
    if (!has_scaled) {
      seurat_obj <- ScaleData(seurat_obj, features = features, verbose = FALSE)
    }

    # Memory budget guard for PCA
    mem_check <- memory_budget_guard(seurat_obj)
    if (!mem_check$ok) {
      warning("Memory budget exceeded: ", mem_check$message, " ", mem_check$suggestion)
    }

    # Run PCA (with approx fallback for version compatibility)
    use_approx <- isTRUE(params$use_irlba)
    pca_obj <- tryCatch(
      RunPCA(
        seurat_obj,
        features = features,
        npcs = params$n_pcs,
        rev.pca = FALSE,
        weight.by.var = TRUE,
        verbose = FALSE,
        reduction.name = "pca",
        approx = use_approx
      ),
      error = function(e) {
        # Fallback: retry without approx if it failed
        if (use_approx) {
          warning("Approximate PCA failed, retrying with exact PCA: ", e$message)
          RunPCA(seurat_obj, features = features, npcs = params$n_pcs,
                 verbose = FALSE, reduction.name = "pca")
        } else {
          stop(e)
        }
      }
    )
    
    # Calculate additional statistics
    stdev <- pca_obj@reductions$pca@stdev
    total_var <- sum(stdev^2)
    var_explained <- stdev^2 / total_var * 100
    cumvar_explained <- cumsum(var_explained)
    
    list(
      success = TRUE,
      object = pca_obj,
      n_pcs = params$n_pcs,
      var_explained = var_explained,
      cumvar_explained = cumvar_explained,
      stdev = stdev,
      features_used = features,
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      n_pcs = 0,
      var_explained = NULL,
      cumvar_explained = NULL,
      stdev = NULL,
      features_used = NULL,
      error = e$message
    )
  })
}

safely_run_nonlinear_reduction <- function(seurat_obj, params) {
  tryCatch({
    set.seed(params$reduction_seed)

    # For very large datasets, notify user about expected time
    n_cells <- ncol(seurat_obj)
    if (n_cells > 100000) {
      message("Large dataset: ", format(n_cells, big.mark = ","), " cells. UMAP/tSNE may take several minutes.")
    }

    method <- params$reduction_method
    input_reduction <- params$reduction_input

    # Clamp dims to what the input reduction actually has
    avail_dims <- ncol(Embeddings(seurat_obj, input_reduction))
    n_dims <- min(params$dims_to_use, avail_dims)
    dims <- 1:n_dims

    if (method == "umap") {
      # Enhanced UMAP with additional parameters
      umap_config <- list(
        n_neighbors = params$n_neighbors,
        min_dist = params$min_dist,
        spread = params$spread,
        metric = params$umap_metric,
        n_components = 2,
        return.model = params$return_model
      )
      
      # Add optional parameters if specified
      if (!is.null(params$a) && !is.na(params$a)) umap_config$a <- params$a
      if (!is.null(params$b) && !is.na(params$b)) umap_config$b <- params$b
      if (params$densmap) umap_config$densmap <- TRUE
      
      reduced_obj <- RunUMAP(
        seurat_obj,
        reduction = input_reduction,
        dims = dims,
        umap.method = "uwot",
        reduction.name = "umap",
        verbose = FALSE,
        n.neighbors = umap_config$n_neighbors,
        min.dist = umap_config$min_dist,
        spread = umap_config$spread,
        metric = umap_config$metric,
        n.components = 2,
        return.model = umap_config$return.model
      )
      
    } else if (method == "tsne") {
      # Enhanced t-SNE
      reduced_obj <- RunTSNE(
        seurat_obj,
        reduction = input_reduction,
        dims = dims,
        perplexity = params$perplexity,
        theta = params$theta,
        max_iter = params$max_iter,
        eta = params$learning_rate,
        dim.embed = 2,
        reduction.name = "tsne",
        check_duplicates = FALSE,
        verbose = FALSE
      )
    } else {
      stop("Unknown reduction method: ", method)
    }
    
    # Calculate quality metrics
    embedding <- Embeddings(reduced_obj, method)
    quality_metrics <- calculate_embedding_quality(seurat_obj, embedding, input_reduction, dims)
    
    list(
      success = TRUE,
      object = reduced_obj,
      method = method,
      input_reduction = input_reduction,
      dims_used = dims,
      parameters = params,
      quality_metrics = quality_metrics,
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      method = params$reduction_method,
      input_reduction = params$reduction_input,
      dims_used = NULL,
      parameters = params,
      quality_metrics = NULL,
      error = e$message
    )
  })
}

safely_run_3d_embedding <- function(seurat_obj, params) {
  tryCatch({
    set.seed(params$seed_3d)
    
    method <- params$method3d
    input_reduction <- params$input_reduction_3d
    dims <- 1:params$dims3d
    
    if (method == "umap") {
      reduced_obj <- RunUMAP(
        seurat_obj,
        reduction = input_reduction,
        dims = dims,
        umap.method = "uwot",
        n.neighbors = params$n_neighbors_3d,
        min.dist = params$min_dist_3d,
        n.components = 3,
        reduction.name = paste0(method, "3d"),
        verbose = FALSE
      )
    } else if (method == "tsne") {
      reduced_obj <- RunTSNE(
        seurat_obj,
        reduction = input_reduction,
        dims = dims,
        perplexity = params$perplexity_3d,
        theta = params$theta_3d,
        dim.embed = 3,
        reduction.name = paste0(method, "3d"),
        check_duplicates = FALSE,
        verbose = FALSE
      )
    } else if (method == "pca") {
      # Use first 3 PCs for 3D visualization
      reduced_obj <- seurat_obj
      if (!"pca" %in% names(reduced_obj@reductions)) {
        # Ensure scale.data exists before RunPCA
        has_scaled <- tryCatch({
          sd <- safe_get_assay_data(reduced_obj, "scale.data")
          !is.null(sd) && nrow(sd) > 0
        }, error = function(e) FALSE)
        if (!has_scaled) {
          features <- VariableFeatures(reduced_obj)
          if (length(features) == 0) features <- rownames(reduced_obj)
          reduced_obj <- ScaleData(reduced_obj, features = features, verbose = FALSE)
        }
        reduced_obj <- RunPCA(reduced_obj, npcs = max(3, params$dims3d), verbose = FALSE)
      }
    } else {
      stop("Unknown 3D method: ", method)
    }
    
    list(
      success = TRUE,
      object = reduced_obj,
      method = method,
      input_reduction = input_reduction,
      dims_used = dims,
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      method = params$method3d,
      input_reduction = params$input_reduction_3d,
      dims_used = NULL,
      error = e$message
    )
  })
}

safely_run_diffusion_maps <- function(seurat_obj, params) {
  tryCatch({
    if (!requireNamespace("destiny", quietly = TRUE)) {
      stop("Package 'destiny' required for diffusion maps")
    }
    
    # Get features to use
    if (params$dm_features == "hvg") {
      features <- VariableFeatures(seurat_obj)
      if (length(features) == 0) {
        features <- rownames(seurat_obj)
      }
    } else {
      features <- rownames(seurat_obj)
    }
    
    # Subset to top variable genes for performance
    if (length(features) > 3000) {
      if (params$dm_features == "hvg") {
        features <- head(features, 3000)
      } else {
        # Select most variable genes
        hvg_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000, verbose = FALSE)
        features <- VariableFeatures(hvg_obj)
      }
    }
    
    # Issue 4.3: For SCT objects, switch to RNA assay log-normalized data for DiffusionMap
    # DiffusionMap expects log-normalized counts; SCT residuals can contain negative values
    # and have different distributional properties that are not suitable for diffusion distance calculations.
    dm_assay <- DefaultAssay(seurat_obj)
    if (tolower(dm_assay) == "sct") {
      if ("RNA" %in% names(seurat_obj@assays)) {
        message("DiffusionMap: DefaultAssay is SCT. Switching to RNA assay log-normalized data. ",
                "SCT residuals are not suitable for diffusion distance computation because they ",
                "can contain negative values and have non-standard variance structure.")
        dm_assay <- "RNA"
      } else {
        warning("DefaultAssay is SCT but RNA assay not found. Proceeding with SCT data.")
      }
    }

    # Get expression matrix from the appropriate assay
    expr_matrix <- as.matrix(safe_get_assay_data(seurat_obj, "data", assay = dm_assay)[features, ])

    # Create diffusion map
    dm_params <- list(
      k = params$dm_k,
      n_eigs = params$dm_ndc,
      distance = params$dm_distance,
      density_norm = params$dm_density_norm
    )
    
    if (!is.null(params$dm_sigma) && !is.na(params$dm_sigma)) {
      dm_params$sigma <- params$dm_sigma
    }
    
    dm <- do.call(destiny::DiffusionMap, c(list(data = t(expr_matrix)), dm_params))
    
    # Extract embedding and create Seurat reduction
    embedding <- as.matrix(destiny::eigenvectors(dm))
    colnames(embedding) <- paste0("DC_", 1:ncol(embedding))
    rownames(embedding) <- colnames(seurat_obj)
    
    # Add to Seurat object
    dm_reduction <- CreateDimReducObject(
      embeddings = embedding,
      key = "DC_",
      assay = DefaultAssay(seurat_obj),
      stdev = destiny::eigenvalues(dm)
    )
    seurat_obj[["diffmap"]] <- dm_reduction
    
    # Extract additional results
    eigenvalues <- destiny::eigenvalues(dm)
    
    list(
      success = TRUE,
      object = seurat_obj,
      dm = dm,
      embedding = embedding,
      eigenvalues = eigenvalues,
      n_components = params$dm_ndc,
      features_used = features,
      parameters = params,
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      dm = NULL,
      embedding = NULL,
      eigenvalues = NULL,
      n_components = 0,
      features_used = NULL,
      parameters = params,
      error = e$message
    )
  })
}

safely_compute_dpt <- function(seurat_obj, dm_results, root_column) {
  tryCatch({
    if (!requireNamespace("destiny", quietly = TRUE)) {
      stop("Package 'destiny' required for diffusion pseudotime")
    }
    
    # Get root cells
    root_groups <- unique(seurat_obj@meta.data[[root_column]])
    if (length(root_groups) == 0) {
      stop("No groups found in root column")
    }
    
    # Use first group as root (could be made interactive)
    root_group <- root_groups[1]
    root_cells <- which(seurat_obj@meta.data[[root_column]] == root_group)
    
    if (length(root_cells) == 0) {
      stop("No root cells found")
    }

    # Validate root cell index is within diffusion map bounds
    n_dm_cells <- nrow(dm_results$dm)
    root_idx <- root_cells[1]
    if (root_idx > n_dm_cells) {
      # Diffusion map may have been computed on a subset — find matching cell
      dm_cells <- rownames(dm_results$dm@eigenvectors)
      if (!is.null(dm_cells)) {
        obj_cells <- colnames(seurat_obj)[root_cells]
        match_idx <- which(dm_cells %in% obj_cells)
        if (length(match_idx) > 0) {
          root_idx <- match_idx[1]
        } else {
          root_idx <- 1L
          warning("Root cell not found in diffusion map. Using first cell as root.")
        }
      } else {
        root_idx <- min(root_idx, n_dm_cells)
      }
    }

    # Compute DPT — tips requires integer index into diffusion map
    dpt <- destiny::DPT(dm_results$dm, tips = as.integer(root_idx))
    
    # Extract pseudotime
    pseudotime <- destiny::dpt(dpt)[, 1]  # First branch
    names(pseudotime) <- colnames(seurat_obj)
    
    # Add to Seurat metadata
    seurat_obj$diffusion_pseudotime <- pseudotime
    
    list(
      success = TRUE,
      object = seurat_obj,
      dpt = dpt,
      pseudotime = pseudotime,
      root_group = root_group,
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      object = NULL,
      dpt = NULL,
      pseudotime = NULL,
      root_group = NULL,
      error = e$message
    )
  })
}

# Quality metrics calculation
calculate_embedding_quality <- function(seurat_obj, embedding, input_reduction, dims) {
  tryCatch({
    # Get original high-dimensional data
    original_data <- Embeddings(seurat_obj, input_reduction)[, dims]

    # Calculate preservation of local neighborhoods
    n_neighbors <- min(30, nrow(embedding) - 1)

    # Issue 4.4: Use stratified sampling by cluster for quality assessment
    local_preservation <- calculate_local_preservation(
      original_data, embedding, n_neighbors, seurat_obj = seurat_obj)

    list(
      local_preservation = local_preservation,
      n_neighbors_used = n_neighbors,
      embedding_variance = apply(embedding, 2, var)
    )

  }, error = function(e) {
    list(
      local_preservation = NA,
      n_neighbors_used = NA,
      embedding_variance = c(NA, NA)
    )
  })
}

calculate_local_preservation <- function(original_data, embedding, k = 30, seurat_obj = NULL) {
  # Simplified local preservation calculation
  # This measures how well local neighborhoods are preserved

  if (nrow(original_data) < k + 1) {
    return(NA)
  }

  # Issue 4.4: Stratified sampling by cluster ensures small clusters are represented
  set.seed(42)
  sample_idx <- NULL
  if (!is.null(seurat_obj)) {
    # Try to find a cluster column
    cluster_col <- NULL
    for (cand in c("seurat_clusters", "cluster", "Cluster", "cell_type", "CellType")) {
      if (cand %in% colnames(seurat_obj@meta.data)) {
        cluster_col <- cand
        break
      }
    }
    if (!is.null(cluster_col)) {
      clusters <- seurat_obj@meta.data[[cluster_col]]
      sample_idx <- unlist(lapply(unique(clusters), function(cl) {
        cl_idx <- which(clusters == cl)
        n_take <- min(50, length(cl_idx))
        sample(cl_idx, n_take)
      }))
    }
  }
  # Fallback to random sampling if no clusters found
  if (is.null(sample_idx)) {
    n_sample <- min(500, nrow(original_data))
    sample_idx <- sample(nrow(original_data), n_sample)
  }

  orig_sample <- original_data[sample_idx, ]
  embed_sample <- embedding[sample_idx, ]
  
  # Calculate average preservation (simplified)
  tryCatch({
    orig_dist <- as.matrix(dist(orig_sample))
    embed_dist <- as.matrix(dist(embed_sample))
    
    # For each point, find k nearest neighbors in original space
    # and measure how many are still neighbors in embedding
    preservation_scores <- sapply(1:nrow(orig_sample), function(i) {
      orig_neighbors <- order(orig_dist[i, ])[2:(k+1)]  # Exclude self
      embed_neighbors <- order(embed_dist[i, ])[2:(k+1)]
      
      length(intersect(orig_neighbors, embed_neighbors)) / k
    })
    
    mean(preservation_scores)
    
  }, error = function(e) {
    NA
  })
}

# Display creation functions
create_enhanced_elbow_plot <- function(pca_results) {
  pc_data <- data.frame(
    PC = 1:length(pca_results$var_explained),
    Variance_Explained = pca_results$var_explained,
    Cumulative_Variance = pca_results$cumvar_explained
  )

  # Issue 4.1: Detect elbow point for annotation
  elbow_pc <- detect_elbow_point(pca_results$stdev)

  p1 <- plot_ly(
    data = pc_data,
    x = ~PC,
    y = ~Variance_Explained,
    type = "scatter",
    mode = "lines+markers",
    name = "Individual",
    line = list(color = "blue"),
    marker = list(color = "blue")
  )

  p1 <- p1 %>% add_trace(
    y = ~Cumulative_Variance,
    name = "Cumulative",
    line = list(color = "red"),
    marker = list(color = "red")
  )

  p1 %>% layout(
    autosize = TRUE,
    title = paste0("PCA Elbow Plot (suggested elbow: PC ", elbow_pc, ")"),
    xaxis = list(title = "Principal Component"),
    yaxis = list(title = "Variance Explained (%)"),
    shapes = list(
      list(
        type = "line",
        x0 = elbow_pc, x1 = elbow_pc,
        y0 = 0, y1 = 1,
        yref = "paper",
        line = list(color = "green", width = 2, dash = "dash")
      )
    ),
    annotations = list(
      list(
        x = elbow_pc,
        y = 1,
        yref = "paper",
        text = paste0("Elbow: PC ", elbow_pc),
        showarrow = TRUE,
        arrowhead = 2,
        ax = 40,
        ay = -30,
        font = list(color = "green", size = 12)
      )
    )
  )
}

create_pc_stats_table <- function(pca_results) {
  pc_stats <- data.frame(
    PC = paste0("PC", 1:length(pca_results$var_explained)),
    Variance_Explained = round(pca_results$var_explained, 2),
    Cumulative_Variance = round(pca_results$cumvar_explained, 2),
    Standard_Deviation = round(pca_results$stdev, 3)
  )
  
  DT::datatable(head(pc_stats, 20), 
                options = list(dom = 't', pageLength = 15, scrollY = "300px"))
}

# (PC heatmaps are now rendered as individual interactive plotly panels in the server)

create_pca_scatter_plot <- function(seurat_obj, pc_x, pc_y, color_by) {
  embedding <- Embeddings(seurat_obj, "pca")
  
  pc_x_num <- as.numeric(gsub("PC_", "", pc_x))
  pc_y_num <- as.numeric(gsub("PC_", "", pc_y))
  
  plot_data <- data.frame(
    PC_X = embedding[, pc_x_num],
    PC_Y = embedding[, pc_y_num],
    Cell = rownames(embedding)
  )
  
  # Add color information
  if (color_by != "none" && color_by %in% colnames(seurat_obj@meta.data)) {
    plot_data$Color <- seurat_obj@meta.data[[color_by]]
    
    p <- plot_ly(
      data = plot_data,
      x = ~PC_X,
      y = ~PC_Y,
      color = ~Color,
      type = "scatter",
      mode = "markers",
      marker = list(size = 4, opacity = 0.8),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", pc_x, ": %{x}<br>", pc_y, ": %{y}<extra></extra>")
    )
  } else {
    p <- plot_ly(
      data = plot_data,
      x = ~PC_X,
      y = ~PC_Y,
      type = "scatter",
      mode = "markers",
      marker = list(size = 4, opacity = 0.8, color = "blue"),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", pc_x, ": %{x}<br>", pc_y, ": %{y}<extra></extra>")
    )
  }
  
  p %>% layout(
    autosize = TRUE,
    title = paste(pc_x, "vs", pc_y),
    xaxis = list(title = pc_x),
    yaxis = list(title = pc_y)
  )
}

create_loadings_plot <- function(seurat_obj, loading_pc, top_n, loading_type) {
  loadings <- Loadings(seurat_obj, reduction = "pca")
  pc_num <- as.numeric(gsub("PC_", "", loading_pc))
  
  if (pc_num > ncol(loadings)) {
    return(plotly_empty() %>% layout(title = "PC not available"))
  }
  
  pc_loadings <- loadings[, pc_num]
  names(pc_loadings) <- rownames(loadings)
  
  # Filter based on type
  if (loading_type == "pos") {
    pc_loadings <- pc_loadings[pc_loadings > 0]
    pc_loadings <- sort(pc_loadings, decreasing = TRUE)
  } else if (loading_type == "neg") {
    pc_loadings <- pc_loadings[pc_loadings < 0]
    pc_loadings <- sort(pc_loadings, decreasing = FALSE)
  } else {
    pc_loadings <- pc_loadings[order(abs(pc_loadings), decreasing = TRUE)]
  }
  
  # Get top loadings
  top_loadings <- head(pc_loadings, top_n)
  
  loading_data <- data.frame(
    Gene = names(top_loadings),
    Loading = as.numeric(top_loadings)
  )
  loading_data$Gene <- factor(loading_data$Gene, levels = rev(loading_data$Gene))
  
  plot_ly(
    data = loading_data,
    x = ~Loading,
    y = ~Gene,
    type = "bar",
    orientation = "h",
    marker = list(color = ifelse(loading_data$Loading > 0, "red", "blue"))
  ) %>% layout(
    autosize = TRUE,
    title = paste("Top Loadings for", loading_pc),
    xaxis = list(title = "Loading"),
    yaxis = list(title = "Gene")
  )
}

create_enhanced_reduction_plot <- function(seurat_obj, reduction_results, input_params) {
  method <- reduction_results$method
  embedding <- Embeddings(seurat_obj, method)
  
  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Cell = rownames(embedding)
  )
  
  # Add color information
  color_by <- input_params$color_by_reduction
  if (color_by != "none" && color_by %in% colnames(seurat_obj@meta.data)) {
    plot_data$Color <- seurat_obj@meta.data[[color_by]]
    
    # Determine if continuous or discrete
    if (is.numeric(plot_data$Color)) {
      p <- plot_ly(
        data = plot_data,
        x = ~Dim1,
        y = ~Dim2,
        color = ~Color,
        type = "scatter",
        mode = "markers",
        marker = list(size = input_params$point_size, opacity = input_params$point_alpha),
        text = ~Cell,
        hovertemplate = paste0("Cell: %{text}<br>", toupper(method), "1: %{x}<br>", toupper(method), "2: %{y}<extra></extra>")
      ) %>% colorbar(title = color_by)
    } else {
      # Discrete coloring
      colors <- get_color_palette(length(unique(plot_data$Color)))
      p <- plot_ly(
        data = plot_data,
        x = ~Dim1,
        y = ~Dim2,
        color = ~Color,
        colors = colors,
        type = "scatter",
        mode = "markers",
        marker = list(size = input_params$point_size, opacity = input_params$point_alpha),
        text = ~Cell,
        hovertemplate = paste0("Cell: %{text}<br>", toupper(method), "1: %{x}<br>", toupper(method), "2: %{y}<extra></extra>")
      )
    }
  } else {
    p <- plot_ly(
      data = plot_data,
      x = ~Dim1,
      y = ~Dim2,
      type = "scatter",
      mode = "markers",
      marker = list(size = input_params$point_size, opacity = input_params$point_alpha, color = "blue"),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", toupper(method), "1: %{x}<br>", toupper(method), "2: %{y}<extra></extra>")
    )
  }
  
  p %>% layout(
    autosize = TRUE,
    title = paste(toupper(method), "Embedding"),
    xaxis = list(title = paste0(toupper(method), "_1")),
    yaxis = list(title = paste0(toupper(method), "_2")),
    showlegend = input_params$show_legend
  )
}

create_reduction_quality_table <- function(reduction_results) {
  quality_metrics <- reduction_results$quality_metrics
  
  quality_df <- data.frame(
    Metric = c("Local Preservation", "Embedding Variance (Dim1)", "Embedding Variance (Dim2)"),
    Value = c(
      round(quality_metrics$local_preservation, 3),
      round(quality_metrics$embedding_variance[1], 3),
      round(quality_metrics$embedding_variance[2], 3)
    )
  )
  
  DT::datatable(quality_df, options = list(dom = 't'))
}

create_reduction_params_table <- function(reduction_results, input_params) {
  method <- reduction_results$method
  
  if (method == "umap") {
    params_df <- data.frame(
      Parameter = c("Method", "Neighbors", "Min Distance", "Spread", "Metric"),
      Value = c("UMAP", input_params$n_neighbors, input_params$min_dist, 
               input_params$spread, input_params$umap_metric)
    )
  } else if (method == "tsne") {
    params_df <- data.frame(
      Parameter = c("Method", "Perplexity", "Theta", "Max Iterations", "Learning Rate"),
      Value = c("t-SNE", input_params$perplexity, input_params$theta, 
               input_params$max_iter, input_params$learning_rate)
    )
  } else {
    params_df <- data.frame(
      Parameter = c("Method", "Input Reduction", "Dimensions"),
      Value = c(toupper(method), reduction_results$input_reduction, length(reduction_results$dims_used))
    )
  }
  
  DT::datatable(params_df, options = list(dom = 't'))
}

create_density_plot <- function(seurat_obj, method) {
  if (!method %in% names(seurat_obj@reductions)) {
    return(plotly_empty() %>% layout(title = "Reduction not available"))
  }
  
  embedding <- Embeddings(seurat_obj, method)[, 1:2]
  
  plot_ly(
    x = embedding[, 1],
    y = embedding[, 2],
    type = "histogram2dcontour",
    colorscale = "Viridis"
  ) %>% layout(
    autosize = TRUE,
    title = paste(toupper(method), "Density Plot"),
    xaxis = list(title = paste0(toupper(method), "_1")),
    yaxis = list(title = paste0(toupper(method), "_2"))
  )
}

create_3d_plot <- function(seurat_obj, embedding_results, input_params) {
  method <- embedding_results$method
  reduction_name <- if (method == "pca") "pca" else paste0(method, "3d")
  
  if (!reduction_name %in% names(seurat_obj@reductions)) {
    return(plotly_empty() %>% layout(title = "3D embedding not available"))
  }
  
  embedding <- Embeddings(seurat_obj, reduction_name)[, 1:3]
  
  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Dim3 = embedding[, 3],
    Cell = rownames(embedding)
  )
  
  # Add color information
  color_by <- input_params$color_by_3d
  if (color_by != "none" && color_by %in% colnames(seurat_obj@meta.data)) {
    plot_data$Color <- seurat_obj@meta.data[[color_by]]
    
    p <- plot_ly(
      data = plot_data,
      x = ~Dim1,
      y = ~Dim2,
      z = ~Dim3,
      color = ~Color,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = input_params$point_size_3d, opacity = input_params$point_opacity_3d),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", 
                            toupper(method), "1: %{x}<br>", 
                            toupper(method), "2: %{y}<br>", 
                            toupper(method), "3: %{z}<extra></extra>")
    )
  } else {
    p <- plot_ly(
      data = plot_data,
      x = ~Dim1,
      y = ~Dim2,
      z = ~Dim3,
      type = "scatter3d",
      mode = "markers",
      marker = list(size = input_params$point_size_3d, opacity = input_params$point_opacity_3d, color = "blue"),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", 
                            toupper(method), "1: %{x}<br>", 
                            toupper(method), "2: %{y}<br>", 
                            toupper(method), "3: %{z}<extra></extra>")
    )
  }
  
  p %>% layout(
    autosize = TRUE,
    title = paste("3D", toupper(method), "Embedding"),
    scene = list(
      xaxis = list(title = paste0(toupper(method), "_1"), showgrid = input_params$show_axes),
      yaxis = list(title = paste0(toupper(method), "_2"), showgrid = input_params$show_axes),
      zaxis = list(title = paste0(toupper(method), "_3"), showgrid = input_params$show_axes),
      camera = list(projection = list(type = input_params$camera_projection))
    )
  )
}

create_diffusion_plot <- function(seurat_obj, dm_x, dm_y, color_by) {
  if (!"diffmap" %in% names(seurat_obj@reductions)) {
    return(plotly_empty() %>% layout(title = "Diffusion map not available"))
  }
  
  embedding <- Embeddings(seurat_obj, "diffmap")
  
  x_num <- as.numeric(gsub("DC_", "", dm_x))
  y_num <- as.numeric(gsub("DC_", "", dm_y))
  
  plot_data <- data.frame(
    DC_X = embedding[, x_num],
    DC_Y = embedding[, y_num],
    Cell = rownames(embedding)
  )
  
  # Add color information
  if (color_by != "none" && color_by %in% colnames(seurat_obj@meta.data)) {
    plot_data$Color <- seurat_obj@meta.data[[color_by]]
    
    p <- plot_ly(
      data = plot_data,
      x = ~DC_X,
      y = ~DC_Y,
      color = ~Color,
      type = "scatter",
      mode = "markers",
      marker = list(size = 4, opacity = 0.8),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", dm_x, ": %{x}<br>", dm_y, ": %{y}<extra></extra>")
    )
  } else {
    p <- plot_ly(
      data = plot_data,
      x = ~DC_X,
      y = ~DC_Y,
      type = "scatter",
      mode = "markers",
      marker = list(size = 4, opacity = 0.8, color = "blue"),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", dm_x, ": %{x}<br>", dm_y, ": %{y}<extra></extra>")
    )
  }
  
  p %>% layout(
    autosize = TRUE,
    title = paste("Diffusion Map:", dm_x, "vs", dm_y),
    xaxis = list(title = dm_x),
    yaxis = list(title = dm_y)
  )
}

create_eigenvalue_plot <- function(dm_results) {
  eigenvals <- dm_results$eigenvalues
  
  eigenval_data <- data.frame(
    Component = 1:length(eigenvals),
    Eigenvalue = eigenvals
  )
  
  plot_ly(
    data = eigenval_data,
    x = ~Component,
    y = ~Eigenvalue,
    type = "scatter",
    mode = "lines+markers",
    marker = list(color = "blue"),
    line = list(color = "blue")
  ) %>% layout(
    autosize = TRUE,
    title = "Diffusion Map Eigenvalues",
    xaxis = list(title = "Diffusion Component"),
    yaxis = list(title = "Eigenvalue", type = "log")
  )
}

create_eigenvalue_table <- function(dm_results) {
  eigenvals <- dm_results$eigenvalues
  
  eigenval_df <- data.frame(
    Component = paste0("DC_", 1:length(eigenvals)),
    Eigenvalue = round(eigenvals, 6),
    Log_Eigenvalue = round(log10(eigenvals), 3)
  )
  
  DT::datatable(eigenval_df, options = list(dom = 't', pageLength = 10))
}

create_dpt_plot <- function(seurat_obj, dpt_results) {
  if (!"diffmap" %in% names(seurat_obj@reductions)) {
    return(plotly_empty() %>% layout(title = "Diffusion map not available"))
  }
  
  embedding <- Embeddings(seurat_obj, "diffmap")[, 1:2]
  
  plot_data <- data.frame(
    DC1 = embedding[, 1],
    DC2 = embedding[, 2],
    Pseudotime = dpt_results$pseudotime,
    Cell = rownames(embedding)
  )
  
  plot_ly(
    data = plot_data,
    x = ~DC1,
    y = ~DC2,
    color = ~Pseudotime,
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.8),
    text = ~Cell,
    hovertemplate = "Cell: %{text}<br>DC1: %{x}<br>DC2: %{y}<br>Pseudotime: %{marker.color}<extra></extra>"
  ) %>% layout(
    autosize = TRUE,
    title = "Diffusion Pseudotime",
    xaxis = list(title = "DC1"),
    yaxis = list(title = "DC2")
  ) %>% colorbar(title = "Pseudotime")
}
