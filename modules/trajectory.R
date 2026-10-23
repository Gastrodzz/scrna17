# Enhanced Trajectory Analysis Module (Comprehensive pseudotime and lineage inference)
# modules/trajectory.R

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(slingshot)
  library(destiny)
  library(plotly)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(RColorBrewer)
  library(future)
  library(promises)
})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

.get_best_2d <- function(obj) {
  reds <- names(obj@reductions)
  if ("umap" %in% reds) return("umap")
  if ("tsne" %in% reds) return("tsne")
  if ("pca"  %in% reds) return("pca")
  NULL
}

.get_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    base <- RColorBrewer::brewer.pal(min(12, max(3, n)), "Set3")
    if (n > length(base)) grDevices::colorRampPalette(base)(n) else base[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, "Dynamic")
  }
}

.scale01 <- function(x) {
  if (all(is.na(x))) return(x)
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

.as_sce_with_reduction <- function(obj, red = "umap") {
  if (!(red %in% names(obj@reductions))) {
    stop(sprintf("Reduction '%s' not found in Seurat object.", red))
  }
  # Join split layers before conversion to avoid Seurat v5 issues
  for (assay_name in names(obj@assays)) {
    a <- obj[[assay_name]]
    if (inherits(a, "Assay5") && length(SeuratObject::Layers(a)) > 1) {
      tryCatch({ obj[[assay_name]] <- JoinLayers(a) }, error = function(e) NULL)
    }
  }
  sce <- as.SingleCellExperiment(obj)
  reducedDims(sce)[[toupper(red)]] <- Embeddings(obj, red)
  sce
}

.as_cds <- function(obj) {
  if (requireNamespace("SeuratWrappers", quietly = TRUE)) {
    return(SeuratWrappers::as.cell_data_set(obj))
  }
  expr <- safe_get_assay_data(obj, "counts")
  gene_meta <- data.frame(gene_short_name = rownames(expr), row.names = rownames(expr))
  cell_meta <- obj@meta.data
  monocle3::new_cell_data_set(
    expression_data = expr,
    cell_metadata = cell_meta,
    gene_metadata = gene_meta
  )
}

# =============================================================================
# MAIN TRAJECTORY UI
# =============================================================================

trajectoryUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Advanced Trajectory Analysis"),
          p("Comprehensive pseudotime analysis using multiple state-of-the-art methods with validation and comparison."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("trajectory"),

          # Global settings
          fluidRow(
            column(3,
              selectInput(ns("cluster_col"), "Cluster Column", choices = NULL)
            ),
            column(3,
              selectInput(ns("root_cluster"), "Root Cluster", choices = NULL)
            ),
            column(3,
              selectInput(ns("reduction_for_traj"), "Embedding",
                         choices = c("PCA (Recommended)" = "pca", "UMAP (visualization only)" = "umap", "t-SNE (visualization only)" = "tsne"), selected = "pca")
            ),
            column(3,
              numericInput(ns("n_dims_traj"), "Dimensions", 
                          value = 30, min = 5, max = 100)
            )
          ),
          
          # Advanced options
          fluidRow(
            column(3,
              checkboxInput(ns("write_to_meta"), "Save to metadata", TRUE)
            ),
            column(3,
              checkboxInput(ns("parallel_traj"), "Use parallel processing", TRUE)
            ),
            column(3,
              numericInput(ns("min_branch_len"), "Min branch length", 
                          value = 0.01, min = 0.001, max = 0.1, step = 0.01)
            ),
            column(3,
              actionButton(ns("refresh_choices"), "Refresh",
                          class = "btn-secondary")
            )
          ),

          # AI-powered trajectory tools
          fluidRow(
            column(12,
              div(style = "margin-bottom: 8px;",
                actionButton(ns("ai_root_suggest"), "AI Suggest Root Cell",
                           icon = icon("robot"), class = "ai-btn btn-sm"),
                actionButton(ns("ai_trajectory_narrative"), "AI Interpret Trajectory",
                           icon = icon("robot"), class = "ai-btn btn-sm")
              ),
              uiOutput(ns("ai_thread_root")),
              uiOutput(ns("ai_thread_trajectory"))
            )
          )
        )
      )
    ),

    # Method-specific panels
    tabsetPanel(
      id = ns("trajectory_methods"),
      
      # =========================== SLINGSHOT ===========================
      tabPanel("Slingshot",            icon = icon("route"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("Slingshot Configuration"),
              fluidRow(
                column(3,
                  numericInput(ns("sling_shrink"), "Curve shrinkage", 
                             value = 1.0, min = 0, max = 1, step = 0.1)
                ),
                column(3,
                  selectInput(ns("sling_lineage_mode"), "Lineage combination",
                            choices = c("First lineage" = "first",
                                      "Average lineages" = "average",
                                      "Principal curve" = "principal"),
                            selected = "principal")
                ),
                column(3,
                  checkboxInput(ns("sling_extend"), "Extend curves", TRUE)
                ),
                column(3,
                  numericInput(ns("sling_thresh"), "Convergence threshold",
                             value = 0.001, min = 0.0001, max = 0.01, step = 0.0001)
                )
              ),
              fluidRow(
                column(6,
                  checkboxInput(ns("sling_show_curves"), "Show curves", TRUE)
                ),
                column(6,
                  actionButton(ns("run_sling"), "Run Slingshot", 
                             class = "btn-primary", width = "100%")
                )
              )
            )
          )
        ),
        
        # Results
        fluidRow(
          column(9,
            withSpinner(plotlyOutput(ns("sling_plot"), height = "600px"))
          ),
          column(3,
            h5("Slingshot Results"),
            verbatimTextOutput(ns("sling_summary")),
            br(),
            h5("Download"),
            downloadButton(ns("dl_sling_pt"), "Pseudotime CSV", 
                          class = "btn-success", width = "100%"),
            br(), br(),
            downloadButton(ns("dl_sling_curves"), "Curves Data", 
                          class = "btn-info", width = "100%")
          )
        )
      ),

      # =========================== MONOCLE 3 ===========================
      tabPanel("Monocle 3",            icon = icon("leaf"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("Monocle 3 Configuration"),
              fluidRow(
                column(3,
                  numericInput(ns("cds_dim"), "Preprocessing dims", 
                             value = 50, min = 10, max = 100)
                ),
                column(3,
                  selectInput(ns("cds_reduce_method"), "Reduction method",
                            choices = c("UMAP", "tSNE", "PCA"), selected = "UMAP")
                ),
                column(3,
                  selectInput(ns("cds_norm_method"), "Normalization",
                            choices = c("log", "size_factors"), selected = "log")
                ),
                column(3,
                  checkboxInput(ns("cds_use_gene_short_name"), "Use gene short names", TRUE)
                )
              ),
              fluidRow(
                column(3,
                  checkboxInput(ns("cds_cluster_cells"), "Cluster cells", TRUE)
                ),
                column(3,
                  numericInput(ns("cds_resolution"), "Clustering resolution",
                             value = 0.01, min = 0.001, max = 0.5, step = 0.005),
                  helpText("Default 0.01. Lower values = fewer partitions. Range 0.01-0.05 recommended.", style="font-size:10px;color:#888;")
                ),
                column(3,
                  checkboxInput(ns("cds_close_loop"), "Close loop", FALSE),
                  checkboxInput(ns("cds_use_partition"), "Use partitions", TRUE)
                ),
                column(3,
                  actionButton(ns("run_cds"), "Run Monocle 3",
                             class = "btn-primary", width = "100%")
                )
              )
            )
          )
        ),
        
        # Results  
        fluidRow(
          column(9,
            withSpinner(plotlyOutput(ns("cds_plot"), height = "600px"))
          ),
          column(3,
            h5("Monocle 3 Results"),
            verbatimTextOutput(ns("cds_summary")),
            br(),
            h5("Download"),
            downloadButton(ns("dl_cds_pt"), "Pseudotime CSV", 
                          class = "btn-success", width = "100%"),
            br(), br(),
            downloadButton(ns("dl_cds_graph"), "Principal Graph", 
                          class = "btn-info", width = "100%")
          )
        )
      ),

      # ======================== DIFFUSION MAPS ==========================
      tabPanel("Diffusion Pseudotime", icon = icon("water"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("Diffusion Map Configuration"),
              fluidRow(
                column(3,
                  numericInput(ns("dpt_k"), "k-nearest neighbors", 
                             value = 30, min = 5, max = 200)
                ),
                column(3,
                  numericInput(ns("dpt_ndc"), "Diffusion components", 
                             value = 20, min = 5, max = 50)
                ),
                column(3,
                  numericInput(ns("dpt_topgenes"), "Max features for DiffusionMap",
                             value = 2000, min = 100, max = 10000)
                ),
                column(3,
                  selectInput(ns("dpt_distance"), "Distance metric",
                            choices = c("euclidean", "cosine", "correlation"),
                            selected = "euclidean")
                )
              ),
              fluidRow(
                column(3,
                  numericInput(ns("dpt_sigma"), "Sigma (bandwidth)", 
                             value = NULL, min = 0.1, max = 10)
                ),
                column(3,
                  checkboxInput(ns("dpt_density_norm"), "Density normalization", TRUE)
                ),
                column(3,
                  checkboxInput(ns("dpt_auto_sigma"), "Auto-estimate sigma", TRUE)
                ),
                column(3,
                  actionButton(ns("run_dpt"), "Run DPT", 
                             class = "btn-primary", width = "100%")
                )
              )
            )
          )
        ),
        
        # Results
        fluidRow(
          column(9,
            withSpinner(plotlyOutput(ns("dpt_plot"), height = "600px"))
          ),
          column(3,
            h5("DPT Results"),
            verbatimTextOutput(ns("dpt_summary")),
            br(),
            h5("Download"),
            downloadButton(ns("dl_dpt_pt"), "Pseudotime CSV", 
                          class = "btn-success", width = "100%"),
            br(), br(),
            downloadButton(ns("dl_dpt_dm"), "Diffusion Map", 
                          class = "btn-info", width = "100%")
          )
        )
      ),

      # ======================== COMPARISON ===========================
      tabPanel("Method Comparison",    icon = icon("columns"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("Compare Trajectory Methods"),
              fluidRow(
                column(4,
                  checkboxGroupInput(ns("methods_to_compare"), "Methods to compare",
                    choices = c("Slingshot" = "slingshot",
                              "Monocle 3" = "monocle3", 
                              "Diffusion Pseudotime" = "dpt"),
                    selected = c("slingshot", "monocle3"))
                ),
                column(4,
                  selectInput(ns("comparison_metric"), "Comparison metric",
                    choices = c("Correlation" = "correlation",
                              "Rank correlation" = "rank_correlation",
                              "Concordance" = "concordance"),
                    selected = "correlation")
                ),
                column(4,
                  div(style = "margin-top: 25px;",
                    actionButton(ns("run_comparison"), "Compare Methods", 
                               class = "btn-warning", width = "100%")
                  )
                )
              )
            )
          )
        ),
        
        # Comparison results
        fluidRow(
          column(6,
            h5("Correlation Matrix"),
            withSpinner(plotlyOutput(ns("comparison_heatmap"), height = "400px"))
          ),
          column(6,
            h5("Pseudotime Distributions"),
            withSpinner(plotlyOutput(ns("comparison_violin"), height = "400px"))
          )
        ),
        br(),
        fluidRow(
          column(10,
            h5("Side-by-side Comparison"),
            withSpinner(plotlyOutput(ns("comparison_sidebyside"), height = "500px"))
          ),
          column(2,
            h5("Download"),
            downloadButton(ns("dl_comparison_csv"), "Comparison CSV",
                          class = "btn-success", width = "100%")
          )
        )
      ),

      # ======================= GENE DYNAMICS ========================
      tabPanel("Gene Dynamics",        icon = icon("chart-area"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("Gene Expression Dynamics"),
              fluidRow(
                column(4,
                  textAreaInput(ns("dynamics_genes"), "Genes to analyze",
                    placeholder = "CD34, KIT, GATA1, HBB",
                    rows = 3)
                ),
                column(4,
                  selectInput(ns("dynamics_method"), "Pseudotime method",
                    choices = c("Slingshot" = "slingshot",
                              "Monocle 3" = "monocle3",
                              "DPT" = "dpt"),
                    selected = "slingshot")
                ),
                column(4,
                  div(style = "margin-top: 25px;",
                    actionButton(ns("analyze_dynamics"), "Analyze Dynamics", 
                               class = "btn-success", width = "100%")
                  )
                )
              ),
              fluidRow(
                column(3,
                  numericInput(ns("n_knots"), "Smoothing knots",
                             value = 6, min = 3, max = 15)
                ),
                column(3,
                  checkboxInput(ns("show_confidence"), "Show confidence bands", TRUE)
                ),
                column(3,
                  selectInput(ns("smooth_method"), "Smoothing method",
                    choices = c("loess", "gam", "lm"), selected = "loess")
                ),
                column(3,
                  numericInput(ns("smooth_span"), "Loess span",
                    value = 0.3, min = 0.1, max = 1.0, step = 0.05),
                  checkboxInput(ns("normalize_dynamics"), "Normalize expression", TRUE)
                )
              ),
              fluidRow(
                column(3,
                  checkboxInput(ns("split_by_lineage"), "Split by lineage", FALSE)
                ),
                column(3,
                  selectInput(ns("lineage_selection"), "Lineages to show",
                    choices = NULL, multiple = TRUE)
                ),
                column(3,
                  downloadButton(ns("dl_dynamics_csv"), "Download Dynamics CSV",
                                class = "btn-success", width = "100%")
                )
              )
            )
          )
        ),

        # Gene dynamics results
        fluidRow(
          column(12,
            withSpinner(plotlyOutput(ns("dynamics_plot"), height = "600px"))
          )
        )
      ),

      # ===================== BRANCH ANALYSIS ========================
      tabPanel("Branch Analysis",      icon = icon("code-branch"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("Lineage Branch Analysis"),
              fluidRow(
                column(4,
                  selectInput(ns("branch_method"), "Method for branching",
                    choices = c("Slingshot" = "slingshot", "Monocle 3" = "monocle3"),
                    selected = "slingshot")
                ),
                column(4,
                  numericInput(ns("branch_threshold"), "Branch point threshold",
                             value = 0.5, min = 0.1, max = 1, step = 0.1)
                ),
                column(4,
                  actionButton(ns("analyze_branches"), "Analyze Branches", 
                             class = "btn-info", width = "100%")
                )
              )
            )
          )
        ),
        
        # Branch analysis results
        fluidRow(
          column(6,
            h5("Branch Structure"),
            withSpinner(plotlyOutput(ns("branch_structure"), height = "400px"))
          ),
          column(6,
            h5("Branch Statistics"),
            withSpinner(DT::dataTableOutput(ns("branch_stats"))),
            br(),
            downloadButton(ns("dl_branch_csv"), "Download Branch CSV",
                          class = "btn-success", width = "100%")
          )
        )
      ),

      # ===================== RNA VELOCITY ========================
      tabPanel("RNA Velocity",         icon = icon("angle-double-right"),
        br(),
        fluidRow(
          column(12,
            wellPanel(
              h5("RNA Velocity Analysis"),
              p("Estimate RNA velocity from spliced/unspliced counts or load pre-computed velocity data."),
              helpText("For RNA velocity analysis, pre-compute velocities using scVelo or velociraptor and import the results. Required: spliced/unspliced count matrices in the @misc slot or a .loom file with velocity layers."),
              fluidRow(
                column(4,
                  selectInput(ns("velocity_source"), "Velocity data source",
                    choices = c("Upload .h5ad file" = "h5ad",
                              "Seurat object misc slot" = "misc",
                              "Compute from spliced/unspliced" = "compute"),
                    selected = "compute")
                ),
                column(4,
                  conditionalPanel(
                    condition = sprintf("input['%s'] == 'h5ad'", ns("velocity_source")),
                    fileInput(ns("velocity_h5ad"), "Upload .h5ad with velocity",
                             accept = c(".h5ad"))
                  )
                ),
                column(4,
                  selectInput(ns("velocity_method"), "Velocity calculation method",
                    choices = c("Steady-state ratio" = "steady_state",
                              "Dynamical model (reticulate)" = "dynamical"),
                    selected = "steady_state")
                )
              ),
              fluidRow(
                column(3,
                  selectInput(ns("velocity_embedding"), "Embedding for velocity plot",
                    choices = c("umap", "tsne", "pca"), selected = "umap")
                ),
                column(3,
                  numericInput(ns("velocity_arrow_scale"), "Arrow scale",
                             value = 1.0, min = 0.1, max = 5, step = 0.1)
                ),
                column(3,
                  numericInput(ns("velocity_grid_n"), "Grid resolution",
                             value = 30, min = 10, max = 80, step = 5)
                ),
                column(3,
                  div(style = "margin-top: 25px;",
                    actionButton(ns("run_velocity"), "Compute Velocity",
                               class = "btn-primary", width = "100%")
                  )
                )
              )
            )
          )
        ),

        # Velocity results
        fluidRow(
          column(8,
            h5("Velocity Stream Plot"),
            withSpinner(plotlyOutput(ns("velocity_stream_plot"), height = "600px"))
          ),
          column(4,
            h5("Velocity Summary"),
            verbatimTextOutput(ns("velocity_summary")),
            br(),
            h5("Velocity Confidence"),
            withSpinner(plotlyOutput(ns("velocity_confidence_plot"), height = "250px")),
            br(),
            h5("Velocity Magnitude"),
            withSpinner(plotlyOutput(ns("velocity_magnitude_plot"), height = "250px")),
            br(),
            downloadButton(ns("dl_velocity_csv"), "Download Velocity CSV", class = "btn-success", width = "100%")
          )
        )
      )
    )
  )
}

# =============================================================================
# TRAJECTORY SERVER
# =============================================================================

trajectoryServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Reactive values for storing results
    sling_res <- reactiveVal(NULL)
    cds_res <- reactiveVal(NULL)
    dpt_res <- reactiveVal(NULL)
    comparison_res <- reactiveVal(NULL)
    dynamics_res <- reactiveVal(NULL)
    branch_res <- reactiveVal(NULL)
    velocity_res <- reactiveVal(NULL)

    # ====================== UI UPDATES ======================
    update_choices <- function() {
      req(values$active_dataset)
      tryCatch({
        md <- values$active_dataset@meta.data

        # Cluster columns: character/factor columns + anything matching cluster patterns
        cand <- colnames(md)[sapply(md, function(x) is.factor(x) || is.character(x))]
        # Also include integer/numeric columns that look like cluster IDs
        cluster_pat <- grep("cluster|leiden|louvain|slm|seurat_clusters", colnames(md),
                            value = TRUE, ignore.case = TRUE)
        cand <- sort(unique(union(cand, cluster_pat)))

        if (length(cand) == 0) cand <- colnames(md)[1]

        # Pick best default
        default_col <- if ("seurat_clusters" %in% cand) "seurat_clusters"
                       else if ("cell_type_consensus" %in% cand) "cell_type_consensus"
                       else cand[1]

        updateSelectInput(session, "cluster_col", choices = cand, selected = default_col)

        # Immediately populate root_cluster from the default column
        if (default_col %in% colnames(md)) {
          levs <- sort(unique(as.character(md[[default_col]])))
          updateSelectInput(session, "root_cluster", choices = levs,
                            selected = if (length(levs)) levs[1] else NULL)
        }

        # Available reductions
        reds <- names(values$active_dataset@reductions)
        if (length(reds) > 0) {
          best <- .get_best_2d(values$active_dataset)
          if (is.null(best)) best <- reds[1]
          updateSelectInput(session, "reduction_for_traj", choices = reds, selected = best)
        }
      }, error = function(e) {
        message("Trajectory update_choices error: ", e$message)
        showNotification(paste("Trajectory setup error:", e$message), type = "warning", duration = 5)
      })
    }

    observe({
      req(values$active_dataset)
      update_choices()
    })

    # Update root_cluster choices when cluster column changes
    observeEvent(input$cluster_col, {
      req(values$active_dataset, input$cluster_col)
      md <- values$active_dataset@meta.data
      if (input$cluster_col %in% colnames(md)) {
        levs <- sort(unique(as.character(md[[input$cluster_col]])))
        updateSelectInput(session, "root_cluster", choices = levs,
                          selected = if (length(levs)) levs[1] else NULL)
      }
    })

    # ====================== SLINGSHOT ======================
    observeEvent(input$run_sling, {
      req(values$active_dataset, input$cluster_col, input$root_cluster)

      if (tolower(input$reduction_for_traj) == "umap") {
        showNotification(
          "WARNING: UMAP distorts distances non-linearly. Trajectory results on UMAP are visualization convenience only and should not be interpreted quantitatively. Use PCA for publication-grade analysis.",
          type = "warning", duration = 12
        )
      }

      if (input$parallel_traj) {
        tryCatch(
          plan(multisession, workers = min(4, future::availableCores() - 1)),
          error = function(e) plan(sequential)
        )
      }

      withProgress(message = "Running Slingshot trajectory analysis...", value = 0, {
        on.exit({ plan(sequential); gc() }, add = TRUE)

        result <- safely_run_slingshot(values$active_dataset, input)

        if (result$success) {
          sling_res(result)

          # Save to metadata if requested
          if (input$write_to_meta) {
            obj <- values$active_dataset
            obj$pseudotime_slingshot <- result$pseudotime
            obj$slingshot_lineage <- result$lineage_assignment
            values$active_dataset <- obj
          }

          showNotification("Slingshot analysis completed", type = "message")
        } else {
          showNotification(paste("Slingshot failed:", result$error), type = "error")
        }
      })
    })

    # ====================== MONOCLE 3 ======================
    observeEvent(input$run_cds, {
      req(values$active_dataset, input$cluster_col, input$root_cluster)

      if (input$parallel_traj) {
        tryCatch(
          plan(multisession, workers = min(4, future::availableCores() - 1)),
          error = function(e) plan(sequential)
        )
      }

      withProgress(message = "Running Monocle 3 trajectory analysis...", value = 0, {
        on.exit({ plan(sequential); gc() }, add = TRUE)

        result <- safely_run_monocle3(values$active_dataset, input)

        if (result$success) {
          cds_res(result)

          # Save to metadata if requested
          if (input$write_to_meta) {
            obj <- values$active_dataset
            obj$pseudotime_monocle3 <- result$pseudotime
            if (!is.null(result$pseudotime_raw)) {
              obj$pseudotime_raw <- result$pseudotime_raw
            }
            obj$monocle3_partition <- result$partition
            values$active_dataset <- obj
          }

          showNotification("Monocle 3 analysis completed", type = "message")
        } else {
          showNotification(paste("Monocle 3 failed:", result$error), type = "error")
        }
      })
    })

    # ====================== DIFFUSION PSEUDOTIME ======================
    observeEvent(input$run_dpt, {
      req(values$active_dataset, input$cluster_col, input$root_cluster)

      withProgress(message = "Computing Diffusion Pseudotime...", value = 0, {
        
        result <- safely_run_dpt(values$active_dataset, input)
        
        if (result$success) {
          dpt_res(result)
          
          # Save to metadata if requested
          if (input$write_to_meta) {
            obj <- values$active_dataset
            obj$pseudotime_dpt <- result$pseudotime
            obj$diffusion_component_1 <- result$dc1
            obj$diffusion_component_2 <- result$dc2
            values$active_dataset <- obj
          }
          
          showNotification("Diffusion pseudotime completed", type = "message")
        } else {
          showNotification(paste("DPT failed:", result$error), type = "error")
        }
      })
    })

    # ====================== METHOD COMPARISON ======================
    observeEvent(input$run_comparison, {
      req(length(input$methods_to_compare) >= 2)
      
      available_results <- list()
      if ("slingshot" %in% input$methods_to_compare && !is.null(sling_res())) {
        available_results$slingshot <- sling_res()$pseudotime
      }
      if ("monocle3" %in% input$methods_to_compare && !is.null(cds_res())) {
        available_results$monocle3 <- cds_res()$pseudotime
      }
      if ("dpt" %in% input$methods_to_compare && !is.null(dpt_res())) {
        available_results$dpt <- dpt_res()$pseudotime
      }
      
      if (length(available_results) < 2) {
        showNotification("Need at least 2 completed trajectory analyses to compare", type = "error")
        return()
      }
      
      withProgress(message = "Comparing trajectory methods...", {
        result <- compare_trajectory_methods(available_results, input$comparison_metric)
        comparison_res(result)
        showNotification("Method comparison completed", type = "message")
      })
    })

    # ====================== LINEAGE SELECTION UPDATE ======================
    # Update lineage choices when slingshot results become available
    observeEvent(sling_res(), {
      res <- sling_res()
      if (!is.null(res) && res$success && res$n_lineages > 1) {
        lineage_choices <- paste0("Lineage ", seq_len(res$n_lineages))
        updateSelectInput(session, "lineage_selection",
                          choices = lineage_choices,
                          selected = lineage_choices)
      }
    })

    # ====================== GENE DYNAMICS ======================
    observeEvent(input$analyze_dynamics, {
      req(values$active_dataset, nzchar(input$dynamics_genes))

      # Get selected pseudotime
      selected_method <- input$dynamics_method
      pseudotime_data <- NULL
      lineage_data    <- NULL

      if (selected_method == "slingshot" && !is.null(sling_res())) {
        pseudotime_data <- sling_res()$pseudotime
        # Provide per-lineage pseudotime if requested
        if (isTRUE(input$split_by_lineage) && !is.null(sling_res()$sce)) {
          lineage_data <- list(
            assignment = sling_res()$lineage_assignment,
            n_lineages = sling_res()$n_lineages,
            sce = sling_res()$sce,
            selected = input$lineage_selection
          )
        }
      } else if (selected_method == "monocle3" && !is.null(cds_res())) {
        pseudotime_data <- cds_res()$pseudotime
      } else if (selected_method == "dpt" && !is.null(dpt_res())) {
        pseudotime_data <- dpt_res()$pseudotime
      }

      if (is.null(pseudotime_data)) {
        showNotification(paste("No", selected_method, "results available"), type = "error")
        return()
      }

      withProgress(message = "Analyzing gene dynamics...", {
        result <- analyze_gene_dynamics(values$active_dataset, pseudotime_data, input,
                                        lineage_data = lineage_data)
        dynamics_res(result)
        showNotification("Gene dynamics analysis completed", type = "message")
      })
    })

    # ====================== BRANCH ANALYSIS ======================
    observeEvent(input$analyze_branches, {
      req(values$active_dataset)
      
      # Get branching data from selected method
      branch_data <- NULL
      if (input$branch_method == "slingshot" && !is.null(sling_res())) {
        branch_data <- sling_res()
      } else if (input$branch_method == "monocle3" && !is.null(cds_res())) {
        branch_data <- cds_res()
      }
      
      if (is.null(branch_data)) {
        showNotification(paste("No", input$branch_method, "results available"), type = "error")
        return()
      }
      
      withProgress(message = "Analyzing lineage branches...", {
        result <- analyze_lineage_branches(branch_data, input)
        branch_res(result)
        showNotification("Branch analysis completed", type = "message")
      })
    })

    # ====================== AI: ROOT CELL SUGGESTION ======================
    observeEvent(input$ai_root_suggest, {
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
      md  <- obj@meta.data

      # Get annotation columns
      annot_cols <- grep("cell_type|annotation|label", colnames(md),
                         value = TRUE, ignore.case = TRUE)

      # Get cluster info
      cluster_cols <- grep("cluster|leiden|louvain|seurat_clusters", colnames(md),
                           value = TRUE, ignore.case = TRUE)
      all_label_cols <- unique(c(annot_cols, cluster_cols))

      # Collect available cell type / cluster labels
      label_info <- if (length(all_label_cols) > 0) {
        paste(vapply(all_label_cols, function(col) {
          levs <- sort(unique(as.character(md[[col]])))
          paste0(col, ": [", paste(levs, collapse = ", "), "]")
        }, character(1)), collapse = "\n")
      } else {
        "No annotation or cluster columns found"
      }

      tissue  <- if (!is.null(values$tissue_global) && nzchar(values$tissue_global)) values$tissue_global else "unknown tissue"
      species <- if (!is.null(values$species_global) && nzchar(values$species_global)) values$species_global else "unknown species"

      prompt <- paste0(
        "I'm performing trajectory analysis on ", tissue, " data from ", species, ". ",
        "The following cell types/clusters are present:\n", label_info, "\n\n",
        "Based on developmental biology and the tissue context, which cluster or cell type ",
        "should be the root/starting point for trajectory inference? ",
        "Explain your reasoning based on known differentiation hierarchies."
      )

      output$ai_thread_root <- renderUI({ ai_thread_loading("AI analyzing root cell candidates...") })
      result <- call_ai("root_suggestion", prompt, values)

      if (result$ok) {
        output$ai_thread_root <- renderUI({
          ai_thread_message(
            content_html = .mdtohtml(result$text),
            title = "AI Root Cell Suggestion",
            model = get_model_for_task("root_suggestion", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out
          )
        })
      } else {
        output$ai_thread_root <- renderUI({ ai_thread_error(result$error) })
      }
    })

    # ====================== AI: TRAJECTORY NARRATIVE ======================
    observeEvent(input$ai_trajectory_narrative, {
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
      md  <- obj@meta.data

      # Check for pseudotime columns
      pt_cols <- grep("pseudotime|dpt|monocle", colnames(md),
                      value = TRUE, ignore.case = TRUE)

      if (length(pt_cols) == 0) {
        showNotification("No pseudotime data found in metadata. Please run a trajectory analysis first.",
                         type = "warning", duration = 8)
        return()
      }

      # Get annotation/cluster labels
      annot_cols <- grep("cell_type|annotation|label", colnames(md),
                         value = TRUE, ignore.case = TRUE)
      cluster_cols <- grep("cluster|leiden|louvain|seurat_clusters", colnames(md),
                           value = TRUE, ignore.case = TRUE)
      label_col <- if (length(annot_cols) > 0) annot_cols[1] else if (length(cluster_cols) > 0) cluster_cols[1] else NULL

      # Compute cell type ordering by mean pseudotime
      ordering_info <- ""
      if (!is.null(label_col)) {
        labels <- as.character(md[[label_col]])
        for (ptc in pt_cols) {
          pt_vals <- as.numeric(md[[ptc]])
          if (all(is.na(pt_vals))) next
          mean_pt <- tapply(pt_vals, labels, mean, na.rm = TRUE)
          mean_pt <- sort(mean_pt)
          ordering_info <- paste0(
            ordering_info,
            "Pseudotime column '", ptc, "' - cell type ordering by mean pseudotime:\n",
            paste(paste0("  ", seq_along(mean_pt), ". ", names(mean_pt),
                         " (mean PT: ", round(mean_pt, 3), ")"), collapse = "\n"),
            "\n\n"
          )
        }
      }

      # Check for branch point info
      branch_info <- ""
      n_sling_lineages <- if (!is.null(sling_res()) && sling_res()$success) sling_res()$n_lineages else 0
      n_cds_partitions <- if (!is.null(cds_res()) && cds_res()$success) cds_res()$n_partitions else 0
      if (n_sling_lineages > 1) {
        branch_info <- paste0(branch_info, "Slingshot detected ", n_sling_lineages, " lineages (branch points present). ")
      }
      if (n_cds_partitions > 1) {
        branch_info <- paste0(branch_info, "Monocle 3 detected ", n_cds_partitions, " partitions. ")
      }

      tissue  <- if (!is.null(values$tissue_global) && nzchar(values$tissue_global)) values$tissue_global else "unknown tissue"
      species <- if (!is.null(values$species_global) && nzchar(values$species_global)) values$species_global else "unknown species"

      prompt <- paste0(
        "I have performed trajectory/pseudotime analysis on ", tissue, " data from ", species, ". ",
        "Please interpret the following trajectory results:\n\n",
        ordering_info,
        if (nzchar(branch_info)) paste0("Branch information: ", branch_info, "\n\n") else "",
        "Please provide:\n",
        "1. A biological narrative of the cell type transitions along the trajectory\n",
        "2. Whether the ordering matches expected differentiation hierarchies for this tissue\n",
        "3. Any surprising or unexpected orderings and possible explanations\n",
        "4. The biological significance of the trajectory and potential implications\n",
        "5. Suggested validation experiments or follow-up analyses"
      )

      output$ai_thread_trajectory <- renderUI({ ai_thread_loading("AI interpreting trajectory...") })
      result <- call_ai("trajectory_narrative", prompt, values)

      if (result$ok) {
        output$ai_thread_trajectory <- renderUI({
          ai_thread_message(
            content_html = .mdtohtml(result$text),
            title = "AI Trajectory Interpretation",
            model = get_model_for_task("trajectory_narrative", values),
            tokens_in = result$tokens_in,
            tokens_out = result$tokens_out
          )
        })
      } else {
        output$ai_thread_trajectory <- renderUI({ ai_thread_error(result$error) })
      }
    })

    # ====================== AI HELPER: GENERIC Q&A ======================
    ai_helper_server("trajectory", values, "Trajectory Analysis", context_fn = function() {
      obj <- values$active_dataset
      if (is.null(obj) || !inherits(obj, "Seurat")) return("")
      parts <- character(0)
      md <- obj@meta.data

      # Pseudotime availability
      pt_cols <- grep("pseudotime|dpt|monocle", colnames(md),
                      value = TRUE, ignore.case = TRUE)
      if (length(pt_cols) > 0) {
        parts <- c(parts, paste0("Pseudotime columns available: ", paste(pt_cols, collapse = ", ")))
      } else {
        parts <- c(parts, "No pseudotime computed yet.")
      }

      # Methods with results
      methods_run <- character(0)
      if (!is.null(sling_res()) && sling_res()$success) {
        methods_run <- c(methods_run, paste0("Slingshot (", sling_res()$n_lineages, " lineages)"))
      }
      if (!is.null(cds_res()) && cds_res()$success) {
        methods_run <- c(methods_run, paste0("Monocle 3 (", cds_res()$n_partitions, " partitions)"))
      }
      if (!is.null(dpt_res()) && dpt_res()$success) {
        methods_run <- c(methods_run, paste0("DPT (", dpt_res()$n_components, " components)"))
      }
      if (length(methods_run) > 0) {
        parts <- c(parts, paste0("Trajectory methods completed: ", paste(methods_run, collapse = ", ")))
      }

      # Cell type labels
      annot_cols <- grep("cell_type|annotation|label", colnames(md),
                         value = TRUE, ignore.case = TRUE)
      cluster_cols <- grep("cluster|leiden|louvain|seurat_clusters", colnames(md),
                           value = TRUE, ignore.case = TRUE)
      label_cols <- unique(c(annot_cols, cluster_cols))
      if (length(label_cols) > 0) {
        label_col <- label_cols[1]
        levs <- sort(unique(as.character(md[[label_col]])))
        parts <- c(parts, paste0("Cell types/clusters (", label_col, "): ",
                                  paste(levs, collapse = ", ")))
      }

      paste(parts, collapse = "\n")
    })

    # ====================== OUTPUT RENDERERS ======================

    .traj_error_plot <- function(e) {
      plotly_empty() %>% layout(
        title = "Plot error",
        annotations = list(list(text = paste("Error:", e$message),
          showarrow = FALSE, x = 0.5, y = 0.5, xref = "paper", yref = "paper"))
      )
    }

    # Slingshot plot
    output$sling_plot <- renderPlotly({
      req(sling_res(), values$active_dataset)
      p <- tryCatch(
        create_slingshot_plot(values$active_dataset, sling_res(), input),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Slingshot trajectory with pseudotime",
                       label = "Slingshot Trajectory", source_module = "Trajectory")
      }
      p
    })

    output$sling_summary <- renderText({
      req(sling_res())
      result <- sling_res()
      paste(
        "Method: Slingshot",
        paste("Lineages:", result$n_lineages),
        paste("Cells with pseudotime:", sum(!is.na(result$pseudotime))),
        paste("Mean pseudotime:", round(mean(result$pseudotime, na.rm = TRUE), 3)),
        sep = "\n"
      )
    })

    # Monocle 3 plot
    output$cds_plot <- renderPlotly({
      req(cds_res(), values$active_dataset)
      p <- tryCatch(
        create_monocle3_plot(values$active_dataset, cds_res(), input),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Monocle3 trajectory with pseudotime",
                       label = "Monocle3 Trajectory", source_module = "Trajectory")
      }
      p
    })

    output$cds_summary <- renderText({
      req(cds_res())
      result <- cds_res()
      paste(
        "Method: Monocle 3",
        paste("Partitions:", result$n_partitions),
        paste("Cells with pseudotime:", sum(!is.na(result$pseudotime))),
        paste("Mean pseudotime:", round(mean(result$pseudotime, na.rm = TRUE), 3)),
        sep = "\n"
      )
    })

    # DPT plot
    output$dpt_plot <- renderPlotly({
      req(dpt_res(), values$active_dataset)
      p <- tryCatch(
        create_dpt_plot(values$active_dataset, dpt_res(), input),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Diffusion pseudotime trajectory",
                       label = "DPT Trajectory", source_module = "Trajectory")
      }
      p
    })

    output$dpt_summary <- renderText({
      req(dpt_res())
      result <- dpt_res()
      paste(
        "Method: Diffusion Pseudotime",
        paste("Diffusion components:", result$n_components),
        paste("Cells analyzed:", length(result$pseudotime)),
        paste("Mean pseudotime:", round(mean(result$pseudotime, na.rm = TRUE), 3)),
        sep = "\n"
      )
    })

    # Comparison plots
    output$comparison_heatmap <- renderPlotly({
      req(comparison_res())
      p <- tryCatch(
        create_comparison_heatmap(comparison_res()),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Trajectory method comparison heatmap",
                       label = "Trajectory Comparison Heatmap", source_module = "Trajectory")
      }
      p
    })

    output$comparison_violin <- renderPlotly({
      req(comparison_res())
      p <- tryCatch(
        create_comparison_violin(comparison_res()),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Trajectory method comparison violin",
                       label = "Trajectory Comparison Violin", source_module = "Trajectory")
      }
      p
    })

    output$comparison_sidebyside <- renderPlotly({
      req(comparison_res(), values$active_dataset)
      p <- tryCatch(
        create_sidebyside_comparison(values$active_dataset, comparison_res(), input),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Side-by-side trajectory comparison",
                       label = "Trajectory Side-by-Side", source_module = "Trajectory")
      }
      p
    })

    # Gene dynamics plot
    output$dynamics_plot <- renderPlotly({
      req(dynamics_res())
      p <- tryCatch(
        create_gene_dynamics_plot(dynamics_res(), input),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Gene expression dynamics along trajectory",
                       label = "Gene Dynamics", source_module = "Trajectory")
      }
      p
    })

    # Branch analysis outputs
    output$branch_structure <- renderPlotly({
      req(branch_res())
      p <- tryCatch(
        create_branch_structure_plot(branch_res()),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Trajectory branch structure",
                       label = "Branch Structure", source_module = "Trajectory")
      }
      p
    })

    output$branch_stats <- DT::renderDataTable({
      req(branch_res())
      tbl <- branch_res()$branch_statistics
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Branch Statistics", source_module = "Trajectory")
      }
      tbl
    }, options = list(dom = 't', pageLength = 10))

    # ====================== DOWNLOAD HANDLERS ======================
    
    output$dl_sling_pt <- downloadHandler(
      filename = function() paste0("slingshot_pseudotime_", Sys.Date(), ".csv"),
      content = function(file) {
        req(sling_res())
        result <- sling_res()
        df <- data.frame(
          Cell = names(result$pseudotime),
          Pseudotime = result$pseudotime,
          Lineage = result$lineage_assignment
        )
        write.csv(df, file, row.names = FALSE)
      }
    )

    output$dl_sling_curves <- downloadHandler(
      filename = function() paste0("slingshot_curves_", Sys.Date(), ".csv"),
      content = function(file) {
        req(sling_res())
        if (!is.null(sling_res()$curves_data)) {
          write.csv(sling_res()$curves_data, file, row.names = FALSE)
        } else {
          showNotification("No curve data available for download.", type = "warning")
          write.csv(data.frame(Note = "No curve data available"), file, row.names = FALSE)
        }
      }
    )

    output$dl_cds_pt <- downloadHandler(
      filename = function() paste0("monocle3_pseudotime_", Sys.Date(), ".csv"),
      content = function(file) {
        req(cds_res())
        result <- cds_res()
        df <- data.frame(
          Cell = names(result$pseudotime),
          Pseudotime = result$pseudotime,
          Partition = result$partition
        )
        write.csv(df, file, row.names = FALSE)
      }
    )

    output$dl_dpt_pt <- downloadHandler(
      filename = function() paste0("dpt_pseudotime_", Sys.Date(), ".csv"),
      content = function(file) {
        req(dpt_res())
        result <- dpt_res()
        df <- data.frame(
          Cell = names(result$pseudotime),
          Pseudotime = result$pseudotime,
          DC1 = result$dc1,
          DC2 = result$dc2
        )
        write.csv(df, file, row.names = FALSE)
      }
    )

    output$dl_comparison_csv <- downloadHandler(
      filename = function() paste0("trajectory_comparison_", Sys.Date(), ".csv"),
      content = function(file) {
        req(comparison_res())
        result <- comparison_res()
        # Export the pseudotime matrix with cell identifiers
        df <- as.data.frame(result$pseudotime_matrix)
        df$Cell <- result$common_cells
        write.csv(df, file, row.names = FALSE)
      }
    )

    output$dl_dynamics_csv <- downloadHandler(
      filename = function() paste0("gene_dynamics_", Sys.Date(), ".csv"),
      content = function(file) {
        req(dynamics_res())
        write.csv(dynamics_res()$dynamics_data, file, row.names = FALSE)
      }
    )

    output$dl_branch_csv <- downloadHandler(
      filename = function() paste0("branch_statistics_", Sys.Date(), ".csv"),
      content = function(file) {
        req(branch_res())
        result <- branch_res()
        # Write branch statistics; append DE results if available
        write.csv(result$branch_statistics, file, row.names = FALSE)
        if (!is.null(result$branch_de)) {
          de_file <- sub("\\.csv$", "_de.csv", file)
          write.csv(result$branch_de, de_file, row.names = FALSE)
        }
      }
    )

    # ====================== RNA VELOCITY ======================
    observeEvent(input$run_velocity, {
      req(values$active_dataset)

      withProgress(message = "Computing RNA velocity...", value = 0, {
        result <- compute_rna_velocity(values$active_dataset, input)

        if (result$success) {
          velocity_res(result)
          showNotification("RNA velocity analysis completed", type = "message")
        } else {
          showNotification(paste("Velocity failed:", result$error), type = "error")
        }
      })
    })

    output$velocity_stream_plot <- renderPlotly({
      req(velocity_res(), values$active_dataset)
      p <- tryCatch(
        create_velocity_stream_plot(values$active_dataset, velocity_res(), input),
        error = .traj_error_plot
      )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "RNA velocity stream plot",
                       label = "RNA Velocity Stream", source_module = "Trajectory")
      }
      p
    })

    output$velocity_summary <- renderText({
      req(velocity_res())
      res <- velocity_res()
      paste(
        "Method: RNA Velocity",
        paste("Source:", res$source),
        paste("Calculation:", res$calc_method),
        paste("Cells analyzed:", length(res$magnitude)),
        paste("Mean magnitude:", round(mean(res$magnitude, na.rm = TRUE), 4)),
        paste("Mean confidence:", round(mean(res$confidence, na.rm = TRUE), 4)),
        sep = "\n"
      )
    })

    output$velocity_confidence_plot <- renderPlotly({
      req(velocity_res())
      p <- tryCatch({
        res <- velocity_res()
        plot_ly(x = res$confidence, type = "histogram",
                marker = list(color = "steelblue"),
                hovertemplate = "Confidence: %{x:.3f}<br>Count: %{y}<extra></extra>") %>%
          layout(
            title = "Velocity Confidence Distribution",
            xaxis = list(title = "Confidence Score"),
            yaxis = list(title = "Cells"),
            autosize = TRUE
          )
      }, error = .traj_error_plot)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Velocity confidence distribution",
                       label = "Velocity Confidence", source_module = "Trajectory")
      }
      p
    })

    output$velocity_magnitude_plot <- renderPlotly({
      req(velocity_res())
      p <- tryCatch({
        res <- velocity_res()
        plot_ly(x = res$magnitude, type = "histogram",
                marker = list(color = "firebrick"),
                hovertemplate = "Magnitude: %{x:.4f}<br>Count: %{y}<extra></extra>") %>%
          layout(
            title = "Velocity Magnitude Distribution",
            xaxis = list(title = "Magnitude"),
            yaxis = list(title = "Cells"),
            autosize = TRUE
          )
      }, error = .traj_error_plot)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Trajectory", figure_type = "Velocity magnitude distribution",
                       label = "Velocity Magnitude", source_module = "Trajectory")
      }
      p
    })

    # Velocity download handler
    output$dl_velocity_csv <- downloadHandler(
      filename = function() paste0("rna_velocity_", Sys.Date(), ".csv"),
      content = function(file) {
        req(velocity_res())
        res <- velocity_res()
        df <- data.frame(
          Cell = res$cells,
          Arrow_DX = res$arrow_dx,
          Arrow_DY = res$arrow_dy,
          Magnitude = res$magnitude,
          Confidence = res$confidence,
          Embed_X = res$embedding[, 1],
          Embed_Y = res$embedding[, 2]
        )
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}

# =============================================================================
# CORE ANALYSIS FUNCTIONS
# =============================================================================

safely_run_slingshot <- function(seurat_obj, params) {
  tryCatch({
    if (!requireNamespace("slingshot", quietly = TRUE)) {
      stop("slingshot package required")
    }
    
    red <- tolower(params$reduction_for_traj)
    if (!(red %in% names(seurat_obj@reductions))) {
      stop(sprintf("Reduction '%s' not found", red))
    }

    # Convert to SCE
    sce <- .as_sce_with_reduction(seurat_obj, red)
    
    # Set cluster labels
    colData(sce)$clusters <- as.factor(seurat_obj@meta.data[[params$cluster_col]])
    
    # Run slingshot (set seed for reproducibility)
    set.seed(42)
    extend_val <- if (isTRUE(params$sling_extend)) "y" else "n"
    sce <- slingshot(
      sce,
      clusterLabels = "clusters",
      reducedDim = toupper(red),
      start.clus = params$root_cluster,
      shrink = params$sling_shrink,
      extend = extend_val,
      thresh = params$sling_thresh
    )

    # Extract pseudotime
    pst <- slingPseudotime(sce)
    
    if (is.null(dim(pst))) {
      pseudotime <- pst
      lineage_assignment <- rep(1, length(pst))
    } else {
      if (params$sling_lineage_mode == "first") {
        pseudotime <- pst[, 1]
        lineage_assignment <- rep(1, nrow(pst))
      } else if (params$sling_lineage_mode == "average") {
        pseudotime <- rowMeans(pst, na.rm = TRUE)
        lineage_assignment <- apply(pst, 1, function(x) which.min(x))
      } else {
        # Principal curve - use first lineage but could be enhanced
        pseudotime <- pst[, 1]
        lineage_assignment <- apply(pst, 1, function(x) which.min(x))
      }
    }
    
    # Store raw pseudotime before rescaling
    seurat_obj[["pseudotime_raw"]] <- pseudotime
    pseudotime <- .scale01(pseudotime)
    names(pseudotime) <- colnames(sce)
    names(lineage_assignment) <- colnames(sce)

    # Extract curve data
    curves_data <- NULL
    if (params$sling_show_curves) {
      curves <- slingCurves(sce)
      if (length(curves) > 0) {
        curve_list <- list()
        for (i in seq_along(curves)) {
          curve_df <- as.data.frame(curves[[i]]$s)
          curve_df$lineage <- i
          curve_list[[i]] <- curve_df
        }
        curves_data <- do.call(rbind, curve_list)
      }
    }

    list(
      success = TRUE,
      pseudotime = pseudotime,
      lineage_assignment = lineage_assignment,
      n_lineages = if (is.null(dim(pst))) 1L else ncol(pst),
      sce = sce,
      curves_data = curves_data,
      method = "slingshot",
      reduction = red,
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      pseudotime = NULL,
      lineage_assignment = NULL,
      n_lineages = 0,
      sce = NULL,
      curves_data = NULL,
      method = "slingshot",
      reduction = NULL,
      error = e$message
    )
  })
}

safely_run_monocle3 <- function(seurat_obj, params) {
  tryCatch({
    if (!requireNamespace("monocle3", quietly = TRUE)) {
      stop("monocle3 package required")
    }

    # Convert to CDS
    cds <- .as_cds(seurat_obj)

    # Detect if data is already log-normalized (avoid double normalization)
    norm_method <- params$cds_norm_method
    sample_rows <- min(100, nrow(cds))
    sample_cols <- min(100, ncol(cds))
    sample_vals <- monocle3::exprs(cds)[seq_len(sample_rows), seq_len(sample_cols)]
    if (max(sample_vals, na.rm = TRUE) < 20) {
      # Values < 20 typical for log1p-normalized data; skip re-normalization
      norm_method <- "none"
      showNotification("Data appears pre-normalized. Skipping Monocle3 internal normalization.", type = "message")
    }

    # Use HVGs from Seurat if available
    hvgs <- VariableFeatures(seurat_obj)
    use_genes <- if (length(hvgs) > 0) hvgs else NULL

    # Guard num_dim against dataset dimensions
    num_dim <- min(params$cds_dim, ncol(cds) - 1, nrow(cds) - 1)

    # Preprocess
    cds <- monocle3::preprocess_cds(
      cds,
      num_dim = num_dim,
      method = "PCA",
      norm_method = norm_method,
      use_genes = use_genes,
      verbose = FALSE
    )

    # Reduce dimension
    cds <- monocle3::reduce_dimension(
      cds,
      reduction_method = params$cds_reduce_method,
      verbose = FALSE
    )

    if (params$cds_reduce_method == "UMAP") {
      showNotification("WARNING: UMAP distorts distances. Monocle3 trajectory topology is approximate.", type = "warning", duration = 10)
    }

    # Cluster cells (required for learn_graph partitions)
    set.seed(42)
    cds <- monocle3::cluster_cells(
      cds,
      reduction_method = params$cds_reduce_method,
      resolution = if (isTRUE(params$cds_cluster_cells)) params$cds_resolution else 1e-4,
      verbose = FALSE
    )

    # Learn graph
    use_partition <- if (!is.null(params$cds_use_partition)) isTRUE(params$cds_use_partition) else TRUE
    min_branch <- if (!is.null(params$min_branch_len)) params$min_branch_len else 10
    cds <- monocle3::learn_graph(
      cds,
      use_partition = use_partition,
      close_loop = params$cds_close_loop,
      minimal_branch_len = min_branch,
      verbose = FALSE
    )
    
    # Get root cells
    root_cells <- rownames(colData(cds))[
      seurat_obj@meta.data[[params$cluster_col]] == params$root_cluster
    ]
    
    if (length(root_cells) == 0) {
      stop("No cells found in root cluster")
    }
    
    # Order cells
    cds <- monocle3::order_cells(cds, root_cells = root_cells)
    
    # Extract pseudotime
    pseudotime <- monocle3::pseudotime(cds)
    pseudotime <- as.numeric(pseudotime)
    # Handle Inf pseudotime from unrooted partitions
    pseudotime[is.infinite(pseudotime)] <- NA_real_
    # Store raw pseudotime before rescaling (returned in result list for write-back)
    pseudotime_raw <- pseudotime
    pseudotime <- .scale01(pseudotime)
    names(pseudotime) <- colnames(cds)
    names(pseudotime_raw) <- colnames(cds)

    # Extract partition info
    if (params$cds_cluster_cells) {
      partition <- monocle3::partitions(cds)
    } else {
      partition <- rep(1, ncol(cds))
      names(partition) <- colnames(cds)
    }

    list(
      success = TRUE,
      pseudotime = pseudotime,
      pseudotime_raw = pseudotime_raw,
      partition = partition,
      n_partitions = length(unique(partition)),
      cds = cds,
      method = "monocle3",
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      pseudotime = NULL,
      partition = NULL,
      n_partitions = 0,
      cds = NULL,
      method = "monocle3",
      error = e$message
    )
  })
}

safely_run_dpt <- function(seurat_obj, params) {
  tryCatch({
    if (!requireNamespace("destiny", quietly = TRUE)) {
      stop("destiny package required")
    }

    # Get features
    if (params$dpt_topgenes > 0) {
      hvgs <- VariableFeatures(seurat_obj)
      if (length(hvgs) == 0) {
        seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = params$dpt_topgenes, verbose = FALSE)
        hvgs <- VariableFeatures(seurat_obj)
      }
      genes <- head(hvgs, params$dpt_topgenes)
    } else {
      genes <- rownames(seurat_obj)
    }

    # Sketch gate: downsample for diffusion map (very expensive O(n^2))
    dpt_sketched <- FALSE
    n_cells <- ncol(seurat_obj)
    if (n_cells > 30000) {
      message("Large dataset (", format(n_cells, big.mark = ","), " cells). Downsampling to 30,000 for diffusion map.")
      full_obj <- seurat_obj
      seurat_obj <- sketch_downsample(seurat_obj, max_cells = 30000L, method = "geometric")
      dpt_sketched <- TRUE
    }

    # Memory guard before dense matrix conversion
    n_elements <- length(genes) * ncol(seurat_obj)
    if (n_elements > 1.2e9) {
      stop("Expression matrix too large (", format(n_elements, scientific = TRUE),
           " elements). Reduce number of genes or cells.")
    }

    # Get expression matrix
    expr_matrix <- as.matrix(safe_get_assay_data(seurat_obj, "data")[genes, ])

    # Create diffusion map
    dm_params <- list(
      k = params$dpt_k,
      n_eigs = params$dpt_ndc,
      distance = params$dpt_distance,
      density_norm = params$dpt_density_norm
    )

    # Use PCA pre-reduction when gene count is high (avoids destiny warning)
    if (length(genes) > 500) {
      dm_params$n_pcs <- min(50L, length(genes) - 1L, ncol(seurat_obj) - 1L)
    }

    if (!is.null(params$dpt_sigma) && !is.na(params$dpt_sigma) && !params$dpt_auto_sigma) {
      dm_params$sigma <- params$dpt_sigma
    }

    dm <- suppressWarnings(do.call(destiny::DiffusionMap, c(list(data = t(expr_matrix)), dm_params)))

    # Find root cell
    root_cells <- rownames(seurat_obj@meta.data)[
      seurat_obj@meta.data[[params$cluster_col]] == params$root_cluster
    ]
    
    if (length(root_cells) == 0) {
      stop("No cells found in root cluster")
    }
    
    # Use medoid of root cluster in diffusion space
    dm_coords <- destiny::eigenvectors(dm)
    root_coords <- dm_coords[root_cells, 1:2, drop = FALSE]
    root_medoid <- root_cells[which.min(rowSums((root_coords - colMeans(root_coords))^2))]

    # Compute DPT
    dpt <- destiny::DPT(dm, tips = root_medoid)
    
    # Extract pseudotime (first branch)
    pseudotime <- destiny::dpt(dpt)[, 1]
    # Store raw pseudotime before rescaling
    seurat_obj[["pseudotime_raw"]] <- as.numeric(pseudotime)
    pseudotime <- .scale01(as.numeric(pseudotime))
    names(pseudotime) <- rownames(dm_coords)

    # Extract diffusion components
    dc1 <- dm_coords[, 1]
    dc2 <- dm_coords[, 2]

    # Project pseudotime to full dataset if sketched
    if (dpt_sketched) {
      full_pseudotime <- rep(NA_real_, ncol(full_obj))
      names(full_pseudotime) <- colnames(full_obj)
      # Copy computed values for sketched cells
      full_pseudotime[names(pseudotime)] <- pseudotime
      # Project to remaining cells via nearest-neighbor in PCA space
      unseen <- setdiff(colnames(full_obj), names(pseudotime))
      if (length(unseen) > 0 && "pca" %in% names(full_obj@reductions)) {
        pca_all <- Embeddings(full_obj, "pca")
        pca_seen <- pca_all[names(pseudotime), , drop = FALSE]
        pca_unseen <- pca_all[unseen, , drop = FALSE]
        # Process in chunks to limit memory
        chunk_size <- 5000L
        for (start in seq(1, length(unseen), by = chunk_size)) {
          end <- min(start + chunk_size - 1, length(unseen))
          chunk_cells <- unseen[start:end]
          dists <- as.matrix(dist(rbind(pca_unseen[chunk_cells, , drop = FALSE], pca_seen)))
          n_chunk <- length(chunk_cells)
          n_seen <- nrow(pca_seen)
          dist_block <- dists[1:n_chunk, (n_chunk + 1):(n_chunk + n_seen), drop = FALSE]
          nn_idx <- max.col(-dist_block, ties.method = "first")
          full_pseudotime[chunk_cells] <- pseudotime[nn_idx]
        }
      } else {
        # Fallback: interpolate with median pseudotime per cluster
        if (!is.null(params$cluster_col) && params$cluster_col %in% colnames(full_obj@meta.data)) {
          cl_seen <- full_obj@meta.data[names(pseudotime), params$cluster_col]
          cl_medians <- tapply(pseudotime, cl_seen, median, na.rm = TRUE)
          cl_unseen <- full_obj@meta.data[unseen, params$cluster_col]
          full_pseudotime[unseen] <- unname(cl_medians[as.character(cl_unseen)])
        }
      }
      pseudotime <- full_pseudotime
      # Update coords for full dataset
      dc1_full <- rep(NA_real_, ncol(full_obj))
      dc2_full <- rep(NA_real_, ncol(full_obj))
      names(dc1_full) <- colnames(full_obj)
      names(dc2_full) <- colnames(full_obj)
      dc1_full[names(dc1)] <- dc1
      dc2_full[names(dc2)] <- dc2
      dc1 <- dc1_full
      dc2 <- dc2_full
    }

    list(
      success = TRUE,
      pseudotime = pseudotime,
      dc1 = dc1,
      dc2 = dc2,
      n_components = params$dpt_ndc,
      dm = dm,
      dpt = dpt,
      genes_used = genes,
      method = "dpt",
      error = NULL
    )
    
  }, error = function(e) {
    list(
      success = FALSE,
      pseudotime = NULL,
      dc1 = NULL,
      dc2 = NULL,
      n_components = 0,
      dm = NULL,
      dpt = NULL,
      genes_used = NULL,
      method = "dpt",
      error = e$message
    )
  })
}

# =============================================================================
# COMPARISON AND ANALYSIS FUNCTIONS
# =============================================================================

compare_trajectory_methods <- function(results_list, metric = "correlation") {
  # Create comparison matrix
  methods <- names(results_list)
  n_methods <- length(methods)
  
  # Find common cells
  common_cells <- Reduce(intersect, lapply(results_list, names))
  
  if (length(common_cells) < 10) {
    stop("Not enough common cells between methods")
  }
  
  # Extract pseudotimes for common cells
  pt_matrix <- do.call(cbind, lapply(results_list, function(x) x[common_cells]))
  colnames(pt_matrix) <- methods
  
  # Compute correlation matrix
  if (metric == "correlation") {
    cor_matrix <- cor(pt_matrix, use = "complete.obs")
  } else if (metric == "rank_correlation") {
    cor_matrix <- cor(pt_matrix, use = "complete.obs", method = "spearman")
  } else {
    # Concordance - simplified implementation
    cor_matrix <- cor(pt_matrix, use = "complete.obs")
  }
  
  list(
    correlation_matrix = cor_matrix,
    pseudotime_matrix = pt_matrix,
    common_cells = common_cells,
    methods = methods
  )
}

analyze_gene_dynamics <- function(seurat_obj, pseudotime, params,
                                  lineage_data = NULL) {
  # Parse genes
  genes <- unique(trimws(unlist(strsplit(params$dynamics_genes, ","))))
  genes <- genes[genes %in% rownames(seurat_obj)]

  if (length(genes) == 0) {
    stop("No valid genes specified")
  }

  # Get expression data
  expr_data <- safe_get_assay_data(seurat_obj, "data")[genes, names(pseudotime), drop = FALSE]

  # Determine whether to split by lineage
  split_by_lineage <- !is.null(lineage_data) && isTRUE(params$split_by_lineage)

  # Create data frame for analysis
  dynamics_data <- list()

  if (split_by_lineage) {
    # Per-lineage pseudotime from slingshot SCE
    pst_mat <- NULL
    tryCatch({
      pst_mat <- slingPseudotime(lineage_data$sce)
    }, error = function(e) NULL)

    # Determine which lineages the user selected
    selected_ids <- seq_len(lineage_data$n_lineages)
    if (!is.null(lineage_data$selected) && length(lineage_data$selected) > 0) {
      # Parse "Lineage 1", "Lineage 2", ... to numeric ids
      parsed <- as.integer(gsub("\\D", "", lineage_data$selected))
      parsed <- parsed[!is.na(parsed) & parsed >= 1 & parsed <= lineage_data$n_lineages]
      if (length(parsed) > 0) selected_ids <- parsed
    }

    for (gene in genes) {
      gene_expr_all <- as.numeric(expr_data[gene, ])
      names(gene_expr_all) <- names(pseudotime)

      for (lid in selected_ids) {
        lineage_label <- paste0("Lineage_", lid)

        if (!is.null(pst_mat) && ncol(pst_mat) >= lid) {
          # Use lineage-specific pseudotime
          lin_pt <- pst_mat[, lid]
          names(lin_pt) <- rownames(pst_mat)
          valid <- !is.na(lin_pt)
          cells <- names(lin_pt)[valid]
          pt_vals <- .scale01(lin_pt[valid])
        } else {
          # Fall back: use assignment-based partitioning
          cells <- names(lineage_data$assignment)[lineage_data$assignment == lid]
          cells <- intersect(cells, names(pseudotime))
          pt_vals <- pseudotime[cells]
        }

        if (length(cells) < 3) next

        g_expr <- gene_expr_all[cells]
        if (params$normalize_dynamics) {
          g_expr <- scale(g_expr)[, 1]
        }

        gene_df <- data.frame(
          Gene       = gene,
          Pseudotime = pt_vals,
          Expression = g_expr,
          Cell       = cells,
          Lineage    = lineage_label,
          stringsAsFactors = FALSE
        )
        gene_df <- gene_df[!is.na(gene_df$Pseudotime), ]
        dynamics_data[[paste0(gene, "_", lineage_label)]] <- gene_df
      }
    }
  } else {
    for (gene in genes) {
      gene_expr <- as.numeric(expr_data[gene, ])
      if (params$normalize_dynamics) {
        gene_expr <- scale(gene_expr)[, 1]
      }

      gene_df <- data.frame(
        Gene       = gene,
        Pseudotime = pseudotime,
        Expression = gene_expr,
        Cell       = names(pseudotime),
        Lineage    = "All",
        stringsAsFactors = FALSE
      )

      # Remove cells with missing pseudotime
      gene_df <- gene_df[!is.na(gene_df$Pseudotime), ]

      dynamics_data[[gene]] <- gene_df
    }
  }

  combined_data <- do.call(rbind, dynamics_data)
  rownames(combined_data) <- NULL

  list(
    dynamics_data    = combined_data,
    genes            = genes,
    method           = params$dynamics_method,
    split_by_lineage = split_by_lineage
  )
}

analyze_lineage_branches <- function(trajectory_data, params) {

  # Initialize all output variables to avoid scope issues

  branch_stats <- NULL
  branch_de <- NULL
  branch_probs <- NULL

  if (trajectory_data$method == "slingshot") {
    # Extract lineage information from Slingshot results
    n_lineages <- trajectory_data$n_lineages
    lineage_ids <- trajectory_data$lineage_assignment
    pseudotime  <- trajectory_data$pseudotime

    cell_counts    <- as.numeric(table(lineage_ids))
    mean_pt        <- tapply(pseudotime, lineage_ids, mean, na.rm = TRUE)
    sd_pt          <- tapply(pseudotime, lineage_ids, sd, na.rm = TRUE)
    median_pt      <- tapply(pseudotime, lineage_ids, median, na.rm = TRUE)
    min_pt         <- tapply(pseudotime, lineage_ids, min, na.rm = TRUE)
    max_pt         <- tapply(pseudotime, lineage_ids, max, na.rm = TRUE)

    branch_stats <- data.frame(
      Branch          = paste("Lineage", seq_len(n_lineages)),
      N_Cells         = cell_counts,
      Mean_Pseudotime = as.numeric(mean_pt),
      SD_Pseudotime   = as.numeric(sd_pt),
      Median_Pseudotime = as.numeric(median_pt),
      Min_Pseudotime  = as.numeric(min_pt),
      Max_Pseudotime  = as.numeric(max_pt),
      stringsAsFactors = FALSE
    )

    # --- Compute branch probabilities for each cell ---
    # Use slingshot curve weights when the SCE object is available
    branch_probs <- NULL
    if (!is.null(trajectory_data$sce)) {
      tryCatch({
        pst_mat <- slingPseudotime(trajectory_data$sce)
        if (!is.null(dim(pst_mat)) && ncol(pst_mat) > 1) {
          # Weight = inverse of pseudotime-based distance, NA means cell not on lineage
          weights <- !is.na(pst_mat)
          row_sums <- rowSums(weights)
          row_sums[row_sums == 0] <- 1
          branch_probs <- sweep(weights * 1.0, 1, row_sums, "/")
          colnames(branch_probs) <- paste0("Prob_Lineage_", seq_len(ncol(branch_probs)))
          rownames(branch_probs) <- rownames(pst_mat)
        }
      }, error = function(e) NULL)
    }

    # --- Differential expression between branches ---
    branch_de <- NULL
    if (n_lineages >= 2 && !is.null(trajectory_data$sce)) {
      tryCatch({
        sce <- trajectory_data$sce
        expr_mat <- as.matrix(assay(sce, "logcounts"))
        de_pairs <- list()

        for (i in 1:(n_lineages - 1)) {
          for (j in (i + 1):n_lineages) {
            cells_i <- names(lineage_ids)[lineage_ids == i]
            cells_j <- names(lineage_ids)[lineage_ids == j]

            if (length(cells_i) < 3 || length(cells_j) < 3) next

            # Compute mean expression per branch
            mean_i <- rowMeans(expr_mat[, cells_i, drop = FALSE])
            mean_j <- rowMeans(expr_mat[, cells_j, drop = FALSE])
            log2fc <- mean_i - mean_j  # already log-transformed data

            # Wilcoxon rank-sum test on top variable genes for speed
            var_genes <- names(sort(apply(expr_mat[, c(cells_i, cells_j)], 1, var),
                                    decreasing = TRUE))[1:min(500, nrow(expr_mat))]

            pvals <- vapply(var_genes, function(g) {
              tryCatch(
                wilcox.test(expr_mat[g, cells_i], expr_mat[g, cells_j])$p.value,
                error = function(e) NA_real_
              )
            }, numeric(1))

            padj <- p.adjust(pvals, method = "BH")

            de_df <- data.frame(
              Gene     = var_genes,
              Log2FC   = log2fc[var_genes],
              PValue   = pvals,
              FDR      = padj,
              Comparison = paste0("Lineage_", i, "_vs_", j),
              stringsAsFactors = FALSE
            )
            de_df <- de_df[order(de_df$FDR), ]
            de_pairs[[paste0(i, "_vs_", j)]] <- de_df
          }
        }

        if (length(de_pairs) > 0) {
          branch_de <- do.call(rbind, de_pairs)
          rownames(branch_de) <- NULL
        }
      }, error = function(e) NULL)
    }

  } else if (trajectory_data$method == "monocle3") {
    # Monocle 3 branch analysis using the principal graph
    pseudotime <- trajectory_data$pseudotime

    if (!is.null(trajectory_data$cds)) {
      tryCatch({
        cds <- trajectory_data$cds
        pr_graph <- monocle3::principal_graph(cds)

        # Identify branch points (degree > 2) in the principal graph
        if (!is.null(pr_graph) && length(pr_graph) > 0) {
          g <- pr_graph[[1]]
          if (requireNamespace("igraph", quietly = TRUE)) {
            deg <- igraph::degree(g)
            branch_nodes <- names(deg)[deg > 2]
            leaf_nodes   <- names(deg)[deg == 1]
            n_branches   <- max(1, length(leaf_nodes))

            # Partition info
            partitions <- monocle3::partitions(cds)
            part_counts <- table(partitions)

            branch_stats <- data.frame(
              Branch          = paste("Partition", names(part_counts)),
              N_Cells         = as.numeric(part_counts),
              Mean_Pseudotime = tapply(pseudotime, partitions, mean, na.rm = TRUE),
              SD_Pseudotime   = tapply(pseudotime, partitions, sd, na.rm = TRUE),
              Median_Pseudotime = tapply(pseudotime, partitions, median, na.rm = TRUE),
              Min_Pseudotime  = tapply(pseudotime, partitions, min, na.rm = TRUE),
              Max_Pseudotime  = tapply(pseudotime, partitions, max, na.rm = TRUE),
              stringsAsFactors = FALSE
            )
          } else {
            stop("igraph required for Monocle 3 graph analysis")
          }
        } else {
          stop("No principal graph found")
        }
      }, error = function(e) {
        branch_stats <<- data.frame(
          Branch          = "Main trajectory (graph analysis failed)",
          N_Cells         = length(pseudotime),
          Mean_Pseudotime = mean(pseudotime, na.rm = TRUE),
          SD_Pseudotime   = sd(pseudotime, na.rm = TRUE),
          Median_Pseudotime = median(pseudotime, na.rm = TRUE),
          Min_Pseudotime  = min(pseudotime, na.rm = TRUE),
          Max_Pseudotime  = max(pseudotime, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    } else {
      branch_stats <- data.frame(
        Branch          = "Main trajectory",
        N_Cells         = length(pseudotime),
        Mean_Pseudotime = mean(pseudotime, na.rm = TRUE),
        SD_Pseudotime   = sd(pseudotime, na.rm = TRUE),
        Median_Pseudotime = median(pseudotime, na.rm = TRUE),
        Min_Pseudotime  = min(pseudotime, na.rm = TRUE),
        Max_Pseudotime  = max(pseudotime, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }

  list(
    branch_statistics  = branch_stats,
    branch_de          = branch_de,
    branch_probs       = branch_probs,
    method             = trajectory_data$method
  )
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

create_slingshot_plot <- function(seurat_obj, sling_results, params) {
  red <- sling_results$reduction
  embedding <- Embeddings(seurat_obj, red)[, 1:2]

  # Safely map pseudotime and lineage, handling potential name mismatches
  cell_names <- rownames(embedding)
  pt_vals <- sling_results$pseudotime[cell_names]
  lin_vals <- sling_results$lineage_assignment[cell_names]
  # Fall back to positional mapping if names don't match
  if (all(is.na(pt_vals)) && length(sling_results$pseudotime) == length(cell_names)) {
    pt_vals <- sling_results$pseudotime
    lin_vals <- sling_results$lineage_assignment
  }

  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Pseudotime = pt_vals,
    Lineage = factor(lin_vals),
    CellID = cell_names,
    stringsAsFactors = FALSE
  )

  # Create base plot
  p <- plot_ly(
    data = plot_data,
    x = ~Dim1,
    y = ~Dim2,
    color = ~Pseudotime,
    colors = viridis::viridis(100),
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.8),
    text = ~CellID,
    hovertemplate = paste0(
      "Cell: %{text}<br>",
      toupper(red), "1: %{x}<br>",
      toupper(red), "2: %{y}<br>",
      "Pseudotime: %{marker.color:.3f}<extra></extra>"
    )
  ) %>% colorbar(title = "Pseudotime")

  # Add curves if available
  if (!is.null(sling_results$curves_data) && isTRUE(params$sling_show_curves)) {
    curve_data <- sling_results$curves_data
    if (ncol(curve_data) >= 3) {
      for (lineage in unique(curve_data$lineage)) {
        ld <- curve_data[curve_data$lineage == lineage, ]
        p <- p %>% add_trace(
          x = ld[, 1],
          y = ld[, 2],
          type = "scatter",
          mode = "lines",
          line = list(width = 3, color = "black"),
          showlegend = FALSE,
          hoverinfo = "none",
          inherit = FALSE
        )
      }
    }
  }

  p %>% layout(
    title = "Slingshot Trajectory Analysis",
    xaxis = list(title = paste0(toupper(red), "_1")),
    yaxis = list(title = paste0(toupper(red), "_2"))
  )
}

create_monocle3_plot <- function(seurat_obj, cds_results, params) {
  red <- tolower(params$reduction_for_traj)
  embedding <- Embeddings(seurat_obj, red)[, 1:2]
  cell_names <- rownames(embedding)

  pt_vals <- cds_results$pseudotime[cell_names]
  part_vals <- cds_results$partition[cell_names]
  if (all(is.na(pt_vals)) && length(cds_results$pseudotime) == length(cell_names)) {
    pt_vals <- cds_results$pseudotime
    part_vals <- cds_results$partition
  }

  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Pseudotime = pt_vals,
    Partition = factor(part_vals),
    CellID = cell_names,
    stringsAsFactors = FALSE
  )

  plot_ly(
    data = plot_data,
    x = ~Dim1,
    y = ~Dim2,
    color = ~Pseudotime,
    colors = viridis::plasma(100),
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.8),
    text = ~CellID,
    hovertemplate = paste0(
      "Cell: %{text}<br>",
      toupper(red), "1: %{x}<br>",
      toupper(red), "2: %{y}<br>",
      "Pseudotime: %{marker.color:.3f}<br>",
      "Partition: %{customdata}<extra></extra>"
    ),
    customdata = ~Partition
  ) %>%
    colorbar(title = "Pseudotime") %>%
    layout(
      title = "Monocle 3 Trajectory Analysis",
      xaxis = list(title = paste0(toupper(red), "_1")),
      yaxis = list(title = paste0(toupper(red), "_2"))
    )
}

create_dpt_plot <- function(seurat_obj, dpt_results, params) {
  plot_data <- data.frame(
    DC1 = dpt_results$dc1,
    DC2 = dpt_results$dc2,
    Pseudotime = dpt_results$pseudotime,
    CellID = if (!is.null(names(dpt_results$pseudotime))) names(dpt_results$pseudotime)
             else paste0("cell_", seq_along(dpt_results$pseudotime)),
    stringsAsFactors = FALSE
  )

  plot_ly(
    data = plot_data,
    x = ~DC1,
    y = ~DC2,
    color = ~Pseudotime,
    colors = viridis::cividis(100),
    type = "scatter",
    mode = "markers",
    marker = list(size = 4, opacity = 0.8),
    text = ~CellID,
    hovertemplate = "Cell: %{text}<br>DC1: %{x}<br>DC2: %{y}<br>Pseudotime: %{marker.color:.3f}<extra></extra>"
  ) %>%
    colorbar(title = "Pseudotime") %>%
    layout(
      title = "Diffusion Pseudotime Analysis",
      xaxis = list(title = "Diffusion Component 1"),
      yaxis = list(title = "Diffusion Component 2")
    )
}

create_comparison_heatmap <- function(comparison_results) {
  cor_matrix <- comparison_results$correlation_matrix
  
  # Convert to long format for plotting
  cor_long <- expand.grid(
    Method1 = rownames(cor_matrix),
    Method2 = colnames(cor_matrix)
  )
  cor_long$Correlation <- as.vector(cor_matrix)
  
  plot_ly(
    data = cor_long,
    x = ~Method1,
    y = ~Method2,
    z = ~Correlation,
    type = "heatmap",
    colorscale = "RdBu",
    hovertemplate = "Methods: %{x} vs %{y}<br>Correlation: %{z:.3f}<extra></extra>"
  ) %>%
    layout(
      title = "Method Correlation Matrix",
      xaxis = list(title = "Method"),
      yaxis = list(title = "Method")
    )
}

create_comparison_violin <- function(comparison_results) {
  pt_matrix <- comparison_results$pseudotime_matrix
  
  # Convert to long format
  violin_data <- tidyr::pivot_longer(
    as.data.frame(pt_matrix),
    cols = everything(),
    names_to = "Method",
    values_to = "Pseudotime"
  )
  
  plot_ly(
    data = violin_data,
    x = ~Method,
    y = ~Pseudotime,
    color = ~Method,
    type = "violin",
    box = list(visible = TRUE)
  ) %>%
    layout(
      title = "Pseudotime Distributions by Method",
      xaxis = list(title = "Method"),
      yaxis = list(title = "Pseudotime")
    )
}

create_sidebyside_comparison <- function(seurat_obj, comparison_results, params) {
  tryCatch({
    pt_matrix <- comparison_results$pseudotime_matrix
    methods   <- comparison_results$methods
    common_cells <- comparison_results$common_cells

    # Get embedding coordinates
    red <- tolower(params$reduction_for_traj)
    embedding <- Embeddings(seurat_obj, red)[common_cells, 1:2]
    dim_names <- colnames(embedding)

    plot_list <- lapply(seq_along(methods), function(i) {
      method <- methods[i]
      pt_vals <- pt_matrix[, method]

      df <- data.frame(
        Dim1 = embedding[, 1],
        Dim2 = embedding[, 2],
        Pseudotime = pt_vals
      )

      plot_ly(
        data = df,
        x = ~Dim1,
        y = ~Dim2,
        color = ~Pseudotime,
        colors = viridis::viridis(256),
        type = "scatter",
        mode = "markers",
        marker = list(size = 4, opacity = 0.7),
        text = ~paste0("Pseudotime: ", round(Pseudotime, 3)),
        hoverinfo = "text",
        showlegend = (i == 1)
      ) %>%
        layout(
          annotations = list(
            list(
              text = method,
              xref = "paper", yref = "paper",
              x = 0.5, y = 1.05,
              showarrow = FALSE,
              font = list(size = 14, weight = "bold")
            )
          ),
          xaxis = list(title = dim_names[1]),
          yaxis = list(title = dim_names[2])
        )
    })

    plotly::subplot(plot_list, nrows = 1, shareY = TRUE, titleX = TRUE, titleY = TRUE) %>%
      layout(
        title = list(text = "Side-by-Side Pseudotime Comparison", y = 0.98),
        coloraxis = list(colorscale = "Viridis", colorbar = list(title = "Pseudotime"))
      )
  }, error = function(e) {
    plotly_empty() %>%
      layout(
        title = "Side-by-side comparison error",
        annotations = list(
          text = paste("Error:", e$message),
          showarrow = FALSE, x = 0.5, y = 0.5
        )
      )
  })
}

create_gene_dynamics_plot <- function(dynamics_results, params) {
  plot_data <- dynamics_results$dynamics_data
  split_mode <- isTRUE(dynamics_results$split_by_lineage)

  # Determine unique groups: gene only, or gene x lineage
  if (split_mode && "Lineage" %in% colnames(plot_data)) {
    group_keys <- unique(plot_data[, c("Gene", "Lineage"), drop = FALSE])
  } else {
    group_keys <- data.frame(
      Gene    = unique(plot_data$Gene),
      Lineage = "All",
      stringsAsFactors = FALSE
    )
  }

  # Build a colour palette across all groups
  n_groups <- nrow(group_keys)
  palette  <- .get_palette(n_groups)

  p <- plot_ly()

  for (gi in seq_len(n_groups)) {
    gene    <- group_keys$Gene[gi]
    lineage <- group_keys$Lineage[gi]
    grp_col <- palette[gi]

    if (split_mode && "Lineage" %in% colnames(plot_data)) {
      gene_data <- plot_data[plot_data$Gene == gene & plot_data$Lineage == lineage, ]
    } else {
      gene_data <- plot_data[plot_data$Gene == gene, ]
    }
    if (nrow(gene_data) < 5) next

    trace_label <- if (split_mode) paste0(gene, " (", lineage, ")") else gene

    # Add raw points
    p <- p %>% add_trace(
      data = gene_data,
      x = ~Pseudotime,
      y = ~Expression,
      type = "scatter",
      mode = "markers",
      name = paste(trace_label, "data"),
      marker = list(size = 3, opacity = 0.4, color = grp_col),
      showlegend = FALSE
    )

    # Fit smooth curve with confidence bands
    pred_grid <- data.frame(
      Pseudotime = seq(min(gene_data$Pseudotime), max(gene_data$Pseudotime), length.out = 100)
    )

    smooth_pred <- NULL

    loess_span <- if (!is.null(params$smooth_span) && is.numeric(params$smooth_span)) params$smooth_span else 0.3
    gam_knots <- if (!is.null(params$n_knots) && is.numeric(params$n_knots)) as.integer(params$n_knots) else 6L

    if (params$smooth_method == "loess") {
      smooth_fit <- loess(Expression ~ Pseudotime, data = gene_data, span = loess_span)
      raw_pred <- predict(smooth_fit, pred_grid$Pseudotime, se = TRUE)
      smooth_pred <- data.frame(
        Pseudotime = pred_grid$Pseudotime,
        Expression = raw_pred$fit,
        SE = raw_pred$se.fit
      )
    } else if (params$smooth_method == "gam") {
      if (!requireNamespace("mgcv", quietly = TRUE)) {
        smooth_fit <- loess(Expression ~ Pseudotime, data = gene_data, span = loess_span)
        raw_pred <- predict(smooth_fit, pred_grid$Pseudotime, se = TRUE)
        smooth_pred <- data.frame(
          Pseudotime = pred_grid$Pseudotime,
          Expression = raw_pred$fit,
          SE = raw_pred$se.fit
        )
      } else {
        smooth_fit <- mgcv::gam(Expression ~ s(Pseudotime, k = gam_knots), data = gene_data)
        raw_pred <- mgcv::predict.gam(smooth_fit, newdata = pred_grid, se.fit = TRUE)
        smooth_pred <- data.frame(
          Pseudotime = pred_grid$Pseudotime,
          Expression = raw_pred$fit,
          SE = raw_pred$se.fit
        )
      }
    } else if (params$smooth_method == "lm") {
      smooth_fit <- lm(Expression ~ poly(Pseudotime, 3), data = gene_data)
      raw_pred <- predict(smooth_fit, newdata = pred_grid, se.fit = TRUE)
      smooth_pred <- data.frame(
        Pseudotime = pred_grid$Pseudotime,
        Expression = raw_pred$fit,
        SE = raw_pred$se.fit
      )
    }

    if (!is.null(smooth_pred)) {
      # Add smooth curve trace
      p <- p %>% add_trace(
        data = smooth_pred,
        x = ~Pseudotime,
        y = ~Expression,
        type = "scatter",
        mode = "lines",
        name = trace_label,
        line = list(width = 3, color = grp_col)
      )

      # Add confidence bands if requested
      if (isTRUE(params$show_confidence) && !all(is.na(smooth_pred$SE))) {
        upper <- smooth_pred$Expression + 1.96 * smooth_pred$SE
        lower <- smooth_pred$Expression - 1.96 * smooth_pred$SE

        rgb_vals <- paste(col2rgb(grp_col), collapse = ",")

        p <- p %>% add_trace(
          x = smooth_pred$Pseudotime,
          y = upper,
          type = "scatter",
          mode = "lines",
          line = list(width = 0),
          showlegend = FALSE,
          hoverinfo = "skip",
          name = paste(trace_label, "upper")
        ) %>% add_trace(
          x = smooth_pred$Pseudotime,
          y = lower,
          type = "scatter",
          mode = "lines",
          line = list(width = 0),
          fill = "tonexty",
          fillcolor = paste0("rgba(", rgb_vals, ",0.2)"),
          showlegend = FALSE,
          hoverinfo = "skip",
          name = paste(trace_label, "CI")
        )
      }
    }
  }

  title_text <- if (split_mode) {
    "Gene Expression Dynamics Along Pseudotime (Per-Lineage)"
  } else {
    "Gene Expression Dynamics Along Pseudotime"
  }

  p %>% layout(
    title = title_text,
    xaxis = list(title = "Pseudotime"),
    yaxis = list(title = "Expression Level")
  )
}

create_branch_structure_plot <- function(branch_results) {
  tryCatch({
    branch_stats <- branch_results$branch_statistics
    method <- branch_results$method

    if (method == "slingshot" && nrow(branch_stats) > 1) {
      # Tree-like visualization for slingshot lineages
      n_branches <- nrow(branch_stats)

      # Build a simple tree: root node at top, branches fan out below
      # Root node
      node_x <- c(0.5)
      node_y <- c(1.0)
      node_labels <- c("Root")
      node_sizes  <- c(20)
      node_colors <- c(0)

      # Branch nodes arranged in a fan below the root
      branch_x <- seq(0.1, 0.9, length.out = n_branches)
      branch_y <- rep(0.3, n_branches)

      cell_counts <- as.numeric(branch_stats$N_Cells)
      mean_pts    <- as.numeric(branch_stats$Mean_Pseudotime)

      # Normalize node sizes: scale cell counts to marker sizes [15, 60]
      size_range <- range(cell_counts)
      if (diff(size_range) == 0) {
        node_sz <- rep(35, n_branches)
      } else {
        node_sz <- 15 + 45 * (cell_counts - size_range[1]) / diff(size_range)
      }

      node_x      <- c(node_x, branch_x)
      node_y      <- c(node_y, branch_y)
      node_labels <- c(node_labels, as.character(branch_stats$Branch))
      node_sizes  <- c(node_sizes, node_sz)
      node_colors <- c(0, mean_pts)

      # Edges from root to each branch
      edge_x <- as.vector(rbind(rep(0.5, n_branches), branch_x, rep(NA, n_branches)))
      edge_y <- as.vector(rbind(rep(1.0, n_branches), branch_y, rep(NA, n_branches)))

      hover_text <- c(
        "Root (branch point)",
        paste0(branch_stats$Branch,
               "<br>Cells: ", cell_counts,
               "<br>Mean Pseudotime: ", round(mean_pts, 3))
      )

      p <- plot_ly() %>%
        add_trace(
          x = edge_x, y = edge_y,
          type = "scatter", mode = "lines",
          line = list(color = "grey", width = 2),
          hoverinfo = "none", showlegend = FALSE
        ) %>%
        add_trace(
          x = node_x, y = node_y,
          type = "scatter", mode = "markers+text",
          marker = list(
            size = node_sizes,
            color = node_colors,
            colorscale = "Viridis",
            colorbar = list(title = "Mean\nPseudotime"),
            line = list(color = "black", width = 1)
          ),
          text = node_labels,
          textposition = "bottom center",
          hovertext = hover_text,
          hoverinfo = "text",
          showlegend = FALSE
        ) %>%
        layout(
          title = "Lineage Branch Structure",
          xaxis = list(visible = FALSE, range = c(-0.05, 1.05)),
          yaxis = list(visible = FALSE, range = c(0.0, 1.15)),
          plot_bgcolor = "white"
        )

      p

    } else {
      # Bar chart with cells per branch and mean pseudotime on secondary axis
      cell_counts <- as.numeric(branch_stats$N_Cells)
      mean_pts    <- as.numeric(branch_stats$Mean_Pseudotime)

      plot_ly() %>%
        add_trace(
          x = branch_stats$Branch,
          y = cell_counts,
          type = "bar",
          name = "Number of Cells",
          marker = list(color = "steelblue"),
          hovertemplate = "%{x}<br>Cells: %{y}<extra></extra>"
        ) %>%
        add_trace(
          x = branch_stats$Branch,
          y = mean_pts,
          type = "scatter",
          mode = "markers+lines",
          name = "Mean Pseudotime",
          yaxis = "y2",
          marker = list(color = "firebrick", size = 10),
          line = list(color = "firebrick", width = 2),
          hovertemplate = "%{x}<br>Mean Pseudotime: %{y:.3f}<extra></extra>"
        ) %>%
        layout(
          title = "Branch Statistics",
          xaxis = list(title = "Branch"),
          yaxis = list(title = "Number of Cells", side = "left"),
          yaxis2 = list(
            title = "Mean Pseudotime",
            side = "right",
            overlaying = "y",
            showgrid = FALSE
          ),
          legend = list(x = 0.7, y = 1.1, orientation = "h")
        )
    }
  }, error = function(e) {
    plotly_empty() %>%
      layout(
        title = "Branch structure error",
        annotations = list(
          text = paste("Error:", e$message),
          showarrow = FALSE, x = 0.5, y = 0.5
        )
      )
  })
}

# =============================================================================
# RNA VELOCITY FUNCTIONS
# =============================================================================

compute_rna_velocity <- function(seurat_obj, params) {
  tryCatch({
    source_type <- params$velocity_source
    calc_method <- params$velocity_method

    spliced   <- NULL
    unspliced <- NULL
    velocity  <- NULL

    # ----- Source: Seurat misc slot -----
    if (source_type == "misc") {
      misc <- seurat_obj@misc
      if (!is.null(misc$velocity)) {
        velocity  <- misc$velocity
        spliced   <- misc$spliced
        unspliced <- misc$unspliced
        # Check provenance
        if (is.null(misc$velocity_metadata)) {
          showNotification(
            "Velocity data found but no provenance metadata (misc$velocity_metadata). Results may not be reproducible. Consider storing source tool and parameters.",
            type = "warning", duration = 10
          )
        }
      } else {
        stop("No velocity data found in Seurat object @misc slot. ",
             "Expected 'velocity', 'spliced', and 'unspliced' matrices.")
      }
    }

    # ----- Source: upload h5ad via reticulate -----
    if (source_type == "h5ad") {
      if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("reticulate package required to read .h5ad files")
      }
      h5ad_path <- params$velocity_h5ad$datapath
      if (is.null(h5ad_path)) {
        stop("No .h5ad file uploaded")
      }

      sc <- reticulate::import("scanpy", convert = FALSE)
      adata <- sc$read_h5ad(h5ad_path)

      # Extract velocity matrix if present
      if ("velocity" %in% reticulate::py_to_r(adata$layers$keys())) {
        velocity <- reticulate::py_to_r(adata$layers[["velocity"]])
        rownames(velocity) <- reticulate::py_to_r(adata$obs_names$to_list())
        colnames(velocity) <- reticulate::py_to_r(adata$var_names$to_list())
        velocity <- t(velocity)  # genes x cells
      }
      if ("spliced" %in% reticulate::py_to_r(adata$layers$keys())) {
        spliced <- reticulate::py_to_r(adata$layers[["spliced"]])
        rownames(spliced) <- reticulate::py_to_r(adata$obs_names$to_list())
        colnames(spliced) <- reticulate::py_to_r(adata$var_names$to_list())
        spliced <- t(spliced)
      }
      if ("unspliced" %in% reticulate::py_to_r(adata$layers$keys())) {
        unspliced <- reticulate::py_to_r(adata$layers[["unspliced"]])
        rownames(unspliced) <- reticulate::py_to_r(adata$obs_names$to_list())
        colnames(unspliced) <- reticulate::py_to_r(adata$var_names$to_list())
        unspliced <- t(unspliced)
      }
    }

    # ----- Source: compute from spliced/unspliced assays -----
    if (source_type == "compute") {
      assay_names <- Assays(seurat_obj)
      # Look for spliced / unspliced assays
      if ("spliced" %in% assay_names && "unspliced" %in% assay_names) {
        spliced   <- safe_get_assay_data(seurat_obj, "counts", assay = "spliced")
        unspliced <- safe_get_assay_data(seurat_obj, "counts", assay = "unspliced")
      } else if (!is.null(seurat_obj@misc$spliced) && !is.null(seurat_obj@misc$unspliced)) {
        spliced   <- seurat_obj@misc$spliced
        unspliced <- seurat_obj@misc$unspliced
      } else {
        stop("Spliced and unspliced count matrices not found. ",
             "Provide them as separate assays ('spliced', 'unspliced') or in @misc.")
      }
    }

    # ----- Compute velocity if not already provided -----
    if (is.null(velocity) && !is.null(spliced) && !is.null(unspliced)) {
      # Align genes and cells
      common_genes <- intersect(rownames(spliced), rownames(unspliced))
      common_cells <- intersect(colnames(spliced), colnames(unspliced))
      common_cells <- intersect(common_cells, colnames(seurat_obj))

      if (length(common_genes) < 10 || length(common_cells) < 10) {
        stop("Not enough shared genes/cells between spliced and unspliced matrices")
      }

      spliced   <- as.matrix(spliced[common_genes, common_cells])
      unspliced <- as.matrix(unspliced[common_genes, common_cells])

      if (calc_method == "steady_state") {
        # Steady-state model: velocity ~ unspliced - gamma * spliced
        # Estimate gamma per gene as slope of unspliced vs spliced
        gamma <- vapply(seq_len(nrow(spliced)), function(i) {
          s_i <- spliced[i, ]
          u_i <- unspliced[i, ]
          if (sd(s_i) == 0) return(0)
          coef(lm(u_i ~ s_i + 0))[1]
        }, numeric(1))
        names(gamma) <- common_genes

        # Velocity = unspliced - gamma * spliced
        velocity <- unspliced - sweep(spliced, 1, gamma, "*")

      } else {
        # Dynamical model: attempt via reticulate + scvelo
        if (requireNamespace("reticulate", quietly = TRUE)) {
          tryCatch({
            scv <- reticulate::import("scvelo", convert = FALSE)
            sc  <- reticulate::import("scanpy", convert = FALSE)
            np  <- reticulate::import("numpy", convert = FALSE)

            # Build a minimal AnnData from spliced/unspliced
            X_r <- t(spliced)
            adata <- sc$AnnData(X = np$array(X_r))
            adata$layers[["spliced"]]   <- np$array(t(spliced))
            adata$layers[["unspliced"]] <- np$array(t(unspliced))

            scv$pp$filter_and_normalize(adata, min_shared_counts = as.integer(5))
            scv$pp$moments(adata, n_pcs = as.integer(30), n_neighbors = as.integer(30))
            scv$tl$recover_dynamics(adata)
            scv$tl$velocity(adata, mode = "dynamical")

            velocity <- reticulate::py_to_r(adata$layers[["velocity"]])
            velocity <- t(velocity)
            rownames(velocity) <- common_genes[seq_len(nrow(velocity))]
            colnames(velocity) <- common_cells[seq_len(ncol(velocity))]
          }, error = function(e) {
            warning("scVelo dynamical model failed, falling back to steady-state: ", e$message)
            gamma <- vapply(seq_len(nrow(spliced)), function(i) {
              s_i <- spliced[i, ]
              u_i <- unspliced[i, ]
              if (sd(s_i) == 0) return(0)
              coef(lm(u_i ~ s_i + 0))[1]
            }, numeric(1))
            velocity <<- unspliced - sweep(spliced, 1, gamma, "*")
          })
        } else {
          warning("reticulate not available, falling back to steady-state model")
          gamma <- vapply(seq_len(nrow(spliced)), function(i) {
            s_i <- spliced[i, ]
            u_i <- unspliced[i, ]
            if (sd(s_i) == 0) return(0)
            coef(lm(u_i ~ s_i + 0))[1]
          }, numeric(1))
          velocity <- unspliced - sweep(spliced, 1, gamma, "*")
        }
      }
    }

    if (is.null(velocity)) {
      stop("Could not obtain or compute velocity matrix")
    }

    # ----- Project velocity onto embedding -----
    cells <- intersect(colnames(velocity), colnames(seurat_obj))
    red   <- tolower(params$velocity_embedding)
    if (!(red %in% names(seurat_obj@reductions))) {
      red <- .get_best_2d(seurat_obj)
      if (is.null(red)) stop("No 2D embedding found")
    }

    embedding <- Embeddings(seurat_obj, red)[cells, 1:2]
    vel_mat   <- as.matrix(velocity[, cells])

    # Project velocity onto embedding using correlation-based approach
    # For each cell, compute correlation of velocity vector with
    # direction to neighbors in embedding space
    n_cells  <- length(cells)
    arrow_dx <- numeric(n_cells)
    arrow_dy <- numeric(n_cells)

    # Use PCA loadings if available for projection, else use kNN approach
    pca_embed <- NULL
    if ("pca" %in% names(seurat_obj@reductions)) {
      pca_embed <- Embeddings(seurat_obj, "pca")[cells, , drop = FALSE]
    }

    # kNN-based projection onto 2D embedding
    k_nn <- min(30, n_cells - 1)

    # Pre-fetch expression data once (avoid repeated subsetting in loop)
    expr_data_all <- as.matrix(safe_get_assay_data(seurat_obj, "data")[
      rownames(vel_mat), cells, drop = FALSE
    ])

    for (ci in seq_len(n_cells)) {
      if (ci %% 100 == 0) {
        incProgress(0.5 / n_cells, detail = sprintf("Projecting cell %d/%d", ci, n_cells))
      }

      cell_embed <- embedding[ci, ]
      diffs <- sweep(embedding, 2, cell_embed, "-")
      dists <- sqrt(rowSums(diffs^2))
      dists[ci] <- Inf
      nn_idx <- order(dists)[1:k_nn]

      # Velocity correlation with expression differences for neighbors
      vel_cell <- vel_mat[, ci]
      expr_diffs <- expr_data_all[, nn_idx, drop = FALSE] -
                    as.numeric(expr_data_all[, ci])

      # Cosine similarity between velocity and expression difference
      cos_sim <- apply(expr_diffs, 2, function(d) {
        norm_v <- sqrt(sum(vel_cell^2))
        norm_d <- sqrt(sum(d^2))
        if (norm_v == 0 || norm_d == 0) return(0)
        sum(vel_cell * d) / (norm_v * norm_d)
      })

      # Weight neighbors by transition probability (softmax of cosine sim)
      weights <- exp(cos_sim) / sum(exp(cos_sim))

      # Weighted average displacement in embedding space
      displacements <- diffs[nn_idx, , drop = FALSE]
      arrow_dx[ci] <- sum(weights * displacements[, 1])
      arrow_dy[ci] <- sum(weights * displacements[, 2])
    }

    rm(expr_data_all)
    gc()

    # Compute confidence and magnitude per cell
    magnitude <- sqrt(arrow_dx^2 + arrow_dy^2)
    # Confidence: consistency of velocity direction among neighbors
    confidence <- vapply(seq_len(n_cells), function(ci) {
      cell_embed <- embedding[ci, ]
      diffs <- sweep(embedding, 2, cell_embed, "-")
      dists <- sqrt(rowSums(diffs^2))
      dists[ci] <- Inf
      nn_idx <- order(dists)[1:k_nn]

      nn_dx <- arrow_dx[nn_idx]
      nn_dy <- arrow_dy[nn_idx]

      # Average cosine similarity of arrow directions with this cell's arrow
      my_dir <- c(arrow_dx[ci], arrow_dy[ci])
      my_norm <- sqrt(sum(my_dir^2))
      if (my_norm == 0) return(0)

      cos_sims <- vapply(seq_along(nn_idx), function(ni) {
        nn_dir <- c(nn_dx[ni], nn_dy[ni])
        nn_norm <- sqrt(sum(nn_dir^2))
        if (nn_norm == 0) return(0)
        sum(my_dir * nn_dir) / (my_norm * nn_norm)
      }, numeric(1))

      mean(cos_sims, na.rm = TRUE)
    }, numeric(1))

    # Clip confidence to [0, 1]
    confidence <- pmax(0, pmin(1, (confidence + 1) / 2))

    names(magnitude)  <- cells
    names(confidence) <- cells

    list(
      success     = TRUE,
      arrow_dx    = arrow_dx,
      arrow_dy    = arrow_dy,
      embedding   = embedding,
      cells       = cells,
      magnitude   = magnitude,
      confidence  = confidence,
      source      = source_type,
      calc_method = calc_method,
      reduction   = red,
      error       = NULL
    )

  }, error = function(e) {
    list(
      success     = FALSE,
      arrow_dx    = NULL,
      arrow_dy    = NULL,
      embedding   = NULL,
      cells       = NULL,
      magnitude   = NULL,
      confidence  = NULL,
      source      = NULL,
      calc_method = NULL,
      reduction   = NULL,
      error       = e$message
    )
  })
}

create_velocity_stream_plot <- function(seurat_obj, velocity_results, params) {
  tryCatch({
    embedding  <- velocity_results$embedding
    arrow_dx   <- velocity_results$arrow_dx
    arrow_dy   <- velocity_results$arrow_dy
    cells      <- velocity_results$cells
    magnitude  <- velocity_results$magnitude
    confidence <- velocity_results$confidence
    red        <- velocity_results$reduction

    arrow_scale <- if (!is.null(params$velocity_arrow_scale)) params$velocity_arrow_scale else 1.0
    grid_n      <- if (!is.null(params$velocity_grid_n)) params$velocity_grid_n else 30

    # Grid-based arrow aggregation for stream-like visualization
    x_range <- range(embedding[, 1])
    y_range <- range(embedding[, 2])
    x_breaks <- seq(x_range[1], x_range[2], length.out = grid_n + 1)
    y_breaks <- seq(y_range[1], y_range[2], length.out = grid_n + 1)

    grid_cx <- (x_breaks[-1] + x_breaks[-length(x_breaks)]) / 2
    grid_cy <- (y_breaks[-1] + y_breaks[-length(y_breaks)]) / 2

    # Assign cells to grid bins
    x_bin <- findInterval(embedding[, 1], x_breaks, all.inside = TRUE)
    y_bin <- findInterval(embedding[, 2], y_breaks, all.inside = TRUE)

    # Aggregate arrows per grid cell
    grid_arrows <- list()
    for (i in seq_len(grid_n)) {
      for (j in seq_len(grid_n)) {
        mask <- x_bin == i & y_bin == j
        if (sum(mask) < 2) next
        mean_dx <- mean(arrow_dx[mask])
        mean_dy <- mean(arrow_dy[mask])
        mean_mag <- sqrt(mean_dx^2 + mean_dy^2)
        if (mean_mag < 1e-10) next

        grid_arrows[[length(grid_arrows) + 1]] <- list(
          x = grid_cx[i], y = grid_cy[j],
          dx = mean_dx, dy = mean_dy,
          mag = mean_mag,
          n_cells = sum(mask)
        )
      }
    }

    # Normalize arrow lengths for display
    if (length(grid_arrows) > 0) {
      all_mags <- vapply(grid_arrows, function(a) a$mag, numeric(1))
      max_mag  <- max(all_mags)
      step_x   <- diff(x_range) / grid_n
      step_y   <- diff(y_range) / grid_n
      base_len <- min(step_x, step_y) * 0.8 * arrow_scale
    }

    # Get cluster colors for background points
    cluster_col <- NULL
    if (!is.null(params$cluster_col) && params$cluster_col %in% colnames(seurat_obj@meta.data)) {
      cluster_col <- seurat_obj@meta.data[cells, params$cluster_col]
    }

    plot_df <- data.frame(
      x = embedding[, 1],
      y = embedding[, 2],
      confidence = confidence,
      magnitude  = magnitude,
      cell = cells
    )
    if (!is.null(cluster_col)) {
      plot_df$cluster <- as.character(cluster_col)
    }

    # Create base scatter plot colored by confidence
    if (!is.null(cluster_col)) {
      p <- plot_ly(
        data = plot_df,
        x = ~x, y = ~y,
        color = ~cluster,
        colors = .get_palette(length(unique(plot_df$cluster))),
        type = "scatter",
        mode = "markers",
        marker = list(size = 4, opacity = 0.5),
        text = ~paste0("Cell: ", cell, "<br>Cluster: ", cluster,
                        "<br>Confidence: ", round(confidence, 3),
                        "<br>Magnitude: ", round(magnitude, 4)),
        hoverinfo = "text"
      )
    } else {
      p <- plot_ly(
        data = plot_df,
        x = ~x, y = ~y,
        color = ~confidence,
        colors = viridis::viridis(100),
        type = "scatter",
        mode = "markers",
        marker = list(size = 4, opacity = 0.5),
        text = ~paste0("Cell: ", cell,
                        "<br>Confidence: ", round(confidence, 3),
                        "<br>Magnitude: ", round(magnitude, 4)),
        hoverinfo = "text"
      )
    }

    # Add grid-based velocity arrows using annotations
    arrow_annotations <- list()
    for (ai in seq_along(grid_arrows)) {
      a <- grid_arrows[[ai]]
      # Normalize direction and scale
      norm_factor <- base_len / max_mag
      dx_scaled <- a$dx * norm_factor
      dy_scaled <- a$dy * norm_factor

      arrow_annotations[[ai]] <- list(
        x = a$x + dx_scaled,
        y = a$y + dy_scaled,
        ax = a$x,
        ay = a$y,
        xref = "x", yref = "y",
        axref = "x", ayref = "y",
        showarrow = TRUE,
        arrowhead = 2,
        arrowsize = 1,
        arrowwidth = 1.5,
        arrowcolor = "rgba(0,0,0,0.6)"
      )
    }

    p <- p %>% layout(
      title = paste("RNA Velocity Stream Plot (", toupper(red), ")"),
      xaxis = list(title = paste0(toupper(red), "_1")),
      yaxis = list(title = paste0(toupper(red), "_2")),
      annotations = arrow_annotations
    )

    p

  }, error = function(e) {
    plotly_empty() %>%
      layout(
        title = "Velocity stream plot error",
        annotations = list(
          text = paste("Error:", e$message),
          showarrow = FALSE, x = 0.5, y = 0.5
        )
      )
  })
}
