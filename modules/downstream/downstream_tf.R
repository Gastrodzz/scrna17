# Downstream TF & GRN — Transcription Factor Activity + Gene Regulatory Networks
# modules/downstream/downstream_tf.R

# Built-in TF regulon database (DOROTHEA-like)
# Tries decoupleR's CollecTRI network first; falls back to built-in ~20 TFs
.get_tf_regulon_db <- function(species = "human") {

  # Attempt to load decoupleR regulons (CollecTRI — comprehensive, peer-reviewed)
  if (requireNamespace("decoupleR", quietly = TRUE)) {
    tryCatch({
      organism <- if (species == "mouse") "mouse" else "human"
      net <- decoupleR::get_collectri(organism = organism, split_complexes = FALSE)
      if (is.data.frame(net) && nrow(net) > 100) {
        # Standardize column names to match built-in format
        colnames(net)[colnames(net) == "source"] <- "tf"
        colnames(net)[colnames(net) == "mor"]    <- "weight"
        net <- net[, c("tf", "target", "weight")]
        net$weight[is.na(net$weight)] <- 1
        message("[TF] Loaded ", nrow(net), " TF-target pairs from decoupleR CollecTRI (", organism, ")")
        return(net)
      }
    }, error = function(e) {
      message("[TF] decoupleR available but CollecTRI load failed: ", e$message,
              ". Falling back to built-in database.")
    })
  }

  # MINIMAL FALLBACK DATABASE - NOT PUBLICATION-GRADE
  # Only used when decoupleR::get_collectri() is unavailable
  # Contains ~20 well-known TFs for exploratory analysis only
  showNotification(
    "decoupleR not available. Using minimal fallback TF database (~20 TFs). Results are exploratory only. Install decoupleR for publication-grade analysis: BiocManager::install('decoupleR')",
    type = "warning", duration = 15
  )

  # Each entry: TF, target gene, weight (+1 activation, -1 repression)
  regulons <- list(
    TP53 = list(
      targets = c("CDKN1A","MDM2","BAX","BBC3","PUMA","GADD45A","TIGAR","SESN1",
                   "DDB2","RRM2B","FDXR","TP53I3","ZMAT3","PHLDA3","GDF15"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    MYC = list(
      targets = c("LDHA","CDK4","NCL","NPM1","RPS6","RPL11","EIF4E","SRM","ODC1",
                   "NOP56","MRTO4","BYSL","PES1","POLR1B","UBF"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    STAT1 = list(
      targets = c("IRF1","GBP1","IFIT1","IFIT2","IFIT3","ISG15","MX1","OAS1",
                   "STAT2","TAP1","CXCL10","IDO1","IFITM1","BST2","IFI35"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    STAT3 = list(
      targets = c("BCL2L1","MCL1","CCND1","MYC","VEGFA","HIF1A","MMP9","SOCS3",
                   "BCL2","BIRC5","IL6","IL10","FOS","JUNB","PIM1"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    NRF2 = list(
      targets = c("NQO1","HMOX1","GCLM","GCLC","TXNRD1","AKR1C1","AKR1B10",
                   "GPX2","SRXN1","SLC7A11","FTH1","FTL","ME1","G6PD","TKT"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    HIF1A = list(
      targets = c("VEGFA","SLC2A1","LDHA","PDK1","PGK1","ENO1","PKM","BNIP3",
                   "CA9","LOX","EGLN3","ADM","EPO","HMOX1","SLC2A3"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    NFKB1 = list(
      targets = c("TNF","IL6","IL8","CCL2","ICAM1","VCAM1","BIRC3","BCL2A1",
                   "NFKBIA","TNFAIP3","CXCL1","CXCL2","IRF4","RELB","TRAF1"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    E2F1 = list(
      targets = c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","CDC6","CDT1",
                   "PCNA","RRM2","TK1","TYMS","CCNE1","CDC25A","DHFR"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    SOX2 = list(
      targets = c("POU5F1","NANOG","UTF1","FGF4","FOXD3","ZFP42","SALL4",
                   "DPPA4","TDGF1","LEFTY1","LEFTY2","GDF3","KLF4","LIN28A","PODXL"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    NANOG = list(
      targets = c("POU5F1","SOX2","ESRRB","SALL4","DPPA4","TCL1A","UTF1",
                   "ZFP42","TDGF1","GDF3","FOXD3","FGF4","LEFTY1","LEFTY2","NR5A2"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    POU5F1 = list(
      targets = c("SOX2","NANOG","UTF1","FGF4","DPPA4","SALL4","FOXD3",
                   "ZFP42","TDGF1","LEFTY1","GDF3","KLF4","ESRRB","CDX2","HAND1"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1)
    ),
    GATA1 = list(
      targets = c("HBB","HBA1","HBA2","ALAS2","SLC4A1","GYPA","KLF1","ANK1",
                   "EPB42","EPOR","TAL1","GATA2","NFE2","AHSP","HEMGN"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,-1,1,1,1)
    ),
    GATA3 = list(
      targets = c("IL4","IL5","IL13","IL10","CCR4","CXCR4","STAT6","MAF",
                   "BATF","IRF4","PTGDR2","IL1RL1","PPARG","TBX21","IFNG"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1)
    ),
    FOXP3 = list(
      targets = c("IL2RA","CTLA4","TNFRSF18","IKZF2","LRRC32","ENTPD1","NT5E",
                   "TIGIT","IL10","TGFB1","IL2","IFNG","IL4","IL17A","SATB1"),
      weights = c(1,1,1,1,1,1,1,1,1,1,-1,-1,-1,-1,-1)
    ),
    TBX21 = list(
      targets = c("IFNG","CXCR3","IL12RB2","IL18R1","CCL3","CCL4","GZMB","PRF1",
                   "EOMES","RUNX3","STAT4","ZEB2","BHLHE40","IL4","GATA3"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1)
    ),
    RORC = list(
      targets = c("IL17A","IL17F","IL22","IL23R","CCR6","IL1R1","AHR","BATF",
                   "IRF4","RORA","IL26","CCL20","STAT3","FOXP3","TBX21"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,-1,-1)
    ),
    BCL6 = list(
      targets = c("BCL2","CXCR5","PDCD1","ICOS","IL21","SH2D1A","CD200",
                   "ASCL2","TCF7","TOX2","MAF","BATF","PRDM1","IL2RA","CCR7"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,-1,-1,-1)
    ),
    IRF4 = list(
      targets = c("PRDM1","XBP1","IL10","CCL17","MMP9","SLAMF7","CD138",
                   "IRF8","BATF","SPI1","BCL6","CIITA","PAX5","EBF1","AICDA"),
      weights = c(1,1,1,1,1,1,1,-1,1,1,-1,-1,-1,-1,1)
    ),
    SPI1 = list(
      targets = c("CSF1R","ITGAM","CD68","LYZ","CTSS","TYROBP","FCER1G",
                   "C1QA","C1QB","C1QC","AIF1","MS4A7","FCGR1A","CD14","MAFB"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    ),
    CEBPA = list(
      targets = c("CSF3R","SLC22A4","LTF","MPO","ELANE","PRTN3","AZU1",
                   "CTSG","MS4A3","HP","CEBPE","GFI1","KLF5","ID1","CEBPD"),
      weights = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
    )
  )

  # Build data frame
  tf_list <- list()
  for (tf_name in names(regulons)) {
    reg <- regulons[[tf_name]]
    n_targets <- length(reg$targets)
    tf_list[[length(tf_list) + 1]] <- data.frame(
      tf = rep(tf_name, n_targets),
      target = reg$targets,
      weight = reg$weights,
      stringsAsFactors = FALSE
    )
  }
  tf_db <- do.call(rbind, tf_list)

  if (species == "mouse") {
    tf_db$tf <- .to_mouse_gene(tf_db$tf)
    tf_db$target <- .to_mouse_gene(tf_db$target)
  }

  attr(tf_db, "is_fallback") <- TRUE
  return(tf_db)
}

tfActivityUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Transcription Factor Activity Analysis"),
          p("Estimate TF activity from target gene expression using a built-in TF-target regulon database (DOROTHEA-like)."),

          fluidRow(
            column(3, selectInput(ns("tf_groupby"), "Cluster / group column", choices = NULL)),
            column(3, selectInput(ns("tf_species"), "Species",
              choices = c("Human" = "human", "Mouse" = "mouse"), selected = "human")),
            column(3, numericInput(ns("tf_min_targets"), "Min targets found",
              value = 5, min = 1, max = 20, step = 1)),
            column(3, selectInput(ns("tf_select"), "Select TFs (empty = all)",
              choices = NULL, multiple = TRUE))
          ),

          fluidRow(
            column(3, selectInput(ns("tf_score_method"), "TF Scoring Method",
              choices = c(
                "decoupleR ULM (Recommended)" = "ulm",
                "decoupleR Weighted Mean" = "wmean",
                "decoupleR Weighted Sum" = "wsum",
                "decoupleR VIPER" = "viper"
              ), selected = "ulm")),
            column(3, numericInput(ns("tf_top_n"), "Top TFs to display",
              value = 20, min = 5, max = 50, step = 5)),
            column(3, selectInput(ns("tf_scale"), "Scale heatmap",
              choices = c("Row (per TF)" = "row", "Column (per cluster)" = "column", "None" = "none"),
              selected = "row")),
            column(3, br(), actionButton(ns("run_tf"), "Run TF Activity Analysis",
              class = "btn-primary btn-lg", width = "100%"))
          )
        )
      )
    ),

    fluidRow(
      column(12,
        uiOutput(ns("tf_results_ui"))
      )
    )
  )
}

tfActivityServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    tf_results <- reactiveVal(NULL)
    tf_cell_scores <- reactiveVal(NULL)
    show_tf_results <- reactiveVal(FALSE)

    output$tf_results_ui <- renderUI({
      req(isTRUE(show_tf_results()))
      ns <- session$ns
      wellPanel(
        h4("TF Activity Results"),

        tabsetPanel(
          tabPanel("Activity Heatmap",
            fluidRow(
              column(10, plotlyOutput(ns("tf_heatmap"), height = "650px")),
              column(2,
                h5("Heatmap Options"),
                checkboxInput(ns("tf_hm_cluster_rows"), "Cluster TFs", TRUE),
                checkboxInput(ns("tf_hm_cluster_cols"), "Cluster groups", TRUE),
                br(),
                downloadButton(ns("download_tf_heatmap"), "Download Plot",
                  class = "btn-success", width = "100%")
              )
            )
          ),

          tabPanel("Violin Plots",
            fluidRow(
              column(3, selectInput(ns("tf_violin_select"), "Select TF",
                choices = NULL)),
              column(9, plotlyOutput(ns("tf_violin"), height = "500px"))
            ),
            br(),
            downloadButton(ns("download_tf_violin"), "Download Violin",
              class = "btn-success")
          ),

          tabPanel("UMAP Activity",
            fluidRow(
              column(3, selectInput(ns("tf_umap_select"), "Select TF",
                choices = NULL)),
              column(9, plotlyOutput(ns("tf_umap"), height = "600px"))
            ),
            br(),
            downloadButton(ns("download_tf_umap"), "Download UMAP",
              class = "btn-success")
          ),

          tabPanel("Results Table",
            DT::dataTableOutput(ns("tf_table")),
            br(),
            fluidRow(
              column(6,
                downloadButton(ns("download_tf_csv"), "Download All Results (CSV)",
                  class = "btn-success", width = "100%")),
              column(6,
                downloadButton(ns("download_tf_cell_csv"), "Download Per-Cell Scores (CSV)",
                  class = "btn-info", width = "100%"))
            )
          )
        )
      )
    })

    # Populate UI inputs
    observe({
      req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "tf_groupby", choices = sort(cat_cols))
    })

    observe({
      req(input$tf_species)
      tf_db <- .get_tf_regulon_db(input$tf_species)
      all_tfs <- sort(unique(tf_db$tf))
      updateSelectInput(session, "tf_select", choices = all_tfs, selected = character(0))
    })

    # Main TF activity analysis
    observeEvent(input$run_tf, {
      req(values$active_dataset, input$tf_groupby)

      if (!input$tf_groupby %in% colnames(values$active_dataset@meta.data)) {
        showNotification("Selected group column not found in metadata", type = "error")
        return()
      }

      withProgress(message = "Computing TF activity scores...", value = 0, {
        tryCatch({
          obj <- values$active_dataset
          group_col <- input$tf_groupby
          species <- input$tf_species
          min_targets <- input$tf_min_targets
          score_method <- input$tf_score_method
          selected_tfs <- input$tf_select

          incProgress(0.1, detail = "Loading TF regulon database...")
          tf_db <- .get_tf_regulon_db(species)

          # Prominent warning if using fallback (non-publication-grade) database
          if (isTRUE(attr(tf_db, "is_fallback"))) {
            showNotification(
              HTML(paste0(
                "<strong>WARNING: MINIMAL FALLBACK DATABASE IN USE</strong><br>",
                "Only ~20 TFs available. Results are NOT publication-grade.<br>",
                "Install decoupleR for comprehensive analysis:<br>",
                "<code>BiocManager::install('decoupleR')</code>"
              )),
              type = "warning", duration = 30
            )
          }

          # Get expression data
          assay_name <- DefaultAssay(obj)
          expr_mat <- safe_get_assay_data(obj, "data", assay = assay_name)
          all_genes <- rownames(expr_mat)
          groups <- as.character(obj@meta.data[[group_col]])
          clusters <- sort(unique(groups))

          # Filter TF database to genes present in dataset
          tf_db <- tf_db[tf_db$target %in% all_genes, ]

          # Count targets per TF and filter
          tf_target_counts <- table(tf_db$tf)
          valid_tfs <- names(tf_target_counts)[tf_target_counts >= min_targets]

          if (length(selected_tfs) > 0) {
            valid_tfs <- intersect(valid_tfs, selected_tfs)
          }

          if (length(valid_tfs) == 0) {
            showNotification("No TFs found with enough target genes in the dataset", type = "error")
            return()
          }

          tf_db <- tf_db[tf_db$tf %in% valid_tfs, ]

          incProgress(0.3, detail = "Computing TF activity...")

          target_genes <- unique(tf_db$target)

          # Subsample cells if matrix would be too large
          max_cells_tf <- 5000
          if (ncol(expr_mat) > max_cells_tf) {
            set.seed(42)
            # Stratified sampling: min(50, cluster_size) cells per cluster
            sampled_cells <- unlist(lapply(split(seq_along(groups), groups), function(idx) {
              n_sample <- min(50, length(idx))
              sample(idx, n_sample)
            }))
            if (length(sampled_cells) > max_cells_tf) {
              sampled_cells <- sample(sampled_cells, max_cells_tf)
            }
            target_mat <- as.matrix(expr_mat[target_genes, sampled_cells, drop = FALSE])
            groups <- groups[sampled_cells]
            clusters <- sort(unique(groups))
          } else {
            target_mat <- as.matrix(expr_mat[target_genes, , drop = FALSE])
          }

          tf_names <- sort(valid_tfs)

          incProgress(0.3, detail = "Computing per-cell and cluster activity...")

          # Cell-level activity (kept for UMAP overlay)
          # Try decoupleR's run_ulm / run_wmean / run_wsum first (preferred),
          # then fall back to VIPER, then to custom manual scoring.
          decoupler_cell_succeeded <- FALSE
          if (score_method %in% c("ulm", "wmean", "wsum") &&
              requireNamespace("decoupleR", quietly = TRUE)) {
            tryCatch({
              net_for_dc <- tf_db
              colnames(net_for_dc) <- c("source", "target", "mor")
              net_for_dc$mor[is.na(net_for_dc$mor)] <- 1

              dc_run_fn <- switch(score_method,
                ulm   = decoupleR::run_ulm,
                wmean = decoupleR::run_wmean,
                wsum  = decoupleR::run_wsum
              )

              dc_cell_res <- dc_run_fn(
                mat     = as.matrix(target_mat),
                net     = net_for_dc,
                .source = "source",
                .target = "target",
                .mor    = "mor",
                verbose = FALSE
              )

              # Pick the best statistic: normalized variant if available
              dc_stat <- switch(score_method,
                ulm   = "ulm",
                wmean = "norm_wmean",
                wsum  = "norm_wsum"
              )
              avail_stats <- unique(dc_cell_res$statistic)
              if (!dc_stat %in% avail_stats) dc_stat <- avail_stats[1]

              dc_filt <- dc_cell_res[dc_cell_res$statistic == dc_stat, ]

              dc_wide_df <- tidyr::pivot_wider(
                dc_filt[, c("source", "condition", "score")],
                names_from = "condition", values_from = "score", values_fill = 0
              )
              dc_mat <- as.matrix(dc_wide_df[, -1, drop = FALSE])
              rownames(dc_mat) <- dc_wide_df$source

              # Align to expected TFs and cells
              common_tfs <- intersect(tf_names, rownames(dc_mat))
              common_cells <- intersect(colnames(target_mat), colnames(dc_mat))
              if (length(common_tfs) > 0 && length(common_cells) > 0) {
                cell_activity <- dc_mat[common_tfs, common_cells, drop = FALSE]
                # Pad any missing TFs with zeros
                if (length(setdiff(tf_names, common_tfs)) > 0) {
                  missing_tfs <- setdiff(tf_names, common_tfs)
                  pad <- matrix(0, nrow = length(missing_tfs), ncol = ncol(cell_activity),
                                dimnames = list(missing_tfs, colnames(cell_activity)))
                  cell_activity <- rbind(cell_activity, pad)
                }
                cell_activity <- cell_activity[tf_names, , drop = FALSE]
                decoupler_cell_succeeded <- TRUE
                showNotification(
                  paste0("Using decoupleR::", score_method, " for TF scoring (",
                         length(common_tfs), " TFs)"),
                  type = "message", duration = 8
                )
              }
            }, error = function(e) {
              showNotification(
                paste0("decoupleR::", score_method, " failed, falling back to custom method: ", e$message),
                type = "warning", duration = 10
              )
            })
          }

          viper_succeeded <- FALSE
          if (!decoupler_cell_succeeded && score_method == "viper" &&
              requireNamespace("decoupleR", quietly = TRUE)) {
            # Use decoupleR's run_viper for proper VIPER scoring
            tryCatch({
              net_for_viper <- tf_db
              colnames(net_for_viper) <- c("source", "target", "mor")
              viper_res <- decoupleR::run_viper(
                mat = as.matrix(target_mat),
                net = net_for_viper,
                .source = "source",
                .target = "target",
                .mor = "mor",
                verbose = FALSE
              )
              # Extract activity matrix
              viper_wide <- viper_res %>%
                dplyr::filter(statistic == "viper") %>%
                dplyr::select(source, condition, score) %>%
                tidyr::pivot_wider(names_from = condition, values_from = score)
              viper_mat <- as.matrix(viper_wide[,-1])
              rownames(viper_mat) <- viper_wide$source
              cell_activity <- viper_mat[intersect(tf_names, rownames(viper_mat)), , drop=FALSE]
              viper_succeeded <- TRUE
              showNotification("Using decoupleR::VIPER for TF scoring", type = "message", duration = 8)
            }, error = function(e) {
              showNotification(paste("VIPER failed, falling back to weighted mean:", e$message), type="warning")
              score_method <<- "wmean"  # fallback
            })
          }

          if (!decoupler_cell_succeeded && !viper_succeeded) {
            # Manual scoring fallback: z-score per gene, then weighted mean/sum
            showNotification("Using custom weighted scoring (decoupleR not available)", type = "message", duration = 8)
            gene_means <- rowMeans(target_mat)
            gene_sds <- apply(target_mat, 1, sd)
            gene_sds[gene_sds == 0] <- 1
            z_mat <- (target_mat - gene_means) / gene_sds

            n_cells <- ncol(z_mat)
            cell_activity <- matrix(0, nrow = length(tf_names), ncol = n_cells,
                                    dimnames = list(tf_names, colnames(z_mat)))
            for (tf in tf_names) {
              tf_targets <- tf_db[tf_db$tf == tf, ]
              sub_z <- z_mat[tf_targets$target, , drop = FALSE]
              weighted_expr <- sweep(sub_z, 1, tf_targets$weight, `*`)
              if (score_method %in% c("wmean", "ulm")) {
                cell_activity[tf, ] <- colMeans(weighted_expr)
              } else {
                cell_activity[tf, ] <- colSums(weighted_expr)
              }
            }
          }

          incProgress(0.2, detail = "Aggregating per cluster...")

          # Cluster-level activity: when decoupleR/VIPER succeeded at cell level,
          # aggregate the cell-level scores by cluster mean. Otherwise, use the
          # custom method (per-cluster mean expression -> weighted TF score -> z-score).
          if (decoupler_cell_succeeded || viper_succeeded) {
            # Aggregate cell-level decoupleR/VIPER scores to cluster means
            cell_names <- colnames(cell_activity)
            cell_groups <- groups[match(cell_names, colnames(target_mat))]
            # Handle case where cell_activity has fewer cells
            if (any(is.na(cell_groups))) {
              valid_idx <- !is.na(cell_groups)
              cell_activity_for_agg <- cell_activity[, valid_idx, drop = FALSE]
              cell_groups <- cell_groups[valid_idx]
            } else {
              cell_activity_for_agg <- cell_activity
            }
            cluster_activity <- matrix(0, nrow = nrow(cell_activity_for_agg),
                                       ncol = length(clusters),
                                       dimnames = list(rownames(cell_activity_for_agg), clusters))
            for (cl in clusters) {
              idx <- which(cell_groups == cl)
              if (length(idx) > 1) {
                cluster_activity[, cl] <- rowMeans(cell_activity_for_agg[, idx, drop = FALSE])
              } else if (length(idx) == 1) {
                cluster_activity[, cl] <- cell_activity_for_agg[, idx]
              }
            }
          } else {
            # Fallback: custom cluster-level scoring from expression
            cluster_means <- matrix(0, nrow = length(target_genes), ncol = length(clusters),
                                    dimnames = list(target_genes, clusters))
            for (cl in clusters) {
              idx <- which(groups == cl)
              if (length(idx) > 1) {
                cluster_means[, cl] <- rowMeans(target_mat[, idx, drop = FALSE])
              } else {
                cluster_means[, cl] <- as.numeric(target_mat[, idx, drop = FALSE])
              }
            }

            cluster_activity <- matrix(0, nrow = length(tf_names), ncol = length(clusters),
                                       dimnames = list(tf_names, clusters))
            for (tf in tf_names) {
              tf_targets <- tf_db[tf_db$tf == tf, ]
              sub_expr <- cluster_means[tf_targets$target, , drop = FALSE]
              weighted_expr <- sweep(sub_expr, 1, tf_targets$weight, `*`)
              if (score_method %in% c("wmean", "ulm")) {
                cluster_activity[tf, ] <- colMeans(weighted_expr)
              } else {
                cluster_activity[tf, ] <- colSums(weighted_expr)
              }
            }

            # Z-score cluster-level activity row-wise (per TF across clusters)
            cluster_activity <- t(scale(t(cluster_activity)))
            cluster_activity[is.na(cluster_activity)] <- 0
          }

          # Build results table
          res_list <- list()
          for (tf in tf_names) {
            for (cl in clusters) {
              res_list[[length(res_list) + 1]] <- data.frame(
                TF = tf,
                Cluster = cl,
                Activity = round(cluster_activity[tf, cl], 4),
                Num_Targets = sum(tf_db$tf == tf),
                stringsAsFactors = FALSE
              )
            }
          }
          results_df <- do.call(rbind, res_list)

          tf_results(results_df)
          values$last_tf_results <- results_df
          tf_cell_scores(cell_activity)
          show_tf_results(TRUE)

          # Update violin / UMAP selectors
          updateSelectInput(session, "tf_violin_select", choices = tf_names, selected = tf_names[1])
          updateSelectInput(session, "tf_umap_select", choices = tf_names, selected = tf_names[1])

          incProgress(0.1, detail = "Done!")
          showNotification(paste0("TF activity computed for ", length(tf_names), " TFs across ",
                                  length(clusters), " clusters"), type = "message")
          gc()

        }, error = function(e) {
          showNotification(paste("TF activity error:", e$message), type = "error")
        })
      })
    })

    # --- Heatmap ---
    output$tf_heatmap <- renderPlotly({
      tryCatch({
        req(tf_results())
        df <- tf_results()

        # Pivot to matrix: TF x Cluster
        wide <- tidyr::pivot_wider(df, id_cols = TF, names_from = Cluster,
                                   values_from = Activity)
        mat <- as.matrix(wide[, -1])
        rownames(mat) <- wide$TF

        # Scale if requested
        scale_opt <- input$tf_scale
        if (scale_opt == "row") {
          row_sds <- apply(mat, 1, sd)
          row_sds[row_sds == 0] <- 1
          mat <- (mat - rowMeans(mat)) / row_sds
        } else if (scale_opt == "column") {
          col_sds <- apply(mat, 2, sd)
          col_sds[col_sds == 0] <- 1
          mat <- t((t(mat) - colMeans(mat)) / col_sds)
        }

        # Optionally cluster rows/columns
        top_n <- min(input$tf_top_n, nrow(mat))
        # Select top TFs by variance
        tf_var <- apply(mat, 1, var)
        top_tfs <- names(sort(tf_var, decreasing = TRUE))[seq_len(top_n)]
        mat <- mat[top_tfs, , drop = FALSE]

        if (input$tf_hm_cluster_rows && nrow(mat) > 2) {
          row_ord <- tryCatch(hclust(dist(mat))$order, error = function(e) seq_len(nrow(mat)))
          mat <- mat[row_ord, , drop = FALSE]
        }
        if (input$tf_hm_cluster_cols && ncol(mat) > 2) {
          col_ord <- tryCatch(hclust(dist(t(mat)))$order, error = function(e) seq_len(ncol(mat)))
          mat <- mat[, col_ord, drop = FALSE]
        }

        p <- plot_ly(x = colnames(mat), y = rownames(mat), z = mat,
                type = "heatmap", colors = colorRamp(c("navy", "white", "firebrick3")),
                hovertemplate = "TF: %{y}<br>Cluster: %{x}<br>Activity: %{z:.3f}<extra></extra>") %>%
          layout(title = "TF Activity Heatmap",
                 xaxis = list(title = "Cluster", tickangle = -45),
                 yaxis = list(title = "Transcription Factor"),
                 margin = list(l = 100, b = 80))
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "TF Activity", figure_type = "TF Heatmap",
                         label = "TF Activity Heatmap", source_module = "TF Activity")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })

    # --- Violin Plot ---
    output$tf_violin <- renderPlotly({
      tryCatch({
        req(tf_cell_scores(), input$tf_violin_select, values$active_dataset)
        cell_scores <- tf_cell_scores()
        tf_sel <- input$tf_violin_select
        req(tf_sel %in% rownames(cell_scores))

        obj <- values$active_dataset
        group_col <- input$tf_groupby
        groups <- as.character(obj@meta.data[[group_col]])
        scores <- cell_scores[tf_sel, ]

        plot_df <- data.frame(
          Group = groups,
          Activity = scores,
          stringsAsFactors = FALSE
        )

        clusters <- sort(unique(plot_df$Group))
        colors <- get_color_palette(length(clusters))

        p <- plot_ly()
        for (i in seq_along(clusters)) {
          cl <- clusters[i]
          vals <- plot_df$Activity[plot_df$Group == cl]
          p <- p %>% add_trace(y = vals, x = rep(cl, length(vals)),
                               type = "violin", name = cl,
                               box = list(visible = TRUE),
                               meanline = list(visible = TRUE),
                               fillcolor = colors[i],
                               line = list(color = colors[i]),
                               opacity = 0.7)
        }
        p %>% layout(title = paste("TF Activity:", tf_sel),
                     xaxis = list(title = "Cluster"),
                     yaxis = list(title = "Activity Score"),
                     showlegend = FALSE,
                     autosize = TRUE)
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })

    # --- UMAP colored by TF activity ---
    output$tf_umap <- renderPlotly({
      tryCatch({
        req(tf_cell_scores(), input$tf_umap_select, values$active_dataset)
        cell_scores <- tf_cell_scores()
        tf_sel <- input$tf_umap_select
        req(tf_sel %in% rownames(cell_scores))

        obj <- values$active_dataset
        red <- best_reduction(obj)
        req(!is.null(red))

        coords <- Embeddings(obj, reduction = red)[, 1:2]
        scores <- cell_scores[tf_sel, colnames(coords)]

        plot_df <- data.frame(
          Dim1 = coords[, 1],
          Dim2 = coords[, 2],
          Activity = scores,
          stringsAsFactors = FALSE
        )

        p <- plot_ly(plot_df, x = ~Dim1, y = ~Dim2, color = ~Activity,
                type = "scattergl", mode = "markers",
                marker = list(size = 3, opacity = 0.7),
                colors = c("lightgrey", "gold", "firebrick"),
                hovertemplate = paste0(
                  red, "_1: %{x:.2f}<br>",
                  red, "_2: %{y:.2f}<br>",
                  "Activity: %{marker.color:.3f}<extra></extra>")) %>%
          layout(title = paste("TF Activity:", tf_sel, "on", toupper(red)),
                 xaxis = list(title = paste0(toupper(red), "_1")),
                 yaxis = list(title = paste0(toupper(red), "_2")))
        if (exists("capture_umap_snapshot", mode = "function"))
          capture_umap_snapshot(values, p, "TF Activity", coloring = tf_sel, method = red)
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })

    # --- Results Table ---
    output$tf_table <- DT::renderDT({
      req(tf_results())
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tf_results(), description = "TF Activity",
                            source_module = "Downstream")
      }
      DT::datatable(tf_results(), filter = "top", rownames = FALSE,
                    options = list(pageLength = 20, scrollX = TRUE,
                                   order = list(list(2, "desc"))))
    })

    # --- Downloads ---
    output$download_tf_csv <- downloadHandler(
      filename = function() paste0("tf_activity_", Sys.Date(), ".csv"),
      content = function(file) {
        req(tf_results())
        write.csv(tf_results(), file, row.names = FALSE)
      }
    )

    output$download_tf_cell_csv <- downloadHandler(
      filename = function() paste0("tf_cell_scores_", Sys.Date(), ".csv"),
      content = function(file) {
        req(tf_cell_scores())
        out <- as.data.frame(t(tf_cell_scores()))
        out$Cell <- rownames(out)
        out <- out[, c("Cell", setdiff(colnames(out), "Cell"))]
        write.csv(out, file, row.names = FALSE)
      }
    )

    output$download_tf_heatmap <- downloadHandler(
      filename = function() paste0("tf_heatmap_", Sys.Date(), ".html"),
      content = function(file) {
        p <- tryCatch(plotly_empty(), error = function(e) plotly_empty())
        htmlwidgets::saveWidget(p, file)
      }
    )

    output$download_tf_violin <- downloadHandler(
      filename = function() paste0("tf_violin_", Sys.Date(), ".html"),
      content = function(file) {
        p <- tryCatch(plotly_empty(), error = function(e) plotly_empty())
        htmlwidgets::saveWidget(p, file)
      }
    )

    output$download_tf_umap <- downloadHandler(
      filename = function() paste0("tf_umap_", Sys.Date(), ".html"),
      content = function(file) {
        p <- tryCatch(plotly_empty(), error = function(e) plotly_empty())
        htmlwidgets::saveWidget(p, file)
      }
    )
  })
}


# ---------------------------------------------------------------------------
# decoupleR-based TF activity scoring (preferred when decoupleR is available)
# Uses run_ulm / run_wmean / run_wsum directly on the expression matrix.
# Returns the same format as .run_tf_activity_auto():
#   list(results = data.frame(TF, Cluster, Activity, Num_Targets),
#        cluster_activity = matrix[TF x Cluster])
# ---------------------------------------------------------------------------
run_tf_activity_decoupler <- function(seurat_obj, group_col, species = "human",
                                      method = c("ulm", "wmean", "wsum"),
                                      min_targets = 5, max_cells = 5000) {
  method <- match.arg(method)

  if (!requireNamespace("decoupleR", quietly = TRUE)) {
    stop("decoupleR package is not installed.")
  }

  # 1. Load the CollecTRI regulon network via decoupleR
  organism <- if (species == "mouse") "mouse" else "human"
  net <- tryCatch(
    decoupleR::get_collectri(organism = organism, split_complexes = FALSE),
    error = function(e) {
      stop("Failed to load CollecTRI network: ", e$message)
    }
  )
  if (!is.data.frame(net) || nrow(net) < 100) {
    stop("CollecTRI network too small or invalid (", nrow(net), " rows).")
  }

  # Ensure standard column names expected by decoupleR: source, target, mor
  if (!"source" %in% colnames(net) && "tf" %in% colnames(net)) {
    colnames(net)[colnames(net) == "tf"] <- "source"
  }
  if (!"mor" %in% colnames(net) && "weight" %in% colnames(net)) {
    colnames(net)[colnames(net) == "weight"] <- "mor"
  }
  net$mor[is.na(net$mor)] <- 1

  # 2. Extract expression matrix from Seurat object
  assay_name <- DefaultAssay(seurat_obj)
  expr_mat <- safe_get_assay_data(seurat_obj, "data", assay = assay_name)
  all_genes <- rownames(expr_mat)
  groups <- as.character(seurat_obj@meta.data[[group_col]])
  clusters <- sort(unique(groups))

  # Filter network to genes present in dataset
  net <- net[net$target %in% all_genes, ]
  tf_target_counts <- table(net$source)
  valid_tfs <- names(tf_target_counts)[tf_target_counts >= min_targets]
  if (length(valid_tfs) == 0) {
    stop("No TFs found with >= ", min_targets, " target genes in the dataset.")
  }
  net <- net[net$source %in% valid_tfs, ]

  # 3. Subsample cells if dataset is very large (for performance)
  if (ncol(expr_mat) > max_cells) {
    set.seed(42)
    sampled_cells <- unlist(lapply(split(seq_along(groups), groups), function(idx) {
      n_sample <- min(50, length(idx))
      sample(idx, n_sample)
    }))
    if (length(sampled_cells) > max_cells) {
      sampled_cells <- sample(sampled_cells, max_cells)
    }
    expr_dense <- as.matrix(expr_mat[, sampled_cells, drop = FALSE])
    groups_sub <- groups[sampled_cells]
  } else {
    expr_dense <- as.matrix(expr_mat)
    groups_sub <- groups
  }
  clusters <- sort(unique(groups_sub))

  # 4. Run decoupleR method on the full expression matrix (cell-level scores)
  run_fn <- switch(method,
    ulm   = decoupleR::run_ulm,
    wmean = decoupleR::run_wmean,
    wsum  = decoupleR::run_wsum
  )

  dc_res <- tryCatch(
    run_fn(
      mat     = expr_dense,
      net     = net,
      .source = "source",
      .target = "target",
      .mor    = "mor",
      verbose = FALSE
    ),
    error = function(e) {
      stop("decoupleR::", method, " failed: ", e$message)
    }
  )

  # 5. Extract the primary statistic from decoupleR results
  #    run_ulm returns statistic "ulm", run_wmean returns "norm_wmean",
  #    run_wsum returns "norm_wsum". Pick the normalized variant if present.
  stat_name <- switch(method,
    ulm   = "ulm",
    wmean = "norm_wmean",
    wsum  = "norm_wsum"
  )
  # Fallback: if the expected statistic is not present, use the first available
  available_stats <- unique(dc_res$statistic)
  if (!stat_name %in% available_stats) {
    stat_name <- available_stats[1]
  }

  dc_filtered <- dc_res[dc_res$statistic == stat_name, ]

  # 6. Aggregate cell-level scores to cluster-level (mean per cluster)
  #    dc_filtered has columns: source (TF), condition (cell barcode), score, p_value
  tf_names <- sort(intersect(valid_tfs, unique(dc_filtered$source)))
  if (length(tf_names) == 0) {
    stop("decoupleR returned no TF scores.")
  }

  # Build a TF x cell matrix from decoupleR output
  dc_wide <- tryCatch({
    wide_df <- dc_filtered[, c("source", "condition", "score")]
    wide_df <- tidyr::pivot_wider(wide_df, names_from = "condition",
                                  values_from = "score", values_fill = 0)
    mat <- as.matrix(wide_df[, -1, drop = FALSE])
    rownames(mat) <- wide_df$source
    mat
  }, error = function(e) {
    stop("Failed to reshape decoupleR results: ", e$message)
  })

  # Ensure we only use cells that exist in our (possibly subsampled) groups
  common_cells <- intersect(colnames(dc_wide), colnames(expr_dense))
  dc_wide <- dc_wide[, common_cells, drop = FALSE]
  groups_matched <- groups_sub[match(common_cells, colnames(expr_dense))]

  # Compute cluster-level mean activity
  cluster_activity <- matrix(0, nrow = length(tf_names), ncol = length(clusters),
                             dimnames = list(tf_names, clusters))
  for (cl in clusters) {
    idx <- which(groups_matched == cl)
    if (length(idx) > 1) {
      cluster_activity[, cl] <- rowMeans(dc_wide[tf_names, idx, drop = FALSE])
    } else if (length(idx) == 1) {
      cluster_activity[, cl] <- dc_wide[tf_names, idx]
    }
  }

  # 7. Extract p-values per TF (median p-value across cells, if available)
  has_pval <- "p_value" %in% colnames(dc_filtered)

  # 8. Build results data.frame in the standard format
  tf_db_counts <- table(net$source)
  res_list <- list()
  for (tf in tf_names) {
    for (cl in clusters) {
      pval <- NA_real_
      if (has_pval) {
        tf_cells <- common_cells[groups_matched == cl]
        pvals <- dc_filtered$p_value[dc_filtered$source == tf &
                                      dc_filtered$condition %in% tf_cells]
        if (length(pvals) > 0) pval <- round(median(pvals, na.rm = TRUE), 6)
      }
      res_list[[length(res_list) + 1]] <- data.frame(
        TF = tf,
        Cluster = cl,
        Activity = round(cluster_activity[tf, cl], 4),
        P_value = pval,
        Num_Targets = as.integer(tf_db_counts[tf]),
        stringsAsFactors = FALSE
      )
    }
  }
  results_df <- do.call(rbind, res_list)

  message("[TF] decoupleR::", method, " completed for ", length(tf_names),
          " TFs across ", length(clusters), " clusters (", organism, ")")
  list(results = results_df, cluster_activity = cluster_activity)
}

# Standalone TF activity: returns list(results=df, cluster_activity=matrix)
# Tries decoupleR first (ulm/wmean/wsum), then falls back to custom method.
.run_tf_activity_auto <- function(obj, group_col, species = "human",
                                   min_targets = 5, score_method = "wmean",
                                   max_cells = 5000) {

  # --- Attempt 1: decoupleR-based scoring (preferred) ---
  if (requireNamespace("decoupleR", quietly = TRUE)) {
    # Map legacy score_method names to decoupleR method names
    dc_method <- switch(score_method,
      ulm   = "ulm",
      wmean = "wmean",
      wsum  = "wsum",
      viper = "ulm",   # viper not directly available via this path; use ulm
      "wmean"           # default
    )
    dc_result <- tryCatch({
      res <- run_tf_activity_decoupler(
        seurat_obj  = obj,
        group_col   = group_col,
        species     = species,
        method      = dc_method,
        min_targets = min_targets,
        max_cells   = max_cells
      )
      message("[TF] Using decoupleR::", dc_method, " for TF activity scoring.")
      res
    }, error = function(e) {
      message("[TF] decoupleR scoring failed: ", e$message, ". Falling back to custom method.")
      NULL
    })
    if (!is.null(dc_result)) return(dc_result)
  } else {
    message("[TF] decoupleR not installed. Using custom weighted scoring method.")
  }

  # --- Attempt 2: Custom weighted scoring (fallback) ---
  tryCatch({
    tf_db <- .get_tf_regulon_db(species)
    assay_name <- DefaultAssay(obj)
    expr_mat <- safe_get_assay_data(obj, "data", assay = assay_name)
    all_genes <- rownames(expr_mat)
    groups <- as.character(obj@meta.data[[group_col]])
    clusters <- sort(unique(groups))

    tf_db <- tf_db[tf_db$target %in% all_genes, ]
    tf_target_counts <- table(tf_db$tf)
    valid_tfs <- names(tf_target_counts)[tf_target_counts >= min_targets]
    if (length(valid_tfs) == 0) return(NULL)
    tf_db <- tf_db[tf_db$tf %in% valid_tfs, ]

    target_genes <- unique(tf_db$target)
    if (ncol(expr_mat) > max_cells) {
      set.seed(42)
      # Stratified sampling: min(50, cluster_size) cells per cluster
      sampled_cells <- unlist(lapply(split(seq_along(groups), groups), function(idx) {
        n_sample <- min(50, length(idx))
        sample(idx, n_sample)
      }))
      if (length(sampled_cells) > max_cells) {
        sampled_cells <- sample(sampled_cells, max_cells)
      }
      target_mat <- as.matrix(expr_mat[target_genes, sampled_cells, drop = FALSE])
      groups <- groups[sampled_cells]
      clusters <- sort(unique(groups))
    } else {
      target_mat <- as.matrix(expr_mat[target_genes, , drop = FALSE])
    }

    # Compute per-cluster mean expression first, then weighted TF activity,
    # then z-score at the cluster level (per TF across clusters).
    # This avoids the dataset-composition confound of z-scoring per-cell first
    # (adding/removing a cluster would change every gene's z-score).
    tf_names <- sort(valid_tfs)

    # Step 1: Compute per-cluster mean expression for target genes
    cluster_means <- matrix(0, nrow = length(target_genes), ncol = length(clusters),
                            dimnames = list(target_genes, clusters))
    for (cl in clusters) {
      idx <- which(groups == cl)
      if (length(idx) > 1) {
        cluster_means[, cl] <- rowMeans(target_mat[, idx, drop = FALSE])
      } else {
        cluster_means[, cl] <- as.numeric(target_mat[, idx, drop = FALSE])
      }
    }

    # Step 2: Compute weighted TF activity per cluster (on raw mean expression)
    cluster_activity <- matrix(0, nrow = length(tf_names), ncol = length(clusters),
                               dimnames = list(tf_names, clusters))
    for (tf in tf_names) {
      tf_targets <- tf_db[tf_db$tf == tf, ]
      sub_expr <- cluster_means[tf_targets$target, , drop = FALSE]
      weighted_expr <- tf_targets$weight * sub_expr
      if (score_method == "wmean") {
        cluster_activity[tf, ] <- colMeans(weighted_expr)
      } else {
        cluster_activity[tf, ] <- colSums(weighted_expr)
      }
    }

    # Step 3: Z-score cluster-level activity row-wise (per TF across clusters)
    # This is how decoupleR::run_wmean is designed to be interpreted
    cluster_activity <- t(scale(t(cluster_activity)))
    cluster_activity[is.na(cluster_activity)] <- 0

    res_list <- list()
    for (tf in tf_names) {
      for (cl in clusters) {
        res_list[[length(res_list) + 1]] <- data.frame(
          TF = tf, Cluster = cl,
          Activity = round(cluster_activity[tf, cl], 4),
          Num_Targets = sum(tf_db$tf == tf),
          stringsAsFactors = FALSE)
      }
    }
    results_df <- do.call(rbind, res_list)
    message("[TF] Custom weighted scoring completed for ", length(tf_names), " TFs.")
    list(results = results_df, cluster_activity = cluster_activity)
  }, error = function(e) {
    message("Auto TF activity error: ", e$message)
    NULL
  })
}

# Standalone GRN: returns list(graph, nodes, edges, modules, layout)
.run_grn_auto <- function(obj, group_col, n_genes = 200, cor_method = "spearman",
                           cor_thresh = 0.3, min_degree = 2, community = "louvain") {
  tryCatch({
    if (!requireNamespace("igraph", quietly = TRUE)) return(NULL)
    assay_name <- DefaultAssay(obj)
    var_features <- tryCatch({
      vf <- VariableFeatures(obj)
      if (length(vf) == 0) {
        obj_tmp <- FindVariableFeatures(obj, nfeatures = n_genes)
        VariableFeatures(obj_tmp)
      } else vf
    }, error = function(e) {
      set.seed(42)
      expr_mat <- safe_get_assay_data(obj, "data", assay = assay_name)
      gene_var <- apply(as.matrix(expr_mat[sample(rownames(expr_mat),
                        min(nrow(expr_mat), 5000)), ]), 1, var)
      names(sort(gene_var, decreasing = TRUE))[seq_len(min(n_genes, length(gene_var)))]
    })

    sel_genes <- head(var_features, n_genes)
    if (length(sel_genes) < 10) return(NULL)

    expr_mat <- safe_get_assay_data(obj, "data", assay = assay_name)
    sel_genes <- sel_genes[sel_genes %in% rownames(expr_mat)]

    max_cells_cor <- 5000
    if (ncol(expr_mat) > max_cells_cor) {
      set.seed(42)
      cell_idx <- sample(ncol(expr_mat), max_cells_cor)
      sub_mat <- as.matrix(t(expr_mat[sel_genes, cell_idx]))
    } else {
      sub_mat <- as.matrix(t(expr_mat[sel_genes, ]))
    }

    cor_mat <- cor(sub_mat, method = cor_method, use = "pairwise.complete.obs")
    cor_mat[is.na(cor_mat)] <- 0
    diag(cor_mat) <- 0
    adj_mat <- (abs(cor_mat) >= cor_thresh) * 1

    g <- igraph::graph_from_adjacency_matrix(adj_mat, mode = "undirected",
                                              diag = FALSE, weighted = NULL)
    edge_list <- igraph::as_edgelist(g)
    if (nrow(edge_list) > 0) {
      edge_weights <- sapply(seq_len(nrow(edge_list)), function(i)
        cor_mat[edge_list[i, 1], edge_list[i, 2]])
      igraph::E(g)$weight <- edge_weights
      igraph::E(g)$abs_weight <- abs(edge_weights)
    }

    degrees <- igraph::degree(g)
    keep_nodes <- names(degrees)[degrees >= min_degree]
    if (length(keep_nodes) < 5) return(NULL)
    g <- igraph::induced_subgraph(g, keep_nodes)

    # Use absolute weights for community detection (walktrap/louvain require non-negative weights)
    if (!is.null(igraph::E(g)$weight) && any(igraph::E(g)$weight < 0)) {
      igraph::E(g)$weight <- abs(igraph::E(g)$weight)
    }

    communities_result <- tryCatch({
      switch(community,
        "louvain" = igraph::cluster_louvain(g),
        "walktrap" = igraph::cluster_walktrap(g),
        "fast_greedy" = igraph::cluster_fast_greedy(g),
        "label_prop" = igraph::cluster_label_prop(g),
        igraph::cluster_louvain(g))
    }, error = function(e) igraph::cluster_louvain(g))

    membership <- igraph::membership(communities_result)
    node_degree <- igraph::degree(g)
    node_betweenness <- igraph::betweenness(g, normalized = TRUE)
    node_closeness <- tryCatch(igraph::closeness(g, normalized = TRUE),
                               error = function(e) rep(0, length(igraph::V(g))))

    node_df <- data.frame(
      Gene = names(igraph::V(g)), Degree = node_degree,
      Betweenness = round(node_betweenness, 4),
      Closeness = round(node_closeness, 4),
      Module = as.integer(membership), stringsAsFactors = FALSE)
    node_df <- node_df[order(-node_df$Degree), ]

    el <- igraph::as_edgelist(g)
    edge_df <- data.frame(
      Gene1 = el[, 1], Gene2 = el[, 2],
      Correlation = round(igraph::E(g)$weight, 4),
      Abs_Correlation = round(igraph::E(g)$abs_weight, 4),
      stringsAsFactors = FALSE)

    module_df <- node_df %>%
      group_by(Module) %>%
      summarise(N_Genes = n(), Hub_Gene = Gene[which.max(Degree)],
                Mean_Degree = round(mean(Degree), 2),
                Top_Genes = paste(head(Gene[order(-Degree)], 5), collapse = ", "),
                .groups = "drop")

    lay <- tryCatch(igraph::layout_with_fr(g),
                    error = function(e) igraph::layout_randomly(g))
    rownames(lay) <- names(igraph::V(g))

    list(graph = g, nodes = node_df, edges = edge_df, modules = module_df, layout = lay)
  }, error = function(e) {
    message("Auto GRN error: ", e$message)
    NULL
  })
}

# TF activity heatmap (standalone)
.plot_tf_heatmap_auto <- function(tf_results, top_n = 20, scale_opt = "row") {
  tryCatch({
    if (is.null(tf_results) || !is.data.frame(tf_results) || nrow(tf_results) == 0)
      return(plotly_empty() %>% layout(title = list(text = "No TF activity results", font = list(size = 14))))

    wide <- tidyr::pivot_wider(tf_results, id_cols = TF, names_from = Cluster,
                               values_from = Activity)
    mat <- as.matrix(wide[, -1])
    rownames(mat) <- wide$TF
    if (scale_opt == "row") {
      row_sds <- apply(mat, 1, sd); row_sds[row_sds == 0] <- 1
      mat <- (mat - rowMeans(mat)) / row_sds
    }
    top_n <- min(top_n, nrow(mat))
    tf_var <- apply(mat, 1, var)
    top_tfs <- names(sort(tf_var, decreasing = TRUE))[seq_len(top_n)]
    mat <- mat[top_tfs, , drop = FALSE]
    if (nrow(mat) > 2) {
      ro <- tryCatch(hclust(dist(mat))$order, error = function(e) seq_len(nrow(mat)))
      mat <- mat[ro, , drop = FALSE]
    }
    if (ncol(mat) > 2) {
      co <- tryCatch(hclust(dist(t(mat)))$order, error = function(e) seq_len(ncol(mat)))
      mat <- mat[, co, , drop = FALSE]
    }
    max_abs <- max(abs(mat), na.rm = TRUE)

    plot_ly(x = colnames(mat), y = rownames(mat), z = mat, type = "heatmap",
            zmin = -max_abs, zmax = max_abs,
            colorscale = list(c(0,"#2166ac"), c(0.2,"#67a9cf"), c(0.4,"#d1e5f0"),
                              c(0.5,"#f7f7f7"), c(0.6,"#fddbc7"),
                              c(0.8,"#ef8a62"), c(1,"#b2182b")),
            colorbar = list(title = list(text = "Activity\nZ-score", font = list(size = 10)),
                            len = 0.55, thickness = 12),
            hovertemplate = paste0("<b>%{y}</b><br>Cluster: %{x}<br>",
                                  "Activity Z-score: %{z:.3f}<extra></extra>")) %>%
      layout(
        title = list(
          text = paste0("Transcription Factor Activity<br>",
                        "<span style='font-size:11px;color:#777'>Top ", top_n,
                        " most variable TFs (row-scaled z-scores)</span>"),
          font = list(size = 14, family = "Arial")),
        xaxis = list(title = "Cell Type / Cluster", tickangle = -45,
                     tickfont = list(size = 11), titlefont = list(size = 12)),
        yaxis = list(title = "", tickfont = list(size = 10)),
        margin = list(l = 100, b = 100, t = 55, r = 40),
        paper_bgcolor = "white", plot_bgcolor = "white")
  }, error = function(e) {
    plotly_empty() %>% layout(title = list(text = paste("TF error:", e$message), font = list(size = 12)))
  })
}

# GRN network plot (standalone)
.plot_grn_network_auto <- function(grn_result, max_edges = 500, show_labels = TRUE) {
  tryCatch({
    if (is.null(grn_result))
      return(plotly_empty() %>% layout(title = list(text = "No GRN results", font = list(size = 14))))

    g <- grn_result$graph; node_df <- grn_result$nodes; lay <- grn_result$layout
    node_x <- lay[, 1]; node_y <- lay[, 2]
    names(node_x) <- names(igraph::V(g))
    names(node_y) <- names(igraph::V(g))

    modules <- node_df$Module[match(names(igraph::V(g)), node_df$Gene)]
    unique_modules <- sort(unique(modules))
    n_modules <- length(unique_modules)
    mod_colors <- get_color_palette(max(n_modules, 3))

    degrees <- igraph::degree(g)
    size_range <- c(5, 16)
    if (max(degrees) > min(degrees)) {
      node_sizes <- size_range[1] + (degrees - min(degrees)) /
        (max(degrees) - min(degrees)) * (size_range[2] - size_range[1])
    } else {
      node_sizes <- rep(mean(size_range), length(degrees))
    }

    # Only label top 30% hub genes to reduce clutter
    node_names <- names(igraph::V(g))
    deg_thresh <- quantile(degrees, 0.7)
    label_text <- ifelse(degrees >= deg_thresh & show_labels, node_names, "")

    p <- plot_ly()
    el <- igraph::as_edgelist(g)
    edge_weights <- igraph::E(g)$abs_weight
    if (nrow(el) > max_edges) {
      top_idx <- order(edge_weights, decreasing = TRUE)[seq_len(max_edges)]
      el <- el[top_idx, , drop = FALSE]
    }
    if (nrow(el) > 0) {
      edge_x <- c(); edge_y <- c()
      for (i in seq_len(nrow(el))) {
        edge_x <- c(edge_x, node_x[el[i, 1]], node_x[el[i, 2]], NA)
        edge_y <- c(edge_y, node_y[el[i, 1]], node_y[el[i, 2]], NA)
      }
      p <- p %>% add_trace(x = edge_x, y = edge_y, type = "scattergl", mode = "lines",
        line = list(width = 0.4, color = "rgba(140,140,140,0.18)"),
        hoverinfo = "none", showlegend = FALSE)
    }

    hover_text <- paste0(
      "<b>", node_names, "</b>",
      "<br>Module: ", modules,
      "<br>Degree: ", degrees[node_names],
      "<br>Betweenness: ", round(node_df$Betweenness[match(node_names, node_df$Gene)], 3))

    # Plot nodes grouped by module for legend
    for (m in unique_modules) {
      mask <- modules == m
      if (!any(mask)) next
      hub_gene <- node_df$Gene[node_df$Module == m][which.max(node_df$Degree[node_df$Module == m])]
      mod_label <- paste0("Module ", m, " (hub: ", hub_gene, ")")
      p <- p %>% add_trace(
        x = node_x[node_names[mask]], y = node_y[node_names[mask]],
        type = "scattergl", mode = "markers+text",
        marker = list(size = node_sizes[mask],
                      color = mod_colors[match(m, unique_modules)],
                      line = list(width = 0.8, color = "white"), opacity = 0.88),
        text = label_text[mask], textposition = "top center",
        textfont = list(size = 9, color = "#222", family = "Arial"),
        hovertext = hover_text[mask], hoverinfo = "text",
        name = mod_label, legendgroup = paste("M", m))
    }

    n_v <- igraph::vcount(g); n_e <- igraph::ecount(g)
    # Module summary annotation
    mod_summary <- paste(sapply(head(sort(unique_modules), 5), function(m) {
      hub <- node_df$Gene[node_df$Module == m][which.max(node_df$Degree[node_df$Module == m])]
      n <- sum(modules == m)
      paste0("M", m, ": ", hub, " (", n, " genes)")
    }), collapse = " | ")

    p %>% layout(
      title = list(
        text = paste0("Gene Regulatory Network<br>",
                      "<span style='font-size:11px;color:#777'>",
                      n_v, " genes \u2022 ", n_e, " edges \u2022 ",
                      n_modules, " modules</span>"),
        font = list(size = 14, family = "Arial")),
      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
      legend = list(orientation = "h", x = 0.5, xanchor = "center",
                    y = -0.02, font = list(size = 10)),
      plot_bgcolor = "#fafafa", paper_bgcolor = "white",
      margin = list(t = 55, b = 50, l = 15, r = 15), autosize = TRUE)
  }, error = function(e) {
    plotly_empty() %>% layout(title = list(text = paste("GRN error:", e$message), font = list(size = 12)))
  })
}
