# Downstream CCC — Cell-Cell Communication
# modules/downstream/downstream_ccc.R

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE)) message("igraph recommended for network plots")
})

# ---------------------------------------------------------------------------
# L-R Database Fallback Chain: LIANA -> OmnipathR -> Built-in
# ---------------------------------------------------------------------------

# Primary: fetch L-R pairs from LIANA's curated multi-database resource
.get_lr_database_liana <- function(species = "human") {

  # ---- Attempt 1: LIANA select_resource() ----
  if (requireNamespace("liana", quietly = TRUE)) {
    lr <- tryCatch({
      message("[CCC] Fetching L-R database from LIANA...")
      res <- liana::select_resource("Consensus")
      if (is.null(res) || nrow(res) == 0) stop("Empty LIANA resource")

      # Standardize column names — LIANA uses source_genesymbol / target_genesymbol
      col_map <- list(
        ligand  = intersect(c("source_genesymbol", "ligand", "ligand_complex"), colnames(res))[1],
        receptor = intersect(c("target_genesymbol", "receptor", "receptor_complex"), colnames(res))[1]
      )
      if (is.na(col_map$ligand) || is.na(col_map$receptor)) stop("Unexpected LIANA column names")

      out <- data.frame(
        ligand   = as.character(res[[col_map$ligand]]),
        receptor = as.character(res[[col_map$receptor]]),
        pathway  = if ("category" %in% colnames(res)) as.character(res$category) else "LIANA",
        stringsAsFactors = FALSE
      )
      out <- out[!is.na(out$ligand) & !is.na(out$receptor) &
                   nchar(out$ligand) > 0 & nchar(out$receptor) > 0, ]
      out <- unique(out)

      if (species == "mouse") {
        out$ligand   <- .to_mouse_gene(out$ligand)
        out$receptor <- .to_mouse_gene(out$receptor)
      }

      message("[CCC] Loaded ", nrow(out), " L-R pairs from LIANA Consensus resource")
      out
    }, error = function(e) {
      message("[CCC] LIANA select_resource() failed: ", e$message)
      NULL
    })
    if (!is.null(lr)) return(lr)
  }

  # ---- Attempt 2: OmnipathR intercell network ----
  if (requireNamespace("OmnipathR", quietly = TRUE)) {
    lr <- tryCatch({
      message("[CCC] Fetching L-R pairs from OmnipathR...")
      net <- OmnipathR::import_intercell_network(
        transmitter_param = list(categories = "ligand"),
        receiver_param    = list(categories = "receptor")
      )
      if (is.null(net) || nrow(net) == 0) stop("Empty OmnipathR network")

      # OmnipathR returns source_genesymbol / target_genesymbol
      src_col <- intersect(c("source_genesymbol", "genesymbol_intercell_source"), colnames(net))[1]
      tgt_col <- intersect(c("target_genesymbol", "genesymbol_intercell_target"), colnames(net))[1]
      if (is.na(src_col) || is.na(tgt_col)) stop("Unexpected OmnipathR column names")

      out <- data.frame(
        ligand   = as.character(net[[src_col]]),
        receptor = as.character(net[[tgt_col]]),
        pathway  = "OmniPath",
        stringsAsFactors = FALSE
      )
      out <- out[!is.na(out$ligand) & !is.na(out$receptor) &
                   nchar(out$ligand) > 0 & nchar(out$receptor) > 0, ]
      out <- unique(out)

      if (species == "mouse") {
        out$ligand   <- .to_mouse_gene(out$ligand)
        out$receptor <- .to_mouse_gene(out$receptor)
      }

      message("[CCC] Loaded ", nrow(out), " L-R pairs from OmnipathR")
      out
    }, error = function(e) {
      message("[CCC] OmnipathR import failed: ", e$message)
      NULL
    })
    if (!is.null(lr)) return(lr)
  }

  # ---- Attempt 3: Built-in minimal database ----
  message("[CCC] Using minimal built-in L-R database (~85 pairs). ",
          "Install 'liana' or 'OmnipathR' for comprehensive analysis.")
  .get_lr_database(species)
}

# Built-in MINIMAL ligand-receptor database for human and mouse
# This is a curated fallback with ~85 key pairs covering major signaling pathways.
# For publication-grade analysis, install LIANA (>1,000 curated pairs) or OmnipathR.
.get_lr_database <- function(species = "human") {
  # Core human L-R pairs (ligand, receptor, pathway)
  human_lr <- data.frame(
    ligand = c(
      # Immune checkpoints
      "CD274", "CD274", "PDCD1LG2", "CD80", "CD86", "LGALS9", "PVR",
      # Cytokines
      "IL2", "IL2", "IL6", "IL6", "IL10", "IL10", "TNF", "TNF",
      "IFNG", "IFNG", "TGFB1", "TGFB1", "TGFB2", "TGFB2",
      "IL1B", "IL1B", "IL4", "IL4", "IL12A", "IL12B", "IL15", "IL17A",
      "IL21", "IL33", "IL18", "CSF1", "CSF2",
      # Chemokines
      "CXCL12", "CXCL12", "CCL2", "CCL5", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL13",
      # Growth factors
      "VEGFA", "VEGFA", "EGF", "AREG", "FGF2", "FGF2", "PDGFA", "PDGFB",
      "HGF", "BMP2", "BMP4", "WNT3A", "WNT5A",
      # Notch
      "DLL1", "DLL4", "JAG1", "JAG1",
      # Cell adhesion
      "ICAM1", "VCAM1", "CDH1", "CDH2",
      # TNF superfamily
      "TNFSF10", "FASLG", "CD40LG", "TNFSF4", "TNFSF9",
      # Costimulatory
      "CD28", "ICOS",
      # MHC / Antigen presentation
      "HLA-A", "HLA-B", "HLA-DRA", "B2M",
      # Semaphorins / Plexins
      "SEMA3A", "SEMA4D",
      # Galectins
      "LGALS1", "LGALS3",
      # Complement
      "C3", "C5",
      # Hedgehog
      "SHH", "IHH",
      # Ephrins
      "EFNA1", "EFNB1", "EFNB2"
    ),
    receptor = c(
      # Immune checkpoints
      "PDCD1", "CD80", "PDCD1", "CTLA4", "CTLA4", "HAVCR2", "TIGIT",
      # Cytokines
      "IL2RA", "IL2RB", "IL6R", "IL6ST", "IL10RA", "IL10RB", "TNFRSF1A", "TNFRSF1B",
      "IFNGR1", "IFNGR2", "TGFBR1", "TGFBR2", "TGFBR1", "TGFBR2",
      "IL1R1", "IL1R2", "IL4R", "IL13RA1", "IL12RB1", "IL12RB2", "IL15RA", "IL17RA",
      "IL21R", "IL1RL1", "IL18R1", "CSF1R", "CSF2RA",
      # Chemokines
      "CXCR4", "ACKR3", "CCR2", "CCR5", "CCR7", "CCR7", "CXCR3", "CXCR3", "CXCR5",
      # Growth factors
      "FLT1", "KDR", "EGFR", "EGFR", "FGFR1", "FGFR2", "PDGFRA", "PDGFRB",
      "MET", "BMPR1A", "BMPR1A", "FZD1", "ROR2",
      # Notch
      "NOTCH1", "NOTCH1", "NOTCH1", "NOTCH2",
      # Cell adhesion
      "ITGAL", "ITGA4", "CDH1", "CDH2",
      # TNF superfamily
      "TNFRSF10A", "FAS", "CD40", "TNFRSF4", "TNFRSF9",
      # Costimulatory
      "CD80", "ICOSLG",
      # MHC / Antigen presentation
      "CD8A", "CD8A", "CD4", "CD3D",
      # Semaphorins
      "NRP1", "CD72",
      # Galectins
      "CD44", "LGALS3BP",
      # Complement
      "C3AR1", "C5AR1",
      # Hedgehog
      "PTCH1", "PTCH1",
      # Ephrins
      "EPHA2", "EPHB2", "EPHB4"
    ),
    pathway = c(
      rep("Immune Checkpoint", 7),
      rep("Cytokine", 27),
      rep("Chemokine", 9),
      rep("Growth Factor", 13),
      rep("Notch", 4),
      rep("Cell Adhesion", 4),
      rep("TNF Superfamily", 5),
      rep("Costimulatory", 2),
      rep("Antigen Presentation", 4),
      rep("Semaphorin", 2),
      rep("Galectin", 2),
      rep("Complement", 2),
      rep("Hedgehog", 2),
      rep("Ephrin", 3)
    ),
    stringsAsFactors = FALSE
  )

  if (species == "mouse") {
    # Convert to mouse orthologs (capitalize first letter only)
    human_lr$ligand <- .to_mouse_gene(human_lr$ligand)
    human_lr$receptor <- .to_mouse_gene(human_lr$receptor)
  }

  return(human_lr)
}

cellCommUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Cell-Cell Communication Analysis"),
          p("Infer cell-cell communication from ligand-receptor co-expression across cell types."),

          fluidRow(
            column(3, selectInput(ns("cc_groupby"), "Cell type column", choices = NULL)),
            column(3, selectInput(ns("cc_species"), "Species",
              choices = c("Human" = "human", "Mouse" = "mouse"), selected = "human")),
            column(3, numericInput(ns("cc_min_expr"), "Min % cells expressing",
              value = 0.1, min = 0, max = 1, step = 0.05)),
            column(3, numericInput(ns("cc_pval_cutoff"), "P-value cutoff",
              value = 0.05, min = 0.001, max = 0.25, step = 0.01))
          ),

          fluidRow(
            column(3, numericInput(ns("cc_min_score"), "Min communication score",
              value = 0.1, min = 0, max = 5, step = 0.05)),
            column(3, numericInput(ns("cc_nperm"), "Permutations",
              value = 1000, min = 500, max = 10000, step = 50)),
            column(3, selectInput(ns("cc_pathways"), "Filter by pathway",
              choices = c("All"), multiple = TRUE, selected = "All")),
            column(3, numericInput(ns("cc_top_n"), "Top interactions to show",
              value = 30, min = 5, max = 200, step = 5))
          ),

          helpText("Estimated time: ~1 min per 1000 permutations for typical datasets. Minimum 500 for reliable p-values."),

          actionButton(ns("run_cc"), "Run Communication Analysis",
            class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),

    fluidRow(
      column(12,
        uiOutput(ns("cc_results_ui"))
      )
    )
  )
}

cellCommServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    cc_results <- reactiveVal(NULL)
    show_cc_results <- reactiveVal(FALSE)

    output$cc_results_ui <- renderUI({
      req(isTRUE(show_cc_results()))
      ns <- session$ns
      wellPanel(
        h4("Communication Results"),

        tabsetPanel(
          tabPanel("Chord Diagram",
            fluidRow(
              column(9, plotlyOutput(ns("cc_chord"), height = "650px")),
              column(3,
                h5("Display Options"),
                sliderInput(ns("cc_chord_min"), "Min score to display",
                  min = 0, max = 2, value = 0.1, step = 0.05),
                checkboxInput(ns("cc_chord_labels"), "Show labels", TRUE),
                br(),
                downloadButton(ns("download_cc_chord"), "Download Plot",
                  class = "btn-success", width = "100%")
              )
            )
          ),

          tabPanel("L-R Heatmap",
            fluidRow(
              column(9, plotlyOutput(ns("cc_heatmap"), height = "650px")),
              column(3,
                h5("Heatmap Options"),
                selectInput(ns("cc_heatmap_pathway"), "Pathway filter",
                  choices = c("All"), selected = "All"),
                numericInput(ns("cc_heatmap_top"), "Top pairs", value = 30, min = 5, max = 100),
                br(),
                downloadButton(ns("download_cc_heatmap"), "Download Heatmap",
                  class = "btn-success", width = "100%")
              )
            )
          ),

          tabPanel("Network",
            fluidRow(
              column(9, plotlyOutput(ns("cc_network"), height = "650px")),
              column(3,
                h5("Network Options"),
                selectInput(ns("cc_net_pathway"), "Pathway", choices = c("All")),
                sliderInput(ns("cc_net_min_score"), "Min edge score",
                  min = 0, max = 2, value = 0.2, step = 0.05),
                br(),
                downloadButton(ns("download_cc_network"), "Download Network",
                  class = "btn-success", width = "100%")
              )
            )
          ),

          tabPanel("Top Pathways",
            fluidRow(
              column(12, plotlyOutput(ns("cc_pathway_bar"), height = "500px"))
            ),
            br(),
            fluidRow(
              column(12, DT::dataTableOutput(ns("cc_pathway_table")))
            )
          ),

          tabPanel("Results Table",
            DT::dataTableOutput(ns("cc_table")),
            br(),
            fluidRow(
              column(6,
                downloadButton(ns("download_cc_csv"), "Download All Results (CSV)",
                  class = "btn-success", width = "100%")),
              column(6,
                downloadButton(ns("download_cc_sig"), "Download Significant (CSV)",
                  class = "btn-info", width = "100%"))
            )
          )
        )
      )
    })

    # Populate UI
    observe({
      req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "cc_groupby", choices = sort(cat_cols))
    })

    # Main analysis
    observeEvent(input$run_cc, {
      req(values$active_dataset, input$cc_groupby)

      if (!input$cc_groupby %in% colnames(values$active_dataset@meta.data)) {
        showNotification("Selected column not found in metadata", type = "error")
        return()
      }

      # Check for biological replication
      sample_col <- grep("orig.ident|sample|donor|batch", colnames(values$active_dataset@meta.data), value = TRUE, ignore.case = TRUE)[1]
      n_samples <- if (!is.na(sample_col)) length(unique(values$active_dataset@meta.data[[sample_col]])) else 1
      if (n_samples <= 1) {
        showNotification(
          "WARNING: Single-sample CCC analysis without biological replication is exploratory only and not suitable for publication. Consider adding biological replicates.",
          type = "warning", duration = 15
        )
      }

      withProgress(message = "Analyzing cell-cell communication...", value = 0, {
        tryCatch({
          obj <- values$active_dataset
          group_col <- input$cc_groupby
          species <- input$cc_species
          min_expr <- input$cc_min_expr
          nperm <- max(input$cc_nperm, 500)

          # ---- Try LIANA first for full analysis ----
          liana_result <- NULL
          if (requireNamespace("liana", quietly = TRUE)) {
            incProgress(0.1, detail = "Running LIANA analysis...")
            showNotification("LIANA detected - using multi-method consensus analysis",
                             type = "message", duration = 5)
            liana_result <- tryCatch({
              .run_ccc_liana(obj, group_col, species = species,
                            min_score = input$cc_min_score)
            }, error = function(e) {
              showNotification(paste("LIANA failed, falling back to built-in method:", e$message),
                               type = "warning", duration = 8)
              NULL
            })
          }

          if (!is.null(liana_result) && nrow(liana_result) > 0) {
            # LIANA succeeded — use its results directly
            res_df <- liana_result
            incProgress(0.9, detail = "LIANA analysis complete!")
            showNotification(
              paste0("LIANA analysis complete! ", nrow(res_df), " interactions found"),
              type = "message", duration = 8
            )
          } else {
            # ---- Fallback: built-in permutation method ----
            if (requireNamespace("liana", quietly = TRUE)) {
              showNotification("LIANA returned no results, using built-in method",
                               type = "warning", duration = 5)
            }

          incProgress(0.1, detail = "Loading L-R database...")
          # Use enhanced fallback chain: LIANA resource -> OmnipathR -> built-in
          lr_db <- .get_lr_database_liana(species)
          db_size <- nrow(lr_db)
          showNotification(
            paste0("Using L-R database with ", db_size, " pairs"),
            type = "message", duration = 5
          )

          incProgress(0.2, detail = "Computing expression fractions...")

          # Get expression data
          assay_name <- DefaultAssay(obj)
          expr_mat <- safe_get_assay_data(obj, "data", assay = assay_name)
          groups <- as.character(obj@meta.data[[group_col]])
          cell_types <- sort(unique(groups))

          if (length(cell_types) < 2) {
            showNotification("Need at least 2 cell types for communication analysis", type = "error")
            return()
          }

          # Filter L-R pairs to those present in dataset
          all_genes <- rownames(expr_mat)
          lr_db <- lr_db[lr_db$ligand %in% all_genes & lr_db$receptor %in% all_genes, ]

          if (nrow(lr_db) == 0) {
            showNotification("No ligand-receptor pairs found in the dataset", type = "error")
            return()
          }

          incProgress(0.3, detail = "Computing communication scores...")

          # Compute mean expression and fraction per cell type
          # Reverse log-normalization to get linear-scale expression for CCC scoring
          # L-R interaction scores should be computed on linear-scale, not log-space
          # Library size normalization is inherited from the upstream NormalizeData step
          # Expression values are in counts-per-10K (CP10K) units after expm1() reversal,
          # providing a common reference scale across cell types regardless of library size
          ct_mean <- list()
          ct_frac <- list()
          for (ct in cell_types) {
            cells_idx <- which(groups == ct)
            sub_mat <- expr_mat[, cells_idx, drop = FALSE]
            ct_mean[[ct]] <- Matrix::rowMeans(expm1(sub_mat))
            ct_frac[[ct]] <- Matrix::rowMeans(sub_mat > 0)
          }

          # Compute communication scores
          results_list <- list()
          total_pairs <- nrow(lr_db) * length(cell_types) * length(cell_types)
          counter <- 0

          for (i in seq_len(nrow(lr_db))) {
            lig <- lr_db$ligand[i]
            rec <- lr_db$receptor[i]
            pw <- lr_db$pathway[i]

            for (src in cell_types) {
              for (tgt in cell_types) {
                counter <- counter + 1

                lig_expr <- ct_mean[[src]][lig]
                rec_expr <- ct_mean[[tgt]][rec]
                lig_frac <- ct_frac[[src]][lig]
                rec_frac <- ct_frac[[tgt]][rec]

                # Skip if expression is below threshold
                if (lig_frac < min_expr || rec_frac < min_expr) next

                # Communication score = product of fraction expressing and mean expression
                score <- (lig_frac * lig_expr) * (rec_frac * rec_expr)

                if (score <= 0) next

                results_list[[length(results_list) + 1]] <- data.frame(
                  source = src,
                  target = tgt,
                  ligand = lig,
                  receptor = rec,
                  pathway = pw,
                  lr_pair = paste0(lig, "_", rec),
                  score = score,
                  ligand_expr = lig_expr,
                  receptor_expr = rec_expr,
                  ligand_frac = lig_frac,
                  receptor_frac = rec_frac,
                  stringsAsFactors = FALSE
                )
              }
            }
          }

          if (length(results_list) == 0) {
            showNotification("No significant communication found with current thresholds", type = "warning")
            return()
          }

          res_df <- do.call(rbind, results_list)

          incProgress(0.6, detail = "Computing p-values by permutation...")

          # Permutation test for significance
          pvals <- rep(1, nrow(res_df))
          if (nperm > 0 && nrow(res_df) > 0) {
            set.seed(42)
            # Pre-compute permuted scores
            perm_scores <- matrix(0, nrow = nrow(res_df), ncol = nperm)
            for (p in seq_len(nperm)) {
              perm_groups <- sample(groups)
              perm_mean <- list()
              for (ct in cell_types) {
                cells_idx <- which(perm_groups == ct)
                if (length(cells_idx) > 0) {
                  perm_mean[[ct]] <- Matrix::rowMeans(expm1(expr_mat[, cells_idx, drop = FALSE]))
                } else {
                  perm_mean[[ct]] <- rep(0, nrow(expr_mat))
                  names(perm_mean[[ct]]) <- rownames(expr_mat)
                }
              }
              for (r in seq_len(nrow(res_df))) {
                perm_lig_expr <- perm_mean[[res_df$source[r]]][res_df$ligand[r]]
                perm_rec_expr <- perm_mean[[res_df$target[r]]][res_df$receptor[r]]
                ps <- (res_df$ligand_frac[r] * perm_lig_expr) * (res_df$receptor_frac[r] * perm_rec_expr)
                perm_scores[r, p] <- ifelse(is.na(ps), 0, ps)
              }
            }
            for (r in seq_len(nrow(res_df))) {
              pvals[r] <- (sum(perm_scores[r, ] >= res_df$score[r]) + 1) / (nperm + 1)
            }
          }

          res_df$pval <- pvals
          res_df$pval_adj <- p.adjust(pvals, method = "BH")

          # Filter by score and p-value
          res_df <- res_df %>%
            filter(score >= input$cc_min_score) %>%
            arrange(pval_adj, desc(score))

          incProgress(0.9, detail = "Finalizing results...")

          } # end of else (built-in fallback path)

          cc_results(res_df)
          values$last_cc_results <- res_df
          show_cc_results(TRUE)

          # Update pathway filter choices
          avail_pathways <- sort(unique(res_df$pathway))
          updateSelectInput(session, "cc_pathways", choices = c("All", avail_pathways), selected = "All")
          updateSelectInput(session, "cc_heatmap_pathway", choices = c("All", avail_pathways), selected = "All")
          updateSelectInput(session, "cc_net_pathway", choices = c("All", avail_pathways), selected = "All")

          n_sig <- sum(res_df$pval_adj < input$cc_pval_cutoff, na.rm = TRUE)
          showNotification(paste("Communication analysis complete!", nrow(res_df),
            "interactions,", n_sig, "significant"), type = "message")

          incProgress(1, detail = "Done!")
          gc()

        }, error = function(e) {
          showNotification(paste("Communication analysis failed:", e$message), type = "error")
        })
      })
    })

    # Filtered results helper
    .get_filtered_cc <- function() {
      req(cc_results())
      df <- cc_results()
      if (!is.null(input$cc_pathways) && !("All" %in% input$cc_pathways)) {
        df <- df %>% filter(pathway %in% input$cc_pathways)
      }
      df
    }

    # Chord diagram (circular network via plotly)
    output$cc_chord <- renderPlotly({
      tryCatch({
        df <- .get_filtered_cc()
        req(nrow(df) > 0)

        min_score <- input$cc_chord_min
        agg <- df %>%
          filter(score >= min_score, pval_adj < input$cc_pval_cutoff) %>%
          group_by(source, target) %>%
          summarise(total_score = sum(score), n_pairs = n(), .groups = "drop")

        if (nrow(agg) == 0) {
          return(plotly_empty() %>% layout(title = "No interactions pass current filters"))
        }

        cell_types <- sort(unique(c(agg$source, agg$target)))
        n_ct <- length(cell_types)
        angles <- seq(0, 2 * pi, length.out = n_ct + 1)[1:n_ct]
        node_x <- cos(angles)
        node_y <- sin(angles)
        names(node_x) <- cell_types
        names(node_y) <- cell_types

        colors <- get_color_palette(n_ct)

        # Build edge traces
        p <- plot_ly()

        max_score <- max(agg$total_score)
        for (i in seq_len(nrow(agg))) {
          src <- agg$source[i]
          tgt <- agg$target[i]
          width <- 1 + 4 * (agg$total_score[i] / max_score)
          src_idx <- which(cell_types == src)
          p <- p %>% add_trace(
            x = c(node_x[src], node_x[tgt]),
            y = c(node_y[src], node_y[tgt]),
            type = "scatter", mode = "lines",
            line = list(width = width, color = colors[src_idx]),
            opacity = 0.5,
            hoverinfo = "text",
            text = paste0(src, " -> ", tgt, "\nScore: ", round(agg$total_score[i], 3),
                          "\nPairs: ", agg$n_pairs[i]),
            showlegend = FALSE
          )
        }

        # Add nodes
        hover_text <- sapply(cell_types, function(ct) {
          out <- sum(agg$total_score[agg$source == ct], na.rm = TRUE)
          inc <- sum(agg$total_score[agg$target == ct], na.rm = TRUE)
          paste0(ct, "\nOutgoing: ", round(out, 2), "\nIncoming: ", round(inc, 2))
        })

        p <- p %>% add_trace(
          x = node_x, y = node_y,
          type = "scatter", mode = "markers+text",
          marker = list(size = 20, color = colors, line = list(color = "black", width = 1.5)),
          text = if (input$cc_chord_labels) cell_types else rep("", n_ct),
          textposition = "top center",
          textfont = list(size = 10),
          hovertext = hover_text,
          hoverinfo = "text",
          showlegend = FALSE
        )

        p <- p %>% layout(
          title = "Cell-Cell Communication Network",
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = "",
                       scaleanchor = "x", scaleratio = 1),
          plot_bgcolor = "white",
          autosize = TRUE
        )
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Cell Communication", figure_type = "CCC Chord",
                         label = "CCC Chord Diagram", source_module = "Cell Communication")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Chord error:", e$message))
      })
    })

    # L-R Heatmap
    output$cc_heatmap <- renderPlotly({
      tryCatch({
        df <- cc_results()
        req(df, nrow(df) > 0)

        pw_filter <- input$cc_heatmap_pathway
        if (!is.null(pw_filter) && pw_filter != "All") {
          df <- df %>% filter(pathway == pw_filter)
        }

        # Show top pairs
        top_pairs <- df %>%
          filter(pval_adj < input$cc_pval_cutoff) %>%
          arrange(pval_adj, desc(score)) %>%
          head(input$cc_heatmap_top)

        if (nrow(top_pairs) == 0) {
          return(plotly_empty() %>% layout(title = "No significant pairs to display"))
        }

        # Build a matrix: rows = LR pairs, cols = source->target
        top_pairs$interaction <- paste0(top_pairs$source, " -> ", top_pairs$target)
        top_pairs$pair_label <- paste0(top_pairs$ligand, " | ", top_pairs$receptor)

        hm_data <- top_pairs %>%
          select(pair_label, interaction, score) %>%
          pivot_wider(names_from = interaction, values_from = score, values_fill = 0)

        mat <- as.matrix(hm_data[, -1, drop = FALSE])
        rownames(mat) <- hm_data$pair_label
        if (ncol(mat) == 0 || nrow(mat) == 0) {
          return(plotly_empty() %>% layout(title = "Not enough data for heatmap"))
        }

        p <- plot_ly(
          x = colnames(mat),
          y = rownames(mat),
          z = mat,
          type = "heatmap",
          colorscale = "Viridis",
          hovertemplate = "Pair: %{y}<br>Interaction: %{x}<br>Score: %{z:.3f}<extra></extra>"
        ) %>% layout(
          title = "Ligand-Receptor Expression Heatmap",
          xaxis = list(title = "Source -> Target", tickangle = 45),
          yaxis = list(title = "L-R Pair"),
          margin = list(l = 150, b = 120),
          autosize = TRUE
        )
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Cell Communication", figure_type = "CCC Heatmap",
                         label = "L-R Expression Heatmap", source_module = "Cell Communication")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Heatmap error:", e$message))
      })
    })

    # Network visualization
    output$cc_network <- renderPlotly({
      tryCatch({
        df <- cc_results()
        req(df, nrow(df) > 0)

        pw_filter <- input$cc_net_pathway
        if (!is.null(pw_filter) && pw_filter != "All") {
          df <- df %>% filter(pathway == pw_filter)
        }

        df <- df %>% filter(score >= input$cc_net_min_score, pval_adj < input$cc_pval_cutoff)

        if (nrow(df) == 0) {
          return(plotly_empty() %>% layout(title = "No interactions pass current filters"))
        }

        # Aggregate by source-target
        agg <- df %>%
          group_by(source, target) %>%
          summarise(total_score = sum(score), n_pairs = n(),
                    top_pairs = paste(head(lr_pair, 3), collapse = ", "), .groups = "drop")

        # Build igraph
        if (!requireNamespace("igraph", quietly = TRUE)) {
          return(plotly_empty() %>% layout(title = "igraph package required for network visualization"))
        }

        nodes <- sort(unique(c(agg$source, agg$target)))
        g <- igraph::graph_from_data_frame(agg, directed = TRUE, vertices = data.frame(name = nodes))
        layout_coords <- igraph::layout_with_fr(g)

        node_x <- layout_coords[, 1]
        node_y <- layout_coords[, 2]
        names(node_x) <- nodes
        names(node_y) <- nodes

        # Compute node strength
        node_strength <- sapply(nodes, function(n) {
          sum(agg$total_score[agg$source == n | agg$target == n])
        })
        node_size <- 10 + 25 * (node_strength - min(node_strength)) /
                     max(1, max(node_strength) - min(node_strength))

        colors <- get_color_palette(length(nodes))

        p <- plot_ly()

        # Edges
        max_score <- max(agg$total_score)
        for (i in seq_len(nrow(agg))) {
          src <- agg$source[i]
          tgt <- agg$target[i]
          w <- 0.5 + 3 * (agg$total_score[i] / max_score)
          p <- p %>% add_trace(
            x = c(node_x[src], node_x[tgt]),
            y = c(node_y[src], node_y[tgt]),
            type = "scatter", mode = "lines",
            line = list(width = w, color = "rgba(150,150,150,0.6)"),
            hoverinfo = "text",
            text = paste0(src, " -> ", tgt, "\nScore: ", round(agg$total_score[i], 3),
                          "\nPairs: ", agg$n_pairs[i], "\nTop: ", agg$top_pairs[i]),
            showlegend = FALSE
          )
        }

        # Nodes
        p <- p %>% add_trace(
          x = node_x, y = node_y,
          type = "scatter", mode = "markers+text",
          marker = list(size = node_size, color = colors,
                        line = list(color = "black", width = 1)),
          text = nodes, textposition = "top center",
          textfont = list(size = 10),
          hovertext = paste0(nodes, "\nStrength: ", round(node_strength, 2)),
          hoverinfo = "text",
          showlegend = FALSE
        )

        p %>% layout(
          title = paste("Signaling Network -", ifelse(pw_filter == "All", "All Pathways", pw_filter)),
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE, title = ""),
          plot_bgcolor = "white",
          autosize = TRUE
        )
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Network error:", e$message))
      })
    })

    # Top pathways bar chart
    output$cc_pathway_bar <- renderPlotly({
      tryCatch({
        df <- cc_results()
        req(df, nrow(df) > 0)

        pw_summary <- df %>%
          filter(pval_adj < input$cc_pval_cutoff) %>%
          group_by(pathway) %>%
          summarise(
            n_interactions = n(),
            mean_score = mean(score),
            total_score = sum(score),
            .groups = "drop"
          ) %>%
          arrange(desc(total_score))

        if (nrow(pw_summary) == 0) {
          return(plotly_empty() %>% layout(title = "No significant pathways"))
        }

        # Defactorize to prevent plotly range() error on factors
        pw_summary$pathway <- as.character(pw_summary$pathway)
        pw_summary$n_interactions <- as.numeric(pw_summary$n_interactions)
        pw_summary$total_score <- as.numeric(pw_summary$total_score)
        pw_top <- head(pw_summary, 25)

        plot_ly(
          data = pw_top,
          x = ~pathway,
          y = ~total_score,
          color = ~n_interactions,
          type = "bar",
          hovertemplate = "Pathway: %{x}<br>Total Score: %{y:.3f}<br>Interactions: %{marker.color}<extra></extra>"
        ) %>% layout(
          title = "Top Signaling Pathways",
          xaxis = list(title = "Pathway",
                       categoryorder = "array",
                       categoryarray = pw_top$pathway[order(pw_top$total_score)]),
          yaxis = list(title = "Total Communication Score"),
          autosize = TRUE
        )
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Pathway plot error:", e$message))
      })
    })

    # Pathway summary table
    output$cc_pathway_table <- DT::renderDataTable({
      req(cc_results())
      tryCatch({
        cc_results() %>%
          filter(pval_adj < input$cc_pval_cutoff) %>%
          group_by(pathway) %>%
          summarise(
            n_interactions = n(),
            n_lr_pairs = n_distinct(lr_pair),
            mean_score = round(mean(score), 4),
            max_score = round(max(score), 4),
            mean_pval = signif(mean(pval_adj), 3),
            .groups = "drop"
          ) %>%
          arrange(desc(n_interactions))
      }, error = function(e) data.frame(Error = e$message))
    }, options = list(pageLength = 15, scrollX = TRUE))

    # Results table
    output$cc_table <- DT::renderDataTable({
      req(cc_results())
      tryCatch({
        tbl <- cc_results() %>%
          mutate(
            score = round(as.numeric(score), 4),
            ligand_expr = round(as.numeric(ligand_expr), 3),
            receptor_expr = round(as.numeric(receptor_expr), 3),
            ligand_frac = round(as.numeric(ligand_frac), 3),
            receptor_frac = round(as.numeric(receptor_frac), 3),
            pval = signif(as.numeric(pval), 3),
            pval_adj = signif(as.numeric(pval_adj), 3)
          )
        if (exists("register_csv_export", mode = "function")) {
          register_csv_export(values, data = tbl, description = "CCC Interactions",
                              source_module = "Downstream")
        }
        tbl
      }, error = function(e) data.frame(Error = e$message))
    }, options = list(pageLength = 25, scrollX = TRUE, scrollY = "400px"),
    filter = "top")

    # Downloads
    output$download_cc_csv <- downloadHandler(
      filename = function() paste0("cell_communication_all_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          write.csv(cc_results(), file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Download error:", e$message), type = "error")
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )

    output$download_cc_sig <- downloadHandler(
      filename = function() paste0("cell_communication_significant_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          sig <- cc_results() %>% filter(pval_adj < input$cc_pval_cutoff)
          write.csv(sig, file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Download error:", e$message), type = "error")
          write.csv(data.frame(Error = e$message), file, row.names = FALSE)
        })
      }
    )

    output$download_cc_chord <- downloadHandler(
      filename = function() paste0("cc_chord_", Sys.Date(), ".html"),
      content = function(file) {
        p <- tryCatch({
          df <- .get_filtered_cc()
          agg <- df %>%
            filter(score >= input$cc_chord_min, pval_adj < input$cc_pval_cutoff) %>%
            group_by(source, target) %>%
            summarise(total_score = sum(score), .groups = "drop")
          plotly_empty() %>% layout(title = "Chord Diagram Export")
        }, error = function(e) plotly_empty())
        htmlwidgets::saveWidget(p, file)
      }
    )

    output$download_cc_heatmap <- downloadHandler(
      filename = function() paste0("cc_heatmap_", Sys.Date(), ".html"),
      content = function(file) {
        p <- tryCatch(plotly_empty(), error = function(e) plotly_empty())
        htmlwidgets::saveWidget(p, file)
      }
    )

    output$download_cc_network <- downloadHandler(
      filename = function() paste0("cc_network_", Sys.Date(), ".html"),
      content = function(file) {
        p <- tryCatch(plotly_empty(), error = function(e) plotly_empty())
        htmlwidgets::saveWidget(p, file)
      }
    )
  })
}

# ---------------------------------------------------------------------------
# LIANA-powered CCC: uses liana::liana_wrap() when available
# Returns results in the same format as .run_ccc_auto() for downstream compatibility
# ---------------------------------------------------------------------------
.run_ccc_liana <- function(obj, group_col, species = "human",
                           min_score = 0.1) {
  tryCatch({
    if (!requireNamespace("liana", quietly = TRUE)) {
      stop("LIANA package not available")
    }

    message("[CCC] Running LIANA-powered cell-cell communication analysis...")

    # LIANA expects Idents to be set to the grouping variable
    Seurat::Idents(obj) <- group_col

    cell_types <- sort(unique(as.character(obj@meta.data[[group_col]])))
    if (length(cell_types) < 2) {
      stop("Need at least 2 cell types for communication analysis")
    }

    # Run LIANA with multiple methods for consensus scoring
    liana_res <- liana::liana_wrap(
      obj,
      method = c("natmi", "cellphonedb", "sca"),
      resource = "Consensus",
      idents_col = group_col
    )

    # Aggregate results across methods
    liana_agg <- liana::liana_aggregate(liana_res)

    if (is.null(liana_agg) || nrow(liana_agg) == 0) {
      stop("LIANA returned no results")
    }

    # Standardize to match .run_ccc_auto() output format
    # LIANA aggregate columns: source, target, ligand.complex, receptor.complex,
    #   aggregate_rank, natmi, cellphonedb, sca, ...
    lig_col <- intersect(c("ligand.complex", "ligand"), colnames(liana_agg))[1]
    rec_col <- intersect(c("receptor.complex", "receptor"), colnames(liana_agg))[1]
    src_col <- intersect(c("source", "sender"), colnames(liana_agg))[1]
    tgt_col <- intersect(c("target", "receiver"), colnames(liana_agg))[1]

    if (is.na(lig_col) || is.na(rec_col) || is.na(src_col) || is.na(tgt_col)) {
      stop("Unexpected LIANA output column names: ", paste(colnames(liana_agg), collapse = ", "))
    }

    # Use aggregate_rank as score (lower = better, so invert for consistency)
    score_col <- if ("aggregate_rank" %in% colnames(liana_agg)) "aggregate_rank" else NULL

    res_df <- data.frame(
      source   = as.character(liana_agg[[src_col]]),
      target   = as.character(liana_agg[[tgt_col]]),
      ligand   = as.character(liana_agg[[lig_col]]),
      receptor = as.character(liana_agg[[rec_col]]),
      pathway  = "LIANA",
      lr_pair  = paste0(liana_agg[[lig_col]], "_", liana_agg[[rec_col]]),
      score    = if (!is.null(score_col)) (1 - liana_agg[[score_col]]) else 1,
      ligand_expr  = NA_real_,
      receptor_expr = NA_real_,
      ligand_frac  = NA_real_,
      receptor_frac = NA_real_,
      stringsAsFactors = FALSE
    )

    # LIANA's aggregate_rank already functions as a p-value (lower = more significant)
    if (!is.null(score_col)) {
      res_df$pval     <- liana_agg[[score_col]]
      res_df$pval_adj <- p.adjust(liana_agg[[score_col]], method = "BH")
    } else {
      res_df$pval     <- NA_real_
      res_df$pval_adj <- NA_real_
    }

    # Filter by minimum score
    res_df <- res_df[res_df$score >= min_score, ]
    res_df <- res_df[order(res_df$pval_adj, -res_df$score), ]

    message("[CCC] LIANA analysis complete: ", nrow(res_df), " interactions found")
    res_df
  }, error = function(e) {
    message("[CCC] LIANA CCC failed: ", e$message)
    NULL
  })
}

# Standalone CCC: returns data.frame of interactions (built-in permutation method)
.run_ccc_auto <- function(obj, group_col, species = "human",
                           min_expr = 0.1, nperm = 1000, min_score = 0.1) {
  tryCatch({
    # Use enhanced fallback chain: LIANA resource -> OmnipathR -> built-in
    lr_db <- .get_lr_database_liana(species)
    assay_name <- DefaultAssay(obj)
    expr_mat <- safe_get_assay_data(obj, "data", assay = assay_name)
    groups <- as.character(obj@meta.data[[group_col]])
    cell_types <- sort(unique(groups))

    if (length(cell_types) < 2) return(NULL)

    # Enforce minimum permutations for reliable p-values
    nperm <- max(nperm, 500)

    all_genes <- rownames(expr_mat)
    lr_db <- lr_db[lr_db$ligand %in% all_genes & lr_db$receptor %in% all_genes, ]
    if (nrow(lr_db) == 0) return(NULL)

    # Reverse log-normalization to get linear-scale expression for CCC scoring
    # L-R interaction scores should be computed on linear-scale, not log-space
    # Library size normalization is inherited from the upstream NormalizeData step
    ct_mean <- list(); ct_frac <- list()
    for (ct in cell_types) {
      cells_idx <- which(groups == ct)
      sub_mat <- expr_mat[, cells_idx, drop = FALSE]
      ct_mean[[ct]] <- Matrix::rowMeans(expm1(sub_mat))
      ct_frac[[ct]] <- Matrix::rowMeans(sub_mat > 0)
    }

    results_list <- list()
    for (i in seq_len(nrow(lr_db))) {
      lig <- lr_db$ligand[i]; rec <- lr_db$receptor[i]; pw <- lr_db$pathway[i]
      for (src in cell_types) {
        for (tgt in cell_types) {
          lig_frac <- ct_frac[[src]][lig]; rec_frac <- ct_frac[[tgt]][rec]
          if (lig_frac < min_expr || rec_frac < min_expr) next
          lig_expr <- ct_mean[[src]][lig]; rec_expr <- ct_mean[[tgt]][rec]
          score <- (lig_frac * lig_expr) * (rec_frac * rec_expr)
          if (score <= 0) next
          results_list[[length(results_list) + 1]] <- data.frame(
            source = src, target = tgt, ligand = lig, receptor = rec,
            pathway = pw, lr_pair = paste0(lig, "_", rec), score = score,
            ligand_expr = lig_expr, receptor_expr = rec_expr,
            ligand_frac = lig_frac, receptor_frac = rec_frac,
            stringsAsFactors = FALSE)
        }
      }
    }
    if (length(results_list) == 0) return(NULL)
    res_df <- do.call(rbind, results_list)

    # Permutation test (also uses expm1 for consistency)
    pvals <- rep(1, nrow(res_df))
    if (nperm > 0 && nrow(res_df) > 0) {
      set.seed(42)
      perm_scores <- matrix(0, nrow = nrow(res_df), ncol = nperm)
      for (p in seq_len(nperm)) {
        perm_groups <- sample(groups)
        perm_mean <- list()
        for (ct in cell_types) {
          cells_idx <- which(perm_groups == ct)
          if (length(cells_idx) > 0) {
            perm_mean[[ct]] <- Matrix::rowMeans(expm1(expr_mat[, cells_idx, drop = FALSE]))
          } else {
            perm_mean[[ct]] <- rep(0, nrow(expr_mat))
            names(perm_mean[[ct]]) <- rownames(expr_mat)
          }
        }
        for (r in seq_len(nrow(res_df))) {
          perm_lig_expr <- perm_mean[[res_df$source[r]]][res_df$ligand[r]]
          perm_rec_expr <- perm_mean[[res_df$target[r]]][res_df$receptor[r]]
          ps <- (res_df$ligand_frac[r] * perm_lig_expr) * (res_df$receptor_frac[r] * perm_rec_expr)
          perm_scores[r, p] <- ifelse(is.na(ps), 0, ps)
        }
      }
      for (r in seq_len(nrow(res_df))) {
        pvals[r] <- (sum(perm_scores[r, ] >= res_df$score[r]) + 1) / (nperm + 1)
      }
    }
    res_df$pval <- pvals
    res_df$pval_adj <- p.adjust(pvals, method = "BH")
    res_df <- res_df %>% filter(score >= min_score) %>% arrange(pval_adj, desc(score))
    res_df
  }, error = function(e) {
    message("Auto CCC error: ", e$message)
    NULL
  })
}

# CCC chord diagram (standalone)
.plot_ccc_chord_auto <- function(cc_df, min_score = 0.1, pval_cutoff = 0.05) {
  tryCatch({
    if (is.null(cc_df) || nrow(cc_df) == 0)
      return(plotly_empty() %>% layout(title = list(text = "No CCC results available", font = list(size = 14))))

    agg <- cc_df %>%
      filter(score >= min_score, pval_adj < pval_cutoff) %>%
      group_by(source, target) %>%
      summarise(total_score = sum(score), n_pairs = n(), .groups = "drop")
    if (nrow(agg) == 0)
      return(plotly_empty() %>% layout(title = list(
        text = "No significant interactions (adj. p < 0.05)", font = list(size = 14))))

    cell_types <- sort(unique(c(agg$source, agg$target)))
    n_ct <- length(cell_types)
    angles <- seq(0, 2 * pi, length.out = n_ct + 1)[1:n_ct]
    node_x <- cos(angles); node_y <- sin(angles)
    names(node_x) <- cell_types; names(node_y) <- cell_types

    # Node size proportional to total interaction strength
    node_strength <- sapply(cell_types, function(ct)
      sum(agg$total_score[agg$source == ct | agg$target == ct], na.rm = TRUE))
    max_str <- max(node_strength, na.rm = TRUE)
    node_sizes <- 14 + 16 * (node_strength / max_str)

    colors <- get_color_palette(n_ct)
    names(colors) <- cell_types
    p <- plot_ly()
    max_score <- max(agg$total_score, na.rm = TRUE)

    # Curved edges (spline through offset midpoint)
    for (i in seq_len(nrow(agg))) {
      src <- agg$source[i]; tgt <- agg$target[i]
      w <- 0.5 + 3.5 * (agg$total_score[i] / max_score)
      alpha <- 0.25 + 0.35 * (agg$total_score[i] / max_score)
      src_idx <- which(cell_types == src)
      col_rgb <- col2rgb(colors[src_idx])
      lc <- sprintf("rgba(%d,%d,%d,%.2f)", col_rgb[1], col_rgb[2], col_rgb[3], alpha)
      mx <- (node_x[src] + node_x[tgt]) / 2 * 0.35
      my <- (node_y[src] + node_y[tgt]) / 2 * 0.35
      p <- p %>% add_trace(
        x = c(node_x[src], mx, node_x[tgt]),
        y = c(node_y[src], my, node_y[tgt]),
        type = "scatter", mode = "lines",
        line = list(width = w, color = lc, shape = "spline"),
        hoverinfo = "text",
        text = paste0(src, " \u2192 ", tgt,
                      "\nScore: ", round(agg$total_score[i], 3),
                      "\nL-R pairs: ", agg$n_pairs[i]),
        showlegend = FALSE)
    }

    hover_text <- sapply(cell_types, function(ct) {
      out_n <- sum(agg$n_pairs[agg$source == ct], na.rm = TRUE)
      in_n <- sum(agg$n_pairs[agg$target == ct], na.rm = TRUE)
      out_s <- round(sum(agg$total_score[agg$source == ct], na.rm = TRUE), 2)
      in_s <- round(sum(agg$total_score[agg$target == ct], na.rm = TRUE), 2)
      paste0("<b>", ct, "</b>",
             "<br>Outgoing: ", out_n, " pairs (score: ", out_s, ")",
             "<br>Incoming: ", in_n, " pairs (score: ", in_s, ")")
    })
    p <- p %>% add_trace(
      x = node_x, y = node_y, type = "scatter", mode = "markers+text",
      marker = list(size = node_sizes, color = colors,
                    line = list(color = "white", width = 2.5), opacity = 0.92),
      text = cell_types, textposition = "top center",
      textfont = list(size = 12, color = "#222", family = "Arial", weight = "bold"),
      hovertext = hover_text, hoverinfo = "text", showlegend = FALSE)

    n_sig <- nrow(agg); total_pairs <- sum(agg$n_pairs)
    p %>% layout(
      title = list(
        text = paste0("Cell-Cell Communication Network<br>",
                      "<span style='font-size:11px;color:#777'>",
                      n_sig, " significant axes \u2022 ", total_pairs, " L-R pairs</span>"),
        font = list(size = 14, family = "Arial")),
      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                   title = "", range = c(-1.45, 1.45)),
      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE,
                   title = "", scaleanchor = "x", scaleratio = 1, range = c(-1.45, 1.45)),
      plot_bgcolor = "white", paper_bgcolor = "white",
      margin = list(t = 55, b = 15, l = 15, r = 15), autosize = TRUE)
  }, error = function(e) {
    plotly_empty() %>% layout(title = list(text = paste("CCC error:", e$message), font = list(size = 12)))
  })
}

# CCC heatmap (standalone)
.plot_ccc_heatmap_auto <- function(cc_df, top_n = 30, pval_cutoff = 0.05) {
  tryCatch({
    if (is.null(cc_df) || nrow(cc_df) == 0)
      return(plotly_empty() %>% layout(title = list(text = "No CCC results", font = list(size = 14))))

    sig <- cc_df %>% filter(pval_adj < pval_cutoff)
    if (nrow(sig) == 0)
      return(plotly_empty() %>% layout(title = list(
        text = "No significant L-R pairs (adj. p < 0.05)", font = list(size = 14))))

    top_pairs <- sig %>% arrange(pval_adj, desc(score)) %>% head(top_n)
    top_pairs$interaction <- paste0(top_pairs$source, " \u2192 ", top_pairs$target)
    top_pairs$pair_label <- paste0(top_pairs$ligand, " - ", top_pairs$receptor)
    hm_data <- top_pairs %>%
      select(pair_label, interaction, score) %>%
      tidyr::pivot_wider(names_from = interaction, values_from = score, values_fill = 0)
    mat <- as.matrix(hm_data[, -1, drop = FALSE])
    rownames(mat) <- hm_data$pair_label
    if (ncol(mat) == 0 || nrow(mat) == 0) {
      return(plotly_empty() %>% layout(title = list(text = "Not enough data for heatmap", font = list(size = 14))))
    }
    if (nrow(mat) > 2) {
      ro <- tryCatch(hclust(dist(mat))$order, error = function(e) seq_len(nrow(mat)))
      mat <- mat[ro, , drop = FALSE]
    }
    if (ncol(mat) > 2) {
      co <- tryCatch(hclust(dist(t(mat)))$order, error = function(e) seq_len(ncol(mat)))
      mat <- mat[, co, drop = FALSE]
    }

    plot_ly(x = colnames(mat), y = rownames(mat), z = mat, type = "heatmap",
            colorscale = list(c(0,"#f7fbff"), c(0.15,"#deebf7"), c(0.3,"#c6dbef"),
                              c(0.5,"#6baed6"), c(0.7,"#3182bd"),
                              c(0.85,"#2171b5"), c(1,"#08306b")),
            colorbar = list(title = list(text = "L-R\nScore", font = list(size = 10)),
                            len = 0.5, thickness = 12),
            hovertemplate = paste0("<b>%{y}</b><br>Axis: %{x}<br>",
                                  "Score: %{z:.3f}<extra></extra>")) %>%
      layout(
        title = list(
          text = paste0("Top Ligand-Receptor Interactions<br>",
                        "<span style='font-size:11px;color:#777'>",
                        nrow(sig), " significant pairs (adj. p < 0.05)</span>"),
          font = list(size = 14, family = "Arial")),
        xaxis = list(title = "", tickangle = 45, tickfont = list(size = 10)),
        yaxis = list(title = "", tickfont = list(size = 10)),
        margin = list(l = 170, b = 140, t = 55, r = 40),
        paper_bgcolor = "white", plot_bgcolor = "white", autosize = TRUE)
  }, error = function(e) {
    plotly_empty() %>% layout(title = list(text = paste("Heatmap error:", e$message), font = list(size = 12)))
  })
}
