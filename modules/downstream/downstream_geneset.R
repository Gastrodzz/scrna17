# ──────────────────────────────────────────────────────────────────────────────
# downstream_geneset.R — Gene Set / Module Scoring
#
# Provides AddModuleScore-based gene set scoring with support for:
#   - Custom gene sets (comma-separated)
#   - Hallmark / canonical pathway sets via msigdbr
#   - UCell scoring (when available)
#   - Visualization: violin, feature, ridge plots
# ──────────────────────────────────────────────────────────────────────────────

geneSetScoringUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Gene Set / Module Scoring"),
          p("Score cells for expression of custom or curated gene sets using Seurat's AddModuleScore."),

          fluidRow(
            column(6,
              selectInput(ns("gs_source"), "Gene Set Source",
                choices = c(
                  "Custom gene list" = "custom",
                  "MSigDB Hallmark" = "hallmark",
                  "MSigDB C2 Canonical Pathways" = "c2",
                  "MSigDB C5 GO Biological Process" = "c5_bp",
                  "MSigDB C7 Immunologic" = "c7"
                ),
                selected = "custom")
            ),
            column(6,
              selectInput(ns("gs_species"), "Species",
                choices = c("Human" = "Homo sapiens", "Mouse" = "Mus musculus"),
                selected = "Homo sapiens")
            )
          ),

          # Custom gene list input
          conditionalPanel(
            condition = "input.gs_source == 'custom'",
            ns = ns,
            textAreaInput(ns("custom_genes"), "Enter gene names (comma or newline separated)",
              value = "", rows = 4, width = "100%",
              placeholder = "e.g. CD3D, CD3E, CD4, CD8A, FOXP3"),
            textInput(ns("custom_set_name"), "Gene set name", value = "CustomScore")
          ),

          # MSigDB pathway selector
          conditionalPanel(
            condition = "input.gs_source != 'custom'",
            ns = ns,
            uiOutput(ns("msigdb_pathway_selector"))
          ),

          fluidRow(
            column(4,
              selectInput(ns("scoring_method"), "Scoring Method",
                choices = c("AddModuleScore (Seurat)" = "seurat"),
                selected = "seurat")
            ),
            column(4,
              numericInput(ns("ctrl_genes"), "Control genes per bin",
                value = 100, min = 10, max = 500, step = 10)
            ),
            column(4,
              selectInput(ns("group_by"), "Group cells by",
                choices = c("seurat_clusters"), selected = "seurat_clusters")
            )
          ),

          actionButton(ns("run_scoring"), "Run Gene Set Scoring",
            class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),

    fluidRow(
      column(12, uiOutput(ns("gs_results_ui")))
    )
  )
}

geneSetScoringServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    scored_sets <- reactiveVal(list())

    # Update group_by choices when dataset changes
    observe({
      req(values$active_dataset)
      md_cols <- colnames(values$active_dataset@meta.data)
      factor_cols <- md_cols[sapply(values$active_dataset@meta.data, function(x)
        is.factor(x) || is.character(x))]
      updateSelectInput(session, "group_by", choices = factor_cols,
        selected = if ("seurat_clusters" %in% factor_cols) "seurat_clusters" else factor_cols[1])
    })

    # Add UCell if available
    observe({
      choices <- c("AddModuleScore (Seurat)" = "seurat")
      if (requireNamespace("UCell", quietly = TRUE)) {
        choices <- c(choices, "UCell" = "ucell")
      }
      updateSelectInput(session, "scoring_method", choices = choices)
    })

    # MSigDB pathway selector
    output$msigdb_pathway_selector <- renderUI({
      req(input$gs_source != "custom")
      ns <- session$ns

      if (!requireNamespace("msigdbr", quietly = TRUE)) {
        return(p("msigdbr package required. Install with: install.packages('msigdbr')",
                 style = "color: red;"))
      }

      category <- switch(input$gs_source,
        "hallmark" = "H",
        "c2" = "C2",
        "c5_bp" = "C5",
        "c7" = "C7"
      )

      subcategory <- if (input$gs_source == "c5_bp") "GO:BP" else NULL

      tryCatch({
        gene_sets <- msigdbr::msigdbr(species = input$gs_species, category = category,
                                       subcategory = subcategory)
        pathway_names <- sort(unique(gene_sets$gs_name))

        selectizeInput(ns("selected_pathways"), "Select gene sets (max 10)",
          choices = pathway_names, selected = head(pathway_names, 1),
          multiple = TRUE, options = list(maxItems = 10,
            placeholder = "Type to search pathways..."))
      }, error = function(e) {
        p(paste("Failed to load MSigDB gene sets:", e$message), style = "color: red;")
      })
    })

    # Run scoring
    observeEvent(input$run_scoring, {
      req(values$active_dataset)

      tryCatch({
        withProgress(message = "Scoring gene sets...", value = 0, {
          obj <- values$active_dataset
          gene_sets_list <- list()

          if (input$gs_source == "custom") {
            raw <- input$custom_genes
            if (is.null(raw) || nchar(trimws(raw)) == 0) {
              showNotification("Please enter gene names", type = "error")
              return()
            }
            genes <- trimws(unlist(strsplit(raw, "[,\n\r]+")))
            genes <- genes[nchar(genes) > 0]
            present <- intersect(genes, rownames(obj))

            if (length(present) < 2) {
              showNotification(
                paste0("Only ", length(present), " of ", length(genes),
                       " genes found in dataset. Need at least 2."),
                type = "error")
              return()
            }

            set_name <- trimws(input$custom_set_name)
            if (nchar(set_name) == 0) set_name <- "CustomScore"
            gene_sets_list[[set_name]] <- present

          } else {
            # MSigDB gene sets
            req(input$selected_pathways)
            if (!requireNamespace("msigdbr", quietly = TRUE)) {
              showNotification("msigdbr package required", type = "error")
              return()
            }

            category <- switch(input$gs_source,
              "hallmark" = "H", "c2" = "C2", "c5_bp" = "C5", "c7" = "C7"
            )
            subcategory <- if (input$gs_source == "c5_bp") "GO:BP" else NULL

            all_sets <- msigdbr::msigdbr(species = input$gs_species,
                                          category = category,
                                          subcategory = subcategory)

            for (pw in input$selected_pathways) {
              pw_genes <- all_sets$gene_symbol[all_sets$gs_name == pw]
              present <- intersect(pw_genes, rownames(obj))
              if (length(present) >= 2) {
                gene_sets_list[[pw]] <- present
              }
            }

            if (length(gene_sets_list) == 0) {
              showNotification("No gene sets had >= 2 matching genes in dataset", type = "error")
              return()
            }
          }

          incProgress(0.3, detail = paste("Scoring", length(gene_sets_list), "gene sets..."))

          # Score each gene set
          score_names <- c()
          for (i in seq_along(gene_sets_list)) {
            set_name <- names(gene_sets_list)[i]
            genes <- gene_sets_list[[set_name]]

            # Clean name for column
            col_name <- gsub("[^A-Za-z0-9_]", "_", set_name)
            col_name <- paste0("GS_", col_name)

            if (input$scoring_method == "ucell" && requireNamespace("UCell", quietly = TRUE)) {
              features_list <- list(genes)
              names(features_list) <- col_name
              set.seed(42)
              obj <- UCell::AddModuleScore_UCell(obj, features = features_list)
            } else {
              # Avoid name collision: if col_name already exists, make unique
              existing_cols <- colnames(obj@meta.data)
              safe_name <- col_name
              if (safe_name %in% existing_cols) {
                safe_name <- make.unique(c(existing_cols, safe_name), sep = "_")
                safe_name <- safe_name[length(safe_name)]
              }
              set.seed(42)
              obj <- AddModuleScore(obj, features = list(genes),
                                     ctrl = input$ctrl_genes, name = safe_name)
              # Seurat appends "1" to name
              actual_col <- paste0(safe_name, "1")
              if (actual_col %in% colnames(obj@meta.data)) {
                obj@meta.data[[col_name]] <- obj@meta.data[[actual_col]]
                obj@meta.data[[actual_col]] <- NULL
              }
            }

            score_names <- c(score_names, col_name)
            incProgress(0.3 + 0.5 * (i / length(gene_sets_list)))
          }

          obj@misc$module_score_seed <- 42
          values$active_dataset <- obj

          scored_sets(list(
            score_names = score_names,
            gene_sets = gene_sets_list,
            method = input$scoring_method,
            group_by = input$group_by
          ))

          incProgress(1, detail = "Done!")
          showNotification(
            paste("Scored", length(score_names), "gene set(s). Columns added to metadata."),
            type = "message")
        })
      }, error = function(e) {
        showNotification(paste("Gene set scoring failed:", e$message), type = "error")
      })
    })

    # Results UI
    output$gs_results_ui <- renderUI({
      req(scored_sets(), length(scored_sets()$score_names) > 0)
      ns <- session$ns

      score_names <- scored_sets()$score_names
      wellPanel(
        h4("Gene Set Scoring Results"),
        fluidRow(
          column(4,
            selectInput(ns("display_score"), "Select score to display",
              choices = score_names, selected = score_names[1])
          ),
          column(4,
            selectInput(ns("plot_type"), "Plot type",
              choices = c("Violin" = "violin", "Feature (UMAP)" = "feature",
                          "Ridge" = "ridge", "Box" = "box"),
              selected = "violin")
          ),
          column(4,
            selectInput(ns("display_group_by"), "Group by",
              choices = colnames(values$active_dataset@meta.data),
              selected = scored_sets()$group_by)
          )
        ),
        plotlyOutput(ns("gs_plot"), height = "500px"),
        br(),
        h5("Score Summary per Group"),
        DT::dataTableOutput(ns("gs_summary_table"))
      )
    })

    # Main plot
    output$gs_plot <- renderPlotly({
      req(values$active_dataset, input$display_score)
      obj <- values$active_dataset
      score_col <- input$display_score
      req(score_col %in% colnames(obj@meta.data))

      group_col <- input$display_group_by %||% "seurat_clusters"
      if (!group_col %in% colnames(obj@meta.data)) group_col <- "seurat_clusters"

      md <- obj@meta.data
      md$group <- as.factor(md[[group_col]])
      md$score <- md[[score_col]]

      if (input$plot_type == "violin") {
        p <- ggplot(md, aes(x = group, y = score, fill = group)) +
          geom_violin(scale = "width", trim = TRUE) +
          geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.5) +
          theme_minimal() +
          labs(title = score_col, x = group_col, y = "Module Score") +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
        ggplotly(p, tooltip = c("x", "y"))

      } else if (input$plot_type == "feature") {
        red <- find_best_reduction(obj)
        if (is.null(red)) {
          return(plotly_empty() %>% layout(title = "No UMAP/tSNE reduction found"))
        }
        emb <- Embeddings(obj, reduction = red)
        md$dim1 <- emb[, 1]
        md$dim2 <- emb[, 2]
        plot_ly(md, x = ~dim1, y = ~dim2, color = ~score,
                type = "scatter", mode = "markers",
                marker = list(size = 3, opacity = 0.6),
                colors = viridis::viridis(100)) %>%
          layout(title = paste(score_col, "on", toupper(red)),
                 xaxis = list(title = paste0(toupper(red), "_1")),
                 yaxis = list(title = paste0(toupper(red), "_2")))

      } else if (input$plot_type == "ridge") {
        p <- ggplot(md, aes(x = score, y = group, fill = group)) +
          ggridges::geom_density_ridges(scale = 1.2, alpha = 0.7) +
          theme_minimal() +
          labs(title = score_col, x = "Module Score", y = group_col) +
          theme(legend.position = "none")
        tryCatch(
          ggplotly(p),
          error = function(e) {
            # ggridges may not convert to plotly cleanly — fall back to violin
            p2 <- ggplot(md, aes(x = group, y = score, fill = group)) +
              geom_violin(scale = "width") +
              theme_minimal() + labs(title = score_col) +
              theme(legend.position = "none")
            ggplotly(p2)
          }
        )

      } else {
        # Box plot
        p <- ggplot(md, aes(x = group, y = score, fill = group)) +
          geom_boxplot(outlier.size = 0.5) +
          theme_minimal() +
          labs(title = score_col, x = group_col, y = "Module Score") +
          theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
        ggplotly(p, tooltip = c("x", "y"))
      }
    })

    # Summary table
    output$gs_summary_table <- DT::renderDataTable({
      req(values$active_dataset, input$display_score)
      obj <- values$active_dataset
      score_col <- input$display_score
      group_col <- input$display_group_by %||% "seurat_clusters"
      req(score_col %in% colnames(obj@meta.data))
      if (!group_col %in% colnames(obj@meta.data)) group_col <- "seurat_clusters"

      md <- obj@meta.data
      summary_df <- md %>%
        dplyr::group_by(Group = .data[[group_col]]) %>%
        dplyr::summarise(
          N_cells = dplyr::n(),
          Mean = round(mean(.data[[score_col]], na.rm = TRUE), 4),
          Median = round(median(.data[[score_col]], na.rm = TRUE), 4),
          SD = round(sd(.data[[score_col]], na.rm = TRUE), 4),
          Min = round(min(.data[[score_col]], na.rm = TRUE), 4),
          Max = round(max(.data[[score_col]], na.rm = TRUE), 4),
          .groups = "drop"
        )

      # Compute pairwise Wilcoxon p-values (each cluster vs all others)
      # then apply global BH correction across ALL comparisons
      groups <- unique(md[[group_col]])
      if (length(groups) >= 2) {
        all_pvals <- c()
        pval_labels <- c()
        for (g in groups) {
          in_group <- md[[score_col]][md[[group_col]] == g]
          out_group <- md[[score_col]][md[[group_col]] != g]
          pv <- tryCatch(
            wilcox.test(in_group, out_group)$p.value,
            error = function(e) NA_real_
          )
          all_pvals <- c(all_pvals, pv)
          pval_labels <- c(pval_labels, as.character(g))
        }
        # Global BH correction across all cluster comparisons
        adj_pvals <- p.adjust(all_pvals, method = "BH")
        pval_df <- data.frame(
          Group = pval_labels,
          P_value = signif(all_pvals, 4),
          P_adj_BH = signif(adj_pvals, 4),
          stringsAsFactors = FALSE
        )
        summary_df <- dplyr::left_join(summary_df,
          pval_df, by = "Group")
      }

      DT::datatable(summary_df, options = list(pageLength = 20, dom = 'tip'),
                     rownames = FALSE)
    })
  })
}
