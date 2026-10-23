# ──────────────────────────────────────────────────────────────────────────────
# downstream_ai.R — AI Analysis tab + Claude Tool Use (function-calling) agent
# Part of the downstream/ module family.  Sourced by downstream_main.R.
#
# Dependencies (loaded before this file):
#   ai_assistant.R  → .mdtohtml(), call_claude_api(), call_ai(),
#                      ai_thread_message/loading/error, AI_EXPERT_PERSONA
#   downstream_enrichment.R → .run_pathway_enrichment_auto(), create_ora_dot_plot(), create_gsea_plot()
#   downstream_ccc.R        → .run_ccc_auto(), .plot_ccc_chord_auto(), .plot_ccc_heatmap_auto()
#   downstream_tf.R         → .run_tf_activity_auto(), .run_grn_auto(),
#                              .plot_tf_heatmap_auto(), .plot_grn_network_auto()
# ──────────────────────────────────────────────────────────────────────────────

# ── Helper: status badge (used only in AI module) ────────────────────────────
.make_badge <- function(label, is_ready, detail) {
  cls <- if (is_ready) "ai-status-badge ai-badge-ready" else "ai-status-badge ai-badge-empty"
  icon_name <- if (is_ready) "check-circle" else "circle"
  tags$span(class = cls,
    icon(icon_name, style = paste0("font-size:10px;color:",
                                    if (is_ready) "#4caf50" else "#ff9800", ";")),
    paste0(" ", label, ": ", detail)
  )
}

# ==============================================================================
# Claude Tool Use (function-calling) architecture
# ==============================================================================

# ── DOWNSTREAM_TOOLS: JSON-schema definitions for Claude's tools parameter ───
DOWNSTREAM_TOOLS <- list(
  list(
    name = "get_dataset_summary",
    description = paste0(
      "Get an overview of the loaded scRNA-seq dataset including cell count, ",
      "gene count, available metadata columns, clustering info, and reductions."),
    input_schema = list(type = "object", properties = list(), required = list())
  ),
  list(
    name = "run_differential_expression",
    description = paste0(
      "Run differential expression analysis between two groups of cells. ",
      "Returns top significant genes with log2FC and adjusted p-values."),
    input_schema = list(
      type = "object",
      properties = list(
        group_column = list(type = "string",
          description = "Metadata column containing group labels (e.g. 'seurat_clusters', 'cell_type')"),
        group1 = list(type = "string",
          description = "Name of the first group (ident.1)"),
        group2 = list(type = "string",
          description = "Name of the second group (ident.2). Use 'all_others' to compare group1 vs all remaining cells."),
        test_method = list(type = "string",
          enum = list("wilcox", "MAST", "DESeq2", "t"),
          description = "Statistical test to use. Default: wilcox"),
        max_genes = list(type = "integer",
          description = "Maximum number of top genes to return. Default: 50")
      ),
      required = list("group_column", "group1")
    )
  ),
  list(
    name = "run_pathway_enrichment",
    description = paste0(
      "Run pathway enrichment analysis (ORA and/or GSEA) on DE results. ",
      "Requires DE results to be available."),
    input_schema = list(
      type = "object",
      properties = list(
        species = list(type = "string",
          enum = list("Homo sapiens", "Mus musculus"),
          description = "Species for gene set database. Default: Homo sapiens"),
        database = list(type = "string",
          enum = list("H", "KEGG", "REACTOME", "C5", "WP"),
          description = "Gene set database. H=Hallmark, C5=GO BP. Default: H"),
        method = list(type = "string",
          enum = list("ORA", "GSEA", "BOTH"),
          description = "Enrichment method. Default: BOTH")
      ),
      required = list()
    )
  ),
  list(
    name = "run_ccc_analysis",
    description = paste0(
      "Run cell-cell communication analysis to identify significant ",
      "ligand-receptor interactions between cell populations."),
    input_schema = list(
      type = "object",
      properties = list(
        group_column = list(type = "string",
          description = "Metadata column for cell type/cluster labels"),
        species = list(type = "string",
          enum = list("human", "mouse"),
          description = "Species. Default: human"),
        min_score = list(type = "number",
          description = "Minimum interaction score threshold. Default: 0.1")
      ),
      required = list("group_column")
    )
  ),
  list(
    name = "run_tf_activity",
    description = "Score transcription factor activity across cell clusters using regulon databases.",
    input_schema = list(
      type = "object",
      properties = list(
        group_column = list(type = "string",
          description = "Metadata column for cell clusters"),
        species = list(type = "string",
          enum = list("human", "mouse"),
          description = "Species. Default: human")
      ),
      required = list("group_column")
    )
  ),
  list(
    name = "run_grn_analysis",
    description = paste0(
      "Build a gene regulatory network from highly variable genes using ",
      "correlation-based methods with community detection."),
    input_schema = list(
      type = "object",
      properties = list(
        group_column = list(type = "string",
          description = "Metadata column for cell clusters"),
        n_genes = list(type = "integer",
          description = "Number of top variable genes to use. Default: 200"),
        cor_threshold = list(type = "number",
          description = "Correlation threshold for edge inclusion. Default: 0.3")
      ),
      required = list("group_column")
    )
  ),
  list(
    name = "get_gene_expression",
    description = "Get expression statistics for specific genes across clusters or conditions.",
    input_schema = list(
      type = "object",
      properties = list(
        genes = list(type = "array", items = list(type = "string"),
          description = "List of gene names to query"),
        group_column = list(type = "string",
          description = "Metadata column for grouping. Default: seurat_clusters")
      ),
      required = list("genes")
    )
  ),
  list(
    name = "summarise_available_results",
    description = paste0(
      "List all analysis results currently available in the session ",
      "(DE, enrichment, CCC, TF, GRN) with summary statistics."),
    input_schema = list(type = "object", properties = list(), required = list())
  ),
  list(
    name = "get_metadata_columns",
    description = paste0(
      "List all metadata columns with their types and example values. ",
      "Useful for discovering available grouping variables, annotations, and scores."),
    input_schema = list(type = "object", properties = list(), required = list())
  ),
  list(
    name = "run_gene_set_scoring",
    description = paste0(
      "Score cells for a custom gene set using Seurat's AddModuleScore. ",
      "Adds a score column to metadata and returns per-group summary statistics."),
    input_schema = list(
      type = "object",
      properties = list(
        genes = list(type = "array", items = list(type = "string"),
          description = "Gene names to score"),
        set_name = list(type = "string",
          description = "Name for the gene set (used as column name). Default: CustomScore"),
        group_column = list(type = "string",
          description = "Metadata column for per-group summary. Default: seurat_clusters")
      ),
      required = list("genes")
    )
  ),
  list(
    name = "get_umap_data",
    description = paste0(
      "Get UMAP (or tSNE/PCA) embedding coordinates for all cells along with ",
      "a selected metadata column. Useful for understanding spatial distribution."),
    input_schema = list(
      type = "object",
      properties = list(
        color_by = list(type = "string",
          description = "Metadata column to include. Default: seurat_clusters"),
        n_sample = list(type = "integer",
          description = "Number of cells to sample (for large datasets). Default: 500")
      ),
      required = list()
    )
  )
)

# ── callClaudeWithTools: API call with tools parameter ────────────────────────
callClaudeWithTools <- function(api_key, model, messages, tools,
                                max_tokens = 4096L, timeout = 120,
                                system_prompt = NULL) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    return(list(ok = FALSE, stop_reason = "error", content = list(),
                text = NULL, tool_calls = list(),
                tokens_in = NA, tokens_out = NA,
                error = "httr2 package required"))
  }

  tryCatch({
    # Handle thinking models: strip suffix, add thinking config, enforce min tokens
    is_thinking <- grepl("thinking", model, fixed = TRUE)
    api_model <- if (is_thinking) sub("-thinking$", "", model) else model

    if (is_thinking) {
      max_tokens <- max(max_tokens, 16384L)
    }

    body <- list(
      model    = api_model,
      max_tokens = as.integer(max_tokens),
      messages = messages,
      tools    = tools
    )

    if (is_thinking) {
      body$thinking <- list(
        type = "enabled",
        budget_tokens = as.integer(max_tokens * 0.6)
      )
    }

    if (!is.null(system_prompt) && nzchar(system_prompt)) {
      body$system <- system_prompt
    }

    req <- httr2::request("https://api.anthropic.com/v1/messages") |>
      httr2::req_headers(
        "x-api-key"         = api_key,
        "anthropic-version" = "2023-06-01",
        "content-type"      = "application/json"
      ) |>
      httr2::req_body_json(body) |>
      httr2::req_timeout(timeout)

    # Add beta header for extended thinking models
    if (is_thinking) {
      req <- req |> httr2::req_headers("anthropic-beta" = "interleaved-thinking-2025-05-14")
    }

    # Retry loop with exponential backoff for transient errors
    max_retries <- 3L
    for (attempt in seq_len(max_retries + 1L)) {
      resp_result <- tryCatch({
        resp <- httr2::req_perform(req)
        list(ok = TRUE, resp = resp)
      }, error = function(e) {
        retryable <- grepl("429|rate.limit|503|529|overloaded", e$message, ignore.case = TRUE)
        list(ok = FALSE, error = e$message, retryable = retryable)
      })

      if (resp_result$ok || !isTRUE(resp_result$retryable) || attempt > max_retries) break
      Sys.sleep(2^attempt + runif(1))
    }
    if (!resp_result$ok) stop(resp_result$error)
    resp <- resp_result$resp

    resp_body <- httr2::resp_body_json(resp)

    stop_reason <- resp_body$stop_reason %||% "end_turn"
    content     <- resp_body$content %||% list()
    tokens_in   <- tryCatch(resp_body$usage$input_tokens,  error = function(e) NA)
    tokens_out  <- tryCatch(resp_body$usage$output_tokens, error = function(e) NA)

    text_blocks <- list()
    tool_calls  <- list()
    for (block in content) {
      if (identical(block$type, "text")) {
        text_blocks <- c(text_blocks, list(block))
      } else if (identical(block$type, "tool_use")) {
        tool_calls <- c(tool_calls, list(list(
          id    = block$id,
          name  = block$name,
          input = block$input
        )))
      }
    }

    combined_text <- paste(
      sapply(text_blocks, function(b) b$text %||% ""), collapse = "\n")

    list(ok = TRUE, stop_reason = stop_reason, content = content,
         text = combined_text, tool_calls = tool_calls,
         tokens_in = tokens_in, tokens_out = tokens_out, error = NULL)

  }, error = function(e) {
    list(ok = FALSE, stop_reason = "error", content = list(),
         text = NULL, tool_calls = list(),
         tokens_in = NA, tokens_out = NA,
         error = paste("API error:", e$message))
  })
}

# ── executeTool: dispatch tool calls to R functions ───────────────────────────
# Returns a structured list: list(ok = TRUE/FALSE, result = "text", is_error = logical)
# The is_error flag is passed to Claude's tool_result so it knows a tool failed
# and can adapt its plan (e.g., skip enrichment if DE returned 0 genes).
executeTool <- function(tool_name, tool_input, values) {
  tryCatch({
    switch(tool_name,

      "get_dataset_summary" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        meta <- obj@meta.data
        n_cells <- ncol(obj); n_genes <- nrow(obj)
        reds <- tryCatch(names(obj@reductions), error = function(e) character(0))
        cat_cols <- colnames(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
        num_cols <- colnames(meta)[sapply(meta, is.numeric)]
        cluster_info <- ""
        for (col in head(cat_cols, 5)) {
          tbl <- sort(table(meta[[col]]), decreasing = TRUE)
          cluster_info <- paste0(cluster_info, "\n  ", col, " (", length(tbl), " groups): ",
            paste(head(paste0(names(tbl), "=", as.integer(tbl)), 10), collapse = ", "))
        }
        list(ok = TRUE, result = paste0(
          "Dataset: ", n_cells, " cells x ", n_genes, " genes\n",
          "Assay: ", DefaultAssay(obj), "\n",
          "Reductions: ", paste(reds, collapse = ", "), "\n",
          "Categorical columns: ", paste(cat_cols, collapse = ", "), "\n",
          "Numeric columns: ", paste(num_cols, collapse = ", "), "\n",
          "Group distributions:", cluster_info))
      },

      "run_differential_expression" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        group_col <- tool_input$group_column
        group1 <- tool_input$group1
        group2 <- tool_input$group2 %||% NULL
        test_method <- tool_input$test_method %||% "wilcox"
        max_genes <- tool_input$max_genes %||% 50
        if (!group_col %in% colnames(obj@meta.data))
          return(list(ok = FALSE, result = paste0("Column '", group_col, "' not found.")))
        if (identical(group2, "all_others")) group2 <- NULL

        # Use the UNIFIED DE engine — same as manual UI path
        result <- tryCatch(
          run_de_analysis(
            obj = obj,
            group_col = group_col,
            ident1 = group1,
            ident2 = group2,
            test_method = test_method,
            analysis_type = "singlecell",
            seed = 42L,
            min_pct_prefilter = 0.01,
            min_pct_display = 0,
            logfc_threshold = 0,
            only_pos = FALSE
          ),
          error = function(e) list(data = NULL, error = e$message)
        )

        if (!is.null(result$error) || is.null(result$data)) {
          return(list(ok = FALSE, result = paste("DE failed:", result$error %||% "No results")))
        }

        de <- result$data
        sig <- de[!is.na(de$p_val_adj) & de$p_val_adj < 0.05, ]
        sig <- sig[order(-abs(sig$avg_log2FC)), ]
        top <- head(sig, max_genes)
        values$last_de_results <- de
        list(ok = TRUE, result = paste0(
          "DE: ", group1, " vs ", if (is.null(group2)) "all others" else group2, "\n",
          "Test: ", test_method, " | Total: ", nrow(de), " | Sig: ", nrow(sig), "\n\n",
          paste(paste0("  ", top$gene, " | log2FC=", round(top$avg_log2FC, 3),
                       " | padj=", signif(top$p_val_adj, 3)), collapse = "\n")))
      },

      "run_pathway_enrichment" = {
        de <- values$last_de_results
        if (is.null(de) || !is.data.frame(de) || nrow(de) == 0)
          return(list(ok = FALSE, result = "No DE results. Run DE first."))
        species  <- tool_input$species  %||% "Homo sapiens"
        database <- tool_input$database %||% "H"
        method   <- tool_input$method   %||% "BOTH"
        pw_res <- tryCatch(
          .run_pathway_enrichment_auto(de, species = species,
                                       database = database, method = method),
          error = function(e) list(ok = FALSE, result = paste("Enrichment failed:", e$message)))
        if (is.list(pw_res) && !is.null(pw_res$ok)) return(pw_res)
        if (is.null(pw_res)) return(list(ok = FALSE, result = "No enrichment results."))
        values$last_enrich_results <- pw_res
        parts <- c()
        if (!is.null(pw_res$ora$up) && nrow(pw_res$ora$up) > 0) {
          top <- head(pw_res$ora$up, 15)
          parts <- c(parts, paste0("ORA Up (", nrow(pw_res$ora$up), "):\n",
            paste(paste0("  ", seq_len(nrow(top)), ". ", as.character(top$Description),
                         " (padj=", signif(as.numeric(as.character(top$p.adjust)), 3), ")"),
                  collapse = "\n")))
        }
        if (!is.null(pw_res$ora$down) && nrow(pw_res$ora$down) > 0) {
          top <- head(pw_res$ora$down, 15)
          parts <- c(parts, paste0("ORA Down (", nrow(pw_res$ora$down), "):\n",
            paste(paste0("  ", seq_len(nrow(top)), ". ", as.character(top$Description),
                         " (padj=", signif(as.numeric(as.character(top$p.adjust)), 3), ")"),
                  collapse = "\n")))
        }
        if (!is.null(pw_res$gsea) && nrow(pw_res$gsea) > 0) {
          top <- head(pw_res$gsea[order(-abs(pw_res$gsea$NES)), ], 15)
          parts <- c(parts, paste0("GSEA (", nrow(pw_res$gsea), "):\n",
            paste(paste0("  ", seq_len(nrow(top)), ". ", as.character(top$pathway),
                         " (NES=", round(top$NES, 2), ", padj=", signif(top$padj, 3), ")"),
                  collapse = "\n")))
        }
        list(ok = TRUE, result = paste(parts, collapse = "\n\n"))
      },

      "run_ccc_analysis" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        group_col <- tool_input$group_column
        species   <- tool_input$species   %||% "human"
        min_score <- tool_input$min_score  %||% 0.1
        if (!group_col %in% colnames(obj@meta.data))
          return(list(ok = FALSE, result = paste0("Column '", group_col, "' not found.")))
        cc_res <- tryCatch(
          .run_ccc_auto(obj, group_col, species = species,
                        min_expr = 0.1, nperm = 50, min_score = min_score),
          error = function(e) list(ok = FALSE, result = paste("CCC failed:", e$message)))
        if (is.list(cc_res) && !is.null(cc_res$ok)) return(cc_res)
        if (is.null(cc_res) || !is.data.frame(cc_res) || nrow(cc_res) == 0)
          return(list(ok = FALSE, result = "No significant interactions found."))
        values$last_cc_results <- cc_res
        sig <- cc_res[cc_res$pval_adj < 0.05, ]
        top <- head(sig[order(sig$pval_adj), ], 25)
        pw_counts <- sort(table(sig$pathway), decreasing = TRUE)
        list(ok = TRUE, result = paste0(
          "CCC: ", nrow(cc_res), " total, ", nrow(sig), " significant\n",
          "Top pathways: ", paste(head(paste0(names(pw_counts), " (", as.integer(pw_counts), ")"), 8), collapse = ", "), "\n\n",
          paste(paste0("  ", top$source, "->", top$target, ": ", top$ligand, "-", top$receptor,
                       " (", top$pathway, ", score=", round(top$score, 3), ")"), collapse = "\n")))
      },

      "run_tf_activity" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        group_col <- tool_input$group_column
        species   <- tool_input$species %||% "human"
        if (!group_col %in% colnames(obj@meta.data))
          return(list(ok = FALSE, result = paste0("Column '", group_col, "' not found.")))
        tf_res <- tryCatch(
          .run_tf_activity_auto(obj, group_col, species = species,
                                min_targets = 5, score_method = "wmean", max_cells = 5000),
          error = function(e) list(ok = FALSE, result = paste("TF failed:", e$message)))
        if (is.list(tf_res) && !is.null(tf_res$ok)) return(tf_res)
        if (is.null(tf_res)) return(list(ok = FALSE, result = "No TF results."))
        values$last_tf_results <- tf_res$results
        df <- tf_res$results
        top <- head(df[order(-abs(df$Activity)), ], 20)
        list(ok = TRUE, result = paste0(
          "TF Activity: ", length(unique(df$TF)), " TFs, ", length(unique(df$Cluster)), " clusters\n\n",
          paste(paste0("  ", top$TF, " in ", top$Cluster, ": ", round(top$Activity, 3)), collapse = "\n")))
      },

      "run_grn_analysis" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        group_col  <- tool_input$group_column
        n_genes    <- tool_input$n_genes       %||% 200
        cor_thresh <- tool_input$cor_threshold  %||% 0.3
        if (!group_col %in% colnames(obj@meta.data))
          return(list(ok = FALSE, result = paste0("Column '", group_col, "' not found.")))
        grn_res <- tryCatch(
          .run_grn_auto(obj, group_col, n_genes = n_genes,
                        cor_method = "spearman", cor_thresh = cor_thresh,
                        min_degree = 2, community = "louvain"),
          error = function(e) list(ok = FALSE, result = paste("GRN failed:", e$message)))
        if (is.list(grn_res) && !is.null(grn_res$ok)) return(grn_res)
        if (is.null(grn_res)) return(list(ok = FALSE, result = "No GRN results."))
        mods <- grn_res$modules
        list(ok = TRUE, result = paste0(
          "GRN: ", igraph::vcount(grn_res$graph), " genes, ",
          igraph::ecount(grn_res$graph), " edges, ", nrow(mods), " modules\n\n",
          paste(paste0("  Module ", mods$Module, ": hub=", mods$Hub_Gene,
                       " (", mods$N_Genes, " genes: ", mods$Top_Genes, ")"), collapse = "\n")))
      },

      "get_gene_expression" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        genes     <- tool_input$genes
        group_col <- tool_input$group_column %||% "seurat_clusters"
        if (!group_col %in% colnames(obj@meta.data))
          return(list(ok = FALSE, result = paste0("Column '", group_col, "' not found.")))
        avail  <- intersect(genes, rownames(obj))
        missing <- setdiff(genes, rownames(obj))
        if (length(avail) == 0)
          return(list(ok = FALSE, result = paste0("None found. Missing: ", paste(missing, collapse = ", "))))
        expr_data <- tryCatch(Seurat::FetchData(obj, vars = c(avail, group_col)), error = function(e) NULL)
        if (is.null(expr_data)) return(list(ok = FALSE, result = "Failed to fetch expression."))
        parts <- c()
        for (gene in avail) {
          stats <- tapply(expr_data[[gene]], expr_data[[group_col]], function(x) {
            x <- x[is.finite(x)]
            c(mean = round(mean(x), 3), pct = round(mean(x > 0) * 100, 1))
          })
          stats_text <- paste(sapply(names(stats), function(g) {
            paste0("    ", g, ": mean=", stats[[g]]["mean"], ", %expr=", stats[[g]]["pct"])
          }), collapse = "\n")
          parts <- c(parts, paste0("  ", gene, ":\n", stats_text))
        }
        list(ok = TRUE, result = paste0("Expression across ", group_col, ":\n",
          paste(parts, collapse = "\n\n"),
          if (length(missing) > 0) paste0("\n\nNot found: ", paste(missing, collapse = ", ")) else ""))
      },

      "get_metadata_columns" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        md <- obj@meta.data
        parts <- c()
        for (col in colnames(md)) {
          vals <- md[[col]]
          if (is.factor(vals) || is.character(vals)) {
            tbl <- sort(table(vals), decreasing = TRUE)
            parts <- c(parts, paste0("  ", col, " [categorical, ", length(tbl), " levels]: ",
              paste(head(paste0(names(tbl), "(", as.integer(tbl), ")"), 6), collapse = ", ")))
          } else if (is.numeric(vals)) {
            parts <- c(parts, paste0("  ", col, " [numeric]: range=",
              round(min(vals, na.rm = TRUE), 3), " to ", round(max(vals, na.rm = TRUE), 3),
              ", mean=", round(mean(vals, na.rm = TRUE), 3)))
          } else {
            parts <- c(parts, paste0("  ", col, " [", class(vals)[1], "]"))
          }
        }
        list(ok = TRUE, result = paste0("Metadata columns (", ncol(md), "):\n",
          paste(parts, collapse = "\n")))
      },

      "run_gene_set_scoring" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        genes <- tool_input$genes
        set_name <- tool_input$set_name %||% "CustomScore"
        group_col <- tool_input$group_column %||% "seurat_clusters"

        present <- intersect(genes, rownames(obj))
        if (length(present) < 2)
          return(list(ok = FALSE, result = paste0("Only ", length(present),
            " of ", length(genes), " genes found. Need >= 2.")))

        col_name <- paste0("GS_", gsub("[^A-Za-z0-9_]", "_", set_name))
        obj <- AddModuleScore(obj, features = list(present), ctrl = 100, name = col_name)
        actual_col <- paste0(col_name, "1")
        if (actual_col %in% colnames(obj@meta.data)) {
          obj@meta.data[[col_name]] <- obj@meta.data[[actual_col]]
          obj@meta.data[[actual_col]] <- NULL
        }
        values$active_dataset <- obj

        if (!group_col %in% colnames(obj@meta.data)) group_col <- "seurat_clusters"
        stats <- tapply(obj@meta.data[[col_name]], obj@meta.data[[group_col]], function(x) {
          c(mean = round(mean(x, na.rm = TRUE), 4), sd = round(sd(x, na.rm = TRUE), 4))
        })
        stats_text <- paste(sapply(names(stats), function(g) {
          paste0("  ", g, ": mean=", stats[[g]]["mean"], ", sd=", stats[[g]]["sd"])
        }), collapse = "\n")
        list(ok = TRUE, result = paste0(
          "Scored ", length(present), " genes as '", col_name, "'\n",
          "Missing: ", paste(setdiff(genes, rownames(obj)), collapse = ", "), "\n\n",
          "Scores by ", group_col, ":\n", stats_text))
      },

      "get_umap_data" = {
        obj <- values$active_dataset
        if (is.null(obj)) return(list(ok = FALSE, result = "No dataset loaded."))
        color_by <- tool_input$color_by %||% "seurat_clusters"
        n_sample <- tool_input$n_sample %||% 500L

        red <- find_best_reduction(obj)
        if (is.null(red)) return(list(ok = FALSE, result = "No UMAP/tSNE/PCA reduction found."))
        emb <- Embeddings(obj, reduction = red)
        if (!color_by %in% colnames(obj@meta.data)) color_by <- "seurat_clusters"

        # Stratified sampling by cluster to ensure small clusters are represented
        n <- min(n_sample, nrow(emb))
        set.seed(42)
        cluster_col <- NULL
        for (cc in c("seurat_clusters", "cluster", color_by)) {
          if (cc %in% colnames(obj@meta.data)) { cluster_col <- cc; break }
        }
        if (!is.null(cluster_col)) {
          clusters <- obj@meta.data[[cluster_col]]
          idx <- unlist(lapply(unique(clusters), function(cl) {
            cl_idx <- which(clusters == cl)
            sample(cl_idx, min(ceiling(n * length(cl_idx) / nrow(emb)), length(cl_idx)))
          }))
          if (length(idx) > n) idx <- sample(idx, n)
        } else {
          idx <- sample(seq_len(nrow(emb)), n)
        }
        df <- data.frame(
          cell = rownames(emb)[idx],
          dim1 = round(emb[idx, 1], 3),
          dim2 = round(emb[idx, 2], 3),
          group = obj@meta.data[[color_by]][idx]
        )
        centroids <- aggregate(cbind(dim1, dim2) ~ group, data = df, FUN = mean)
        list(ok = TRUE, result = paste0(
          "Embedding: ", toupper(red), " (", nrow(emb), " cells, sampled ", n, ")\n",
          "Colored by: ", color_by, "\n\n",
          "Group centroids:\n",
          paste(paste0("  ", centroids$group, ": (", round(centroids$dim1, 2), ", ",
                       round(centroids$dim2, 2), ")"), collapse = "\n")))
      },

      "summarise_available_results" = {
        parts <- c()
        has_data <- !is.null(values$active_dataset)
        parts <- c(parts, paste0("Dataset: ",
          if (has_data) paste0(ncol(values$active_dataset), " cells x ", nrow(values$active_dataset), " genes")
          else "Not loaded"))
        has_de <- !is.null(values$last_de_results) && is.data.frame(values$last_de_results)
        if (has_de) {
          sig <- sum(values$last_de_results$p_val_adj < 0.05, na.rm = TRUE)
          parts <- c(parts, paste0("DE: ", nrow(values$last_de_results), " genes, ", sig, " significant"))
        } else parts <- c(parts, "DE: Not run")
        has_pw <- !is.null(values$last_enrich_results)
        if (has_pw) {
          n_ora <- 0
          if (!is.null(values$last_enrich_results$ora$up)) n_ora <- n_ora + nrow(values$last_enrich_results$ora$up)
          if (!is.null(values$last_enrich_results$ora$down)) n_ora <- n_ora + nrow(values$last_enrich_results$ora$down)
          n_gsea <- if (!is.null(values$last_enrich_results$gsea)) nrow(values$last_enrich_results$gsea) else 0
          parts <- c(parts, paste0("Enrichment: ", n_ora, " ORA + ", n_gsea, " GSEA"))
        } else parts <- c(parts, "Enrichment: Not run")
        has_cc <- !is.null(tryCatch(values$last_cc_results, error = function(e) NULL)) &&
                  is.data.frame(tryCatch(values$last_cc_results, error = function(e) NULL))
        if (has_cc) parts <- c(parts, paste0("CCC: ", nrow(values$last_cc_results), " interactions"))
        else parts <- c(parts, "CCC: Not run")
        has_tf <- !is.null(tryCatch(values$last_tf_results, error = function(e) NULL)) &&
                  is.data.frame(tryCatch(values$last_tf_results, error = function(e) NULL))
        if (has_tf) parts <- c(parts, paste0("TF: ", length(unique(values$last_tf_results$TF)), " TFs"))
        else parts <- c(parts, "TF: Not run")
        list(ok = TRUE, result = paste(parts, collapse = "\n"))
      },

      # default
      list(ok = FALSE, result = paste0("Unknown tool: ", tool_name))
    )
  }, error = function(e) {
    list(ok = FALSE, result = paste("ERROR:", tool_name, "failed:", e$message),
         is_error = TRUE)
  })
}

# ── runAgenticAnalysis: agentic loop with tool calls ──────────────────────────
runAgenticAnalysis <- function(api_key, model, user_prompt, values,
                               max_iterations = 10, max_tokens = 8192L,
                               timeout = 120, progress_fn = NULL) {
  system_prompt <- paste0(
    if (exists("AI_EXPERT_PERSONA")) AI_EXPERT_PERSONA else "",
    "\n\nYou are an AI assistant integrated into a scRNA-seq analysis platform. ",
    "You have access to tools that can run analyses on the user's loaded dataset. ",
    "Use the tools to explore the data and answer the user's question. ",
    "After running analyses, provide a comprehensive biological interpretation.\n\n",
    "Strategy:\n",
    "1. First check what data and results are available (summarise_available_results)\n",
    "2. Run additional analyses as needed to answer the question\n",
    "3. After gathering all needed data, provide your interpretation\n"
  )

  messages <- list(list(role = "user", content = user_prompt))
  total_tokens_in  <- 0
  total_tokens_out <- 0
  tool_log  <- list()
  iteration <- 0
  accumulated_text <- ""   # accumulates the final answer across max_tokens continuations

  while (iteration < max_iterations) {
    iteration <- iteration + 1
    if (!is.null(progress_fn)) {
      progress_fn(iteration / max_iterations,
                  detail = paste0("Agent step ", iteration, "/", max_iterations, "..."))
    }

    result <- callClaudeWithTools(
      api_key = api_key, model = model, messages = messages,
      tools = DOWNSTREAM_TOOLS, max_tokens = max_tokens,
      timeout = timeout, system_prompt = system_prompt
    )

    if (!result$ok) {
      return(list(ok = FALSE, text = NULL, tool_log = tool_log,
                  tokens_in = total_tokens_in, tokens_out = total_tokens_out,
                  error = result$error))
    }

    total_tokens_in  <- total_tokens_in  + (result$tokens_in  %||% 0)
    total_tokens_out <- total_tokens_out + (result$tokens_out %||% 0)

    # Final response — no more tool calls
    if (length(result$tool_calls) == 0) {
      accumulated_text <- paste0(accumulated_text, result$text %||% "")
      # If the answer was cut off by the token limit, continue it instead of
      # returning a truncated response — the user always gets the full answer.
      if (identical(result$stop_reason, "max_tokens") && iteration < max_iterations) {
        messages <- c(messages, list(
          list(role = "assistant", content = result$text),
          list(role = "user", content = paste(
            "Continue the previous answer from exactly where it stopped.",
            "Do not repeat anything already written and do not add a preamble."))
        ))
        next
      }
      return(list(ok = TRUE, text = accumulated_text, tool_log = tool_log,
                  tokens_in = total_tokens_in, tokens_out = total_tokens_out,
                  error = NULL))
    }

    # Append assistant message (with tool_use blocks) to conversation
    messages <- c(messages, list(list(role = "assistant", content = result$content)))

    # Execute each tool call and build tool_result blocks
    tool_result_blocks <- list()
    for (tc in result$tool_calls) {
      tool_output <- executeTool(tc$name, tc$input, values)
      tool_log <- c(tool_log, list(list(
        tool = tc$name, input = tc$input,
        success = isTRUE(tool_output$ok),
        result_preview = substr(tool_output$result, 1, 200)
      )))
      # Truncate large tool results to avoid overwhelming context
      result_text <- tool_output$result
      if (nchar(result_text) > 4000) {
        result_text <- paste0(
          substr(result_text, 1, 3800),
          "\n\n... [truncated — ", nchar(tool_output$result) - 3800, " chars omitted] ..."
        )
      }
      tool_result_block <- list(
        type        = "tool_result",
        tool_use_id = tc$id,
        content     = result_text
      )
      # Signal tool failure to Claude so it can adapt its plan
      if (isTRUE(tool_output$is_error) || !isTRUE(tool_output$ok)) {
        tool_result_block$is_error <- TRUE
      }
      tool_result_blocks <- c(tool_result_blocks, list(tool_result_block))
    }

    # Append tool results as user message
    messages <- c(messages, list(list(role = "user", content = tool_result_blocks)))
  }

  # Max iterations reached
  list(ok = TRUE,
    text = paste0(result$text %||% "",
      "\n\n*[Agent reached maximum iterations (", max_iterations, "). Analysis may be incomplete.]*"),
    tool_log = tool_log,
    tokens_in = total_tokens_in, tokens_out = total_tokens_out,
    error = NULL)
}

# ==============================================================================
# Existing AI Analysis UI + Server (from original downstream.R)
# ==============================================================================
aiAnalysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    # CSS moved to www/ai.css

    # ---- Row 1: API Setup ----
    fluidRow(
      column(12,
        div(class = "ai-section-card",
          fluidRow(
            column(12,
              h4(style = "margin-top:0;", icon("robot"), " AI-Powered Analysis Assistant"),
              p(style = "color: #666; margin-bottom: 12px;",
                "Use Claude AI to interpret your scRNA-seq results. API key and model are set in ",
                tags$b("Data Management > AI Config"), ".")
            )
          ),
          fluidRow(
            column(4,
              selectInput(ns("ai_mode"), "Analysis Mode", width = "100%",
                choices = c(
                  "Interpret Results" = "interpret",
                  "Generate Hypotheses" = "hypotheses",
                  "Full Report" = "report",
                  "Compare Clusters" = "compare",
                  "Methods Section" = "methods",
                  "Figure Legends" = "legends",
                  "Literature Connections" = "literature",
                  "Marker Interpretation" = "marker_interp",
                  "Drug Target Analysis" = "drug_targets",
                  "Custom Question" = "custom"
                ), selected = "interpret")
            ),
            column(4,
              numericInput(ns("ai_max_tokens"), "Max Response Length",
                value = 8192, min = 1024, max = 16384, step = 1024, width = "100%")
            ),
            column(4,
              div(style = "margin-top: 25px; font-size: 13px; color: #666;",
                uiOutput(ns("ai_key_status_badge"))
              )
            )
          )
        )
      )
    ),

    # ---- Row 2: Data Context Status ----
    fluidRow(
      column(12,
        div(class = "ai-section-card",
          fluidRow(
            column(6,
              h5(style = "margin-top:0;", "Available Data Context"),
              uiOutput(ns("ai_context_status"))
            ),
            column(4,
              h5(style = "margin-top:0;", "Include in Prompt"),
              checkboxGroupInput(ns("ai_context"), NULL,
                choices = c(
                  "Dataset overview (cells, genes, clusters)" = "summary",
                  "Differential expression (top DE genes)" = "de",
                  "Pathway enrichment (ORA & GSEA)" = "pathway",
                  "Cell-cell communication (L-R pairs)" = "ccc",
                  "TF activity scores" = "tf",
                  "Cell type annotations" = "annotations",
                  "Cluster markers (top 5 per cluster)" = "markers"
                ), selected = c("summary", "de", "pathway"),
                inline = FALSE)
            ),
            column(2,
              h5(style = "margin-top:0;", "Annotation Column"),
              selectInput(ns("ai_annot_col"), NULL, choices = NULL, width = "100%"),
              helpText("Used for annotation context", style = "font-size: 11px;")
            )
          )
        )
      )
    ),

    # ---- Row 3: Cluster comparison (conditional) ----
    conditionalPanel(
      condition = "input.ai_mode == 'compare'", ns = ns,
      fluidRow(
        column(12,
          div(class = "ai-section-card",
            h5(style = "margin-top:0;", "Cluster Comparison Settings"),
            fluidRow(
              column(4, selectInput(ns("ai_compare_col"), "Grouping Column",
                choices = NULL, width = "100%")),
              column(4, selectInput(ns("ai_compare_g1"), "Group 1",
                choices = NULL, width = "100%")),
              column(4, selectInput(ns("ai_compare_g2"), "Group 2",
                choices = NULL, width = "100%"))
            )
          )
        )
      )
    ),

    # ---- Row 4: Custom query (conditional) ----
    conditionalPanel(
      condition = "input.ai_mode == 'custom'", ns = ns,
      fluidRow(
        column(12,
          div(class = "ai-section-card",
            textAreaInput(ns("ai_custom_query"), "Your Question:",
              placeholder = paste0(
                "Examples:\n",
                "- What biological processes explain the upregulation of inflammatory genes in cluster 3?\n",
                "- Are there signatures of epithelial-mesenchymal transition in my data?\n",
                "- What are the key differences between the T cell and macrophage populations?\n",
                "- Suggest validation experiments for the top DE genes"
              ), rows = 4, width = "100%")
          )
        )
      )
    ),

    # ---- Row 5: Action buttons ----
    fluidRow(
      column(12,
        div(class = "ai-section-card",
          fluidRow(
            column(2, actionButton(ns("ai_run"), "Generate",
              class = "btn-success", width = "100%", icon = icon("paper-plane"))),
            column(2, actionButton(ns("ai_followup_btn"), "Follow Up",
              class = "btn-primary", width = "100%", icon = icon("reply"))),
            column(2, actionButton(ns("ai_preview"), "Preview Prompt",
              class = "btn-default", width = "100%", icon = icon("eye"))),
            column(2, downloadButton(ns("ai_download_html"), "HTML Report",
              class = "btn-info", width = "100%")),
            column(2, downloadButton(ns("ai_download_md"), "Markdown",
              class = "btn-default", width = "100%")),
            column(2, actionButton(ns("ai_clear"), "Clear All",
              class = "btn-warning", width = "100%", icon = icon("trash")))
          ),
          # Follow-up input (hidden until clicked)
          conditionalPanel(
            condition = paste0("input['", ns("ai_followup_btn"), "'] > 0"),
            div(style = "margin-top: 10px;",
              fluidRow(
                column(10,
                  textAreaInput(ns("ai_followup_text"), NULL,
                    placeholder = "Ask a follow-up question about the analysis above...",
                    rows = 2, width = "100%")
                ),
                column(2,
                  actionButton(ns("ai_send_followup"), "Send",
                    class = "btn-primary", width = "100%",
                    style = "margin-top: 25px;", icon = icon("paper-plane"))
                )
              )
            )
          )
        )
      )
    ),

    # ---- Row 6: Prompt preview (collapsible) ----
    uiOutput(ns("ai_prompt_preview_area")),

    # ---- Row 7: Chat conversation ----
    fluidRow(
      column(12,
        uiOutput(ns("ai_chat_area"))
      )
    ),

    # ---- Row 8: Auto-Run All Analyses ----
    fluidRow(column(12,
      div(class = "ai-section-card",
        style = "border-left: 4px solid #4caf50;",
        div(style = "display: flex; align-items: center; margin-bottom: 8px;",
          icon("rocket", style = "font-size: 18px; color: #4caf50; margin-right: 8px;"),
          h5(style = "margin: 0; font-weight: 600;", "Run All Downstream Analyses")
        ),
        p(style = "color: #666; font-size: 12px; margin-bottom: 12px;",
          "Automatically run Pathway Enrichment, Cell-Cell Communication, TF Activity, and Gene Regulatory Network analyses. ",
          "Results are displayed below, added to the Figure Gallery, and included as AI context."),
        fluidRow(
          column(3, selectInput(ns("ai_auto_groupcol"), "Cell Type / Cluster Column", choices = NULL)),
          column(2, selectInput(ns("ai_auto_species"), "Species",
            choices = c("Human" = "human", "Mouse" = "mouse"), selected = "human")),
          column(2, selectInput(ns("ai_auto_pw_db"), "Pathway DB",
            choices = c("Hallmark" = "H", "KEGG" = "KEGG", "Reactome" = "REACTOME",
                        "GO BP" = "C5", "WikiPathways" = "WP"), selected = "H")),
          column(2, actionButton(ns("ai_run_all"), "Run All Analyses",
            class = "btn-success", width = "100%",
            style = "font-weight: 600; padding: 8px 12px; font-size: 14px;",
            icon = icon("play-circle"))),
          column(3, uiOutput(ns("ai_auto_status")))
        )
      )
    )),

    # ---- Row 9: Auto-generated figures (6-panel grid) ----
    fluidRow(column(12,
      uiOutput(ns("ai_auto_figures_ui"))
    )),

    # ---- Row 10: Pathway figures (auto-generated after AI response) ----
    fluidRow(
      column(12,
        div(style = "margin-top: 10px;",
          checkboxInput(ns("show_pathway_figs"), "Show Pathway Figures (from enrichment results)", FALSE)
        ),
        conditionalPanel(
          condition = paste0("input['", ns("show_pathway_figs"), "']"),
          div(class = "ai-section-card",
            h5("Auto-Generated Pathway Figures"),
            fluidRow(
              column(6, plotlyOutput(ns("ai_ora_dot_plot"), height = "450px")),
              column(6, plotlyOutput(ns("ai_gsea_plot"), height = "450px"))
            )
          )
        )
      )
    )
  )
}

aiAnalysisServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Initialize cross-module AI conversation log
    if (is.null(values$ai_conversation_log)) {
      values$ai_conversation_log <- list()
    }

    # Reactive state
    chat_history <- reactiveVal(list())      # list of list(role, content, timestamp, tokens)
    api_messages <- reactiveVal(list())      # messages for API multi-turn
    show_preview <- reactiveVal(FALSE)
    is_processing <- reactiveVal(FALSE)
    ai_auto_results <- reactiveVal(NULL)     # list with pathway, ccc, tf, grn results
    ai_auto_status_info <- reactiveVal(NULL) # list of status per analysis

    # ---- API key status badge ----
    output$ai_key_status_badge <- renderUI({
      key <- values$ai_api_key
      if (!is.null(key) && nzchar(trimws(key))) {
        model <- values$ai_model %||% get_default_model(values$ai_api_key)
        model_label <- if (grepl("opus", model, ignore.case = TRUE)) "Opus"
                       else if (grepl("haiku", model, ignore.case = TRUE)) "Haiku"
                       else "Sonnet"
        tagList(
          icon("check-circle", style = "color: #4caf50;"),
          tags$span(" API key set | Model: ", tags$b(model_label))
        )
      } else {
        tagList(
          icon("exclamation-triangle", style = "color: #ff9800;"),
          tags$span(" No API key — set in ", tags$b("Data Management > AI Config"))
        )
      }
    })

    # ---- Populate auto-run group column selector ----
    observe({
      req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) {
        (is.factor(x) || is.character(x)) && length(unique(x)) >= 2 && length(unique(x)) <= 100
      })]
      annot_like <- grep("cell_type|annotation|label|singler|consensus|seurat_clusters|cluster",
                         cat_cols, value = TRUE, ignore.case = TRUE)
      others <- setdiff(cat_cols, annot_like)
      ordered <- c(sort(annot_like), sort(others))
      sel <- if (length(annot_like) > 0) annot_like[1] else if (length(ordered) > 0) ordered[1] else NULL
      updateSelectInput(session, "ai_auto_groupcol", choices = ordered, selected = sel)
    })

    # ---- Run All Analyses handler ----
    observeEvent(input$ai_run_all, {
      req(values$active_dataset, input$ai_auto_groupcol)

      obj <- values$active_dataset
      group_col <- input$ai_auto_groupcol
      species_map <- c("human" = "Homo sapiens", "mouse" = "Mus musculus")
      species_code <- input$ai_auto_species
      species_full <- species_map[[species_code]]
      pw_db <- input$ai_auto_pw_db

      if (!group_col %in% colnames(obj@meta.data)) {
        showNotification("Selected column not found in metadata", type = "error")
        return()
      }

      # Check for DE results (needed for pathway enrichment)
      has_de <- !is.null(values$last_de_results) &&
                is.data.frame(values$last_de_results) && nrow(values$last_de_results) > 0

      status <- list(pathway = "skipped", ccc = "pending", tf = "pending", grn = "pending")
      results <- list()

      withProgress(message = "Running all downstream analyses...", value = 0, {

        # 1. Pathway Enrichment
        if (has_de) {
          status$pathway <- "running"
          ai_auto_status_info(status)
          incProgress(0.05, detail = "Running pathway enrichment...")
          tryCatch({
            pw_res <- .run_pathway_enrichment_auto(values$last_de_results,
                                                    species = species_full,
                                                    database = pw_db, method = "BOTH")
            if (!is.null(pw_res)) {
              results$pathway <- pw_res
              values$last_enrich_results <- pw_res
              status$pathway <- "done"
            } else {
              status$pathway <- "no results"
            }
          }, error = function(e) {
            status$pathway <<- paste("error:", e$message)
          })
        } else {
          status$pathway <- "skipped (no DE)"
        }
        ai_auto_status_info(status)
        incProgress(0.2, detail = "Pathway done. Running CCC...")

        # 2. Cell-Cell Communication
        status$ccc <- "running"
        ai_auto_status_info(status)
        tryCatch({
          cc_res <- .run_ccc_auto(obj, group_col, species = species_code,
                                   min_expr = 0.1, nperm = 50, min_score = 0.1)
          if (!is.null(cc_res) && is.data.frame(cc_res) && nrow(cc_res) > 0) {
            results$ccc <- cc_res
            values$last_cc_results <- cc_res
            status$ccc <- "done"
          } else {
            status$ccc <- "no results"
          }
        }, error = function(e) {
          status$ccc <<- paste("error:", e$message)
        })
        ai_auto_status_info(status)
        incProgress(0.25, detail = "CCC done. Running TF activity...")

        # 3. TF Activity
        status$tf <- "running"
        ai_auto_status_info(status)
        tryCatch({
          tf_res <- .run_tf_activity_auto(obj, group_col, species = species_code,
                                           min_targets = 5, score_method = "wmean",
                                           max_cells = 5000)
          if (!is.null(tf_res)) {
            results$tf <- tf_res
            values$last_tf_results <- tf_res$results
            status$tf <- "done"
          } else {
            status$tf <- "no results"
          }
        }, error = function(e) {
          status$tf <<- paste("error:", e$message)
        })
        ai_auto_status_info(status)
        incProgress(0.25, detail = "TF done. Running gene network...")

        # 4. Gene Regulatory Network
        status$grn <- "running"
        ai_auto_status_info(status)
        tryCatch({
          grn_res <- .run_grn_auto(obj, group_col, n_genes = 200,
                                    cor_method = "spearman", cor_thresh = 0.3,
                                    min_degree = 2, community = "louvain")
          if (!is.null(grn_res)) {
            results$grn <- grn_res
            status$grn <- "done"
          } else {
            status$grn <- "no results"
          }
        }, error = function(e) {
          status$grn <<- paste("error:", e$message)
        })
        ai_auto_status_info(status)
        incProgress(0.2, detail = "All analyses complete!")
      })

      ai_auto_results(results)

      # Auto-check relevant context checkboxes
      ctx_to_check <- c("summary", "de")
      if (status$pathway == "done") ctx_to_check <- c(ctx_to_check, "pathway")
      if (status$ccc == "done") ctx_to_check <- c(ctx_to_check, "ccc")
      if (status$tf == "done") ctx_to_check <- c(ctx_to_check, "tf")
      ctx_to_check <- c(ctx_to_check, "markers")
      updateCheckboxGroupInput(session, "ai_context", selected = ctx_to_check)

      n_done <- sum(c(status$pathway, status$ccc, status$tf, status$grn) == "done")
      showNotification(
        paste0("Auto-run complete! ", n_done, "/4 analyses succeeded."),
        type = if (n_done >= 3) "message" else "warning", duration = 6)
      gc()
    })

    # ---- Auto-run status badges ----
    output$ai_auto_status <- renderUI({
      status <- ai_auto_status_info()
      if (is.null(status)) {
        return(div(style = "padding-top: 8px; color: #999; font-size: 12px;",
                   icon("info-circle"), " Click 'Run All' to start analyses"))
      }
      badge_fn <- function(label, val) {
        is_err <- grepl("^error:", val)
        color <- if (is_err) "#f44336" else switch(val,
          "done" = "#4caf50", "running" = "#2196f3",
          "pending" = "#bdbdbd", "skipped" = "#ff9800",
          "skipped (no DE)" = "#ff9800",
          "no results" = "#ff9800", "#f44336")
        icon_nm <- if (is_err) "times-circle" else switch(val,
          "done" = "check-circle", "running" = "spinner",
          "pending" = "circle", "skipped" = "forward",
          "skipped (no DE)" = "forward",
          "no results" = "exclamation-triangle", "times-circle")
        spin_class <- if (val == "running") " fa-spin" else ""
        tags$div(style = "display: inline-block; margin: 2px 3px;",
          tags$span(
            style = paste0("display: inline-flex; align-items: center; padding: 3px 8px; ",
                           "border-radius: 12px; font-size: 11px; font-weight: 500; ",
                           "background:", color, "15; color:", color, ";"),
            tags$i(class = paste0("fa fa-", icon_nm, spin_class),
                   style = "font-size: 10px; margin-right: 4px;"),
            label
          )
        )
      }
      div(style = "padding-top: 4px;",
        badge_fn("Pathway", status$pathway),
        badge_fn("CCC", status$ccc),
        badge_fn("TF", status$tf),
        badge_fn("GRN", status$grn)
      )
    })

    # ---- Auto-run figures UI (conditionally rendered) ----
    output$ai_auto_figures_ui <- renderUI({
      res <- ai_auto_results()
      if (is.null(res) || length(res) == 0) return(NULL)
      ns <- session$ns

      has_pw <- !is.null(res$pathway)
      has_ccc <- !is.null(res$ccc)
      has_tf <- !is.null(res$tf)
      has_grn <- !is.null(res$grn)
      if (!has_pw && !has_ccc && !has_tf && !has_grn) return(NULL)

      # -- Summary statistics --
      stats <- list()
      if (has_pw) {
        n_ora <- 0
        if (!is.null(res$pathway$ora$up)) n_ora <- n_ora + nrow(res$pathway$ora$up)
        if (!is.null(res$pathway$ora$down)) n_ora <- n_ora + nrow(res$pathway$ora$down)
        n_gsea <- if (!is.null(res$pathway$gsea)) nrow(res$pathway$gsea) else 0
        stats[["pw"]] <- paste0(n_ora, " ORA + ", n_gsea, " GSEA pathways")
      }
      if (has_ccc) {
        n_sig <- sum(res$ccc$pval_adj < 0.05, na.rm = TRUE)
        n_axes <- length(unique(paste0(res$ccc$source, res$ccc$target)))
        stats[["ccc"]] <- paste0(n_sig, " significant L-R interactions across ", n_axes, " axes")
      }
      if (has_tf) {
        n_tfs <- length(unique(res$tf$results$TF))
        stats[["tf"]] <- paste0(n_tfs, " transcription factors scored")
      }
      if (has_grn) {
        stats[["grn"]] <- paste0(igraph::vcount(res$grn$graph), " genes, ",
                                 igraph::ecount(res$grn$graph), " edges, ",
                                 nrow(res$grn$modules), " modules")
      }

      summary_items <- lapply(stats, function(s) {
        tags$span(style = "margin-right: 16px; color: #444; font-size: 12px;",
          tags$i(class = "fa fa-check-circle", style = "color:#4caf50; margin-right: 4px;"), s)
      })
      summary_bar <- div(
        style = paste0("background: linear-gradient(135deg, #e8f5e9, #f1f8e9); ",
                       "border-radius: 6px; padding: 10px 16px; margin-bottom: 16px; ",
                       "border-left: 4px solid #4caf50;"),
        do.call(tagList, summary_items)
      )

      # -- Figure sections --
      fig_card <- function(...) {
        div(style = paste0("border: 1px solid #e0e0e0; border-radius: 6px; ",
                           "padding: 4px; background: white;"), ...)
      }
      section_hdr <- function(icon_name, label, color) {
        tags$h6(style = paste0("color:", color, "; margin: 14px 0 8px 0; ",
                               "font-weight: 600; font-size: 13px;"),
          tags$i(class = paste0("fa fa-", icon_name), style = "margin-right: 6px;"), label)
      }

      sections <- list()
      if (has_pw) {
        sections[[length(sections) + 1]] <- tagList(
          section_hdr("dna", "Pathway Enrichment", "#1976d2"),
          fluidRow(
            column(6, fig_card(plotlyOutput(ns("ai_auto_ora"), height = "380px"))),
            column(6, fig_card(plotlyOutput(ns("ai_auto_gsea"), height = "380px")))
          )
        )
      }
      if (has_ccc) {
        sections[[length(sections) + 1]] <- tagList(
          section_hdr("project-diagram", "Cell-Cell Communication", "#388e3c"),
          fluidRow(
            column(6, fig_card(plotlyOutput(ns("ai_auto_ccc_chord"), height = "420px"))),
            column(6, fig_card(plotlyOutput(ns("ai_auto_ccc_heatmap"), height = "420px")))
          )
        )
      }
      if (has_tf || has_grn) {
        tf_grn_cols <- list()
        if (has_tf) tf_grn_cols[[length(tf_grn_cols) + 1]] <-
          column(6, fig_card(plotlyOutput(ns("ai_auto_tf_heatmap"), height = "420px")))
        if (has_grn) tf_grn_cols[[length(tf_grn_cols) + 1]] <-
          column(6, fig_card(plotlyOutput(ns("ai_auto_grn_network"), height = "420px")))
        sections[[length(sections) + 1]] <- tagList(
          section_hdr("network-wired", "Regulatory Analysis", "#e64a19"),
          do.call(fluidRow, tf_grn_cols)
        )
      }

      div(style = paste0("background: #fafafa; border-radius: 10px; padding: 18px; ",
                         "margin-top: 10px; border: 1px solid #e0e0e0;"),
        div(style = "display: flex; align-items: center; margin-bottom: 10px;",
          tags$i(class = "fa fa-chart-area",
                 style = "font-size: 18px; color: #1976d2; margin-right: 8px;"),
          tags$h5(style = "margin: 0; font-weight: 600; color: #333;", "Analysis Results")),
        summary_bar,
        tagList(sections),
        # AI Interpret All button
        div(style = paste0("margin-top: 16px; padding-top: 14px; ",
                           "border-top: 1px solid #e0e0e0; text-align: center;"),
          actionButton(ns("ai_interpret_all_auto"), "AI Interpret All Results",
            icon = icon("brain"), class = "ai-btn",
            style = "padding: 10px 28px; font-size: 14px; border-radius: 10px;"),
          tags$p(style = "margin: 6px 0 0 0; font-size: 11px; color: #999;",
            "Generate a comprehensive biological interpretation of all analyses")
        ),
        uiOutput(ns("ai_auto_interpretation"))
      )
    })

    # ---- 6 Auto-run figure renderers ----
    output$ai_auto_ora <- renderPlotly({
      res <- ai_auto_results()
      req(res, res$pathway)
      p <- tryCatch({
        enrich <- res$pathway
        if (!is.null(enrich$ora$up) && is.data.frame(enrich$ora$up) && nrow(enrich$ora$up) > 0) {
          create_ora_dot_plot(enrich$ora$up, "ORA - Upregulated Pathways")
        } else if (!is.null(enrich$ora$down) && is.data.frame(enrich$ora$down) && nrow(enrich$ora$down) > 0) {
          create_ora_dot_plot(enrich$ora$down, "ORA - Downregulated Pathways")
        } else {
          plotly_empty() %>% layout(title = "No ORA results")
        }
      }, error = function(e) plotly_empty() %>% layout(title = paste("ORA error:", e$message)))
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Auto-Run", figure_type = "ORA Dot Plot",
                       label = "Auto-Run ORA", source_module = "AI Analysis")
      }
      p
    })

    output$ai_auto_gsea <- renderPlotly({
      res <- ai_auto_results()
      req(res, res$pathway)
      p <- tryCatch({
        enrich <- res$pathway
        if (!is.null(enrich$gsea) && is.data.frame(enrich$gsea) && nrow(enrich$gsea) > 0) {
          create_gsea_plot(enrich$gsea)
        } else {
          plotly_empty() %>% layout(title = "No GSEA results")
        }
      }, error = function(e) plotly_empty() %>% layout(title = paste("GSEA error:", e$message)))
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Auto-Run", figure_type = "GSEA Plot",
                       label = "Auto-Run GSEA", source_module = "AI Analysis")
      }
      p
    })

    output$ai_auto_ccc_chord <- renderPlotly({
      res <- ai_auto_results()
      req(res, res$ccc)
      p <- .plot_ccc_chord_auto(res$ccc)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Auto-Run", figure_type = "CCC Chord",
                       label = "Auto-Run CCC Network", source_module = "AI Analysis")
      }
      p
    })

    output$ai_auto_ccc_heatmap <- renderPlotly({
      res <- ai_auto_results()
      req(res, res$ccc)
      p <- .plot_ccc_heatmap_auto(res$ccc)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Auto-Run", figure_type = "CCC Heatmap",
                       label = "Auto-Run L-R Heatmap", source_module = "AI Analysis")
      }
      p
    })

    output$ai_auto_tf_heatmap <- renderPlotly({
      res <- ai_auto_results()
      req(res, res$tf)
      p <- .plot_tf_heatmap_auto(res$tf$results)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Auto-Run", figure_type = "TF Heatmap",
                       label = "Auto-Run TF Activity", source_module = "AI Analysis")
      }
      p
    })

    output$ai_auto_grn_network <- renderPlotly({
      res <- ai_auto_results()
      req(res, res$grn)
      p <- .plot_grn_network_auto(res$grn)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Auto-Run", figure_type = "GRN Network",
                       label = "Auto-Run Gene Network", source_module = "AI Analysis")
      }
      p
    })

    # ---- AI Interpret All Auto-Run Results ----
    observeEvent(input$ai_interpret_all_auto, {
      res <- ai_auto_results()
      req(res)

      api_key <- trimws(values$ai_api_key %||% "")
      if (!nzchar(api_key)) {
        showNotification("Set your API key in AI Config first.", type = "error")
        return()
      }

      output$ai_auto_interpretation <- renderUI({
        ai_thread_loading("AI is interpreting all results...")
      })

      # Build comprehensive context from all auto-run results
      ctx <- c()

      # DE context
      if (!is.null(values$last_de_results) && is.data.frame(values$last_de_results)) {
        de <- values$last_de_results
        sig <- de[!is.na(de$p_val_adj) & de$p_val_adj < 0.05, ]
        up_g <- head(sig$gene[sig$avg_log2FC > 0][order(-sig$avg_log2FC[sig$avg_log2FC > 0])], 25)
        dn_g <- head(sig$gene[sig$avg_log2FC < 0][order(sig$avg_log2FC[sig$avg_log2FC < 0])], 25)
        ctx <- c(ctx, paste0(
          "## Differential Expression\n",
          "Total significant: ", nrow(sig), " genes\n",
          "Top upregulated: ", paste(up_g, collapse = ", "), "\n",
          "Top downregulated: ", paste(dn_g, collapse = ", ")))
      }

      # Pathway context
      if (!is.null(res$pathway)) {
        pw_text <- ""
        if (!is.null(res$pathway$ora$up) && nrow(res$pathway$ora$up) > 0) {
          top_pw <- head(as.character(res$pathway$ora$up$Description), 15)
          pw_text <- paste0(pw_text, "ORA Up: ", paste(top_pw, collapse = "; "), "\n")
        }
        if (!is.null(res$pathway$ora$down) && nrow(res$pathway$ora$down) > 0) {
          top_pw <- head(as.character(res$pathway$ora$down$Description), 15)
          pw_text <- paste0(pw_text, "ORA Down: ", paste(top_pw, collapse = "; "), "\n")
        }
        if (!is.null(res$pathway$gsea) && nrow(res$pathway$gsea) > 0) {
          gsea <- res$pathway$gsea
          top_gsea <- head(gsea[order(-abs(gsea$NES)), ], 15)
          pw_text <- paste0(pw_text, "GSEA top: ",
            paste(paste0(as.character(top_gsea$pathway), " (NES=", round(top_gsea$NES, 2), ")"),
                  collapse = "; "), "\n")
        }
        ctx <- c(ctx, paste0("## Pathway Enrichment\n", pw_text))
      }

      # CCC context
      if (!is.null(res$ccc) && is.data.frame(res$ccc) && nrow(res$ccc) > 0) {
        sig_ccc <- res$ccc[res$ccc$pval_adj < 0.05, ]
        top_int <- head(sig_ccc[order(sig_ccc$pval_adj), ], 20)
        ccc_text <- paste(paste0(top_int$source, "->", top_int$target, ": ",
                                 top_int$ligand, "-", top_int$receptor,
                                 " (", top_int$pathway, ", score=", round(top_int$score, 3), ")"),
                          collapse = "\n")
        pw_counts <- sort(table(sig_ccc$pathway), decreasing = TRUE)
        ctx <- c(ctx, paste0(
          "## Cell-Cell Communication\n",
          "Significant interactions: ", nrow(sig_ccc), "\n",
          "Top pathways: ", paste(paste0(names(head(pw_counts, 10)), " (", as.integer(head(pw_counts, 10)), ")"), collapse = ", "), "\n",
          "Top interactions:\n", ccc_text))
      }

      # TF context
      if (!is.null(res$tf) && !is.null(res$tf$results)) {
        tf_df <- res$tf$results
        top_tfs <- head(tf_df[order(-abs(tf_df$Activity)), ], 20)
        ctx <- c(ctx, paste0(
          "## Transcription Factor Activity\n",
          "Total TFs analysed: ", length(unique(tf_df$TF)), "\n",
          "Top variable TFs: ",
          paste(paste0(top_tfs$TF, " (", top_tfs$Cluster, ": ", round(top_tfs$Activity, 3), ")"),
                collapse = "; ")))
      }

      # GRN context
      if (!is.null(res$grn)) {
        n_v <- igraph::vcount(res$grn$graph)
        n_e <- igraph::ecount(res$grn$graph)
        mods <- res$grn$modules
        mod_text <- paste(paste0("Module ", mods$Module, ": hub=", mods$Hub_Gene,
                                 " (", mods$N_Genes, " genes: ", mods$Top_Genes, ")"),
                          collapse = "\n")
        ctx <- c(ctx, paste0(
          "## Gene Regulatory Network\n",
          n_v, " genes, ", n_e, " edges, ", nrow(mods), " modules\n",
          mod_text))
      }

      context_text <- paste(ctx, collapse = "\n\n")

      # Build prompt
      expert_persona <- if (exists("AI_EXPERT_PERSONA")) AI_EXPERT_PERSONA else ""
      prompt <- paste0(
        expert_persona, "\n\n",
        "Below are comprehensive results from an automated scRNA-seq downstream analysis:\n\n",
        "---\n", context_text, "\n---\n\n",
        "## Task: Integrated Biological Interpretation\n\n",
        "Provide a comprehensive, publication-quality biological interpretation that integrates ALL the analyses above.\n\n",
        "Structure your response as:\n\n",
        "## Summary\n",
        "2-3 sentence executive summary of the most important biological findings.\n\n",
        "## Key Biological Signals\n",
        "- Dominant biological processes, activated vs repressed gene programs\n",
        "- Connect DE genes to pathway results and known biology\n\n",
        "## Cell Communication & Signaling\n",
        "- Key ligand-receptor axes and their biological significance\n",
        "- Which cell populations are communicating and through what mechanisms\n\n",
        "## Regulatory Landscape\n",
        "- Master transcription factors driving observed expression programs\n",
        "- Network modules and hub genes: what do they regulate?\n\n",
        "## Integrated Model\n",
        "- How do all findings connect into a coherent biological narrative?\n",
        "- Proposed mechanistic model linking DE, pathways, communication, and TF activity\n\n",
        "## Recommended Validation Experiments\n",
        "- 3-5 specific, actionable experiments to validate key findings\n",
        "- Include technique, expected result, and priority\n\n",
        "Use specific gene names and pathway names from the data. Be precise and evidence-based."
      )

      ai_model <- values$ai_model %||% get_default_model(values$ai_api_key)

      withProgress(message = "AI interpreting all results...", value = 0.3, {
        result <- tryCatch({
          if (exists("call_ai", mode = "function")) {
            call_ai("comprehensive_interpretation", prompt, values)
          } else {
            call_claude_api(api_key, ai_model,
              list(list(role = "user", content = prompt)),
              max_tokens = 16384L, timeout = 300)
          }
        }, error = function(e) {
          list(ok = FALSE, error = e$message)
        })
      })

      if (isTRUE(result$ok)) {
        html_content <- .mdtohtml(result$text)
        output$ai_auto_interpretation <- renderUI({
          if (exists("ai_thread_message", mode = "function")) {
            ai_thread_message(
              content_html = html_content,
              title = "Integrated Analysis Interpretation",
              model = ai_model,
              tokens_in = result$tokens_in,
              tokens_out = result$tokens_out
            )
          } else {
            div(class = "ai-thread-content", style = "padding: 16px;", HTML(html_content))
          }
        })
      } else {
        err_msg <- result$error %||% "Unknown error during AI interpretation"
        output$ai_auto_interpretation <- renderUI({
          if (exists("ai_thread_error", mode = "function")) {
            ai_thread_error(err_msg)
          } else {
            div(style = "color: #c62828; padding: 12px;", icon("exclamation-triangle"), err_msg)
          }
        })
      }
    })

    # ---- Populate comparison selectors ----
    observe({
      req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) {
        (is.factor(x) || is.character(x)) && length(unique(x)) >= 2 && length(unique(x)) <= 100
      })]
      updateSelectInput(session, "ai_compare_col", choices = sort(cat_cols))
    })

    observeEvent(input$ai_compare_col, {
      req(values$active_dataset, input$ai_compare_col)
      grps <- sort(unique(as.character(values$active_dataset@meta.data[[input$ai_compare_col]])))
      updateSelectInput(session, "ai_compare_g1", choices = grps, selected = grps[1])
      updateSelectInput(session, "ai_compare_g2", choices = grps,
        selected = if (length(grps) >= 2) grps[2] else grps[1])
    })

    # ---- Populate annotation column selector ----
    observe({
      req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) {
        (is.factor(x) || is.character(x)) && length(unique(x)) >= 2 && length(unique(x)) <= 200
      })]
      # Prioritize annotation-like columns
      annot_like <- grep("cell_type|annotation|label|singler|consensus", cat_cols,
                         value = TRUE, ignore.case = TRUE)
      others <- setdiff(cat_cols, annot_like)
      ordered <- c(sort(annot_like), sort(others))
      sel <- if (length(annot_like) > 0) annot_like[1] else if (length(ordered) > 0) ordered[1] else NULL
      updateSelectInput(session, "ai_annot_col", choices = ordered, selected = sel)
    })

    # ---- Context status badges ----
    output$ai_context_status <- renderUI({
      badges <- list()

      has_data <- !is.null(values$active_dataset)
      badges[[1]] <- .make_badge("Dataset", has_data,
        if (has_data) paste0(ncol(values$active_dataset), " cells") else "Not loaded")

      has_de <- !is.null(values$last_de_results) &&
                is.data.frame(values$last_de_results) && nrow(values$last_de_results) > 0
      badges[[2]] <- .make_badge("DE", has_de,
        if (has_de) paste0(nrow(values$last_de_results), " genes") else "Not run")

      has_pw <- !is.null(values$last_enrich_results)
      badges[[3]] <- .make_badge("Pathways", has_pw,
        if (has_pw) "Available" else "Not run")

      has_cc <- !is.null(tryCatch(values$last_cc_results, error = function(e) NULL)) &&
                is.data.frame(tryCatch(values$last_cc_results, error = function(e) NULL))
      badges[[4]] <- .make_badge("CCC", has_cc,
        if (has_cc) paste0(nrow(values$last_cc_results), " interactions") else "Not run")

      has_tf <- !is.null(tryCatch(values$last_tf_results, error = function(e) NULL)) &&
                is.data.frame(tryCatch(values$last_tf_results, error = function(e) NULL))
      badges[[5]] <- .make_badge("TF Activity", has_tf,
        if (has_tf) "Available" else "Not run")

      has_annot <- has_data && any(grepl("cell_type|celltype|annotation|label",
        colnames(values$active_dataset@meta.data), ignore.case = TRUE))
      badges[[6]] <- .make_badge("Annotations", has_annot,
        if (has_annot) "Available" else "Not found")

      tagList(badges)
    })

    # ---- Build comprehensive data context ----
    build_full_context <- function() {
      ctx_parts <- list()
      selected <- input$ai_context

      # 1. Dataset Summary
      if ("summary" %in% selected && !is.null(values$active_dataset)) {
        tryCatch({
          obj <- values$active_dataset
          meta <- obj@meta.data
          n_cells <- ncol(obj)
          n_genes <- nrow(obj)
          assay <- DefaultAssay(obj)

          # Available reductions
          reductions <- tryCatch(names(obj@reductions), error = function(e) character(0))

          # Numeric metadata summary
          num_cols <- colnames(meta)[sapply(meta, is.numeric)]
          num_summary <- ""
          for (col in head(num_cols, 6)) {
            vals <- meta[[col]]
            vals <- vals[is.finite(vals)]
            if (length(vals) > 0) {
              num_summary <- paste0(num_summary,
                "\n  - ", col, ": median=", round(median(vals), 1),
                ", mean=", round(mean(vals), 1),
                ", range=[", round(min(vals), 1), ", ", round(max(vals), 1), "]")
            }
          }

          # Categorical metadata
          cat_cols <- colnames(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
          cat_summary <- ""
          for (col in head(cat_cols, 5)) {
            tbl <- sort(table(meta[[col]]), decreasing = TRUE)
            cat_summary <- paste0(cat_summary,
              "\n  - ", col, " (", length(tbl), " levels): ",
              paste(head(paste0(names(tbl), "=", as.integer(tbl)), 15), collapse = ", "),
              if (length(tbl) > 15) "..." else "")
          }

          ctx_parts$summary <- paste0(
            "## Dataset Summary\n",
            "- Total cells: ", format(n_cells, big.mark = ","), "\n",
            "- Total genes: ", format(n_genes, big.mark = ","), "\n",
            "- Default assay: ", assay, "\n",
            "- Available reductions: ", paste(reductions, collapse = ", "), "\n",
            "- Numeric metadata:", num_summary, "\n",
            "- Categorical metadata:", cat_summary, "\n"
          )
        }, error = function(e) NULL)
      }

      # 2. DE Results
      if ("de" %in% selected && !is.null(values$last_de_results)) {
        tryCatch({
          de <- values$last_de_results
          if (is.data.frame(de) && nrow(de) > 0) {
            # Determine column names flexibly
            gene_col <- intersect(c("gene", "Gene", "feature", "row.names"), colnames(de))[1]
            fc_col <- intersect(c("avg_log2FC", "logFC", "log2FoldChange"), colnames(de))[1]
            padj_col <- intersect(c("p_val_adj", "padj", "FDR", "adj.P.Val"), colnames(de))[1]

            if (!is.na(gene_col) && !is.na(fc_col)) {
              sig <- de[!is.na(de[[padj_col]]) & de[[padj_col]] < 0.05, ]
              up_genes <- sig[sig[[fc_col]] > 0, ]
              up_genes <- up_genes[order(-up_genes[[fc_col]]), ]
              down_genes <- sig[sig[[fc_col]] < 0, ]
              down_genes <- down_genes[order(down_genes[[fc_col]]), ]

              format_gene <- function(row) {
                paste0("  ", row[[gene_col]],
                       " (log2FC=", round(row[[fc_col]], 2),
                       ", padj=", signif(row[[padj_col]], 3), ")")
              }

              up_text <- if (nrow(up_genes) > 0) {
                paste(apply(head(up_genes, 30), 1, format_gene), collapse = "\n")
              } else "  None"

              down_text <- if (nrow(down_genes) > 0) {
                paste(apply(head(down_genes, 30), 1, format_gene), collapse = "\n")
              } else "  None"

              ctx_parts$de <- paste0(
                "## Differential Expression Results\n",
                "- Total genes tested: ", format(nrow(de), big.mark = ","), "\n",
                "- Significant (padj < 0.05): ", nrow(sig), "\n",
                "- Upregulated: ", nrow(up_genes), "\n",
                "- Downregulated: ", nrow(down_genes), "\n\n",
                "### Top 30 Upregulated Genes\n", up_text, "\n\n",
                "### Top 30 Downregulated Genes\n", down_text, "\n"
              )
            }
          }
        }, error = function(e) NULL)
      }

      # 3. Pathway Enrichment
      if ("pathway" %in% selected && !is.null(values$last_enrich_results)) {
        tryCatch({
          enrich <- values$last_enrich_results
          pw_text <- ""

          if (!is.null(enrich$ora$up) && is.data.frame(enrich$ora$up) && nrow(enrich$ora$up) > 0) {
            top <- head(enrich$ora$up, 20)
            desc_col <- intersect(c("Description", "ID", "pathway"), colnames(top))[1]
            padj_col <- intersect(c("p.adjust", "padj", "qvalue"), colnames(top))[1]
            count_col <- intersect(c("Count", "count", "size"), colnames(top))[1]
            if (!is.na(desc_col) && !is.na(padj_col)) {
              pw_text <- paste0(pw_text, "### ORA - Upregulated Pathways\n")
              for (i in seq_len(nrow(top))) {
                cnt <- if (!is.na(count_col)) paste0(", genes=", top[[count_col]][i]) else ""
                ratio <- if ("GeneRatio" %in% colnames(top)) paste0(", ratio=", top$GeneRatio[i]) else ""
                pw_text <- paste0(pw_text, "  ", i, ". ", top[[desc_col]][i],
                  " (padj=", signif(as.numeric(as.character(top[[padj_col]][i])), 3), cnt, ratio, ")\n")
              }
              pw_text <- paste0(pw_text, "\n")
            }
          }

          if (!is.null(enrich$ora$down) && is.data.frame(enrich$ora$down) && nrow(enrich$ora$down) > 0) {
            top <- head(enrich$ora$down, 20)
            desc_col <- intersect(c("Description", "ID", "pathway"), colnames(top))[1]
            padj_col <- intersect(c("p.adjust", "padj", "qvalue"), colnames(top))[1]
            if (!is.na(desc_col) && !is.na(padj_col)) {
              pw_text <- paste0(pw_text, "### ORA - Downregulated Pathways\n")
              for (i in seq_len(nrow(top))) {
                pw_text <- paste0(pw_text, "  ", i, ". ", top[[desc_col]][i],
                  " (padj=", signif(as.numeric(as.character(top[[padj_col]][i])), 3), ")\n")
              }
              pw_text <- paste0(pw_text, "\n")
            }
          }

          if (!is.null(enrich$gsea) && is.data.frame(enrich$gsea) && nrow(enrich$gsea) > 0) {
            top <- head(enrich$gsea, 20)
            nm_col <- intersect(c("pathway", "Description", "ID"), colnames(top))[1]
            if (!is.na(nm_col) && "NES" %in% colnames(top)) {
              pw_text <- paste0(pw_text, "### GSEA Results (ranked by |NES|)\n")
              top <- top[order(-abs(top$NES)), ]
              for (i in seq_len(nrow(top))) {
                direction <- if (top$NES[i] > 0) "UP" else "DOWN"
                pw_text <- paste0(pw_text, "  ", i, ". ", top[[nm_col]][i],
                  " (NES=", round(top$NES[i], 2),
                  ", padj=", signif(top$padj[i], 3),
                  ", direction=", direction, ")\n")
              }
              pw_text <- paste0(pw_text, "\n")
            }
          }

          if (nzchar(pw_text)) {
            ctx_parts$pathway <- paste0("## Pathway & Gene Set Enrichment\n", pw_text)
          }
        }, error = function(e) NULL)
      }

      # 4. Cell-Cell Communication
      if ("ccc" %in% selected) {
        tryCatch({
          cc_res <- values$last_cc_results
          if (!is.null(cc_res) && is.data.frame(cc_res) && nrow(cc_res) > 0) {
            sig <- cc_res[!is.na(cc_res$pval_adj) & cc_res$pval_adj < 0.05, ]
            sig <- sig[order(sig$pval_adj), ]

            # Pathway summary
            pw_counts <- sort(table(sig$pathway), decreasing = TRUE)
            # Source-target summary
            sig$pair <- paste0(sig$source, " -> ", sig$target)
            pair_counts <- sort(table(sig$pair), decreasing = TRUE)

            top_int <- head(sig, 25)
            int_text <- paste(paste0("  ", top_int$source, " -> ", top_int$target,
              ": ", top_int$ligand, "-", top_int$receptor,
              " (", top_int$pathway, ", score=", round(top_int$score, 3),
              ", padj=", signif(top_int$pval_adj, 3), ")"), collapse = "\n")

            ctx_parts$ccc <- paste0(
              "## Cell-Cell Communication\n",
              "- Total interactions tested: ", nrow(cc_res), "\n",
              "- Significant (padj < 0.05): ", nrow(sig), "\n",
              "- Pathway distribution: ", paste(head(paste0(names(pw_counts), "=", as.integer(pw_counts)), 8), collapse = ", "), "\n",
              "- Most active pairs: ", paste(head(paste0(names(pair_counts), " (", as.integer(pair_counts), ")"), 6), collapse = ", "), "\n\n",
              "### Top 25 Significant Interactions\n", int_text, "\n"
            )
          }
        }, error = function(e) NULL)
      }

      # 5. TF Activity
      if ("tf" %in% selected) {
        tryCatch({
          tf_res <- values$last_tf_results
          if (!is.null(tf_res) && is.data.frame(tf_res) && nrow(tf_res) > 0) {
            # Top activating TFs per cluster
            clusters <- unique(tf_res$Cluster)
            tf_text <- ""
            for (cl in head(sort(clusters), 10)) {
              cl_data <- tf_res[tf_res$Cluster == cl, ]
              cl_data <- cl_data[order(-cl_data$Activity), ]
              top_act <- head(cl_data, 5)
              bot_act <- tail(cl_data[order(cl_data$Activity), ], 3)
              tf_text <- paste0(tf_text, "  - Cluster ", cl,
                ": activated=[", paste(paste0(top_act$TF, "(", round(top_act$Activity, 2), ")"), collapse = ", "),
                "], repressed=[", paste(paste0(bot_act$TF, "(", round(bot_act$Activity, 2), ")"), collapse = ", "),
                "]\n")
            }
            ctx_parts$tf <- paste0(
              "## Transcription Factor Activity\n",
              "- TFs analyzed: ", length(unique(tf_res$TF)), "\n",
              "- Clusters: ", length(clusters), "\n",
              "### Per-Cluster TF Activity (top activated & repressed)\n", tf_text
            )
          }
        }, error = function(e) NULL)
      }

      # 6. Cell type annotations (use selected annotation column if available)
      if ("annotations" %in% selected && !is.null(values$active_dataset)) {
        tryCatch({
          meta <- values$active_dataset@meta.data
          # Prioritize user-selected column
          sel_col <- input$ai_annot_col
          annot_cols <- grep("cell_type|celltype|annotation|label|singler|ai_predicted",
            colnames(meta), ignore.case = TRUE, value = TRUE)
          if (!is.null(sel_col) && nzchar(sel_col) && sel_col %in% colnames(meta)) {
            annot_cols <- unique(c(sel_col, annot_cols))
          }
          if (length(annot_cols) > 0) {
            annot_text <- ""
            for (col in head(annot_cols, 3)) {
              tbl <- sort(table(meta[[col]]), decreasing = TRUE)
              annot_text <- paste0(annot_text, "### ", col, "\n",
                paste(paste0("  - ", names(tbl), ": ", as.integer(tbl), " cells (",
                  round(as.integer(tbl) / sum(tbl) * 100, 1), "%)"), collapse = "\n"), "\n\n")
            }
            ctx_parts$annotations <- paste0("## Cell Type Annotations\n", annot_text)
          }
        }, error = function(e) NULL)
      }

      # 7. Cluster markers (use cached if available, otherwise compute fast)
      if ("markers" %in% selected && !is.null(values$active_dataset)) {
        tryCatch({
          obj <- values$active_dataset
          meta <- obj@meta.data

          # Check for cached markers first
          markers <- tryCatch(values$last_cluster_markers, error = function(e) NULL)

          if (is.null(markers) || !is.data.frame(markers) || nrow(markers) == 0) {
            # Compute markers with aggressive speed settings
            cluster_col <- intersect(c("seurat_clusters", "cluster", "Cluster"),
              colnames(meta))[1]
            if (!is.na(cluster_col) && length(unique(meta[[cluster_col]])) >= 2) {
              Idents(obj) <- cluster_col
              markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.3,
                logfc.threshold = 0.7, max.cells.per.ident = 300,
                test.use = "wilcox", verbose = FALSE)
              # Cache for future use
              values$last_cluster_markers <- markers
            }
          }

          if (!is.null(markers) && is.data.frame(markers) && nrow(markers) > 0) {
            top_markers <- markers %>%
              group_by(cluster) %>%
              slice_max(order_by = avg_log2FC, n = 5) %>%
              ungroup()

            mk_text <- ""
            for (cl in sort(unique(as.character(top_markers$cluster)))) {
              cl_mk <- top_markers[as.character(top_markers$cluster) == cl, ]
              mk_text <- paste0(mk_text, "  - Cluster ", cl, ": ",
                paste(paste0(cl_mk$gene, "(FC=", round(cl_mk$avg_log2FC, 1), ")"),
                  collapse = ", "), "\n")
            }
            ctx_parts$markers <- paste0(
              "## Cluster Marker Genes (top 5 per cluster)\n", mk_text)
          }
        }, error = function(e) NULL)
      }

      return(ctx_parts)
    }

    # ---- Build prompt based on mode ----
    build_prompt <- function() {
      ctx <- build_full_context()
      mode <- input$ai_mode

      context_text <- paste(unlist(ctx), collapse = "\n\n")
      if (!nzchar(trimws(context_text))) {
        context_text <- "[No analysis data available. Please load a dataset and run analyses first.]"
      }

      # Use the shared AI_EXPERT_PERSONA if available, otherwise inline version
      expert_persona <- if (exists("AI_EXPERT_PERSONA")) AI_EXPERT_PERSONA else paste0(
        "You are a senior principal investigator and expert computational biologist with 15+ ",
        "years of experience in single-cell genomics. You are recognised as a leading authority ",
        "in scRNA-seq analysis, multi-omics integration, and translational bioinformatics.")
      system_preamble <- paste0(
        expert_persona, "\n\n",
        "Below are results from a scRNA-seq analysis pipeline:\n\n",
        "---\n", context_text, "\n---\n\n"
      )

      prompt <- switch(mode,
        "interpret" = paste0(system_preamble,
          "## Task: Comprehensive Biological Interpretation\n\n",
          "Provide a detailed biological interpretation of these results. Structure your response as:\n\n",
          "### 1. Executive Summary\n",
          "2-3 sentences summarizing the most important biological findings.\n\n",
          "### 2. Key Biological Signals\n",
          "- What are the dominant biological processes revealed by DE and pathway analysis?\n",
          "- Which gene programs are activated vs repressed?\n",
          "- Are there known disease signatures or cell state transitions?\n\n",
          "### 3. Cell Population Characterization\n",
          "- What do expression patterns suggest about cell identities?\n",
          "- Are there activated, stressed, proliferating, or transitional states?\n",
          "- How do clusters relate to known cell types?\n\n",
          "### 4. Signaling & Communication\n",
          "- What signaling axes are active?\n",
          "- How do cell populations communicate?\n",
          "- Which pathways coordinate the observed changes?\n\n",
          "### 5. Regulatory Landscape\n",
          "- Which transcription factors drive the observed expression programs?\n",
          "- Are there master regulators or regulatory hubs?\n\n",
          "### 6. Clinical/Biological Relevance\n",
          "- How do these findings relate to known biology and disease?\n",
          "- What is the potential translational significance?\n\n",
          "### 7. Limitations & Caveats\n",
          "- What are the analytical limitations?\n",
          "- What alternative interpretations exist?\n\n",
          "Use specific gene names and pathway names from the data. Be precise and evidence-based."
        ),

        "hypotheses" = paste0(system_preamble,
          "## Task: Generate Testable Biological Hypotheses\n\n",
          "Based on these scRNA-seq analysis results, generate rigorous, testable hypotheses.\n\n",
          "### Primary Hypotheses (5 hypotheses, strongest evidence)\n",
          "For EACH hypothesis provide:\n",
          "1. **Hypothesis Statement**: Clear, falsifiable statement\n",
          "2. **Supporting Evidence**: Specific genes, pathways, scores from the data above\n",
          "3. **Mechanism**: Proposed biological mechanism\n",
          "4. **Experimental Test**: Specific experiment (wet lab) to validate\n",
          "   - Technique (qPCR, Western blot, CRISPR KO, flow cytometry, etc.)\n",
          "   - Expected result if hypothesis is true\n",
          "   - Expected result if hypothesis is false\n",
          "5. **Computational Validation**: Additional bioinformatics analysis to strengthen evidence\n",
          "6. **Priority Score**: High/Medium/Low based on impact and feasibility\n\n",
          "### Secondary Hypotheses (3-5 hypotheses, require more evidence)\n",
          "Same format but briefer.\n\n",
          "### Integrated Model\n",
          "Propose a coherent biological model that connects 2-3 of the primary hypotheses.\n",
          "Include a text description of how the key findings integrate.\n\n",
          "Be specific: use actual gene names and pathway names from the data."
        ),

        "report" = paste0(system_preamble,
          "## Task: Publication-Quality Analysis Report\n\n",
          "Generate a comprehensive, publication-quality analysis report suitable for ",
          "inclusion in a research paper supplementary material or lab notebook. Use formal scientific language.\n\n",
          "### 1. Abstract (150-200 words)\n\n",
          "### 2. Dataset Overview\n",
          "- Sample characteristics, cell numbers, quality metrics\n",
          "- Clustering resolution and cell population composition\n\n",
          "### 3. Differential Expression Analysis\n",
          "- Summary of DE analysis results\n",
          "- Key upregulated genes: biological function, known roles, significance\n",
          "- Key downregulated genes: biological function, known roles, significance\n",
          "- Functional categorization (immune response, metabolism, signaling, etc.)\n\n",
          "### 4. Pathway & Gene Set Enrichment\n",
          "- Significantly enriched pathways and their biological meaning\n",
          "- Pathway crosstalk and interconnections\n",
          "- Comparison between ORA and GSEA results\n\n",
          "### 5. Cell-Cell Communication Analysis\n",
          "- Overview of significant ligand-receptor interactions\n",
          "- Key communication axes between cell populations\n",
          "- Biological significance of top interactions\n\n",
          "### 6. Transcription Factor Activity\n",
          "- Key TFs driving expression programs\n",
          "- Regulatory circuits and master regulators\n\n",
          "### 7. Integrated Biological Model\n",
          "- How all analyses converge into a coherent biological narrative\n",
          "- Proposed mechanisms and regulatory circuits\n",
          "- Comparison with published findings in the field\n\n",
          "### 8. Key Findings & Conclusions\n",
          "- Numbered list of 5-8 key findings\n",
          "- Each finding supported by specific data\n\n",
          "### 9. Future Directions\n",
          "- Recommended follow-up experiments\n",
          "- Additional computational analyses\n",
          "- Potential clinical/translational implications\n\n",
          "### 10. Methods Summary\n",
          "- Brief description of analytical approaches and tools\n\n",
          "Use markdown tables where helpful. Be thorough and precise."
        ),

        "compare" = {
          g1 <- input$ai_compare_g1
          g2 <- input$ai_compare_g2
          col <- input$ai_compare_col

          # Extract comparison-specific data
          compare_ctx <- ""
          tryCatch({
            if (!is.null(values$active_dataset) && !is.null(col)) {
              obj <- values$active_dataset
              meta <- obj@meta.data
              g1_cells <- sum(meta[[col]] == g1, na.rm = TRUE)
              g2_cells <- sum(meta[[col]] == g2, na.rm = TRUE)
              compare_ctx <- paste0(
                "\n## Comparison: ", g1, " vs ", g2, " (column: ", col, ")\n",
                "- ", g1, ": ", g1_cells, " cells\n",
                "- ", g2, ": ", g2_cells, " cells\n")
            }
          }, error = function(e) NULL)

          paste0(system_preamble, compare_ctx,
            "\n## Task: Cluster/Condition Comparison\n\n",
            "Compare **", g1, "** vs **", g2, "** populations:\n\n",
            "### 1. Key Differences\n",
            "- Which genes are specifically up/down in each group?\n",
            "- What biological processes distinguish them?\n\n",
            "### 2. Shared Features\n",
            "- What programs are active in both populations?\n\n",
            "### 3. Cell State Interpretation\n",
            "- What biological states do these groups represent?\n",
            "- Are they different stages of the same process or distinct cell types?\n\n",
            "### 4. Regulatory Differences\n",
            "- Which TFs/pathways are differentially active?\n\n",
            "### 5. Communication Differences\n",
            "- How do signaling patterns differ?\n\n",
            "### 6. Biological Significance\n",
            "- What is the functional implication of these differences?\n\n",
            "Base your analysis on the DE genes, pathways, and other data provided above."
          )
        },

        "methods" = paste0(system_preamble,
          "## Task: Methods Section for Publication\n\n",
          "Write a detailed Methods section suitable for a peer-reviewed publication. ",
          "Include:\n\n",
          "### Single-cell RNA-seq Data Processing\n",
          "- Quality control criteria and filtering\n",
          "- Normalization and scaling approach\n",
          "- Feature selection method\n\n",
          "### Dimensionality Reduction and Clustering\n",
          "- PCA, UMAP/t-SNE parameters\n",
          "- Clustering algorithm and resolution\n\n",
          "### Differential Expression Analysis\n",
          "- Statistical test used\n",
          "- Multiple testing correction\n",
          "- Significance thresholds\n\n",
          "### Pathway and Gene Set Enrichment\n",
          "- Databases used (MSigDB, KEGG, Reactome, etc.)\n",
          "- ORA and GSEA methods and parameters\n\n",
          "### Cell-Cell Communication\n",
          "- Ligand-receptor database source\n",
          "- Scoring method and permutation testing\n\n",
          "### Transcription Factor Activity\n",
          "- TF regulon database and scoring method\n\n",
          "### Software and Versions\n",
          "- R version, Seurat, and other key packages\n\n",
          "Write in third person, past tense, formal scientific style. ",
          "Infer reasonable parameters from the dataset characteristics."
        ),

        "legends" = paste0(system_preamble,
          "## Task: Generate Figure Legends\n\n",
          "Write publication-quality figure legends for common scRNA-seq figures. ",
          "For each figure below, write a detailed legend that includes:\n",
          "- Panel labels (a, b, c, ...)\n",
          "- What is shown (plot type, axes, colors)\n",
          "- Key observations to highlight\n",
          "- Statistical details where relevant\n\n",
          "Generate legends for:\n\n",
          "**Figure 1: Dataset Overview**\n",
          "(a) UMAP of all cells colored by cluster\n",
          "(b) UMAP colored by sample/condition\n",
          "(c) Cell count bar plot per cluster\n\n",
          "**Figure 2: Quality Control**\n",
          "(a) Violin plots of nFeature, nCount, percent.mt per cluster\n",
          "(b) Scatter plot nCount vs nFeature colored by percent.mt\n\n",
          "**Figure 3: Differential Expression**\n",
          "(a) Volcano plot of DE genes\n",
          "(b) Heatmap of top DE genes across clusters\n\n",
          "**Figure 4: Pathway Analysis**\n",
          "(a) Dot plot of top enriched pathways\n",
          "(b) GSEA enrichment plot for top pathways\n\n",
          "**Figure 5: Cell Communication**\n",
          "(a) Chord diagram of significant interactions\n",
          "(b) Heatmap of ligand-receptor interaction scores\n\n",
          "Use the actual gene names and pathway names from the data above."
        ),

        "literature" = paste0(system_preamble,
          "## Task: Literature Connections & Context\n\n",
          "Connect these scRNA-seq findings to published literature:\n\n",
          "### 1. Gene-Level Literature Context\n",
          "For the top 15 DE genes, provide:\n",
          "- Known biological function\n",
          "- Relevant publications and disease associations\n",
          "- Known regulatory mechanisms\n\n",
          "### 2. Pathway-Level Context\n",
          "For the top enriched pathways:\n",
          "- Role in relevant biological processes\n",
          "- Known cross-talk with other enriched pathways\n",
          "- Therapeutic targets within these pathways\n\n",
          "### 3. Cell Type/State Context\n",
          "- How do the observed cell populations compare to published single-cell atlases?\n",
          "- Are there known cell state transitions or differentiation trajectories?\n\n",
          "### 4. Disease Relevance\n",
          "- Which findings are associated with specific diseases?\n",
          "- Are there known drug targets among the DE genes?\n",
          "- Relevant clinical trials or therapeutic implications\n\n",
          "### 5. Methodological Context\n",
          "- How do these results compare to similar published analyses?\n",
          "- Any notable discrepancies or novel findings?\n\n",
          "Cite specific studies and databases where possible (e.g., PMID numbers, database names)."
        ),

        "marker_interp" = paste0(system_preamble,
          "## Task: Marker Gene Interpretation for Cell Type Identification\n\n",
          "Analyze the marker genes from each cluster to identify cell types:\n\n",
          "### For each cluster:\n",
          "1. **Top markers**: List the most informative marker genes\n",
          "2. **Cell type prediction**: Most likely cell type based on markers\n",
          "3. **Confidence**: How confident is the assignment? (High/Medium/Low)\n",
          "4. **Key evidence**: Which specific markers support the assignment?\n",
          "5. **Alternative interpretations**: Could these markers indicate a different cell type?\n",
          "6. **Cell state**: Is this a resting, activated, proliferating, or transitional state?\n\n",
          "### Summary Table\n",
          "Provide a markdown table: | Cluster | Cell Type | Confidence | Key Markers | State |\n\n",
          "### Novel or Rare Populations\n",
          "Identify any clusters that may represent:\n",
          "- Novel cell states\n",
          "- Transitional/intermediate populations\n",
          "- Doublets or low-quality clusters\n\n",
          "Use actual marker genes from the data above."
        ),

        "drug_targets" = paste0(system_preamble,
          "## Task: Drug Target Analysis from scRNA-seq Data\n\n",
          "Identify potential drug targets from the differential expression and pathway results:\n\n",
          "### 1. Druggable Gene Targets\n",
          "For each significantly DE gene that is a known drug target:\n",
          "- Gene name and expression change\n",
          "- Known drugs targeting this gene/protein\n",
          "- Drug class (small molecule, antibody, etc.)\n",
          "- FDA approval status\n",
          "- Relevance to the observed biology\n\n",
          "### 2. Druggable Pathway Targets\n",
          "For enriched pathways with known therapeutic interventions:\n",
          "- Pathway name\n",
          "- Available drugs/inhibitors\n",
          "- Clinical stage\n\n",
          "### 3. Surface Markers for Targeting\n",
          "- Identify upregulated surface proteins (CD markers, receptors)\n",
          "- Potential for antibody-drug conjugates or CAR-T therapy\n\n",
          "### 4. Combination Therapy Opportunities\n",
          "- Based on pathway crosstalk, suggest rational drug combinations\n\n",
          "### 5. Priority Target List\n",
          "Provide a ranked table: | Rank | Gene | Drug | Mechanism | Evidence | Feasibility |\n\n",
          "Be specific about drug names and cite known clinical applications."
        ),

        "custom" = paste0(system_preamble,
          "## User Question\n",
          if (is.null(input$ai_custom_query)) "" else input$ai_custom_query,
          "\n\n",
          "Please provide a thorough, scientifically rigorous answer based on the scRNA-seq data above. ",
          "Structure your response with clear headers. ",
          "Reference specific genes, pathways, and cell types from the data. ",
          "Suggest follow-up analyses and experiments where relevant.\n"
        )
      )

      return(prompt)
    }

    # ---- Preview prompt ----
    observeEvent(input$ai_preview, {
      show_preview(!show_preview())
    })

    output$ai_prompt_preview_area <- renderUI({
      if (!show_preview()) return(NULL)

      prompt <- tryCatch(build_prompt(), error = function(e)
        paste("Error building prompt:", e$message))

      fluidRow(column(12,
        div(class = "ai-section-card",
          style = "background: #fafafa;",
          div(style = "display:flex; justify-content:space-between; align-items:center;",
            h5(style = "margin:0;", "Prompt Preview (editable before sending)"),
            span(style = "color:#888; font-size:11px;",
              paste(nchar(prompt), "characters"))
          ),
          textAreaInput(ns("ai_prompt_edit"), NULL,
            value = prompt, rows = 15, width = "100%"),
          p(style = "color: #888; font-size: 11px; margin-top: 5px;",
            "Edit the prompt above then click 'Generate' to send it.")
        )
      ))
    })

    # ---- Clear ----
    observeEvent(input$ai_clear, {
      chat_history(list())
      api_messages(list())
      show_preview(FALSE)
    })

    # ---- Main generate ----
    observeEvent(input$ai_run, {
      api_key <- values$ai_api_key
      if (is.null(api_key) || !nzchar(trimws(api_key))) {
        showNotification("Please set your Claude API key in Data Management > AI Config first.", type = "error", duration = 8)
        return()
      }

      if (input$ai_mode == "custom" &&
          (is.null(input$ai_custom_query) || !nzchar(trimws(input$ai_custom_query)))) {
        showNotification("Please enter your question in the text box above.", type = "error")
        return()
      }

      if (is.null(values$active_dataset)) {
        showNotification("Please load a dataset before running AI analysis.", type = "warning")
        return()
      }

      is_processing(TRUE)

      # Build prompt (use edited version if preview is open)
      prompt <- if (show_preview() && !is.null(input$ai_prompt_edit) && nzchar(input$ai_prompt_edit)) {
        input$ai_prompt_edit
      } else {
        tryCatch(build_prompt(), error = function(e) {
          showNotification(paste("Error building prompt:", e$message), type = "error")
          is_processing(FALSE)
          return(NULL)
        })
      }
      if (is.null(prompt)) return()

      withProgress(message = "AI Analysis in progress...", value = 0, {
        incProgress(0.1, detail = "Preparing request...")

        # Build messages for API
        msgs <- list(list(role = "user", content = prompt))

        ai_model <- values$ai_model %||% get_default_model(values$ai_api_key)
        incProgress(0.2, detail = paste("Calling", ai_model, "..."))
        result <- call_claude_api(
          api_key = trimws(api_key),
          model = ai_model,
          messages = msgs,
          max_tokens = as.integer(input$ai_max_tokens),
          timeout = 240
        )

        if (!result$ok) {
          showNotification(result$error, type = "error", duration = 12)
          is_processing(FALSE)
          return()
        }

        incProgress(0.9, detail = "Rendering response...")

        # Store in chat history
        mode_label <- switch(input$ai_mode,
          "interpret" = "Interpret Results",
          "hypotheses" = "Generate Hypotheses",
          "report" = "Full Report",
          "compare" = "Compare Clusters",
          "methods" = "Methods Section",
          "legends" = "Figure Legends",
          "literature" = "Literature Connections",
          "marker_interp" = "Marker Interpretation",
          "drug_targets" = "Drug Target Analysis",
          "custom" = "Custom Query",
          input$ai_mode)

        user_msg <- list(
          role = "user",
          content = paste0("[", mode_label, "] ", if (input$ai_mode == "custom") input$ai_custom_query else mode_label),
          timestamp = Sys.time(),
          tokens = result$tokens_in
        )
        assistant_msg <- list(
          role = "assistant",
          content = result$text,
          timestamp = Sys.time(),
          tokens = result$tokens_out
        )

        # Update state
        hist <- chat_history()
        hist[[length(hist) + 1]] <- user_msg
        hist[[length(hist) + 1]] <- assistant_msg
        chat_history(hist)

        # Store API messages for multi-turn
        api_msgs <- list(
          list(role = "user", content = prompt),
          list(role = "assistant", content = result$text)
        )
        api_messages(api_msgs)

        show_preview(FALSE)
        is_processing(FALSE)

        showNotification(
          paste0("AI analysis complete! (", result$tokens_out, " tokens)"),
          type = "message", duration = 5)

        # Auto-enable pathway figures if enrichment results exist
        if (!is.null(values$last_enrich_results)) {
          updateCheckboxInput(session, "show_pathway_figs", value = TRUE)
        }

        incProgress(1, detail = "Done!")
      })
    })

    # ---- AI Pathway figure renderers ----
    output$ai_ora_dot_plot <- renderPlotly({
      req(input$show_pathway_figs)
      enrich <- values$last_enrich_results
      req(enrich)
      p <- tryCatch({
        if (!is.null(enrich$ora$up) && is.data.frame(enrich$ora$up) && nrow(enrich$ora$up) > 0) {
          create_ora_dot_plot(enrich$ora$up, "ORA - Upregulated Pathways")
        } else if (!is.null(enrich$ora$down) && is.data.frame(enrich$ora$down) && nrow(enrich$ora$down) > 0) {
          create_ora_dot_plot(enrich$ora$down, "ORA - Downregulated Pathways")
        } else {
          plot_ly() %>% layout(title = "No ORA results available")
        }
      }, error = function(e) {
        plot_ly() %>% layout(title = paste("ORA plot error:", e$message))
      })
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Analysis", figure_type = "ORA Dot Plot",
                       label = "ORA Dot Plot", source_module = "AI Analysis")
      }
      p
    })

    output$ai_gsea_plot <- renderPlotly({
      req(input$show_pathway_figs)
      enrich <- values$last_enrich_results
      req(enrich)
      p <- tryCatch({
        if (!is.null(enrich$gsea) && is.data.frame(enrich$gsea) && nrow(enrich$gsea) > 0) {
          create_gsea_plot(enrich$gsea)
        } else {
          plot_ly() %>% layout(title = "No GSEA results available")
        }
      }, error = function(e) {
        plot_ly() %>% layout(title = paste("GSEA plot error:", e$message))
      })
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "AI Analysis", figure_type = "GSEA Plot",
                       label = "GSEA Enrichment", source_module = "AI Analysis")
      }
      p
    })

    # ---- Follow-up ----
    observeEvent(input$ai_send_followup, {
      api_key <- values$ai_api_key
      followup <- input$ai_followup_text

      if (is.null(api_key) || !nzchar(trimws(api_key))) {
        showNotification("Please set your Claude API key in Data Management > AI Config first.", type = "error")
        return()
      }
      if (is.null(followup) || !nzchar(trimws(followup))) {
        showNotification("Please type a follow-up question.", type = "error")
        return()
      }
      if (length(api_messages()) == 0) {
        showNotification("Run an initial analysis first before sending follow-ups.", type = "warning")
        return()
      }

      is_processing(TRUE)

      withProgress(message = "Processing follow-up...", value = 0, {
        incProgress(0.2, detail = "Sending follow-up question...")

        # Append follow-up to existing conversation
        msgs <- api_messages()
        msgs[[length(msgs) + 1]] <- list(role = "user", content = followup)

        ai_model <- values$ai_model %||% get_default_model(values$ai_api_key)
        result <- call_claude_api(
          api_key = trimws(api_key),
          model = ai_model,
          messages = msgs,
          max_tokens = as.integer(input$ai_max_tokens),
          timeout = 240
        )

        if (!result$ok) {
          showNotification(result$error, type = "error", duration = 12)
          is_processing(FALSE)
          return()
        }

        incProgress(0.9, detail = "Rendering response...")

        # Update chat history
        hist <- chat_history()
        hist[[length(hist) + 1]] <- list(
          role = "user", content = followup,
          timestamp = Sys.time(), tokens = result$tokens_in)
        hist[[length(hist) + 1]] <- list(
          role = "assistant", content = result$text,
          timestamp = Sys.time(), tokens = result$tokens_out)
        chat_history(hist)

        # Update API messages
        msgs[[length(msgs) + 1]] <- list(role = "assistant", content = result$text)
        api_messages(msgs)

        # Log to cross-module conversation history
        values$ai_conversation_log <- c(values$ai_conversation_log, list(
          list(role = "user", content = followup, timestamp = Sys.time(), module = "ai_analysis"),
          list(role = "assistant", content = result$text, timestamp = Sys.time(), module = "ai_analysis")
        ))

        # Clear follow-up input
        updateTextAreaInput(session, "ai_followup_text", value = "")

        is_processing(FALSE)
        showNotification("Follow-up complete!", type = "message", duration = 3)
        incProgress(1, detail = "Done!")
      })
    })

    # ---- Render chat area ----
    output$ai_chat_area <- renderUI({
      hist <- chat_history()
      if (length(hist) == 0) {
        return(fluidRow(column(12,
          div(class = "ai-section-card",
            style = "text-align: center; padding: 40px; color: #999;",
            icon("robot", style = "font-size: 48px; margin-bottom: 15px; color: #ccc;"),
            h4(style = "color: #888;", "No analysis yet"),
            p("Select an analysis mode and click 'Generate' to start."),
            p(style = "font-size: 12px;",
              "Tip: Run DE, pathway, and CCC analyses first to provide richer context for the AI.")
          )
        )))
      }

      # Render each message
      msg_elements <- lapply(seq_along(hist), function(i) {
        msg <- hist[[i]]
        if (msg$role == "user") {
          div(class = "ai-chat-msg ai-chat-user",
            div(style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:5px;",
              span(style = "font-weight:600; color:#1976d2;", "You"),
              span(style = "color:#999; font-size:11px;",
                if (!is.null(msg$timestamp)) format(msg$timestamp, "%H:%M:%S") else "")
            ),
            p(style = "margin:0; white-space: pre-wrap;", msg$content)
          )
        } else {
          html_content <- tryCatch(.mdtohtml(msg$content), error = function(e) {
            gsub("\n", "<br>", msg$content)
          })

          token_info <- if (!is.null(msg$tokens) && !is.na(msg$tokens)) {
            paste0(" (", format(msg$tokens, big.mark = ","), " tokens)")
          } else ""

          div(class = "ai-chat-msg ai-chat-assistant",
            div(style = paste0("display:flex; align-items:center; gap:8px; padding:8px 14px; ",
                               "background:linear-gradient(135deg,#667eea 0%,#764ba2 100%); ",
                               "border-radius:10px 10px 0 0; margin:-12px -14px 12px -14px;"),
              tags$i(class = "fa fa-robot", style = "color:#fff; font-size:14px;"),
              span(style = "font-weight:700; color:#fff; font-size:13px;", "AI Assistant"),
              span(style = "color:rgba(255,255,255,0.75); font-size:11px; margin-left:auto;",
                paste0(if (!is.null(msg$timestamp)) format(msg$timestamp, "%H:%M:%S") else "",
                       token_info))
            ),
            div(class = "ai-thread-content", HTML(html_content))
          )
        }
      })

      fluidRow(column(12,
        div(class = "ai-section-card",
          style = "max-height: 75vh; overflow-y: auto;",
          tagList(msg_elements)
        )
      ))
    })

    # ---- Download HTML report ----
    output$ai_download_html <- downloadHandler(
      filename = function() {
        paste0("scRNA_AI_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".html")
      },
      content = function(file) {
        tryCatch({
          hist <- chat_history()
          if (length(hist) == 0) {
            writeLines("<html><body><p>No analysis results to export.</p></body></html>", file)
            return()
          }

          # Collect all assistant messages
          assistant_msgs <- hist[sapply(hist, function(m) m$role == "assistant")]
          if (length(assistant_msgs) == 0) {
            writeLines("<html><body><p>No AI responses to export.</p></body></html>", file)
            return()
          }

          # Build HTML body
          body_parts <- lapply(seq_along(hist), function(i) {
            msg <- hist[[i]]
            if (msg$role == "user") {
              paste0('<div style="background:#e3f2fd;padding:10px 15px;margin:10px 0;border-radius:8px;border-left:4px solid #1976d2;">',
                '<p style="font-weight:600;color:#1976d2;margin:0 0 5px 0;">You</p>',
                '<p style="margin:0;">', htmltools::htmlEscape(msg$content), '</p></div>')
            } else {
              html_body <- tryCatch(.mdtohtml(msg$content),
                error = function(e) gsub("\n", "<br>", htmltools::htmlEscape(msg$content)))
              paste0('<div style="background:#f5f5f5;padding:15px 20px;margin:10px 0;border-radius:8px;border-left:4px solid #4caf50;">',
                '<p style="font-weight:600;color:#2e7d32;margin:0 0 8px 0;">AI Assistant</p>',
                html_body, '</div>')
            }
          })

          html_doc <- paste0(
            '<!DOCTYPE html><html><head><meta charset="utf-8">',
            '<title>scRNA-seq AI Analysis Report</title>',
            '<style>',
            'body{font-family:"Segoe UI",Helvetica,Arial,sans-serif;max-width:1000px;margin:0 auto;padding:30px;line-height:1.65;color:#333;}',
            'h1{color:#1a237e;border-bottom:3px solid #3c8dbc;padding-bottom:10px;}',
            'h2{color:#283593;margin-top:25px;border-bottom:1px solid #e0e0e0;padding-bottom:5px;}',
            'h3{color:#1565c0;margin-top:20px;}',
            'h4{color:#1976d2;}',
            'table{border-collapse:collapse;width:100%;margin:15px 0;}',
            'th{background:#3c8dbc;color:white;padding:8px 12px;text-align:left;border:1px solid #2980b9;}',
            'td{padding:8px 12px;border:1px solid #ddd;}',
            'tr:nth-child(even){background:#f8f9fa;}',
            'code{background:#f4f4f4;padding:2px 6px;border-radius:3px;font-size:0.9em;}',
            'blockquote{border-left:4px solid #3c8dbc;margin:15px 0;padding:10px 20px;background:#f8f9fa;}',
            'ul,ol{margin:5px 0;}',
            'li{margin:3px 0;}',
            '@media print{body{max-width:100%;} .no-print{display:none;}}',
            '</style></head><body>',
            '<h1>scRNA-seq AI Analysis Report</h1>',
            '<p style="color:#888;">Generated: ', format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
            ' | Model: ', input$ai_model, '</p><hr>',
            paste(body_parts, collapse = "\n"),
            '<hr><p style="color:#aaa;font-size:11px;text-align:center;">',
            'Generated by scRNA4 AI Analysis Assistant</p>',
            '</body></html>'
          )

          writeLines(html_doc, file)
        }, error = function(e) {
          writeLines(paste("<html><body><p>Error generating report:", e$message, "</p></body></html>"), file)
        })
      }
    )

    # ---- Download Markdown ----
    output$ai_download_md <- downloadHandler(
      filename = function() {
        paste0("scRNA_AI_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".md")
      },
      content = function(file) {
        tryCatch({
          hist <- chat_history()
          if (length(hist) == 0) {
            writeLines("# No analysis results to export", file)
            return()
          }

          md_parts <- c(
            paste0("# scRNA-seq AI Analysis Report\n"),
            paste0("*Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                   " | Model: ", input$ai_model, "*\n"),
            "---\n"
          )

          for (msg in hist) {
            if (msg$role == "user") {
              md_parts <- c(md_parts, paste0("\n## Query\n\n> ", msg$content, "\n"))
            } else {
              md_parts <- c(md_parts, paste0("\n", msg$content, "\n\n---\n"))
            }
          }

          writeLines(paste(md_parts, collapse = "\n"), file)
        }, error = function(e) {
          writeLines(paste("# Error generating report\n\n", e$message), file)
        })
      }
    )
  })
}
