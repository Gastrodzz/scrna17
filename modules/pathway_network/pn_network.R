# modules/pathway_network/pn_network.R
# Network analysis engines: Gene Regulatory Networks (WGCNA + correlation),
# co-expression, pathway-network (Jaccard) construction + visualizations.
#
# L-R interaction logic from old network_analysis.R is intentionally NOT
# salvaged — modules/downstream/downstream_ccc.R (CellChat) is the
# authoritative CCC implementation.
#
# Depends on: modules/core_robustness.R (safe_get_assay_data),
#             modules/pathway_network/pn_core.R (sanitiser).

# ---------------------------------------------------------------------------
# GRN — top-level dispatcher (WGCNA or correlation).
# ---------------------------------------------------------------------------

#' @param seurat_obj Seurat object.
#' @param params list with grn_method, grn_cell_type ("all" or a cell-type label),
#'   grn_top_genes, grn_min_cells, plus method-specific tuning fields.
#' @return list(success, network, modules, method, genes_used, error).
build_gene_regulatory_network <- function(seurat_obj, params) {
  tryCatch({
    obj_subset <- if (!identical(params$grn_cell_type, "all") &&
                      !is.null(params$grn_cell_type) &&
                      nzchar(params$grn_cell_type)) {
      ct_col <- params$grn_cell_type_col %||% .pn_detect_cell_type_col(seurat_obj)
      if (is.null(ct_col)) stop("Could not auto-detect a cell-type column.")
      keep <- rownames(seurat_obj@meta.data)[
        as.character(seurat_obj@meta.data[[ct_col]]) == params$grn_cell_type]
      if (length(keep) < 10)
        stop("Cell-type subset has fewer than 10 cells (n=", length(keep), ").")
      subset(seurat_obj, cells = keep)
    } else {
      seurat_obj
    }

    if (length(Seurat::VariableFeatures(obj_subset)) == 0) {
      obj_subset <- Seurat::FindVariableFeatures(
        obj_subset, nfeatures = params$grn_top_genes %||% 2000L, verbose = FALSE)
    }
    hvgs <- utils::head(Seurat::VariableFeatures(obj_subset),
                        params$grn_top_genes %||% 2000L)

    expr <- safe_get_assay_data(obj_subset, "data")[hvgs, , drop = FALSE]
    keep_genes <- rownames(expr)[
      Matrix::rowSums(expr > 0) >= (params$grn_min_cells %||% 10L)]
    if (length(keep_genes) < 20)
      stop("Too few genes pass min_cells filter (n=", length(keep_genes), ").")
    expr <- expr[keep_genes, , drop = FALSE]

    method <- params$grn_method %||% "wgcna"
    res <- switch(method,
      "wgcna"       = .pn_run_wgcna(expr, params, obj_subset),
      "correlation" = .pn_run_correlation_grn(expr, params),
      stop("Unsupported GRN method: ", method,
           ". Use 'wgcna' or 'correlation'.")
    )

    list(success = TRUE, network = res$network, modules = res$modules,
         MEs = res$MEs %||% NULL, power = res$power %||% NA,
         method = method, genes_used = keep_genes, error = NULL)
  }, error = function(e) {
    list(success = FALSE, network = NULL, modules = NULL,
         method = params$grn_method %||% "wgcna",
         genes_used = NULL, error = conditionMessage(e))
  })
}

# Local detector — pn_core uses the de_core helpers via a different name; we
# duplicate a tiny one here so pn_network is self-contained.
.pn_detect_cell_type_col <- function(obj) {
  cands <- c("cell_type", "celltype", "cell.type", "annotation",
             "seurat_clusters", "cluster")
  found <- intersect(cands, colnames(obj@meta.data))
  if (length(found) > 0) return(found[1])
  cat_cols <- colnames(obj@meta.data)[vapply(obj@meta.data,
    function(x) is.factor(x) || is.character(x), logical(1))]
  if (length(cat_cols) > 0) cat_cols[1] else NULL
}

# ---------------------------------------------------------------------------
# WGCNA — pseudobulk-aggregated, memory-bounded.
# ---------------------------------------------------------------------------
.pn_run_wgcna <- function(expr, params, seurat_obj = NULL) {
  if (!requireNamespace("WGCNA", quietly = TRUE))
    stop("WGCNA not installed. Install: install.packages('WGCNA').")

  # Cap genes for memory safety (TOM is genes^2 dense)
  max_genes <- min(nrow(expr), params$wgcna_max_genes %||% 2000L)
  if (nrow(expr) > max_genes) {
    gene_var <- apply(expr, 1, stats::var)
    top <- names(sort(gene_var, decreasing = TRUE))[seq_len(max_genes)]
    expr <- expr[top, , drop = FALSE]
  }
  est_mem_gb <- (nrow(expr)^2 * 8) / (1024^3)
  if (est_mem_gb > 4)
    stop(sprintf("Estimated TOM memory: %.1f GB (max 4 GB). Reduce gene count.",
                 est_mem_gb))

  # Pseudobulk aggregation: WGCNA assumes sample independence — single-cell
  # data violates this. Aggregate by sample x cluster when possible.
  if (!is.null(seurat_obj)) {
    cluster_col <- grep("seurat_clusters|cluster", colnames(seurat_obj@meta.data),
                        value = TRUE, ignore.case = TRUE)[1]
    sample_col  <- grep("orig.ident|sample|donor|batch",
                        colnames(seurat_obj@meta.data),
                        value = TRUE, ignore.case = TRUE)[1]
    if (is.na(sample_col) || length(unique(seurat_obj@meta.data[[sample_col]])) < 3) {
      stop("WGCNA requires ≥3 biological replicates (pseudobulk). ",
           "Single-cell-level WGCNA produces spurious modules. ",
           "Consider hdWGCNA for true single-cell network inference.")
    }
    if (!is.na(cluster_col) && !is.na(sample_col)) {
      meta <- seurat_obj@meta.data
      groups <- paste(meta[[sample_col]], meta[[cluster_col]], sep = "_")
      full_expr <- safe_get_assay_data(seurat_obj, "data")
      # Restrict to genes already chosen
      full_expr <- full_expr[intersect(rownames(full_expr), rownames(expr)), ]
      uniq <- unique(groups)
      pb <- vapply(uniq, function(g) {
        cells <- which(groups == g)
        if (length(cells) > 1) Matrix::rowMeans(full_expr[, cells, drop = FALSE])
        else as.numeric(full_expr[, cells])
      }, numeric(nrow(full_expr)))
      colnames(pb) <- uniq
      expr <- pb
    }
  }

  expr_t <- t(as.matrix(expr))   # WGCNA: samples=rows

  if (isTRUE(params$wgcna_detect_outliers)) {
    gsg <- WGCNA::goodSamplesGenes(expr_t, verbose = 0)
    if (!gsg$allOK) expr_t <- expr_t[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
  }

  power <- params$wgcna_power %||% 0L
  sft <- NULL
  if (power == 0) {
    sft <- WGCNA::pickSoftThreshold(expr_t, powerVector = 1:20, verbose = 0)
    power <- sft$powerEstimate %||% 6L
    if (is.na(power)) power <- 6L
  }

  tom_type <- params$wgcna_tom_type %||% "signed"
  tom_file <- tempfile(pattern = "wgcna_tom_", fileext = "")
  net <- WGCNA::blockwiseModules(
    expr_t, power = power, TOMType = tom_type,
    minModuleSize = params$wgcna_min_module_size %||% 30L,
    reassignThreshold = 0,
    mergeCutHeight = params$wgcna_merge_threshold %||% 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    saveTOMs = TRUE, saveTOMFileBase = tom_file, verbose = 0
  )

  tom_filename <- paste0(tom_file, "-block.1.RData")
  TOM <- tryCatch({ load(tom_filename); TOM },
                  error = function(e) WGCNA::TOMsimilarityFromExpr(
                    expr_t, power = power, TOMType = tom_type))
  tryCatch(file.remove(tom_filename), error = function(e) NULL)

  threshold <- stats::quantile(TOM[upper.tri(TOM)], 0.95)
  pairs <- which(TOM > threshold, arr.ind = TRUE)
  pairs <- pairs[pairs[, 1] < pairs[, 2], , drop = FALSE]
  if (nrow(pairs) > 5000) {
    w <- TOM[pairs]
    pairs <- pairs[order(w, decreasing = TRUE)[seq_len(5000)], , drop = FALSE]
  }
  edges <- data.frame(
    from = colnames(expr_t)[pairs[, 1]],
    to   = colnames(expr_t)[pairs[, 2]],
    weight = TOM[pairs],
    stringsAsFactors = FALSE
  )
  rm(TOM, expr_t); gc()

  net_g <- igraph::graph_from_data_frame(edges, directed = FALSE)
  modules_df <- data.frame(gene = names(net$colors),
                           module = net$colors,
                           stringsAsFactors = FALSE)

  list(network = net_g, modules = modules_df, MEs = net$MEs,
       power = power, sft = sft)
}

# ---------------------------------------------------------------------------
# Correlation network — Spearman / Pearson / Kendall.
# ---------------------------------------------------------------------------
.pn_run_correlation_grn <- function(expr, params) {
  # Cap genes
  if (nrow(expr) > 2000L) {
    gv <- apply(expr, 1, stats::var)
    expr <- expr[names(sort(gv, decreasing = TRUE))[seq_len(2000L)], , drop = FALSE]
  }
  # Cap cells
  if (ncol(expr) > 5000L) {
    set.seed(as.integer(params$seed %||% 42L))
    expr <- expr[, sort(sample.int(ncol(expr), 5000L)), drop = FALSE]
  }
  cor_mat <- stats::cor(t(as.matrix(expr)),
                        method = params$cor_method %||% "spearman")
  n <- ncol(expr)
  cor_sq <- pmin(cor_mat^2, 1 - 1e-10)
  t_stat <- cor_mat * sqrt((n - 2) / (1 - cor_sq))
  p <- 2 * stats::pt(-abs(t_stat), n - 2)

  sig <- (abs(cor_mat) > (params$cor_threshold %||% 0.7)) &
         (p < (params$cor_p_threshold %||% 0.05))
  sig[is.na(sig)] <- FALSE
  diag(sig) <- FALSE
  idx <- which(sig, arr.ind = TRUE)
  idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
  edges <- data.frame(
    from = rownames(cor_mat)[idx[, 1]],
    to   = colnames(cor_mat)[idx[, 2]],
    weight = abs(cor_mat[idx]),
    correlation = cor_mat[idx],
    p_value = p[idx],
    stringsAsFactors = FALSE
  )
  rm(cor_mat, p, t_stat); gc()
  list(network = igraph::graph_from_data_frame(edges, directed = FALSE),
       modules = NULL)
}

# ---------------------------------------------------------------------------
# Co-expression — user-supplied gene list, pairwise correlation.
# ---------------------------------------------------------------------------
build_coexpression_network <- function(seurat_obj, params) {
  tryCatch({
    genes <- unique(trimws(unlist(strsplit(
      params$coexp_genes %||% "", "[,\\s]+"))))
    genes <- genes[nzchar(genes) & genes %in% rownames(seurat_obj)]
    if (length(genes) < 2)
      stop("Need at least 2 valid genes (found ", length(genes), ").")

    expr <- as.matrix(safe_get_assay_data(seurat_obj, "data")[genes, , drop = FALSE])
    cor_mat <- stats::cor(t(expr), method = params$coexp_method %||% "spearman")
    cor_mat[abs(cor_mat) < (params$coexp_threshold %||% 0.3)] <- 0
    diag(cor_mat) <- 0

    idx <- which(cor_mat != 0, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    if (nrow(idx) > (params$coexp_top_n %||% 200L)) {
      vals <- abs(cor_mat[idx])
      idx <- idx[order(vals, decreasing = TRUE)[seq_len(params$coexp_top_n %||% 200L)], , drop = FALSE]
    }
    edges <- data.frame(
      from = rownames(cor_mat)[idx[, 1]],
      to   = colnames(cor_mat)[idx[, 2]],
      correlation = cor_mat[idx],
      stringsAsFactors = FALSE
    )
    list(success = TRUE,
         network = igraph::graph_from_data_frame(edges, directed = FALSE),
         edges = edges, correlation_matrix = cor_mat,
         genes = genes, error = NULL)
  }, error = function(e) {
    list(success = FALSE, network = NULL, edges = NULL,
         correlation_matrix = NULL, genes = NULL, error = conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Pathway network — Jaccard overlap between enriched pathways.
# ---------------------------------------------------------------------------
build_pathway_network <- function(enrichment_results, params) {
  tryCatch({
    df <- if (!is.null(enrichment_results$ora)) {
      do.call(rbind, c(
        list(if (!is.null(enrichment_results$ora$up))   enrichment_results$ora$up),
        list(if (!is.null(enrichment_results$ora$down)) enrichment_results$ora$down)
      ))
    } else {
      stop("No ORA results in enrichment_results.")
    }
    if (is.null(df) || nrow(df) == 0)
      stop("No enriched pathways available.")
    if (!"geneID" %in% colnames(df))
      stop("ORA result has no geneID column — cannot compute Jaccard.")

    # Limit to top-N for performance
    top_n <- params$pathway_top_n %||% 30L
    df <- utils::head(df, top_n)

    pw_genes <- lapply(seq_len(nrow(df)), function(i)
      unlist(strsplit(as.character(df$geneID[i]), "[/]")))
    names(pw_genes) <- df$ID

    n <- length(pw_genes)
    if (n < 2) stop("Need at least 2 pathways for a network.")
    jacc <- matrix(0, nrow = n, ncol = n,
                   dimnames = list(names(pw_genes), names(pw_genes)))
    overlap_min <- params$pathway_overlap_threshold %||% 3L
    jacc_min    <- params$pathway_jaccard_threshold %||% 0.1
    for (i in seq_len(n - 1)) {
      gi <- pw_genes[[i]]
      for (j in (i + 1):n) {
        gj <- pw_genes[[j]]
        inter <- length(intersect(gi, gj))
        if (inter < overlap_min) next
        uni <- length(union(gi, gj))
        jc <- if (uni > 0) inter / uni else 0
        if (jc >= jacc_min) jacc[i, j] <- jacc[j, i] <- jc
      }
    }
    idx <- which(jacc > 0, arr.ind = TRUE)
    idx <- idx[idx[, 1] < idx[, 2], , drop = FALSE]
    if (nrow(idx) == 0)
      stop("No pathway pairs cleared the overlap / Jaccard thresholds.")
    edges <- data.frame(
      from = rownames(jacc)[idx[, 1]],
      to   = colnames(jacc)[idx[, 2]],
      jaccard_index = jacc[idx],
      stringsAsFactors = FALSE
    )
    g <- igraph::graph_from_data_frame(edges, directed = FALSE)
    comm <- igraph::cluster_louvain(g)
    clusters <- data.frame(
      pathway = names(igraph::membership(comm)),
      cluster = as.numeric(igraph::membership(comm)),
      stringsAsFactors = FALSE
    )
    list(success = TRUE, network = g, network_data = edges,
         clusters = clusters, pathway_genes = pw_genes, error = NULL)
  }, error = function(e) {
    list(success = FALSE, network = NULL, network_data = NULL,
         clusters = NULL, pathway_genes = NULL, error = conditionMessage(e))
  })
}

# ---------------------------------------------------------------------------
# Visualizations — visNetwork and plotly helpers.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# visNetwork physics: stabilize once, then freeze. Without this the simulation
# keeps shaking the layout indefinitely. The pattern below runs ~200 physics
# iterations to settle the node positions then disables physics so the layout
# stays put while the user pans/zooms.
# ---------------------------------------------------------------------------
.pn_static_physics <- function(net) {
  net |>
    visNetwork::visPhysics(
      stabilization = list(enabled = TRUE, iterations = 200,
                           updateInterval = 25, onlyDynamicEdges = FALSE,
                           fit = TRUE),
      enabled = FALSE) |>
    visNetwork::visEdges(smooth = list(enabled = FALSE)) |>
    visNetwork::visInteraction(navigationButtons = TRUE,
                               dragNodes = TRUE,
                               zoomView = TRUE,
                               dragView = TRUE)
}

#' visNetwork rendering of a GRN (or any igraph result with $network).
render_grn_network <- function(grn_results) {
  if (is.null(grn_results) || is.null(grn_results$network) ||
      igraph::vcount(grn_results$network) == 0) {
    return(visNetwork::visNetwork(
      data.frame(id = 1, label = "No network — try different settings"),
      data.frame(from = numeric(0), to = numeric(0))))
  }
  g <- grn_results$network
  deg <- igraph::degree(g)
  bet <- igraph::betweenness(g)
  nodes <- data.frame(
    id = igraph::V(g)$name,
    label = igraph::V(g)$name,
    size  = scales::rescale(deg, to = c(10, 50)),
    color = scales::rescale(bet, to = c(0.3, 1)),
    stringsAsFactors = FALSE
  )
  edges <- igraph::as_data_frame(g, what = "edges")
  edges$color <- "#9aa3af"; edges$width <- 2
  visNetwork::visNetwork(nodes, edges) |>
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
    visNetwork::visLayout(randomSeed = 123) |>
    .pn_static_physics()
}

render_coexp_network <- function(coexp_results) {
  if (is.null(coexp_results) || is.null(coexp_results$network) ||
      igraph::vcount(coexp_results$network) == 0) {
    return(visNetwork::visNetwork(
      data.frame(id = 1, label = "No co-expression network"),
      data.frame(from = numeric(0), to = numeric(0))))
  }
  g <- coexp_results$network
  deg <- igraph::degree(g)
  nodes <- data.frame(
    id = igraph::V(g)$name, label = igraph::V(g)$name,
    size = scales::rescale(deg, to = c(15, 40)),
    color = "#91D1C2", stringsAsFactors = FALSE)
  edges <- igraph::as_data_frame(g, what = "edges")
  edges$width <- scales::rescale(abs(edges$correlation), to = c(1, 8))
  edges$color <- ifelse(edges$correlation > 0, "#E64B35", "#3C5488")
  visNetwork::visNetwork(nodes, edges) |>
    visNetwork::visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) |>
    visNetwork::visLayout(randomSeed = 123) |>
    .pn_static_physics()
}

render_correlation_matrix <- function(coexp_results) {
  if (is.null(coexp_results) || is.null(coexp_results$correlation_matrix))
    return(plotly::plotly_empty())
  m <- coexp_results$correlation_matrix
  plotly::plot_ly(z = ~m, x = colnames(m), y = rownames(m),
                  type = "heatmap", colorscale = "RdBu", zmid = 0) |>
    plotly::layout(title = "Gene correlation matrix")
}

render_pathway_network <- function(pathway_results) {
  if (is.null(pathway_results) || is.null(pathway_results$network) ||
      igraph::vcount(pathway_results$network) == 0) {
    return(visNetwork::visNetwork(
      data.frame(id = 1, label = "No pathway network"),
      data.frame(from = numeric(0), to = numeric(0))))
  }
  g <- pathway_results$network
  comm <- igraph::cluster_louvain(g)
  nodes <- data.frame(
    id = igraph::V(g)$name, label = igraph::V(g)$name,
    group = as.character(igraph::membership(comm)), size = 25,
    stringsAsFactors = FALSE)
  edges <- igraph::as_data_frame(g, what = "edges")
  edges$width <- scales::rescale(edges$jaccard_index, to = c(1, 8))
  visNetwork::visNetwork(nodes, edges) |>
    visNetwork::visOptions(highlightNearest = TRUE) |>
    visNetwork::visLayout(randomSeed = 123) |>
    visNetwork::visLegend() |>
    .pn_static_physics()
}

# Null-coalescing operator
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
