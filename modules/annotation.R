# =============================================================================
# Cell Type Annotation Module (v3.0)
# modules/annotation.R
#
# Automatic annotation via ScType (Ianevski et al., Nat Comms 2022)
# Uses local ScTypeDB_full.xlsx marker database
# AI annotation via Claude API (separate tab)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(ggplot2)
  library(R6)
  library(openxlsx)
})

# ---- Null-coalescing operator (if not yet defined) --------------------------
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

# =============================================================================
# Package availability checks
# =============================================================================

# clustree availability
.clustree_available <- requireNamespace("clustree", quietly = TRUE)

# =============================================================================
# Marker database metadata
# =============================================================================

MARKER_DB_VERSION <- list(
  version  = "2.0.0",
  source   = "Curated from literature + PanglaoDB",
  date     = "2025-06-01",
  note     = "Heuristic marker database with weighted positive/negative markers"
)

# =============================================================================
# Utility helpers
# =============================================================================

# Local alias for core_robustness helpers (guard with exists)
.safe_default_assay <- function(obj) {
  if (exists("safe_default_assay", mode = "function")) safe_default_assay(obj)
  else tryCatch(DefaultAssay(obj), error = function(e) "RNA")
}

# Case-insensitive gene mapping: upper-case key -> original gene name in object
.case_insensitive_map <- function(genes_in_obj) {
  setNames(genes_in_obj, toupper(genes_in_obj))
}

# Min-max normalize each row of a matrix to [0, 1]
.row_minmax01 <- function(m) {
  if (is.null(dim(m))) return(m)
  rmin <- apply(m, 1, min, na.rm = TRUE)
  rmax <- apply(m, 1, max, na.rm = TRUE)
  rng <- pmax(1e-9, rmax - rmin)
  sweep(sweep(m, 1, rmin, FUN = "-"), 1, rng, FUN = "/")
}

# Rank-based normalization: replace values with quantile ranks per row
.row_rank_normalize <- function(m) {
  if (is.null(dim(m))) return(m)
  t(apply(m, 1, function(x) {
    r <- rank(x, ties.method = "average")
    (r - 1) / max(1, length(r) - 1)
  }))
}

# Detect species from gene names
# Species detection now lives in modules/species_utils.R (single source of truth).
# Thin alias keeps existing call sites working; returns just the species label.
.detect_species <- function(genes) {
  sp <- detect_species(genes)$species
  if (sp == "unknown") "human" else sp   # preserve prior human-default behaviour
}

# Sanitize a cell-type label so cluster names stay clean ("just the cell type").
# Strips anything a scorer/AI may append: trailing "(p<0.01)", "score: 0.8",
# "p value ...", extra lines, stray brackets. Vectorised.
.clean_celltype_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- "Unknown"
  x <- sub("[\r\n].*$", "", x)                                  # first line only
  x <- gsub("\\s*[\\(\\[\\{].*$", "", x, perl = TRUE)           # drop trailing (...)/[...]
  # drop trailing score/p-value/confidence fragments (with or without punctuation)
  x <- gsub("(?i)\\s*[,;:|/-]?\\s*(p[ ._-]?val(ue)?|score|adj\\.?p|padj|fdr|q[ ._-]?val(ue)?|conf(idence)?|es[ ._-]?max)\\b.*$",
            "", x, perl = TRUE)
  x <- gsub("\\s+", " ", x)                                     # collapse whitespace
  x <- trimws(x)
  x[!nzchar(x)] <- "Unknown"
  x
}

# Clean the verbose cell-type names used by whole-organism atlases so a cluster
# shows just the cell type. Format differs per atlas:
#   MCA  "B cell_Rps27rt high"                    -> "B cell"        (marker suffix after "_")
#   HCL  "AT2.cell.Adult.Lung2."                  -> "AT2 cell"      (dot-separated stage/tissue)
#        "Vascular.endothelial.cell_VWF.high..."  -> "Vascular endothelial cell"
#   celldex "Epithelial_cells"                    -> "Epithelial cells" (underscore = space)
.clean_atlas_label <- function(x, atlas = c("celldex", "mca", "hcl")) {
  atlas <- match.arg(atlas)
  x <- as.character(x); x[is.na(x)] <- "Unknown"
  if (atlas == "mca") {
    x <- sub("_.*$", "", x)                                     # drop "_Marker high"
    x <- sub("\\s*\\(.*$", "", x)                               # drop "(Tissue)" suffix
  } else if (atlas == "hcl") {
    x <- sub("_.*$", "", x)                                     # drop "_Marker high"
    x <- sub("\\s*\\(.*$", "", x)                               # drop "(Tissue)" suffix
    x <- gsub("\\.", " ", x)                                    # dots -> spaces
    x <- gsub("(?i)\\s+(Adult|Fetal|Neonat\\w*|Newborn|Embryo\\w*|Male|Female)\\b.*$",
              "", x, perl = TRUE)                               # drop stage/sex + trailing
    x <- gsub("\\s+[A-Z][a-z]+[0-9]+\\s*$", "", x, perl = TRUE) # drop trailing Tissue<digit>
  } else {                                                      # celldex
    x <- gsub("_", " ", x)                                      # B_cell -> B cell
  }
  x <- gsub("\\s+", " ", x); x <- trimws(x)
  x[!nzchar(x)] <- "Unknown"
  x
}

# Quick color palette
.get_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    base <- RColorBrewer::brewer.pal(min(12, max(3, n)), "Set3")
    if (n > length(base)) grDevices::colorRampPalette(base)(n) else base[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, "Dynamic")
  }
}

# Best available reduction
.best_reduction <- function(obj, prefer_3d = FALSE) {
  reds <- names(obj@reductions)
  if (prefer_3d && "umap3d" %in% reds) return("umap3d")
  for (r in c("umap", "tsne", "pca")) {
    if (r %in% reds) return(r)
  }
  if (length(reds) > 0) return(reds[1])
  NULL
}

# Safe prepare_for_de wrapper
.safe_prepare_for_de <- function(obj) {
  if (exists("prepare_for_de", mode = "function")) {
    tryCatch(prepare_for_de(obj), error = function(e) obj)
  } else {
    obj
  }
}

# Safe should_sketch wrapper
.safe_should_sketch <- function(values, budget = 50000L) {
  if (exists("should_sketch", mode = "function")) {
    tryCatch(should_sketch(values, budget = budget), error = function(e) FALSE)
  } else {
    FALSE
  }
}

# Safe sketch_downsample wrapper
.safe_sketch_downsample <- function(obj, max_cells = 50000L, method = "geometric") {
  if (exists("sketch_downsample", mode = "function")) {
    tryCatch(sketch_downsample(obj, max_cells = max_cells, method = method),
             error = function(e) obj)
  } else {
    obj
  }
}

# Safe get assay data wrapper
.safe_get_assay_data <- function(obj, layer = "data", assay = NULL) {
  if (exists("safe_get_assay_data", mode = "function")) {
    safe_get_assay_data(obj, layer, assay = assay)
  } else {
    tryCatch({
      a <- assay %||% DefaultAssay(obj)
      GetAssayData(obj, layer = layer, assay = a)
    }, error = function(e) {
      tryCatch(GetAssayData(obj, slot = layer), error = function(e2) NULL)
    })
  }
}

# =============================================================================
# AI prompt builder (exported helper)
# =============================================================================

build_ai_annotation_prompt_structured <- function(cluster_markers_df, tissue, species) {
  clusters <- unique(cluster_markers_df$cluster)
  marker_lines <- sapply(clusters, function(cl) {
    sub_df <- cluster_markers_df[cluster_markers_df$cluster == cl &
                                   !is.na(cluster_markers_df$p_val_adj) &
                                   cluster_markers_df$p_val_adj < 0.05, ]
    sub_df <- sub_df[order(-sub_df$avg_log2FC), ]
    top <- head(sub_df, 15)
    paste0("Cluster ", cl, ": ",
           paste0(top$gene, "(FC=", round(top$avg_log2FC, 2), ")"),
           collapse = ", ")
  })

  paste0(
    "You are annotating scRNA-seq clusters from ", species,
    " ", tissue, " tissue.\n\n",
    "Respond with ONLY a valid JSON array. No markdown, no explanation.\n",
    "Each element must have exactly these fields:\n",
    '{"cluster":"<id>","cell_type":"<name>","subtype":"<name or null>",',
    '"confidence":"high|medium|low",',
    '"key_markers":["GENE1","GENE2","GENE3"],',
    '"notes":"<one sentence>"}\n\n',
    "Clusters and their top markers:\n",
    paste(marker_lines, collapse = "\n"), "\n\n",
    "IMPORTANT: Output ONLY the JSON array starting with [ and ending with ]."
  )
}

# =============================================================================
# Enhanced marker-based annotator (R6 class)
# =============================================================================

EnhancedCellTypeAnnotator <- R6::R6Class(
  "EnhancedCellTypeAnnotator",
  public = list(
    marker_database = NULL,

    initialize = function(marker_db = NULL) {
      self$marker_database <- marker_db %||% (if (exists("CELL_TYPE_MARKERS")) CELL_TYPE_MARKERS else NULL)
      if (is.null(self$marker_database)) stop("Marker database is NULL.")
    },

    # ------------------------------------------------------------------
    # Multi-signal scoring engine
    # Returns: list(scores, confidence, weighted_scores, pct_scores,
    #               specificity_scores)
    # ------------------------------------------------------------------
    score_cell_types_enhanced = function(seurat_obj, assay = NULL) {
      assay <- assay %||% .safe_default_assay(seurat_obj)
      expr <- .safe_get_assay_data(seurat_obj, "data", assay = assay)
      if (is.null(expr)) stop("Cannot retrieve expression data from assay '", assay, "'")
      if (!inherits(expr, "dgCMatrix")) expr <- as(as.matrix(expr), "dgCMatrix")

      genes_map <- .case_insensitive_map(rownames(expr))
      cell_types <- names(self$marker_database)

      n_cells <- ncol(expr)
      n_ct    <- length(cell_types)
      ws   <- matrix(0, nrow = n_cells, ncol = n_ct,
                      dimnames = list(colnames(expr), cell_types))
      pcs  <- ws   # detection breadth
      spec <- ws   # specificity
      neg  <- ws   # negative marker penalty

      # --- Gene specificity: inverse of how many cell types share each marker ---
      all_markers <- unique(toupper(unlist(lapply(self$marker_database, function(m) names(m)))))
      gene_ct_count <- setNames(integer(length(all_markers)), all_markers)
      for (ct in cell_types) {
        mk <- self$marker_database[[ct]]
        if (is.null(mk)) next
        for (g in toupper(names(mk))) {
          if (g %in% names(gene_ct_count)) gene_ct_count[g] <- gene_ct_count[g] + 1L
        }
      }

      # --- Per-gene global statistics ---
      gene_means <- Matrix::rowMeans(expr)
      gene_sds   <- sqrt(pmax(0, Matrix::rowMeans(expr^2) - gene_means^2))

      total_matches    <- 0L
      n_markers_per_ct <- integer(n_ct)

      for (i in seq_along(cell_types)) {
        ct <- cell_types[i]
        mk <- self$marker_database[[ct]]
        if (is.null(mk) || !length(mk)) next

        names(mk) <- toupper(names(mk))
        keep <- names(mk) %in% names(genes_map)
        if (!any(keep)) next

        idx_genes <- genes_map[names(mk)[keep]]
        w         <- as.numeric(mk[keep])
        mk_names  <- names(mk)[keep]
        n_markers_per_ct[i] <- sum(keep)
        total_matches       <- total_matches + length(w)

        # BUG FIX #4: Relaxed detection filter — use 1% instead of 5%
        # and never filter out the LAST remaining markers
        detection_rate <- Matrix::rowMeans(expr[idx_genes, , drop = FALSE] > 0)
        # Adaptive threshold: 1% for rare types (few markers), 2% for well-covered types
        det_threshold <- if (length(idx_genes) <= 5) 0.005 else 0.01
        reliable <- detection_rate >= det_threshold
        if (sum(reliable) >= 2) {
          idx_genes <- idx_genes[reliable]
          w <- w[reliable]
          mk_names <- mk_names[reliable]
          n_markers_per_ct[i] <- sum(reliable)
        }
        # If fewer than 2 pass the filter, keep all — do not discard

        sub <- expr[idx_genes, , drop = FALSE]

        # === SIGNAL 1: Weighted expression score (with soft thresholding) ===
        g_mu  <- gene_means[idx_genes]
        g_sd  <- pmax(gene_sds[idx_genes], 0.1)
        g_median <- apply(sub, 1, function(x) stats::median(x))
        g_mad    <- pmax(0.1, apply(sub, 1, function(x) stats::mad(x, constant = 1.4826)))
        thresh   <- g_median + 1.0 * g_mad

        if (nrow(sub) == 1) {
          expr_vals <- as.numeric(sub)
          above <- pmin(1, pmax(0, (expr_vals - thresh[1] * 0.5) / pmax(0.1, thresh[1])))
          ws[, i] <- expr_vals * w * above
        } else {
          thresh_mat <- matrix(thresh, nrow = length(thresh), ncol = n_cells)
          above <- pmin(1, pmax(0, (sub - thresh_mat * 0.5) / pmax(0.1, thresh_mat)))
          dampened <- sub * above
          ws[, i] <- as.numeric(crossprod(w, dampened))
        }

        # === SIGNAL 2: Detection breadth ===
        pos_mask <- w > 0
        if (any(pos_mask)) {
          bin_pos <- sub[pos_mask, , drop = FALSE] > 0
          w_pos   <- abs(w[pos_mask])
          denom_pos <- sum(w_pos)
          if (sum(pos_mask) == 1) {
            pcs[, i] <- as.numeric(bin_pos) * w_pos / denom_pos
          } else {
            pcs[, i] <- as.numeric(crossprod(w_pos, bin_pos)) / denom_pos
          }
        }

        # === SIGNAL 3: Specificity-weighted score ===
        spec_w <- ifelse(w > 0,
          w * 1 / pmax(1, gene_ct_count[mk_names]),
          0)
        if (nrow(sub) == 1) {
          spec[, i] <- as.numeric(sub > 0) * spec_w
        } else {
          bin_all <- sub > 0
          spec[, i] <- as.numeric(crossprod(spec_w, bin_all))
        }

        # === SIGNAL 4: Negative marker penalty ===
        neg_mask <- w < 0
        if (any(neg_mask)) {
          neg_sub <- sub[neg_mask, , drop = FALSE]
          neg_w   <- abs(w[neg_mask])
          if (sum(neg_mask) == 1) {
            neg[, i] <- as.numeric(neg_sub > 0) * neg_w
          } else {
            neg[, i] <- as.numeric(crossprod(neg_w, neg_sub > 0))
          }
        }
      }

      if (total_matches == 0) stop("No marker genes matched the dataset.")

      # --- Normalize each signal ---
      ws_n   <- .row_minmax01(ws)
      pcs_n  <- .row_minmax01(pcs)
      spec_n <- .row_minmax01(spec)

      # BUG FIX #4: Also apply rank-based normalization and blend
      ws_rank   <- .row_rank_normalize(ws)
      pcs_rank  <- .row_rank_normalize(pcs)
      spec_rank <- .row_rank_normalize(spec)

      # Blend min-max and rank normalization (50/50)
      ws_blend   <- 0.5 * ws_n + 0.5 * ws_rank
      pcs_blend  <- 0.5 * pcs_n + 0.5 * pcs_rank
      spec_blend <- 0.5 * spec_n + 0.5 * spec_rank

      neg_n <- neg
      for (i in seq_len(n_ct)) {
        mx <- max(neg[, i], na.rm = TRUE)
        if (mx > 0) neg_n[, i] <- neg[, i] / mx
      }

      # --- Combined score ---
      # Expression (35%) + Detection breadth (25%) + Specificity (25%) + Rank-boost (5%) - Negative (10%)
      combined <- 0.35 * ws_blend + 0.25 * pcs_blend + 0.25 * spec_blend +
                  0.05 * ws_rank - 0.10 * neg_n
      combined[combined < 0] <- 0

      # --- Coverage penalty ---
      max_markers <- max(n_markers_per_ct, na.rm = TRUE)
      for (i in seq_len(n_ct)) {
        if (n_markers_per_ct[i] < 2) {
          combined[, i] <- combined[, i] * (n_markers_per_ct[i] / 2) * 0.4
        } else if (n_markers_per_ct[i] < 3) {
          combined[, i] <- combined[, i] * 0.7
        } else if (n_markers_per_ct[i] < max_markers * 0.3) {
          combined[, i] <- combined[, i] * 0.85
        }
      }

      combined <- .row_minmax01(combined)

      # --- Confidence scoring (entropy-based) ---
      best <- apply(combined, 1, max, na.rm = TRUE)
      ord3 <- t(apply(combined, 1, function(x) {
        s <- sort(x, decreasing = TRUE)
        c(s[1], if (length(s) > 1) s[2] else 0, if (length(s) > 2) s[3] else 0)
      }))
      margin_1_2 <- ord3[, 1] - ord3[, 2]

      top3_sum <- ord3[, 1] + ord3[, 2] + ord3[, 3] + 1e-10
      p1 <- ord3[, 1] / top3_sum
      p2 <- ord3[, 2] / top3_sum
      p3 <- ord3[, 3] / top3_sum
      entropy <- -(p1 * log2(p1 + 1e-10) + p2 * log2(p2 + 1e-10) + p3 * log2(p3 + 1e-10))
      entropy_norm <- 1 - pmin(1, entropy / 1.585)

      conf <- 0.35 * best + 0.35 * margin_1_2 + 0.30 * entropy_norm

      list(
        scores = combined,
        confidence = data.frame(
          best_score = best,
          score_diff = margin_1_2,
          overall    = conf,
          row.names  = rownames(combined)
        ),
        weighted_scores     = ws,
        pct_scores          = pcs,
        specificity_scores  = spec
      )
    },

    # ------------------------------------------------------------------
    # Cell-level annotation
    # ------------------------------------------------------------------
    annotate_cells_with_confidence = function(seurat_obj, scores, confidence,
                                              min_conf = 0.25,
                                              min_raw_score = NULL,
                                              column = "cell_type_enhanced") {
      if (is.null(min_raw_score)) {
        all_max <- apply(scores, 1, max)
        nonzero <- all_max[all_max > 0]
        min_raw_score <- if (length(nonzero) > 10) quantile(nonzero, 0.25) else 0.1
      }
      score_range <- max(scores) - min(scores)
      if (score_range < 0.05) {
        min_raw_score <- max(scores) + 1
      }
      labs <- unname(colnames(scores)[max.col(scores, ties.method = "first")])
      raw_max <- apply(scores, 1, max)
      labs[raw_max < min_raw_score] <- "Unknown"
      labs[confidence$overall < min_conf] <- "Unknown"
      seurat_obj[[column]] <- unname(labs)
      seurat_obj[["annotation_score"]]      <- unname(apply(scores, 1, max))
      seurat_obj[["annotation_confidence"]] <- unname(confidence$overall)
      seurat_obj[["score_difference"]]      <- unname(confidence$score_diff)
      seurat_obj
    },

    # ------------------------------------------------------------------
    # Aggregate cell-level scores by cluster
    # ------------------------------------------------------------------
    aggregate_scores_by_cluster = function(scores, confidence, cluster_labels) {
      cluster_labels <- as.character(cluster_labels)
      clusters <- sort(unique(cluster_labels))
      n_ct <- ncol(scores)

      agg_scores <- matrix(0, nrow = length(clusters), ncol = n_ct,
                            dimnames = list(clusters, colnames(scores)))
      for (cl in clusters) {
        mask <- cluster_labels == cl
        if (sum(mask) == 1) {
          agg_scores[cl, ] <- scores[mask, ]
        } else {
          agg_scores[cl, ] <- colMeans(scores[mask, , drop = FALSE], na.rm = TRUE)
        }
      }

      agg_scores <- .row_minmax01(agg_scores)

      best <- apply(agg_scores, 1, max, na.rm = TRUE)
      ord3 <- t(apply(agg_scores, 1, function(x) {
        s <- sort(x, decreasing = TRUE)
        c(s[1], if (length(s) > 1) s[2] else 0, if (length(s) > 2) s[3] else 0)
      }))
      margin_1_2 <- ord3[, 1] - ord3[, 2]
      top3_sum <- ord3[, 1] + ord3[, 2] + ord3[, 3] + 1e-10
      p1 <- ord3[, 1] / top3_sum
      p2 <- ord3[, 2] / top3_sum
      p3 <- ord3[, 3] / top3_sum
      entropy <- -(p1 * log2(p1 + 1e-10) + p2 * log2(p2 + 1e-10) + p3 * log2(p3 + 1e-10))
      entropy_norm <- 1 - pmin(1, entropy / 1.585)
      conf <- 0.35 * best + 0.35 * margin_1_2 + 0.30 * entropy_norm

      agg_confidence <- data.frame(
        best_score = best,
        score_diff = margin_1_2,
        overall    = conf,
        row.names  = clusters
      )

      best_types <- colnames(agg_scores)[max.col(agg_scores, ties.method = "first")]
      names(best_types) <- clusters

      list(
        cluster_scores     = agg_scores,
        cluster_confidence = agg_confidence,
        cluster_labels     = best_types
      )
    },

    # ------------------------------------------------------------------
    # Cluster-level annotation
    # ------------------------------------------------------------------
    annotate_clusters_with_confidence = function(seurat_obj, scores, confidence,
                                                  cluster_labels,
                                                  min_conf = 0.25,
                                                  min_raw_score = NULL,
                                                  column = "cell_type_enhanced") {
      if (is.null(min_raw_score)) {
        all_max <- apply(scores, 1, max)
        nonzero <- all_max[all_max > 0]
        min_raw_score <- if (length(nonzero) > 10) quantile(nonzero, 0.25) else 0.1
      }
      score_range <- max(scores) - min(scores)
      if (score_range < 0.05) {
        min_raw_score <- max(scores) + 1
      }
      agg <- self$aggregate_scores_by_cluster(scores, confidence, cluster_labels)
      cluster_labels <- as.character(cluster_labels)

      cell_labels     <- unname(agg$cluster_labels[cluster_labels])
      cell_conf       <- unname(agg$cluster_confidence$overall[
        match(cluster_labels, rownames(agg$cluster_confidence))
      ])
      cell_best_score <- unname(agg$cluster_confidence$best_score[
        match(cluster_labels, rownames(agg$cluster_confidence))
      ])
      cell_score_diff <- unname(agg$cluster_confidence$score_diff[
        match(cluster_labels, rownames(agg$cluster_confidence))
      ])

      cell_labels[cell_best_score < min_raw_score] <- "Unknown"
      cell_labels[cell_conf < min_conf]            <- "Unknown"

      seurat_obj[[column]]                  <- cell_labels
      seurat_obj[["annotation_score"]]      <- cell_best_score
      seurat_obj[["annotation_confidence"]] <- cell_conf
      seurat_obj[["score_difference"]]      <- cell_score_diff
      seurat_obj
    }
  )
)

# =============================================================================
# Contract validation & job helpers
# =============================================================================

validate_annotation_contract <- function(values, obj) {
  if (!exists("validate_seurat_contract", mode = "function")) return(TRUE)
  contract <- validate_seurat_contract(
    obj,
    stage        = "annotation",
    min_cells    = 20,
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
        module  = "annotation",
        action  = "contract_validation_failed",
        level   = "error",
        details = list(errors = contract$errors, warnings = contract$warnings,
                       suggestions = contract$suggestions)
      )
    }
    return(FALSE)
  }
  if (length(contract$warnings) > 0) {
    showNotification(paste("Contract warnings:", paste(contract$warnings, collapse = " | ")),
                     type = "warning", duration = 8)
  }
  TRUE
}

begin_annotation_job <- function(values, action, meta = list()) {
  if (!exists("start_job", mode = "function")) return(NULL)
  start_job(values, module = "annotation", action = action, meta = meta)
}

end_annotation_job <- function(values, job_id, status = "done", message = NULL, details = list()) {
  if (is.null(job_id) || !exists("finish_job", mode = "function")) return(invisible(NULL))
  finish_job(values, job_id = job_id, status = status, message = message, details = details)
}

checkpoint_annotation <- function(values, action, obj, params = list(), method = NULL) {
  if (!exists("register_dataset_checkpoint", mode = "function")) return(invisible(NULL))
  register_dataset_checkpoint(
    values,
    module  = "annotation",
    action  = action,
    obj     = obj,
    params  = params,
    methods = list(method = method)
  )
}

# =============================================================================
# Annotation disagreement summary (exported helper)
# =============================================================================

annotation_disagreement_summary <- function(meta_df) {
  if (is.null(meta_df) || !"annotation_disagreement" %in% colnames(meta_df)) {
    total <- if (is.null(meta_df)) 0 else nrow(meta_df)
    return(data.frame(
      Metric = c("Total cells", "Disagreement cells", "Disagreement rate (%)"),
      Value  = c(total, 0, 0),
      stringsAsFactors = FALSE
    ))
  }

  total <- nrow(meta_df)
  disagree_vec <- as.logical(meta_df$annotation_disagreement)
  disagree_vec[is.na(disagree_vec)] <- FALSE
  disag <- sum(disagree_vec)
  rate  <- if (total > 0) round(100 * disag / total, 2) else 0
  low_conf <- if ("annotation_confidence" %in% colnames(meta_df)) {
    sum(meta_df$annotation_confidence < 0.25, na.rm = TRUE)
  } else {
    NA_integer_
  }

  data.frame(
    Metric = c("Total cells", "Disagreement cells", "Disagreement rate (%)",
               "Low-confidence cells (<0.25)"),
    Value  = c(total, disag, rate, low_conf),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Shared: FindAllMarkers for annotation (used by both ScType and PanglaoDB)
# =============================================================================

.find_markers_for_annotation <- function(seurat_obj, cluster_col = "seurat_clusters") {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      cluster_col <- "seurat_clusters"
    } else {
      seurat_obj$seurat_clusters <- Idents(seurat_obj)
      cluster_col <- "seurat_clusters"
    }
  }
  Idents(seurat_obj) <- cluster_col
  clusters <- sort(unique(Idents(seurat_obj)))
  if (length(clusters) < 2) {
    stop("Need at least 2 clusters to find markers, found ", length(clusters), ".")
  }

  marker_obj <- seurat_obj
  max_cells <- 200L
  if (ncol(marker_obj) > max_cells * length(clusters)) {
    set.seed(42)
    cells_keep <- c()
    for (cl in clusters) {
      cl_cells <- WhichCells(marker_obj, idents = cl)
      if (length(cl_cells) > max_cells) cl_cells <- sample(cl_cells, max_cells)
      cells_keep <- c(cells_keep, cl_cells)
    }
    marker_obj <- subset(marker_obj, cells = cells_keep)
  }
  marker_obj <- .safe_prepare_for_de(marker_obj)

  markers_df <- tryCatch(
    FindAllMarkers(marker_obj, only.pos = TRUE, min.pct = 0.1,
                   logfc.threshold = 0.25, verbose = FALSE),
    error = function(e) {
      tryCatch(
        FindAllMarkers(marker_obj, only.pos = TRUE, min.pct = 0.01,
                       logfc.threshold = 0.1, verbose = FALSE),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(markers_df) || nrow(markers_df) == 0) {
    stop("FindAllMarkers returned no markers. Ensure data is normalized and clusters exist.")
  }
  markers_df
}

# =============================================================================
# ScType: Automated Cell Type Annotation
# Ianevski et al., Nature Communications 13, 1246 (2022)
# Marker-based approach: finds cluster markers -> matches against ScTypeDB
# =============================================================================

.sctype_db_path <- function() {
  candidates <- c(
    file.path(getwd(), "ScTypeDB_full.xlsx"),
    file.path(dirname(getwd()), "ScTypeDB_full.xlsx"),
    file.path(getwd(), "sc-type-master", "ScTypeDB_full.xlsx")
  )
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NULL
}

.load_sctype_db <- function(db_path = NULL) {
  if (is.null(db_path)) db_path <- .sctype_db_path()
  if (is.null(db_path) || !file.exists(db_path)) {
    stop("ScTypeDB_full.xlsx not found. Place it in the app root directory.")
  }
  db <- openxlsx::read.xlsx(db_path)
  db$markers_pos <- lapply(seq_len(nrow(db)), function(i) {
    raw <- gsub(" ", "", as.character(db$geneSymbolmore1[i]))
    m <- unlist(strsplit(raw, ","))
    m <- toupper(m[m != "" & m != "NA" & !is.na(m)])
    unique(sort(m))
  })
  db$markers_neg <- lapply(seq_len(nrow(db)), function(i) {
    raw <- gsub(" ", "", as.character(db$geneSymbolmore2[i]))
    m <- unlist(strsplit(raw, ","))
    m <- toupper(m[m != "" & m != "NA" & !is.na(m)])
    unique(sort(m))
  })
  # Pre-compute marker specificity: how many cell types each marker appears in
  all_pos <- unlist(db$markers_pos)
  marker_freq <- table(all_pos)
  n_types <- length(unique(db$cellName))
  db$.marker_specificity <- lapply(db$markers_pos, function(genes) {
    if (length(genes) == 0) return(numeric(0))
    freqs <- as.numeric(marker_freq[genes])
    freqs[is.na(freqs)] <- 1
    # Higher weight for markers unique to fewer cell types
    setNames(1 / sqrt(pmax(freqs, 1)), genes)
  })
  db
}

# Improved scoring: FC-weighted overlap * marker specificity, negative penalty, Jaccard normalization
.score_cluster_vs_sctype <- function(cluster_markers_upper, cluster_fc, db) {
  n_cluster <- length(cluster_markers_upper)
  scores <- vapply(seq_len(nrow(db)), function(i) {
    pos <- db$markers_pos[[i]]
    neg <- db$markers_neg[[i]]
    spec <- db$.marker_specificity[[i]]
    n_pos <- length(pos)
    if (n_pos == 0) return(0)

    pos_hits <- intersect(cluster_markers_upper, pos)
    n_hits <- length(pos_hits)
    if (n_hits == 0) return(0)

    # Weighted score: log2FC * specificity for each matching marker
    fc_vals <- cluster_fc[pos_hits]
    sp_vals <- spec[pos_hits]
    score_pos <- sum(fc_vals * sp_vals, na.rm = TRUE)

    # Negative penalty: cluster markers that appear in the cell type's negative list
    neg_hits <- intersect(cluster_markers_upper, neg)
    score_neg <- sum(cluster_fc[neg_hits], na.rm = TRUE) * 0.5

    # Normalize: geometric mean of coverage fraction and raw score
    coverage <- n_hits / n_pos  # fraction of DB markers found
    raw <- (score_pos - score_neg) / sqrt(n_pos)
    raw * (1 + coverage)
  }, numeric(1))

  n_overlap <- vapply(seq_len(nrow(db)), function(i) {
    length(intersect(cluster_markers_upper, db$markers_pos[[i]]))
  }, integer(1))

  data.frame(
    cellName     = db$cellName,
    tissue       = db$tissueType,
    score        = scores,
    n_overlap    = n_overlap,
    n_db_markers = vapply(db$markers_pos, length, integer(1)),
    stringsAsFactors = FALSE
  )
}

# Map a free-text tissue to a ScTypeDB tissueType (improves accuracy by scoring
# only the relevant tissue's cell types). Returns NULL when no confident match.
.match_sctype_tissue <- function(tissue, db_tissues) {
  if (is.null(tissue) || !nzchar(trimws(tissue))) return(NULL)
  t <- tolower(trimws(tissue))
  aliases <- list(
    "Immune system" = c("pbmc","blood","immune","bone marrow","lymph","peripheral","leukocyte","wbc","spleen"),
    "Brain"     = c("brain","cortex","neuron","cerebr","hippocamp","cns","neural"),
    "Lung"      = c("lung","pulmonary","airway","alveol","bronch"),
    "Liver"     = c("liver","hepat"),
    "Kidney"    = c("kidney","renal","nephron"),
    "Heart"     = c("heart","cardiac","myocard"),
    "Pancreas"  = c("pancrea","islet"),
    "Intestine" = c("intestin","gut","colon","bowel","ileum"),
    "Stomach"   = c("stomach","gastric"),
    "Muscle"    = c("muscle","skeletal","myo"),
    "Eye"       = c("eye","retina","ocular"),
    "Placenta"  = c("placenta","trophoblast"),
    "Thymus"    = c("thymus","thymic"),
    "Spleen"    = c("spleen","splenic"),
    "Adrenal"   = c("adrenal")
  )
  for (dt in db_tissues)                     # direct substring either way
    if (grepl(tolower(dt), t, fixed = TRUE) || grepl(t, tolower(dt), fixed = TRUE)) return(dt)
  for (canon in names(aliases))              # alias match
    if (canon %in% db_tissues &&
        any(vapply(aliases[[canon]], function(a) grepl(a, t, fixed = TRUE), logical(1))))
      return(canon)
  NULL
}

run_sctype_annotation <- function(seurat_obj, db_path = NULL, cluster_col = "seurat_clusters",
                                   top_n_markers = 30, markers_df = NULL, tissue = NULL) {
  db <- .load_sctype_db(db_path)

  # Tissue-restricted scoring (ScType's intended workflow) when a tissue is known.
  # Falls back to scoring across all tissues when there is no confident match.
  matched_tissue <- .match_sctype_tissue(tissue, unique(db$tissueType))
  if (!is.null(matched_tissue)) {
    db_f <- db[db$tissueType == matched_tissue, , drop = FALSE]
    if (nrow(db_f) >= 3) {
      db <- db_f
      message("[ScType] Restricted scoring to tissue: ", matched_tissue,
              " (", nrow(db), " cell types)")
    }
  }

  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      cluster_col <- "seurat_clusters"
    } else {
      seurat_obj$seurat_clusters <- Idents(seurat_obj)
      cluster_col <- "seurat_clusters"
    }
  }
  Idents(seurat_obj) <- cluster_col
  clusters <- sort(unique(Idents(seurat_obj)))

  if (is.null(markers_df)) {
    markers_df <- .find_markers_for_annotation(seurat_obj, cluster_col)
  }

  all_results <- list()
  cluster_assignments <- data.frame(
    cluster = character(0), type = character(0), tissue = character(0),
    scores = numeric(0), n_overlap = integer(0), ncells = integer(0),
    stringsAsFactors = FALSE
  )

  for (cl in clusters) {
    cl_char <- as.character(cl)
    cl_markers <- markers_df[markers_df$cluster == cl_char, , drop = FALSE]
    ncells_cl <- sum(Idents(seurat_obj) == cl)
    if (nrow(cl_markers) == 0) {
      cluster_assignments <- rbind(cluster_assignments, data.frame(
        cluster = cl_char, type = "Unknown", tissue = NA_character_,
        scores = 0, n_overlap = 0L, ncells = ncells_cl,
        stringsAsFactors = FALSE))
      next
    }
    cl_markers <- cl_markers[order(cl_markers$avg_log2FC, decreasing = TRUE), ]
    cl_markers <- utils::head(cl_markers, top_n_markers)

    gene_names <- toupper(cl_markers$gene)
    gene_fc    <- setNames(cl_markers$avg_log2FC, gene_names)

    scored <- .score_cluster_vs_sctype(gene_names, gene_fc, db)
    scored <- scored[order(scored$score, decreasing = TRUE), ]
    all_results[[cl_char]] <- utils::head(scored, 10)
    best <- scored[1, ]

    cluster_assignments <- rbind(cluster_assignments, data.frame(
      cluster = cl_char, type = best$cellName, tissue = best$tissue,
      scores = best$score, n_overlap = best$n_overlap, ncells = ncells_cl,
      stringsAsFactors = FALSE))
  }

  cluster_assignments$type[
    cluster_assignments$n_overlap < 2 | cluster_assignments$scores <= 0
  ] <- "Unknown"

  label_map  <- setNames(cluster_assignments$type, cluster_assignments$cluster)
  tissue_map <- setNames(cluster_assignments$tissue, cluster_assignments$cluster)
  score_map  <- setNames(cluster_assignments$scores, cluster_assignments$cluster)
  cl_vec     <- as.character(seurat_obj@meta.data[[cluster_col]])

  seurat_obj@meta.data[["cell_type_sctype"]] <- unname(label_map[cl_vec])
  seurat_obj@meta.data[["sctype_tissue"]]    <- unname(tissue_map[cl_vec])
  seurat_obj@meta.data[["sctype_score"]]     <- unname(as.numeric(score_map[cl_vec]))

  list(obj = seurat_obj, markers = markers_df,
       cluster_results = all_results, cluster_assignments = cluster_assignments)
}

# =============================================================================
# PanglaoDB: Automated Cell Type Annotation
# Uses local PanglaoDB_markers TSV with ubiquitousness, sensitivity, specificity
# =============================================================================

.panglao_db_path <- function() {
  candidates <- c(
    file.path(getwd(), "PanglaoDB", "PanglaoDB_markers_27_Mar_2020.tsv"),
    file.path(getwd(), "PanglaoDB_markers_27_Mar_2020.tsv"),
    file.path(dirname(getwd()), "PanglaoDB", "PanglaoDB_markers_27_Mar_2020.tsv"),
    Sys.glob(file.path(getwd(), "PanglaoDB", "PanglaoDB_markers*.tsv")),
    Sys.glob(file.path(getwd(), "PanglaoDB_markers*.tsv"))
  )
  candidates <- candidates[nchar(candidates) > 0]
  for (p in candidates) {
    if (file.exists(p)) return(normalizePath(p, winslash = "/"))
  }
  NULL
}

.load_panglao_db <- function(db_path = NULL, species = "human") {
  if (is.null(db_path)) db_path <- .panglao_db_path()
  if (is.null(db_path) || !file.exists(db_path)) {
    stop("PanglaoDB markers TSV not found. Place PanglaoDB_markers_27_Mar_2020.tsv in the app directory.")
  }
  db <- read.delim(db_path, stringsAsFactors = FALSE, check.names = FALSE)

  # Standardize column names
  colnames(db) <- gsub("^\\s+|\\s+$", "", colnames(db))
  col_map <- c(
    "species" = "species",
    "official gene symbol" = "gene",
    "cell type" = "cell_type",
    "ubiquitousness index" = "ubiquitousness",
    "canonical marker" = "canonical",
    "organ" = "organ",
    "germ layer" = "germ_layer",
    "sensitivity_human" = "sensitivity_human",
    "sensitivity_mouse" = "sensitivity_mouse",
    "specificity_human" = "specificity_human",
    "specificity_mouse" = "specificity_mouse",
    "gene type" = "gene_type"
  )
  for (orig in names(col_map)) {
    idx <- which(tolower(colnames(db)) == tolower(orig))
    if (length(idx) == 1) colnames(db)[idx] <- col_map[[orig]]
  }

  # Filter by species
  sp_code <- if (tolower(species) == "mouse") "Mm" else "Hs"
  db <- db[grepl(sp_code, db$species, fixed = TRUE), ]

  # Ensure numeric columns
  for (col in c("ubiquitousness", "sensitivity_human", "sensitivity_mouse",
                "specificity_human", "specificity_mouse")) {
    if (col %in% colnames(db)) {
      db[[col]] <- suppressWarnings(as.numeric(db[[col]]))
      db[[col]][is.na(db[[col]])] <- 0
    }
  }
  if (!"canonical" %in% colnames(db)) db$canonical <- 1
  db$canonical <- suppressWarnings(as.numeric(db$canonical))
  db$canonical[is.na(db$canonical)] <- 0

  db$gene <- toupper(db$gene)

  # Aggregate markers per cell type
  cell_types <- unique(db$cell_type)
  ct_list <- lapply(cell_types, function(ct) {
    sub <- db[db$cell_type == ct, , drop = FALSE]
    organs <- unique(na.omit(sub$organ))
    organs <- organs[organs != "NA" & organs != ""]

    sens_col <- if (tolower(species) == "mouse") "sensitivity_mouse" else "sensitivity_human"
    spec_col <- if (tolower(species) == "mouse") "specificity_mouse" else "specificity_human"

    list(
      cell_type = ct,
      organ     = paste(organs, collapse = "; "),
      genes     = sub$gene,
      ubiq      = sub$ubiquitousness,
      canonical = sub$canonical,
      sensitivity = if (sens_col %in% colnames(sub)) sub[[sens_col]] else rep(0, nrow(sub)),
      specificity = if (spec_col %in% colnames(sub)) sub[[spec_col]] else rep(0, nrow(sub))
    )
  })
  names(ct_list) <- cell_types
  ct_list
}

# PanglaoDB scoring: for each cluster's markers, score against all PanglaoDB cell types
# Uses: FC weighting * (1 - ubiquitousness) * (1 + specificity), canonical markers get boost
.score_cluster_vs_panglao <- function(cluster_markers_upper, cluster_fc, ct_list) {
  results <- lapply(ct_list, function(ct) {
    db_genes  <- ct$genes
    ubiq      <- ct$ubiq
    canonical <- ct$canonical
    spec      <- ct$specificity
    n_markers <- length(db_genes)
    if (n_markers == 0) {
      return(data.frame(cell_type = ct$cell_type, organ = ct$organ,
                        score = 0, n_overlap = 0L, n_db_markers = 0L,
                        stringsAsFactors = FALSE))
    }

    hits_idx <- which(db_genes %in% cluster_markers_upper)
    n_hits   <- length(hits_idx)
    if (n_hits == 0) {
      return(data.frame(cell_type = ct$cell_type, organ = ct$organ,
                        score = 0, n_overlap = 0L, n_db_markers = n_markers,
                        stringsAsFactors = FALSE))
    }

    hit_genes <- db_genes[hits_idx]
    hit_fc    <- cluster_fc[hit_genes]
    hit_fc[is.na(hit_fc)] <- 0.25

    # Weights per matching marker
    # Penalize ubiquitous genes (they're less informative)
    ubiq_weight <- 1 - pmin(ubiq[hits_idx], 0.9)
    # Boost specific markers
    spec_weight <- 1 + spec[hits_idx]
    # Canonical markers get extra weight
    canon_weight <- ifelse(canonical[hits_idx] == 1, 1.5, 1.0)

    weighted_score <- sum(hit_fc * ubiq_weight * spec_weight * canon_weight, na.rm = TRUE)

    # Normalize: score / sqrt(n_db_markers) * (1 + coverage)
    coverage <- n_hits / n_markers
    final_score <- (weighted_score / sqrt(n_markers)) * (1 + coverage)

    data.frame(cell_type = ct$cell_type, organ = ct$organ,
               score = final_score, n_overlap = as.integer(n_hits),
               n_db_markers = as.integer(n_markers),
               stringsAsFactors = FALSE)
  })
  do.call(rbind, results)
}

run_panglao_annotation <- function(seurat_obj, db_path = NULL, cluster_col = "seurat_clusters",
                                    top_n_markers = 30, markers_df = NULL, species = NULL) {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
      cluster_col <- "seurat_clusters"
    } else {
      seurat_obj$seurat_clusters <- Idents(seurat_obj)
      cluster_col <- "seurat_clusters"
    }
  }
  Idents(seurat_obj) <- cluster_col
  clusters <- sort(unique(Idents(seurat_obj)))

  # Auto-detect species
  if (is.null(species)) {
    species <- .detect_species(rownames(seurat_obj))
  }

  ct_list <- .load_panglao_db(db_path, species = species)

  if (is.null(markers_df)) {
    markers_df <- .find_markers_for_annotation(seurat_obj, cluster_col)
  }

  all_results <- list()
  cluster_assignments <- data.frame(
    cluster = character(0), type = character(0), organ = character(0),
    scores = numeric(0), n_overlap = integer(0), ncells = integer(0),
    stringsAsFactors = FALSE
  )

  for (cl in clusters) {
    cl_char <- as.character(cl)
    cl_markers <- markers_df[markers_df$cluster == cl_char, , drop = FALSE]
    ncells_cl <- sum(Idents(seurat_obj) == cl)
    if (nrow(cl_markers) == 0) {
      cluster_assignments <- rbind(cluster_assignments, data.frame(
        cluster = cl_char, type = "Unknown", organ = NA_character_,
        scores = 0, n_overlap = 0L, ncells = ncells_cl,
        stringsAsFactors = FALSE))
      next
    }

    cl_markers <- cl_markers[order(cl_markers$avg_log2FC, decreasing = TRUE), ]
    cl_markers <- utils::head(cl_markers, top_n_markers)

    gene_names <- toupper(cl_markers$gene)
    gene_fc    <- setNames(cl_markers$avg_log2FC, gene_names)

    scored <- .score_cluster_vs_panglao(gene_names, gene_fc, ct_list)
    scored <- scored[order(scored$score, decreasing = TRUE), ]
    all_results[[cl_char]] <- utils::head(scored, 10)
    best <- scored[1, ]

    cluster_assignments <- rbind(cluster_assignments, data.frame(
      cluster = cl_char, type = best$cell_type, organ = best$organ,
      scores = best$score, n_overlap = best$n_overlap, ncells = ncells_cl,
      stringsAsFactors = FALSE))
  }

  cluster_assignments$type[
    cluster_assignments$n_overlap < 2 | cluster_assignments$scores <= 0
  ] <- "Unknown"

  label_map <- setNames(cluster_assignments$type, cluster_assignments$cluster)
  organ_map <- setNames(cluster_assignments$organ, cluster_assignments$cluster)
  score_map <- setNames(cluster_assignments$scores, cluster_assignments$cluster)
  cl_vec    <- as.character(seurat_obj@meta.data[[cluster_col]])

  seurat_obj@meta.data[["cell_type_panglao"]] <- unname(label_map[cl_vec])
  seurat_obj@meta.data[["panglao_organ"]]     <- unname(organ_map[cl_vec])
  seurat_obj@meta.data[["panglao_score"]]     <- unname(as.numeric(score_map[cl_vec]))

  list(obj = seurat_obj, markers = markers_df,
       cluster_results = all_results, cluster_assignments = cluster_assignments,
       species = species)
}

# =============================================================================
# clustifyr: Reference-correlation Cell Type Annotation (celldex references)
# =============================================================================
# Methodologically independent of ScType/PanglaoDB (marker-overlap on
# FindAllMarkers). clustifyr correlates each CLUSTER's average expression
# (pseudobulk) against curated reference cell-type profiles (Spearman on shared
# variable genes) and assigns each cluster its best-correlated type. Operating at
# cluster level makes it robust: every cluster is scored independently, so it will
# not collapse all clusters onto one label the way per-cell classifiers can.

# Human/mouse reference loaders (celldex; fetched via ExperimentHub, cached locally)
.REFCORR_REFS <- list(
  "HumanPrimaryCellAtlas"        = list(species = "human", fn = function() celldex::HumanPrimaryCellAtlasData()),
  "BlueprintEncode"              = list(species = "human", fn = function() celldex::BlueprintEncodeData()),
  "MonacoImmune"                 = list(species = "human", fn = function() celldex::MonacoImmuneData()),
  "DatabaseImmuneCellExpression" = list(species = "human", fn = function() celldex::DatabaseImmuneCellExpressionData()),
  "NovershternHematopoietic"     = list(species = "human", fn = function() celldex::NovershternHematopoieticData()),
  "MouseRNAseq"                  = list(species = "mouse", fn = function() celldex::MouseRNAseqData()),
  "ImmGen"                       = list(species = "mouse", fn = function() celldex::ImmGenData()),
  # Han-lab whole-organism atlases (genes x cell-type matrices; very granular)
  "HumanCellLandscape (HCL)"     = list(species = "human", type = "matrix",
    fn = function() { e <- new.env(); utils::data("ref.expr", package = "scHCL", envir = e); as.matrix(get("ref.expr", e)) }),
  "MouseCellAtlas (MCA)"         = list(species = "mouse", type = "matrix",
    fn = function() { e <- new.env(); utils::data("ref.expr", package = "scMCA", envir = e); as.matrix(get("ref.expr", e)) })
)

#' References available for a given species (named vector display -> key)
clustifyr_refs_for_species <- function(species) {
  sp <- if (exists("normalize_species")) normalize_species(species) else species
  if (!sp %in% c("human", "mouse")) sp <- "human"
  keys <- names(.REFCORR_REFS)[vapply(.REFCORR_REFS, function(x) x$species == sp, logical(1))]
  stats::setNames(keys, keys)
}

.default_clustifyr_ref <- function(species) {
  sp <- if (exists("normalize_species")) normalize_species(species) else species
  if (identical(sp, "mouse")) "MouseRNAseq" else "HumanPrimaryCellAtlas"
}

#' Build a genes x cell-type average matrix from a celldex reference.
.celldex_ref_matrix <- function(ref, label_field = "label.main") {
  labs <- SummarizedExperiment::colData(ref)[[label_field]]
  if (is.null(labs)) labs <- SummarizedExperiment::colData(ref)[["label.main"]]
  labs <- as.character(labs)
  mat  <- as.matrix(SummarizedExperiment::assay(ref))
  uniq <- unique(labs)
  refm <- vapply(uniq, function(l) rowMeans(mat[, labs == l, drop = FALSE]),
                 numeric(nrow(mat)))
  colnames(refm) <- uniq
  rownames(refm) <- rownames(mat)
  refm
}

#' Run clustifyr reference-correlation annotation.
#' Stores per-cell `cell_type_clustifyr` + `clustifyr_score` (cluster-uniform).
#' @return list(obj, ref, n_labels, cluster_labels)
run_clustifyr_annotation <- function(seurat_obj, species = "human", ref_name = NULL,
                                     cluster_col = "seurat_clusters", assay = NULL,
                                     label_field = "label.main", query_genes = NULL) {
  for (p in c("clustifyr", "celldex", "SummarizedExperiment"))
    if (!requireNamespace(p, quietly = TRUE))
      stop("Package '", p, "' is required for clustifyr annotation. Install it from Bioconductor.")

  if (is.null(assay) || !nzchar(assay)) assay <- SeuratObject::DefaultAssay(seurat_obj)

  # Ensure a normalized "data" layer is available
  a <- seurat_obj[[assay]]
  if (inherits(a, "Assay5")) {
    lyrs <- SeuratObject::Layers(a)
    if (any(grepl("^data", lyrs)) && !("data" %in% lyrs))
      seurat_obj[[assay]] <- SeuratObject::JoinLayers(seurat_obj[[assay]])
  }
  logdata <- tryCatch(
    SeuratObject::GetAssayData(seurat_obj, assay = assay, layer = "data"),
    error = function(e) tryCatch(
      SeuratObject::GetAssayData(seurat_obj, assay = assay, slot = "data"),
      error = function(e2) NULL)
  )
  if (is.null(logdata) || nrow(logdata) == 0 ||
      isTRUE(all.equal(max(logdata), min(logdata))))
    stop("Normalized data not found. Run Normalization (Preprocessing) before clustifyr.")

  if (!cluster_col %in% colnames(seurat_obj@meta.data))
    stop("Cluster column '", cluster_col, "' not found. Run Clustering first.")
  clusters <- as.character(seurat_obj@meta.data[[cluster_col]])
  names(clusters) <- rownames(seurat_obj@meta.data)

  # Resolve + load reference; build genes x cell-type matrix.
  # celldex refs are SummarizedExperiments (averaged by label); Han-lab atlases
  # (HCL/MCA) are already genes x cell-type matrices.
  if (is.null(ref_name) || !nzchar(ref_name) || is.null(.REFCORR_REFS[[ref_name]]))
    ref_name <- .default_clustifyr_ref(species)
  entry   <- .REFCORR_REFS[[ref_name]]
  raw     <- entry$fn()
  ref_mat <- if (identical(entry$type, "matrix")) as.matrix(raw)
             else .celldex_ref_matrix(raw, label_field)

  # Query genes: focus the correlation on DISCRIMINATIVE genes for higher
  # per-cluster accuracy — cluster markers (from FindAllMarkers) + variable
  # features, intersected with genes shared with the reference. Fall back to all
  # shared genes when too few discriminative genes overlap the reference.
  shared <- intersect(rownames(logdata), rownames(ref_mat))
  candidate <- character(0)
  if (!is.null(query_genes)) candidate <- union(candidate, as.character(query_genes))
  vf <- tryCatch(SeuratObject::VariableFeatures(seurat_obj), error = function(e) NULL)
  if (!is.null(vf) && length(vf) > 0) candidate <- union(candidate, as.character(vf))
  qgenes <- intersect(candidate, shared)
  if (length(qgenes) < 50) qgenes <- shared    # not enough discriminative overlap
  if (length(qgenes) < 20)
    stop("Only ", length(qgenes), " shared genes between data and reference. ",
         "Check species/gene-symbol format (symbols vs Ensembl IDs).")

  cells <- intersect(colnames(logdata), names(clusters))
  meta  <- data.frame(cluster = clusters[cells], row.names = cells)

  cor_mat <- clustifyr::clustify(
    input       = as.matrix(logdata[, cells, drop = FALSE]),
    metadata    = meta,
    ref_mat     = ref_mat,
    cluster_col = "cluster",
    query_genes = qgenes
  )
  calls <- clustifyr::cor_to_call(cor_mat)          # cluster, type, r
  cl_ids <- as.character(calls[[1]])
  lab_map   <- stats::setNames(as.character(calls$type), cl_ids)
  score_map <- stats::setNames(as.numeric(calls$r),      cl_ids)
  # Treat low-confidence / unassigned calls as Unknown
  lab_map[is.na(lab_map) | grepl("^r<|unassigned", lab_map, ignore.case = TRUE)] <- "Unknown"
  # Clean the reference's verbose cell-type names to just the cell type.
  # (.clean_atlas_label drops names, so re-attach the cluster ids afterwards.)
  atlas <- if (grepl("MCA|MouseCellAtlas", ref_name)) "mca"
           else if (grepl("HCL|HumanCellLandscape", ref_name)) "hcl"
           else "celldex"
  keep_unknown <- lab_map == "Unknown"
  lab_map <- stats::setNames(.clean_atlas_label(lab_map, atlas), names(lab_map))
  lab_map[keep_unknown] <- "Unknown"

  all_cells <- rownames(seurat_obj@meta.data)
  cl_all    <- as.character(seurat_obj@meta.data[[cluster_col]])
  seurat_obj@meta.data[["cell_type_clustifyr"]] <- unname(lab_map[cl_all])
  seurat_obj@meta.data[["clustifyr_score"]]     <- unname(as.numeric(score_map[cl_all]))

  list(obj = seurat_obj, ref = ref_name,
       n_labels = length(setdiff(unique(lab_map), "Unknown")),
       cluster_labels = lab_map)
}

# =============================================================================
# General N-method consensus (majority vote; ties broken by score)
# =============================================================================
# Inputs are per-cell character vectors, cluster-uniform for all three methods
# (ScType, PanglaoDB, clustifyr all label at cluster level). Returns consensus
# label, source, score, and a disagreement flag.

.build_annotation_consensus <- function(label_list, score_list) {
  meth <- names(label_list)
  n <- length(label_list[[1]])
  labs <- do.call(cbind, lapply(label_list, function(x) { x[is.na(x)] <- "Unknown"; as.character(x) }))
  scos <- do.call(cbind, lapply(score_list, function(x) suppressWarnings(as.numeric(x))))
  colnames(labs) <- meth; colnames(scos) <- meth

  consensus <- character(n); source <- character(n)
  score <- numeric(n); disagree <- logical(n)

  # Operate on unique label/score rows for speed (cluster-uniform data collapses well)
  key <- apply(cbind(labs, round(scos, 4)), 1, paste, collapse = "|")
  for (k in unique(key)) {
    idx <- which(key == k)
    i <- idx[1]
    l <- labs[i, ]; s <- scos[i, ]
    known <- l != "Unknown" & !is.na(l)
    if (!any(known)) {
      consensus[idx] <- "Unknown"; source[idx] <- "none"
      score[idx] <- NA_real_; disagree[idx] <- NA
      next
    }
    lk <- l[known]; sk <- s[known]; mk <- meth[known]
    tab <- sort(table(lk), decreasing = TRUE)
    top <- names(tab)[1]
    if (length(tab) > 1 && tab[1] == tab[2]) top <- lk[which.max(sk)]  # tie -> top score
    agree <- mk[lk == top]
    consensus[idx] <- top
    # "all" only when every method (not just the voting ones) agrees
    source[idx]    <- if (length(agree) == length(meth)) "all" else paste(agree, collapse = "+")
    score[idx]     <- suppressWarnings(max(sk[lk == top], na.rm = TRUE))
    disagree[idx]  <- length(unique(lk)) > 1
  }
  score[!is.finite(score)] <- NA_real_
  list(consensus = consensus, source = source, score = score, disagreement = disagree)
}

# =============================================================================
# clustree integration helper (BUG FIX #5)
# =============================================================================

run_clustree_analysis <- function(seurat_obj, resolutions = seq(0.1, 2.0, by = 0.1),
                                   n_dims = 30, prefix = NULL) {
  if (!.clustree_available) {
    stop("clustree package not available. Install: install.packages('clustree')")
  }
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("ggraph package required for clustree. Install: install.packages('ggraph')")
  }

  # Ensure neighbors exist
  if (length(grep("_snn_", names(seurat_obj@graphs))) == 0) {
    n_pcs <- min(n_dims, ncol(Embeddings(seurat_obj, "pca")))
    k_param <- min(max(20L, floor(sqrt(ncol(seurat_obj)) / 2)), ncol(seurat_obj) - 1L)
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca",
                                 dims = 1:n_pcs, k.param = k_param, verbose = FALSE)
  }

  # Run clustering at each resolution
  for (res in resolutions) {
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
  }

  # Auto-detect prefix from the Seurat object's meta.data (e.g. "RNA_snn_res.",
  # "SCT_snn_res.", "integrated_snn_res.") instead of hard-coding "RNA_snn_res.".
  md <- seurat_obj@meta.data
  all_res <- grep("_snn_res\\.[0-9]", colnames(md), value = TRUE)
  if (is.null(prefix) || !any(startsWith(all_res, prefix))) {
    get_pref <- function(x) sub("[0-9]+(\\.[0-9]+)?$", "", x)
    prefs    <- vapply(all_res, get_pref, character(1))
    tbl      <- sort(table(prefs), decreasing = TRUE)
    prefix   <- if (length(tbl)) names(tbl)[1] else "RNA_snn_res."
  }
  res_cols <- all_res[startsWith(all_res, prefix)]

  # Drop resolution columns with <2 unique non-NA values
  res_cols <- res_cols[vapply(res_cols, function(cc) {
    u <- unique(md[[cc]][!is.na(md[[cc]])])
    length(u) >= 2
  }, logical(1))]
  if (length(res_cols) < 2) {
    stop("Clustree needs >=2 resolution columns with >=2 clusters each.")
  }

  # Pass a plain factor data.frame (not the Seurat object) to clustree to
  # avoid S4 / S7 dispatch issues with newer ggplot2 / ggraph versions.
  md_cl <- md[, res_cols, drop = FALSE]
  for (cc in res_cols) md_cl[[cc]] <- as.factor(md_cl[[cc]])

  # Build clustree and then strip the ggraph-custom edge guides. ggplot2 >= 4.0
  # introduced an S7 guide dispatch that cannot find ggraph's "edge_colourbar"
  # unless ggraph is attached at *print* time; suppressing these guides avoids
  # the "Unknown guide: edge_colourbar" render crash.
  p <- clustree::clustree(md_cl, prefix = prefix) +
    ggplot2::guides(edge_colour = "none",
                    edge_alpha  = "none",
                    edge_width  = "none") +
    ggplot2::theme(legend.position = "right") +
    ggplot2::ggtitle("Clustree: Resolution Exploration")

  # Materialize the plot (force evaluation) so downstream renderPlot
  # does not re-evaluate any lazy promises from the calling scope.
  force(p)

  list(obj = seurat_obj, plot = p, resolutions = resolutions, prefix = prefix)
}

# =============================================================================
# UI: Automatic Annotation (ScType)
# =============================================================================

automaticAnnotationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          if (exists("build_method_card_ui", mode = "function"))
            build_method_card_ui("annotation"),
          h4("Automatic Cell Type Annotation"),
          helpText("Finds top markers per cluster (FindAllMarkers), then matches them ",
                   "against marker databases. No tissue selection needed."),

          fluidRow(
            column(3,
              selectInput(ns("annotation_method"), "Method",
                choices = c("ScType" = "sctype",
                            "PanglaoDB" = "panglao",
                            "clustifyr (reference-based)" = "clustifyr",
                            "ScType + PanglaoDB (consensus)" = "both",
                            "All three (consensus)" = "all"),
                selected = "all"
              ),
              conditionalPanel(
                condition = "input.annotation_method == 'clustifyr' || input.annotation_method == 'all'",
                ns = ns,
                uiOutput(ns("clustifyr_ref_ui"))
              ),
              uiOutput(ns("group_column_ui"))
            ),
            column(3,
              numericInput(ns("top_n_markers"), "Top markers per cluster",
                value = 30, min = 5, max = 100, step = 5),
              helpText(textOutput(ns("db_status_text")), style = "font-size: 11px;")
            ),
            column(3,
              actionButton(ns("run_annotation"), "Run Annotation",
                class = "btn-primary btn-lg", width = "100%"),
              br(), br(),
              helpText(textOutput(ns("species_text")), style = "font-size: 11px;")
            ),
            column(3,
              checkboxInput(ns("show_clustree"),
                            "Show clustree resolution explorer", FALSE),
              conditionalPanel(
                condition = "input.show_clustree",
                ns = ns,
                sliderInput(ns("clustree_res_range"), "Resolution Range",
                  min = 0.1, max = 3.0, value = c(0.2, 1.5), step = 0.1
                ),
                actionButton(ns("run_clustree"), "Generate clustree",
                  class = "btn-info", width = "100%")
              )
            )
          )
        )
      )
    ),

    # Clustree output
    fluidRow(
      column(12, uiOutput(ns("clustree_output")))
    ),

    # Annotation results
    fluidRow(
      column(12, uiOutput(ns("annotation_results")))
    )
  )
}

# =============================================================================
# Server: Automatic Annotation
# =============================================================================

automaticAnnotationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    sctype_result  <- reactiveVal(NULL)
    panglao_result <- reactiveVal(NULL)
    clustifyr_result <- reactiveVal(NULL)

    # DB status text
    output$db_status_text <- renderText({
      sc_ok <- tryCatch({
        p <- .sctype_db_path()
        !is.null(p) && file.exists(p)
      }, error = function(e) FALSE)
      pg_ok <- tryCatch({
        p <- .panglao_db_path()
        !is.null(p) && file.exists(p)
      }, error = function(e) FALSE)
      paste0("ScType: ", if (sc_ok) "OK" else "missing",
             " | PanglaoDB: ", if (pg_ok) "OK" else "missing")
    })

    # Species auto-detect
    output$species_text <- renderText({
      req(values$active_dataset)
      sp <- .detect_species(rownames(values$active_dataset))
      paste0("Detected species: ", sp)
    })

    # clustifyr reference selector (species-aware)
    output$clustifyr_ref_ui <- renderUI({
      req(values$active_dataset)
      sp    <- .detect_species(rownames(values$active_dataset))
      refs  <- clustifyr_refs_for_species(sp)
      selectInput(session$ns("clustifyr_ref"),
                  "clustifyr reference",
                  choices  = refs,
                  selected = .default_clustifyr_ref(sp))
    })

    # Group column selector
    output$group_column_ui <- renderUI({
      req(values$active_dataset)
      md   <- values$active_dataset@meta.data
      cand <- colnames(md)[sapply(md, function(x) is.factor(x) || is.character(x))]
      cand <- union(cand, grep("cluster|leiden|louvain|slm|res\\.",
                                colnames(md), value = TRUE, ignore.case = TRUE))
      selectInput(session$ns("groupby"), "Cluster column",
                  choices = sort(unique(cand)), selected = "seurat_clusters")
    })

    # ---------------------------------------------------------------
    # clustree handler
    # ---------------------------------------------------------------
    # Store the clustree plot in a reactiveVal so renderPlot does not
    # depend on an outer observer scope (which was causing
    # "object '<x>' not found" evaluation errors).
    clustree_plot_rv <- reactiveVal(NULL)

    observeEvent(input$run_clustree, {
      req(values$active_dataset)
      if (!.clustree_available) {
        showNotification("clustree package not installed. Run install.packages('clustree')",
                         type = "error", duration = 10)
        return(NULL)
      }

      withProgress(message = "Running clustree analysis...", value = 0, {
        tryCatch({
          obj <- values$active_dataset
          obj <- .safe_prepare_for_de(obj)

          if (!"pca" %in% names(obj@reductions)) {
            incProgress(0.2, detail = "Running PCA...")
            n_hvgs <- if (length(VariableFeatures(obj)) > 0) length(VariableFeatures(obj)) else 2000L
            obj <- FindVariableFeatures(obj, nfeatures = n_hvgs, verbose = FALSE)
            if (DefaultAssay(obj) != "SCT") obj <- ScaleData(obj, verbose = FALSE)
            obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
          }

          res_range <- input$clustree_res_range
          resolutions <- seq(res_range[1], res_range[2], by = 0.1)

          incProgress(0.4, detail = "Clustering at multiple resolutions...")
          result <- run_clustree_analysis(obj, resolutions = resolutions)

          values$active_dataset <- result$obj
          clustree_plot_rv(result$plot)

          output$clustree_output <- renderUI({
            wellPanel(
              h5("Clustree: Resolution Exploration"),
              plotOutput(session$ns("clustree_plot"), height = "600px")
            )
          })

          incProgress(1.0)
          showNotification("Clustree analysis complete.", type = "message")
        }, error = function(e) {
          showNotification(paste("Clustree failed:", conditionMessage(e)),
                           type = "error", duration = 10)
        })
      })
    })

    output$clustree_plot <- renderPlot({
      p <- clustree_plot_rv()
      if (is.null(p)) {
        plot.new(); title("Run clustree to see the resolution tree.")
        return(invisible(NULL))
      }
      tryCatch(print(p), error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Clustree render error:", conditionMessage(e)), cex = 1)
      })
    })

    # ---------------------------------------------------------------
    # Core annotation runner
    # ---------------------------------------------------------------
    observeEvent(input$run_annotation, {
      req(values$active_dataset)
      obj <- values$active_dataset

      if (ncol(obj) < 100) {
        showNotification(
          "Warning: Dataset has fewer than 100 cells. Results may be unreliable.",
          type = "warning", duration = 10
        )
      }

      method      <- input$annotation_method %||% "both"
      cluster_col <- input$groupby %||% "seurat_clusters"
      top_n       <- input$top_n_markers %||% 30L

      job_id <- begin_annotation_job(values, "run_annotation", meta = list(
        method = method, cluster_col = cluster_col, top_n = top_n
      ))

      if (!validate_annotation_contract(values, values$active_dataset)) {
        end_annotation_job(values, job_id, status = "failed",
                           message = "contract validation failed")
        return()
      }

      tryCatch({
        obj <- values$active_dataset
        obj <- .safe_prepare_for_de(obj)

        withProgress(message = "Running annotation...", value = 0, {

          # Step 1: Find markers once (shared by both methods)
          incProgress(0.1, detail = "Finding cluster markers (FindAllMarkers)...")
          markers_df <- .find_markers_for_annotation(obj, cluster_col)

          run_sctype  <- method %in% c("sctype", "both", "all")
          run_panglao <- method %in% c("panglao", "both", "all")
          run_clustifyr <- method %in% c("clustifyr", "all")

          sc_res <- NULL
          pg_res <- NULL
          cf_res <- NULL

          # Step 2: ScType
          if (run_sctype) {
            sc_db <- .sctype_db_path()
            if (is.null(sc_db) || !file.exists(sc_db)) {
              showNotification("ScTypeDB_full.xlsx not found, skipping ScType.",
                               type = "warning", duration = 8)
              run_sctype <- FALSE
            } else {
              incProgress(0.3, detail = "Scoring against ScTypeDB...")
              sc_res <- run_sctype_annotation(
                obj, db_path = sc_db, cluster_col = cluster_col,
                top_n_markers = top_n, markers_df = markers_df,
                tissue = values$tissue_global %||% input$tissue_context %||% NULL
              )
              obj <- sc_res$obj
              sctype_result(sc_res)
            }
          }

          # Step 3: PanglaoDB
          if (run_panglao) {
            pg_db <- .panglao_db_path()
            if (is.null(pg_db) || !file.exists(pg_db)) {
              showNotification("PanglaoDB markers TSV not found, skipping PanglaoDB.",
                               type = "warning", duration = 8)
              run_panglao <- FALSE
            } else {
              incProgress(0.5, detail = "Scoring against PanglaoDB...")
              pg_res <- run_panglao_annotation(
                obj, db_path = pg_db, cluster_col = cluster_col,
                top_n_markers = top_n, markers_df = markers_df
              )
              obj <- pg_res$obj
              panglao_result(pg_res)
            }
          }

          # Step 3.5: clustifyr (reference-based, methodologically independent)
          if (run_clustifyr) {
            incProgress(0.6, detail = "Running clustifyr (reference-based)...")
            cf_res <- tryCatch(
              run_clustifyr_annotation(
                obj,
                species     = .detect_species(rownames(obj)),
                ref_name    = input$clustifyr_ref %||% NULL,
                cluster_col = cluster_col,
                query_genes = if (!is.null(markers_df) && "gene" %in% colnames(markers_df))
                                unique(as.character(markers_df$gene)) else NULL
              ),
              error = function(e) {
                showNotification(paste("clustifyr skipped:", conditionMessage(e)),
                                 type = "warning", duration = 10)
                NULL
              }
            )
            if (!is.null(cf_res)) {
              obj <- cf_res$obj
              clustifyr_result(cf_res)
            } else {
              run_clustifyr <- FALSE
            }
          }

          incProgress(0.75, detail = "Building consensus...")

          # Step 4: Consensus across whichever methods ran (majority vote, ties by
          # score). All three methods label at cluster level, so votes are directly
          # comparable.
          label_list <- list(); score_list <- list()
          if (run_sctype) {
            label_list$ScType <- as.character(obj@meta.data[["cell_type_sctype"]])
            score_list$ScType <- as.numeric(obj@meta.data[["sctype_score"]])
          }
          if (run_panglao) {
            label_list$PanglaoDB <- as.character(obj@meta.data[["cell_type_panglao"]])
            score_list$PanglaoDB <- as.numeric(obj@meta.data[["panglao_score"]])
          }
          if (run_clustifyr) {
            # clustifyr labels are already cluster-level (uniform within a cluster)
            label_list$clustifyr <- as.character(obj@meta.data[["cell_type_clustifyr"]])
            score_list$clustifyr <- as.numeric(obj@meta.data[["clustifyr_score"]])
          }

          if (length(label_list) >= 2) {
            cons <- .build_annotation_consensus(label_list, score_list)
            obj@meta.data[["cell_type_consensus"]]         <- cons$consensus
            obj@meta.data[["annotation_consensus_source"]] <- cons$source
            obj@meta.data[["annotation_score"]]            <- cons$score
            obj@meta.data[["annotation_disagreement"]]     <- cons$disagreement
          } else if (length(label_list) == 1) {
            only <- names(label_list)[1]
            obj@meta.data[["cell_type_consensus"]]         <- label_list[[1]]
            obj@meta.data[["annotation_consensus_source"]] <- tolower(only)
            obj@meta.data[["annotation_score"]]            <- score_list[[1]]
            obj@meta.data[["annotation_disagreement"]]     <- NA
          }
          obj@meta.data[["annotation_confidence"]] <- NA_real_
          obj@meta.data[["score_difference"]]       <- NA_real_

          # Keep all cluster labels clean ("just the cell type") — strip any
          # score/p-value/parenthetical text a scorer or the AI may have appended.
          for (lbl_col in c("cell_type_consensus", "cell_type_sctype",
                            "cell_type_panglao", "cell_type_clustifyr",
                            "cell_type_enhanced")) {
            if (lbl_col %in% colnames(obj@meta.data))
              obj@meta.data[[lbl_col]] <- .clean_celltype_label(obj@meta.data[[lbl_col]])
          }

          n_types <- length(unique(
            obj@meta.data[["cell_type_consensus"]][
              obj@meta.data[["cell_type_consensus"]] != "Unknown"]))

          methods_used <- c()
          if (run_sctype)  methods_used <- c(methods_used, "ScType")
          if (run_panglao) methods_used <- c(methods_used, "PanglaoDB")
          if (run_clustifyr) methods_used <- c(methods_used, paste0("clustifyr(", cf_res$ref, ")"))

          showNotification(
            paste0(paste(methods_used, collapse = " + "), ": ",
                   n_types, " cell types identified"),
            type = "message", duration = 8
          )

          values$active_dataset <- obj
          checkpoint_annotation(
            values,
            action = "run_annotation",
            obj    = obj,
            params = list(
              method      = method,
              cluster_col = cluster_col,
              top_n       = top_n,
              groupby     = input$groupby
            ),
            method = method
          )
          end_annotation_job(values, job_id, status = "done",
                             message = "Annotation completed")
          values$annotation_done <- TRUE
          incProgress(1)
        })

        showNotification("Annotation complete.", type = "message")
        gc()
      }, error = function(e) {
        end_annotation_job(values, job_id, status = "error", message = e$message)
        showNotification(paste("Annotation failed:", e$message),
                         type = "error", duration = 15)
      })
    })

    # ---------------------------------------------------------------
    # Results UI
    # ---------------------------------------------------------------
    output$annotation_results <- renderUI({
      req(values$active_dataset)
      md_cols  <- colnames(values$active_dataset@meta.data)
      ann_cols <- c("cell_type_consensus", "cell_type_sctype",
                    "cell_type_panglao", "cell_type_clustifyr")
      req(any(ann_cols %in% md_cols))

      ns <- session$ns
      has_3d_umap <- !is.null(values$active_dataset) &&
        "umap3d" %in% names(values$active_dataset@reductions)

      # Per-method comparison plots for whichever methods produced labels
      method_plots <- list()
      if ("cell_type_sctype"  %in% md_cols) method_plots[["ann_dimplot_sctype"]]  <- "ScType"
      if ("cell_type_panglao" %in% md_cols) method_plots[["ann_dimplot_panglao"]] <- "PanglaoDB"
      if ("cell_type_clustifyr" %in% md_cols) method_plots[["ann_dimplot_clustifyr"]] <- "clustifyr"
      show_comparison <- length(method_plots) >= 2
      col_w <- if (length(method_plots) >= 3) 4 else 6

      wellPanel(
        h4("Annotation Results"),
        div(style = "margin-bottom: 25px;",
          fluidRow(
            column(6, plotlyOutput(ns("ann_dimplot"), height = "520px")),
            column(6, DT::DTOutput(ns("ann_table")))
          )
        ),
        if (show_comparison) {
          div(style = "margin-bottom: 25px;",
            h5("Method Comparison"),
            fluidRow(
              lapply(names(method_plots), function(oid)
                column(col_w, plotlyOutput(ns(oid), height = "420px")))
            )
          )
        },
        if (has_3d_umap) {
          div(style = "margin-bottom: 25px;",
            h5("3D UMAP - Annotated Clusters"),
            plotlyOutput(ns("ann_3d_umap"), height = "600px")
          )
        } else {
          div(style = "margin-bottom: 15px;",
            actionButton(ns("run_3d_umap"), "Generate 3D UMAP for Annotation View",
                         class = "btn-info btn-sm"),
            helpText("Runs UMAP with 3 components for 3D visualization.",
                     style = "font-size: 11px;")
          )
        },
        downloadButton(ns("dl_annotations"), "Download annotations (CSV)")
      )
    })

    # ---------------------------------------------------------------
    # 2D DimPlot
    # ---------------------------------------------------------------
    # Helper: make a 2D annotation plotly
    .make_ann_dimplot <- function(obj, labcol, title_suffix = "") {
      red <- .best_reduction(obj)
      if (is.null(red)) return(plotly_empty() %>% layout(title = "No reduction"))
      if (!(labcol %in% colnames(obj@meta.data)) || all(is.na(obj@meta.data[[labcol]])))
        return(plotly_empty() %>% layout(title = paste("No", title_suffix, "annotations")))

      emb  <- Embeddings(obj, red)[, 1:2, drop = FALSE]
      labs <- obj@meta.data[[labcol]]
      n_total <- nrow(emb)
      idx <- if (n_total > 80000 && exists("adaptive_plot_index", mode = "function")) {
        adaptive_plot_index(n_total, max_points = 80000)
      } else seq_len(n_total)

      plot_ly(
        x = emb[idx, 1], y = emb[idx, 2],
        type = "scatter", mode = "markers",
        marker = list(size = 3, opacity = 0.9),
        color = labs[idx], colors = .get_palette(length(unique(labs))),
        text = colnames(obj)[idx],
        hovertemplate = paste0(
          "Cell: %{text}<br>", toupper(red), "1: %{x:.2f}<br>",
          toupper(red), "2: %{y:.2f}<br>",
          "Label: ", labs[idx], "<extra></extra>")
      ) %>% layout(
        title = paste(toupper(red), "-", if (nchar(title_suffix) > 0) title_suffix else labcol),
        xaxis = list(title = paste0(toupper(red), "1")),
        yaxis = list(title = paste0(toupper(red), "2")))
    }

    output$ann_dimplot <- renderPlotly({
      req(values$active_dataset)
      obj <- values$active_dataset
      labcol <- NULL
      for (cand in c("cell_type_consensus", "cell_type_sctype", "cell_type_panglao", "cell_type_clustifyr")) {
        if (cand %in% colnames(obj@meta.data) && !all(is.na(obj@meta.data[[cand]]))) {
          labcol <- cand; break
        }
      }
      if (is.null(labcol)) return(plotly_empty() %>% layout(title = "No annotations available"))
      p <- .make_ann_dimplot(obj, labcol, "Consensus")
      if (exists("capture_umap_snapshot", mode = "function"))
        capture_umap_snapshot(values, p, "Annotation", coloring = labcol,
                              method = .best_reduction(obj))
      if (exists("capture_figure", mode = "function"))
        capture_figure(values, p, "Annotation",
                       figure_type = "DimPlot colored by annotation",
                       label = "Annotation DimPlot", source_module = "Annotation")
      p
    })

    # ---------------------------------------------------------------
    # Comparison DimPlots (when both methods ran)
    # ---------------------------------------------------------------
    output$ann_dimplot_sctype <- renderPlotly({
      req(values$active_dataset)
      .make_ann_dimplot(values$active_dataset, "cell_type_sctype", "ScType")
    })

    output$ann_dimplot_panglao <- renderPlotly({
      req(values$active_dataset)
      .make_ann_dimplot(values$active_dataset, "cell_type_panglao", "PanglaoDB")
    })

    output$ann_dimplot_clustifyr <- renderPlotly({
      req(values$active_dataset)
      req("cell_type_clustifyr" %in% colnames(values$active_dataset@meta.data))
      .make_ann_dimplot(values$active_dataset, "cell_type_clustifyr", "clustifyr")
    })

    # ---------------------------------------------------------------
    # 3D UMAP
    # ---------------------------------------------------------------
    observeEvent(input$run_3d_umap, {
      req(values$active_dataset)
      obj <- values$active_dataset
      withProgress(message = "Computing 3D UMAP...", {
        tryCatch({
          if (!"pca" %in% names(obj@reductions)) {
            incProgress(0.2, detail = "Running PCA first...")
            n_hvgs <- if (length(VariableFeatures(obj)) > 0) length(VariableFeatures(obj)) else 2000L
            obj <- FindVariableFeatures(obj, nfeatures = n_hvgs, verbose = FALSE)
            if (DefaultAssay(obj) != "SCT") obj <- ScaleData(obj, verbose = FALSE)
            obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
          }
          n_dims <- min(30, ncol(Embeddings(obj, "pca")))
          incProgress(0.5, detail = "Running 3D UMAP...")
          obj <- RunUMAP(obj, dims = 1:n_dims, n.components = 3L,
                         umap.method = "uwot", reduction.name = "umap3d", verbose = FALSE)
          values$active_dataset <- obj
          showNotification("3D UMAP computed.", type = "message")
        }, error = function(e) {
          showNotification(paste("3D UMAP failed:", e$message), type = "error")
        })
      })
    })

    output$ann_3d_umap <- renderPlotly({
      req(values$active_dataset)
      obj <- values$active_dataset
      req("umap3d" %in% names(obj@reductions))
      emb <- Embeddings(obj, "umap3d")
      req(ncol(emb) >= 3)

      labcol <- NULL
      for (cand in c("cell_type_consensus", "cell_type_sctype", "cell_type_panglao", "cell_type_clustifyr")) {
        if (cand %in% colnames(obj@meta.data) && !all(is.na(obj@meta.data[[cand]]))) {
          labcol <- cand
          break
        }
      }
      if (is.null(labcol)) return(plotly_empty())
      labs <- obj@meta.data[[labcol]]

      n_cells <- nrow(emb)
      idx <- if (n_cells > 50000 && exists("adaptive_plot_index", mode = "function")) {
        adaptive_plot_index(n_cells, max_points = 50000)
      } else {
        seq_len(n_cells)
      }

      plot_ly(
        x = emb[idx, 1], y = emb[idx, 2], z = emb[idx, 3],
        type = "scatter3d", mode = "markers",
        marker = list(size = 2, opacity = 0.8),
        color  = labs[idx],
        colors = .get_palette(length(unique(labs))),
        text   = colnames(obj)[idx],
        hovertemplate = paste0(
          "Cell: %{text}<br>UMAP1: %{x:.2f}<br>UMAP2: %{y:.2f}<br>UMAP3: %{z:.2f}<br>",
          "Label: ", labs[idx], "<extra></extra>"
        )
      ) %>%
        layout(
          title = "3D UMAP - Annotated Clusters",
          scene = list(
            xaxis = list(title = "UMAP1"),
            yaxis = list(title = "UMAP2"),
            zaxis = list(title = "UMAP3")
          )
        )
    })

    # ---------------------------------------------------------------
    # Results table
    # ---------------------------------------------------------------
    output$ann_table <- DT::renderDT({
      req(values$active_dataset)

      sc_res <- sctype_result()
      pg_res <- panglao_result()

      # Build a combined table per cluster
      if (!is.null(sc_res) || !is.null(pg_res)) {
        sc_tbl <- if (!is.null(sc_res) && !is.null(sc_res$cluster_assignments)) {
          d <- sc_res$cluster_assignments[, c("cluster", "type", "tissue", "scores", "n_overlap"), drop = FALSE]
          colnames(d) <- c("Cluster", "ScType", "ScType_tissue", "ScType_score", "ScType_overlap")
          d$ScType_score <- round(d$ScType_score, 3)
          d
        } else NULL

        pg_tbl <- if (!is.null(pg_res) && !is.null(pg_res$cluster_assignments)) {
          d <- pg_res$cluster_assignments[, c("cluster", "type", "organ", "scores", "n_overlap"), drop = FALSE]
          colnames(d) <- c("Cluster", "PanglaoDB", "PanglaoDB_organ", "PanglaoDB_score", "PanglaoDB_overlap")
          d$PanglaoDB_score <- round(d$PanglaoDB_score, 3)
          d
        } else NULL

        if (!is.null(sc_tbl) && !is.null(pg_tbl)) {
          tbl <- merge(sc_tbl, pg_tbl, by = "Cluster", all = TRUE)
          # Add consensus column
          md <- values$active_dataset@meta.data
          if ("cell_type_consensus" %in% colnames(md)) {
            cl_col <- if ("seurat_clusters" %in% colnames(md)) "seurat_clusters" else colnames(md)[1]
            cons <- unique(data.frame(
              Cluster = as.character(md[[cl_col]]),
              Consensus = md[["cell_type_consensus"]],
              stringsAsFactors = FALSE
            ))
            tbl <- merge(tbl, cons, by = "Cluster", all.x = TRUE)
          }
        } else {
          tbl <- sc_tbl %||% pg_tbl
        }

        if (exists("register_csv_export", mode = "function")) {
          register_csv_export(values, data = tbl, description = "Annotation Results",
                              source_module = "Annotation")
        }
        DT::datatable(tbl, rownames = FALSE,
                      options = list(pageLength = 25, scrollX = TRUE),
                      caption = "Cell type assignments per cluster")
      } else {
        df <- values$active_dataset@meta.data
        keep <- c("seurat_clusters", "cell_type_consensus", "cell_type_sctype",
                  "cell_type_panglao", "cell_type_clustifyr",
                  "sctype_score", "panglao_score", "clustifyr_score")
        keep <- keep[keep %in% colnames(df)]
        if (length(keep) == 0) return(DT::datatable(data.frame(Message = "Run annotation first")))
        tbl <- unique(df[, keep, drop = FALSE])
        DT::datatable(tbl, rownames = FALSE,
                      options = list(pageLength = 25, scrollX = TRUE))
      }
    })

    output$dl_annotations <- downloadHandler(
      filename = function() paste0("annotations_", Sys.Date(), ".csv"),
      content = function(file) {
        df <- values$active_dataset@meta.data
        keep <- c("cell_type_consensus", "cell_type_sctype", "cell_type_panglao",
                  "cell_type_clustifyr",
                  "sctype_tissue", "panglao_organ", "sctype_score", "panglao_score",
                  "clustifyr_score", "annotation_consensus_source", "seurat_clusters")
        keep <- keep[keep %in% colnames(df)]
        write.csv(df[, keep, drop = FALSE], file, row.names = TRUE)
      }
    )
  })
}

# =============================================================================
# UI: Marker Validation
# =============================================================================

markerValidationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Marker / Annotation Validation"),
          fluidRow(
            column(4,
              textInput(ns("marker_genes"), "Markers (comma-separated)",
                        placeholder = "CD3D, CD4, MS4A1")
            ),
            column(4,
              selectInput(ns("plot_type"), "Plot",
                choices = c("DotPlot", "Heatmap", "Violin", "Feature Plot (UMAP)"),
                selected = "DotPlot")
            ),
            column(4,
              selectInput(ns("groupby_val"), "Group by", choices = NULL)
            )
          ),
          fluidRow(
            column(4,
              checkboxInput(ns("order_by_mean"),
                            "Order groups by mean marker expression", TRUE)
            ),
            column(4,
              conditionalPanel(
                condition = "input.plot_type == 'Feature Plot (UMAP)'",
                ns = ns,
                numericInput(ns("feature_pt_size"), "Point size",
                             value = 1.5, min = 0.5, max = 5, step = 0.5)
              )
            ),
            column(4,
              actionButton(ns("run_validation"), "Run Validation",
                           class = "btn-primary", icon = icon("check"))
            )
          )
        )
      )
    ),
    fluidRow(
      column(12, uiOutput(ns("validation_output")))
    )
  )
}

# =============================================================================
# Server: Marker Validation
# =============================================================================

markerValidationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    observe({
      req(values$active_dataset)
      meta_cols   <- colnames(values$active_dataset@meta.data)
      factor_cols <- meta_cols[sapply(values$active_dataset@meta.data,
                                       function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "groupby_val",
                        choices = c(factor_cols, meta_cols),
                        selected = if ("cell_type_consensus" %in% meta_cols)
                                     "cell_type_consensus"
                                   else "seurat_clusters")
    })

    observeEvent(input$run_validation, {
      req(values$active_dataset)
      genes <- trimws(unlist(strsplit(input$marker_genes, ",")))
      genes <- genes[genes %in% rownames(values$active_dataset)]
      if (!length(genes)) {
        showNotification("No valid markers found in dataset.", type = "error")
        return(NULL)
      }

      # Dynamic output UI
      output$validation_output <- renderUI({
        if (input$plot_type == "Feature Plot (UMAP)") {
          panels <- lapply(seq_along(genes), function(i) {
            column(6, plotlyOutput(session$ns(paste0("feature_gene_", i)), height = "400px"))
          })
          all_panels <- c(
            list(column(6, plotlyOutput(session$ns("feature_clusters"), height = "400px"))),
            panels
          )
          tagList(
            fluidRow(do.call(tagList, all_panels)),
            br(),
            downloadButton(session$ns("dl_avg_expr_table"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(session$ns("avg_expr_table"))
          )
        } else {
          tagList(
            plotlyOutput(session$ns("val_plot"), height = "650px"),
            br(),
            downloadButton(session$ns("dl_avg_expr_table"), "CSV", class = "btn-xs btn-default"),
            DT::dataTableOutput(session$ns("avg_expr_table"))
          )
        }
      })

      # Helper: prepare object with joined layers
      .prep_obj <- function(obj) {
        assay_name <- DefaultAssay(obj)
        a <- obj[[assay_name]]
        if (inherits(a, "Assay5") && length(SeuratObject::Layers(a)) > 1) {
          tryCatch({ obj[[assay_name]] <- JoinLayers(a) }, error = function(e) NULL)
          a <- obj[[assay_name]]
        }
        if (inherits(a, "Assay5") && !"data" %in% Layers(a) && "counts" %in% Layers(a)) {
          obj <- NormalizeData(obj, verbose = FALSE)
        }
        obj
      }

      # Helper: best reduction
      .get_reduction <- function(obj) {
        .best_reduction(obj)
      }

      # Cluster reference plot for feature mode
      output$feature_clusters <- renderPlotly({
        req(values$active_dataset)
        tryCatch({
          obj <- values$active_dataset
          red <- .get_reduction(obj)
          if (is.null(red)) return(plotly_empty() %>% layout(title = "No reduction available"))

          emb    <- Embeddings(obj, red)[, 1:2]
          groups <- as.character(obj@meta.data[[input$groupby_val]])
          n_groups <- length(unique(groups))
          pal <- .get_palette(n_groups)

          p <- plot_ly(
            x = emb[, 1], y = emb[, 2],
            color = groups, colors = pal,
            type = "scatter", mode = "markers",
            marker = list(size = input$feature_pt_size, opacity = 0.8),
            hoverinfo = "text",
            text = paste0("Cell: ", colnames(obj), "\n", input$groupby_val, ": ", groups)
          ) %>% layout(
            title = list(text = paste("Clusters:", input$groupby_val), font = list(size = 13)),
            xaxis = list(title = paste0(toupper(red), "_1")),
            yaxis = list(title = paste0(toupper(red), "_2")),
            autosize = TRUE,
            legend = list(font = list(size = 9), itemsizing = "constant")
          )
          if (exists("capture_figure", mode = "function")) {
            capture_figure(values, p, "Annotation",
                           figure_type = "Cluster reference plot for validation",
                           label = "Feature Clusters", source_module = "Annotation")
          }
          p
        }, error = function(e) plotly_empty() %>% layout(title = paste("Error:", e$message)))
      })

      # Individual gene feature plots
      for (gi in seq_along(genes)) {
        local({
          gene_i <- genes[gi]
          output[[paste0("feature_gene_", gi)]] <- renderPlotly({
            req(values$active_dataset)
            tryCatch({
              obj <- .prep_obj(values$active_dataset)
              red <- .get_reduction(obj)
              if (is.null(red)) return(plotly_empty() %>% layout(title = "No reduction available"))

              emb       <- Embeddings(obj, red)[, 1:2]
              expr_vals <- as.numeric(.safe_get_assay_data(obj, "data")[gene_i, ])

              ord <- order(expr_vals)
              plot_ly(
                x = emb[ord, 1], y = emb[ord, 2],
                type = "scatter", mode = "markers",
                marker = list(
                  size = input$feature_pt_size,
                  opacity = 0.85,
                  color = expr_vals[ord],
                  colorscale = list(
                    c(0, "#D3D3D3"), c(0.15, "#C6DBEF"), c(0.4, "#6BAED6"),
                    c(0.65, "#2171B5"), c(0.85, "#08519C"), c(1, "#6A0DAD")
                  ),
                  colorbar = list(title = "Expr", len = 0.6),
                  showscale = TRUE
                ),
                hoverinfo = "text",
                text = paste0("Cell: ", colnames(obj)[ord], "\n",
                              gene_i, ": ", round(expr_vals[ord], 2))
              ) %>% layout(
                title = list(text = gene_i, font = list(size = 14, color = "#2C3E50")),
                xaxis = list(title = paste0(toupper(red), "_1")),
                yaxis = list(title = paste0(toupper(red), "_2")),
                autosize = TRUE
              )
            }, error = function(e) plotly_empty() %>% layout(title = paste(gene_i, "error:", e$message)))
          })
        })
      }

      # Standard plot (DotPlot / Heatmap / Violin)
      output$val_plot <- renderPlotly({
        req(values$active_dataset)
        if (input$plot_type == "Feature Plot (UMAP)") return(NULL)
        tryCatch({
          obj        <- .prep_obj(values$active_dataset)
          assay_name <- DefaultAssay(obj)

          if (input$plot_type == "DotPlot") {
            p <- DotPlot(obj, features = genes, group.by = input$groupby_val) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1))
            p <- ggplotly(p)
          } else if (input$plot_type == "Heatmap") {
            avg <- AverageExpression(obj, features = genes, group.by = input$groupby_val,
                                     assays = assay_name)
            avg_data <- avg[[assay_name]]
            if (is.null(avg_data)) avg_data <- avg[[1]]
            avg_data <- as.matrix(avg_data)
            valid_genes <- genes[genes %in% rownames(avg_data)]
            avg_data <- avg_data[valid_genes, , drop = FALSE]
            df <- as.data.frame(avg_data)
            df$gene <- rownames(df)
            df <- tidyr::pivot_longer(df, cols = -gene, names_to = "group", values_to = "expr")
            p <- ggplot(df, aes(x = group, y = gene, fill = expr)) +
              geom_tile() + scale_fill_viridis_c() +
              theme_minimal() +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
              labs(x = "Group", y = "Gene", fill = "Expr")
            p <- ggplotly(p)
          } else {
            # Violin
            p <- VlnPlot(obj, features = genes, group.by = input$groupby_val, pt.size = 0)
            p <- ggplotly(p)
          }
          if (exists("capture_figure", mode = "function")) {
            capture_figure(values, p, "Annotation",
                           figure_type = paste(input$plot_type, "validation plot"),
                           label = "Validation Plot", source_module = "Annotation")
          }
          p
        }, error = function(e) {
          plotly_empty() %>% layout(title = paste("Validation plot error:", e$message))
        })
      })

      # Average expression table
      output$avg_expr_table <- DT::renderDT({
        tryCatch({
          obj <- .prep_obj(values$active_dataset)
          assay_name <- DefaultAssay(obj)
          avg <- AverageExpression(obj, features = genes, group.by = input$groupby_val,
                                   assays = assay_name)
          avg_data <- avg[[assay_name]]
          if (is.null(avg_data)) avg_data <- avg[[1]]
          avg_data <- as.matrix(avg_data)
          if (isTRUE(input$order_by_mean)) {
            ord <- order(colMeans(avg_data, na.rm = TRUE), decreasing = TRUE)
            avg_data <- avg_data[, ord, drop = FALSE]
          }
          DT::datatable(round(avg_data, 3), options = list(pageLength = 10, scrollX = TRUE))
        }, error = function(e) {
          DT::datatable(data.frame(Error = e$message))
        })
      })

      output$dl_avg_expr_table <- downloadHandler(
        filename = function() paste0("avg_expression_", Sys.Date(), ".csv"),
        content = function(file) {
          tryCatch({
            obj <- .prep_obj(values$active_dataset)
            cur_genes <- unique(trimws(unlist(strsplit(input$marker_genes, ","))))
            cur_genes <- cur_genes[cur_genes %in% rownames(obj)]
            avg <- AverageExpression(obj, features = cur_genes, group.by = input$groupby_val,
                                     assays = DefaultAssay(obj))
            avg_data <- avg[[DefaultAssay(obj)]]
            if (is.null(avg_data)) avg_data <- avg[[1]]
            write.csv(round(as.matrix(avg_data), 3), file, row.names = TRUE)
          }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
        }
      )
    })
  })
}

# =============================================================================
# UI: Manual Annotation (Rename / Merge)
# =============================================================================

manualAnnotationUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Manual Annotation (Rename/Merge)"),
    fluidRow(
      column(4, selectInput(ns("source_column"), "Source labels", choices = NULL)),
      column(4, textInput(ns("from_label"), "From label")),
      column(4, textInput(ns("to_label"), "To (new) label"))
    ),
    fluidRow(
      column(3, actionButton(ns("rename_one"), "Rename label", class = "btn-warning")),
      column(3, actionButton(ns("merge_into"), "Merge into label", class = "btn-danger")),
      column(3, actionButton(ns("undo_last"), "Undo Last", class = "btn-default")),
      column(3, downloadButton(ns("dl_manual_csv"), "Export CSV"))
    ),
    br(),
    helpText("Tip: You can write results back to 'cell_type_consensus' for downstream use."),
    hr(),
    fluidRow(
      column(6, uiOutput(ns("label_summary_ui"))),
      column(6, plotlyOutput(ns("manual_dimplot"), height = "400px"))
    )
  )
}

# =============================================================================
# Server: Manual Annotation
# =============================================================================

manualAnnotationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    undo_stack <- reactiveVal(list())

    observe({
      req(values$active_dataset)
      meta_cols <- colnames(values$active_dataset@meta.data)
      updateSelectInput(session, "source_column", choices = meta_cols,
                        selected = if ("cell_type_consensus" %in% meta_cols)
                                     "cell_type_consensus"
                                   else meta_cols[1])
    })

    # Save current state for undo
    .save_undo <- function() {
      if (is.null(values$active_dataset)) return()
      col <- input$source_column
      if (is.null(col) || !col %in% colnames(values$active_dataset@meta.data)) return()
      old_vals <- as.character(values$active_dataset@meta.data[[col]])
      stack <- undo_stack()
      stack[[length(stack) + 1]] <- list(col = col, vals = old_vals)
      # Keep max 10 undo steps
      if (length(stack) > 10) stack <- stack[(length(stack) - 9):length(stack)]
      undo_stack(stack)
    }

    observeEvent(input$rename_one, {
      req(values$active_dataset, input$source_column,
          nzchar(input$from_label), nzchar(input$to_label))
      .save_undo()
      v <- as.character(values$active_dataset@meta.data[[input$source_column]])
      n_changed <- sum(v == input$from_label)
      v[v == input$from_label] <- input$to_label
      values$active_dataset@meta.data[[input$source_column]] <- v
      showNotification(
        sprintf("Renamed '%s' -> '%s' (%d cells).", input$from_label, input$to_label, n_changed),
        type = "message"
      )
    })

    observeEvent(input$merge_into, {
      req(values$active_dataset, input$source_column,
          nzchar(input$from_label), nzchar(input$to_label))
      .save_undo()
      v <- as.character(values$active_dataset@meta.data[[input$source_column]])
      n_changed <- sum(v == input$from_label)
      v[v == input$from_label] <- input$to_label
      values$active_dataset@meta.data[[input$source_column]] <- v
      showNotification(
        sprintf("Merged '%s' into '%s' (%d cells).", input$from_label, input$to_label, n_changed),
        type = "warning"
      )
    })

    observeEvent(input$undo_last, {
      stack <- undo_stack()
      if (length(stack) == 0) {
        showNotification("Nothing to undo.", type = "warning")
        return()
      }
      last <- stack[[length(stack)]]
      stack[[length(stack)]] <- NULL
      undo_stack(stack)
      values$active_dataset@meta.data[[last$col]] <- last$vals
      showNotification("Undo successful.", type = "message")
    })

    # Label summary
    output$label_summary_ui <- renderUI({
      req(values$active_dataset, input$source_column)
      col <- input$source_column
      if (!col %in% colnames(values$active_dataset@meta.data)) return(NULL)
      v <- as.character(values$active_dataset@meta.data[[col]])
      tbl <- sort(table(v), decreasing = TRUE)
      df <- data.frame(Label = names(tbl), Count = as.integer(tbl), stringsAsFactors = FALSE)
      tags$div(
        h5("Current Labels"),
        DT::renderDT(DT::datatable(df, options = list(pageLength = 15, dom = "tip"), rownames = FALSE))
      )
    })

    # Live preview dimplot
    output$manual_dimplot <- renderPlotly({
      req(values$active_dataset, input$source_column)
      obj <- values$active_dataset
      col <- input$source_column
      if (!col %in% colnames(obj@meta.data)) return(plotly_empty())
      red <- .best_reduction(obj)
      if (is.null(red)) return(plotly_empty())

      emb  <- Embeddings(obj, red)[, 1:2, drop = FALSE]
      labs <- as.character(obj@meta.data[[col]])

      plot_ly(
        x = emb[, 1], y = emb[, 2],
        type = "scatter", mode = "markers",
        marker = list(size = 2, opacity = 0.7),
        color = labs, colors = .get_palette(length(unique(labs))),
        hoverinfo = "text",
        text = paste0(labs, "\n", colnames(obj))
      ) %>% layout(
        title = paste(toupper(red), "-", col),
        xaxis = list(title = paste0(toupper(red), "1")),
        yaxis = list(title = paste0(toupper(red), "2")),
        showlegend = TRUE
      )
    })

    # Download handler
    output$dl_manual_csv <- downloadHandler(
      filename = function() paste0("manual_annotations_", Sys.Date(), ".csv"),
      content = function(file) {
        req(values$active_dataset, input$source_column)
        df <- data.frame(
          cell_barcode = colnames(values$active_dataset),
          label = as.character(values$active_dataset@meta.data[[input$source_column]]),
          stringsAsFactors = FALSE
        )
        colnames(df)[2] <- input$source_column
        write.csv(df, file, row.names = FALSE)
      }
    )
  })
}

# =============================================================================
# UI: Cluster Editor (rename, recolor, merge with live preview)
# =============================================================================

clusterEditorUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Cluster Editor"),
    fluidRow(
      column(6, selectInput(ns("editor_column"), "Annotation Column", choices = NULL)),
      column(6, selectInput(ns("editor_reduction"), "Preview Reduction", choices = NULL))
    ),
    uiOutput(ns("cluster_rows")),
    fluidRow(
      column(4, actionButton(ns("apply_changes"), "Apply All Changes",
                             class = "btn-primary", width = "100%")),
      column(4, actionButton(ns("merge_selected"), "Merge Checked",
                             class = "btn-warning", width = "100%")),
      column(4, textInput(ns("merge_target"), "Merge into name", placeholder = "New merged label"))
    ),
    hr(),
    h5("Preview"),
    plotlyOutput(ns("editor_dimplot"), height = "500px"),
    br(),
    fluidRow(
      column(4, downloadButton(ns("dl_editor_csv"), "Download Labels CSV")),
      column(4, downloadButton(ns("dl_editor_colors"), "Download Color Map CSV")),
      column(4, actionButton(ns("reset_colors"), "Reset Colors", class = "btn-default"))
    )
  )
}

# =============================================================================
# Server: Cluster Editor
# =============================================================================

clusterEditorServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    custom_colors <- reactiveVal(list())

    observe({
      req(values$active_dataset)
      md <- colnames(values$active_dataset@meta.data)
      char_cols <- md[sapply(values$active_dataset@meta.data,
                              function(x) is.factor(x) || is.character(x))]
      if (length(char_cols) == 0) char_cols <- md[1]
      selected <- if ("cell_type_consensus" %in% char_cols) "cell_type_consensus"
                  else if ("seurat_clusters" %in% char_cols) "seurat_clusters"
                  else char_cols[1]
      updateSelectInput(session, "editor_column", choices = char_cols, selected = selected)

      reds <- names(values$active_dataset@reductions)
      best <- if ("umap" %in% reds) "umap" else if ("tsne" %in% reds) "tsne" else reds[1]
      updateSelectInput(session, "editor_reduction", choices = reds, selected = best)
    })

    # Dynamic rows: one per cluster
    output$cluster_rows <- renderUI({
      req(values$active_dataset, input$editor_column)
      col <- input$editor_column
      if (!col %in% colnames(values$active_dataset@meta.data)) return(NULL)
      labels  <- sort(unique(as.character(values$active_dataset@meta.data[[col]])))
      n       <- length(labels)
      palette <- .get_palette(n)
      cc      <- isolate(custom_colors())

      rows <- lapply(seq_along(labels), function(i) {
        lbl         <- labels[i]
        default_col <- if (!is.null(cc[[lbl]])) cc[[lbl]] else palette[i]
        fluidRow(
          column(1, checkboxInput(ns(paste0("merge_", i)), label = NULL, value = FALSE)),
          column(2, tags$span(style = "line-height:34px; font-weight:bold;", lbl)),
          column(4, textInput(ns(paste0("rename_", i)), label = NULL, value = lbl, placeholder = lbl)),
          column(3, tags$input(
            type = "color", id = ns(paste0("color_", i)),
            value = default_col,
            style = "width:60px;height:34px;padding:2px;cursor:pointer;border:1px solid #ccc;border-radius:4px;",
            onchange = sprintf("Shiny.setInputValue('%s', this.value)", ns(paste0("color_", i)))
          )),
          column(2, tags$small(style = "line-height:34px;",
            paste0("n=", sum(values$active_dataset@meta.data[[col]] == lbl))))
        )
      })
      tagList(
        fluidRow(
          column(1, tags$b("Sel")), column(2, tags$b("Current")),
          column(4, tags$b("New Name")), column(3, tags$b("Color")),
          column(2, tags$b("Cells"))
        ),
        rows
      )
    })

    # Apply changes
    observeEvent(input$apply_changes, {
      req(values$active_dataset, input$editor_column)
      col    <- input$editor_column
      labels <- sort(unique(as.character(values$active_dataset@meta.data[[col]])))
      v      <- as.character(values$active_dataset@meta.data[[col]])
      new_colors <- list()

      for (i in seq_along(labels)) {
        new_name <- input[[paste0("rename_", i)]]
        new_col  <- input[[paste0("color_", i)]]
        if (!is.null(new_name) && nzchar(new_name) && new_name != labels[i]) {
          v[v == labels[i]] <- new_name
          if (!is.null(new_col)) new_colors[[new_name]] <- new_col
        } else {
          if (!is.null(new_col)) new_colors[[labels[i]]] <- new_col
        }
      }

      values$active_dataset@meta.data[[col]] <- v
      custom_colors(new_colors)
      if (!is.null(values$cluster_colors)) {
        values$cluster_colors <- unlist(new_colors)
      }
      showNotification("Cluster names and colors updated.", type = "message")
    })

    # Merge selected
    observeEvent(input$merge_selected, {
      req(values$active_dataset, input$editor_column, nzchar(input$merge_target))
      col    <- input$editor_column
      labels <- sort(unique(as.character(values$active_dataset@meta.data[[col]])))
      v      <- as.character(values$active_dataset@meta.data[[col]])
      merged <- 0

      for (i in seq_along(labels)) {
        if (isTRUE(input[[paste0("merge_", i)]])) {
          v[v == labels[i]] <- input$merge_target
          merged <- merged + 1
        }
      }

      if (merged > 0) {
        values$active_dataset@meta.data[[col]] <- v
        showNotification(sprintf("Merged %d clusters into '%s'.", merged, input$merge_target),
                         type = "message")
      } else {
        showNotification("No clusters selected for merging.", type = "warning")
      }
    })

    # Reset colors
    observeEvent(input$reset_colors, {
      custom_colors(list())
      showNotification("Colors reset to defaults.", type = "message")
    })

    # Live UMAP preview
    output$editor_dimplot <- renderPlotly({
      req(values$active_dataset, input$editor_column, input$editor_reduction)
      tryCatch({
        obj <- values$active_dataset
        col <- input$editor_column
        red <- input$editor_reduction
        if (!red %in% names(obj@reductions)) return(plotly_empty())
        if (!col %in% colnames(obj@meta.data)) return(plotly_empty())

        embedding   <- Embeddings(obj, red)[, 1:2]
        labels      <- as.character(obj@meta.data[[col]])
        unique_labs <- sort(unique(labels))
        n           <- length(unique_labs)

        cc      <- custom_colors()
        palette <- .get_palette(n)
        colors  <- setNames(palette, unique_labs)
        for (lbl in names(cc)) {
          if (lbl %in% names(colors)) colors[[lbl]] <- cc[[lbl]]
        }

        df <- data.frame(
          Dim1 = embedding[, 1], Dim2 = embedding[, 2],
          Label = labels, stringsAsFactors = FALSE
        )

        p <- plot_ly()
        for (lbl in unique_labs) {
          sub_df <- df[df$Label == lbl, ]
          p <- p %>% add_trace(
            data = sub_df, x = ~Dim1, y = ~Dim2,
            type = "scatter", mode = "markers",
            marker = list(size = 3, color = colors[[lbl]], opacity = 0.7),
            name = lbl,
            hoverinfo = "text",
            text = ~paste0(lbl, "\n", round(Dim1, 2), ", ", round(Dim2, 2))
          )
        }

        p <- p %>% layout(
          title = paste(toupper(red), "- Cluster Editor Preview"),
          xaxis = list(title = colnames(embedding)[1]),
          yaxis = list(title = colnames(embedding)[2]),
          autosize = TRUE
        )
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Annotation",
                         figure_type = "Cluster editor preview dimplot",
                         label = "Editor DimPlot", source_module = "Annotation")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Preview error:", e$message))
      })
    })

    # Download labels CSV
    output$dl_editor_csv <- downloadHandler(
      filename = function() paste0("cluster_labels_", Sys.Date(), ".csv"),
      content = function(file) {
        req(values$active_dataset, input$editor_column)
        df <- data.frame(
          cell = colnames(values$active_dataset),
          label = as.character(values$active_dataset@meta.data[[input$editor_column]]),
          stringsAsFactors = FALSE
        )
        colnames(df)[2] <- input$editor_column
        write.csv(df, file, row.names = FALSE)
      }
    )

    # Download color map CSV
    output$dl_editor_colors <- downloadHandler(
      filename = function() paste0("color_map_", Sys.Date(), ".csv"),
      content = function(file) {
        cc <- custom_colors()
        if (length(cc) == 0) {
          write.csv(data.frame(label = character(), color = character()), file, row.names = FALSE)
        } else {
          write.csv(data.frame(label = names(cc), color = unlist(cc)), file, row.names = FALSE)
        }
      }
    )
  })
}

# =============================================================================
# Helper: Find top marker genes per cluster (exported)
# =============================================================================

find_top_markers_per_cluster <- function(seurat_obj, cluster_col, top_n = 20) {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column '", cluster_col, "' not found. Run clustering first.")
  }

  Idents(seurat_obj) <- cluster_col
  clusters <- sort(unique(Idents(seurat_obj)))

  if (length(clusters) < 2) {
    stop("Need at least 2 clusters to find markers, but found ", length(clusters), ".")
  }

  # Downsample before joining layers
  max_cells_per_cluster <- 200
  n_cells <- ncol(seurat_obj)
  if (n_cells > max_cells_per_cluster * length(clusters)) {
    set.seed(42)
    cells_to_keep <- c()
    for (cl in clusters) {
      cl_cells <- WhichCells(seurat_obj, idents = cl)
      if (length(cl_cells) > max_cells_per_cluster) {
        cl_cells <- sample(cl_cells, max_cells_per_cluster)
      }
      cells_to_keep <- c(cells_to_keep, cl_cells)
    }
    marker_obj <- subset(seurat_obj, cells = cells_to_keep)
  } else {
    marker_obj <- seurat_obj
  }

  marker_obj <- .safe_prepare_for_de(marker_obj)
  markers_list <- list()

  # Attempt FindAllMarkers
  markers_df <- tryCatch(
    FindAllMarkers(marker_obj, only.pos = TRUE, min.pct = 0.1,
                   logfc.threshold = 0.25, verbose = FALSE),
    warning = function(w) { invokeRestart("muffleWarning") },
    error = function(e) NULL
  )

  if (is.null(markers_df) || nrow(markers_df) == 0) {
    markers_df <- tryCatch(
      FindAllMarkers(marker_obj, only.pos = TRUE, min.pct = 0.01,
                     logfc.threshold = 0.1, verbose = FALSE),
      warning = function(w) { invokeRestart("muffleWarning") },
      error = function(e) NULL
    )
  }

  if (!is.null(markers_df) && nrow(markers_df) > 0) {
    for (cl in clusters) {
      cl_char <- as.character(cl)
      sub <- markers_df[markers_df$cluster == cl_char, , drop = FALSE]
      if (nrow(sub) == 0) next
      sub <- sub[order(sub$avg_log2FC, decreasing = TRUE), ]
      sub <- utils::head(sub, top_n)
      genes_with_fc <- setNames(sub$avg_log2FC, sub$gene)
      markers_list[[cl_char]] <- genes_with_fc
    }
  } else {
    # Fallback: average expression fold change
    avg_expr <- tryCatch(
      AverageExpression(marker_obj, assays = .safe_default_assay(marker_obj), layer = "data"),
      error = function(e) AverageExpression(marker_obj, assays = .safe_default_assay(marker_obj))
    )
    assay_name <- names(avg_expr)[1]
    avg_mat <- as.matrix(avg_expr[[assay_name]])
    global_mean <- rowMeans(avg_mat)
    global_mean[global_mean == 0] <- 1e-9

    for (cl in clusters) {
      cl_char <- as.character(cl)
      if (!(cl_char %in% colnames(avg_mat))) next
      fc <- avg_mat[, cl_char] / global_mean
      fc <- sort(fc, decreasing = TRUE)
      fc <- utils::head(fc, top_n)
      markers_list[[cl_char]] <- log2(fc)
    }
  }

  markers_list
}

# =============================================================================
# Helper: Parse AI response JSON (exported)
# =============================================================================

parse_ai_response <- function(ai_text) {
  json_str <- paste(ai_text, collapse = "\n")

  # Strip markdown code fences
  json_str <- sub("^.*```(?:json)?\\s*\\n?", "", json_str, perl = TRUE)
  json_str <- sub("\\s*```\\s*$", "", json_str, perl = TRUE)

  # Find JSON array start
  arr_start <- regexpr("\\[", json_str)
  if (arr_start > 0) {
    json_str <- substring(json_str, arr_start)
  }

  # Handle truncated JSON
  n_open_bracket  <- nchar(gsub("[^\\[]", "", json_str))
  n_close_bracket <- nchar(gsub("[^\\]]", "", json_str))
  n_open_brace    <- nchar(gsub("[^{]", "", json_str))
  n_close_brace   <- nchar(gsub("[^}]", "", json_str))

  if (n_open_bracket > n_close_bracket || n_open_brace > n_close_brace) {
    last_complete <- regexpr(".*\\}", json_str, perl = TRUE)
    if (last_complete > 0) {
      json_str <- regmatches(json_str, last_complete)
      json_str <- sub(",\\s*$", "", json_str)
      json_str <- paste0(json_str, "]")
    }
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(json_str, simplifyDataFrame = TRUE),
    error = function(e) {
      stop("Failed to parse AI response as JSON. Raw response:\n", substr(ai_text, 1, 500))
    }
  )

  if (!is.data.frame(parsed)) {
    if (is.list(parsed)) {
      parsed <- do.call(rbind, lapply(parsed, as.data.frame, stringsAsFactors = FALSE))
    } else {
      stop("AI response did not produce a valid data.frame.")
    }
  }

  expected_cols <- c("cluster", "cell_type", "confidence", "reasoning", "alternative", "key_markers")
  for (col in expected_cols) {
    if (!(col %in% colnames(parsed))) {
      parsed[[col]] <- NA
    }
  }

  # Convert key_markers list column to string
  if ("key_markers" %in% colnames(parsed) && is.list(parsed$key_markers)) {
    parsed$key_markers <- vapply(parsed$key_markers, function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(x, collapse = ", ")
    }, character(1))
  }

  parsed$cluster    <- as.character(parsed$cluster)
  parsed$confidence <- as.numeric(parsed$confidence)
  parsed
}

# =============================================================================
# UI: AI-Powered Annotation (Claude API)
# =============================================================================

aiAnnotationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(icon("robot"), " AI Annotation"),
          p("Use Claude AI to identify cell types from top marker genes per cluster. ",
            "API key and model are configured in ",
            tags$b("Data Management > AI Config"), "."),

          fluidRow(
            column(3,
              selectInput(ns("ai_cluster_col"), "Cluster Column", choices = NULL)
            ),
            column(3,
              selectInput(ns("ai_species"), "Species",
                choices = c("Human" = "human", "Mouse" = "mouse",
                            "Rat" = "rat", "Zebrafish" = "zebrafish"),
                selected = "human")
            ),
            column(3,
              textInput(ns("ai_annotation_name"), "Annotation Column Name",
                value = "cell_type_ai", placeholder = "e.g., cell_type_ai_v2")
            ),
            column(3,
              textInput(ns("tissue_context"), "Tissue Context",
                placeholder = "e.g., PBMC, brain, tumor, liver")
            )
          ),

          fluidRow(
            column(3,
              numericInput(ns("top_n_genes"), "Top N Genes per Cluster",
                value = 20, min = 5, max = 50, step = 5)
            ),
            column(3,
              checkboxInput(ns("include_marker_context"),
                            "Include known marker context", TRUE),
              checkboxInput(ns("validate_with_markers"),
                            "Cross-validate with marker DB", TRUE)
            ),
            column(3,
              checkboxInput(ns("recluster_for_annotation"),
                            "Re-cluster before annotation", FALSE),
              conditionalPanel(
                condition = "input.recluster_for_annotation",
                ns = ns,
                sliderInput(ns("annotation_resolution"), "Clustering Resolution",
                  min = 0.1, max = 2.0, value = 0.5, step = 0.1)
              )
            ),
            column(3,
              div(style = "margin-top: 5px;",
                actionButton(ns("run_ai_annotation"), "Run AI Annotation",
                  icon = icon("robot"), class = "btn-primary btn-lg", width = "100%")
              )
            )
          ),

          # Conversation history toggle
          fluidRow(
            column(12,
              checkboxInput(ns("show_conversation"), "Show AI conversation details", FALSE),
              conditionalPanel(
                condition = "input.show_conversation",
                ns = ns,
                wellPanel(
                  style = "max-height: 300px; overflow-y: auto; background: #f9f9f9;",
                  h5("AI Conversation History"),
                  uiOutput(ns("ai_conversation_log"))
                )
              )
            )
          )
        )
      )
    ),

    fluidRow(
      column(12,
        uiOutput(ns("ai_results_ui"))
      )
    )
  )
}

# =============================================================================
# Server: AI-Powered Annotation (Claude API)
# =============================================================================

aiAnnotationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    ai_results     <- reactiveVal(NULL)
    last_annot_col <- reactiveVal("cell_type_ai")
    ai_conv_log    <- reactiveVal(list())

    output$ai_results_ui <- renderUI({
      req(ai_results())
      ns <- session$ns
      wellPanel(
        h4("AI Annotation Results"),

        tabsetPanel(
          tabPanel("Predictions",
            div(style = "margin-bottom: 25px;",
              fluidRow(
                column(6,
                  plotlyOutput(ns("ai_dimplot"), height = "550px")
                ),
                column(6,
                  DT::DTOutput(ns("ai_results_table"))
                )
              )
            )
          ),
          tabPanel("AI Reasoning",
            div(style = "margin-bottom: 25px;",
              fluidRow(
                column(12,
                  downloadButton(ns("dl_ai_reasoning_table"), "CSV", class = "btn-xs btn-default"),
                  DT::DTOutput(ns("ai_reasoning_table"))
                )
              )
            )
          ),
          tabPanel("Comparison",
            div(style = "margin-bottom: 25px;",
              fluidRow(
                column(6,
                  h5("AI vs ScType"),
                  plotlyOutput(ns("ai_comparison_plot"), height = "450px")
                ),
                column(6,
                  h5("Agreement Summary"),
                  downloadButton(ns("dl_ai_comparison_table"), "CSV", class = "btn-xs btn-default"),
                  DT::DTOutput(ns("ai_comparison_table"))
                )
              )
            )
          ),
          tabPanel("Marker Heatmap",
            div(style = "margin-bottom: 25px;",
              plotlyOutput(ns("ai_marker_heatmap"), height = "600px")
            )
          ),
          tabPanel("3D View",
            div(style = "margin-bottom: 25px;",
              plotlyOutput(ns("ai_3d_umap"), height = "600px"),
              helpText("Requires 3D UMAP reduction ('umap3d') to be computed.")
            )
          )
        ),

        br(),
        fluidRow(
          column(6,
            downloadButton(ns("dl_ai_annotations"), "Download AI Annotations (CSV)",
              class = "btn-success", width = "100%")
          ),
          column(6,
            actionButton(ns("apply_ai_labels"), "Apply AI Labels as Consensus",
              class = "btn-warning", width = "100%")
          )
        )
      )
    })

    # Update cluster column choices
    observe({
      req(values$active_dataset)
      md   <- values$active_dataset@meta.data
      cand <- colnames(md)[sapply(md, function(x) is.factor(x) || is.character(x) || is.integer(x))]
      cand <- union(cand, grep("cluster|leiden|louvain|slm|res\\.",
                                colnames(md), value = TRUE, ignore.case = TRUE))
      cand <- sort(unique(cand))
      sel  <- if ("seurat_clusters" %in% cand) "seurat_clusters" else cand[1]
      updateSelectInput(session, "ai_cluster_col", choices = cand, selected = sel)
    })

    # AI conversation log
    output$ai_conversation_log <- renderUI({
      log <- ai_conv_log()
      if (length(log) == 0) return(helpText("No conversation yet."))
      tags$div(
        lapply(log, function(entry) {
          tags$div(
            style = paste0("margin-bottom: 8px; padding: 6px; border-left: 3px solid ",
                            if (entry$role == "user") "#007bff" else "#28a745", ";"),
            tags$strong(entry$role), ": ",
            tags$span(style = "white-space: pre-wrap; font-size: 11px;",
                       substr(entry$content, 1, 500),
                       if (nchar(entry$content) > 500) "..." else "")
          )
        })
      )
    })

    # ---------------------------------------------------------------
    # Run AI annotation
    # ---------------------------------------------------------------
    observeEvent(input$run_ai_annotation, {
      api_key <- values$ai_api_key
      if (is.null(api_key) || !nzchar(trimws(api_key))) {
        showNotification(
          "Please set your Claude API key in Data Management > AI Config first.",
          type = "error", duration = 8
        )
        return(NULL)
      }

      if (!requireNamespace("httr2", quietly = TRUE)) {
        showNotification(
          "The 'httr2' package is required. Install: install.packages('httr2')",
          type = "error", duration = 10
        )
        return(NULL)
      }
      if (!requireNamespace("jsonlite", quietly = TRUE)) {
        showNotification(
          "The 'jsonlite' package is required. Install: install.packages('jsonlite')",
          type = "error", duration = 10
        )
        return(NULL)
      }

      req(values$active_dataset)
      obj         <- values$active_dataset
      cluster_col <- input$ai_cluster_col
      species     <- input$ai_species
      tissue      <- input$tissue_context
      top_n       <- input$top_n_genes
      model       <- values$ai_model %||% get_default_model(values$ai_api_key)

      if (is.null(cluster_col) || !(cluster_col %in% colnames(obj@meta.data))) {
        showNotification("Selected cluster column not found.", type = "error")
        return(NULL)
      }

      withProgress(message = "AI Annotation in progress...", value = 0, {

        # Step 0: Prepare
        obj <- .safe_prepare_for_de(obj)

        # Optional re-clustering
        if (isTRUE(input$recluster_for_annotation)) {
          incProgress(0.05, detail = "Re-clustering for annotation...")
          tryCatch({
            res <- input$annotation_resolution

            if (!"pca" %in% names(obj@reductions)) {
              n_hvgs <- if (length(VariableFeatures(obj)) > 0) length(VariableFeatures(obj)) else 2000L
              obj <- FindVariableFeatures(obj, nfeatures = n_hvgs, verbose = FALSE)
              if (DefaultAssay(obj) != "SCT") obj <- ScaleData(obj, verbose = FALSE)
              max_pcs <- max(2L, min(ncol(obj) - 1L, nrow(obj) - 1L, 50L))
              n_pcs <- min(max(10L, as.integer(values$nPCsSelected %||% 30L)), max_pcs)
              obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
            }

            n_dims  <- min(30L, ncol(Embeddings(obj, "pca")))
            k_param <- min(max(20L, floor(sqrt(ncol(obj)) / 2)), ncol(obj) - 1L)

            if (exists("run_standard_graph_clustering", mode = "function")) {
              reclustered <- run_standard_graph_clustering(
                obj, reduction = "pca", dims = n_dims, k_param = k_param,
                resolution = res, algorithm = "leiden",
                cluster_name = "seurat_clusters",
                graph_prefix = "annotation_recluster",
                compute_quality = FALSE, random_seed = 42L
              )
              obj         <- reclustered$object
              cluster_col <- reclustered$context$cluster_name
            } else {
              obj <- FindNeighbors(obj, reduction = "pca", dims = seq_len(n_dims),
                                   k.param = k_param, verbose = FALSE)
              obj <- FindClusters(obj, resolution = res, random.seed = 42L, verbose = FALSE)
              cluster_col <- "seurat_clusters"
            }

            showNotification(
              paste("Re-clustered at resolution", res),
              type = "message"
            )
          }, error = function(e) {
            showNotification(paste("Re-clustering failed:", e$message), type = "warning")
          })
        }

        # Step 1: Find markers
        incProgress(0.1, detail = "Finding top marker genes per cluster...")

        markers_list <- tryCatch(
          find_top_markers_per_cluster(obj, cluster_col, top_n),
          error = function(e) {
            showNotification(paste("Error finding markers:", e$message),
                             type = "error", duration = 10)
            NULL
          }
        )

        if (is.null(markers_list) || length(markers_list) == 0) {
          showNotification("No markers found. Cannot proceed.", type = "error")
          return(NULL)
        }

        # Step 2: Build prompt
        incProgress(0.3, detail = "Building prompt for Claude...")

        tissue_str <- if (nzchar(trimws(tissue))) tissue else "unknown tissue"
        include_marker_ctx <- isTRUE(input$include_marker_context)

        marker_text_parts <- vapply(names(markers_list), function(cl) {
          genes_fc  <- markers_list[[cl]]
          gene_strs <- paste0(names(genes_fc), " (avg_log2FC=", round(genes_fc, 2), ")")
          paste0("Cluster ", cl, ": ", paste(gene_strs, collapse = ", "))
        }, character(1))

        marker_text <- paste(marker_text_parts, collapse = "\n")

        # Marker database context
        marker_db_context <- ""
        if (include_marker_ctx && exists("CELL_TYPE_MARKERS") && length(CELL_TYPE_MARKERS) > 0) {
          all_cluster_genes <- unique(toupper(unlist(lapply(markers_list, names))))
          relevant_types <- character(0)
          for (ct_name in names(CELL_TYPE_MARKERS)) {
            ct_genes <- toupper(names(CELL_TYPE_MARKERS[[ct_name]]))
            overlap  <- length(intersect(all_cluster_genes, ct_genes))
            if (overlap >= 2) relevant_types <- c(relevant_types, ct_name)
          }
          if (length(relevant_types) > 0) {
            relevant_types <- head(relevant_types, 30)
            db_lines <- vapply(relevant_types, function(ct) {
              mk     <- CELL_TYPE_MARKERS[[ct]]
              top_mk <- head(sort(mk, decreasing = TRUE), 8)
              paste0("  ", ct, ": ", paste(names(top_mk), collapse = ", "))
            }, character(1))
            marker_db_context <- paste0(
              "\n\nKnown cell type marker reference (curated database):\n",
              paste(db_lines, collapse = "\n"), "\n"
            )
          }
        }

        prompt <- paste0(
          "You are an expert single-cell biologist specializing in cell type identification.\n\n",
          "## Experiment Context\n",
          "- Species: ", species, "\n",
          "- Tissue/Sample: ", tissue_str, "\n",
          "- Number of clusters: ", length(markers_list), "\n\n",
          "## Top Differentially Expressed Genes Per Cluster\n",
          marker_text, "\n",
          marker_db_context, "\n",
          "## Instructions\n",
          "For each cluster, identify the most likely cell type. Consider:\n",
          "1. Canonical marker gene expression patterns\n",
          "2. The tissue context (expected cell types in ", tissue_str, ")\n",
          "3. Co-expression patterns and marker combinations\n",
          "4. Relative fold changes (higher FC = more reliable)\n\n",
          "Return your response as a JSON array with objects containing:\n",
          "- \"cluster\": the cluster ID (string)\n",
          "- \"cell_type\": predicted cell type (standard nomenclature)\n",
          "- \"confidence\": number 0-1\n",
          "- \"reasoning\": concise explanation citing key markers\n",
          "- \"alternative\": alternative cell type if confidence < 0.8, else null\n",
          "- \"key_markers\": array of top 5-8 supporting markers\n\n",
          "Return ONLY the JSON array. No markdown, no code fences."
        )

        # Save to conversation log
        log <- ai_conv_log()
        log[[length(log) + 1]] <- list(role = "user", content = prompt)
        ai_conv_log(log)

        # Step 3: Call Claude API
        incProgress(0.5, detail = "Calling Claude API...")

        api_result <- if (exists("call_claude_api", mode = "function")) {
          call_claude_api(
            api_key    = api_key,
            model      = model,
            messages   = list(list(role = "user", content = prompt)),
            max_tokens = 8192L,
            timeout    = 120
          )
        } else {
          tryCatch({
            req <- httr2::request("https://api.anthropic.com/v1/messages") |>
              httr2::req_headers(
                "x-api-key" = api_key,
                "anthropic-version" = "2023-06-01",
                "content-type" = "application/json"
              ) |>
              httr2::req_body_json(list(
                model = model, max_tokens = 8192L,
                messages = list(list(role = "user", content = prompt))
              )) |>
              httr2::req_timeout(120)
            resp <- httr2::req_perform(req)
            body <- httr2::resp_body_json(resp)
            list(ok = TRUE, text = body$content[[1]]$text)
          }, error = function(e) {
            list(ok = FALSE, error = e$message)
          })
        }

        api_response <- if (isTRUE(api_result$ok)) {
          api_result$text
        } else {
          showNotification(paste("API error:", api_result$error), type = "error", duration = 15)
          NULL
        }

        if (is.null(api_response)) return(NULL)

        # Save response to conversation log
        log <- ai_conv_log()
        log[[length(log) + 1]] <- list(role = "assistant", content = api_response)
        ai_conv_log(log)

        # Step 4: Parse response
        incProgress(0.8, detail = "Parsing AI response...")

        parse_fn <- if (exists("parse_ai_annotation_response", mode = "function")) {
          parse_ai_annotation_response
        } else {
          parse_ai_response
        }
        parsed <- tryCatch(
          parse_fn(api_response),
          error = function(e) {
            showNotification(paste("Parse error:", e$message), type = "error", duration = 15)
            NULL
          }
        )

        if (is.null(parsed) || nrow(parsed) == 0) {
          showNotification("AI returned no valid predictions.", type = "error")
          return(NULL)
        }

        # Step 5: Cross-validate with marker DB
        if (isTRUE(input$validate_with_markers) && exists("CELL_TYPE_MARKERS") &&
            length(CELL_TYPE_MARKERS) > 0) {
          incProgress(0.85, detail = "Cross-validating with marker database...")

          ai_scoring_obj <- obj
          if (.safe_should_sketch(values, budget = 50000L)) {
            ai_scoring_obj <- .safe_sketch_downsample(obj, max_cells = 50000L, method = "geometric")
          }

          annotator_inst <- EnhancedCellTypeAnnotator$new(CELL_TYPE_MARKERS)
          marker_scores <- tryCatch(
            annotator_inst$score_cell_types_enhanced(
              ai_scoring_obj, assay = .safe_default_assay(ai_scoring_obj)
            ),
            error = function(e) NULL
          )

          if (!is.null(marker_scores)) {
            clusters_vec      <- as.character(ai_scoring_obj@meta.data[[cluster_col]])
            marker_validation <- character(nrow(parsed))
            marker_agreement  <- numeric(nrow(parsed))

            for (i in seq_len(nrow(parsed))) {
              cl       <- parsed$cluster[i]
              cl_cells <- which(clusters_vec == cl)
              if (length(cl_cells) > 0) {
                cl_scores      <- colMeans(marker_scores$scores[cl_cells, , drop = FALSE])
                top_marker_type <- names(which.max(cl_scores))
                marker_validation[i] <- top_marker_type

                ai_pred <- tolower(parsed$cell_type[i])
                mk_pred <- tolower(top_marker_type)
                if (!is.na(ai_pred) && !is.na(mk_pred) &&
                    (ai_pred == mk_pred ||
                     grepl(ai_pred, mk_pred, fixed = TRUE) ||
                     grepl(mk_pred, ai_pred, fixed = TRUE))) {
                  marker_agreement[i] <- 1.0
                } else {
                  lineage_terms <- c("t cell", "b cell", "nk", "monocyte", "macrophage",
                                     "dendritic", "neutrophil", "epithelial", "endothelial",
                                     "fibroblast", "stem", "progenitor", "plasma",
                                     "mast", "erythro", "platelet", "astrocyte",
                                     "oligodendrocyte", "microglia", "neuron")
                  shared <- any(vapply(lineage_terms, function(term) {
                    grepl(term, ai_pred) && grepl(term, mk_pred)
                  }, logical(1)))
                  marker_agreement[i] <- if (shared) 0.5 else 0.0
                }
              } else {
                marker_validation[i] <- NA_character_
                marker_agreement[i]  <- NA_real_
              }
            }

            parsed$marker_db_match     <- marker_validation
            parsed$marker_agreement    <- marker_agreement
            parsed$adjusted_confidence <- pmin(1.0,
              parsed$confidence * 0.7 + marker_agreement * 0.3
            )
          }
        }

        # Step 6: Map predictions to cells
        incProgress(0.9, detail = "Applying AI labels to cells...")

        annot_name <- make.names(input$ai_annotation_name)
        if (!nzchar(annot_name)) annot_name <- "cell_type_ai"
        conf_name <- paste0(annot_name, "_confidence")

        cluster_to_type <- setNames(parsed$cell_type, parsed$cluster)
        cell_clusters   <- as.character(obj@meta.data[[cluster_col]])
        ai_labels       <- unname(cluster_to_type[cell_clusters])
        ai_labels[is.na(ai_labels)] <- "Unknown"
        obj@meta.data[[annot_name]] <- ai_labels

        if ("adjusted_confidence" %in% colnames(parsed)) {
          cluster_to_conf <- setNames(parsed$adjusted_confidence, parsed$cluster)
          obj@meta.data[[conf_name]] <- as.numeric(unname(cluster_to_conf[cell_clusters]))
        }

        values$active_dataset <- obj
        last_annot_col(annot_name)
        ai_results(parsed)

        incProgress(1, detail = "Done!")
      })

      showNotification(
        paste0("AI annotation complete! Column: '", last_annot_col(), "'"),
        type = "message"
      )
    })

    # ---------------------------------------------------------------
    # AI DimPlot
    # ---------------------------------------------------------------
    output$ai_dimplot <- renderPlotly({
      req(values$active_dataset)
      acol <- last_annot_col()
      req(acol %in% colnames(values$active_dataset@meta.data))

      obj <- values$active_dataset
      red <- .best_reduction(obj)
      if (is.null(red)) return(plotly_empty())

      emb  <- Embeddings(obj, red)[, 1:2, drop = FALSE]
      labs <- obj@meta.data[[acol]]

      p <- plot_ly(
        x = emb[, 1], y = emb[, 2],
        type = "scatter", mode = "markers",
        marker = list(size = 3, opacity = 0.9),
        color = labs, colors = .get_palette(length(unique(labs))),
        text = paste0("Cell: ", colnames(obj), "\nType: ", labs),
        hovertemplate = paste0(
          "%{text}<br>",
          toupper(red), "1: %{x:.2f}<br>",
          toupper(red), "2: %{y:.2f}<extra></extra>"
        )
      ) %>%
        layout(
          title = paste(toupper(red), "colored by", acol),
          xaxis = list(title = paste0(toupper(red), "1")),
          yaxis = list(title = paste0(toupper(red), "2"))
        )
      if (exists("capture_umap_snapshot", mode = "function"))
        capture_umap_snapshot(values, p, "Annotation (AI)", coloring = acol, method = red)
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Annotation",
                       figure_type = "AI annotation dimplot",
                       label = "AI Annotation DimPlot", source_module = "Annotation")
      }
      p
    })

    # 3D UMAP for AI
    output$ai_3d_umap <- renderPlotly({
      req(values$active_dataset)
      obj <- values$active_dataset
      acol <- last_annot_col()
      if (!"umap3d" %in% names(obj@reductions) || !acol %in% colnames(obj@meta.data)) {
        return(plotly_empty() %>% layout(title = "3D UMAP not available. Generate it first."))
      }
      emb <- Embeddings(obj, "umap3d")
      if (ncol(emb) < 3) return(plotly_empty())

      labs <- obj@meta.data[[acol]]
      n_cells <- nrow(emb)
      idx <- if (n_cells > 50000 && exists("adaptive_plot_index", mode = "function")) {
        adaptive_plot_index(n_cells, max_points = 50000)
      } else {
        seq_len(n_cells)
      }

      plot_ly(
        x = emb[idx, 1], y = emb[idx, 2], z = emb[idx, 3],
        type = "scatter3d", mode = "markers",
        marker = list(size = 2, opacity = 0.8),
        color = labs[idx], colors = .get_palette(length(unique(labs))),
        text = colnames(obj)[idx],
        hovertemplate = paste0(
          "Cell: %{text}<br>UMAP1: %{x:.2f}<br>UMAP2: %{y:.2f}<br>UMAP3: %{z:.2f}<br>",
          "Type: ", labs[idx], "<extra></extra>"
        )
      ) %>%
        layout(
          title = paste("3D UMAP -", acol),
          scene = list(
            xaxis = list(title = "UMAP1"),
            yaxis = list(title = "UMAP2"),
            zaxis = list(title = "UMAP3")
          )
        )
    })

    # ---------------------------------------------------------------
    # Results table
    # ---------------------------------------------------------------
    output$ai_results_table <- DT::renderDT({
      req(ai_results())
      res <- ai_results()
      display_cols <- c("cluster", "cell_type", "confidence")
      if ("adjusted_confidence" %in% colnames(res)) display_cols <- c(display_cols, "adjusted_confidence")
      if ("marker_db_match" %in% colnames(res))     display_cols <- c(display_cols, "marker_db_match")
      if ("marker_agreement" %in% colnames(res))    display_cols <- c(display_cols, "marker_agreement")
      display_cols <- c(display_cols, "alternative")
      if ("key_markers" %in% colnames(res))          display_cols <- c(display_cols, "key_markers")
      display_cols <- display_cols[display_cols %in% colnames(res)]

      df <- res[, display_cols, drop = FALSE]
      col_names <- c("Cluster", "AI Prediction", "AI Confidence")
      if ("adjusted_confidence" %in% display_cols) col_names <- c(col_names, "Adj. Confidence")
      if ("marker_db_match" %in% display_cols)     col_names <- c(col_names, "Marker DB Match")
      if ("marker_agreement" %in% display_cols)    col_names <- c(col_names, "Agreement")
      col_names <- c(col_names, "Alternative")
      if ("key_markers" %in% display_cols)          col_names <- c(col_names, "Key Markers")
      colnames(df) <- col_names

      dt <- DT::datatable(df,
        options = list(pageLength = 25, scrollX = TRUE, dom = "tip"),
        rownames = FALSE
      ) %>% DT::formatRound("AI Confidence", digits = 2)

      if ("Adj. Confidence" %in% col_names) {
        dt <- dt %>% DT::formatRound("Adj. Confidence", digits = 2)
      }
      if ("Agreement" %in% col_names) {
        dt <- dt %>% DT::formatStyle("Agreement",
          backgroundColor = DT::styleInterval(c(0.25, 0.75), c("#f8d7da", "#fff3cd", "#d4edda"))
        )
      }
      dt <- dt %>% DT::formatStyle("AI Confidence",
        background = DT::styleColorBar(c(0, 1), "lightgreen"),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = df,
                            description = "AI Annotation Results",
                            source_module = "Annotation")
      }
      dt
    })

    # Reasoning table
    output$ai_reasoning_table <- DT::renderDT({
      req(ai_results())
      res  <- ai_results()
      keep <- intersect(c("cluster", "cell_type", "confidence", "reasoning", "alternative"),
                         colnames(res))
      df <- res[, keep, drop = FALSE]
      nice_names <- c(cluster = "Cluster", cell_type = "Cell Type", confidence = "Confidence",
                      reasoning = "Reasoning", alternative = "Alternative")
      colnames(df) <- nice_names[colnames(df)]
      DT::datatable(df,
        options = list(pageLength = 20, scrollX = TRUE, dom = "tip",
          columnDefs = list(list(width = "40%",
                                 targets = which(colnames(df) == "Reasoning") - 1L))),
        rownames = FALSE
      )
    })

    # Comparison plot
    output$ai_comparison_plot <- renderPlotly({
      req(values$active_dataset)
      obj <- values$active_dataset
      md  <- obj@meta.data

      compare_col <- NULL
      for (candidate in c("cell_type_consensus", "cell_type_sctype", "cell_type_panglao", "cell_type_clustifyr", "cell_type_enhanced")) {
        if (candidate %in% colnames(md) && !all(is.na(md[[candidate]]))) {
          compare_col <- candidate
          break
        }
      }

      acol <- last_annot_col()
      if (is.null(compare_col) || !(acol %in% colnames(md))) {
        return(plot_ly() %>%
          layout(title = "No ScType annotation available for comparison",
                 xaxis = list(visible = FALSE), yaxis = list(visible = FALSE)))
      }

      red <- .best_reduction(obj)
      if (is.null(red)) return(plotly_empty())
      emb <- Embeddings(obj, red)[, 1:2, drop = FALSE]

      ai_labs     <- md[[acol]]
      marker_labs <- md[[compare_col]]

      p <- plot_ly(
        x = emb[, 1], y = emb[, 2],
        type = "scatter", mode = "markers",
        marker = list(size = 2, opacity = 0.7),
        color = ai_labs,
        colors = .get_palette(length(unique(ai_labs))),
        text = paste0("AI: ", ai_labs, "\n", compare_col, ": ", marker_labs),
        hovertemplate = "%{text}<extra></extra>",
        showlegend = TRUE
      ) %>%
        layout(
          title = paste(acol, "vs", compare_col),
          xaxis = list(title = paste0(toupper(red), "1")),
          yaxis = list(title = paste0(toupper(red), "2"))
        )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Annotation",
                       figure_type = "AI vs marker-based comparison",
                       label = "AI Comparison Plot", source_module = "Annotation")
      }
      p
    })

    # Comparison table
    output$ai_comparison_table <- DT::renderDT({
      req(values$active_dataset)
      obj <- values$active_dataset
      md  <- obj@meta.data

      compare_col <- NULL
      for (candidate in c("cell_type_consensus", "cell_type_sctype", "cell_type_panglao", "cell_type_clustifyr", "cell_type_enhanced")) {
        if (candidate %in% colnames(md) && !all(is.na(md[[candidate]]))) {
          compare_col <- candidate
          break
        }
      }

      acol <- last_annot_col()
      if (is.null(compare_col) || !(acol %in% colnames(md))) {
        return(DT::datatable(
          data.frame(Message = "No ScType annotation available for comparison."),
          rownames = FALSE
        ))
      }

      cluster_col <- input$ai_cluster_col
      if (is.null(cluster_col) || !(cluster_col %in% colnames(md))) {
        cluster_col <- "seurat_clusters"
      }

      clusters <- sort(unique(as.character(md[[cluster_col]])))
      comparison <- data.frame(
        Cluster    = clusters,
        AI_Type    = NA_character_,
        Marker_Type = NA_character_,
        Agreement  = NA_character_,
        stringsAsFactors = FALSE
      )

      for (i in seq_along(clusters)) {
        cl   <- clusters[i]
        mask <- as.character(md[[cluster_col]]) == cl
        ai_tab <- sort(table(md[[acol]][mask]), decreasing = TRUE)
        mk_tab <- sort(table(md[[compare_col]][mask]), decreasing = TRUE)
        ai_top <- if (length(ai_tab) > 0) names(ai_tab)[1] else "Unknown"
        mk_top <- if (length(mk_tab) > 0) names(mk_tab)[1] else "Unknown"
        comparison$AI_Type[i]    <- ai_top
        comparison$Marker_Type[i] <- mk_top
        comparison$Agreement[i]  <- if (tolower(ai_top) == tolower(mk_top)) "Yes" else "No"
      }

      colnames(comparison) <- c("Cluster", "AI Prediction",
                                 paste0("Marker (", compare_col, ")"), "Agree?")

      DT::datatable(comparison,
        options = list(pageLength = 20, scrollX = TRUE, dom = "tip"),
        rownames = FALSE
      ) %>%
        DT::formatStyle("Agree?",
          backgroundColor = DT::styleEqual(c("Yes", "No"), c("#d4edda", "#f8d7da"))
        )
    })

    # Marker heatmap
    output$ai_marker_heatmap <- renderPlotly({
      req(ai_results(), values$active_dataset)
      res  <- ai_results()
      acol <- last_annot_col()
      req(acol %in% colnames(values$active_dataset@meta.data))

      obj <- values$active_dataset

      all_markers <- character(0)
      if ("key_markers" %in% colnames(res)) {
        for (i in seq_len(nrow(res))) {
          mk <- res$key_markers[i]
          if (is.na(mk) || !nzchar(mk)) next
          genes <- trimws(unlist(strsplit(as.character(mk), ",")))
          genes <- genes[nzchar(genes)]
          all_markers <- c(all_markers, genes)
        }
      }
      all_markers <- unique(all_markers)

      if (length(all_markers) < 2) {
        return(plot_ly() %>% layout(title = "Not enough markers for heatmap"))
      }

      avail <- intersect(all_markers, rownames(obj))
      if (length(avail) < 2) {
        return(plot_ly() %>% layout(title = "Marker genes not found in expression data"))
      }
      avail <- head(avail, 50)

      expr <- tryCatch({
        mat <- .safe_get_assay_data(obj, "data", .safe_default_assay(obj))
        as.matrix(mat[avail, , drop = FALSE])
      }, error = function(e) NULL)
      if (is.null(expr)) {
        return(plot_ly() %>% layout(title = "Could not retrieve expression data"))
      }

      cell_types <- as.character(obj@meta.data[[acol]])
      unique_ct  <- sort(unique(cell_types))
      mean_mat   <- matrix(0, nrow = length(avail), ncol = length(unique_ct),
                            dimnames = list(avail, unique_ct))
      for (ct in unique_ct) {
        idx <- which(cell_types == ct)
        if (length(idx) > 1) {
          mean_mat[, ct] <- rowMeans(expr[, idx, drop = FALSE])
        } else if (length(idx) == 1) {
          mean_mat[, ct] <- expr[, idx]
        }
      }

      scaled <- t(apply(mean_mat, 1, function(x) {
        if (sd(x) == 0) return(rep(0, length(x)))
        (x - mean(x)) / sd(x)
      }))

      p <- plot_ly(
        x = colnames(scaled), y = rownames(scaled),
        z = scaled, type = "heatmap",
        colorscale = list(c(0, "#313695"), c(0.5, "#ffffbf"), c(1, "#a50026")),
        hovertemplate = "Gene: %{y}<br>Cell Type: %{x}<br>Z-score: %{z:.2f}<extra></extra>"
      ) %>%
        layout(
          title = "Marker Gene Expression by Cell Type (AI)",
          xaxis = list(title = "", tickangle = -45),
          yaxis = list(title = "", autorange = "reversed"),
          margin = list(l = 100, b = 120)
        )
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Annotation",
                       figure_type = "AI marker heatmap",
                       label = "AI Marker Heatmap", source_module = "Annotation")
      }
      p
    })

    # ---------------------------------------------------------------
    # Downloads
    # ---------------------------------------------------------------
    output$dl_ai_annotations <- downloadHandler(
      filename = function() paste0("ai_annotations_", Sys.Date(), ".csv"),
      content = function(file) {
        req(ai_results())
        results_df <- ai_results()
        acol <- last_annot_col()
        if (!is.null(values$active_dataset) && acol %in% colnames(values$active_dataset@meta.data)) {
          cell_df <- data.frame(
            cell_barcode = colnames(values$active_dataset),
            annotation   = values$active_dataset@meta.data[[acol]],
            stringsAsFactors = FALSE
          )
          colnames(cell_df)[2] <- acol
          write.csv(results_df, file, row.names = FALSE)
          summary_lines <- readLines(file)
          writeLines(c(
            "# === Cluster-Level Summary ===",
            summary_lines,
            "",
            "# === Cell-Level Annotations ===",
            utils::capture.output(write.csv(cell_df, row.names = FALSE))
          ), file)
        } else {
          write.csv(results_df, file, row.names = FALSE)
        }
      }
    )

    output$dl_ai_reasoning_table <- downloadHandler(
      filename = function() paste0("ai_reasoning_", Sys.Date(), ".csv"),
      content = function(file) {
        req(ai_results())
        res  <- ai_results()
        keep <- intersect(c("cluster", "cell_type", "confidence", "reasoning", "alternative"),
                           colnames(res))
        write.csv(res[, keep, drop = FALSE], file, row.names = FALSE)
      }
    )

    output$dl_ai_comparison_table <- downloadHandler(
      filename = function() paste0("ai_comparison_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          obj <- values$active_dataset
          md  <- obj@meta.data
          compare_col <- NULL
          for (candidate in c("cell_type_consensus", "cell_type_sctype", "cell_type_panglao", "cell_type_clustifyr", "cell_type_enhanced")) {
            if (candidate %in% colnames(md) && !all(is.na(md[[candidate]]))) {
              compare_col <- candidate
              break
            }
          }
          acol <- last_annot_col()
          if (is.null(compare_col) || !(acol %in% colnames(md))) {
            write.csv(data.frame(Message = "No comparison available"), file, row.names = FALSE)
          } else {
            cluster_col <- if (!is.null(input$ai_cluster_col) &&
                               isTRUE(input$ai_cluster_col %in% colnames(md)))
                             input$ai_cluster_col
                           else "seurat_clusters"
            clusters <- sort(unique(as.character(md[[cluster_col]])))
            comp <- data.frame(Cluster = clusters, AI_Type = NA_character_,
                               Marker_Type = NA_character_, Agreement = NA_character_,
                               stringsAsFactors = FALSE)
            for (i in seq_along(clusters)) {
              mask <- as.character(md[[cluster_col]]) == clusters[i]
              ai_tab <- sort(table(md[[acol]][mask]), decreasing = TRUE)
              mk_tab <- sort(table(md[[compare_col]][mask]), decreasing = TRUE)
              comp$AI_Type[i]     <- if (length(ai_tab) > 0) names(ai_tab)[1] else "Unknown"
              comp$Marker_Type[i] <- if (length(mk_tab) > 0) names(mk_tab)[1] else "Unknown"
              comp$Agreement[i]   <- if (tolower(comp$AI_Type[i]) == tolower(comp$Marker_Type[i]))
                                       "Yes" else "No"
            }
            write.csv(comp, file, row.names = FALSE)
          }
        }, error = function(e) write.csv(data.frame(Error = e$message), file, row.names = FALSE))
      }
    )

    # Apply AI labels as consensus
    observeEvent(input$apply_ai_labels, {
      req(values$active_dataset)
      acol <- last_annot_col()
      req(acol %in% colnames(values$active_dataset@meta.data))

      obj <- values$active_dataset
      obj@meta.data[["cell_type_consensus"]] <- obj@meta.data[[acol]]
      values$active_dataset <- obj

      showNotification(
        paste0("'", acol, "' applied to 'cell_type_consensus' column."),
        type = "message"
      )
    })

  })
}

# =============================================================================
# END OF MODULE
# =============================================================================
