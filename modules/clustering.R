# ============================================================================
# Enhanced Graph-Based Clustering Module
# modules/clustering.R
# Comprehensive clustering engine for scRNA-seq analysis
# Supports: Leiden, Louvain, SLM, Multilevel algorithms
# Features: Single resolution, Resolution sweep, Optimal resolution,
#           Sub-clustering, Quality metrics, Clustree, 3D visualization,
#           Split-by-batch, Highlight UMAP, AI integration
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(plotly)
  library(dplyr)
  library(rlang)
  library(igraph)
  library(ggplot2)
  library(viridis)
  library(cluster)
  library(Matrix)
})

# get_color_palette(), best_reduction(), safe_default_assay(),
# detect_cell_type_col() are defined in core_robustness.R and sourced before
# this module. We guard every external dependency with exists() so the module
# can still load standalone during development or testing.


# =============================================================================
# 1. compute_silhouette_safe
#    Computes the average silhouette width for a set of cluster labels.
#    Downsamples when the dataset is large to keep runtime bounded (the
#    distance matrix is O(n^2) in memory). Returns NA_real_ on any error
#    so callers never need to handle exceptions.
# =============================================================================
compute_silhouette_safe <- function(embeddings, labels, max_cells = 5000L, seed = 42L) {
  # Validate inputs early

  if (is.null(embeddings) || is.null(labels)) return(NA_real_)
  if (!is.matrix(embeddings) && !is.data.frame(embeddings)) {
    embeddings <- tryCatch(as.matrix(embeddings), error = function(e) NULL)
    if (is.null(embeddings)) return(NA_real_)
  }

  n <- nrow(embeddings)
  if (n < 3L) return(NA_real_)

  labels <- as.character(labels)
  if (length(labels) != n) return(NA_real_)

  # Silhouette requires at least two distinct clusters
  unique_labels <- unique(labels)
  if (length(unique_labels) < 2L) return(NA_real_)

  # Downsample preserving proportional representation when dataset is large.
  # The distance matrix for silhouette is O(n^2) so we cap at max_cells.
  max_cells <- max(10L, as.integer(max_cells))
  if (n > max_cells) {
    set.seed(seed)
    idx <- sample_cells_stratified(labels, max_cells = max_cells, seed = seed)
    embeddings <- embeddings[idx, , drop = FALSE]
    labels <- labels[idx]
  }

  # Convert labels to integer factor codes (cluster::silhouette requires integer)
  label_int <- as.integer(as.factor(labels))

  sil <- tryCatch({
    d <- dist(embeddings)
    cluster::silhouette(label_int, d)
  }, error = function(e) {
    # Common failure: all cells in one cluster after downsampling, or
    # dist() runs out of memory on very wide embeddings.
    NULL
  })

  if (is.null(sil)) return(NA_real_)

  avg <- mean(sil[, 3], na.rm = TRUE)
  if (is.finite(avg)) avg else NA_real_
}


# =============================================================================
# 2. validate_clustering_contract
#    Checks that the Seurat object meets minimum requirements before clustering.
#    Delegates to the platform-wide validate_seurat_contract() when available
#    (from core_robustness.R). When the validator is absent, performs lightweight
#    sanity checks directly. Returns TRUE if the object passes, FALSE otherwise,
#    showing user-facing notifications for failures and warnings.
# =============================================================================
validate_clustering_contract <- function(values, obj, require_reductions = NULL) {
  # If the platform contract validator exists, use it for a thorough check
  if (exists("validate_seurat_contract", mode = "function")) {
    contract <- validate_seurat_contract(
      obj,
      stage           = "clustering",
      require_reductions = require_reductions,
      min_cells       = 20,
      min_features    = 50
    )
    if (!isTRUE(contract$ok)) {
      msg <- if (exists("format_contract_feedback", mode = "function")) {
        format_contract_feedback(contract)
      } else {
        paste(contract$errors, collapse = " | ")
      }
      showNotification(msg, type = "error", duration = 12)
      # Audit log the failure if the audit system is available
      if (exists("append_audit_log", mode = "function")) {
        append_audit_log(
          values,
          module  = "clustering",
          action  = "contract_validation_failed",
          level   = "error",
          details = list(
            errors      = contract$errors,
            warnings    = contract$warnings,
            suggestions = contract$suggestions
          )
        )
      }
      return(FALSE)
    }
    # Surface warnings even when the contract passes
    if (length(contract$warnings) > 0) {
      showNotification(
        paste("Contract warnings:", paste(contract$warnings, collapse = " | ")),
        type = "warning", duration = 8
      )
    }
    return(TRUE)
  }

  # Fallback lightweight validation when core_robustness.R is not loaded
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    showNotification("No valid Seurat object found.", type = "error", duration = 10)
    return(FALSE)
  }
  if (ncol(obj) < 20) {
    showNotification(
      paste0("Dataset has only ", ncol(obj), " cells. At least 20 are needed for clustering."),
      type = "error", duration = 10
    )
    return(FALSE)
  }
  if (!is.null(require_reductions)) {
    avail <- names(obj@reductions)
    missing <- setdiff(require_reductions, avail)
    if (length(missing) > 0) {
      showNotification(
        paste0("Missing required reductions: ", paste(missing, collapse = ", ")),
        type = "error", duration = 10
      )
      return(FALSE)
    }
  }
  TRUE
}


# =============================================================================
# 3. begin_clustering_job
#    Registers a new job with the platform's job tracker. The job tracker
#    (start_job from core_robustness.R) is optional — when absent we return
#    NULL silently. The caller stores the returned job_id and passes it to
#    end_clustering_job() when the operation finishes.
# =============================================================================
begin_clustering_job <- function(values, action, meta = list()) {
  if (!exists("start_job", mode = "function")) return(NULL)
  tryCatch(
    start_job(values, module = "clustering", action = action, meta = meta),
    error = function(e) {
      # Job tracking is non-critical; never let it block clustering
      message("[Clustering] Job tracker start_job failed: ", conditionMessage(e))
      NULL
    }
  )
}


# =============================================================================
# 4. end_clustering_job
#    Marks a previously started job as finished. Gracefully no-ops when
#    job_id is NULL (job tracking was unavailable) or when finish_job is
#    not defined.
# =============================================================================
end_clustering_job <- function(values, job_id, status = "done", message = NULL, details = list()) {
  if (is.null(job_id) || !exists("finish_job", mode = "function")) return(invisible(NULL))
  tryCatch(
    finish_job(values, job_id = job_id, status = status, message = message, details = details),
    error = function(e) {
      message("[Clustering] Job tracker finish_job failed: ", conditionMessage(e))
      invisible(NULL)
    }
  )
}


# =============================================================================
# 5. checkpoint_clustering
#    Registers a dataset checkpoint in the provenance system so users can
#    revert to this state later. Wraps register_dataset_checkpoint() from
#    core_robustness.R; silently no-ops when the function is absent.
# =============================================================================
checkpoint_clustering <- function(values, action, obj, params = list(), method = NULL) {
  if (!exists("register_dataset_checkpoint", mode = "function")) return(invisible(NULL))
  tryCatch(
    register_dataset_checkpoint(
      values,
      module  = "clustering",
      action  = action,
      obj     = obj,
      params  = params,
      methods = list(method = method)
    ),
    error = function(e) {
      message("[Clustering] Checkpoint registration failed: ", conditionMessage(e))
      invisible(NULL)
    }
  )
}


# =============================================================================
# 6. default_clustering_quality
#    Returns a fresh quality-metrics structure with all fields initialized to
#    their empty/NA defaults. Used to reset quality state before a new
#    clustering run and to ensure consistent downstream access patterns.
# =============================================================================
default_clustering_quality <- function() {
  list(
    silhouette_scores         = NULL,
    silhouette_cluster_labels = NULL,
    silhouette_avg            = NA_real_,
    modularity                = NA_real_,
    calinski_score            = NA_real_
  )
}


# =============================================================================
# 7. notify_clustering
#    Thin wrapper around Shiny's showNotification. Exists so helper functions
#    that live outside the server context can still emit UI notifications
#    without having to check availability every time. Falls back to message()
#    when running outside a Shiny session (e.g., during testing).
# =============================================================================
notify_clustering <- function(message, type = "warning", duration = 8) {
  if (exists("showNotification", mode = "function")) {
    tryCatch(
      showNotification(message, type = type, duration = duration),
      error = function(e) {
        # Outside a Shiny session showNotification will error; fall back to console
        base::message("[Clustering] ", type, ": ", message)
      }
    )
  } else {
    base::message("[Clustering] ", type, ": ", message)
  }
  invisible(message)
}


# =============================================================================
# 8. leiden_backend_available
#    Pure R check for whether igraph::cluster_leiden is available. Since
#    igraph >= 2.0, cluster_leiden is included and works without any Python
#    dependency. We verify that the function is actually callable rather than
#    just checking package version, because some minimal igraph installs may
#    omit it.
# =============================================================================
leiden_backend_available <- function() {
  tryCatch({
    if (!requireNamespace("igraph", quietly = TRUE)) return(FALSE)
    ver <- utils::packageVersion("igraph")
    if (ver < "2.0.0") {
      message("igraph < 2.0 detected — Leiden unavailable, using Louvain fallback")
      return(FALSE)
    }
    is.function(igraph::cluster_leiden)
  }, error = function(e) FALSE)
}


# =============================================================================
# 9. sanitize_cluster_column_name
#    Cleans a user-supplied column name so it is safe for use in Seurat
#    meta.data (which is just a data.frame). Trims whitespace, replaces
#    non-alphanumeric characters with underscores, and falls back to a
#    default when the result is empty.
# =============================================================================
sanitize_cluster_column_name <- function(name, fallback = "seurat_clusters") {
  out <- as.character(name %||% fallback)
  out <- trimws(out)
  if (!nzchar(out)) out <- fallback
  # Replace any run of non-alphanumeric / non-underscore characters

  out <- gsub("[^A-Za-z0-9_]+", "_", out)
  # Strip leading/trailing underscores that may result from replacement
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) out <- fallback
  out
}


# =============================================================================
# 10. format_resolution_id
#     Formats a numeric resolution value into a clean string suitable for use
#     in column names (e.g., 0.800 -> "0.8", 1.000 -> "1"). Trailing zeros
#     after the decimal point are stripped for readability while preserving
#     enough precision to distinguish close resolutions.
# =============================================================================
format_resolution_id <- function(resolution, digits = 3L) {
  res_num <- suppressWarnings(as.numeric(resolution))
  if (is.na(res_num)) return("0.8")
  out <- formatC(res_num, format = "f", digits = as.integer(digits))
  # Strip trailing zeros: "0.800" -> "0.8", "1.000" -> "1"
  out <- sub("0+$", "", out)
  out <- sub("\\.$", "", out)
  out
}


# =============================================================================
# 11. sample_cells_stratified
#     Stratified downsampling that preserves proportional cluster representation.
#     Given a label vector, returns a set of indices (sorted) no larger than
#     max_cells. Each group receives a base quota of floor(max_cells / n_groups)
#     cells; remaining slots are filled randomly from leftover cells. This
#     ensures tiny clusters are always represented in the sample.
# =============================================================================
sample_cells_stratified <- function(labels, max_cells = 3000L, seed = 42L) {
  labels <- as.character(labels)
  n <- length(labels)
  max_cells <- max(2L, as.integer(max_cells))
  if (n <= max_cells) return(seq_len(n))

  groups <- split(seq_len(n), labels)
  n_groups <- max(1L, length(groups))
  # Each group gets at least 2 cells (if available) so we never lose a cluster
  base_quota <- max(2L, floor(max_cells / n_groups))

  set.seed(seed)
  selected <- unlist(lapply(groups, function(idx) {
    keep <- min(length(idx), base_quota)
    if (length(idx) <= keep) idx else sample(idx, keep)
  }), use.names = FALSE)

  # Fill any remaining capacity with randomly drawn leftover cells
  remaining <- setdiff(seq_len(n), selected)
  slots_left <- max_cells - length(selected)
  if (slots_left > 0L && length(remaining) > 0L) {
    selected <- c(selected, sample(remaining, min(slots_left, length(remaining))))
  }

  sort(unique(selected))
}


# =============================================================================
# 12. select_clustering_reduction
#     Selects the best high-dimensional reduction for graph construction.
#     Preference order: PCA > Harmony > integrated > LSI > ICA (all >2D).
#     If allow_2d is FALSE (default), reductions with only 2 columns (like
#     UMAP/tSNE) are excluded because they are poor inputs for SNN graph
#     construction — they discard too much variance.
#
#     Returns the name of the chosen reduction, or NULL if none qualifies.
# =============================================================================
select_clustering_reduction <- function(seurat_obj, requested = NULL, allow_2d = FALSE) {
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) return(NULL)

  reductions <- names(seurat_obj@reductions)
  if (length(reductions) == 0) return(NULL)

  # Canonical preference order for clustering (high-dimensional first)
  high_dim_pref <- c("pca", "harmony", "integrated", "integrated.dr", "lsi", "ica")

  # Helper: checks that a reduction exists and has sufficient dimensionality

  is_valid_reduction <- function(red) {
    if (!red %in% reductions) return(FALSE)
    n_dims <- tryCatch(ncol(Embeddings(seurat_obj, red)), error = function(e) 0L)
    isTRUE(allow_2d) || n_dims > 2L
  }

  # Honour an explicit user request first
  requested <- requested %||% NULL
  if (!is.null(requested) && nzchar(requested) && is_valid_reduction(requested)) {
    return(requested)
  }

  # Walk the preference list
  for (red in high_dim_pref) {
    if (is_valid_reduction(red)) return(red)
  }

  # Last resort: any reduction that passes the dimensionality check
  candidates <- reductions[vapply(reductions, is_valid_reduction, logical(1))]
  if (length(candidates) > 0) candidates[1] else NULL
}


# =============================================================================
# 13. resolve_graph_clustering_context
#     Central validation & normalization of all clustering parameters. Produces
#     a fully resolved "context" list that downstream code can use without
#     further checking. Handles:
#       - Reduction selection (delegates to select_clustering_reduction)
#       - Dimension count capping to available dimensions
#       - k-neighbor clamping to [5, n_cells-1]
#       - Resolution validation (must be > 0)
#       - Algorithm mapping to Seurat's integer codes:
#           leiden = 4, louvain = 1, slm = 3, multilevel = 2
#       - Automatic Leiden -> Louvain fallback if igraph::cluster_leiden
#         is unavailable
#       - NN method selection: rann for <= 5000 cells, annoy for larger
#       - Graph naming (nn_name, graph_name) based on graph_prefix
# =============================================================================
resolve_graph_clustering_context <- function(seurat_obj,
                                             reduction    = NULL,
                                             dims         = 30L,
                                             k_param      = 20L,
                                             resolution   = 0.8,
                                             algorithm    = "leiden",
                                             cluster_name = "seurat_clusters",
                                             prune_snn    = 1/15,
                                             graph_type   = "snn",
                                             graph_prefix = "clustering",
                                             random_seed  = 42L) {
  # --- Object validation ---
  if (is.null(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("A valid Seurat object is required for clustering.")
  }
  n_cells <- ncol(seurat_obj)
  if (n_cells < 3L) {
    stop("Need at least 3 cells to run graph-based clustering.")
  }

  # --- Graph type ---
  graph_type_key <- tolower(as.character(graph_type %||% "snn"))
  if (graph_type_key == "knn") graph_type_key <- "dense_snn"
  if (!graph_type_key %in% c("snn", "dense_snn")) graph_type_key <- "snn"

  # --- Reduction ---
  reduction_use <- select_clustering_reduction(
    seurat_obj, requested = reduction, allow_2d = FALSE
  )
  if (is.null(reduction_use)) {
    stop("No suitable high-dimensional reduction found. Run PCA or Harmony first.")
  }

  # --- Dimensions ---
  available_dims <- ncol(Embeddings(seurat_obj, reduction_use))
  dims_use <- suppressWarnings(as.integer(dims %||% min(30L, available_dims)))
  if (is.na(dims_use) || dims_use < 2L) dims_use <- min(30L, available_dims)
  dims_use <- min(max(2L, dims_use), available_dims)

  # --- k-nearest neighbors ---
  k_use <- suppressWarnings(as.integer(k_param %||% 20L))
  if (is.na(k_use) || k_use < 2L) k_use <- 20L
  k_use <- min(max(5L, k_use), n_cells - 1L)

  # Warn about unusually high k values that tend to merge communities
  if (k_use > 50L) {
    notify_clustering(
      paste0("k-neighbors=", k_use, " is unusually high (standard: 20). ",
             "This may produce a single giant cluster. Consider lowering to 20-30."),
      type = "warning", duration = 10
    )
  }

  # --- Resolution ---
  resolution_use <- suppressWarnings(as.numeric(resolution %||% 0.8))
  if (is.na(resolution_use) || resolution_use <= 0) resolution_use <- 0.8

  # --- SNN pruning ---
  prune_use <- suppressWarnings(as.numeric(prune_snn %||% (1 / 15)))
  if (is.na(prune_use)) prune_use <- 1 / 15
  # Dense SNN: no pruning (prune = 0). Standard SNN: clamp to [0, 1].
  prune_use <- if (graph_type_key == "dense_snn") 0 else min(max(prune_use, 0), 1)

  # --- Random seed ---
  seed_use <- suppressWarnings(as.integer(random_seed %||% 42L))
  if (is.na(seed_use)) seed_use <- 42L

  # --- Algorithm ---
  algo_key <- tolower(as.character(algorithm %||% "louvain"))
  if (!algo_key %in% c("leiden", "louvain", "slm", "multilevel")) {
    algo_key <- "louvain"
  }
  # Leiden requires igraph >= 2.0 with cluster_leiden
  if (algo_key == "leiden" && !leiden_backend_available()) {
    stop("Leiden algorithm not available: igraph::cluster_leiden missing. ",
         "Update igraph: install.packages('igraph') or install the leiden package: install.packages('leiden')")
  }

  # --- Graph naming ---
  graph_prefix <- sanitize_cluster_column_name(graph_prefix, fallback = "clustering")
  nn_name    <- paste0(graph_prefix, "_nn")
  graph_name <- paste0(graph_prefix, "_snn")

  # --- NN method selection ---
  # rann is exact and fast for small datasets; annoy uses approximate
  # nearest neighbors with random projection trees and scales better.
  nn_method <- if (n_cells <= 5000L) "rann" else "annoy"
  # More trees = more accurate but slower; use fewer for very large datasets
  n_trees   <- if (n_cells >= 100000L) 25L else 50L

  list(
    reduction      = reduction_use,
    dims           = dims_use,
    available_dims = available_dims,
    k_param        = k_use,
    resolution     = resolution_use,
    cluster_name   = sanitize_cluster_column_name(cluster_name, fallback = "seurat_clusters"),
    prune_snn      = prune_use,
    graph_type     = graph_type_key,
    graph_prefix   = graph_prefix,
    nn_name        = nn_name,
    graph_name     = graph_name,
    random_seed    = seed_use,
    algorithm_key  = algo_key,
    algorithm_code = switch(algo_key,
      "leiden"     = 4L,
      "louvain"    = 1L,
      "slm"        = 3L,
      "multilevel" = 2L,
      1L  # default safety net
    ),
    nn_method      = nn_method,
    n_trees        = n_trees
  )
}


# =============================================================================
# 14. find_clusters_safe
#     CRITICAL wrapper around Seurat::FindClusters.
#
#     For algorithm == 4 (Leiden):
#       Sets leiden_method = "igraph" and leiden_objective_function = "modularity"
#       to guarantee the pure-R igraph backend is used. This avoids any Python
#       leidenalg dependency. Also enforces random.seed >= 1 as required by
#       Seurat's internal RunLeiden.
#
#     For all other algorithms (Louvain=1, Multilevel=2, SLM=3):
#       Delegates directly to FindClusters with no extra parameters needed —
#       these all use igraph's C-level implementations and are always available.
#
#     Returns the Seurat object with cluster assignments in the active ident
#     and meta.data. Errors propagate to the caller (no internal tryCatch)
#     because the caller needs to handle them in context.
# =============================================================================
find_clusters_safe <- function(object, graph.name, resolution, algorithm,
                               random.seed = 42L, verbose = FALSE) {
  algo <- as.integer(algorithm)

  if (algo == 4L) {
    # Leiden via igraph — pure R, no Python dependency
    message("[Clustering] FindClusters Leiden via igraph (pure R)")
    return(FindClusters(
      object,
      graph.name                = graph.name,
      resolution                = resolution,
      algorithm                 = 4L,
      leiden_method             = "igraph",
      leiden_objective_function  = "modularity",
      random.seed               = max(1L, as.integer(random.seed)),
      verbose                   = verbose
    ))
  }

  # Louvain (1), Multilevel (2), SLM (3) — all pure igraph, no Python
  message("[Clustering] FindClusters algorithm=", algo, " (igraph)")
  FindClusters(
    object,
    graph.name  = graph.name,
    resolution  = resolution,
    algorithm   = algo,
    random.seed = as.integer(random.seed),
    verbose     = verbose
  )
}


# =============================================================================
# 15. record_clustering_run
#     Stores metadata about the clustering run inside the Seurat object itself
#     (in @misc) so that downstream modules (annotation, export, AI) can
#     inspect what parameters were used. Two slots are maintained:
#       - clustering_last_run:  the most recent run (overwritten each time)
#       - clustering_history:   append-only list of all runs in this session
#     Returns the modified Seurat object (must be captured by the caller).
# =============================================================================
record_clustering_run <- function(seurat_obj, context, mode = "single", extra = list()) {
  run_info <- c(
    list(
      mode           = mode,
      algorithm      = context$algorithm_key,
      algorithm_code = context$algorithm_code,
      reduction      = context$reduction,
      dims           = context$dims,
      k_param        = context$k_param,
      resolution     = context$resolution,
      cluster_name   = context$cluster_name,
      graph_name     = context$graph_name,
      nn_name        = context$nn_name,
      graph_type     = context$graph_type,
      prune_snn      = context$prune_snn,
      random_seed    = context$random_seed,
      timestamp      = as.character(Sys.time())
    ),
    extra
  )

  seurat_obj@misc$clustering_last_run <- run_info

  # Append to history (create the list if this is the first run)
  history <- seurat_obj@misc$clustering_history %||% list()
  seurat_obj@misc$clustering_history <- c(history, list(run_info))

  seurat_obj
}


# =============================================================================
# 16. build_resolution_grid
#     Constructs a numeric vector of resolution values for sweep mode. The
#     grid is capped at max_points to prevent excessively long runs — when
#     the natural grid exceeds the cap, we switch to a linearly spaced grid
#     of exactly max_points values and warn the user.
# =============================================================================
build_resolution_grid <- function(res_min, res_max, res_step, max_points = 20L) {
  res_min  <- as.numeric(res_min  %||% 0.1)
  res_max  <- as.numeric(res_max  %||% 2.0)
  res_step <- as.numeric(res_step %||% 0.1)

  if (is.na(res_min) || is.na(res_max) || is.na(res_step) || res_step <= 0) {
    stop("Resolution sweep parameters are invalid. Ensure min, max, and step are positive numbers.")
  }
  # Swap if user provided them in reverse order
  if (res_max < res_min) {
    tmp     <- res_min
    res_min <- res_max
    res_max <- tmp
  }

  grid <- seq(res_min, res_max, by = res_step)

  if (length(grid) > max_points) {
    notify_clustering(
      paste0("Resolution sweep capped at ", max_points,
             " grid points (requested ", length(grid),
             ") to keep runtime manageable."),
      type = "warning",
      duration = 6
    )
    grid <- seq(res_min, res_max, length.out = max_points)
  }

  unique(round(grid, 3))
}


# =============================================================================
# 17. get_color_palette
#     Generates a vector of n distinguishable colors. Uses RColorBrewer's Set3
#     palette as the base (12 colours); when more are needed, interpolates via
#     colorRampPalette. Falls back to hcl.colors("Dynamic") when RColorBrewer
#     is not installed. Guarded with exists() so it does not collide with the
#     canonical version in core_robustness.R.
# =============================================================================
if (!exists("get_color_palette", mode = "function")) {
  get_color_palette <- function(n) {
    n <- max(1L, as.integer(n))
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      base <- RColorBrewer::brewer.pal(min(12, max(3, n)), "Set3")
      if (n > length(base)) {
        grDevices::colorRampPalette(base)(n)
      } else {
        base[seq_len(n)]
      }
    } else {
      grDevices::hcl.colors(n, "Dynamic")
    }
  }
}


# =============================================================================
# 18. best_reduction
#     Returns the name of the "best" 2D reduction available for visualization.
#     Preference: UMAP > tSNE > PCA > first available. Guarded so the
#     canonical definition in core_robustness.R takes precedence when loaded.
# =============================================================================
if (!exists("best_reduction", mode = "function")) {
  best_reduction <- function(obj, prefer = c("umap", "tsne", "pca")) {
    if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
    reds <- names(obj@reductions)
    if (length(reds) == 0) return(NULL)
    for (nm in prefer) {
      if (nm %in% reds) return(nm)
    }
    reds[1]
  }
}


# =============================================================================
# =============================================================================
#
#  SECTION A:  graphClusteringUI
#  SECTION B:  graphClusteringServer
#  SECTION C:  subclusteringUI
#  SECTION D:  subclusteringServer
#
# =============================================================================
# =============================================================================


# =============================================================================
# SECTION A: graphClusteringUI
# =============================================================================

graphClusteringUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Advanced Graph-Based Clustering"),
          p("Identify cell clusters using state-of-the-art graph-based algorithms with comprehensive validation."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("clustering"),

          # Algorithm selection
          fluidRow(
            column(4,
              selectInput(ns("cluster_algorithm"), "Clustering Algorithm",
                choices = c(
                  "Leiden (Recommended)" = "leiden",
                  "Louvain" = "louvain",
                  "SLM (Smart Local Moving)" = "slm",
                  "Louvain + Multilevel Refinement" = "multilevel"
                ),
                selected = "leiden"),
              helpText("Leiden guarantees well-connected communities; recommended for publication.")
            ),
            column(4,
              selectInput(ns("reduction_cluster"), "Input Reduction", choices = NULL)
            ),
            column(4,
              selectInput(ns("resolution_mode"), "Resolution Strategy",
                choices = c(
                  "Single Resolution" = "single",
                  "Resolution Sweep" = "sweep",
                  "Optimal Resolution" = "optimal"
                ),
                selected = "single")
            )
          ),

          # --- Single resolution ---
          conditionalPanel(
            condition = "input.resolution_mode == 'single'",
            ns = ns,
            fluidRow(
              column(12,
                sliderInput(ns("resolution"), "Resolution",
                  min = 0.1, max = 3, value = 0.8, step = 0.1),
                helpText("Tip: For datasets < 30K cells, resolution 0.3-0.5 typically yields biologically meaningful clusters.")
              )
            )
          ),

          # --- Sweep mode ---
          conditionalPanel(
            condition = "input.resolution_mode == 'sweep'",
            ns = ns,
            fluidRow(
              column(4,
                numericInput(ns("res_min"), "Min Resolution", value = 0.1, min = 0.01, max = 1)
              ),
              column(4,
                numericInput(ns("res_max"), "Max Resolution", value = 2.0, min = 0.5, max = 5)
              ),
              column(4,
                numericInput(ns("res_step"), "Step Size", value = 0.1, min = 0.05, max = 0.5)
              )
            )
          ),

          # --- Optimal mode ---
          conditionalPanel(
            condition = "input.resolution_mode == 'optimal'",
            ns = ns,
            fluidRow(
              column(6,
                selectInput(ns("optimization_metric"), "Optimization Metric",
                  choices = c(
                    "Silhouette Score" = "silhouette",
                    "Modularity" = "modularity",
                    "Calinski-Harabasz" = "calinski"
                  ),
                  selected = "silhouette")
              ),
              column(6,
                numericInput(ns("target_clusters"), "Target # Clusters (optional)",
                  value = NA, min = 2, max = 50)
              )
            )
          ),

          # Graph construction parameters
          fluidRow(
            column(3,
              numericInput(ns("dims_cluster"), "Dimensions to Use",
                value = 30, min = 5, max = 100)
            ),
            column(3,
              numericInput(ns("n_neighbors_cluster"), "k-neighbors",
                value = 20, min = 5, max = 200)
            ),
            column(3,
              selectInput(ns("graph_type"), "Graph Type",
                choices = c("Pruned SNN" = "snn", "Dense SNN" = "dense_snn"),
                selected = "snn")
            ),
            column(3,
              numericInput(ns("prune_snn"), "Prune SNN",
                value = round(1/15, 4), min = 0, max = 1, step = 0.01)
            )
          ),

          # Advanced options
          fluidRow(
            column(4,
              checkboxInput(ns("compute_silhouette"), "Compute silhouette scores", TRUE)
            ),
            column(4,
              checkboxInput(ns("find_markers"), "Find cluster markers", TRUE)
            ),
            column(4,
              textInput(ns("cluster_name"), "Cluster column name", value = "seurat_clusters")
            )
          ),

          # AI buttons
          fluidRow(
            column(12,
              div(style = "margin-bottom: 8px;",
                actionButton(ns("ai_resolution_advisor"), "AI Suggest Resolution",
                  icon = icon("robot"), class = "ai-btn btn-sm")
              ),
              uiOutput(ns("ai_thread_resolution"))
            )
          ),

          actionButton(ns("run_clustering"), "Run Clustering",
            class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),

    # Results display
    fluidRow(
      column(12,
        uiOutput(ns("clustering_results_ui"))
      )
    )
  )
}


# =============================================================================
# SECTION B: graphClusteringServer
# =============================================================================

graphClusteringServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # -------------------------------------------------------------------
    # 1. Export marker test to parent so subclustering can inherit
    # -------------------------------------------------------------------
    observe({
      session$userData$parent_marker_test <- input$marker_test
    })

    # -------------------------------------------------------------------
    # 2. Leiden availability check
    # -------------------------------------------------------------------
    # Leiden availability is checked at runtime inside resolve_graph_clustering_context().
    # We do NOT disable the Leiden dropdown at startup because the igraph DLL can
    # report unavailability transiently during initial library loading on Windows.
    # If cluster_leiden is truly absent, the runtime fallback to Louvain + a
    # user notification will handle it gracefully.

    # -------------------------------------------------------------------
    # 3. Adaptive k-neighbors based on dataset size
    # -------------------------------------------------------------------
    observe({
      req(values$active_dataset)
      n_cells <- ncol(values$active_dataset)
      adaptive_k <- max(10L, min(30L, round(log10(n_cells) * 7)))
      updateNumericInput(session, "n_neighbors_cluster", value = adaptive_k)
    })

    # -------------------------------------------------------------------
    # 4. Local reactive values
    # -------------------------------------------------------------------
    show_clustering_results  <- reactiveVal(FALSE)
    clustering_results       <- reactiveVal(NULL)
    resolution_sweep_results <- reactiveVal(NULL)
    cluster_markers          <- reactiveVal(NULL)

    # -------------------------------------------------------------------
    # 5. Results UI renderer
    # -------------------------------------------------------------------
    output$clustering_results_ui <- renderUI({
      req(isTRUE(show_clustering_results()))
      ns <- session$ns

      # Pre-compute choices from current dataset to avoid timing issues
      obj <- values$active_dataset
      vis_choices   <- if (!is.null(obj) && inherits(obj, "Seurat")) names(obj@reductions) else c("umap", "pca")
      vis_selected  <- if (!is.null(obj) && inherits(obj, "Seurat")) {
        red <- best_reduction(obj)
        if (!is.null(red)) red else vis_choices[1]
      } else vis_choices[1]

      color_choices <- if (!is.null(obj) && inherits(obj, "Seurat")) {
        cname <- input$cluster_name %||% "seurat_clusters"
        c("Clusters" = cname, colnames(obj@meta.data))
      } else c("seurat_clusters")

      marker_cluster_choices <- tryCatch({
        res <- clustering_results()
        if (!is.null(res) && !is.null(res$cluster_labels)) sort(unique(res$cluster_labels)) else NULL
      }, error = function(e) NULL)

      wellPanel(
        h4("Clustering Results"),
        tabsetPanel(
          # --- Tab 1: Cluster Visualization ---
          tabPanel("Cluster Visualization",
            fluidRow(
              column(8,
                plotlyOutput(ns("cluster_plot"), height = "600px")
              ),
              column(4,
                h5("Visualization Options"),
                selectInput(ns("vis_reduction"), "Reduction",
                  choices = vis_choices, selected = vis_selected),
                selectInput(ns("color_by_cluster"), "Color by",
                  choices = color_choices, selected = color_choices[1]),
                sliderInput(ns("point_size_cluster"), "Point size",
                  min = 0.1, max = 5, value = 1, step = 0.1),
                sliderInput(ns("point_alpha_cluster"), "Opacity",
                  min = 0.1, max = 1, value = 0.8, step = 0.1),
                checkboxInput(ns("show_cluster_centers"), "Show cluster centers", FALSE),
                br(),
                downloadButton(ns("download_cluster_plot"), "Save Plot", class = "btn-success")
              )
            )
          ),

          # --- Tab 2: Cluster Statistics ---
          tabPanel("Cluster Statistics",
            fluidRow(
              column(6,
                h5("Cluster Summary"),
                downloadButton(ns("dl_cluster_summary"), "CSV", class = "btn-xs btn-default"),
                DT::dataTableOutput(ns("cluster_summary"))
              ),
              column(6,
                h5("Size Distribution"),
                plotlyOutput(ns("cluster_distribution"), height = "400px")
              )
            )
          ),

          # --- Tab 3: Quality Metrics ---
          tabPanel("Quality Metrics",
            fluidRow(
              column(6,
                h5("Clustering Quality"),
                downloadButton(ns("dl_quality_metrics_table"), "CSV", class = "btn-xs btn-default"),
                DT::dataTableOutput(ns("quality_metrics_table"))
              ),
              column(6,
                h5("Silhouette Analysis"),
                plotlyOutput(ns("silhouette_plot"), height = "400px")
              )
            )
          ),

          # --- Tab 4: Cluster Markers ---
          tabPanel("Cluster Markers",
            fluidRow(
              column(4,
                h5("Marker Options"),
                numericInput(ns("top_markers"), "Top markers per cluster",
                  value = 10, min = 5, max = 50),
                selectInput(ns("marker_test"), "Statistical test",
                  choices = c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR"),
                  selected = "wilcox"),
                actionButton(ns("find_all_markers"), "Find All Cluster Markers",
                  class = "btn-primary", width = "100%"),
                hr(),
                h5("Single Cluster"),
                selectInput(ns("marker_cluster"), "Select cluster",
                  choices = marker_cluster_choices),
                actionButton(ns("find_cluster_markers"), "Find Markers for Selected",
                  class = "btn-info", width = "100%"),
                hr(),
                actionButton(ns("ai_interpret_markers"), "AI Interpret Markers",
                  icon = icon("robot"), class = "ai-btn", width = "100%")
              ),
              column(8,
                downloadButton(ns("dl_cluster_markers_table"), "CSV", class = "btn-xs btn-default"),
                DT::dataTableOutput(ns("cluster_markers_table")),
                uiOutput(ns("ai_marker_interpretation"))
              )
            )
          ),

          # --- Tab 5: Resolution Comparison ---
          tabPanel("Resolution Comparison",
            conditionalPanel(
              condition = "input.resolution_mode == 'sweep'",
              ns = ns,
              fluidRow(
                column(6,
                  h5("Resolution Sweep Results"),
                  downloadButton(ns("dl_resolution_sweep_table"), "CSV", class = "btn-xs btn-default"),
                  DT::dataTableOutput(ns("resolution_sweep_table"))
                ),
                column(6,
                  h5("Quality vs Resolution"),
                  plotlyOutput(ns("resolution_quality_plot"), height = "400px")
                )
              )
            )
          ),

          # --- Tab 6: Clustree (Stability) ---
          tabPanel("Clustree (Stability)",
            fluidRow(
              column(12,
                actionButton(ns("run_clustree"), "Generate Clustree",
                  class = "btn-info", icon = icon("sitemap")),
                helpText("Clustree shows how clusters split/merge across resolutions. Requires a resolution sweep first or multiple clustering columns."),
                plotOutput(ns("clustree_plot"), height = "700px")
              )
            )
          ),

          # --- Tab 7: 3D View ---
          tabPanel("3D View",
            fluidRow(
              column(4,
                actionButton(ns("compute_3d_umap"), "Compute 3D UMAP",
                  class = "btn-info", icon = icon("cube")),
                selectInput(ns("color_by_3d"), "Color by", choices = NULL)
              ),
              column(8,
                plotlyOutput(ns("cluster_plot_3d"), height = "600px")
              )
            )
          )
        )
      )
    })

    # -------------------------------------------------------------------
    # 6. Input updaters
    # -------------------------------------------------------------------

    # Update reduction choices when dataset changes
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      reductions <- names(obj@reductions)
      if (length(reductions) == 0) return()

      # For clustering input: prefer PCA > Harmony
      cluster_red <- if ("pca" %in% reductions) "pca"
                     else if ("harmony" %in% reductions) "harmony"
                     else if (exists("best_reduction", mode = "function")) best_reduction(obj)
                     else reductions[1]

      updateSelectInput(session, "reduction_cluster",
        choices = reductions, selected = cluster_red)
    })

    # Update marker cluster dropdown when clustering results arrive
    observeEvent(clustering_results(), {
      req(clustering_results())
      clusters <- unique(clustering_results()$cluster_labels)
      if (length(clusters) > 0) {
        updateSelectInput(session, "marker_cluster", choices = sort(clusters))
      }
    })

    # -------------------------------------------------------------------
    # 7. Main clustering execution
    # -------------------------------------------------------------------
    observeEvent(input$run_clustering, {
      req(values$active_dataset)

      job_id <- begin_clustering_job(values, "run_clustering", meta = list(
        algorithm       = input$cluster_algorithm,
        reduction       = input$reduction_cluster,
        resolution_mode = input$resolution_mode
      ))

      tryCatch({
        obj <- values$active_dataset
        selected_reduction <- input$reduction_cluster

        # --- Guard: dataset loaded ---
        if (is.null(obj) || !inherits(obj, "Seurat")) {
          showNotification("No valid Seurat object loaded.", type = "error", duration = 8)
          end_clustering_job(values, job_id, status = "failed", message = "no dataset")
          return()
        }

        # --- Guard: reduction selected / auto-select ---
        if (is.null(selected_reduction) || !nzchar(selected_reduction)) {
          auto_red <- select_clustering_reduction(obj)
          if (is.null(auto_red)) {
            showNotification("No suitable reduction found. Run PCA or Harmony first.", type = "error", duration = 8)
            end_clustering_job(values, job_id, status = "failed", message = "no reduction")
            return()
          }
          selected_reduction <- auto_red
          updateSelectInput(session, "reduction_cluster", selected = auto_red)
          showNotification(paste("Auto-selected reduction:", auto_red), type = "warning", duration = 6)
        }

        if (!selected_reduction %in% names(obj@reductions)) {
          fallback <- select_clustering_reduction(obj)
          if (!is.null(fallback)) {
            selected_reduction <- fallback
            updateSelectInput(session, "reduction_cluster", selected = fallback)
            showNotification(paste("Selected reduction missing; using fallback:", fallback), type = "warning", duration = 8)
          } else {
            showNotification("Selected reduction not found. Run dimensionality reduction first.", type = "error", duration = 8)
            end_clustering_job(values, job_id, status = "failed", message = "missing reduction")
            return()
          }
        }

        # --- Guard: block UMAP/tSNE clustering ---
        if (selected_reduction %in% c("umap", "tsne", "umap3d", "tsne3d")) {
          showNotification(
            paste0("Clustering on ", toupper(selected_reduction), " is methodologically incorrect. ",
                   "These 2D embeddings discard most biological variance. Use PCA or Harmony instead."),
            type = "error", duration = 15
          )
          better_red <- intersect(c("pca", "harmony", "integrated"), names(obj@reductions))[1]
          if (!is.na(better_red)) updateSelectInput(session, "reduction_cluster", selected = better_red)
          end_clustering_job(values, job_id, status = "failed", message = "invalid reduction for clustering")
          return()
        }

        # --- Contract validation ---
        if (!validate_clustering_contract(values, obj, require_reductions = selected_reduction)) {
          end_clustering_job(values, job_id, status = "failed", message = "contract validation failed")
          return()
        }

        # --- Dimension validation and sanitization ---
        avail_dims <- ncol(Embeddings(obj, selected_reduction))
        if (avail_dims < 2) {
          showNotification(
            paste("Reduction '", selected_reduction, "' has only", avail_dims, "dimension(s). Need at least 2."),
            type = "error", duration = 8)
          end_clustering_job(values, job_id, status = "failed", message = "insufficient dimensions")
          return()
        }

        dims_use <- min(input$dims_cluster %||% 30L, avail_dims)
        if ((input$dims_cluster %||% 30L) > avail_dims) {
          showNotification(
            paste("Requested", input$dims_cluster, "dims but", selected_reduction, "only has", avail_dims,
                  ". Will use", avail_dims, "dims."),
            type = "warning", duration = 6)
        }

        withProgress(message = "Running clustering analysis...", value = 0, {

          # --- Runtime budget guard ---
          if (exists("runtime_budget_guard", mode = "function")) {
            budget <- runtime_budget_guard(obj, task = "clustering")
            if (!isTRUE(budget$ok)) {
              showNotification(paste(budget$message, budget$suggestion %||% ""), type = "error", duration = 10)
              end_clustering_job(values, job_id, status = "failed", message = budget$message %||% "budget exceeded")
              return()
            }
          }

          incProgress(0.2, detail = "Building neighbor graph...")

          # --- Build params from input, sanitize NA values ---
          params <- reactiveValuesToList(input)
          params$reduction_cluster <- selected_reduction
          params$dims_cluster      <- dims_use
          # Sanitize NA values in critical params
          params$n_neighbors_cluster <- params$n_neighbors_cluster %||% 20L
          params$prune_snn           <- params$prune_snn %||% round(1/15, 4)
          params$cluster_name        <- params$cluster_name %||% "seurat_clusters"
          params$resolution          <- params$resolution %||% 0.8

          # --- Dispatch to run mode ---
          result <- tryCatch({
            if (input$resolution_mode == "single") {
              run_single_resolution_clustering(obj, params)
            } else if (input$resolution_mode == "sweep") {
              run_resolution_sweep_clustering(obj, params)
            } else if (input$resolution_mode == "optimal") {
              run_optimal_resolution_clustering(obj, params)
            } else {
              run_single_resolution_clustering(obj, params)
            }
          }, error = function(e) {
            list(success = FALSE, error = conditionMessage(e))
          })

          # --- Auto-retry with k=20 if single cluster detected and k was higher ---
          if (!isTRUE(result$success) &&
              !is.null(result$error) &&
              grepl("only 1 cluster", result$error, ignore.case = TRUE) &&
              (params$n_neighbors_cluster %||% 20L) > 20L) {
            showNotification("Only 1 cluster found with current k. Retrying with k=20...", type = "warning", duration = 5)
            params$n_neighbors_cluster <- 20L
            result <- tryCatch({
              if (input$resolution_mode == "single") {
                run_single_resolution_clustering(obj, params)
              } else if (input$resolution_mode == "sweep") {
                run_resolution_sweep_clustering(obj, params)
              } else if (input$resolution_mode == "optimal") {
                run_optimal_resolution_clustering(obj, params)
              } else {
                run_single_resolution_clustering(obj, params)
              }
            }, error = function(e) {
              list(success = FALSE, error = conditionMessage(e))
            })
          }

          # --- On success ---
          if (isTRUE(result$success)) {
            values$active_dataset  <- result$object
            values$clustering_done <- TRUE
            clustering_results(result)

            if (input$resolution_mode == "sweep" && !is.null(result$sweep_results)) {
              resolution_sweep_results(result$sweep_results)
            }

            show_clustering_results(TRUE)

            # Checkpoint
            checkpoint_clustering(
              values,
              action = paste0("run_clustering_", input$resolution_mode),
              obj    = result$object,
              params = list(
                algorithm       = input$cluster_algorithm,
                reduction       = selected_reduction,
                dims            = dims_use,
                neighbors       = params$n_neighbors_cluster,
                resolution_mode = input$resolution_mode,
                resolution      = if (input$resolution_mode == "single") params$resolution else NA
              ),
              method = paste0("Seurat::FindNeighbors + Seurat::FindClusters (", input$cluster_algorithm, ")")
            )

            end_clustering_job(values, job_id, status = "done", message = "clustering completed")

            n_clusters <- length(unique(result$cluster_labels))
            showNotification(paste("Clustering complete!", n_clusters, "clusters found."), type = "message")
          } else {
            # --- On failure ---
            end_clustering_job(values, job_id, status = "failed", message = result$error %||% "unknown error")
            showNotification(paste("Clustering failed:", result$error %||% "Unknown error"), type = "error")
          }

          incProgress(1)
        })

      }, error = function(e) {
        # --- Outer tryCatch for uncaught errors ---
        end_clustering_job(values, job_id, status = "error", message = conditionMessage(e))
        showNotification(paste("Clustering error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })

    # -------------------------------------------------------------------
    # 8. Find markers handlers
    # -------------------------------------------------------------------

    observeEvent(input$find_all_markers, {
      req(values$active_dataset, clustering_results())
      tryCatch({
        withProgress(message = "Finding markers for all clusters...", value = 0, {
          result <- if (exists("find_all_cluster_markers", mode = "function")) {
            find_all_cluster_markers(values$active_dataset, input)
          } else {
            list(success = FALSE, error = "find_all_cluster_markers not available")
          }
          if (isTRUE(result$success)) {
            cluster_markers(result$markers)
            n_clust <- length(unique(result$markers$cluster))
            showNotification(paste("Markers found for", n_clust, "clusters."), type = "message")
          } else {
            showNotification(paste("Marker finding failed:", result$error %||% "Unknown error"), type = "error")
          }
        })
      }, error = function(e) {
        showNotification(paste("Marker finding error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })

    observeEvent(input$find_cluster_markers, {
      req(values$active_dataset, input$marker_cluster)
      tryCatch({
        withProgress(message = "Finding cluster markers...", value = 0, {
          result <- if (exists("find_cluster_specific_markers", mode = "function")) {
            find_cluster_specific_markers(values$active_dataset, input)
          } else {
            list(success = FALSE, error = "find_cluster_specific_markers not available")
          }
          if (isTRUE(result$success)) {
            cluster_markers(result$markers)
            showNotification("Markers found.", type = "message")
          } else {
            showNotification(paste("Marker finding failed:", result$error %||% "Unknown error"), type = "error")
          }
        })
      }, error = function(e) {
        showNotification(paste("Marker finding error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })

    # -------------------------------------------------------------------
    # 9. AI marker interpretation
    # -------------------------------------------------------------------
    observeEvent(input$ai_interpret_markers, {
      # Guard API key
      if (is.null(values$ai_api_key) || !nzchar(trimws(values$ai_api_key %||% ""))) {
        showNotification("Please set your Claude API key in Data Management > AI Config first.",
          type = "warning", duration = 8)
        return()
      }
      markers <- cluster_markers()
      if (is.null(markers) || nrow(markers) == 0) {
        showNotification("Find cluster markers first, then interpret.", type = "warning", duration = 6)
        return()
      }

      # Loading state
      output$ai_marker_interpretation <- renderUI({
        if (exists("ai_thread_loading", mode = "function")) {
          ai_thread_loading("AI interpreting cluster markers...")
        } else {
          div(icon("spinner", class = "fa-spin"), "AI interpreting cluster markers...")
        }
      })

      # Build marker summary per cluster
      clusters <- unique(markers$cluster)
      marker_summary <- paste(lapply(clusters, function(cl) {
        cl_markers <- markers[markers$cluster == cl, ]
        top_genes  <- head(cl_markers$gene, 20)
        fc_vals    <- tryCatch(head(round(cl_markers$avg_log2FC, 2), 20), error = function(e) rep(NA, length(top_genes)))
        paste0("Cluster ", cl, ": ", paste(paste0(top_genes, " (FC=", fc_vals, ")"), collapse = ", "))
      }), collapse = "\n")

      tissue  <- if (!is.null(values$tissue_global) && nzchar(values$tissue_global %||% "")) values$tissue_global else "unknown"
      species <- if (!is.null(values$species_global) && nzchar(values$species_global %||% "")) values$species_global else "unknown"

      prompt_text <- paste0(
        "## Task: Interpret Cluster Markers\n\n",
        "Below are the top marker genes per cluster from a scRNA-seq clustering analysis.\n",
        "Tissue: ", tissue, " | Species: ", species, "\n\n",
        marker_summary, "\n\n",
        "For each cluster, provide:\n",
        "1. **Predicted Cell Type** - the most likely cell identity based on these markers\n",
        "2. **Key Evidence** - which markers support this identity and why\n",
        "3. **Biological Notes** - any notable activation states, subtypes, or functional programs\n\n",
        "Be concise but specific. Use standard cell type nomenclature."
      )

      prompt <- if (exists("build_expert_prompt", mode = "function")) {
        build_expert_prompt(prompt_text, values = values, include_project = TRUE, include_dataset = TRUE)
      } else {
        prompt_text
      }

      result <- tryCatch({
        if (exists("call_ai", mode = "function")) {
          withProgress(message = "AI interpreting markers...", value = 0.5, {
            call_ai("de_interpretation", prompt, values)
          })
        } else {
          list(ok = FALSE, error = "AI assistant (call_ai) not available.")
        }
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })

      if (isTRUE(result$ok)) {
        output$ai_marker_interpretation <- renderUI({
          content_html <- if (exists(".mdtohtml", mode = "function")) .mdtohtml(result$text) else result$text
          model_name   <- if (exists("get_model_for_task", mode = "function")) get_model_for_task("de_interpretation", values) else NULL
          if (exists("ai_thread_message", mode = "function")) {
            ai_thread_message(
              content_html = content_html,
              title        = "AI Marker Interpretation",
              model        = model_name,
              tokens_in    = result$tokens_in %||% NA,
              tokens_out   = result$tokens_out %||% NA
            )
          } else {
            div(class = "ai-result", HTML(content_html))
          }
        })
      } else {
        output$ai_marker_interpretation <- renderUI({
          if (exists("ai_thread_error", mode = "function")) {
            ai_thread_error(result$error %||% "Unknown AI error")
          } else {
            div(class = "text-danger", result$error %||% "Unknown AI error")
          }
        })
      }
    })

    # -------------------------------------------------------------------
    # 10. Result renderers
    # -------------------------------------------------------------------

    output$cluster_plot <- renderPlotly({
      req(clustering_results(), values$active_dataset)
      tryCatch({
        p <- if (exists("create_enhanced_cluster_plot", mode = "function")) {
          create_enhanced_cluster_plot(values$active_dataset, clustering_results(), input)
        } else {
          plotly_empty() %>% layout(title = "create_enhanced_cluster_plot not available")
        }
        if (exists("capture_umap_snapshot", mode = "function")) {
          capture_umap_snapshot(values, p, "Clustering",
            coloring = "clusters",
            method   = tryCatch(clustering_results()$algorithm, error = function(e) ""))
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type   = "Cluster UMAP/tSNE embedding colored by cluster identity",
            label         = "Cluster Plot",
            source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    output$cluster_summary <- DT::renderDataTable({
      req(clustering_results())
      tbl <- if (exists("create_cluster_summary_table", mode = "function")) {
        create_cluster_summary_table(clustering_results())
      } else {
        data.frame(Cluster = unique(clustering_results()$cluster_labels),
                   Count = as.numeric(table(clustering_results()$cluster_labels)))
      }
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Cluster Summary", source_module = "Clustering")
      }
      tbl
    }, options = list(pageLength = 15, scrollX = TRUE))

    output$cluster_distribution <- renderPlotly({
      req(clustering_results())
      tryCatch({
        p <- if (exists("create_cluster_distribution_plot", mode = "function")) {
          create_cluster_distribution_plot(clustering_results())
        } else {
          labels <- clustering_results()$cluster_labels
          counts <- as.data.frame(table(Cluster = labels))
          plot_ly(counts, x = ~Cluster, y = ~Freq, type = "bar") %>%
            layout(title = "Cluster Size Distribution", autosize = TRUE)
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type = "Bar chart of cell counts per cluster",
            label = "Cluster Distribution", source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Distribution plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    output$quality_metrics_table <- DT::renderDataTable({
      req(clustering_results())
      tbl <- if (exists("create_quality_metrics_table", mode = "function")) {
        create_quality_metrics_table(clustering_results())
      } else {
        res <- clustering_results()
        data.frame(
          Metric = c("Silhouette (avg)", "Modularity", "Calinski-Harabasz"),
          Value  = c(
            round(res$silhouette_avg %||% NA_real_, 4),
            round(res$modularity %||% NA_real_, 4),
            round(res$calinski_score %||% NA_real_, 4)
          )
        )
      }
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Cluster Quality Metrics", source_module = "Clustering")
      }
      tbl
    }, options = list(dom = "t"))

    output$silhouette_plot <- renderPlotly({
      req(clustering_results())
      if (!is.null(clustering_results()$silhouette_scores)) {
        tryCatch({
          p <- if (exists("create_silhouette_plot", mode = "function")) {
            create_silhouette_plot(clustering_results())
          } else {
            plotly_empty() %>% layout(title = "create_silhouette_plot not available", autosize = TRUE)
          }
          if (exists("capture_figure", mode = "function")) {
            capture_figure(values, p, "Clustering",
              figure_type = "Silhouette score distribution per cluster",
              label = "Silhouette Plot", source_module = "Clustering")
          }
          p
        }, error = function(e) {
          plotly_empty() %>% layout(title = paste("Silhouette plot error:", conditionMessage(e)), autosize = TRUE)
        })
      } else {
        plotly_empty() %>% layout(title = "Silhouette analysis not computed", autosize = TRUE)
      }
    })

    output$cluster_markers_table <- DT::renderDataTable({
      req(cluster_markers())
      tbl <- cluster_markers()
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Cluster Markers", source_module = "Clustering")
      }
      tbl
    }, options = list(pageLength = 15, scrollX = TRUE))

    output$resolution_sweep_table <- DT::renderDataTable({
      req(resolution_sweep_results())
      tbl <- resolution_sweep_results()
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Resolution Sweep", source_module = "Clustering")
      }
      tbl
    }, options = list(pageLength = 10, scrollX = TRUE))

    output$resolution_quality_plot <- renderPlotly({
      req(resolution_sweep_results())
      tryCatch({
        p <- if (exists("create_resolution_quality_plot", mode = "function")) {
          create_resolution_quality_plot(resolution_sweep_results())
        } else {
          plotly_empty() %>% layout(title = "create_resolution_quality_plot not available", autosize = TRUE)
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type = "Line plot of clustering quality metrics across resolutions",
            label = "Resolution Quality", source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Quality plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    # -------------------------------------------------------------------
    # 11. Download handlers
    # -------------------------------------------------------------------

    output$download_cluster_plot <- downloadHandler(
      filename = function() paste0("cluster_plot_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- if (exists("create_enhanced_cluster_plot", mode = "function")) {
            create_enhanced_cluster_plot(values$active_dataset, clustering_results(), input)
          } else {
            plotly_empty()
          }
          htmlwidgets::saveWidget(p, file)
        }, error = function(e) {
          showNotification(paste("Export error:", conditionMessage(e)), type = "error")
          htmlwidgets::saveWidget(plotly_empty(), file)
        })
      }
    )

    output$dl_cluster_summary <- downloadHandler(
      filename = function() paste0("cluster_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(clustering_results())
        tryCatch({
          tbl <- if (exists("create_cluster_summary_table", mode = "function")) {
            create_cluster_summary_table(clustering_results())
          } else {
            data.frame(Cluster = unique(clustering_results()$cluster_labels))
          }
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = conditionMessage(e)), file, row.names = FALSE))
      }
    )

    output$dl_quality_metrics_table <- downloadHandler(
      filename = function() paste0("quality_metrics_", Sys.Date(), ".csv"),
      content = function(file) {
        req(clustering_results())
        tryCatch({
          tbl <- if (exists("create_quality_metrics_table", mode = "function")) {
            create_quality_metrics_table(clustering_results())
          } else {
            data.frame(Metric = "N/A", Value = NA)
          }
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = conditionMessage(e)), file, row.names = FALSE))
      }
    )

    output$dl_cluster_markers_table <- downloadHandler(
      filename = function() paste0("cluster_markers_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          req(cluster_markers())
          write.csv(cluster_markers(), file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Download error:", conditionMessage(e)), type = "error")
          write.csv(data.frame(Error = conditionMessage(e)), file, row.names = FALSE)
        })
      }
    )

    output$dl_resolution_sweep_table <- downloadHandler(
      filename = function() paste0("resolution_sweep_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          req(resolution_sweep_results())
          write.csv(resolution_sweep_results(), file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Download error:", conditionMessage(e)), type = "error")
          write.csv(data.frame(Error = conditionMessage(e)), file, row.names = FALSE)
        })
      }
    )

    # -------------------------------------------------------------------
    # 12. AI resolution advisor
    # -------------------------------------------------------------------
    observeEvent(input$ai_resolution_advisor, {
      if (is.null(values$ai_api_key) || !nzchar(trimws(values$ai_api_key %||% ""))) {
        showNotification("Please set your Claude API key in the Data Management tab first.",
          type = "warning", duration = 8)
        return()
      }
      if (is.null(values$active_dataset)) {
        showNotification("No active dataset loaded.", type = "warning", duration = 6)
        return()
      }

      obj     <- values$active_dataset
      n_cells <- ncol(obj)
      tissue  <- if (!is.null(values$tissue_global) && nzchar(values$tissue_global %||% "")) values$tissue_global else "unknown tissue"
      species <- if (!is.null(values$species_global) && nzchar(values$species_global %||% "")) values$species_global else "unknown species"

      # Collect existing cluster columns and their cluster counts
      md <- obj@meta.data
      cluster_cols <- grep("cluster|leiden|louvain|seurat_clusters", colnames(md),
                           value = TRUE, ignore.case = TRUE)
      cluster_info <- if (length(cluster_cols) > 0) {
        paste(vapply(cluster_cols, function(col) {
          n_clust <- length(unique(md[[col]]))
          paste0(col, ": ", n_clust, " clusters")
        }, character(1)), collapse = "; ")
      } else {
        "none yet"
      }

      prompt <- paste0(
        "I have ", format(n_cells, big.mark = ","), " cells from ", tissue, " (", species, "). ",
        "Current clustering columns and cluster counts: ", cluster_info, ". ",
        "Suggest an optimal resolution for Leiden/Louvain clustering. ",
        "Consider expected cell type diversity for this tissue. ",
        "Return the suggested resolution (numeric) and brief reasoning."
      )

      # Show loading state
      output$ai_thread_resolution <- renderUI({
        if (exists("ai_thread_loading", mode = "function")) {
          ai_thread_loading("AI analyzing clustering resolution...")
        } else {
          div(icon("spinner", class = "fa-spin"), "AI analyzing clustering resolution...")
        }
      })

      result <- tryCatch({
        if (exists("call_ai", mode = "function")) {
          call_ai("resolution_advisor", prompt, values)
        } else {
          list(ok = FALSE, error = "AI assistant (call_ai) not available.")
        }
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })

      if (isTRUE(result$ok)) {
        # Extract suggested resolution number from AI response
        res_num <- suppressWarnings(
          as.numeric(regmatches(result$text, regexpr("\\b(0\\.\\d+|[12]\\.\\d+|[12])\\b", result$text)))
        )
        ai_res_val <- if (length(res_num) > 0 && !is.na(res_num[1]) && res_num[1] >= 0.1 && res_num[1] <= 3.0) {
          res_num[1]
        } else {
          NULL
        }

        apply_btn <- if (!is.null(ai_res_val)) {
          actionButton(ns("apply_ai_resolution"), paste0("Apply Resolution ", ai_res_val),
            class = "ai-thread-apply-btn", icon = icon("check"))
        } else {
          NULL
        }

        output$ai_thread_resolution <- renderUI({
          content_html <- if (exists(".mdtohtml", mode = "function")) .mdtohtml(result$text) else result$text
          model_name   <- if (exists("get_model_for_task", mode = "function")) get_model_for_task("resolution_advisor", values) else NULL
          if (exists("ai_thread_message", mode = "function")) {
            ai_thread_message(
              content_html = content_html,
              title        = "AI Resolution Suggestion",
              model        = model_name,
              tokens_in    = result$tokens_in %||% NA,
              tokens_out   = result$tokens_out %||% NA,
              apply_btn    = apply_btn
            )
          } else {
            tagList(div(class = "ai-result", HTML(content_html)), apply_btn)
          }
        })

        values$.ai_resolution_suggestion <- ai_res_val
      } else {
        output$ai_thread_resolution <- renderUI({
          if (exists("ai_thread_error", mode = "function")) {
            ai_thread_error(result$error %||% "Unknown AI error")
          } else {
            div(class = "text-danger", result$error %||% "Unknown AI error")
          }
        })
      }
    })

    # -------------------------------------------------------------------
    # 13. Apply AI resolution
    # -------------------------------------------------------------------
    observeEvent(input$apply_ai_resolution, {
      res_val <- values$.ai_resolution_suggestion
      req(!is.null(res_val))
      updateSliderInput(session, "resolution", value = res_val)
      showNotification(paste0("Resolution set to ", res_val, "."), type = "message", duration = 4)
    })

    # -------------------------------------------------------------------
    # 14. AI helper server
    # -------------------------------------------------------------------
    if (exists("ai_helper_server", mode = "function")) {
      ai_helper_server("graphCluster", values, "Clustering", context_fn = function() {
        obj <- values$active_dataset
        if (is.null(obj) || !inherits(obj, "Seurat")) return("")
        parts <- character(0)
        n_cells <- ncol(obj)
        parts <- c(parts, paste0("Total cells: ", format(n_cells, big.mark = ",")))

        md <- obj@meta.data
        cluster_cols <- grep("cluster|leiden|louvain|seurat_clusters", colnames(md),
                             value = TRUE, ignore.case = TRUE)
        if (length(cluster_cols) > 0) {
          for (col in cluster_cols) {
            n_clust  <- length(unique(md[[col]]))
            sizes    <- sort(table(md[[col]]), decreasing = TRUE)
            top_sizes <- head(sizes, 5)
            size_str  <- paste0(names(top_sizes), " (n=", top_sizes, ")", collapse = ", ")
            parts <- c(parts, paste0("Column '", col, "': ", n_clust, " clusters. Largest: ", size_str))
          }
        } else {
          parts <- c(parts, "No clustering columns found yet.")
        }

        reds <- names(obj@reductions)
        if (length(reds) > 0) {
          parts <- c(parts, paste0("Available reductions: ", paste(reds, collapse = ", ")))
        }
        paste(parts, collapse = "\n")
      })
    }

    # -------------------------------------------------------------------
    # 15. Clustree integration
    # -------------------------------------------------------------------
    # Build clustree plot eagerly on button click, store in reactiveVal.
    # renderPlot reads only from the reactiveVal so there is no lazy
    # dependency on outer observer-scope variables.
    clustree_plot_rv <- reactiveVal(NULL)

    observeEvent(input$run_clustree, {
      req(values$active_dataset)
      if (!requireNamespace("clustree", quietly = TRUE)) {
        showNotification("Package 'clustree' is not installed.", type = "error", duration = 8)
        return()
      }
      if (!requireNamespace("ggraph", quietly = TRUE)) {
        showNotification("Package 'ggraph' is required by clustree.", type = "error", duration = 8)
        return()
      }

      obj <- values$active_dataset
      md  <- obj@meta.data

      cand <- unique(c(
        grep("_res\\.[0-9]", colnames(md), value = TRUE),
        grep("^cluster_res_[0-9]", colnames(md), value = TRUE)
      ))
      cand <- cand[vapply(cand, function(cc) {
        v <- md[[cc]]
        if (is.null(v)) return(FALSE)
        length(unique(v[!is.na(v)])) >= 2
      }, logical(1))]

      if (length(cand) < 2) {
        showNotification(
          "Clustree requires at least 2 resolution columns with >=2 clusters each. Run a resolution sweep first.",
          type = "warning", duration = 8)
        return()
      }

      get_prefix <- function(x) sub("[0-9]+(\\.[0-9]+)?$", "", x)
      prefs <- vapply(cand, get_prefix, character(1))
      tbl   <- sort(table(prefs), decreasing = TRUE)
      chosen_prefix <- names(tbl)[1]
      if (is.na(chosen_prefix) || tbl[[1]] < 2) {
        showNotification(
          "Could not find >=2 resolution columns sharing a common prefix.",
          type = "warning", duration = 8)
        return()
      }
      res_cols <- cand[prefs == chosen_prefix]

      md_cl <- md[, res_cols, drop = FALSE]
      for (cc in res_cols) md_cl[[cc]] <- as.factor(md_cl[[cc]])

      p <- tryCatch(
        clustree::clustree(md_cl, prefix = chosen_prefix) +
          # Suppress ggraph's custom "edge_colourbar" guide (breaks on
          # ggplot2 >= 4.0 S7 guide dispatch) — the tree is readable without it.
          ggplot2::guides(edge_colour = "none",
                          edge_alpha  = "none",
                          edge_width  = "none") +
          ggplot2::theme(legend.position = "right"),
        error = function(e) e
      )
      if (inherits(p, "error")) {
        showNotification(paste("Clustree error:", conditionMessage(p)),
                         type = "error", duration = 10)
        clustree_plot_rv(NULL)
        return()
      }
      force(p)
      clustree_plot_rv(p)
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

    # -------------------------------------------------------------------
    # 16. 3D UMAP
    # -------------------------------------------------------------------

    # Compute 3D UMAP
    observeEvent(input$compute_3d_umap, {
      req(values$active_dataset)
      obj <- values$active_dataset

      # Need a reduction to build 3D UMAP from
      red_source <- select_clustering_reduction(obj)
      if (is.null(red_source)) {
        showNotification("No suitable reduction for 3D UMAP. Run PCA first.", type = "error", duration = 8)
        return()
      }

      tryCatch({
        withProgress(message = "Computing 3D UMAP...", value = 0.3, {
          n_dims <- min(30L, ncol(Embeddings(obj, red_source)))
          obj <- RunUMAP(obj,
            reduction      = red_source,
            dims           = 1:n_dims,
            n.components   = 3L,
            reduction.name = "umap3d",
            reduction.key  = "UMAP3D_",
            verbose        = FALSE
          )
          values$active_dataset <- obj
          incProgress(0.7, detail = "Done")
        })
        showNotification("3D UMAP computed.", type = "message", duration = 4)
      }, error = function(e) {
        showNotification(paste("3D UMAP error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })

    # 3D UMAP plot
    output$cluster_plot_3d <- renderPlotly({
      req(values$active_dataset)
      obj <- values$active_dataset
      if (!"umap3d" %in% names(obj@reductions)) {
        return(plotly_empty() %>% layout(title = "Click 'Compute 3D UMAP' first", autosize = TRUE))
      }

      tryCatch({
        emb <- Embeddings(obj, "umap3d")
        md  <- obj@meta.data

        color_col <- input$color_by_3d %||% "seurat_clusters"
        if (!color_col %in% colnames(md)) {
          color_col <- colnames(md)[1]
        }

        color_vals <- as.character(md[[color_col]])
        unique_vals <- sort(unique(color_vals))
        colors <- if (exists("get_color_palette", mode = "function")) {
          get_color_palette(length(unique_vals))
        } else {
          grDevices::hcl.colors(length(unique_vals), "Dynamic")
        }

        plot_ly(
          x = emb[, 1], y = emb[, 2], z = emb[, 3],
          color = color_vals,
          colors = colors,
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 2, opacity = 0.7),
          text = paste0("Cell: ", rownames(emb), "\n", color_col, ": ", color_vals),
          hoverinfo = "text"
        ) %>% layout(
          title = paste("3D UMAP colored by", color_col),
          scene = list(
            xaxis = list(title = "UMAP3D_1"),
            yaxis = list(title = "UMAP3D_2"),
            zaxis = list(title = "UMAP3D_3")
          ),
          autosize = TRUE
        )
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("3D plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    # Update color_by_3d choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      md_cols <- colnames(obj@meta.data)
      cat_cols <- md_cols[vapply(obj@meta.data, function(x) is.factor(x) || is.character(x), logical(1))]
      cluster_name <- input$cluster_name %||% "seurat_clusters"
      choices <- unique(c(cluster_name, cat_cols))
      selected <- if (cluster_name %in% choices) cluster_name else choices[1]
      updateSelectInput(session, "color_by_3d", choices = choices, selected = selected)
    })

  })
}


# =============================================================================
# SECTION C: subclusteringUI
# =============================================================================

subclusteringUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Sub-clustering Analysis"),
          p("Perform focused clustering on specific cell populations to reveal finer structure."),

          fluidRow(
            column(4,
              selectInput(ns("parent_cluster_col"), "Parent clustering", choices = NULL)
            ),
            column(4,
              uiOutput(ns("clusters_to_subcluster_ui"))
            ),
            column(4,
              textInput(ns("subcluster_name"), "Sub-cluster column name", value = "subclusters")
            )
          ),

          # Sub-clustering parameters
          fluidRow(
            column(3,
              numericInput(ns("sub_resolution"), "Resolution", value = 0.5, min = 0.1, max = 2, step = 0.1)
            ),
            column(3,
              selectInput(ns("sub_algorithm"), "Algorithm",
                choices = c("Leiden" = "leiden", "Louvain" = "louvain", "SLM" = "slm"),
                selected = "leiden")
            ),
            column(3,
              numericInput(ns("sub_dims"), "Dimensions", value = 30, min = 5, max = 100)
            ),
            column(3,
              numericInput(ns("sub_k"), "k-neighbors", value = 20, min = 5, max = 200)
            )
          ),

          # Resolution sweep option
          fluidRow(
            column(4,
              checkboxInput(ns("sub_resolution_sweep"), "Auto-optimize resolution", FALSE)
            ),
            column(4,
              conditionalPanel(
                condition = paste0("input['", ns("sub_resolution_sweep"), "']"),
                numericInput(ns("sweep_res_min"), "Min resolution", value = 0.1, min = 0.05, max = 1, step = 0.05)
              )
            ),
            column(4,
              conditionalPanel(
                condition = paste0("input['", ns("sub_resolution_sweep"), "']"),
                numericInput(ns("sweep_res_max"), "Max resolution", value = 1.5, min = 0.5, max = 3, step = 0.1)
              )
            )
          ),

          # Advanced options
          fluidRow(
            column(3,
              checkboxInput(ns("recompute_hvg"), "Recompute HVGs on subset", TRUE)
            ),
            column(3,
              checkboxInput(ns("recompute_pca"), "Recompute PCA on subset", TRUE)
            ),
            column(3,
              checkboxInput(ns("recompute_umap"), "Recompute UMAP on subset", TRUE)
            ),
            column(3,
              selectInput(ns("split_by_subcluster"), "Split visualization by",
                choices = NULL)
            )
          ),

          actionButton(ns("run_subclustering"), "Run Sub-clustering",
            class = "btn-primary btn-lg", width = "100%")
        )
      )
    ),

    # Results
    fluidRow(
      column(12,
        uiOutput(ns("sub_results_ui"))
      )
    )
  )
}


# =============================================================================
# SECTION D: subclusteringServer
# =============================================================================

subclusteringServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # -------------------------------------------------------------------
    # 1. Local reactive values
    # -------------------------------------------------------------------
    show_sub_results      <- reactiveVal(FALSE)
    subclustering_results <- reactiveVal(NULL)
    sub_markers           <- reactiveVal(NULL)
    sub_sweep_results     <- reactiveVal(NULL)

    # -------------------------------------------------------------------
    # 2. Results UI renderer
    # -------------------------------------------------------------------
    output$sub_results_ui <- renderUI({
      req(isTRUE(show_sub_results()))
      ns <- session$ns

      # Pre-compute choices
      obj <- values$active_dataset
      vis_reds <- if (!is.null(obj) && inherits(obj, "Seurat")) names(obj@reductions) else c("umap", "pca")
      vis_red_sel <- if ("umap" %in% vis_reds) "umap" else vis_reds[1]

      sub_res <- subclustering_results()
      sub_labels <- if (!is.null(sub_res) && !is.null(sub_res$subcluster_labels)) {
        sort(unique(na.omit(sub_res$subcluster_labels)))
      } else NULL

      color_choices <- if (!is.null(obj) && inherits(obj, "Seurat")) {
        md_cols <- colnames(obj@meta.data)
        cat_cols <- md_cols[vapply(obj@meta.data, function(x) is.factor(x) || is.character(x), logical(1))]
        c("Clusters" = "subclusters", cat_cols)
      } else c("subclusters")

      wellPanel(
        h4("Sub-clustering Results"),
        tabsetPanel(

          # --- Tab 1: Sub-cluster Plot ---
          tabPanel("Sub-cluster Plot",
            fluidRow(
              column(8,
                plotlyOutput(ns("sub_plot"), height = "600px")
              ),
              column(4,
                wellPanel(
                  h5("Visualization Options"),
                  selectInput(ns("sub_vis_reduction"), "Reduction",
                    choices = vis_reds, selected = vis_red_sel),
                  selectInput(ns("sub_color_by"), "Color by",
                    choices = color_choices, selected = "subclusters"),
                  sliderInput(ns("sub_point_size"), "Point size",
                    min = 0.5, max = 6, value = 1.5, step = 0.5),
                  sliderInput(ns("sub_point_alpha"), "Opacity",
                    min = 0.1, max = 1, value = 0.7, step = 0.1),
                  checkboxInput(ns("sub_show_centers"), "Show cluster centers", TRUE)
                ),
                h5("Sub-cluster Summary"),
                downloadButton(ns("dl_sub_summary_table"), "CSV", class = "btn-xs btn-default"),
                DT::dataTableOutput(ns("sub_summary_table"))
              )
            )
          ),

          # --- Tab 2: Comparison ---
          tabPanel("Comparison",
            fluidRow(
              column(6,
                h5("Before Sub-clustering"),
                plotlyOutput(ns("before_subclustering"), height = "400px")
              ),
              column(6,
                h5("After Sub-clustering"),
                plotlyOutput(ns("after_subclustering"), height = "400px")
              )
            )
          ),

          # --- Tab 3: Highlight on Full UMAP ---
          tabPanel("Highlight on Full UMAP",
            fluidRow(
              column(12,
                plotlyOutput(ns("highlight_umap"), height = "600px"),
                helpText("Grey cells are background; colored cells are from the subclustered population.")
              )
            )
          ),

          # --- Tab 4: Sub-cluster Markers ---
          tabPanel("Sub-cluster Markers",
            fluidRow(
              column(4,
                selectInput(ns("subcluster_for_markers"), "Select sub-cluster",
                  choices = sub_labels),
                selectInput(ns("sub_marker_mode"), "Comparison mode",
                  choices = c(
                    "vs parent cluster" = "vs_parent",
                    "vs all subclusters" = "vs_all",
                    "pairwise vs siblings" = "pairwise"
                  ),
                  selected = "vs_parent"),
                actionButton(ns("find_sub_markers"), "Find Markers", class = "btn-info")
              ),
              column(8,
                downloadButton(ns("dl_sub_markers_table"), "CSV", class = "btn-xs btn-default"),
                DT::dataTableOutput(ns("sub_markers_table"))
              )
            )
          ),

          # --- Tab 5: Composition ---
          tabPanel("Composition",
            plotlyOutput(ns("sub_composition_plot"), height = "500px")
          ),

          # --- Tab 6: Split by Batch ---
          tabPanel("Split by Batch",
            plotlyOutput(ns("split_by_batch_plot"), height = "500px")
          ),

          # --- Tab 7: Resolution Sweep ---
          tabPanel("Resolution Sweep",
            uiOutput(ns("sweep_results_ui"))
          )
        )
      )
    })

    # -------------------------------------------------------------------
    # 3. Input updaters
    # -------------------------------------------------------------------

    # Parent cluster column choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      meta_cols <- colnames(obj@meta.data)
      cluster_cols <- meta_cols[vapply(obj@meta.data, function(x) is.factor(x) || is.character(x), logical(1))]
      cluster_cols <- union(cluster_cols, grep("cluster", meta_cols, value = TRUE, ignore.case = TRUE))
      if (length(cluster_cols) > 0) {
        updateSelectInput(session, "parent_cluster_col", choices = sort(cluster_cols))
      }
    })

    # Clusters to subcluster UI (dynamic based on parent column)
    output$clusters_to_subcluster_ui <- renderUI({
      req(values$active_dataset, input$parent_cluster_col)
      if (!input$parent_cluster_col %in% colnames(values$active_dataset@meta.data)) {
        return(NULL)
      }
      clusters <- unique(values$active_dataset@meta.data[[input$parent_cluster_col]])
      clusters <- clusters[!is.na(clusters)]
      selectizeInput(ns("clusters_to_subcluster"), "Clusters to sub-cluster",
        choices = sort(clusters), multiple = TRUE)
    })

    # Split-by choices (batch/condition columns)
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      md_cols <- colnames(obj@meta.data)
      cat_cols <- md_cols[vapply(obj@meta.data, function(x) {
        (is.factor(x) || is.character(x)) && length(unique(x)) > 1 && length(unique(x)) < 50
      }, logical(1))]
      split_choices <- c("(none)" = "", cat_cols)
      updateSelectInput(session, "split_by_subcluster", choices = split_choices)
    })

    # Update sub-cluster choices for markers and vis controls
    observeEvent(subclustering_results(), {
      req(subclustering_results())
      subclusters <- unique(subclustering_results()$subcluster_labels)
      subclusters <- subclusters[!is.na(subclusters)]
      updateSelectInput(session, "subcluster_for_markers", choices = sort(subclusters))

      # Update visualization controls
      obj <- values$active_dataset
      if (!is.null(obj) && inherits(obj, "Seurat")) {
        reds <- names(obj@reductions)
        updateSelectInput(session, "sub_vis_reduction", choices = reds,
          selected = if ("umap" %in% reds) "umap" else reds[1])

        md_cols  <- colnames(obj@meta.data)
        cat_cols <- md_cols[vapply(obj@meta.data, function(x) is.factor(x) || is.character(x), logical(1))]
        color_choices <- c("Clusters" = "subclusters", cat_cols)
        updateSelectInput(session, "sub_color_by", choices = color_choices, selected = "subclusters")
      }
    })

    # -------------------------------------------------------------------
    # 4. Main subclustering execution
    # -------------------------------------------------------------------
    observeEvent(input$run_subclustering, ignoreInit = TRUE, {
      req(values$active_dataset, input$parent_cluster_col, input$clusters_to_subcluster)

      if (length(input$clusters_to_subcluster) == 0) {
        showNotification("Select at least one cluster to sub-cluster.", type = "error")
        return()
      }

      tryCatch({
        withProgress(message = "Running sub-clustering...", value = 0, {
          obj <- values$active_dataset

          # If resolution sweep is enabled, find optimal resolution first
          if (isTRUE(input$sub_resolution_sweep)) {
            incProgress(0.2, detail = "Running resolution sweep...")

            parent_clusters <- input$clusters_to_subcluster
            cells_mask <- obj@meta.data[[input$parent_cluster_col]] %in% parent_clusters
            temp_sub <- tryCatch(
              subset(obj, cells = colnames(obj)[cells_mask]),
              error = function(e) NULL
            )

            if (is.null(temp_sub) || ncol(temp_sub) < 10) {
              showNotification("Too few cells in selected clusters for sweep.", type = "error", duration = 8)
              return()
            }

            # Re-process subset if requested
            if (isTRUE(input$recompute_hvg)) {
              temp_sub <- tryCatch(
                FindVariableFeatures(temp_sub, selection.method = "vst", nfeatures = 2000, verbose = FALSE),
                error = function(e) temp_sub)
            }
            if (isTRUE(input$recompute_pca)) {
              temp_sub <- tryCatch({
                temp_sub <- ScaleData(temp_sub, verbose = FALSE)
                RunPCA(temp_sub, npcs = max(30, input$sub_dims %||% 30L), verbose = FALSE)
              }, error = function(e) temp_sub)
            }

            sub_dims_use <- min(input$sub_dims %||% 30L,
                                tryCatch(ncol(Embeddings(temp_sub, "pca")), error = function(e) 30L))

            temp_sub <- tryCatch({
              set.seed(42)
              FindNeighbors(temp_sub, dims = 1:sub_dims_use, k.param = input$sub_k %||% 20L, verbose = FALSE)
            }, error = function(e) temp_sub)

            algorithm <- switch(input$sub_algorithm %||% "leiden",
              "leiden" = 4L, "louvain" = 1L, "slm" = 3L, 4L)

            sweep <- tryCatch({
              if (exists("run_subclustering_sweep", mode = "function")) {
                run_subclustering_sweep(
                  temp_sub,
                  dims      = sub_dims_use,
                  k_param   = input$sub_k %||% 20L,
                  algorithm = algorithm,
                  res_min   = input$sweep_res_min %||% 0.1,
                  res_max   = input$sweep_res_max %||% 1.5
                )
              } else {
                list(best_resolution = 0.5, sweep_data = NULL)
              }
            }, error = function(e) {
              showNotification(paste("Sweep failed:", conditionMessage(e)), type = "warning", duration = 6)
              list(best_resolution = 0.5, sweep_data = NULL)
            })

            sub_sweep_results(sweep$sweep_data)

            # Use optimal resolution
            params <- reactiveValuesToList(input)
            params$sub_resolution <- sweep$best_resolution

            showNotification(
              paste0("Optimal resolution: ", round(sweep$best_resolution, 2)),
              type = "message", duration = 5)

            incProgress(0.3, detail = "Clustering at optimal resolution...")
            result <- if (exists("run_subclustering_analysis", mode = "function")) {
              run_subclustering_analysis(obj, params)
            } else {
              list(success = FALSE, error = "run_subclustering_analysis not available")
            }

          } else {
            # Standard subclustering without sweep
            result <- if (exists("run_subclustering_analysis", mode = "function")) {
              run_subclustering_analysis(obj, input)
            } else {
              list(success = FALSE, error = "run_subclustering_analysis not available")
            }
          }

          if (isTRUE(result$success)) {
            values$active_dataset <- result$object
            subclustering_results(result)
            show_sub_results(TRUE)

            n_subclusters <- length(unique(na.omit(result$subcluster_labels)))
            showNotification(paste("Sub-clustering complete!", n_subclusters, "sub-clusters found."), type = "message")
          } else {
            showNotification(paste("Sub-clustering failed:", result$error %||% "Unknown error"), type = "error")
          }

          incProgress(1)
        })
      }, error = function(e) {
        showNotification(paste("Sub-clustering error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })

    # -------------------------------------------------------------------
    # 5. Find sub-cluster markers
    # -------------------------------------------------------------------
    observeEvent(input$find_sub_markers, ignoreInit = TRUE, {
      req(values$active_dataset, input$subcluster_for_markers, subclustering_results())
      tryCatch({
        withProgress(message = "Finding sub-cluster markers...", value = 0, {
          # Inherit marker test from parent module if available
          test_use    <- session$userData$parent_marker_test %||% "wilcox"
          marker_mode <- input$sub_marker_mode %||% "vs_all"

          result <- if (exists("find_subcluster_markers", mode = "function")) {
            find_subcluster_markers(
              values$active_dataset,
              input$subcluster_for_markers,
              input$subcluster_name %||% "subclusters",
              mode     = marker_mode,
              test_use = test_use
            )
          } else {
            list(success = FALSE, error = "find_subcluster_markers not available")
          }

          if (isTRUE(result$success)) {
            sub_markers(result$markers)
            showNotification("Sub-cluster markers found.", type = "message")
          } else {
            showNotification(paste("Sub-cluster marker finding failed:", result$error %||% "Unknown error"), type = "error")
          }
        })
      }, error = function(e) {
        showNotification(paste("Sub-cluster marker error:", conditionMessage(e)), type = "error", duration = 10)
      })
    })

    # -------------------------------------------------------------------
    # 6. Result renderers
    # -------------------------------------------------------------------

    output$sub_plot <- renderPlotly({
      req(subclustering_results(), values$active_dataset)
      tryCatch({
        p <- if (exists("create_subclustering_plot", mode = "function")) {
          create_subclustering_plot(values$active_dataset, subclustering_results(), input)
        } else {
          plotly_empty() %>% layout(title = "create_subclustering_plot not available", autosize = TRUE)
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type   = "Sub-cluster embedding colored by sub-cluster identity",
            label         = "Sub-cluster Plot",
            source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Sub-cluster plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    output$sub_summary_table <- DT::renderDataTable({
      req(subclustering_results())
      tryCatch({
        if (exists("create_subcluster_summary_table", mode = "function")) {
          create_subcluster_summary_table(subclustering_results())
        } else {
          labels <- na.omit(subclustering_results()$subcluster_labels)
          data.frame(Subcluster = sort(unique(labels)), Count = as.numeric(table(labels)))
        }
      }, error = function(e) {
        data.frame(Error = conditionMessage(e))
      })
    }, options = list(dom = "t", pageLength = 10))

    output$before_subclustering <- renderPlotly({
      req(subclustering_results(), values$active_dataset)
      tryCatch({
        p <- if (exists("create_before_subclustering_plot", mode = "function")) {
          create_before_subclustering_plot(values$active_dataset, subclustering_results(), input)
        } else {
          plotly_empty() %>% layout(title = "create_before_subclustering_plot not available", autosize = TRUE)
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type   = "Embedding showing original clusters before sub-clustering",
            label         = "Before Subclustering",
            source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Before plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    output$after_subclustering <- renderPlotly({
      req(subclustering_results(), values$active_dataset)
      tryCatch({
        p <- if (exists("create_after_subclustering_plot", mode = "function")) {
          create_after_subclustering_plot(values$active_dataset, subclustering_results(), input)
        } else {
          plotly_empty() %>% layout(title = "create_after_subclustering_plot not available", autosize = TRUE)
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type   = "Embedding showing sub-clusters after sub-clustering",
            label         = "After Subclustering",
            source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("After plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    output$sub_markers_table <- DT::renderDataTable({
      req(sub_markers())
      tbl <- sub_markers()
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Sub-cluster Markers", source_module = "Clustering")
      }
      tbl
    }, options = list(pageLength = 15, scrollX = TRUE))

    # -------------------------------------------------------------------
    # 7. Highlight UMAP renderer
    # -------------------------------------------------------------------
    output$highlight_umap <- renderPlotly({
      req(subclustering_results(), values$active_dataset)
      tryCatch({
        obj <- values$active_dataset
        sub_res <- subclustering_results()

        # Determine reduction for full UMAP
        vis_red <- if ("umap" %in% names(obj@reductions)) "umap"
                   else if ("tsne" %in% names(obj@reductions)) "tsne"
                   else if (exists("best_reduction", mode = "function")) best_reduction(obj)
                   else names(obj@reductions)[1]

        if (is.null(vis_red) || !vis_red %in% names(obj@reductions)) {
          return(plotly_empty() %>% layout(title = "No 2D reduction available for highlight", autosize = TRUE))
        }

        emb <- Embeddings(obj, vis_red)
        md  <- obj@meta.data

        # Identify subclustered cells
        subcluster_col <- sub_res$subcluster_col_name %||% (input$subcluster_name %||% "subclusters")
        if (subcluster_col %in% colnames(md)) {
          is_sub <- !is.na(md[[subcluster_col]])
        } else if (!is.null(sub_res$subclustered_cells)) {
          is_sub <- rownames(md) %in% sub_res$subclustered_cells
        } else {
          is_sub <- rep(FALSE, nrow(md))
        }

        bg_cells  <- !is_sub
        sub_cells <- is_sub

        # Grey background
        p <- plot_ly() %>%
          add_trace(
            x = emb[bg_cells, 1], y = emb[bg_cells, 2],
            type = "scatter", mode = "markers",
            marker = list(color = "#CCCCCC", size = 2, opacity = 0.3),
            name = "Background",
            hoverinfo = "skip"
          )

        # Colored subclustered cells
        if (any(sub_cells)) {
          sub_labels <- if (subcluster_col %in% colnames(md)) {
            as.character(md[[subcluster_col]][sub_cells])
          } else {
            rep("subclustered", sum(sub_cells))
          }

          unique_labels <- sort(unique(sub_labels))
          colors <- if (exists("get_color_palette", mode = "function")) {
            get_color_palette(length(unique_labels))
          } else {
            grDevices::hcl.colors(length(unique_labels), "Dynamic")
          }

          for (i in seq_along(unique_labels)) {
            lbl <- unique_labels[i]
            mask <- sub_labels == lbl
            p <- p %>% add_trace(
              x = emb[which(sub_cells)[mask], 1],
              y = emb[which(sub_cells)[mask], 2],
              type = "scatter", mode = "markers",
              marker = list(color = colors[i], size = 4, opacity = 0.8),
              name = lbl,
              text = paste0("Sub-cluster: ", lbl),
              hoverinfo = "text"
            )
          }
        }

        axis_label <- toupper(vis_red)
        p %>% layout(
          title = "Sub-clustered Cells Highlighted on Full Dataset",
          xaxis = list(title = paste0(axis_label, "_1")),
          yaxis = list(title = paste0(axis_label, "_2")),
          autosize = TRUE
        )
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Highlight UMAP error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    # -------------------------------------------------------------------
    # 8. Split by batch renderer
    # -------------------------------------------------------------------
    output$split_by_batch_plot <- renderPlotly({
      req(subclustering_results(), values$active_dataset)
      split_col <- input$split_by_subcluster
      if (is.null(split_col) || !nzchar(split_col)) {
        return(plotly_empty() %>% layout(
          title = "Select a 'Split visualization by' variable in the panel above.",
          autosize = TRUE))
      }

      tryCatch({
        obj <- values$active_dataset
        md  <- obj@meta.data
        sub_res <- subclustering_results()

        if (!split_col %in% colnames(md)) {
          return(plotly_empty() %>% layout(title = paste("Column", split_col, "not found"), autosize = TRUE))
        }

        # Determine reduction
        vis_red <- if ("umap" %in% names(obj@reductions)) "umap"
                   else if ("tsne" %in% names(obj@reductions)) "tsne"
                   else if (exists("best_reduction", mode = "function")) best_reduction(obj)
                   else names(obj@reductions)[1]

        if (is.null(vis_red)) {
          return(plotly_empty() %>% layout(title = "No 2D reduction available", autosize = TRUE))
        }

        emb <- Embeddings(obj, vis_red)
        subcluster_col <- sub_res$subcluster_col_name %||% (input$subcluster_name %||% "subclusters")

        # Filter to subclustered cells only
        if (subcluster_col %in% colnames(md)) {
          mask <- !is.na(md[[subcluster_col]])
        } else if (!is.null(sub_res$subclustered_cells)) {
          mask <- rownames(md) %in% sub_res$subclustered_cells
        } else {
          mask <- rep(TRUE, nrow(md))
        }

        sub_emb    <- emb[mask, , drop = FALSE]
        sub_md     <- md[mask, , drop = FALSE]
        split_vals <- as.character(sub_md[[split_col]])
        sub_labels <- if (subcluster_col %in% colnames(sub_md)) as.character(sub_md[[subcluster_col]]) else rep("sub", sum(mask))

        groups <- sort(unique(split_vals))
        n_groups <- length(groups)
        if (n_groups == 0 || n_groups > 20) {
          return(plotly_empty() %>% layout(
            title = paste("Split variable has", n_groups, "groups (need 1-20)."),
            autosize = TRUE))
        }

        unique_labels <- sort(unique(sub_labels))
        colors <- if (exists("get_color_palette", mode = "function")) {
          get_color_palette(length(unique_labels))
        } else {
          grDevices::hcl.colors(length(unique_labels), "Dynamic")
        }
        color_map <- setNames(colors, unique_labels)

        n_cols <- min(n_groups, 3)
        n_rows <- ceiling(n_groups / n_cols)

        p <- plotly::subplot(
          lapply(seq_along(groups), function(i) {
            grp <- groups[i]
            grp_mask <- split_vals == grp
            plot_ly(
              x = sub_emb[grp_mask, 1], y = sub_emb[grp_mask, 2],
              color = sub_labels[grp_mask],
              colors = color_map,
              type = "scatter", mode = "markers",
              marker = list(size = 3, opacity = 0.7),
              showlegend = (i == 1),
              text = paste0("Sub-cluster: ", sub_labels[grp_mask]),
              hoverinfo = "text"
            ) %>% layout(annotations = list(
              list(text = grp, x = 0.5, y = 1.05, xref = "paper", yref = "paper",
                   showarrow = FALSE, font = list(size = 12))
            ))
          }),
          nrows = n_rows,
          shareX = TRUE, shareY = TRUE,
          titleX = FALSE, titleY = FALSE
        ) %>% layout(
          title = paste("Sub-clusters split by", split_col),
          autosize = TRUE
        )
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Split plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    # -------------------------------------------------------------------
    # 9. Download handlers
    # -------------------------------------------------------------------

    output$dl_sub_summary_table <- downloadHandler(
      filename = function() paste0("subcluster_summary_", Sys.Date(), ".csv"),
      content = function(file) {
        req(subclustering_results())
        tryCatch({
          tbl <- if (exists("create_subcluster_summary_table", mode = "function")) {
            create_subcluster_summary_table(subclustering_results())
          } else {
            labels <- na.omit(subclustering_results()$subcluster_labels)
            data.frame(Subcluster = sort(unique(labels)), Count = as.numeric(table(labels)))
          }
          write.csv(tbl, file, row.names = FALSE)
        }, error = function(e) write.csv(data.frame(Error = conditionMessage(e)), file, row.names = FALSE))
      }
    )

    output$dl_sub_markers_table <- downloadHandler(
      filename = function() paste0("subcluster_markers_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          req(sub_markers())
          write.csv(sub_markers(), file, row.names = FALSE)
        }, error = function(e) {
          write.csv(data.frame(Error = conditionMessage(e)), file, row.names = FALSE)
        })
      }
    )

    # -------------------------------------------------------------------
    # 10. Composition plot renderer
    # -------------------------------------------------------------------
    output$sub_composition_plot <- renderPlotly({
      req(subclustering_results(), values$active_dataset)
      tryCatch({
        p <- if (exists("create_subcluster_composition_plot", mode = "function")) {
          create_subcluster_composition_plot(values$active_dataset, subclustering_results())
        } else {
          # Fallback: simple stacked bar
          sub_res <- subclustering_results()
          obj <- values$active_dataset
          md  <- obj@meta.data
          subcluster_col <- sub_res$subcluster_col_name %||% (input$subcluster_name %||% "subclusters")
          parent_col     <- input$parent_cluster_col

          if (subcluster_col %in% colnames(md) && parent_col %in% colnames(md)) {
            comp <- as.data.frame(table(
              Subcluster = md[[subcluster_col]],
              Parent     = md[[parent_col]]
            ))
            comp <- comp[comp$Freq > 0, ]
            plot_ly(comp, x = ~Subcluster, y = ~Freq, color = ~Parent, type = "bar") %>%
              layout(barmode = "stack", title = "Sub-cluster Composition", autosize = TRUE)
          } else {
            plotly_empty() %>% layout(title = "Composition data not available", autosize = TRUE)
          }
        }
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Clustering",
            figure_type   = "Stacked bar showing parent cluster composition per sub-cluster",
            label         = "Sub-cluster Composition",
            source_module = "Clustering")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Composition plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    # Resolution sweep results tab
    output$sweep_results_ui <- renderUI({
      sweep_data <- sub_sweep_results()
      ns <- session$ns
      if (is.null(sweep_data)) {
        return(div(style = "padding: 20px; color: #888;",
          p("No resolution sweep data available."),
          p("Enable 'Auto-optimize resolution' and re-run sub-clustering to see sweep results.")))
      }
      tagList(
        plotlyOutput(ns("sweep_quality_plot"), height = "400px"),
        DT::dataTableOutput(ns("sweep_data_table"))
      )
    })

    output$sweep_quality_plot <- renderPlotly({
      req(sub_sweep_results())
      tryCatch({
        if (exists("create_resolution_quality_plot", mode = "function")) {
          create_resolution_quality_plot(sub_sweep_results())
        } else {
          plotly_empty() %>% layout(title = "create_resolution_quality_plot not available", autosize = TRUE)
        }
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Sweep plot error:", conditionMessage(e)), autosize = TRUE)
      })
    })

    output$sweep_data_table <- DT::renderDataTable({
      req(sub_sweep_results())
      DT::datatable(sub_sweep_results(), options = list(dom = "t", pageLength = 20))
    })

  })
}



# ============================================================================
# BACKEND HELPER FUNCTIONS
# Called by graphClusteringServer / subclusteringServer UI/Server sections.
# Appended after line 2785.
# ============================================================================


# =============================================================================
# 1. build_neighbor_graph
#    Constructs the SNN (shared nearest-neighbor) graph that is the foundation
#    of all graph-based clustering algorithms. Validates inputs via
#    resolve_graph_clustering_context() and calls Seurat::FindNeighbors.
# =============================================================================
build_neighbor_graph <- function(seurat_obj, params, graph_prefix = "clustering") {

  # Resolve naming conventions — UI inputs may use either naming scheme
  reduction   <- params$reduction_cluster %||% params$reduction %||% NULL
  dims        <- params$dims_cluster      %||% params$dims      %||% 30L
  k_param     <- params$n_neighbors_cluster %||% params$k_param %||% 20L
  resolution  <- params$resolution        %||% 0.8
  algorithm   <- params$cluster_algorithm %||% params$algorithm %||% "leiden"
  cluster_name <- params$cluster_name     %||% "seurat_clusters"
  prune_snn   <- params$prune_snn        %||% 1/15
  graph_type  <- params$graph_type        %||% "snn"
  random_seed <- params$random_seed       %||% 42L

  # Central context resolution — validates and normalizes everything
  context <- resolve_graph_clustering_context(
    seurat_obj   = seurat_obj,
    reduction    = reduction,
    dims         = as.integer(dims),
    k_param      = as.integer(k_param),
    resolution   = as.numeric(resolution),
    algorithm    = algorithm,
    cluster_name = cluster_name,
    prune_snn    = as.numeric(prune_snn),
    graph_type   = graph_type,
    graph_prefix = graph_prefix,
    random_seed  = as.integer(random_seed)
  )

  n_cells <- ncol(seurat_obj)

  # Inform user about large-dataset processing
  if (n_cells > 100000L) {
    message("[Clustering] Large dataset detected (", format(n_cells, big.mark = ","),
            " cells). Graph construction may take several minutes.")
    notify_clustering(
      paste0("Building neighbor graph for ", format(n_cells, big.mark = ","),
             " cells. This may take a few minutes..."),
      type = "message", duration = 10
    )
  }

  # Build the neighbor graph
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction   = context$reduction,
    dims        = seq_len(context$dims),
    k.param     = context$k_param,
    compute.SNN = TRUE,
    prune.SNN   = context$prune_snn,
    nn.method   = context$nn_method,
    n.trees     = context$n_trees,
    graph.name  = c(context$nn_name, context$graph_name),
    verbose     = FALSE
  )

  # Validate graph was actually created
  available_graphs <- names(seurat_obj@graphs)
  if (!context$graph_name %in% available_graphs) {
    stop("Graph construction failed: '", context$graph_name,
         "' not found in object. Available graphs: ",
         paste(available_graphs, collapse = ", "))
  }

  list(object = seurat_obj, context = context)
}


# =============================================================================
# 2. run_single_resolution_clustering
#    Complete pipeline for clustering at a single resolution: build graph,
#    find clusters, compute quality metrics, record provenance.
# =============================================================================
run_single_resolution_clustering <- function(seurat_obj, params) {
  tryCatch({
    # Build the neighbor graph
    graph_result <- build_neighbor_graph(seurat_obj, params)
    seurat_obj <- graph_result$object
    context    <- graph_result$context

    # Find clusters using the safe wrapper (NEVER call FindClusters directly)
    seurat_obj <- find_clusters_safe(
      object      = seurat_obj,
      graph.name  = context$graph_name,
      resolution  = context$resolution,
      algorithm   = context$algorithm_code,
      random.seed = context$random_seed,
      verbose     = FALSE
    )

    # Extract cluster labels
    cluster_labels <- as.character(Idents(seurat_obj))
    seurat_obj[[context$cluster_name]] <- cluster_labels
    n_clusters <- length(unique(cluster_labels))

    # Compute quality metrics
    quality <- compute_clustering_quality(
      seurat_obj         = seurat_obj,
      reduction          = context$reduction,
      dims               = context$dims,
      cluster_labels     = cluster_labels,
      metrics            = c("silhouette", "modularity"),
      graph_name         = context$graph_name,
      silhouette_max_cells = params$silhouette_max_cells %||% 5000L,
      seed               = context$random_seed
    )

    # Record this run in the object's provenance
    seurat_obj <- record_clustering_run(seurat_obj, context, mode = "single", extra = list(
      n_clusters     = n_clusters,
      silhouette_avg = quality$silhouette_avg,
      modularity     = quality$modularity
    ))

    list(
      success                  = TRUE,
      object                   = seurat_obj,
      cluster_labels           = cluster_labels,
      n_clusters               = n_clusters,
      resolution               = context$resolution,
      algorithm                = context$algorithm_key,
      silhouette_scores        = quality$silhouette_scores,
      silhouette_cluster_labels = quality$silhouette_cluster_labels,
      silhouette_avg           = quality$silhouette_avg,
      modularity               = quality$modularity,
      error                    = NULL
    )
  }, error = function(e) {
    list(
      success                  = FALSE,
      object                   = NULL,
      cluster_labels           = NULL,
      n_clusters               = NULL,
      resolution               = params$resolution %||% 0.8,
      algorithm                = params$cluster_algorithm %||% params$algorithm %||% "leiden",
      silhouette_scores        = NULL,
      silhouette_cluster_labels = NULL,
      silhouette_avg           = NA_real_,
      modularity               = NA_real_,
      error                    = conditionMessage(e)
    )
  })
}


# =============================================================================
# 3. run_resolution_sweep_clustering
#    Builds the neighbor graph once, then tests multiple resolutions to find
#    the best one. Returns a data.frame of results and the final clustering
#    at the optimal resolution.
# =============================================================================
run_resolution_sweep_clustering <- function(seurat_obj, params) {
  tryCatch({
    # Build graph once
    graph_result <- build_neighbor_graph(seurat_obj, params)
    seurat_obj <- graph_result$object
    context    <- graph_result$context

    # Build resolution grid
    resolutions <- build_resolution_grid(
      res_min  = params$res_min  %||% 0.1,
      res_max  = params$res_max  %||% 2.0,
      res_step = params$res_step %||% 0.1
    )

    # Cap silhouette sample size based on dataset size
    n_cells <- ncol(seurat_obj)
    sil_max <- min(
      params$silhouette_max_cells %||% 5000L,
      if (n_cells > 50000L) 3000L else if (n_cells > 20000L) 4000L else 5000L
    )

    # Storage for sweep results
    sweep_rows <- vector("list", length(resolutions))

    for (i in seq_along(resolutions)) {
      res <- resolutions[i]
      col_name <- paste0("seurat_clusters_res.", format_resolution_id(res))

      res_result <- tryCatch({
        tmp_obj <- find_clusters_safe(
          object      = seurat_obj,
          graph.name  = context$graph_name,
          resolution  = res,
          algorithm   = context$algorithm_code,
          random.seed = context$random_seed,
          verbose     = FALSE
        )
        labels <- as.character(Idents(tmp_obj))
        # Store in metadata for later reference
        seurat_obj[[col_name]] <- labels
        n_clust <- length(unique(labels))

        # Compute quality only for non-trivial clusterings
        if (n_clust > 1L) {
          q <- compute_clustering_quality(
            seurat_obj         = seurat_obj,
            reduction          = context$reduction,
            dims               = context$dims,
            cluster_labels     = labels,
            metrics            = c("silhouette", "modularity"),
            graph_name         = context$graph_name,
            silhouette_max_cells = sil_max,
            seed               = context$random_seed
          )
          list(n = n_clust, sil = q$silhouette_avg, mod = q$modularity)
        } else {
          list(n = n_clust, sil = NA_real_, mod = NA_real_)
        }
      }, error = function(e) {
        message("[Clustering] Sweep failed at resolution ", res, ": ", conditionMessage(e))
        list(n = NA_integer_, sil = NA_real_, mod = NA_real_)
      })

      sweep_rows[[i]] <- data.frame(
        Resolution       = res,
        N_Clusters       = res_result$n,
        Silhouette_Score = res_result$sil,
        Modularity       = res_result$mod,
        stringsAsFactors = FALSE
      )
    }

    sweep_data <- do.call(rbind, sweep_rows)

    # Find best resolution: prefer silhouette, fall back to modularity
    valid_sil <- sweep_data[is.finite(sweep_data$Silhouette_Score), ]
    if (nrow(valid_sil) > 0) {
      best_idx <- which.max(valid_sil$Silhouette_Score)
      best_res <- valid_sil$Resolution[best_idx]
    } else {
      valid_mod <- sweep_data[is.finite(sweep_data$Modularity), ]
      if (nrow(valid_mod) > 0) {
        best_idx <- which.max(valid_mod$Modularity)
        best_res <- valid_mod$Resolution[best_idx]
      } else {
        best_res <- resolutions[ceiling(length(resolutions) / 2)]
      }
    }

    # Run final clustering at best resolution
    seurat_obj <- find_clusters_safe(
      object      = seurat_obj,
      graph.name  = context$graph_name,
      resolution  = best_res,
      algorithm   = context$algorithm_code,
      random.seed = context$random_seed,
      verbose     = FALSE
    )

    cluster_labels <- as.character(Idents(seurat_obj))
    seurat_obj[[context$cluster_name]] <- cluster_labels
    n_clusters <- length(unique(cluster_labels))

    # Final quality metrics
    quality <- compute_clustering_quality(
      seurat_obj         = seurat_obj,
      reduction          = context$reduction,
      dims               = context$dims,
      cluster_labels     = cluster_labels,
      metrics            = c("silhouette", "modularity"),
      graph_name         = context$graph_name,
      silhouette_max_cells = sil_max,
      seed               = context$random_seed
    )

    # Record run
    seurat_obj <- record_clustering_run(seurat_obj, context, mode = "sweep", extra = list(
      n_clusters     = n_clusters,
      best_resolution = best_res,
      silhouette_avg = quality$silhouette_avg,
      modularity     = quality$modularity,
      n_resolutions_tested = length(resolutions)
    ))

    # Override context resolution to the best one found
    context$resolution <- best_res

    list(
      success                  = TRUE,
      object                   = seurat_obj,
      cluster_labels           = cluster_labels,
      n_clusters               = n_clusters,
      resolution               = best_res,
      algorithm                = context$algorithm_key,
      silhouette_scores        = quality$silhouette_scores,
      silhouette_cluster_labels = quality$silhouette_cluster_labels,
      silhouette_avg           = quality$silhouette_avg,
      modularity               = quality$modularity,
      sweep_results            = sweep_data,
      error                    = NULL
    )
  }, error = function(e) {
    list(
      success                  = FALSE,
      object                   = NULL,
      cluster_labels           = NULL,
      n_clusters               = NULL,
      resolution               = NULL,
      algorithm                = params$cluster_algorithm %||% params$algorithm %||% "leiden",
      silhouette_scores        = NULL,
      silhouette_cluster_labels = NULL,
      silhouette_avg           = NA_real_,
      modularity               = NA_real_,
      sweep_results            = NULL,
      error                    = conditionMessage(e)
    )
  })
}


# =============================================================================
# 4. run_optimal_resolution_clustering
#    Two-pass optimization: coarse grid search followed by fine refinement
#    around the best coarse point. Supports silhouette, modularity, and
#    Calinski-Harabasz as optimization metrics.
# =============================================================================
run_optimal_resolution_clustering <- function(seurat_obj, params) {
  tryCatch({
    # Build graph once
    graph_result <- build_neighbor_graph(seurat_obj, params)
    seurat_obj <- graph_result$object
    context    <- graph_result$context

    optimization_metric <- tolower(params$optimization_metric %||% "silhouette")
    if (!optimization_metric %in% c("silhouette", "modularity", "calinski")) {
      optimization_metric <- "silhouette"
    }

    target_clusters <- suppressWarnings(as.integer(params$target_clusters %||% 0L))
    if (is.na(target_clusters) || target_clusters < 0L) target_clusters <- 0L

    n_cells <- ncol(seurat_obj)
    sil_max <- min(
      params$silhouette_max_cells %||% 5000L,
      if (n_cells > 50000L) 3000L else if (n_cells > 20000L) 4000L else 5000L
    )

    # Inner function: evaluate a single resolution and return a quality score
    evaluate_resolution <- function(res) {
      tryCatch({
        tmp_obj <- find_clusters_safe(
          object      = seurat_obj,
          graph.name  = context$graph_name,
          resolution  = res,
          algorithm   = context$algorithm_code,
          random.seed = context$random_seed,
          verbose     = FALSE
        )
        labels <- as.character(Idents(tmp_obj))
        n_clust <- length(unique(labels))

        if (n_clust <= 1L) return(list(score = -Inf, n_clusters = n_clust))

        q <- compute_clustering_quality(
          seurat_obj         = seurat_obj,
          reduction          = context$reduction,
          dims               = context$dims,
          cluster_labels     = labels,
          metrics            = optimization_metric,
          graph_name         = context$graph_name,
          silhouette_max_cells = sil_max,
          seed               = context$random_seed
        )

        score <- switch(optimization_metric,
          "silhouette" = q$silhouette_avg %||% NA_real_,
          "modularity" = q$modularity %||% NA_real_,
          "calinski"   = q$calinski_score %||% NA_real_,
          NA_real_
        )
        if (!is.finite(score)) score <- -Inf

        # Apply penalty for deviating from target cluster count
        if (target_clusters > 0L && is.finite(score)) {
          deviation <- abs(n_clust - target_clusters) / max(1L, target_clusters)
          penalty <- 1 - min(1, deviation * 0.5)
          score <- score * penalty
        }

        list(score = score, n_clusters = n_clust)
      }, error = function(e) {
        list(score = -Inf, n_clusters = NA_integer_)
      })
    }

    # Phase 1: Coarse grid search (0.1 to 3.0, ~14 points)
    coarse_grid <- seq(0.1, 3.0, length.out = 14)
    coarse_scores <- vapply(coarse_grid, function(res) {
      evaluate_resolution(res)$score
    }, numeric(1))

    best_coarse_idx <- which.max(coarse_scores)
    best_coarse_res <- coarse_grid[best_coarse_idx]

    # Phase 2: Fine refinement around the best coarse resolution (+-0.2, step 0.05)
    fine_min <- max(0.01, best_coarse_res - 0.2)
    fine_max <- best_coarse_res + 0.2
    fine_grid <- seq(fine_min, fine_max, by = 0.05)
    # Remove points already tested in the coarse grid
    fine_grid <- fine_grid[!fine_grid %in% coarse_grid]

    if (length(fine_grid) > 0) {
      fine_scores <- vapply(fine_grid, function(res) {
        evaluate_resolution(res)$score
      }, numeric(1))

      all_resolutions <- c(coarse_grid, fine_grid)
      all_scores      <- c(coarse_scores, fine_scores)
    } else {
      all_resolutions <- coarse_grid
      all_scores      <- coarse_scores
    }

    # Select overall best
    best_idx <- which.max(all_scores)
    optimal_res <- all_resolutions[best_idx]
    optimal_score <- all_scores[best_idx]

    # Run final clustering at optimal resolution
    seurat_obj <- find_clusters_safe(
      object      = seurat_obj,
      graph.name  = context$graph_name,
      resolution  = optimal_res,
      algorithm   = context$algorithm_code,
      random.seed = context$random_seed,
      verbose     = FALSE
    )

    cluster_labels <- as.character(Idents(seurat_obj))
    seurat_obj[[context$cluster_name]] <- cluster_labels
    n_clusters <- length(unique(cluster_labels))

    # Full quality at optimal resolution
    quality <- compute_clustering_quality(
      seurat_obj         = seurat_obj,
      reduction          = context$reduction,
      dims               = context$dims,
      cluster_labels     = cluster_labels,
      metrics            = c("silhouette", "modularity"),
      graph_name         = context$graph_name,
      silhouette_max_cells = sil_max,
      seed               = context$random_seed
    )

    # Record run
    context$resolution <- optimal_res
    seurat_obj <- record_clustering_run(seurat_obj, context, mode = "optimal", extra = list(
      n_clusters          = n_clusters,
      optimization_metric = optimization_metric,
      optimization_score  = optimal_score,
      silhouette_avg      = quality$silhouette_avg,
      modularity          = quality$modularity,
      n_resolutions_tested = length(all_resolutions)
    ))

    list(
      success                  = TRUE,
      object                   = seurat_obj,
      cluster_labels           = cluster_labels,
      n_clusters               = n_clusters,
      resolution               = optimal_res,
      algorithm                = context$algorithm_key,
      silhouette_scores        = quality$silhouette_scores,
      silhouette_cluster_labels = quality$silhouette_cluster_labels,
      silhouette_avg           = quality$silhouette_avg,
      modularity               = quality$modularity,
      optimization_score       = optimal_score,
      error                    = NULL
    )
  }, error = function(e) {
    list(
      success                  = FALSE,
      object                   = NULL,
      cluster_labels           = NULL,
      n_clusters               = NULL,
      resolution               = NULL,
      algorithm                = params$cluster_algorithm %||% params$algorithm %||% "leiden",
      silhouette_scores        = NULL,
      silhouette_cluster_labels = NULL,
      silhouette_avg           = NA_real_,
      modularity               = NA_real_,
      optimization_score       = NA_real_,
      error                    = conditionMessage(e)
    )
  })
}


# =============================================================================
# 5. run_standard_graph_clustering
#    Convenience wrapper used by subclustering and other modules that need
#    a quick clustering pipeline without the full params machinery.
# =============================================================================
run_standard_graph_clustering <- function(seurat_obj,
                                          reduction          = NULL,
                                          dims               = 30L,
                                          k_param            = 20L,
                                          resolution         = 0.8,
                                          algorithm          = "leiden",
                                          cluster_name       = "seurat_clusters",
                                          prune_snn          = 1/15,
                                          graph_type         = "snn",
                                          graph_prefix       = "clustering",
                                          compute_quality    = TRUE,
                                          quality_metrics    = c("silhouette", "modularity"),
                                          silhouette_max_cells = 5000L,
                                          random_seed        = 42L) {

  # Build a params list mimicking the UI input structure
  params <- list(
    reduction    = reduction,
    dims         = dims,
    k_param      = k_param,
    resolution   = resolution,
    algorithm    = algorithm,
    cluster_name = cluster_name,
    prune_snn    = prune_snn,
    graph_type   = graph_type,
    random_seed  = random_seed,
    silhouette_max_cells = silhouette_max_cells
  )

  # Build neighbor graph
  graph_result <- build_neighbor_graph(seurat_obj, params, graph_prefix = graph_prefix)
  seurat_obj <- graph_result$object
  context    <- graph_result$context

  # Find clusters via safe wrapper
  seurat_obj <- find_clusters_safe(
    object      = seurat_obj,
    graph.name  = context$graph_name,
    resolution  = context$resolution,
    algorithm   = context$algorithm_code,
    random.seed = context$random_seed,
    verbose     = FALSE
  )

  cluster_labels <- as.character(Idents(seurat_obj))
  seurat_obj[[context$cluster_name]] <- cluster_labels

  # Compute quality metrics if requested
  quality <- default_clustering_quality()
  if (isTRUE(compute_quality) && length(unique(cluster_labels)) > 1L) {
    quality <- compute_clustering_quality(
      seurat_obj         = seurat_obj,
      reduction          = context$reduction,
      dims               = context$dims,
      cluster_labels     = cluster_labels,
      metrics            = quality_metrics,
      graph_name         = context$graph_name,
      silhouette_max_cells = silhouette_max_cells,
      seed               = context$random_seed
    )
  }

  list(
    object         = seurat_obj,
    cluster_labels = cluster_labels,
    context        = context,
    quality        = quality
  )
}


# =============================================================================
# 6. compute_clustering_quality
#    Computes one or more quality metrics for a given clustering. All metrics
#    are computed inside tryCatch so that partial failures do not block the
#    rest of the pipeline.
# =============================================================================
compute_clustering_quality <- function(seurat_obj,
                                       reduction,
                                       dims,
                                       cluster_labels,
                                       metrics            = c("silhouette", "modularity"),
                                       graph_name         = NULL,
                                       silhouette_max_cells = 5000L,
                                       modularity_max_cells = 50000L,
                                       seed               = 42L) {
  tryCatch({
    result <- default_clustering_quality()

    # Validate inputs
    if (is.null(seurat_obj) || is.null(cluster_labels)) return(result)
    cluster_labels <- as.character(cluster_labels)
    n_cells <- ncol(seurat_obj)
    if (n_cells < 3L || length(unique(cluster_labels)) < 2L) return(result)

    # Get embedding matrix and clamp dimensions
    embedding <- tryCatch(
      Embeddings(seurat_obj, reduction),
      error = function(e) NULL
    )
    if (is.null(embedding)) return(result)

    available_dims <- ncol(embedding)
    dims_use <- min(as.integer(dims %||% 30L), available_dims)
    dims_use <- max(2L, dims_use)
    embedding <- embedding[, seq_len(dims_use), drop = FALSE]

    metrics <- tolower(metrics)

    # --- Silhouette ---
    if ("silhouette" %in% metrics) {
      tryCatch({
        sil_max <- max(10L, as.integer(silhouette_max_cells %||% 5000L))
        if (n_cells > sil_max) {
          idx <- sample_cells_stratified(cluster_labels, max_cells = sil_max, seed = seed)
          sil_embedding <- embedding[idx, , drop = FALSE]
          sil_labels    <- cluster_labels[idx]
        } else {
          sil_embedding <- embedding
          sil_labels    <- cluster_labels
          idx           <- seq_len(n_cells)
        }

        # Need at least 2 clusters in sample
        if (length(unique(sil_labels)) >= 2L) {
          label_int <- as.integer(as.factor(sil_labels))
          d <- dist(sil_embedding)
          sil <- cluster::silhouette(label_int, d)
          result$silhouette_scores        <- sil[, 3]
          result$silhouette_cluster_labels <- sil_labels
          result$silhouette_avg           <- mean(sil[, 3], na.rm = TRUE)
          if (!is.finite(result$silhouette_avg)) result$silhouette_avg <- NA_real_
        }
      }, error = function(e) {
        message("[Clustering] Silhouette computation failed: ", conditionMessage(e))
      })
    }

    # --- Calinski-Harabasz ---
    if ("calinski" %in% metrics) {
      tryCatch({
        if (requireNamespace("fpc", quietly = TRUE) &&
            exists("calinhara", where = asNamespace("fpc"), mode = "function")) {
          label_int <- as.integer(as.factor(cluster_labels))
          result$calinski_score <- fpc::calinhara(embedding, label_int)
          if (!is.finite(result$calinski_score)) result$calinski_score <- NA_real_
        } else {
          message("[Clustering] fpc::calinhara not available; skipping Calinski-Harabasz.")
        }
      }, error = function(e) {
        message("[Clustering] Calinski-Harabasz computation failed: ", conditionMessage(e))
      })
    }

    # --- Modularity ---
    if ("modularity" %in% metrics) {
      tryCatch({
        mod_max <- as.integer(modularity_max_cells %||% 25000L)
        if (n_cells > mod_max) {
          message("[Clustering] Skipping modularity: dataset (", n_cells,
                  " cells) exceeds modularity_max_cells (", mod_max, ").")
        } else {
          # Resolve graph name
          gname <- graph_name %||% NULL
          if (is.null(gname) || !gname %in% names(seurat_obj@graphs)) {
            # Try to find any SNN graph
            avail_graphs <- names(seurat_obj@graphs)
            snn_graphs <- avail_graphs[grepl("snn", avail_graphs, ignore.case = TRUE)]
            gname <- if (length(snn_graphs) > 0) snn_graphs[1] else NULL
          }

          if (!is.null(gname)) {
            graph_mat <- seurat_obj@graphs[[gname]]
            # Convert to sparse matrix then to igraph
            sparse_mat <- methods::as(graph_mat, "dgCMatrix")
            ig <- igraph::graph_from_adjacency_matrix(
              sparse_mat, mode = "undirected", weighted = TRUE, diag = FALSE
            )
            membership <- as.integer(as.factor(cluster_labels))
            result$modularity <- igraph::modularity(ig, membership)
            if (!is.finite(result$modularity)) result$modularity <- NA_real_
          } else {
            message("[Clustering] No SNN graph found for modularity computation.")
          }
        }
      }, error = function(e) {
        message("[Clustering] Modularity computation failed: ", conditionMessage(e))
      })
    }

    result
  }, error = function(e) {
    message("[Clustering] Quality computation failed: ", conditionMessage(e))
    default_clustering_quality()
  })
}


# =============================================================================
# 7. find_cluster_specific_markers
#    Finds DE markers for a single cluster versus all other cells.
# =============================================================================
find_cluster_specific_markers <- function(seurat_obj, params) {
  tryCatch({
    cluster_col <- params$cluster_name %||% "seurat_clusters"
    marker_cluster <- params$marker_cluster
    marker_test <- params$marker_test %||% "wilcox"

    # Validate cluster column exists
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop("Cluster column '", cluster_col, "' not found in metadata. ",
           "Available: ", paste(head(colnames(seurat_obj@meta.data), 10), collapse = ", "))
    }

    Idents(seurat_obj) <- cluster_col

    # Handle Seurat v5 split layers
    if (exists("JoinLayers", mode = "function")) {
      tryCatch({
        layers <- Layers(seurat_obj)
        if (length(layers) > 1 && any(grepl("^data\\.", layers))) {
          seurat_obj <- JoinLayers(seurat_obj)
        }
      }, error = function(e) {
        message("[Clustering] JoinLayers skipped: ", conditionMessage(e))
      })
    }

    # Ensure data layer exists
    tryCatch({
      assay <- DefaultAssay(seurat_obj)
      layers <- Layers(seurat_obj, assay = assay)
      if (!"data" %in% layers && "counts" %in% layers) {
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
      }
    }, error = function(e) NULL)

    markers <- FindMarkers(
      seurat_obj,
      ident.1    = marker_cluster,
      test.use   = marker_test,
      verbose    = FALSE
    )

    # Ensure numeric columns (guard against factors)
    numeric_cols <- c("avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
    for (col in intersect(numeric_cols, colnames(markers))) {
      if (is.factor(markers[[col]])) {
        markers[[col]] <- as.numeric(as.character(markers[[col]]))
      }
    }

    list(success = TRUE, markers = markers, error = NULL)
  }, error = function(e) {
    list(success = FALSE, markers = NULL, error = conditionMessage(e))
  })
}


# =============================================================================
# 8. find_all_cluster_markers
#    Finds DE markers for all clusters. Includes progressive fallback with
#    lower thresholds and a final AverageExpression-based fallback.
# =============================================================================
find_all_cluster_markers <- function(seurat_obj, params) {
  tryCatch({
    cluster_col <- params$cluster_name %||% "seurat_clusters"
    marker_test <- params$marker_test  %||% "wilcox"
    top_n       <- as.integer(params$top_n_markers %||% 10L)

    # Validate cluster column
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop("Cluster column '", cluster_col, "' not found in metadata.")
    }

    Idents(seurat_obj) <- cluster_col
    clusters <- levels(Idents(seurat_obj))

    # Downsample to max 200 cells per cluster for speed
    max_per_cluster <- 200L
    cell_counts <- table(Idents(seurat_obj))
    needs_downsample <- any(cell_counts > max_per_cluster)

    if (needs_downsample) {
      set.seed(params$random_seed %||% 42L)
      cells_keep <- unlist(lapply(names(cell_counts), function(cl) {
        cl_cells <- WhichCells(seurat_obj, idents = cl)
        if (length(cl_cells) > max_per_cluster) {
          sample(cl_cells, max_per_cluster)
        } else {
          cl_cells
        }
      }))
      work_obj <- subset(seurat_obj, cells = cells_keep)
    } else {
      work_obj <- seurat_obj
    }

    # Handle Seurat v5 split layers
    if (exists("JoinLayers", mode = "function")) {
      tryCatch({
        layers <- Layers(work_obj)
        if (length(layers) > 1 && any(grepl("^data\\.", layers))) {
          work_obj <- JoinLayers(work_obj)
        }
      }, error = function(e) {
        message("[Clustering] JoinLayers skipped: ", conditionMessage(e))
      })
    }

    # Ensure data layer
    tryCatch({
      assay <- DefaultAssay(work_obj)
      layers <- Layers(work_obj, assay = assay)
      if (!"data" %in% layers && "counts" %in% layers) {
        work_obj <- NormalizeData(work_obj, verbose = FALSE)
      }
    }, error = function(e) NULL)

    # Progressive attempts with lower thresholds
    attempts <- list(
      list(min.pct = 0.25, logfc.threshold = 0.25),
      list(min.pct = 0.10, logfc.threshold = 0.10),
      list(min.pct = 0.05, logfc.threshold = 0.05)
    )

    markers <- NULL

    for (attempt in attempts) {
      markers <- tryCatch({
        FindAllMarkers(
          work_obj,
          test.use         = marker_test,
          min.pct          = attempt$min.pct,
          logfc.threshold  = attempt$logfc.threshold,
          only.pos         = TRUE,
          verbose          = FALSE
        )
      }, error = function(e) {
        message("[Clustering] FindAllMarkers attempt failed (min.pct=",
                attempt$min.pct, "): ", conditionMessage(e))
        NULL
      })

      if (!is.null(markers) && nrow(markers) > 0) break
    }

    # Final fallback: AverageExpression fold-change method
    if (is.null(markers) || nrow(markers) == 0) {
      message("[Clustering] All FindAllMarkers attempts failed. Using AverageExpression fallback.")
      markers <- tryCatch({
        avg_expr <- AverageExpression(work_obj, assays = DefaultAssay(work_obj), verbose = FALSE)
        expr_mat <- avg_expr[[DefaultAssay(work_obj)]]

        if (is.null(expr_mat) || ncol(expr_mat) < 2) {
          stop("AverageExpression returned insufficient data.")
        }

        # For each cluster, compute log2 fold change vs. mean of all others
        all_markers <- lapply(colnames(expr_mat), function(cl) {
          this_cl <- expr_mat[, cl]
          other_cls <- rowMeans(expr_mat[, setdiff(colnames(expr_mat), cl), drop = FALSE])

          # Pseudo-count to avoid log(0)
          lfc <- log2((this_cl + 1e-9) / (other_cls + 1e-9))
          top_genes <- head(order(lfc, decreasing = TRUE), top_n)

          data.frame(
            gene         = rownames(expr_mat)[top_genes],
            cluster      = cl,
            avg_log2FC   = lfc[top_genes],
            p_val        = NA_real_,
            p_val_adj    = NA_real_,
            pct.1        = NA_real_,
            pct.2        = NA_real_,
            stringsAsFactors = FALSE
          )
        })
        do.call(rbind, all_markers)
      }, error = function(e) {
        message("[Clustering] AverageExpression fallback failed: ", conditionMessage(e))
        NULL
      })
    }

    if (is.null(markers) || nrow(markers) == 0) {
      return(list(success = FALSE, markers = NULL, error = "No markers found with any method."))
    }

    # Ensure numeric columns
    numeric_cols <- c("avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
    for (col in intersect(numeric_cols, colnames(markers))) {
      if (is.factor(markers[[col]])) {
        markers[[col]] <- as.numeric(as.character(markers[[col]]))
      }
    }

    # Return top N per cluster
    if ("cluster" %in% colnames(markers)) {
      markers <- markers %>%
        dplyr::group_by(cluster) %>%
        dplyr::slice_head(n = top_n) %>%
        dplyr::ungroup() %>%
        as.data.frame()
    }

    list(success = TRUE, markers = markers, error = NULL)
  }, error = function(e) {
    list(success = FALSE, markers = NULL, error = conditionMessage(e))
  })
}


# =============================================================================
# 9. prepare_subclustering_subset
#    Prepares a subset Seurat object for subclustering by optionally
#    recomputing HVGs, PCA, and UMAP on the subset.
# =============================================================================
prepare_subclustering_subset <- function(sub_obj,
                                         recompute_hvg  = TRUE,
                                         recompute_pca  = TRUE,
                                         recompute_umap = TRUE,
                                         sub_dims       = 20L,
                                         random_seed    = 42L) {

  # Use prepare_for_de if available (from core_robustness.R)
  if (exists("prepare_for_de", mode = "function")) {
    sub_obj <- tryCatch(
      prepare_for_de(sub_obj),
      error = function(e) {
        message("[Clustering] prepare_for_de failed: ", conditionMessage(e))
        sub_obj
      }
    )
  }

  reduction_use <- "pca"
  dims_use <- as.integer(sub_dims %||% 20L)
  if (is.na(dims_use) || dims_use < 2L) dims_use <- 20L

  # Recompute variable features on the subset
  if (isTRUE(recompute_hvg)) {
    sub_obj <- tryCatch({
      FindVariableFeatures(sub_obj, selection.method = "vst",
                           nfeatures = 2000, verbose = FALSE)
    }, error = function(e) {
      message("[Clustering] FindVariableFeatures on subset failed: ", conditionMessage(e))
      sub_obj
    })
  }

  # Recompute PCA on the subset
  if (isTRUE(recompute_pca)) {
    sub_obj <- tryCatch({
      sub_obj <- ScaleData(sub_obj, verbose = FALSE)
      n_features <- length(VariableFeatures(sub_obj))
      npcs <- min(dims_use, max(5L, n_features - 1L), ncol(sub_obj) - 1L)
      sub_obj <- RunPCA(sub_obj, npcs = npcs, verbose = FALSE, seed.use = random_seed)
      # Clamp dims to what is actually available
      dims_use <<- min(dims_use, ncol(Embeddings(sub_obj, "pca")))
      sub_obj
    }, error = function(e) {
      message("[Clustering] PCA on subset failed: ", conditionMessage(e))
      # Fall back to existing reduction if available
      reduction_use <<- select_clustering_reduction(sub_obj, allow_2d = FALSE) %||% "pca"
      if (reduction_use %in% names(sub_obj@reductions)) {
        dims_use <<- min(dims_use, ncol(Embeddings(sub_obj, reduction_use)))
      }
      sub_obj
    })
  } else {
    # Use existing reduction
    reduction_use <- select_clustering_reduction(sub_obj, allow_2d = FALSE) %||% "pca"
    if (reduction_use %in% names(sub_obj@reductions)) {
      dims_use <- min(dims_use, ncol(Embeddings(sub_obj, reduction_use)))
    }
  }

  # Recompute UMAP for visualization
  if (isTRUE(recompute_umap)) {
    sub_obj <- tryCatch({
      RunUMAP(sub_obj, dims = seq_len(dims_use), reduction = reduction_use,
              verbose = FALSE, seed.use = random_seed)
    }, error = function(e) {
      message("[Clustering] UMAP on subset failed: ", conditionMessage(e))
      sub_obj
    })
  }

  list(object = sub_obj, reduction = reduction_use, dims = dims_use)
}


# =============================================================================
# 10. run_subclustering_analysis
#     Complete subclustering pipeline: subsets to selected clusters,
#     recomputes embeddings, clusters, and maps sub-clusters back to the
#     parent object with "parent.sub" naming.
# =============================================================================
run_subclustering_analysis <- function(seurat_obj, params) {
  tryCatch({
    cluster_col      <- params$cluster_name     %||% params$sub_cluster_col %||% "seurat_clusters"
    selected_clusters <- params$selected_clusters %||% params$sub_selected_clusters
    recompute_hvg    <- params$sub_recompute_hvg  %||% TRUE
    recompute_pca    <- params$sub_recompute_pca  %||% TRUE
    recompute_umap   <- params$sub_recompute_umap %||% TRUE
    sub_dims         <- params$sub_dims           %||% 20L
    sub_resolution   <- params$sub_resolution     %||% 0.5
    sub_algorithm    <- params$sub_algorithm      %||% params$cluster_algorithm %||% "leiden"
    sub_k_param      <- params$sub_k_param        %||% 15L
    random_seed      <- params$random_seed        %||% 42L

    # Validate cluster column
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop("Cluster column '", cluster_col, "' not found in metadata.")
    }

    # Validate selected clusters
    available <- unique(as.character(seurat_obj@meta.data[[cluster_col]]))
    if (is.null(selected_clusters) || length(selected_clusters) == 0) {
      stop("No clusters selected for subclustering.")
    }
    missing <- setdiff(as.character(selected_clusters), available)
    if (length(missing) > 0) {
      notify_clustering(
        paste0("Clusters not found (skipped): ", paste(missing, collapse = ", ")),
        type = "warning", duration = 8
      )
    }
    valid_clusters <- intersect(as.character(selected_clusters), available)
    if (length(valid_clusters) == 0) {
      stop("None of the selected clusters exist in the data.")
    }

    # Subset to selected clusters
    cells_keep <- WhichCells(seurat_obj, expression = !!sym(cluster_col) %in% valid_clusters)
    if (length(cells_keep) < 10L) {
      stop("Too few cells (", length(cells_keep), ") in selected clusters for subclustering.")
    }
    sub_obj <- subset(seurat_obj, cells = cells_keep)

    # Prepare the subset (recompute HVG/PCA/UMAP)
    prep <- prepare_subclustering_subset(
      sub_obj        = sub_obj,
      recompute_hvg  = isTRUE(recompute_hvg),
      recompute_pca  = isTRUE(recompute_pca),
      recompute_umap = isTRUE(recompute_umap),
      sub_dims       = as.integer(sub_dims),
      random_seed    = as.integer(random_seed)
    )
    sub_obj <- prep$object

    # Run clustering on the subset
    sub_cluster_result <- run_standard_graph_clustering(
      seurat_obj         = sub_obj,
      reduction          = prep$reduction,
      dims               = prep$dims,
      k_param            = as.integer(sub_k_param),
      resolution         = as.numeric(sub_resolution),
      algorithm          = sub_algorithm,
      cluster_name       = "sub_clusters",
      graph_prefix       = "subcluster",
      compute_quality    = TRUE,
      silhouette_max_cells = params$silhouette_max_cells %||% 5000L,
      random_seed        = as.integer(random_seed)
    )

    sub_obj <- sub_cluster_result$object
    sub_labels_raw <- as.character(sub_cluster_result$cluster_labels)

    # Map sub-clusters: "parent.sub" format
    parent_labels <- as.character(sub_obj@meta.data[[cluster_col]])
    subcluster_labels <- paste0(parent_labels, ".", sub_labels_raw)
    sub_obj[["subcluster_id"]] <- subcluster_labels
    n_subclusters <- length(unique(subcluster_labels))

    # Map back to the full object
    full_subcluster <- rep(NA_character_, ncol(seurat_obj))
    names(full_subcluster) <- colnames(seurat_obj)
    full_subcluster[colnames(sub_obj)] <- subcluster_labels
    seurat_obj[["subcluster_id"]] <- full_subcluster

    # Store provenance
    seurat_obj@misc$subclustering_provenance <- list(
      parent_column     = cluster_col,
      parent_clusters   = valid_clusters,
      n_cells_subset    = ncol(sub_obj),
      n_subclusters     = n_subclusters,
      sub_resolution    = sub_resolution,
      sub_algorithm     = sub_algorithm,
      sub_dims          = prep$dims,
      sub_reduction     = prep$reduction,
      recomputed_hvg    = isTRUE(recompute_hvg),
      recomputed_pca    = isTRUE(recompute_pca),
      recomputed_umap   = isTRUE(recompute_umap),
      timestamp         = as.character(Sys.time())
    )

    list(
      success           = TRUE,
      object            = seurat_obj,
      subcluster_labels = subcluster_labels,
      sub_object        = sub_obj,
      parent_clusters   = valid_clusters,
      n_subclusters     = n_subclusters,
      error             = NULL
    )
  }, error = function(e) {
    list(
      success           = FALSE,
      object            = NULL,
      subcluster_labels = NULL,
      sub_object        = NULL,
      parent_clusters   = NULL,
      n_subclusters     = NULL,
      error             = conditionMessage(e)
    )
  })
}


# =============================================================================
# 11. run_subclustering_sweep
#     Tests multiple resolutions for subclustering. Builds the graph once,
#     then evaluates each resolution.
# =============================================================================
run_subclustering_sweep <- function(sub_obj,
                                    reduction = "pca",
                                    dims      = 20L,
                                    k_param   = 15L,
                                    algorithm = "leiden",
                                    res_min   = 0.1,
                                    res_max   = 1.5,
                                    res_step  = 0.1) {
  tryCatch({
    # Build graph once
    params <- list(
      reduction   = reduction,
      dims        = dims,
      k_param     = k_param,
      algorithm   = algorithm,
      cluster_name = "sweep_temp",
      random_seed = 42L
    )
    graph_result <- build_neighbor_graph(sub_obj, params, graph_prefix = "subsweep")
    sub_obj <- graph_result$object
    context <- graph_result$context

    resolutions <- build_resolution_grid(res_min, res_max, res_step)

    sweep_rows <- vector("list", length(resolutions))

    for (i in seq_along(resolutions)) {
      res <- resolutions[i]
      res_result <- tryCatch({
        tmp_obj <- find_clusters_safe(
          object      = sub_obj,
          graph.name  = context$graph_name,
          resolution  = res,
          algorithm   = context$algorithm_code,
          random.seed = context$random_seed,
          verbose     = FALSE
        )
        labels <- as.character(Idents(tmp_obj))
        n_clust <- length(unique(labels))

        sil <- NA_real_
        if (n_clust > 1L) {
          sil <- compute_silhouette_safe(
            Embeddings(sub_obj, context$reduction)[, seq_len(context$dims), drop = FALSE],
            labels,
            max_cells = min(5000L, ncol(sub_obj)),
            seed = context$random_seed
          )
        }

        list(n = n_clust, sil = sil)
      }, error = function(e) {
        list(n = NA_integer_, sil = NA_real_)
      })

      sweep_rows[[i]] <- data.frame(
        Resolution       = res,
        N_Clusters       = res_result$n,
        Silhouette_Score = res_result$sil,
        stringsAsFactors = FALSE
      )
    }

    sweep_data <- do.call(rbind, sweep_rows)

    # Best resolution by silhouette
    valid <- sweep_data[is.finite(sweep_data$Silhouette_Score), ]
    if (nrow(valid) > 0) {
      best_res <- valid$Resolution[which.max(valid$Silhouette_Score)]
    } else {
      best_res <- resolutions[ceiling(length(resolutions) / 2)]
    }

    list(sweep_data = sweep_data, best_resolution = best_res)
  }, error = function(e) {
    message("[Clustering] Subclustering sweep failed: ", conditionMessage(e))
    list(sweep_data = NULL, best_resolution = NULL)
  })
}


# =============================================================================
# 12. find_subcluster_markers
#     Finds markers for a specific sub-cluster using different comparison
#     modes: "vs_parent" (siblings), "vs_all" (all subclustered cells), or
#     "pairwise" (each sibling individually).
# =============================================================================
find_subcluster_markers <- function(seurat_obj,
                                    subcluster_id,
                                    subcluster_column = "subcluster_id",
                                    mode              = "vs_parent",
                                    test_use          = "wilcox") {
  tryCatch({
    # Validate column
    if (!subcluster_column %in% colnames(seurat_obj@meta.data)) {
      stop("Subcluster column '", subcluster_column, "' not found in metadata.")
    }

    all_labels <- as.character(seurat_obj@meta.data[[subcluster_column]])
    if (!subcluster_id %in% all_labels) {
      stop("Subcluster '", subcluster_id, "' not found in column '", subcluster_column, "'.")
    }

    # Filter to non-NA cells only
    valid_cells <- !is.na(all_labels)
    if (sum(valid_cells) < 10) {
      stop("Too few cells with subcluster labels for marker analysis.")
    }

    work_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[valid_cells])
    Idents(work_obj) <- subcluster_column

    # Handle Seurat v5 split layers
    if (exists("JoinLayers", mode = "function")) {
      tryCatch({
        layers <- Layers(work_obj)
        if (length(layers) > 1 && any(grepl("^data\\.", layers))) {
          work_obj <- JoinLayers(work_obj)
        }
      }, error = function(e) {
        message("[Clustering] JoinLayers skipped: ", conditionMessage(e))
      })
    }

    # Ensure data layer
    tryCatch({
      assay <- DefaultAssay(work_obj)
      layers <- Layers(work_obj, assay = assay)
      if (!"data" %in% layers && "counts" %in% layers) {
        work_obj <- NormalizeData(work_obj, verbose = FALSE)
      }
    }, error = function(e) NULL)

    # Parse parent cluster from "parent.sub" format
    parent <- sub("\\..*$", "", subcluster_id)
    all_subclusters <- unique(as.character(Idents(work_obj)))
    siblings <- all_subclusters[grepl(
      paste0("^", gsub("([.\\\\^$|?*+(){}\\[\\]])", "\\\\\\1", parent), "\\."),
      all_subclusters
    )]
    siblings <- setdiff(siblings, subcluster_id)

    mode <- tolower(mode %||% "vs_parent")

    if (mode == "vs_parent") {
      # Compare to sibling sub-clusters sharing same parent
      if (length(siblings) == 0) {
        # No siblings: fall back to vs_all
        markers <- FindMarkers(work_obj, ident.1 = subcluster_id,
                               test.use = test_use, verbose = FALSE)
      } else {
        markers <- FindMarkers(work_obj, ident.1 = subcluster_id, ident.2 = siblings,
                               test.use = test_use, verbose = FALSE)
      }
    } else if (mode == "vs_all") {
      # Compare to all other subclustered cells
      markers <- FindMarkers(work_obj, ident.1 = subcluster_id,
                             test.use = test_use, verbose = FALSE)
    } else if (mode == "pairwise") {
      # Compare to each sibling individually
      if (length(siblings) == 0) {
        markers <- FindMarkers(work_obj, ident.1 = subcluster_id,
                               test.use = test_use, verbose = FALSE)
        markers$comparison <- "vs_all"
      } else {
        markers_list <- lapply(siblings, function(sib) {
          tryCatch({
            m <- FindMarkers(work_obj, ident.1 = subcluster_id, ident.2 = sib,
                             test.use = test_use, verbose = FALSE)
            m$gene <- rownames(m)
            m$comparison <- paste0("vs_", sib)
            m
          }, error = function(e) NULL)
        })
        markers_list <- markers_list[!vapply(markers_list, is.null, logical(1))]
        if (length(markers_list) == 0) {
          stop("All pairwise comparisons failed.")
        }
        markers <- do.call(rbind, markers_list)
      }
    } else {
      stop("Unknown marker mode: '", mode, "'. Use 'vs_parent', 'vs_all', or 'pairwise'.")
    }

    # Ensure numeric columns
    numeric_cols <- c("avg_log2FC", "p_val", "p_val_adj", "pct.1", "pct.2")
    for (col in intersect(numeric_cols, colnames(markers))) {
      if (is.factor(markers[[col]])) {
        markers[[col]] <- as.numeric(as.character(markers[[col]]))
      }
    }

    list(success = TRUE, markers = markers, error = NULL)
  }, error = function(e) {
    list(success = FALSE, markers = NULL, error = conditionMessage(e))
  })
}


# =============================================================================
# 13. create_enhanced_cluster_plot
#     Interactive plotly scatter plot of cells colored by cluster assignment
#     with optional cluster center labels.
# =============================================================================
create_enhanced_cluster_plot <- function(seurat_obj, clustering_results, input_params) {
  tryCatch({
    # Determine reduction for visualization
    vis_reduction <- input_params$vis_reduction %||% NULL
    if (is.null(vis_reduction) || !vis_reduction %in% names(seurat_obj@reductions)) {
      vis_reduction <- best_reduction(seurat_obj) %||% "umap"
    }
    if (!vis_reduction %in% names(seurat_obj@reductions)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No reduction available for visualization",
                              autosize = TRUE))
    }

    embeddings <- Embeddings(seurat_obj, vis_reduction)
    dim1 <- embeddings[, 1]
    dim2 <- embeddings[, 2]

    # Determine coloring
    color_by <- input_params$color_by_cluster %||% "cluster"
    if (color_by == "cluster" || is.null(color_by)) {
      color_labels <- as.character(
        clustering_results$cluster_labels %||% Idents(seurat_obj)
      )
    } else if (color_by %in% colnames(seurat_obj@meta.data)) {
      color_labels <- as.character(seurat_obj@meta.data[[color_by]])
    } else {
      color_labels <- as.character(Idents(seurat_obj))
    }

    # Ensure lengths match
    if (length(color_labels) != nrow(embeddings)) {
      color_labels <- as.character(Idents(seurat_obj))
    }

    unique_labels <- sort(unique(color_labels))
    colors <- get_color_palette(length(unique_labels))
    color_map <- setNames(colors, unique_labels)

    # Plot parameters — UI names these point_size_cluster / point_alpha_cluster
    point_size  <- as.numeric(input_params$point_size_cluster %||% input_params$point_size %||% 1.5)
    point_alpha <- as.numeric(input_params$point_alpha_cluster %||% input_params$point_alpha %||% 0.7)
    show_centers <- isTRUE(input_params$show_cluster_centers %||% input_params$show_centers %||% FALSE)

    # Axis labels
    red_upper <- toupper(vis_reduction)
    x_label <- paste0(red_upper, "_1")
    y_label <- paste0(red_upper, "_2")

    # Build plotly scatter
    cell_names <- colnames(seurat_obj)

    p <- plotly::plot_ly() %>%
      plotly::layout(
        xaxis = list(title = x_label),
        yaxis = list(title = y_label),
        title = paste0("Clusters (", vis_reduction, ")"),
        legend = list(title = list(text = "Cluster")),
        autosize = TRUE
      )

    # Add traces per cluster for proper legend grouping
    for (cl in unique_labels) {
      mask <- color_labels == cl
      p <- p %>% plotly::add_trace(
        x = dim1[mask],
        y = dim2[mask],
        type = "scatter",
        mode = "markers",
        marker = list(
          size    = point_size,
          color   = color_map[cl],
          opacity = point_alpha
        ),
        name = cl,
        text = paste0("Cell: ", cell_names[mask],
                       "<br>", x_label, ": ", round(dim1[mask], 2),
                       "<br>", y_label, ": ", round(dim2[mask], 2),
                       "<br>Cluster: ", cl),
        hovertemplate = paste0(
          "Cell: %{text}<extra>", cl, "</extra>"
        ),
        hoverinfo = "text"
      )
    }

    # Optional cluster centers
    if (isTRUE(show_centers)) {
      center_df <- data.frame(
        x     = tapply(dim1, color_labels, mean, na.rm = TRUE),
        y     = tapply(dim2, color_labels, mean, na.rm = TRUE),
        label = names(tapply(dim1, color_labels, mean, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )

      p <- p %>%
        plotly::add_markers(
          data = center_df,
          x = ~x, y = ~y,
          marker = list(size = 10, color = "black", symbol = "diamond",
                        line = list(color = "white", width = 1)),
          showlegend = FALSE,
          hoverinfo = "skip"
        ) %>%
        plotly::add_text(
          data = center_df,
          x = ~x, y = ~y, text = ~label,
          textfont = list(size = 11, color = "black"),
          showlegend = FALSE,
          hoverinfo = "skip"
        )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 14. create_cluster_summary_table
#     Summary data.frame with cluster ID, cell count, and percentage.
# =============================================================================
create_cluster_summary_table <- function(clustering_results) {
  tryCatch({
    labels <- as.character(clustering_results$cluster_labels)
    if (is.null(labels) || length(labels) == 0) {
      return(data.frame(Cluster = character(), Cell_Count = integer(),
                        Percentage = numeric(), stringsAsFactors = FALSE))
    }

    counts <- sort(table(labels), decreasing = TRUE)
    total  <- sum(counts)

    data.frame(
      Cluster    = names(counts),
      Cell_Count = as.integer(counts),
      Percentage = round(as.numeric(counts) / total * 100, 1),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(Cluster = character(), Cell_Count = integer(),
               Percentage = numeric(), stringsAsFactors = FALSE)
  })
}


# =============================================================================
# 15. create_cluster_distribution_plot
#     Bar chart of cluster sizes.
# =============================================================================
create_cluster_distribution_plot <- function(clustering_results) {
  tryCatch({
    labels <- as.character(clustering_results$cluster_labels)
    if (is.null(labels) || length(labels) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No clustering data available", autosize = TRUE))
    }

    counts <- table(labels)
    # Sort by cluster label (numeric if possible)
    sorted_names <- tryCatch(
      names(counts)[order(as.numeric(names(counts)))],
      warning = function(w) sort(names(counts)),
      error   = function(e) sort(names(counts))
    )
    counts <- counts[sorted_names]

    n_clusters <- length(counts)
    colors <- get_color_palette(n_clusters)

    plotly::plot_ly(
      x      = names(counts),
      y      = as.integer(counts),
      type   = "bar",
      marker = list(color = colors),
      text   = paste0(names(counts), ": ", as.integer(counts), " cells"),
      hoverinfo = "text"
    ) %>%
      plotly::layout(
        title  = "Cluster Size Distribution",
        xaxis  = list(title = "Cluster", type = "category"),
        yaxis  = list(title = "Number of Cells"),
        autosize = TRUE
      )
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Distribution plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 16. create_quality_metrics_table
#     Summary table of clustering quality metrics.
# =============================================================================
create_quality_metrics_table <- function(clustering_results) {
  tryCatch({
    n_clusters <- clustering_results$n_clusters %||% NA
    sil_avg    <- clustering_results$silhouette_avg %||% NA_real_
    mod_val    <- clustering_results$modularity %||% NA_real_
    res_val    <- clustering_results$resolution %||% NA_real_

    # Format nicely
    format_val <- function(x, digits = 4) {
      if (is.null(x) || !is.finite(x)) return("N/A")
      formatC(x, format = "f", digits = digits)
    }

    data.frame(
      Metric = c("Number of Clusters",
                  "Average Silhouette Score",
                  "Modularity",
                  "Resolution"),
      Value  = c(as.character(n_clusters %||% "N/A"),
                 format_val(sil_avg),
                 format_val(mod_val),
                 format_val(res_val, 3)),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(Metric = "Error", Value = conditionMessage(e),
               stringsAsFactors = FALSE)
  })
}


# =============================================================================
# 17. create_silhouette_plot
#     Histogram of silhouette scores per cluster with a vertical average line.
# =============================================================================
create_silhouette_plot <- function(clustering_results) {
  tryCatch({
    sil_scores <- clustering_results$silhouette_scores
    sil_labels <- clustering_results$silhouette_cluster_labels
    sil_avg    <- clustering_results$silhouette_avg %||% NA_real_

    if (is.null(sil_scores) || length(sil_scores) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No silhouette data available", autosize = TRUE))
    }

    # Guard against length mismatch
    if (is.null(sil_labels) || length(sil_labels) != length(sil_scores)) {
      # Fall back to a single histogram
      sil_labels <- rep("All", length(sil_scores))
    }

    unique_labels <- sort(unique(sil_labels))
    colors <- get_color_palette(length(unique_labels))
    color_map <- setNames(colors, unique_labels)

    p <- plotly::plot_ly()

    for (cl in unique_labels) {
      mask <- sil_labels == cl
      p <- p %>% plotly::add_histogram(
        x          = sil_scores[mask],
        name       = cl,
        marker     = list(color = color_map[cl]),
        opacity    = 0.6,
        nbinsx     = 30
      )
    }

    p <- p %>% plotly::layout(
      barmode  = "overlay",
      title    = "Silhouette Score Distribution",
      xaxis    = list(title = "Silhouette Score", range = c(-1, 1)),
      yaxis    = list(title = "Count"),
      autosize = TRUE
    )

    # Add average line
    if (is.finite(sil_avg)) {
      p <- p %>% plotly::layout(
        shapes = list(
          list(
            type = "line",
            x0 = sil_avg, x1 = sil_avg,
            y0 = 0, y1 = 1,
            yref = "paper",
            line = list(color = "red", dash = "dash", width = 2)
          )
        ),
        annotations = list(
          list(
            x = sil_avg, y = 1, yref = "paper",
            text = paste0("Avg: ", round(sil_avg, 3)),
            showarrow = FALSE,
            font = list(color = "red", size = 12),
            yanchor = "bottom"
          )
        )
      )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Silhouette plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 18. create_resolution_quality_plot
#     Dual-line chart showing Silhouette (blue) and Modularity (red) across
#     the resolution sweep.
# =============================================================================
create_resolution_quality_plot <- function(sweep_results) {
  tryCatch({
    if (is.null(sweep_results) || !is.data.frame(sweep_results) ||
        nrow(sweep_results) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No sweep results available", autosize = TRUE))
    }

    p <- plotly::plot_ly(data = sweep_results) %>%
      plotly::layout(
        title  = "Resolution vs Quality Metrics",
        xaxis  = list(title = "Resolution"),
        yaxis  = list(title = "Score"),
        legend = list(x = 0.02, y = 0.98),
        autosize = TRUE
      )

    # Silhouette line (blue)
    if ("Silhouette_Score" %in% colnames(sweep_results)) {
      valid_sil <- is.finite(sweep_results$Silhouette_Score)
      if (any(valid_sil)) {
        p <- p %>% plotly::add_trace(
          x    = sweep_results$Resolution[valid_sil],
          y    = sweep_results$Silhouette_Score[valid_sil],
          type = "scatter",
          mode = "lines+markers",
          name = "Silhouette",
          line = list(color = "blue", width = 2),
          marker = list(color = "blue", size = 6)
        )
      }
    }

    # Modularity line (red)
    if ("Modularity" %in% colnames(sweep_results)) {
      valid_mod <- is.finite(sweep_results$Modularity)
      if (any(valid_mod)) {
        p <- p %>% plotly::add_trace(
          x    = sweep_results$Resolution[valid_mod],
          y    = sweep_results$Modularity[valid_mod],
          type = "scatter",
          mode = "lines+markers",
          name = "Modularity",
          line = list(color = "red", width = 2),
          marker = list(color = "red", size = 6)
        )
      }
    }

    # Add N_Clusters as secondary info on hover
    if ("N_Clusters" %in% colnames(sweep_results)) {
      p <- p %>% plotly::add_trace(
        x    = sweep_results$Resolution,
        y    = rep(0, nrow(sweep_results)),
        type = "scatter",
        mode = "markers",
        marker = list(size = 0.1, opacity = 0),
        text = paste0("k=", sweep_results$N_Clusters),
        hoverinfo = "text+x",
        name = "N Clusters",
        showlegend = FALSE
      )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Quality plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 19. create_subclustering_plot
#     Interactive scatter of subclustered cells colored by sub-cluster ID.
# =============================================================================
create_subclustering_plot <- function(seurat_obj, subclustering_results, input_params) {
  tryCatch({
    vis_reduction <- input_params$sub_vis_reduction %||% NULL
    if (is.null(vis_reduction) || !vis_reduction %in% names(seurat_obj@reductions)) {
      vis_reduction <- best_reduction(seurat_obj) %||% "umap"
    }
    if (!vis_reduction %in% names(seurat_obj@reductions)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No reduction available", autosize = TRUE))
    }

    point_size  <- as.numeric(input_params$sub_point_size  %||% 1.5)
    point_alpha <- as.numeric(input_params$sub_point_alpha %||% 0.7)
    show_centers <- isTRUE(input_params$sub_show_centers %||% TRUE)
    color_by    <- input_params$sub_color_by %||% "subcluster"

    # Get subcluster labels and filter to non-NA
    subcluster_col <- "subcluster_id"
    if (!subcluster_col %in% colnames(seurat_obj@meta.data)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustering results found",
                              autosize = TRUE))
    }

    all_labels <- seurat_obj@meta.data[[subcluster_col]]
    valid_mask <- !is.na(all_labels)
    if (sum(valid_mask) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustered cells", autosize = TRUE))
    }

    embeddings <- Embeddings(seurat_obj, vis_reduction)
    dim1 <- embeddings[valid_mask, 1]
    dim2 <- embeddings[valid_mask, 2]
    cell_names <- colnames(seurat_obj)[valid_mask]

    if (color_by == "subcluster" || is.null(color_by)) {
      color_labels <- as.character(all_labels[valid_mask])
    } else if (color_by %in% colnames(seurat_obj@meta.data)) {
      color_labels <- as.character(seurat_obj@meta.data[[color_by]][valid_mask])
    } else {
      color_labels <- as.character(all_labels[valid_mask])
    }

    unique_labels <- sort(unique(color_labels))
    colors <- get_color_palette(length(unique_labels))
    color_map <- setNames(colors, unique_labels)

    red_upper <- toupper(vis_reduction)
    x_label <- paste0(red_upper, "_1")
    y_label <- paste0(red_upper, "_2")

    p <- plotly::plot_ly() %>%
      plotly::layout(
        title  = paste0("Sub-clusters (", vis_reduction, ")"),
        xaxis  = list(title = x_label),
        yaxis  = list(title = y_label),
        legend = list(title = list(text = "Sub-cluster")),
        autosize = TRUE
      )

    for (cl in unique_labels) {
      mask <- color_labels == cl
      p <- p %>% plotly::add_trace(
        x = dim1[mask],
        y = dim2[mask],
        type = "scatter",
        mode = "markers",
        marker = list(size = point_size, color = color_map[cl],
                      opacity = point_alpha),
        name = cl,
        text = paste0("Cell: ", cell_names[mask],
                       "<br>", x_label, ": ", round(dim1[mask], 2),
                       "<br>", y_label, ": ", round(dim2[mask], 2),
                       "<br>Sub-cluster: ", color_labels[mask]),
        hoverinfo = "text"
      )
    }

    if (isTRUE(show_centers)) {
      center_df <- data.frame(
        x     = tapply(dim1, color_labels, mean, na.rm = TRUE),
        y     = tapply(dim2, color_labels, mean, na.rm = TRUE),
        label = names(tapply(dim1, color_labels, mean, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )

      p <- p %>%
        plotly::add_markers(
          data = center_df,
          x = ~x, y = ~y,
          marker = list(size = 10, color = "black", symbol = "diamond",
                        line = list(color = "white", width = 1)),
          showlegend = FALSE, hoverinfo = "skip"
        ) %>%
        plotly::add_text(
          data = center_df,
          x = ~x, y = ~y, text = ~label,
          textfont = list(size = 10, color = "black"),
          showlegend = FALSE, hoverinfo = "skip"
        )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Subclustering plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 20. create_subcluster_summary_table
#     Summary table with sub-cluster ID, cell count, and percentage.
# =============================================================================
create_subcluster_summary_table <- function(subclustering_results) {
  tryCatch({
    labels <- as.character(subclustering_results$subcluster_labels)
    if (is.null(labels) || length(labels) == 0) {
      return(data.frame(SubCluster = character(), Cell_Count = integer(),
                        Percentage = numeric(), stringsAsFactors = FALSE))
    }

    counts <- sort(table(labels), decreasing = TRUE)
    total  <- sum(counts)

    data.frame(
      SubCluster = names(counts),
      Cell_Count = as.integer(counts),
      Percentage = round(as.numeric(counts) / total * 100, 1),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(SubCluster = character(), Cell_Count = integer(),
               Percentage = numeric(), stringsAsFactors = FALSE)
  })
}


# =============================================================================
# 21. create_before_subclustering_plot
#     Shows selected cells colored by their original parent cluster assignment
#     (before subclustering).
# =============================================================================
create_before_subclustering_plot <- function(seurat_obj, subclustering_results,
                                             input_params) {
  tryCatch({
    vis_reduction <- input_params$sub_vis_reduction %||% NULL
    if (is.null(vis_reduction) || !vis_reduction %in% names(seurat_obj@reductions)) {
      vis_reduction <- best_reduction(seurat_obj) %||% "umap"
    }
    if (!vis_reduction %in% names(seurat_obj@reductions)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No reduction available", autosize = TRUE))
    }

    # Identify subclustered cells
    subcluster_col <- "subcluster_id"
    if (!subcluster_col %in% colnames(seurat_obj@meta.data)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustering results found",
                              autosize = TRUE))
    }

    all_labels <- seurat_obj@meta.data[[subcluster_col]]
    valid_mask <- !is.na(all_labels)
    if (sum(valid_mask) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustered cells", autosize = TRUE))
    }

    embeddings <- Embeddings(seurat_obj, vis_reduction)
    dim1 <- embeddings[valid_mask, 1]
    dim2 <- embeddings[valid_mask, 2]
    cell_names <- colnames(seurat_obj)[valid_mask]

    # Get parent cluster labels
    provenance <- seurat_obj@misc$subclustering_provenance
    parent_col <- provenance$parent_column %||% "seurat_clusters"

    if (parent_col %in% colnames(seurat_obj@meta.data)) {
      parent_labels <- as.character(seurat_obj@meta.data[[parent_col]][valid_mask])
    } else {
      # Extract from subcluster ID (parent.sub format)
      parent_labels <- sub("\\..*$", "", as.character(all_labels[valid_mask]))
    }

    unique_parents <- sort(unique(parent_labels))
    colors <- get_color_palette(length(unique_parents))
    color_map <- setNames(colors, unique_parents)

    red_upper <- toupper(vis_reduction)
    x_label <- paste0(red_upper, "_1")
    y_label <- paste0(red_upper, "_2")

    p <- plotly::plot_ly() %>%
      plotly::layout(
        title  = "Before Subclustering (Parent Clusters)",
        xaxis  = list(title = x_label),
        yaxis  = list(title = y_label),
        legend = list(title = list(text = "Parent Cluster")),
        autosize = TRUE
      )

    point_size  <- as.numeric(input_params$sub_point_size  %||% 1.5)
    point_alpha <- as.numeric(input_params$sub_point_alpha %||% 0.7)

    for (cl in unique_parents) {
      mask <- parent_labels == cl
      p <- p %>% plotly::add_trace(
        x = dim1[mask],
        y = dim2[mask],
        type = "scatter",
        mode = "markers",
        marker = list(size = point_size, color = color_map[cl],
                      opacity = point_alpha),
        name = cl,
        text = paste0("Cell: ", cell_names[mask],
                       "<br>Parent: ", cl),
        hoverinfo = "text"
      )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Before plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 22. create_after_subclustering_plot
#     Shows selected cells colored by their sub-cluster assignment
#     (after subclustering).
# =============================================================================
create_after_subclustering_plot <- function(seurat_obj, subclustering_results,
                                            input_params) {
  tryCatch({
    vis_reduction <- input_params$sub_vis_reduction %||% NULL
    if (is.null(vis_reduction) || !vis_reduction %in% names(seurat_obj@reductions)) {
      vis_reduction <- best_reduction(seurat_obj) %||% "umap"
    }
    if (!vis_reduction %in% names(seurat_obj@reductions)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No reduction available", autosize = TRUE))
    }

    subcluster_col <- "subcluster_id"
    if (!subcluster_col %in% colnames(seurat_obj@meta.data)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustering results found",
                              autosize = TRUE))
    }

    all_labels <- seurat_obj@meta.data[[subcluster_col]]
    valid_mask <- !is.na(all_labels)
    if (sum(valid_mask) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustered cells", autosize = TRUE))
    }

    embeddings <- Embeddings(seurat_obj, vis_reduction)
    dim1 <- embeddings[valid_mask, 1]
    dim2 <- embeddings[valid_mask, 2]
    cell_names <- colnames(seurat_obj)[valid_mask]
    sub_labels <- as.character(all_labels[valid_mask])

    unique_labels <- sort(unique(sub_labels))
    colors <- get_color_palette(length(unique_labels))
    color_map <- setNames(colors, unique_labels)

    red_upper <- toupper(vis_reduction)
    x_label <- paste0(red_upper, "_1")
    y_label <- paste0(red_upper, "_2")

    point_size  <- as.numeric(input_params$sub_point_size  %||% 1.5)
    point_alpha <- as.numeric(input_params$sub_point_alpha %||% 0.7)
    show_centers <- isTRUE(input_params$sub_show_centers %||% TRUE)

    p <- plotly::plot_ly() %>%
      plotly::layout(
        title  = "After Subclustering (Sub-cluster Assignment)",
        xaxis  = list(title = x_label),
        yaxis  = list(title = y_label),
        legend = list(title = list(text = "Sub-cluster")),
        autosize = TRUE
      )

    for (cl in unique_labels) {
      mask <- sub_labels == cl
      p <- p %>% plotly::add_trace(
        x = dim1[mask],
        y = dim2[mask],
        type = "scatter",
        mode = "markers",
        marker = list(size = point_size, color = color_map[cl],
                      opacity = point_alpha),
        name = cl,
        text = paste0("Cell: ", cell_names[mask],
                       "<br>Sub-cluster: ", cl),
        hoverinfo = "text"
      )
    }

    if (isTRUE(show_centers)) {
      center_df <- data.frame(
        x     = tapply(dim1, sub_labels, mean, na.rm = TRUE),
        y     = tapply(dim2, sub_labels, mean, na.rm = TRUE),
        label = names(tapply(dim1, sub_labels, mean, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )

      p <- p %>%
        plotly::add_markers(
          data = center_df,
          x = ~x, y = ~y,
          marker = list(size = 10, color = "black", symbol = "diamond",
                        line = list(color = "white", width = 1)),
          showlegend = FALSE, hoverinfo = "skip"
        ) %>%
        plotly::add_text(
          data = center_df,
          x = ~x, y = ~y, text = ~label,
          textfont = list(size = 10, color = "black"),
          showlegend = FALSE, hoverinfo = "skip"
        )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("After plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}


# =============================================================================
# 23. create_subcluster_composition_plot
#     Stacked bar chart showing parent cluster composition per sub-cluster.
# =============================================================================
create_subcluster_composition_plot <- function(seurat_obj, subclustering_results) {
  tryCatch({
    subcluster_col <- "subcluster_id"
    if (!subcluster_col %in% colnames(seurat_obj@meta.data)) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustering results found",
                              autosize = TRUE))
    }

    all_labels <- seurat_obj@meta.data[[subcluster_col]]
    valid_mask <- !is.na(all_labels)
    if (sum(valid_mask) == 0) {
      return(plotly::plotly_empty() %>%
               plotly::layout(title = "No subclustered cells", autosize = TRUE))
    }

    sub_labels <- as.character(all_labels[valid_mask])

    # Get parent labels
    provenance <- seurat_obj@misc$subclustering_provenance
    parent_col <- provenance$parent_column %||% "seurat_clusters"

    if (parent_col %in% colnames(seurat_obj@meta.data)) {
      parent_labels <- as.character(seurat_obj@meta.data[[parent_col]][valid_mask])
    } else {
      parent_labels <- sub("\\..*$", "", sub_labels)
    }

    # Cross-tabulation
    comp_table <- table(SubCluster = sub_labels, Parent = parent_labels)
    comp_df <- as.data.frame(comp_table, stringsAsFactors = FALSE)
    colnames(comp_df) <- c("SubCluster", "Parent", "Count")

    # Sort sub-clusters
    sub_order <- tryCatch(
      sort(unique(comp_df$SubCluster)[order(
        as.numeric(gsub("^.*\\.", "", unique(comp_df$SubCluster)))
      )]),
      error = function(e) sort(unique(comp_df$SubCluster))
    )

    unique_parents <- sort(unique(comp_df$Parent))
    colors <- get_color_palette(length(unique_parents))
    color_map <- setNames(colors, unique_parents)

    p <- plotly::plot_ly() %>%
      plotly::layout(
        barmode = "stack",
        title   = "Sub-cluster Composition by Parent Cluster",
        xaxis   = list(title = "Sub-cluster", categoryorder = "array",
                        categoryarray = sub_order),
        yaxis   = list(title = "Number of Cells"),
        legend  = list(title = list(text = "Parent Cluster")),
        autosize = TRUE
      )

    for (parent in unique_parents) {
      parent_data <- comp_df[comp_df$Parent == parent, ]
      p <- p %>% plotly::add_trace(
        x      = parent_data$SubCluster,
        y      = parent_data$Count,
        type   = "bar",
        name   = parent,
        marker = list(color = color_map[parent]),
        text   = paste0(parent, ": ", parent_data$Count),
        hoverinfo = "text"
      )
    }

    p
  }, error = function(e) {
    plotly::plotly_empty() %>%
      plotly::layout(title = paste("Composition plot error:", conditionMessage(e)),
                     autosize = TRUE)
  })
}
