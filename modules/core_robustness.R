# Core robustness, provenance, and validation helpers
# modules/core_robustness.R

suppressPackageStartupMessages({
  library(Seurat)
  library(shiny)
})

if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}

.state_get <- function(values, name, default = NULL) {
  val <- NULL
  if (inherits(values, "reactivevalues")) {
    val <- shiny::isolate(values[[name]])
  } else if (is.environment(values) || is.list(values)) {
    val <- values[[name]]
  }
  if (is.null(val)) default else val
}

.state_set <- function(values, name, value) {
  values[[name]] <- value
  invisible(value)
}

.safe_now_iso <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
}

.json_encode <- function(x, pretty = FALSE) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(x, auto_unbox = TRUE, null = "null", pretty = pretty)
  } else {
    paste(capture.output(str(x)), collapse = " ")
  }
}

.safe_layer_matrix <- function(obj, assay = NULL, preferred_layers = c("data", "counts")) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
  if (is.null(assay)) assay <- DefaultAssay(obj)
  if (is.null(assay) || !assay %in% names(obj@assays)) assay <- names(obj@assays)[1]
  if (is.null(assay)) return(NULL)

  assay_obj <- obj[[assay]]
  if (inherits(assay_obj, "Assay5")) {
    layers <- SeuratObject::Layers(assay_obj)
    for (layer_name in preferred_layers) {
      if (layer_name %in% layers) {
        return(SeuratObject::LayerData(assay_obj, layer = layer_name))
      }
      split_layers <- grep(paste0("^", layer_name, "\\."), layers, value = TRUE)
      if (length(split_layers) > 0) {
        joined <- tryCatch(SeuratObject::JoinLayers(assay_obj), error = function(e) NULL)
        if (!is.null(joined) && layer_name %in% SeuratObject::Layers(joined)) {
          return(SeuratObject::LayerData(joined, layer = layer_name))
        }
      }
    }
    return(NULL)
  }

  for (layer_name in preferred_layers) {
    mat <- tryCatch(SeuratObject::GetAssayData(obj, assay = assay, layer = layer_name), error = function(e) NULL)
    if (!is.null(mat)) return(mat)
  }
  NULL
}

# =============================================================================
# Canonical utility functions (sourced once, used by all modules)
# =============================================================================

# Seurat v4/v5 compatible data accessor (handles split layers in v5)
safe_get_assay_data <- function(obj, layer_name = "data", assay = NULL, default = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(obj)
  if (!assay %in% names(obj@assays)) assay <- names(obj@assays)[1]

  assay_obj <- obj[[assay]]

  if (inherits(assay_obj, "Assay5")) {
    lyrs <- Layers(assay_obj)
    if (layer_name %in% lyrs) {
      return(LayerData(assay_obj, layer = layer_name))
    }
    split_lyrs <- grep(paste0("^", layer_name, "\\."), lyrs, value = TRUE)
    if (length(split_lyrs) > 0) {
      joined_assay <- JoinLayers(assay_obj)
      return(LayerData(joined_assay, layer = layer_name))
    }
    if (!is.null(default)) return(default)
    stop("Layer '", layer_name, "' not found in assay '", assay,
         "'. Available: ", paste(lyrs, collapse = ", "))
  }

  tryCatch(
    GetAssayData(obj, layer = layer_name),
    error = function(e) {
      if (!is.null(default)) return(default)
      stop(e)
    }
  )
}

# Prepare a Seurat object for DE / data access (joins split layers, normalizes)
prepare_for_de <- function(obj, assay = NULL) {
  if (is.null(assay)) assay <- DefaultAssay(obj)
  if (!assay %in% names(obj@assays)) assay <- names(obj@assays)[1]
  assay_obj <- obj[[assay]]
  if (inherits(assay_obj, "Assay5")) {
    lyrs <- Layers(assay_obj)
    if (any(grepl("\\.", lyrs))) {
      obj[[assay]] <- JoinLayers(obj[[assay]])
      lyrs <- Layers(obj[[assay]])
    }
    if (!"data" %in% lyrs && "counts" %in% lyrs) {
      obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
    }
  }
  obj
}

compute_seurat_signature <- function(obj, sample_cells = 256L, sample_features = 256L) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NA_character_)

  sample_cells <- max(1L, as.integer(sample_cells))
  sample_features <- max(1L, as.integer(sample_features))

  # Lightweight deterministic signature to detect stale/reactive mismatch.
  dims <- c(n_features = nrow(obj), n_cells = ncol(obj))
  assays <- sort(names(obj@assays))
  reductions <- sort(names(obj@reductions))
  meta_cols <- sort(colnames(obj@meta.data))
  default_assay <- DefaultAssay(obj) %||% assays[1] %||% ""

  cell_ids <- sort(colnames(obj))
  feature_ids <- sort(rownames(obj))
  cell_sample <- utils::head(cell_ids, sample_cells)
  feature_sample <- utils::head(feature_ids, sample_features)

  mat <- .safe_layer_matrix(obj, assay = default_assay, preferred_layers = c("data", "counts"))
  expr_digest <- "none"
  if (!is.null(mat) && length(cell_sample) > 0 && length(feature_sample) > 0) {
    sub <- tryCatch(mat[feature_sample, cell_sample, drop = FALSE], error = function(e) NULL)
    if (!is.null(sub)) {
      vals <- as.numeric(sub)
      vals <- vals[is.finite(vals)]
      if (length(vals) == 0) vals <- 0
      sparse_summary <- list(
        sum = as.numeric(sum(vals)),
        mean = as.numeric(mean(vals)),
        nonzero = as.numeric(sum(vals != 0)),
        q = as.numeric(stats::quantile(vals, probs = c(0, 0.5, 0.9, 0.99, 1), na.rm = TRUE))
      )
      expr_digest <- .json_encode(sparse_summary, pretty = FALSE)
    }
  }

  signature_input <- list(
    dims = dims,
    assays = assays,
    reductions = reductions,
    meta_cols = meta_cols,
    default_assay = default_assay,
    cells = cell_sample,
    features = feature_sample,
    expr = expr_digest
  )
  payload <- .json_encode(signature_input, pretty = FALSE)

  if (requireNamespace("digest", quietly = TRUE)) {
    digest::digest(payload, algo = "xxhash64")
  } else {
    paste0("sig_", as.integer(sum(utf8ToInt(payload)) %% 1000000007L))
  }
}

append_audit_log <- function(values, module, action, level = "info", details = list()) {
  if (is.null(values$analysis_audit)) {
    values$analysis_audit <- data.frame(
      timestamp = character(0),
      module = character(0),
      action = character(0),
      level = character(0),
      details = character(0),
      stringsAsFactors = FALSE
    )
  }

  row <- data.frame(
    timestamp = .safe_now_iso(),
    module = as.character(module %||% ""),
    action = as.character(action %||% ""),
    level = as.character(level %||% "info"),
    details = as.character(.json_encode(details, pretty = FALSE)),
    stringsAsFactors = FALSE
  )
  values$analysis_audit <- rbind(values$analysis_audit, row)
  invisible(row)
}

start_job <- function(values, module, action, meta = list()) {
  job_id <- paste0("job_", format(as.numeric(Sys.time()) * 1000, scientific = FALSE, digits = 15), "_", sample.int(9999, 1))
  row <- list(
    job_id = job_id,
    module = module,
    action = action,
    status = "running",
    progress = 0,
    message = "started",
    started_at = .safe_now_iso(),
    finished_at = NA_character_,
    meta = meta
  )
  values$job_queue <- c(values$job_queue %||% list(), list(row))
  append_audit_log(values, module = module, action = paste0("job_start:", action), details = list(job_id = job_id, meta = meta))
  job_id
}

update_job <- function(values, job_id, progress = NULL, message = NULL, status = NULL) {
  jobs <- values$job_queue %||% list()
  if (length(jobs) == 0) return(invisible(NULL))

  idx <- which(vapply(jobs, function(x) identical(x$job_id, job_id), logical(1)))
  if (length(idx) == 0) return(invisible(NULL))
  i <- idx[1]

  if (!is.null(progress)) jobs[[i]]$progress <- max(0, min(100, as.numeric(progress)))
  if (!is.null(message)) jobs[[i]]$message <- as.character(message)
  if (!is.null(status)) jobs[[i]]$status <- as.character(status)
  values$job_queue <- jobs
  invisible(jobs[[i]])
}

finish_job <- function(values, job_id, status = "done", message = NULL, details = list()) {
  jobs <- values$job_queue %||% list()
  if (length(jobs) == 0) return(invisible(NULL))
  idx <- which(vapply(jobs, function(x) identical(x$job_id, job_id), logical(1)))
  if (length(idx) == 0) return(invisible(NULL))
  i <- idx[1]

  jobs[[i]]$status <- as.character(status %||% "done")
  jobs[[i]]$progress <- if (jobs[[i]]$status %in% c("done", "success")) 100 else jobs[[i]]$progress
  jobs[[i]]$message <- as.character(message %||% jobs[[i]]$message)
  jobs[[i]]$finished_at <- .safe_now_iso()
  values$job_queue <- jobs

  append_audit_log(
    values,
    module = jobs[[i]]$module,
    action = paste0("job_finish:", jobs[[i]]$action),
    level = if (jobs[[i]]$status %in% c("failed", "error")) "error" else "info",
    details = c(list(job_id = job_id, status = jobs[[i]]$status), details)
  )
  invisible(jobs[[i]])
}

job_queue_table <- function(values, limit = 50L) {
  jobs <- values$job_queue %||% list()
  if (length(jobs) == 0) {
    return(data.frame(
      job_id = character(0),
      module = character(0),
      action = character(0),
      status = character(0),
      progress = numeric(0),
      message = character(0),
      started_at = character(0),
      finished_at = character(0),
      stringsAsFactors = FALSE
    ))
  }

  df <- do.call(rbind, lapply(jobs, function(j) {
    data.frame(
      job_id = as.character(j$job_id %||% ""),
      module = as.character(j$module %||% ""),
      action = as.character(j$action %||% ""),
      status = as.character(j$status %||% ""),
      progress = as.numeric(j$progress %||% 0),
      message = as.character(j$message %||% ""),
      started_at = as.character(j$started_at %||% ""),
      finished_at = as.character(j$finished_at %||% ""),
      stringsAsFactors = FALSE
    )
  }))

  utils::tail(df, as.integer(limit))
}

find_best_reduction <- function(obj, preferred = c("umap", "tsne", "pca", "harmony", "integrated")) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
  reds <- names(obj@reductions)
  if (length(reds) == 0) return(NULL)
  for (nm in preferred) {
    if (nm %in% reds) return(nm)
  }
  reds[1]
}

# Alias used by many modules — delegates to find_best_reduction
best_reduction <- function(obj, prefer = c("umap", "tsne", "pca")) {
  find_best_reduction(obj, preferred = prefer)
}

# Canonical color palette (consolidated from clustering, visualization, dim reduction)
get_color_palette <- function(n) {
  n <- max(1L, as.integer(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    base <- RColorBrewer::brewer.pal(min(12, max(3, n)), "Set3")
    if (n > length(base)) grDevices::colorRampPalette(base)(n) else base[seq_len(n)]
  } else {
    grDevices::hcl.colors(n, "Dynamic")
  }
}

# Safe default assay accessor (consolidated from visualization, dim reduction, annotation)
safe_default_assay <- function(obj) {
  a <- DefaultAssay(obj)
  if (!is.null(a) && a %in% names(obj@assays)) return(a)
  names(obj@assays)[1]
}

# Detect best cell type column from metadata
detect_cell_type_col <- function(seurat_obj) {
  md <- colnames(seurat_obj@meta.data)
  candidates <- c("cell_type_consensus", "cell_type_ai", "cell_type_enhanced",
                   "cell_type_singler", "cell_type", "seurat_clusters",
                   "celltype", "CellType")
  for (col in candidates) {
    if (col %in% md) return(col)
  }
  char_cols <- md[sapply(seurat_obj@meta.data, function(x) is.factor(x) || is.character(x))]
  if (length(char_cols) > 0) return(char_cols[1])
  md[1]
}

adaptive_plot_index <- function(n_cells, max_points = 50000L, seed = 42L) {
  n_cells <- as.integer(n_cells %||% 0L)
  if (n_cells <= 0) return(integer(0))
  max_points <- max(1000L, as.integer(max_points))
  if (n_cells <= max_points) return(seq_len(n_cells))
  set.seed(seed)
  sort(sample.int(n_cells, max_points))
}

update_scalable_mode <- function(values, obj, threshold_cells = 100000L, threshold_features = 30000L) {
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    values$scalable_mode <- FALSE
    return(FALSE)
  }
  large <- ncol(obj) >= as.integer(threshold_cells) || nrow(obj) >= as.integer(threshold_features)
  values$scalable_mode <- isTRUE(large)
  large
}

runtime_budget_guard <- function(obj, task = "analysis", max_cells = 1000000L) {
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    return(list(ok = FALSE, message = "No active Seurat object found.", suggestion = "Load or select a dataset first."))
  }
  n_cells <- ncol(obj)
  if (n_cells > as.integer(max_cells)) {
    return(list(
      ok = FALSE,
      message = paste0("Task '", task, "' exceeds runtime budget for ", format(max_cells, big.mark = ","), " cells."),
      suggestion = "Use scalable mode, sketch/downsample, or run this step offline."
    ))
  }
  list(ok = TRUE, message = "within budget", suggestion = NA_character_)
}

memory_budget_guard <- function(obj, max_matrix_elements = 1.2e9) {
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    return(list(ok = FALSE, message = "No active dataset.", suggestion = "Load a dataset first."))
  }
  estimated <- as.numeric(nrow(obj)) * as.numeric(ncol(obj))
  if (is.finite(estimated) && estimated > max_matrix_elements) {
    return(list(
      ok = FALSE,
      message = paste0("Estimated dense matrix size is too large (", format(estimated, scientific = TRUE), " elements)."),
      suggestion = "Use sparse operations, fewer features, or sampling before this step."
    ))
  }
  list(ok = TRUE, message = "within budget", suggestion = NA_character_)
}

validate_seurat_contract <- function(
  obj,
  stage = "generic",
  require_reductions = NULL,
  require_metadata = NULL,
  min_cells = 20L,
  min_features = 50L
) {
  errors <- character(0)
  warnings <- character(0)
  suggestions <- character(0)

  if (is.null(obj)) {
    errors <- c(errors, "No active dataset found.")
    suggestions <- c(suggestions, "Upload or select a dataset before running this module.")
    return(list(ok = FALSE, errors = errors, warnings = warnings, suggestions = suggestions))
  }
  if (!inherits(obj, "Seurat")) {
    errors <- c(errors, "Active object is not a Seurat object.")
    suggestions <- c(suggestions, "Reload the dataset using a Seurat-compatible input format.")
    return(list(ok = FALSE, errors = errors, warnings = warnings, suggestions = suggestions))
  }

  if (ncol(obj) < as.integer(min_cells)) {
    errors <- c(errors, paste0("Dataset has ", ncol(obj), " cells; minimum required is ", min_cells, "."))
  }
  if (nrow(obj) < as.integer(min_features)) {
    errors <- c(errors, paste0("Dataset has ", nrow(obj), " features; minimum required is ", min_features, "."))
  }

  assay <- DefaultAssay(obj)
  if (is.null(assay) || !(assay %in% names(obj@assays))) {
    errors <- c(errors, "Default assay is missing or invalid.")
  }

  if (!is.null(require_metadata) && length(require_metadata) > 0) {
    missing_meta <- setdiff(require_metadata, colnames(obj@meta.data))
    if (length(missing_meta) > 0) {
      errors <- c(errors, paste0("Missing metadata columns: ", paste(missing_meta, collapse = ", "), "."))
      suggestions <- c(suggestions, "Run earlier preprocessing/QC steps to create required metadata columns.")
    }
  }

  reds <- names(obj@reductions)
  if (!is.null(require_reductions) && length(require_reductions) > 0) {
    if (!any(require_reductions %in% reds)) {
      errors <- c(errors, paste0("Missing required reduction(s): ", paste(require_reductions, collapse = " or "), "."))
      suggestions <- c(suggestions, "Run PCA/UMAP/tSNE before this step.")
    }
  }

  stage <- tolower(stage %||% "generic")

  if (stage %in% c("preprocessing", "preprocess")) {
    has_counts <- !is.null(.safe_layer_matrix(obj, assay = assay, preferred_layers = c("counts")))
    if (!has_counts) {
      warnings <- c(warnings, "Counts layer not found in default assay.")
      suggestions <- c(suggestions, "Ensure raw counts are available for normalization.")
    }
  }
  if (stage %in% c("dimensionality_reduction", "dimred")) {
    has_norm <- !is.null(.safe_layer_matrix(obj, assay = assay, preferred_layers = c("data")))
    if (!has_norm) {
      warnings <- c(warnings, "Normalized data layer not found.")
      suggestions <- c(suggestions, "Run normalization before dimensionality reduction for best results.")
    }
  }
  if (stage %in% c("clustering", "cluster")) {
    if (!any(c("pca", "harmony", "integrated", "umap", "tsne") %in% reds)) {
      errors <- c(errors, "Clustering requires an embedding/reduction (e.g., PCA).")
      suggestions <- c(suggestions, "Run dimensionality reduction first.")
    }
  }
  if (stage %in% c("annotation", "annotate")) {
    cluster_cols <- grep("cluster|seurat_clusters|leiden|louvain", colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
    if (length(cluster_cols) == 0) {
      warnings <- c(warnings, "No cluster-like metadata column found.")
      suggestions <- c(suggestions, "Run clustering first for cluster-level annotation or switch annotation level to cells.")
    }
  }
  if (stage %in% c("trajectory")) {
    if (!any(c("pca", "umap", "tsne", "harmony", "integrated") %in% reds)) {
      errors <- c(errors, "Trajectory analysis requires a dimensionality reduction (e.g., PCA, UMAP).")
      suggestions <- c(suggestions, "Run dimensionality reduction before trajectory analysis.")
    }
    cluster_cols <- grep("cluster|seurat_clusters|leiden|louvain", colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
    if (length(cluster_cols) == 0) {
      warnings <- c(warnings, "No cluster column found. Trajectory analysis works best with clustered data.")
      suggestions <- c(suggestions, "Run clustering before trajectory analysis.")
    }
  }
  if (stage %in% c("downstream", "de")) {
    # Check for grouping columns suitable for DE
    md <- obj@meta.data
    cat_cols <- colnames(md)[vapply(md, function(x) {
      (is.factor(x) || is.character(x)) && length(unique(x)) >= 2
    }, logical(1))]
    if (length(cat_cols) == 0) {
      errors <- c(errors, "No categorical grouping column found for differential expression.")
      suggestions <- c(suggestions, "Run clustering or add group labels before DE analysis.")
    }
  }
  if (stage %in% c("export")) {
    # Minimal checks - just object exists and has cells
    if (ncol(obj) == 0) {
      errors <- c(errors, "Dataset has 0 cells. Nothing to export.")
    }
  }

  list(
    ok = length(errors) == 0,
    errors = errors,
    warnings = unique(warnings),
    suggestions = unique(suggestions)
  )
}

format_contract_feedback <- function(contract_result) {
  if (is.null(contract_result)) return("Contract validation failed with unknown error.")
  out <- character(0)
  if (length(contract_result$errors) > 0) {
    out <- c(out, paste0("Errors: ", paste(contract_result$errors, collapse = " | ")))
  }
  if (length(contract_result$warnings) > 0) {
    out <- c(out, paste0("Warnings: ", paste(contract_result$warnings, collapse = " | ")))
  }
  if (length(contract_result$suggestions) > 0) {
    out <- c(out, paste0("Recovery: ", paste(contract_result$suggestions, collapse = " | ")))
  }
  paste(out, collapse = "\n")
}

# =============================================================================
# Fallback / recovery helpers
# =============================================================================

# Find best available reduction, with fallback chain
ensure_reduction <- function(obj, needed = "umap",
                             fallback_chain = c("umap", "tsne", "pca", "harmony", "integrated")) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
  reds <- names(obj@reductions)
  if (length(reds) == 0) return(NULL)
  if (needed %in% reds) return(needed)
  for (r in fallback_chain) {
    if (r %in% reds) return(r)
  }
  reds[1]
}

# Find best available grouping column
ensure_grouping_column <- function(obj, preferred = "seurat_clusters") {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
  md <- obj@meta.data
  if (preferred %in% colnames(md)) return(preferred)

  # Try common cluster/annotation columns
  candidates <- c("seurat_clusters", "cell_type_consensus", "cell_type_enhanced",
                   "cell_type_singler", "cell_type_ai")
  for (cand in candidates) {
    if (cand %in% colnames(md)) return(cand)
  }

  # Try pattern-based search
  cluster_cols <- grep("cluster|leiden|louvain|cell_type|annotation",
                       colnames(md), value = TRUE, ignore.case = TRUE)
  if (length(cluster_cols) > 0) return(cluster_cols[1])
  NULL
}

# Try multiple column name candidates, return the first that exists
safe_get_metadata_column <- function(obj, candidates) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)
  for (col in candidates) {
    if (col %in% colnames(obj@meta.data)) return(col)
  }
  NULL
}

# =============================================================================
# Scalable mode helpers (sketch/downsample for large datasets)
# =============================================================================

#' Downsample a Seurat object for expensive operations
#'
#' @param obj Seurat object
#' @param max_cells Maximum cells to keep
#' @param method "random" or "geometric" (geometric sketch preserves rare populations)
#' @param seed Random seed
#' @return Subsetted Seurat object
sketch_downsample <- function(obj, max_cells = 50000L, method = "random", seed = 42L) {
  if (is.null(obj) || !inherits(obj, "Seurat")) return(obj)
  n <- ncol(obj)
  max_cells <- as.integer(max_cells)
  if (n <= max_cells) return(obj)

  set.seed(seed)
  if (method == "geometric" && "pca" %in% names(obj@reductions)) {
    # Geometric sketching: sample more from sparse regions in PCA space
    pca <- Embeddings(obj, "pca")[, 1:min(10, ncol(Embeddings(obj, "pca"))), drop = FALSE]
    # Use k-means to identify cluster centers and sample proportionally
    k <- min(max_cells %/% 10, n %/% 5, 200)
    if (k < 2) k <- 2
    km <- tryCatch(kmeans(pca, centers = k, nstart = 3, iter.max = 20),
                   error = function(e) NULL)
    if (!is.null(km)) {
      # Sample proportionally from each k-means cluster
      cluster_ids <- km$cluster
      cells_per_cluster <- table(cluster_ids)
      target_per_cluster <- pmax(1L, round(max_cells * cells_per_cluster / n))
      selected <- unlist(lapply(seq_along(cells_per_cluster), function(i) {
        cl_cells <- which(cluster_ids == names(cells_per_cluster)[i])
        sample(cl_cells, min(target_per_cluster[i], length(cl_cells)))
      }))
      selected <- sort(unique(selected))
      if (length(selected) > max_cells) selected <- sample(selected, max_cells)
      return(obj[, selected])
    }
  }

  # Random downsampling fallback
  idx <- sort(sample.int(n, max_cells))
  obj[, idx]
}

#' Check if sketching/downsampling is needed for a given task
#'
#' @param values reactiveValues (checks scalable_mode and active_dataset)
#' @param budget Maximum cells for this operation
#' @return logical
should_sketch <- function(values, budget = 50000L) {
  obj <- values$active_dataset
  if (is.null(obj) || !inherits(obj, "Seurat")) return(FALSE)
  ncol(obj) > as.integer(budget)
}

register_dataset_checkpoint <- function(values, module, action, obj = NULL, params = list(), methods = list(), record_unchanged = FALSE) {
  obj <- obj %||% values$active_dataset
  if (is.null(obj) || !inherits(obj, "Seurat")) return(NULL)

  sig <- compute_seurat_signature(obj)
  old_sig <- values$dataset_signature %||% NA_character_
  changed <- !identical(sig, old_sig)

  if (is.null(values$dataset_version)) values$dataset_version <- 0L
  if (changed) values$dataset_version <- as.integer(values$dataset_version) + 1L

  if (!changed && !isTRUE(record_unchanged)) {
    return(invisible(NULL))
  }

  entry <- list(
    version = as.integer(values$dataset_version),
    signature = sig,
    changed = changed,
    module = module %||% "",
    action = action %||% "",
    timestamp = .safe_now_iso(),
    n_cells = ncol(obj),
    n_features = nrow(obj),
    params = params,
    methods = methods
  )
  values$dataset_signature <- sig
  values$dataset_history <- c(values$dataset_history %||% list(), list(entry))

  append_audit_log(
    values,
    module = module,
    action = paste0("checkpoint:", action),
    details = list(version = values$dataset_version, signature = sig, changed = changed, params = params, methods = methods)
  )
  invisible(entry)
}

is_dataset_stale <- function(values, signature) {
  if (is.null(signature) || is.na(signature)) return(FALSE)
  !identical(as.character(signature), as.character(values$dataset_signature %||% NA_character_))
}

build_method_inventory <- function() {
  data.frame(
    module = c(
      "data_input",
      "quality_control", "quality_control", "quality_control", "quality_control",
      "preprocessing", "preprocessing", "preprocessing", "preprocessing", "preprocessing",
      "dimensionality_reduction", "dimensionality_reduction", "dimensionality_reduction", "dimensionality_reduction",
      "clustering", "clustering",
      "annotation", "annotation", "annotation",
      "downstream", "downstream", "downstream", "downstream", "downstream",
      "trajectory", "trajectory",
      "pathway_network", "pathway_network",
      "export"
    ),
    method = c(
      "Seurat object ingestion",
      "DoubletFinder", "scDblFinder", "SoupX", "decontX",
      "SCTransform", "Seurat CCA integration", "Seurat RPCA integration", "Harmony", "MNN (fastMNN)",
      "PCA (irlba)", "UMAP (uwot)", "t-SNE (Rtsne)", "DiffusionMap (destiny)",
      "Leiden clustering", "Louvain/SLM clustering",
      "Marker-based scoring", "SingleR reference mapping", "Consensus annotation",
      "Wilcoxon/MAST DE", "Pseudobulk DESeq2", "Pseudobulk edgeR", "Gene set enrichment (fgsea/clusterProfiler)", "Cell-cell communication",
      "Diffusion pseudotime (DPT)", "Monocle3 trajectory",
      "WGCNA co-expression networks", "TF activity inference",
      "RMarkdown + JSON provenance export"
    ),
    package = c(
      "Seurat",
      "DoubletFinder", "scDblFinder", "SoupX", "celda",
      "Seurat/sctransform", "Seurat", "Seurat", "harmony", "batchelor",
      "Seurat/irlba", "Seurat/uwot", "Seurat/Rtsne", "destiny",
      "Seurat/igraph", "Seurat/igraph",
      "internal marker DB", "SingleR/celldex", "internal consensus",
      "Seurat/MAST", "DESeq2", "edgeR", "fgsea/clusterProfiler/msigdbr", "internal/CellChat",
      "destiny", "monocle3",
      "WGCNA", "internal (DoRothEA-style)",
      "rmarkdown/jsonlite"
    ),
    category = c(
      "canonical",
      "acceptable_alternative", "canonical", "canonical", "acceptable_alternative",
      "canonical", "canonical", "canonical", "canonical", "acceptable_alternative",
      "canonical", "canonical", "canonical", "acceptable_alternative",
      "canonical", "acceptable_alternative",
      "acceptable_alternative", "canonical", "canonical",
      "canonical", "canonical", "canonical", "canonical", "acceptable_alternative",
      "canonical", "acceptable_alternative",
      "canonical", "acceptable_alternative",
      "canonical"
    ),
    citation = c(
      "Hao et al., Cell 2021. doi:10.1016/j.cell.2021.04.048",
      "McGinnis et al., Cell Syst 2019. doi:10.1016/j.cels.2019.03.003",
      "Germain et al., F1000Research 2021. doi:10.12688/f1000research.73600.2",
      "Young & Beber, GigaScience 2020. doi:10.1093/gigascience/giaa151",
      "Yang et al., Genome Biol 2020. doi:10.1186/s13059-020-1950-6",
      "Hafemeister & Satija, Genome Biol 2019. doi:10.1186/s13059-019-1874-1",
      "Stuart et al., Cell 2019. doi:10.1016/j.cell.2019.05.031",
      "Stuart et al., Cell 2019. doi:10.1016/j.cell.2019.05.031",
      "Korsunsky et al., Nat Methods 2019. doi:10.1038/s41592-019-0619-0",
      "Haghverdi et al., Nat Biotechnol 2018. doi:10.1038/nbt.4091",
      "Satija et al., Nat Biotechnol 2015. doi:10.1038/nbt.3192",
      "McInnes et al., JOSS 2018. doi:10.21105/joss.00861",
      "van der Maaten & Hinton, JMLR 2008.",
      "Angerer et al., Bioinformatics 2016. doi:10.1093/bioinformatics/btv715",
      "Traag et al., Sci Rep 2019. doi:10.1038/s41598-019-41695-z",
      "Blondel et al., J Stat Mech 2008. doi:10.1088/1742-5468/2008/10/P10008",
      "Marker-enrichment heuristic (literature-curated DB)",
      "Aran et al., Nat Immunol 2019. doi:10.1038/s41590-018-0276-y",
      "Internal consensus rule set",
      "Finak et al., Genome Biol 2015. doi:10.1186/s13059-015-0844-5",
      "Love et al., Genome Biol 2014. doi:10.1186/s13059-014-0550-8",
      "Robinson et al., Bioinformatics 2010. doi:10.1093/bioinformatics/btp616",
      "Korotkevich et al., bioRxiv 2021. doi:10.1101/060012",
      "Jin et al., Nat Commun 2021. doi:10.1038/s41467-021-21246-9",
      "Haghverdi et al., Nat Methods 2016. doi:10.1038/nmeth.3971",
      "Cao et al., Nature 2019. doi:10.1038/s41586-019-0969-x",
      "Langfelder & Horvath, BMC Bioinf 2008. doi:10.1186/1471-2105-9-559",
      "Garcia-Alonso et al., Genome Res 2019. doi:10.1101/gr.240663.118",
      "Xie et al., R Journal 2018."
    ),
    limitations = c(
      "Input schema heterogeneity requires strict validation.",
      "Sensitive to pK/expected rate settings.", "Requires Bioc stack and conversion to SCE.", "Best with raw droplet matrix and clusters.", "May over-correct in low-contamination samples.",
      "Can be memory-intensive on very large data.", "Requires overlapping features/cells across batches.", "Faster than CCA but may miss nonlinear batch effects.", "Can over-correct biology if misconfigured.", "Assumes batch effects are linear; requires shared biology.",
      "Linear representation only.", "Stochastic embedding; seed-dependent.", "Perplexity-sensitive; poor at preserving global structure.", "Expensive for very large matrices.",
      "Resolution-sensitive; depends on graph construction.", "Can be less stable than Leiden.",
      "Depends on marker DB quality.", "Reference-dependent and species-sensitive.", "Disagreement requires explicit reporting.",
      "Assumes sufficient cells per group.", "Requires sample-level replication (min 2 per group).", "Requires sample-level replication (min 2 per group).", "Results depend on gene set database quality.", "Requires sufficient cell numbers per type; database-dependent.",
      "Trajectory assumptions may not fit all systems.", "Requires sufficient trajectory signal; sensitive to root selection.",
      "O(genes^2) memory; limited to ~2,000 genes.", "Requires curated TF-target database; limited species coverage.",
      "Report quality depends on recorded provenance."
    ),
    rationale = c(
      "Canonical ingestion ecosystem for scRNA-seq analysis.",
      "Widely used doublet detector.", "Bioconductor-standard doublet detector.", "Ambient RNA removal with profile estimation.", "Bayesian decontamination without raw matrix requirement.",
      "Academic default normalization workflow.", "Canonical multi-dataset integration.", "Scalable integration using reciprocal PCA.", "Widely accepted fast batch correction.", "Mutual nearest neighbors for batch correction.",
      "Standard first linear reduction.", "Standard nonlinear visualization/manifold tool.", "Classic nonlinear embedding for visualization.", "Common trajectory manifold alternative.",
      "Current default graph clustering standard.", "Fallback clustering alternatives.",
      "Fast interpretable baseline.", "Reference-based robust annotation.", "Improves robustness vs single method.",
      "Standard single-cell DE family.", "Study-level inference with negative binomial model.", "Study-level inference with quasi-likelihood framework.", "Pathway-level summarization of DE results.", "Infer intercellular signaling from expression data.",
      "Provides diffusion-based trajectory inference.", "Graph-based trajectory learning.",
      "Identifies gene co-expression modules.", "Infer TF regulatory activity from target expression.",
      "Supports reproducibility and peer review."
    ),
    stringsAsFactors = FALSE
  )
}

get_method_card <- function(module, inventory = NULL) {
  inventory <- inventory %||% build_method_inventory()
  mod <- tolower(module %||% "")
  inventory[tolower(inventory$module) == mod, , drop = FALSE]
}

build_method_card_ui <- function(module, inventory = NULL) {
  card <- get_method_card(module, inventory)
  if (nrow(card) == 0) return(NULL)

  rows <- lapply(seq_len(nrow(card)), function(i) {
    r <- card[i, , drop = FALSE]
    tags$tr(
      tags$td(r$method),
      tags$td(r$package),
      tags$td(r$category),
      tags$td(r$citation)
    )
  })

  tags$div(
    class = "method-card",
    style = "margin: 10px 0; padding: 10px; border: 1px solid #d9e2ef; border-radius: 8px; background: #f9fbff;",
    tags$strong(paste("Method Card -", module)),
    tags$table(
      style = "width:100%; margin-top:8px; font-size:12px;",
      tags$thead(
        tags$tr(tags$th("Method"), tags$th("Package"), tags$th("Class"), tags$th("Citation"))
      ),
      tags$tbody(rows)
    )
  )
}

collect_environment_metadata <- function() {
  si <- utils::sessionInfo()
  loaded <- names(si$otherPkgs %||% list())
  pkg_versions <- lapply(loaded, function(p) {
    list(package = p, version = as.character(utils::packageVersion(p)))
  })
  list(
    timestamp = .safe_now_iso(),
    r_version = R.version.string,
    platform = R.version$platform,
    loaded_packages = pkg_versions
  )
}

create_provenance_bundle <- function(values) {
  list(
    generated_at = .safe_now_iso(),
    dataset_version = as.integer(.state_get(values, "dataset_version", 0L)),
    dataset_signature = as.character(.state_get(values, "dataset_signature", NA_character_)),
    dataset_history = .state_get(values, "dataset_history", list()),
    audit_log = .state_get(values, "analysis_audit", data.frame()),
    method_inventory = .state_get(values, "method_inventory", build_method_inventory()),
    environment = collect_environment_metadata()
  )
}

write_provenance_bundle <- function(values, file) {
  bundle <- create_provenance_bundle(values)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(bundle, path = file, pretty = TRUE, auto_unbox = TRUE, na = "null")
  } else {
    dput(bundle, file = file)
  }
  invisible(file)
}

init_robustness_state <- function(values) {
  .state_set(values, "dataset_version", as.integer(.state_get(values, "dataset_version", 0L)))
  .state_set(values, "dataset_signature", as.character(.state_get(values, "dataset_signature", NA_character_)))
  .state_set(values, "dataset_history", .state_get(values, "dataset_history", list()))
  .state_set(values, "analysis_audit", .state_get(values, "analysis_audit", data.frame(
    timestamp = character(0),
    module = character(0),
    action = character(0),
    level = character(0),
    details = character(0),
    stringsAsFactors = FALSE
  )))
  .state_set(values, "job_queue", .state_get(values, "job_queue", list()))
  .state_set(values, "scalable_mode", isTRUE(.state_get(values, "scalable_mode", FALSE)))
  .state_set(values, "method_inventory", .state_get(values, "method_inventory", build_method_inventory()))
  invisible(TRUE)
}
