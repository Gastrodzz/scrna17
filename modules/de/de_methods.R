# modules/de/de_methods.R
# DE engines and dispatcher. Sourced before de_server.R.
# Depends on: modules/de/de_core.R, modules/core_robustness.R (prepare_for_de).

# --- Constants ---------------------------------------------------------------
.SC_METHODS <- c("wilcox", "MAST", "t", "bimod", "LR", "negbinom")
.PB_METHODS <- c("DESeq2", "edgeR", "limma-voom")
.MIN_CELLS_SC <- 10L
.MIN_REPS_PB <- 2L
.METHOD_PKG <- list(
  MAST       = "MAST",
  DESeq2     = "DESeq2",
  edgeR      = "edgeR",
  `limma-voom` = "limma",
  negbinom   = "MASS",
  LR         = "stats",
  wilcox     = NULL, t = NULL, bimod = NULL
)

# --- Pre-flight validation ---------------------------------------------------
preflight_de <- function(obj, group_col, ident1, ident2 = NULL,
                         analysis_type = c("singlecell", "pseudobulk"),
                         test_method, sample_col = NULL,
                         min_cells = .MIN_CELLS_SC,
                         min_reps  = .MIN_REPS_PB) {
  analysis_type <- match.arg(analysis_type)
  if (!group_col %in% colnames(obj@meta.data))
    stop("Group column '", group_col, "' not found in metadata.")
  if (!ident1 %in% as.character(obj@meta.data[[group_col]]))
    stop("ident1 '", ident1, "' not found in column '", group_col, "'.")

  # Method-specific package check
  pkg <- .METHOD_PKG[[test_method]]
  if (!is.null(pkg) && !requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' required for test '", test_method,
         "'. Install: BiocManager::install('", pkg, "').")

  # Cell-count check (single-cell)
  if (analysis_type == "singlecell") {
    n1 <- sum(obj@meta.data[[group_col]] == ident1)
    if (n1 < min_cells)
      stop("Group '", ident1, "' has fewer than ", min_cells,
           " cells (n=", n1, "). Use pseudobulk or pick another group.")
    if (!is.null(ident2) && !ident2 %in% c("rest", "all", "")) {
      n2 <- sum(obj@meta.data[[group_col]] == ident2)
      if (n2 < min_cells)
        stop("Group '", ident2, "' has fewer than ", min_cells, " cells (n=", n2, ").")
    }
  }

  # Replicate check (pseudobulk)
  if (analysis_type == "pseudobulk") {
    if (is.null(sample_col) || !nzchar(sample_col))
      stop("sample_col is required for pseudobulk DE.")
    if (!sample_col %in% colnames(obj@meta.data))
      stop("Sample column '", sample_col, "' not found in metadata.")

    md <- obj@meta.data
    s1 <- unique(md[[sample_col]][md[[group_col]] == ident1])
    if (length(s1) < min_reps)
      stop("Pseudobulk requires at least ", min_reps, " replicates per group; ",
           "group '", ident1, "' has ", length(s1), ".")
    if (!is.null(ident2) && !ident2 %in% c("rest", "all", "")) {
      s2 <- unique(md[[sample_col]][md[[group_col]] == ident2])
      if (length(s2) < min_reps)
        stop("Pseudobulk requires at least ", min_reps, " replicates per group; ",
             "group '", ident2, "' has ", length(s2), ".")
    }
  }

  invisible(TRUE)
}

# --- Output normalisation ----------------------------------------------------

#' Normalise an engine's raw output to the canonical DE schema.
#' @param raw data.frame with at least p_val, avg_log2FC, pct.1, pct.2 and rownames=gene.
normalize_de_output <- function(raw, method, seed,
                                n_cells_group1, n_cells_group2) {
  if (is.null(raw) || nrow(raw) == 0) {
    return(.empty_de_df())
  }
  cn <- colnames(raw)
  fc_col <- if ("avg_log2FC" %in% cn) "avg_log2FC"
            else if ("avg_logFC" %in% cn) "avg_logFC"
            else if ("logFC" %in% cn)     "logFC"
            else NA_character_
  p_col   <- if ("p_val"   %in% cn) "p_val"
             else if ("PValue" %in% cn) "PValue"
             else if ("pvalue" %in% cn) "pvalue"
             else NA_character_
  padj_col <- if ("p_val_adj" %in% cn) "p_val_adj"
              else if ("FDR" %in% cn)  "FDR"
              else if ("padj" %in% cn) "padj"
              else NA_character_

  gene <- if ("gene" %in% cn) as.character(raw$gene) else rownames(raw)
  avg_log2FC <- if (!is.na(fc_col)) as.numeric(raw[[fc_col]]) else rep(NA_real_, nrow(raw))
  p_val      <- if (!is.na(p_col))   as.numeric(raw[[p_col]])   else rep(NA_real_, nrow(raw))
  p_val_adj  <- if (!is.na(padj_col)) as.numeric(raw[[padj_col]])
                else stats::p.adjust(p_val, method = "BH")
  pct.1 <- if ("pct.1" %in% cn) as.numeric(raw$pct.1) else rep(NA_real_, nrow(raw))
  pct.2 <- if ("pct.2" %in% cn) as.numeric(raw$pct.2) else rep(NA_real_, nrow(raw))

  AUC <- ifelse(is.na(pct.1) | is.na(pct.2), NA_real_, (pct.1 + (1 - pct.2)) / 2)
  cohens_d <- avg_log2FC / log2(exp(1))   # log2 -> natural log

  data.frame(
    gene = gene,
    avg_log2FC = avg_log2FC,
    p_val = p_val,
    p_val_adj = p_val_adj,
    pct.1 = pct.1,
    pct.2 = pct.2,
    AUC = AUC,
    cohens_d = cohens_d,
    method = rep(method, length(gene)),
    seed = rep(as.integer(seed), length(gene)),
    n_cells_group1 = rep(as.integer(n_cells_group1), length(gene)),
    n_cells_group2 = rep(as.integer(n_cells_group2), length(gene)),
    stringsAsFactors = FALSE
  )
}

.empty_de_df <- function() {
  data.frame(
    gene = character(0), avg_log2FC = numeric(0),
    p_val = numeric(0), p_val_adj = numeric(0),
    pct.1 = numeric(0), pct.2 = numeric(0),
    AUC = numeric(0), cohens_d = numeric(0),
    method = character(0), seed = integer(0),
    n_cells_group1 = integer(0), n_cells_group2 = integer(0),
    stringsAsFactors = FALSE
  )
}

# --- Single-cell engine ------------------------------------------------------

#' Single-cell DE via Seurat::FindMarkers.
#' Methods: wilcox, MAST, t, bimod, LR, negbinom.
#' Always tests on the FULL dataset — downsampling is forbidden by design
#' (per user requirement: every cell counts).
run_single_cell_de <- function(obj, group_col, ident1, ident2 = NULL,
                               test_method = "wilcox",
                               seed = 42L,
                               latent_vars = NULL,
                               min_pct_prefilter = 0.01,
                               logfc_threshold = 0,
                               only_pos = FALSE,
                               min_pct_display = 0) {
  if (!test_method %in% .SC_METHODS)
    stop("Unsupported single-cell method: ", test_method,
         ". Supported: ", paste(.SC_METHODS, collapse = ", "))

  obj <- prepare_for_de(obj)
  Seurat::Idents(obj) <- group_col
  set.seed(seed)

  # MAST: auto-add CDR (cellular detection rate) as a latent variable
  if (test_method == "MAST") {
    expr <- safe_get_assay_data(obj, "data")
    obj$CDR <- as.numeric(scale(colMeans(expr > 0)))
    latent_vars <- unique(c(latent_vars, "CDR"))
  }

  ident2_arg <- if (is.null(ident2) || ident2 %in% c("rest", "")) NULL else ident2
  cells1 <- which(as.character(obj@meta.data[[group_col]]) == ident1)
  cells2 <- if (is.null(ident2_arg)) {
    which(as.character(obj@meta.data[[group_col]]) != ident1)
  } else {
    which(as.character(obj@meta.data[[group_col]]) == ident2_arg)
  }
  n1 <- length(cells1); n2 <- length(cells2)

  raw <- Seurat::FindMarkers(
    obj, ident.1 = ident1, ident.2 = ident2_arg,
    test.use = test_method,
    min.pct = min_pct_prefilter,
    logfc.threshold = 0,
    min.diff.pct = -Inf,
    only.pos = FALSE,
    latent.vars = latent_vars,
    max.cells.per.ident = Inf,        # full dataset always — no downsampling
    random.seed = seed,
    verbose = FALSE
  )
  raw$gene <- rownames(raw)

  out <- normalize_de_output(raw, method = test_method, seed = seed,
                             n_cells_group1 = n1, n_cells_group2 = n2)

  # Display-only post-filters
  if (isTRUE(only_pos))     out <- out[out$avg_log2FC > 0, , drop = FALSE]
  if (logfc_threshold > 0)  out <- out[abs(out$avg_log2FC) >= logfc_threshold, , drop = FALSE]
  if (min_pct_display > 0)
    out <- out[pmax(out$pct.1, out$pct.2, na.rm = TRUE) >= min_pct_display, , drop = FALSE]

  out <- out[order(out$p_val_adj, -abs(out$avg_log2FC)), ]

  attr(out, "de_audit") <- build_de_audit(
    seed = seed, test = test_method, group_col = group_col,
    ident1 = ident1, ident2 = ident2,
    n_cells_total = ncol(obj), n_group1 = n1, n_group2 = n2,
    latent_vars = latent_vars,
    n_genes_tested = nrow(out)
  )
  out
}

# --- Top-level dispatcher ----------------------------------------------------

#' Top-level DE dispatcher used by the Shiny server.
#' Always tests on the full dataset — no downsampling, no max-cell caps.
#' Returns: list(data = data.frame|NULL, error = character|NULL).
run_de_analysis <- function(obj,
                            group_col, ident1, ident2 = NULL,
                            test_method = "wilcox",
                            analysis_type = c("singlecell", "pseudobulk"),
                            sample_col = NULL,
                            subset_col = NULL, subset_val = NULL,
                            latent_vars = NULL,
                            seed = 42L,
                            min_pct_prefilter = 0.01,
                            min_pct_display   = 0.10,
                            logfc_threshold   = 0,
                            only_pos = FALSE) {
  analysis_type <- match.arg(analysis_type)

  tryCatch({
    # 1. Optional cell-type subset
    if (!is.null(subset_col) && nzchar(subset_col) &&
        !is.null(subset_val) && nzchar(subset_val)) {
      keep <- which(as.character(obj@meta.data[[subset_col]]) == subset_val)
      if (length(keep) < .MIN_CELLS_SC)
        stop("Too few cells (", length(keep), ") for '", subset_val,
             "' in '", subset_col, "'. Need at least ", .MIN_CELLS_SC, ".")
      obj <- subset(obj, cells = keep)
    }

    # 2. Pre-flight (will throw on bad inputs)
    preflight_de(obj, group_col = group_col, ident1 = ident1, ident2 = ident2,
                 analysis_type = analysis_type, test_method = test_method,
                 sample_col = sample_col)

    # 3. Route
    data <- if (analysis_type == "pseudobulk") {
      run_pseudobulk_de(obj, group_col = group_col,
                        ident1 = ident1, ident2 = ident2,
                        sample_col = sample_col, method = test_method,
                        seed = seed)
    } else {
      run_single_cell_de(obj, group_col = group_col,
                         ident1 = ident1, ident2 = ident2,
                         test_method = test_method, seed = seed,
                         latent_vars = latent_vars,
                         min_pct_prefilter = min_pct_prefilter,
                         min_pct_display   = min_pct_display,
                         logfc_threshold   = logfc_threshold,
                         only_pos = only_pos)
    }
    list(data = data, error = NULL)
  }, error = function(e) list(data = NULL, error = conditionMessage(e)))
}

# --- Pseudobulk engine -------------------------------------------------------

#' Pseudobulk DE: aggregate counts per sample, run DESeq2/edgeR/limma-voom.
#' @param method one of "DESeq2", "edgeR", "limma-voom".
run_pseudobulk_de <- function(obj, group_col, ident1, ident2 = NULL,
                              sample_col, method = "DESeq2", seed = 42L) {
  if (!method %in% .PB_METHODS)
    stop("Unsupported pseudobulk method: ", method,
         ". Supported: ", paste(.PB_METHODS, collapse = ", "))

  obj <- prepare_for_de(obj)
  set.seed(seed)

  md <- obj@meta.data
  cells1 <- rownames(md)[md[[group_col]] == ident1]
  ident2_label <- if (is.null(ident2) || ident2 %in% c("rest", "all", ""))
                    "rest" else ident2
  cells2 <- if (ident2_label == "rest")
              rownames(md)[md[[group_col]] != ident1]
            else
              rownames(md)[md[[group_col]] == ident2_label]

  n1_cells <- length(cells1); n2_cells <- length(cells2)
  s1 <- unique(md[cells1, sample_col]); s2 <- unique(md[cells2, sample_col])

  counts_mat <- safe_get_assay_data(obj, "counts")
  agg <- function(samples, group) {
    do.call(cbind, lapply(samples, function(s) {
      cells <- rownames(md)[md[[sample_col]] == s &
                            ((group == ident1 & md[[group_col]] == ident1) |
                             (group != ident1 & md[[group_col]] != ident1))]
      if (length(cells) == 0) return(NULL)
      Matrix::rowSums(counts_mat[, cells, drop = FALSE])
    }))
  }
  m1 <- agg(s1, ident1); m2 <- agg(s2, "other")
  if (is.null(m1) || is.null(m2))
    return(.empty_de_df_with_audit(method, seed, n1_cells, n2_cells,
                                   group_col, ident1, ident2,
                                   ncol(obj)))
  bulk <- cbind(m1, m2)
  colnames(bulk) <- c(s1, s2)
  group_factor <- factor(c(rep(ident1, length(s1)),
                           rep(ident2_label, length(s2))),
                         levels = c(ident2_label, ident1))
  col_data <- data.frame(group = group_factor, row.names = colnames(bulk))

  raw <- switch(method,
    "DESeq2"     = .pb_deseq2(bulk, col_data, ident1, ident2_label),
    "edgeR"      = .pb_edger (bulk, col_data, ident1, ident2_label),
    "limma-voom" = .pb_voom  (bulk, col_data, ident1, ident2_label)
  )

  out <- normalize_de_output(raw, method = method, seed = seed,
                             n_cells_group1 = n1_cells,
                             n_cells_group2 = n2_cells)
  out <- out[!is.na(out$p_val_adj), , drop = FALSE]
  out <- out[order(out$p_val_adj, -abs(out$avg_log2FC)), ]

  attr(out, "de_audit") <- build_de_audit(
    seed = seed, test = method, group_col = group_col,
    ident1 = ident1, ident2 = ident2,
    n_cells_total = ncol(obj),
    n_group1 = n1_cells, n_group2 = n2_cells,
    latent_vars = NULL,
    n_genes_tested = nrow(out)
  )
  out
}

.empty_de_df_with_audit <- function(method, seed, n1, n2,
                                    group_col, ident1, ident2, n_total) {
  out <- .empty_de_df()
  attr(out, "de_audit") <- build_de_audit(
    seed = seed, test = method, group_col = group_col,
    ident1 = ident1, ident2 = ident2,
    n_cells_total = n_total, n_group1 = n1, n_group2 = n2,
    latent_vars = NULL,
    n_genes_tested = 0L
  )
  out
}

.pb_deseq2 <- function(bulk, col_data, ident1, ident2_label) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(bulk), colData = col_data, design = ~group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(
    dds, contrast = c("group", ident1, ident2_label)))
  res$gene <- rownames(res)
  names(res)[names(res) == "log2FoldChange"] <- "avg_log2FC"
  names(res)[names(res) == "padj"]           <- "p_val_adj"
  names(res)[names(res) == "pvalue"]         <- "p_val"
  res
}

.pb_edger <- function(bulk, col_data, ident1, ident2_label) {
  y <- edgeR::DGEList(counts = round(bulk), group = col_data$group)
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group, data = col_data)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit)
  res <- as.data.frame(edgeR::topTags(qlf, n = Inf, sort.by = "none"))
  res$gene <- rownames(res)
  names(res)[names(res) == "logFC"]  <- "avg_log2FC"
  names(res)[names(res) == "FDR"]    <- "p_val_adj"
  names(res)[names(res) == "PValue"] <- "p_val"
  res
}

.pb_voom <- function(bulk, col_data, ident1, ident2_label) {
  y <- edgeR::DGEList(counts = round(bulk), group = col_data$group)
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group, data = col_data)
  v <- limma::voom(y, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  coef_name <- paste0("group", ident1)
  res <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  names(res)[names(res) == "logFC"]    <- "avg_log2FC"
  names(res)[names(res) == "P.Value"]  <- "p_val"
  names(res)[names(res) == "adj.P.Val"] <- "p_val_adj"
  res
}

#' Build the audit-log list attached to every DE result.
#' By design, every DE run uses the FULL dataset — no downsampling or max-cell
#' fields are reported because they are constants ("none" / Inf) for every run.
build_de_audit <- function(seed, test, group_col, ident1, ident2,
                           n_cells_total, n_group1, n_group2,
                           latent_vars, n_genes_tested) {
  list(
    seed           = as.integer(seed),
    test           = test,
    group_col      = group_col,
    ident1         = ident1,
    ident2         = if (is.null(ident2)) NA_character_ else ident2,
    n_cells_total  = as.integer(n_cells_total),
    n_group1       = as.integer(n_group1),
    n_group2       = as.integer(n_group2),
    latent_vars    = if (is.null(latent_vars)) character(0) else as.character(latent_vars),
    n_genes_tested = as.integer(n_genes_tested),
    timestamp      = Sys.time(),
    seurat_version = as.character(utils::packageVersion("Seurat")),
    r_version      = R.version$version.string
  )
}
