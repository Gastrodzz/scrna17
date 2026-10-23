# =============================================================================
# Downstream Core — Shared Helpers
# modules/downstream/downstream_core.R
#
# Shared utilities used by downstream_ccc, downstream_tf, downstream_geneset,
# downstream_composition, and downstream_ai sub-modules. Also referenced by
# modules/pathway_network/pn_*.R via .sanitize_enrichment_df / safe_renderPlotly.
# (DE has moved to modules/de/. Pathway enrichment has moved to modules/pathway_network/.)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(ggplot2)
  library(future)
  library(promises)
})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R
# .mdtohtml(), call_claude_api(), call_ai(), build_expert_prompt(),
# ai_thread_message(), ai_thread_loading(), ai_thread_error(),
# get_model_for_task(), AI_EXPERT_PERSONA are defined in ai_assistant.R

# --- Shared Data Accessors ---

#' Get the active Seurat object from values
getActiveSeuratObj <- function(values) {
  values$active_dataset
}

#' Get the most recent DE results from values
getActiveDE <- function(values) {
  values$last_de_results
}

#' Set a deterministic analysis seed
buildAnalysisSeed <- function(seed = NULL) {
  s <- if (!is.null(seed)) as.integer(seed) else 42L
  set.seed(s)
  s
}

# --- Data Sanitisation Helpers ---

#' Strip ALL factors from a data frame so plotly never sees one.
#' plotly_build() calls range() during serialization which errors on factors.
.defactorize <- function(df) {
  for (col in colnames(df)) {
    if (is.factor(df[[col]])) {
      df[[col]] <- as.character(df[[col]])
    }
  }
  # Also drop list columns (e.g. leadingEdge) which plotly can't serialize
  list_cols <- sapply(df, is.list)
  if (any(list_cols)) df <- df[, !list_cols, drop = FALSE]
  df
}

#' Sanitize enrichment data frame — convert ALL factor columns to character/numeric
.sanitize_enrichment_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)
  num_cols <- c("pvalue", "p.adjust", "qvalue", "qvalues", "Count", "NES",
                "ES", "nMoreExtreme", "size", "padj", "pval", "log2err")
  for (col in intersect(num_cols, colnames(df))) {
    df[[col]] <- as.numeric(as.character(df[[col]]))
  }
  char_cols <- c("ID", "Description", "GeneRatio", "BgRatio", "geneID",
                 "Direction", "pathway", "leadingEdge")
  for (col in intersect(char_cols, colnames(df))) {
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
  }
  for (col in colnames(df)) {
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
  }
  df
}

#' Safely extract a column from a data frame
.safe_col <- function(df, col, n = NULL) {
  tryCatch({
    if (!col %in% colnames(df)) return(character(0))
    vals <- df[[col]]
    if (!is.null(n)) vals <- head(vals, n)
    vals
  }, error = function(e) character(0))
}

#' Register a result in values for cross-module access
registerResult <- function(values, slot, data) {
  values[[slot]] <- data
}

#' Safe plotly renderer — strips factors from all trace data before serialization
#' to prevent plotly_build() calling range() on factors (which throws an error).
#' Use this in place of renderPlotly() wherever factor columns may be present.
safe_renderPlotly <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) expr <- substitute(expr)
  renderPlotly({
    tryCatch({
      p <- eval(expr, envir = env)
      if (inherits(p, "plotly")) {
        # Defactorize all trace data
        if (!is.null(p$x$data)) {
          p$x$data <- lapply(p$x$data, function(trace) {
            for (nm in names(trace)) {
              if (is.factor(trace[[nm]])) trace[[nm]] <- as.character(trace[[nm]])
            }
            trace
          })
        }
      }
      p
    }, error = function(e) {
      plotly_empty() %>% layout(title = list(text = paste("Plot error:", e$message)))
    })
  })
}

#' Convert human gene names to mouse orthologs
#' Fallback chain: homologene → biomaRt::getLDS() → case conversion (last resort)
.to_mouse_gene <- function(human_genes) {
  if (length(human_genes) == 0) return(character(0))

  # --- Strategy 1: homologene (fast, offline) ---
  if (requireNamespace("homologene", quietly = TRUE)) {
    tryCatch({
      mapped <- homologene::homologene(human_genes, inTax = 9606, outTax = 10090)
      if (is.data.frame(mapped) && nrow(mapped) > 0) {
        result <- human_genes
        idx <- match(human_genes, mapped[[1]])
        valid <- !is.na(idx)
        result[valid] <- mapped[[2]][idx[valid]]
        n_mapped <- sum(valid)
        if (n_mapped > 0) {
          message(sprintf("[ortholog] homologene mapped %d/%d genes to mouse",
                          n_mapped, length(human_genes)))
          return(result)
        }
      }
    }, error = function(e) {
      message("[ortholog] homologene failed: ", e$message)
    })
  }

  # --- Strategy 2: biomaRt::getLDS() (comprehensive, requires internet) ---
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    tryCatch({
      human_mart <- biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                                          mirror = "useast")
      mouse_mart <- biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",
                                          mirror = "useast")
      orth <- biomaRt::getLDS(
        attributes = "hgnc_symbol", filters = "hgnc_symbol",
        values = human_genes, mart = human_mart,
        attributesL = "mgi_symbol", martL = mouse_mart,
        uniqueRows = TRUE
      )
      if (is.data.frame(orth) && nrow(orth) > 0) {
        # Take first match per human gene (some have multiple orthologs)
        orth <- orth[!duplicated(orth[[1]]), ]
        result <- human_genes
        idx <- match(human_genes, orth[[1]])
        valid <- !is.na(idx)
        result[valid] <- orth[[2]][idx[valid]]
        message(sprintf("[ortholog] biomaRt mapped %d/%d genes to mouse",
                        sum(valid), length(human_genes)))
        return(result)
      }
    }, error = function(e) {
      message("[ortholog] biomaRt getLDS failed: ", e$message,
              ". Falling back to case conversion.")
    })
  }

  # --- Strategy 3: Case conversion (last resort, imprecise) ---
  # Standard mouse nomenclature: first letter uppercase, rest lowercase
  # This is WRONG for many gene families (Hla, mt-, etc.) but better than nothing
  warning("Neither homologene nor biomaRt available. Using case conversion for mouse genes ",
          "(inaccurate for HLA, miRNA, etc.). Install: install.packages('homologene') or ",
          "BiocManager::install('biomaRt')")
  paste0(toupper(substr(human_genes, 1, 1)),
         tolower(substr(human_genes, 2, nchar(human_genes))))
}
