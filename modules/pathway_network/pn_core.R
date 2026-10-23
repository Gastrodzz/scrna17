# modules/pathway_network/pn_core.R
# Pure helpers for pathway & network analysis. No Shiny imports.
# Species, gene-ID conversion, ranked-list construction, sanitisation.

# ---------------------------------------------------------------------------
# Sole source of truth for the "no species set yet" fallback. Centralizing it
# here means the literal exists in exactly one place; everywhere else uses
# pn_resolve_species(values). When the fallback fires, the helper also returns
# a flag the caller can use to surface a notification.
# ---------------------------------------------------------------------------
PN_DEFAULT_SPECIES <- "Homo sapiens"

#' Resolve the species to use for an enrichment / network analysis.
#' @param values shiny reactiveValues (must be accessed in a reactive context).
#' @return list(species = character, fell_back = logical).
pn_resolve_species <- function(values) {
  if (!is.null(values$species_global) && nzchar(values$species_global)) {
    return(list(species = values$species_global, fell_back = FALSE))
  }
  list(species = PN_DEFAULT_SPECIES, fell_back = TRUE)
}

# ---------------------------------------------------------------------------
# Species mapping: accept Latin OR lowercase names, normalise once.
# Mirrors and extends modules/de/de_core.R helpers.
# ---------------------------------------------------------------------------

#' Map a species name to its OrgDb package name. Latin or lowercase.
pn_orgdb_for_species <- function(species) {
  key <- .pn_normalize_species_key(species)
  map <- c(
    human     = "org.Hs.eg.db",
    mouse     = "org.Mm.eg.db",
    rat       = "org.Rn.eg.db",
    zebrafish = "org.Dr.eg.db"
  )
  pkg <- map[key]
  if (is.na(pkg)) stop("Unsupported species: ", species,
                       ". Supported: human, mouse, rat, zebrafish.")
  unname(pkg)
}

#' Map a species name to its KEGG 3-letter organism code. Latin or lowercase.
pn_kegg_code_for_species <- function(species) {
  key <- .pn_normalize_species_key(species)
  map <- c(human = "hsa", mouse = "mmu", rat = "rno", zebrafish = "dre")
  code <- map[key]
  if (is.na(code)) stop("No KEGG code for species: ", species)
  unname(code)
}

#' Map a species name to a Reactome organism string ("human", "mouse", "rat", "zebrafish").
pn_reactome_for_species <- function(species) {
  key <- .pn_normalize_species_key(species)
  if (is.na(key)) stop("Unsupported species: ", species)
  key
}

#' Map a species name to msigdbr's expected species string (Latin form).
pn_msigdb_for_species <- function(species) {
  key <- .pn_normalize_species_key(species)
  map <- c(
    human     = "Homo sapiens",
    mouse     = "Mus musculus",
    rat       = "Rattus norvegicus",
    zebrafish = "Danio rerio"
  )
  out <- map[key]
  if (is.na(out)) stop("No MSigDB mapping for species: ", species)
  unname(out)
}

#' Map species name to NCBI taxonomy ID (used by STRINGdb).
pn_taxid_for_species <- function(species) {
  key <- .pn_normalize_species_key(species)
  map <- c(human = 9606L, mouse = 10090L, rat = 10116L, zebrafish = 7955L)
  out <- map[key]
  if (is.na(out)) stop("No NCBI tax id for species: ", species)
  unname(out)
}

.pn_normalize_species_key <- function(species) {
  if (is.null(species) || !nzchar(species)) return(NA_character_)
  s <- tolower(trimws(species))
  if (s %in% c("human", "homo sapiens", "h. sapiens", "hs"))      return("human")
  if (s %in% c("mouse", "mus musculus", "m. musculus", "mm"))     return("mouse")
  if (s %in% c("rat", "rattus norvegicus", "r. norvegicus", "rn")) return("rat")
  if (s %in% c("zebrafish", "danio rerio", "d. rerio", "dr"))     return("zebrafish")
  NA_character_
}

# ---------------------------------------------------------------------------
# Gene ID conversion: SYMBOL <-> ENTREZID (and ENSEMBL on demand).
# Single canonical wrapper. All ORA/GSEA paths go through this.
# ---------------------------------------------------------------------------

#' Convert gene IDs via clusterProfiler::bitr() with the right OrgDb.
#' @param genes character vector of input IDs.
#' @param from / to keytype names accepted by AnnotationDbi.
#' @param species species name (Latin or lowercase).
#' @return data.frame with columns <from> and <to>; NAs dropped.
convert_gene_ids <- function(genes, from = "SYMBOL", to = "ENTREZID",
                             species = "human") {
  if (length(genes) == 0)
    return(data.frame(stats::setNames(list(character(0), character(0)),
                                      c(from, to)), stringsAsFactors = FALSE))
  orgdb <- pn_orgdb_for_species(species)
  if (!requireNamespace(orgdb, quietly = TRUE))
    stop("OrgDb package '", orgdb, "' not installed. Install: BiocManager::install('", orgdb, "').")
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("clusterProfiler not installed. Install: BiocManager::install('clusterProfiler').")
  res <- tryCatch(
    suppressMessages(suppressWarnings(
      clusterProfiler::bitr(genes, fromType = from, toType = to,
                            OrgDb = orgdb, drop = TRUE))),
    error = function(e) data.frame(stats::setNames(list(character(0), character(0)),
                                                   c(from, to)), stringsAsFactors = FALSE)
  )
  unique(res)
}

#' Convenience: SYMBOL -> ENTREZID, returning a character vector keyed by SYMBOL.
symbols_to_entrez <- function(genes, species) {
  m <- convert_gene_ids(genes, from = "SYMBOL", to = "ENTREZID", species = species)
  if (nrow(m) == 0) return(character(0))
  out <- m$ENTREZID
  names(out) <- m$SYMBOL
  out
}

# ---------------------------------------------------------------------------
# GSEA ranked-list construction. Always uses NOMINAL p_val (not p_val_adj),
# per Subramanian 2005 + Love 2014: BH adjustment compresses the lower tail
# non-uniformly, distorting NES.
# ---------------------------------------------------------------------------

#' Build a sorted, named numeric vector for GSEA from a DE result.
#' @param de_result data.frame with columns gene, avg_log2FC, p_val (and optional p_val_adj).
#' @param metric one of "signed_neglog10p", "log2fc", "t_statistic".
#' @param seed integer; used to add tiny tie-breaking jitter.
build_ranked_list <- function(de_result,
                              metric = c("signed_neglog10p", "log2fc", "t_statistic"),
                              seed = 42L) {
  metric <- match.arg(metric)
  df <- de_result[!is.na(de_result$avg_log2FC), , drop = FALSE]
  if ("p_val" %in% colnames(df)) {
    df <- df[!is.na(df$p_val), , drop = FALSE]
    p_col <- "p_val"
  } else if ("p_val_adj" %in% colnames(df)) {
    df <- df[!is.na(df$p_val_adj), , drop = FALSE]
    p_col <- "p_val_adj"
  } else {
    p_col <- NA_character_
  }

  scores <- switch(metric,
    signed_neglog10p = {
      if (is.na(p_col)) sign(df$avg_log2FC)   # fallback: just direction
      else sign(df$avg_log2FC) * -log10(pmax(df[[p_col]], 1e-300))
    },
    log2fc       = df$avg_log2FC,
    t_statistic  = if ("t_stat" %in% colnames(df)) as.numeric(df$t_stat)
                   else if ("statistic" %in% colnames(df)) as.numeric(df$statistic)
                   else sign(df$avg_log2FC) * -log10(pmax(df[[p_col]], 1e-300))
  )

  names(scores) <- df$gene
  scores <- scores[!is.na(scores)]
  scores <- sort(scores, decreasing = TRUE)

  # Tiny jitter to break ties (prevents fgsea ties warning, deterministic via seed)
  set.seed(as.integer(seed))
  scores <- scores + stats::runif(length(scores), -1e-6, 1e-6)
  sort(scores, decreasing = TRUE)
}

# ---------------------------------------------------------------------------
# Sanitise enrichment data frames so plotly's range() / serialize never sees
# factors or list columns. Lifted (and trimmed) from downstream_core.R.
# ---------------------------------------------------------------------------

pn_sanitize_enrichment_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(df)
  num_cols <- c("pvalue", "p.adjust", "qvalue", "qvalues", "Count", "NES",
                "ES", "nMoreExtreme", "size", "padj", "pval", "log2err")
  for (col in intersect(num_cols, colnames(df))) {
    df[[col]] <- as.numeric(as.character(df[[col]]))
  }
  for (col in colnames(df)) {
    if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
  }
  list_cols <- vapply(df, is.list, logical(1))
  if (any(list_cols)) {
    for (col in colnames(df)[list_cols]) {
      df[[col]] <- vapply(df[[col]], function(x)
        paste(utils::head(x, 20), collapse = ", "), character(1))
    }
  }
  df
}

# ---------------------------------------------------------------------------
# Background gene universe: spec §3 says always derive from the active dataset,
# i.e. rownames(values$active_dataset). Provide a helper so callers don't
# touch the Seurat object directly.
# ---------------------------------------------------------------------------

#' Background gene universe for ORA. Always rownames(seurat_obj) when present.
#' Falls back to NULL if no Seurat object is given (callers must override).
pn_background_universe <- function(seurat_obj = NULL, custom = NULL) {
  if (!is.null(custom) && length(custom) > 0) return(unique(as.character(custom)))
  if (!is.null(seurat_obj)) return(unique(rownames(seurat_obj)))
  NULL
}
