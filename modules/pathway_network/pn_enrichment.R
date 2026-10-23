# modules/pathway_network/pn_enrichment.R
# ORA + GSEA engines for pathway enrichment.
# Depends on: modules/pathway_network/pn_core.R (sourced first by pn_server.R).

# ---------------------------------------------------------------------------
# Package readiness — lazy load on first use, not at source time.
# ---------------------------------------------------------------------------
.pn_ensure_pkgs <- function(method = c("ORA", "GSEA", "BOTH")) {
  method <- match.arg(method)
  pkgs <- c("clusterProfiler")
  if (method %in% c("GSEA", "BOTH")) pkgs <- c(pkgs, "fgsea")
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE))
      stop(sprintf("Package '%s' required. Install: BiocManager::install('%s')", p, p))
    if (!paste0("package:", p) %in% search()) {
      suppressPackageStartupMessages(library(p, character.only = TRUE))
    }
  }
}

# ---------------------------------------------------------------------------
# Build a gene-set list (term -> character vector of gene symbols) for the
# requested database. For databases that have a live API (KEGG, REACTOME, GO,
# WikiPathways), return a sentinel list and let the engines take the live path.
# For MSigDB collections (H, C2, C5, C6), return a real named list.
# ---------------------------------------------------------------------------

#' @param database one of: "H","C2","C5","C6","KEGG","REACTOME","GO_BP","GO_MF",
#'   "GO_CC","WP","REACTOME_PA","STRING","MSIGDB_KEGG","MSIGDB_REACTOME","MSIGDB_WP"
#' @param species lowercase or Latin
load_gene_sets <- function(database, species) {
  if (!requireNamespace("msigdbr", quietly = TRUE) &&
      database %in% c("H","C2","C5","C6","MSIGDB_KEGG","MSIGDB_REACTOME","MSIGDB_WP")) {
    stop("msigdbr not installed. Install: BiocManager::install('msigdbr').")
  }

  # Live-API databases — do not load gene sets here, return a sentinel.
  if (database %in% c("KEGG", "REACTOME", "REACTOME_PA",
                      "GO_BP", "GO_MF", "GO_CC", "WP", "STRING")) {
    return(list(.live_api = TRUE, database = database, species = species))
  }

  msig_species <- pn_msigdb_for_species(species)

  # msigdbr 10.0+ renamed "category"/"subcategory" -> "collection"/"subcollection"
  # and split "CP:KEGG" into "CP:KEGG_LEGACY" + "CP:KEGG_MEDICUS". We accept the
  # new arg names if available, fall back to the old for older installs, and
  # merge both KEGG variants when MSIGDB_KEGG is requested.
  .msig_query <- function(coll, sub = NULL) {
    fn <- msigdbr::msigdbr
    if ("collection" %in% names(formals(fn))) {
      if (is.null(sub)) fn(species = msig_species, collection = coll)
      else              fn(species = msig_species, collection = coll, subcollection = sub)
    } else {
      if (is.null(sub)) fn(species = msig_species, category = coll)
      else              fn(species = msig_species, category = coll, subcategory = sub)
    }
  }
  msig_data <- if (database == "MSIGDB_KEGG") {
    # KEGG was split — combine legacy + medicus to recover the old behaviour
    a <- tryCatch(.msig_query("C2", "CP:KEGG_LEGACY"),  error = function(e) NULL)
    b <- tryCatch(.msig_query("C2", "CP:KEGG_MEDICUS"), error = function(e) NULL)
    do.call(rbind, Filter(Negate(is.null), list(a, b)))
  } else if (database == "MSIGDB_REACTOME") {
    .msig_query("C2", "CP:REACTOME")
  } else if (database == "MSIGDB_WP") {
    .msig_query("C2", "CP:WIKIPATHWAYS")
  } else {
    # H, C2, C5, C6 — top-level collections
    .msig_query(database)
  }
  if (nrow(msig_data) == 0) return(NULL)

  gene_sets <- split(msig_data$gene_symbol, msig_data$gs_name)
  # Filter very small / very large sets to keep statistics reliable
  gene_sets[lengths(gene_sets) >= 10 & lengths(gene_sets) <= 500]
}

# ---------------------------------------------------------------------------
# Gene-list preparation: build (up_genes, down_genes, ranked, universe) from
# a DE result data.frame.
# ---------------------------------------------------------------------------

#' @param de_result data.frame with at least gene, avg_log2FC, p_val, p_val_adj.
#' @param ora_pval threshold on adj p for ORA gene selection.
#' @param ora_lfc absolute log2FC threshold for ORA gene selection.
#' @param gsea_metric ranking metric for GSEA.
#' @param universe character vector of background genes (rownames(seurat_obj)).
#' @param seed integer.
prepare_gene_lists <- function(de_result, ora_pval = 0.05, ora_lfc = 0.25,
                               gsea_metric = "signed_neglog10p",
                               universe = NULL, seed = 42L) {
  d <- de_result
  d$avg_log2FC <- as.numeric(as.character(d$avg_log2FC))
  d$p_val_adj  <- as.numeric(as.character(d$p_val_adj))
  if ("p_val" %in% colnames(d)) d$p_val <- as.numeric(as.character(d$p_val))
  d$gene <- as.character(d$gene)

  up_genes <- d$gene[!is.na(d$p_val_adj) &
                     !is.na(d$avg_log2FC) &
                     d$p_val_adj < ora_pval &
                     d$avg_log2FC >  ora_lfc]
  down_genes <- d$gene[!is.na(d$p_val_adj) &
                       !is.na(d$avg_log2FC) &
                       d$p_val_adj < ora_pval &
                       d$avg_log2FC < -ora_lfc]
  ranked <- build_ranked_list(d, metric = gsea_metric, seed = seed)
  bg <- if (!is.null(universe) && length(universe) > 0) unique(universe)
        else unique(d$gene)

  list(up_genes = up_genes, down_genes = down_genes,
       ranked_genes = ranked, universe = bg)
}

# ---------------------------------------------------------------------------
# ORA — Over-Representation Analysis.
# Dispatches by database: live (clusterProfiler::enrich*) or MSigDB (enricher).
# ---------------------------------------------------------------------------

#' @param gene_lists list from prepare_gene_lists().
#' @param gene_sets either (a) a sentinel list from load_gene_sets() for live-API DBs
#'                  or (b) a named list of gene sets for MSigDB.
#' @param params list with ora_pval_cutoff, max_pathways.
#' @return list(up = data.frame|NULL, down = data.frame|NULL).
run_ora_analysis <- function(gene_lists, gene_sets, params) {
  .pn_ensure_pkgs("ORA")
  pcut <- params$ora_pval_cutoff %||% 0.05

  # ------ Live-API path ------
  if (is.list(gene_sets) && isTRUE(gene_sets$.live_api)) {
    db <- gene_sets$database
    sp <- gene_sets$species
    .one_dir <- function(genes, direction) {
      if (length(genes) == 0) return(NULL)
      out <- .pn_run_live_ora(genes, db, sp, pcut)
      if (is.null(out) || nrow(out) == 0) return(NULL)
      out$Direction <- direction
      out
    }
    return(list(up   = .one_dir(gene_lists$up_genes, "Up"),
                down = .one_dir(gene_lists$down_genes, "Down")))
  }

  # ------ MSigDB path (gene_sets is a named list of vectors) ------
  if (is.null(gene_sets) || length(gene_sets) == 0)
    return(list(up = NULL, down = NULL,
                error = "No gene sets loaded for the selected database."))

  term2gene <- do.call(rbind, lapply(names(gene_sets), function(term)
    data.frame(term = term, gene = gene_sets[[term]], stringsAsFactors = FALSE)))

  .one_dir <- function(genes, direction) {
    if (length(genes) == 0) return(NULL)
    res <- tryCatch(
      clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = term2gene,
        universe = gene_lists$universe,
        pAdjustMethod = "BH",
        pvalueCutoff = pcut,
        qvalueCutoff = pcut),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    out <- pn_sanitize_enrichment_df(as.data.frame(res))
    out$Direction <- direction
    utils::head(out, params$max_pathways %||% 50)
  }
  list(up   = .one_dir(gene_lists$up_genes,   "Up"),
       down = .one_dir(gene_lists$down_genes, "Down"))
}

#' Internal: live-API ORA dispatch for KEGG/REACTOME/GO/WP/STRING/REACTOME_PA.
.pn_run_live_ora <- function(genes, db, species, pcut) {
  if (db == "STRING") {
    if (!requireNamespace("STRINGdb", quietly = TRUE))
      stop("STRINGdb not installed. Install: BiocManager::install('STRINGdb').")
    taxid <- pn_taxid_for_species(species)
    sdb <- STRINGdb::STRINGdb$new(version = "12.0", species = taxid,
                                  score_threshold = 400)
    mapped <- sdb$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
    if (nrow(mapped) < 3) return(NULL)
    raw <- sdb$get_enrichment(mapped$STRING_id)
    if (is.null(raw) || nrow(raw) == 0) return(NULL)
    return(.pn_string_to_ora_df(raw))
  }

  entrez <- symbols_to_entrez(genes, species)
  if (length(entrez) == 0) return(NULL)

  res <- if (db == "KEGG") {
    clusterProfiler::enrichKEGG(
      gene = unname(entrez), organism = pn_kegg_code_for_species(species),
      pvalueCutoff = pcut, pAdjustMethod = "BH")
  } else if (db %in% c("REACTOME", "REACTOME_PA")) {
    if (!requireNamespace("ReactomePA", quietly = TRUE))
      stop("ReactomePA not installed. Install: BiocManager::install('ReactomePA').")
    ReactomePA::enrichPathway(
      gene = unname(entrez), organism = pn_reactome_for_species(species),
      pvalueCutoff = pcut, pAdjustMethod = "BH", readable = TRUE)
  } else if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
    ont <- sub("GO_", "", db)
    clusterProfiler::enrichGO(
      gene = unname(entrez), OrgDb = pn_orgdb_for_species(species), ont = ont,
      pvalueCutoff = pcut, pAdjustMethod = "BH", readable = TRUE)
  } else if (db == "WP") {
    clusterProfiler::enrichWP(
      gene = unname(entrez), organism = pn_msigdb_for_species(species),
      pvalueCutoff = pcut, pAdjustMethod = "BH")
  } else {
    return(NULL)
  }

  if (is.null(res) || nrow(res) == 0) return(NULL)
  pn_sanitize_enrichment_df(as.data.frame(res))
}

.pn_string_to_ora_df <- function(s) {
  if (is.null(s) || nrow(s) == 0) return(NULL)
  bg_total <- if ("number_of_genes_in_background" %in% colnames(s))
    max(as.integer(s$number_of_genes_in_background), na.rm = TRUE) else 20000L
  data.frame(
    ID          = as.character(s$term),
    Description = if ("description" %in% colnames(s)) as.character(s$description) else as.character(s$term),
    GeneRatio   = paste0(s$number_of_genes, "/",
                         if ("number_of_genes_in_background" %in% colnames(s))
                           s$number_of_genes_in_background else s$number_of_genes),
    BgRatio     = paste0(if ("number_of_genes_in_background" %in% colnames(s))
                           s$number_of_genes_in_background else s$number_of_genes, "/", bg_total),
    pvalue      = as.numeric(s$p_value),
    p.adjust    = as.numeric(if ("fdr" %in% colnames(s)) s$fdr else s$p_value),
    qvalue      = as.numeric(if ("fdr" %in% colnames(s)) s$fdr else s$p_value),
    Count       = as.integer(s$number_of_genes),
    geneID      = if ("preferredNames" %in% colnames(s)) as.character(s$preferredNames) else "",
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# GSEA — Gene Set Enrichment Analysis.
# Dispatches by database: live (cP::gse*) or MSigDB (fgsea on gene-set list).
# ---------------------------------------------------------------------------

#' @param ranked_genes named numeric vector (sorted descending).
#' @param gene_sets sentinel list (live API) OR named list of gene sets.
#' @param params list with gsea_min_size, gsea_max_size, gsea_nperm,
#'        gsea_pval_cutoff, max_pathways, seed.
#' @return data.frame or NULL.
run_gsea_analysis <- function(ranked_genes, gene_sets, params) {
  .pn_ensure_pkgs("GSEA")
  pcut  <- params$gsea_pval_cutoff %||% 0.05
  nperm <- max(params$gsea_nperm %||% 1000L, 1000L)
  seed  <- as.integer(params$seed %||% 42L)

  # ------ Live-API path ------
  if (is.list(gene_sets) && isTRUE(gene_sets$.live_api)) {
    db <- gene_sets$database
    sp <- gene_sets$species
    return(.pn_run_live_gsea(ranked_genes, db, sp, pcut, seed,
                             params$max_pathways %||% 50))
  }

  if (is.null(gene_sets) || length(gene_sets) == 0) return(NULL)

  # fgsea scoreType: "std" assumes both signs in stats; if all stats are >=0
  # (rare but happens when DE returns only positive log2FC) we must use "pos",
  # otherwise fgsea silently drops pathways with NA p-values.
  has_neg <- any(ranked_genes < 0, na.rm = TRUE)
  has_pos <- any(ranked_genes > 0, na.rm = TRUE)
  score_type <- if (has_neg && has_pos) "std"
                else if (has_pos) "pos"
                else "neg"

  set.seed(seed)
  res <- tryCatch(
    fgsea::fgsea(pathways = gene_sets, stats = ranked_genes,
                 minSize = params$gsea_min_size %||% 15L,
                 maxSize = params$gsea_max_size %||% 500L,
                 nPermSimple = nperm,
                 scoreType = score_type),
    error = function(e) NULL
  )
  if (is.null(res) || nrow(res) == 0) return(NULL)

  df <- as.data.frame(res)
  df <- pn_sanitize_enrichment_df(df)
  if ("padj" %in% colnames(df))
    df <- df[order(df$padj), , drop = FALSE]
  utils::head(df, params$max_pathways %||% 50)
}

#' Internal: live-API GSEA dispatch.
.pn_run_live_gsea <- function(ranked_genes, db, species, pcut, seed, max_pathways) {
  syms <- names(ranked_genes)
  m <- convert_gene_ids(syms, "SYMBOL", "ENTREZID", species = species)
  if (nrow(m) == 0) return(NULL)
  entrez_ranks <- ranked_genes[m$SYMBOL]
  names(entrez_ranks) <- m$ENTREZID
  entrez_ranks <- sort(entrez_ranks[!is.na(entrez_ranks)], decreasing = TRUE)
  if (length(entrez_ranks) < 10) return(NULL)

  set.seed(seed)
  res <- if (db == "KEGG") {
    clusterProfiler::gseKEGG(
      geneList = entrez_ranks, organism = pn_kegg_code_for_species(species),
      pvalueCutoff = pcut, pAdjustMethod = "BH",
      verbose = FALSE, seed = TRUE)
  } else if (db %in% c("REACTOME", "REACTOME_PA")) {
    if (!requireNamespace("ReactomePA", quietly = TRUE))
      stop("ReactomePA not installed. Install: BiocManager::install('ReactomePA').")
    ReactomePA::gsePathway(
      geneList = entrez_ranks, organism = pn_reactome_for_species(species),
      pvalueCutoff = pcut, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE)
  } else if (db %in% c("GO_BP", "GO_MF", "GO_CC")) {
    ont <- sub("GO_", "", db)
    clusterProfiler::gseGO(
      geneList = entrez_ranks, OrgDb = pn_orgdb_for_species(species), ont = ont,
      pvalueCutoff = pcut, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE)
  } else if (db == "WP") {
    clusterProfiler::gseWP(
      geneList = entrez_ranks, organism = pn_msigdb_for_species(species),
      pvalueCutoff = pcut, pAdjustMethod = "BH", verbose = FALSE, seed = TRUE)
  } else {
    return(NULL)
  }
  if (is.null(res) || nrow(res) == 0) return(NULL)

  df <- as.data.frame(res)
  # Standardise column names so downstream renderers work for both fgsea and cP
  if ("Description" %in% colnames(df) && !"pathway" %in% colnames(df))
    df$pathway <- df$Description
  if ("pvalue" %in% colnames(df) && !"pval" %in% colnames(df))
    df$pval <- df$pvalue
  if ("p.adjust" %in% colnames(df) && !"padj" %in% colnames(df))
    df$padj <- df$p.adjust
  df <- pn_sanitize_enrichment_df(df)
  if ("padj" %in% colnames(df))
    df <- df[order(df$padj), , drop = FALSE]
  utils::head(df, max_pathways)
}

# ---------------------------------------------------------------------------
# Audit metadata for an enrichment run.
# ---------------------------------------------------------------------------
build_enrichment_audit <- function(database, method, species,
                                   n_input, n_background, seed,
                                   n_results) {
  list(
    database     = database,
    method       = method,
    species      = species,
    n_input      = as.integer(n_input),
    n_background = as.integer(n_background),
    seed         = as.integer(seed),
    n_results    = as.integer(n_results),
    timestamp    = Sys.time(),
    r_version    = R.version$version.string
  )
}

# ---------------------------------------------------------------------------
# Comparative enrichment — run ORA / GSEA on N DE results side by side.
# Re-uses run_ora_analysis / run_gsea_analysis per condition, then merges
# results into long-form data frames with a "condition" column. This is the
# open-source equivalent of clusterProfiler::compareCluster() and IPA's
# "Comparison Analysis" feature.
# ---------------------------------------------------------------------------

#' @param de_runs named list. Each element is a list with $data (DE data.frame
#'   with columns gene, avg_log2FC, p_val_adj, optional p_val) and $audit.
#'   The list names become the "condition" labels (e.g. "ko vs wt · osteo").
#' @param database one of the keys accepted by load_gene_sets().
#' @param species species name (Latin or lowercase).
#' @param method "ORA", "GSEA", or "BOTH".
#' @param params list with: ora_pval, ora_lfc, ora_pval_cutoff, gsea_metric,
#'   gsea_min_size, gsea_max_size, gsea_nperm, gsea_pval_cutoff, max_pathways,
#'   universe (vector or NULL), seed.
#' @return list(per_condition=<named list>, combined_ora=<df|NULL>,
#'              combined_gsea=<df|NULL>, sig_pathways=<named list per condition>,
#'              audit=<list>).
run_compare_enrichment <- function(de_runs, database, species,
                                   method = c("ORA", "GSEA", "BOTH"),
                                   params = list()) {
  method <- match.arg(method)
  if (length(de_runs) < 2)
    stop("run_compare_enrichment needs at least 2 DE runs (got ",
         length(de_runs), ").")

  defaults <- list(
    ora_pval = 0.05, ora_lfc = 0.25, ora_pval_cutoff = 0.05,
    gsea_metric = "signed_neglog10p",
    gsea_min_size = 15L, gsea_max_size = 500L, gsea_nperm = 1000L,
    gsea_pval_cutoff = 0.05, max_pathways = 30L,
    universe = NULL, seed = 42L
  )
  for (k in names(defaults))
    if (is.null(params[[k]])) params[[k]] <- defaults[[k]]

  # Load gene sets once (live-API path returns a sentinel that is harmless
  # to share across calls; MSigDB path returns a list that's reused).
  gs <- load_gene_sets(database, species)

  per_condition <- vector("list", length(de_runs))
  names(per_condition) <- names(de_runs)

  for (cond in names(de_runs)) {
    df <- de_runs[[cond]]$data
    if (is.null(df) || nrow(df) == 0) {
      per_condition[[cond]] <- list(ora = NULL, gsea = NULL,
                                     n_genes_input = 0L, error = "empty DE")
      next
    }
    gl <- prepare_gene_lists(
      df,
      ora_pval = params$ora_pval, ora_lfc = params$ora_lfc,
      gsea_metric = params$gsea_metric,
      universe = params$universe, seed = params$seed)

    out <- list(ora = NULL, gsea = NULL,
                n_genes_input = length(gl$up_genes) + length(gl$down_genes))
    if (method %in% c("ORA", "BOTH")) {
      out$ora <- tryCatch(
        run_ora_analysis(gl, gs, params),
        error = function(e) list(up = NULL, down = NULL,
                                  error = conditionMessage(e)))
    }
    if (method %in% c("GSEA", "BOTH")) {
      out$gsea <- tryCatch(
        run_gsea_analysis(gl$ranked_genes, gs, params),
        error = function(e) NULL)
    }
    per_condition[[cond]] <- out
  }

  combined_ora  <- .combine_ora_across_conditions(per_condition)
  combined_gsea <- .combine_gsea_across_conditions(per_condition)
  sig_pathways  <- .significant_pathways_per_condition(
    per_condition, ora_padj_cutoff = params$ora_pval_cutoff,
    gsea_padj_cutoff = params$gsea_pval_cutoff)

  list(
    per_condition = per_condition,
    combined_ora  = combined_ora,
    combined_gsea = combined_gsea,
    sig_pathways  = sig_pathways,
    audit         = build_compare_audit(
                       conditions = names(de_runs),
                       database = database, method = method,
                       species = species,
                       n_genes_per_condition = vapply(per_condition,
                         function(p) as.integer(p$n_genes_input %||% 0L),
                         integer(1)),
                       seed = as.integer(params$seed))
  )
}

# Merge per-condition ORA results into a long-form data.frame.
.combine_ora_across_conditions <- function(per_condition) {
  parts <- list()
  for (cond in names(per_condition)) {
    o <- per_condition[[cond]]$ora
    if (is.null(o)) next
    if (!is.null(o$up) && nrow(o$up) > 0) {
      x <- o$up; x$condition <- cond; parts[[length(parts) + 1L]] <- x
    }
    if (!is.null(o$down) && nrow(o$down) > 0) {
      x <- o$down; x$condition <- cond; parts[[length(parts) + 1L]] <- x
    }
  }
  if (length(parts) == 0) return(NULL)
  # Robust rbind: harmonize columns across conditions
  all_cols <- unique(unlist(lapply(parts, colnames)))
  parts <- lapply(parts, function(d) {
    miss <- setdiff(all_cols, colnames(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, parts)
}

# Merge per-condition GSEA results into a long-form data.frame.
.combine_gsea_across_conditions <- function(per_condition) {
  parts <- list()
  for (cond in names(per_condition)) {
    g <- per_condition[[cond]]$gsea
    if (is.null(g) || nrow(g) == 0) next
    x <- g; x$condition <- cond; parts[[length(parts) + 1L]] <- x
  }
  if (length(parts) == 0) return(NULL)
  all_cols <- unique(unlist(lapply(parts, colnames)))
  parts <- lapply(parts, function(d) {
    miss <- setdiff(all_cols, colnames(d))
    for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
  do.call(rbind, parts)
}

# Per-condition sets of significant pathway IDs (used by Venn/UpSet helpers).
.significant_pathways_per_condition <- function(per_condition,
                                                ora_padj_cutoff = 0.05,
                                                gsea_padj_cutoff = 0.05) {
  out <- list()
  for (cond in names(per_condition)) {
    p <- per_condition[[cond]]
    ora_ids <- character(0); gsea_ids <- character(0)
    if (!is.null(p$ora$up) && "p.adjust" %in% colnames(p$ora$up))
      ora_ids <- c(ora_ids, p$ora$up$ID[p$ora$up$p.adjust < ora_padj_cutoff])
    if (!is.null(p$ora$down) && "p.adjust" %in% colnames(p$ora$down))
      ora_ids <- c(ora_ids, p$ora$down$ID[p$ora$down$p.adjust < ora_padj_cutoff])
    if (!is.null(p$gsea) && nrow(p$gsea) > 0) {
      padj_col <- if ("padj" %in% colnames(p$gsea)) "padj"
                  else if ("p.adjust" %in% colnames(p$gsea)) "p.adjust"
                  else NA_character_
      id_col   <- if ("pathway" %in% colnames(p$gsea)) "pathway"
                  else if ("ID" %in% colnames(p$gsea)) "ID"
                  else if ("Description" %in% colnames(p$gsea)) "Description"
                  else NA_character_
      if (!is.na(padj_col) && !is.na(id_col))
        gsea_ids <- p$gsea[[id_col]][as.numeric(p$gsea[[padj_col]]) <
                                       gsea_padj_cutoff]
    }
    out[[cond]] <- list(
      ora_ids  = unique(stats::na.omit(as.character(ora_ids))),
      gsea_ids = unique(stats::na.omit(as.character(gsea_ids))),
      union    = unique(stats::na.omit(as.character(c(ora_ids, gsea_ids))))
    )
  }
  out
}

#' Audit metadata for a comparative enrichment run.
build_compare_audit <- function(conditions, database, method, species,
                                n_genes_per_condition, seed) {
  list(
    n_conditions          = length(conditions),
    conditions            = as.character(conditions),
    database              = database,
    method                = method,
    species               = species,
    n_genes_per_condition = as.integer(n_genes_per_condition),
    seed                  = as.integer(seed),
    timestamp             = Sys.time(),
    r_version             = R.version$version.string
  )
}

# Null-coalescing operator (used by ORA/GSEA functions above).
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a)) a else b
}
