# Verify the gallery push helpers write to values$figure_gallery and
# values$csv_library when the umap_gallery helpers are available.
suppressPackageStartupMessages({
  library(shiny)
  source("../../modules/umap_gallery.R")
  source("../../modules/core_robustness.R")
  source("../../modules/de/de_core.R")
  source("../../modules/de/de_methods.R")
  source("../../modules/de/de_server.R")
  source("../../modules/pathway_network/pn_core.R")
  source("../../modules/pathway_network/pn_enrichment.R")
  source("../../modules/pathway_network/pn_network.R")
  source("../../modules/pathway_network/pn_server.R")
})

# Reusable fixture: tiny DE result + run_compare_enrichment output
.fixture <- function() {
  set.seed(11)
  ifn <- c("STAT1","IRF1","JAK1","ISG15","MX1","OAS1","IFIT1","CXCL10","HLA-A",
           "B2M","TAP1","PSMB8","IRF7","IRF9","RSAD2","BST2","CCL5","CXCL11",
           "TAP2","PSMB9")
  emt <- c("VIM","SNAI1","ZEB1","TWIST1","FN1","CDH2","MMP2","TGFB1","COL1A1",
           "ACTA2","TAGLN","SPARC","CTGF","BMP2","SMAD2","SMAD3","LOX",
           "ITGB1","LOXL2","POSTN")
  mk <- function(g) {
    df <- data.frame(
      gene = g,
      avg_log2FC = runif(length(g), 1, 4),
      p_val      = 10^(-runif(length(g), 5, 30)),
      p_val_adj  = 10^(-runif(length(g), 4, 25)),
      pct.1      = runif(length(g), 0.2, 0.9),
      pct.2      = runif(length(g), 0.0, 0.3))
    attr(df, "de_audit") <- list(
      test = "wilcox", group_col = "genotype",
      ident1 = "ko", ident2 = "wt",
      n_cells_total = 200L, n_group1 = 100L, n_group2 = 100L,
      latent_vars = character(0),
      n_genes_tested = nrow(df),
      timestamp = Sys.time(),
      seurat_version = "5.4.0", r_version = R.version$version.string)
    df
  }
  de1 <- mk(ifn); de2 <- mk(emt)
  bg <- unique(c(ifn, emt, paste0("BG", 1:5000)))
  list(de1 = de1, de2 = de2, universe = bg,
       runs = list(
         "Epithelial: ko vs wt"  = list(data = de1, audit = attr(de1, "de_audit")),
         "Mesenchymal: ko vs wt" = list(data = de2, audit = attr(de2, "de_audit"))))
}

test_that(".de_push_to_gallery writes figures + CSVs into values", {
  fx <- .fixture()
  v <- reactiveValues()
  fake_input <- list(pval_cutoff = 0.05, logfc_threshold = 0.10,
                     min_pct_display = 0.10, label_top_n = 20)
  isolate(.de_push_to_gallery(v, fx$de1, attr(fx$de1, "de_audit"),
                              fake_input, label = "ko vs wt · wilcox"))
  isolate({
    expect_true(!is.null(v$figure_gallery))
    expect_gte(length(v$figure_gallery), 3)   # volcano + MA + venn (no Seurat -> no heatmap)
    fig_labels <- vapply(v$figure_gallery, function(g) g$label, character(1))
    expect_true(any(grepl("Volcano", fig_labels)))
    expect_true(any(grepl("MA plot|Expression overlap", fig_labels)))
    expect_true(!is.null(v$csv_library))
    csv_descs <- vapply(v$csv_library, function(g) g$description, character(1))
    expect_true(any(grepl("DE results", csv_descs)))
    expect_true(any(grepl("audit log", csv_descs)))
  })
})

test_that(".pn_push_enrichment_to_gallery writes when MSigDB ORA+GSEA runs", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("fgsea")
  fx <- .fixture()
  gl <- prepare_gene_lists(fx$de1, ora_pval = 0.05, ora_lfc = 0.25,
                           universe = fx$universe, seed = 42L)
  gs <- load_gene_sets("H", "human")
  ora  <- run_ora_analysis(gl, gs,  list(ora_pval_cutoff = 0.5, max_pathways = 10))
  gsea <- run_gsea_analysis(gl$ranked_genes, gs,
            list(gsea_min_size = 5, gsea_max_size = 500,
                 gsea_nperm = 1000, max_pathways = 10, seed = 42L))
  v <- reactiveValues()
  isolate(.pn_push_enrichment_to_gallery(v,
    list(ora = ora, gsea = gsea), "H", "BOTH", 20L))
  isolate({
    expect_true(!is.null(v$figure_gallery))
    expect_gte(length(v$figure_gallery), 2)
    fig_labels <- vapply(v$figure_gallery, function(g) g$label, character(1))
    expect_true(any(grepl("dot plot", fig_labels)))
    expect_true(any(grepl("bar plot", fig_labels)))
  })
})

test_that(".pn_push_compare_to_gallery writes dotplot + overlap + 3 CSVs", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("fgsea")
  fx <- .fixture()
  res <- run_compare_enrichment(de_runs = fx$runs,
    database = "H", species = "human", method = "BOTH",
    params = list(ora_pval_cutoff = 0.5, gsea_pval_cutoff = 0.5,
                  gsea_min_size = 5, gsea_max_size = 500, gsea_nperm = 1000,
                  max_pathways = 20, universe = fx$universe, seed = 42L))
  v <- reactiveValues()
  isolate(.pn_push_compare_to_gallery(v, res, "H", "BOTH", 20L))
  isolate({
    fig_labels <- vapply(v$figure_gallery %||% list(),
                         function(g) g$label, character(1))
    csv_descs <- vapply(v$csv_library %||% list(),
                        function(g) g$description, character(1))
    expect_true(any(grepl("Comparative dotplot",  fig_labels)))
    expect_true(any(grepl("Pathway overlap",      fig_labels)))
    expect_true(any(grepl("Comparative ORA",      csv_descs)))
    expect_true(any(grepl("Comparative GSEA",     csv_descs)))
    expect_true(any(grepl("Pathway-by-condition", csv_descs)))
  })
})

test_that("Hard-coded species literal exists in exactly one place (pn_core.R)", {
  files <- list.files(c("../../modules/de", "../../modules/pathway_network"),
                      pattern = "\\.R$", full.names = TRUE, recursive = FALSE)
  hits <- character(0)
  for (f in files) {
    if (basename(f) %in% c("pn_core.R", "de_core.R")) next
    txt <- readLines(f, warn = FALSE)
    # Look for "Homo sapiens" or "Mus musculus" as code, not in comments/doc lines
    code_lines <- txt[!grepl("^\\s*#'?", txt)]
    if (any(grepl('"Homo sapiens"|"Mus musculus"', code_lines)))
      hits <- c(hits, f)
  }
  expect_length(hits, 0)
})

test_that("Hex-color string literals do not appear in plotly call bodies (only in palette constants)", {
  for (f in c("../../modules/de/de_server.R",
              "../../modules/pathway_network/pn_server.R")) {
    txt <- readLines(f, warn = FALSE)
    # Strip palette definition + comments + rgba tuples for transparent backgrounds
    bad <- grep('color = "#[0-9A-Fa-f]{6}"', txt, value = TRUE)
    # Should be zero — palette tokens are referenced as .DE_PAL$xxx / .PN_PAL$xxx
    expect_length(bad, 0)
  }
})

test_that("pn_resolve_species returns fell_back=TRUE when species_global is unset", {
  v <- reactiveValues()
  isolate({
    r <- pn_resolve_species(v)
    expect_equal(r$species, PN_DEFAULT_SPECIES)
    expect_true(r$fell_back)
    v$species_global <- "Mus musculus"
    r2 <- pn_resolve_species(v)
    expect_equal(r2$species, "Mus musculus")
    expect_false(r2$fell_back)
  })
})
