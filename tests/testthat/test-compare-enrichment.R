source("../../modules/pathway_network/pn_core.R")
source("../../modules/pathway_network/pn_enrichment.R")

# Build a tiny well-mocked 2-condition DE run set with REAL human gene names so
# MSigDB Hallmark gives a sensible signal we can assert on.
.build_two_de_runs <- function() {
  set.seed(2026)
  # Condition 1 — IFN-γ response
  ifn <- c("STAT1","IRF1","STAT2","JAK1","JAK2","ISG15","MX1","OAS1","OAS2",
           "OAS3","IFIT1","IFIT2","IFIT3","CXCL10","CXCL9","HLA-A","HLA-B",
           "HLA-C","B2M","TAP1","TAP2","PSMB8","PSMB9","CCL5","CXCL11",
           "IRF7","IRF9","DDX60","RSAD2","XAF1","BST2","SAMD9","HERC6",
           "TRIM21","TRIM22","USP18")
  de1 <- data.frame(gene = ifn,
                    avg_log2FC = runif(length(ifn), 1, 4),
                    p_val      = 10^(-runif(length(ifn), 5, 30)),
                    p_val_adj  = 10^(-runif(length(ifn), 4, 25)))
  # Condition 2 — EMT
  emt <- c("VIM","SNAI1","SNAI2","ZEB1","ZEB2","TWIST1","FN1","CDH2",
           "MMP2","MMP9","TGFB1","TGFB2","SERPINE1","COL1A1","COL1A2",
           "COL3A1","ACTA2","TAGLN","S100A4","SPARC","POSTN","FAP","LOX",
           "LOXL2","ITGB1","ITGA5","TIMP1","CTGF","CXCL12","BMP2",
           "SMAD2","SMAD3")
  de2 <- data.frame(gene = emt,
                    avg_log2FC = runif(length(emt), 1, 4),
                    p_val      = 10^(-runif(length(emt), 5, 30)),
                    p_val_adj  = 10^(-runif(length(emt), 4, 25)))
  bg <- unique(c(ifn, emt, paste0("BG", 1:5000)))
  list(
    runs = list(
      "Epithelial: ko vs wt"  = list(data = de1,
        audit = list(test = "wilcox", ident1 = "ko", ident2 = "wt")),
      "Mesenchymal: ko vs wt" = list(data = de2,
        audit = list(test = "wilcox", ident1 = "ko", ident2 = "wt"))
    ),
    universe = bg
  )
}

test_that("run_compare_enrichment errors with fewer than 2 runs", {
  fixture <- .build_two_de_runs()
  expect_error(
    run_compare_enrichment(de_runs = fixture$runs[1],
                           database = "H", species = "human",
                           method = "ORA", params = list(universe = fixture$universe)),
    "at least 2"
  )
})

test_that("run_compare_enrichment(MSigDB Hallmark, ORA) returns per-condition + combined", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("clusterProfiler")
  fixture <- .build_two_de_runs()
  res <- run_compare_enrichment(
    de_runs = fixture$runs, database = "H", species = "human",
    method = "ORA",
    params = list(ora_pval = 0.05, ora_lfc = 0.25,
                  ora_pval_cutoff = 0.5, max_pathways = 20,
                  universe = fixture$universe, seed = 42L))
  expect_named(res,
    c("per_condition", "combined_ora", "combined_gsea",
      "sig_pathways", "audit"))
  expect_equal(length(res$per_condition), 2)
  expect_true(!is.null(res$combined_ora))
  expect_true("condition" %in% colnames(res$combined_ora))
  expect_setequal(unique(res$combined_ora$condition), names(fixture$runs))
})

test_that("run_compare_enrichment(MSigDB Hallmark, GSEA) recovers expected biology", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("fgsea")
  fixture <- .build_two_de_runs()
  res <- run_compare_enrichment(
    de_runs = fixture$runs, database = "H", species = "human",
    method = "GSEA",
    params = list(gsea_min_size = 5, gsea_max_size = 500,
                  gsea_nperm = 1000, gsea_pval_cutoff = 0.5,
                  max_pathways = 20, universe = fixture$universe, seed = 42L))
  expect_true(!is.null(res$combined_gsea))
  # Condition labels carried through
  expect_setequal(unique(res$combined_gsea$condition), names(fixture$runs))
  # Per-condition pathway sets exist and overlap reasonably with the truth
  by_cond <- split(res$combined_gsea, res$combined_gsea$condition)
  pw_col <- if ("pathway" %in% colnames(res$combined_gsea)) "pathway"
            else "Description"
  ifn_hits <- grepl("INTERFERON|IL6_JAK_STAT3|ALLOGRAFT",
                    by_cond[[1]][[pw_col]])
  emt_hits <- grepl("EPITHELIAL_MESENCHYMAL|COAGULATION|MYOGENESIS",
                    by_cond[[2]][[pw_col]])
  expect_true(any(ifn_hits))
  expect_true(any(emt_hits))
})

test_that("build_compare_audit captures the right metadata", {
  a <- build_compare_audit(
    conditions = c("A", "B", "C"),
    database = "H", method = "BOTH", species = "human",
    n_genes_per_condition = c(120L, 80L, 50L),
    seed = 42L)
  expect_equal(a$n_conditions, 3L)
  expect_equal(a$conditions, c("A", "B", "C"))
  expect_equal(a$database, "H")
  expect_equal(a$method, "BOTH")
  expect_equal(a$seed, 42L)
  expect_equal(a$n_genes_per_condition, c(120L, 80L, 50L))
  expect_s3_class(a$timestamp, "POSIXct")
})

test_that("sig_pathways structure is correct (Venn input)", {
  skip_if_not_installed("msigdbr")
  skip_if_not_installed("clusterProfiler")
  fixture <- .build_two_de_runs()
  res <- run_compare_enrichment(
    de_runs = fixture$runs, database = "H", species = "human",
    method = "BOTH",
    params = list(ora_pval_cutoff = 0.5, gsea_pval_cutoff = 0.5,
                  gsea_min_size = 5, gsea_max_size = 500,
                  gsea_nperm = 1000, max_pathways = 20,
                  universe = fixture$universe, seed = 42L))
  expect_equal(length(res$sig_pathways), 2)
  for (cond in names(res$sig_pathways)) {
    s <- res$sig_pathways[[cond]]
    expect_named(s, c("ora_ids", "gsea_ids", "union"))
    expect_true(is.character(s$ora_ids))
    expect_true(is.character(s$gsea_ids))
    expect_true(is.character(s$union))
  }
})

test_that("combine helpers harmonise columns even when conditions differ in schema", {
  per_cond <- list(
    A = list(ora = list(
      up = data.frame(ID = "P1", Description = "x", p.adjust = 0.01,
                      Direction = "Up", stringsAsFactors = FALSE),
      down = NULL)),
    B = list(ora = list(
      up = data.frame(ID = "P2", Description = "y", p.adjust = 0.02,
                      Direction = "Up", extra_col = "z",
                      stringsAsFactors = FALSE),
      down = NULL))
  )
  combined <- .combine_ora_across_conditions(per_cond)
  expect_true(!is.null(combined))
  expect_equal(nrow(combined), 2)
  expect_true("extra_col" %in% colnames(combined))
  expect_true("condition" %in% colnames(combined))
  expect_setequal(unique(combined$condition), c("A", "B"))
})
