source("../../modules/de/de_core.R")
source("../../modules/de/de_methods.R")
source("../../modules/core_robustness.R")  # prepare_for_de, safe_get_assay_data

test_that("run_single_cell_de(wilcox) returns canonical schema with non-zero rows", {
  obj <- make_mock_seurat()
  res <- run_single_cell_de(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    test_method = "wilcox", seed = 42L,
    latent_vars = NULL, min_pct_prefilter = 0.01, logfc_threshold = 0
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2",
                    "AUC","cohens_d","method","seed",
                    "n_cells_group1","n_cells_group2") %in% colnames(res)))
  expect_gt(nrow(res), 0)
  expect_equal(unique(res$method), "wilcox")
  expect_equal(unique(res$seed), 42L)
})

test_that("run_single_cell_de finds known DE genes from the mock", {
  obj <- make_mock_seurat()
  res <- run_single_cell_de(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    test_method = "wilcox", seed = 42L,
    latent_vars = NULL, min_pct_prefilter = 0.01, logfc_threshold = 0
  )
  truth <- attr(obj, "true_de_up_in_A")
  top10 <- head(res$gene[order(res$p_val_adj)], 10)
  expect_true(length(intersect(truth, top10)) >= 4)   # at least 4/5 recovered
})

test_that("run_single_cell_de is reproducible across two calls with same seed", {
  obj <- make_mock_seurat()
  a <- run_single_cell_de(obj, "group", "group_A", "group_B",
                          test_method = "wilcox", seed = 7L,
                          latent_vars = NULL,
                          min_pct_prefilter = 0.01, logfc_threshold = 0)
  b <- run_single_cell_de(obj, "group", "group_A", "group_B",
                          test_method = "wilcox", seed = 7L,
                          latent_vars = NULL,
                          min_pct_prefilter = 0.01, logfc_threshold = 0)
  expect_equal(a$gene, b$gene)
  expect_equal(a$p_val, b$p_val)
  expect_equal(a$avg_log2FC, b$avg_log2FC)
})

test_that("run_single_cell_de attaches a complete audit log", {
  obj <- make_mock_seurat()
  res <- run_single_cell_de(obj, "group", "group_A", "group_B",
                            test_method = "wilcox", seed = 42L,
                            latent_vars = NULL,
                            min_pct_prefilter = 0.01, logfc_threshold = 0)
  audit <- attr(res, "de_audit")
  expect_type(audit, "list")
  expect_equal(audit$test, "wilcox")
  expect_equal(audit$ident1, "group_A")
  expect_equal(audit$n_group1, 100L)
  expect_equal(audit$n_group2, 100L)
})

test_that("run_single_cell_de(MAST) auto-injects CDR into latent_vars", {
  skip_if_not_installed("MAST")
  obj <- make_mock_seurat()
  res <- run_single_cell_de(obj, "group", "group_A", "group_B",
                            test_method = "MAST", seed = 42L,
                            latent_vars = NULL,
                            min_pct_prefilter = 0.01, logfc_threshold = 0)
  audit <- attr(res, "de_audit")
  expect_true("CDR" %in% audit$latent_vars)
})

test_that("run_pseudobulk_de(DESeq2) returns canonical schema", {
  skip_if_not_installed("DESeq2")
  obj <- make_mock_seurat()
  res <- run_pseudobulk_de(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    sample_col = "sample", method = "DESeq2", seed = 42L
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("gene","avg_log2FC","p_val","p_val_adj","method","seed",
                    "n_cells_group1","n_cells_group2") %in% colnames(res)))
  expect_equal(unique(res$method), "DESeq2")
  truth <- attr(obj, "true_de_up_in_A")
  top10 <- head(res$gene[order(res$p_val_adj)], 10)
  expect_true(length(intersect(truth, top10)) >= 3)
})

test_that("run_pseudobulk_de(edgeR) returns canonical schema", {
  skip_if_not_installed("edgeR")
  obj <- make_mock_seurat()
  res <- run_pseudobulk_de(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    sample_col = "sample", method = "edgeR", seed = 42L
  )
  expect_equal(unique(res$method), "edgeR")
  expect_gt(nrow(res), 0)
})

test_that("run_pseudobulk_de(limma-voom) returns canonical schema", {
  skip_if_not_installed("limma")
  skip_if_not_installed("edgeR")   # voom uses edgeR::DGEList()
  obj <- make_mock_seurat()
  res <- run_pseudobulk_de(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    sample_col = "sample", method = "limma-voom", seed = 42L
  )
  expect_equal(unique(res$method), "limma-voom")
  expect_gt(nrow(res), 0)
})

test_that("run_pseudobulk_de attaches an audit log with cell counts", {
  skip_if_not_installed("DESeq2")
  obj <- make_mock_seurat()
  res <- run_pseudobulk_de(obj, "group", "group_A", "group_B",
                           sample_col = "sample", method = "DESeq2", seed = 42L)
  audit <- attr(res, "de_audit")
  expect_equal(audit$test, "DESeq2")
  expect_equal(audit$n_group1, 100L)   # cell counts
  expect_equal(audit$n_group2, 100L)
})
