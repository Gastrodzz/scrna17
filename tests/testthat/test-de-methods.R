source("../../modules/de/de_core.R")
source("../../modules/de/de_methods.R")
source("../../modules/core_robustness.R")  # prepare_for_de, safe_get_assay_data

test_that("preflight_de errors when group_col missing from metadata", {
  obj <- make_mock_seurat()
  expect_error(
    preflight_de(obj, group_col = "nope", ident1 = "group_A",
                 analysis_type = "singlecell", test_method = "wilcox"),
    "not found in metadata"
  )
})

test_that("preflight_de errors when ident1 has fewer than 10 cells (singlecell)", {
  obj <- make_mock_seurat()
  obj$tiny_group <- "X"
  obj$tiny_group[1:5] <- "Y"   # only 5 cells in 'Y'
  expect_error(
    preflight_de(obj, group_col = "tiny_group", ident1 = "Y",
                 analysis_type = "singlecell", test_method = "wilcox"),
    "fewer than 10 cells"
  )
})

test_that("preflight_de errors when sample_col missing for pseudobulk", {
  obj <- make_mock_seurat()
  expect_error(
    preflight_de(obj, group_col = "group", ident1 = "group_A",
                 analysis_type = "pseudobulk", test_method = "DESeq2",
                 sample_col = "no_such_col"),
    "Sample column"
  )
})

test_that("preflight_de errors when fewer than 2 replicates per group (pseudobulk)", {
  obj <- make_mock_seurat()
  obj$one_rep <- "S1"   # every cell belongs to a single sample
  expect_error(
    preflight_de(obj, group_col = "group", ident1 = "group_A",
                 analysis_type = "pseudobulk", test_method = "DESeq2",
                 sample_col = "one_rep"),
    "at least 2 replicates"
  )
})

test_that("preflight_de passes on a well-formed mock object (singlecell)", {
  obj <- make_mock_seurat()
  expect_silent(
    preflight_de(obj, group_col = "group", ident1 = "group_A",
                 analysis_type = "singlecell", test_method = "wilcox")
  )
})

test_that("preflight_de passes on a well-formed mock object (pseudobulk)", {
  obj <- make_mock_seurat()
  expect_silent(
    preflight_de(obj, group_col = "group", ident1 = "group_A",
                 analysis_type = "pseudobulk", test_method = "DESeq2",
                 sample_col = "sample")
  )
})

test_that("normalize_de_output coerces names and adds effect sizes", {
  raw <- data.frame(
    p_val = c(0.01, 0.5),
    avg_log2FC = c(2, -0.1),
    pct.1 = c(0.8, 0.2),
    pct.2 = c(0.1, 0.3),
    p_val_adj = c(0.05, 0.7),
    row.names = c("Gene001", "Gene042")
  )
  out <- normalize_de_output(
    raw, method = "wilcox", seed = 42L,
    n_cells_group1 = 100L, n_cells_group2 = 100L
  )
  expect_true(all(c("gene", "avg_log2FC", "p_val", "p_val_adj",
                    "pct.1", "pct.2", "AUC", "cohens_d", "method",
                    "seed", "n_cells_group1", "n_cells_group2")
                  %in% colnames(out)))
  expect_equal(out$gene, c("Gene001", "Gene042"))
  expect_equal(out$method, c("wilcox", "wilcox"))
  expect_equal(out$seed, c(42L, 42L))
  expect_equal(out$n_cells_group1, c(100L, 100L))
  # AUC = (pct.1 + (1 - pct.2)) / 2
  expect_equal(out$AUC, c((0.8 + 0.9) / 2, (0.2 + 0.7) / 2))
})

test_that("run_de_analysis routes singlecell+wilcox correctly", {
  obj <- make_mock_seurat()
  out <- run_de_analysis(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    test_method = "wilcox", analysis_type = "singlecell", seed = 42L
  )
  expect_null(out$error)
  expect_s3_class(out$data, "data.frame")
  expect_equal(unique(out$data$method), "wilcox")
})

test_that("run_de_analysis routes pseudobulk+DESeq2 correctly", {
  skip_if_not_installed("DESeq2")
  obj <- make_mock_seurat()
  out <- run_de_analysis(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    test_method = "DESeq2", analysis_type = "pseudobulk",
    sample_col = "sample", seed = 42L
  )
  expect_null(out$error)
  expect_equal(unique(out$data$method), "DESeq2")
})

test_that("run_de_analysis returns a structured error on bad inputs", {
  obj <- make_mock_seurat()
  out <- run_de_analysis(
    obj, group_col = "nope", ident1 = "group_A",
    test_method = "wilcox", analysis_type = "singlecell"
  )
  expect_null(out$data)
  expect_match(out$error, "not found in metadata")
})

test_that("run_de_analysis subsets by cell-type column when subset_col is set", {
  obj <- make_mock_seurat()
  obj$ct <- ifelse(seq_len(ncol(obj)) %% 2 == 0, "ct_alpha", "ct_beta")
  out <- run_de_analysis(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    test_method = "wilcox", analysis_type = "singlecell",
    subset_col = "ct", subset_val = "ct_alpha", seed = 42L
  )
  expect_null(out$error)
  audit <- attr(out$data, "de_audit")
  expect_lt(audit$n_cells_total, 200L)   # subset reduced cell count
})

test_that("build_de_audit captures full reproducibility metadata (no downsample fields)", {
  audit <- build_de_audit(
    seed = 42L, test = "wilcox", group_col = "group",
    ident1 = "group_A", ident2 = "group_B",
    n_cells_total = 200L, n_group1 = 100L, n_group2 = 100L,
    latent_vars = NULL, n_genes_tested = 100L
  )
  expect_named(audit,
    c("seed", "test", "group_col", "ident1", "ident2",
      "n_cells_total", "n_group1", "n_group2",
      "latent_vars", "n_genes_tested",
      "timestamp", "seurat_version", "r_version"),
    ignore.order = TRUE
  )
  # Downsampling fields must NOT appear — DE always tests on the full dataset.
  expect_false("downsampled" %in% names(audit))
  expect_false("max_cells"   %in% names(audit))
  expect_equal(audit$test, "wilcox")
  expect_s3_class(audit$timestamp, "POSIXct")
})
