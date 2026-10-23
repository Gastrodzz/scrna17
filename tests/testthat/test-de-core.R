source("../../modules/de/de_core.R")

test_that("orgdb_for_species accepts Latin and lowercase forms", {
  expect_equal(orgdb_for_species("Homo sapiens"), "org.Hs.eg.db")
  expect_equal(orgdb_for_species("human"),         "org.Hs.eg.db")
  expect_equal(orgdb_for_species("Mus musculus"),  "org.Mm.eg.db")
  expect_equal(orgdb_for_species("MOUSE"),         "org.Mm.eg.db")
  expect_equal(orgdb_for_species("Rattus norvegicus"), "org.Rn.eg.db")
  expect_equal(orgdb_for_species("zebrafish"),     "org.Dr.eg.db")
})

test_that("orgdb_for_species errors with informative message on unknown species", {
  expect_error(orgdb_for_species("frog"), "Unsupported species")
})

test_that("kegg_code_for_species accepts Latin and lowercase forms", {
  expect_equal(kegg_code_for_species("Homo sapiens"), "hsa")
  expect_equal(kegg_code_for_species("mouse"),        "mmu")
  expect_equal(kegg_code_for_species("Rattus norvegicus"), "rno")
  expect_equal(kegg_code_for_species("Danio rerio"),  "dre")
})

test_that("detect_cell_type_col prefers canonical names then falls back", {
  obj <- make_mock_seurat()
  # No canonical name present yet -> falls back to the first categorical column
  # (Seurat puts 'orig.ident' first in metadata).
  expect_equal(detect_cell_type_col(obj), "orig.ident")
  # Adding a canonical name flips the result.
  obj$cell_type <- obj$group
  expect_equal(detect_cell_type_col(obj), "cell_type")
})

test_that("detect_sample_col prefers canonical names (orig.ident > sample)", {
  obj <- make_mock_seurat()
  # Seurat::CreateSeuratObject auto-creates orig.ident; it beats 'sample'.
  expect_equal(detect_sample_col(obj), "orig.ident")
  obj@meta.data$orig.ident <- NULL
  expect_equal(detect_sample_col(obj), "sample")
})

test_that("detect_sample_col returns NULL when no candidate exists", {
  obj <- make_mock_seurat()
  obj@meta.data$sample <- NULL
  obj@meta.data$orig.ident <- NULL
  expect_null(detect_sample_col(obj))
})
