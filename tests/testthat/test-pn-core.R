source("../../modules/pathway_network/pn_core.R")

test_that("pn_orgdb_for_species accepts Latin or lowercase", {
  expect_equal(pn_orgdb_for_species("Homo sapiens"), "org.Hs.eg.db")
  expect_equal(pn_orgdb_for_species("mouse"),        "org.Mm.eg.db")
  expect_equal(pn_orgdb_for_species("Rattus norvegicus"), "org.Rn.eg.db")
  expect_equal(pn_orgdb_for_species("zebrafish"),    "org.Dr.eg.db")
})

test_that("pn_kegg_code_for_species returns 3-letter codes for 4 species", {
  expect_equal(pn_kegg_code_for_species("human"),    "hsa")
  expect_equal(pn_kegg_code_for_species("Mus musculus"), "mmu")
  expect_equal(pn_kegg_code_for_species("rat"),      "rno")
  expect_equal(pn_kegg_code_for_species("Danio rerio"), "dre")
})

test_that("pn_msigdb_for_species returns Latin names", {
  expect_equal(pn_msigdb_for_species("human"), "Homo sapiens")
  expect_equal(pn_msigdb_for_species("Mouse"), "Mus musculus")
})

test_that("pn_taxid_for_species returns NCBI ids", {
  expect_equal(pn_taxid_for_species("human"), 9606L)
  expect_equal(pn_taxid_for_species("mouse"), 10090L)
})

test_that("pn_reactome_for_species returns lowercase canonical key", {
  expect_equal(pn_reactome_for_species("Homo sapiens"), "human")
  expect_equal(pn_reactome_for_species("MOUSE"),        "mouse")
})

test_that("build_ranked_list returns a sorted named numeric vector", {
  de <- data.frame(
    gene = paste0("G", 1:5),
    avg_log2FC = c(2, -1, 3, -0.5, 0.1),
    p_val = c(0.001, 0.01, 0.0001, 0.05, 0.5)
  )
  ranks <- build_ranked_list(de, metric = "signed_neglog10p", seed = 42L)
  expect_type(ranks, "double")
  expect_equal(length(ranks), 5)
  expect_true(all(diff(ranks) <= 0))   # sorted descending
  expect_setequal(names(ranks), de$gene)
})

test_that("build_ranked_list metric=log2fc uses fold change directly", {
  de <- data.frame(gene = c("A","B","C"), avg_log2FC = c(-1, 2, 0.5),
                   p_val = c(0.5, 0.001, 0.01))
  ranks <- build_ranked_list(de, metric = "log2fc", seed = 1L)
  # Top should be 'B' (FC=2), bottom 'A' (FC=-1)
  expect_equal(names(ranks)[1], "B")
  expect_equal(names(ranks)[length(ranks)], "A")
})

test_that("pn_sanitize_enrichment_df coerces factors and collapses list columns", {
  df <- data.frame(
    ID = factor(c("A","B")),
    pvalue = c(0.01, 0.02),
    Direction = c("Up","Down"),
    stringsAsFactors = TRUE
  )
  df$leadingEdge <- list(c("g1","g2","g3"), c("g4","g5"))
  out <- pn_sanitize_enrichment_df(df)
  expect_false(is.factor(out$ID))
  expect_type(out$leadingEdge, "character")
  expect_equal(out$leadingEdge[1], "g1, g2, g3")
})

test_that("pn_background_universe prefers custom, falls back to seurat rownames", {
  fake_obj <- list()
  attr(fake_obj, "fake") <- TRUE
  fake_obj <- structure(list(),
                        class = "Seurat",
                        fake = TRUE)
  # Use a real-looking object: just a list with a rownames<- method override is hard;
  # instead test the custom-override branch and the NULL-fallback branch.
  expect_equal(pn_background_universe(custom = c("x","y","y","z")), c("x","y","z"))
  expect_null(pn_background_universe())
})
