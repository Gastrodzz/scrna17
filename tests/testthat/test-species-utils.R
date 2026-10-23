source(file.path("..", "..", "modules", "species_utils.R"))

test_that("human symbols detect as human with MT pattern ^MT-", {
  genes <- c("MT-ND1","MT-CO1","PTPRC","CD3D","CD8A","MS4A1","NKG7","EPCAM","GAPDH")
  r <- detect_species(genes)
  expect_equal(r$species, "human")
  expect_equal(r$mt_pattern, "^MT-")
})

test_that("mouse symbols detect as mouse with MT pattern ^mt-", {
  genes <- c("mt-Nd1","mt-Co1","Ptprc","Cd3d","Cd8a","Ms4a1","Nkg7","Epcam","Gapdh")
  r <- detect_species(genes)
  expect_equal(r$species, "mouse")
  expect_equal(r$mt_pattern, "^mt-")
})

test_that("ensembl id matrix is flagged unknown/low", {
  genes <- paste0("ENSG000000", sprintf("%05d", 1:50))
  r <- detect_species(genes)
  expect_equal(r$species, "unknown")
  expect_match(r$note, "Ensembl")
})

test_that("normalize_species maps latin and short forms", {
  expect_equal(normalize_species("Homo sapiens"), "human")
  expect_equal(normalize_species("Mus musculus"), "mouse")
  expect_equal(normalize_species("mouse"), "mouse")
  expect_equal(normalize_species("Rattus norvegicus"), "rat")
  expect_equal(normalize_species("Other"), "")
})

test_that("species_patterns returns human patterns for unknown/blank", {
  expect_equal(species_patterns("")$mt, "^MT-")
  expect_equal(species_patterns("Mus musculus")$mt, "^mt-")
})
