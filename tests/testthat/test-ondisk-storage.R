suppressWarnings(suppressMessages({
  library(Seurat); library(Matrix)
}))
source(file.path("..", "..", "modules", "ondisk_storage.R"))

make_obj <- function(nc = 60, ng = 300) {
  set.seed(1)
  m <- Matrix::rsparsematrix(ng, nc, density = 0.15)
  m@x <- abs(round(m@x * 10)) + 1
  rownames(m) <- paste0("Gene", seq_len(ng))
  colnames(m) <- paste0("Cell", seq_len(nc))
  suppressWarnings(CreateSeuratObject(counts = m))
}

test_that("small object is unchanged by maybe_store_ondisk (below threshold)", {
  skip_if_not(requireNamespace("BPCells", quietly = TRUE), "BPCells not installed")
  so <- make_obj(nc = 60)
  out <- maybe_store_ondisk(so, dir = tempfile("od_"), threshold = 50000L)
  expect_false(is_ondisk(out))                       # stays in-memory
  expect_equal(ncol(out), 60)
})

test_that("store_seurat_ondisk moves counts to disk and preserves values", {
  skip_if_not(requireNamespace("BPCells", quietly = TRUE), "BPCells not installed")
  so <- make_obj(nc = 60, ng = 300)
  orig <- as.matrix(SeuratObject::LayerData(so, layer = "counts", assay = "RNA"))
  dir <- tempfile("od_")
  out <- store_seurat_ondisk(so, dir = dir)
  expect_true(is_ondisk(out))                        # now on disk
  on_disk <- as.matrix(SeuratObject::LayerData(out, layer = "counts", assay = "RNA"))
  expect_true(all(on_disk == orig))                  # values identical
})

test_that("normalization works on an on-disk object", {
  skip_if_not(requireNamespace("BPCells", quietly = TRUE), "BPCells not installed")
  so <- store_seurat_ondisk(make_obj(nc = 60), dir = tempfile("od_"))
  so <- suppressWarnings(NormalizeData(so, verbose = FALSE))
  expect_true("data" %in% SeuratObject::Layers(so[["RNA"]]))
})

test_that("converting an already on-disk object is a no-op (idempotent)", {
  skip_if_not(requireNamespace("BPCells", quietly = TRUE), "BPCells not installed")
  so <- store_seurat_ondisk(make_obj(nc = 60), dir = tempfile("od_"))
  again <- store_seurat_ondisk(so, dir = tempfile("od2_"))
  expect_true(is_ondisk(again))
})
