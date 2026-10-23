# =============================================================================
# ondisk_storage.R — On-disk (BPCells) matrix storage to cut RAM
# =============================================================================
# Large scRNA-seq datasets keep their whole count matrix in RAM. BPCells stores
# the counts on disk (memory-mapped) and streams them, so RAM drops from
# O(whole matrix) to O(working set). Seurat v5 operates directly on the on-disk
# matrix. This module is written to be SAFE: every entry point degrades to a
# no-op (returns the object unchanged) if BPCells is missing or anything fails,
# so it can never break the in-memory path that already works.
# =============================================================================

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

#' Is a Seurat object's counts layer already backed on disk (BPCells)?
is_ondisk <- function(obj, assay = NULL) {
  if (is.null(obj)) return(FALSE)
  if (is.null(assay)) assay <- tryCatch(SeuratObject::DefaultAssay(obj), error = function(e) "RNA")
  m <- tryCatch(SeuratObject::LayerData(obj, layer = "counts", assay = assay),
                error = function(e) NULL)
  !is.null(m) && inherits(m, "IterableMatrix")
}

#' Number of cells in a Seurat object (cheap; used for the auto-threshold).
.n_cells <- function(obj) tryCatch(ncol(obj), error = function(e) 0L)

#' Convert a Seurat object's counts to on-disk BPCells storage.
#'
#' @param obj   Seurat object.
#' @param dir   Persistent directory to hold the on-disk matrix. MUST be a stable
#'   path (e.g. the project's session folder) because a saved object references
#'   it by path — a temp dir would be lost on restart.
#' @param assay Assay to convert (default = active assay).
#' @return the object with an on-disk counts layer, or the original object
#'   unchanged if BPCells is unavailable / conversion fails / already on disk.
store_seurat_ondisk <- function(obj, dir, assay = NULL, overwrite = TRUE) {
  if (is.null(obj)) return(obj)
  if (!requireNamespace("BPCells", quietly = TRUE)) return(obj)
  if (is.null(assay)) assay <- tryCatch(SeuratObject::DefaultAssay(obj),
                                        error = function(e) "RNA")
  if (is_ondisk(obj, assay)) return(obj)  # already on disk

  counts <- tryCatch(SeuratObject::LayerData(obj, layer = "counts", assay = assay),
                     error = function(e) NULL)
  if (is.null(counts) || !inherits(counts, c("dgCMatrix", "dgTMatrix", "matrix", "Matrix")))
    return(obj)

  cdir <- file.path(dir, paste0("bpcells_", assay))
  tryCatch({
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
    if (dir.exists(cdir)) {
      if (isTRUE(overwrite)) unlink(cdir, recursive = TRUE) else return(obj)
    }
    BPCells::write_matrix_dir(mat = counts, dir = cdir)
    obj[[assay]]$counts <- BPCells::open_matrix_dir(cdir)
    gc(FALSE)
    obj
  }, error = function(e) {
    warning("[on-disk] BPCells conversion failed (staying in-memory): ",
            conditionMessage(e))
    obj
  })
}

#' Apply on-disk storage automatically when a dataset is large enough to matter.
#' Returns the object unchanged for small datasets or when BPCells is missing.
#'
#' @param threshold cell count at/above which on-disk storage kicks in.
maybe_store_ondisk <- function(obj, dir, threshold = 50000L, assay = NULL) {
  if (is.null(obj)) return(obj)
  if (!requireNamespace("BPCells", quietly = TRUE)) return(obj)
  if (.n_cells(obj) < threshold) return(obj)
  store_seurat_ondisk(obj, dir = dir, assay = assay)
}

#' Is on-disk storage enabled? OFF unless env SCRNA_ONDISK is set (=1/true).
#' Kept off by default so the fully-in-RAM path (and its self-contained session
#' RDS save/restore) is unchanged unless a deployer opts in.
ondisk_enabled <- function() {
  v <- tolower(Sys.getenv("SCRNA_ONDISK", ""))
  nzchar(v) && v %in% c("1", "true", "yes", "on")
}

#' Convert a freshly loaded dataset to on-disk storage IFF enabled (env-gated).
#' Safe no-op when disabled, BPCells missing, small, or on any failure.
apply_ondisk_if_enabled <- function(obj, dir) {
  if (!ondisk_enabled()) return(obj)
  thr <- suppressWarnings(as.integer(Sys.getenv("SCRNA_ONDISK_THRESHOLD", "50000")))
  if (is.na(thr) || thr < 1L) thr <- 50000L
  maybe_store_ondisk(obj, dir = dir, threshold = thr)
}
