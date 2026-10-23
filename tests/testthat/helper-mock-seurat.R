# tests/testthat/helper-mock-seurat.R
# Builds a tiny deterministic Seurat v5 object for DE testing.
# Returns: 200 cells, 100 genes, 2 groups, 6 samples (3/group), 5 truly-DE genes.
make_mock_seurat <- function(seed = 1L) {
  set.seed(seed)
  n_cells_per_group <- 100L
  n_genes <- 100L
  n_samples_per_group <- 3L

  # Cells: 100 in group_A across samples A1/A2/A3, 100 in group_B across B1/B2/B3
  cell_ids <- paste0("cell_", seq_len(2 * n_cells_per_group))
  group <- rep(c("group_A", "group_B"), each = n_cells_per_group)
  sample <- c(
    rep(paste0("A", seq_len(n_samples_per_group)), length.out = n_cells_per_group),
    rep(paste0("B", seq_len(n_samples_per_group)), length.out = n_cells_per_group)
  )

  gene_ids <- paste0("Gene", sprintf("%03d", seq_len(n_genes)))
  de_up_in_A <- gene_ids[1:5]   # higher in group_A

  # Counts: Poisson background; DE genes get +5 mean shift in group_A
  counts <- matrix(rpois(n_genes * 2 * n_cells_per_group, lambda = 2),
                   nrow = n_genes, dimnames = list(gene_ids, cell_ids))
  shift <- matrix(0, nrow = n_genes, ncol = 2 * n_cells_per_group,
                  dimnames = list(gene_ids, cell_ids))
  shift[de_up_in_A, group == "group_A"] <- rpois(5 * n_cells_per_group, lambda = 8)
  counts <- counts + shift

  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj$group <- group
  obj$sample <- sample
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  attr(obj, "true_de_up_in_A") <- de_up_in_A
  obj
}
