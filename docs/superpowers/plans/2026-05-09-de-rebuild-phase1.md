# DE Backend Rebuild — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the 2 227-line monolithic `modules/downstream/downstream_de.R` with a 4-file modular `modules/de/` backend that runs all 9 supported DE methods (6 single-cell + 3 pseudobulk) reproducibly, with audit metadata, on a 25 k-cell dataset.

**Architecture:** One dispatcher (`run_de_analysis`) calls one of two engines (`run_single_cell_de` via Seurat `FindMarkers`; `run_pseudobulk_de` via aggregated counts → DESeq2/edgeR/limma-voom). All methods return the same data-frame schema and an attached `attr(result, "de_audit")` list. Pre-flight validation hard-stops on missing replicates or empty groups. Seed is propagated from UI to every randomized step. UI is split-screen (340 px config + flex results) per spec §7. AI controls and unrelated sub-tabs are removed from the Downstream Analysis page.

**Tech Stack:** R 4.x, Seurat v5 (layer-safe), Bioconductor (MAST, DESeq2, edgeR, limma), Shiny + shinydashboard, plotly, DT, testthat.

**Project context (read before starting):**
- This project is **not a git repository** (`git rev-parse` confirms). Replace any “commit” checkpoint with a logical save-point — write a one-line note to `docs/superpowers/plans/_progress.md` so the executor can resume.
- Working directory is `/Users/midou/Desktop/scrna15`. All paths below are relative to it.
- Spec lives at [docs/superpowers/specs/2026-05-09-de-pathway-network-rebuild-design.md](docs/superpowers/specs/2026-05-09-de-pathway-network-rebuild-design.md) — sections referenced by `§N`.
- `values$active_dataset` is a Seurat v5 object; `values$species_global` is the **Latin name** (`"Homo sapiens"`, `"Mus musculus"`, `"Rattus norvegicus"`, `"Danio rerio"`); `values$last_de_results` is the canonical DE slot read by `modules/export.R` (8 sites) — **must continue to be written by this module**.
- Helpers `prepare_for_de()` and `safe_get_assay_data()` live in [modules/core_robustness.R](modules/core_robustness.R); reuse them, do not redefine.
- During Phase 1 the existing `pathwayServer`, `cellCommServer`, `tfActivityServer`, `grnServer`, `geneSetScoringServer`, `compositionServer`, `aiAnalysisServer` files stay on disk untouched. Only their UI/server **calls** disappear from `app.R`. They will be re-mounted under the new `pathway_network` module in Phase 2/3.

**Out of scope (later phases):**
- New `modules/pathway_network/` directory — Phase 2.
- Deletion of `modules/network_analysis.R` and merging GRN/co-expression/pathway-network — Phase 3.
- Reintroducing AI Auto-DE — deferred indefinitely per user instruction.

---

## File Structure

**New files (created by this plan):**

| Path | Responsibility |
|---|---|
| `modules/de/de_core.R` | Pure helpers: species → OrgDb / KEGG code, metadata-column detection, output normalization, audit-log builder. No Shiny. |
| `modules/de/de_methods.R` | Pre-flight validators + dispatcher `run_de_analysis()` + engines `run_single_cell_de()` / `run_pseudobulk_de()`. No Shiny. |
| `modules/de/de_ui.R` | `deUI(id)` — split-screen config (340 px) + results panel. Method pills, threshold inputs, Audit Log tab. |
| `modules/de/de_server.R` | `deServer(id, values)` — reactive wiring, output renderers, sources `de_core.R` + `de_methods.R` at the top. |
| `tests/testthat/helper-mock-seurat.R` | testthat helper: builds a tiny deterministic Seurat object with two groups (≥ 3 replicates each) and known DE genes. |
| `tests/testthat/test-de-core.R` | Unit tests for `de_core.R` helpers. |
| `tests/testthat/test-de-methods.R` | Unit tests for pre-flight validation, dispatcher routing, output schema, audit log, reproducibility. |
| `tests/testthat/test-de-engines.R` | End-to-end tests calling `run_de_analysis()` once per method on the mock object. |
| `docs/superpowers/plans/_progress.md` | Single-line progress notes (replaces git commits since repo is not under git). |

**Modified files:**

| Path | Change |
|---|---|
| [app.R](app.R) | (1) Add `safe_source("modules/de/de_server.R")` to the `modules_loaded` block (line ~88). (2) Replace the `tabsetPanel` body of the Downstream tab (lines ~709-750) with a single `deUI("de")`. (3) Remove `pathwayUI`, `cellCommUI`, `tfActivityUI`, `grnUI`, `geneSetScoringUI`, `compositionUI`, `aiAnalysisUI` calls from that tab. (4) In the lazy `.downstream_init` observer (lines ~1782-1795), keep only `deServer("de", values)` — drop the other server calls. |
| [modules/downstream/downstream_main.R](modules/downstream/downstream_main.R) | Remove `"downstream_de.R"` from the `sub_modules` vector (line ~55). |
| [modules/downstream/downstream_de.R](modules/downstream/downstream_de.R) | **Deleted in final task.** |

**Untouched files (sourced verbatim by the new module):**
- [modules/core_robustness.R](modules/core_robustness.R) — provides `prepare_for_de()`, `safe_get_assay_data()`.
- [modules/downstream/downstream_core.R](modules/downstream/downstream_core.R) — provides `.sanitize_enrichment_df()` (not used here, leave alone).

---

## Conventions for Every Task

- **TDD discipline:** Write the failing test, run it to confirm the failure mode, write the minimal code to pass, run again to confirm green, then checkpoint. Tests live in `tests/testthat/`.
- **Run tests with:** `Rscript -e "testthat::test_file('tests/testthat/test-<name>.R', reporter='summary')"` from the project root. Expected output for a passing file ends with `[ FAIL 0 | WARN 0 | SKIP 0 | PASS N ]`.
- **Run the full suite with:** `Rscript -e "testthat::test_dir('tests/testthat', reporter='summary')"`.
- **Checkpoint syntax (replaces git commit):** append a single dated line to `docs/superpowers/plans/_progress.md`, e.g. `2026-05-09 task 3 done — pre-flight validators tested`.
- **No hard-coded values** outside UI defaults — every column name, threshold, species code, file path is derived (spec §3 + §9.3).
- **Seurat v5 only.** Always go through `prepare_for_de()` before reading data. Never touch `@data` / `@counts` directly.
- **Latin-OR-lowercase species:** every helper that takes a species string must accept either `"Homo sapiens"` or `"human"` (case-insensitive). The mapping helper normalises internally.

---

## Task 1: Test fixtures — mock Seurat builder

**Files:**
- Create: `tests/testthat/helper-mock-seurat.R`

testthat auto-sources files matching `helper-*.R` before any test file. The fixture must produce a deterministic Seurat v5 object with two groups (`group_A`, `group_B`), six samples (three per group — meets pseudobulk’s ≥ 2 replicates requirement), and a known set of 5 strongly DE genes.

- [ ] **Step 1: Write the helper**

```r
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
  background <- setdiff(gene_ids, de_up_in_A)

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
```

- [ ] **Step 2: Sanity-run the helper**

Run: `Rscript -e "source('tests/testthat/helper-mock-seurat.R'); o <- make_mock_seurat(); cat(ncol(o), 'cells,', nrow(o), 'genes,', length(unique(o\$sample)), 'samples\n')"`
Expected: `200 cells, 100 genes, 6 samples`

- [ ] **Step 3: Checkpoint**

Append to `docs/superpowers/plans/_progress.md`:
`2026-05-09 task 1 done — mock seurat builder ready`

---

## Task 2: `de_core.R` — species + metadata helpers (TDD)

**Files:**
- Create: `modules/de/de_core.R`
- Test: `tests/testthat/test-de-core.R`

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-de-core.R
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
  expect_equal(detect_cell_type_col(obj), "group")  # 'group' is the only categorical
  obj$cell_type <- obj$group
  expect_equal(detect_cell_type_col(obj), "cell_type")
})

test_that("detect_sample_col prefers canonical names", {
  obj <- make_mock_seurat()
  expect_equal(detect_sample_col(obj), "sample")
  obj$orig.ident <- obj$sample
  expect_equal(detect_sample_col(obj), "orig.ident")
})

test_that("detect_sample_col returns NULL when no candidate exists", {
  obj <- make_mock_seurat()
  obj@meta.data$sample <- NULL
  obj@meta.data$orig.ident <- NULL
  expect_null(detect_sample_col(obj))
})
```

- [ ] **Step 2: Run tests; confirm they fail with `could not find function`**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-core.R', reporter='summary')"`
Expected: every test errors with `could not find function "orgdb_for_species"` (or similar). 0 PASS.

- [ ] **Step 3: Write the minimal implementation**

```r
# modules/de/de_core.R
# Pure helpers shared by de_methods.R, de_server.R. No Shiny imports.

#' Map a species name (Latin or lowercase) to its OrgDb package name.
#' @param species character(1). Examples: "Homo sapiens", "human", "MOUSE".
#' @return character(1). Package name like "org.Hs.eg.db".
orgdb_for_species <- function(species) {
  key <- .normalize_species_key(species)
  map <- c(
    human     = "org.Hs.eg.db",
    mouse     = "org.Mm.eg.db",
    rat       = "org.Rn.eg.db",
    zebrafish = "org.Dr.eg.db"
  )
  pkg <- map[key]
  if (is.na(pkg)) stop("Unsupported species: ", species,
                       ". Supported: human, mouse, rat, zebrafish.")
  unname(pkg)
}

#' Map a species name to its KEGG 3-letter organism code.
kegg_code_for_species <- function(species) {
  key <- .normalize_species_key(species)
  map <- c(human = "hsa", mouse = "mmu", rat = "rno", zebrafish = "dre")
  code <- map[key]
  if (is.na(code)) stop("No KEGG code for species: ", species)
  unname(code)
}

#' Internal: collapse Latin or English species name to a canonical lowercase key.
.normalize_species_key <- function(species) {
  if (is.null(species) || !nzchar(species)) return(NA_character_)
  s <- tolower(trimws(species))
  if (s %in% c("human", "homo sapiens", "h. sapiens", "hs"))      return("human")
  if (s %in% c("mouse", "mus musculus", "m. musculus", "mm"))     return("mouse")
  if (s %in% c("rat", "rattus norvegicus", "r. norvegicus", "rn")) return("rat")
  if (s %in% c("zebrafish", "danio rerio", "d. rerio", "dr"))     return("zebrafish")
  NA_character_
}

#' Pick the most plausible cell-type / cluster column from a Seurat object.
#' @return character(1) or NULL.
detect_cell_type_col <- function(obj) {
  candidates <- c("cell_type", "celltype", "cell.type", "annotation",
                  "seurat_clusters", "cluster")
  cols <- colnames(obj@meta.data)
  found <- intersect(candidates, cols)
  if (length(found) > 0) return(found[1])
  cat_cols <- cols[vapply(obj@meta.data,
                          function(x) is.factor(x) || is.character(x),
                          logical(1))]
  if (length(cat_cols) > 0) cat_cols[1] else NULL
}

#' Pick the most plausible sample / replicate column. Returns NULL if none found.
detect_sample_col <- function(obj) {
  candidates <- c("orig.ident", "sample", "sample_id", "donor",
                  "patient", "replicate", "batch")
  found <- intersect(candidates, colnames(obj@meta.data))
  if (length(found) > 0) found[1] else NULL
}
```

- [ ] **Step 4: Run tests; confirm green**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-core.R', reporter='summary')"`
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 8 ]`.

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 2 done — de_core.R helpers + tests green`

---

## Task 3: `de_methods.R` — pre-flight validators (TDD)

**Files:**
- Create: `modules/de/de_methods.R` (start with the validators)
- Test: `tests/testthat/test-de-methods.R`

The validators stop bad inputs **before** any expensive work runs. Per spec §4.3:
- `sample_col` exists in metadata (pseudobulk only).
- ≥ 2 replicates per group (pseudobulk only).
- ≥ 10 cells per group (single-cell).
- `data` layer present — call `NormalizeData()` if only `counts`.
- Method package is installed.

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-de-methods.R
source("../../modules/de/de_core.R")
source("../../modules/de/de_methods.R")

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
```

- [ ] **Step 2: Run; confirm failures**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-methods.R', reporter='summary')"`
Expected: 6 errors `could not find function "preflight_de"`.

- [ ] **Step 3: Write the minimal implementation**

```r
# modules/de/de_methods.R
# DE engines and dispatcher. Sourced before de_server.R.
# Depends on: modules/de/de_core.R, modules/core_robustness.R (prepare_for_de).

# --- Constants ---------------------------------------------------------------
.SC_METHODS <- c("wilcox", "MAST", "t", "bimod", "LR", "negbinom")
.PB_METHODS <- c("DESeq2", "edgeR", "limma-voom")
.MIN_CELLS_SC <- 10L
.MIN_REPS_PB <- 2L
.METHOD_PKG <- list(
  MAST       = "MAST",
  DESeq2     = "DESeq2",
  edgeR      = "edgeR",
  `limma-voom` = "limma",
  negbinom   = "MASS",
  LR         = "stats",
  wilcox     = NULL, t = NULL, bimod = NULL
)

# --- Pre-flight validation ---------------------------------------------------
preflight_de <- function(obj, group_col, ident1, ident2 = NULL,
                         analysis_type = c("singlecell", "pseudobulk"),
                         test_method, sample_col = NULL,
                         min_cells = .MIN_CELLS_SC,
                         min_reps  = .MIN_REPS_PB) {
  analysis_type <- match.arg(analysis_type)
  if (!group_col %in% colnames(obj@meta.data))
    stop("Group column '", group_col, "' not found in metadata.")
  if (!ident1 %in% as.character(obj@meta.data[[group_col]]))
    stop("ident1 '", ident1, "' not found in column '", group_col, "'.")

  # Method-specific package check
  pkg <- .METHOD_PKG[[test_method]]
  if (!is.null(pkg) && !requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' required for test '", test_method,
         "'. Install: BiocManager::install('", pkg, "').")

  # Cell-count check (single-cell)
  if (analysis_type == "singlecell") {
    n1 <- sum(obj@meta.data[[group_col]] == ident1)
    if (n1 < min_cells)
      stop("Group '", ident1, "' has fewer than ", min_cells,
           " cells (n=", n1, "). Use pseudobulk or pick another group.")
    if (!is.null(ident2) && !ident2 %in% c("rest", "all", "")) {
      n2 <- sum(obj@meta.data[[group_col]] == ident2)
      if (n2 < min_cells)
        stop("Group '", ident2, "' has fewer than ", min_cells, " cells (n=", n2, ").")
    }
  }

  # Replicate check (pseudobulk)
  if (analysis_type == "pseudobulk") {
    if (is.null(sample_col) || !nzchar(sample_col))
      stop("sample_col is required for pseudobulk DE.")
    if (!sample_col %in% colnames(obj@meta.data))
      stop("Sample column '", sample_col, "' not found in metadata.")

    md <- obj@meta.data
    s1 <- unique(md[[sample_col]][md[[group_col]] == ident1])
    if (length(s1) < min_reps)
      stop("Pseudobulk requires at least ", min_reps, " replicates per group; ",
           "group '", ident1, "' has ", length(s1), ".")
    if (!is.null(ident2) && !ident2 %in% c("rest", "all", "")) {
      s2 <- unique(md[[sample_col]][md[[group_col]] == ident2])
      if (length(s2) < min_reps)
        stop("Pseudobulk requires at least ", min_reps, " replicates per group; ",
             "group '", ident2, "' has ", length(s2), ".")
    }
  }

  invisible(TRUE)
}
```

- [ ] **Step 4: Run; confirm green**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-methods.R', reporter='summary')"`
Expected: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]`.

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 3 done — pre-flight validators tested`

---

## Task 4: Output normalisation + audit log (TDD)

**Files:**
- Modify: `modules/de/de_methods.R` (append helpers)
- Modify: `tests/testthat/test-de-methods.R` (append tests)

Goal: every engine returns the schema specified in §4.1: `gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, AUC, cohens_d, method, seed, n_cells_group1, n_cells_group2`. The audit list (§4.4) is attached as `attr(result, "de_audit")`.

- [ ] **Step 1: Append failing tests**

```r
# Append to tests/testthat/test-de-methods.R
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

test_that("build_de_audit captures full reproducibility metadata", {
  audit <- build_de_audit(
    seed = 42L, test = "wilcox", group_col = "group",
    ident1 = "group_A", ident2 = "group_B",
    n_cells_total = 200L, n_group1 = 100L, n_group2 = 100L,
    downsampled = FALSE, max_cells = Inf, latent_vars = NULL,
    n_genes_tested = 100L
  )
  expect_named(audit,
    c("seed", "test", "group_col", "ident1", "ident2",
      "n_cells_total", "n_group1", "n_group2",
      "downsampled", "max_cells", "latent_vars", "n_genes_tested",
      "timestamp", "seurat_version", "r_version"),
    ignore.order = TRUE
  )
  expect_equal(audit$test, "wilcox")
  expect_s3_class(audit$timestamp, "POSIXct")
})
```

- [ ] **Step 2: Run; confirm failures**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-methods.R', reporter='summary')"`
Expected: 2 new failures (`could not find function "normalize_de_output"` / `"build_de_audit"`).

- [ ] **Step 3: Append the implementation**

```r
# Append to modules/de/de_methods.R

#' Normalise an engine's raw output to the canonical DE schema.
#' @param raw data.frame with at least p_val, avg_log2FC, pct.1, pct.2 and rownames=gene.
normalize_de_output <- function(raw, method, seed,
                                n_cells_group1, n_cells_group2) {
  if (is.null(raw) || nrow(raw) == 0) {
    return(.empty_de_df())
  }
  cn <- colnames(raw)
  fc_col <- if ("avg_log2FC" %in% cn) "avg_log2FC"
            else if ("avg_logFC" %in% cn) "avg_logFC"
            else if ("logFC" %in% cn)     "logFC"
            else NA_character_
  p_col   <- if ("p_val"   %in% cn) "p_val"
             else if ("PValue" %in% cn) "PValue"
             else if ("pvalue" %in% cn) "pvalue"
             else NA_character_
  padj_col <- if ("p_val_adj" %in% cn) "p_val_adj"
              else if ("FDR" %in% cn)  "FDR"
              else if ("padj" %in% cn) "padj"
              else NA_character_

  gene <- if ("gene" %in% cn) as.character(raw$gene) else rownames(raw)
  avg_log2FC <- if (!is.na(fc_col)) as.numeric(raw[[fc_col]]) else rep(NA_real_, nrow(raw))
  p_val      <- if (!is.na(p_col))   as.numeric(raw[[p_col]])   else rep(NA_real_, nrow(raw))
  p_val_adj  <- if (!is.na(padj_col)) as.numeric(raw[[padj_col]])
                else stats::p.adjust(p_val, method = "BH")
  pct.1 <- if ("pct.1" %in% cn) as.numeric(raw$pct.1) else rep(NA_real_, nrow(raw))
  pct.2 <- if ("pct.2" %in% cn) as.numeric(raw$pct.2) else rep(NA_real_, nrow(raw))

  AUC <- ifelse(is.na(pct.1) | is.na(pct.2), NA_real_, (pct.1 + (1 - pct.2)) / 2)
  cohens_d <- avg_log2FC / log2(exp(1))   # log2 → natural log

  data.frame(
    gene = gene,
    avg_log2FC = avg_log2FC,
    p_val = p_val,
    p_val_adj = p_val_adj,
    pct.1 = pct.1,
    pct.2 = pct.2,
    AUC = AUC,
    cohens_d = cohens_d,
    method = rep(method, length(gene)),
    seed = rep(as.integer(seed), length(gene)),
    n_cells_group1 = rep(as.integer(n_cells_group1), length(gene)),
    n_cells_group2 = rep(as.integer(n_cells_group2), length(gene)),
    stringsAsFactors = FALSE
  )
}

.empty_de_df <- function() {
  data.frame(
    gene = character(0), avg_log2FC = numeric(0),
    p_val = numeric(0), p_val_adj = numeric(0),
    pct.1 = numeric(0), pct.2 = numeric(0),
    AUC = numeric(0), cohens_d = numeric(0),
    method = character(0), seed = integer(0),
    n_cells_group1 = integer(0), n_cells_group2 = integer(0),
    stringsAsFactors = FALSE
  )
}

#' Build the audit-log list attached to every DE result.
build_de_audit <- function(seed, test, group_col, ident1, ident2,
                           n_cells_total, n_group1, n_group2,
                           downsampled, max_cells, latent_vars,
                           n_genes_tested) {
  list(
    seed           = as.integer(seed),
    test           = test,
    group_col      = group_col,
    ident1         = ident1,
    ident2         = if (is.null(ident2)) NA_character_ else ident2,
    n_cells_total  = as.integer(n_cells_total),
    n_group1       = as.integer(n_group1),
    n_group2       = as.integer(n_group2),
    downsampled    = isTRUE(downsampled),
    max_cells      = if (is.finite(max_cells)) as.integer(max_cells) else Inf,
    latent_vars    = if (is.null(latent_vars)) character(0) else as.character(latent_vars),
    n_genes_tested = as.integer(n_genes_tested),
    timestamp      = Sys.time(),
    seurat_version = as.character(utils::packageVersion("Seurat")),
    r_version      = R.version$version.string
  )
}
```

- [ ] **Step 4: Run; confirm green**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-methods.R', reporter='summary')"`
Expected: 8 PASS total (6 from Task 3 + 2 new).

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 4 done — output normaliser + audit log tested`

---

## Task 5: Single-cell engine (TDD)

**Files:**
- Modify: `modules/de/de_methods.R` (append `run_single_cell_de`)
- Create: `tests/testthat/test-de-engines.R`

The single-cell engine wraps `Seurat::FindMarkers()` for the 6 supported methods. Per spec §4.1: full dataset by default, deterministic seed, MAST gets CDR added to `latent.vars` automatically. Returns the canonical schema via `normalize_de_output()`.

- [ ] **Step 1: Write the failing tests**

```r
# tests/testthat/test-de-engines.R
source("../../modules/de/de_core.R")
source("../../modules/de/de_methods.R")
source("../../modules/core_robustness.R")  # prepare_for_de, safe_get_assay_data

test_that("run_single_cell_de(wilcox) returns canonical schema with non-zero rows", {
  obj <- make_mock_seurat()
  res <- run_single_cell_de(
    obj, group_col = "group", ident1 = "group_A", ident2 = "group_B",
    test_method = "wilcox", seed = 42L, max_cells_per_ident = Inf,
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
    test_method = "wilcox", seed = 42L, max_cells_per_ident = Inf,
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
                          max_cells_per_ident = Inf, latent_vars = NULL,
                          min_pct_prefilter = 0.01, logfc_threshold = 0)
  b <- run_single_cell_de(obj, "group", "group_A", "group_B",
                          test_method = "wilcox", seed = 7L,
                          max_cells_per_ident = Inf, latent_vars = NULL,
                          min_pct_prefilter = 0.01, logfc_threshold = 0)
  expect_equal(a$gene, b$gene)
  expect_equal(a$p_val, b$p_val)
  expect_equal(a$avg_log2FC, b$avg_log2FC)
})

test_that("run_single_cell_de attaches a complete audit log", {
  obj <- make_mock_seurat()
  res <- run_single_cell_de(obj, "group", "group_A", "group_B",
                            test_method = "wilcox", seed = 42L,
                            max_cells_per_ident = Inf, latent_vars = NULL,
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
                            max_cells_per_ident = Inf, latent_vars = NULL,
                            min_pct_prefilter = 0.01, logfc_threshold = 0)
  audit <- attr(res, "de_audit")
  expect_true("CDR" %in% audit$latent_vars)
})
```

- [ ] **Step 2: Run; confirm failures**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-engines.R', reporter='summary')"`
Expected: 5 failures `could not find function "run_single_cell_de"`.

- [ ] **Step 3: Append the implementation**

```r
# Append to modules/de/de_methods.R

#' Single-cell DE via Seurat::FindMarkers.
#' Methods: wilcox, MAST, t, bimod, LR, negbinom.
run_single_cell_de <- function(obj, group_col, ident1, ident2 = NULL,
                               test_method = "wilcox",
                               seed = 42L,
                               max_cells_per_ident = Inf,
                               latent_vars = NULL,
                               min_pct_prefilter = 0.01,
                               logfc_threshold = 0,
                               only_pos = FALSE,
                               min_pct_display = 0) {
  if (!test_method %in% .SC_METHODS)
    stop("Unsupported single-cell method: ", test_method,
         ". Supported: ", paste(.SC_METHODS, collapse = ", "))

  obj <- prepare_for_de(obj)
  Seurat::Idents(obj) <- group_col
  set.seed(seed)

  # MAST: auto-add CDR (cellular detection rate) as a latent variable
  if (test_method == "MAST") {
    expr <- safe_get_assay_data(obj, "data")
    obj$CDR <- as.numeric(scale(colMeans(expr > 0)))
    latent_vars <- unique(c(latent_vars, "CDR"))
  }

  ident2_arg <- if (is.null(ident2) || ident2 %in% c("rest", "")) NULL else ident2
  cells1 <- which(as.character(obj@meta.data[[group_col]]) == ident1)
  cells2 <- if (is.null(ident2_arg)) {
    which(as.character(obj@meta.data[[group_col]]) != ident1)
  } else {
    which(as.character(obj@meta.data[[group_col]]) == ident2_arg)
  }
  n1 <- length(cells1); n2 <- length(cells2)

  raw <- Seurat::FindMarkers(
    obj, ident.1 = ident1, ident.2 = ident2_arg,
    test.use = test_method,
    min.pct = min_pct_prefilter,
    logfc.threshold = 0,
    min.diff.pct = -Inf,
    only.pos = FALSE,
    latent.vars = latent_vars,
    max.cells.per.ident = max_cells_per_ident,
    random.seed = seed,
    verbose = FALSE
  )
  raw$gene <- rownames(raw)

  out <- normalize_de_output(raw, method = test_method, seed = seed,
                             n_cells_group1 = n1, n_cells_group2 = n2)

  # Display-only post-filters
  if (isTRUE(only_pos))     out <- out[out$avg_log2FC > 0, , drop = FALSE]
  if (logfc_threshold > 0)  out <- out[abs(out$avg_log2FC) >= logfc_threshold, , drop = FALSE]
  if (min_pct_display > 0)
    out <- out[pmax(out$pct.1, out$pct.2, na.rm = TRUE) >= min_pct_display, , drop = FALSE]

  out <- out[order(out$p_val_adj, -abs(out$avg_log2FC)), ]

  attr(out, "de_audit") <- build_de_audit(
    seed = seed, test = test_method, group_col = group_col,
    ident1 = ident1, ident2 = ident2,
    n_cells_total = ncol(obj), n_group1 = n1, n_group2 = n2,
    downsampled = is.finite(max_cells_per_ident),
    max_cells = max_cells_per_ident,
    latent_vars = latent_vars,
    n_genes_tested = nrow(out)
  )
  out
}
```

- [ ] **Step 4: Run; confirm green**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-engines.R', reporter='summary')"`
Expected: 5 PASS (or 4 PASS + 1 SKIP if MAST is not installed; that is acceptable).

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 5 done — single-cell engine 6 methods tested`

---

## Task 6: Pseudobulk engine — DESeq2 + edgeR + limma-voom (TDD)

**Files:**
- Modify: `modules/de/de_methods.R` (append `run_pseudobulk_de` and helpers)
- Modify: `tests/testthat/test-de-engines.R` (append tests)

Aggregate counts per `sample_col`, run the chosen tool, normalise output. **`limma-voom` is new** — Seurat doesn’t expose it; we implement it directly with `limma::voom()` + `lmFit()` + `topTable()`.

- [ ] **Step 1: Append failing tests**

```r
# Append to tests/testthat/test-de-engines.R
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

test_that("run_pseudobulk_de attaches an audit log with sample counts", {
  skip_if_not_installed("DESeq2")
  obj <- make_mock_seurat()
  res <- run_pseudobulk_de(obj, "group", "group_A", "group_B",
                           sample_col = "sample", method = "DESeq2", seed = 42L)
  audit <- attr(res, "de_audit")
  expect_equal(audit$test, "DESeq2")
  expect_equal(audit$n_group1, 100L)   # cell counts
  expect_equal(audit$n_group2, 100L)
})
```

- [ ] **Step 2: Run; confirm failures**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-engines.R', reporter='summary')"`
Expected: 4 new failures (`could not find function "run_pseudobulk_de"`).

- [ ] **Step 3: Append the implementation**

```r
# Append to modules/de/de_methods.R

#' Pseudobulk DE: aggregate counts per sample, run DESeq2/edgeR/limma-voom.
#' @param method one of "DESeq2", "edgeR", "limma-voom".
run_pseudobulk_de <- function(obj, group_col, ident1, ident2 = NULL,
                              sample_col, method = "DESeq2", seed = 42L) {
  if (!method %in% .PB_METHODS)
    stop("Unsupported pseudobulk method: ", method,
         ". Supported: ", paste(.PB_METHODS, collapse = ", "))

  obj <- prepare_for_de(obj)
  set.seed(seed)

  md <- obj@meta.data
  cells1 <- rownames(md)[md[[group_col]] == ident1]
  ident2_label <- if (is.null(ident2) || ident2 %in% c("rest", "all", ""))
                    "rest" else ident2
  cells2 <- if (ident2_label == "rest")
              rownames(md)[md[[group_col]] != ident1]
            else
              rownames(md)[md[[group_col]] == ident2_label]

  n1_cells <- length(cells1); n2_cells <- length(cells2)
  s1 <- unique(md[cells1, sample_col]); s2 <- unique(md[cells2, sample_col])

  counts_mat <- safe_get_assay_data(obj, "counts")
  agg <- function(samples, group) {
    do.call(cbind, lapply(samples, function(s) {
      cells <- rownames(md)[md[[sample_col]] == s &
                            ((group == ident1 & md[[group_col]] == ident1) |
                             (group != ident1 & md[[group_col]] != ident1))]
      if (length(cells) == 0) return(NULL)
      Matrix::rowSums(counts_mat[, cells, drop = FALSE])
    }))
  }
  m1 <- agg(s1, ident1); m2 <- agg(s2, "other")
  if (is.null(m1) || is.null(m2))
    return(.empty_de_df_with_audit(method, seed, n1_cells, n2_cells,
                                   group_col, ident1, ident2,
                                   ncol(obj)))
  bulk <- cbind(m1, m2)
  colnames(bulk) <- c(s1, s2)
  group_factor <- factor(c(rep(ident1, length(s1)),
                           rep(ident2_label, length(s2))),
                         levels = c(ident2_label, ident1))
  col_data <- data.frame(group = group_factor, row.names = colnames(bulk))

  raw <- switch(method,
    "DESeq2"     = .pb_deseq2(bulk, col_data, ident1, ident2_label),
    "edgeR"      = .pb_edger (bulk, col_data, ident1, ident2_label),
    "limma-voom" = .pb_voom  (bulk, col_data, ident1, ident2_label)
  )

  out <- normalize_de_output(raw, method = method, seed = seed,
                             n_cells_group1 = n1_cells,
                             n_cells_group2 = n2_cells)
  out <- out[!is.na(out$p_val_adj), , drop = FALSE]
  out <- out[order(out$p_val_adj, -abs(out$avg_log2FC)), ]

  attr(out, "de_audit") <- build_de_audit(
    seed = seed, test = method, group_col = group_col,
    ident1 = ident1, ident2 = ident2,
    n_cells_total = ncol(obj),
    n_group1 = n1_cells, n_group2 = n2_cells,
    downsampled = FALSE, max_cells = Inf, latent_vars = NULL,
    n_genes_tested = nrow(out)
  )
  out
}

.empty_de_df_with_audit <- function(method, seed, n1, n2,
                                    group_col, ident1, ident2, n_total) {
  out <- .empty_de_df()
  attr(out, "de_audit") <- build_de_audit(
    seed = seed, test = method, group_col = group_col,
    ident1 = ident1, ident2 = ident2,
    n_cells_total = n_total, n_group1 = n1, n_group2 = n2,
    downsampled = FALSE, max_cells = Inf, latent_vars = NULL,
    n_genes_tested = 0L
  )
  out
}

.pb_deseq2 <- function(bulk, col_data, ident1, ident2_label) {
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(bulk), colData = col_data, design = ~group)
  dds <- DESeq2::DESeq(dds, quiet = TRUE)
  res <- as.data.frame(DESeq2::results(
    dds, contrast = c("group", ident1, ident2_label)))
  res$gene <- rownames(res)
  names(res)[names(res) == "log2FoldChange"] <- "avg_log2FC"
  names(res)[names(res) == "padj"]           <- "p_val_adj"
  names(res)[names(res) == "pvalue"]         <- "p_val"
  res
}

.pb_edger <- function(bulk, col_data, ident1, ident2_label) {
  y <- edgeR::DGEList(counts = round(bulk), group = col_data$group)
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group, data = col_data)
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  qlf <- edgeR::glmQLFTest(fit)
  res <- as.data.frame(edgeR::topTags(qlf, n = Inf, sort.by = "none"))
  res$gene <- rownames(res)
  names(res)[names(res) == "logFC"]  <- "avg_log2FC"
  names(res)[names(res) == "FDR"]    <- "p_val_adj"
  names(res)[names(res) == "PValue"] <- "p_val"
  res
}

.pb_voom <- function(bulk, col_data, ident1, ident2_label) {
  y <- edgeR::DGEList(counts = round(bulk), group = col_data$group)
  y <- edgeR::calcNormFactors(y)
  design <- stats::model.matrix(~group, data = col_data)
  v <- limma::voom(y, design = design, plot = FALSE)
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  coef_name <- paste0("group", ident1)
  res <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$gene <- rownames(res)
  names(res)[names(res) == "logFC"]    <- "avg_log2FC"
  names(res)[names(res) == "P.Value"]  <- "p_val"
  names(res)[names(res) == "adj.P.Val"] <- "p_val_adj"
  res
}
```

- [ ] **Step 4: Run; confirm green**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-engines.R', reporter='summary')"`
Expected: 9 PASS total (5 SC + 4 PB), or PASS+SKIP for any missing optional package. **No failures.**

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 6 done — pseudobulk engine 3 methods tested`

---

## Task 7: Top-level dispatcher `run_de_analysis()` (TDD)

**Files:**
- Modify: `modules/de/de_methods.R` (append dispatcher)
- Modify: `tests/testthat/test-de-methods.R` (append routing tests)

This is the single entry point used by the server. It runs `preflight_de()`, then routes to the right engine. It is the only function `de_server.R` calls.

- [ ] **Step 1: Append failing tests**

```r
# Append to tests/testthat/test-de-methods.R
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
```

- [ ] **Step 2: Run; confirm failures**

Run: `Rscript -e "testthat::test_file('tests/testthat/test-de-methods.R', reporter='summary')"`
Expected: 4 new failures.

- [ ] **Step 3: Append the implementation**

```r
# Append to modules/de/de_methods.R

#' Top-level DE dispatcher used by the Shiny server.
#' Returns: list(data = data.frame|NULL, error = character|NULL).
run_de_analysis <- function(obj,
                            group_col, ident1, ident2 = NULL,
                            test_method = "wilcox",
                            analysis_type = c("singlecell", "pseudobulk"),
                            sample_col = NULL,
                            subset_col = NULL, subset_val = NULL,
                            latent_vars = NULL,
                            seed = 42L,
                            max_cells_per_ident = Inf,
                            min_pct_prefilter = 0.01,
                            min_pct_display   = 0.10,
                            logfc_threshold   = 0,
                            only_pos = FALSE) {
  analysis_type <- match.arg(analysis_type)

  tryCatch({
    # 1. Optional cell-type subset
    if (!is.null(subset_col) && nzchar(subset_col) &&
        !is.null(subset_val) && nzchar(subset_val)) {
      keep <- which(as.character(obj@meta.data[[subset_col]]) == subset_val)
      if (length(keep) < .MIN_CELLS_SC)
        stop("Too few cells (", length(keep), ") for '", subset_val,
             "' in '", subset_col, "'. Need at least ", .MIN_CELLS_SC, ".")
      obj <- subset(obj, cells = keep)
    }

    # 2. Pre-flight (will throw on bad inputs)
    preflight_de(obj, group_col = group_col, ident1 = ident1, ident2 = ident2,
                 analysis_type = analysis_type, test_method = test_method,
                 sample_col = sample_col)

    # 3. Route
    data <- if (analysis_type == "pseudobulk") {
      run_pseudobulk_de(obj, group_col = group_col,
                        ident1 = ident1, ident2 = ident2,
                        sample_col = sample_col, method = test_method,
                        seed = seed)
    } else {
      run_single_cell_de(obj, group_col = group_col,
                         ident1 = ident1, ident2 = ident2,
                         test_method = test_method, seed = seed,
                         max_cells_per_ident = max_cells_per_ident,
                         latent_vars = latent_vars,
                         min_pct_prefilter = min_pct_prefilter,
                         min_pct_display   = min_pct_display,
                         logfc_threshold   = logfc_threshold,
                         only_pos = only_pos)
    }
    list(data = data, error = NULL)
  }, error = function(e) list(data = NULL, error = conditionMessage(e)))
}
```

- [ ] **Step 4: Run; confirm full backend test suite green**

Run: `Rscript -e "testthat::test_dir('tests/testthat', reporter='summary')"`
Expected: all PASS (with optional SKIPs for uninstalled packages). Total ≥ 22 passing tests.

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 7 done — dispatcher + full backend suite green`

---

## Task 8: `de_ui.R` — split-screen layout

**Files:**
- Create: `modules/de/de_ui.R`

Split-screen per spec §7: left 340 px config, right flex results. Method choice via pill buttons (radioGroupButtons). Audit Log tab in results. **No AI buttons or Auto-DE wellPanel.**

There is no UI test in this task — the test is functional verification at Task 13. The goal here is to lay down a UI that compiles and namespaces cleanly.

- [ ] **Step 1: Write the UI**

```r
# modules/de/de_ui.R
# UI for Differential Expression. Sourced by de_server.R (which is sourced from app.R).

deUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::tags$style(shiny::HTML(sprintf("
      .de-split-%s        { display: flex; gap: 16px; align-items: stretch; }
      .de-config-%s       { flex: 0 0 340px; max-height: 85vh; overflow-y: auto;
                            padding: 12px; border: 1px solid #e3e8ef;
                            border-radius: 8px; background: #fafbfc; }
      .de-results-%s      { flex: 1 1 auto; min-width: 0; }
      .de-section-title   { font-weight: 600; font-size: 13px;
                            text-transform: uppercase; letter-spacing: 0.5px;
                            color: #0D1526; margin-top: 12px; margin-bottom: 6px; }
      .de-help            { font-size: 10px; color: #888; margin-top: 2px; }
    ", id, id, id))),

    shiny::div(class = paste0("de-split-", id),

      # ============ CONFIG (left) ============
      shiny::div(class = paste0("de-config-", id),
        shiny::h4("Differential Expression",
                  style = "margin-top:0; color:#2563EB;"),

        shiny::div(class = "de-section-title", "Comparison"),
        shiny::selectInput(ns("groupby"), "Group by", choices = NULL),
        shiny::selectInput(ns("ident1"),  "Group 1 (test)",     choices = NULL),
        shiny::selectInput(ns("ident2"),  "Group 2 (reference)",
                           choices = c("rest", "all")),

        shiny::div(class = "de-section-title", "Method"),
        shinyWidgets::radioGroupButtons(
          ns("analysis_type"),
          label = NULL,
          choices = c("Single-cell" = "singlecell",
                      "Pseudobulk"  = "pseudobulk"),
          selected = "singlecell", justified = TRUE, size = "sm"
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'singlecell'", ns("analysis_type")),
          shinyWidgets::radioGroupButtons(
            ns("sc_method"),
            label = NULL,
            choices = c("Wilcoxon" = "wilcox", "MAST" = "MAST",
                        "t-test" = "t", "Bimod" = "bimod",
                        "LR" = "LR", "Negbinom" = "negbinom"),
            selected = "wilcox", justified = FALSE, size = "sm",
            individual = TRUE
          )
        ),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] == 'pseudobulk'", ns("analysis_type")),
          shinyWidgets::radioGroupButtons(
            ns("pb_method"),
            label = NULL,
            choices = c("DESeq2" = "DESeq2", "edgeR" = "edgeR",
                        "limma-voom" = "limma-voom"),
            selected = "DESeq2", justified = TRUE, size = "sm"
          ),
          shiny::selectInput(ns("sample_col"),
                             "Sample / replicate column",
                             choices = NULL),
          shiny::div(class = "de-help",
                     "Pseudobulk requires ≥2 replicates per group.")
        ),

        shiny::div(class = "de-section-title", "Subset (optional)"),
        shiny::selectInput(ns("subset_col"),
                           "Cell-type column",
                           choices = c("(none)" = "")),
        shiny::conditionalPanel(
          condition = sprintf("input['%s'] != ''", ns("subset_col")),
          shiny::selectInput(ns("subset_val"), "Cell type", choices = NULL)
        ),

        shiny::div(class = "de-section-title", "Covariates (single-cell)"),
        shiny::selectInput(ns("latent_vars"),
                           "Latent variables", choices = NULL, multiple = TRUE),
        shiny::div(class = "de-help",
                   "MAST: CDR is added automatically."),

        shiny::div(class = "de-section-title", "Display thresholds"),
        shiny::numericInput(ns("min_pct_display"),
                            "Min % expressing", value = 0.10,
                            min = 0, max = 1, step = 0.01),
        shiny::numericInput(ns("logfc_threshold"),
                            "|log2FC| display", value = 0.10,
                            min = 0, max = 5, step = 0.05),
        shiny::numericInput(ns("pval_cutoff"),
                            "Adj p cutoff (highlight)", value = 0.05,
                            min = 1e-6, max = 0.5, step = 0.01),
        shiny::checkboxInput(ns("only_pos"),
                             "Only positive markers", value = FALSE),

        shiny::div(class = "de-section-title", "Reproducibility"),
        shiny::numericInput(ns("seed"),
                            "Random seed", value = 42, min = 1,
                            max = 99999, step = 1),
        shiny::numericInput(ns("max_cells_per_ident"),
                            "Max cells per group",
                            value = Inf, min = 100, step = 1000),
        shiny::div(class = "de-help",
                   "Inf = use all cells (recommended)."),

        shiny::br(),
        shiny::actionButton(ns("run"), "Run DE",
          class = "btn-primary btn-lg",
          style = "width:100%; background:#2563EB; border-color:#2563EB;",
          icon = shiny::icon("play"))
      ),

      # ============ RESULTS (right) ============
      shiny::div(class = paste0("de-results-", id),
        shiny::uiOutput(ns("results_ui"))
      )
    )
  )
}
```

- [ ] **Step 2: Sanity-load the file**

Run: `Rscript -e "source('modules/de/de_ui.R'); cat('OK: deUI =', exists('deUI'), '\n')"`
Expected: `OK: deUI = TRUE`

- [ ] **Step 3: Checkpoint**

Append: `2026-05-09 task 8 done — de_ui.R compiles`

---

## Task 9: `de_server.R` — reactive wiring

**Files:**
- Create: `modules/de/de_server.R`

Server sources `de_core.R`, `de_methods.R`, and `de_ui.R` at the top so that `app.R` only needs to source one file. Implements the reactive graph from spec §4.5: dropdown population, run button, six output renderers (volcano, MA, table, heatmap, gene plots, audit), and writes `values$last_de_results` for downstream consumers.

- [ ] **Step 1: Write the server**

```r
# modules/de/de_server.R
# Sourced by app.R. Sources core_robustness.R, de_core.R, de_methods.R, de_ui.R.

local({
  this_dir <- tryCatch({
    f <- sys.frame(1)$ofile
    if (!is.null(f) && nzchar(f)) normalizePath(dirname(f), mustWork = TRUE)
    else file.path(getwd(), "modules", "de")
  }, error = function(e) file.path(getwd(), "modules", "de"))

  src <- function(rel, base = this_dir) {
    p <- normalizePath(file.path(base, rel), mustWork = FALSE)
    if (file.exists(p)) source(p, local = FALSE)
  }
  # core_robustness.R provides prepare_for_de + safe_get_assay_data
  src("../core_robustness.R", base = dirname(this_dir))
  src("de_core.R")
  src("de_methods.R")
  src("de_ui.R")
})

deServer <- function(id, values) {
  shiny::moduleServer(id, function(input, output, session) {

    de_results       <- shiny::reactiveVal(NULL)
    de_audit         <- shiny::reactiveVal(NULL)
    show_results     <- shiny::reactiveVal(FALSE)

    # ---- Populate dropdowns when active_dataset changes ----
    shiny::observe({
      shiny::req(values$active_dataset)
      meta <- values$active_dataset@meta.data
      cat_cols <- colnames(meta)[vapply(meta,
                    function(x) is.factor(x) || is.character(x), logical(1))]
      num_cols <- colnames(meta)[vapply(meta, is.numeric, logical(1))]

      shiny::updateSelectInput(session, "groupby",   choices = sort(cat_cols))
      shiny::updateSelectInput(session, "subset_col",
                               choices = c("(none)" = "", sort(cat_cols)))
      shiny::updateSelectInput(session, "sample_col",
                               choices = sort(cat_cols),
                               selected = detect_sample_col(values$active_dataset))
      shiny::updateSelectInput(session, "latent_vars", choices = num_cols)
    })

    shiny::observeEvent(input$groupby, {
      shiny::req(values$active_dataset, input$groupby)
      meta <- values$active_dataset@meta.data
      if (input$groupby %in% colnames(meta)) {
        groups <- sort(unique(as.character(meta[[input$groupby]])))
        groups <- groups[!is.na(groups)]
        shiny::updateSelectInput(session, "ident1", choices = groups)
        shiny::updateSelectInput(session, "ident2",
          choices = c("rest", "all", groups))
      }
    })

    shiny::observeEvent(input$subset_col, {
      shiny::req(values$active_dataset)
      col <- input$subset_col
      if (!is.null(col) && nzchar(col) &&
          col %in% colnames(values$active_dataset@meta.data)) {
        vals <- sort(unique(as.character(
          values$active_dataset@meta.data[[col]])))
        vals <- vals[!is.na(vals)]
        shiny::updateSelectInput(session, "subset_val", choices = vals)
      } else {
        shiny::updateSelectInput(session, "subset_val", choices = character(0))
      }
    })

    # ---- Run DE ----
    shiny::observeEvent(input$run, {
      shiny::req(values$active_dataset, input$groupby, input$ident1)

      method <- if (input$analysis_type == "pseudobulk") input$pb_method
                else input$sc_method
      max_cells <- if (is.null(input$max_cells_per_ident) ||
                       is.na(input$max_cells_per_ident) ||
                       !is.finite(input$max_cells_per_ident)) Inf
                   else input$max_cells_per_ident

      shiny::withProgress(
        message = paste("Running", method, "..."), value = 0.1, {
        result <- run_de_analysis(
          obj = values$active_dataset,
          group_col = input$groupby,
          ident1 = input$ident1,
          ident2 = input$ident2,
          test_method = method,
          analysis_type = input$analysis_type,
          sample_col = if (input$analysis_type == "pseudobulk")
                         input$sample_col else NULL,
          subset_col = if (nzchar(input$subset_col)) input$subset_col else NULL,
          subset_val = if (!is.null(input$subset_val) &&
                           nzchar(input$subset_val)) input$subset_val else NULL,
          latent_vars = if (length(input$latent_vars) > 0)
                          input$latent_vars else NULL,
          seed = as.integer(input$seed %||% 42L),
          max_cells_per_ident = max_cells,
          min_pct_display = input$min_pct_display %||% 0.10,
          logfc_threshold = input$logfc_threshold %||% 0,
          only_pos = isTRUE(input$only_pos)
        )
        shiny::incProgress(0.9)

        if (!is.null(result$data)) {
          de_results(result$data)
          de_audit(attr(result$data, "de_audit"))
          values$last_de_results <- result$data
          show_results(TRUE)
          n_sig <- sum(result$data$p_val_adj < (input$pval_cutoff %||% 0.05),
                       na.rm = TRUE)
          shiny::showNotification(
            sprintf("DE complete → %d significant genes (adj p < %.3f).",
                    n_sig, input$pval_cutoff %||% 0.05),
            type = "message", duration = 5)
        } else {
          show_results(FALSE)
          shiny::showNotification(
            paste("DE failed:", result$error),
            type = "error", duration = 10)
        }
      })
    })

    # ---- Results UI scaffold ----
    output$results_ui <- shiny::renderUI({
      shiny::req(isTRUE(show_results()))
      ns <- session$ns
      shiny::tabsetPanel(
        shiny::tabPanel("Volcano",
          shiny::br(), plotly::plotlyOutput(ns("volcano"), height = "560px")),
        shiny::tabPanel("MA Plot",
          shiny::br(), plotly::plotlyOutput(ns("ma"), height = "560px")),
        shiny::tabPanel("Results Table",
          shiny::br(), DT::DTOutput(ns("table")),
          shiny::br(),
          shiny::downloadButton(ns("download_csv"),
                                "Download all (CSV)", class = "btn-success")),
        shiny::tabPanel("Heatmap",
          shiny::br(), plotly::plotlyOutput(ns("heatmap"), height = "640px")),
        shiny::tabPanel("Gene plots",
          shiny::br(),
          shiny::textAreaInput(ns("genes_to_plot"),
                               "Genes (comma-separated)",
                               rows = 2,
                               placeholder = "e.g. Gene001, Gene002"),
          shiny::actionButton(ns("plot_genes"), "Plot",
                              class = "btn-primary"),
          shiny::br(), shiny::br(),
          plotly::plotlyOutput(ns("gene_plots"), height = "560px")),
        shiny::tabPanel("Audit Log",
          shiny::br(), DT::DTOutput(ns("audit_table")))
      )
    })

    # ---- Renderers ----
    output$volcano <- plotly::renderPlotly({
      shiny::req(de_results())
      .render_volcano(de_results(), input)
    })
    output$ma <- plotly::renderPlotly({
      shiny::req(de_results())
      .render_ma(de_results(), input)
    })
    output$table <- DT::renderDT({
      shiny::req(de_results())
      df <- de_results()
      keep_cols <- c("gene","avg_log2FC","p_val","p_val_adj","pct.1","pct.2",
                     "AUC","cohens_d","method")
      keep_cols <- intersect(keep_cols, colnames(df))
      DT::datatable(
        df[, keep_cols, drop = FALSE],
        options = list(pageLength = 25, scrollX = TRUE),
        rownames = FALSE) |>
        DT::formatSignif(c("p_val","p_val_adj"), 3) |>
        DT::formatRound(c("avg_log2FC","pct.1","pct.2","AUC","cohens_d"), 3)
    })
    output$heatmap <- plotly::renderPlotly({
      shiny::req(de_results(), values$active_dataset)
      .render_heatmap(de_results(), values$active_dataset, input)
    })

    shiny::observeEvent(input$plot_genes, {
      output$gene_plots <- plotly::renderPlotly({
        shiny::req(values$active_dataset, input$genes_to_plot)
        genes <- trimws(strsplit(input$genes_to_plot, ",")[[1]])
        genes <- genes[nzchar(genes)]
        genes <- intersect(genes, rownames(values$active_dataset))
        shiny::validate(shiny::need(length(genes) > 0,
                                    "No valid genes found in dataset."))
        .render_gene_box(values$active_dataset, genes, input$groupby)
      })
    })

    output$audit_table <- DT::renderDT({
      shiny::req(de_audit())
      a <- de_audit()
      df <- data.frame(
        Field = names(a),
        Value = vapply(a, function(v) {
          if (length(v) == 0) "(none)"
          else if (length(v) > 1) paste(v, collapse = ", ")
          else as.character(v)
        }, character(1)),
        stringsAsFactors = FALSE
      )
      DT::datatable(df, options = list(pageLength = 20, dom = "t"),
                    rownames = FALSE)
    })

    output$download_csv <- shiny::downloadHandler(
      filename = function()
        sprintf("de_%s_%s.csv",
                attr(de_results(), "de_audit")$test %||% "result",
                format(Sys.time(), "%Y%m%d_%H%M%S")),
      content = function(file) utils::write.csv(
        de_results(), file, row.names = FALSE)
    )
  })
}

# ---- Plot helpers (file-local) ----------------------------------------------
.render_volcano <- function(df, input) {
  pcut <- input$pval_cutoff %||% 0.05
  fcut <- input$logfc_threshold %||% 0.10
  df <- df[!is.na(df$avg_log2FC) & !is.na(df$p_val_adj), ]
  df$neglog10p <- -log10(pmax(df$p_val_adj, 1e-300))
  df$status <- ifelse(df$p_val_adj < pcut & df$avg_log2FC >  fcut, "Up",
               ifelse(df$p_val_adj < pcut & df$avg_log2FC < -fcut, "Down",
                                                                   "NS"))
  pal <- c(Up = "#E64B35", Down = "#3C5488", NS = "#9aa3af")
  top <- head(df[order(df$p_val_adj), "gene"],
              input$label_top_n %||% 20)
  df$label <- ifelse(df$gene %in% top, df$gene, "")
  plotly::plot_ly(df, x = ~avg_log2FC, y = ~neglog10p,
                  type = "scatter", mode = "markers+text",
                  color = ~status, colors = pal,
                  text = ~label, textposition = "top center",
                  hoverinfo = "text",
                  hovertext = ~paste0(gene, "<br>log2FC=", round(avg_log2FC, 2),
                                      "<br>adj.p=", signif(p_val_adj, 3))) |>
    plotly::layout(
      title = "Volcano",
      xaxis = list(title = "avg log2 fold change"),
      yaxis = list(title = "-log10(adj p)"))
}

.render_ma <- function(df, input) {
  df <- df[!is.na(df$avg_log2FC) & !is.na(df$pct.1) & !is.na(df$pct.2), ]
  df$mean_expr <- (df$pct.1 + df$pct.2) / 2
  pcut <- input$pval_cutoff %||% 0.05
  df$status <- ifelse(df$p_val_adj < pcut, "sig", "ns")
  plotly::plot_ly(df, x = ~mean_expr, y = ~avg_log2FC,
                  type = "scatter", mode = "markers",
                  color = ~status,
                  colors = c(sig = "#E64B35", ns = "#9aa3af"),
                  text = ~gene,
                  hoverinfo = "text") |>
    plotly::layout(
      title = "MA plot",
      xaxis = list(title = "mean detection rate"),
      yaxis = list(title = "avg log2 fold change"))
}

.render_heatmap <- function(df, obj, input) {
  top_n <- as.integer(input$label_top_n %||% 20)
  pcut  <- input$pval_cutoff %||% 0.05
  top <- head(df[df$p_val_adj < pcut, "gene"], top_n)
  if (length(top) < 2) {
    return(plotly::plotly_empty() |>
           plotly::layout(title = "Not enough significant genes for heatmap"))
  }
  expr <- safe_get_assay_data(obj, "data")[top, , drop = FALSE]
  group <- as.character(obj@meta.data[[input$groupby]])
  ord <- order(group)
  expr <- expr[, ord, drop = FALSE]; group <- group[ord]
  plotly::plot_ly(z = as.matrix(expr), x = colnames(expr), y = rownames(expr),
                  type = "heatmap", colors = viridis::viridis(50)) |>
    plotly::layout(title = sprintf("Top %d significant genes", length(top)))
}

.render_gene_box <- function(obj, genes, group_col) {
  expr <- safe_get_assay_data(obj, "data")[genes, , drop = FALSE]
  long <- do.call(rbind, lapply(genes, function(g) {
    data.frame(gene = g,
               expr = as.numeric(expr[g, ]),
               group = as.character(obj@meta.data[[group_col]]))
  }))
  plotly::plot_ly(long, x = ~group, y = ~expr,
                  color = ~gene, type = "box") |>
    plotly::layout(boxmode = "group",
                   yaxis = list(title = "log-normalised expression"),
                   xaxis = list(title = group_col))
}

# Null-coalescing operator (used heavily above)
`%||%` <- function(a, b) if (!is.null(a)) a else b
```

- [ ] **Step 2: Sanity-load**

Run: `Rscript -e "source('modules/de/de_server.R'); cat('OK: deUI=', exists('deUI'), 'deServer=', exists('deServer'), '\n')"`
Expected: `OK: deUI= TRUE deServer= TRUE`

- [ ] **Step 3: Checkpoint**

Append: `2026-05-09 task 9 done — de_server.R compiles, sources work`

---

## Task 10: Wire `app.R` — source new module + replace tab body

**Files:**
- Modify: `app.R` (3 specific edits)

- [ ] **Step 1: Edit 1 — register the new module loader**

In `app.R`, find the `modules_loaded <- list(` block (line ~73) and the entry `downstream = safe_source("modules/downstream/downstream_main.R"),`. **Add a sibling line right after it:**

```r
  de = safe_source("modules/de/de_server.R"),
```

The block now contains both `downstream` (still loads CCC, TF, geneset, composition, ai, core; no longer loads downstream_de.R after Task 12) and `de` (loads new DE module).

- [ ] **Step 2: Edit 2 — replace the Downstream tab body**

In `app.R`, find the `tabItem( tabName = "downstream", ...)` block starting at line ~700. Replace the entire `box(...)` content (currently containing the `tabsetPanel` with 8 sub-tabs) with:

```r
        fluidRow(
          box(
            title = div(icon("microscope"), "Differential Expression"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            if (modules_loaded$de) deUI("de")
            else p("DE module not available. Check R console for details.")
          )
        )
```

The closing `tabItem` braces remain.

- [ ] **Step 3: Edit 3 — slim the lazy server init for downstream**

Find the `.downstream_init` observer in `app.R` (line ~1782). Replace its body:

```r
  .downstream_init <- shiny::reactiveVal(FALSE)
  shiny::observe({
    shiny::req(input$tabs == "downstream", !.downstream_init())
    if (modules_loaded$de) {
      deServer("de", values)
    }
    .downstream_init(TRUE)
  })
```

The other server functions (`pathwayServer`, `cellCommServer`, `tfActivityServer`, `grnServer`, `geneSetScoringServer`, `compositionServer`, `aiAnalysisServer`) are **removed from this observer**. They will be re-mounted under the `pathway_network` module in Phase 2. Their source files stay on disk untouched.

- [ ] **Step 4: Smoke-load app.R syntactically**

Run: `Rscript -e "parse('app.R'); cat('app.R parses OK\n')"`
Expected: `app.R parses OK`. Any syntax error must be fixed before continuing.

- [ ] **Step 5: Checkpoint**

Append: `2026-05-09 task 10 done — app.R wired to new DE module`

---

## Task 11: Update `downstream_main.R` to drop old DE source

**Files:**
- Modify: [modules/downstream/downstream_main.R](modules/downstream/downstream_main.R)

- [ ] **Step 1: Edit the `sub_modules` vector**

Replace:

```r
  sub_modules <- c(
    "downstream_core.R",
    "downstream_de.R",
    "downstream_enrichment.R",
    "downstream_ccc.R",
    "downstream_tf.R",
    "downstream_geneset.R",
    "downstream_composition.R",
    "downstream_ai.R"
  )
```

With:

```r
  sub_modules <- c(
    "downstream_core.R",
    "downstream_enrichment.R",
    "downstream_ccc.R",
    "downstream_tf.R",
    "downstream_geneset.R",
    "downstream_composition.R",
    "downstream_ai.R"
  )
```

(Only `downstream_de.R` removed; the rest stay because their source files still exist and may be re-mounted in Phase 2/3. `downstream_enrichment.R` and `downstream_ai.R` will be removed in later phases.)

- [ ] **Step 2: Sanity-source**

Run: `Rscript -e "source('modules/downstream/downstream_main.R'); cat('downstream_main loads OK\n')"`
Expected: `downstream_main loads OK` (warnings about Shiny context are fine; errors are not).

- [ ] **Step 3: Checkpoint**

Append: `2026-05-09 task 11 done — downstream_main.R cleaned up`

---

## Task 12: Delete old `downstream_de.R`

**Files:**
- Delete: [modules/downstream/downstream_de.R](modules/downstream/downstream_de.R)

- [ ] **Step 1: Confirm nothing in the codebase still references it**

Run: `grep -rn "downstream_de.R\|run_pseudobulk_de_api\|run_pseudobulk_de_edger\|create_enhanced_volcano_plot\|create_enhanced_ma_plot\|create_de_summary_table\|create_de_distribution_plot\|create_de_gene_expression_plots" modules/ app.R 2>/dev/null`
Expected: **no output** (everything has been removed in earlier tasks).

If any line is printed, the corresponding file/code must be cleaned up before deletion. Common offenders: lingering source statements, helper imports.

- [ ] **Step 2: Delete the file**

Run: `rm modules/downstream/downstream_de.R`

- [ ] **Step 3: Re-run the test suite — must still pass**

Run: `Rscript -e "testthat::test_dir('tests/testthat', reporter='summary')"`
Expected: same green result as Task 7 step 4.

- [ ] **Step 4: Checkpoint**

Append: `2026-05-09 task 12 done — old downstream_de.R deleted, suite green`

---

## Task 13: Functional verification

**Files:** none modified.

This task is the manual UI/UX acceptance check that closes Phase 1. The executor runs the full app, exercises every method, and confirms the spec §10 success criteria for Phase 1.

- [ ] **Step 1: Boot the app**

Run: `Rscript -e "shiny::runApp('.', launch.browser = FALSE, port = 4781)"` (foreground; open `http://127.0.0.1:4781` in a browser).
Expected: app loads, "Downstream Analysis" sidebar item appears.

- [ ] **Step 2: Navigate to Downstream Analysis tab**

Confirm:
- Page shows ONE area: "Differential Expression" (no Pathway / CCC / TF / GRN / Geneset / Composition / AI sub-tabs).
- Layout is split-screen (left config 340 px, right results placeholder).
- No "AI Recommend Test", "AI Interpret DE", or "AI Auto-Run DE" panels are visible anywhere.

- [ ] **Step 3: Run all 6 single-cell methods**

Load any project / dataset that has at least one categorical metadata column with two well-populated levels (≥ 50 cells each). For each of `wilcox`, `MAST`, `t`, `bimod`, `LR`, `negbinom`:
1. Pick `analysis_type = singlecell`, choose the method pill, set ident1 / ident2.
2. Click Run DE.
3. Verify: results panel populates, Volcano + MA + Table + Audit Log tabs all render without error.
4. Open the Audit Log tab → confirm `seed`, `test`, `n_group1`, `n_group2`, `seurat_version`, `r_version` are all populated (no `NA`).

- [ ] **Step 4: Run all 3 pseudobulk methods**

Use the same dataset. Switch `analysis_type = pseudobulk`, pick a sample column with ≥ 2 replicates per group. For each of `DESeq2`, `edgeR`, `limma-voom`:
1. Pick the pill, click Run DE.
2. Verify results render and the Audit Log shows the right test name.
3. Check that the Results Table `method` column matches the chosen method.

- [ ] **Step 5: Reproducibility check**

1. Run `wilcox` with seed 42, copy the top-10 gene list (table sorted by adj p).
2. Run `wilcox` again with seed 42 — top-10 must be identical (same genes, same p-values to 3 sig figs).
3. Run with seed 7 — at least the order or values should differ if any randomized step depended on the seed.

- [ ] **Step 6: Pre-flight error surfacing**

1. Try pseudobulk with a sample column that has only 1 sample per group → error notification appears: "Pseudobulk requires at least 2 replicates per group..." (no crash).
2. Try a `group_col` with one level that has < 10 cells → error notification: "fewer than 10 cells".

- [ ] **Step 7: Cross-module hand-off**

1. After a successful DE run, open the Export tab.
2. Confirm Export still recognises DE results (export.R reads `values$last_de_results`; the data frame must be present).

- [ ] **Step 8: Spec §10 success-criteria checklist**

Tick off:
- [ ] §10.1 — all 9 methods ran on a 25 k-cell dataset (or as close as the test dataset allows; document the actual cell count).
- [ ] §10.2 — bit-for-bit identical results across two same-seed runs (Step 5).
- [ ] §10.5 — no code duplication (`run_de_analysis()` exists only in `modules/de/de_methods.R`; verify `grep -rn "^run_de_analysis " modules/`).
- [ ] §10.7 (DE part) — `modules/downstream/downstream_de.R` no longer exists.
- [ ] §10.9 — dead-code checklist passes (run the orphan-detection greps from spec §9.2 and confirm no DE-related orphans).

- [ ] **Step 9: Final checkpoint**

Append: `2026-05-09 task 13 done — Phase 1 acceptance verified, ready for Phase 2`

---

## Self-Review Notes

**Spec coverage map (Phase 1 only):**

| Spec section | Covered by |
|---|---|
| §2.2 sidebar — Downstream = DE only | Task 10 step 2 |
| §2.3 new files `de/{core,methods,ui,server}.R` | Tasks 2, 3-7, 8, 9 |
| §2.3 delete old `downstream_de.R` | Task 12 |
| §2.4 update `downstream_main.R` | Task 11 |
| §2.4 `app.R` sources `de_server.R` directly | Task 10 step 1 |
| §3 no hard-coding — species, columns, thresholds | Task 2 (helpers) + Task 9 (UI defaults only) |
| §4.1 design principles — single dispatcher, full dataset, deterministic, Seurat v5 safe, standardised output | Tasks 4, 5, 6, 7 |
| §4.2 9 methods (6 SC + 3 PB) | Task 5 (SC) + Task 6 (PB inc. limma-voom) |
| §4.3 pre-flight validation | Task 3 |
| §4.4 audit log via `attr(result, "de_audit")` | Tasks 4, 5, 6 |
| §4.5 reactive wiring | Task 9 |
| §7 split-screen UI, no AI buttons, audit log tab | Task 8 |
| §9.1 + §9.3 dead-code / dynamic-only | Task 13 step 8 |
| §10 success criteria (DE-relevant items) | Task 13 |

**No gaps detected.** Every Phase 1 spec requirement maps to at least one task.

**Type consistency:** dispatcher signature `run_de_analysis(obj, group_col, ident1, ident2, test_method, analysis_type, sample_col, subset_col, subset_val, latent_vars, seed, max_cells_per_ident, min_pct_prefilter, min_pct_display, logfc_threshold, only_pos)` is the same in Tasks 7 (definition) and 9 (call site). Engine signatures match between Tasks 5/6 (definitions) and 7 (calls). Output schema column list `gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, AUC, cohens_d, method, seed, n_cells_group1, n_cells_group2` is identical in Tasks 4 (definition), 5/6 (engines), 7 (tests), 9 (renderers).

**Placeholders:** none — every step has the actual code, command, or expected output.
