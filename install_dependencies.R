# =============================================================================
# scRNA7 - Dependency Installation Script
# Run this script to install all required packages on a new machine
# Usage: Rscript install_dependencies.R
# =============================================================================

cat("=== scRNA-seq Platform Dependency Installer (with Auth) ===\n\n")

# Helper: install CRAN package if not already installed
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# Helper: install Bioconductor package if not already installed
install_bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s (Bioconductor)...\n", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
  } else {
    cat(sprintf("  %s already installed\n", pkg))
  }
}

# ---- Step 1: CRAN Packages ----
cat("\n--- Step 1: Installing CRAN packages ---\n")

cran_packages <- c(
  # Shiny framework
  "shiny", "shinydashboard", "shinyWidgets", "shinyjs",
  "shinycssloaders", "shinyFiles", "waiter",
  # Interactive output
  "DT", "plotly", "htmlwidgets", "visNetwork",
  # Data manipulation
  "dplyr", "tidyr", "Matrix", "data.table", "readxl", "rlang",
  # Visualization
  "ggplot2", "viridis", "RColorBrewer", "patchwork", "circlize", "corrplot", "scales", "ggridges",
  # Dimensionality reduction
  "uwot", "Rtsne", "irlba",
  # Clustering & networks
  "igraph", "cluster", "factoextra",
  # WGCNA
  "WGCNA",
  # Async / parallelism
  "future", "promises", "later",
  # Utilities
  "R6", "base64enc", "digest", "htmltools", "jsonlite", "markdown",
  "openxlsx", "rmarkdown", "knitr", "zip",
  # ── NEW: Authentication system (DBI + SQLite + crypto) ─────────────────────
  "DBI",        # database interface abstraction (required by auth)
  "RSQLite",    # SQLite backend for user/session/project database
  "openssl",    # cryptographic random token generation
  "bcrypt",     # secure password hashing (argon2/bcrypt; falls back to digest)
  # ────────────────────────────────────────────────────────────────────────────
  # API calls
  "httr2",
  # Graphics export
  "webshot2",
  # Batch correction
  "harmony",
  # Cluster validation
  "fpc",
  # Cluster stability visualization
  "clustree", "ggraph",
  # Gene symbol correction (used by scType)
  "HGNChelper",
  # Leiden algorithm (community detection)
  "leiden",
  # KNN propagation (DoubletFinder subsampling)
  "FNN",
  # System / platform
  "unix",
  # Misc
  "mgcv", "reticulate", "remotes", "readxl"
)

for (pkg in cran_packages) {
  tryCatch(install_if_missing(pkg), error = function(e) {
    cat(sprintf("  WARNING: Failed to install %s: %s\n", pkg, e$message))
  })
}

# ---- Step 2: BiocManager ----
cat("\n--- Step 2: Setting up BiocManager ---\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# ---- Step 3: Core Bioconductor Packages ----
cat("\n--- Step 3: Installing Bioconductor packages ---\n")

bioc_core <- c(
  # Seurat ecosystem
  "Seurat", "SeuratObject",
  # Single-cell data structures
  "SingleCellExperiment", "SummarizedExperiment", "BiocGenerics", "S4Vectors",
  # Normalization & preprocessing
  "sctransform", "glmGamPoi",
  # Batch correction
  "batchelor", "sva", "BiocSingular",
  # Cell type annotation (clustifyr = reference-correlation; celldex = references)
  "clustifyr", "celldex",
  # OrgDbs — required by modules/pathway_network/ pn_orgdb_for_species()
  # for the 4 supported species (Latin-or-lowercase normalised in pn_core.R).
  "org.Hs.eg.db",   # human
  "org.Mm.eg.db",   # mouse
  "org.Rn.eg.db",   # rat        (NEW for Phase 2 — required by pn_core)
  "org.Dr.eg.db",   # zebrafish  (NEW for Phase 2 — required by pn_core)
  "AnnotationDbi",
  # Differential expression — required by modules/de/de_methods.R
  # All 9 supported methods need these (limma drives both DESeq2/edgeR contrast
  # construction and the new limma-voom pseudobulk path).
  "limma", "DESeq2", "edgeR", "MAST",
  # Pathway enrichment — required by modules/pathway_network/pn_enrichment.R
  # for both the live-API path (clusterProfiler::enrich*/gse* + ReactomePA) and
  # the MSigDB path (msigdbr + fgsea).
  "clusterProfiler", "enrichplot", "fgsea", "msigdbr", "ReactomePA",
  # STRING PPI enrichment — required for the "STRING PPI" database in pn_ui.R
  "STRINGdb",       # NEW for Phase 2
  # Trajectory analysis
  "slingshot", "destiny",
  # Visualization
  "ComplexHeatmap"
)

for (pkg in bioc_core) {
  tryCatch(install_bioc_if_missing(pkg), error = function(e) {
    cat(sprintf("  WARNING: Failed to install %s: %s\n", pkg, e$message))
  })
}

# ---- Step 4: Optional Bioconductor Packages ----
cat("\n--- Step 4: Installing optional packages (failures are OK) ---\n")

bioc_optional <- c(
  "SeuratWrappers",
  "monocle3",
  "scDblFinder",
  "SoupX",
  "celda",           # decontX ambient RNA removal
  "zellkonverter",
  "hdf5r",
  "BPCells",         # on-disk (memory-mapped) matrices — low-RAM large datasets
  "AUCell",          # SCENIC / gene set scoring
  "UCell",           # gene set scoring (AddModuleScore_UCell)
  "speckle",         # cell composition / propeller testing
  "decoupleR",       # TF activity inference (downstream_tf.R)
  "OmnipathR",       # CollecTRI regulons + L-R pairs (downstream_tf.R / ccc.R)
  "liana",           # CCC consensus method (downstream_ccc.R, optional)
  "CellChat",        # CCC inference (downstream_ccc.R, optional)
  "escape",          # gene set enrichment in single-cell
  "scrapper",        # required for SingleR >= 2.0
  "loomR"            # Loom file format support
)

# Optional CRAN packages (loaded on demand by pathway_network sub-modules).
cran_optional <- c(
  "ggupset",         # upset-style overlap plots for enrichment results
  "homologene",      # offline ortholog mapping (downstream_core.R fallback)
  "entropy",         # clustifyr dependency
  "shinythemes"      # required to install scHCL / scMCA (Han-lab atlases)
)
for (pkg in cran_optional) {
  tryCatch(install_if_missing(pkg), error = function(e) {
    cat(sprintf("  SKIPPED (optional): %s - %s\n", pkg, e$message))
  })
}

for (pkg in bioc_optional) {
  tryCatch(install_bioc_if_missing(pkg), error = function(e) {
    cat(sprintf("  SKIPPED (optional): %s - %s\n", pkg, e$message))
  })
}

# ---- Step 5: GitHub Packages ----
cat("\n--- Step 5: Installing GitHub packages ---\n")

github_packages <- list(
  "DoubletFinder" = "chris-mcginnis-ucsf/DoubletFinder",
  "SeuratDisk"    = "mojaveazure/seurat-disk",
  # Han-lab reference atlases used as clustifyr references (annotation.R):
  # scHCL = Human Cell Landscape, scMCA = Mouse Cell Atlas. Need shinythemes.
  "scHCL"         = "ggjlab/scHCL",
  "scMCA"         = "ggjlab/scMCA"
)

for (pkg_name in names(github_packages)) {
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    cat(sprintf("Installing %s from GitHub...\n", pkg_name))
    tryCatch(
      remotes::install_github(github_packages[[pkg_name]], quiet = TRUE, upgrade = "never"),
      error = function(e) cat(sprintf("  SKIPPED (optional): %s - %s\n", pkg_name, e$message))
    )
  } else {
    cat(sprintf("  %s already installed\n", pkg_name))
  }
}

# ---- Step 6: Verification ----
cat("\n\n=== Verification ===\n")

all_packages <- c(cran_packages, bioc_core, bioc_optional, names(github_packages))
# Deduplicate
all_packages <- unique(all_packages)
installed <- vapply(all_packages, function(p) requireNamespace(p, quietly = TRUE), logical(1))

n_installed <- sum(installed)
n_total <- length(all_packages)
n_missing <- n_total - n_installed

cat(sprintf("Installed: %d / %d packages\n", n_installed, n_total))

if (n_missing > 0) {
  cat("\nMissing packages (optional ones may be OK to skip):\n")
  for (pkg in all_packages[!installed]) {
    is_optional <- pkg %in% c(bioc_optional, names(github_packages))
    cat(sprintf("  - %s%s\n", pkg, if (is_optional) " (optional)" else " ** REQUIRED **"))
  }
}

cat("\n=== Done! Run the app with: shiny::runApp('.') ===\n")
cat("\nFirst launch will create:  data/scrna4_auth.db  (SQLite user database)\n")
cat("Project RDS files saved to: data/user_sessions/<user_id>/<project_id>/analysis.rds\n\n")
