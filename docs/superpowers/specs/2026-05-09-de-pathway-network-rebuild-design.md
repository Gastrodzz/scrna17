# Design Spec: DE & Pathway/Network Analysis Rebuild
**Date:** 2026-05-09  
**Status:** Approved  
**Scope:** Phase 1 (DE) + Phase 2 (Enrichment) + Phase 3 (Network consolidation)

---

## 1. Problem Statement

The existing scRNA-seq platform has three critical failures:

1. **DE is broken** — only MAST runs but silently analyzes <50% of cells on large datasets (25k+ cells). All other methods (Wilcoxon, DESeq2, edgeR, limma-voom, t-test, Bimod, LR, Negbinom) crash or produce wrong results. Results are not reproducible across runs.

2. **Enrichment is mis-wired** — `clusterProfiler`, `msigdbr`, `ReactomePA`, `fgsea` are loaded but not correctly connected to their databases. Gene universes are wrong, gene ID conversion is missing, results are inaccurate and not reproducible.

3. **Duplication across modules** — GRN logic exists in both `downstream_tf.R` and `network_analysis.R`. L-R interactions exist in both `downstream_ccc.R` and `network_analysis.R`. Pathway network visualization depends on `values$last_enrich_results` which is often NULL. AI buttons are scattered everywhere.

---

## 2. Approved Architecture

### 2.1 New Sidebar Structure

| Sidebar Item | Sub-tabs |
|---|---|
| **Downstream Analysis** | Differential Analysis (only) |
| **Pathway & Network Analysis** | Pathway Enrichment, Pathway Network, Gene Regulatory Networks, Co-expression, Cell-Cell Communication, TF Activity, Gene Set Scoring, Composition |

### 2.2 New File Structure

```
modules/
  de/
    de_core.R          # shared helpers: layer joining, seed mgmt, audit log
    de_methods.R       # all 9 method engines (single-cell + pseudobulk)
    de_ui.R            # UI only — split-screen layout
    de_server.R        # server: reactive wiring, output renderers
  pathway_network/
    pn_core.R          # shared helpers: gene ID conversion, background sets
    pn_enrichment.R    # ORA + GSEA engines, all database wiring
    pn_network.R       # GRN (consolidated), co-expression
    pn_ui.R            # UI — all 8 sub-tabs
    pn_server.R        # server: reactive wiring, output renderers
```

### 2.3 Files to DELETE / MODIFY

**Delete entirely:**
- `modules/downstream/downstream_de.R` → replaced by `modules/de/`
- `modules/downstream/downstream_enrichment.R` → replaced by `modules/pathway_network/pn_enrichment.R`
- `modules/network_analysis.R` → all surviving content moved to `modules/pathway_network/`

**Modify in-place (keep file, edit content):**
- `modules/downstream/downstream_tf.R` → delete `build_gene_regulatory_network()` and all GRN-related functions; keep TF activity scoring only
- `modules/downstream/downstream_main.R` → remove `downstream_de.R` and `downstream_enrichment.R` from source list; `modules/de/de_server.R` is sourced by `app.R` directly

**Files that stay untouched (re-routed only):**
- `modules/downstream/downstream_ccc.R` — sourced by `pn_server.R`
- `modules/downstream/downstream_geneset.R` — sourced by `pn_server.R`
- `modules/downstream/downstream_composition.R` — sourced by `pn_server.R`
- `modules/downstream/downstream_ai.R` — left as-is for now

### 2.4 Files to UPDATE

- `modules/downstream/downstream_main.R` → source `modules/de/de_server.R` instead of `downstream_de.R`; remove `downstream_enrichment.R`
- `app.R` → add `pathway_network` module source, rename sidebar entry, remove old `network` module

---

## 3. Global Engineering Constraint — No Hard-Coding

This applies to **every file** in the rebuild. Any value that could vary across datasets, species, machines, or user sessions must be derived at runtime, not written as a literal.

### What must never be hard-coded

| Category | Bad (hard-coded) | Good (derived) |
|---|---|---|
| Column names | `obj@meta.data$seurat_clusters` | `obj@meta.data[[group_col]]` where `group_col` is detected or user-selected |
| Cell type labels | `"T cells"`, `"Endothelial"` | Read from `unique(obj@meta.data[[cell_type_col]])` |
| Gene names | `c("CD3D", "CD8A")` | Read from `rownames(obj)` or user input |
| Species / OrgDb | `"org.Hs.eg.db"` literal | `orgdb_for_species(values$species_global)` helper that maps species → OrgDb |
| Thresholds | `pval < 0.05`, `logfc > 0.25` | Always from UI inputs: `input$pval_cutoff`, `input$logfc_threshold` |
| Memory limits | `max_genes <- 2000` | Derived from available RAM: `floor(memuse::Sys.meminfo()$totalram * 0.3 / 8 / ncol(expr))` or exposed as UI parameter |
| Cell count caps | `if (ncol > 50000)` | Warning threshold exposed as app-level config or UI input; never silent |
| File paths | `"modules/downstream/..."` | Resolved relative to `this_file` or app root, same strategy as `downstream_main.R` |
| Package names | `library("MAST")` inline | Always via `requireNamespace()` check + informative error if missing |
| Permutation counts | `nPermSimple = 1000` | From `input$gsea_nperm` |
| Random seeds | `set.seed(42)` baked in | Always `set.seed(input$de_seed)` or `set.seed(input$enrich_seed)` |
| Assay names | `"RNA"` | `DefaultAssay(obj)` |
| Layer names | `"data"`, `"counts"` | Checked via `Layers(obj[[assay]])` before access |
| Number of top genes shown | `head(..., 20)` | From `input$label_top_n` or `input$max_pathways` |
| KEGG organism code | `"hsa"` | Derived from species: `kegg_code_for_species(values$species_global)` |
| Reactome organism | `"human"` | Derived from `values$species_global` |
| Background gene set | any fixed list | `rownames(values$active_dataset)` always, with custom override |

### Helper functions required in `pn_core.R` and `de_core.R`

```r
# Maps species name → OrgDb package name
orgdb_for_species <- function(species) {
  map <- c(human="org.Hs.eg.db", mouse="org.Mm.eg.db",
           rat="org.Rn.eg.db", zebrafish="org.Dr.eg.db")
  pkg <- map[tolower(species)]
  if (is.na(pkg)) stop("Unsupported species: ", species)
  pkg
}

# Maps species name → KEGG organism code
kegg_code_for_species <- function(species) {
  map <- c(human="hsa", mouse="mmu", rat="rno", zebrafish="dre")
  code <- map[tolower(species)]
  if (is.na(code)) stop("No KEGG code for species: ", species)
  code
}

# Detects the best cell-type / cluster column from Seurat metadata
detect_cell_type_col <- function(obj) {
  candidates <- c("cell_type", "celltype", "cell.type", "annotation",
                  "seurat_clusters", "cluster")
  found <- intersect(candidates, colnames(obj@meta.data))
  if (length(found) > 0) found[1] else colnames(obj@meta.data)[1]
}

# Detects sample / replicate column
detect_sample_col <- function(obj) {
  candidates <- c("orig.ident", "sample", "sample_id", "donor",
                  "patient", "replicate", "batch")
  found <- intersect(candidates, colnames(obj@meta.data))
  if (length(found) > 0) found[1] else NULL
}
```

---

## 4. Phase 1 — DE Backend Rebuild

### 4.1 Core Design Principles (in addition to §3)

- **One entry point**: `run_de_analysis()` in `de_methods.R` is the single dispatcher used by all UI paths. No parallel code paths.
- **Full dataset by default**: `max_cells_per_ident = Inf`. Subsampling is an explicit user choice, never silent.
- **Deterministic**: fixed seed passed to every randomized operation. Seed stored in audit metadata.
- **Seurat v5 safe**: always call `prepare_for_de()` which joins split layers, ensures `data` slot exists.
- **Standardized output**: every method returns the same column schema: `gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2, AUC, cohens_d, method, seed, n_cells_group1, n_cells_group2`.

### 4.2 Supported Methods

**Single-cell (via `FindMarkers`):**

| Method | Key parameter | Notes |
|---|---|---|
| `wilcox` | — | Fast, non-parametric. Default for exploratory. |
| `MAST` | `latent.vars = c("CDR", user_vars)` | CDR auto-computed and added. Never duplicated. |
| `t` | — | Parametric, assumes normality. |
| `bimod` | — | Likelihood ratio for bimodal distributions. |
| `LR` | `latent.vars` | Logistic regression with covariates. |
| `negbinom` | — | Negative binomial GLM. Requires MASS. |

**Pseudobulk (count aggregation per sample):**

| Method | Package | Notes |
|---|---|---|
| `DESeq2` | `DESeq2` | Requires ≥2 replicates per group. Round counts. Wald test. |
| `edgeR` | `edgeR` | Quasi-likelihood F-test. TMM normalisation. |
| `limma-voom` | `limma` | voom precision weights. Best for unbalanced designs. |

### 4.3 Pre-flight Validation (run before any method)

1. Check `sample_col` exists in metadata (pseudobulk only)
2. Check ≥2 replicates per group (pseudobulk only) — hard stop with clear message
3. Check ≥10 cells per group (single-cell) — hard stop
4. Check `data` layer populated — call `NormalizeData()` if counts-only
5. Estimate runtime based on cell count + method; show notification if >1 min
6. If `max_cells_per_ident` is set, downsample with `set.seed(seed)` and flag result

### 4.4 Audit Log

Every result carries an `attr(result, "de_audit")` list:

```r
list(
  seed           = seed,
  test           = test_method,
  group_col      = group_col,
  ident1         = ident1,
  ident2         = ident2,
  n_cells_total  = ncol(obj),
  n_group1       = length(cells1),
  n_group2       = length(cells2),
  downsampled    = is.finite(max_cells_per_ident),
  max_cells      = max_cells_per_ident,
  latent_vars    = latent_vars,
  n_genes_tested = nrow(result),
  timestamp      = Sys.time(),
  seurat_version = packageVersion("Seurat"),
  r_version      = R.version$version.string
)
```

This is rendered in the "Audit Log" result tab.

### 4.5 Reactive Wiring (de_server.R)

```
values$active_dataset
  └─► observe() → populate de_groupby, de_ident1, de_ident2, covariates
  
input$run_de
  └─► run_de_analysis() [in de_methods.R]
        └─► de_results(result)
              └─► show_de_results(TRUE)
                    ├─► output$volcano_plot
                    ├─► output$ma_plot
                    ├─► output$de_table
                    ├─► output$de_heatmap
                    ├─► output$de_gene_plots
                    └─► output$de_audit_table
                    
input$send_to_enrichment
  └─► values$de_results_for_enrichment ← de_results()
  └─► updateTabItems(session, "tabs", "pathway_network")
```

---

## 5. Phase 2 — Pathway & Network Enrichment Backend

### 5.1 Core Design Principles

- **One gene ID conversion step**: `convert_gene_ids()` in `pn_core.R` handles SYMBOL→ENTREZID→ENSEMBL. Called once, cached in the result.
- **Background gene set**: always the full set of detected genes from `rownames(values$active_dataset)`. Custom background is an override.
- **Per-database dispatch**: each database has its own tested wrapper function. No generic catch-all.
- **Species-aware**: OrgDb package selected dynamically: `org.Hs.eg.db` (human), `org.Mm.eg.db` (mouse), `org.Rn.eg.db` (rat), `org.Dr.eg.db` (zebrafish).

### 5.2 Database Wiring

**ORA — Over-Representation Analysis:**

| Database | Package | Function | Gene ID |
|---|---|---|---|
| GO BP/MF/CC | `clusterProfiler` | `enrichGO()` | SYMBOL (keyType="SYMBOL") |
| KEGG | `clusterProfiler` | `enrichKEGG()` | ENTREZID (auto-convert) |
| Reactome | `ReactomePA` | `enrichPathway()` | ENTREZID (auto-convert) |
| WikiPathways | `clusterProfiler` | `enrichWP()` | ENTREZID (auto-convert) |
| MSigDB H/C2/C5/C6 | `msigdbr` + `clusterProfiler` | `enricher()` | SYMBOL |

**GSEA — Gene Set Enrichment Analysis:**

| Database | Package | Function | Ranking metric |
|---|---|---|---|
| GO BP/MF/CC | `clusterProfiler` | `gseGO()` | signed −log10(p)×sign(FC) |
| KEGG | `clusterProfiler` | `gseKEGG()` | signed −log10(p)×sign(FC) |
| Reactome | `ReactomePA` | `gsePathway()` | signed −log10(p)×sign(FC) |
| MSigDB all | `fgsea` | `fgsea()` | signed −log10(p)×sign(FC) |

### 5.3 Gene ID Conversion (`pn_core.R`)

```r
convert_gene_ids <- function(genes, from = "SYMBOL", to = "ENTREZID", species = "human") {
  orgdb <- switch(species,
    human    = "org.Hs.eg.db",
    mouse    = "org.Mm.eg.db",
    rat      = "org.Rn.eg.db",
    zebrafish = "org.Dr.eg.db"
  )
  requireNamespace(orgdb, quietly = TRUE)
  result <- clusterProfiler::bitr(genes, fromType = from, toType = to,
                                   OrgDb = orgdb, drop = TRUE)
  # Return named vector: ENTREZID → SYMBOL for traceability
  result
}
```

### 5.4 Ranked Gene List Construction (GSEA)

```r
build_ranked_list <- function(de_result, metric = "signed_neglog10p") {
  df <- de_result[!is.na(de_result$p_val) & !is.na(de_result$avg_log2FC), ]
  ranks <- switch(metric,
    signed_neglog10p = -log10(pmax(df$p_val, 1e-300)) * sign(df$avg_log2FC),
    log2fc           = df$avg_log2FC,
    t_statistic      = if ("t" %in% colnames(df)) df$t else df$avg_log2FC
  )
  names(ranks) <- df$gene
  sort(ranks, decreasing = TRUE)
}
```

### 5.5 Reproducibility

- `fgsea`: always called with `set.seed(seed)` + `nPermSimple = nperm`
- `clusterProfiler GSEA`: `seed = TRUE` parameter
- All results carry audit metadata: databases used, n_genes_input, n_background, species, seed, timestamp

### 5.6 Reactive Wiring (pn_server.R)

```
values$de_results_for_enrichment   (set by DE "Send" button)
  └─► observe() → auto-populate gene input, show "DE results loaded" badge

values$active_dataset
  └─► observe() → set background gene universe = rownames(obj)

input$run_enrichment
  └─► run_ora_analysis()   [pn_enrichment.R]  → enrichment_results$ora
  └─► run_gsea_analysis()  [pn_enrichment.R]  → enrichment_results$gsea
        └─► enrichment_results() ← combined list
              ├─► output$dotplot
              ├─► output$barplot
              ├─► output$gsea_plot
              ├─► output$enrichment_map   (Jaccard pathway network, rendered by pn_network.R)
              ├─► output$upset_plot
              └─► output$enrich_table
              
input$view_pathway_network
  └─► updateTabsetPanel to "Pathway Network" sub-tab
  └─► values$last_enrich_results ← enrichment_results()   [for pn_network.R]
```

---

## 6. Phase 3 — Network Consolidation

### 5.1 What Moves Where

| Feature | Old location | New location | Action |
|---|---|---|---|
| GRN (WGCNA + correlation) | `network_analysis.R` + `downstream_tf.R` | `pn_network.R` | Keep `network_analysis.R` WGCNA (has pseudobulk aggregation + memory guards). Delete GRN block from `downstream_tf.R`. |
| Co-expression | `network_analysis.R` | `pn_network.R` | Move verbatim |
| L-R Interactions | `network_analysis.R` | — | Delete entirely — `downstream_ccc.R` (CellChat) is the authoritative CCC implementation |
| Pathway Networks | `network_analysis.R` | `pn_server.R` "Pathway Network" sub-tab | Move viz logic; input now comes from `values$last_enrich_results` set by `pn_server.R` in same module |
| GRN logic in TF module | `downstream_tf.R` | — | Delete `build_gene_regulatory_network()` block; keep only TF activity scoring functions |
| CCC | `downstream_ccc.R` | `pn_server.R` "CCC" sub-tab | File stays; `pn_server.R` sources it and calls `cccUI(ns("ccc"))` + `cccServer("ccc", values)` |
| TF Activity | `downstream_tf.R` | `pn_server.R` "TF Activity" sub-tab | File stays (GRN block removed); `pn_server.R` sources it and calls `tfUI(ns("tf"))` + `tfServer("tf", values)` |
| Gene Set Scoring | `downstream_geneset.R` | `pn_server.R` "Gene Set Scoring" sub-tab | File stays; `pn_server.R` sources and calls its UI/server |
| Composition | `downstream_composition.R` | `pn_server.R` "Composition" sub-tab | File stays; `pn_server.R` sources and calls its UI/server |

### 5.2 values$ Reactive Contracts

All cross-module communication goes through the shared `values` reactiveValues object:

```r
# Set by DE module:
values$de_results_for_enrichment   # data.frame with standard DE schema
values$de_audit                    # audit metadata list

# Set by Enrichment module:
values$last_enrich_results         # list(ora=..., gsea=...) for pathway network tab

# Read by all modules:
values$active_dataset              # Seurat object
values$species_global              # "human" | "mouse" | "rat" | "zebrafish"
values$tissue_global               # free-text tissue name
```

---

## 7. UI Design (Approved)

- **Design tokens**: use existing `app.css` variables (`--primary: #2563EB`, `--font-sans: Inter`, dark sidebar `#0D1526`)
- **Layout**: split-screen — scrollable config panel (left, 340px) + scrollable results panel (right, flex-1)
- **Downstream Analysis tab**: contains only "Differential Analysis" sub-tab
- **Pathway & Network Analysis tab**: contains 8 sub-tabs (Pathway Enrichment, Pathway Network, GRN, Co-expression, CCC, TF Activity, Gene Set Scoring, Composition)
- **No AI buttons** in DE or Enrichment (removed; to be re-added in a later sprint)
- **Method selection**: pill buttons (not dropdowns) for method choice
- **Audit Log tab**: always present in DE results, shows full reproducibility metadata

---

## 8. Package Dependencies

### Required (must be installed)
```r
# DE
"Seurat", "MAST", "DESeq2", "edgeR", "limma", "MASS"

# Enrichment
"clusterProfiler", "msigdbr", "fgsea", "ReactomePA",
"org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db", "org.Dr.eg.db",
"enrichplot", "ggplot2", "ggupset"

# Network
"igraph", "visNetwork", "WGCNA", "corrplot"
```

### Checked at runtime (optional, degrade gracefully)
```r
"decoupleR"   # TF activity (falls back to built-in regulon db)
"circlize"    # chord diagrams
"piano"       # alternative GSEA
```

---

## 9. Code Hygiene — Dead Code & Dynamic-Only Policy

This applies to **every file touched during the rebuild**, not just new files.

### 9.1 Dead Code Removal Rules

Every phase must pass this checklist before moving to the next:

| Check | Rule |
|---|---|
| Commented-out code blocks | Delete entirely — git history is the undo. No `# old code` blocks left behind |
| Unused function definitions | Delete any function not called by at least one active code path |
| Orphaned `output$` renderers | Delete any `output$xxx <- render*({})` with no matching `xxxOutput()` in UI |
| Orphaned `observeEvent` / `observe` | Delete any observer referencing a non-existent `input$` or `values$` slot |
| Unused `reactive()` / `reactiveVal()` | Delete if not consumed downstream |
| Unreachable `if/else` branches | Delete branches whose condition can never be TRUE given current UI |
| Duplicate function definitions | If same function name defined in two files, delete the inferior copy |
| Leftover `print()` / `message()` debug statements | Delete all — use `shiny::showNotification()` for user-facing feedback only |
| Unused `library()` / `requireNamespace()` calls | Remove packages no longer needed after refactor |
| Empty `tryCatch` catch blocks | Replace with a real error message or remove entirely |
| Fallback hardcoded gene/pathway lists | Remove all built-in fallback lists (e.g. the minimal L-R pairs in `load_lr_database()`) — if a package isn't installed, fail loudly with install instructions |

### 9.2 Orphan Detection — Opus 4.7 Must Run These After Each Phase

```bash
# All output$ IDs defined in server files
grep -rn "output\$" modules/ | sed 's/.*output\$\([a-zA-Z_]*\).*/\1/' | sort -u

# All *Output() calls in UI files (what's actually rendered)
grep -rn "Output(" modules/ | grep -v "^Binary" | sort -u

# Any output ID in server with no UI match = orphan → delete

# Functions defined but never called
grep -rn "^[a-zA-Z_][a-zA-Z0-9_.]*[[:space:]]*<-[[:space:]]*function" modules/
# Cross-reference against call sites — any with 0 call sites → delete
```

### 9.3 Dynamic-Only Final Checklist (per file, before commit)

- [ ] No `$seurat_clusters`, `$cell_type`, `$orig.ident` as `$literal` — always `[[col_var]]`
- [ ] No species strings (`"human"`, `"hsa"`, `"org.Hs.eg.db"`) outside the mapping helpers in `de_core.R` / `pn_core.R`
- [ ] No gene name literals (`c("CD3D", "CD8A", ...)`) anywhere except test fixtures
- [ ] No threshold literals (`0.05`, `0.25`, `2000`) outside UI default values in `*_ui.R` files
- [ ] No `set.seed(42)` — always `set.seed(input$xxx_seed)` from UI
- [ ] No `library(pkg)` inside module files — only `requireNamespace()` with informative `stop()` on failure
- [ ] No `head(x, 20)` or fixed display counts — always from UI input
- [ ] No hardcoded file paths — all resolved via `file.path()` relative to app root
- [ ] No `DefaultAssay(obj) == "RNA"` string comparisons — always `DefaultAssay(obj)` dynamically

---

## 10. Success Criteria

1. All 9 DE methods run to completion on a 25k-cell dataset without error
2. Results are bit-for-bit identical across two runs with the same seed
3. ORA and GSEA run correctly against all 8 databases for human and mouse
4. Gene universe is always set to all detected genes (never default to database genes)
5. No code duplication — each analysis function exists in exactly one file
6. `values$last_enrich_results` is always populated before the Pathway Network tab renders
7. All old files listed in §2.3 are deleted from the codebase
8. `install_dependencies.R` updated with all required packages
9. §9.1 dead-code checklist passes for every modified file
10. §9.3 dynamic-only checklist passes for every modified file
