# scRNA15 — Production Roadmap

**Date:** 2026-05-11
**Status:** Approved (user decisions captured 2026-05-11)
**Author:** Architecture review
**Scope:** End-to-end plan to take the platform from "runs locally" to "production multi-tenant scientific platform."

## Locked-in decisions (2026-05-11)

| # | Decision | Implication |
|---|---|---|
| 1 | **Deployment**: Docker-first, supports either a 128 GB-RAM VPS or GCP (Cloud Run / GKE). User will pick later. | Plan builds artefacts that work on both; tilt config to "scale-vertically on one box" by default. |
| 2 | **AI model**: user-chosen per session. Opus, Sonnet, Haiku selectable in the AI panel. | Add a model picker in P3-1; cost-attribution logged per call. |
| 3 | **Spatial transcriptomics**: deferred to v2. | P5-7 removed from this roadmap. |
| 4 | **CellTypist via reticulate**: OK. | P4-3 keeps the CellTypist row; reticulate becomes a hard dep. |
| 5 | **Data retention**: **1 month** of inactivity → auto-archive. | Adds a cron + archival path (see new P7-3). |
| 6 | **Tier plans** with dataset-count caps: **Free = 1**, **Basic = 5**, **Standard = 10**, **Pro = 50**. | Adds a full plans/quotas subsystem — promoted from P7-4 to its own phase (now **Phase 4.5**) because it touches auth, UI, persistence, and admin tooling. |

---

## Executive summary

The platform is **functionally complete and architecturally sound**. The gaps preventing production deployment are *operational*, not *scientific* — concurrency, persistence, async cancellation, RAM ceilings, and unified AI. This document organizes those gaps into 7 phases over ~3 months. Every item has: problem statement, scientific/engineering solution, files affected, acceptance criteria, and effort. The plan is biased toward **reproducibility, robustness, and academic rigor** over feature breadth.

**Effort estimates use this scale:**
- S = ½–1 day · M = 2–4 days · L = 1–2 weeks · XL = 3+ weeks

---

## Current-state inventory

| Layer | Status |
|---|---|
| Auth (SQLite + bcrypt + 8-char min + WAL) | ✅ baseline; missing rate-limiting + 2FA |
| Project persistence (RDS in `data/user_sessions/<uid>/<pid>/analysis.rds`) | ⚠ no locking, no schema versioning, no metadata/UI-state save |
| Modules (19, ~50 700 LOC) | ✅ modular; ❌ 14/19 untested |
| Async (`ExtendedTask` for AI only) | ⚠ 77 other long-running ops are blocking |
| Memory (all in-RAM; `future.globals.maxSize=16 GB`) | ❌ collapses past ~50–80k cells |
| Gallery (figures + CSVs) | ⚠ CSVs registered but no download buttons in UI |
| Annotation methods | ⚠ ScType + PanglaoDB; SingleR not in dropdown; no scoring exposure |
| Trajectory | ⚠ Slingshot + Monocle 3 + DPT + velocity stubs; 110 KB monolithic |
| CCC | ⚠ CellChat path slow on >25k cells; no async, no progress cancel |
| Subclustering | ⚠ implemented but not first-class (no hierarchy tracking, no merge-back) |
| AI features | ⚠ 6 separate AI entry points, no shared context, no tool-use loop |
| Deployment artifacts | ❌ none (no Dockerfile, renv.lock, CI, .gitignore, README) |

---

# PHASE 0 — Pre-production blockers (Week 1)

> The work that has to happen before anyone other than the developer runs the app. Pulled forward from the May-10 review because it's required by everything downstream.

## P0-1 · Environment freeze: `renv` + Dockerfile
**Problem:** `install_dependencies.R` installs latest of everything → the msigdbr `category→collection` rename broke us at runtime.

**Solution:**
- `renv::init(bare=TRUE)` → `renv::snapshot()` → commit `renv.lock`
- Multi-stage Dockerfile on `rocker/r-ver:4.5.0` that runs `renv::restore()` at build time
- `.dockerignore` excludes `data/`, `omnipathr-log/`, `sessions/`

**Files:** new `renv.lock`, `Dockerfile`, `.dockerignore`, `.gitignore`

**Acceptance:** `docker build -t scrna15:0.1.0 . && docker run -p 3838:3838 scrna15:0.1.0` boots to login page in <90 s.

**Effort:** S

## P0-2 · `.gitignore` + `README.md`
**Problem:** `data/scrna4_auth.db` and `data/user_sessions/*/analysis.rds` are inside the project root.

**Files:** new `.gitignore` excluding `data/`, `*.rds`, `*.db`, `.Rhistory`, `.RData`, `.DS_Store`, `omnipathr-log/`, `sessions/`, `renv/library/`; new `README.md` with quickstart + module map.

**Acceptance:** `git status` after a fresh clone shows zero data/secret files.

**Effort:** S

## P0-3 · GitHub Actions CI
**Problem:** No regression net. The sprintf-percent bug shipped because headless smoke tests didn't call UI functions.

**Solution:** `.github/workflows/ci.yml`:
```yaml
jobs:
  lint:    parse() every R file; testthat::test_dir; check for inline hex / hard-coded species (regression guards already exist in tests)
  ui-render: source app.R + as.character(deUI('de')) + as.character(pnUI('pn')) → fail on any sprintf-style error
  docker:  docker build to validate Dockerfile
```

**Acceptance:** CI green on a clean push; CI red when a literal `"Homo sapiens"` is added outside `pn_core.R`.

**Effort:** S

## P0-4 · Sanitized error handler + structured logging
**Problem:** `options(shiny.error = ...)` shows raw R errors in `showNotification` → XSS risk and unprofessional UX.

**Solution:**
- Add `logger` to deps; init in `app.R` start with rotating file appender to `logs/app.log`
- Wrap global handler: log full trace server-side with a UUID; show the user `"Internal error — please report ID {uuid}"`
- Map known error patterns → friendly messages (table-driven in `core_robustness.R`)

**Files:** `app.R` (error block), new `modules/util/logger_setup.R`, `core_robustness.R` (error-map table).

**Acceptance:** Triggering `stop("<script>alert(1)</script>")` from any module logs the trace and shows a sanitized message in UI.

**Effort:** M

## P0-5 · Login rate-limiting + lockout
**Problem:** `authenticate_user()` has no attempt cap.

**Solution:** new `login_attempts(username, ip, ts)` table; 5 failures in 15 min → reject for 15 min with timing-safe response. Log to `logs/auth.log`.

**Files:** `modules/auth/auth_db.R`, `modules/auth/auth_module.R`.

**Acceptance:** Automated test (with `httr2`) — sixth failure returns 429-equivalent + delay.

**Effort:** M

---

# PHASE 1 — Persistence + resume + cancel (Weeks 2–3)

> The four user-visible asks: save full UI state, resume where the user left off, cancel without crash, run two projects in parallel.

## P1-1 · Full project state snapshot (schema-versioned)
**Problem:** `save_project_state()` saves data objects but **not UI parameters, DE history, enrichment history, AI conversations, last-viewed tab, or thresholds**. User reopens project → starts from scratch.

**Solution:** expand the save bundle with a versioned schema:

```r
project_bundle <- list(
  schema_version    = 3L,
  app_version       = packageVersion(...),
  saved_at          = Sys.time(),

  # Data
  seurat_list       = ...,
  active_dataset    = ...,                   # large — see P2-1 for on-disk option

  # Analysis state
  de_history        = values$de_history,
  enrich_history    = values$enrich_history,    # NEW — mirror DE history
  cmp_history       = values$cmp_history,       # NEW
  grn_results       = values$grn_results_last,
  trajectory_state  = ...,
  ccc_results       = ...,

  # UI state (what was selected, what thresholds, which tab)
  ui_state = list(
    active_tab           = input$tabs,
    de_inputs            = .snapshot_de_inputs(input),
    enrich_inputs        = .snapshot_enrich_inputs(input),
    cmp_inputs           = .snapshot_cmp_inputs(input),
    network_inputs       = ...,
    annotation_inputs    = ...,
    selected_genes       = ...,
    custom_palettes      = values$cluster_colors,
    last_viewed_figures  = ...
  ),

  # AI
  ai_conversations      = values$ai_conversations,
  ai_context_summary    = values$ai_context_summary,

  # Reproducibility
  analysis_audit        = values$analysis_audit,
  session_info          = .capture_session_info(),

  # Resume pointer
  resume_state = list(
    last_completed_step  = "differential_expression",
    in_progress_step     = NULL,
    progress_percent     = NULL
  )
)
```

Plus a migration function `.migrate_bundle(b)` that handles schema v1 → v2 → v3.

**Files:** `modules/auth/auth_db.R` (save/load), `app.R` (UI-state restore on session start), new `modules/util/state_snapshot.R`.

**Acceptance:** User runs DE with custom thresholds → logs out → logs back in → DE tab pre-populated with same thresholds + same results visible. Schema version reported in audit log.

**Effort:** L

## P1-2 · "Resume" banner + step tracker
**Problem:** Even after P1-1, the user doesn't know where they were.

**Solution:**
- Track the pipeline graph in `values$analysis_steps`:
  `data_loaded → qc_done → preprocessing_done → dimred_done → clustering_done → annotation_done → de_done → enrichment_done`
- On project reopen, show a banner at the top of the main panel: *"Welcome back. Last completed: Pathway Enrichment (epithelial cluster). Continue?"* with a button that updates `input$tabs` to the next step
- Visual breadcrumb in the sidebar showing completed steps in green, in-progress in amber, pending in gray

**Files:** new `modules/util/resume_banner.R`, `app.R`.

**Acceptance:** Banner appears with correct last-step text after reload; click forwards to the right tab; dismiss-once-per-session option.

**Effort:** M

## P1-3 · Cancellable async tasks — kill the "crash on stop" problem
**Problem:** 77 blocking `withProgress` calls. User has no way to abort a 30-minute DE run without killing the whole app.

**Solution:** introduce a uniform task abstraction in `modules/util/task_manager.R`:

```r
# Task launcher
run_cancellable_task(values, task_id, label, fn, on_success, on_error, on_cancel) {
  # Creates an ExtendedTask, registers in values$task_queue
  # Each task carries: id, label, status (queued/running/done/cancelled/error),
  # start/end, progress%, kill_flag, result, owner_project_id
}

# UI: status dock at the bottom of every page
taskDockUI()  →  pill list of active tasks with × button each
```

Every heavy analysis (DE, enrichment, comparative, GRN, WGCNA, trajectory, CCC, integration) wraps its core call:

```r
# Before:
withProgress({ run_de_analysis(...) })

# After:
run_cancellable_task(values, "de", label = .auto_de_label(...),
  fn = function(progress, kill_flag) {
    run_de_analysis(..., on_progress = progress, kill_flag = kill_flag)
  },
  on_success = function(res) { de_results(res$data); ... },
  on_cancel  = function() showNotification("DE cancelled", type = "warning")
)
```

Engines accept a `kill_flag` reactive and check it at safe points (per-gene loop in DE; per-pathway loop in fgsea wrapper; per-block in WGCNA). On `kill_flag()` → clean return with `attr(result, "cancelled") = TRUE`.

**Files:** new `modules/util/task_manager.R`, modifications to: `de_methods.R`, `pn_enrichment.R`, `pn_network.R`, `downstream_ccc.R`, `downstream_tf.R`, `preprocessing.R` (integration), `clustering.R`, `trajectory.R`.

**Acceptance:** Start a DE run → click Cancel mid-run → result panel shows "Cancelled at gene 1432 / 5400", app remains responsive, no error notification.

**Effort:** L

## P1-4 · Multi-project queue + concurrent execution
**Problem:** Running an analysis in project A blocks project B's session.

**Solution:**
- Per-Shiny-session ID is already isolated; the bottleneck is the **R worker count**. Configure `simple_scheduler 4` in `shiny-server.conf` (or `mirai::daemons(4)` for self-hosted) so 4 user sessions run truly in parallel.
- Within a session, queue tasks: `values$task_queue` is a FIFO; only N=2 tasks run concurrently per project (configurable). Surplus tasks show "Queued" status.
- Persist queue state to project bundle (P1-1) so a killed Shiny worker can recover the queue on next load.
- **Hard locks**: when a task is operating on `values$active_dataset`, it acquires an in-memory lock; other tasks queue. Cross-session locks via SQLite advisory row in `data/scrna4_auth.db`.

**Files:** new `modules/util/task_manager.R` (extended from P1-3), `modules/auth/auth_db.R` (advisory locks).

**Acceptance:** Two browser tabs on different projects can run DE simultaneously; one user's load doesn't block the other.

**Effort:** L

## P1-5 · Gallery CSV download button (clean fix)
**Problem:** [modules/umap_gallery.R:517-540](modules/umap_gallery.R) renders `csv_table` with `Description / Module / Rows / Columns / Time / File` — but no per-row Download button.

**Solution:**
- Add an Action column rendered as a Shiny actionButton via DT's `escape = FALSE` + JS callback
- Wire a `downloadHandler` per CSV id that reads `lib[[id]]$filename` and streams it back
- Bulk action: "Download all as ZIP" (already have `zip` in deps)
- Also add: re-name, delete, push-to-S3 (P5 placeholder)

**Files:** `modules/umap_gallery.R` (CSV Library tab body + server).

**Acceptance:** Every row has a working Download button; "Download all" produces a ZIP with sensible filenames.

**Effort:** S

---

# PHASE 2 — Memory + scale for 100k+ cells (Weeks 4–5)

> Right now everything is in-RAM Seurat. 100k-cell datasets OOM on typical 32 GB machines.

## P2-1 · On-disk storage with BPCells
**Problem:** A 200k-cell sparse count matrix is ~6 GB in memory; with `data` + `scale.data` layers and PCA it's 15+ GB. Most academic machines don't have it.

**Solution:** integrate Seurat v5's **BPCells** backend.

- New "Storage mode" radio in Data Upload: `Auto (recommended)` | `In-memory` | `On-disk (BPCells)`
- Auto rule: file size × 3 > 50 % of `memuse::Sys.meminfo()$totalram` → use BPCells
- On-disk storage path: `data/user_sessions/<uid>/<pid>/bpcells/<dataset>/`
- All downstream modules already use `LayerData()` / `safe_get_assay_data()` → BPCells is transparent to them. Verify each module.
- Update `prepare_for_de()` to handle BPCells layers (already supports Seurat v5 layered assays — verify counts → data conversion works on-disk).
- Add explicit DESeq2/edgeR pseudobulk path that **materializes only the per-sample aggregate** (small matrix) and never the full counts.

**Files:** `modules/data_input.R`, `modules/preprocessing.R`, `modules/core_robustness.R` (extend `safe_get_assay_data`), all engine helpers.

**Acceptance:** Load a 200k-cell dataset on a 16 GB machine; QC + UMAP + clustering complete without OOM; peak RAM <8 GB.

**Effort:** L

## P2-2 · Memory pre-flight checks
**Problem:** No guard before heavy ops. WGCNA on 5000 genes = 200 MB TOM but we don't check available RAM before allocating.

**Solution:** add `mem_preflight(required_gb, op_name)` helper in `core_robustness.R`. Every heavy op (integration, WGCNA TOM, pseudobulk DESeq2, fgsea on giant gene sets) calls it. If insufficient: warn the user, offer downsample / on-disk / smaller-gene-set fallback.

**Acceptance:** Trying to run WGCNA on 8000 genes with 4 GB free → clear error message *"Estimated TOM memory 0.5 GB exceeds 40 % of free RAM. Try ≤4000 genes."*

**Effort:** S

## P2-3 · Async + chunked compute for the 8 heaviest paths
**Problem:** Even on-disk, some ops are slow. The user wants progress + cancel.

**Specifically convert to async + chunked:**
1. Integration (CCA / RPCA / Harmony) — already supports BPCells via Seurat v5
2. WGCNA TOM (chunk genes per block)
3. fgsea on large gene sets (already vectorized; just wrap in promise)
4. CCC CellChat (P4-4 below)
5. Pseudobulk DESeq2 (already fast; just async)
6. Live-API enrichment (KEGG / Reactome HTTP — async with timeout)
7. WGCNA `pickSoftThreshold` (slow on big data)
8. SingleR scoring (per-cell loop)

All use the `task_manager` from P1-3.

**Effort:** L (1 day per path × 8)

## P2-4 · Disk-backed AI conversation index (RAG prep)
**Problem:** AI context can grow large; in-memory keeps eating session RAM.

**Solution:** store AI conversations in SQLite (`ai_conversations` table) with one row per turn; embed via local sentence-transformer (`sentence-transformers/all-MiniLM-L6-v2` via reticulate) and store vectors in `duckdb`. Used by the unified AI in Phase 3.

**Effort:** M

---

# PHASE 3 — Unified AI + tool use + RAG (Weeks 6–7)

> Six scattered AI buttons → one assistant that knows context, can act on the data, and answers questions with retrieved evidence.

## P3-1 · Single AI assistant entry point (user-chosen model)
**Problem:** Six AI features (AI auto-DE, AI test recommender, AI pathway interpreter, AI annotation, AI plot builder, AI assistant) — each fires independently, with no shared context.

**Solution:** consolidate into one persistent panel.

- New floating "AI Assistant" button (bottom-right of every page)
- Click → slide-in chat panel (right side, 420 px)
- Chat is persistent across tabs and across project reloads (P1-1 saves history)
- The same chat handles every previous AI function via **tool use** (P3-2)
- Top bar of the panel: project name · current cell count · last analysis stage (gives the AI context for free)
- **Model picker** in the panel header — user selects per-session: `Opus` (best quality, slow, expensive), `Sonnet` (default, balanced), `Haiku` (fast, cheap). Selection saved to project bundle.
- **Token + cost meter** below the picker — running total of tokens/USD this session (read from Anthropic response headers); resets per session; useful for capped academic budgets
- **Free-tier fallback**: if the user has no API key and is on the Free plan (see Phase 4.5), the AI panel is replaced with an upgrade prompt — no silent failure

**Files:** redesign `modules/ai_assistant.R` as the **only** AI entry point; delete the AI button blocks in `annotation.R`, `visualization.R`, `export.R` (already mostly done for DE/enrichment).

**Effort:** M

## P3-2 · Anthropic tool-use loop
**Problem:** AI can only generate text; it can't act.

**Solution:** implement a tool-use agent loop using Anthropic's Tool Use API.

Tools the AI can call:

| Tool | What it does | Result |
|---|---|---|
| `get_dataset_summary()` | reads `values$active_dataset` | n_cells, n_genes, clusters, species, ... |
| `list_de_runs()` | reads `values$de_history` | named list with summaries |
| `run_de(params)` | calls `run_de_analysis()` | saves to history → returns top genes |
| `run_enrichment(params)` | calls `run_ora/gsea_analysis` | returns top pathways |
| `update_input(input_id, value)` | calls `updateXxxInput()` server-side | confirms; UI changes |
| `switch_tab(tab_id)` | `updateTabItems()` | confirms |
| `lookup_gene(gene)` | reads expression matrix + DE results | per-cluster expr + DE stats |
| `search_pubmed(query)` | (if PubMed MCP available) | top 5 papers |
| `query_pathway(name)` | reads gene sets via msigdbr | member genes + description |
| `save_to_gallery(plot_id)` | calls `capture_figure()` | confirms |

The loop continues until the AI emits a final user-facing answer.

**Files:** rewrite `modules/ai_assistant.R::call_claude_api` to handle `tool_use` responses; new `modules/ai_assistant_tools.R` with each tool's R impl.

**Acceptance:** User: "What are the top markers of cluster 5 and which pathways are they in?" → AI calls `list_de_runs` → finds nothing → calls `run_de(group_by=cluster, ident1=5, ident2=rest)` → calls `run_enrichment(database=H)` → reports the result.

**Effort:** L

## P3-3 · RAG over project artifacts
**Problem:** The AI doesn't know what's in this project's DE/enrichment results without explicit calls.

**Solution:** index project artifacts as embeddings.

- On every successful analysis run, generate a textual summary ("DE of mutant vs wt within osteoblast cluster found Hes1 down (FC=-2.1, padj=1e-8), ...") and embed it
- Store embeddings + summaries in `data/user_sessions/<uid>/<pid>/rag.duckdb`
- On every AI query, do top-K cosine search and prepend matches to the system prompt
- Hybrid retrieval: also include the latest 3 audit logs verbatim regardless of similarity
- Embedding model: `BAAI/bge-small-en-v1.5` via reticulate (or Anthropic embedding API if Anthropic ships one by then)

**Files:** new `modules/util/rag.R`, integration hooks in DE/enrichment/network run handlers.

**Acceptance:** User asks "Did we already test Hes1?" without specifying which run — the AI retrieves the right DE result via RAG and answers.

**Effort:** L

## P3-4 · Persistent conversation + provenance
**Problem:** Conversations disappear on session end.

**Solution:**
- `ai_conversations` table in SQLite (one row per turn: id, project_id, role, content, tool_calls, tool_results, ts)
- New "AI History" tab in the project dashboard
- Every assistant message gets a provenance footer linking to the analyses it cited (clickable → opens the figure in the gallery)
- Export conversation as Markdown

**Effort:** M

---

# PHASE 4 — Analysis depth (Weeks 8–10)

> The user-flagged scientific gaps: subclustering, trajectory, annotation, CCC.

## P4-1 · Subclustering rebuild — first-class hierarchical clustering
**Problem:** [clustering.R:2004](modules/clustering.R) has a subclustering panel but it's flat — no parent/child tracking, no merge-back, no annotation per level.

**Solution:** new module `modules/subcluster/`:

- **Cluster tree data structure**: `values$cluster_tree` is a tree with each node = `{cluster_id, parent_id, cells, n_cells, label, label_method, depth}`
- **Subcluster a node**: select node → re-run FindVariableFeatures + PCA + Leiden on subset → child nodes appended
- **Merge-back labels**: child labels propagate to a hierarchical column `cluster_L0`, `cluster_L1`, ... in `obj@meta.data`
- **Tree visualization** using `collapsibleTree` (small dep) — clickable node opens the cluster's per-node QC + UMAP highlight
- **Auto-annotation per node**: runs the annotation pipeline (ScType + PanglaoDB + optional SingleR) on each subcluster
- **Reproducibility**: every subcluster operation saved to `values$subcluster_history` with seed, resolution, n_PCs
- **Persistence**: tree saved to project bundle (P1-1)

**Scientific contract:**
- Variable feature selection is **re-done on the subset** (vital — global HVGs may be dominated by other lineages)
- PCA is **re-computed** on the subset (also vital)
- Optional: integration is re-applied (Harmony per subcluster) for cross-condition subclustering

**Files:** new `modules/subcluster/{subcluster_core,subcluster_ui,subcluster_server}.R`; integration in `clustering.R` to delegate.

**Acceptance:** Click "T cells" → "Subcluster" → 5 child clusters appear → annotate → see "CD4 naive / CD8 effector / Treg / NK / gd-T" propagated to the main UMAP under `cluster_L1`.

**Effort:** L

## P4-2 · Trajectory rebuild — modular, audited, async
**Problem:** [trajectory.R](modules/trajectory.R) is 110 KB monolithic.

**Solution:** split into `modules/trajectory/`:

```
trajectory_core.R      # Seurat ↔ SingleCellExperiment ↔ CellDataSet converters; root-cell helpers
trajectory_methods.R   # 4 engines: slingshot, monocle3, dpt (destiny), velocity (velociraptor)
trajectory_ui.R        # split-screen UI matching DE/PN style
trajectory_server.R    # async wiring + gallery integration
```

**Scientific contract:**
- Every method outputs to a canonical schema: `data.frame(cell, pseudotime, lineage, prob)` + optional `trajectory_graph` (igraph) + audit log
- Methods are run via the same `task_manager` as DE → cancellable, gallery-wired
- **RNA velocity via `velociraptor`** (Bioconductor R port of scVelo) — closes the Fig 7 gap from the Notch2 paper
- **Root-cell selection**: by cluster picker (default), by gene-expression-top-N (e.g. "cells highest in CD34"), or by manual cell list
- **Branch comparison**: when ≥2 lineages, automatic per-branch DE + enrichment via the existing engines (a button in the trajectory tab → populates DE history with lineage labels)

**Files:** new directory `modules/trajectory/`, delete `modules/trajectory.R` after migration.

**Acceptance:** Run Slingshot + Monocle 3 on the same dataset → side-by-side pseudotime distributions per cluster; cancel mid-run works.

**Effort:** XL

## P4-3 · Annotation expansion + SingleR with full scoring
**Problem:** Only ScType + PanglaoDB exposed. SingleR is in deps but not wired. PanglaoDB scoring is naive (top-marker overlap, no statistical test).

**Solution:** annotation method dropdown expands to:

| Method | Backend | Output | Use case |
|---|---|---|---|
| **ScType** (current) | local Excel DB | per-cluster top label + score | Fast, no internet |
| **PanglaoDB** (current) | local TSV | per-cluster top label + score | Tissue-specific |
| **SingleR (full)** | `SingleR::SingleR()` + `celldex` references | per-cell label + delta score + top genes per label + heatmap | Per-cell resolution |
| **Azimuth** | `Azimuth::RunAzimuth()` | per-cell label + L1/L2/L3 hierarchy | Reference atlas mapping |
| **scGate** | `scGate::scGate()` | per-cell binary gates (marker-based) | Strict marker-based purity |
| **CellTypist** (optional) | reticulate → Python | per-cell label + confidence | Cross-tissue reference |
| **Consensus** | majority vote across selected methods | per-cell label + agreement % | Robust |

**SingleR scoring exposure:**
- Show the delta score histogram per cluster (`scater::plotScoreHeatmap`-equivalent in plotly)
- Show per-cluster top-3 candidate labels with their scores (not just the winner)
- Flag low-confidence assignments (delta < 0.1) for manual review

**Use the existing PanglaoDB asset**: the directory `PanglaoDB/` has `PanglaoDB_markers_27_Mar_2020.tsv` — currently loaded but only used for overlap counting. Add a hypergeometric test for marker enrichment per cluster, with multiple-testing correction across clusters.

**Files:** `modules/annotation.R` (extend dropdown + scoring renderers), new `modules/cell_type_markers.R` additions for SingleR / Azimuth / scGate adapters.

**Effort:** L

## P4-4 · CCC speed-up + cancellation
**Problem:** CellChat on 25k+ cells takes 5–30 min and can't be cancelled.

**Solution:**

1. **Async + cancel** via `task_manager` (P1-3). User clicks Cancel → CellChat permutation loop breaks cleanly.
2. **L-R pair pre-filter**: drop pairs where neither ligand nor receptor is expressed >5% in any cluster before permutation. Typically cuts pair count by 80%.
3. **Per-pair parallelism**: CellChat already supports `future::plan(multisession)` — set workers = available cores − 1.
4. **Subsample option** (explicit): "Use up to N cells per cluster". Defaults to Inf; user sets 1000 for exploration.
5. **Cache by dataset signature**: if user re-runs with same params on same `values$dataset_signature`, return cached result.
6. **Alternative fast engine — NICHES or scLink**: add a "Fast (no permutations)" method that just scores L-R co-expression without significance testing. Useful for exploration before the full run.

**Files:** `modules/downstream/downstream_ccc.R` extensive.

**Acceptance:** CellChat on the 27k-cell test dataset completes in <8 min with default params; cancel-mid-run works; cache returns in <1 s on re-run.

**Effort:** L

---

# PHASE 4.5 — Tier plans, quotas, retention (Week 10–11, parallel-track)

> User-confirmed business model: 4 plan tiers gated by dataset count, 1-month inactivity retention. Promoted from a small line-item in P7-4 to its own phase because it touches **auth, persistence, every upload path, every analysis run, an admin dashboard, and a cron**.

## P45-1 · Plans + quota schema
**Problem:** No notion of plan tier in the auth DB. Every user has unlimited datasets.

**Solution:** schema additions to `data/scrna4_auth.db`:

```sql
-- Plans (seeded, immutable after first migration)
CREATE TABLE IF NOT EXISTS plans (
  id            TEXT PRIMARY KEY,            -- "free", "basic", "standard", "pro"
  name          TEXT NOT NULL,
  max_datasets  INTEGER NOT NULL,             -- 1, 5, 10, 50
  max_ram_gb    INTEGER NOT NULL,             -- per-run cap; see P45-3
  max_jobs_day  INTEGER NOT NULL,             -- daily analysis cap
  ai_calls_day  INTEGER NOT NULL,             -- AI calls per day
  retention_days INTEGER NOT NULL DEFAULT 30, -- archive after N days inactivity
  priority      INTEGER NOT NULL,             -- queue priority (lower = first)
  price_monthly_usd REAL DEFAULT 0,
  description   TEXT
);

INSERT OR IGNORE INTO plans VALUES
  ('free',     'Free',     1,  16,  5,   20,  30, 100, 0,     '1 dataset, 5 jobs/day'),
  ('basic',    'Basic',    5,  32,  20,  100, 30, 50,  9.99,  '5 datasets'),
  ('standard', 'Standard', 10, 64,  50,  300, 30, 20,  19.99, '10 datasets'),
  ('pro',      'Pro',      50, 128, 200, 1000,30, 10,  49.99, '50 datasets, top priority');

ALTER TABLE users ADD COLUMN plan_id TEXT NOT NULL DEFAULT 'free' REFERENCES plans(id);
ALTER TABLE users ADD COLUMN plan_started TIMESTAMP DEFAULT CURRENT_TIMESTAMP;
ALTER TABLE users ADD COLUMN plan_expires TIMESTAMP DEFAULT NULL;     -- NULL = no expiry (free)
ALTER TABLE users ADD COLUMN jobs_today INTEGER NOT NULL DEFAULT 0;   -- reset by cron
ALTER TABLE users ADD COLUMN ai_calls_today INTEGER NOT NULL DEFAULT 0;
ALTER TABLE users ADD COLUMN last_quota_reset DATE DEFAULT CURRENT_DATE;

ALTER TABLE projects ADD COLUMN last_accessed TIMESTAMP DEFAULT CURRENT_TIMESTAMP;
ALTER TABLE projects ADD COLUMN archived_at   TIMESTAMP DEFAULT NULL;  -- non-NULL = archived
ALTER TABLE projects ADD COLUMN size_bytes    INTEGER DEFAULT NULL;    -- updated on save
```

All new columns have defaults so the migration is backward-compatible. Migration runs once on first boot of `>= 1.0.0`.

**Files:** new `modules/auth/auth_db_migrations.R` with a list of versioned migrations; `auth_db.R` calls the migrator at init.

**Acceptance:** Migration runs idempotently on existing DB; `SELECT plan_id, max_datasets FROM users JOIN plans ON users.plan_id = plans.id LIMIT 5` works.

**Effort:** S

## P45-2 · Quota enforcement (the actual gates)
**Problem:** No checks anywhere.

**Solution:** four enforcement points, all using a single `check_quota(user_id, what, n)` helper in `auth_db.R`:

1. **Dataset count at upload**:
   - `data_input.R` upload handler → before adding to `values$seurat_list`, call `check_quota(uid, "datasets", current_count + 1)`
   - On fail: notification *"Free plan allows 1 dataset. Upgrade to Basic (5) or higher to add more."* + Upgrade button
2. **Daily job count at run**:
   - `task_manager.R` `run_cancellable_task` → before queuing, increment `users.jobs_today`; if `> plan.max_jobs_day` → reject with notification
   - Cron resets daily at midnight UTC (see P45-5)
3. **AI calls per day**:
   - `ai_assistant.R::call_claude_api` → increment + check `users.ai_calls_today`; if `>=` plan cap → block with "AI daily limit reached"
4. **Per-run RAM cap**:
   - Extended `mem_preflight()` from P2-2 → also rejects if estimated RAM > `plan.max_ram_gb`

All four gates write to `logs/quota.log` on every block (for support + abuse detection).

**Files:** `modules/auth/auth_db.R` (helper), `modules/data_input.R`, `modules/util/task_manager.R`, `modules/ai_assistant.R`, `modules/core_robustness.R`.

**Acceptance:** A Free user with 1 dataset tries to upload a 2nd → blocked with helpful message. A Basic user runs the 21st job in 24h → blocked. Both events appear in `quota.log`.

**Effort:** M

## P45-3 · Account / plan management UI
**Problem:** Users need to see their plan + usage + upgrade.

**Solution:** new tab "Account" in the project dashboard:

- **Header card**: current plan name + price + renewal date
- **Usage strip**: 4 progress bars (datasets / jobs today / AI calls today / disk used) — each color-coded green/amber/red based on percent of cap
- **Plan-comparison table** showing all 4 tiers with their limits, current plan highlighted, "Upgrade" / "Downgrade" / "Cancel" buttons
- **Billing**: deferred to a later phase OR integrated via Stripe (`stripr` R wrapper). For now, "Contact us to upgrade" goes to an admin queue.
- **Admin path**: superuser flag in `users` table → admin can set plan manually via an admin dashboard (P45-6)

**Files:** new `modules/account/` directory with `account_ui.R`, `account_server.R`.

**Effort:** M

## P45-4 · 1-month inactivity archival
**Problem:** User said "1 month storage" — projects untouched for 30 days should auto-archive.

**Solution:**
- A **nightly cron** (system cron or `cronR` in the app) runs `archive_inactive_projects()`:
  - `SELECT * FROM projects WHERE archived_at IS NULL AND (julianday('now') - julianday(last_accessed)) > 30`
  - For each: gzip the project bundle → move to `data/archives/<uid>/<pid>.rds.gz`; set `archived_at = now()`
  - Send notification email to user (P7-5) 7 days before archival
- **Restore path**: archived projects still listed on the project dashboard with a "Restore" button (within 30 days of archival → free; after 60 days → premium plans only). After 90 days → permanent deletion.
- Every `load_project_state()` updates `last_accessed = CURRENT_TIMESTAMP` (so opening keeps it fresh).

**Files:** new `modules/util/retention_cron.R`; modifications to `auth_db.R` (touch on load) + `project_module.R` (Restore button).

**Acceptance:** Test: create a project, manually `UPDATE projects SET last_accessed = datetime('now','-31 days')`, run the cron → project is archived, file gzipped, restorable.

**Effort:** M

## P45-5 · Daily counter reset cron
**Problem:** `users.jobs_today` and `users.ai_calls_today` need to reset every midnight UTC.

**Solution:** part of the same cron as P45-4. Single SQL: `UPDATE users SET jobs_today=0, ai_calls_today=0, last_quota_reset=CURRENT_DATE WHERE last_quota_reset < CURRENT_DATE`.

Also a safety check on every `check_quota` call: if `users.last_quota_reset < CURRENT_DATE`, do an inline reset before checking (handles cases where cron hasn't fired yet).

**Effort:** S

## P45-6 · Admin dashboard
**Problem:** No way to manage users / plans / abuse.

**Solution:** new sidebar item "Admin" (visible only to `users.is_admin = 1`):
- Users table: filter by plan, see last login, jobs today, datasets count
- Set plan + plan_expires per user
- Force-archive a project
- View quota.log + auth.log tails
- Active sessions list with "Force logout"
- AI usage report: tokens + USD by user (read from logged Anthropic API responses)

**Effort:** L

---

# PHASE 5 — Scientific capabilities (Weeks 12–13)

> New features that take the platform from "good" to "comprehensive."

## P5-1 · Batch-effect validation (kBET / LISI / silhouette)
**Problem:** Integration runs but never reports *whether it worked*.

**Solution:** new "Integration QC" sub-tab in Preprocessing:
- **kBET** (rejection rate by batch)
- **iLISI / cLISI** (integration vs cell-type LISI)
- **Silhouette** before / after integration, faceted by batch and by cell type
- Side-by-side UMAP colored by batch, with quantitative score below

Auto-runs after integration; results saved to project bundle.

**Effort:** M

## P5-2 · Differential abundance testing (propeller + miloR)
**Problem:** "35–45 % drop in osteoblast cluster" is a percentage, not a p-value.

**Solution:** new "Differential Abundance" sub-tab under Pathway & Network → Composition:
- `speckle::propeller()` — beta-regression on cluster proportions (already in deps)
- `miloR::testNhoods()` (optional) — neighborhood-level (sub-cluster-resolution) abundance test
- Output: per-cluster log-fold-change + p-value, volcano-style plot

**Effort:** M

## P5-3 · Cell-cycle + sex-chromosome QC scoring
**Problem:** No automatic check for cell-cycle confounders or sex mismatches.

**Solution:** add to QC module:
- `Seurat::CellCycleScoring()` with built-in gene lists (Tirosh 2016) — auto-runs after normalization
- Sex check via XIST + Y-chromosome gene scoring
- Both visualized on UMAP

**Effort:** S

## P5-4 · Methods report generator
**Problem:** Users have to write the methods section manually.

**Solution:** new button "Generate Methods Report" in Export module:
- Walks `values$analysis_audit` chronologically
- Renders to Markdown + Word using the saved templates
- Includes: package versions (already captured), seed values, parameter values, sample counts at each filter step
- Output drops into Figure Gallery as a downloadable .docx + .md

**Effort:** M

## P5-5 · Public-dataset import (GEO + CellxGene)
**Problem:** Users have to download data manually.

**Solution:** new "Import from Public Repository" tab in Data Upload:
- Enter GSE accession → uses `GEOquery` to fetch SOFT + supplementary files
- Enter CellxGene Census ID → uses `cellxgene.census` (Python via reticulate or `cellxgene.census.r` if available)
- Auto-detects 10x vs h5ad format

**Effort:** M

## P5-6 · h5ad / loom export for Python interop
**Problem:** Users want to take results to scanpy / scVI.

**Solution:** Export module adds:
- `zellkonverter::writeH5AD()` of the current Seurat object
- `SeuratDisk::SaveLoom()` (already in deps)
- All `values$de_history` entries exported as `obs/varm` annotations

**Effort:** S

## P5-7 · Spatial transcriptomics support — **deferred to v2.0**
Not in scope for this roadmap per user decision (2026-05-11). Architectural placeholders only:
- The data-input module keeps `.h5ad` import working so the path for `anndata` + spatial metadata stays open
- New sidebar items can be added without core changes — when v2 starts, add a `modules/spatial/` directory mirroring `modules/de/` shape

---

# PHASE 6 — Testing + observability (continuous from Week 2)

## P6-1 · Test the 5 highest-stakes modules
Current test coverage is ~10%. Add tests in this order:

1. `modules/auth/auth_db.R` — login attempts, session token validation, project save/load roundtrip with versioning (after P1-1)
2. `modules/data_input.R` — fixture-based tests for 10x, h5ad, h5seurat, rds, csv
3. `modules/preprocessing.R` — integration with 2-sample mock
4. `modules/clustering.R` + new `modules/subcluster/` — resolution sweep, tree integrity
5. `modules/export.R` — full project save → load roundtrip equality

**Acceptance:** 60 %+ statement coverage on these 5 modules via `covr`.

**Effort:** L

## P6-2 · `shinytest2` E2E smoke tests
Add `shinytest2` integration tests that boot the app, log in as a test user, load the demo dataset, run a DE comparison, run enrichment, and screenshot key views. Run in CI.

**Effort:** M

## P6-3 · Structured logging + Sentry-style error tracking
Build on P0-4. Add Sentry (or `errlogger` self-hosted) for production error tracking.

**Effort:** S after P0-4

---

# PHASE 7 — Deployment + ops (Week 13+)

## P7-1 · Self-hosted Docker + nginx + HTTPS
- Docker compose with `shiny-server`, `nginx` (reverse proxy + TLS termination), `letsencrypt`
- Volume mounts for `data/`, `logs/`
- `shiny-server.conf` with `simple_scheduler 4`

**Effort:** M

## P7-2 · CI/CD pipeline to production
- Tag-based release: `git tag v1.0.0` → CI builds image → pushes to ghcr.io → server pulls
- Blue/green via two `docker-compose` stacks

**Effort:** M

## P7-3 · Backup + disaster recovery
- Nightly `sqlite_backup` of `data/scrna4_auth.db` → S3 with 30-day retention
- Weekly rsync of `data/user_sessions/` → S3
- Documented restore procedure

**Effort:** M

## P7-4 · Per-user resource quotas
- `users.quota_ram_gb`, `users.quota_disk_gb`, `users.quota_jobs_per_day`
- Enforce in `task_manager` before queuing
- Admin dashboard to adjust

**Effort:** M

## P7-5 · Telemetry (anonymized) + admin dashboard
- Per-module usage counts (anonymized)
- Average runtime per analysis
- Error rate
- Active users
- Admin-only sidebar item

**Effort:** L

---

# Implementation order — recommended sprints

| Sprint | Phase | Weeks | Outcome |
|---|---|---|---|
| Sprint 0 | P0 (all 5 items) | Week 1 | Deployable container; CI green; secure auth |
| Sprint 1 | P1-1 + P1-3 + P1-5 | Weeks 2–3 | Full state save; cancellable tasks; gallery downloads work |
| Sprint 2 | P1-2 + P1-4 + P2-2 | Week 4 | Resume banner; multi-project parallelism; memory pre-flight |
| Sprint 3 | P2-1 + P2-3 | Week 5 | BPCells on-disk; 8 paths async |
| Sprint 4 | P3 (all) | Weeks 6–7 | Unified AI + tool use + RAG (with model picker) |
| Sprint 5 | P4-1 + P4-3 | Weeks 8–9 | Subclustering + annotation rebuilds (incl. CellTypist) |
| Sprint 6 | P4-2 + P4-4 + P45 | Weeks 10–11 | Trajectory + CCC + **tier plans / quotas / 30-day retention / admin** |
| Sprint 7 | P5-1..P5-6 | Week 12 | Sci-quality features (batch QC, propeller, methods report, GEO/CellxGene import, h5ad export) |
| Sprint 8 | P6 + P7 (deployment) | Week 13 | E2E tests, Docker prod, backups, deploy decision (VPS vs GCP) |

After Sprint 8, the platform is **production**. Phase 4.5 (tier system) runs **parallel** in Sprint 6 because it's mostly auth/DB/admin work and doesn't conflict with the analysis rebuilds.

---

# Acceptance criteria summary

For each phase, the definition-of-done. We do not advance to the next phase until these pass.

| Phase | DoD |
|---|---|
| P0 | App boots from `docker run`; CI green; login lockout works; raw errors don't reach UI |
| P1 | Full project state survives logout/login including UI; all heavy ops cancellable; gallery CSVs downloadable |
| P2 | 200k-cell dataset processes on 16 GB RAM; 8 heavy paths run async with progress + cancel |
| P3 | Single AI assistant with model picker; calls tools; retrieves project context via RAG; persists conversation; per-call cost meter |
| P4 | Subclustering produces hierarchical tree with audit; trajectory split into modules with velocity working; annotation has SingleR with scoring + CellTypist + Azimuth; CCC <8 min on 27k cells |
| **P4.5** | **DB migration adds plan/quota columns; dataset-count + daily-job + AI + RAM gates enforced at 4 points; Account UI shows usage strip; 30-day inactivity cron archives projects; admin can set plans** |
| P5 | Methods report auto-generates Word doc; propeller p-values for differential abundance; cell-cycle scored; GEO/CellxGene import works; h5ad export works |
| P6 | 60%+ coverage on auth/data_input/preprocessing/clustering/export; shinytest2 E2E in CI |
| P7 | Self-hosted Docker stack with TLS + nightly backups; deployment chosen between 128 GB VPS or GCP; quota enforcement integrated with billing path |

---

# Risk register

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| BPCells API changes in Seurat v5 minor releases | M | H | Pin Seurat version via renv; abstract storage via `safe_get_assay_data()` |
| Tool-use loop infinite-runs the AI | M | M | Hard cap 10 tool calls per turn; cancel button on AI panel |
| Async `future_promise` workers lose state | L | H | Persist task queue to SQLite (P1-4) |
| Anthropic API rate limits during RAG indexing | M | L | Local embedding model (P3-3); throttle to 5 req/s |
| 100k+ cell integrations exceed disk budget | M | H | Quota enforcement P7-4; explicit user warning |
| Schema migration of saved projects (P1-1) breaks old data | M | H | Versioned migration functions with regression tests on stored fixtures |
| Scientific reproducibility regression | L | H | Audit log + regression tests on the mock dataset for every engine |

---

# Why this plan is "scientific + reproducible + robust"

- **Every engine writes an audit log** (already in place) with: package versions, seed, parameters, cell counts, timestamps, R version — extended in P1-1 to also capture UI state and dataset signatures
- **Schema versioning** on project bundles (P1-1) makes saved analyses portable across app upgrades
- **Pinned dependencies** via renv (P0-1) make any analysis re-runnable years later
- **CI regression guards** (P0-3) protect against silent semantic drift
- **Methods report generator** (P5-4) makes the science publishable directly from the app
- **Two-method consensus** as a UX pattern across annotation, integration QC, and trajectory — surfaces method disagreement instead of hiding it
- **Cancellable tasks** (P1-3) preserve work even on user abort — no "lost session" silent data loss
- **On-disk storage** (P2-1) makes large-dataset analyses tractable on modest academic hardware
- **Tool-using AI with provenance** (P3-2 + P3-4) means every AI claim is linked to the underlying analysis — no hallucinated conclusions

---

# Resolved questions (2026-05-11)

| # | Question | Decision | Status |
|---|---|---|---|
| 1 | Deployment | Docker-first; user will choose between 128 GB VPS or GCP later | ✅ Plan supports both |
| 2 | AI model | User-chosen at runtime (Opus/Sonnet/Haiku picker) | ✅ Built into P3-1 |
| 3 | Spatial transcriptomics | Deferred to v2 | ✅ Removed from scope |
| 4 | CellTypist | OK via reticulate | ✅ Included in P4-3 |
| 5 | Data retention | 1 month inactivity → archive | ✅ Built into P45-4 |
| 6 | Quotas | 4 tiers: Free=1 / Basic=5 / Standard=10 / Pro=50 datasets | ✅ Built into Phase 4.5 |

# Still open (low impact — can decide later)

| # | Question | Default if no answer |
|---|---|---|
| 7 | Concurrent-user budget | Defaults to `simple_scheduler 4` (good for 4–8 users on the 128 GB VPS); revisit when usage data exists |
| 8 | GDPR / IRB / encryption-at-rest | Default: SQLite + filesystem only on first launch; add at-rest encryption (P7-3) when first protected-data user signs up |
| 9 | Billing integration | Default: admin manually sets plans in admin dashboard (P45-6). Stripe added in Sprint 9+ when needed. |

---

**Total estimated effort:** ~13 weeks for one full-time engineer to reach the end of Sprint 8 (production).
**Minimum for "stop the bleeding":** Sprint 0 only (1 week) — gets deployable + secure baseline.
**Recommended order to commit:** Sprint 0 → Sprint 1 → Sprint 3 → Sprint 4 → Sprint 6 (the tier system is needed early because every other feature needs to respect quotas). Sprints 2, 5, 7, 8 fill in around them.

