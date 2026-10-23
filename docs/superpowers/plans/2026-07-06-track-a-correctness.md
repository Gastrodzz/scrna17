# Track A — Correctness & Robustness Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix five correctness/robustness issues in the scRNA-seq Shiny platform: auto-updating model list, species auto-detection (+ a silent mito bug), cell-cycle consolidation, removal of the redundant "Ask AI" buttons, and an edit path for project/AI settings (fixing the "stuck" bug).

**Architecture:** Introduce two new self-contained utility files (`ai_models.R`, `species_utils.R`) as single sources of truth, then rewire existing modules to consume them and remove dead/duplicated code. Pure logic is covered by `testthat`; Shiny wiring is verified by running the app locally.

**Tech Stack:** R, Shiny, shinydashboard, Seurat v5, httr2, DBI/RSQLite, testthat.

**Spec:** `docs/superpowers/specs/2026-07-06-track-a-correctness-design.md`

---

## Notes for the implementer

- **This project is NOT a git repo.** "Commit" steps are written as **checkpoints**: run the tests / launch the app and confirm the stated result before moving on. If you want history, run `git init` first and replace each checkpoint with a real commit.
- **Launch the app** with: `R -e "shiny::runApp('.', port=7788, launch.browser=FALSE)"` (or open `app.R` in RStudio and click *Run App*). The login page appears first; register a throwaway account to reach the dashboard.
- **Run unit tests** with: `R -e "testthat::test_dir('tests/testthat')"` (or `testthat::test_file('tests/testthat/test-ai-models.R')` for one file).
- `%||%` (null-coalesce) is defined in `app.R:115`, **after** modules are sourced, so new files must not rely on it at source time — each new file defines its own guard.

---

## Task 1: Create `modules/ai_models.R` (central model registry)

**Files:**
- Create: `modules/ai_models.R`
- Test: `tests/testthat/test-ai-models.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-ai-models.R`:

```r
source(file.path("..", "..", "modules", "ai_models.R"))

test_that("fallback list is returned when no api key", {
  ch <- list_claude_models(api_key = NULL, include_thinking = FALSE)
  expect_true(is.character(ch))
  expect_true(length(ch) >= 3)
  expect_true(all(grepl("^claude-", unname(ch))))
  expect_true(!is.null(names(ch)))          # display names present
})

test_that("thinking variants are appended for opus/sonnet only", {
  ch <- list_claude_models(api_key = NULL, include_thinking = TRUE, force = TRUE)
  think <- unname(ch)[grepl("-thinking$", unname(ch))]
  expect_true(length(think) >= 1)
  expect_true(all(grepl("opus|sonnet", think)))   # no haiku-thinking
  expect_false(any(grepl("haiku.*-thinking", think)))
})

test_that("get_default_model prefers a non-thinking sonnet", {
  m <- get_default_model(api_key = NULL)
  expect_true(is.character(m) && length(m) == 1)
  expect_true(grepl("sonnet", m))
  expect_false(grepl("-thinking$", m))
})

test_that("cache is reused within ttl and bypassed with force", {
  a <- list_claude_models(api_key = NULL, force = TRUE)
  b <- list_claude_models(api_key = NULL, force = FALSE)
  expect_identical(a, b)
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `R -e "testthat::test_file('tests/testthat/test-ai-models.R')"`
Expected: FAIL — `could not find function "list_claude_models"`.

- [ ] **Step 3: Create `modules/ai_models.R`**

```r
# modules/ai_models.R
# Central registry for Claude models. Live-fetches from Anthropic's /v1/models
# endpoint when an API key is available, caches the result, and falls back to a
# curated static list otherwise. Single source of truth for every model dropdown
# and default-model resolution in the app.

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# --- Curated fallback (newest -> oldest); used when no key / offline / API error
AI_MODELS_FALLBACK <- c(
  "Claude Opus 4.6 (Most capable)"  = "claude-opus-4-6",
  "Claude Sonnet 4.6 (Recommended)" = "claude-sonnet-4-6",
  "Claude Sonnet 4.5"               = "claude-sonnet-4-5-20250929",
  "Claude Haiku 4.5 (Fastest)"      = "claude-haiku-4-5-20251001"
)
AI_MODELS_DEFAULT <- "claude-sonnet-4-6"

# Families that support extended thinking (used to synthesize "-thinking" variants)
.AI_THINKING_FAMILIES <- c("opus", "sonnet")

# Module-level cache
.ai_models_cache <- new.env(parent = emptyenv())
.ai_models_cache$data <- NULL
.ai_models_cache$ts   <- 0

#' List available Claude models as a named vector (display_name -> id).
#' @param api_key Anthropic key; when supplied, live-fetches from /v1/models.
#' @param force Bypass the cache.
#' @param ttl Cache lifetime in seconds.
#' @param include_thinking Append synthesized "-thinking" display variants.
list_claude_models <- function(api_key = NULL, force = FALSE, ttl = 3600,
                               include_thinking = TRUE) {
  now <- as.numeric(Sys.time())
  if (!force && !is.null(.ai_models_cache$data) &&
      (now - .ai_models_cache$ts) < ttl) {
    return(.ai_models_cache$data)
  }

  live <- NULL
  if (!is.null(api_key) && nzchar(trimws(api_key)) &&
      requireNamespace("httr2", quietly = TRUE)) {
    live <- tryCatch({
      resp <- httr2::request("https://api.anthropic.com/v1/models") |>
        httr2::req_headers(
          "x-api-key"         = trimws(api_key),
          "anthropic-version" = "2023-06-01"
        ) |>
        httr2::req_url_query(limit = 100) |>
        httr2::req_timeout(20) |>
        httr2::req_perform()
      body <- httr2::resp_body_json(resp)
      data <- body$data
      if (is.null(data) || length(data) == 0) stop("empty model list")
      ids  <- vapply(data, function(m) m$id %||% NA_character_, character(1))
      disp <- vapply(data, function(m) m$display_name %||% m$id %||% NA_character_, character(1))
      keep <- !is.na(ids) & nzchar(ids) & grepl("^claude-", ids)
      ids <- ids[keep]; disp <- disp[keep]
      if (length(ids) == 0) stop("no claude models returned")
      stats::setNames(ids, disp)
    }, error = function(e) NULL)
  }

  base <- if (!is.null(live)) live else AI_MODELS_FALLBACK

  if (isTRUE(include_thinking)) {
    thinking <- character(0)
    for (i in seq_along(base)) {
      id <- base[[i]]; nm <- names(base)[i]
      if (any(vapply(.AI_THINKING_FAMILIES,
                     function(f) grepl(f, id, ignore.case = TRUE), logical(1)))) {
        thinking[paste0(nm, " — Extended Thinking")] <- paste0(id, "-thinking")
      }
    }
    base <- c(base, thinking)
  }

  # Cache live results for the full ttl; cache fallback only briefly so that a
  # later key entry re-triggers a live fetch soon.
  .ai_models_cache$data <- base
  .ai_models_cache$ts   <- if (!is.null(live)) now else (now - ttl + 60)
  base
}

#' Named choices vector for selectInput (display_name -> id).
ai_model_choices <- function(api_key = NULL, force = FALSE) {
  list_claude_models(api_key = api_key, force = force)
}

#' Resolve a sensible default model id (non-thinking sonnet preferred).
get_default_model <- function(api_key = NULL) {
  ch <- tryCatch(list_claude_models(api_key = api_key, include_thinking = FALSE),
                 error = function(e) AI_MODELS_FALLBACK)
  ids <- unname(ch)
  sonnet <- ids[grepl("sonnet", ids, ignore.case = TRUE) & !grepl("-thinking$", ids)]
  if (length(sonnet) > 0) return(sonnet[1])
  if (length(ids) > 0) return(ids[1])
  AI_MODELS_DEFAULT
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `R -e "testthat::test_file('tests/testthat/test-ai-models.R')"`
Expected: PASS (4 tests).

- [ ] **Step 5: Checkpoint** — tests green. Proceed.

---

## Task 2: Wire `ai_models.R` into the app and remove hardcoded model lists

**Files:**
- Modify: `app.R:71-92` (source order), `app.R:1177`, `app.R:1122-1153` (`.reset_values`)
- Modify: `modules/ai_assistant.R:15-30`, `:341-347`, `:349-378`, `:722-731`
- Modify: `modules/auth/project_module.R:212-221`
- Modify: `modules/annotation.R:2807`, `modules/downstream/downstream_ai.R` (4 defaults)

- [ ] **Step 1: Source `ai_models.R` before `ai_assistant.R`**

In `app.R`, inside the `modules_loaded <- list(...)` block, add a line as the FIRST entry (before `core_robustness`):

```r
modules_loaded <- list(
  # Central registries first (single source of truth, no deps)
  ai_models = safe_source("modules/ai_models.R"),
  # Utility modules
  core_robustness = safe_source("modules/core_robustness.R"),
  ai_assistant = safe_source("modules/ai_assistant.R"),
  ...
```

- [ ] **Step 2: Replace the hardcoded `AI_MODEL_CONFIG` list**

In `modules/ai_assistant.R`, replace lines 15-30 (`AI_MODEL_CONFIG <- list(...)`) with:

```r
AI_MODEL_CONFIG <- list(
  # available_models/default_model now resolve dynamically via ai_models.R.
  # Kept for backward-compat: these accessors read the registry.
  default_model    = tryCatch(get_default_model(), error = function(e) "claude-sonnet-4-6"),
  available_models = tryCatch(as.list(ai_model_choices()), error = function(e) list()),
  api_timeout      = 240,
  context_budget   = 6000
)
```

- [ ] **Step 3: Generalize `get_maxtokens` model checks**

In `modules/ai_assistant.R`, replace lines 375-377:

```r
  if (grepl("opus", model, ignore.case = TRUE))        pmin(as.integer(base * 2L), 32000L)
  else if (grepl("haiku", model, ignore.case = TRUE))  pmin(base, 8192L)
  else                                                 pmin(base, 16384L)
```

- [ ] **Step 4: Use registry default in `get_model_for_task`**

In `modules/ai_assistant.R:343`, replace:

```r
  if (is.null(user_model) || !nzchar(user_model)) user_model <- get_default_model(values$ai_api_key)
```

- [ ] **Step 5: Make `aiConfigUI` model dropdown dynamic**

In `modules/ai_assistant.R`, replace the `selectInput(ns("model"), ...)` block (lines 723-731) with a placeholder that the server fills:

```r
        selectInput(ns("model"), "Default Model", width = "100%",
          choices = ai_model_choices(), selected = get_default_model()),
```

- [ ] **Step 6: Make the project-create dropdown dynamic**

In `modules/auth/project_module.R`, replace the `selectInput(ns("new_ai_model"), ...)` block (lines 214-220) with:

```r
                  selectInput(ns("new_ai_model"), label = NULL,
                              choices = ai_model_choices(),
                              selected = get_default_model(),
                              selectize = FALSE, width = "100%")
```

Then, still in `project_module.R` server, add an observer (place it right after the `observeEvent(input$btn_gen_summary, {...})` block, before `observeEvent(input$confirm_create, ...)`) so pasting a key refreshes the live list:

```r
    # Refresh model list live once a key is entered (debounced)
    observeEvent(input$new_ai_api_key, {
      key <- trimws(input$new_ai_api_key %||% "")
      if (nzchar(key)) {
        ch <- tryCatch(ai_model_choices(api_key = key, force = TRUE),
                       error = function(e) NULL)
        if (!is.null(ch)) {
          sel <- isolate(input$new_ai_model) %||% get_default_model(key)
          updateSelectInput(session, "new_ai_model", choices = ch,
                            selected = if (sel %in% ch) sel else get_default_model(key))
        }
      }
    }, ignoreInit = TRUE)
```

- [ ] **Step 7: Replace remaining hardcoded defaults**

Replace each of these literal defaults with `get_default_model()`:
- `app.R:1177` → `ai_model = get_default_model(),`
- `app.R:943` (`values$ai_model <- auth_state$... %||% "claude-sonnet-4-6"`) → `%||% get_default_model()`
- `app.R` `.reset_values` (add near the ai_cache reset, ~line 1152): `values$ai_model <- get_default_model()`
- `modules/annotation.R:2807` → `model <- values$ai_model %||% get_default_model(values$ai_api_key)`
- `modules/auth/project_module.R:386` and `:421` (`%||% "claude-sonnet-4-6"`) → `%||% get_default_model()`
- `modules/downstream/downstream_ai.R:999,1513,2282,2426` → replace `%||% "claude-sonnet-4-5-20250929"` with `%||% get_default_model(values$ai_api_key)`

Find them all with:
```bash
grep -rn '"claude-sonnet-4-6"\|"claude-sonnet-4-5-20250929"\|"claude-opus-4-6"\|"claude-haiku-4-5-20251001"' --include="*.R" app.R modules/ | grep -v "modules/ai_models.R"
```
Expected after edits: only `modules/ai_models.R` (the fallback list) should still contain literal model IDs.

- [ ] **Step 8: Checkpoint** — launch the app, create a project WITHOUT a key: the model dropdown shows the fallback list. Paste a valid key: the dropdown repopulates from the live `/v1/models` response (a newly released model appears without code changes). Confirm no errors in the R console.

---

## Task 3: Create `modules/species_utils.R` (shared species detection)

**Files:**
- Create: `modules/species_utils.R`
- Test: `tests/testthat/test-species-utils.R`

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-species-utils.R`:

```r
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
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `R -e "testthat::test_file('tests/testthat/test-species-utils.R')"`
Expected: FAIL — `could not find function "detect_species"`.

- [ ] **Step 3: Create `modules/species_utils.R`**

```r
# modules/species_utils.R
# Shared species detection + gene-pattern resolution. One source of truth for
# mitochondrial/ribosomal regexes and species labels, used by QC, preprocessing,
# and annotation so patterns never drift between modules.

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

.SPECIES_MARKERS <- list(
  human = c("PTPRC","EPCAM","PECAM1","CD19","CD3D","CD4","CD8A","FOXP3",
            "COL1A1","KRT19","MS4A1","GATA3","NKG7","S100A4"),
  mouse = c("Ptprc","Epcam","Pecam1","Cd19","Cd3d","Cd4","Cd8a","Foxp3",
            "Col1a1","Krt19","Ms4a1","Gata3","Nkg7","S100a4")
)

.SPECIES_PATTERNS <- list(
  human     = list(mt = "^MT-",  ribo = "^RP[SL]\\d"),
  mouse     = list(mt = "^mt-",  ribo = "^Rp[sl]\\d"),
  rat       = list(mt = "^Mt-",  ribo = "^Rp[sl]\\d"),
  zebrafish = list(mt = "^mt-",  ribo = "^rp[sl]\\d")
)

#' Map a free-text / latin species label to the app's canonical key.
normalize_species <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return("")
  s <- tolower(trimws(x))
  if (grepl("sapiens|human|homo", s))        return("human")
  if (grepl("musculus|mouse|\\bmus\\b", s))  return("mouse")
  if (grepl("norvegicus|\\brat\\b|rattus", s)) return("rat")
  if (grepl("rerio|zebrafish|danio", s))     return("zebrafish")
  if (s %in% c("human","mouse","rat","zebrafish")) return(s)
  ""
}

#' Detect species from gene identifiers.
#' @return list(species, mt_pattern, ribo_pattern, confidence, note)
detect_species <- function(genes) {
  genes <- as.character(genes)
  if (length(genes) == 0) {
    return(list(species = "unknown", mt_pattern = "^MT-", ribo_pattern = "^RP[SL]\\d",
                confidence = "none", note = "No genes supplied."))
  }
  if (mean(grepl("^ENS[A-Z]*G\\d", genes)) > 0.5) {
    return(list(species = "unknown", mt_pattern = "^MT-", ribo_pattern = "^RP[SL]\\d",
                confidence = "low",
                note = "Matrix uses Ensembl gene IDs; convert to symbols for reliable MT/ribo detection."))
  }
  n_mt_h <- sum(grepl("^MT-", genes)); n_mt_m <- sum(grepl("^mt-", genes))
  n_h <- sum(genes %in% .SPECIES_MARKERS$human)
  n_m <- sum(genes %in% .SPECIES_MARKERS$mouse)
  h_score <- n_mt_h * 2 + n_h
  m_score <- n_mt_m * 2 + n_m
  sp <- if (h_score == 0 && m_score == 0) "unknown"
        else if (h_score >= m_score) "human" else "mouse"
  conf <- if (sp == "unknown") "none"
          else if (abs(h_score - m_score) >= 6) "high"
          else if (abs(h_score - m_score) >= 2) "medium" else "low"
  pat <- .SPECIES_PATTERNS[[sp]] %||% .SPECIES_PATTERNS$human
  list(species = sp, mt_pattern = pat$mt, ribo_pattern = pat$ribo,
       confidence = conf, note = sprintf("human_score=%d, mouse_score=%d", h_score, m_score))
}

#' MT/ribo patterns for a known species key (free-text accepted).
species_patterns <- function(species_key) {
  .SPECIES_PATTERNS[[normalize_species(species_key)]] %||% .SPECIES_PATTERNS$human
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `R -e "testthat::test_file('tests/testthat/test-species-utils.R')"`
Expected: PASS (4 tests).

- [ ] **Step 5: Source it in `app.R`** — add to `modules_loaded` right after `ai_models`:

```r
  species_utils = safe_source("modules/species_utils.R"),
```

- [ ] **Step 6: Checkpoint** — tests green, app still launches.

---

## Task 4: Fix the silent mito bug + wire species auto-detection

**Files:**
- Modify: `app.R:947-966` (the auto-populate loop), add a data-load observer
- Modify: `modules/quality_control.R:200-217` (MT compute)
- Modify: `modules/annotation.R:78-97` (remove the local `.detect_species`, defer to shared)

- [ ] **Step 1: Stop pushing latin species into the `human/mouse` dropdown**

In `app.R`, the loop at lines 951-960 pushes raw `sp` (e.g. `"Homo sapiens"`) into `qcMetrics-species_qc` whose choices are `"human"/"mouse"`. Change the species-select loop to normalize the value for the QC/species dropdowns. Replace the block:

```r
    sp <- values$species_global
    ts <- values$tissue_global
    sp_norm <- normalize_species(sp)   # "" when unknown/other
    # species selects (choices use human/mouse vocabulary)
    for (full_id in c("dataUpload-species", "qcMetrics-species_qc",
                      "cellCycle-species_cc", "cellCycle-cc_species",
                      "aiAnn-ai_species", "preproc-cc_species",
                      "enrichment-species", "geneset-gs_species",
                      "ccc-cc_species", "tf-tf_species",
                      "lr-lr_species", "autoAnn-ai_auto_species")) {
      val <- if (nzchar(sp_norm)) sp_norm else sp
      tryCatch(updateSelectInput(session, full_id, selected = val),
               error = function(e) NULL)
    }
```

> Note: if any of those selects use latin choices rather than `human/mouse`, `updateSelectInput` with a non-matching value is a no-op (safe). This change fixes `species_qc` specifically without breaking others.

- [ ] **Step 2: Make QC MT computation self-correcting**

In `modules/quality_control.R`, replace lines 200-217 (the `mt_pattern <- switch(...)` through the `else { ... showNotification("No mitochondrial genes found with pattern") }` block) with:

```r
          # Resolve MT pattern from the dropdown, then verify against the data.
          mt_pattern <- switch(input$species_qc,
            "human"  = "^MT-",
            "mouse"  = "^mt-",
            "custom" = input$custom_mt_pattern,
            "^MT-"
          )
          mt_genes <- grep(mt_pattern, rownames(obj), value = TRUE)

          # Self-correct: if the chosen pattern finds nothing, auto-detect from data
          if (length(mt_genes) == 0) {
            det <- detect_species(rownames(obj))
            if (det$species != "unknown" && det$mt_pattern != mt_pattern) {
              mt_pattern <- det$mt_pattern
              mt_genes   <- grep(mt_pattern, rownames(obj), value = TRUE)
              if (length(mt_genes) > 0) {
                showNotification(sprintf(
                  "Species dropdown found no MT genes; auto-detected %s (pattern '%s', %d genes).",
                  det$species, mt_pattern, length(mt_genes)),
                  type = "warning", duration = 8)
              }
            }
          }

          if (length(mt_genes) > 0) {
            obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
            message("Found ", length(mt_genes), " mitochondrial genes")
          } else {
            tmp <- rep(0, ncol(obj)); names(tmp) <- colnames(obj)
            obj[["percent.mt"]] <- tmp
            showNotification(
              "No mitochondrial genes found. If this is expected (pre-filtered matrix), ignore; otherwise check gene naming (symbols vs Ensembl IDs).",
              type = "warning", duration = 10)
          }
```

- [ ] **Step 3: Auto-detect on data load (UX hint)**

In `app.R`, add a new observer just after the `observeEvent(auth_state$current_project_id, ...)` block (after line 1006). This sets the QC species dropdown from the actual matrix when a dataset becomes active:

```r
  # Auto-detect species from the active dataset and pre-set QC controls
  observeEvent(values$active_dataset, {
    req(values$active_dataset)
    det <- tryCatch(detect_species(rownames(values$active_dataset)),
                    error = function(e) NULL)
    if (is.null(det) || det$species == "unknown") return()
    sel <- if (det$species %in% c("human","mouse")) det$species else "custom"
    tryCatch(updateSelectInput(session, "qcMetrics-species_qc", selected = sel),
             error = function(e) NULL)
    if (sel == "custom")
      tryCatch(updateTextInput(session, "qcMetrics-custom_mt_pattern",
                               value = det$mt_pattern), error = function(e) NULL)
    showNotification(sprintf("Auto-detected species: %s (%s confidence).",
                             det$species, det$confidence),
                     type = "message", duration = 5)
  }, ignoreInit = TRUE)
```

- [ ] **Step 4: De-duplicate `.detect_species` in annotation**

In `modules/annotation.R`, delete the local `.detect_species <- function(genes) {...}` (lines 78-97) and add, at the top of the file (after any `library()`/`source()` lines), a thin alias so existing call sites keep working:

```r
# Species detection now lives in modules/species_utils.R (single source of truth)
.detect_species <- function(genes) detect_species(genes)$species
```

> `species_utils.R` is sourced before `annotation.R` in `app.R`, so `detect_species` is available.

- [ ] **Step 5: Checkpoint** — launch the app, load a **mouse** dataset (or any matrix with `mt-` genes). Go to QC → Compute metrics. Confirm: (a) a "Auto-detected species: mouse" notification appears, (b) `percent.mt` is non-zero (check the QC violin/scatter), (c) filtering on MT works. Re-run tests: `R -e "testthat::test_dir('tests/testthat')"` — all green.

---

## Task 5: Remove the redundant "Ask AI" buttons

**Files:**
- Modify: `modules/visualization.R`, `modules/trajectory.R`, `modules/export.R`,
  `modules/preprocessing.R`, `modules/quality_control.R`,
  `modules/dimensionality_reduction.R`, `modules/clustering.R`
- Modify: `modules/ai_assistant.R` (delete `ai_helper_button_ui` + `ai_helper_server` defs, ~1276-1351)

- [ ] **Step 1: Locate every usage**

Run:
```bash
grep -rn "ai_helper_button_ui\|ai_helper_server\|ai_helper_panel\|ai_helper_question\|ai_helper_ask\|ai_helper_thread" --include="*.R" modules/ | grep -v "modules/ai_assistant.R"
```
This lists the UI insertions and server registrations to remove (visualization:86/571, trajectory:137/1030, export:273/492, preprocessing:365, quality_control:441, dimensionality_reduction:130/594, clustering:855/1786).

- [ ] **Step 2: Remove each UI button call**

For each `ai_helper_button_ui(ns, "<Module>")` occurrence, delete the whole call. Where it sits inside a `tagList(...)`/comma-separated UI list (e.g. `visualization.R:86`, `clustering.R:855`), also remove the trailing/leading comma so the surrounding `tagList`/`div` stays valid. Where it is guarded by `if (exists("ai_helper_button_ui", ...))`, remove the entire `if (...)` expression. Where it is UNguarded (`trajectory.R:137`, `dimensionality_reduction.R:130`), remove that line.

- [ ] **Step 3: Remove each server registration**

For each `ai_helper_server("<id>", values, "<Module>", context_fn = function() {...})` occurrence — including any wrapping `if (exists("ai_helper_server", mode = "function")) { ... }` — delete the entire block. These end at the closing `})` of the `ai_helper_server(...)` call.

- [ ] **Step 4: Delete the now-dead definitions**

In `modules/ai_assistant.R`, delete the section `# Cross-Module AI Helper Button` including both `ai_helper_button_ui <- function(...) {...}` and `ai_helper_server <- function(...) {...}` (approximately lines 1273-1351; confirm the exact end is the closing `}` of `ai_helper_server`). Keep everything after it intact.

- [ ] **Step 5: Verify nothing references the removed helpers**

Run:
```bash
grep -rn "ai_helper" --include="*.R" modules/ app.R
```
Expected: **no output**.

- [ ] **Step 6: Checkpoint** — launch the app; open Visualization, Trajectory, Export, Preprocessing, QC, Dim-Reduction, Clustering tabs one by one. Each renders with **no "Ask AI" button** and **no errors** in the console. The dedicated AI features (annotation, QC threshold advisor, DE advisor, downstream interpretation) still work.

---

## Task 6: Consolidate cell cycle + add it to the QC report

**Files:**
- Investigate/Modify: `modules/quality_control.R` (Cell Cycle sub-tab + `qcReport`)
- Investigate: `modules/preprocessing.R` (cell-cycle scoring is the canonical compute)
- Modify: `app.R:1683` (QC module registrations) if the QC-Metrics area still computes CC

- [ ] **Step 1: Map the current cell-cycle surface**

Run:
```bash
grep -n "CellCycleScoring\|S.Score\|G2M.Score\|Phase\|cell_cycle\|cellCycle\|cc_" modules/quality_control.R
grep -n "CellCycleScoring\|S.Score\|G2M.Score\|Phase\|cc_" modules/preprocessing.R
sed -n '460,520p' app.R   # QC tabItem UI — find the Cell Cycle sub-tabPanel
```
Decide: (a) confirm the **QC Metrics** compute path does NOT run `CellCycleScoring` (the comment at `quality_control.R:234` says it moved); (b) confirm the dedicated **Cell Cycle** sub-tab exists; (c) confirm Preprocessing computes `S.Score/G2M.Score/Phase`.

- [ ] **Step 2: Ensure scoring is Preprocessing-only; Cell Cycle tab consumes results**

If the dedicated Cell Cycle sub-tab (`cellCycleServer`) recomputes scoring on raw data, change it to: if `all(c("S.Score","G2M.Score","Phase") %in% colnames(obj@meta.data))`, visualize those; otherwise show an actionable message instead of erroring:

```r
      if (!all(c("S.Score","G2M.Score") %in% colnames(obj@meta.data))) {
        showNotification(
          "Cell-cycle scores not found. Run 'Cell Cycle Scoring' in the Preprocessing tab (requires normalized data), then return here.",
          type = "warning", duration = 8)
        return()
      }
```

- [ ] **Step 3: Add a Cell-Cycle section to the QC report**

Locate the `qcReport` render (search `output$` in the `qcReportServer` function). Add a Cell-Cycle block to the report body. Insert this helper near the top of `qcReportServer` and call it where the report HTML/markdown is assembled:

```r
    cc_report_section <- function(obj) {
      md <- obj@meta.data
      if (!all(c("S.Score","G2M.Score","Phase") %in% colnames(md))) {
        return(tags$div(class = "qc-report-section",
          tags$h4("Cell Cycle"),
          tags$p(style = "color:#888;",
                 "Not scored. Run Cell Cycle Scoring in Preprocessing to include this in the report.")))
      }
      ph <- table(md$Phase)
      tags$div(class = "qc-report-section",
        tags$h4("Cell Cycle"),
        tags$ul(
          tags$li(sprintf("G1: %d cells (%.1f%%)",  ph["G1"]  %||% 0L, 100 * (ph["G1"]  %||% 0L)/nrow(md))),
          tags$li(sprintf("S: %d cells (%.1f%%)",   ph["S"]   %||% 0L, 100 * (ph["S"]   %||% 0L)/nrow(md))),
          tags$li(sprintf("G2M: %d cells (%.1f%%)", ph["G2M"] %||% 0L, 100 * (ph["G2M"] %||% 0L)/nrow(md)))
        ),
        tags$p(sprintf("Mean S.Score = %.3f, Mean G2M.Score = %.3f",
                       mean(md$S.Score), mean(md$G2M.Score)))
      )
    }
```

> `%||%` is available at runtime (defined in `app.R:115`). Wire `cc_report_section(obj)` into the report's `tagList`/UI assembly next to the doublet/ambient summaries.

- [ ] **Step 4: Checkpoint** — launch the app; run normalization + Cell Cycle Scoring in Preprocessing; open the dedicated Cell Cycle tab (renders from stored scores, no error); open the QC report and confirm the Cell Cycle section shows phase counts. Try the Cell Cycle tab BEFORE scoring — confirm it shows the guidance message, not an error.

---

## Task 7: Dashboard "Edit Project" modal

**Files:**
- Modify: `modules/auth/project_module.R` (UI: add edit modal + pencil button; server: prefill + save)

- [ ] **Step 1: Add a pencil (edit) button to each project card**

In `projectDashboardServer`'s `output$projects_grid` card actions (the `tags$div(class = "proj-card-actions", ...)` block around lines 543-566), add an Edit button before the Delete button:

```r
              tags$button(
                class   = "btn-proj-edit",
                type    = "button",
                title   = "Edit project",
                onclick = sprintf(
                  "Shiny.setInputValue('%s', %d, {priority:'event'})",
                  ns("edit_project"), p$id
                ),
                tags$i(class = "fa fa-pen")
              ),
```

- [ ] **Step 2: Add the Edit modal markup**

In `projectDashboardUI`, after the Create modal `tags$div(id = ns("modal_create"), ...)` block (before the inline-JS `tags$script`), add an Edit modal. It reuses the create fields with `edit_` ids:

```r
    # ── Edit Project modal ───────────────────────────────────────────────────
    tags$div(id = ns("modal_edit"), class = "proj-modal-overlay proj-modal--hidden",
      tags$div(class = "proj-modal proj-modal--wide",
        tags$div(class = "proj-modal-header",
          tags$div(class = "proj-modal-title",
            tags$div(class = "proj-modal-title-icon", tags$i(class = "fa fa-pen")),
            tags$div(tags$h4("Edit project"),
                     tags$p(class = "proj-modal-subtitle",
                            "Update any field, including the AI API key and model."))),
          tags$button(class = "proj-modal-close", type = "button",
                      onclick = sprintf("projHideModal('%s')", ns("modal_edit")),
                      `aria-label` = "Close", HTML("&times;"))
        ),
        tags$div(class = "proj-modal-body",
          tags$div(class = "proj-form-field",
            tags$label("Project name"),
            textInput(ns("edit_name"), label = NULL, width = "100%")),
          tags$div(class = "proj-form-field",
            tags$label("Description"),
            textAreaInput(ns("edit_desc"), label = NULL, rows = 2, width = "100%")),
          tags$div(class = "proj-grid-2",
            tags$div(class = "proj-form-field",
              tags$label("Species"),
              selectInput(ns("edit_species"), label = NULL,
                choices = list("— Not specified —" = "",
                               "Homo sapiens" = "Homo sapiens",
                               "Mus musculus" = "Mus musculus",
                               "Rattus norvegicus" = "Rattus norvegicus",
                               "Danio rerio" = "Danio rerio",
                               "Other" = "Other"),
                selectize = FALSE, width = "100%")),
            tags$div(class = "proj-form-field",
              tags$label("Tissue / organ"),
              textInput(ns("edit_tissue"), label = NULL, width = "100%"))),
          tags$div(class = "proj-form-field",
            tags$label("Experimental conditions"),
            textInput(ns("edit_conditions"), label = NULL, width = "100%")),
          tags$hr(),
          tags$div(class = "proj-grid-2",
            tags$div(class = "proj-form-field",
              tags$label("Claude API key"),
              passwordInput(ns("edit_ai_api_key"), label = NULL,
                            placeholder = "sk-ant-… (leave blank to keep current)",
                            width = "100%")),
            tags$div(class = "proj-form-field",
              tags$label("Model"),
              selectInput(ns("edit_ai_model"), label = NULL,
                          choices = ai_model_choices(), selectize = FALSE, width = "100%"))),
          uiOutput(ns("edit_alert"))
        ),
        tags$div(class = "proj-modal-footer",
          tags$div(class = "proj-modal-footer-actions",
            tags$button(type = "button", class = "btn-proj-cancel",
                        onclick = sprintf("projHideModal('%s')", ns("modal_edit")),
                        "Cancel"),
            actionButton(ns("confirm_edit"),
                         tagList(tags$i(class = "fa fa-check"), " Save changes"),
                         class = "btn-proj-confirm")))
      )
    ),
```

- [ ] **Step 3: Prefill + open the modal**

In `projectDashboardServer`, add (near the other observers):

```r
    edit_msg <- reactiveVal(NULL)
    output$edit_alert <- renderUI({
      m <- edit_msg(); if (is.null(m)) return(NULL)
      tags$div(class = paste0("auth-alert auth-alert--", m$type), m$text)
    })
    editing_id <- reactiveVal(NULL)

    observeEvent(input$edit_project, {
      pid  <- as.integer(input$edit_project)
      proj <- get_project(pid, auth_state$user_id)
      if (nrow(proj) == 0) { showNotification("Project not found.", type = "error"); return() }
      editing_id(pid)
      updateTextInput(session, "edit_name",  value = proj$name[1])
      updateTextAreaInput(session, "edit_desc", value = proj$description[1] %||% "")
      updateSelectInput(session, "edit_species", selected = proj$species[1] %||% "")
      updateTextInput(session, "edit_tissue", value = proj$tissue[1] %||% "")
      updateTextInput(session, "edit_conditions", value = proj$conditions[1] %||% "")
      # refresh model list against the saved key, preselect saved model
      key <- proj$ai_api_key[1] %||% ""
      ch  <- tryCatch(ai_model_choices(api_key = if (nzchar(key)) key else NULL),
                      error = function(e) ai_model_choices())
      sel <- proj$ai_model[1] %||% get_default_model()
      updateSelectInput(session, "edit_ai_model", choices = ch,
                        selected = if (sel %in% ch) sel else get_default_model())
      updateTextInput(session, "edit_ai_api_key", value = "")  # never echo the key
      edit_msg(NULL)
      shinyjs::runjs(sprintf("projShowModal('%s')", ns("modal_edit")))
    })
```

- [ ] **Step 4: Save changes (and sync if currently open)**

```r
    observeEvent(input$confirm_edit, {
      pid <- editing_id(); req(pid)
      name <- trimws(input$edit_name %||% "")
      if (!nzchar(name)) { edit_msg(list(type = "error", text = "Project name is required.")); return() }

      fields <- list(
        name        = name,
        description = input$edit_desc %||% "",
        species     = input$edit_species %||% "",
        tissue      = input$edit_tissue %||% "",
        conditions  = input$edit_conditions %||% "",
        ai_model    = input$edit_ai_model %||% get_default_model()
      )
      # Only overwrite the key when a new one is typed
      new_key <- trimws(input$edit_ai_api_key %||% "")
      if (nzchar(new_key)) fields$ai_api_key <- new_key

      ok <- do.call(update_project_metadata,
                    c(list(project_id = pid, user_id = auth_state$user_id), fields))
      if (isTRUE(ok)) {
        shinyjs::runjs(sprintf("projHideModal('%s')", ns("modal_edit")))
        refresh_trigger(refresh_trigger() + 1)
        # If editing the currently-open project, sync live state
        if (identical(auth_state$current_project_id, pid)) {
          auth_state$current_project_name    <- name
          auth_state$current_project_species <- fields$species
          auth_state$current_project_tissue  <- fields$tissue
          auth_state$current_project_ai_model <- fields$ai_model
          if (!is.null(fields$ai_api_key))
            auth_state$current_project_ai_api_key <- fields$ai_api_key
        }
        showNotification("Project updated.", type = "message", duration = 4)
      } else {
        edit_msg(list(type = "error", text = "Update failed. Check the logs."))
      }
    })
```

- [ ] **Step 5: Add minimal CSS for the edit button**

In the project stylesheet (find with `grep -rln "btn-proj-delete" www/`), add a `.btn-proj-edit` rule mirroring `.btn-proj-delete` (same size, neutral/blue hover). If no dedicated file, add an inline `tags$style` in `projectDashboardUI`:

```r
    tags$style(HTML(".btn-proj-edit{background:transparent;border:1px solid #d0d7de;color:#57606a;border-radius:6px;padding:6px 10px;cursor:pointer;margin-right:6px}.btn-proj-edit:hover{background:#0d6efd;color:#fff;border-color:#0d6efd}")),
```

- [ ] **Step 6: Checkpoint** — launch app; create a project with NO key; click the pencil → modal prefilled; enter a key + pick a model → Save; reopen the modal → confirms model saved and key field blank (not echoed). Open the project and confirm an AI feature now works. **This closes the "forgot the key" trap from the dashboard side.**

---

## Task 8: In-app "AI Settings" panel (re-enable `aiConfigServer`)

**Files:**
- Modify: `app.R:236-320` (sidebar menu — add item), `app.R` tabItems (add tab body), `app.R:1661` (remove `FALSE &&`)
- Modify: `modules/ai_assistant.R` `aiConfigServer` (persist to DB; hydrate from `values`)

- [ ] **Step 1: Add the sidebar menu item**

In `app.R` `sidebarMenu(...)` (after the `Export & Report` `menuItem`, ~line 317), add:

```r
      menuItem("AI Settings",
        tabName = "ai_settings",
        icon = icon("robot"))
```

- [ ] **Step 2: Add the tab body**

In the `tabItems(...)` block (after the `export` tabItem, ~line 805), add:

```r
      tabItem(tabName = "ai_settings",
        fluidRow(column(12, h3("AI Assistant Settings"),
          p("Set or change the Claude API key and model for this project. Changes are saved to the project."),
          aiConfigUI("aiConfig")))
      ),
```

- [ ] **Step 3: Re-enable the server (remove the dead `FALSE &&`)**

In `app.R:1661`, replace:

```r
  if (isTRUE(modules_loaded$ai_assistant) && exists("aiConfigServer", mode = "function")) {
    aiConfigServer("aiConfig", values, auth_state)
  }
```

- [ ] **Step 4: Hydrate the panel from `values` and persist changes to the DB**

In `modules/ai_assistant.R`, change the `aiConfigServer` signature and add hydration + DB persistence. Replace `aiConfigServer <- function(id, values) {` with:

```r
aiConfigServer <- function(id, values, auth_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Hydrate controls from current project state when the panel appears
    observe({
      updateTextInput(session, "api_key", value = values$ai_api_key %||% "")
    })
    observeEvent(TRUE, {
      ch  <- tryCatch(ai_model_choices(api_key = values$ai_api_key),
                      error = function(e) ai_model_choices())
      sel <- values$ai_model %||% get_default_model(values$ai_api_key)
      updateSelectInput(session, "model", choices = ch,
                        selected = if (sel %in% ch) sel else get_default_model())
    }, once = TRUE)
```

Then, inside the existing `observeEvent(input$api_key, {...})` and `observeEvent(input$model, {...})` handlers, ALSO persist to the DB when a project is open. Append to the `input$api_key` observer body:

```r
      if (!is.null(auth_state) && !is.null(auth_state$current_project_id) &&
          nzchar(trimws(input$api_key %||% ""))) {
        update_project_metadata(auth_state$current_project_id, auth_state$user_id,
                                ai_api_key = trimws(input$api_key))
        # refresh the live model list against the new key
        ch <- tryCatch(ai_model_choices(api_key = input$api_key, force = TRUE),
                       error = function(e) NULL)
        if (!is.null(ch)) updateSelectInput(session, "model", choices = ch,
                                            selected = isolate(input$model))
      }
```

Append to the `input$model` observer body:

```r
      if (!is.null(auth_state) && !is.null(auth_state$current_project_id)) {
        update_project_metadata(auth_state$current_project_id, auth_state$user_id,
                                ai_model = input$model)
      }
```

- [ ] **Step 5: Checkpoint** — launch app; open a project with NO key; go to **AI Settings**; paste a key + choose a model; see "API key configured"; click **Test API Connection** → success. Switch to another tab and use an AI feature → works. Reopen the project later → the key/model persisted. **This closes the "forgot the key" trap from the in-app side.**

---

## Task 9: Fix credential leak on reset + de-freeze project save

**Files:**
- Modify: `app.R:1117-1153` (`.reset_values`)
- Modify: `app.R:1039-1077` (switch project) and `app.R:1080-1114` (logout) save blocks

- [ ] **Step 1: Clear AI credentials in `.reset_values`**

In `app.R` `.reset_values()`, add (near `values$ai_cache <- list()`):

```r
    values$ai_api_key        <- NULL
    values$ai_model          <- get_default_model()
    values$ai_custom_persona <- ""
    values$ai_project_summary <- ""
```

- [ ] **Step 2: Wrap the switch-project save in a visible progress modal**

In `observeEvent(input$sidebar_switch_project, {...})`, wrap the `save_project_state(pid, uid, session_data)` call so the UI shows "Saving…" instead of appearing frozen. Replace the `tryCatch({ ... save_project_state(...) }, error = ...)` with:

```r
    withProgress(message = "Saving project…", value = 0.5, {
      tryCatch({
        session_data <- list(
          seurat_list = values$seurat_list, active_dataset = values$active_dataset,
          active_name = values$active_name, cluster_colors = values$cluster_colors,
          analysis_flags = list(
            data_loaded = values$data_loaded, qc_done = values$qc_done,
            normalization_done = values$normalization_done,
            preprocessing_done = values$preprocessing_done,
            features_selected = values$features_selected,
            scaling_done = values$scaling_done, dimred_done = values$dimred_done,
            clustering_done = values$clustering_done, annotation_done = values$annotation_done,
            trajectory_done = values$trajectory_done),
          dataset_version = values$dataset_version, dataset_signature = values$dataset_signature,
          dataset_history = values$dataset_history, analysis_audit = values$analysis_audit,
          job_queue = values$job_queue, method_inventory = values$method_inventory)
        setProgress(0.7, detail = "Writing to disk…")
        save_project_state(pid, uid, session_data)
        setProgress(1, detail = "Done")
      }, error = function(e) showNotification(paste("Save failed:", conditionMessage(e)),
                                              type = "error", duration = 8))
    })
```

Apply the same `withProgress` wrap to the identical save block in `observeEvent(input$sidebar_logout, {...})`.

- [ ] **Step 3: Checkpoint** — launch app; load a dataset; click **Switch Project**. Confirm a "Saving… / Writing to disk…" progress bar shows (no frozen-looking UI), the dashboard returns, and reopening restores state. Confirm the previous project's key does NOT leak into a different project opened next (open a keyless project → AI Settings shows "No API key set").

---

## Final verification

- [ ] **Run the full test suite:** `R -e "testthat::test_dir('tests/testthat')"` — all green (incl. the two new files; the pre-existing DE/gallery/pn tests still pass).
- [ ] **Manual smoke test (end-to-end):** register → create keyless project → open → load mouse data → QC shows non-zero `percent.mt` (Task 4) → no "Ask AI" buttons anywhere (Task 5) → add key via AI Settings (Task 8) → run an AI feature → normalize + score cell cycle → QC report shows Cell Cycle section (Task 6) → edit project name/model from dashboard pencil (Task 7) → switch project shows Saving progress (Task 9) → new keyless project has no leaked key (Task 9).
- [ ] **Grep guard:** `grep -rn "ai_helper" --include="*.R" modules/ app.R` → no output; `grep -rn '"claude-' --include="*.R" app.R modules/ | grep -v modules/ai_models.R` → no output.

---

## Self-review (author checklist — completed)

- **Spec coverage:** A1→Tasks 1-2; A2→Tasks 3-4; A3→Task 6; A4→Task 5; A5→Tasks 7,8,9. All spec sections mapped.
- **Type consistency:** `list_claude_models`/`ai_model_choices`/`get_default_model`/`detect_species`/`normalize_species`/`species_patterns`/`update_project_metadata` used with the signatures they are defined with.
- **Placeholders:** none — Task 6 contains explicit investigation steps (reading specified line ranges) because the exact cell-cycle/report code must be confirmed against the live file; all new code to add is provided in full.
