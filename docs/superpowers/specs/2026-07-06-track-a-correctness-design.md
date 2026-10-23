# Track A — Correctness & Robustness Design

**Date:** 2026-07-06
**Scope:** scRNA-seq R/Shiny platform — Track A (correctness bugs + robustness, low risk, high ROI)
**Status:** Approved (design), pending implementation plan

Track A covers five issues from the code review, chosen because they fix real
correctness bugs and remove dead/confusing paths without touching the heavy
analytical engines. Security (API-key encryption, password hashing) is **Track B**
and is explicitly out of scope here. New annotation methods, downstream hardening,
and AI re-architecture are **Track C**.

---

## A1. Central model registry (live fetch + static fallback)

**Problem.** Model IDs are hardcoded in 5 places with drifted defaults
(`claude-sonnet-4-6` vs `claude-sonnet-4-5-20250929`). No call to Anthropic's
`GET /v1/models`, so new models never appear automatically.

**Design.**
- New file `modules/ai_models.R`, sourced **before** all AI-consuming code.
- `list_claude_models(api_key = NULL, force = FALSE)`:
  - When `api_key` is present, `GET https://api.anthropic.com/v1/models` with
    headers `x-api-key`, `anthropic-version: 2023-06-01`.
  - Parse `data[].id` + `data[].display_name`; cache in a module-level environment
    with a timestamp; reuse if `< 3600s` old (unless `force`).
  - On no key / offline / any error → return curated `AI_MODELS_FALLBACK`.
- Extended thinking stays working: `/v1/models` returns only real IDs, so generate
  `…-thinking` display variants programmatically for Opus/Sonnet-family models.
  `call_claude_api()` already strips the `-thinking` suffix and sets the thinking
  request body — unchanged.
- Public helpers: `get_default_model(api_key)`, `ai_model_choices(api_key)`
  (named vector for `selectInput`).
- Replace all hardcoded copies:
  - `ai_assistant.R` `AI_MODEL_CONFIG$available_models` + `default_model`
  - `ai_assistant.R` `aiConfigUI` dropdown → dynamic
  - `project_module.R` create dropdown → dynamic (fallback when key blank;
    refresh on key paste)
  - `app.R:1177` + `.reset_values` default → `get_default_model()`
  - `annotation.R` / `downstream_ai.R` `-4-5-` defaults → `%||% get_default_model()`
  - `get_maxtokens()` hardcoded `== "claude-opus-4-6"` → `grepl("opus"/"haiku")`

**Interface.** Every model dropdown and every "default model" fallback in the app
resolves through `ai_models.R`. Consumers never hardcode an ID again.

---

## A2. Species auto-detection + silent mito bug

**Problem.** `.detect_species()` exists only in annotation. QC uses a manual
`human/mouse/custom` dropdown. `app.R:952` pushes the project species value
(`"Homo sapiens"`) into that dropdown whose choices are `"human"/"mouse"` →
no match → silent fallback to `"human"`. For a mouse dataset, MT pattern `^MT-`
finds 0 genes → `percent.mt = 0` for all cells → mito filtering silently disabled.

**Design.**
- Hoist detection into shared `modules/species_utils.R`:
  - `detect_species(genes)` → `list(species, mt_pattern, ribo_pattern, confidence)`.
    Recognizes human/mouse robustly; flags Ensembl-ID matrices (warn); minimal
    rat/zebrafish support.
  - `normalize_species(x)` maps `"Homo sapiens"→"human"`, `"Mus musculus"→"mouse"`,
    etc.
- On data load: auto-run detection, set the QC species select + MT pattern, show
  an "auto-detected: human (high confidence)" hint. If project species and detected
  species disagree, **warn** rather than silently default.
- Reuse the same util in annotation, preprocessing, and any species-dependent
  module (single source of truth). Ribo pattern made species-aware; note MRPL/MRPS.

---

## A3. Cell-cycle consolidation

**Problem.** Cell-cycle scoring lives in both the QC tab (`cellCycleServer`) and
Preprocessing; it needs normalized data, so "run it straight away" in raw QC fails.
Comment at `quality_control.R:234` says it "moved to Preprocessing" but the QC tab
was never removed.

**Design.**
- Remove cell-cycle UI/compute from the QC **Metrics** area. Scoring stays in
  **Preprocessing** (post-normalization): `S.Score / G2M.Score / Phase`.
- Keep **one dedicated Cell Cycle tab** (current sidebar position) that visualizes
  the scores; if scoring hasn't run, guide the user to Preprocessing rather than
  erroring on raw data.
- Add a **Cell-Cycle section to the final QC report**, and improve the report:
  structured sections, before/after cell counts, thresholds used, doublet + ambient
  + CC summary. Exact wiring of `cellCycleServer` vs. Preprocessing `cc_*` inputs to
  be de-duplicated during implementation.

---

## A4. Remove the "Ask AI" buttons

**Problem.** Generic `ai_helper_button_ui()` / `ai_helper_server()` "ask about your
data" boxes are wired into 8 modules, redundant with the main assistant, and
sometimes called without an `exists()` guard (trajectory, dim-reduction) so a
failed AI load errors those tabs.

**Design.**
- Remove `ai_helper_button_ui()` / `ai_helper_server()` from visualization,
  trajectory, export, preprocessing, quality_control, dimensionality_reduction,
  clustering; delete the now-dead definitions in `ai_assistant.R`.
- Keep purposeful, context-specific AI features (annotation, QC threshold advisor,
  DE test advisor, downstream interpretation).

---

## A5. Edit path + the "stuck" bug (both routes)

**Problem.** There is no edit-project UI (dashboard is create/open/delete only) and
the in-app AI config panel is dead code (`app.R:1661 if (FALSE …)`). A user who
forgot the API key has no path to add it. "Switch Project" also runs a synchronous
`saveRDS` of the whole Seurat object → single-thread freeze that reads as "stuck".
`.reset_values()` never clears `ai_*`, leaking credentials across projects.

**Design.**
- **Dashboard "Edit Project"**: pencil action on each card → edit modal (reuses the
  create-modal layout), prefilled via `get_project()`, saved through the existing
  `update_project_metadata()`. If the edited project is currently open, sync
  `auth_state` + `values` so changes take effect immediately.
- **In-app "AI Settings" panel**: re-enable a cleaned-up `aiConfigServer` in the
  sidebar to set/change API key + model live mid-analysis; writes to `values` and
  persists via `update_project_metadata()`. Masked current key + "Test key" button.
- **Fix the freeze/stuck**:
  - Root cause (no edit path) resolved by the two routes above.
  - `.reset_values()` also clears `ai_api_key/ai_model/ai_custom_persona/ai_project_summary`.
  - Wrap switch/logout `save_project_state` (large synchronous `saveRDS`) in a
    visible progress modal + re-entrancy guard so it reads as "Saving…" not frozen.
  - Full background/async save is an optional stretch (heavier, riskier) — flagged,
    not committed.

**Scope guard.** API-key encryption + password hashing hardening are Track B. A5
does not worsen the plaintext situation; Track B encrypts the column afterward.

---

## Cross-cutting notes

- **Testing.** User can run the app locally. Each change gets a concrete manual test
  step (e.g., load a mouse dataset → confirm `percent.mt` is non-zero; open project
  without key → add via AI Settings → confirm an AI call works).
- **Sourcing order.** `ai_models.R` and `species_utils.R` must load before modules
  that consume them; verify in `app.R`'s source block.
- **No git repo.** Spec saved to disk; not committed (project is not under git).
