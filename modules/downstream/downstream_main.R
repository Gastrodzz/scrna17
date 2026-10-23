# ──────────────────────────────────────────────────────────────────────────────
# downstream_main.R — Orchestrator for all downstream analysis sub-modules.
# app.R sources only this file; it then sources each sub-module in order.
#
# Module hierarchy:
#   downstream_core.R        → shared helpers & data accessors
#   downstream_ccc.R         → Cell-Cell Communication (consumed by pathway_network)
#   downstream_tf.R          → TF Activity (GRN removed; consumed by pathway_network)
#   downstream_geneset.R     → Gene Set / Module Scoring (consumed by pathway_network)
#   downstream_composition.R → Cell composition (consumed by pathway_network)
#   downstream_ai.R          → AI Analysis tab (currently unmounted; will return)
#
# Differential Expression now lives in modules/de/ (sourced by app.R directly).
# Pathway & Network Analysis lives in modules/pathway_network/ and re-uses
# downstream_ccc.R, downstream_tf.R (TF activity only), downstream_geneset.R,
# downstream_composition.R via its own server.
# ──────────────────────────────────────────────────────────────────────────────

local({
  # Robust module directory resolution with multiple fallback strategies
  module_dir <- tryCatch({
    # Strategy 1: Relative to this file's location (works when sourced directly)
    this_file <- sys.frame(1)$ofile
    if (!is.null(this_file) && nzchar(this_file)) {
      normalizePath(file.path(dirname(this_file), "downstream"), mustWork = TRUE)
    } else {
      stop("ofile not available")
    }
  }, error = function(e) {
    # Strategy 2: Try well-known relative paths from the working directory
    candidates <- c(
      "modules/downstream",
      file.path(getwd(), "modules", "downstream"),
      # Handle case where app.R is in a subdirectory
      file.path(dirname(getwd()), "modules", "downstream")
    )
    # Also check relative to app root if shiny is running
    if (requireNamespace("shiny", quietly = TRUE)) {
      app_dir <- tryCatch(shiny::getShinyOption("appDir", NULL), error = function(e) NULL)
      if (!is.null(app_dir)) {
        candidates <- c(file.path(app_dir, "modules", "downstream"), candidates)
      }
    }
    for (cand in candidates) {
      resolved <- tryCatch(normalizePath(cand, mustWork = TRUE), error = function(e) NULL)
      if (!is.null(resolved) && dir.exists(resolved)) return(resolved)
    }
    stop("Could not locate modules/downstream directory. Searched: ",
         paste(candidates, collapse = ", "))
  })

  if (!dir.exists(module_dir)) {
    stop("[downstream] Module directory not found: ", module_dir,
         ". Ensure modules/downstream/ exists relative to the app root.")
  }

  sub_modules <- c(
    "downstream_core.R",
    "downstream_ccc.R",
    "downstream_tf.R",
    "downstream_geneset.R",
    "downstream_composition.R",
    "downstream_ai.R"
  )

  for (mod in sub_modules) {
    mod_path <- file.path(module_dir, mod)
    if (file.exists(mod_path)) {
      tryCatch(
        source(mod_path, local = FALSE),
        error = function(e) {
          message("[downstream] Failed to source ", mod, ": ", e$message)
        }
      )
    } else {
      message("[downstream] Sub-module not found: ", mod_path)
    }
  }
})
