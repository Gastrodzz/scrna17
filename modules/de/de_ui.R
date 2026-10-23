# modules/de/de_ui.R
# UI for Differential Expression. Sourced by de_server.R (which is sourced from app.R).
#
# Design system: minimalism + scientific dashboard.
# - Uses app.css design tokens (--primary, --surface, --text-*, --radius, --t)
# - SVG icons via shiny::icon (Font Awesome) — no emojis
# - 175ms transitions, visible focus rings, 4.5:1 contrast
# - Empty-state placeholder before first run; status strip with dataset context
# - Single CTA per panel ("Run DE") — green (success) for primary action

deUI <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    .de_styles(id),

    shiny::div(class = paste0("de-shell-", id),

      # ============ CONFIG (left) ============
      shiny::div(class = paste0("de-config-", id),

        # Header
        shiny::div(class = "de-header",
          shiny::div(class = "de-header-icon", shiny::icon("flask")),
          shiny::div(
            shiny::div(class = "de-title", "Differential Expression"),
            shiny::div(class = "de-subtitle", "Find genes that differ between groups")
          )
        ),

        # Dataset status strip
        shiny::uiOutput(ns("dataset_status")),

        # ── Section: comparison ─────────────────────────────────────────────
        .de_section("Comparison", "filter",
          shiny::selectInput(ns("groupby"), "Group by",
            choices = NULL, width = "100%"),
          shiny::selectInput(ns("ident1"),  "Group 1 (test)",
            choices = NULL, width = "100%"),
          shiny::selectInput(ns("ident2"),  "Group 2 (reference)",
            choices = c("rest", "all"), width = "100%")
        ),

        # ── Section: method ─────────────────────────────────────────────────
        .de_section("Method", "calculator",
          shinyWidgets::radioGroupButtons(
            ns("analysis_type"), label = NULL,
            choices = c("Single-cell" = "singlecell",
                        "Pseudobulk"  = "pseudobulk"),
            selected = "singlecell", justified = TRUE, size = "sm",
            status = "primary"
          ),
          shiny::div(class = "de-method-help",
            shiny::icon("info-circle"),
            " Single-cell: every cell as a sample. Pseudobulk: aggregate counts per replicate."
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'singlecell'", ns("analysis_type")),
            shiny::div(class = "de-pill-grid",
              shinyWidgets::radioGroupButtons(
                ns("sc_method"), label = NULL,
                choices = c("Wilcoxon" = "wilcox", "MAST" = "MAST",
                            "t-test"  = "t",      "Bimod" = "bimod",
                            "LR"      = "LR",     "Negbinom" = "negbinom"),
                selected = "wilcox", size = "sm", individual = TRUE,
                status = "primary-outline"
              )
            )
          ),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] == 'pseudobulk'", ns("analysis_type")),
            shinyWidgets::radioGroupButtons(
              ns("pb_method"), label = NULL,
              choices = c("DESeq2" = "DESeq2", "edgeR" = "edgeR",
                          "limma-voom" = "limma-voom"),
              selected = "DESeq2", justified = TRUE, size = "sm",
              status = "primary-outline"
            ),
            shiny::div(style = "margin-top: 8px;",
              shiny::selectInput(ns("sample_col"),
                "Sample / replicate column", choices = NULL, width = "100%")),
            shiny::div(class = "de-callout de-callout-warning",
              shiny::icon("exclamation-triangle"),
              shiny::span("Pseudobulk requires ≥2 replicates per group."))
          )
        ),

        # ── Section: optional cell-type subset ──────────────────────────────
        .de_section("Subset (optional)", "layer-group",
          shiny::selectInput(ns("subset_col"), "Cell-type column",
            choices = c("(none)" = ""), width = "100%"),
          shiny::conditionalPanel(
            condition = sprintf("input['%s'] != ''", ns("subset_col")),
            shiny::selectInput(ns("subset_val"), "Cell type",
              choices = NULL, width = "100%"))
        ),

        # ── Section: covariates ─────────────────────────────────────────────
        .de_section("Covariates", "sliders-h",
          shiny::selectInput(ns("latent_vars"), "Latent variables (single-cell)",
            choices = NULL, multiple = TRUE, width = "100%"),
          shiny::div(class = "de-help",
            "MAST: cellular detection rate (CDR) is added automatically.")
        ),

        # ── Section: display thresholds ─────────────────────────────────────
        .de_section("Display thresholds", "eye",
          shiny::div(class = "de-grid-2",
            shiny::numericInput(ns("min_pct_display"), "Min % expressing",
              value = 0.10, min = 0, max = 1, step = 0.01),
            shiny::numericInput(ns("logfc_threshold"), "|log2FC| display",
              value = 0.10, min = 0, max = 5, step = 0.05),
            shiny::numericInput(ns("pval_cutoff"), "Adj p (highlight)",
              value = 0.05, min = 1e-6, max = 0.5, step = 0.01),
            shiny::numericInput(ns("label_top_n"), "Top genes shown",
              value = 20, min = 5, max = 200, step = 5)
          ),
          shiny::checkboxInput(ns("only_pos"), "Only positive markers",
            value = FALSE)
        ),

        # ── Section: reproducibility ────────────────────────────────────────
        .de_section("Reproducibility", "shield-alt",
          shiny::numericInput(ns("seed"), "Random seed",
            value = 42, min = 1, max = 99999, step = 1, width = "100%"),
          shiny::div(class = "de-callout de-callout-info",
            shiny::icon("check-circle"),
            shiny::span("DE always tests on the ", shiny::strong("full dataset"),
                        " — no downsampling."))
        ),

        # ── Section: run label ──────────────────────────────────────────────
        .de_section("Save to history", "bookmark",
          shiny::textInput(ns("run_label"),
            "Label (optional — auto-generated if blank)",
            value = "", width = "100%",
            placeholder = "e.g. KO vs WT (osteoblast)"),
          shiny::div(class = "de-help",
            "Each successful run is saved to the in-session history (max 50). ",
            "Pick any 2+ in ", shiny::strong("Comparative Enrichment"),
            " to enrich side-by-side.")
        ),

        # CTA + per-cluster loop
        shiny::div(class = "de-cta-wrap",
          shiny::actionButton(ns("run"), "Run DE Analysis",
            class = "de-cta", icon = shiny::icon("play")),
          shiny::br(), shiny::br(),
          shiny::actionButton(ns("run_per_cluster"),
            "Run DE in every cell type (loop)",
            class = "de-cta-secondary", icon = shiny::icon("repeat")),
          shiny::div(class = "de-help",
            "Loops the comparison above across every value of the ",
            shiny::strong("Subset by cell-type column"),
            " — populates the history with one run per cell type.")
        )
      ),

      # ============ RESULTS (right) ============
      shiny::div(class = paste0("de-results-", id),
        shiny::uiOutput(ns("results_ui"))
      )
    )
  )
}

# ---------------------------------------------------------------------------
# Style block — scoped via the module id so two instances don't collide.
# Uses app.css design tokens: --primary, --surface, --text-*, --radius, --t.
# ---------------------------------------------------------------------------
.de_styles <- function(id) {
  # Use gsub('__ID__', id, ...) instead of sprintf to keep CSS percent signs
  # literal (sprintf would treat 'width: 100%;' as a format specifier).
  css <- "
.de-shell-__ID__ {
  display: flex; gap: 18px; align-items: stretch;
  font-family: var(--font-sans);
  color: var(--text-primary);
}
.de-config-__ID__ {
  flex: 0 0 360px; max-height: 88vh; overflow-y: auto;
  padding: 18px;
  background: var(--surface);
  border: 1px solid var(--border);
  border-radius: var(--radius-lg);
  box-shadow: var(--shadow-sm);
}
.de-config-__ID__::-webkit-scrollbar { width: 8px; }
.de-config-__ID__::-webkit-scrollbar-thumb {
  background: var(--border); border-radius: 4px;
}
.de-results-__ID__ { flex: 1 1 auto; min-width: 0; }

.de-config-__ID__ .de-header {
  display: flex; align-items: center; gap: 12px;
  padding-bottom: 14px; margin-bottom: 14px;
  border-bottom: 1px solid var(--border);
}
.de-config-__ID__ .de-header-icon {
  width: 38px; height: 38px; border-radius: var(--radius);
  background: var(--primary-light); color: var(--primary);
  display: flex; align-items: center; justify-content: center;
  font-size: 18px;
}
.de-config-__ID__ .de-title {
  font-weight: 700; font-size: 16px; color: var(--text-primary);
  letter-spacing: -0.01em;
}
.de-config-__ID__ .de-subtitle {
  font-size: 12px; color: var(--text-secondary); margin-top: 2px;
}

.de-config-__ID__ .de-section {
  margin: 18px 0 4px; padding-bottom: 14px;
  border-bottom: 1px dashed var(--border-light);
}
.de-config-__ID__ .de-section:last-of-type { border-bottom: none; }
.de-config-__ID__ .de-section-title {
  display: flex; align-items: center; gap: 8px;
  font-weight: 600; font-size: 11px; letter-spacing: 0.06em;
  text-transform: uppercase; color: var(--text-secondary);
  margin-bottom: 10px;
}
.de-config-__ID__ .de-section-title .icon {
  color: var(--primary); width: 14px; text-align: center;
}

.de-config-__ID__ .form-control,
.de-config-__ID__ .selectize-input {
  border-radius: var(--radius-sm) !important;
  border: 1px solid var(--border) !important;
  font-size: 13px;
  transition: border-color var(--t), box-shadow var(--t);
}
.de-config-__ID__ .form-control:focus,
.de-config-__ID__ .selectize-input.focus {
  border-color: var(--primary) !important;
  box-shadow: 0 0 0 3px var(--border-focus) !important;
  outline: none;
}
.de-config-__ID__ label {
  font-size: 12px; font-weight: 500; color: var(--text-secondary);
  margin-bottom: 4px;
}

.de-config-__ID__ .de-help {
  font-size: 11px; color: var(--text-muted); margin-top: 6px;
  line-height: 1.5;
}
.de-config-__ID__ .de-method-help {
  font-size: 11px; color: var(--text-secondary); margin: 8px 0;
  padding: 8px 10px; background: var(--primary-light);
  border-radius: var(--radius-sm);
  display: flex; gap: 8px; align-items: flex-start;
}
.de-config-__ID__ .de-method-help i { color: var(--primary); margin-top: 2px; }

.de-config-__ID__ .de-pill-grid { margin-top: 10px; }
.de-config-__ID__ .btn-group-justified { width: 100%; }

.de-config-__ID__ .de-grid-2 {
  display: grid; grid-template-columns: 1fr 1fr; gap: 10px;
}

.de-config-__ID__ .de-callout {
  display: flex; gap: 8px; align-items: flex-start;
  padding: 8px 10px; border-radius: var(--radius-sm);
  font-size: 11px; line-height: 1.5; margin-top: 8px;
}
.de-config-__ID__ .de-callout-warning {
  background: var(--warning-light); color: var(--warning);
}
.de-config-__ID__ .de-callout-warning i { margin-top: 2px; }
.de-config-__ID__ .de-callout-warning span { color: var(--text-primary); }
.de-config-__ID__ .de-callout-info {
  background: var(--success-light); color: var(--success);
}
.de-config-__ID__ .de-callout-info i { margin-top: 2px; }
.de-config-__ID__ .de-callout-info span { color: var(--text-primary); }

/* Dataset status strip */
.de-config-__ID__ .de-dataset-strip {
  display: flex; gap: 10px; align-items: center;
  padding: 10px 12px; border-radius: var(--radius-sm);
  background: var(--bg); border: 1px solid var(--border-light);
  font-size: 11px; color: var(--text-secondary); margin-bottom: 10px;
}
.de-config-__ID__ .de-dataset-strip .dot {
  width: 8px; height: 8px; border-radius: 50%;
  background: var(--success); flex: 0 0 auto;
}
.de-config-__ID__ .de-dataset-strip .dot.muted { background: var(--text-muted); }
.de-config-__ID__ .de-dataset-strip .meta { color: var(--text-primary); font-weight: 500; }

/* Run button — single CTA */
.de-config-__ID__ .de-cta-wrap { margin-top: 18px; }
.de-config-__ID__ .de-cta {
  width: 100%; padding: 12px 16px;
  background: var(--success); color: #fff;
  border: none; border-radius: var(--radius);
  font-weight: 600; font-size: 14px; letter-spacing: 0.01em;
  cursor: pointer;
  box-shadow: 0 1px 2px rgba(22,163,74,0.18);
  transition: background var(--t), box-shadow var(--t), transform var(--t);
}
.de-config-__ID__ .de-cta:hover {
  background: #15803D;
  box-shadow: 0 3px 8px rgba(22,163,74,0.25);
}
.de-config-__ID__ .de-cta:active { transform: translateY(1px); }
.de-config-__ID__ .de-cta:focus-visible {
  outline: none;
  box-shadow: 0 0 0 3px rgba(22,163,74,0.30);
}
.de-config-__ID__ .de-cta i { margin-right: 8px; }

.de-config-__ID__ .de-cta-secondary {
  width: 100%; padding: 10px 14px;
  background: transparent; color: var(--primary);
  border: 1px solid var(--primary); border-radius: var(--radius);
  font-weight: 600; font-size: 13px; cursor: pointer;
  transition: background var(--t), color var(--t);
}
.de-config-__ID__ .de-cta-secondary:hover {
  background: var(--primary-light);
}
.de-config-__ID__ .de-cta-secondary:focus-visible {
  outline: none;
  box-shadow: 0 0 0 3px var(--border-focus);
}
.de-config-__ID__ .de-cta-secondary i { margin-right: 8px; }
.de-config-__ID__ .de-cta-secondary:disabled {
  opacity: 0.45; cursor: not-allowed;
}

/* Results panel */
.de-results-__ID__ .de-empty {
  display: flex; flex-direction: column; align-items: center;
  justify-content: center; min-height: 60vh;
  padding: 40px; text-align: center;
  background: var(--surface); border: 1px dashed var(--border);
  border-radius: var(--radius-lg);
}
.de-results-__ID__ .de-empty .icon {
  width: 64px; height: 64px; border-radius: 50%;
  background: var(--primary-light); color: var(--primary);
  display: flex; align-items: center; justify-content: center;
  font-size: 26px; margin-bottom: 16px;
}
.de-results-__ID__ .de-empty h4 {
  margin: 0 0 8px; color: var(--text-primary); font-weight: 600;
}
.de-results-__ID__ .de-empty p {
  color: var(--text-secondary); font-size: 13px; max-width: 420px;
  line-height: 1.55; margin: 0;
}
.de-results-__ID__ .de-results-card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: var(--radius-lg); padding: 16px;
  box-shadow: var(--shadow-sm);
}
.de-results-__ID__ .nav-tabs > li > a {
  font-size: 13px; font-weight: 500; color: var(--text-secondary);
  border-radius: 0; transition: color var(--t), border-color var(--t);
}
.de-results-__ID__ .nav-tabs > li.active > a {
  color: var(--primary); border-bottom: 2px solid var(--primary);
}

/* ---- Volcano summary panel (cards + venn + stats) -------------------- */
.de-results-__ID__ .de-summary-section {
  margin-top: 22px;
  display: flex; flex-direction: column; gap: 16px;
}
.de-results-__ID__ .de-section-h {
  display: flex; align-items: center; gap: 8px;
  font-weight: 600; font-size: 12px; letter-spacing: 0.06em;
  text-transform: uppercase; color: var(--text-secondary);
  margin: 0 0 8px;
}
.de-results-__ID__ .de-section-h .icon { color: var(--primary); }

.de-results-__ID__ .de-metrics {
  display: grid; grid-template-columns: repeat(4, 1fr); gap: 12px;
}
.de-results-__ID__ .de-metric-card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: var(--radius); padding: 14px 16px;
  display: flex; flex-direction: column; gap: 4px;
  position: relative; overflow: hidden;
  transition: box-shadow var(--t), transform var(--t);
}
.de-results-__ID__ .de-metric-card:hover {
  box-shadow: var(--shadow);
}
.de-results-__ID__ .de-metric-card .label {
  font-size: 11px; color: var(--text-secondary);
  text-transform: uppercase; letter-spacing: 0.04em; font-weight: 600;
}
.de-results-__ID__ .de-metric-card .value {
  font-size: 26px; font-weight: 700; color: var(--text-primary);
  line-height: 1.1; letter-spacing: -0.02em;
  font-variant-numeric: tabular-nums;
}
.de-results-__ID__ .de-metric-card .sub {
  font-size: 11px; color: var(--text-muted); margin-top: 2px;
}
.de-results-__ID__ .de-metric-card .accent {
  position: absolute; left: 0; top: 0; bottom: 0; width: 3px;
}
.de-results-__ID__ .de-metric-card.up   .value { color: var(--danger); }
.de-results-__ID__ .de-metric-card.up   .accent { background: var(--danger); }
.de-results-__ID__ .de-metric-card.down .value { color: var(--primary); }
.de-results-__ID__ .de-metric-card.down .accent { background: var(--primary); }
.de-results-__ID__ .de-metric-card.neutral .accent { background: var(--info); }
.de-results-__ID__ .de-metric-card.success .accent { background: var(--success); }

.de-results-__ID__ .de-split-row {
  display: grid; grid-template-columns: 1.1fr 1fr; gap: 16px;
}
.de-results-__ID__ .de-card {
  background: var(--surface); border: 1px solid var(--border);
  border-radius: var(--radius); padding: 16px;
}
.de-results-__ID__ .de-card-title {
  font-size: 13px; font-weight: 600; color: var(--text-primary);
  margin: 0 0 12px;
}

.de-results-__ID__ .de-cond-block {
  display: flex; flex-direction: column; gap: 10px;
}
.de-results-__ID__ .de-cond-row {
  display: flex; gap: 12px; align-items: center;
  padding: 12px 14px;
  background: var(--bg); border: 1px solid var(--border-light);
  border-radius: var(--radius-sm);
}
.de-results-__ID__ .de-cond-row .badge-grp {
  width: 32px; height: 32px; border-radius: 50%;
  display: flex; align-items: center; justify-content: center;
  font-weight: 700; font-size: 13px; color: #fff;
  flex: 0 0 auto;
}
.de-results-__ID__ .de-cond-row .badge-grp.g1 { background: var(--danger); }
.de-results-__ID__ .de-cond-row .badge-grp.g2 { background: var(--primary); }
.de-results-__ID__ .de-cond-row .meta { display: flex; flex-direction: column; }
.de-results-__ID__ .de-cond-row .name {
  font-size: 14px; font-weight: 600; color: var(--text-primary);
}
.de-results-__ID__ .de-cond-row .sub {
  font-size: 11px; color: var(--text-secondary);
}
.de-results-__ID__ .de-cond-vs {
  text-align: center; font-size: 11px; color: var(--text-muted);
  text-transform: uppercase; letter-spacing: 0.1em; font-weight: 600;
  margin: -4px 0;
}
.de-results-__ID__ .de-cond-meta {
  display: grid; grid-template-columns: auto 1fr; gap: 6px 12px;
  margin-top: 10px; font-size: 12px;
}
.de-results-__ID__ .de-cond-meta dt {
  color: var(--text-secondary); font-weight: 500;
}
.de-results-__ID__ .de-cond-meta dd {
  margin: 0; color: var(--text-primary); font-weight: 500;
  font-variant-numeric: tabular-nums;
}

.de-results-__ID__ .de-stat-table {
  width: 100%; border-collapse: collapse; font-size: 12px;
}
.de-results-__ID__ .de-stat-table td {
  padding: 8px 10px; border-bottom: 1px solid var(--border-light);
}
.de-results-__ID__ .de-stat-table tr:last-child td { border-bottom: none; }
.de-results-__ID__ .de-stat-table td:first-child {
  color: var(--text-secondary); width: 55%;
}
.de-results-__ID__ .de-stat-table td:last-child {
  font-weight: 600; color: var(--text-primary);
  font-variant-numeric: tabular-nums;
  font-family: var(--font-mono);
  text-align: right;
}
.de-results-__ID__ .de-stat-table tr:hover td {
  background: var(--primary-light);
}

@media (max-width: 1100px) {
  .de-results-__ID__ .de-metrics { grid-template-columns: repeat(2, 1fr); }
  .de-results-__ID__ .de-split-row { grid-template-columns: 1fr; }
}
"
  shiny::tags$style(shiny::HTML(gsub("__ID__", id, css, fixed = TRUE)))
}

# ---------------------------------------------------------------------------
# Section helper — uniform "[icon] TITLE" stripe + slot for inputs.
# ---------------------------------------------------------------------------
.de_section <- function(title, icon, ...) {
  shiny::div(class = "de-section",
    shiny::div(class = "de-section-title",
      shiny::span(class = "icon", shiny::icon(icon)),
      shiny::span(title)
    ),
    ...
  )
}
