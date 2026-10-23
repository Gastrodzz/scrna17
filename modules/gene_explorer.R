# =============================================================================
# Gene Explorer Module (Trailmaker-style)
# -----------------------------------------------------------------------------
# Three-panel layout (UMAP + Cell sets/Annotate / Genes+DE) and bottom DE
# heatmap. Lets the user build cell sets by lasso, assign annotation labels
# per set, and commit a consolidated `user_annotation` factor that all
# downstream modules can consume.
#
# Metadata contract:
#   - `gex_<set_name>`         logical per-cell mask per cell set
#   - `gene_explorer_sets`     factor: last-write-wins label per cell set
#   - `user_annotation`        factor: user-assigned cluster labels (kept
#                              up-to-date whenever "Apply annotations" runs)
# =============================================================================

# --- helpers -----------------------------------------------------------------

.ge_safe_default_assay <- function(obj) {
  if (exists("safe_default_assay", mode = "function"))
    tryCatch(safe_default_assay(obj), error = function(e) "RNA")
  else tryCatch(DefaultAssay(obj), error = function(e) "RNA")
}

.ge_safe_get_assay_data <- function(obj, layer = "data", assay = NULL) {
  if (exists("safe_get_assay_data", mode = "function")) {
    tryCatch(safe_get_assay_data(obj, layer, assay = assay),
             error = function(e) NULL)
  } else {
    tryCatch({
      a <- assay %||% DefaultAssay(obj)
      GetAssayData(obj, layer = layer, assay = a)
    }, error = function(e) {
      tryCatch(GetAssayData(obj, slot = layer), error = function(e2) NULL)
    })
  }
}

.ge_prepare_for_de <- function(obj) {
  if (exists("prepare_for_de", mode = "function"))
    tryCatch(prepare_for_de(obj), error = function(e) obj)
  else obj
}

.ge_best_reduction <- function(obj) {
  reds <- names(obj@reductions)
  for (r in c("umap", "tsne", "pca")) if (r %in% reds) return(r)
  if (length(reds) > 0) reds[1] else NULL
}

.ge_clip <- function(x, lo = 0.0, hi = 0.99) {
  if (!length(x)) return(x)
  qs <- stats::quantile(x, c(lo, hi), na.rm = TRUE)
  # Sparse genes (expressed in <(1-hi)*100% of cells) make the upper quantile
  # collapse onto the lower one — e.g. both 0 — which would flatten every cell
  # to a single value and hide the signal. Fall back to the true max so the
  # positive cells stay visible.
  if (!is.finite(qs[2]) || qs[2] <= qs[1]) {
    mx <- suppressWarnings(max(x, na.rm = TRUE))
    if (!is.finite(mx) || mx <= qs[1]) return(x)
    qs[2] <- mx
  }
  pmin(pmax(x, qs[1]), qs[2])
}

# Build a (canonical -> assays) registry of all features across assays once per
# session. Used for robust, case-insensitive gene resolution.
.ge_build_feature_index <- function(obj) {
  assays <- names(obj@assays)
  idx <- list()  # key = toupper(gene) -> data.frame(gene, assay)
  for (a in assays) {
    feats <- rownames(obj[[a]])
    if (!length(feats)) next
    keys <- toupper(feats)
    for (i in seq_along(keys)) {
      k <- keys[i]
      if (is.null(idx[[k]])) {
        idx[[k]] <- list(gene = feats[i], assays = a)
      } else {
        idx[[k]]$assays <- unique(c(idx[[k]]$assays, a))
      }
    }
  }
  idx
}

# Resolve a user-entered gene symbol to (canonical_name, assay). Falls back
# gracefully when exact case doesn't match.
.ge_resolve_gene <- function(feat_index, gene, prefer_assay = NULL) {
  if (is.null(gene) || !nzchar(gene)) return(NULL)
  hit <- feat_index[[toupper(gene)]]
  if (is.null(hit)) return(NULL)
  assay <- if (!is.null(prefer_assay) && prefer_assay %in% hit$assays)
             prefer_assay
           else hit$assays[1]
  list(gene = hit$gene, assay = assay)
}

# Fetch log-normalized expression for a gene over a set of cells, robustly.
.ge_get_expr <- function(obj, feat_index, gene, cells, default_assay) {
  info <- .ge_resolve_gene(feat_index, gene, prefer_assay = default_assay)
  if (is.null(info)) return(NULL)
  dat <- .ge_safe_get_assay_data(obj, "data", info$assay)
  if (is.null(dat) || !info$gene %in% rownames(dat)) {
    # fallback: normalize counts on the fly for this gene
    counts <- .ge_safe_get_assay_data(obj, "counts", info$assay)
    if (is.null(counts) || !info$gene %in% rownames(counts)) return(NULL)
    keep <- intersect(cells, colnames(counts))
    if (!length(keep)) return(NULL)
    sub <- counts[info$gene, keep, drop = TRUE]
    libsz <- Matrix::colSums(counts[, keep, drop = FALSE])
    vals <- as.numeric(log1p(1e4 * sub / pmax(as.numeric(libsz), 1)))
    names(vals) <- keep
    return(list(values = vals, gene = info$gene, assay = info$assay,
                source = "counts(on-the-fly)"))
  }
  keep <- intersect(cells, colnames(dat))
  if (!length(keep)) return(NULL)
  vals <- as.numeric(dat[info$gene, keep, drop = TRUE])
  names(vals) <- keep
  list(values = vals, gene = info$gene, assay = info$assay, source = "data")
}

.ge_prefix <- "gex_"
.ge_consolidated_col <- "gene_explorer_sets"
.ge_annot_col <- "user_annotation"

# Stable qualitative palette suitable for overlaying cell sets on the UMAP.
.ge_set_palette <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#F39B7F", "#3C5488", "#8491B4",
  "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#FFA500", "#7B68EE",
  "#2CA02C", "#9467BD", "#17BECF", "#E377C2"
)

# =============================================================================
# UI
# =============================================================================

geneExplorerUI <- function(id) {
  ns <- NS(id)

  tagList(
    tags$style(HTML("
      /* ===================== Gene Explorer v2 ===================== */
      .ge-wrap { --ge-border:#e5e7eb; --ge-border-strong:#d1d5db;
                 --ge-bg:#ffffff; --ge-bg-soft:#f8fafc;
                 --ge-text:#111827; --ge-text-muted:#6b7280;
                 --ge-accent:#4DBBD5; --ge-accent-deep:#1f6f87;
                 --ge-success:#16a34a; --ge-danger:#dc2626;
                 --ge-purple:#6f42c1;
                 --ge-radius:10px; --ge-radius-sm:6px;
                 --ge-shadow:0 1px 2px rgba(16,24,40,.04); }

      .ge-wrap .ge-hero {
        display:flex; align-items:center; justify-content:space-between;
        padding:14px 18px 14px;
        background: linear-gradient(135deg,#f8fafc 0%, #eef4fb 100%);
        border:1px solid var(--ge-border); border-radius: var(--ge-radius);
        margin-bottom:14px;
      }
      .ge-wrap .ge-hero-left  { display:flex; align-items:center; gap:14px; }
      .ge-wrap .ge-hero-icon  { width:38px; height:38px; border-radius:10px;
                                 background: linear-gradient(135deg,#4DBBD5,#6f42c1);
                                 color:#fff; display:flex; align-items:center;
                                 justify-content:center; font-size:16px;
                                 box-shadow:0 4px 10px rgba(77,187,213,.28); }
      .ge-wrap .ge-hero-title { font-weight:700; font-size:16px; color:var(--ge-text);
                                 line-height:1.2; }
      .ge-wrap .ge-hero-sub   { font-size:12px; color:var(--ge-text-muted);
                                 margin-top:2px; }
      .ge-wrap .ge-steps { display:flex; gap:6px; align-items:center;
                           font-size:11.5px; color:var(--ge-text-muted); }
      .ge-wrap .ge-step  { padding:4px 10px; border-radius:999px;
                           background:#fff; border:1px solid var(--ge-border);
                           font-weight:600; color:#374151; }
      .ge-wrap .ge-step .num { color:var(--ge-accent); margin-right:5px; }
      .ge-wrap .ge-step-sep { color:#cbd5e1; }

      /* Panels */
      .ge-wrap .ge-panel {
        border:1px solid var(--ge-border); border-radius: var(--ge-radius);
        background: var(--ge-bg); padding:14px; height:100%;
        box-shadow: var(--ge-shadow);
      }
      .ge-wrap .ge-panel-header {
        display:flex; justify-content:space-between; align-items:center;
        margin-bottom:12px; padding-bottom:10px;
        border-bottom:1px solid var(--ge-border);
      }
      .ge-wrap .ge-panel-title {
        font-weight:700; font-size:13px; color:var(--ge-text);
        letter-spacing:0.02em; display:flex; align-items:center; gap:8px;
      }
      .ge-wrap .ge-panel-title .fa { color: var(--ge-accent); font-size:13px; }

      /* Input labels */
      .ge-wrap .ge-field { margin-bottom:12px; }
      .ge-wrap .ge-field > label,
      .ge-wrap .ge-field-label {
        display:block; font-size:11.5px; font-weight:600;
        color:#374151; margin-bottom:4px;
        text-transform:uppercase; letter-spacing:0.04em;
      }
      .ge-wrap .shiny-input-container { margin-bottom:0 !important; width:100%; }

      /* Controls row */
      .ge-wrap .ge-controls { display:grid; gap:10px; margin-bottom:10px; }
      .ge-wrap .ge-controls.c-4 { grid-template-columns: repeat(4,minmax(0,1fr)); }
      .ge-wrap .ge-controls.c-3 { grid-template-columns: repeat(3,minmax(0,1fr)); }
      .ge-wrap .ge-controls.c-2 { grid-template-columns: repeat(2,minmax(0,1fr)); }
      @media (max-width: 900px) {
        .ge-wrap .ge-controls.c-4 { grid-template-columns: repeat(2,minmax(0,1fr)); }
      }

      /* Cell set row */
      .ge-wrap .ge-set-row {
        padding:8px 10px; border-radius: var(--ge-radius-sm); margin-bottom:5px;
        background: var(--ge-bg-soft); display:flex; align-items:center;
        justify-content:space-between; border:1px solid transparent;
        transition: background 150ms ease, border-color 150ms ease;
      }
      .ge-wrap .ge-set-row:hover {
        background:#eef4fb; border-color:#dbe9fa;
      }
      .ge-wrap .ge-set-row.active {
        background: linear-gradient(90deg,#dbe9fa 0%, #eaf5fb 100%);
        border-color:#bcdaf0; box-shadow: inset 3px 0 0 var(--ge-accent);
      }
      .ge-wrap .ge-set-left { display:flex; align-items:center; min-width:0; }
      .ge-wrap .ge-set-swatch {
        width:12px; height:12px; border-radius:50%;
        display:inline-block; margin-right:10px; flex-shrink:0;
        box-shadow: inset 0 0 0 1px rgba(0,0,0,.1);
      }
      .ge-wrap .ge-set-name {
        font-weight:600; font-size:13px; color:#111827;
        overflow:hidden; text-overflow:ellipsis; white-space:nowrap;
      }
      .ge-wrap .ge-set-count {
        color:#6b7280; font-size:11px; margin-left:8px;
        background:#fff; border:1px solid var(--ge-border);
        padding:1px 6px; border-radius:999px; font-weight:500;
      }
      .ge-wrap .ge-set-actions { display:flex; align-items:center; gap:2px; flex-shrink:0; }
      .ge-wrap .ge-set-actions a {
        width:26px; height:26px; display:inline-flex;
        align-items:center; justify-content:center;
        color:#6b7280; font-size:12px; border-radius:6px;
        cursor:pointer; transition: background 150ms ease, color 150ms ease;
      }
      .ge-wrap .ge-set-actions a:hover { background:#fff; color:#111827; text-decoration:none; }
      .ge-wrap .ge-set-actions a.ge-rm:hover { background:#fef2f2; color: var(--ge-danger); }
      .ge-wrap .ge-set-actions a.ge-hide.is-hidden { opacity:.45; }

      /* Annotate rows — vertical stack so the input has full width */
      .ge-wrap .ge-annot-row {
        padding:10px 0 12px;
        border-bottom:1px solid #f1f5f9;
      }
      .ge-wrap .ge-annot-row:last-child { border-bottom:none; padding-bottom:4px; }
      .ge-wrap .ge-annot-head {
        display:flex; align-items:center; gap:8px; margin-bottom:6px;
      }
      .ge-wrap .ge-annot-name {
        font-weight:600; font-size:12.5px; color:#111827;
        overflow:hidden; text-overflow:ellipsis; white-space:nowrap;
        flex:1;
      }
      .ge-wrap .ge-annot-row .form-control,
      .ge-wrap .ge-annot-row input[type='text'] {
        height:34px; font-size:13px;
        border:1px solid #d1d5db; border-radius:6px;
        background:#fff;
      }
      .ge-wrap .ge-annot-row .form-control:focus,
      .ge-wrap .ge-annot-row input[type='text']:focus {
        border-color: var(--ge-accent);
        box-shadow: 0 0 0 3px rgba(77,187,213,0.15);
        outline:none;
      }
      .ge-wrap .ge-annot-empty {
        display:flex; align-items:center; gap:8px;
        padding:14px 12px; margin:6px 0;
        background: var(--ge-bg-soft);
        border:1px dashed var(--ge-border-strong);
        border-radius: var(--ge-radius-sm);
        color: var(--ge-text-muted); font-size:12px;
      }
      .ge-wrap .ge-annot-empty .fa { color: var(--ge-accent); }

      /* Small badge / info pill */
      .ge-wrap .ge-badge {
        display:inline-block; padding:2px 8px; border-radius:999px;
        background:#eef4fb; color: var(--ge-accent-deep);
        font-size:11px; font-weight:600;
      }
      .ge-wrap .ge-info-pill {
        background:#eef4fb; color: var(--ge-accent-deep);
        border-radius:999px; padding:3px 10px;
        font-size:11.5px; font-weight:600;
        display:inline-flex; align-items:center; gap:5px;
      }
      .ge-wrap .ge-section-hint {
        font-size:11.5px; color: var(--ge-text-muted);
        margin-top:8px; line-height:1.4;
      }

      /* Toolbars */
      .ge-wrap .ge-toolbar {
        display:flex; align-items:center; gap:6px; flex-wrap:wrap;
      }
      .ge-wrap .ge-toolbar .btn { padding:5px 10px; font-size:12px; }

      /* Tabs polish */
      .ge-wrap .nav-tabs {
        border-bottom:1px solid var(--ge-border); margin-bottom:10px;
      }
      .ge-wrap .nav-tabs > li > a {
        padding:8px 12px; font-size:12.5px; font-weight:600;
        color:var(--ge-text-muted); border:none; border-radius:0;
        cursor:pointer;
      }
      .ge-wrap .nav-tabs > li > a:hover { background:transparent; color: var(--ge-text); }
      .ge-wrap .nav-tabs > li.active > a,
      .ge-wrap .nav-tabs > li.active > a:hover,
      .ge-wrap .nav-tabs > li.active > a:focus {
        color: var(--ge-accent-deep); background:transparent;
        border:none; border-bottom:2px solid var(--ge-accent);
      }

      /* Action button polish — cursor-pointer everywhere */
      .ge-wrap .btn { cursor:pointer; transition: all 150ms ease; }
      .ge-wrap .btn-primary { background: var(--ge-accent); border-color: var(--ge-accent); }
      .ge-wrap .btn-primary:hover { background: var(--ge-accent-deep); border-color: var(--ge-accent-deep); }
      .ge-wrap .btn-success { background: var(--ge-success); border-color: var(--ge-success); }
      .ge-wrap .btn-danger  { background: var(--ge-danger);  border-color: var(--ge-danger); }

      /* Row-of-actions wrapper with equal-height alignment */
      .ge-wrap .ge-row-actions { display:flex; gap:8px; align-items:stretch; }
      .ge-wrap .ge-row-actions > * { flex:1; }
      .ge-wrap .ge-row-actions .btn { width:100%; }

      /* DE heatmap controls — one clean grid */
      .ge-wrap .ge-de-controls {
        display:grid;
        grid-template-columns: repeat(6, minmax(0,1fr));
        gap:10px; align-items:end;
      }
      @media (max-width: 1100px) {
        .ge-wrap .ge-de-controls { grid-template-columns: repeat(3, minmax(0,1fr)); }
      }
      @media (max-width: 600px) {
        .ge-wrap .ge-de-controls { grid-template-columns: repeat(2, minmax(0,1fr)); }
      }
      .ge-wrap .ge-de-controls .btn { width:100%; }

      /* Export buttons stack */
      .ge-wrap .ge-stack { display:flex; flex-direction:column; gap:8px; }
      .ge-wrap .ge-stack .btn,
      .ge-wrap .ge-stack .shiny-download-link { width:100%; }

      /* Inline radio group visual */
      .ge-wrap .shiny-options-group { margin-top:2px; }
      .ge-wrap .radio-inline { padding-left:22px; margin-right:10px; font-size:12px; }

      /* DT tables inside gene-explorer panels — keep pagination visible */
      .ge-wrap .dataTables_wrapper {
        font-size:12px;
        display:flex; flex-direction:column;
      }
      .ge-wrap .dataTables_wrapper .dataTables_filter {
        margin-bottom:6px; text-align:left;
      }
      .ge-wrap .dataTables_wrapper .dataTables_filter input {
        width:100%; margin-left:0; border:1px solid var(--ge-border);
        border-radius:6px; padding:5px 8px; font-size:12px;
      }
      .ge-wrap .dataTables_wrapper table.dataTable {
        font-size:12px; margin:0 !important;
      }
      .ge-wrap .dataTables_wrapper table.dataTable thead th {
        padding:6px 8px; font-size:11.5px;
        font-weight:700; color:#374151;
        background:#f9fafb;
      }
      .ge-wrap .dataTables_wrapper table.dataTable tbody td {
        padding:5px 8px;
      }
      .ge-wrap .dataTables_wrapper .dataTables_info,
      .ge-wrap .dataTables_wrapper .dataTables_paginate {
        padding:8px 0 2px !important; font-size:11.5px;
      }
      .ge-wrap .dataTables_wrapper .dataTables_paginate .paginate_button {
        padding:3px 8px; margin:0 1px;
      }
      .ge-wrap .dataTables_scrollBody { border-bottom:1px solid var(--ge-border); }

      /* Right-panel tab content grows — table gets remaining space, pagination below */
      .ge-wrap .ge-right-panel .tab-content,
      .ge-wrap .ge-right-panel .tab-pane { overflow:visible; }
      .ge-wrap .ge-right-panel .dataTables_wrapper { overflow:visible; }

      /* Panels — equal heights for the three top panels */
      .ge-wrap .ge-top-row > [class*='col-'] { display:flex; }
      .ge-wrap .ge-top-row .ge-panel { width:100%; display:flex; flex-direction:column; }
    ")),

    div(class = "ge-wrap",

      # ───────── Hero header ─────────────────────────────────────────
      div(class = "ge-hero",
        div(class = "ge-hero-left",
          div(class = "ge-hero-icon", icon("compass")),
          div(
            div(class = "ge-hero-title", "Gene Explorer"),
            div(class = "ge-hero-sub",
                "Interactive lasso selection, cell-set annotation and on-the-fly differential expression.")
          )
        ),
        div(class = "ge-steps",
          tags$span(class = "ge-step", tags$span(class = "num", "1"), "Lasso"),
          tags$span(class = "ge-step-sep", HTML("&rarr;")),
          tags$span(class = "ge-step", tags$span(class = "num", "2"), "Cell set"),
          tags$span(class = "ge-step-sep", HTML("&rarr;")),
          tags$span(class = "ge-step", tags$span(class = "num", "3"), "Annotate"),
          tags$span(class = "ge-step-sep", HTML("&rarr;")),
          tags$span(class = "ge-step", tags$span(class = "num", "4"), "Save")
        )
      ),

      div(class = "ge-top-row",
      fluidRow(
        # ================== LEFT: UMAP =====================
        column(
          width = 6,
          div(class = "ge-panel",
            div(class = "ge-panel-header",
              div(class = "ge-panel-title", icon("map"), " UMAP"),
              div(class = "ge-info-pill",
                  icon("mouse-pointer"),
                  textOutput(ns("selection_info"), inline = TRUE))
            ),
            div(class = "ge-controls c-4",
              div(class = "ge-field",
                tags$label(`for` = ns("color_mode"), "Color by"),
                selectInput(ns("color_mode"), NULL,
                            choices  = c("Cluster"         = "cluster",
                                         "Gene expression" = "gene",
                                         "Module score"    = "module",
                                         "Cell sets"       = "sets",
                                         "User annotation" = "user_annot"),
                            selected = "cluster", width = "100%")
              ),
              div(class = "ge-field",
                tags$label(`for` = ns("group_by"), "Group / idents"),
                selectInput(ns("group_by"), NULL, choices = NULL, width = "100%")
              ),
              div(class = "ge-field",
                tags$label(`for` = ns("max_points"), "Max points"),
                numericInput(ns("max_points"), NULL, value = 50000L,
                             min = 5000L, max = 300000L, step = 5000L,
                             width = "100%")
              ),
              div(class = "ge-field",
                tags$label(`for` = ns("point_size"), "Point size"),
                selectInput(ns("point_size"), NULL,
                            choices  = c("Small" = 3, "Medium" = 5, "Large" = 7),
                            selected = 5, width = "100%")
              )
            ),
            plotlyOutput(ns("umap"), height = "520px"),
            tags$div(style = "height:12px;"),
            div(class = "ge-field",
              tags$label(`for` = ns("set_name"), "Name the current lasso selection"),
              div(class = "ge-row-actions",
                div(style = "flex:3;",
                    textInput(ns("set_name"), NULL,
                              placeholder = "e.g., cluster_of_interest",
                              width = "100%")),
                div(style = "flex:1.2;",
                    actionButton(ns("save_set"), "Add cell set",
                                 icon = icon("plus"), class = "btn-primary",
                                 width = "100%")),
                div(style = "flex:0.8;",
                    actionButton(ns("clear_lasso"), "Clear",
                                 icon = icon("eraser"), class = "btn-default",
                                 width = "100%"))
              )
            )
          )
        ),

        # ================== MIDDLE: Cell sets / Annotate ====
        column(
          width = 3,
          div(class = "ge-panel",
            div(class = "ge-panel-header",
              div(class = "ge-panel-title", icon("layer-group"), " Cell sets")
            ),
            tabsetPanel(
              id = ns("mid_tabs"), type = "tabs",

              tabPanel("Cell sets",
                tags$div(style = "height:10px;"),
                uiOutput(ns("cell_sets_ui")),
                tags$hr(style = "margin:10px 0;"),
                div(class = "ge-toolbar",
                  actionButton(ns("show_all_sets"), "Show all",
                               icon = icon("eye"), class = "btn-default btn-sm"),
                  actionButton(ns("hide_all_sets"), "Hide all",
                               icon = icon("eye-slash"), class = "btn-default btn-sm"),
                  actionButton(ns("clear_sets"), "Delete all",
                               icon = icon("trash"), class = "btn-danger btn-sm")
                ),
                div(class = "ge-section-hint",
                    "Lasso on the UMAP, name the selection, then click 'Add cell set'.")
              ),

              tabPanel("Annotate",
                tags$div(style = "height:10px;"),
                uiOutput(ns("annot_rows_ui")),
                tags$hr(style = "margin:10px 0;"),
                div(class = "ge-field",
                  tags$label(`for` = ns("base_annot_col"), "Base labels"),
                  selectizeInput(ns("base_annot_col"), NULL,
                                 choices = NULL, width = "100%",
                                 options = list(dropdownParent = "body"))
                ),
                actionButton(ns("apply_annot"),
                             "Apply annotations",
                             icon  = icon("check-circle"),
                             class = "btn-success",
                             width = "100%"),
                div(class = "ge-section-hint",
                    paste0("Cells outside any cell set keep their 'Base label'. ",
                           "Result is written to metadata column '",
                           .ge_annot_col, "'."))
              ),

              tabPanel("Export",
                tags$div(style = "height:10px;"),
                div(class = "ge-stack",
                  actionButton(ns("commit_to_meta"),
                               "Commit cell sets as metadata",
                               icon = icon("save"), class = "btn-primary"),
                  actionButton(ns("promote_set"),
                               "Promote active set to Idents()",
                               icon = icon("flag"), class = "btn-default"),
                  downloadButton(ns("dl_sets_csv"),
                                 "Download cell sets (CSV)",
                                 class = "btn-default")
                ),
                div(class = "ge-section-hint",
                    paste0("Commit writes logical columns '", .ge_prefix, "*' ",
                           "plus a consolidated factor '",
                           .ge_consolidated_col, "' so DE / composition / ",
                           "pathway modules can use them."))
              )
            )
          )
        ),

        # ================== RIGHT: Genes / DE ================
        column(
          width = 3,
          div(class = "ge-panel ge-right-panel",
            div(class = "ge-panel-header",
              div(class = "ge-panel-title", icon("dna"), " Genes")
            ),
            tabsetPanel(
              id = ns("right_tabs"), type = "tabs",

              tabPanel("Gene list",
                tags$div(style = "height:10px;"),
                div(class = "ge-field",
                  tags$label(`for` = ns("genes"), "Selected genes"),
                  selectizeInput(ns("genes"), NULL, choices = NULL,
                                 multiple = TRUE, width = "100%",
                                 options = list(
                                   placeholder = "Type to search genes (or click a row below)...",
                                   plugins = list("remove_button"),
                                   maxItems = 50L,
                                   closeAfterSelect = TRUE
                                 ))
                ),
                div(class = "ge-controls c-2",
                  div(class = "ge-field",
                    tags$label("Sort by"),
                    radioButtons(ns("gene_sort"), NULL,
                                 choices  = c("Dispersion"  = "dispersion",
                                              "Mean expr"   = "mean",
                                              "A \u2192 Z"  = "alpha"),
                                 selected = "dispersion", inline = TRUE)
                  ),
                  div(class = "ge-field",
                    tags$label("Filter"),
                    checkboxInput(ns("only_hvg"), "HVG only", value = FALSE)
                  )
                ),
                DT::DTOutput(ns("gene_table")),
                div(class = "ge-section-hint",
                    textOutput(ns("gene_status"), inline = TRUE))
              ),

              tabPanel("Differential expression",
                tags$div(style = "height:10px;"),
                div(style = "font-size:12px; color:#374151; margin-bottom:8px;",
                    textOutput(ns("de_context"), inline = TRUE)),
                DT::DTOutput(ns("de_table"), height = "440px")
              )
            )
          )
        )
      )),

      tags$div(style = "height:14px;"),

      # ================== BOTTOM: DE heatmap =================
      fluidRow(
        column(12,
          div(class = "ge-panel",
            div(class = "ge-panel-header",
              div(class = "ge-panel-title", icon("th"),
                  " Heatmap — top DE genes"),
              div(class = "ge-info-pill",
                  icon("wave-square"),
                  textOutput(ns("de_status"), inline = TRUE))
            ),
            div(class = "ge-de-controls",
              div(class = "ge-field",
                tags$label(`for` = ns("de_top_n"), "Top N"),
                numericInput(ns("de_top_n"), NULL, value = 20L,
                             min = 5L, max = 100L, step = 5L, width = "100%")
              ),
              div(class = "ge-field",
                tags$label(`for` = ns("de_scale"), "Scale"),
                selectInput(ns("de_scale"), NULL,
                            choices  = c("Row z-score" = "row", "None" = "none"),
                            selected = "row", width = "100%")
              ),
              div(class = "ge-field",
                tags$label(`for` = ns("de_max_cells"), "Max cells / ident"),
                numericInput(ns("de_max_cells"), NULL, value = 1500L,
                             min = 200L, max = 5000L, step = 100L, width = "100%")
              ),
              div(class = "ge-field",
                tags$label(`for` = ns("de_mode"), "Comparison"),
                selectInput(ns("de_mode"), NULL,
                            choices  = c("Selection vs rest" = "vs_rest",
                                         "All cell sets"     = "all_sets"),
                            selected = "vs_rest", width = "100%")
              ),
              div(class = "ge-field",
                tags$label(HTML("&nbsp;")),
                actionButton(ns("compute_de"), "Compute",
                             icon = icon("play"), class = "btn-primary")
              ),
              div(class = "ge-field",
                tags$label(HTML("&nbsp;")),
                downloadButton(ns("dl_de_csv"), "Download DE",
                               class = "btn-default")
              )
            ),
            tags$hr(style = "margin:12px 0;"),
            plotOutput(ns("de_heatmap"), height = "520px")
          )
        )
      )
    )
  )
}

# =============================================================================
# Server
# =============================================================================

geneExplorerServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

    rv <- reactiveValues(
      cell_sets        = list(),   # name -> character vector of cell ids
      set_colors       = character(0),
      set_visible      = logical(0),
      set_labels       = character(0),   # name -> user-typed label
      active_set       = NULL,
      lasso_cells      = character(0),
      de_markers       = NULL,
      de_context_text  = "Not computed yet.",
      gene_status_text = "",
      feat_index       = NULL,
      last_obj_sig     = NULL
    )

    has_obj <- reactive({
      !is.null(values$active_dataset) && inherits(values$active_dataset, "Seurat")
    })

    # Build / refresh feature index when the object changes
    observe({
      req(has_obj())
      obj <- values$active_dataset
      sig <- paste(ncol(obj), nrow(obj), paste(names(obj@assays), collapse = ","), sep = ":")
      if (!identical(rv$last_obj_sig, sig)) {
        rv$feat_index <- .ge_build_feature_index(obj)
        rv$last_obj_sig <- sig
      }
    })

    # --- group_by choices ---------------------------------------------------
    observe({
      req(has_obj())
      meta <- values$active_dataset@meta.data
      cand <- names(meta)[vapply(meta, function(x)
        is.factor(x) || is.character(x) ||
        (is.numeric(x) && length(unique(x)) < 50),
        logical(1))]
      preferred <- intersect(
        c("seurat_clusters", "cell_type_consensus", "cell_type",
          "predicted.celltype", "ident", .ge_annot_col),
        cand)
      sel <- if (length(preferred)) preferred[1] else cand[1]
      updateSelectInput(session, "group_by", choices = cand,
                        selected = isolate(input$group_by) %||% sel)
      updateSelectizeInput(session, "base_annot_col",
                           choices = c("(none)" = "", cand),
                           selected = isolate(input$base_annot_col) %||% sel)
    })

    # --- server-side gene selectize (all features across all assays) --------
    observe({
      req(has_obj(), !is.null(rv$feat_index))
      # show canonical gene names from the default assay first, then add others
      obj <- values$active_dataset
      da <- .ge_safe_default_assay(obj)
      genes <- rownames(obj[[da]])
      others <- setdiff(
        unlist(lapply(rv$feat_index, function(x) x$gene), use.names = FALSE),
        genes)
      all_genes <- c(genes, others)
      updateSelectizeInput(session, "genes",
                          choices = all_genes, selected = character(0),
                          server = TRUE)
    })

    # --- gene stats table (cached per object signature) --------------------
    gene_stats <- reactive({
      req(has_obj())
      obj <- values$active_dataset
      assay <- .ge_safe_default_assay(obj)
      genes <- rownames(obj[[assay]])

      disp <- rep(NA_real_, length(genes))
      mean_expr <- rep(NA_real_, length(genes))
      is_hvg <- rep(FALSE, length(genes))

      meta_feat <- tryCatch(obj[[assay]]@meta.features,
                            error = function(e) NULL)
      if (!is.null(meta_feat) && nrow(meta_feat) == length(genes)) {
        if ("vst.variance.standardized" %in% colnames(meta_feat)) {
          disp <- as.numeric(meta_feat$vst.variance.standardized)
        } else if ("mvp.dispersion.scaled" %in% colnames(meta_feat)) {
          disp <- as.numeric(meta_feat$mvp.dispersion.scaled)
        }
        if ("vst.mean" %in% colnames(meta_feat)) {
          mean_expr <- as.numeric(meta_feat$vst.mean)
        } else if ("mean" %in% colnames(meta_feat)) {
          mean_expr <- as.numeric(meta_feat$mean)
        }
        if ("vst.variable" %in% colnames(meta_feat)) {
          is_hvg <- as.logical(meta_feat$vst.variable)
        }
      }

      hvg <- tryCatch(VariableFeatures(obj), error = function(e) character(0))
      if (length(hvg)) is_hvg <- is_hvg | (genes %in% hvg)

      data.frame(Gene = genes,
                 Dispersion = round(disp, 3),
                 MeanExpr   = round(mean_expr, 3),
                 HVG        = is_hvg,
                 stringsAsFactors = FALSE)
    })

    gene_stats_sorted <- reactive({
      df <- gene_stats()
      if (isTRUE(input$only_hvg)) df <- df[df$HVG, , drop = FALSE]
      sort_by <- input$gene_sort %||% "dispersion"
      if (sort_by == "dispersion") {
        df <- df[order(-df$Dispersion, na.last = TRUE), , drop = FALSE]
      } else if (sort_by == "mean") {
        df <- df[order(-df$MeanExpr, na.last = TRUE), , drop = FALSE]
      } else {
        df <- df[order(df$Gene), , drop = FALSE]
      }
      df
    })

    output$gene_table <- DT::renderDT({
      df <- gene_stats_sorted()
      DT::datatable(
        df[, c("Gene", "Dispersion", "MeanExpr")],
        rownames = FALSE,
        selection = "single",
        escape = TRUE,
        options = list(
          pageLength      = 10L,
          lengthChange    = FALSE,
          scrollY         = "300px",
          scrollCollapse  = TRUE,
          scrollX         = FALSE,
          dom             = "ftip",
          searchHighlight = TRUE,
          columnDefs      = list(list(className = "dt-center", targets = 1:2))
        )
      )
    }, server = TRUE)

    observeEvent(input$gene_table_rows_selected, {
      req(input$gene_table_rows_selected)
      df <- gene_stats_sorted()
      if (nrow(df) < input$gene_table_rows_selected) return()
      g <- df$Gene[input$gene_table_rows_selected]
      cur <- input$genes %||% character(0)
      new <- unique(c(cur, g))
      updateSelectizeInput(session, "genes", selected = new)
      if (length(new) == 1L)
        updateSelectInput(session, "color_mode", selected = "gene")
      else if (length(new) >= 2L)
        updateSelectInput(session, "color_mode", selected = "module")
    })

    # --- UMAP data frame (downsampled) --------------------------------------
    umap_df <- reactive({
      req(has_obj())
      obj <- values$active_dataset
      red <- .ge_best_reduction(obj)
      req(!is.null(red))
      emb <- Embeddings(obj, reduction = red)[, 1:2, drop = FALSE]
      colnames(emb) <- c("Dim1", "Dim2")
      df <- data.frame(
        Cell = rownames(emb),
        Dim1 = emb[, 1],
        Dim2 = emb[, 2],
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      gb <- input$group_by
      if (!is.null(gb) && gb %in% colnames(obj@meta.data)) {
        df$Group <- as.character(obj@meta.data[[gb]])
      } else df$Group <- "all"
      max_n <- max(5000L, as.integer(input$max_points %||% 50000L))
      if (nrow(df) > max_n) {
        set.seed(42)
        idx <- sample.int(nrow(df), max_n)
        df <- df[sort(idx), , drop = FALSE]
      }
      df
    })

    # --- color / overlay ----------------------------------------------------
    color_vec <- reactive({
      req(has_obj())
      df <- umap_df(); req(nrow(df) > 0)
      mode <- input$color_mode %||% "cluster"
      genes <- input$genes %||% character(0)
      obj <- values$active_dataset
      assay <- .ge_safe_default_assay(obj)

      # Cell-set overlay
      if (mode == "sets") {
        if (length(rv$cell_sets) == 0L)
          return(list(kind = "cat", values = rep("(no sets)", nrow(df)),
                      title = "Cell sets", palette = c(`(no sets)` = "#bdc3c7")))
        # last-write-wins for overlapping sets
        label <- rep("other", nrow(df))
        names(label) <- df$Cell
        for (nm in names(rv$cell_sets)) {
          if (isFALSE(rv$set_visible[[nm]])) next
          hit <- intersect(rv$cell_sets[[nm]], df$Cell)
          if (length(hit)) label[hit] <- nm
        }
        pal <- c(rv$set_colors, other = "#dfe4ea")
        return(list(kind = "cat", values = unname(label),
                    title = "Cell sets", palette = pal))
      }

      if (mode == "user_annot" && .ge_annot_col %in% colnames(obj@meta.data)) {
        vals <- as.character(obj@meta.data[[.ge_annot_col]])
        names(vals) <- colnames(obj)
        out <- vals[df$Cell]; out[is.na(out)] <- "unassigned"
        return(list(kind = "cat", values = unname(out),
                    title = .ge_annot_col, palette = NULL))
      }

      if (mode == "cluster" || length(genes) == 0L) {
        return(list(kind = "cat", values = df$Group,
                    title = input$group_by %||% "group", palette = NULL))
      }

      # gene / module
      fi <- rv$feat_index
      present <- character(0); missing_g <- character(0)
      for (g in genes) {
        info <- .ge_resolve_gene(fi, g, prefer_assay = assay)
        if (is.null(info)) missing_g <- c(missing_g, g)
        else present <- c(present, info$gene)
      }
      present <- unique(present)

      if (length(missing_g)) {
        rv$gene_status_text <- paste0(
          "Missing genes: ", paste(missing_g, collapse = ", "),
          " (not found in any assay)."
        )
      } else if (length(present)) {
        rv$gene_status_text <- paste0(
          "Using ", length(present), " gene(s) on assay '", assay, "'.")
      } else {
        rv$gene_status_text <- ""
      }

      if (!length(present)) {
        return(list(kind = "cat", values = df$Group,
                    title = input$group_by %||% "group", palette = NULL))
      }

      # single-gene expression (fast path)
      if (mode == "gene" || (mode == "module" && length(present) == 1L)) {
        expr <- .ge_get_expr(obj, fi, present[1], df$Cell, assay)
        if (is.null(expr)) {
          return(list(kind = "cat", values = df$Group,
                      title = input$group_by %||% "group", palette = NULL))
        }
        out <- rep(NA_real_, nrow(df))
        out[match(names(expr$values), df$Cell)] <- as.numeric(expr$values)
        out <- .ge_clip(out, 0, 0.99)
        return(list(kind = "num", values = out,
                    title = paste0("expr: ", present[1])))
      }

      # module score: compute mean of z-scored rows over selected cells
      # (fast, no Seurat::AddModuleScore call needed)
      mat <- NULL
      for (g in present) {
        e <- .ge_get_expr(obj, fi, g, df$Cell, assay)
        if (is.null(e)) next
        v <- rep(NA_real_, nrow(df))
        v[match(names(e$values), df$Cell)] <- as.numeric(e$values)
        mat <- rbind(mat, v)
      }
      if (is.null(mat) || !nrow(mat)) {
        return(list(kind = "cat", values = df$Group,
                    title = input$group_by %||% "group", palette = NULL))
      }
      mu <- rowMeans(mat, na.rm = TRUE)
      sdv <- apply(mat, 1, stats::sd, na.rm = TRUE); sdv[sdv == 0 | is.na(sdv)] <- 1
      z <- sweep(sweep(mat, 1, mu, "-"), 1, sdv, "/")
      score <- colMeans(z, na.rm = TRUE)
      score <- .ge_clip(score, 0.01, 0.99)
      list(kind = "num", values = as.numeric(score),
           title = paste0("module (", length(present), " genes)"))
    })

    output$gene_status <- renderText({ rv$gene_status_text })

    # --- UMAP plot ----------------------------------------------------------
    output$umap <- renderPlotly({
      req(has_obj())
      df <- umap_df(); req(nrow(df) > 0)
      cv <- color_vec()
      use_webgl <- nrow(df) > 10000L
      scatter_type <- if (use_webgl) "scattergl" else "scatter"
      psize <- as.numeric(input$point_size %||% 4)

      plt <- if (cv$kind == "cat") {
        grps <- unique(cv$values)
        # put "other" last
        if ("other" %in% grps) grps <- c(setdiff(grps, "other"), "other")
        if (!is.null(cv$palette)) {
          pal <- cv$palette
          miss <- setdiff(grps, names(pal))
          if (length(miss)) {
            extra <- grDevices::hcl.colors(length(miss), "Dynamic")
            names(extra) <- miss
            pal <- c(pal, extra)
          }
        } else {
          pal_base <- if (exists("get_enhanced_color_palette", mode = "function"))
                        get_enhanced_color_palette(max(3L, length(grps)),
                                                    "qualitative")
                      else grDevices::hcl.colors(max(3L, length(grps)), "Dynamic")
          pal <- setNames(pal_base[seq_along(grps)], grps)
        }
        p <- plot_ly(type = scatter_type, mode = "markers",
                     source = ns("umap_src"))
        for (g in grps) {
          sub <- df[cv$values == g, , drop = FALSE]
          if (!nrow(sub)) next
          p <- add_trace(p, x = sub$Dim1, y = sub$Dim2,
                         type = scatter_type, mode = "markers",
                         name = g, text = sub$Cell, customdata = sub$Cell,
                         marker = list(size = psize,
                                       opacity = if (identical(g, "other")) 0.25 else 0.8,
                                       color = pal[[g]],
                                       line = list(width = 0)),
                         hovertemplate = paste0("Cell: %{text}<br>",
                                                cv$title, ": %{data.name}<extra></extra>"))
        }
        p
      } else {
        plot_ly(data = df,
                x = ~Dim1, y = ~Dim2,
                type = scatter_type, mode = "markers",
                text = ~Cell, customdata = ~Cell,
                source = ns("umap_src"),
                marker = list(size = psize, opacity = 0.85,
                              color = cv$values,
                              colorscale = "Viridis",
                              showscale = TRUE,
                              colorbar = list(title = cv$title,
                                              thickness = 10, len = 0.6),
                              line = list(width = 0)),
                hovertemplate = "Cell: %{text}<br>val: %{marker.color:.2f}<extra></extra>")
      }

      plt %>%
        layout(title = list(text = paste0("UMAP \u2014 ", cv$title),
                            font = list(size = 13)),
               dragmode = "lasso",
               xaxis = list(title = "UMAP_1", zeroline = FALSE,
                            showgrid = FALSE),
               yaxis = list(title = "UMAP_2", zeroline = FALSE,
                            showgrid = FALSE),
               plot_bgcolor = "#fafbfc", paper_bgcolor = "#ffffff",
               legend = list(itemsizing = "constant", font = list(size = 11))) %>%
        config(displaylogo = FALSE,
               modeBarButtonsToAdd = list("lasso2d", "select2d"),
               toImageButtonOptions = list(format = "png", scale = 2)) %>%
        event_register("plotly_selected") %>%
        event_register("plotly_deselect")
    })

    # --- lasso selection ----------------------------------------------------
    observe({
      sel <- event_data("plotly_selected", source = ns("umap_src"))
      # Lassoing an empty area returns a list (not a data.frame), so
      # nrow(sel) is NULL and `!nrow(sel)` would error with "invalid
      # argument type". NROW() handles vectors, lists, and data.frames.
      if (is.null(sel) || NROW(sel) == 0L) return()
      cd <- sel$customdata
      if (is.null(cd)) return()
      cells <- as.character(cd)
      cells <- cells[!is.na(cells) & nzchar(cells)]
      rv$lasso_cells <- unique(cells)
    })

    observe({
      ev <- event_data("plotly_deselect", source = ns("umap_src"))
      if (!is.null(ev)) rv$lasso_cells <- character(0)
    })

    observeEvent(input$clear_lasso, {
      rv$lasso_cells <- character(0)
      showNotification("Cleared lasso selection.", type = "message",
                       duration = 2)
    })

    output$selection_info <- renderText({
      n <- length(rv$lasso_cells); ns_total <- length(rv$cell_sets)
      paste0(n, " cells lassoed  \u2022  ", ns_total, " cell set(s)")
    })

    # --- save cell set ------------------------------------------------------
    observeEvent(input$save_set, {
      nm <- trimws(input$set_name %||% "")
      if (!nzchar(nm)) {
        showNotification("Please name the cell set first.",
                         type = "warning"); return()
      }
      if (length(rv$lasso_cells) == 0L) {
        showNotification("Lasso cells on the UMAP first.",
                         type = "warning"); return()
      }
      nm <- gsub("[^A-Za-z0-9_.-]+", "_", nm)
      if (nm %in% names(rv$cell_sets))
        nm <- paste0(nm, "_", format(Sys.time(), "%H%M%S"))

      sets <- rv$cell_sets
      sets[[nm]] <- rv$lasso_cells
      rv$cell_sets <- sets

      # color
      n <- length(sets)
      col <- .ge_set_palette[((n - 1L) %% length(.ge_set_palette)) + 1L]
      rv$set_colors[[nm]] <- col
      rv$set_visible[[nm]] <- TRUE
      rv$set_labels[[nm]] <- ""
      rv$active_set <- nm

      updateTextInput(session, "set_name", value = "")
      showNotification(paste0("Saved cell set '", nm, "' (",
                              length(rv$lasso_cells), " cells)."),
                       type = "message", duration = 3)
    })

    # --- cell sets list UI --------------------------------------------------
    output$cell_sets_ui <- renderUI({
      sets <- rv$cell_sets
      if (length(sets) == 0L) {
        return(div(style = "color:#95a5a6; font-size:12px; font-style:italic; padding:8px 0;",
                   "No cell sets yet. Lasso cells, give them a name, and click 'Add cell set'."))
      }
      tagList(lapply(names(sets), function(nm) {
        sw <- rv$set_colors[[nm]] %||% "#888"
        visible <- isTRUE(rv$set_visible[[nm]])
        is_active <- identical(rv$active_set, nm)
        div(class = paste("ge-set-row", if (is_active) "active" else ""),
          div(class = "ge-set-left",
            tags$span(class = "ge-set-swatch",
                      style = paste0("background:", sw, ";")),
            actionLink(ns(paste0("pick_", nm)), nm,
                       class = "ge-set-name",
                       style = if (!visible) "opacity:0.45;" else ""),
            tags$span(class = "ge-set-count", length(sets[[nm]]))
          ),
          div(class = "ge-set-actions",
            actionLink(ns(paste0("vis_", nm)), label = NULL,
                       icon = icon(if (visible) "eye" else "eye-slash"),
                       class = paste("ge-hide", if (visible) "" else "is-hidden"),
                       title = "Toggle visibility"),
            actionLink(ns(paste0("rn_", nm)), label = NULL,
                       icon = icon("edit"), title = "Rename"),
            actionLink(ns(paste0("rm_", nm)), label = NULL,
                       icon = icon("times"), class = "ge-rm",
                       title = "Delete")
          )
        )
      }))
    })

    # --- annotate rows UI ---------------------------------------------------
    # IMPORTANT: only depend on the set NAMES (rv$cell_sets names + set_colors).
    # Reading rv$set_labels[[nm]] reactively would rebuild this UI on every
    # keystroke, destroying and recreating the textInput -> flicker / cursor jump.
    # Swatch colour and previous label are read via isolate().
    output$annot_rows_ui <- renderUI({
      nms <- names(rv$cell_sets)
      if (!length(nms)) {
        return(div(class = "ge-annot-empty",
                   icon("layer-group"),
                   tags$span("No cell sets yet. Create one first.")))
      }
      colors <- isolate(rv$set_colors)
      labels <- isolate(rv$set_labels)
      tagList(lapply(nms, function(nm) {
        sw <- colors[[nm]] %||% "#888"
        div(class = "ge-annot-row",
          div(class = "ge-annot-head",
            tags$span(class = "ge-set-swatch",
                      style = paste0("background:", sw, ";")),
            tags$span(class = "ge-annot-name", nm)
          ),
          textInput(ns(paste0("lbl_", nm)), NULL,
                    value = labels[[nm]] %||% "",
                    placeholder = "cell type label, e.g. 'B cell'",
                    width = "100%")
        )
      }))
    })

    # --- register per-set observers once ------------------------------------
    registered_obs <- reactiveValues(names = character(0))
    observe({
      sets <- rv$cell_sets
      new_names <- setdiff(names(sets), registered_obs$names)
      if (!length(new_names)) return()
      for (nm in new_names) {
        local({
          this <- nm

          observeEvent(input[[paste0("pick_", this)]], {
            if (!this %in% names(rv$cell_sets)) return()
            rv$active_set <- this
            rv$lasso_cells <- rv$cell_sets[[this]]
            showNotification(paste0("Active: ", this),
                             type = "message", duration = 1.5)
          }, ignoreInit = TRUE)

          observeEvent(input[[paste0("vis_", this)]], {
            if (!this %in% names(rv$cell_sets)) return()
            rv$set_visible[[this]] <- !isTRUE(rv$set_visible[[this]])
          }, ignoreInit = TRUE)

          observeEvent(input[[paste0("rn_", this)]], {
            if (!this %in% names(rv$cell_sets)) return()
            showModal(modalDialog(
              title = paste0("Rename '", this, "'"),
              textInput(ns("rename_to"), "New name", value = this),
              footer = tagList(
                modalButton("Cancel"),
                actionButton(ns("rename_confirm"), "Rename",
                             class = "btn-primary")
              ), size = "s", easyClose = TRUE
            ))
            rv$active_set <- this
          }, ignoreInit = TRUE)

          observeEvent(input[[paste0("rm_", this)]], {
            sets2 <- rv$cell_sets; sets2[[this]] <- NULL
            rv$cell_sets <- sets2
            rv$set_visible[[this]] <- NULL
            rv$set_colors[[this]] <- NULL
            rv$set_labels[[this]] <- NULL
            if (identical(rv$active_set, this)) rv$active_set <- NULL
            showNotification(paste0("Removed: ", this),
                             type = "message", duration = 1.5)
          }, ignoreInit = TRUE)

          observeEvent(input[[paste0("lbl_", this)]], {
            if (!this %in% names(rv$cell_sets)) return()
            rv$set_labels[[this]] <- input[[paste0("lbl_", this)]]
          }, ignoreInit = TRUE)
        })
      }
      registered_obs$names <- c(registered_obs$names, new_names)
    })

    observeEvent(input$rename_confirm, {
      nm_old <- rv$active_set
      nm_new <- trimws(input$rename_to %||% "")
      if (is.null(nm_old) || !nzchar(nm_new)) { removeModal(); return() }
      nm_new <- gsub("[^A-Za-z0-9_.-]+", "_", nm_new)
      if (nm_new %in% names(rv$cell_sets) && !identical(nm_new, nm_old)) {
        showNotification("Name already exists.", type = "warning")
        return()
      }
      sets <- rv$cell_sets
      nms <- names(sets); nms[nms == nm_old] <- nm_new
      names(sets) <- nms
      rv$cell_sets <- sets
      # migrate parallel maps
      for (fld in c("set_colors", "set_visible", "set_labels")) {
        vec <- rv[[fld]]
        if (nm_old %in% names(vec)) {
          vec[[nm_new]] <- vec[[nm_old]]; vec[[nm_old]] <- NULL
          rv[[fld]] <- vec
        }
      }
      if (identical(rv$active_set, nm_old)) rv$active_set <- nm_new
      removeModal()
      showNotification(paste0("Renamed to '", nm_new, "'."),
                       type = "message", duration = 2)
    })

    # --- show/hide all, delete all -----------------------------------------
    observeEvent(input$show_all_sets, {
      for (nm in names(rv$cell_sets)) rv$set_visible[[nm]] <- TRUE
    })
    observeEvent(input$hide_all_sets, {
      for (nm in names(rv$cell_sets)) rv$set_visible[[nm]] <- FALSE
    })
    observeEvent(input$clear_sets, {
      rv$cell_sets <- list(); rv$set_visible <- logical(0)
      rv$set_colors <- character(0); rv$set_labels <- character(0)
      rv$active_set <- NULL; rv$lasso_cells <- character(0)
      showNotification("Cleared all cell sets.", type = "message")
    })

    # --- commit sets to metadata -------------------------------------------
    observeEvent(input$commit_to_meta, {
      req(has_obj())
      sets <- rv$cell_sets
      if (!length(sets)) {
        showNotification("No cell sets to commit.", type = "warning"); return()
      }
      obj <- values$active_dataset
      all_cells <- colnames(obj)
      for (nm in names(sets)) {
        obj@meta.data[[paste0(.ge_prefix, nm)]] <- all_cells %in% sets[[nm]]
      }
      consol <- rep(NA_character_, length(all_cells))
      names(consol) <- all_cells
      for (nm in names(sets)) {
        hit <- intersect(sets[[nm]], all_cells); consol[hit] <- nm
      }
      consol[is.na(consol)] <- "other"
      obj@meta.data[[.ge_consolidated_col]] <- factor(
        consol, levels = c(names(sets), "other"))
      values$active_dataset <- obj
      showNotification(
        paste0("Committed ", length(sets), " set(s) as metadata columns '",
               .ge_prefix, "*' and '", .ge_consolidated_col, "'."),
        type = "message", duration = 5)
    })

    # --- apply annotations (final goal) ------------------------------------
    observeEvent(input$apply_annot, {
      req(has_obj())
      sets <- rv$cell_sets
      if (!length(sets)) {
        showNotification("No cell sets. Create some first.",
                         type = "warning"); return()
      }
      labels <- vapply(names(sets), function(nm) {
        v <- rv$set_labels[[nm]] %||% ""
        trimws(v)
      }, character(1))
      unlabeled <- names(sets)[!nzchar(labels)]
      if (length(unlabeled)) {
        showNotification(paste0("Please label all cell sets first. Missing: ",
                                paste(unlabeled, collapse = ", ")),
                         type = "warning", duration = 6)
        return()
      }

      obj <- values$active_dataset
      all_cells <- colnames(obj)

      base_col <- input$base_annot_col %||% ""
      if (nzchar(base_col) && base_col %in% colnames(obj@meta.data)) {
        out <- as.character(obj@meta.data[[base_col]])
      } else {
        out <- rep("unassigned", length(all_cells))
      }
      names(out) <- all_cells

      # apply set labels — last-write-wins
      for (nm in names(sets)) {
        hit <- intersect(sets[[nm]], all_cells)
        out[hit] <- labels[[nm]]
      }
      out[is.na(out) | !nzchar(out)] <- "unassigned"

      lvls <- c(unique(labels),
                setdiff(unique(out), c(unique(labels), "unassigned")),
                "unassigned")
      obj@meta.data[[.ge_annot_col]] <- factor(unname(out), levels = lvls)

      values$active_dataset <- obj

      # switch color mode to show the result
      updateSelectInput(session, "color_mode", selected = "user_annot")

      showNotification(
        paste0("Wrote '", .ge_annot_col, "' to metadata (",
               length(unique(labels)), " labels, base='",
               if (nzchar(base_col)) base_col else "unassigned", "')."),
        type = "message", duration = 5)
    })

    # --- promote -----------------------------------------------------------
    observeEvent(input$promote_set, {
      req(has_obj())
      nm <- rv$active_set
      if (is.null(nm) || !nm %in% names(rv$cell_sets)) {
        showNotification("Pick a cell set first.", type = "warning"); return()
      }
      obj <- values$active_dataset
      cells <- intersect(rv$cell_sets[[nm]], colnames(obj))
      if (!length(cells)) {
        showNotification("Empty set.", type = "warning"); return()
      }
      new_i <- as.character(Idents(obj)); names(new_i) <- colnames(obj)
      new_i[cells] <- nm
      Idents(obj) <- factor(new_i)
      values$active_dataset <- obj
      showNotification(paste0("Promoted '", nm, "' to Idents."),
                       type = "message")
    })

    # --- download sets CSV --------------------------------------------------
    output$dl_sets_csv <- downloadHandler(
      filename = function() paste0("gene_explorer_sets_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"),
                                   ".csv"),
      content = function(file) {
        sets <- rv$cell_sets
        if (!length(sets)) {
          utils::write.csv(data.frame(note = "No sets."), file,
                           row.names = FALSE); return()
        }
        rows <- do.call(rbind, lapply(names(sets), function(nm) {
          data.frame(set = nm,
                     label = rv$set_labels[[nm]] %||% "",
                     cell = sets[[nm]],
                     stringsAsFactors = FALSE)
        }))
        utils::write.csv(rows, file, row.names = FALSE)
      }
    )

    # --- DE computation -----------------------------------------------------
    de_result <- eventReactive(input$compute_de, {
      req(has_obj())

      obj <- .ge_prepare_for_de(values$active_dataset)
      all_cells <- colnames(obj)
      mode <- input$de_mode %||% "vs_rest"

      if (mode == "all_sets") {
        if (length(rv$cell_sets) < 2L) {
          showNotification("Need at least 2 cell sets for 'All cell sets' mode.",
                           type = "warning"); return(NULL)
        }
        grp <- rep(NA_character_, length(all_cells))
        names(grp) <- all_cells
        for (nm in names(rv$cell_sets)) {
          hit <- intersect(rv$cell_sets[[nm]], all_cells); grp[hit] <- nm
        }
        keep <- !is.na(grp)
        obj2 <- subset(obj, cells = all_cells[keep])
        Idents(obj2) <- factor(grp[keep], levels = names(rv$cell_sets))

        mk <- tryCatch(
          withProgress(message = "DE: all cell sets (one-vs-rest per set)...",
                       value = 0.1, {
            incProgress(0.3, detail = "FindAllMarkers...")
            m <- Seurat::FindAllMarkers(
              obj2, test.use = "wilcox", only.pos = TRUE,
              min.pct = 0.1, logfc.threshold = 0.25,
              max.cells.per.ident = as.integer(input$de_max_cells %||% 1500L),
              verbose = FALSE)
            incProgress(0.9, detail = "Done"); m
          }),
          error = function(e) {
            showNotification(paste("DE failed:", e$message), type = "error",
                             duration = 6); NULL
          })
        if (!is.null(mk) && nrow(mk) > 0)
          mk <- mk[order(mk$cluster, -mk$avg_log2FC), , drop = FALSE]
        rv$de_markers <- mk
        rv$de_context_text <- if (is.null(mk)) "DE failed."
                              else paste0("All cell sets: ", length(rv$cell_sets),
                                          " groups, ", nrow(mk), " markers.")
        return(list(mode = "all_sets", markers = mk, groups = grp[keep]))
      }

      # selection vs rest
      sel <- rv$lasso_cells
      if (!is.null(rv$active_set) && rv$active_set %in% names(rv$cell_sets)) {
        sel <- rv$cell_sets[[rv$active_set]]
      }
      if (length(sel) < 10L) {
        showNotification("Select \u226510 cells (lasso or pick a set).",
                         type = "warning"); return(NULL)
      }
      if (length(sel) > 5000L) {
        showNotification(paste0("Selection is large (", length(sel),
                                "); capping to 5000 for DE."),
                         type = "warning", duration = 4)
        set.seed(42); sel <- sample(sel, 5000L)
      }
      sel <- intersect(sel, all_cells); rest <- setdiff(all_cells, sel)
      if (!length(sel) || !length(rest)) {
        showNotification("Empty selection.", type = "error"); return(NULL)
      }

      res <- tryCatch(
        withProgress(message = "Computing DE (selection vs rest)...",
                     value = 0.1, {
          incProgress(0.2, detail = "FindMarkers...")
          m <- Seurat::FindMarkers(
            obj, ident.1 = sel, ident.2 = rest,
            test.use = "wilcox", only.pos = TRUE,
            min.pct = 0.1, logfc.threshold = 0.25,
            max.cells.per.ident = as.integer(input$de_max_cells %||% 1500L),
            verbose = FALSE)
          incProgress(0.8, detail = "Done")
          m$gene <- rownames(m); m
        }),
        error = function(e) {
          showNotification(paste("DE failed:", e$message), type = "error",
                           duration = 6); NULL
        })
      if (!is.null(res) && nrow(res) > 0)
        res <- res[order(-res$avg_log2FC), , drop = FALSE]
      rv$de_markers <- res
      rv$de_context_text <-
        paste0("Selection vs rest: ", length(sel), " vs ", length(rest),
               " cells, ", if (is.null(res)) 0 else nrow(res), " markers.")
      list(mode = "vs_rest", markers = res, selected = sel)
    }, ignoreInit = TRUE)

    output$de_context <- renderText({ rv$de_context_text })
    output$de_status  <- renderText({
      mk <- rv$de_markers
      if (is.null(mk) || !nrow(mk)) "No DE yet"
      else paste0(nrow(mk), " markers")
    })

    # --- DE table ----------------------------------------------------------
    output$de_table <- DT::renderDT({
      mk <- rv$de_markers
      if (is.null(mk) || !nrow(mk)) {
        return(DT::datatable(data.frame(Info = "Click 'Compute' to run DE."),
                             rownames = FALSE,
                             options = list(dom = "t")))
      }
      show <- mk
      if ("gene" %in% colnames(show)) {
        cols <- intersect(c("gene", "cluster", "avg_log2FC", "pct.1", "pct.2",
                            "p_val", "p_val_adj"), colnames(show))
      } else {
        cols <- intersect(c("cluster", "avg_log2FC", "pct.1", "pct.2",
                            "p_val", "p_val_adj"), colnames(show))
        show$gene <- rownames(show)
        cols <- unique(c("gene", cols))
      }
      show <- show[, cols, drop = FALSE]
      for (nc in intersect(c("avg_log2FC", "pct.1", "pct.2"), colnames(show)))
        show[[nc]] <- round(show[[nc]], 3)
      for (nc in intersect(c("p_val", "p_val_adj"), colnames(show)))
        show[[nc]] <- signif(show[[nc]], 3)

      DT::datatable(
        show,
        rownames = FALSE,
        selection = "single",
        filter = "top",
        options = list(
          pageLength = 15L,
          scrollY = "370px",
          scrollX = TRUE,
          dom = "ftip",
          searchHighlight = TRUE
        )
      )
    }, server = TRUE)

    # clicking a DE row adds the gene to the selection and colors by it
    observeEvent(input$de_table_rows_selected, {
      req(input$de_table_rows_selected)
      mk <- rv$de_markers; req(!is.null(mk), nrow(mk) > 0)
      g <- if ("gene" %in% colnames(mk)) mk$gene[input$de_table_rows_selected]
           else rownames(mk)[input$de_table_rows_selected]
      cur <- input$genes %||% character(0)
      updateSelectizeInput(session, "genes",
                           selected = unique(c(cur, g)))
      updateSelectInput(session, "color_mode", selected = "gene")
    })

    # --- DE heatmap --------------------------------------------------------
    output$de_heatmap <- renderPlot({
      dr <- de_result()
      req(dr, !is.null(dr$markers), nrow(dr$markers) > 0)
      obj <- .ge_prepare_for_de(values$active_dataset)
      top_n <- max(3L, as.integer(input$de_top_n %||% 20L))

      if (identical(dr$mode, "all_sets")) {
        mk <- dr$markers
        top_per <- do.call(rbind, lapply(split(mk, mk$cluster), function(sub) {
          utils::head(sub[order(-sub$avg_log2FC), , drop = FALSE], top_n)
        }))
        top_genes <- unique(top_per$gene)
        cells_groups <- dr$groups
        cells_use <- names(cells_groups)
        if (length(cells_use) > 6000L) {
          set.seed(42); idx <- sample.int(length(cells_use), 6000L)
          cells_use <- cells_use[idx]
          cells_groups <- cells_groups[idx]
        }
        grp <- factor(cells_groups, levels = unique(cells_groups))
      } else {
        mk <- dr$markers
        top_genes <- utils::head(mk$gene, top_n)
        sel <- intersect(dr$selected, colnames(obj))
        rest <- setdiff(colnames(obj), sel)
        set.seed(42)
        n_rest <- min(length(rest), max(200L, length(sel)))
        rest_sub <- sample(rest, n_rest)
        cells_use <- c(sel, rest_sub)
        grp <- factor(c(rep("Selection", length(sel)),
                        rep("Rest", length(rest_sub))),
                      levels = c("Selection", "Rest"))
      }

      dat <- .ge_safe_get_assay_data(obj, "data",
                                      .ge_safe_default_assay(obj))
      req(!is.null(dat))
      top_genes <- intersect(top_genes, rownames(dat))
      req(length(top_genes) >= 1L)

      m <- as.matrix(dat[top_genes, cells_use, drop = FALSE])

      if (identical(input$de_scale, "row")) {
        mu <- rowMeans(m, na.rm = TRUE)
        sdv <- apply(m, 1, stats::sd, na.rm = TRUE); sdv[sdv == 0 | is.na(sdv)] <- 1
        m <- sweep(sweep(m, 1, mu, "-"), 1, sdv, "/")
        m[is.na(m)] <- 0
        m[m >  3] <-  3; m[m < -3] <- -3
      }

      if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stats::heatmap(m, Rowv = NA, Colv = NA, scale = "none",
                       labCol = "", main = "Top DE genes")
        return()
      }

      col_fun <- if (requireNamespace("circlize", quietly = TRUE)) {
        if (identical(input$de_scale, "row"))
          circlize::colorRamp2(c(-2, 0, 2),
                               c("#440154", "#21908C", "#FDE725"))
        else
          circlize::colorRamp2(range(m, na.rm = TRUE),
                               c("#1a0033", "#FDE725"))
      } else NULL

      if (identical(dr$mode, "all_sets")) {
        lvls <- levels(grp)
        grp_cols <- setNames(
          .ge_set_palette[((seq_along(lvls) - 1L) %%
                             length(.ge_set_palette)) + 1L],
          lvls)
      } else {
        grp_cols <- c(Selection = "#E64B35", Rest = "#95a5a6")
      }

      ann <- ComplexHeatmap::HeatmapAnnotation(
        Group = grp,
        col = list(Group = grp_cols),
        annotation_legend_param = list(Group = list(title = "Group"))
      )

      ht <- ComplexHeatmap::Heatmap(
        m,
        name = if (identical(input$de_scale, "row")) "z-score" else "expr",
        col = col_fun,
        cluster_rows = nrow(m) > 1,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        column_split = grp,
        top_annotation = ann,
        row_names_gp = grid::gpar(fontsize = 9),
        column_title_gp = grid::gpar(fontsize = 11, fontface = "bold"),
        heatmap_legend_param = list(direction = "horizontal")
      )
      ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom",
                           annotation_legend_side = "bottom",
                           merge_legend = TRUE)
    })

    output$dl_de_csv <- downloadHandler(
      filename = function() paste0("gene_explorer_de_",
                                   format(Sys.time(), "%Y%m%d_%H%M%S"),
                                   ".csv"),
      content = function(file) {
        mk <- rv$de_markers
        if (is.null(mk) || !nrow(mk)) {
          utils::write.csv(data.frame(note = "No DE yet."), file,
                           row.names = FALSE); return()
        }
        utils::write.csv(mk, file, row.names = FALSE)
      }
    )
  })
}
