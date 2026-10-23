# Professional Single-Cell RNA-seq Analysis Platform (Enhanced)
# app.R - Main Application with Enhanced Features and Error Handling

# =============================================================================
# SYSTEM CONFIGURATION & SETUP
# =============================================================================

# Remove R's vector memory limit (default 16GB can cause SCTransform to fail)
if (.Platform$OS.type != "windows") {
  tryCatch(
    unix::rlimit_as(cur = Inf, max = Inf),
    error = function(e) NULL
  )
}
# Also set via env-style option (works on macOS R.app / RStudio)
Sys.setenv("R_MAX_VSIZE" = "128Gb")

# Set file upload limit FIRST (before any library calls).
# Production mode (set env SCRNA_PRODUCTION=1 on the server) sanitises error
# messages shown to users; development shows full errors for debugging.
.scrna_production <- nzchar(Sys.getenv("SCRNA_PRODUCTION"))
options(
  shiny.maxRequestSize = 8000 * 1024^2,          # 8 GB upload limit (large h5ad/10x)
  shiny.sanitize.errors = .scrna_production,      # hide internals from users in prod
  future.globals.maxSize = 16 * 1024^3            # 16 GB — Seurat passes full objects to workers
)

# Configure parallel processing - use sequential to avoid memory duplication
# SCTransform and other heavy operations copy the full object to each worker,
# which can exhaust system memory and cause crashes/freezes
library(future)
plan(sequential)

# Core packages
suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(shinyjs)
  library(shinycssloaders)
  library(shinyFiles)
  library(DT)
  library(plotly)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  # Load Bioconductor S4 infrastructure before Seurat to avoid
  # "as.SimpleList not found" errors with Seurat v5 layer operations
  for (.pkg in c("BiocGenerics", "S4Vectors", "SummarizedExperiment")) {
    if (requireNamespace(.pkg, quietly = TRUE)) library(.pkg, character.only = TRUE)
  }
  library(Seurat)
  library(Matrix)
  library(promises)
  library(waiter)
})

# =============================================================================
# SOURCE MODULES
# =============================================================================

# Helper function to safely source modules
safe_source <- function(file) {
  tryCatch({
    source(file)
    TRUE
  }, error = function(e) {
    warning(paste("Failed to load module:", file, "-", e$message))
    FALSE
  })
}

# Load all modules with error handling
modules_loaded <- list(
  # Central registries first (single source of truth, no deps)
  ai_models = safe_source("modules/ai_models.R"),
  species_utils = safe_source("modules/species_utils.R"),
  ondisk_storage = safe_source("modules/ondisk_storage.R"),
  # Utility modules first (define helpers used by other modules)
  core_robustness = safe_source("modules/core_robustness.R"),
  ai_assistant = safe_source("modules/ai_assistant.R"),
  plot_popup = safe_source("modules/plot_popup.R"),
  umap_gallery = safe_source("modules/umap_gallery.R"),
  # Data and analysis modules
  markers = safe_source("modules/cell_type_markers.R"),
  data_input = safe_source("modules/data_input.R"),
  quality_control = safe_source("modules/quality_control.R"),
  preprocessing = safe_source("modules/preprocessing.R"),
  dimensionality_reduction = safe_source("modules/dimensionality_reduction.R"),
  clustering = safe_source("modules/clustering.R"),
  annotation = safe_source("modules/annotation.R"),
  gene_explorer = safe_source("modules/gene_explorer.R"),
  trajectory = safe_source("modules/trajectory.R"),
  visualization = safe_source("modules/visualization.R"),
  downstream = safe_source("modules/downstream/downstream_main.R"),
  de = safe_source("modules/de/de_server.R"),
  pn = safe_source("modules/pathway_network/pn_server.R"),
  export = safe_source("modules/export.R")
)

# Check if all modules loaded successfully
failed_modules <- names(modules_loaded)[!unlist(modules_loaded)]
if (length(failed_modules) > 0) {
  warning(paste("Failed to load modules:", paste(failed_modules, collapse = ", ")))
}

# =============================================================================
# AUTH & PROJECT MODULES
# =============================================================================
safe_source("modules/auth/auth_crypto.R")
safe_source("modules/auth/auth_db.R")
safe_source("modules/auth/auth_module.R")
safe_source("modules/auth/project_module.R")

# Initialise SQLite auth database (creates tables on first run)
tryCatch(init_auth_db(), error = function(e) warning("Auth DB init: ", e$message))

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Null coalescing operator
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Enhanced color palette function
get_enhanced_color_palette <- function(n, type = "qualitative") {
  n <- max(1L, as.integer(n))
  
  if (type == "qualitative") {
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      if (n <= 3) {
        RColorBrewer::brewer.pal(3, "Set2")[1:n]
      } else if (n <= 8) {
        RColorBrewer::brewer.pal(n, "Set2")
      } else if (n <= 12) {
        RColorBrewer::brewer.pal(n, "Set3")
      } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
      }
    } else {
      grDevices::hcl.colors(n, "Dynamic")
    }
  } else if (type == "sequential") {
    if (requireNamespace("viridis", quietly = TRUE)) {
      viridis::viridis(n)
    } else {
      grDevices::hcl.colors(n, "Viridis")
    }
  } else {
    RColorBrewer::brewer.pal(min(11, max(3, n)), "RdYlBu")
  }
}

# Memory usage monitoring
get_memory_usage <- function() {
  tryCatch({
    mem <- gc(verbose = FALSE)
    used_mb <- sum(mem[, 2])  # Mb used
    paste(round(used_mb), "MB R memory used")
  }, error = function(e) "Memory info unavailable")
}

# Session validation
validate_session <- function(values) {
  issues <- character(0)
  
  if (is.null(values$active_dataset)) {
    issues <- c(issues, "No active dataset")
  }
  
  if (length(values$seurat_list) == 0) {
    issues <- c(issues, "No datasets loaded")
  }
  
  if (length(issues) > 0) {
    return(list(valid = FALSE, issues = issues))
  } else {
    return(list(valid = TRUE, issues = character(0)))
  }
}

# =============================================================================
# ENHANCED UI
# =============================================================================

main_ui_def <- dashboardPage(
  title = "scRNA-seq Analysis Platform",
  skin = "blue",

  # ============================= HEADER =====================================
  dashboardHeader(
    title = "",
    titleWidth = 0,

    # Centered platform title in navbar
    tags$li(class = "dropdown navbar-brand-wrapper",
      tags$a(class = "navbar-brand-center", href = "#",
        tags$i(class = "fa fa-dna navbar-brand-icon"),
        tags$span("scRNA-seq", class = "brand-text-main"),
        tags$span("Platform", class = "brand-text-sub")
      )
    ),

    # Notification dropdown
    dropdownMenu(
      type = "notifications",
      headerText = "System Status",
      icon = icon("bell"),

      notificationItem(
        text = "Memory usage",
        icon = icon("memory"),
        status = "info"
      ),

      notificationItem(
        text = "Session validation",
        icon = icon("check-circle"),
        status = "success"
      )
    ),

    # Task dropdown
    dropdownMenu(
      type = "tasks",
      headerText = "Analysis Progress",
      icon = icon("tasks"),

      taskItem(
        text = "Data Loading",
        value = 0,
        color = "blue"
      )
    )
  ),

  # ============================= SIDEBAR ====================================
  dashboardSidebar(
    width = 300,

    # User / project identity banner (rendered in server after login)
    uiOutput("sidebarUserBanner"),

    sidebarMenu(
      id = "tabs",
      
      # Analysis workflow menu items
      menuItem("Data Upload & Management",
        tabName = "data",
        icon = icon("database"),
        badgeLabel = "Start",
        badgeColor = "green"
      ),

      menuItem("Quality Control",
        tabName = "qc",
        icon = icon("check-circle"),
        badgeLabel = "QC",
        badgeColor = "blue"
      ),

      menuItem("Preprocessing",
        tabName = "preprocess",
        icon = icon("cogs"),
        badgeLabel = "Prep",
        badgeColor = "orange"
      ),

      menuItem("Dimensionality Reduction",
        tabName = "dimred",
        icon = icon("compress"),
        badgeLabel = "DimRed",
        badgeColor = "purple"
      ),

      menuItem("Clustering",
        tabName = "cluster",
        icon = icon("project-diagram"),
        badgeLabel = "Cluster",
        badgeColor = "red"
      ),

      menuItem("Cell Type Annotation",
        tabName = "annotate",
        icon = icon("tags"),
        badgeLabel = "Annotate",
        badgeColor = "teal"
      ),

      menuItem("Visualization",
        tabName = "visualize",
        icon = icon("chart-line"),
        badgeLabel = "Viz",
        badgeColor = "lime"
      ),

      menuItem("Downstream Analysis",
        tabName = "downstream",
        icon = icon("microscope"),
        badgeLabel = "DE",
        badgeColor = "navy"
      ),

      menuItem("Trajectory Analysis",
        tabName = "trajectory",
        icon = icon("route"),
        badgeLabel = "Trajectory",
        badgeColor = "yellow"
      ),

      menuItem("Pathway & Network Analysis",
        tabName = "pathway_network",
        icon = icon("project-diagram"),
        badgeLabel = "Pathway",
        badgeColor = "olive"
      ),

      menuItem("Comparison View",
        tabName = "comparison",
        icon = icon("columns"),
        badgeLabel = "Compare",
        badgeColor = "aqua"
      ),

      menuItem("Export & Report",
        tabName = "export",
        icon = icon("file-export"),
        badgeLabel = "Export",
        badgeColor = "maroon"
      ),

      hr(),

      # Session management section
      div(class = "session-controls sidebar-section",
        tags$h5(class = "sidebar-section-title", icon("sliders-h"), "Session"),
        div(style = "display: flex; gap: 6px; margin-bottom: 6px;",
          actionButton("saveSession", "Save Session",
            icon = icon("cloud-upload-alt"),
            class = "btn-success btn-sm",
            style = "flex: 1; padding: 6px 10px; font-size: 12px; box-sizing: border-box;"
          ),
          actionButton("loadSession", "Load Session",
            icon = icon("cloud-download-alt"),
            class = "btn-info btn-sm",
            style = "flex: 1; padding: 6px 10px; font-size: 12px; box-sizing: border-box;"
          )
        ),
        actionButton("resetSession", "Reset",
          icon = icon("redo-alt"),
          class = "btn-outline-light btn-sm",
          style = "display: block; width: 50%; margin: 4px auto 0; padding: 4px 8px; font-size: 11px; box-sizing: border-box; text-align: center;",
          onclick = "return confirm('Reset your entire session? All unsaved work will be lost.');"
        )
      ),

      hr(),

      # Active dataset info
      div(class = "sidebar-section sidebar-card",
        tags$h5(class = "sidebar-section-title", icon("database"), "Active Dataset"),
        uiOutput("datasetInfo")
      ),

      hr(),

      # Analysis progress
      div(class = "sidebar-section",
        tags$h5(class = "sidebar-section-title", icon("chart-line"), "Progress"),
        uiOutput("progressInfo")
      ),

      hr(),

      # System info
      div(class = "sidebar-section",
        tags$h5(class = "sidebar-section-title", icon("microchip"), "System"),
        uiOutput("systemInfo"),
        br(),
        uiOutput("scalableInfo"),
        br(),
        uiOutput("jobsInfo")
      )
    )
  ),

  # ============================== BODY ======================================
  dashboardBody(
    # Error handling
    div(id = "error-container"),
    
    # Loading screen
    div(
      id = "loading-screen",
      style = "position: fixed; top: 0; left: 0; width: 100%; height: 100%; 
              background: rgba(255,255,255,0.9); z-index: 9999; display: none;",
      div(
        style = "position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);
                text-align: center;",
        icon("spinner", class = "fa-spin fa-3x"),
        h3("Loading...", style = "margin-top: 20px;")
      )
    ),

    # Main content tabs
    tabItems(

      # ========================== DATA MANAGEMENT ==========================
      tabItem(
        tabName = "data",

        # Welcome message for new users (controlled by JavaScript show/hide)
        div(
          id = "welcome-panel",
          style = "display: none;",
          fluidRow(
            column(12,
              div(class = "welcome-card",
                div(class = "welcome-card-header",
                  icon("dna", class = "fa-lg"), "Welcome to scRNA-seq Analysis Platform"
                ),
                div(class = "welcome-card-body",
                  p(icon("upload", class = "welcome-icon"), "Start by uploading your single-cell RNA sequencing data below."),
                  p(icon("file-alt", class = "welcome-icon"),
                    tags$strong("Supported formats: "),
                    "10X H5, H5AD, Seurat RDS, CSV/TSV, Excel, and more."
                  )
                )
              )
            )
          )
        ),

        fluidRow(
          box(
            title = div(icon("database"), "Data Upload & Management"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$data_input) {
              tabsetPanel(
                id = "data_tabs",
                tabPanel("Upload Files",
                  icon = icon("upload"),
                  br(),
                  dataUploadUI("dataUpload")
                ),
                tabPanel("Batch Upload",
                  icon = icon("layer-group"),
                  br(),
                  batchUploadUI("batchUpload")
                ),
                tabPanel("Manage Datasets",
                  icon = icon("folder"),
                  br(),
                  dataManagementUI("dataManage")
                ),
                tabPanel("Data Validation",
                  icon = icon("check-double"),
                  br(),
                  dataValidationUI("dataValidate")
                )
              )
            } else {
              p("Data input module not available. Check R console for details.")
            }
          )
        )
      ),

      # ========================= QUALITY CONTROL =========================
      tabItem(
        tabName = "qc",
        fluidRow(
          box(
            title = div(icon("check-circle"), "Advanced Quality Control"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$quality_control) {
              tabsetPanel(
                tabPanel("QC Metrics",
                  icon = icon("chart-bar"),
                  br(),
                  qcMetricsUI("qcMetrics")
                ),
                tabPanel("Advanced Filtering",
                  icon = icon("filter"),
                  br(),
                  qcFilteringUI("qcFilter")
                ),
                tabPanel("Doublet Detection",
                  icon = icon("users"),
                  br(),
                  doubletDetectionUI("doublet")
                ),
                tabPanel("Ambient RNA",
                  icon = icon("flask"),
                  br(),
                  ambientRNAUI("ambient")
                ),
                tabPanel("Cell Cycle",
                  icon = icon("sync"),
                  br(),
                  cellCycleUI("cellCycle")
                ),
                tabPanel("Batch Effects",
                  icon = icon("layer-group"),
                  br(),
                  batchEffectsUI("batchEffects")
                ),
                tabPanel("QC Report",
                  icon = icon("file-alt"),
                  br(),
                  qcReportUI("qcReport")
                )
              )
            } else {
              p("Quality control module not available. Check R console for details.")
            }
          )
        )
      ),

      # ========================= PREPROCESSING ===========================
      tabItem(
        tabName = "preprocess",
        fluidRow(
          box(
            title = div(icon("cogs"), "Preprocessing & Integration"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$preprocessing) {
              tabsetPanel(
                tabPanel("Normalization",
                  icon = icon("balance-scale"),
                  br(),
                  normalizationUI("normalize")
                ),
                tabPanel("Feature Selection",
                  icon = icon("crosshairs"),
                  br(),
                  featureSelectionUI("features")
                ),
                tabPanel("Scaling",
                  icon = icon("ruler"),
                  br(),
                  scalingUI("scale")
                ),
                tabPanel("Integration",
                  icon = icon("puzzle-piece"),
                  br(),
                  integrationUI("integrate")
                )
              )
            } else {
              p("Preprocessing module not available. Check R console for details.")
            }
          )
        )
      ),

      # =================== DIMENSIONALITY REDUCTION ===================
      tabItem(
        tabName = "dimred",
        fluidRow(
          box(
            title = div(icon("compress"), "Dimensionality Reduction"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$dimensionality_reduction) {
              tabsetPanel(
                tabPanel("PCA",
                  icon = icon("chart-line"),
                  br(),
                  pcaUI("pca")
                ),
                tabPanel("UMAP / t-SNE",
                  icon = icon("map"),
                  br(),
                  umapTsneUI("umapTsne")
                ),
                tabPanel("3D Embeddings",
                  icon = icon("cube"),
                  br(),
                  embedding3dUI("embed3d")
                ),
                tabPanel("Diffusion Maps",
                  icon = icon("water"),
                  br(),
                  diffusionUI("diffusion")
                )
              )
            } else {
              p("Dimensionality reduction module not available. Check R console for details.")
            }
          )
        )
      ),

      # ========================== CLUSTERING ============================
      tabItem(
        tabName = "cluster",
        fluidRow(
          box(
            title = div(icon("project-diagram"), "Advanced Clustering Analysis"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$clustering) {
              tabsetPanel(
                tabPanel("Graph Clustering",
                  icon = icon("share-alt"),
                  br(),
                  graphClusteringUI("graphCluster")
                ),
                tabPanel("Sub-clustering",
                  icon = icon("search-plus"),
                  br(),
                  subclusteringUI("subCluster")
                )
              )
            } else {
              p("Clustering module not available. Check R console for details.")
            }
          )
        )
      ),

      # ========================= ANNOTATION =============================
      tabItem(
        tabName = "annotate",
        fluidRow(
          box(
            title = div(icon("tags"), "Cell Type Annotation"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$annotation) {
              tabsetPanel(
                tabPanel("Automatic",
                  icon = icon("magic"),
                  br(),
                  automaticAnnotationUI("autoAnn")
                ),
                tabPanel("Validation",
                  icon = icon("check"),
                  br(),
                  markerValidationUI("markValid")
                ),
                tabPanel("Manual",
                  icon = icon("edit"),
                  br(),
                  manualAnnotationUI("manualAnn"),
                  hr(),
                  clusterEditorUI("clusterEditor")
                ),
                tabPanel("AI Annotation",
                  icon = icon("robot"),
                  br(),
                  aiAnnotationUI("aiAnn")
                ),
                tabPanel("Gene Explorer",
                  icon = icon("compass"),
                  br(),
                  if (isTRUE(modules_loaded$gene_explorer)) {
                    geneExplorerUI("geneExp")
                  } else {
                    p("Gene Explorer module not available. Check R console for details.")
                  }
                )
              )
            } else {
              p("Annotation module not available. Check R console for details.")
            }
          )
        )
      ),

      # ========================= TRAJECTORY =============================
      tabItem(
        tabName = "trajectory",
        fluidRow(
          box(
            title = div(icon("route"), "Trajectory Analysis"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$trajectory) {
              trajectoryUI("trajectory")
            } else {
              p("Trajectory module not available. Check R console for details.")
            }
          )
        )
      ),

      # ========================= DOWNSTREAM =============================
      tabItem(
        tabName = "downstream",
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
      ),

      # ======================== VISUALIZATION ===========================
      tabItem(
        tabName = "visualize",
        fluidRow(
          box(
            title = div(icon("chart-line"), "Advanced Visualization Suite"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$visualization) {
              tabsetPanel(
                tabPanel("AI Plot Builder",
                  icon = icon("robot"),
                  br(),
                  aiPlotBuilderUI("viz")
                ),
                tabPanel("Interactive Plots",
                  icon = icon("chart-bar"),
                  br(),
                  interactivePlotsUI("viz")
                ),
                tabPanel("Expression Explorer",
                  icon = icon("search"),
                  br(),
                  expressionExplorerUI("exprExplorer")
                ),
                tabPanel("Custom Heatmap",
                  icon = icon("th"),
                  br(),
                  heatmapBuilderUI("heatmapBuilder")
                ),
                tabPanel("Cluster Comparison",
                  icon = icon("columns"),
                  br(),
                  clusterComparisonUI("clusterComp")
                ),
                tabPanel("Publication Figures",
                  icon = icon("image"),
                  br(),
                  figurePanelUI("figPanel")
                ),
                if (modules_loaded$umap_gallery) tabPanel("Figure & Export Gallery",
                  icon = icon("images"),
                  br(),
                  umapGalleryUI("umapGallery")
                )
              )
            } else {
              p("Visualization module not available. Check R console for details.")
            }
          )
        )
      ),

      # ====================== COMPARISON VIEW ==========================
      tabItem(
        tabName = "comparison",
        fluidRow(
          box(
            title = div(icon("columns"), "Split-Screen Comparison View"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$visualization) {
              comparisonViewUI("comparisonView")
            } else {
              p("Visualization module not available. Check R console for details.")
            }
          )
        )
      ),

      # ================== PATHWAY & NETWORK ANALYSIS ====================
      tabItem(
        tabName = "pathway_network",
        fluidRow(
          box(
            title = div(icon("project-diagram"), "Pathway & Network Analysis"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            if (modules_loaded$pn) pnUI("pn")
            else p("Pathway & Network module not available. Check R console for details.")
          )
        )
      ),

      # =========================== EXPORT ===============================
      tabItem(
        tabName = "export",
        fluidRow(
          box(
            title = div(icon("file-export"), "Export & Reporting"),
            status = "primary",
            solidHeader = TRUE,
            width = 12,

            if (modules_loaded$export) {
              dataExportUI("exporter")
            } else {
              p("Export module not available. Check R console for details.")
            }
          )
        )
      )
    ),

    # Console log panel
    div(
      id = "console-panel",
      style = "position: fixed; bottom: 0; left: 230px; right: 0; z-index: 950;",
      div(
        style = "background: #2c3e50; color: #ecf0f1; cursor: pointer; padding: 6px 15px; font-size: 12px;",
        onclick = "$('#console-output').slideToggle(200);",
        icon("terminal"), " Console Log ",
        span(id = "console-toggle", style = "float: right;", icon("chevron-up"))
      ),
      div(
        id = "console-output",
        style = "display: none; background: #1a252f; color: #2ecc71; font-family: 'Courier New', monospace;
                 font-size: 11px; max-height: 200px; overflow-y: auto; padding: 10px; white-space: pre-wrap;",
        verbatimTextOutput("consoleLog")
      )
    )
  )
)

# =============================================================================
# ROOT UI  — Auth / Project overlay  +  Main analysis dashboard
# =============================================================================

ui <- tagList(
  # Global head (CSS + shinyjs + waiter)
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = NA),
    tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),
    tags$link(rel = "stylesheet", href = "app.css"),
    tags$link(rel = "stylesheet", href = "ai.css"),
    tags$link(rel = "stylesheet", href = "auth.css")
  ),
  useShinyjs(),
  useWaiter(),
  if (exists("plotPopupJS", mode = "function")) plotPopupJS(),

  # ── Auth / Project selection overlay (shown when not in analysis mode) ──
  div(id = "page_auth_projects",
    uiOutput("auth_project_ui")
  ),

  # ── Main analysis dashboard (hidden until user opens a project) ─────────
  div(id = "page_main", style = "display:none;",
    main_ui_def
  )
)

# =============================================================================
# ENHANCED SERVER
# =============================================================================

server <- function(input, output, session) {

  # ======================== AUTH STATE =================================

  # Tracks login / project selection across the whole session
  auth_state <- reactiveValues(
    logged_in                         = FALSE,
    user_id                           = NULL,
    username                          = NULL,
    email                             = NULL,
    full_name                         = NULL,
    session_token                     = NULL,
    current_project_id                = NULL,
    current_project_name              = NULL,
    current_project_species           = "",
    current_project_tissue            = "",
    current_project_conditions        = "",
    current_project_description       = "",
    current_project_biological_question = "",
    current_project_experiment_design = "",
    current_project_background        = ""
  )

  # ── Dynamic page UI (auth page / project dashboard / nothing) ────────────
  output$auth_project_ui <- renderUI({
    if (!isTRUE(auth_state$logged_in)) {
      authPageUI("auth")
    } else if (is.null(auth_state$current_project_id)) {
      projectDashboardUI("projects")
    } else {
      NULL   # main app is shown via shinyjs
    }
  })

  # ── Page visibility controller ────────────────────────────────────────────
  observe({
    in_main_app <- isTRUE(auth_state$logged_in) && !is.null(auth_state$current_project_id)
    if (in_main_app) {
      shinyjs::hide("page_auth_projects")
      shinyjs::show("page_main")
    } else {
      shinyjs::show("page_auth_projects")
      shinyjs::hide("page_main")
    }
  })

  # ── Auth module server (handles login / register) ─────────────────────────
  authPageServer("auth", auth_state)

  # ── Project dashboard server (handles create / open / delete) ────────────
  projectDashboardServer("projects", auth_state)

  # ── On project open: load saved state + populate values ──────────────────
  observeEvent(auth_state$current_project_id, {
    req(auth_state$current_project_id)
    pid <- auth_state$current_project_id
    uid <- auth_state$user_id

    # Stable, persistent per-project directory for optional on-disk (BPCells)
    # matrix storage — must live with the session RDS so a saved object's
    # on-disk pointer stays valid across restarts.
    values$ondisk_dir <- file.path(getwd(), "data", "user_sessions",
                                   as.character(uid), as.character(pid), "ondisk")

    # Restore project context into values
    values$species_global       <- auth_state$current_project_species        %||% ""
    values$tissue_global        <- auth_state$current_project_tissue         %||% ""
    values$conditions_global    <- auth_state$current_project_conditions     %||% ""
    values$biological_question  <- auth_state$current_project_biological_question %||% ""
    values$experiment_design    <- auth_state$current_project_experiment_design   %||% ""
    values$project_background   <- auth_state$current_project_background     %||% ""
    # AI assistant per-project config (used by ai_assistant.R + downstream)
    values$ai_api_key           <- auth_state$current_project_ai_api_key     %||% ""
    values$ai_model             <- auth_state$current_project_ai_model       %||% get_default_model()
    values$ai_custom_persona    <- auth_state$current_project_ai_persona     %||% ""
    values$ai_project_summary   <- auth_state$current_project_ai_summary     %||% ""

    # Auto-populate per-module species / tissue inputs (module namespaces)
    # so users don't have to re-enter what they gave during project creation.
    sp <- values$species_global
    ts <- values$tissue_global
    sp_norm <- normalize_species(sp)   # "" when unknown/other
    # species selects (some use human/mouse vocabulary, others latin — push the
    # canonical value where it matches; updateSelectInput is a no-op otherwise)
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
    # tissue text inputs (incl. Advanced QC Filtering's "Tissue Context for AI")
    for (full_id in c("qcMetrics-tissue_context", "qcFilter-tissue_context",
                      "aiAnn-tissue_context", "autoAnn-tissue")) {
      tryCatch(updateTextInput(session, full_id, value = ts),
               error = function(e) NULL)
    }

    # Load saved analysis state (if any)
    saved <- tryCatch(load_project_state(pid, uid), error = function(e) NULL)
    if (!is.null(saved)) {
      withProgress(message = paste0("Loading project: ", auth_state$current_project_name), value = 0, {
        incProgress(0.3, detail = "Restoring datasets…")
        if (!is.null(saved$seurat_list))    values$seurat_list    <- saved$seurat_list
        if (!is.null(saved$active_dataset)) values$active_dataset <- saved$active_dataset
        if (!is.null(saved$active_name))    values$active_name    <- saved$active_name
        if (!is.null(saved$cluster_colors)) values$cluster_colors <- saved$cluster_colors

        incProgress(0.3, detail = "Restoring analysis flags…")
        if (!is.null(saved$analysis_flags)) {
          for (flag in names(saved$analysis_flags)) {
            values[[flag]] <- saved$analysis_flags[[flag]]
          }
        }
        if (!is.null(saved$dataset_version))   values$dataset_version   <- saved$dataset_version
        if (!is.null(saved$dataset_signature))  values$dataset_signature  <- saved$dataset_signature
        if (!is.null(saved$dataset_history))    values$dataset_history    <- saved$dataset_history
        if (!is.null(saved$analysis_audit))     values$analysis_audit     <- saved$analysis_audit
        if (!is.null(saved$job_queue))          values$job_queue          <- saved$job_queue
        if (!is.null(saved$method_inventory))   values$method_inventory   <- saved$method_inventory

        # Restore generated figures / tables / gallery snapshots
        incProgress(0.1, detail = "Restoring figures & tables…")
        if (!is.null(saved$figure_gallery))    values$figure_gallery    <- saved$figure_gallery
        if (!is.null(saved$csv_library))       values$csv_library       <- saved$csv_library
        if (!is.null(saved$umap_gallery))      values$umap_gallery      <- saved$umap_gallery
        if (!is.null(saved$cached_embeddings)) values$cached_embeddings <- saved$cached_embeddings

        if (length(values$seurat_list) > 0) {
          updateSelectInput(session, "active_dataset_select",
                            choices  = names(values$seurat_list),
                            selected = values$active_name)
        }
        incProgress(0.4, detail = "Done!")
        showNotification(paste0("Project \"", auth_state$current_project_name, "\" loaded."),
                         type = "message", duration = 4)
      })
    } else {
      showNotification(
        paste0("Project \"", auth_state$current_project_name, "\" opened. Upload your data to begin."),
        type = "message", duration = 5
      )
    }
  }, ignoreInit = TRUE)

  # ── Auto-detect species from the active dataset and pre-set QC controls ────
  observeEvent(values$active_dataset, {
    req(values$active_dataset)
    det <- tryCatch(detect_species(rownames(values$active_dataset)),
                    error = function(e) NULL)
    if (is.null(det) || det$species == "unknown") return()
    sel <- if (det$species %in% c("human", "mouse")) det$species else "custom"
    tryCatch(updateSelectInput(session, "qcMetrics-species_qc", selected = sel),
             error = function(e) NULL)
    if (sel == "custom")
      tryCatch(updateTextInput(session, "qcMetrics-custom_mt_pattern",
                               value = det$mt_pattern), error = function(e) NULL)
    showNotification(sprintf("Auto-detected species: %s (%s confidence).",
                             det$species, det$confidence),
                     type = "message", duration = 5)
  }, ignoreInit = TRUE)

  # ── Sidebar user / project banner ─────────────────────────────────────────
  output$sidebarUserBanner <- renderUI({
    if (!isTRUE(auth_state$logged_in) || is.null(auth_state$current_project_id)) return(NULL)

    uname    <- auth_state$full_name %||% auth_state$username %||% "User"
    initials <- paste(toupper(substr(strsplit(uname, " ")[[1]], 1, 1)), collapse = "")
    initials <- substr(initials, 1, 2)
    proj_name <- auth_state$current_project_name %||% "—"

    tags$div(class = "sidebar-user-banner",
      tags$div(class = "sub-banner-user",
        tags$div(class = "sidebar-avatar", initials),
        tags$div(class = "sidebar-user-text",
          tags$span(class = "sidebar-user-name",  uname),
          tags$span(class = "sidebar-user-role",  "Researcher")
        )
      ),
      tags$div(class = "sidebar-project-name",
        tags$i(class = "fa fa-folder-open"), " ",
        tags$strong(proj_name)
      ),
      actionButton("sidebar_switch_project", "Switch Project",
                   icon  = icon("exchange-alt"),
                   class = "btn-sidebar-switch"),
      actionButton("sidebar_logout", "Sign Out",
                   icon  = icon("sign-out-alt"),
                   class = "btn-sidebar-logout")
    )
  })

  # ── Switch project: save current state then return to project dashboard ───
  observeEvent(input$sidebar_switch_project, {
    req(auth_state$current_project_id)
    pid <- auth_state$current_project_id
    uid <- auth_state$user_id

    # Auto-save the full analysis state (datasets + figures + tables) with visible
    # progress, then return to the dashboard.
    .persist_state(show_progress = TRUE)

    # Reset values and go back to project dashboard
    .reset_values()
    auth_state$current_project_id   <- NULL
    auth_state$current_project_name <- NULL
  })

  # ── Sidebar logout ────────────────────────────────────────────────────────
  observeEvent(input$sidebar_logout, {
    # Save the full analysis state (datasets + figures + tables) before logging out
    .persist_state(show_progress = TRUE)
    invalidate_session_token(auth_state$session_token)
    .reset_values()
    auth_state$logged_in          <- FALSE
    auth_state$user_id            <- NULL
    auth_state$username           <- NULL
    auth_state$full_name          <- NULL
    auth_state$session_token      <- NULL
    auth_state$current_project_id <- NULL
    auth_state$current_project_name <- NULL
  })

  # ── Helper: collect the FULL analysis state for persistence ──────────────
  # Single source of truth for what gets saved (datasets, flags, provenance AND
  # the generated figures / tables / gallery snapshots that were previously lost).
  .collect_session_data <- function() {
    list(
      seurat_list    = values$seurat_list,
      active_dataset = values$active_dataset,
      active_name    = values$active_name,
      cluster_colors = values$cluster_colors,
      analysis_flags = list(
        data_loaded        = values$data_loaded,
        qc_done            = values$qc_done,
        normalization_done = values$normalization_done,
        preprocessing_done = values$preprocessing_done,
        features_selected  = values$features_selected,
        scaling_done       = values$scaling_done,
        dimred_done        = values$dimred_done,
        clustering_done    = values$clustering_done,
        annotation_done    = values$annotation_done,
        trajectory_done    = values$trajectory_done
      ),
      dataset_version   = values$dataset_version,
      dataset_signature = values$dataset_signature,
      dataset_history   = values$dataset_history,
      analysis_audit    = values$analysis_audit,
      job_queue         = values$job_queue,
      method_inventory  = values$method_inventory,
      # Generated outputs — figures, tables, gallery snapshots
      figure_gallery    = values$figure_gallery,
      csv_library       = values$csv_library,
      umap_gallery      = values$umap_gallery,
      cached_embeddings = values$cached_embeddings
    )
  }

  # Signature that changes whenever there is something new worth saving.
  .state_signature <- function() {
    paste(
      values$dataset_version   %||% 0L,
      length(values$seurat_list),
      length(values$figure_gallery),
      length(values$csv_library),
      length(values$umap_gallery),
      sep = "_"
    )
  }

  # Persist the current state to the open project (safe no-op if none open).
  .persist_state <- function(show_progress = FALSE) {
    pid <- auth_state$current_project_id
    uid <- auth_state$user_id
    if (is.null(pid) || is.null(uid)) return(invisible(FALSE))
    do_save <- function() {
      tryCatch({
        save_project_state(pid, uid, .collect_session_data())
        last_saved_signature(.state_signature())
        TRUE
      }, error = function(e) {
        showNotification(paste("Save failed:", conditionMessage(e)),
                         type = "error", duration = 8)
        FALSE
      })
    }
    if (isTRUE(show_progress)) {
      withProgress(message = "Saving project…", value = 0.4, {
        setProgress(0.7, detail = "Writing to disk…")
        ok <- do_save()
        setProgress(1, detail = "Done")
        ok
      })
    } else do_save()
  }

  # Tracks the last-persisted state signature (for autosave change detection).
  last_saved_signature <- reactiveVal("")

  # ── Helper: reset all analysis values ────────────────────────────────────
  .reset_values <- function() {
    values$seurat_list        <- list()
    values$active_dataset     <- NULL
    values$active_name        <- NULL
    values$data_loaded        <- FALSE
    values$qc_done            <- FALSE
    values$normalization_done <- FALSE
    values$preprocessing_done <- FALSE
    values$features_selected  <- FALSE
    values$scaling_done       <- FALSE
    values$dimred_done        <- FALSE
    values$clustering_done    <- FALSE
    values$annotation_done    <- FALSE
    values$trajectory_done    <- FALSE
    values$last_de_results    <- NULL
    values$last_enrich_results <- NULL
    values$error_count        <- 0
    values$analysis_start_time <- Sys.time()
    values$dataset_version    <- 0L
    values$dataset_signature  <- NA_character_
    values$dataset_history    <- list()
    values$analysis_audit     <- values$analysis_audit[0, , drop = FALSE]
    values$job_queue          <- list()
    values$scalable_mode      <- FALSE
    values$project_background <- ""
    values$experiment_design  <- ""
    values$biological_question <- ""
    values$species_global     <- ""
    values$tissue_global      <- ""
    values$conditions_global  <- ""
    values$cluster_colors     <- NULL
    values$umap_gallery       <- list()
    values$figure_gallery     <- list()
    values$csv_library        <- list()
    values$cached_embeddings  <- list()
    values$ai_cache           <- list()
    # Clear per-project AI credentials so they never leak into the next project
    values$ai_api_key         <- NULL
    values$ai_model           <- get_default_model()
    values$ai_custom_persona  <- ""
    values$ai_project_summary <- ""
  }

  # ======================== REACTIVE VALUES ============================

  # Central state shared across modules
  values <- reactiveValues(
    seurat_list = list(),
    active_dataset = NULL,
    active_name = NULL,
    
    # Analysis status flags
    data_loaded = FALSE,
    qc_done = FALSE,
    normalization_done = FALSE,
    preprocessing_done = FALSE,
    features_selected = FALSE,
    scaling_done = FALSE,
    dimred_done = FALSE,
    clustering_done = FALSE,
    annotation_done = FALSE,
    trajectory_done = FALSE,
    
    # AI assistant (global config)
    ai_api_key = NULL,
    ai_model = get_default_model(),
    ai_action_queue = list(),

    # Project context (for AI prompts)
    project_background = "",
    experiment_design = "",
    biological_question = "",
    species_global = "",
    tissue_global = "",
    conditions_global = "",

    # Results storage
    last_de_results = NULL,
    last_enrich_results = NULL,
    last_qc_metrics = NULL,
    
    # System info
    memory_usage = "",
    analysis_start_time = Sys.time(),
    
    # Error tracking
    last_error = NULL,
    error_count = 0,

    # UMAP gallery snapshots
    umap_gallery = list(),

    # Universal figure gallery
    figure_gallery = list(),

    # CSV export library
    csv_library = list(),

    # Custom cluster colors from cluster editor
    cluster_colors = NULL,

    # AI response cache (keyed by dataset_signature::task_type::prompt_digest)
    ai_cache = list(),

    # Cached UMAP/reduction embedding data frames
    cached_embeddings = list(),

    # Robustness/provenance state
    dataset_version = 0L,
    dataset_signature = NA_character_,
    dataset_history = list(),
    analysis_audit = data.frame(
      timestamp = character(0),
      module = character(0),
      action = character(0),
      level = character(0),
      details = character(0),
      stringsAsFactors = FALSE
    ),
    job_queue = list(),
    scalable_mode = FALSE,
    method_inventory = NULL
  )

  if (exists("init_robustness_state", mode = "function")) {
    init_robustness_state(values)
  }

  # ===================== ERROR HANDLING ================================
  
  # Global error handler — an uncaught error in one operation must never take the
  # session down. Recover the parallel plan (a heavy op may have left workers) and
  # reclaim memory so the app stays responsive after an aborted/failed task.
  options(shiny.error = function() {
    values$error_count <- (values$error_count %||% 0) + 1
    values$last_error  <- geterrmessage()
    tryCatch(future::plan(future::sequential), error = function(e) NULL)
    tryCatch(gc(FALSE), error = function(e) NULL)
    showNotification(
      if (isTRUE(.scrna_production))
        "Something went wrong with that operation. It was cancelled safely — you can try again."
      else paste("Error occurred:", values$last_error),
      type = "error", duration = 10
    )
  })

  # ===================== SIDEBAR UPDATES ===============================
  
  # Dataset information display
  output$datasetInfo <- renderUI({
    if (is.null(values$active_dataset)) {
      return(
        div(
          style = "text-align: center; opacity: 0.7;",
          icon("exclamation-triangle"),
          br(),
          "No dataset loaded",
          br(),
          tags$small("Upload data to get started")
        )
      )
    }
    
    obj <- values$active_dataset
    
    tagList(
      div(
        style = "display: flex; align-items: center; margin-bottom: 8px;",
        icon("database", style = "margin-right: 8px;"),
        strong(values$active_name)
      ),
      
      div(
        style = "font-size: 12px; line-height: 1.4;",
        span(icon("circle", style = "color: #2ecc71; margin-right: 5px;"), 
             paste(format(ncol(obj), big.mark = ","), "cells")),
        br(),
        span(icon("circle", style = "color: #3498db; margin-right: 5px;"), 
             paste(format(nrow(obj), big.mark = ","), "genes")),
        br(),
        span(icon("circle", style = "color: #e74c3c; margin-right: 5px;"), 
             paste("Assays:", paste(names(obj@assays), collapse = ", "))),
        br(),
        
        if (length(obj@reductions) > 0) {
          span(icon("circle", style = "color: #f39c12; margin-right: 5px;"), 
               paste("Reductions:", paste(names(obj@reductions), collapse = ", ")))
        }
      )
    )
  })

  # Analysis progress indicators (enhanced stepper with active/completed/pending states)
  output$progressInfo <- renderUI({
    steps <- list(
      list(name = "Data Loading", flag = values$data_loaded),
      list(name = "Quality Control", flag = values$qc_done),
      list(name = "Normalization", flag = values$normalization_done),
      list(name = "Feature Selection", flag = values$features_selected),
      list(name = "Scaling", flag = values$scaling_done),
      list(name = "Preprocessing", flag = values$preprocessing_done),
      list(name = "Dim Reduction", flag = values$dimred_done),
      list(name = "Clustering", flag = values$clustering_done),
      list(name = "Annotation", flag = values$annotation_done)
    )

    step_items <- lapply(seq_along(steps), function(i) {
      step <- steps[[i]]
      if (isTRUE(step$flag)) {
        # Completed
        div(style = "padding: 4px 0; font-size: 12px; color: #2ecc71;",
          icon("check-circle", style = "margin-right: 6px;"),
          step$name
        )
      } else {
        # Pending
        div(style = "padding: 4px 0; font-size: 12px; color: rgba(255,255,255,0.5);",
          icon("circle", class = "far", style = "margin-right: 6px;"),
          step$name
        )
      }
    })

    do.call(tagList, step_items)
  })

  # System information
  output$systemInfo <- renderUI({
    tagList(
      div(
        style = "font-size: 11px; line-height: 1.3;",
        span(icon("memory"), " ", values$memory_usage),
        br(),
        span(icon("clock"), " ", 
             format(difftime(Sys.time(), values$analysis_start_time), digits = 2)),
        br(),
        if (values$error_count > 0) {
          span(icon("exclamation-triangle"), " ", values$error_count, " errors")
        }
      )
    )
  })

  output$scalableInfo <- renderUI({
    mode <- if (isTRUE(values$scalable_mode)) "ON" else "OFF"
    color <- if (isTRUE(values$scalable_mode)) "#2ecc71" else "#f39c12"
    ver <- as.integer(values$dataset_version %||% 0L)
    tagList(
      div(
        style = "font-size: 11px; line-height: 1.3;",
        span(icon("tachometer-alt"), " Scalable mode: ",
          span(style = paste0("color:", color, "; font-weight:600;"), mode)
        ),
        br(),
        span(icon("fingerprint"), " Dataset version: ", ver)
      )
    )
  })

  output$jobsInfo <- renderUI({
    jobs <- values$job_queue %||% list()
    active_jobs <- sum(vapply(jobs, function(x) identical(x$status, "running"), logical(1)))
    total_jobs <- length(jobs)
    last_status <- if (total_jobs > 0) as.character(tail(jobs, 1)[[1]]$status %||% "unknown") else "none"
    tagList(
      div(
        style = "font-size: 11px; line-height: 1.3;",
        span(icon("tasks"), " Jobs: ", active_jobs, " running / ", total_jobs, " total"),
        br(),
        span(icon("history"), " Last status: ", last_status)
      )
    )
  })

  # Update memory usage periodically
  observe({
    invalidateLater(30000)  # Update every 30 seconds
    values$memory_usage <- get_memory_usage()
  })

  # ===================== SESSION MANAGEMENT ============================
  
  # ── Helper: build serialisable session_data from current values ─────────────
  .build_session_data <- function() {
    .prep_for_save <- function(obj) {
      if (is.null(obj)) return(NULL)
      tryCatch({
        for (assay_name in names(obj@assays)) {
          a <- obj[[assay_name]]
          if (inherits(a, "Assay5") && length(SeuratObject::Layers(a)) > 1)
            obj[[assay_name]] <- JoinLayers(a)
        }
        obj
      }, error = function(e) obj)
    }
    list(
      seurat_list      = lapply(values$seurat_list, .prep_for_save),
      active_dataset   = .prep_for_save(values$active_dataset),
      active_name      = values$active_name,
      cluster_colors   = values$cluster_colors,
      dataset_version  = values$dataset_version,
      dataset_signature = values$dataset_signature,
      dataset_history  = values$dataset_history,
      analysis_audit   = values$analysis_audit,
      job_queue        = values$job_queue,
      method_inventory = values$method_inventory,
      analysis_flags   = list(
        data_loaded        = values$data_loaded,
        qc_done            = values$qc_done,
        normalization_done = values$normalization_done,
        preprocessing_done = values$preprocessing_done,
        features_selected  = values$features_selected,
        scaling_done       = values$scaling_done,
        dimred_done        = values$dimred_done,
        clustering_done    = values$clustering_done,
        annotation_done    = values$annotation_done,
        trajectory_done    = values$trajectory_done
      ),
      project_context  = list(
        project_background  = values$project_background,
        experiment_design   = values$experiment_design,
        biological_question = values$biological_question,
        species_global      = values$species_global,
        tissue_global       = values$tissue_global,
        conditions_global   = values$conditions_global
      ),
      timestamp = Sys.time()
    )
  }

  observeEvent(input$saveSession, {
    # ── Project-aware save ────────────────────────────────────────────────────
    if (!is.null(auth_state$current_project_id)) {
      withProgress(message = "Saving project…", value = 0, {
        tryCatch({
          incProgress(0.3, detail = "Preparing data…")
          sd <- .build_session_data()
          incProgress(0.4, detail = "Writing to project storage…")
          result <- save_project_state(auth_state$current_project_id,
                                       auth_state$user_id, sd)
          gc()
          incProgress(0.3, detail = "Done!")
          if (isTRUE(result$success)) {
            showNotification(paste0("Project \"", auth_state$current_project_name,
                                    "\" saved."), type = "message", duration = 4)
          } else {
            showNotification(paste("Save failed:", result$message), type = "error")
          }
        }, error = function(e) {
          showNotification(paste("Failed to save:", e$message), type = "error", duration = 10)
        })
      })
    } else {
      # ── Legacy file-system save (no project context) ──────────────────────
      showModal(modalDialog(
        title = "Save Session",
        textInput("sessionName", "Session name:",
                  value = paste0("session_", format(Sys.Date(), "%Y%m%d"))),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("confirmSave", "Save", class = "btn-primary")
        )
      ))
    }
  })

  observeEvent(input$confirmSave, {
    if (nzchar(input$sessionName %||% "")) {
      removeModal()
      withProgress(message = "Saving session...", value = 0, {
        tryCatch({
          incProgress(0.3, detail = "Preparing data...")
          sd <- .build_session_data()
          incProgress(0.4, detail = "Writing to disk...")
          sessions_dir <- file.path(getwd(), "sessions")
          if (!dir.exists(sessions_dir)) dir.create(sessions_dir, recursive = TRUE)
          filepath <- file.path(sessions_dir, paste0(input$sessionName, ".rds"))
          saveRDS(sd, filepath, compress = "gzip")
          rm(sd); gc()
          incProgress(0.3, detail = "Done!")
          showNotification(paste("Session saved:", filepath), type = "message")
        }, error = function(e) {
          showNotification(paste("Failed to save session:", e$message), type = "error", duration = 10)
        })
      })
    }
  })
  
  observeEvent(input$loadSession, {
    showModal(modalDialog(
      title = "Load Session",
      fileInput("sessionFile", "Choose session file (.rds):", accept = ".rds"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirmLoad", "Load", class = "btn-primary")
      )
    ))
  })
  
  observeEvent(input$confirmLoad, {
    req(input$sessionFile)
    removeModal()

    withProgress(message = "Loading session...", value = 0, {
      tryCatch({
        incProgress(0.2, detail = "Reading file...")
        session_data <- readRDS(input$sessionFile$datapath)

        incProgress(0.3, detail = "Restoring datasets...")
        # Restore core data
        if (!is.null(session_data$seurat_list)) {
          values$seurat_list <- session_data$seurat_list
        }
        if (!is.null(session_data$active_dataset)) {
          values$active_dataset <- session_data$active_dataset
        }
        if (!is.null(session_data$active_name)) {
          values$active_name <- session_data$active_name
        }
        if (!is.null(session_data$cluster_colors)) {
          values$cluster_colors <- session_data$cluster_colors
        }
        if (!is.null(session_data$dataset_version)) {
          values$dataset_version <- session_data$dataset_version
        }
        if (!is.null(session_data$dataset_signature)) {
          values$dataset_signature <- session_data$dataset_signature
        }
        if (!is.null(session_data$dataset_history)) {
          values$dataset_history <- session_data$dataset_history
        }
        if (!is.null(session_data$analysis_audit)) {
          values$analysis_audit <- session_data$analysis_audit
        }
        if (!is.null(session_data$job_queue)) {
          values$job_queue <- session_data$job_queue
        }
        if (!is.null(session_data$method_inventory)) {
          values$method_inventory <- session_data$method_inventory
        }

        # Restore project context
        if (!is.null(session_data$project_context)) {
          for (key in names(session_data$project_context)) {
            values[[key]] <- session_data$project_context[[key]] %||% ""
          }
          # Project Setup tab removed — context now comes from project creation.
        }

        incProgress(0.3, detail = "Restoring analysis state...")
        # Restore flags
        if (!is.null(session_data$analysis_flags)) {
          for (flag in names(session_data$analysis_flags)) {
            if (!is.null(session_data$analysis_flags[[flag]])) {
              values[[flag]] <- session_data$analysis_flags[[flag]]
            }
          }
        }

        # Update dataset selector if seurat_list has names
        if (length(values$seurat_list) > 0) {
          updateSelectInput(session, "active_dataset_select",
                            choices = names(values$seurat_list),
                            selected = values$active_name)
        }

        incProgress(0.2, detail = "Done!")
        gc()
        showNotification(paste("Session loaded successfully. Saved:",
                               session_data$timestamp), type = "message")

      }, error = function(e) {
        showNotification(paste("Failed to load session:", e$message), type = "error", duration = 10)
      })
    })
  })
  
  observeEvent(input$resetSession, {
    .reset_values()
    showNotification("Session reset", type = "warning")
  })

  # ===================== CONSOLE LOG CAPTURE ========================

  # Console log capture
  log_buffer <- reactiveVal(character(0))

  # Capture messages
  observe({
    invalidateLater(2000)
    # Capture recent R messages
    msgs <- log_buffer()
    timestamp <- format(Sys.time(), "[%H:%M:%S]")

    # Add system status periodically
    if (length(msgs) == 0 || length(msgs) %% 10 == 0) {
      new_msg <- paste(timestamp, "System ready. Memory:", get_memory_usage())
      msgs <- c(msgs, new_msg)
    }

    # Keep last 100 messages
    if (length(msgs) > 100) msgs <- tail(msgs, 100)
    log_buffer(msgs)
  })

  output$consoleLog <- renderText({
    paste(log_buffer(), collapse = "\n")
  })

  # ===================== MODULE SERVERS =============================
  
  # Update data_loaded flag when dataset changes
  observe({
    if (!is.null(values$active_dataset) && length(values$seurat_list) > 0) {
      values$data_loaded <- TRUE
    } else {
      values$data_loaded <- FALSE
    }
  })

  # Guided workflow prerequisites for major tabs
  last_valid_tab <- reactiveVal("data")
  observeEvent(input$tabs, {
    req(input$tabs)
    requested <- input$tabs
    obj <- values$active_dataset
    has_pca <- !is.null(obj) && inherits(obj, "Seurat") && ("pca" %in% names(obj@reductions))

    gate <- switch(
      requested,
      "qc" = list(ok = isTRUE(values$data_loaded), msg = "Load data before running quality control."),
      "preprocess" = list(ok = isTRUE(values$data_loaded), msg = "Load data before preprocessing."),
      "dimred" = list(ok = isTRUE(values$normalization_done) || isTRUE(values$data_loaded), msg = "Run normalization first for dimensionality reduction."),
      "cluster" = list(ok = has_pca, msg = "Run PCA before clustering."),
      "annotate" = list(ok = isTRUE(values$clustering_done), msg = "Run clustering before annotation."),
      "trajectory" = list(ok = has_pca || isTRUE(values$clustering_done), msg = "Run dimensionality reduction/clustering before trajectory analysis."),
      "downstream" = list(ok = isTRUE(values$clustering_done), msg = "Run clustering before downstream analysis."),
      "visualize" = list(ok = isTRUE(values$data_loaded), msg = "Load data before visualization."),
      "comparison" = list(ok = length(values$seurat_list) >= 1, msg = "Load at least one dataset for comparison."),
      "network" = list(ok = isTRUE(values$clustering_done), msg = "Run clustering/annotation before network analysis."),
      "export" = list(ok = isTRUE(values$data_loaded), msg = "Load data before export."),
      list(ok = TRUE, msg = NULL)
    )

    if (!isTRUE(gate$ok)) {
      showNotification(gate$msg, type = "warning", duration = 6)
      updateTabItems(session, "tabs", selected = last_valid_tab())
    } else {
      last_valid_tab(requested)
    }
  }, ignoreInit = TRUE)

  # AI is configured at project creation and can be edited later from the
  # dashboard's "Edit Project" (pencil) modal, so no in-app AI Settings tab.

  # Project Setup tab removed — biological context + AI setup are collected
  # in the Create Project dialog and loaded into `values` when a project opens.

  # Data Management modules
  if (modules_loaded$data_input) {
    dataUploadServer("dataUpload", values)
    batchUploadServer("batchUpload", values)
    dataManagementServer("dataManage", values)
    dataValidationServer("dataValidate", values)
  }

  # Quality Control modules
  if (modules_loaded$quality_control) {
    qcMetricsServer("qcMetrics", values)
    qcFilteringServer("qcFilter", values)
    doubletDetectionServer("doublet", values)
    ambientRNAServer("ambient", values)
    cellCycleServer("cellCycle", values)
    batchEffectsServer("batchEffects", values)
    qcReportServer("qcReport", values)
  }

  # Preprocessing modules
  if (modules_loaded$preprocessing) {
    normalizationServer("normalize", values)
    featureSelectionServer("features", values)
    scalingServer("scale", values)
    integrationServer("integrate", values)
  }

  # Dimensionality Reduction modules
  if (modules_loaded$dimensionality_reduction) {
    pcaServer("pca", values)
    umapTsneServer("umapTsne", values)
    embedding3dServer("embed3d", values)
    diffusionServer("diffusion", values)
  }

  # Clustering modules
  if (modules_loaded$clustering) {
    graphClusteringServer("graphCluster", values)
    subclusteringServer("subCluster", values)
  }

  # Annotation modules
  if (modules_loaded$annotation) {
    automaticAnnotationServer("autoAnn", values)
    markerValidationServer("markValid", values)
    manualAnnotationServer("manualAnn", values)
    tryCatch(clusterEditorServer("clusterEditor", values), error = function(e) NULL)
    aiAnnotationServer("aiAnn", values)
  }

  # Gene Explorer (Trailmaker-style sub-tab of Annotation)
  if (isTRUE(modules_loaded$gene_explorer)) {
    tryCatch(geneExplorerServer("geneExp", values),
             error = function(e) warning("Gene Explorer init: ", e$message))
  }

  # Trajectory analysis (lazy: only initialize when tab visited)
  .traj_init <- reactiveVal(FALSE)
  observe({
    req(input$tabs == "trajectory", !.traj_init())
    if (modules_loaded$trajectory) {
      trajectoryServer("trajectory", values)
    }
    .traj_init(TRUE)
  })

  # Downstream analysis modules (lazy) — Phase 1: DE only.
  # Pathway/CCC/TF/GRN/Geneset/Composition/AI servers will be re-mounted under
  # the new pathway_network module in Phase 2/3.
  .downstream_init <- reactiveVal(FALSE)
  observe({
    req(input$tabs == "downstream", !.downstream_init())
    if (modules_loaded$de) {
      deServer("de", values)
    }
    .downstream_init(TRUE)
  })

  # Visualization modules (lazy)
  .viz_init <- reactiveVal(FALSE)
  observe({
    req(input$tabs == "visualize", !.viz_init())
    if (modules_loaded$visualization) {
      interactivePlotsServer("viz", values)
      tryCatch(expressionExplorerServer("exprExplorer", values), error = function(e) NULL)
      tryCatch(heatmapBuilderServer("heatmapBuilder", values), error = function(e) NULL)
      tryCatch(clusterComparisonServer("clusterComp", values), error = function(e) NULL)
      tryCatch(figurePanelServer("figPanel", values), error = function(e) NULL)
    }
    .viz_init(TRUE)
  })

  # Comparison View module (lazy)
  .comp_init <- reactiveVal(FALSE)
  observe({
    req(input$tabs == "comparison", !.comp_init())
    if (modules_loaded$visualization) {
      tryCatch(comparisonViewServer("comparisonView", values), error = function(e) NULL)
    }
    .comp_init(TRUE)
  })

  # UMAP Gallery module
  if (isTRUE(modules_loaded$umap_gallery)) {
    tryCatch(umapGalleryServer("umapGallery", values), error = function(e) NULL)
  }

  # Pathway & Network analysis (lazy) — replaces the old network_analysis module.
  .pn_init <- reactiveVal(FALSE)
  observe({
    req(input$tabs == "pathway_network", !.pn_init())
    if (modules_loaded$pn) {
      pnServer("pn", values)
    }
    .pn_init(TRUE)
  })

  # Export modules (lazy)
  .export_init <- reactiveVal(FALSE)
  observe({
    req(input$tabs == "export", !.export_init())
    if (modules_loaded$export) {
      dataExportServer("exporter", values)
    }
    .export_init(TRUE)
  })

  # ===================== DATASET SYNC ==============================

  # Keep seurat_list in sync when modules modify active_dataset
  observe({
    req(values$active_dataset, values$active_name)
    if (!is.null(values$active_name) && values$active_name %in% names(values$seurat_list)) {
      values$seurat_list[[values$active_name]] <- values$active_dataset
    }
  })

  observeEvent(values$active_dataset, {
    obj <- values$active_dataset
    req(obj)
    if (!inherits(obj, "Seurat")) return()

    if (exists("update_scalable_mode", mode = "function")) {
      update_scalable_mode(values, obj)
    }
    if (exists("register_dataset_checkpoint", mode = "function")) {
      register_dataset_checkpoint(
        values,
        module = "app",
        action = "active_dataset_update",
        obj = obj,
        params = list(active_name = values$active_name)
      )
    }
  }, ignoreInit = TRUE)

  # Pre-compute embedding data frames for faster plot rendering
  observeEvent(values$active_dataset, {
    obj <- values$active_dataset
    req(obj)
    if (!inherits(obj, "Seurat")) return()
    reds <- names(obj@reductions)
    cached <- list()
    for (r in reds) {
      tryCatch({
        emb <- Seurat::Embeddings(obj, r)
        if (ncol(emb) >= 2) {
          cached[[r]] <- data.frame(
            Dim1 = emb[, 1],
            Dim2 = emb[, 2],
            Cell = rownames(emb),
            stringsAsFactors = FALSE
          )
        }
      }, error = function(e) NULL)
    }
    values$cached_embeddings <- cached
  }, ignoreInit = TRUE)

  # ===================== FINAL SETUP ================================

  # Show welcome message to new users
  observe({
    if (is.null(values$active_dataset) && length(values$seurat_list) == 0) {
      runjs("$('#welcome-panel').show();")
    } else {
      runjs("$('#welcome-panel').hide();")
    }
  })
  
  # Validate session periodically
  observe({
    invalidateLater(10000)  # Check every 10 seconds
    validation <- validate_session(values)
    if (!validation$valid) {
      # Could show warnings or suggestions here
    }

    if (!is.null(values$active_dataset) && exists("validate_seurat_contract", mode = "function")) {
      contract <- validate_seurat_contract(values$active_dataset, stage = "generic", min_cells = 1, min_features = 1)
      if (!isTRUE(contract$ok) && exists("append_audit_log", mode = "function")) {
        append_audit_log(
          values,
          module = "app",
          action = "periodic_contract_check_failed",
          level = "warning",
          details = list(errors = contract$errors, warnings = contract$warnings)
        )
      }
    }
  })
  
  # ── Periodic autosave ─────────────────────────────────────────────────────
  # Every 3 minutes, if a project is open and something changed since the last
  # save, persist the full state silently. Protects against browser crashes /
  # accidental tab closes without waiting for an explicit Switch/Sign-out.
  observe({
    invalidateLater(180000)  # 3 minutes
    isolate({
      if (is.null(auth_state$current_project_id)) return()
      sig <- .state_signature()
      if (identical(sig, last_saved_signature())) return()  # nothing new
      tryCatch({
        save_project_state(auth_state$current_project_id, auth_state$user_id,
                           .collect_session_data())
        last_saved_signature(sig)
        message("[Autosave] Project ", auth_state$current_project_id, " saved.")
      }, error = function(e) message("[Autosave] failed: ", conditionMessage(e)))
    })
  })

  # ── Clean up + final save on session end (browser close / disconnect) ─────
  session$onSessionEnded(function() {
    plan(sequential)  # Reset parallel processing
    # Persist whatever state exists so a closed tab never loses progress.
    isolate({
      pid <- auth_state$current_project_id
      uid <- auth_state$user_id
      if (!is.null(pid) && !is.null(uid)) {
        tryCatch(save_project_state(pid, uid, .collect_session_data()),
                 error = function(e) NULL)
      }
    })
  })
}

# =============================================================================
# RUN APPLICATION
# =============================================================================

# Add error handling for app startup
tryCatch({
  shinyApp(ui = ui, server = server)
}, error = function(e) {
  cat("Failed to start application:", e$message, "\n")
  stop(e)
})
