# Data Input Module (improved - non-freezing, better batch upload)
# modules/data_input.R

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(readxl)
  library(future)
  library(promises)
})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R

# ============================== UI: Upload ====================================

dataUploadUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Upload Single-Cell Data"),
          helpText(paste("Maximum file size:", round(getOption("shiny.maxRequestSize", 5*1024^2)/(1024^2), 0), "MB per file")),
          
          # Single/Multiple file upload
          fileInput(ns("files"), "Choose files to upload",
            multiple = TRUE,
            accept = c(".h5",".h5ad",".h5seurat",".rds",".csv",".tsv",".txt",
                       ".xlsx",".xls",".mtx",".gz",".zip",".tar",".loom")),
          
          fluidRow(
            column(4,
              selectInput(ns("format"), "File Format",
                choices = c(
                  "Auto-detect"      = "auto",
                  "10X H5"           = "10x_h5", 
                  "10X MTX folder"   = "10x_mtx",
                  "H5AD (AnnData)"   = "h5ad",
                  "H5Seurat"         = "h5seurat",
                  "Seurat RDS"       = "seurat_rds",
                  "CSV/TSV matrix"   = "csv_tsv",
                  "Excel"            = "excel",
                  "Loom"             = "loom",
                  "MTX (sparse trio)"= "mtx"
                )
              )
            ),
            column(4,
              selectInput(ns("species"), "Species",
                choices = c("Human"="human","Mouse"="mouse","Other"="other"),
                selected = "human")
            ),
            column(4,
              div(style="margin-top: 25px;",
                checkboxInput(ns("preview_before_load"), "Preview before loading", FALSE)
              )
            )
          ),
          
          actionButton(ns("upload"), "Load Selected Files",
            class = "btn-primary btn-lg", icon = icon("upload"), width="100%"),
          
          hr(),
          
          # Batch options
          h5("Batch Processing Options"),
          fluidRow(
            column(6,
              checkboxInput(ns("auto_combine"), "Combine all uploaded files", value = FALSE)
            ),
            column(6,
              selectInput(ns("auto_combine_method"), "Combine method",
                choices = c("Merge (append cells)" = "merge",
                            "Harmony Integration" = "harmony",
                            "CCA Integration" = "cca",
                            "RPCA Integration" = "rpca",
                            "FastMNN Integration" = "fastmnn"),
                selected = "merge")
            )
          ),
          textInput(ns("auto_combine_name"), "Combined dataset name", value = "combined_dataset")
        )
      )
    ),
    
    # Progress and status
    fluidRow(
      column(12,
        conditionalPanel(
          condition = sprintf("output['%s']", ns("show_progress")),
          ns = ns,
          wellPanel(
            h4("Upload Progress"),
            verbatimTextOutput(ns("upload_progress")),
            progressBar(ns("progress_bar"), value = 0, striped = TRUE)
          )
        ),
        uiOutput(ns("upload_status"))
      )
    )
  )
}

# ============================ UI: Batch Upload =================================

batchUploadUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Batch File Upload"),
          p("Upload multiple files or a ZIP archive containing your data files."),
          
          tabsetPanel(
            # Multi-file upload
            tabPanel("Multiple Files", icon = icon("copy"),
              br(),
              fileInput(ns("batch_files"), "Select multiple data files",
                multiple = TRUE,
                accept = c(".h5",".h5ad",".h5seurat",".rds",".csv",".tsv",".txt",
                          ".xlsx",".xls",".mtx",".loom")),
              helpText("Hold Ctrl/Cmd to select multiple files from your computer."),
              actionButton(ns("load_batch_files"), "Load All Selected Files",
                          icon = icon("play"), class = "btn-success btn-lg", width = "100%")
            ),
            
            # ZIP archive upload  
            tabPanel("ZIP Archive",    icon = icon("file-archive"),
              br(),
              fileInput(ns("zip_file"), "Upload ZIP archive",
                accept = c(".zip", ".tar", ".tar.gz")),
              helpText("Upload a ZIP or TAR archive containing your data files. All supported files will be extracted and loaded."),
              actionButton(ns("load_zip"), "Extract and Load Archive",
                          icon = icon("file-archive"), class = "btn-info btn-lg", width = "100%")
            ),
            
            # 10X folder structure
            tabPanel("10X Folder",     icon = icon("folder-open"),
              br(),
              p("For 10X Genomics data: upload the 3 required files from your cellranger output."),
              fileInput(ns("tenx_matrix"), "matrix.mtx[.gz]", accept = c(".mtx", ".gz")),
              fileInput(ns("tenx_barcodes"), "barcodes.tsv[.gz]", accept = c(".tsv", ".gz")),
              fileInput(ns("tenx_features"), "features.tsv[.gz] or genes.tsv[.gz]", accept = c(".tsv", ".gz")),
              textInput(ns("tenx_name"), "Dataset name", value = "10X_dataset"),
              actionButton(ns("load_tenx"), "Load 10X Dataset",
                          icon = icon("folder-open"), class = "btn-warning btn-lg", width = "100%")
            )
          )
        )
      )
    ),
    
    fluidRow(
      column(12,
        conditionalPanel(
          condition = sprintf("output['%s']", ns("show_batch_progress")),
          ns = ns,
          wellPanel(
            h4("Batch Processing Progress"),
            verbatimTextOutput(ns("batch_progress")),
            progressBar(ns("batch_progress_bar"), value = 0, striped = TRUE)
          )
        ),
        uiOutput(ns("batch_status"))
      )
    )
  )
}

# ======================== UI: Manage datasets ==================================

dataManagementUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Dataset Management"),

          # Dataset overview
          conditionalPanel(
            condition = "true", # Always show
            DT::dataTableOutput(ns("dataset_table"))
          ),

          br(),

          # Management actions
          fluidRow(
            column(6,
              wellPanel(
                h5("Rename Dataset"),
                selectInput(ns("rename_from"), "Select dataset to rename", choices = NULL),
                textInput(ns("rename_to"), "New name", ""),
                actionButton(ns("do_rename"), "Rename", class = "btn-warning", width = "100%")
              )
            ),
            column(6,
              wellPanel(
                h5("Delete Dataset"),
                selectInput(ns("delete_dataset"), "Select dataset to delete", choices = NULL),
                actionButton(ns("do_delete"), "Delete", class = "btn-danger", width = "100%",
                           onclick = "return confirm('Are you sure you want to delete this dataset?');")
              )
            )
          ),

          hr(),

          # Combining datasets
          wellPanel(
            h5("Combine Selected Datasets"),
            helpText("Select 2 or more datasets in the table above, then configure combination settings below."),
            fluidRow(
              column(4,
                selectInput(ns("combine_method"), "Combine method",
                  choices = c("Merge (append cells)" = "merge",
                              "Harmony Integration" = "harmony",
                              "CCA Integration" = "cca",
                              "RPCA Integration" = "rpca",
                              "FastMNN Integration" = "fastmnn"),
                  selected = "merge")
              ),
              column(4,
                textInput(ns("combine_name"), "Combined dataset name", "combined")
              ),
              column(4,
                div(style="margin-top: 25px;",
                  actionButton(ns("combine_selected"), "Combine Selected",
                             class = "btn-success", width = "100%")
                )
              )
            ),

            hr(),

            # --- Merge options panel (visible when 2+ datasets are selected) ---
            uiOutput(ns("merge_options_panel"))
          ),

          # --- Merge preview panel ---
          uiOutput(ns("merge_preview_panel"))
        )
      )
    )
  )
}

# ======================== UI: Data validation ==================================

dataValidationUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Data Quality Validation"),
          p("Comprehensive validation of loaded datasets to ensure data integrity and compatibility."),
          
          fluidRow(
            column(6,
              selectInput(ns("validate_dataset"), "Select dataset to validate", choices = NULL)
            ),
            column(6,
              div(style="margin-top: 25px;",
                actionButton(ns("run_validation"), "Run Validation",
                           class = "btn-primary", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    fluidRow(
      column(12,
        conditionalPanel(
          condition = sprintf("output['%s']", ns("show_validation_results")),
          ns = ns,
          wellPanel(
            h4("Validation Results"),
            verbatimTextOutput(ns("validation_log")),
            br(),
            DT::dataTableOutput(ns("validation_details"))
          )
        )
      )
    )
  )
}

# ============================== Server: Upload =================================

dataUploadServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    # Progress tracking
    progress_val <- reactiveVal(0)
    progress_msg <- reactiveVal("")
    show_progress <- reactiveVal(FALSE)
    
    output$show_progress <- reactive({ show_progress() })
    outputOptions(output, "show_progress", suspendWhenHidden = FALSE)
    
    output$upload_progress <- renderText({ progress_msg() })
    output$progress_bar <- renderUI({
      progressBar("progress", value = progress_val(), striped = TRUE)
    })
    
    # Main upload handler
    observeEvent(input$upload, {
      req(input$files)

      tryCatch({
      show_progress(TRUE)
      progress_val(0)
      progress_msg("Initializing upload...")
      
      files_df <- input$files
      total_files <- nrow(files_df)
      loaded_datasets <- character(0)
      
      # Process files sequentially to avoid freezing
      for (i in seq_len(total_files)) {
        file_info <- files_df[i, ]
        
        # Update progress
        progress_val(round((i-1)/total_files * 100))
        progress_msg(paste("Processing file", i, "of", total_files, ":", file_info$name))
        
        # Allow UI to update
        Sys.sleep(0.1)
        
        # Load file
        result <- safely_load_file(file_info, input$format, input$species, input$preview_before_load)
        
        if (!is.null(result$data)) {
          # Generate unique name
          base_name <- tools::file_path_sans_ext(file_info$name)
          unique_name <- make_unique_name(names(values$seurat_list), base_name)

          # Optionally move large datasets to on-disk (BPCells) storage to cut RAM.
          # Env-gated (SCRNA_ONDISK) and fully guarded — a no-op / in-memory
          # fallback when disabled, small, or on any failure.
          obj <- apply_ondisk_if_enabled(
            result$data,
            file.path(values$ondisk_dir %||% file.path(tempdir(), "ondisk"), unique_name))

          # Store dataset
          values$seurat_list[[unique_name]] <- obj
          loaded_datasets <- c(loaded_datasets, unique_name)

          # Set as active if first successful load
          if (is.null(values$active_dataset)) {
            values$active_dataset <- obj
            values$active_name <- unique_name
          }

          showNotification(paste("Loaded:", file_info$name, "→", unique_name), type = "message")
        } else {
          showNotification(paste("Failed:", file_info$name, "-", result$error), type = "error")
        }
      }
      
      # Final progress
      progress_val(100)
      progress_msg(paste("Complete! Loaded", length(loaded_datasets), "of", total_files, "files."))
      
      # Auto-combine if requested
      if (input$auto_combine && length(loaded_datasets) >= 2) {
        progress_msg("Combining datasets...")

        combine_name <- make_unique_name(names(values$seurat_list), input$auto_combine_name)

        result <- safely_combine_objects(
          values$seurat_list[loaded_datasets],
          input$auto_combine_method,
          combine_name,
          add_batch_label = TRUE,
          batch_col_name = "batch",
          progress_callback = function(msg, val = NULL) {
            progress_msg(paste("Combining:", msg))
            if (!is.null(val)) progress_val(val)
          }
        )

        if (!is.null(result$data)) {
          obj <- apply_ondisk_if_enabled(
            result$data,
            file.path(values$ondisk_dir %||% file.path(tempdir(), "ondisk"), combine_name))
          values$seurat_list[[combine_name]] <- obj
          values$active_dataset <- obj
          values$active_name <- combine_name
          showNotification(
            paste("Combined", length(loaded_datasets), "datasets ->", combine_name,
                  "(", ncol(result$data), "cells)"),
            type = "message", duration = 8
          )
        } else {
          showNotification(paste("Auto-combine failed:", result$error), type = "error", duration = 10)
        }
      }
      
      # Hide progress after delay
      later::later(function() { show_progress(FALSE) }, 3)
      }, error = function(e) {
        show_progress(FALSE)
        showNotification(paste("Upload error:", e$message), type = "error", duration = 10)
      })
    })
    
    # Upload status display
    output$upload_status <- renderUI({
      if (length(values$seurat_list) > 0) {
        tagList(
          h4("Loaded Datasets"),
          DT::dataTableOutput(session$ns("loaded_datasets_table"))
        )
      }
    })
    
    output$loaded_datasets_table <- DT::renderDataTable({
      req(length(values$seurat_list) > 0)
      
      create_dataset_summary_table(values$seurat_list, values$active_name)
    })
  })
}

# ============================= Server: Batch Upload ========================

batchUploadServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Batch progress tracking
    batch_progress_val <- reactiveVal(0)
    batch_progress_msg <- reactiveVal("")
    show_batch_progress <- reactiveVal(FALSE)

    # Track files loaded in the most recent batch for auto-combine
    last_batch_loaded <- reactiveVal(character(0))

    output$show_batch_progress <- reactive({ show_batch_progress() })
    outputOptions(output, "show_batch_progress", suspendWhenHidden = FALSE)

    output$batch_progress <- renderText({ batch_progress_msg() })
    output$batch_progress_bar <- renderUI({
      progressBar("batch_progress", value = batch_progress_val(), striped = TRUE)
    })

    # Multi-file batch upload
    observeEvent(input$load_batch_files, {
      req(input$batch_files)
      tryCatch({
        show_batch_progress(TRUE)
        batch_progress_val(0)
        batch_progress_msg("Starting batch upload...")

        loaded_names <- process_multiple_files(input$batch_files, values,
                                               batch_progress_val, batch_progress_msg)
        last_batch_loaded(loaded_names)

        if (length(loaded_names) >= 2) {
          batch_progress_msg(paste("Batch complete! Loaded", length(loaded_names),
                                   "datasets. Ready for merging from Data Management tab."))
        }

        later::later(function() { show_batch_progress(FALSE) }, 4)
      }, error = function(e) {
        show_batch_progress(FALSE)
        showNotification(paste("Batch upload error:", e$message), type = "error", duration = 10)
      })
    })

    # ZIP archive upload
    observeEvent(input$load_zip, {
      req(input$zip_file)
      tryCatch({
        show_batch_progress(TRUE)
        batch_progress_val(0)
        batch_progress_msg("Extracting archive...")

        loaded_names <- extract_and_process_archive(input$zip_file$datapath, values,
                                                    batch_progress_val, batch_progress_msg)
        last_batch_loaded(if (is.character(loaded_names)) loaded_names else character(0))

        later::later(function() { show_batch_progress(FALSE) }, 4)
      }, error = function(e) {
        show_batch_progress(FALSE)
        showNotification(paste("Archive upload error:", e$message), type = "error", duration = 10)
      })
    })

    # 10X dataset upload (supports .gz compressed files)
    observeEvent(input$load_tenx, {
      req(input$tenx_matrix, input$tenx_barcodes, input$tenx_features)
      tryCatch({
        show_batch_progress(TRUE)
        batch_progress_val(10)
        batch_progress_msg("Processing 10X dataset...")

        result <- process_10x_files(
          input$tenx_matrix$datapath,
          input$tenx_barcodes$datapath,
          input$tenx_features$datapath,
          input$tenx_name,
          values
        )

        if (result$success) {
          batch_progress_val(100)
          batch_progress_msg("10X dataset loaded successfully.")
          showNotification("10X dataset loaded successfully", type = "message")
        } else {
          batch_progress_msg(paste("10X loading failed:", result$error))
          showNotification(paste("10X loading failed:", result$error), type = "error")
        }

        later::later(function() { show_batch_progress(FALSE) }, 3)
      }, error = function(e) {
        show_batch_progress(FALSE)
        showNotification(paste("10X upload error:", e$message), type = "error", duration = 10)
      })
    })

    # Batch status display
    output$batch_status <- renderUI({
      if (length(values$seurat_list) > 0) {
        tagList(
          h4("All Loaded Datasets"),
          DT::dataTableOutput(session$ns("all_datasets_table"))
        )
      }
    })

    output$all_datasets_table <- DT::renderDataTable({
      req(length(values$seurat_list) > 0)
      create_dataset_summary_table(values$seurat_list, values$active_name)
    })
  })
}

# ======================= Server: Dataset Management ========================

dataManagementServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Update choices when datasets change
    observe({
      dataset_names <- names(values$seurat_list)
      updateSelectInput(session, "rename_from", choices = dataset_names)
      updateSelectInput(session, "delete_dataset", choices = dataset_names)
    })

    # Dataset table
    output$dataset_table <- DT::renderDataTable({
      req(length(values$seurat_list) > 0)
      create_dataset_summary_table(values$seurat_list, values$active_name)
    })

    # Set active dataset by clicking a single row
    observeEvent(input$dataset_table_rows_selected, {
      selected_rows <- input$dataset_table_rows_selected
      if (length(selected_rows) == 1) {
        dataset_names <- names(values$seurat_list)
        selected_name <- dataset_names[selected_rows]
        values$active_dataset <- values$seurat_list[[selected_name]]
        values$active_name <- selected_name
        showNotification(paste("Active dataset:", selected_name), type = "message")
      }
    })

    # Rename dataset
    observeEvent(input$do_rename, {
      req(input$rename_from, input$rename_to)
      old_name <- input$rename_from
      new_name <- trimws(input$rename_to)

      if (nzchar(new_name) && old_name %in% names(values$seurat_list)) {
        # Ensure unique name
        new_name <- make_unique_name(names(values$seurat_list), new_name)

        # Rename
        values$seurat_list[[new_name]] <- values$seurat_list[[old_name]]
        values$seurat_list[[old_name]] <- NULL

        # Update active if needed
        if (identical(values$active_name, old_name)) {
          values$active_name <- new_name
        }

        showNotification(paste("Renamed:", old_name, "->", new_name), type = "message")
        updateTextInput(session, "rename_to", value = "")
      }
    })

    # Delete dataset
    observeEvent(input$do_delete, {
      req(input$delete_dataset)
      dataset_to_delete <- input$delete_dataset

      if (dataset_to_delete %in% names(values$seurat_list)) {
        values$seurat_list[[dataset_to_delete]] <- NULL

        # Update active dataset if deleted
        if (identical(values$active_name, dataset_to_delete)) {
          remaining_names <- names(values$seurat_list)
          if (length(remaining_names) > 0) {
            values$active_dataset <- values$seurat_list[[remaining_names[1]]]
            values$active_name <- remaining_names[1]
          } else {
            values$active_dataset <- NULL
            values$active_name <- NULL
          }
        }

        showNotification(paste("Deleted dataset:", dataset_to_delete), type = "warning")
      }
    })

    # ---- Merge options panel (dynamic UI shown when 2+ datasets selected) ----
    output$merge_options_panel <- renderUI({
      selected_indices <- input$dataset_table_rows_selected
      n_datasets <- length(values$seurat_list)

      # Show merge options when 2+ datasets are selected OR 2+ datasets exist
      if (n_datasets < 2 && (is.null(selected_indices) || length(selected_indices) < 2)) {
        return(NULL)
      }

      ns <- session$ns

      # Collect all metadata columns across selected datasets
      if (!is.null(selected_indices) && length(selected_indices) >= 2) {
        dataset_names <- names(values$seurat_list)
        sel_names <- dataset_names[selected_indices]
        all_meta_cols <- unique(unlist(lapply(sel_names, function(nm) {
          colnames(values$seurat_list[[nm]]@meta.data)
        })))
      } else {
        all_meta_cols <- unique(unlist(lapply(values$seurat_list, function(obj) {
          colnames(obj@meta.data)
        })))
      }

      # Remove columns that are always kept internally
      optional_cols <- setdiff(all_meta_cols, c("orig.ident", "nCount_RNA", "nFeature_RNA"))

      tagList(
        h5("Merge Options"),
        fluidRow(
          column(6,
            checkboxInput(ns("add_batch_label"), "Add batch/sample label column", value = TRUE),
            conditionalPanel(
              condition = paste0("input['", ns("add_batch_label"), "']"),
              textInput(ns("batch_col_name"), "Batch column name", value = "batch")
            )
          ),
          column(6,
            checkboxGroupInput(ns("include_meta_cols"), "Additional metadata columns to include",
              choices = optional_cols,
              selected = optional_cols
            ),
            helpText("Essential columns (orig.ident, nCount_RNA, nFeature_RNA) are always included.")
          )
        )
      )
    })

    # ---- Merge preview panel (shows summary of what merged dataset would look like) ----
    output$merge_preview_panel <- renderUI({
      selected_indices <- input$dataset_table_rows_selected

      if (is.null(selected_indices) || length(selected_indices) < 2) {
        return(NULL)
      }

      ns <- session$ns
      dataset_names <- names(values$seurat_list)
      sel_names <- dataset_names[selected_indices]

      # Build a preview summary
      preview_rows <- lapply(sel_names, function(nm) {
        obj <- values$seurat_list[[nm]]
        data.frame(
          Dataset = nm,
          Cells = ncol(obj),
          Genes = nrow(obj),
          stringsAsFactors = FALSE
        )
      })
      preview_df <- do.call(rbind, preview_rows)

      # Compute union gene count
      all_genes <- unique(unlist(lapply(sel_names, function(nm) {
        rownames(values$seurat_list[[nm]])
      })))
      union_gene_count <- length(all_genes)

      # Compute intersection gene count
      gene_sets <- lapply(sel_names, function(nm) rownames(values$seurat_list[[nm]]))
      intersect_gene_count <- length(Reduce(intersect, gene_sets))

      total_cells <- sum(preview_df$Cells)

      wellPanel(
        style = "background-color: #f0f8ff; border: 1px solid #b0d4f1;",
        h5("Merge Preview"),
        fluidRow(
          column(6,
            tags$table(class = "table table-condensed table-bordered",
              style = "margin-bottom: 10px;",
              tags$thead(
                tags$tr(
                  tags$th("Source Dataset"),
                  tags$th("Cells"),
                  tags$th("Genes")
                )
              ),
              tags$tbody(
                lapply(seq_len(nrow(preview_df)), function(r) {
                  tags$tr(
                    tags$td(preview_df$Dataset[r]),
                    tags$td(format(preview_df$Cells[r], big.mark = ",")),
                    tags$td(format(preview_df$Genes[r], big.mark = ","))
                  )
                }),
                tags$tr(style = "font-weight: bold; background-color: #e8f4fd;",
                  tags$td("MERGED TOTAL"),
                  tags$td(format(total_cells, big.mark = ",")),
                  tags$td(paste0(format(union_gene_count, big.mark = ","), " (union)"))
                )
              )
            )
          ),
          column(6,
            tags$ul(
              tags$li(paste0("Total cells after merge: ", format(total_cells, big.mark = ","))),
              tags$li(paste0("Gene union (all genes): ", format(union_gene_count, big.mark = ","))),
              tags$li(paste0("Gene intersection (shared): ", format(intersect_gene_count, big.mark = ","))),
              tags$li(paste0("Number of datasets: ", length(sel_names))),
              if (!is.null(input$combine_method) && input$combine_method != "merge") {
                tags$li(style = "color: #d9534f;",
                        paste0(toupper(input$combine_method), " integration selected - this may take several minutes for large datasets."))
              }
            )
          )
        )
      )
    })

    # ---- Combine selected datasets ----
    observeEvent(input$combine_selected, {
      selected_indices <- input$dataset_table_rows_selected

      if (length(selected_indices) < 2) {
        showNotification("Select at least 2 datasets in the table to combine", type = "error")
        return()
      }

      dataset_names <- names(values$seurat_list)
      selected_datasets <- values$seurat_list[dataset_names[selected_indices]]

      # Determine the combined name, ensure it is unique
      combine_name <- trimws(input$combine_name)
      if (!nzchar(combine_name)) combine_name <- "combined"
      combine_name <- make_unique_name(names(values$seurat_list), combine_name)

      # Read merge options (use defaults if the dynamic UI hasn't rendered yet)
      add_batch <- if (!is.null(input$add_batch_label)) input$add_batch_label else TRUE
      batch_col <- if (!is.null(input$batch_col_name) && nzchar(input$batch_col_name)) {
        input$batch_col_name
      } else {
        "batch"
      }
      include_cols <- input$include_meta_cols  # may be NULL, which is fine

      withProgress(message = "Combining datasets...", value = 0, {
        progress_cb <- function(msg, val = NULL) {
          if (!is.null(val)) setProgress(value = val / 100, message = msg)
          else setProgress(message = msg)
        }

        result <- tryCatch({
          combined <- combine_seurat_objects(
            selected_datasets,
            method = input$combine_method,
            name = combine_name,
            add_batch_label = add_batch,
            batch_col_name = batch_col,
            include_meta_cols = include_cols,
            progress_callback = progress_cb
          )
          list(data = combined, error = NULL)
        }, error = function(e) {
          list(data = NULL, error = e$message)
        })

        if (!is.null(result$data)) {
          obj <- apply_ondisk_if_enabled(
            result$data,
            file.path(values$ondisk_dir %||% file.path(tempdir(), "ondisk"), combine_name))
          values$seurat_list[[combine_name]] <- obj
          values$active_dataset <- obj
          values$active_name <- combine_name
          showNotification(
            paste("Combined", length(selected_datasets), "datasets ->", combine_name,
                  "(", ncol(result$data), "cells,", nrow(result$data), "genes)"),
            type = "message", duration = 8
          )
        } else {
          showNotification(paste("Combination failed:", result$error), type = "error", duration = 10)
        }
      })
    })
  })
}

# ======================= Server: Data Validation ===========================

dataValidationServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    
    show_validation_results <- reactiveVal(FALSE)
    
    output$show_validation_results <- reactive({ show_validation_results() })
    outputOptions(output, "show_validation_results", suspendWhenHidden = FALSE)
    
    # Update dataset choices
    observe({
      dataset_names <- names(values$seurat_list)
      updateSelectInput(session, "validate_dataset", 
                       choices = c("Active Dataset" = "active", dataset_names),
                       selected = "active")
    })
    
    # Run validation
    observeEvent(input$run_validation, {
      dataset_to_validate <- if (input$validate_dataset == "active") {
        values$active_dataset
      } else {
        values$seurat_list[[input$validate_dataset]]
      }
      
      req(dataset_to_validate)
      
      validation_results <- comprehensive_data_validation(dataset_to_validate)
      
      # Display results
      output$validation_log <- renderText({
        paste(validation_results$messages, collapse = "\n")
      })
      
      output$validation_details <- DT::renderDataTable({
        validation_results$details_table
      }, options = list(pageLength = 15, scrollX = TRUE))
      
      show_validation_results(TRUE)
      
      if (validation_results$overall_status == "PASS") {
        showNotification("Dataset validation passed", type = "message")
      } else {
        showNotification("Dataset validation found issues", type = "warning")
      }
    })
  })
}

# ========================= Helper Functions ==========================

safely_load_file <- function(file_info, format, species, preview) {
  tryCatch({
    fmt <- if (format == "auto") detect_file_format(file_info$datapath, file_info$name) else format
    data <- load_data_by_format(file_info$datapath, fmt, species)
    list(data = data, error = NULL)
  }, error = function(e) {
    list(data = NULL, error = e$message)
  })
}

safely_combine_objects <- function(object_list, method, name, add_batch_label = TRUE,
                                   batch_col_name = "batch", include_meta_cols = NULL,
                                   progress_callback = NULL) {
  tryCatch({
    combined <- combine_seurat_objects(
      object_list, method, name,
      add_batch_label = add_batch_label,
      batch_col_name = batch_col_name,
      include_meta_cols = include_meta_cols,
      progress_callback = progress_callback
    )
    list(data = combined, error = NULL)
  }, error = function(e) {
    list(data = NULL, error = e$message)
  })
}

make_unique_name <- function(existing_names, base_name) {
  if (!base_name %in% existing_names) return(base_name)
  
  counter <- 1
  while (paste0(base_name, "_", counter) %in% existing_names) {
    counter <- counter + 1
  }
  paste0(base_name, "_", counter)
}

create_dataset_summary_table <- function(seurat_list, active_name) {
  if (length(seurat_list) == 0) return(data.frame())
  
  summary_data <- do.call(rbind, lapply(names(seurat_list), function(name) {
    obj <- seurat_list[[name]]
    data.frame(
      Dataset = name,
      Cells = ncol(obj),
      Genes = nrow(obj), 
      Assays = paste(names(obj@assays), collapse = ", "),
      Reductions = paste(names(obj@reductions), collapse = ", "),
      Active = ifelse(name == active_name, "✓", ""),
      stringsAsFactors = FALSE
    )
  }))
  
  DT::datatable(summary_data, 
                options = list(pageLength = 10, scrollX = TRUE),
                selection = "multiple")
}

process_multiple_files <- function(files_df, values, progress_val, progress_msg) {
  total_files <- nrow(files_df)
  loaded_count <- 0
  loaded_names <- character(0)

  for (i in seq_len(total_files)) {
    file_info <- files_df[i, ]
    progress_val(round((i-1)/total_files * 100))
    progress_msg(paste("Loading", i, "of", total_files, ":", file_info$name))

    # Allow the UI to update between files
    Sys.sleep(0.1)

    result <- safely_load_file(file_info, "auto", "human", FALSE)

    if (!is.null(result$data)) {
      base_name <- tools::file_path_sans_ext(file_info$name)
      # Strip double extensions like .mtx.gz -> base stays clean
      base_name <- tools::file_path_sans_ext(base_name)
      if (nchar(base_name) == 0) base_name <- paste0("dataset_", i)
      unique_name <- make_unique_name(names(values$seurat_list), base_name)
      values$seurat_list[[unique_name]] <- result$data
      loaded_count <- loaded_count + 1
      loaded_names <- c(loaded_names, unique_name)

      if (is.null(values$active_dataset)) {
        values$active_dataset <- result$data
        values$active_name <- unique_name
      }

      showNotification(paste("Loaded:", file_info$name, "->", unique_name), type = "message")
    } else {
      showNotification(paste("Failed to load:", file_info$name, "-", result$error), type = "error")
    }
  }

  progress_val(100)
  progress_msg(paste("Batch complete! Loaded", loaded_count, "of", total_files, "files."))

  # Return the names of successfully loaded datasets
  invisible(loaded_names)
}

extract_and_process_archive <- function(archive_path, values, progress_val, progress_msg) {
  temp_dir <- tempdir()
  extract_dir <- file.path(temp_dir, paste0("extracted_", as.numeric(Sys.time())))
  dir.create(extract_dir, recursive = TRUE)

  # Ensure cleanup of extracted files when done
  on.exit({
    tryCatch(unlink(extract_dir, recursive = TRUE), error = function(e) NULL)
  }, add = TRUE)

  # Extract archive
  tryCatch({
    if (grepl("\\.zip$", archive_path, ignore.case = TRUE)) {
      utils::unzip(archive_path, exdir = extract_dir)
    } else if (grepl("\\.(tar|tar\\.gz|tgz)$", archive_path, ignore.case = TRUE)) {
      utils::untar(archive_path, exdir = extract_dir)
    } else {
      showNotification("Unsupported archive format. Use .zip, .tar, or .tar.gz", type = "error")
      return(character(0))
    }
  }, error = function(e) {
    showNotification(paste("Failed to extract archive:", e$message), type = "error")
    return(character(0))
  })

  # Find supported files (including .gz compressed variants)
  all_files <- list.files(extract_dir,
                         pattern = "\\.(h5|h5ad|h5seurat|rds|csv|tsv|txt|xlsx|xls|mtx|loom)(\\.gz)?$",
                         full.names = TRUE, recursive = TRUE, ignore.case = TRUE)

  # Check for 10X folder structures inside the archive (matrix.mtx + barcodes.tsv + features/genes.tsv)
  mtx_files <- grep("matrix\\.mtx", all_files, value = TRUE, ignore.case = TRUE)
  tenx_dirs <- unique(dirname(mtx_files))

  # Process 10X directories first, then remove individual mtx/tsv files from the general list
  loaded_names <- character(0)
  tenx_component_files <- character(0)

  for (tenx_dir in tenx_dirs) {
    dir_files <- list.files(tenx_dir, full.names = TRUE, ignore.case = TRUE)
    has_matrix <- any(grepl("matrix\\.mtx", dir_files, ignore.case = TRUE))
    has_barcodes <- any(grepl("barcodes\\.tsv", dir_files, ignore.case = TRUE))
    has_features <- any(grepl("(features|genes)\\.tsv", dir_files, ignore.case = TRUE))

    if (has_matrix && has_barcodes && has_features) {
      progress_msg(paste("Detected 10X folder structure in:", basename(tenx_dir)))

      result <- tryCatch({
        counts <- Read10X(tenx_dir)
        ds_name <- basename(tenx_dir)
        if (ds_name == "." || nchar(ds_name) == 0) ds_name <- "10X_archive"
        obj <- CreateSeuratObject(counts = counts, project = ds_name, min.cells = 3, min.features = 200)
        list(data = obj, error = NULL)
      }, error = function(e) {
        list(data = NULL, error = e$message)
      })

      if (!is.null(result$data)) {
        ds_name <- basename(tenx_dir)
        if (ds_name == "." || nchar(ds_name) == 0) ds_name <- "10X_archive"
        unique_name <- make_unique_name(names(values$seurat_list), ds_name)
        values$seurat_list[[unique_name]] <- result$data
        loaded_names <- c(loaded_names, unique_name)

        if (is.null(values$active_dataset)) {
          values$active_dataset <- result$data
          values$active_name <- unique_name
        }
        showNotification(paste("Loaded 10X folder:", unique_name), type = "message")
      }

      # Mark these files so we skip them in general processing
      tenx_component_files <- c(tenx_component_files, dir_files)
    }
  }

  # Remove 10X component files from the general list
  all_files <- setdiff(all_files, tenx_component_files)
  # Also remove bare .tsv/.mtx files that were part of 10X dirs
  all_files <- all_files[!all_files %in% list.files(tenx_dirs, full.names = TRUE, recursive = FALSE)]

  if (length(all_files) == 0 && length(loaded_names) == 0) {
    showNotification("No supported files found in archive", type = "error")
    return(character(0))
  }

  # Process remaining files
  for (i in seq_along(all_files)) {
    file_path <- all_files[i]
    file_name <- basename(file_path)
    progress_val(round((i-1)/length(all_files) * 100))
    progress_msg(paste("Processing", i, "of", length(all_files), ":", file_name))

    Sys.sleep(0.1)

    file_info <- list(datapath = file_path, name = file_name)
    result <- safely_load_file(file_info, "auto", "human", FALSE)

    if (!is.null(result$data)) {
      base_name <- tools::file_path_sans_ext(file_name)
      base_name <- tools::file_path_sans_ext(base_name)  # handle double extensions like .mtx.gz
      if (nchar(base_name) == 0) base_name <- paste0("archive_dataset_", i)
      unique_name <- make_unique_name(names(values$seurat_list), base_name)
      values$seurat_list[[unique_name]] <- result$data
      loaded_names <- c(loaded_names, unique_name)

      if (is.null(values$active_dataset)) {
        values$active_dataset <- result$data
        values$active_name <- unique_name
      }
      showNotification(paste("Loaded:", file_name, "->", unique_name), type = "message")
    } else {
      showNotification(paste("Failed:", file_name, "-", result$error), type = "warning")
    }
  }

  progress_val(100)
  progress_msg(paste("Archive processing complete!", length(loaded_names), "datasets loaded."))

  invisible(loaded_names)
}

process_10x_files <- function(matrix_path, barcodes_path, features_path, dataset_name, values) {
  tryCatch({
    # Create temporary directory with 10X structure
    temp_dir <- tempfile()
    dir.create(temp_dir)

    # Ensure cleanup on exit
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

    # Helper: copy file with the correct target name, preserving .gz extension if present
    copy_10x_file <- function(src_path, target_base) {
      src_name <- basename(src_path)
      is_gz <- grepl("\\.gz$", src_name, ignore.case = TRUE)
      target_name <- if (is_gz) paste0(target_base, ".gz") else target_base
      file.copy(src_path, file.path(temp_dir, target_name))
    }

    copy_10x_file(matrix_path, "matrix.mtx")
    copy_10x_file(barcodes_path, "barcodes.tsv")

    # features.tsv or genes.tsv - Seurat Read10X accepts either, but features.tsv is preferred
    feat_src_name <- basename(features_path)
    if (grepl("genes", feat_src_name, ignore.case = TRUE)) {
      copy_10x_file(features_path, "genes.tsv")
    } else {
      copy_10x_file(features_path, "features.tsv")
    }

    # Load using Seurat's Read10X
    counts <- Read10X(temp_dir)

    # If counts is a list (multi-modal), take the Gene Expression component
    if (is.list(counts) && !is.matrix(counts) && !inherits(counts, "dgCMatrix")) {
      if ("Gene Expression" %in% names(counts)) {
        counts <- counts[["Gene Expression"]]
      } else {
        counts <- counts[[1]]
      }
    }

    seurat_obj <- CreateSeuratObject(counts = counts, project = dataset_name,
                                     min.cells = 3, min.features = 200)

    # Store dataset
    unique_name <- make_unique_name(names(values$seurat_list), dataset_name)
    values$seurat_list[[unique_name]] <- seurat_obj

    if (is.null(values$active_dataset)) {
      values$active_dataset <- seurat_obj
      values$active_name <- unique_name
    }

    list(success = TRUE, error = NULL)
  }, error = function(e) {
    list(success = FALSE, error = e$message)
  })
}

comprehensive_data_validation <- function(seurat_obj) {
  messages <- c()
  details <- list()
  status <- "PASS"
  
  # Basic structure checks
  if (nrow(seurat_obj) == 0) {
    messages <- c(messages, "❌ CRITICAL: No genes in dataset")
    status <- "FAIL"
  } else {
    messages <- c(messages, paste("✅ Genes:", nrow(seurat_obj)))
  }
  
  if (ncol(seurat_obj) == 0) {
    messages <- c(messages, "❌ CRITICAL: No cells in dataset")
    status <- "FAIL"
  } else {
    messages <- c(messages, paste("✅ Cells:", ncol(seurat_obj)))
  }
  
  # Assay checks
  assays <- names(seurat_obj@assays)
  messages <- c(messages, paste("✅ Assays:", paste(assays, collapse = ", ")))
  
  if (!"RNA" %in% assays) {
    messages <- c(messages, "⚠️ WARNING: No RNA assay found")
    status <- "WARNING"
  }
  
  # Metadata checks
  required_meta <- c("nCount_RNA", "nFeature_RNA")
  missing_meta <- required_meta[!required_meta %in% colnames(seurat_obj@meta.data)]
  if (length(missing_meta) > 0) {
    messages <- c(messages, paste("⚠️ WARNING: Missing metadata:", paste(missing_meta, collapse = ", ")))
    status <- "WARNING"
  }
  
  # Data quality metrics
  if ("nCount_RNA" %in% colnames(seurat_obj@meta.data)) {
    count_stats <- summary(seurat_obj$nCount_RNA)
    details$nCount_stats <- count_stats
    messages <- c(messages, paste("📊 UMI counts - Median:", round(count_stats["Median"]), 
                                 "Mean:", round(count_stats["Mean"])))
  }
  
  if ("nFeature_RNA" %in% colnames(seurat_obj@meta.data)) {
    feature_stats <- summary(seurat_obj$nFeature_RNA)
    details$nFeature_stats <- feature_stats
    messages <- c(messages, paste("📊 Gene counts - Median:", round(feature_stats["Median"]), 
                                 "Mean:", round(feature_stats["Mean"])))
  }
  
  # Reductions
  reductions <- names(seurat_obj@reductions)
  if (length(reductions) > 0) {
    messages <- c(messages, paste("✅ Reductions:", paste(reductions, collapse = ", ")))
  } else {
    messages <- c(messages, "ℹ️ INFO: No dimensional reductions computed")
  }
  
  # Create details table
  details_df <- data.frame(
    Metric = c("Total Genes", "Total Cells", "Assays", "Reductions", "Required Metadata"),
    Value = c(nrow(seurat_obj), ncol(seurat_obj), 
             paste(assays, collapse = ", "),
             ifelse(length(reductions) > 0, paste(reductions, collapse = ", "), "None"),
             paste(required_meta[required_meta %in% colnames(seurat_obj@meta.data)], collapse = ", ")),
    Status = c("✅", "✅", "✅", 
              ifelse(length(reductions) > 0, "✅", "ℹ️"),
              ifelse(length(missing_meta) == 0, "✅", "⚠️"))
  )
  
  list(
    messages = messages,
    details_table = details_df,
    overall_status = status,
    details = details
  )
}

# File format detection and loading functions (simplified but robust)
detect_file_format <- function(filepath, filename) {
  ext <- tolower(tools::file_ext(filename))
  
  format_map <- list(
    "h5" = "10x_h5",  # Assume 10x format for .h5
    "h5ad" = "h5ad",
    "h5seurat" = "h5seurat", 
    "rds" = "seurat_rds",
    "csv" = "csv_tsv",
    "tsv" = "csv_tsv",
    "txt" = "csv_tsv",
    "xlsx" = "excel",
    "xls" = "excel",
    "loom" = "loom",
    "mtx" = "mtx"
  )
  
  format_map[[ext]] %||% "unknown"
}

validate_and_convert_gene_ids <- function(seurat_obj, species = "human") {
  gene_names <- rownames(seurat_obj)
  n_genes <- length(gene_names)

  # Detect gene ID format
  ensembl_human <- grepl("^ENSG\\d", gene_names)
  ensembl_mouse <- grepl("^ENSMUSG\\d", gene_names)
  numeric_only <- grepl("^\\d+$", gene_names)

  pct_ensembl <- sum(ensembl_human | ensembl_mouse) / n_genes
  pct_numeric <- sum(numeric_only) / n_genes

  if (pct_ensembl > 0.5) {
    # Attempt Ensembl to symbol conversion
    tryCatch({
      mapped <- map_ensembl_to_symbol(gene_names[ensembl_human | ensembl_mouse], species)
      valid_map <- mapped[!is.na(mapped) & mapped != ""]
      if (length(valid_map) > 0) {
        new_names <- gene_names
        new_names[names(valid_map)] <- valid_map
        # Handle duplicates
        dups <- duplicated(new_names)
        if (any(dups)) {
          new_names[dups] <- paste0(new_names[dups], "_dup_", seq_len(sum(dups)))
        }
        rownames(seurat_obj) <- new_names
        message("Converted ", length(valid_map), " Ensembl IDs to gene symbols")
      }
    }, error = function(e) {
      warning("Gene ID conversion failed: ", e$message, ". Keeping original IDs.")
    })
  } else if (pct_numeric > 0.5) {
    warning("Gene IDs appear to be numeric (Entrez IDs). Consider converting to gene symbols using org.Hs.eg.db::mapIds().")
  }

  seurat_obj
}

load_data_by_format <- function(filepath, format, species = "human") {
  switch(format,
    "10x_h5" = {
      if (!requireNamespace("hdf5r", quietly = TRUE)) {
        stop("Package 'hdf5r' required for .h5 files")
      }
      counts <- Read10X_h5(filepath)
      obj <- CreateSeuratObject(counts = counts, project = "scRNA", min.cells = 3, min.features = 200)
      validate_and_convert_gene_ids(obj, species)
    },
    
    "h5seurat" = {
      if (!requireNamespace("SeuratDisk", quietly = TRUE)) {
        stop("Package 'SeuratDisk' required for .h5seurat files")
      }
      SeuratDisk::LoadH5Seurat(filepath)
    },
    
    "seurat_rds" = {
      readRDS(filepath)
    },
    
    "csv_tsv" = {
      # Auto-detect separator
      first_line <- readLines(filepath, n = 1)
      sep <- if (grepl(",", first_line)) "," else "\t"
      
      data <- data.table::fread(filepath, sep = sep, header = TRUE)
      
      # Convert to matrix (genes as rows, cells as columns)
      gene_names <- data[[1]]
      expr_matrix <- as.matrix(data[, -1, drop = FALSE])
      rownames(expr_matrix) <- gene_names
      
      CreateSeuratObject(counts = expr_matrix, project = "scRNA")
    },
    
    "excel" = {
      data <- readxl::read_excel(filepath)
      gene_names <- data[[1]]
      expr_matrix <- as.matrix(data[, -1, drop = FALSE])
      rownames(expr_matrix) <- gene_names
      CreateSeuratObject(counts = expr_matrix, project = "scRNA")
    },
    
    "h5ad" = {
      # Try SeuratDisk first, then zellkonverter
      if (requireNamespace("SeuratDisk", quietly = TRUE)) {
        # Convert h5ad to h5seurat first
        h5seurat_file <- sub("\\.h5ad$", ".h5seurat", filepath)
        SeuratDisk::Convert(filepath, dest = "h5seurat", overwrite = TRUE)
        obj <- SeuratDisk::LoadH5Seurat(h5seurat_file)
        unlink(h5seurat_file)  # cleanup
        obj
      } else if (requireNamespace("zellkonverter", quietly = TRUE) &&
                 requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        sce <- zellkonverter::readH5AD(filepath)
        as.Seurat(sce, counts = "X", data = NULL)
      } else {
        stop("Package 'SeuratDisk' or 'zellkonverter' required for .h5ad files. Install with: remotes::install_github('mojaveazure/seurat-disk') or BiocManager::install('zellkonverter')")
      }
    },

    "loom" = {
      if (requireNamespace("SeuratDisk", quietly = TRUE)) {
        loom_con <- SeuratDisk::Connect(filepath, mode = "r")
        obj <- as.Seurat(loom_con)
        loom_con$close_all()
        obj
      } else if (requireNamespace("loomR", quietly = TRUE)) {
        loom_con <- loomR::connect(filepath, mode = "r")
        expr_matrix <- loom_con[["matrix"]][,]
        gene_names <- loom_con[["row_attrs/Gene"]][]
        cell_names <- loom_con[["col_attrs/CellID"]][]
        loom_con$close_all()
        rownames(expr_matrix) <- gene_names
        colnames(expr_matrix) <- cell_names
        CreateSeuratObject(counts = expr_matrix, project = "scRNA")
      } else {
        stop("Package 'SeuratDisk' or 'loomR' required for .loom files")
      }
    },

    {
      stop("Unsupported file format: ", format)
    }
  )
}

combine_seurat_objects <- function(object_list, method, name, add_batch_label = TRUE,
                                   batch_col_name = "batch", include_meta_cols = NULL,
                                   progress_callback = NULL) {
  if (length(object_list) < 2) stop("Need at least 2 objects to combine")

  dataset_names <- names(object_list)
  if (is.null(dataset_names)) {
    dataset_names <- paste0("dataset_", seq_along(object_list))
  }

  # Helper to report progress
  report_progress <- function(msg, value = NULL) {
    if (is.function(progress_callback)) {
      progress_callback(msg, value)
    }
  }

  # --- Step 1: Add batch / orig.ident metadata to each object ---
  report_progress("Preparing datasets for merge...", 5)
  for (i in seq_along(object_list)) {
    obj <- object_list[[i]]
    obj_name <- dataset_names[i]

    # Always set orig.ident to the dataset name so source is traceable
    obj$orig.ident <- obj_name

    # Add a batch label column if requested
    if (add_batch_label && nzchar(batch_col_name)) {
      obj@meta.data[[batch_col_name]] <- obj_name
    }

    object_list[[i]] <- obj
  }

  # --- Step 2: Ensure all objects share the same gene universe (union) ---
  report_progress("Computing union of gene sets across datasets...", 15)
  all_genes <- unique(unlist(lapply(object_list, rownames)))
  per_obj_genes <- lapply(object_list, rownames)
  genes_differ <- !all(sapply(per_obj_genes, function(g) setequal(g, all_genes)))

  if (genes_differ) {
    # Issue 1.3: Report per-dataset unique genes filled with zeros
    for (i in seq_along(object_list)) {
      ds_genes <- rownames(object_list[[i]])
      unique_to_ds <- setdiff(ds_genes, Reduce(union, per_obj_genes[-i]))
      if (length(unique_to_ds) > 0) {
        warning(length(unique_to_ds), " genes unique to '", dataset_names[i],
                "' will be filled with 0 in other datasets (e.g., ",
                paste(head(unique_to_ds, 3), collapse = ", "), ")")
      }
    }
    report_progress(paste0("Gene sets differ across datasets. Union contains ",
                           length(all_genes), " genes. Padding missing genes with zeros..."), 20)

    for (i in seq_along(object_list)) {
      obj <- object_list[[i]]
      missing_genes <- setdiff(all_genes, rownames(obj))

      if (length(missing_genes) > 0) {
        # Create a sparse zero matrix for the missing genes
        default_assay <- DefaultAssay(obj)
        counts_mat <- safe_get_assay_data(obj, "counts")

        zero_mat <- Matrix::sparseMatrix(
          i = integer(0), j = integer(0),
          x = numeric(0),
          dims = c(length(missing_genes), ncol(obj)),
          dimnames = list(missing_genes, colnames(obj))
        )

        # Combine existing counts with zero-padded missing genes
        combined_counts <- rbind(counts_mat, zero_mat)
        # Reorder to a consistent gene order
        combined_counts <- combined_counts[all_genes, , drop = FALSE]

        # Rebuild the Seurat object preserving metadata
        meta_backup <- obj@meta.data
        obj <- CreateSeuratObject(
          counts = combined_counts,
          project = obj@project.name,
          meta.data = meta_backup
        )
        object_list[[i]] <- obj
      }
    }
  }

  # --- Step 3: Optionally filter metadata columns to carry forward ---
  if (!is.null(include_meta_cols) && length(include_meta_cols) > 0) {
    report_progress("Filtering metadata columns...", 30)
    # Always keep essential columns
    essential_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA")
    if (add_batch_label && nzchar(batch_col_name)) {
      essential_cols <- c(essential_cols, batch_col_name)
    }
    keep_cols <- unique(c(essential_cols, include_meta_cols))

    for (i in seq_along(object_list)) {
      obj <- object_list[[i]]
      available <- intersect(keep_cols, colnames(obj@meta.data))
      obj@meta.data <- obj@meta.data[, available, drop = FALSE]
      object_list[[i]] <- obj
    }
  }

  # --- Step 4: Perform the actual merge / integration ---
  if (method == "merge") {
    report_progress("Merging datasets (simple append)...", 40)

    # Use Seurat merge with add.cell.ids to prevent barcode collisions
    merged <- tryCatch({
      merge(
        x = object_list[[1]],
        y = object_list[-1],
        add.cell.ids = dataset_names,
        project = name
      )
    }, error = function(e) {
      stop("Seurat merge failed: ", e$message,
           "\nTip: ensure all objects are valid Seurat objects with an RNA assay.")
    })

    # Issue 1.4: Check for duplicate barcodes after merge
    barcodes <- colnames(merged)
    n_dup <- sum(duplicated(barcodes))
    if (n_dup > 0) {
      warning(n_dup, " duplicate barcodes detected after merge. ",
              "Cell IDs have been prefixed with dataset names to ensure uniqueness.")
    }

    report_progress("Merge complete.", 100)
    return(merged)

  } else if (method == "harmony") {
    # Harmony integration: merge first, then run Harmony on PCA
    if (!requireNamespace("harmony", quietly = TRUE)) {
      stop("Package 'harmony' required. Install with: install.packages('harmony')")
    }
    report_progress("Merging datasets before Harmony integration...", 30)
    merged <- merge(x = object_list[[1]], y = object_list[-1],
                    add.cell.ids = dataset_names, project = name)
    report_progress("Normalizing merged data...", 50)
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- FindVariableFeatures(merged, verbose = FALSE)
    merged <- ScaleData(merged, verbose = FALSE)
    merged <- RunPCA(merged, verbose = FALSE)
    report_progress("Running Harmony integration...", 70)
    merged <- harmony::RunHarmony(merged, group.by.vars = batch_col, verbose = FALSE)
    report_progress("Harmony integration complete.", 100)
    return(merged)

  } else if (method %in% c("cca", "rpca")) {
    # CCA / RPCA integration via Seurat
    report_progress("Normalizing each dataset...", 30)
    for (i in seq_along(object_list)) {
      object_list[[i]] <- NormalizeData(object_list[[i]], verbose = FALSE)
      object_list[[i]] <- FindVariableFeatures(object_list[[i]], verbose = FALSE)
    }
    report_progress("Selecting integration features...", 50)
    features <- SelectIntegrationFeatures(object_list, nfeatures = 3000)
    if (method == "rpca") {
      for (i in seq_along(object_list)) {
        object_list[[i]] <- ScaleData(object_list[[i]], features = features, verbose = FALSE)
        object_list[[i]] <- RunPCA(object_list[[i]], features = features, verbose = FALSE)
      }
    }
    report_progress("Finding integration anchors...", 65)
    anchors <- FindIntegrationAnchors(
      object.list = object_list,
      anchor.features = features,
      reduction = if (method == "rpca") "rpca" else "cca",
      verbose = FALSE
    )
    report_progress("Integrating data...", 85)
    integrated <- IntegrateData(anchorset = anchors, verbose = FALSE)
    DefaultAssay(integrated) <- "integrated"
    report_progress(paste0(toupper(method), " integration complete."), 100)
    return(integrated)

  } else if (method == "fastmnn") {
    # FastMNN integration via batchelor
    if (!requireNamespace("batchelor", quietly = TRUE)) {
      stop("Package 'batchelor' required. Install with: BiocManager::install('batchelor')")
    }
    report_progress("Merging datasets before FastMNN...", 30)
    merged <- merge(x = object_list[[1]], y = object_list[-1],
                    add.cell.ids = dataset_names, project = name)
    merged <- NormalizeData(merged, verbose = FALSE)
    merged <- FindVariableFeatures(merged, verbose = FALSE)
    report_progress("Running FastMNN...", 60)
    # Convert to SCE, run FastMNN, store corrected embeddings
    sce_list <- lapply(object_list, function(obj) {
      obj <- NormalizeData(obj, verbose = FALSE)
      obj <- FindVariableFeatures(obj, verbose = FALSE)
      as.SingleCellExperiment(obj)
    })
    mnn_out <- do.call(batchelor::fastMNN, c(sce_list, list(d = 50)))
    # Store corrected embeddings in merged object with explicit cell alignment
    mnn_embed <- SingleCellExperiment::reducedDim(mnn_out, "corrected")
    mnn_cells <- colnames(mnn_out)
    if (!is.null(mnn_cells) && !identical(mnn_cells, colnames(merged))) {
      # Align embedding rows to merged object cell order
      mnn_embed <- mnn_embed[match(colnames(merged), mnn_cells), , drop = FALSE]
    }
    rownames(mnn_embed) <- colnames(merged)
    merged[["mnn"]] <- CreateDimReducObject(embeddings = mnn_embed, key = "MNN_", assay = "RNA")
    report_progress("FastMNN integration complete.", 100)
    return(merged)

  } else {
    stop("Unknown combination method: '", method, "'. Supported: 'merge', 'harmony', 'cca', 'rpca', 'fastmnn'.")
  }
}