# Advanced Visualization Module (enhanced with comprehensive plotting options)
# modules/visualization.R

suppressPackageStartupMessages({
  library(Seurat)
  library(plotly)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(htmlwidgets)

})

# safe_get_assay_data() and prepare_for_de() are defined in core_robustness.R

# get_color_palette(), best_reduction(), safe_default_assay(),
# detect_cell_type_col() are defined in core_robustness.R

parse_features <- function(obj, features_str) {
  if (is.null(features_str) || !nzchar(features_str)) return(character(0))
  ft <- unique(trimws(unlist(strsplit(features_str, ","))))
  ft[ft %in% rownames(obj)]
}

# -----------------------------------------------------------------------------
# AI PLOT BUILDER (standalone tab)
# -----------------------------------------------------------------------------

aiPlotBuilderUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          style = "background: linear-gradient(135deg, #f8f6ff 0%, #f0f4ff 100%); border:1px solid #d1c4e9; padding:18px; border-radius: 12px;",
          h4(icon("robot"), " AI Plot Builder", style = "margin-top:0; color: #5e35b1; font-weight: 700;"),
          p(style = "color: #666; font-size: 13px;",
            "Describe any plot in natural language and AI will generate the R code and render it."),
          uiOutput(ns("ai_plot_context_info")),
          fluidRow(
            column(9,
              textInput(ns("nl_plot_request"), NULL,
                        placeholder = "e.g., Violin plot of Tgfb1, Tgfb2, Tgfb3 in SMC cells split by batch",
                        width = "100%")
            ),
            column(3,
              actionButton(ns("ai_generate_plot"), "Generate Plot",
                           icon = icon("magic"), class = "ai-btn btn-sm", width = "100%",
                           style = "margin-top: 0px;")
            )
          ),
          uiOutput(ns("ai_plot_output")),
          uiOutput(ns("ai_plot_code_display"))
        )
      )
    )
  )
}

aiPlotBuilderServer <- function(id, values) {
  # The AI Plot Builder logic is already inside interactivePlotsServer
 # This thin wrapper calls the same module namespace so the existing
 # observeEvent handlers for ai_generate_plot etc. work.
  # (Handled by sharing the same module ID — see app.R server.)
}

# -----------------------------------------------------------------------------
# ENHANCED INTERACTIVE PLOTS UI & SERVER
# -----------------------------------------------------------------------------

interactivePlotsUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Advanced Interactive Visualization Suite"),
          p("Create publication-ready interactive plots with comprehensive customization options."),
          if (exists("build_method_card_ui", mode = "function")) build_method_card_ui("visualization"),

          # Plot type and basic settings
          fluidRow(
            column(3,
              selectInput(ns("plot_type"), "Plot Type",
                choices = c(
                  "Dimensional Reduction" = "dimplot",
                  "Feature Expression" = "feature",
                  "Violin Plot" = "violin",
                  "Dot Plot" = "dotplot",
                  "Heatmap" = "heatmap",
                  "Ridge Plot" = "ridge",
                  "Scatter Plot" = "scatter",
                  "Box Plot" = "boxplot",
                  "Density Plot" = "density",
                  "Volcano Plot" = "volcano",
                  "QC Plot" = "qc",
                  "Correlation Plot" = "correlation"
                ),
                selected = "dimplot")
            ),
            column(3,
              selectInput(ns("reduction"), "Reduction", choices = NULL)
            ),
            column(3,
              selectInput(ns("group_by"), "Group By", choices = NULL)
            ),
            column(3,
              selectInput(ns("split_by"), "Split By", choices = c("None" = "none"))
            )
          ),
          
          # Feature selection (conditional)
          conditionalPanel(
            condition = "input.plot_type == 'feature' || input.plot_type == 'violin' || input.plot_type == 'dotplot' || input.plot_type == 'heatmap' || input.plot_type == 'ridge' || input.plot_type == 'boxplot' || input.plot_type == 'correlation'",
            ns = ns,
            fluidRow(
              column(12,
                textAreaInput(ns("features"), "Features (comma-separated)",
                  placeholder = "CD3D, CD4, MS4A1, PTPRC",
                  rows = 2, width = "100%")
              )
            )
          ),
          
          # Scatter plot specific options
          conditionalPanel(
            condition = "input.plot_type == 'scatter'",
            ns = ns,
            fluidRow(
              column(3, textInput(ns("x_var"), "X variable (gene or meta)", "")),
              column(3, textInput(ns("y_var"), "Y variable (gene or meta)", "")),
              column(3, textInput(ns("color_var"), "Color by (optional)", "")),
              column(3, textInput(ns("size_var"), "Size by (optional)", ""))
            ),
            fluidRow(
              column(3, checkboxInput(ns("log_x"), "Log X", FALSE)),
              column(3, checkboxInput(ns("log_y"), "Log Y", FALSE)),
              column(3, checkboxInput(ns("add_correlation"), "Show correlation", TRUE)),
              column(3, checkboxInput(ns("add_regression"), "Add regression line", FALSE))
            )
          ),
          
          # QC plot options
          conditionalPanel(
            condition = "input.plot_type == 'qc'",
            ns = ns,
            fluidRow(
              column(4,
                selectInput(ns("qc_type"), "QC Plot Type",
                  choices = c(
                    "Metrics Overview" = "overview",
                    "UMI vs Genes" = "umi_genes",
                    "Mitochondrial %" = "mt_percent",
                    "Cell Cycle" = "cell_cycle",
                    "Doublet Scores" = "doublets"
                  ),
                  selected = "overview")
              ),
              column(4,
                selectInput(ns("qc_color_by"), "Color by", choices = NULL)
              ),
              column(4,
                checkboxInput(ns("show_outliers"), "Highlight outliers", TRUE)
              )
            )
          ),
          
          # Volcano plot options
          conditionalPanel(
            condition = "input.plot_type == 'volcano'",
            ns = ns,
            fluidRow(
              column(4,
                numericInput(ns("fc_threshold"), "Log2FC threshold", value = 0.5, min = 0, max = 3, step = 0.1)
              ),
              column(4,
                numericInput(ns("pval_threshold"), "P-value threshold", value = 0.05, min = 0.001, max = 0.1, step = 0.01)
              ),
              column(4,
                textInput(ns("volcano_genes"), "Genes to label (optional)", placeholder = "CD3D, CD4")
              )
            )
          )
        )
      )
    ),
    
    # Advanced customization panel
    fluidRow(
      column(12,
        wellPanel(
          h5("Plot Customization"),
          
          fluidRow(
            column(3,
              selectInput(ns("color_palette"), "Color Palette",
                choices = c(
                  "Default" = "default",
                  "Viridis" = "viridis",
                  "Plasma" = "plasma",
                  "RdYlBu" = "RdYlBu",
                  "Set3" = "Set3",
                  "Dark2" = "Dark2",
                  "Spectral" = "Spectral"
                ),
                selected = "default")
            ),
            column(3,
              sliderInput(ns("point_size"), "Point Size", min = 0.1, max = 5, value = 1, step = 0.1)
            ),
            column(3,
              sliderInput(ns("point_alpha"), "Point Opacity", min = 0.1, max = 1, value = 0.8, step = 0.1)
            ),
            column(3,
              selectInput(ns("theme_choice"), "Theme",
                choices = c("Minimal" = "minimal", "Classic" = "classic", "Dark" = "dark"),
                selected = "minimal")
            )
          ),
          
          fluidRow(
            column(3,
              checkboxInput(ns("show_legend"), "Show Legend", TRUE)
            ),
            column(3,
              checkboxInput(ns("show_axes"), "Show Axes", TRUE)
            ),
            column(3,
              textInput(ns("plot_title"), "Plot Title", "")
            ),
            column(3,
              div(style = "margin-top: 25px;",
                actionButton(ns("generate_plot"), "Generate Plot", 
                           class = "btn-primary btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    
    # Plot output area
    fluidRow(
      column(12,
        conditionalPanel(
          condition = sprintf("output['%s']", ns("show_plot_output")),
          ns = ns,
          wellPanel(
            h4("Visualization Output"),
            
            tabsetPanel(
              tabPanel("Interactive Plot",  icon = icon("chart-line"),
                fluidRow(
                  column(10,
                    uiOutput(ns("plot_output_ui"))
                  ),
                  column(2,
                    h5("Export Options"),
                    br(),
                    downloadButton(ns("download_interactive"), "HTML", class = "btn-success", width = "100%"),
                    br(), br(),
                    downloadButton(ns("download_static"), "PNG", class = "btn-info", width = "100%"),
                    br(), br(),
                    downloadButton(ns("download_pdf"), "PDF", class = "btn-warning", width = "100%"),
                    br(), br(),
                    h6("Dimensions"),
                    numericInput(ns("plot_width"), "Width (in)", value = 8, min = 3, max = 20),
                    numericInput(ns("plot_height"), "Height (in)", value = 6, min = 3, max = 20),
                    numericInput(ns("plot_dpi"), "DPI", value = 300, min = 72, max = 600)
                  )
                )
              ),
              
              tabPanel("Plot Statistics",   icon = icon("table"),
                fluidRow(
                  column(6,
                    h5("Data Summary"),
                    downloadButton(ns("dl_plot_stats_table"), "CSV", class = "btn-xs btn-default"),
                    DT::dataTableOutput(ns("plot_stats_table"))
                  ),
                  column(6,
                    h5("Plot Information"),
                    verbatimTextOutput(ns("plot_info"))
                  )
                )
              ),
              
              tabPanel("Advanced Options",  icon = icon("sliders-h"),
                fluidRow(
                  column(4,
                    h5("Data Filtering"),
                    checkboxInput(ns("filter_data"), "Apply data filters", FALSE),
                    conditionalPanel(
                      condition = "input.filter_data",
                      ns = ns,
                      numericInput(ns("min_cells"), "Min cells per feature", value = 1, min = 1),
                      numericInput(ns("max_cells"), "Max cells to plot", value = 10000, min = 100, max = 50000),
                      checkboxInput(ns("remove_outliers"), "Remove outliers", FALSE)
                    )
                  ),
                  column(4,
                    h5("Style Options"),
                    checkboxInput(ns("use_raster"), "Use rasterization", FALSE),
                    numericInput(ns("raster_dpi"), "Raster DPI", value = 300, min = 72, max = 600),
                    checkboxInput(ns("add_density_contours"), "Add density contours", FALSE),
                    selectInput(ns("contour_color"), "Contour color", 
                              choices = c("black", "white", "red", "blue"), selected = "black")
                  ),
                  column(4,
                    h5("Layout Options"),
                    numericInput(ns("plot_margin"), "Plot margin", value = 1, min = 0, max = 5, step = 0.1),
                    numericInput(ns("legend_size"), "Legend text size", value = 12, min = 8, max = 20),
                    numericInput(ns("axis_text_size"), "Axis text size", value = 10, min = 6, max = 16),
                    checkboxInput(ns("equal_axes"), "Equal axis scaling", TRUE)
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

interactivePlotsServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Debounced slider reactives to prevent flicker during rapid dragging
    point_size_d  <- debounce(reactive(input$point_size), 300)
    point_alpha_d <- debounce(reactive(input$point_alpha), 300)

    show_plot_output <- reactiveVal(FALSE)
    current_plot <- reactiveVal(NULL)
    plot_data_summary <- reactiveVal(NULL)
    
    output$show_plot_output <- reactive({ show_plot_output() })
    outputOptions(output, "show_plot_output", suspendWhenHidden = FALSE)

    # ---- AI Plot Builder: context info showing available columns ----
    output$ai_plot_context_info <- renderUI({
      obj <- values$active_dataset
      if (is.null(obj)) {
        return(div(style = "font-size: 12px; color: #999; margin-bottom: 8px;",
          icon("info-circle"), " Load a dataset to see available columns and genes"))
      }
      meta <- obj@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      num_cols <- colnames(meta)[sapply(meta, is.numeric)]
      reductions <- names(obj@reductions)

      # Build compact chips for columns
      make_chips <- function(cols, color, max_show = 10) {
        shown <- head(cols, max_show)
        chips <- paste(sapply(shown, function(c) {
          paste0("<span style='display:inline-block;background:", color,
                 ";color:#333;padding:1px 7px;border-radius:8px;font-size:11px;margin:1px 2px;'>", c, "</span>")
        }), collapse = "")
        if (length(cols) > max_show) {
          chips <- paste0(chips, "<span style='font-size:11px;color:#999;'> +", length(cols) - max_show, " more</span>")
        }
        chips
      }

      div(style = "font-size: 12px; margin-bottom: 8px; padding: 8px 10px; background: rgba(255,255,255,0.7); border-radius: 8px; border: 1px solid #e8e0f0;",
        tags$div(style = "margin-bottom: 3px;",
          tags$strong(style = "color: #5e35b1;", "Group/Color by: "),
          HTML(make_chips(cat_cols, "#ede7f6", 8))
        ),
        tags$div(style = "margin-bottom: 3px;",
          tags$strong(style = "color: #5e35b1;", "Numeric: "),
          HTML(make_chips(num_cols, "#e8eaf6", 6))
        ),
        if (length(reductions) > 0) tags$div(
          tags$strong(style = "color: #5e35b1;", "Reductions: "),
          HTML(make_chips(reductions, "#e0f2f1", 5))
        ) else NULL,
        tags$div(style = "color: #888; margin-top: 2px;",
          paste0(nrow(obj), " genes | ", ncol(obj), " cells")
        )
      )
    })

    # ---- AI Natural-Language Plot Generator ----
    ai_plot_obj <- reactiveVal(NULL)

    observeEvent(input$ai_generate_plot, {
      if (is.null(values$ai_api_key) || !nzchar(values$ai_api_key)) {
        showNotification("Please set your Claude API key in Data Management > AI Config first.",
                         type = "warning", duration = 6)
        return()
      }
      if (is.null(values$active_dataset)) {
        showNotification("No dataset loaded.", type = "warning", duration = 4)
        return()
      }
      if (!nzchar(input$nl_plot_request)) {
        showNotification("Please describe the plot you want.", type = "warning", duration = 4)
        return()
      }

      output$ai_plot_output <- renderUI({ ai_thread_loading("AI generating your plot...") })
      output$ai_plot_code_display <- renderUI(NULL)

      obj <- values$active_dataset
      available_reductions <- paste(names(obj@reductions), collapse = ", ")
      meta_cols <- paste(colnames(obj@meta.data), collapse = ", ")

      # Show some example values per categorical column for context
      meta <- obj@meta.data
      cat_cols <- colnames(meta)[sapply(meta, function(x) is.factor(x) || is.character(x))]
      col_vals <- paste(sapply(head(cat_cols, 15), function(col) {
        vals <- unique(as.character(meta[[col]]))
        paste0(col, ": [", paste(head(vals, 6), collapse = ", "),
               if (length(vals) > 6) paste0("... (", length(vals), ")") else "", "]")
      }), collapse = "\n")

      gene_names <- rownames(obj)
      sample_genes <- paste(head(gene_names, 100), collapse = ", ")

      prompt <- paste0(
        "The user wants this plot: '", input$nl_plot_request, "'.\n\n",
        "DATASET INFO:\n",
        "- ", ncol(obj), " cells, ", nrow(obj), " genes\n",
        "- Available reductions: ", available_reductions, "\n",
        "- Metadata columns and example values:\n", col_vals, "\n",
        "- Sample gene names: ", sample_genes, "\n\n",
        "IMPORTANT RULES:\n",
        "1. The Seurat object is available as `obj` — do NOT create or load data\n",
        "2. Do NOT use library() calls — all packages are already loaded\n",
        "3. The LAST expression must produce the plot (e.g., the plot call itself)\n",
        "4. Use Seurat functions: DimPlot, FeaturePlot, VlnPlot, DotPlot, DoHeatmap, RidgePlot\n",
        "5. For subsetting cells, use: subset(obj, subset = column_name == 'value')\n",
        "6. If a gene or column does not exist, the code will fail — only use names from the lists above\n",
        "7. Return ONLY the R code inside ```r blocks, nothing else\n"
      )

      withProgress(message = "AI generating plot...", value = 0.3, {
        result <- call_ai("nlplot_codegen", prompt, values)
      })

      if (!result$ok || is.null(result$text)) {
        output$ai_plot_output <- renderUI({ ai_thread_error(result$error %||% "Unknown error") })
        return()
      }

      # Extract code from markdown code block
      code_text <- result$text
      code_match <- regmatches(code_text,
        regexpr("```[rR]?\\s*\\n([\\s\\S]*?)```", code_text, perl = TRUE))
      if (length(code_match) > 0 && nzchar(code_match)) {
        code_text <- sub("^```[rR]?\\s*\\n", "", code_match)
        code_text <- sub("\\s*```$", "", code_text)
      }
      # Remove any library() calls
      code_text <- gsub("library\\([^)]*\\)\\s*\n?", "", code_text)
      code_text <- trimws(code_text)

      # Execute the code with the user's dataset
      plot_result <- tryCatch({
        env <- new.env(parent = globalenv())
        env$obj <- obj
        eval(parse(text = code_text), envir = env)
      }, error = function(e) {
        list(error = e$message)
      })

      if (is.list(plot_result) && !is.null(plot_result$error)) {
        # Show error + code for debugging
        output$ai_plot_output <- renderUI({
          ai_thread_error(paste("Plot execution failed:", plot_result$error))
        })
        output$ai_plot_code_display <- renderUI({
          div(style = "margin-top: 8px;",
            tags$details(
              tags$summary(style = "cursor: pointer; color: #5e35b1; font-weight: 600; font-size: 13px;",
                icon("code"), " Show generated code"),
              tags$pre(style = "background:#1e1e2e; color:#cdd6f4; padding:12px; border-radius:10px; margin-top:6px; border-left: 3px solid #7c4dff;",
                tags$code(code_text))
            )
          )
        })
      } else {
        # Render the plot
        ai_plot_obj(plot_result)

        output$ai_plot_render <- renderPlot({
          p <- ai_plot_obj()
          req(p)
          if (inherits(p, "gg") || inherits(p, "patchwork")) {
            print(p)
          } else if (inherits(p, "Heatmap") || inherits(p, "HeatmapList")) {
            ComplexHeatmap::draw(p)
          } else {
            print(p)
          }
        }, height = 500)

        output$ai_plot_output <- renderUI({
          div(class = "ai-thread-container",
            div(class = "ai-thread-header-bar",
              div(class = "ai-thread-avatar", icon("chart-bar")),
              tags$span(class = "ai-thread-sender", "AI Plot Builder"),
              tags$span(class = "ai-thread-model",
                if (grepl("opus", get_model_for_task("nlplot_codegen", values), ignore.case = TRUE)) "OPUS"
                else if (grepl("haiku", get_model_for_task("nlplot_codegen", values), ignore.case = TRUE)) "HAIKU"
                else "SONNET"),
              tags$span(class = "ai-thread-time", format(Sys.time(), "%H:%M"))
            ),
            div(class = "ai-thread-message",
              div(class = "ai-thread-body",
                div(style = "background: white; border-radius: 10px; padding: 8px; margin-bottom: 8px;",
                  plotOutput(session$ns("ai_plot_render"), height = "500px")
                )
              ),
              div(class = "ai-thread-footer",
                tags$span(class = "ai-thread-tokens",
                  if (!is.na(result$tokens_in) && !is.na(result$tokens_out))
                    tagList(icon("microchip"), paste0(result$tokens_in, " in / ", result$tokens_out, " out"))
                  else ""
                ),
                downloadButton(session$ns("dl_ai_plot"), "Download Plot",
                  class = "ai-thread-apply-btn", style = "padding: 6px 16px;")
              )
            )
          )
        })

        # Show code in collapsible section
        output$ai_plot_code_display <- renderUI({
          div(style = "margin-top: 8px;",
            tags$details(
              tags$summary(style = "cursor: pointer; color: #5e35b1; font-weight: 600; font-size: 13px;",
                icon("code"), " Show generated code"),
              tags$pre(style = "background:#1e1e2e; color:#cdd6f4; padding:12px; border-radius:10px; margin-top:6px; border-left: 3px solid #7c4dff;",
                tags$code(code_text))
            )
          )
        })

        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, plot_result, "AI Plot Builder",
            figure_type = paste("AI-generated:", input$nl_plot_request),
            label = "AI Plot", source_module = "Visualization")
        }
      }
    })

    # Download AI-generated plot
    output$dl_ai_plot <- downloadHandler(
      filename = function() paste0("ai_plot_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
      content = function(file) {
        p <- ai_plot_obj()
        if (!is.null(p)) {
          grDevices::png(file, width = 10, height = 7, units = "in", res = 300)
          tryCatch({
            if (inherits(p, "gg") || inherits(p, "patchwork")) print(p) else print(p)
          }, finally = grDevices::dev.off())
        }
      }
    )

    # ---- AI Helper (generic Q&A modal for Visualization) ----
    if (exists("ai_helper_server", mode = "function")) {
      ai_helper_server("viz", values, "Visualization", context_fn = function() {
        obj <- values$active_dataset
        if (is.null(obj)) return("")
        paste0(
          "Available plot types: DimPlot, FeaturePlot, VlnPlot, DotPlot, DoHeatmap, ",
          "RidgePlot, Scatter, BoxPlot, Density, Volcano, QC, Correlation.\n",
          "Available reductions: ", paste(names(obj@reductions), collapse = ", "), ".\n",
          "Metadata columns: ", paste(colnames(obj@meta.data), collapse = ", "), "."
        )
      })
    }

    # Update UI choices based on active dataset
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      
      # Update reduction choices
      reds <- names(obj@reductions)
      if (length(reds) > 0) {
        best_red <- best_reduction(obj)
        updateSelectInput(session, "reduction", choices = reds, selected = best_red)
      }
      
      # Update metadata choices
      meta_cols <- colnames(obj@meta.data)
      factor_like <- meta_cols[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      
      updateSelectInput(session, "group_by", choices = c("None" = "none", meta_cols))
      updateSelectInput(session, "split_by", choices = c("None" = "none", factor_like))
      updateSelectInput(session, "qc_color_by", choices = c("None" = "none", meta_cols))
    })
    
    # Main plot generation
    observeEvent(input$generate_plot, {
      req(values$active_dataset)
      
      withProgress(message = "Generating plot...", value = 0, {
        obj <- values$active_dataset
        
        incProgress(0.3, detail = "Preparing data...")
        
        # Generate plot based on type
        result <- create_plot_by_type(obj, input)
        
        if (!is.null(result$plot)) {
          current_plot(result$plot)
          plot_data_summary(result$data_summary)
          show_plot_output(TRUE)
          
          showNotification("Plot generated successfully", type = "message")
        } else {
          showNotification(paste("Plot generation failed:", result$error), type = "error")
        }
        
        incProgress(1)
      })
    })
    
    # Dynamic plot output UI based on plot type
    output$plot_output_ui <- renderUI({
      req(current_plot())
      
      if (input$plot_type == "heatmap") {
        # For complex heatmaps, use plotOutput
        plotOutput(session$ns("static_plot"), height = "700px")
      } else {
        # For interactive plots, use plotlyOutput
        plotlyOutput(session$ns("interactive_plot"), height = "700px")
      }
    })
    
    # Render interactive plot
    output$interactive_plot <- renderPlotly({
      req(current_plot())
      tryCatch({
        if (input$plot_type != "heatmap") {
          p <- current_plot() %>% plotly::partial_bundle()
          if (exists("capture_figure", mode = "function")) {
            capture_figure(values, p, "Visualization", figure_type = input$plot_type,
                           label = "Interactive Plot", source_module = "Visualization")
          }
          p
        }
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })
    
    # Render static plot (for heatmaps)
    output$static_plot <- renderPlot({
      req(current_plot())
      if (input$plot_type == "heatmap") {
        current_plot()
      }
    })
    
    # Plot statistics
    output$plot_stats_table <- DT::renderDataTable({
      req(plot_data_summary())
      plot_data_summary()
    }, options = list(dom = 't', pageLength = 10))
    
    # Plot information
    output$plot_info <- renderText({
      req(current_plot(), values$active_dataset)
      
      obj <- values$active_dataset
      info_text <- paste(
        "Plot Type:", input$plot_type,
        "\nDataset:", ncol(obj), "cells,", nrow(obj), "genes",
        "\nActive Assay:", safe_default_assay(obj),
        "\nReductions Available:", paste(names(obj@reductions), collapse = ", "),
        "\nPlot Generated:", Sys.time(),
        sep = ""
      )
      
      # Add plot-specific information
      if (input$plot_type == "feature" && nzchar(input$features)) {
        features <- parse_features(obj, input$features)
        info_text <- paste(info_text, "\nFeatures Plotted:", length(features), sep = "")
      }
      
      if (input$plot_type == "dimplot" && input$reduction != "") {
        info_text <- paste(info_text, "\nReduction Used:", input$reduction, sep = "")
      }
      
      info_text
    })
    
    # Download handlers
    output$dl_plot_stats_table <- downloadHandler(
      filename = function() paste0("plot_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        req(plot_data_summary())
        write.csv(plot_data_summary(), file, row.names = FALSE)
      }
    )

    output$download_interactive <- downloadHandler(
      filename = function() paste0(input$plot_type, "_plot_", Sys.Date(), ".html"),
      content = function(file) {
        if (!is.null(current_plot()) && input$plot_type != "heatmap") {
          htmlwidgets::saveWidget(current_plot(), file)
        }
      }
    )
    
    output$download_static <- downloadHandler(
      filename = function() paste0(input$plot_type, "_plot_", Sys.Date(), ".png"),
      content = function(file) {
        if (!is.null(current_plot())) {
          if (input$plot_type == "heatmap") {
            png(file, width = input$plot_width * input$plot_dpi, 
                height = input$plot_height * input$plot_dpi, res = input$plot_dpi)
            print(current_plot())
            dev.off()
          } else {
            # Convert plotly to static
            plotly_plot <- current_plot()
            {
  ok <- TRUE
  tryCatch({
    plotly::orca(plotly_plot, file = file, width = input$plot_width * input$plot_dpi, height = input$plot_height * input$plot_dpi)
  }, error = function(e) { ok <<- FALSE })
  if (!ok) {
    tmp_html <- tempfile(fileext = ".html")
    htmlwidgets::saveWidget(plotly_plot, tmp_html, selfcontained = TRUE)
    if (requireNamespace("webshot2", quietly = TRUE)) {
      webshot2::webshot(tmp_html, file = file, vwidth = input$plot_width * input$plot_dpi, vheight = input$plot_height * input$plot_dpi)
    } else {
      # Fallback: provide HTML when PNG export isn't available
      file.copy(tmp_html, sub("\\.png$", ".html", file), overwrite = TRUE)
    }
  }
}
}
        }
      }
    )
    
    output$download_pdf <- downloadHandler(
      filename = function() paste0(input$plot_type, "_plot_", Sys.Date(), ".pdf"),
      content = function(file) {
        if (!is.null(current_plot())) {
          if (input$plot_type == "heatmap") {
            pdf(file, width = input$plot_width, height = input$plot_height)
            print(current_plot())
            dev.off()
          } else {
            # For plotly plots, first convert to ggplot if possible, then save as PDF
            showNotification("PDF export for interactive plots requires additional setup", type = "warning")
          }
        }
      }
    )
  })
}

# ============================================================================
# PLOT CREATION FUNCTIONS
# ============================================================================

create_plot_by_type <- function(seurat_obj, params) {
  tryCatch({
    plot_result <- switch(params$plot_type,
      "dimplot" = create_enhanced_dimplot(seurat_obj, params),
      "feature" = create_enhanced_feature_plot(seurat_obj, params),
      "violin" = create_enhanced_violin_plot(seurat_obj, params),
      "dotplot" = create_enhanced_dot_plot(seurat_obj, params),
      "heatmap" = create_enhanced_heatmap(seurat_obj, params),
      "ridge" = create_enhanced_ridge_plot(seurat_obj, params),
      "scatter" = create_enhanced_scatter_plot(seurat_obj, params),
      "boxplot" = create_enhanced_box_plot(seurat_obj, params),
      "density" = create_enhanced_density_plot(seurat_obj, params),
      "volcano" = create_enhanced_volcano_plot(seurat_obj, params),
      "qc" = create_enhanced_qc_plot(seurat_obj, params),
      "correlation" = create_enhanced_correlation_plot(seurat_obj, params),
      list(plot = NULL, data_summary = NULL, error = "Unknown plot type")
    )
    
    plot_result
    
  }, error = function(e) {
    list(plot = NULL, data_summary = NULL, error = e$message)
  })
}

create_enhanced_dimplot <- function(seurat_obj, params) {
  reduction <- params$reduction
  group_by <- if (params$group_by != "none") params$group_by else NULL
  split_by <- if (params$split_by != "none") params$split_by else NULL
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    reduction <- best_reduction(seurat_obj)
  }
  
  if (is.null(reduction)) {
    return(list(plot = NULL, data_summary = NULL, error = "No reduction available"))
  }
  
  # Get embedding
  embedding <- Embeddings(seurat_obj, reduction)[, 1:2]
  
  # Prepare plot data
  plot_data <- data.frame(
    Dim1 = embedding[, 1],
    Dim2 = embedding[, 2],
    Cell = rownames(embedding)
  )

  # Progressive rendering: subsample for large datasets
  n_total <- nrow(plot_data)
  plot_idx <- adaptive_plot_index(n_total, max_points = 50000L)
  subsampled <- length(plot_idx) < n_total
  if (subsampled) {
    plot_data <- plot_data[plot_idx, , drop = FALSE]
  }

  # Add grouping variable
  if (!is.null(group_by) && group_by %in% colnames(seurat_obj@meta.data)) {
    plot_data$Group <- seurat_obj@meta.data[plot_data$Cell, group_by]
    color_var <- ~Group
    colors <- get_enhanced_colors(length(unique(plot_data$Group)), params$color_palette)
  } else {
    color_var <- NULL
    colors <- "blue"
  }
  
  # Create plotly plot
  p <- plot_ly(
    data = plot_data,
    x = ~Dim1,
    y = ~Dim2,
    color = color_var,
    colors = colors,
    type = "scatter",
    mode = "markers",
    marker = list(size = params$point_size * 2, opacity = params$point_alpha),
    text = ~Cell,
    hovertemplate = paste0("Cell: %{text}<br>", toupper(reduction), "1: %{x}<br>", toupper(reduction), "2: %{y}<extra></extra>")
  )
  
  # Apply theme and styling
  p <- apply_plot_theme(p, params, paste(toupper(reduction), "Plot"))

  if (subsampled) {
    p <- p %>% layout(
      annotations = list(list(
        text = paste0("Showing ", format(length(plot_idx), big.mark = ","),
                      " of ", format(n_total, big.mark = ","), " cells"),
        xref = "paper", yref = "paper", x = 1, y = 0,
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 10, color = "#888888")
      ))
    )
  }

  # Data summary
  data_summary <- create_data_summary(plot_data, "Dimensional reduction plot")

  list(plot = p, data_summary = data_summary, error = NULL)
}

create_enhanced_feature_plot <- function(seurat_obj, params) {
  features <- parse_features(seurat_obj, params$features)
  
  if (length(features) == 0) {
    return(list(plot = NULL, data_summary = NULL, error = "No valid features specified"))
  }
  
  reduction <- params$reduction
  if (!reduction %in% names(seurat_obj@reductions)) {
    reduction <- best_reduction(seurat_obj)
  }
  
  if (is.null(reduction)) {
    return(list(plot = NULL, data_summary = NULL, error = "No reduction available"))
  }
  
  # Limit features for performance
  features <- head(features, 9)  # Max 9 features for subplot layout
  
  # Get embedding and expression data
  embedding <- Embeddings(seurat_obj, reduction)[, 1:2]
  expr_data <- safe_get_assay_data(seurat_obj, "data")[features, , drop = FALSE]

  # Progressive rendering: subsample for large datasets
  n_total <- nrow(embedding)
  plot_idx <- adaptive_plot_index(n_total, max_points = 50000L)
  subsampled <- length(plot_idx) < n_total
  if (subsampled) {
    embedding <- embedding[plot_idx, , drop = FALSE]
    expr_data <- expr_data[, plot_idx, drop = FALSE]
  }

  # Create subplots for multiple features
  if (length(features) == 1) {
    # Single feature plot
    plot_data <- data.frame(
      Dim1 = embedding[, 1],
      Dim2 = embedding[, 2],
      Expression = as.numeric(expr_data[1, ]),
      Cell = colnames(expr_data)
    )

    colors <- get_enhanced_colors(100, params$color_palette, continuous = TRUE)

    p <- plot_ly(
      data = plot_data,
      x = ~Dim1,
      y = ~Dim2,
      color = ~Expression,
      colors = colors,
      type = "scatter",
      mode = "markers",
      marker = list(size = params$point_size * 2, opacity = params$point_alpha),
      text = ~Cell,
      hovertemplate = paste0("Cell: %{text}<br>", features[1], ": %{marker.color:.3f}<br>",
                            toupper(reduction), "1: %{x}<br>", toupper(reduction), "2: %{y}<extra></extra>")
    ) %>% colorbar(title = features[1])

    p <- apply_plot_theme(p, params, paste("Expression:", features[1]))

    if (subsampled) {
      p <- p %>% layout(
        annotations = list(list(
          text = paste0("Showing ", format(length(plot_idx), big.mark = ","),
                        " of ", format(n_total, big.mark = ","), " cells"),
          xref = "paper", yref = "paper", x = 1, y = 0,
          xanchor = "right", yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 10, color = "#888888")
        ))
      )
    }

  } else {
    # Multiple feature plots as subplots
    plots <- list()

    for (i in seq_along(features)) {
      feat <- features[i]
      plot_data <- data.frame(
        Dim1 = embedding[, 1],
        Dim2 = embedding[, 2],
        Expression = as.numeric(expr_data[feat, ]),
        Cell = colnames(expr_data)
      )

      colors <- get_enhanced_colors(100, params$color_palette, continuous = TRUE)

      plots[[i]] <- plot_ly(
        data = plot_data,
        x = ~Dim1,
        y = ~Dim2,
        color = ~Expression,
        colors = colors,
        type = "scatter",
        mode = "markers",
        marker = list(size = params$point_size, opacity = params$point_alpha),
        showlegend = FALSE
      ) %>% layout(
        title = feat,
        xaxis = list(title = "", showticklabels = FALSE),
        yaxis = list(title = "", showticklabels = FALSE)
      )
    }

    # Arrange subplots
    n_cols <- ceiling(sqrt(length(features)))
    n_rows <- ceiling(length(features) / n_cols)

    p <- subplot(plots, nrows = n_rows, shareX = TRUE, shareY = TRUE) %>%
      layout(title = "Feature Expression")

    if (subsampled) {
      p <- p %>% layout(
        annotations = list(list(
          text = paste0("Showing ", format(length(plot_idx), big.mark = ","),
                        " of ", format(n_total, big.mark = ","), " cells"),
          xref = "paper", yref = "paper", x = 1, y = 0,
          xanchor = "right", yanchor = "bottom",
          showarrow = FALSE,
          font = list(size = 10, color = "#888888")
        ))
      )
    }
  }

  # Data summary
  data_summary <- create_feature_data_summary(expr_data, features)

  list(plot = p, data_summary = data_summary, error = NULL)
}

create_enhanced_violin_plot <- function(seurat_obj, params) {
  features <- parse_features(seurat_obj, params$features)
  
  if (length(features) == 0) {
    return(list(plot = NULL, data_summary = NULL, error = "No valid features specified"))
  }
  
  group_by <- if (params$group_by != "none") params$group_by else "all_cells"
  
  # Prepare data
  expr_data <- safe_get_assay_data(seurat_obj, "data")[features, , drop = FALSE]

  # Create long format data
  plot_data_list <- list()
  
  for (feat in features) {
    feat_data <- data.frame(
      Expression = as.numeric(expr_data[feat, ]),
      Feature = feat,
      Cell = colnames(seurat_obj)
    )
    
    if (group_by != "all_cells" && group_by %in% colnames(seurat_obj@meta.data)) {
      feat_data$Group <- seurat_obj@meta.data[[group_by]]
    } else {
      feat_data$Group <- "All Cells"
    }
    
    plot_data_list[[feat]] <- feat_data
  }
  
  plot_data <- do.call(rbind, plot_data_list)

  # Progressive rendering: subsample to max 20,000 cells per group
  n_total <- nrow(plot_data)
  max_per_group <- 20000L
  groups <- unique(plot_data$Group)
  if (any(table(plot_data$Group) > max_per_group)) {
    keep <- unlist(lapply(groups, function(g) {
      idx <- which(plot_data$Group == g)
      if (length(idx) > max_per_group) {
        idx[adaptive_plot_index(length(idx), max_points = max_per_group)]
      } else {
        idx
      }
    }))
    plot_data <- plot_data[sort(keep), , drop = FALSE]
  }
  subsampled_violin <- nrow(plot_data) < n_total

  # Create violin plot
  colors <- get_enhanced_colors(length(unique(plot_data$Group)), params$color_palette)

  p <- plot_ly(
    data = plot_data,
    x = ~Group,
    y = ~Expression,
    color = ~Group,
    colors = colors,
    type = "violin",
    box = list(visible = TRUE),
    meanline = list(visible = TRUE)
  ) %>% layout(
    title = "Feature Expression Distribution",
    xaxis = list(title = group_by),
    yaxis = list(title = "Expression Level")
  )

  if (length(features) > 1) {
    p <- p %>% layout(facet_col = ~Feature)
  }

  p <- apply_plot_theme(p, params, "Violin Plot")

  if (subsampled_violin) {
    p <- p %>% layout(
      annotations = list(list(
        text = paste0("Showing ", format(nrow(plot_data), big.mark = ","),
                      " of ", format(n_total, big.mark = ","), " points (max ",
                      format(max_per_group, big.mark = ","), " per group)"),
        xref = "paper", yref = "paper", x = 1, y = 0,
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 10, color = "#888888")
      ))
    )
  }

  # Data summary
  data_summary <- create_violin_data_summary(plot_data)

  list(plot = p, data_summary = data_summary, error = NULL)
}

create_enhanced_dot_plot <- function(seurat_obj, params) {
  features <- parse_features(seurat_obj, params$features)
  
  if (length(features) == 0) {
    return(list(plot = NULL, data_summary = NULL, error = "No valid features specified"))
  }
  
  group_by <- if (params$group_by != "none") params$group_by else NULL
  
  if (is.null(group_by) || !group_by %in% colnames(seurat_obj@meta.data)) {
    # Create artificial grouping
    seurat_obj$temp_group <- "All Cells"
    group_by <- "temp_group"
  }
  
  # Calculate average expression and percent expressed
  groups <- unique(seurat_obj@meta.data[[group_by]])
  
  dot_data_list <- list()
  
  for (group in groups) {
    group_cells <- rownames(seurat_obj@meta.data)[seurat_obj@meta.data[[group_by]] == group]
    group_expr <- safe_get_assay_data(seurat_obj, "data")[features, group_cells, drop = FALSE]
    
    for (feat in features) {
      feat_expr <- as.numeric(group_expr[feat, ])
      
      dot_data_list[[paste(group, feat, sep = "_")]] <- data.frame(
        Group = group,
        Feature = feat,
        Avg_Expression = mean(feat_expr),
        Pct_Expressed = sum(feat_expr > 0) / length(feat_expr) * 100
      )
    }
  }
  
  dot_data <- do.call(rbind, dot_data_list)
  
  # Create dot plot
  p <- plot_ly(
    data = dot_data,
    x = ~Feature,
    y = ~Group,
    color = ~Avg_Expression,
    size = ~Pct_Expressed,
    colors = get_enhanced_colors(100, params$color_palette, continuous = TRUE),
    type = "scatter",
    mode = "markers",
    marker = list(opacity = 0.8, sizemode = 'diameter', sizeref = 2),
    hovertemplate = "Group: %{y}<br>Feature: %{x}<br>Avg Expression: %{marker.color:.3f}<br>% Expressed: %{marker.size:.1f}%<extra></extra>"
  ) %>% layout(
    title = "Dot Plot",
    xaxis = list(title = "Features"),
    yaxis = list(title = "Groups")
  ) %>% colorbar(title = "Avg Expr")
  
  p <- apply_plot_theme(p, params, "Dot Plot")
  
  # Data summary
  data_summary <- create_dot_data_summary(dot_data)
  
  list(plot = p, data_summary = data_summary, error = NULL)
}

create_enhanced_heatmap <- function(seurat_obj, params) {
  features <- parse_features(seurat_obj, params$features)
  
  if (length(features) == 0) {
    return(list(plot = NULL, data_summary = NULL, error = "No valid features specified"))
  }
  
  group_by <- if (params$group_by != "none") params$group_by else NULL

  # Subsample genes for performance if needed (max 100)
  n_genes_total <- length(features)
  if (n_genes_total > 100L) {
    features <- features[adaptive_plot_index(n_genes_total, max_points = 100L)]
  }

  # Get expression data
  expr_data <- safe_get_assay_data(seurat_obj, "data")[features, , drop = FALSE]

  # Subsample cells for performance if needed (max 5,000)
  n_cells_total <- ncol(expr_data)
  cell_idx <- adaptive_plot_index(n_cells_total, max_points = 5000L)
  subsampled_hm <- length(cell_idx) < n_cells_total || n_genes_total > 100L
  if (length(cell_idx) < n_cells_total) {
    expr_data <- expr_data[, cell_idx, drop = FALSE]
  }

  # Convert to matrix for heatmap
  expr_matrix <- as.matrix(expr_data)
  
  # Create annotation if grouping specified
  if (!is.null(group_by) && group_by %in% colnames(seurat_obj@meta.data)) {
    col_annotation <- seurat_obj@meta.data[colnames(expr_matrix), group_by, drop = FALSE]
    names(col_annotation) <- group_by
  } else {
    col_annotation <- NULL
  }
  
  # Create heatmap using ComplexHeatmap
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    # Color function
    if (params$color_palette == "viridis") {
      col_fun <- circlize::colorRamp2(c(0, max(expr_matrix)), viridis::viridis(2))
    } else {
      col_fun <- circlize::colorRamp2(c(0, max(expr_matrix)), c("white", "red"))
    }
    
    # Build titles with subsample info
    col_title <- if (length(cell_idx) < n_cells_total) {
      paste0("Cells (showing ", format(length(cell_idx), big.mark = ","),
             " of ", format(n_cells_total, big.mark = ","), ")")
    } else {
      "Cells"
    }
    row_title <- if (n_genes_total > 100L) {
      paste0("Features (showing ", length(features), " of ", n_genes_total, ")")
    } else {
      "Features"
    }

    # Create heatmap
    ht <- ComplexHeatmap::Heatmap(
      expr_matrix,
      name = "Expression",
      col = col_fun,
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_row_names = length(features) <= 50,
      show_column_names = FALSE,
      column_title = col_title,
      row_title = row_title,
      heatmap_legend_param = list(title = "Expression")
    )
    
    if (!is.null(col_annotation)) {
      ht <- ht + ComplexHeatmap::HeatmapAnnotation(
        df = col_annotation,
        which = "column"
      )
    }
    
    p <- ht
    
  } else {
    # Fallback to simple heatmap
    p <- heatmap(expr_matrix, 
                 main = "Expression Heatmap",
                 xlab = "Cells", ylab = "Features")
  }
  
  # Data summary
  data_summary <- create_heatmap_data_summary(expr_matrix)
  
  list(plot = p, data_summary = data_summary, error = NULL)
}

create_enhanced_scatter_plot <- function(seurat_obj, params) {
  x_var <- params$x_var
  y_var <- params$y_var
  color_var <- params$color_var
  size_var <- params$size_var
  
  if (!nzchar(x_var) || !nzchar(y_var)) {
    return(list(plot = NULL, data_summary = NULL, error = "X and Y variables required"))
  }
  
  # Helper function to get variable data
  get_var_data <- function(var_name) {
    if (var_name %in% rownames(seurat_obj)) {
      as.numeric(safe_get_assay_data(seurat_obj, "data")[var_name, ])
    } else if (var_name %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[var_name]]
    } else {
      NULL
    }
  }
  
  # Get data for variables
  x_data <- get_var_data(x_var)
  y_data <- get_var_data(y_var)
  
  if (is.null(x_data) || is.null(y_data)) {
    return(list(plot = NULL, data_summary = NULL, error = "Variables not found in dataset"))
  }
  
  # Create plot data
  plot_data <- data.frame(
    X = x_data,
    Y = y_data,
    Cell = colnames(seurat_obj)
  )
  
  # Add color variable if specified
  if (nzchar(color_var)) {
    color_data <- get_var_data(color_var)
    if (!is.null(color_data)) {
      plot_data$Color <- color_data
    }
  }
  
  # Add size variable if specified
  if (nzchar(size_var)) {
    size_data <- get_var_data(size_var)
    if (!is.null(size_data) && is.numeric(size_data)) {
      plot_data$Size <- size_data
    }
  }

  # Progressive rendering: subsample for large datasets
  n_total_scatter <- nrow(plot_data)
  scatter_idx <- adaptive_plot_index(n_total_scatter, max_points = 50000L)
  subsampled_scatter <- length(scatter_idx) < n_total_scatter
  if (subsampled_scatter) {
    plot_data <- plot_data[scatter_idx, , drop = FALSE]
  }

  # Build marker list
  marker_opts <- list(opacity = params$point_alpha)
  if (!"Size" %in% colnames(plot_data)) {
    marker_opts$size <- params$point_size * 4
  }

  # Create scatter plot with all options upfront (avoids add_trace type/mode warnings)
  plot_args <- list(
    data = plot_data,
    x = ~X,
    y = ~Y,
    type = "scatter",
    mode = "markers",
    marker = marker_opts,
    text = ~Cell
  )
  if ("Color" %in% colnames(plot_data)) {
    plot_args$color <- ~Color
  }
  if ("Size" %in% colnames(plot_data)) {
    plot_args$size <- ~Size
  }

  p <- do.call(plot_ly, plot_args)
  
  # Apply log transformations if requested
  layout_list <- list(
    title = paste(x_var, "vs", y_var),
    xaxis = list(title = x_var),
    yaxis = list(title = y_var)
  )
  
  if (params$log_x && all(plot_data$X > 0)) {
    layout_list$xaxis$type <- "log"
  }
  
  if (params$log_y && all(plot_data$Y > 0)) {
    layout_list$yaxis$type <- "log"
  }
  
  p <- p %>% layout(layout_list)
  
  # Add correlation if requested
  if (params$add_correlation && is.numeric(plot_data$X) && is.numeric(plot_data$Y)) {
    correlation <- cor(plot_data$X, plot_data$Y, use = "complete.obs")
    p <- p %>% layout(
      title = paste(x_var, "vs", y_var, "(r =", round(correlation, 3), ")")
    )
  }
  
  # Add regression line if requested
  if (params$add_regression && is.numeric(plot_data$X) && is.numeric(plot_data$Y)) {
    lm_model <- lm(Y ~ X, data = plot_data)
    pred_data <- data.frame(X = seq(min(plot_data$X), max(plot_data$X), length.out = 100))
    pred_data$Y <- predict(lm_model, pred_data)
    
    p <- p %>% add_trace(
      data = pred_data,
      x = ~X,
      y = ~Y,
      type = "scatter",
      mode = "lines",
      line = list(color = "red"),
      showlegend = FALSE
    )
  }
  
  p <- apply_plot_theme(p, params, "Scatter Plot")

  if (subsampled_scatter) {
    p <- p %>% layout(
      annotations = list(list(
        text = paste0("Showing ", format(length(scatter_idx), big.mark = ","),
                      " of ", format(n_total_scatter, big.mark = ","), " cells"),
        xref = "paper", yref = "paper", x = 1, y = 0,
        xanchor = "right", yanchor = "bottom",
        showarrow = FALSE,
        font = list(size = 10, color = "#888888")
      ))
    )
  }

  # Data summary
  data_summary <- create_scatter_data_summary(plot_data, x_var, y_var)

  list(plot = p, data_summary = data_summary, error = NULL)
}

create_enhanced_qc_plot <- function(seurat_obj, params) {
  qc_type <- params$qc_type
  
  # Prepare QC data
  meta_data <- seurat_obj@meta.data
  
  # Ensure basic QC metrics exist
  required_metrics <- c("nCount_RNA", "nFeature_RNA")
  missing_metrics <- required_metrics[!required_metrics %in% colnames(meta_data)]
  
  if (length(missing_metrics) > 0) {
    return(list(plot = NULL, data_summary = NULL, 
                error = paste("Missing QC metrics:", paste(missing_metrics, collapse = ", "))))
  }
  
  if (qc_type == "overview") {
    # Multi-panel QC overview
    qc_cols <- intersect(c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo"), 
                        colnames(meta_data))
    
    # Create violin plots for each metric
    plots <- list()
    
    for (i in seq_along(qc_cols)) {
      metric <- qc_cols[i]
      plot_data <- data.frame(
        Value = meta_data[[metric]],
        Metric = metric,
        Group = "All Cells"
      )
      
      plots[[i]] <- plot_ly(
        data = plot_data,
        x = ~Group,
        y = ~Value,
        type = "violin",
        name = metric,
        box = list(visible = TRUE)
      ) %>% layout(
        title = metric,
        showlegend = FALSE
      )
    }
    
    p <- subplot(plots, nrows = 2, shareY = FALSE) %>%
      layout(title = "QC Metrics Overview")
    
  } else if (qc_type == "umi_genes") {
    # UMI vs Genes scatter plot
    plot_data <- data.frame(
      UMI = meta_data$nCount_RNA,
      Genes = meta_data$nFeature_RNA,
      Cell = rownames(meta_data)
    )
    
    if (params$qc_color_by != "none" && params$qc_color_by %in% colnames(meta_data)) {
      plot_data$Color <- meta_data[[params$qc_color_by]]
      color_var <- ~Color
    } else {
      color_var <- NULL
    }
    
    p <- plot_ly(
      data = plot_data,
      x = ~UMI,
      y = ~Genes,
      color = color_var,
      type = "scatter",
      mode = "markers",
      marker = list(size = params$point_size * 2, opacity = params$point_alpha),
      text = ~Cell
    ) %>% layout(
      title = "UMI vs Genes",
      xaxis = list(title = "UMI Count", type = "log"),
      yaxis = list(title = "Gene Count", type = "log")
    )
    
  } else if (qc_type == "mt_percent") {
    # Mitochondrial percentage distribution
    if (!"percent.mt" %in% colnames(meta_data)) {
      return(list(plot = NULL, data_summary = NULL, 
                  error = "Mitochondrial percentage not calculated"))
    }
    
    plot_data <- data.frame(
      MT_Percent = meta_data$percent.mt,
      Cell = rownames(meta_data)
    )
    
    p <- plot_ly(
      data = plot_data,
      x = ~MT_Percent,
      type = "histogram",
      nbinsx = 50
    ) %>% layout(
      title = "Mitochondrial Gene Expression %",
      xaxis = list(title = "Mitochondrial %"),
      yaxis = list(title = "Number of Cells")
    )
    
  } else {
    return(list(plot = NULL, data_summary = NULL, error = "Unknown QC plot type"))
  }
  
  p <- apply_plot_theme(p, params, "QC Plot")
  
  # Data summary
  data_summary <- create_qc_data_summary(meta_data, qc_type)
  
  list(plot = p, data_summary = data_summary, error = NULL)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

get_enhanced_colors <- function(n, palette_name, continuous = FALSE) {
  if (continuous) {
    switch(palette_name,
      "viridis" = viridis::viridis(n),
      "plasma" = viridis::plasma(n),
      "RdYlBu" = RColorBrewer::brewer.pal(min(11, max(3, n)), "RdYlBu"),
      viridis::viridis(n)
    )
  } else {
    switch(palette_name,
      "Set3" = RColorBrewer::brewer.pal(min(12, max(3, n)), "Set3"),
      "Dark2" = RColorBrewer::brewer.pal(min(8, max(3, n)), "Dark2"),
      "Spectral" = RColorBrewer::brewer.pal(min(11, max(3, n)), "Spectral"),
      get_color_palette(n)
    )
  }
}

apply_plot_theme <- function(plot, params, title = "") {
  if (nzchar(params$plot_title)) {
    title <- params$plot_title
  }
  
  theme_list <- list(
    title = title,
    showlegend = params$show_legend,
    font = list(size = 12)
  )
  
  if (params$theme_choice == "dark") {
    theme_list$paper_bgcolor <- "black"
    theme_list$plot_bgcolor <- "black"
    theme_list$font$color <- "white"
  } else if (params$theme_choice == "minimal") {
    theme_list$paper_bgcolor <- "white"
    theme_list$plot_bgcolor <- "white"
  }
  
  plot %>% layout(theme_list)
}

# Data summary creation functions
create_data_summary <- function(plot_data, plot_type) {
  summary_stats <- data.frame(
    Metric = c("Total Points", "Variables", "Plot Type"),
    Value = c(nrow(plot_data), ncol(plot_data), plot_type)
  )
  
  DT::datatable(summary_stats, options = list(dom = 't'))
}

create_feature_data_summary <- function(expr_data, features) {
  summary_list <- list()
  
  for (feat in features) {
    feat_expr <- as.numeric(expr_data[feat, ])
    summary_list[[feat]] <- data.frame(
      Feature = feat,
      Mean_Expression = round(mean(feat_expr), 3),
      Max_Expression = round(max(feat_expr), 3),
      Percent_Expressing = round(sum(feat_expr > 0) / length(feat_expr) * 100, 1),
      stringsAsFactors = FALSE
    )
  }
  
  summary_df <- do.call(rbind, summary_list)
  DT::datatable(summary_df, options = list(dom = 't', pageLength = 10))
}

create_violin_data_summary <- function(plot_data) {
  summary_stats <- plot_data %>%
    group_by(Feature, Group) %>%
    summarise(
      Mean = round(mean(Expression), 3),
      Median = round(median(Expression), 3),
      SD = round(sd(Expression), 3),
      .groups = "drop"
    )
  
  DT::datatable(summary_stats, options = list(dom = 't', pageLength = 10))
}

create_dot_data_summary <- function(dot_data) {
  DT::datatable(dot_data, options = list(dom = 't', pageLength = 15))
}

create_heatmap_data_summary <- function(expr_matrix) {
  summary_df <- data.frame(
    Metric = c("Features", "Cells", "Mean Expression", "Max Expression", "Min Expression"),
    Value = c(nrow(expr_matrix), ncol(expr_matrix), 
             round(mean(expr_matrix), 3), round(max(expr_matrix), 3), round(min(expr_matrix), 3))
  )
  
  DT::datatable(summary_df, options = list(dom = 't'))
}

create_scatter_data_summary <- function(plot_data, x_var, y_var) {
  if (is.numeric(plot_data$X) && is.numeric(plot_data$Y)) {
    correlation <- cor(plot_data$X, plot_data$Y, use = "complete.obs")
    
    summary_df <- data.frame(
      Metric = c("X Variable", "Y Variable", "Data Points", "Correlation"),
      Value = c(x_var, y_var, nrow(plot_data), round(correlation, 3))
    )
  } else {
    summary_df <- data.frame(
      Metric = c("X Variable", "Y Variable", "Data Points"),
      Value = c(x_var, y_var, nrow(plot_data))
    )
  }
  
  DT::datatable(summary_df, options = list(dom = 't'))
}

create_qc_data_summary <- function(meta_data, qc_type) {
  if (qc_type == "overview") {
    qc_metrics <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
    qc_metrics <- qc_metrics[qc_metrics %in% colnames(meta_data)]
    
    summary_list <- list()
    for (metric in qc_metrics) {
      values <- meta_data[[metric]]
      summary_list[[metric]] <- data.frame(
        Metric = metric,
        Mean = round(mean(values, na.rm = TRUE), 2),
        Median = round(median(values, na.rm = TRUE), 2),
        SD = round(sd(values, na.rm = TRUE), 2)
      )
    }
    
    summary_df <- do.call(rbind, summary_list)
  } else {
    summary_df <- data.frame(
      Metric = c("Total Cells", "QC Type"),
      Value = c(nrow(meta_data), qc_type)
    )
  }
  
  DT::datatable(summary_df, options = list(dom = 't'))
}

# Additional plot types (volcano, correlation, box, ridge)
create_enhanced_volcano_plot <- function(seurat_obj, params) {
  tryCatch({
    # --- 1. Obtain differential expression results -------------------------
    de_results <- NULL

    # Check misc$de_results first
    if (!is.null(seurat_obj@misc$de_results)) {
      de_results <- seurat_obj@misc$de_results
    }

    # Check misc$markers (common alternative storage)
    if (is.null(de_results) && !is.null(seurat_obj@misc$markers)) {
      de_results <- seurat_obj@misc$markers
    }

    # If still NULL, run FindAllMarkers with default parameters
    if (is.null(de_results)) {
      # Need an identity class set; use active idents
      if (length(levels(Idents(seurat_obj))) < 2) {
        # Try to set idents from a clustering column
        cluster_cols <- grep("seurat_clusters|cluster|celltype|cell_type",
                             colnames(seurat_obj@meta.data),
                             ignore.case = TRUE, value = TRUE)
        if (length(cluster_cols) > 0) {
          Idents(seurat_obj) <- cluster_cols[1]
        } else {
          return(list(plot = NULL, data_summary = NULL,
                      error = "Volcano plot requires at least 2 identity classes for DE analysis"))
        }
      }
      de_results <- Seurat::FindAllMarkers(seurat_obj,
                                            only.pos = FALSE,
                                            min.pct = 0.1,
                                            logfc.threshold = 0.1,
                                            verbose = FALSE)
    }

    # Validate DE results have required columns
    required_cols <- c("avg_log2FC", "p_val_adj")
    # Handle older Seurat column names
    if ("avg_logFC" %in% colnames(de_results) && !"avg_log2FC" %in% colnames(de_results)) {
      de_results$avg_log2FC <- de_results$avg_logFC
    }
    if (!"gene" %in% colnames(de_results)) {
      if (!is.null(rownames(de_results))) {
        de_results$gene <- rownames(de_results)
      }
    }

    # Ensure numeric columns (guard against factors)
    if ("avg_log2FC" %in% colnames(de_results))
      de_results$avg_log2FC <- as.numeric(as.character(de_results$avg_log2FC))
    if ("p_val_adj" %in% colnames(de_results))
      de_results$p_val_adj <- as.numeric(as.character(de_results$p_val_adj))
    if ("p_val" %in% colnames(de_results))
      de_results$p_val <- as.numeric(as.character(de_results$p_val))

    missing_cols <- required_cols[!required_cols %in% colnames(de_results)]
    if (length(missing_cols) > 0) {
      return(list(plot = NULL, data_summary = NULL,
                  error = paste("DE results missing columns:", paste(missing_cols, collapse = ", "))))
    }

    # --- 2. Prepare plot data ------------------------------------------------
    fc_threshold  <- if (!is.null(params$fc_threshold))  params$fc_threshold  else 0.5
    pval_threshold <- if (!is.null(params$pval_threshold)) params$pval_threshold else 0.05

    plot_data <- data.frame(
      gene        = as.character(de_results$gene),
      log2FC      = de_results$avg_log2FC,
      pval_adj    = de_results$p_val_adj,
      stringsAsFactors = FALSE
    )

    # Replace zero p-values with smallest representable double for -log10
    plot_data$pval_adj[plot_data$pval_adj == 0] <- .Machine$double.xmin
    plot_data$neg_log10_pval <- -log10(plot_data$pval_adj)

    # Classify genes
    plot_data$significance <- "Not Significant"
    plot_data$significance[plot_data$log2FC >=  fc_threshold & plot_data$pval_adj < pval_threshold] <- "Upregulated"
    plot_data$significance[plot_data$log2FC <= -fc_threshold & plot_data$pval_adj < pval_threshold] <- "Downregulated"
    plot_data$significance <- factor(plot_data$significance,
                                     levels = c("Upregulated", "Downregulated", "Not Significant"))

    # Color map
    color_map <- c("Upregulated" = "#E41A1C",
                    "Downregulated" = "#377EB8",
                    "Not Significant" = "#999999")

    # --- 3. Build plotly scatter plot ----------------------------------------
    p <- plot_ly(
      data = plot_data,
      x = ~log2FC,
      y = ~neg_log10_pval,
      color = ~significance,
      colors = color_map,
      type = "scatter",
      mode = "markers",
      marker = list(size = 5, opacity = 0.7),
      text = ~gene,
      hovertemplate = paste0(
        "Gene: %{text}<br>",
        "Log2FC: %{x:.3f}<br>",
        "-log10(p_adj): %{y:.2f}<br>",
        "p_adj: %{customdata:.2e}<extra></extra>"
      ),
      customdata = ~pval_adj
    )

    # Threshold lines
    p <- p %>% layout(
      shapes = list(
        # Vertical FC thresholds
        list(type = "line", x0 =  fc_threshold, x1 =  fc_threshold,
             y0 = 0, y1 = max(plot_data$neg_log10_pval, na.rm = TRUE),
             line = list(dash = "dash", color = "grey50", width = 1)),
        list(type = "line", x0 = -fc_threshold, x1 = -fc_threshold,
             y0 = 0, y1 = max(plot_data$neg_log10_pval, na.rm = TRUE),
             line = list(dash = "dash", color = "grey50", width = 1)),
        # Horizontal p-value threshold
        list(type = "line", x0 = min(plot_data$log2FC, na.rm = TRUE),
             x1 = max(plot_data$log2FC, na.rm = TRUE),
             y0 = -log10(pval_threshold), y1 = -log10(pval_threshold),
             line = list(dash = "dash", color = "grey50", width = 1))
      ),
      xaxis = list(title = "Log2 Fold Change"),
      yaxis = list(title = "-log10(Adjusted P-value)")
    )

    # --- 4. Label specific genes if requested --------------------------------
    if (!is.null(params$volcano_genes) && nzchar(params$volcano_genes)) {
      label_genes <- unique(trimws(unlist(strsplit(params$volcano_genes, ","))))
      label_genes <- label_genes[label_genes %in% plot_data$gene]

      if (length(label_genes) > 0) {
        label_data <- plot_data[plot_data$gene %in% label_genes, , drop = FALSE]
        p <- p %>% add_annotations(
          x = label_data$log2FC,
          y = label_data$neg_log10_pval,
          text = label_data$gene,
          showarrow = TRUE,
          arrowhead = 2,
          arrowsize = 0.5,
          ax = 20,
          ay = -25,
          font = list(size = 10, color = "black")
        )
      }
    }

    p <- apply_plot_theme(p, params, "Volcano Plot")

    # --- 5. Data summary -----------------------------------------------------
    n_up   <- sum(plot_data$significance == "Upregulated")
    n_down <- sum(plot_data$significance == "Downregulated")
    n_ns   <- sum(plot_data$significance == "Not Significant")
    top_up <- plot_data[plot_data$significance == "Upregulated", ]
    top_up <- head(top_up[order(-top_up$log2FC), ], 5)
    top_down <- plot_data[plot_data$significance == "Downregulated", ]
    top_down <- head(top_down[order(top_down$log2FC), ], 5)

    summary_df <- data.frame(
      Metric = c("Total Genes Tested", "Upregulated", "Downregulated",
                  "Not Significant", "FC Threshold", "P-value Threshold",
                  paste0("Top Up: ", top_up$gene),
                  paste0("Top Down: ", top_down$gene)),
      Value = c(nrow(plot_data), n_up, n_down, n_ns,
                fc_threshold, pval_threshold,
                round(top_up$log2FC, 3),
                round(top_down$log2FC, 3)),
      stringsAsFactors = FALSE
    )
    data_summary <- DT::datatable(summary_df, options = list(dom = 't', pageLength = 20))

    list(plot = p, data_summary = data_summary, error = NULL)

  }, error = function(e) {
    list(plot = NULL, data_summary = NULL, error = paste("Volcano plot error:", e$message))
  })
}

create_enhanced_correlation_plot <- function(seurat_obj, params) {
  tryCatch({
    features <- parse_features(seurat_obj, params$features)

    if (length(features) < 2) {
      return(list(plot = NULL, data_summary = NULL,
                  error = "Correlation plot requires at least 2 valid features (comma-separated)"))
    }

    # Cap features at a reasonable number for readability
    features <- head(features, 50)

    # --- 1. Build expression matrix (genes x cells) then transpose ----------
    expr_data <- as.matrix(
      safe_get_assay_data(seurat_obj, "data")[features, , drop = FALSE]
    )
    # rows = genes, cols = cells  ->  transpose so rows = cells, cols = genes
    expr_matrix <- t(expr_data)

    # --- 2. Compute correlation matrix --------------------------------------
    cor_matrix <- cor(expr_matrix, use = "pairwise.complete.obs")

    # --- 3. Two-gene scatter plot shortcut ----------------------------------
    if (length(features) == 2) {
      scatter_data <- data.frame(
        Gene1 = expr_matrix[, 1],
        Gene2 = expr_matrix[, 2],
        Cell  = rownames(expr_matrix)
      )
      r_val <- cor_matrix[1, 2]

      p <- plot_ly(
        data = scatter_data,
        x = ~Gene1,
        y = ~Gene2,
        type = "scatter",
        mode = "markers",
        marker = list(size = 4, opacity = 0.6,
                       color = "steelblue"),
        text = ~Cell,
        hovertemplate = paste0(
          "Cell: %{text}<br>",
          features[1], ": %{x:.3f}<br>",
          features[2], ": %{y:.3f}<extra></extra>"
        )
      ) %>% layout(
        title = paste0(features[1], " vs ", features[2],
                       "  (r = ", round(r_val, 3), ")"),
        xaxis = list(title = features[1]),
        yaxis = list(title = features[2])
      )

      # Add trend line
      if (nrow(scatter_data) > 2) {
        fit <- lm(Gene2 ~ Gene1, data = scatter_data)
        pred_x <- seq(min(scatter_data$Gene1), max(scatter_data$Gene1), length.out = 100)
        pred_y <- predict(fit, newdata = data.frame(Gene1 = pred_x))
        p <- p %>% add_trace(
          x = pred_x, y = pred_y,
          type = "scatter", mode = "lines",
          line = list(color = "red", width = 2, dash = "dash"),
          showlegend = FALSE, inherit = FALSE
        )
      }

      p <- apply_plot_theme(p, params, "Gene Correlation")

      summary_df <- data.frame(
        Metric = c("Gene 1", "Gene 2", "Pearson r", "N cells"),
        Value  = c(features[1], features[2], round(r_val, 4), nrow(scatter_data)),
        stringsAsFactors = FALSE
      )
      data_summary <- DT::datatable(summary_df, options = list(dom = 't'))

      return(list(plot = p, data_summary = data_summary, error = NULL))
    }

    # --- 4. Multi-gene correlation heatmap ----------------------------------
    # Round for display
    cor_text <- matrix(sprintf("%.2f", cor_matrix), nrow = nrow(cor_matrix))

    p <- plot_ly(
      x = features,
      y = features,
      z = cor_matrix,
      type = "heatmap",
      colorscale = list(
        list(0,   "#2166AC"),
        list(0.25, "#67A9CF"),
        list(0.5, "#FFFFFF"),
        list(0.75, "#EF8A62"),
        list(1,   "#B2182B")
      ),
      zmin = -1,
      zmax = 1,
      text = cor_text,
      hovertemplate = paste0(
        "Gene X: %{x}<br>",
        "Gene Y: %{y}<br>",
        "Correlation: %{text}<extra></extra>"
      )
    ) %>% colorbar(title = "Pearson r") %>%
      layout(
        title = "Gene-Gene Correlation Matrix",
        xaxis = list(title = "", tickangle = -45),
        yaxis = list(title = "", autorange = "reversed")
      )

    p <- apply_plot_theme(p, params, "Correlation Heatmap")

    # --- 5. Data summary -----------------------------------------------------
    # Extract upper triangle correlations for summary table
    upper_idx <- which(upper.tri(cor_matrix), arr.ind = TRUE)
    cor_pairs <- data.frame(
      Gene1 = features[upper_idx[, 1]],
      Gene2 = features[upper_idx[, 2]],
      Correlation = round(cor_matrix[upper_idx], 4),
      stringsAsFactors = FALSE
    )
    cor_pairs <- cor_pairs[order(-abs(cor_pairs$Correlation)), ]
    data_summary <- DT::datatable(cor_pairs,
                                   options = list(dom = 't', pageLength = 20))

    list(plot = p, data_summary = data_summary, error = NULL)

  }, error = function(e) {
    list(plot = NULL, data_summary = NULL, error = paste("Correlation plot error:", e$message))
  })
}

create_enhanced_box_plot <- function(seurat_obj, params) {
  tryCatch({
    features <- parse_features(seurat_obj, params$features)

    if (length(features) == 0) {
      return(list(plot = NULL, data_summary = NULL, error = "No valid features specified"))
    }

    # Cap features for layout
    features <- head(features, 12)

    group_by <- if (!is.null(params$group_by) && params$group_by != "none") params$group_by else "all_cells"
    split_by <- if (!is.null(params$split_by) && params$split_by != "none") params$split_by else NULL

    # --- 1. Build long-format data frame -----------------------------------
    expr_data <- safe_get_assay_data(seurat_obj, "data")[features, , drop = FALSE]

    plot_data_list <- list()
    for (feat in features) {
      fd <- data.frame(
        Expression = as.numeric(expr_data[feat, ]),
        Feature    = feat,
        Cell       = colnames(seurat_obj),
        stringsAsFactors = FALSE
      )
      if (group_by != "all_cells" && group_by %in% colnames(seurat_obj@meta.data)) {
        fd$Group <- as.character(seurat_obj@meta.data[[group_by]])
      } else {
        fd$Group <- "All Cells"
      }
      if (!is.null(split_by) && split_by %in% colnames(seurat_obj@meta.data)) {
        fd$Split <- as.character(seurat_obj@meta.data[[split_by]])
      }
      plot_data_list[[feat]] <- fd
    }
    plot_data <- do.call(rbind, plot_data_list)

    colors <- get_enhanced_colors(length(unique(plot_data$Group)), params$color_palette)

    # --- 2. Create plotly box plots ----------------------------------------
    if (length(features) == 1) {
      # Single feature: one panel
      if (!is.null(split_by) && "Split" %in% colnames(plot_data)) {
        # Faceted by split_by
        split_levels <- unique(plot_data$Split)
        plots <- list()
        for (sl in split_levels) {
          sub_data <- plot_data[plot_data$Split == sl, ]
          plots[[sl]] <- plot_ly(
            data = sub_data,
            x = ~Group,
            y = ~Expression,
            color = ~Group,
            colors = colors,
            type = "box",
            boxpoints = "outliers",
            jitter = 0.3,
            pointpos = -1.5,
            marker = list(opacity = 0.5, size = 3),
            showlegend = (sl == split_levels[1])
          ) %>% layout(title = sl)
        }
        n_cols <- min(length(split_levels), 3)
        n_rows <- ceiling(length(split_levels) / n_cols)
        p <- subplot(plots, nrows = n_rows, shareY = TRUE) %>%
          layout(title = paste0(features[1], " (split by ", split_by, ")"))
      } else {
        p <- plot_ly(
          data = plot_data,
          x = ~Group,
          y = ~Expression,
          color = ~Group,
          colors = colors,
          type = "box",
          boxpoints = "outliers",
          jitter = 0.3,
          pointpos = -1.5,
          marker = list(opacity = 0.5, size = 3)
        ) %>% layout(
          title = features[1],
          xaxis = list(title = group_by),
          yaxis = list(title = "Expression Level")
        )
      }
    } else {
      # Multiple features: create subplot per feature
      plots <- list()
      for (i in seq_along(features)) {
        feat <- features[i]
        sub_data <- plot_data[plot_data$Feature == feat, ]
        plots[[i]] <- plot_ly(
          data = sub_data,
          x = ~Group,
          y = ~Expression,
          color = ~Group,
          colors = colors,
          type = "box",
          boxpoints = "outliers",
          jitter = 0.3,
          pointpos = -1.5,
          marker = list(opacity = 0.5, size = 3),
          showlegend = (i == 1)
        ) %>% layout(
          annotations = list(
            text = feat, x = 0.5, y = 1.05,
            xref = "paper", yref = "paper",
            showarrow = FALSE, font = list(size = 12)
          )
        )
      }
      n_cols <- min(length(features), 3)
      n_rows <- ceiling(length(features) / n_cols)
      p <- subplot(plots, nrows = n_rows, shareX = FALSE, shareY = FALSE) %>%
        layout(title = "Gene Expression Box Plots")
    }

    p <- apply_plot_theme(p, params, "Box Plot")

    # --- 3. Summary statistics ----------------------------------------------
    summary_stats <- plot_data %>%
      dplyr::group_by(Feature, Group) %>%
      dplyr::summarise(
        Mean   = round(mean(Expression, na.rm = TRUE), 3),
        Median = round(median(Expression, na.rm = TRUE), 3),
        SD     = round(sd(Expression, na.rm = TRUE), 3),
        Q1     = round(quantile(Expression, 0.25, na.rm = TRUE), 3),
        Q3     = round(quantile(Expression, 0.75, na.rm = TRUE), 3),
        N      = dplyr::n(),
        .groups = "drop"
      )
    data_summary <- DT::datatable(summary_stats,
                                   options = list(dom = 't', pageLength = 20))

    list(plot = p, data_summary = data_summary, error = NULL)

  }, error = function(e) {
    list(plot = NULL, data_summary = NULL, error = paste("Box plot error:", e$message))
  })
}

create_enhanced_ridge_plot <- function(seurat_obj, params) {
  tryCatch({
    features <- parse_features(seurat_obj, params$features)

    if (length(features) == 0) {
      return(list(plot = NULL, data_summary = NULL, error = "No valid features specified"))
    }

    # Cap features for layout
    features <- head(features, 9)

    group_by <- if (!is.null(params$group_by) && params$group_by != "none") params$group_by else "all_cells"

    # --- 1. Build expression data -------------------------------------------
    expr_data <- safe_get_assay_data(seurat_obj, "data")[features, , drop = FALSE]

    if (group_by != "all_cells" && group_by %in% colnames(seurat_obj@meta.data)) {
      groups <- as.character(seurat_obj@meta.data[[group_by]])
    } else {
      groups <- rep("All Cells", ncol(seurat_obj))
    }
    group_levels <- sort(unique(groups))
    colors <- get_enhanced_colors(length(group_levels), params$color_palette)
    names(colors) <- group_levels

    # --- 2. Build ridge plot(s) per feature ---------------------------------
    if (length(features) == 1) {
      feat <- features[1]
      feat_expr <- as.numeric(expr_data[feat, ])

      p <- plot_ly()
      offset_step <- 0
      # Compute a global y-range for offset spacing
      all_densities <- list()
      for (grp in group_levels) {
        vals <- feat_expr[groups == grp]
        if (length(vals[vals > 0]) > 2) {
          d <- density(vals, n = 256, from = min(feat_expr), to = max(feat_expr))
        } else {
          d <- density(vals, n = 256)
        }
        all_densities[[grp]] <- d
      }
      max_density <- max(sapply(all_densities, function(d) max(d$y)))
      offset_gap <- max_density * 0.8

      for (i in seq_along(group_levels)) {
        grp <- group_levels[i]
        d <- all_densities[[grp]]
        offset <- (i - 1) * offset_gap

        p <- p %>% add_trace(
          x = d$x,
          y = d$y + offset,
          type = "scatter",
          mode = "lines",
          fill = "toself",
          fillcolor = paste0(colors[grp], "55"),
          line = list(color = colors[grp], width = 1.5),
          name = grp,
          text = paste0("Group: ", grp, "<br>Expression: ",
                        round(d$x, 3), "<br>Density: ", round(d$y, 4)),
          hoverinfo = "text"
        )
      }

      # Add y-axis tick labels showing group names at their offset positions
      tick_vals <- (seq_along(group_levels) - 1) * offset_gap
      p <- p %>% layout(
        title = paste("Ridge Plot:", feat),
        xaxis = list(title = "Expression Level"),
        yaxis = list(title = "", tickvals = tick_vals,
                     ticktext = group_levels, zeroline = FALSE)
      )

    } else {
      # Multiple features: subplot per feature
      plots <- list()
      for (fi in seq_along(features)) {
        feat <- features[fi]
        feat_expr <- as.numeric(expr_data[feat, ])

        sub_p <- plot_ly()
        all_densities <- list()
        for (grp in group_levels) {
          vals <- feat_expr[groups == grp]
          if (length(vals[vals > 0]) > 2) {
            d <- density(vals, n = 256, from = min(feat_expr), to = max(feat_expr))
          } else {
            d <- density(vals, n = 256)
          }
          all_densities[[grp]] <- d
        }
        max_density <- max(sapply(all_densities, function(d) max(d$y)))
        offset_gap <- max_density * 0.8

        for (i in seq_along(group_levels)) {
          grp <- group_levels[i]
          d <- all_densities[[grp]]
          offset <- (i - 1) * offset_gap

          sub_p <- sub_p %>% add_trace(
            x = d$x,
            y = d$y + offset,
            type = "scatter",
            mode = "lines",
            fill = "toself",
            fillcolor = paste0(colors[grp], "55"),
            line = list(color = colors[grp], width = 1),
            name = grp,
            showlegend = (fi == 1),
            text = paste0("Group: ", grp, "<br>Expression: ",
                          round(d$x, 3), "<br>Density: ", round(d$y, 4)),
            hoverinfo = "text"
          )
        }

        tick_vals <- (seq_along(group_levels) - 1) * offset_gap
        sub_p <- sub_p %>% layout(
          xaxis = list(title = feat),
          yaxis = list(title = "", tickvals = tick_vals,
                       ticktext = group_levels, zeroline = FALSE)
        )
        plots[[fi]] <- sub_p
      }

      n_cols <- min(length(features), 3)
      n_rows <- ceiling(length(features) / n_cols)
      p <- subplot(plots, nrows = n_rows, shareX = FALSE, shareY = FALSE) %>%
        layout(title = "Ridge Plots")
    }

    p <- apply_plot_theme(p, params, "Ridge Plot")

    # --- 3. Summary statistics -----------------------------------------------
    summary_list <- list()
    for (feat in features) {
      feat_expr <- as.numeric(expr_data[feat, ])
      for (grp in group_levels) {
        vals <- feat_expr[groups == grp]
        summary_list[[paste(feat, grp, sep = "__")]] <- data.frame(
          Feature = feat,
          Group   = grp,
          Mean    = round(mean(vals, na.rm = TRUE), 3),
          Median  = round(median(vals, na.rm = TRUE), 3),
          SD      = round(sd(vals, na.rm = TRUE), 3),
          Pct_Expressing = round(sum(vals > 0) / length(vals) * 100, 1),
          N       = length(vals),
          stringsAsFactors = FALSE
        )
      }
    }
    summary_df <- do.call(rbind, summary_list)
    data_summary <- DT::datatable(summary_df,
                                   options = list(dom = 't', pageLength = 20))

    list(plot = p, data_summary = data_summary, error = NULL)

  }, error = function(e) {
    list(plot = NULL, data_summary = NULL, error = paste("Ridge plot error:", e$message))
  })
}

create_enhanced_density_plot <- function(seurat_obj, params) {
  # Create density plots for dimensional reductions
  reduction <- params$reduction
  if (!reduction %in% names(seurat_obj@reductions)) {
    reduction <- best_reduction(seurat_obj)
  }
  
  if (is.null(reduction)) {
    return(list(plot = NULL, data_summary = NULL, error = "No reduction available"))
  }
  
  embedding <- Embeddings(seurat_obj, reduction)[, 1:2]
  
  # Create 2D density plot
  p <- plot_ly(
    x = embedding[, 1],
    y = embedding[, 2],
    type = "histogram2dcontour",
    colorscale = "Viridis"
  ) %>% layout(
    title = paste(toupper(reduction), "Density Plot"),
    xaxis = list(title = paste0(toupper(reduction), "_1")),
    yaxis = list(title = paste0(toupper(reduction), "_2"))
  )
  
  p <- apply_plot_theme(p, params, "Density Plot")

  data_summary <- create_data_summary(data.frame(embedding), "Density plot")

  list(plot = p, data_summary = data_summary, error = NULL)
}

# =============================================================================
# MODULE 1: EXPRESSION EXPLORER
# =============================================================================

expressionExplorerUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Expression Explorer"),
          p("Explore gene expression across clusters, batches, and conditions with full user control."),

          # Gene and metadata selection
          fluidRow(
            column(4,
              textAreaInput(ns("explorer_genes"), "Genes (comma-separated)",
                placeholder = "CD3D, CD4, MS4A1, CD8A, PTPRC", rows = 2, width = "100%")
            ),
            column(4,
              selectInput(ns("explorer_group"), "Group By", choices = NULL),
              selectInput(ns("explorer_split"), "Split By (optional)", choices = c("None" = "none"))
            ),
            column(4,
              selectInput(ns("explorer_plot_type"), "Plot Type",
                choices = c(
                  "Violin Plot" = "violin",
                  "Box Plot" = "box",
                  "Dot Plot" = "dot",
                  "Ridge Plot" = "ridge",
                  "Bar Plot (Mean)" = "bar",
                  "Stacked Violin" = "stacked_violin",
                  "Grouped Bar (Condition x Cluster)" = "grouped_condition_bar"
                ), selected = "violin"),
              selectInput(ns("explorer_assay"), "Assay", choices = NULL)
            )
          ),

          # Subset and filter options
          fluidRow(
            column(3,
              selectInput(ns("subset_column"), "Subset by (optional)", choices = c("None" = "none"))
            ),
            column(3,
              selectInput(ns("subset_values"), "Subset values", choices = NULL, multiple = TRUE)
            ),
            column(3,
              checkboxInput(ns("show_points"), "Show data points", FALSE),
              checkboxInput(ns("log_scale"), "Log scale", FALSE)
            ),
            column(3,
              div(style = "margin-top: 25px;",
                actionButton(ns("run_explorer"), "Generate Plot",
                  class = "btn-primary btn-lg", width = "100%")
              )
            )
          )
        )
      )
    ),
    fluidRow(
      column(10, uiOutput(ns("explorer_plot_container"))),
      column(2,
        h5("Export"),
        numericInput(ns("exp_width"), "Width (in)", 10, min = 4, max = 20),
        numericInput(ns("exp_height"), "Height (in)", 7, min = 3, max = 15),
        numericInput(ns("exp_dpi"), "DPI", 300, min = 72, max = 600),
        downloadButton(ns("dl_explorer_png"), "PNG", class = "btn-success", width = "100%"),
        br(), br(),
        downloadButton(ns("dl_explorer_html"), "HTML", class = "btn-info", width = "100%"),
        br(), br(),
        downloadButton(ns("dl_explorer_csv"), "Data CSV", class = "btn-warning", width = "100%")
      )
    ),
    fluidRow(
      column(12,
        h5("Expression Statistics"),
        downloadButton(ns("dl_explorer_stats"), "CSV", class = "btn-xs btn-default"),
        DT::dataTableOutput(ns("explorer_stats"))
      )
    )
  )
}

expressionExplorerServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    explorer_plot_rv <- reactiveVal(NULL)
    explorer_data_rv <- reactiveVal(NULL)
    explorer_stats_rv <- reactiveVal(NULL)

    # --- Update UI choices when active_dataset changes -----------------------
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset

      meta_cols <- colnames(obj@meta.data)
      factor_cols <- meta_cols[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      assay_names <- names(obj@assays)
      default_assay <- safe_default_assay(obj)

      updateSelectInput(session, "explorer_group",
        choices = c("None" = "none", meta_cols),
        selected = if ("seurat_clusters" %in% meta_cols) "seurat_clusters" else "none")
      updateSelectInput(session, "explorer_split",
        choices = c("None" = "none", factor_cols))
      updateSelectInput(session, "subset_column",
        choices = c("None" = "none", factor_cols))
      updateSelectInput(session, "explorer_assay",
        choices = assay_names, selected = default_assay)
    })

    # --- Update subset_values when subset_column changes ---------------------
    observeEvent(input$subset_column, {
      req(values$active_dataset)
      if (!is.null(input$subset_column) && input$subset_column != "none") {
        vals <- unique(as.character(values$active_dataset@meta.data[[input$subset_column]]))
        updateSelectInput(session, "subset_values", choices = vals, selected = vals)
      } else {
        updateSelectInput(session, "subset_values", choices = character(0))
      }
    })

    # --- Main plot generation ------------------------------------------------
    observeEvent(input$run_explorer, {
      req(values$active_dataset)

      withProgress(message = "Generating expression plot...", value = 0, {
        tryCatch({
          obj <- values$active_dataset

          # Set assay if specified
          if (!is.null(input$explorer_assay) && input$explorer_assay %in% names(obj@assays)) {
            DefaultAssay(obj) <- input$explorer_assay
          }

          # Parse and validate genes
          genes_raw <- unique(trimws(unlist(strsplit(input$explorer_genes, "[,\n]+"))))
          genes_raw <- genes_raw[nzchar(genes_raw)]
          valid_genes <- genes_raw[genes_raw %in% rownames(obj)]

          if (length(valid_genes) == 0) {
            showNotification("No valid genes found in the dataset.", type = "error")
            return(NULL)
          }

          incProgress(0.2, detail = "Subsetting data...")

          # Optional subset
          if (!is.null(input$subset_column) && input$subset_column != "none" &&
              length(input$subset_values) > 0) {
            keep_cells <- rownames(obj@meta.data)[
              as.character(obj@meta.data[[input$subset_column]]) %in% input$subset_values
            ]
            if (length(keep_cells) < 10) {
              showNotification("Subset has fewer than 10 cells. Aborting.", type = "error")
              return(NULL)
            }
            obj <- subset(obj, cells = keep_cells)
          }

          incProgress(0.3, detail = "Extracting expression data...")

          # Get expression data
          expr_mat <- safe_get_assay_data(obj, "data")[valid_genes, , drop = FALSE]

          # Build long-format data frame
          group_col <- if (!is.null(input$explorer_group) && input$explorer_group != "none" &&
                          input$explorer_group %in% colnames(obj@meta.data)) {
            input$explorer_group
          } else {
            NULL
          }

          split_col <- if (!is.null(input$explorer_split) && input$explorer_split != "none" &&
                          input$explorer_split %in% colnames(obj@meta.data)) {
            input$explorer_split
          } else {
            NULL
          }

          long_list <- list()
          for (g in valid_genes) {
            df <- data.frame(
              Cell = colnames(obj),
              Gene = g,
              Expression = as.numeric(expr_mat[g, ]),
              stringsAsFactors = FALSE
            )
            if (!is.null(group_col)) {
              df$Group <- as.character(obj@meta.data[[group_col]])
            } else {
              df$Group <- "All Cells"
            }
            if (!is.null(split_col)) {
              df$Split <- as.character(obj@meta.data[[split_col]])
            }
            long_list[[g]] <- df
          }
          long_df <- do.call(rbind, long_list)
          rownames(long_df) <- NULL

          group_levels <- sort(unique(long_df$Group))
          n_groups <- length(group_levels)
          colors <- get_enhanced_colors(n_groups, "default")

          incProgress(0.5, detail = "Building plot...")

          plot_type <- input$explorer_plot_type

          # ---- VIOLIN ----
          if (plot_type == "violin") {
            if (length(valid_genes) == 1) {
              p <- plot_ly(data = long_df, x = ~Group, y = ~Expression,
                           color = ~Group, colors = colors, type = "violin",
                           box = list(visible = TRUE), meanline = list(visible = TRUE),
                           points = if (input$show_points) "all" else FALSE,
                           jitter = 0.3, pointpos = -1.5,
                           marker = list(size = 2, opacity = 0.4))
              p <- p %>% layout(title = paste("Violin:", valid_genes[1]),
                                xaxis = list(title = if (!is.null(group_col)) group_col else ""),
                                yaxis = list(title = "Expression"))
            } else {
              has_split <- !is.null(split_col) && "Split" %in% colnames(long_df)
              if (has_split) {
                split_levels <- sort(unique(long_df$Split))
                plots <- list()
                idx <- 0
                first_shown <- FALSE
                for (g in valid_genes) {
                  for (sl in split_levels) {
                    idx <- idx + 1
                    sub_df <- long_df[long_df$Gene == g & long_df$Split == sl, ]
                    plots[[idx]] <- plot_ly(data = sub_df, x = ~Group, y = ~Expression,
                                            color = ~Group, colors = colors, type = "violin",
                                            box = list(visible = TRUE), meanline = list(visible = TRUE),
                                            points = if (input$show_points) "all" else FALSE,
                                            jitter = 0.3, pointpos = -1.5,
                                            marker = list(size = 2, opacity = 0.4),
                                            showlegend = !first_shown) %>%
                      layout(annotations = list(text = paste(g, "-", sl), x = 0.5, y = 1.05,
                                                xref = "paper", yref = "paper",
                                                showarrow = FALSE, font = list(size = 11)))
                    first_shown <- TRUE
                  }
                }
                nc <- min(length(split_levels), 4)
                nr <- ceiling(length(plots) / nc)
                p <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE, margin = 0.06) %>%
                  layout(title = paste("Violin Plots split by", split_col))
              } else {
                plots <- list()
                for (i in seq_along(valid_genes)) {
                  g <- valid_genes[i]
                  sub_df <- long_df[long_df$Gene == g, ]
                  plots[[i]] <- plot_ly(data = sub_df, x = ~Group, y = ~Expression,
                                        color = ~Group, colors = colors, type = "violin",
                                        box = list(visible = TRUE), meanline = list(visible = TRUE),
                                        points = if (input$show_points) "all" else FALSE,
                                        jitter = 0.3, pointpos = -1.5,
                                        marker = list(size = 2, opacity = 0.4),
                                        showlegend = (i == 1)) %>%
                    layout(annotations = list(text = g, x = 0.5, y = 1.05,
                                              xref = "paper", yref = "paper",
                                              showarrow = FALSE, font = list(size = 11)))
                }
                nc <- min(length(valid_genes), 3)
                nr <- ceiling(length(valid_genes) / nc)
                p <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE, margin = 0.06) %>%
                  layout(title = "Violin Plots")
              }
            }

          # ---- BOX ----
          } else if (plot_type == "box") {
            if (length(valid_genes) == 1) {
              p <- plot_ly(data = long_df, x = ~Group, y = ~Expression,
                           color = ~Group, colors = colors, type = "box",
                           boxpoints = if (input$show_points) "all" else "outliers",
                           jitter = 0.3, pointpos = -1.5,
                           marker = list(size = 2, opacity = 0.4))
              p <- p %>% layout(title = paste("Box Plot:", valid_genes[1]),
                                xaxis = list(title = if (!is.null(group_col)) group_col else ""),
                                yaxis = list(title = "Expression"))
            } else {
              has_split <- !is.null(split_col) && "Split" %in% colnames(long_df)
              if (has_split) {
                split_levels <- sort(unique(long_df$Split))
                plots <- list()
                idx <- 0
                first_shown <- FALSE
                for (g in valid_genes) {
                  for (sl in split_levels) {
                    idx <- idx + 1
                    sub_df <- long_df[long_df$Gene == g & long_df$Split == sl, ]
                    plots[[idx]] <- plot_ly(data = sub_df, x = ~Group, y = ~Expression,
                                            color = ~Group, colors = colors, type = "box",
                                            boxpoints = if (input$show_points) "all" else "outliers",
                                            jitter = 0.3, pointpos = -1.5,
                                            marker = list(size = 2, opacity = 0.4),
                                            showlegend = !first_shown) %>%
                      layout(annotations = list(text = paste(g, "-", sl), x = 0.5, y = 1.05,
                                                xref = "paper", yref = "paper",
                                                showarrow = FALSE, font = list(size = 11)))
                    first_shown <- TRUE
                  }
                }
                nc <- min(length(split_levels), 4)
                nr <- ceiling(length(plots) / nc)
                p <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE, margin = 0.06) %>%
                  layout(title = paste("Box Plots split by", split_col))
              } else {
                plots <- list()
                for (i in seq_along(valid_genes)) {
                  g <- valid_genes[i]
                  sub_df <- long_df[long_df$Gene == g, ]
                  plots[[i]] <- plot_ly(data = sub_df, x = ~Group, y = ~Expression,
                                        color = ~Group, colors = colors, type = "box",
                                        boxpoints = if (input$show_points) "all" else "outliers",
                                        jitter = 0.3, pointpos = -1.5,
                                        marker = list(size = 2, opacity = 0.4),
                                        showlegend = (i == 1)) %>%
                    layout(annotations = list(text = g, x = 0.5, y = 1.05,
                                              xref = "paper", yref = "paper",
                                              showarrow = FALSE, font = list(size = 11)))
                }
                nc <- min(length(valid_genes), 3)
                nr <- ceiling(length(valid_genes) / nc)
                p <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE, margin = 0.06) %>%
                  layout(title = "Box Plots")
              }
            }

          # ---- DOT ----
          } else if (plot_type == "dot") {
            dot_list <- list()
            for (g in valid_genes) {
              for (grp in group_levels) {
                vals <- long_df$Expression[long_df$Gene == g & long_df$Group == grp]
                dot_list[[paste(g, grp, sep = "___")]] <- data.frame(
                  Gene = g, Group = grp,
                  Avg_Expression = mean(vals),
                  Pct_Expressing = sum(vals > 0) / length(vals) * 100,
                  stringsAsFactors = FALSE
                )
              }
            }
            dot_df <- do.call(rbind, dot_list)
            rownames(dot_df) <- NULL

            p <- plot_ly(data = dot_df, x = ~Gene, y = ~Group,
                         color = ~Avg_Expression, size = ~Pct_Expressing,
                         colors = viridis::viridis(100),
                         type = "scatter", mode = "markers",
                         marker = list(opacity = 0.85, sizemode = "diameter",
                                       sizeref = max(dot_df$Pct_Expressing) / 30),
                         hovertemplate = paste0("Gene: %{x}<br>Group: %{y}<br>",
                                               "Avg Expr: %{marker.color:.3f}<br>",
                                               "%% Expressing: %{marker.size:.1f}%%<extra></extra>")) %>%
              layout(title = "Dot Plot",
                     xaxis = list(title = "", tickangle = -45),
                     yaxis = list(title = "")) %>%
              colorbar(title = "Avg Expr")

          # ---- RIDGE ----
          } else if (plot_type == "ridge") {
            if (length(valid_genes) == 1) {
              feat_expr <- as.numeric(expr_mat[valid_genes[1], ])
              groups_vec <- if (!is.null(group_col)) as.character(obj@meta.data[[group_col]]) else rep("All Cells", ncol(obj))

              p <- plot_ly()
              all_dens <- list()
              for (grp in group_levels) {
                vals <- feat_expr[groups_vec == grp]
                all_dens[[grp]] <- if (length(vals[vals > 0]) > 2) {
                  density(vals, n = 256, from = min(feat_expr), to = max(feat_expr))
                } else {
                  density(vals, n = 256)
                }
              }
              max_d <- max(sapply(all_dens, function(d) max(d$y)))
              gap <- max_d * 0.85

              for (i in seq_along(group_levels)) {
                grp <- group_levels[i]
                d <- all_dens[[grp]]
                offset <- (i - 1) * gap
                p <- p %>% add_trace(
                  x = d$x, y = d$y + offset, type = "scatter", mode = "lines",
                  fill = "toself", fillcolor = paste0(colors[i], "55"),
                  line = list(color = colors[i], width = 1.5), name = grp,
                  hoverinfo = "text",
                  text = paste0("Group: ", grp, "<br>Expr: ", round(d$x, 3),
                                "<br>Density: ", round(d$y, 4))
                )
              }
              tick_vals <- (seq_along(group_levels) - 1) * gap
              p <- p %>% layout(title = paste("Ridge Plot:", valid_genes[1]),
                                xaxis = list(title = "Expression"),
                                yaxis = list(title = "", tickvals = tick_vals,
                                             ticktext = group_levels, zeroline = FALSE))
            } else {
              plots <- list()
              groups_vec <- if (!is.null(group_col)) as.character(obj@meta.data[[group_col]]) else rep("All Cells", ncol(obj))
              for (fi in seq_along(valid_genes)) {
                feat <- valid_genes[fi]
                feat_expr <- as.numeric(expr_mat[feat, ])
                sub_p <- plot_ly()
                all_dens <- list()
                for (grp in group_levels) {
                  vals <- feat_expr[groups_vec == grp]
                  all_dens[[grp]] <- if (length(vals[vals > 0]) > 2) {
                    density(vals, n = 256, from = min(feat_expr), to = max(feat_expr))
                  } else {
                    density(vals, n = 256)
                  }
                }
                max_d <- max(sapply(all_dens, function(d) max(d$y)))
                gap <- max_d * 0.85
                for (i in seq_along(group_levels)) {
                  grp <- group_levels[i]
                  d <- all_dens[[grp]]
                  offset <- (i - 1) * gap
                  sub_p <- sub_p %>% add_trace(
                    x = d$x, y = d$y + offset, type = "scatter", mode = "lines",
                    fill = "toself", fillcolor = paste0(colors[i], "55"),
                    line = list(color = colors[i], width = 1), name = grp,
                    showlegend = (fi == 1), hoverinfo = "text",
                    text = paste0("Group: ", grp, "<br>Expr: ", round(d$x, 3),
                                  "<br>Density: ", round(d$y, 4))
                  )
                }
                tick_vals <- (seq_along(group_levels) - 1) * gap
                sub_p <- sub_p %>% layout(
                  xaxis = list(title = feat),
                  yaxis = list(title = "", tickvals = tick_vals,
                               ticktext = group_levels, zeroline = FALSE))
                plots[[fi]] <- sub_p
              }
              nc <- min(length(valid_genes), 3)
              nr <- ceiling(length(valid_genes) / nc)
              p <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE, margin = 0.06) %>%
                layout(title = "Ridge Plots")
            }

          # ---- BAR (MEAN) ----
          } else if (plot_type == "bar") {
            bar_list <- list()
            for (g in valid_genes) {
              for (grp in group_levels) {
                vals <- long_df$Expression[long_df$Gene == g & long_df$Group == grp]
                n <- length(vals)
                m <- mean(vals)
                se <- if (n > 1) sd(vals) / sqrt(n) else 0
                bar_list[[paste(g, grp, sep = "___")]] <- data.frame(
                  Gene = g, Group = grp, Mean = m, SEM = se, N = n,
                  stringsAsFactors = FALSE
                )
              }
            }
            bar_df <- do.call(rbind, bar_list)
            rownames(bar_df) <- NULL

            if (length(valid_genes) == 1) {
              p <- plot_ly(data = bar_df, x = ~Group, y = ~Mean,
                           color = ~Group, colors = colors, type = "bar",
                           error_y = list(type = "data", array = ~SEM, visible = TRUE),
                           hovertemplate = paste0("Group: %{x}<br>Mean: %{y:.3f}<br>",
                                                  "SEM: %{error_y.array:.3f}<extra></extra>")) %>%
                layout(title = paste("Mean Expression:", valid_genes[1]),
                       xaxis = list(title = if (!is.null(group_col)) group_col else ""),
                       yaxis = list(title = "Mean Expression"))
            } else {
              plots <- list()
              for (i in seq_along(valid_genes)) {
                g <- valid_genes[i]
                sub_df <- bar_df[bar_df$Gene == g, ]
                plots[[i]] <- plot_ly(data = sub_df, x = ~Group, y = ~Mean,
                                      color = ~Group, colors = colors, type = "bar",
                                      error_y = list(type = "data", array = ~SEM, visible = TRUE),
                                      showlegend = (i == 1),
                                      hovertemplate = paste0("Group: %{x}<br>Mean: %{y:.3f}<br>",
                                                             "SEM: %{error_y.array:.3f}<extra></extra>")) %>%
                  layout(annotations = list(text = g, x = 0.5, y = 1.05,
                                            xref = "paper", yref = "paper",
                                            showarrow = FALSE, font = list(size = 11)))
              }
              nc <- min(length(valid_genes), 3)
              nr <- ceiling(length(valid_genes) / nc)
              p <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE, margin = 0.06) %>%
                layout(title = "Mean Expression Bar Plots")
            }

          # ---- STACKED VIOLIN ----
          } else if (plot_type == "stacked_violin") {
            plots <- list()
            for (i in seq_along(valid_genes)) {
              g <- valid_genes[i]
              sub_df <- long_df[long_df$Gene == g, ]
              plots[[i]] <- plot_ly(data = sub_df, x = ~Group, y = ~Expression,
                                    color = ~Group, colors = colors, type = "violin",
                                    box = list(visible = TRUE), meanline = list(visible = TRUE),
                                    points = if (input$show_points) "all" else FALSE,
                                    jitter = 0.3, pointpos = -1.5,
                                    marker = list(size = 2, opacity = 0.4),
                                    showlegend = (i == 1),
                                    scalemode = "width") %>%
                layout(yaxis = list(title = g))
            }
            p <- subplot(plots, nrows = length(valid_genes), shareX = TRUE, shareY = FALSE, margin = 0.04) %>%
              layout(title = "Stacked Violin Plots")

          } else if (plot_type == "grouped_condition_bar") {
            # Grouped bar chart: mean expression + SEM per Gene x Group x Condition
            if (is.null(split_col) || !"Split" %in% colnames(long_df)) {
              showNotification("Grouped bar plot requires 'Split By' to be set (e.g., condition, batch).", type = "warning")
              return(NULL)
            }

            split_levels <- sort(unique(long_df$Split))
            n_splits <- length(split_levels)
            split_colors <- get_enhanced_colors(n_splits, "default")

            if (length(valid_genes) == 1) {
              agg <- long_df %>%
                dplyr::group_by(Group, Split) %>%
                dplyr::summarise(
                  Mean = mean(Expression, na.rm = TRUE),
                  SEM = sd(Expression, na.rm = TRUE) / sqrt(dplyr::n()),
                  .groups = "drop"
                )
              agg$Group <- as.character(agg$Group)
              agg$Split <- as.character(agg$Split)

              p <- plot_ly()
              for (j in seq_along(split_levels)) {
                s <- split_levels[j]
                sub <- agg[agg$Split == s, ]
                p <- p %>% add_trace(
                  data = sub, x = ~Group, y = ~Mean, type = "bar",
                  name = s, marker = list(color = split_colors[j]),
                  error_y = list(type = "data", array = ~SEM, visible = TRUE),
                  hovertemplate = paste0(s, "<br>%{x}<br>Mean: %{y:.3f}<extra></extra>")
                )
              }
              p <- p %>% layout(
                title = paste("Grouped Bar:", valid_genes[1]),
                barmode = "group",
                xaxis = list(title = if (!is.null(group_col)) group_col else "",
                             categoryorder = "array", categoryarray = group_levels),
                yaxis = list(title = "Mean Expression"),
                legend = list(title = list(text = split_col))
              )
            } else {
              # Multi-gene: subplot grid
              plots <- list()
              n_cols_sub <- min(2, length(valid_genes))
              for (i in seq_along(valid_genes)) {
                g <- valid_genes[i]
                sub_df <- long_df[long_df$Gene == g, ]
                agg <- sub_df %>%
                  dplyr::group_by(Group, Split) %>%
                  dplyr::summarise(
                    Mean = mean(Expression, na.rm = TRUE),
                    SEM = sd(Expression, na.rm = TRUE) / sqrt(dplyr::n()),
                    .groups = "drop"
                  )
                agg$Group <- as.character(agg$Group)
                agg$Split <- as.character(agg$Split)

                pi <- plot_ly()
                for (j in seq_along(split_levels)) {
                  s <- split_levels[j]
                  sub <- agg[agg$Split == s, ]
                  pi <- pi %>% add_trace(
                    data = sub, x = ~Group, y = ~Mean, type = "bar",
                    name = s, marker = list(color = split_colors[j]),
                    error_y = list(type = "data", array = ~SEM, visible = TRUE),
                    showlegend = (i == 1),
                    hovertemplate = paste0(g, " | ", s, "<br>%{x}<br>Mean: %{y:.3f}<extra></extra>")
                  )
                }
                pi <- pi %>% layout(
                  barmode = "group",
                  xaxis = list(title = ""),
                  yaxis = list(title = g),
                  annotations = list(list(
                    text = g, xref = "paper", yref = "paper",
                    x = 0.5, y = 1.05, showarrow = FALSE, font = list(size = 13)
                  ))
                )
                plots[[i]] <- pi
              }
              n_rows_sub <- ceiling(length(valid_genes) / n_cols_sub)
              p <- subplot(plots, nrows = n_rows_sub, shareX = TRUE, margin = 0.08) %>%
                layout(title = "Grouped Bar Plot by Condition",
                       barmode = "group")
            }

          } else {
            showNotification("Unknown plot type.", type = "error")
            return(NULL)
          }

          # Log scale
          if (input$log_scale) {
            p <- p %>% layout(yaxis = list(type = "log"))
          }

          explorer_plot_rv(p)
          explorer_data_rv(long_df)

          # Build stats table
          stats_df <- long_df %>%
            dplyr::group_by(Gene, Group) %>%
            dplyr::summarise(
              Mean = round(mean(Expression, na.rm = TRUE), 4),
              Median = round(median(Expression, na.rm = TRUE), 4),
              SD = round(sd(Expression, na.rm = TRUE), 4),
              Pct_Expressing = round(sum(Expression > 0) / dplyr::n() * 100, 1),
              N_Cells = dplyr::n(),
              .groups = "drop"
            )
          explorer_stats_rv(stats_df)

          incProgress(1, detail = "Done!")
          showNotification("Plot generated successfully.", type = "message")

        }, error = function(e) {
          showNotification(paste("Error:", e$message), type = "error")
        })
      })
    })

    # --- Render outputs with dynamic height -----------------------------------

    output$explorer_plot_container <- renderUI({
      n_genes <- length(input$explorer_genes)
      # Account for split-by producing extra panels
      has_split <- !is.null(input$explorer_split) && input$explorer_split != "none"
      n_splits <- 1
      if (has_split && !is.null(values$active_dataset)) {
        tryCatch({
          sc <- input$explorer_split
          if (sc %in% colnames(values$active_dataset@meta.data)) {
            n_splits <- length(unique(values$active_dataset@meta.data[[sc]]))
          }
        }, error = function(e) NULL)
      }
      n_panels <- max(1, n_genes) * n_splits
      nc <- if (has_split) min(n_splits, 4) else min(max(1, n_genes), 3)
      n_rows <- ceiling(n_panels / nc)
      h <- max(650, n_rows * 350)
      plotlyOutput(session$ns("explorer_plot"), height = paste0(h, "px"))
    })

    output$explorer_plot <- renderPlotly({
      req(explorer_plot_rv())
      p <- explorer_plot_rv() %>%
        layout(autosize = TRUE) %>%
        plotly::partial_bundle()
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Visualization", figure_type = "Expression Explorer",
                       label = "Expression Explorer", source_module = "Visualization")
      }
      p
    })

    output$explorer_stats <- DT::renderDT({
      req(explorer_stats_rv())
      tbl <- explorer_stats_rv()
      if (exists("register_csv_export", mode = "function")) {
        register_csv_export(values, data = tbl, description = "Expression Explorer Stats", source_module = "Visualization")
      }
      DT::datatable(tbl, options = list(
        pageLength = 20, scrollX = TRUE, dom = "frtip"
      ))
    })

    # --- Download handlers ---------------------------------------------------
    output$dl_explorer_png <- downloadHandler(
      filename = function() paste0("expression_explorer_", Sys.Date(), ".png"),
      content = function(file) {
        tryCatch({
          p <- explorer_plot_rv()
          if (is.null(p)) return(NULL)
          # Try orca first, then webshot2 fallback
          ok <- TRUE
          tryCatch({
            plotly::orca(p, file = file,
                         width = input$exp_width * input$exp_dpi,
                         height = input$exp_height * input$exp_dpi)
          }, error = function(e) { ok <<- FALSE })
          if (!ok) {
            tmp_html <- tempfile(fileext = ".html")
            htmlwidgets::saveWidget(p, tmp_html, selfcontained = TRUE)
            if (requireNamespace("webshot2", quietly = TRUE)) {
              webshot2::webshot(tmp_html, file = file,
                               vwidth = input$exp_width * input$exp_dpi,
                               vheight = input$exp_height * input$exp_dpi)
            } else {
              showNotification("PNG export requires 'webshot2' package. Exporting HTML instead.",
                               type = "warning")
              file.copy(tmp_html, sub("\\.png$", ".html", file), overwrite = TRUE)
            }
          }
        }, error = function(e) {
          showNotification(paste("PNG export error:", e$message), type = "error")
        })
      }
    )

    output$dl_explorer_html <- downloadHandler(
      filename = function() paste0("expression_explorer_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- explorer_plot_rv()
          if (!is.null(p)) {
            htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
          }
        }, error = function(e) {
          showNotification(paste("HTML export error:", e$message), type = "error")
        })
      }
    )

    output$dl_explorer_csv <- downloadHandler(
      filename = function() paste0("expression_data_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          df <- explorer_data_rv()
          if (!is.null(df)) {
            write.csv(df, file, row.names = FALSE)
          }
        }, error = function(e) {
          showNotification(paste("CSV export error:", e$message), type = "error")
        })
      }
    )

    output$dl_explorer_stats <- downloadHandler(
      filename = function() paste0("expression_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        req(explorer_stats_rv())
        write.csv(explorer_stats_rv(), file, row.names = FALSE)
      }
    )

  })
}


# =============================================================================
# MODULE 2: CUSTOM HEATMAP BUILDER
# =============================================================================

heatmapBuilderUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Custom Heatmap Builder"),
          p("Build publication-ready heatmaps with full control over genes, groups, clustering, and colors."),

          fluidRow(
            column(6,
              selectInput(ns("heatmap_preset"), "Gene Presets",
                choices = c("Custom" = "custom",
                            "Immune Markers" = "immune",
                            "T Cell Subtypes" = "tcell",
                            "Cell Cycle" = "cell_cycle",
                            "Stress/Apoptosis" = "stress",
                            "Metabolism" = "metabolism"),
                selected = "custom"),
              textAreaInput(ns("heatmap_genes"), "Genes (comma-separated or one per line)",
                placeholder = "CD3D, CD4, CD8A, MS4A1, CD14, FCGR3A, NKG7, PPBP", rows = 4, width = "100%"),
              fluidRow(
                column(6,
                  actionButton(ns("find_top_markers"), "Find Top N Markers", class = "btn-info btn-sm", width = "100%")
                ),
                column(6,
                  numericInput(ns("top_n_markers"), "Top N per cluster", 10, min = 3, max = 50)
                )
              ),
              fileInput(ns("heatmap_csv"), "Upload genes from CSV/TXT", accept = c(".csv", ".txt"))
            ),
            column(6,
              selectInput(ns("heatmap_group"), "Group cells by", choices = NULL),
              selectInput(ns("heatmap_assay"), "Assay", choices = NULL),
              selectInput(ns("heatmap_slot"), "Data slot",
                choices = c("Scaled" = "scale.data", "Normalized" = "data", "Counts" = "counts"),
                selected = "scale.data")
            )
          ),

          fluidRow(
            column(3,
              selectInput(ns("heatmap_scale"), "Scaling",
                choices = c("Row (gene)" = "row", "Column (cell)" = "column", "None" = "none"),
                selected = "row")
            ),
            column(3,
              selectInput(ns("heatmap_colors"), "Color scheme",
                choices = c("Blue-White-Red" = "bwr", "Viridis" = "viridis",
                            "Purple-Orange" = "PuOr", "Yellow-Red" = "YlOrRd",
                            "Blue-Yellow" = "by"), selected = "bwr")
            ),
            column(3,
              checkboxInput(ns("cluster_rows"), "Cluster genes", TRUE),
              checkboxInput(ns("cluster_cols"), "Cluster cells", FALSE)
            ),
            column(3,
              checkboxInput(ns("show_gene_names"), "Show gene names", TRUE),
              checkboxInput(ns("average_expression"), "Average per group", TRUE)
            )
          ),

          fluidRow(
            column(4, numericInput(ns("hm_max_cells"), "Max cells to show", 2000, min = 100, max = 10000)),
            column(4, numericInput(ns("hm_fontsize"), "Font size", 10, min = 6, max = 20)),
            column(4,
              actionButton(ns("build_heatmap"), "Build Heatmap",
                class = "btn-danger btn-lg", width = "100%")
            )
          )
        )
      )
    ),
    fluidRow(
      column(12,
        plotOutput(ns("heatmap_output"), height = "700px"),
        br(),
        fluidRow(
          column(4, downloadButton(ns("dl_hm_png"), "PNG", class = "btn-success")),
          column(4, downloadButton(ns("dl_hm_pdf"), "PDF", class = "btn-info")),
          column(4, downloadButton(ns("dl_hm_csv"), "Data CSV", class = "btn-warning"))
        )
      )
    )
  )
}

heatmapBuilderServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    heatmap_obj_rv <- reactiveVal(NULL)
    heatmap_matrix_rv <- reactiveVal(NULL)

    # --- Update choices on dataset change ------------------------------------
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset

      meta_cols <- colnames(obj@meta.data)
      factor_cols <- meta_cols[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      # Also include numeric columns that look like clusters
      cluster_pat <- grep("cluster|leiden|louvain|seurat_clusters", meta_cols, value = TRUE, ignore.case = TRUE)
      all_group_cols <- sort(unique(union(factor_cols, cluster_pat)))
      assay_names <- names(obj@assays)
      default_assay <- safe_default_assay(obj)

      default_sel <- if ("seurat_clusters" %in% all_group_cols) "seurat_clusters" else "none"
      updateSelectInput(session, "heatmap_group",
        choices = c("None" = "none", all_group_cols),
        selected = default_sel)
      updateSelectInput(session, "heatmap_assay",
        choices = assay_names, selected = default_assay)
    })

    # --- Gene presets ---------------------------------------------------------
    observeEvent(input$heatmap_preset, {
      genes <- switch(input$heatmap_preset,
        "immune" = "CD3D, CD3E, CD4, CD8A, CD8B, MS4A1, CD19, CD14, CD16, FCGR3A, NKG7, GNLY, NCAM1, CD34, PPBP, CD68, CSF1R, ITGAX, CD1C, CLEC4C",
        "tcell" = "CD3D, CD3E, CD4, CD8A, CD8B, FOXP3, IL2RA, CTLA4, PDCD1, LAG3, HAVCR2, TOX, TCF7, LEF1, CCR7, SELL, GZMB, PRF1, IFNG, IL17A",
        "cell_cycle" = "MKI67, TOP2A, PCNA, MCM6, CDK1, CCNB1, CCND1, CCNE1, CDK4, RB1, TP53, CDKN1A, AURKA, PLK1, BUB1",
        "stress" = "HSPA1A, HSPA1B, HSP90AA1, DDIT3, ATF4, XBP1, BCL2, BAX, CASP3, CASP9, FOS, JUN, NFE2L2, HIF1A, NFKB1",
        "metabolism" = "HK1, HK2, PKM, LDHA, PDK1, IDH1, IDH2, OGDH, CS, ACO2, SDHA, FH, MDH2, G6PD, TKT",
        NULL
      )
      if (!is.null(genes)) {
        updateTextAreaInput(session, "heatmap_genes", value = genes)
      }
    })

    # --- Find Top Markers per cluster -----------------------------------------
    observeEvent(input$find_top_markers, {
      req(values$active_dataset)
      withProgress(message = "Finding top markers...", value = 0, {
        tryCatch({
          obj <- values$active_dataset
          obj <- prepare_for_de(obj)

          group_col <- input$heatmap_group
          if (!is.null(group_col) && group_col != "none") {
            Idents(obj) <- group_col
          }

          incProgress(0.2, detail = "Running FindAllMarkers...")
          markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
          top_n <- input$top_n_markers
          if (is.null(top_n) || is.na(top_n)) top_n <- 10

          top_genes <- markers %>%
            dplyr::group_by(cluster) %>%
            dplyr::slice_max(order_by = avg_log2FC, n = top_n) %>%
            dplyr::pull(gene) %>% unique()

          updateTextAreaInput(session, "heatmap_genes", value = paste(top_genes, collapse = ", "))
          showNotification(paste("Found", length(top_genes), "top marker genes"), type = "message")
        }, error = function(e) {
          showNotification(paste("Marker finding error:", e$message), type = "error")
        })
      })
    })

    # --- CSV/TXT gene upload --------------------------------------------------
    observeEvent(input$heatmap_csv, {
      req(input$heatmap_csv)
      tryCatch({
        df <- read.csv(input$heatmap_csv$datapath, header = TRUE, stringsAsFactors = FALSE)
        genes <- unique(trimws(unlist(df[, 1])))
        genes <- genes[nzchar(genes)]
        updateTextAreaInput(session, "heatmap_genes", value = paste(genes, collapse = ", "))
        showNotification(paste("Loaded", length(genes), "genes from file"), type = "message")
      }, error = function(e) {
        showNotification(paste("File read error:", e$message), type = "error")
      })
    })

    # --- Build heatmap -------------------------------------------------------
    observeEvent(input$build_heatmap, {
      req(values$active_dataset)

      withProgress(message = "Building heatmap...", value = 0, {
        tryCatch({
          obj <- values$active_dataset

          # Set assay
          if (!is.null(input$heatmap_assay) && input$heatmap_assay %in% names(obj@assays)) {
            DefaultAssay(obj) <- input$heatmap_assay
          }

          # Parse genes
          genes_raw <- unique(trimws(unlist(strsplit(input$heatmap_genes, "[,\n]+"))))
          genes_raw <- genes_raw[nzchar(genes_raw)]
          valid_genes <- genes_raw[genes_raw %in% rownames(obj)]

          if (length(valid_genes) == 0) {
            showNotification("No valid genes found in the dataset.", type = "error")
            return(NULL)
          }

          incProgress(0.2, detail = "Extracting expression...")

          # Get group column
          group_col <- if (!is.null(input$heatmap_group) && input$heatmap_group != "none" &&
                          input$heatmap_group %in% colnames(obj@meta.data)) {
            input$heatmap_group
          } else {
            NULL
          }

          # Get slot data -- handle scale.data carefully
          slot_name <- input$heatmap_slot
          expr_mat <- tryCatch({
            as.matrix(safe_get_assay_data(obj, slot_name)[valid_genes, , drop = FALSE])
          }, error = function(e) {
            # Fallback to data slot if scale.data not available
            showNotification("Requested slot not available, using 'data' slot.", type = "warning")
            as.matrix(safe_get_assay_data(obj, "data")[valid_genes, , drop = FALSE])
          })

          incProgress(0.4, detail = "Processing matrix...")

          if (input$average_expression && !is.null(group_col)) {
            # Compute average expression per group
            groups <- as.character(obj@meta.data[[group_col]])
            group_levels <- sort(unique(groups))
            avg_mat <- matrix(0, nrow = length(valid_genes), ncol = length(group_levels),
                              dimnames = list(valid_genes, group_levels))
            for (grp in group_levels) {
              cells <- colnames(obj)[groups == grp]
              if (length(cells) > 1) {
                avg_mat[, grp] <- rowMeans(expr_mat[, cells, drop = FALSE])
              } else if (length(cells) == 1) {
                avg_mat[, grp] <- expr_mat[, cells]
              }
            }
            expr_mat <- avg_mat
          } else {
            # Subsample cells for performance
            if (ncol(expr_mat) > input$hm_max_cells) {
              set.seed(42)
              sample_idx <- sample(seq_len(ncol(expr_mat)), input$hm_max_cells)
              expr_mat <- expr_mat[, sample_idx, drop = FALSE]
            }
          }

          # Apply scaling
          if (input$heatmap_scale == "row") {
            expr_mat <- t(scale(t(expr_mat)))
            # Clip extreme values
            expr_mat[expr_mat > 3] <- 3
            expr_mat[expr_mat < -3] <- -3
            expr_mat[is.na(expr_mat)] <- 0
          } else if (input$heatmap_scale == "column") {
            expr_mat <- scale(expr_mat)
            expr_mat[expr_mat > 3] <- 3
            expr_mat[expr_mat < -3] <- -3
            expr_mat[is.na(expr_mat)] <- 0
          }

          incProgress(0.6, detail = "Creating heatmap...")

          # Color function
          color_scheme <- input$heatmap_colors
          mat_range <- range(expr_mat, na.rm = TRUE)

          col_fun <- if (color_scheme == "bwr") {
            circlize::colorRamp2(
              c(mat_range[1], mean(mat_range), mat_range[2]),
              c("#2166AC", "#F7F7F7", "#B2182B"))
          } else if (color_scheme == "viridis") {
            circlize::colorRamp2(
              seq(mat_range[1], mat_range[2], length.out = 5),
              viridis::viridis(5))
          } else if (color_scheme == "PuOr") {
            circlize::colorRamp2(
              c(mat_range[1], mean(mat_range), mat_range[2]),
              c("#542788", "#F7F7F7", "#E08214"))
          } else if (color_scheme == "YlOrRd") {
            circlize::colorRamp2(
              seq(mat_range[1], mat_range[2], length.out = 5),
              RColorBrewer::brewer.pal(5, "YlOrRd"))
          } else if (color_scheme == "by") {
            circlize::colorRamp2(
              c(mat_range[1], mean(mat_range), mat_range[2]),
              c("#313695", "#FFFFBF", "#A50026"))
          } else {
            circlize::colorRamp2(
              c(mat_range[1], mean(mat_range), mat_range[2]),
              c("#2166AC", "#F7F7F7", "#B2182B"))
          }

          # Column annotation
          col_ha <- NULL
          if (!is.null(group_col) && !input$average_expression) {
            group_vals <- as.character(obj@meta.data[colnames(expr_mat), group_col])
            n_grps <- length(unique(group_vals))
            grp_colors <- get_enhanced_colors(n_grps, "default")
            names(grp_colors) <- sort(unique(group_vals))
            col_ha <- ComplexHeatmap::HeatmapAnnotation(
              Group = group_vals,
              col = list(Group = grp_colors),
              show_annotation_name = TRUE,
              annotation_name_side = "left"
            )
          }

          # Build the heatmap
          ht <- ComplexHeatmap::Heatmap(
            expr_mat,
            name = "Expression",
            col = col_fun,
            cluster_rows = input$cluster_rows,
            cluster_columns = input$cluster_cols,
            show_row_names = input$show_gene_names && nrow(expr_mat) <= 100,
            show_column_names = ncol(expr_mat) <= 50,
            row_names_gp = grid::gpar(fontsize = input$hm_fontsize),
            column_names_gp = grid::gpar(fontsize = max(6, input$hm_fontsize - 2)),
            top_annotation = col_ha,
            column_title = if (input$average_expression) "Groups" else paste0("Cells (n=", ncol(expr_mat), ")"),
            row_title = paste0("Genes (n=", nrow(expr_mat), ")"),
            heatmap_legend_param = list(title = "Expression", title_gp = grid::gpar(fontsize = input$hm_fontsize))
          )

          heatmap_obj_rv(ht)
          heatmap_matrix_rv(expr_mat)

          incProgress(1, detail = "Done!")
          showNotification("Heatmap built successfully.", type = "message")

        }, error = function(e) {
          showNotification(paste("Heatmap error:", e$message), type = "error")
        })
      })
    })

    # --- Render heatmap ------------------------------------------------------
    output$heatmap_output <- renderPlot({
      req(heatmap_obj_rv())
      ComplexHeatmap::draw(heatmap_obj_rv())
    })

    # Interactive plotly heatmap alternative (avoids ggplotly(ComplexHeatmap) issues)
    output$heatmap_interactive <- renderPlotly({
      req(heatmap_matrix_rv())
      mat <- heatmap_matrix_rv()
      # Limit to manageable size for interactivity
      if (nrow(mat) > 200) mat <- mat[seq_len(200), ]
      if (ncol(mat) > 500) mat <- mat[, seq_len(500)]
      plot_ly(z = as.matrix(mat), x = colnames(mat), y = rownames(mat),
              type = "heatmap", colorscale = "Viridis",
              hovertemplate = "Gene: %{y}<br>Cell: %{x}<br>Expr: %{z:.3f}<extra></extra>") %>%
        layout(title = "Interactive Heatmap",
               xaxis = list(title = "", showticklabels = ncol(mat) <= 50),
               yaxis = list(title = "", autorange = "reversed"))
    })

    # --- Download handlers ---------------------------------------------------
    output$dl_hm_png <- downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".png"),
      content = function(file) {
        tryCatch({
          ht <- heatmap_obj_rv()
          if (!is.null(ht)) {
            png(file, width = 10 * 300, height = 8 * 300, res = 300)
            ComplexHeatmap::draw(ht)
            dev.off()
          }
        }, error = function(e) {
          showNotification(paste("PNG export error:", e$message), type = "error")
        })
      }
    )

    output$dl_hm_pdf <- downloadHandler(
      filename = function() paste0("heatmap_", Sys.Date(), ".pdf"),
      content = function(file) {
        tryCatch({
          ht <- heatmap_obj_rv()
          if (!is.null(ht)) {
            pdf(file, width = 10, height = 8)
            ComplexHeatmap::draw(ht)
            dev.off()
          }
        }, error = function(e) {
          showNotification(paste("PDF export error:", e$message), type = "error")
        })
      }
    )

    output$dl_hm_csv <- downloadHandler(
      filename = function() paste0("heatmap_data_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          mat <- heatmap_matrix_rv()
          if (!is.null(mat)) {
            out <- as.data.frame(mat)
            out <- cbind(Gene = rownames(out), out)
            write.csv(out, file, row.names = FALSE)
          }
        }, error = function(e) {
          showNotification(paste("CSV export error:", e$message), type = "error")
        })
      }
    )

  })
}


# =============================================================================
# MODULE 3: CLUSTER COMPARISON
# =============================================================================

clusterComparisonUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Cluster Comparison"),
          p("Compare gene expression across clusters and batches with statistical testing."),

          fluidRow(
            column(4,
              selectInput(ns("cc_cluster_col"), "Cluster column", choices = NULL),
              selectInput(ns("cc_clusters"), "Select clusters to compare",
                          choices = NULL, multiple = TRUE)
            ),
            column(4,
              textAreaInput(ns("cc_genes"), "Genes to compare (comma-separated)",
                placeholder = "CD3D, CD4, CD8A, MS4A1, NKG7", rows = 2, width = "100%"),
              selectInput(ns("cc_batch_col"), "Batch / condition column (optional)",
                          choices = c("None" = "none"))
            ),
            column(4,
              selectInput(ns("cc_assay"), "Assay", choices = NULL),
              selectInput(ns("cc_reduction"), "Reduction (for UMAP)", choices = NULL),
              actionButton(ns("run_comparison"), "Run Comparison",
                class = "btn-primary btn-lg", width = "100%")
            )
          )
        )
      )
    ),

    # Output tabs
    fluidRow(
      column(12,
        tabsetPanel(
          tabPanel("Violin Comparison",
            uiOutput(ns("cc_violin_container")),
            downloadButton(ns("dl_cc_violin_html"), "HTML", class = "btn-info")
          ),
          tabPanel("Proportions",
            plotlyOutput(ns("cc_proportions"), height = "600px"),
            downloadButton(ns("dl_cc_prop_html"), "HTML", class = "btn-info")
          ),
          tabPanel("Statistical Tests",
            DT::dataTableOutput(ns("cc_stats_table")),
            downloadButton(ns("dl_cc_stats_csv"), "CSV", class = "btn-warning")
          ),
          tabPanel("Split UMAP",
            uiOutput(ns("cc_umap_container")),
            downloadButton(ns("dl_cc_umap_html"), "HTML", class = "btn-info")
          )
        )
      )
    )
  )
}

clusterComparisonServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    cc_violin_rv <- reactiveVal(NULL)
    cc_prop_rv <- reactiveVal(NULL)
    cc_stats_rv <- reactiveVal(NULL)
    cc_umap_rv <- reactiveVal(NULL)

    # --- Update UI choices ---------------------------------------------------
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset

      meta_cols <- colnames(obj@meta.data)
      factor_cols <- meta_cols[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      assay_names <- names(obj@assays)
      reds <- names(obj@reductions)

      updateSelectInput(session, "cc_cluster_col",
        choices = factor_cols,
        selected = if ("seurat_clusters" %in% factor_cols) "seurat_clusters" else factor_cols[1])
      updateSelectInput(session, "cc_batch_col",
        choices = c("None" = "none", factor_cols))
      updateSelectInput(session, "cc_assay",
        choices = assay_names, selected = safe_default_assay(obj))
      updateSelectInput(session, "cc_reduction",
        choices = reds, selected = best_reduction(obj))
    })

    # Update cluster choices when cluster column changes
    observeEvent(input$cc_cluster_col, {
      req(values$active_dataset, input$cc_cluster_col)
      if (input$cc_cluster_col %in% colnames(values$active_dataset@meta.data)) {
        vals <- sort(unique(as.character(values$active_dataset@meta.data[[input$cc_cluster_col]])))
        updateSelectInput(session, "cc_clusters", choices = vals, selected = vals)
      }
    })

    # --- Run comparison ------------------------------------------------------
    observeEvent(input$run_comparison, {
      req(values$active_dataset)

      withProgress(message = "Running cluster comparison...", value = 0, {
        tryCatch({
          obj <- values$active_dataset

          # Set assay
          if (!is.null(input$cc_assay) && input$cc_assay %in% names(obj@assays)) {
            DefaultAssay(obj) <- input$cc_assay
          }

          cluster_col <- input$cc_cluster_col
          selected_clusters <- input$cc_clusters

          if (is.null(cluster_col) || !cluster_col %in% colnames(obj@meta.data)) {
            showNotification("Please select a valid cluster column.", type = "error")
            return(NULL)
          }

          if (is.null(selected_clusters) || length(selected_clusters) < 2) {
            showNotification("Please select at least 2 clusters to compare.", type = "error")
            return(NULL)
          }

          # Subset to selected clusters
          keep_cells <- rownames(obj@meta.data)[
            as.character(obj@meta.data[[cluster_col]]) %in% selected_clusters
          ]
          obj_sub <- subset(obj, cells = keep_cells)

          # Parse genes
          genes_raw <- unique(trimws(unlist(strsplit(input$cc_genes, "[,\n]+"))))
          genes_raw <- genes_raw[nzchar(genes_raw)]
          valid_genes <- genes_raw[genes_raw %in% rownames(obj_sub)]

          batch_col <- if (!is.null(input$cc_batch_col) && input$cc_batch_col != "none" &&
                          input$cc_batch_col %in% colnames(obj_sub@meta.data)) {
            input$cc_batch_col
          } else {
            NULL
          }

          incProgress(0.2, detail = "Building violin plots...")

          # ---- 1. VIOLIN COMPARISON ----
          if (length(valid_genes) > 0) {
            expr_mat <- safe_get_assay_data(obj_sub, "data")[valid_genes, , drop = FALSE]
            clusters <- as.character(obj_sub@meta.data[[cluster_col]])
            cluster_levels <- sort(unique(clusters))
            colors <- get_enhanced_colors(length(cluster_levels), "default")

            if (!is.null(batch_col)) {
              # Side-by-side violins per cluster per batch
              batches <- as.character(obj_sub@meta.data[[batch_col]])
              batch_levels <- sort(unique(batches))

              plots <- list()
              for (gi in seq_along(valid_genes)) {
                g <- valid_genes[gi]
                g_expr <- as.numeric(expr_mat[g, ])
                plot_df <- data.frame(
                  Expression = g_expr, Cluster = clusters, Batch = batches,
                  stringsAsFactors = FALSE
                )
                plots[[gi]] <- plot_ly(data = plot_df, x = ~Cluster, y = ~Expression,
                                       color = ~Batch, type = "violin",
                                       box = list(visible = TRUE),
                                       showlegend = (gi == 1)) %>%
                  layout(annotations = list(text = g, x = 0.5, y = 1.05,
                                            xref = "paper", yref = "paper",
                                            showarrow = FALSE, font = list(size = 11)))
              }
              nc <- min(length(valid_genes), 3)
              nr <- ceiling(length(valid_genes) / nc)
              vp <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE) %>%
                layout(title = paste("Expression by Cluster and", batch_col))
            } else {
              plots <- list()
              for (gi in seq_along(valid_genes)) {
                g <- valid_genes[gi]
                g_expr <- as.numeric(expr_mat[g, ])
                plot_df <- data.frame(Expression = g_expr, Cluster = clusters,
                                      stringsAsFactors = FALSE)
                plots[[gi]] <- plot_ly(data = plot_df, x = ~Cluster, y = ~Expression,
                                       color = ~Cluster, colors = colors, type = "violin",
                                       box = list(visible = TRUE),
                                       showlegend = (gi == 1)) %>%
                  layout(annotations = list(text = g, x = 0.5, y = 1.05,
                                            xref = "paper", yref = "paper",
                                            showarrow = FALSE, font = list(size = 11)))
              }
              nc <- min(length(valid_genes), 3)
              nr <- ceiling(length(valid_genes) / nc)
              vp <- subplot(plots, nrows = nr, shareX = FALSE, shareY = FALSE) %>%
                layout(title = "Expression by Cluster")
            }
            cc_violin_rv(vp)
          } else {
            cc_violin_rv(plotly_empty() %>% layout(title = "No valid genes provided"))
          }

          incProgress(0.4, detail = "Building proportion chart...")

          # ---- 2. PROPORTIONAL BAR CHART ----
          clusters_all <- as.character(obj_sub@meta.data[[cluster_col]])

          if (!is.null(batch_col)) {
            batches_all <- as.character(obj_sub@meta.data[[batch_col]])
            prop_df <- as.data.frame(table(Batch = batches_all, Cluster = clusters_all))
            # Compute proportions within each batch
            prop_df <- prop_df %>%
              dplyr::group_by(Batch) %>%
              dplyr::mutate(Proportion = Freq / sum(Freq) * 100) %>%
              dplyr::ungroup()

            pp <- plot_ly(data = prop_df, x = ~Batch, y = ~Proportion,
                          color = ~Cluster, type = "bar",
                          colors = get_enhanced_colors(length(unique(prop_df$Cluster)), "default"),
                          hovertemplate = paste0("Batch: %{x}<br>Cluster: %{data.name}<br>",
                                                 "Proportion: %{y:.1f}%%<br>Count: %{customdata}<extra></extra>"),
                          customdata = ~Freq) %>%
              layout(barmode = "stack",
                     title = paste("Cell Type Proportions by", batch_col),
                     xaxis = list(title = batch_col),
                     yaxis = list(title = "Proportion (%)"))
          } else {
            count_df <- as.data.frame(table(Cluster = clusters_all))
            count_df$Proportion <- count_df$Freq / sum(count_df$Freq) * 100
            pp <- plot_ly(data = count_df, x = ~Cluster, y = ~Proportion,
                          color = ~Cluster,
                          colors = get_enhanced_colors(length(unique(count_df$Cluster)), "default"),
                          type = "bar",
                          hovertemplate = paste0("Cluster: %{x}<br>Proportion: %{y:.1f}%%<br>",
                                                 "Count: %{customdata}<extra></extra>"),
                          customdata = ~Freq) %>%
              layout(title = "Cell Type Proportions",
                     xaxis = list(title = cluster_col),
                     yaxis = list(title = "Proportion (%)"))
          }
          cc_prop_rv(pp)

          incProgress(0.6, detail = "Running statistical tests...")

          # ---- 3. STATISTICAL COMPARISON TABLE (Wilcoxon) ----
          if (length(valid_genes) > 0 && length(selected_clusters) >= 2) {
            stat_results <- list()
            cluster_pairs <- combn(selected_clusters, 2, simplify = FALSE)

            for (g in valid_genes) {
              g_expr <- as.numeric(expr_mat[g, ])
              for (pair in cluster_pairs) {
                vals1 <- g_expr[clusters == pair[1]]
                vals2 <- g_expr[clusters == pair[2]]

                test_result <- tryCatch({
                  wt <- wilcox.test(vals1, vals2, exact = FALSE)
                  list(p = wt$p.value, stat = wt$statistic)
                }, error = function(e) {
                  list(p = NA, stat = NA)
                })

                stat_results[[paste(g, pair[1], pair[2], sep = "___")]] <- data.frame(
                  Gene = g,
                  Cluster_1 = pair[1],
                  Cluster_2 = pair[2],
                  Mean_1 = round(mean(vals1, na.rm = TRUE), 4),
                  Mean_2 = round(mean(vals2, na.rm = TRUE), 4),
                  Log2FC = round(log2((mean(vals1, na.rm = TRUE) + 1) /
                                       (mean(vals2, na.rm = TRUE) + 1)), 4),
                  Wilcoxon_W = round(as.numeric(test_result$stat), 2),
                  P_value = signif(test_result$p, 4),
                  stringsAsFactors = FALSE
                )
              }
            }
            stats_df <- do.call(rbind, stat_results)
            rownames(stats_df) <- NULL

            # Adjust p-values
            if (nrow(stats_df) > 1) {
              stats_df$P_adjusted <- signif(p.adjust(stats_df$P_value, method = "BH"), 4)
            } else {
              stats_df$P_adjusted <- stats_df$P_value
            }

            stats_df$Significance <- ifelse(stats_df$P_adjusted < 0.001, "***",
                                     ifelse(stats_df$P_adjusted < 0.01, "**",
                                     ifelse(stats_df$P_adjusted < 0.05, "*", "ns")))

            cc_stats_rv(stats_df)
          } else {
            cc_stats_rv(data.frame(Message = "Provide genes and at least 2 clusters for statistical testing."))
          }

          incProgress(0.8, detail = "Building UMAP...")

          # ---- 4. SPLIT UMAP ----
          red_name <- input$cc_reduction
          if (!is.null(red_name) && red_name %in% names(obj_sub@reductions)) {
            emb <- Embeddings(obj_sub, red_name)[, 1:2]
            umap_df <- data.frame(
              Dim1 = emb[, 1], Dim2 = emb[, 2],
              Cluster = as.character(obj_sub@meta.data[[cluster_col]]),
              Cell = colnames(obj_sub),
              stringsAsFactors = FALSE
            )

            umap_colors <- get_enhanced_colors(length(unique(umap_df$Cluster)), "default")

            if (!is.null(batch_col)) {
              umap_df$Batch <- as.character(obj_sub@meta.data[[batch_col]])
              batch_levels <- sort(unique(umap_df$Batch))

              # Adaptive rendering: use scattergl for >50K cells to prevent browser crashes
              scatter_type <- if (nrow(umap_df) > 50000) "scattergl" else "scatter"
              marker_size <- if (nrow(umap_df) > 100000) 1 else if (nrow(umap_df) > 50000) 2 else 3
              marker_opacity <- if (nrow(umap_df) > 100000) 0.3 else 0.7

              umap_plots <- list()
              for (bi in seq_along(batch_levels)) {
                b <- batch_levels[bi]
                sub_df <- umap_df[umap_df$Batch == b, ]
                umap_plots[[bi]] <- plot_ly(data = sub_df, x = ~Dim1, y = ~Dim2,
                                            color = ~Cluster, colors = umap_colors,
                                            type = scatter_type, mode = "markers",
                                            marker = list(size = marker_size, opacity = marker_opacity),
                                            text = ~Cell, showlegend = (bi == 1),
                                            hovertemplate = paste0("Cell: %{text}<br>",
                                                                    "Cluster: %{data.name}<br>",
                                                                    toupper(red_name), "1: %{x:.2f}<br>",
                                                                    toupper(red_name), "2: %{y:.2f}<extra></extra>")) %>%
                  layout(annotations = list(text = b, x = 0.5, y = 1.05,
                                            xref = "paper", yref = "paper",
                                            showarrow = FALSE, font = list(size = 12)))
              }
              nc_u <- min(length(batch_levels), 3)
              nr_u <- ceiling(length(batch_levels) / nc_u)
              up <- subplot(umap_plots, nrows = nr_u, shareX = TRUE, shareY = TRUE) %>%
                layout(title = paste(toupper(red_name), "split by", batch_col))
            } else {
              up <- plot_ly(data = umap_df, x = ~Dim1, y = ~Dim2,
                            color = ~Cluster, colors = umap_colors,
                            type = scatter_type, mode = "markers",
                            marker = list(size = marker_size, opacity = marker_opacity),
                            text = ~Cell,
                            hovertemplate = paste0("Cell: %{text}<br>",
                                                    "Cluster: %{data.name}<br>",
                                                    toupper(red_name), "1: %{x:.2f}<br>",
                                                    toupper(red_name), "2: %{y:.2f}<extra></extra>")) %>%
                layout(title = paste(toupper(red_name), "- Cluster Comparison"),
                       xaxis = list(title = paste0(toupper(red_name), "_1")),
                       yaxis = list(title = paste0(toupper(red_name), "_2")))
            }
            cc_umap_rv(up)
          } else {
            cc_umap_rv(plotly_empty() %>% layout(title = "No dimensional reduction available"))
          }

          incProgress(1, detail = "Done!")
          showNotification("Cluster comparison complete.", type = "message")

        }, error = function(e) {
          showNotification(paste("Comparison error:", e$message), type = "error")
        })
      })
    })

    # --- Render outputs with dynamic heights ---------------------------------

    # Dynamic violin container: height scales with number of genes
    output$cc_violin_container <- renderUI({
      n_genes <- length(input$cc_genes)
      n_rows <- ceiling(max(1, n_genes) / 3)
      h <- max(600, n_rows * 350)
      plotlyOutput(session$ns("cc_violin"), height = paste0(h, "px"))
    })

    output$cc_violin <- renderPlotly({
      req(cc_violin_rv())
      tryCatch({
        p <- cc_violin_rv() %>% layout(autosize = TRUE) %>% plotly::partial_bundle()
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Visualization", figure_type = "Cluster Comparison Violin",
                         label = "Cluster Comparison Violin", source_module = "Visualization")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })

    output$cc_proportions <- renderPlotly({
      req(cc_prop_rv())
      tryCatch({
        p <- cc_prop_rv() %>% layout(autosize = TRUE) %>% plotly::partial_bundle()
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Visualization", figure_type = "Cluster Proportions",
                         label = "Cluster Proportions", source_module = "Visualization")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })

    output$cc_stats_table <- DT::renderDT({
      req(cc_stats_rv())
      DT::datatable(cc_stats_rv(), options = list(
        pageLength = 25, scrollX = TRUE, dom = "frtip"
      ))
    })

    # Dynamic UMAP container: height scales with batch count
    output$cc_umap_container <- renderUI({
      n_batches <- 1
      if (!is.null(input$cc_batch) && nzchar(input$cc_batch) && !is.null(values$active_dataset)) {
        if (input$cc_batch %in% colnames(values$active_dataset@meta.data)) {
          n_batches <- length(unique(values$active_dataset@meta.data[[input$cc_batch]]))
        }
      }
      n_rows <- ceiling(max(1, n_batches) / 3)
      h <- max(600, n_rows * 350)
      plotlyOutput(session$ns("cc_umap"), height = paste0(h, "px"))
    })

    output$cc_umap <- renderPlotly({
      req(cc_umap_rv())
      tryCatch({
        p <- cc_umap_rv() %>% layout(autosize = TRUE) %>% plotly::partial_bundle()
        if (exists("capture_umap_snapshot", mode = "function"))
          capture_umap_snapshot(values, p, "Cluster Comparison", coloring = "clusters", method = "umap")
        if (exists("capture_figure", mode = "function")) {
          capture_figure(values, p, "Visualization", figure_type = "Cluster Comparison UMAP",
                         label = "Cluster Comparison UMAP", source_module = "Visualization")
        }
        p
      }, error = function(e) {
        plotly_empty() %>% layout(title = paste("Error:", e$message))
      })
    })

    # --- Download handlers ---------------------------------------------------
    output$dl_cc_violin_html <- downloadHandler(
      filename = function() paste0("cluster_comparison_violin_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- cc_violin_rv()
          if (!is.null(p)) htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
        }, error = function(e) {
          showNotification(paste("Export error:", e$message), type = "error")
        })
      }
    )

    output$dl_cc_prop_html <- downloadHandler(
      filename = function() paste0("cluster_proportions_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- cc_prop_rv()
          if (!is.null(p)) htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
        }, error = function(e) {
          showNotification(paste("Export error:", e$message), type = "error")
        })
      }
    )

    output$dl_cc_stats_csv <- downloadHandler(
      filename = function() paste0("cluster_comparison_stats_", Sys.Date(), ".csv"),
      content = function(file) {
        tryCatch({
          df <- cc_stats_rv()
          if (!is.null(df)) write.csv(df, file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste("Export error:", e$message), type = "error")
        })
      }
    )

    output$dl_cc_umap_html <- downloadHandler(
      filename = function() paste0("cluster_comparison_umap_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- cc_umap_rv()
          if (!is.null(p)) htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
        }, error = function(e) {
          showNotification(paste("Export error:", e$message), type = "error")
        })
      }
    )

  })
}


# =============================================================================
# MODULE 4: PUBLICATION FIGURE PANEL
# =============================================================================

figurePanelUI <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Publication Figure Panel"),
          p("Create multi-panel publication-ready figures with configurable layout."),

          fluidRow(
            column(3,
              selectInput(ns("fp_layout"), "Panel Layout",
                choices = c("1 x 2" = "1x2", "2 x 1" = "2x1",
                            "2 x 2" = "2x2", "2 x 3" = "2x3",
                            "3 x 2" = "3x2", "1 x 3" = "1x3"),
                selected = "2x2")
            ),
            column(3,
              numericInput(ns("fp_width"), "Figure Width (in)", 12, min = 6, max = 24),
              numericInput(ns("fp_height"), "Figure Height (in)", 10, min = 4, max = 20)
            ),
            column(3,
              numericInput(ns("fp_dpi"), "DPI", 300, min = 72, max = 600),
              numericInput(ns("fp_fontsize"), "Base font size", 12, min = 8, max = 20)
            ),
            column(3,
              selectInput(ns("fp_reduction"), "Reduction", choices = NULL),
              selectInput(ns("fp_group_col"), "Group By", choices = NULL)
            )
          ),

          hr(),
          h5("Panel Configuration"),

          # Panel A
          fluidRow(
            column(6,
              wellPanel(style = "background-color:#f8f9fa;",
                h6("Panel A"),
                selectInput(ns("panel_a_type"), "Plot Type",
                  choices = c("UMAP" = "umap", "Feature Plot" = "feature",
                              "Violin" = "violin", "Dot Plot" = "dotplot",
                              "Bar Chart" = "bar", "Heatmap" = "heatmap"),
                  selected = "umap"),
                textInput(ns("panel_a_genes"), "Genes (if applicable)", placeholder = "CD3D"),
                textInput(ns("panel_a_title"), "Title", "A")
              )
            ),
            column(6,
              wellPanel(style = "background-color:#f8f9fa;",
                h6("Panel B"),
                selectInput(ns("panel_b_type"), "Plot Type",
                  choices = c("UMAP" = "umap", "Feature Plot" = "feature",
                              "Violin" = "violin", "Dot Plot" = "dotplot",
                              "Bar Chart" = "bar", "Heatmap" = "heatmap"),
                  selected = "feature"),
                textInput(ns("panel_b_genes"), "Genes (if applicable)", placeholder = "CD4, CD8A"),
                textInput(ns("panel_b_title"), "Title", "B")
              )
            )
          ),

          # Panel C and D
          fluidRow(
            column(6,
              wellPanel(style = "background-color:#f8f9fa;",
                h6("Panel C"),
                selectInput(ns("panel_c_type"), "Plot Type",
                  choices = c("None" = "none", "UMAP" = "umap", "Feature Plot" = "feature",
                              "Violin" = "violin", "Dot Plot" = "dotplot",
                              "Bar Chart" = "bar", "Heatmap" = "heatmap"),
                  selected = "violin"),
                textInput(ns("panel_c_genes"), "Genes (if applicable)", placeholder = "MS4A1, CD14"),
                textInput(ns("panel_c_title"), "Title", "C")
              )
            ),
            column(6,
              wellPanel(style = "background-color:#f8f9fa;",
                h6("Panel D"),
                selectInput(ns("panel_d_type"), "Plot Type",
                  choices = c("None" = "none", "UMAP" = "umap", "Feature Plot" = "feature",
                              "Violin" = "violin", "Dot Plot" = "dotplot",
                              "Bar Chart" = "bar", "Heatmap" = "heatmap"),
                  selected = "dotplot"),
                textInput(ns("panel_d_genes"), "Genes (if applicable)", placeholder = "CD3D, CD4, CD8A, MS4A1"),
                textInput(ns("panel_d_title"), "Title", "D")
              )
            )
          ),

          # Panel E and F
          fluidRow(
            column(6,
              wellPanel(style = "background-color:#f8f9fa;",
                h6("Panel E"),
                selectInput(ns("panel_e_type"), "Plot Type",
                  choices = c("None" = "none", "UMAP" = "umap", "Feature Plot" = "feature",
                              "Violin" = "violin", "Dot Plot" = "dotplot",
                              "Bar Chart" = "bar", "Heatmap" = "heatmap"),
                  selected = "none"),
                textInput(ns("panel_e_genes"), "Genes (if applicable)", placeholder = ""),
                textInput(ns("panel_e_title"), "Title", "E")
              )
            ),
            column(6,
              wellPanel(style = "background-color:#f8f9fa;",
                h6("Panel F"),
                selectInput(ns("panel_f_type"), "Plot Type",
                  choices = c("None" = "none", "UMAP" = "umap", "Feature Plot" = "feature",
                              "Violin" = "violin", "Dot Plot" = "dotplot",
                              "Bar Chart" = "bar", "Heatmap" = "heatmap"),
                  selected = "none"),
                textInput(ns("panel_f_genes"), "Genes (if applicable)", placeholder = ""),
                textInput(ns("panel_f_title"), "Title", "F")
              )
            )
          ),

          actionButton(ns("build_figure"), "Build Figure Panel",
            class = "btn-success btn-lg", width = "100%")
        )
      )
    ),

    fluidRow(
      column(12,
        plotOutput(ns("figure_output"), height = "800px"),
        br(),
        fluidRow(
          column(4, downloadButton(ns("dl_fig_png"), "PNG", class = "btn-success")),
          column(4, downloadButton(ns("dl_fig_pdf"), "PDF", class = "btn-info")),
          column(4, downloadButton(ns("dl_fig_tiff"), "TIFF", class = "btn-warning"))
        )
      )
    )
  )
}

figurePanelServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    combined_figure_rv <- reactiveVal(NULL)

    # --- Update choices ------------------------------------------------------
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset

      meta_cols <- colnames(obj@meta.data)
      factor_cols <- meta_cols[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      reds <- names(obj@reductions)

      updateSelectInput(session, "fp_reduction",
        choices = reds, selected = best_reduction(obj))
      updateSelectInput(session, "fp_group_col",
        choices = c("None" = "none", factor_cols),
        selected = if ("seurat_clusters" %in% factor_cols) "seurat_clusters" else "none")
    })

    # --- Helper: generate a single ggplot panel ------------------------------
    build_single_panel <- function(obj, panel_type, genes_str, title_label,
                                   reduction, group_col, base_size) {
      if (panel_type == "none") return(NULL)

      tryCatch({
        genes <- unique(trimws(unlist(strsplit(genes_str, "[,\n]+"))))
        genes <- genes[nzchar(genes)]
        valid_genes <- genes[genes %in% rownames(obj)]

        group_by <- if (!is.null(group_col) && group_col != "none" &&
                        group_col %in% colnames(obj@meta.data)) group_col else NULL

        theme_pub <- theme_minimal(base_size = base_size) +
          theme(plot.title = element_text(face = "bold", size = base_size + 2),
                legend.position = "right")

        # ---- UMAP ----
        if (panel_type == "umap") {
          if (is.null(reduction) || !reduction %in% names(obj@reductions)) {
            reduction <- best_reduction(obj)
          }
          if (is.null(reduction)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No reduction") + theme_void())

          emb <- as.data.frame(Embeddings(obj, reduction)[, 1:2])
          colnames(emb) <- c("Dim1", "Dim2")
          if (!is.null(group_by)) {
            emb$Group <- obj@meta.data[[group_by]]
          }

          p <- ggplot(emb, aes(x = Dim1, y = Dim2)) +
            labs(title = title_label,
                 x = paste0(toupper(reduction), "_1"),
                 y = paste0(toupper(reduction), "_2")) +
            theme_pub

          if (!is.null(group_by)) {
            n_grps <- length(unique(emb$Group))
            p <- p + geom_point(aes(color = Group), size = 0.5, alpha = 0.7) +
              scale_color_manual(values = get_enhanced_colors(n_grps, "default"))
          } else {
            p <- p + geom_point(size = 0.5, alpha = 0.7, color = "steelblue")
          }
          return(p)
        }

        # ---- FEATURE PLOT ----
        if (panel_type == "feature") {
          if (length(valid_genes) == 0) {
            return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid genes") + theme_void() +
                     labs(title = title_label))
          }
          feat <- valid_genes[1]
          if (is.null(reduction) || !reduction %in% names(obj@reductions)) {
            reduction <- best_reduction(obj)
          }
          if (is.null(reduction)) return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No reduction") + theme_void())

          emb <- as.data.frame(Embeddings(obj, reduction)[, 1:2])
          colnames(emb) <- c("Dim1", "Dim2")
          emb$Expression <- as.numeric(safe_get_assay_data(obj, "data")[feat, ])

          p <- ggplot(emb, aes(x = Dim1, y = Dim2, color = Expression)) +
            geom_point(size = 0.5, alpha = 0.7) +
            scale_color_viridis_c(option = "C") +
            labs(title = title_label, subtitle = feat,
                 x = paste0(toupper(reduction), "_1"),
                 y = paste0(toupper(reduction), "_2")) +
            theme_pub
          return(p)
        }

        # ---- VIOLIN ----
        if (panel_type == "violin") {
          if (length(valid_genes) == 0) {
            return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid genes") + theme_void() +
                     labs(title = title_label))
          }
          feat <- valid_genes[1]
          expr_vec <- as.numeric(safe_get_assay_data(obj, "data")[feat, ])
          grp <- if (!is.null(group_by)) as.character(obj@meta.data[[group_by]]) else rep("All", ncol(obj))

          df <- data.frame(Expression = expr_vec, Group = grp, stringsAsFactors = FALSE)
          n_grps <- length(unique(df$Group))

          p <- ggplot(df, aes(x = Group, y = Expression, fill = Group)) +
            geom_violin(scale = "width", trim = TRUE) +
            geom_boxplot(width = 0.1, fill = "white", alpha = 0.5) +
            scale_fill_manual(values = get_enhanced_colors(n_grps, "default")) +
            labs(title = title_label, subtitle = feat, x = "", y = "Expression") +
            theme_pub +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
          return(p)
        }

        # ---- DOT PLOT ----
        if (panel_type == "dotplot") {
          if (length(valid_genes) == 0) {
            return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid genes") + theme_void() +
                     labs(title = title_label))
          }
          grp_vec <- if (!is.null(group_by)) as.character(obj@meta.data[[group_by]]) else rep("All", ncol(obj))
          grp_levels <- sort(unique(grp_vec))
          expr_mat <- safe_get_assay_data(obj, "data")[valid_genes, , drop = FALSE]

          dot_list <- list()
          for (g in valid_genes) {
            for (grp in grp_levels) {
              cells <- which(grp_vec == grp)
              vals <- as.numeric(expr_mat[g, cells])
              dot_list[[paste(g, grp)]] <- data.frame(
                Gene = g, Group = grp,
                Avg_Expr = mean(vals),
                Pct_Expr = sum(vals > 0) / length(vals) * 100,
                stringsAsFactors = FALSE
              )
            }
          }
          dot_df <- do.call(rbind, dot_list)

          p <- ggplot(dot_df, aes(x = Gene, y = Group, size = Pct_Expr, color = Avg_Expr)) +
            geom_point() +
            scale_color_viridis_c(option = "C") +
            scale_size_continuous(range = c(1, 8)) +
            labs(title = title_label, x = "", y = "", color = "Avg Expr", size = "% Expr") +
            theme_pub +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          return(p)
        }

        # ---- BAR CHART ----
        if (panel_type == "bar") {
          grp_vec <- if (!is.null(group_by)) as.character(obj@meta.data[[group_by]]) else rep("All", ncol(obj))
          count_df <- as.data.frame(table(Group = grp_vec))
          count_df$Pct <- count_df$Freq / sum(count_df$Freq) * 100
          n_grps <- length(unique(count_df$Group))

          p <- ggplot(count_df, aes(x = Group, y = Pct, fill = Group)) +
            geom_col() +
            scale_fill_manual(values = get_enhanced_colors(n_grps, "default")) +
            labs(title = title_label, x = "", y = "Percentage (%)") +
            theme_pub +
            theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
          return(p)
        }

        # ---- HEATMAP (as ggplot tile) ----
        if (panel_type == "heatmap") {
          if (length(valid_genes) == 0) {
            return(ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No valid genes") + theme_void() +
                     labs(title = title_label))
          }
          grp_vec <- if (!is.null(group_by)) as.character(obj@meta.data[[group_by]]) else rep("All", ncol(obj))
          grp_levels <- sort(unique(grp_vec))
          expr_mat <- as.matrix(safe_get_assay_data(obj, "data")[valid_genes, , drop = FALSE])

          # Average per group
          avg_mat <- matrix(0, nrow = length(valid_genes), ncol = length(grp_levels),
                            dimnames = list(valid_genes, grp_levels))
          for (grp in grp_levels) {
            cells <- which(grp_vec == grp)
            if (length(cells) > 1) avg_mat[, grp] <- rowMeans(expr_mat[, cells, drop = FALSE])
            else if (length(cells) == 1) avg_mat[, grp] <- expr_mat[, cells]
          }
          # Scale rows
          avg_mat <- t(scale(t(avg_mat)))
          avg_mat[is.na(avg_mat)] <- 0

          hm_df <- as.data.frame(as.table(avg_mat))
          colnames(hm_df) <- c("Gene", "Group", "Expression")

          p <- ggplot(hm_df, aes(x = Group, y = Gene, fill = Expression)) +
            geom_tile(color = "white", linewidth = 0.3) +
            scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
            labs(title = title_label, x = "", y = "", fill = "Z-score") +
            theme_pub +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
          return(p)
        }

        # Fallback
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Unknown panel type") + theme_void()

      }, error = function(e) {
        ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Error:", e$message)) +
          theme_void() + labs(title = title_label)
      })
    }

    # --- Build figure --------------------------------------------------------
    observeEvent(input$build_figure, {
      req(values$active_dataset)

      withProgress(message = "Building figure panel...", value = 0, {
        tryCatch({
          obj <- values$active_dataset
          reduction <- input$fp_reduction
          group_col <- input$fp_group_col
          base_size <- input$fp_fontsize

          # Parse layout
          layout_str <- input$fp_layout
          layout_parts <- as.integer(strsplit(layout_str, "x")[[1]])
          n_rows <- layout_parts[1]
          n_cols <- layout_parts[2]
          n_panels <- n_rows * n_cols

          # Collect panel configs
          panel_configs <- list(
            list(type = input$panel_a_type, genes = input$panel_a_genes, title = input$panel_a_title),
            list(type = input$panel_b_type, genes = input$panel_b_genes, title = input$panel_b_title),
            list(type = input$panel_c_type, genes = input$panel_c_genes, title = input$panel_c_title),
            list(type = input$panel_d_type, genes = input$panel_d_genes, title = input$panel_d_title),
            list(type = input$panel_e_type, genes = input$panel_e_genes, title = input$panel_e_title),
            list(type = input$panel_f_type, genes = input$panel_f_genes, title = input$panel_f_title)
          )

          # Take only the needed panels
          panel_configs <- panel_configs[seq_len(min(n_panels, 6))]

          # Build individual panels
          panel_plots <- list()
          active_count <- 0

          for (i in seq_along(panel_configs)) {
            cfg <- panel_configs[[i]]
            incProgress(i / (length(panel_configs) + 1), detail = paste("Building panel", cfg$title))

            if (cfg$type != "none") {
              p <- build_single_panel(obj, cfg$type, cfg$genes, cfg$title,
                                      reduction, group_col, base_size)
              if (!is.null(p)) {
                active_count <- active_count + 1
                panel_plots[[active_count]] <- p
              }
            }
          }

          if (active_count == 0) {
            showNotification("No active panels configured.", type = "error")
            return(NULL)
          }

          # Combine with patchwork
          combined <- patchwork::wrap_plots(panel_plots, ncol = n_cols, nrow = n_rows) +
            patchwork::plot_annotation(
              tag_levels = "A",
              theme = theme(plot.tag = element_text(size = base_size + 4, face = "bold"))
            )

          combined_figure_rv(combined)

          incProgress(1, detail = "Done!")
          showNotification("Figure panel built successfully.", type = "message")

        }, error = function(e) {
          showNotification(paste("Figure panel error:", e$message), type = "error")
        })
      })
    })

    # --- Render figure -------------------------------------------------------
    output$figure_output <- renderPlot({
      req(combined_figure_rv())
      combined_figure_rv()
    }, width = function() input$fp_width * 72,
       height = function() input$fp_height * 72)

    # --- Download handlers ---------------------------------------------------
    output$dl_fig_png <- downloadHandler(
      filename = function() paste0("figure_panel_", Sys.Date(), ".png"),
      content = function(file) {
        tryCatch({
          fig <- combined_figure_rv()
          if (!is.null(fig)) {
            ggsave(file, plot = fig, width = input$fp_width, height = input$fp_height,
                   dpi = input$fp_dpi, units = "in", device = "png")
          }
        }, error = function(e) {
          showNotification(paste("PNG export error:", e$message), type = "error")
        })
      }
    )

    output$dl_fig_pdf <- downloadHandler(
      filename = function() paste0("figure_panel_", Sys.Date(), ".pdf"),
      content = function(file) {
        tryCatch({
          fig <- combined_figure_rv()
          if (!is.null(fig)) {
            ggsave(file, plot = fig, width = input$fp_width, height = input$fp_height,
                   units = "in", device = "pdf")
          }
        }, error = function(e) {
          showNotification(paste("PDF export error:", e$message), type = "error")
        })
      }
    )

    output$dl_fig_tiff <- downloadHandler(
      filename = function() paste0("figure_panel_", Sys.Date(), ".tiff"),
      content = function(file) {
        tryCatch({
          fig <- combined_figure_rv()
          if (!is.null(fig)) {
            ggsave(file, plot = fig, width = input$fp_width, height = input$fp_height,
                   dpi = input$fp_dpi, units = "in", device = "tiff", compression = "lzw")
          }
        }, error = function(e) {
          showNotification(paste("TIFF export error:", e$message), type = "error")
        })
      }
    )

  })
}

# =============================================================================
# COMPARISON VIEW MODULE (Split-Screen Side-by-Side)
# =============================================================================

comparisonViewUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4(icon("columns"), "Split-Screen Comparison"),
          fluidRow(
            column(3, selectInput(ns("comp_meta_col"), "Group by metadata",
              choices = NULL)),
            column(3, selectInput(ns("comp_reduction"), "Reduction",
              choices = NULL)),
            column(3, textInput(ns("comp_gene"), "Feature (gene)", placeholder = "e.g. CD3D")),
            column(3, radioGroupButtons(ns("comp_mode"), "Mode",
              choices = c("Side by Side" = "side", "Overlay" = "overlay"),
              selected = "side", justified = TRUE, size = "sm"))
          ),
          fluidRow(
            column(6, selectInput(ns("comp_left_val"), "Left panel group", choices = NULL)),
            column(6, selectInput(ns("comp_right_val"), "Right panel group", choices = NULL))
          ),
          actionButton(ns("run_compare"), "Compare", class = "btn-primary", icon = icon("play"))
        )
      )
    ),
    fluidRow(
      column(6,
        box(title = "Left Panel", status = "info", solidHeader = TRUE, width = 12,
          plotlyOutput(ns("comp_plot_left"), height = "500px"),
          verbatimTextOutput(ns("comp_info_left"))
        )
      ),
      column(6,
        box(title = "Right Panel", status = "warning", solidHeader = TRUE, width = 12,
          plotlyOutput(ns("comp_plot_right"), height = "500px"),
          verbatimTextOutput(ns("comp_info_right"))
        )
      )
    ),
    fluidRow(
      column(12,
        downloadButton(ns("dl_comp_left"), "Download Left"),
        downloadButton(ns("dl_comp_right"), "Download Right")
      )
    )
  )
}

comparisonViewServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {

    # Update metadata column choices
    observe({
      req(values$active_dataset)
      obj <- values$active_dataset
      meta_cols <- colnames(obj@meta.data)
      factor_cols <- meta_cols[sapply(obj@meta.data, function(x) is.factor(x) || is.character(x))]
      updateSelectInput(session, "comp_meta_col",
        choices = c(factor_cols, meta_cols),
        selected = if ("seurat_clusters" %in% meta_cols) "seurat_clusters" else factor_cols[1])
    })

    # Update reduction choices
    observe({
      req(values$active_dataset)
      reds <- names(values$active_dataset@reductions)
      updateSelectInput(session, "comp_reduction", choices = reds,
        selected = best_reduction(values$active_dataset))
    })

    # Update group value choices when metadata column changes
    observe({
      req(values$active_dataset, input$comp_meta_col)
      obj <- values$active_dataset
      if (input$comp_meta_col %in% colnames(obj@meta.data)) {
        vals <- sort(unique(as.character(obj@meta.data[[input$comp_meta_col]])))
        updateSelectInput(session, "comp_left_val", choices = vals,
          selected = vals[1])
        updateSelectInput(session, "comp_right_val", choices = vals,
          selected = if (length(vals) >= 2) vals[2] else vals[1])
      }
    })

    # Reactive: build a plotly for a subset of cells
    build_comparison_plot <- function(obj, cells, red, gene, label) {
      if (length(cells) == 0) {
        return(plotly_empty() %>% layout(title = paste(label, "- No cells")))
      }
      sub_obj <- obj[, cells]
      embed <- tryCatch(Embeddings(sub_obj, red), error = function(e) NULL)
      if (is.null(embed)) {
        return(plotly_empty() %>% layout(title = paste(label, "- Reduction not found")))
      }
      df <- as.data.frame(embed[, 1:min(2, ncol(embed)), drop = FALSE])
      colnames(df) <- c("Dim1", "Dim2")[1:ncol(df)]
      if (ncol(df) < 2) df$Dim2 <- 0

      if (nzchar(gene) && gene %in% rownames(sub_obj)) {
        expr <- tryCatch(
          as.numeric(safe_get_assay_data(sub_obj, "data")[gene, ]),
          error = function(e) rep(0, ncol(sub_obj))
        )
        df$Expression <- expr
        p <- plot_ly(df, x = ~Dim1, y = ~Dim2, color = ~Expression,
                     colors = viridis::viridis(50),
                     type = "scatter", mode = "markers",
                     marker = list(size = 4, opacity = 0.7),
                     text = ~paste("Expr:", round(Expression, 2))) %>%
          layout(title = paste(label, "-", gene))
      } else {
        clusters <- as.character(sub_obj@meta.data[["seurat_clusters"]])
        if (is.null(clusters)) clusters <- rep("NA", ncol(sub_obj))
        df$Cluster <- clusters
        pal <- get_color_palette(length(unique(df$Cluster)))
        p <- plot_ly(df, x = ~Dim1, y = ~Dim2, color = ~Cluster,
                     colors = pal,
                     type = "scatter", mode = "markers",
                     marker = list(size = 4, opacity = 0.7)) %>%
          layout(title = label)
      }
      p %>% layout(
        xaxis = list(title = paste0(toupper(red), "_1")),
        yaxis = list(title = paste0(toupper(red), "_2"))
      ) %>% plotly::partial_bundle()
    }

    # Store comparison results
    comp_left <- reactiveVal(NULL)
    comp_right <- reactiveVal(NULL)
    comp_info_l <- reactiveVal("")
    comp_info_r <- reactiveVal("")

    observeEvent(input$run_compare, {
      req(values$active_dataset, input$comp_meta_col, input$comp_left_val, input$comp_right_val, input$comp_reduction)
      obj <- values$active_dataset
      meta_col <- input$comp_meta_col
      red <- input$comp_reduction
      gene <- trimws(input$comp_gene)

      left_cells <- colnames(obj)[as.character(obj@meta.data[[meta_col]]) == input$comp_left_val]
      right_cells <- colnames(obj)[as.character(obj@meta.data[[meta_col]]) == input$comp_right_val]

      withProgress(message = "Building comparison...", value = 0, {
        incProgress(0.3, detail = "Left panel...")
        comp_left(build_comparison_plot(obj, left_cells, red, gene, input$comp_left_val))
        comp_info_l(paste0(input$comp_left_val, ": ", length(left_cells), " cells"))

        incProgress(0.3, detail = "Right panel...")
        comp_right(build_comparison_plot(obj, right_cells, red, gene, input$comp_right_val))
        comp_info_r(paste0(input$comp_right_val, ": ", length(right_cells), " cells"))

        incProgress(0.4, detail = "Done")
      })
    })

    output$comp_plot_left <- renderPlotly({
      req(comp_left())
      p <- comp_left()
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Visualization", figure_type = "Comparison Plot",
                       label = "Comparison Plot Left", source_module = "Visualization")
      }
      p
    })

    output$comp_plot_right <- renderPlotly({
      req(comp_right())
      p <- comp_right()
      if (exists("capture_figure", mode = "function")) {
        capture_figure(values, p, "Visualization", figure_type = "Comparison Plot",
                       label = "Comparison Plot Right", source_module = "Visualization")
      }
      p
    })

    output$comp_info_left <- renderText({ comp_info_l() })
    output$comp_info_right <- renderText({ comp_info_r() })

    output$dl_comp_left <- downloadHandler(
      filename = function() paste0("comparison_left_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- comp_left()
          if (!is.null(p)) htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
        }, error = function(e) {
          showNotification(paste("Export error:", e$message), type = "error")
        })
      }
    )

    output$dl_comp_right <- downloadHandler(
      filename = function() paste0("comparison_right_", Sys.Date(), ".html"),
      content = function(file) {
        tryCatch({
          p <- comp_right()
          if (!is.null(p)) htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
        }, error = function(e) {
          showNotification(paste("Export error:", e$message), type = "error")
        })
      }
    )
  })
}