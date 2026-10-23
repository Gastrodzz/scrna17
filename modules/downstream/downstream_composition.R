# ──────────────────────────────────────────────────────────────────────────────
# downstream_composition.R — Cell Composition / Abundance Testing Module
# Tests whether cell type proportions differ between conditions using propeller
# ──────────────────────────────────────────────────────────────────────────────

# =============================================================================
# UI
# =============================================================================

compositionUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Cell Type Composition Analysis"),
          p("Test whether cell type proportions differ between conditions ",
            "using the propeller method (speckle package) or chi-square test."),

          fluidRow(
            column(3,
              selectInput(ns("celltype_col"), "Cell Type Column",
                choices = NULL,
                selected = NULL)
            ),
            column(3,
              selectInput(ns("condition_col"), "Condition/Group Column",
                choices = NULL,
                selected = NULL)
            ),
            column(3,
              selectInput(ns("sample_col"), "Sample/Replicate Column",
                choices = NULL,
                selected = NULL),
              helpText("Required for propeller. Each sample = one biological replicate.",
                       style = "font-size: 11px; color: #888;")
            ),
            column(3,
              div(style = "margin-top: 25px;",
                actionButton(ns("run_composition"), "Run Composition Test",
                  class = "btn-primary", width = "100%", icon = icon("chart-bar"))
              )
            )
          ),

          fluidRow(
            column(12,
              div(style = "margin-top: 8px;",
                actionButton(ns("ai_interpret_composition"),
                  "AI Interpret Composition",
                  icon = icon("robot"), class = "ai-btn btn-sm")
              ),
              uiOutput(ns("ai_thread_composition"))
            )
          )
        )
      )
    ),

    fluidRow(
      column(6,
        wellPanel(
          h5("Proportion Bar Chart"),
          shinycssloaders::withSpinner(
            plotly::plotlyOutput(ns("composition_barplot"), height = "500px")
          )
        )
      ),
      column(6,
        wellPanel(
          h5("Composition Test Results"),
          shinycssloaders::withSpinner(
            DT::dataTableOutput(ns("composition_table"))
          )
        )
      )
    )
  )
}

# =============================================================================
# SERVER
# =============================================================================

compositionServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive storage for results
    composition_results <- reactiveVal(NULL)

    # Update column choices when dataset changes
    observe({
      req(values$active_dataset)
      md <- values$active_dataset@meta.data
      cat_cols <- colnames(md)[sapply(md, function(x) is.factor(x) || is.character(x))]

      # Detect likely cell type columns
      ct_cols <- grep("cell_type|annotation|label|celltype", cat_cols,
                      value = TRUE, ignore.case = TRUE)
      ct_choices <- if (length(ct_cols) > 0) ct_cols else cat_cols

      updateSelectInput(session, "celltype_col", choices = ct_choices,
                        selected = if (length(ct_cols) > 0) ct_cols[1] else NULL)
      updateSelectInput(session, "condition_col", choices = cat_cols)
      updateSelectInput(session, "sample_col", choices = cat_cols)
    })

    # Run composition test
    observeEvent(input$run_composition, {
      req(values$active_dataset, input$celltype_col, input$condition_col, input$sample_col)

      obj <- values$active_dataset
      md <- obj@meta.data

      # Validate columns exist
      for (col in c(input$celltype_col, input$condition_col, input$sample_col)) {
        if (!col %in% colnames(md)) {
          showNotification(paste0("Column '", col, "' not found in metadata."),
                           type = "error")
          return()
        }
      }

      withProgress(message = "Running composition analysis...", value = 0.3, {

        ct <- as.character(md[[input$celltype_col]])
        cond <- as.character(md[[input$condition_col]])
        sample_id <- as.character(md[[input$sample_col]])

        # Compute proportions per sample
        prop_df <- as.data.frame(table(
          CellType = ct, Condition = cond, Sample = sample_id
        ))
        prop_df <- prop_df[prop_df$Freq > 0 | TRUE, ]

        # Try propeller (speckle) first, fall back to chi-square
        incProgress(0.4, detail = "Statistical testing...")
        test_results <- tryCatch({
          if (!requireNamespace("speckle", quietly = TRUE)) {
            stop("speckle not installed")
          }
          res <- speckle::propeller(
            clusters = ct,
            sample = sample_id,
            group = cond
          )
          res$method <- "propeller"
          res
        }, error = function(e) {
          showNotification(
            paste0("propeller unavailable (", e$message,
                   "). Using chi-square test instead."),
            type = "warning", duration = 6
          )
          # Per-cell-type Fisher's exact test with BH correction
          tab <- table(ct, cond)
          ct_names <- rownames(tab)
          fisher_results <- lapply(ct_names, function(ct_name) {
            ct_row <- tab[ct_name, , drop = FALSE]
            rest_row <- colSums(tab[rownames(tab) != ct_name, , drop = FALSE])
            test_mat <- rbind(as.numeric(ct_row), as.numeric(rest_row))
            ft <- tryCatch(fisher.test(test_mat), error = function(e) {
              # Fall back to chi-square for large tables
              chisq.test(test_mat)
            })
            data.frame(
              CellType = ct_name,
              P.Value = ft$p.value,
              stringsAsFactors = FALSE
            )
          })
          fisher_df <- do.call(rbind, fisher_results)
          fisher_df$FDR <- p.adjust(fisher_df$P.Value, method = "BH")
          data.frame(
            BaselineProp.clusters = rowSums(tab) / sum(tab),
            PropRatio = apply(tab, 1, function(x) {
              props <- x / colSums(tab)[names(x)]
              if (length(props) >= 2) max(props) / max(min(props), 1e-6) else 1
            }),
            Tstatistic = NA_real_,
            P.Value = fisher_df$P.Value[match(ct_names, fisher_df$CellType)],
            FDR = fisher_df$FDR[match(ct_names, fisher_df$CellType)],
            method = "fisher_exact",
            row.names = ct_names
          )
        })

        incProgress(0.3, detail = "Building visualizations...")

        # Compute proportions for plotting
        prop_table <- prop.table(table(ct, cond), margin = 2)
        plot_df <- as.data.frame(prop_table)
        colnames(plot_df) <- c("CellType", "Condition", "Proportion")

        composition_results(list(
          test = test_results,
          proportions = plot_df,
          raw_table = table(ct, cond)
        ))

        showNotification("Composition analysis complete.", type = "message")
      })
    })

    # Bar plot
    output$composition_barplot <- plotly::renderPlotly({
      req(composition_results())
      plot_df <- composition_results()$proportions

      p <- ggplot(plot_df, aes(x = Condition, y = Proportion, fill = CellType)) +
        geom_bar(stat = "identity", position = "stack") +
        theme_minimal() +
        labs(x = "Condition", y = "Proportion", fill = "Cell Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      plotly::ggplotly(p, tooltip = c("fill", "y"))
    })

    # Results table
    output$composition_table <- DT::renderDataTable({
      req(composition_results())
      res <- composition_results()$test
      if (is.null(res)) return(NULL)

      # Format numeric columns
      num_cols <- sapply(res, is.numeric)
      for (col in names(which(num_cols))) {
        if (grepl("P.Value|FDR|pval", col, ignore.case = TRUE)) {
          res[[col]] <- signif(res[[col]], 3)
        } else {
          res[[col]] <- round(res[[col]], 4)
        }
      }

      DT::datatable(res,
        options = list(pageLength = 20, scrollX = TRUE, dom = "frtip"),
        rownames = TRUE, class = "compact stripe"
      )
    })

    # AI interpretation
    observeEvent(input$ai_interpret_composition, {
      req(values$ai_api_key, composition_results())
      res <- composition_results()$test

      if (is.null(res) || nrow(res) == 0) {
        showNotification("No composition results to interpret.", type = "warning")
        return()
      }

      output$ai_thread_composition <- renderUI({
        if (exists("ai_thread_loading", mode = "function")) {
          ai_thread_loading("AI interpreting composition changes...")
        } else {
          tags$div(class = "text-muted", "Analyzing...")
        }
      })

      # Build summary of top changed cell types
      if ("FDR" %in% colnames(res)) {
        sig <- res[order(res$FDR), ]
        sig <- head(sig, 10)
      } else {
        sig <- head(res, 10)
      }

      prompt <- paste0(
        "Interpret these cell type composition changes between conditions ",
        "in ", values$tissue_global %||% "tissue", " (",
        values$species_global %||% "species", "):\n\n",
        paste(capture.output(print(sig)), collapse = "\n"),
        "\n\nWhat biological processes might explain these proportion changes? ",
        "Are any changes concerning or expected?"
      )

      result <- tryCatch({
        if (exists("call_ai", mode = "function")) {
          call_ai("de_interpretation", prompt, values = values)
        } else {
          list(ok = FALSE, error = "AI not available")
        }
      }, error = function(e) {
        list(ok = FALSE, error = e$message)
      })

      output$ai_thread_composition <- renderUI({
        if (isTRUE(result$ok)) {
          rendered <- if (exists(".mdtohtml", mode = "function")) {
            .mdtohtml(result$text)
          } else {
            tags$pre(result$text)
          }
          div(class = "ai-response-content", rendered)
        } else {
          div(class = "text-danger",
            icon("exclamation-triangle"),
            paste("AI error:", result$error %||% "Unknown error"))
        }
      })
    })

  })
}
