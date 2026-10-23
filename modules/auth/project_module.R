# =============================================================================
# project_module.R — Project Dashboard Shiny Module
# scRNA-seq Platform
# =============================================================================

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

#' Project Dashboard UI
#' @param id  Module namespace
#' @param username  Display name (reactive-safe via renderUI in caller)
projectDashboardUI <- function(id, username = "Researcher") {
  ns <- NS(id)

  tagList(
    tags$div(class = "proj-dashboard",

      # ── top header bar ──────────────────────────────────────────────────
      tags$div(class = "proj-header",
        # left: logo + brand
        tags$div(class = "proj-header-brand",
          tags$div(class = "proj-header-logo",
            tags$i(class = "fa fa-dna")
          ),
          tags$div(
            tags$span("scRNA-seq Platform", class = "proj-header-title"),
            tags$span("Project Workspace",  class = "proj-header-sub")
          )
        ),

        # right: user chip + logout
        tags$div(class = "proj-header-user",
          tags$div(class = "proj-user-chip",
            tags$div(class = "proj-user-avatar",
              uiOutput(ns("user_initials"))
            ),
            tags$div(class = "proj-user-info",
              tags$span(class = "proj-user-name",  uiOutput(ns("user_display_name"))),
              tags$span(class = "proj-user-email", uiOutput(ns("user_email_display")))
            )
          ),
          actionButton(ns("change_pw_btn"), "Change Password",
                       icon  = icon("key"),
                       class = "btn-proj-logout"),
          actionButton(ns("logout_btn"), "Sign Out",
                       icon  = icon("sign-out-alt"),
                       class = "btn-proj-logout")
        )
      ),

      # ── main content area ────────────────────────────────────────────────
      tags$div(class = "proj-main",

        # section toolbar
        tags$div(class = "proj-toolbar",
          tags$div(class = "proj-toolbar-left",
            tags$h2(class = "proj-section-title",
              tags$i(class = "fa fa-folder-open"), " My Projects"
            ),
            uiOutput(ns("proj_count_badge"))
          ),
          tags$div(class = "proj-toolbar-right",
            actionButton(ns("refresh_btn"), "",
                         icon  = icon("sync-alt"),
                         class = "btn-proj-icon",
                         title = "Refresh"),
            actionButton(ns("new_proj_btn"), tagList(
                           tags$i(class = "fa fa-plus"), " New Project"),
                         class = "btn-proj-new")
          )
        ),

        # projects grid
        uiOutput(ns("projects_grid"))
      )
    ),

    # ── Create Project modal ─────────────────────────────────────────────────
    tags$div(id = ns("modal_create"), class = "proj-modal-overlay proj-modal--hidden",
      tags$div(class = "proj-modal proj-modal--wide",

        tags$div(class = "proj-modal-header",
          tags$div(class = "proj-modal-title",
            tags$div(class = "proj-modal-title-icon",
              tags$i(class = "fa fa-plus")
            ),
            tags$div(
              tags$h4("Create new project"),
              tags$p(class = "proj-modal-subtitle",
                     "Set up your scRNA-seq workspace. All fields except the name are optional and can be edited later.")
            )
          ),
          tags$button(class = "proj-modal-close", type = "button",
                     onclick = sprintf("projHideModal('%s')", ns("modal_create")),
                     `aria-label` = "Close",
                     HTML("&times;"))
        ),

        tags$div(class = "proj-modal-body",

          # =============================================================
          # SECTION 1 — Project basics
          # =============================================================
          tags$section(class = "proj-section",
            tags$header(class = "proj-section-head",
              tags$div(class = "proj-section-num", "1"),
              tags$div(
                tags$h5("Project basics"),
                tags$p("Give the project a clear name and describe the experiment.")
              )
            ),
            tags$div(class = "proj-section-body",
              tags$div(class = "proj-form-field",
                tags$label(`for` = ns("new_name"),
                           "Project name ",
                           tags$span(class = "proj-req-dot", "required")),
                textInput(ns("new_name"), label = NULL,
                          placeholder = "e.g., PBMC_COVID_Study_2025", width = "100%")
              ),
              tags$div(class = "proj-form-field",
                tags$label(`for` = ns("new_desc"), "Description"),
                textAreaInput(ns("new_desc"), label = NULL,
                              placeholder = "Brief description of your experiment…",
                              rows = 2, width = "100%"),
                tags$p(class = "proj-field-hint",
                       "Shown on the project card — keep it short and scannable.")
              ),
              tags$div(class = "proj-grid-2",
                tags$div(class = "proj-form-field",
                  tags$label(`for` = ns("new_species"), "Species"),
                  selectInput(ns("new_species"), label = NULL,
                              choices = list(
                                "— Not specified —"       = "",
                                "Homo sapiens"            = "Homo sapiens",
                                "Mus musculus"            = "Mus musculus",
                                "Rattus norvegicus"       = "Rattus norvegicus",
                                "Danio rerio"             = "Danio rerio",
                                "Drosophila melanogaster" = "Drosophila melanogaster",
                                "Caenorhabditis elegans"  = "Caenorhabditis elegans",
                                "Other"                   = "Other"
                              ),
                              selected = "", selectize = FALSE, width = "100%")
                ),
                tags$div(class = "proj-form-field",
                  tags$label(`for` = ns("new_tissue"), "Tissue / organ"),
                  textInput(ns("new_tissue"), label = NULL,
                            placeholder = "e.g., PBMC, lung, brain", width = "100%")
                )
              ),
              tags$div(class = "proj-form-field",
                tags$label(`for` = ns("new_conditions"), "Experimental conditions"),
                textInput(ns("new_conditions"), label = NULL,
                          placeholder = "e.g., Control vs. Treatment", width = "100%")
              )
            )
          ),

          # =============================================================
          # SECTION 2 — Biological context
          # =============================================================
          tags$section(class = "proj-section",
            tags$header(class = "proj-section-head",
              tags$div(class = "proj-section-num", "2"),
              tags$div(
                tags$h5("Biological context"),
                tags$p("Helps the AI assistant and downstream modules tailor suggestions to your study.")
              )
            ),
            tags$div(class = "proj-section-body",
              tags$div(class = "proj-form-field",
                tags$label(`for` = ns("new_bio_question"), "Biological question"),
                textAreaInput(ns("new_bio_question"), label = NULL,
                              placeholder = "What question is this experiment designed to answer?",
                              rows = 2, width = "100%")
              ),
              tags$div(class = "proj-grid-2",
                tags$div(class = "proj-form-field",
                  tags$label(`for` = ns("new_exp_design"), "Experimental design"),
                  textAreaInput(ns("new_exp_design"), label = NULL,
                                placeholder = "Samples, replicates, comparisons…",
                                rows = 3, width = "100%")
                ),
                tags$div(class = "proj-form-field",
                  tags$label(`for` = ns("new_background"), "Background"),
                  textAreaInput(ns("new_background"), label = NULL,
                                placeholder = "Prior knowledge, hypotheses, key literature…",
                                rows = 3, width = "100%")
                )
              )
            )
          ),

          # =============================================================
          # SECTION 3 — AI assistant setup
          # =============================================================
          tags$section(class = "proj-section",
            tags$header(class = "proj-section-head",
              tags$div(class = "proj-section-num", "3"),
              tags$div(
                tags$h5("AI assistant setup"),
                tags$p("Your API key is stored per-project and is never shared between projects.")
              )
            ),
            tags$div(class = "proj-section-body",
              tags$div(class = "proj-grid-2",
                tags$div(class = "proj-form-field",
                  tags$label(`for` = ns("new_ai_api_key"), "Claude API key"),
                  passwordInput(ns("new_ai_api_key"), label = NULL,
                                placeholder = "sk-ant-…", width = "100%"),
                  tags$p(class = "proj-field-hint",
                         tags$i(class = "fa fa-lock"),
                         " Encrypted at rest in the project database.")
                ),
                tags$div(class = "proj-form-field",
                  tags$label(`for` = ns("new_ai_model"), "Model"),
                  selectInput(ns("new_ai_model"), label = NULL,
                              choices = ai_model_choices(),
                              selected = get_default_model(),
                              selectize = FALSE, width = "100%")
                )
              ),
              tags$div(class = "proj-form-field",
                tags$label(`for` = ns("new_ai_persona"), "Custom persona (optional)"),
                textAreaInput(ns("new_ai_persona"), label = NULL,
                              placeholder = "Extra instructions appended to the default expert persona…",
                              rows = 2, width = "100%")
              ),
              tags$div(class = "proj-form-field",
                tags$div(class = "proj-label-row",
                  tags$label(`for` = ns("new_ai_summary"), "Project summary"),
                  actionButton(ns("btn_gen_summary"),
                               tagList(tags$i(class = "fa fa-wand-magic-sparkles"), " Generate with AI"),
                               class = "btn-proj-generate")
                ),
                textAreaInput(ns("new_ai_summary"), label = NULL,
                              placeholder = "1-2 sentence summary — or click Generate with AI to let Claude draft one from the fields above.",
                              rows = 3, width = "100%"),
                tags$div(class = "proj-ai-summary-status",
                         uiOutput(ns("summary_status")))
              )
            )
          ),

          uiOutput(ns("create_alert"))
        ),

        tags$div(class = "proj-modal-footer",
          tags$p(class = "proj-modal-footer-note",
                 tags$i(class = "fa fa-circle-info"),
                 " You can edit every field later from the project menu."),
          tags$div(class = "proj-modal-footer-actions",
            tags$button(type = "button", class = "btn-proj-cancel",
                       onclick = sprintf("projHideModal('%s')", ns("modal_create")),
                       "Cancel"),
            actionButton(ns("confirm_create"), tagList(
                           tags$i(class = "fa fa-check"), " Create project"),
                         class = "btn-proj-confirm")
          )
        )
      )
    ),

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
                               "Drosophila melanogaster" = "Drosophila melanogaster",
                               "Caenorhabditis elegans" = "Caenorhabditis elegans",
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

    # ── Change Password modal ────────────────────────────────────────────────
    tags$div(id = ns("modal_pw"), class = "proj-modal-overlay proj-modal--hidden",
      tags$div(class = "proj-modal",
        tags$div(class = "proj-modal-header",
          tags$div(class = "proj-modal-title",
            tags$div(class = "proj-modal-title-icon", tags$i(class = "fa fa-key")),
            tags$div(tags$h4("Change password"))),
          tags$button(class = "proj-modal-close", type = "button",
                      onclick = sprintf("projHideModal('%s')", ns("modal_pw")),
                      `aria-label` = "Close", HTML("&times;"))),
        tags$div(class = "proj-modal-body",
          tags$div(class = "proj-form-field",
            tags$label("Current password"),
            passwordInput(ns("pw_current"), label = NULL, width = "100%")),
          tags$div(class = "proj-form-field",
            tags$label("New password"),
            passwordInput(ns("pw_new"), label = NULL, width = "100%")),
          tags$div(class = "proj-form-field",
            tags$label("Confirm new password"),
            passwordInput(ns("pw_confirm"), label = NULL, width = "100%")),
          uiOutput(ns("pw_alert"))),
        tags$div(class = "proj-modal-footer",
          tags$div(class = "proj-modal-footer-actions",
            tags$button(type = "button", class = "btn-proj-cancel",
                        onclick = sprintf("projHideModal('%s')", ns("modal_pw")), "Cancel"),
            actionButton(ns("confirm_pw"),
                         tagList(tags$i(class = "fa fa-check"), " Update"),
                         class = "btn-proj-confirm")))
      )
    ),

    tags$style(HTML(".btn-proj-edit{background:transparent;border:1px solid #d0d7de;color:#57606a;border-radius:6px;padding:6px 10px;cursor:pointer;margin-right:6px}.btn-proj-edit:hover{background:#0d6efd;color:#fff;border-color:#0d6efd}")),

    # ── inline JS helpers ─────────────────────────────────────────────────────
    tags$script(HTML("
window.projShowModal = function(id) {
  var el = document.getElementById(id);
  if (el) { el.classList.remove('proj-modal--hidden'); el.classList.add('proj-modal--visible'); }
};
window.projHideModal = function(id) {
  var el = document.getElementById(id);
  if (el) { el.classList.remove('proj-modal--visible'); el.classList.add('proj-modal--hidden'); }
};
    "))
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

#' Project Dashboard Server
projectDashboardServer <- function(id, auth_state) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ── reactive trigger for refreshing ─────────────────────────────────────
    refresh_trigger <- reactiveVal(0)

    # ── user display ─────────────────────────────────────────────────────────
    output$user_initials <- renderUI({
      name <- auth_state$full_name %||% auth_state$username %||% "?"
      initials <- paste(toupper(substr(strsplit(name, " ")[[1]], 1, 1)), collapse = "")
      initials <- substr(initials, 1, 2)
      tags$span(initials)
    })
    output$user_display_name <- renderUI({
      tags$span(auth_state$full_name %||% auth_state$username %||% "Researcher")
    })
    output$user_email_display <- renderUI({
      tags$span(auth_state$email %||% "")
    })

    # ── load projects ────────────────────────────────────────────────────────
    projects_rv <- reactive({
      refresh_trigger()
      req(auth_state$user_id)
      get_user_projects(auth_state$user_id)
    })

    output$proj_count_badge <- renderUI({
      n <- nrow(projects_rv())
      if (n == 0) return(NULL)
      tags$span(class = "proj-count-badge",
                sprintf("%d project%s", n, if (n == 1) "" else "s"))
    })

    observeEvent(input$refresh_btn, refresh_trigger(refresh_trigger() + 1))

    # ── open create-modal ─────────────────────────────────────────────────────
    observeEvent(input$new_proj_btn, {
      shinyjs::runjs(sprintf("projShowModal('%s')", ns("modal_create")))
    })
    observeEvent(input$new_proj_btn_empty, {
      shinyjs::runjs(sprintf("projShowModal('%s')", ns("modal_create")))
    })

    # create-alert
    create_msg <- reactiveVal(NULL)
    output$create_alert <- renderUI({
      m <- create_msg(); if (is.null(m)) return(NULL)
      tags$div(class = paste0("auth-alert auth-alert--", m$type), m$text)
    })

    # ── AI summary generator ──────────────────────────────────────────────────
    summary_status_msg <- reactiveVal(NULL)
    output$summary_status <- renderUI({
      m <- summary_status_msg(); if (is.null(m)) return(NULL)
      tags$small(style = paste0("color:",
                                if (m$type == "error") "#dc3545" else "#28a745",
                                "; display:block; margin-top:4px;"),
                 m$text)
    })

    observeEvent(input$btn_gen_summary, {
      api_key <- trimws(input$new_ai_api_key %||% "")
      if (!nzchar(api_key)) {
        summary_status_msg(list(type = "error",
                                text = "Enter an API key above before generating."))
        return()
      }
      if (!exists("call_claude_api", mode = "function")) {
        summary_status_msg(list(type = "error",
                                text = "AI helper not loaded."))
        return()
      }

      species    <- trimws(input$new_species %||% "")
      tissue     <- trimws(input$new_tissue  %||% "")
      conditions <- trimws(input$new_conditions %||% "")
      bio_q      <- trimws(input$new_bio_question %||% "")
      exp_d      <- trimws(input$new_exp_design  %||% "")
      bg         <- trimws(input$new_background  %||% "")

      if (!nzchar(species) && !nzchar(tissue) && !nzchar(bio_q)) {
        summary_status_msg(list(type = "error",
                                text = "Fill in at least species, tissue, or biological question first."))
        return()
      }

      summary_status_msg(list(type = "info", text = "Generating…"))

      prompt <- paste0(
        "Produce a concise 1-2 sentence project summary for an scRNA-seq study.\n",
        "Return only the summary, no preamble.\n\n",
        "Species: ",    if (nzchar(species))    species    else "(unspecified)", "\n",
        "Tissue: ",     if (nzchar(tissue))     tissue     else "(unspecified)", "\n",
        "Conditions: ", if (nzchar(conditions)) conditions else "(unspecified)", "\n",
        if (nzchar(bio_q)) paste0("Biological question: ", bio_q, "\n") else "",
        if (nzchar(exp_d)) paste0("Experimental design: ", exp_d, "\n") else "",
        if (nzchar(bg))    paste0("Background: ",         bg,    "\n") else ""
      )

      result <- tryCatch(
        call_claude_api(api_key = api_key,
                        model   = input$new_ai_model %||% get_default_model(),
                        messages = list(list(role = "user", content = prompt))),
        error = function(e) list(success = FALSE, error = conditionMessage(e))
      )

      if (isTRUE(result$ok) && nzchar(result$text %||% "")) {
        updateTextAreaInput(session, "new_ai_summary",
                            value = trimws(result$text))
        summary_status_msg(list(type = "ok", text = "Summary generated."))
      } else {
        summary_status_msg(list(type = "error",
                                text = paste("Generation failed:",
                                             result$error %||% "unknown error")))
      }
    })

    # ── refresh model list live once a key is entered ─────────────────────────
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

    # ── confirm project creation ──────────────────────────────────────────────
    observeEvent(input$confirm_create, {
      name <- trimws(input$new_name %||% "")
      if (nchar(name) == 0) {
        create_msg(list(type = "error", text = "Project name is required."))
        return()
      }

      result <- create_project(
        user_id              = auth_state$user_id,
        name                 = name,
        description          = input$new_desc       %||% "",
        species              = input$new_species    %||% "",
        tissue               = input$new_tissue     %||% "",
        conditions           = input$new_conditions %||% "",
        biological_question  = input$new_bio_question %||% "",
        experiment_design    = input$new_exp_design   %||% "",
        project_background   = input$new_background   %||% "",
        ai_api_key           = input$new_ai_api_key   %||% "",
        ai_model             = input$new_ai_model     %||% get_default_model(),
        ai_custom_persona    = input$new_ai_persona   %||% "",
        ai_project_summary   = input$new_ai_summary   %||% ""
      )

      if (isTRUE(result$success)) {
        create_msg(NULL)
        summary_status_msg(NULL)
        shinyjs::runjs(sprintf("projHideModal('%s')", ns("modal_create")))
        # clear form
        updateTextInput(session, "new_name",            value = "")
        updateTextAreaInput(session, "new_desc",        value = "")
        updateTextInput(session, "new_tissue",          value = "")
        updateTextInput(session, "new_conditions",      value = "")
        updateTextAreaInput(session, "new_bio_question", value = "")
        updateTextAreaInput(session, "new_exp_design",   value = "")
        updateTextAreaInput(session, "new_background",   value = "")
        updateTextInput(session, "new_ai_api_key",      value = "")
        updateTextAreaInput(session, "new_ai_persona",  value = "")
        updateTextAreaInput(session, "new_ai_summary",  value = "")
        refresh_trigger(refresh_trigger() + 1)
        showNotification(paste0("Project \"", name, "\" created!"),
                         type = "message", duration = 4)
      } else {
        create_msg(list(type = "error", text = result$message))
      }
    })

    # ── projects grid ─────────────────────────────────────────────────────────
    output$projects_grid <- renderUI({
      projs <- projects_rv()

      if (nrow(projs) == 0) {
        return(
          tags$div(class = "proj-empty-state",
            tags$div(class = "proj-empty-icon",
              tags$i(class = "fa fa-flask")
            ),
            tags$h3("No projects yet"),
            tags$p("Create your first project to start analysing your single-cell data."),
            actionButton(ns("new_proj_btn_empty"),
                         tagList(tags$i(class = "fa fa-plus"), " Create First Project"),
                         class = "btn-proj-new")
          )
        )
      }

      cards <- lapply(seq_len(nrow(projs)), function(i) {
        p   <- projs[i, ]
        si  <- get_status_info(p$status %||% "empty")

        # species icon
        sp_icon <- if      (grepl("sapiens",  p$species %||% "")) "👤"
                   else if (grepl("musculus", p$species %||% "")) "🐭"
                   else if (grepl("rerio",    p$species %||% "")) "🐟"
                   else "🧬"

        # relative time
        rel_time <- tryCatch({
          dt <- as.POSIXct(p$last_modified)
          d  <- as.numeric(difftime(Sys.time(), dt, units = "secs"))
          if      (d < 60)    "just now"
          else if (d < 3600)  sprintf("%d min ago",  floor(d / 60))
          else if (d < 86400) sprintf("%d hr ago",   floor(d / 3600))
          else                sprintf("%d days ago", floor(d / 86400))
        }, error = function(e) "—")

        # truncate description
        desc_text <- p$description %||% ""
        if (nchar(desc_text) > 80) desc_text <- paste0(substr(desc_text, 1, 78), "…")

        tags$div(class = "proj-card",

          # card header
          tags$div(class = "proj-card-header",
            tags$div(class = "proj-card-species-icon", sp_icon),
            tags$div(class = "proj-card-title-wrap",
              tags$h4(class = "proj-card-name", p$name),
              tags$span(
                class = "proj-status-badge",
                style = paste0(
                  "background:", si$color, "18;",
                  "color:", si$color, ";",
                  "border:1px solid ", si$color, "40;"
                ),
                si$label
              )
            )
          ),

          # description
          tags$p(class = if (nchar(desc_text) == 0) "proj-card-desc proj-card-desc--empty"
                          else "proj-card-desc",
                 if (nchar(desc_text) == 0) "No description" else desc_text),

          # metadata chips
          tags$div(class = "proj-card-chips",
            if (nchar(p$species %||% "") > 0)
              tags$span(class = "proj-chip", tags$i(class = "fa fa-dna"), " ", p$species),
            if (nchar(p$tissue %||% "") > 0)
              tags$span(class = "proj-chip", tags$i(class = "fa fa-microscope"), " ", p$tissue),
            if (nchar(p$conditions %||% "") > 0)
              tags$span(class = "proj-chip", tags$i(class = "fa fa-flask"), " ", p$conditions)
          ),

          # progress bar
          tags$div(class = "proj-progress",
            tags$div(class = "proj-progress-top",
              tags$span("Analysis progress"),
              tags$span(paste0(si$pct, "%"))
            ),
            tags$div(class = "proj-progress-track",
              tags$div(class = "proj-progress-fill",
                style = paste0("width:", si$pct, "%; background:", si$color, ";"))
            )
          ),

          # footer: timestamp + actions
          tags$div(class = "proj-card-footer",
            tags$span(class = "proj-card-time",
              tags$i(class = "fa fa-clock"), " ", rel_time),

            tags$div(class = "proj-card-actions",
              # Open button — sends project id via Shiny.setInputValue
              tags$button(
                class   = "btn-proj-open",
                type    = "button",
                onclick = sprintf(
                  "Shiny.setInputValue('%s', %d, {priority:'event'})",
                  ns("open_project"), p$id
                ),
                tags$i(class = "fa fa-folder-open"), " Open"
              ),
              # Edit button
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
              # Delete button
              tags$button(
                class   = "btn-proj-delete",
                type    = "button",
                title   = "Delete project",
                onclick = sprintf(
                  "if(confirm('Delete project \\u201c%s\\u201d?\\nThis cannot be undone.')) { Shiny.setInputValue('%s', %d, {priority:'event'}); }",
                  gsub("'", "\\'", p$name, fixed = TRUE),
                  ns("delete_project"), p$id
                ),
                tags$i(class = "fa fa-trash-alt")
              )
            )
          )
        )
      })

      tags$div(class = "proj-grid", cards)
    })

    # ── open project ──────────────────────────────────────────────────────────
    observeEvent(input$open_project, {
      pid  <- as.integer(input$open_project)
      proj <- get_project(pid, auth_state$user_id)
      if (nrow(proj) == 0) {
        showNotification("Project not found.", type = "error")
        return()
      }
      # Push all project metadata into auth_state for the main app to consume
      auth_state$current_project_id                    <- pid
      auth_state$current_project_name                  <- proj$name[1]
      auth_state$current_project_species               <- proj$species[1]          %||% ""
      auth_state$current_project_tissue                <- proj$tissue[1]           %||% ""
      auth_state$current_project_conditions            <- proj$conditions[1]       %||% ""
      auth_state$current_project_description           <- proj$description[1]      %||% ""
      auth_state$current_project_biological_question   <- proj$biological_question[1] %||% ""
      auth_state$current_project_experiment_design     <- proj$experiment_design[1]   %||% ""
      auth_state$current_project_background            <- proj$project_background[1]  %||% ""
      # AI assistant per-project config
      auth_state$current_project_ai_api_key            <- proj$ai_api_key[1]           %||% ""
      auth_state$current_project_ai_model              <- proj$ai_model[1]             %||% ""
      auth_state$current_project_ai_persona            <- proj$ai_custom_persona[1]    %||% ""
      auth_state$current_project_ai_summary            <- proj$ai_project_summary[1]   %||% ""
    })

    # ── delete project ────────────────────────────────────────────────────────
    observeEvent(input$delete_project, {
      pid    <- as.integer(input$delete_project)
      result <- delete_project(pid, auth_state$user_id)
      if (isTRUE(result$success)) {
        refresh_trigger(refresh_trigger() + 1)
        showNotification("Project deleted.", type = "warning", duration = 3)
      } else {
        showNotification(result$message, type = "error")
      }
    })

    # ── edit project: prefill + open modal ────────────────────────────────────
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
      # refresh model list against the saved key, preselect the saved model
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

    # ── edit project: save changes (and sync if currently open) ───────────────
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
          auth_state$current_project_name     <- name
          auth_state$current_project_species  <- fields$species
          auth_state$current_project_tissue   <- fields$tissue
          auth_state$current_project_conditions <- fields$conditions
          auth_state$current_project_ai_model <- fields$ai_model
          if (!is.null(fields$ai_api_key))
            auth_state$current_project_ai_api_key <- fields$ai_api_key
        }
        showNotification("Project updated.", type = "message", duration = 4)
      } else {
        edit_msg(list(type = "error", text = "Update failed. Check the logs."))
      }
    })

    # ── change password ───────────────────────────────────────────────────────
    pw_msg <- reactiveVal(NULL)
    output$pw_alert <- renderUI({
      m <- pw_msg(); if (is.null(m)) return(NULL)
      tags$div(class = paste0("auth-alert auth-alert--", m$type), m$text)
    })
    observeEvent(input$change_pw_btn, {
      pw_msg(NULL)
      updateTextInput(session, "pw_current", value = "")
      updateTextInput(session, "pw_new",     value = "")
      updateTextInput(session, "pw_confirm", value = "")
      shinyjs::runjs(sprintf("projShowModal('%s')", ns("modal_pw")))
    })
    observeEvent(input$confirm_pw, {
      if (!identical(input$pw_new, input$pw_confirm)) {
        pw_msg(list(type = "error", text = "New passwords do not match.")); return()
      }
      res <- change_user_password(auth_state$user_id,
                                  input$pw_current %||% "", input$pw_new %||% "")
      if (isTRUE(res$success)) {
        shinyjs::runjs(sprintf("projHideModal('%s')", ns("modal_pw")))
        showNotification("Password updated.", type = "message", duration = 4)
      } else {
        pw_msg(list(type = "error", text = res$message))
      }
    })

    # ── logout ────────────────────────────────────────────────────────────────
    observeEvent(input$logout_btn, {
      invalidate_session_token(auth_state$session_token)
      auth_state$logged_in          <- FALSE
      auth_state$user_id            <- NULL
      auth_state$username           <- NULL
      auth_state$email              <- NULL
      auth_state$full_name          <- NULL
      auth_state$session_token      <- NULL
      auth_state$current_project_id <- NULL
      auth_state$current_project    <- NULL
    })

  })
}
