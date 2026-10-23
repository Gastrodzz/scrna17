# =============================================================================
# auth_module.R — Login / Register Shiny Module
# scRNA-seq Platform
# =============================================================================

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------

#' Full-page Authentication UI (Login + Register tabs)
authPageUI <- function(id) {
  ns <- NS(id)

  tagList(
    # ── page wrapper ──────────────────────────────────────────────────────────
    tags$div(class = "auth-wrapper",

      # animated mesh background blobs
      tags$div(class = "auth-bg",
        tags$div(class = "auth-bg-blob auth-bg-blob-1"),
        tags$div(class = "auth-bg-blob auth-bg-blob-2"),
        tags$div(class = "auth-bg-blob auth-bg-blob-3")
      ),

      # ── centred card ────────────────────────────────────────────────────────
      tags$div(class = "auth-container",

        # logo / platform title
        tags$div(class = "auth-logo",
          tags$div(class = "auth-logo-icon-wrap",
            tags$i(class = "fa fa-dna")
          ),
          tags$div(class = "auth-logo-text",
            tags$h1("scRNA-seq Platform"),
            tags$p("Advanced Single-Cell Analysis Suite")
          )
        ),

        # card
        tags$div(class = "auth-card",

          # ── tab bar ─────────────────────────────────────────────────────────
          tags$div(class = "auth-tabs",
            tags$button(
              id    = ns("tab_login"),
              class = "auth-tab auth-tab--active",
              type  = "button",
              onclick = sprintf("authSwitchTab('%s','login')", ns("")),
              tags$i(class = "fa fa-sign-in-alt"), " Sign In"
            ),
            tags$button(
              id    = ns("tab_register"),
              class = "auth-tab",
              type  = "button",
              onclick = sprintf("authSwitchTab('%s','register')", ns("")),
              tags$i(class = "fa fa-user-plus"), " Register"
            )
          ),

          # ── LOGIN panel ─────────────────────────────────────────────────────
          tags$div(id = ns("panel_login"), class = "auth-panel",

            tags$div(class = "auth-panel-header",
              tags$h3("Welcome back"),
              tags$p("Sign in to continue your research")
            ),

            uiOutput(ns("login_alert")),

            tags$div(class = "auth-field",
              tags$label(tags$i(class = "fa fa-user"), " Username or Email"),
              textInput(ns("login_input"), label = NULL,
                        placeholder = "Enter username or email", width = "100%")
            ),

            tags$div(class = "auth-field",
              tags$label(tags$i(class = "fa fa-lock"), " Password"),
              passwordInput(ns("login_password"), label = NULL,
                            placeholder = "Enter your password", width = "100%")
            ),

            tags$div(class = "auth-actions",
              actionButton(ns("do_login"), label = tagList(
                tags$i(class = "fa fa-sign-in-alt"), " Sign In"
              ), class = "btn-auth-primary", width = "100%")
            ),

            tags$div(class = "auth-switch-text",
              "Don't have an account? ",
              tags$a(href = "#",
                     onclick = sprintf("authSwitchTab('%s','register'); return false;", ns("")),
                     "Create one")
            )
          ),

          # ── REGISTER panel ──────────────────────────────────────────────────
          tags$div(id = ns("panel_register"), class = "auth-panel auth-panel--hidden",

            tags$div(class = "auth-panel-header",
              tags$h3("Create Account"),
              tags$p("Join the single-cell research platform")
            ),

            uiOutput(ns("register_alert")),

            fluidRow(
              column(6,
                tags$div(class = "auth-field",
                  tags$label(tags$i(class = "fa fa-at"), " Username *"),
                  textInput(ns("reg_username"), label = NULL,
                            placeholder = "e.g. ahmed_lab", width = "100%")
                )
              ),
              column(6,
                tags$div(class = "auth-field",
                  tags$label(tags$i(class = "fa fa-id-badge"), " Full Name"),
                  textInput(ns("reg_fullname"), label = NULL,
                            placeholder = "Your name (optional)", width = "100%")
                )
              )
            ),

            tags$div(class = "auth-field",
              tags$label(tags$i(class = "fa fa-envelope"), " Email Address *"),
              textInput(ns("reg_email"), label = NULL,
                        placeholder = "your@email.com", width = "100%")
            ),

            fluidRow(
              column(6,
                tags$div(class = "auth-field",
                  tags$label(tags$i(class = "fa fa-key"), " Password *"),
                  passwordInput(ns("reg_password"), label = NULL,
                                placeholder = "Min. 8 characters", width = "100%")
                )
              ),
              column(6,
                tags$div(class = "auth-field",
                  tags$label(tags$i(class = "fa fa-key"), " Confirm Password *"),
                  passwordInput(ns("reg_confirm"), label = NULL,
                                placeholder = "Repeat password", width = "100%")
                )
              )
            ),

            # password-strength bar
            tags$div(class = "pw-strength-wrap",
              tags$div(id = ns("pw_strength_bar"), class = "pw-strength-bar"),
              tags$span(id = ns("pw_strength_label"), class = "pw-strength-label", "")
            ),

            tags$div(class = "auth-actions",
              actionButton(ns("do_register"), label = tagList(
                tags$i(class = "fa fa-user-plus"), " Create Account"
              ), class = "btn-auth-primary", width = "100%")
            ),

            tags$div(class = "auth-switch-text",
              "Already have an account? ",
              tags$a(href = "#",
                     onclick = sprintf("authSwitchTab('%s','login'); return false;", ns("")),
                     "Sign in")
            )
          )
        ),

        # page footer
        tags$div(class = "auth-footer",
          tags$span(tags$i(class="fa fa-lock"), " Secure local session"),
          tags$span(" · "),
          tags$span(tags$i(class="fa fa-dna"), " Seurat v5"),
          tags$span(" · "),
          tags$span(tags$i(class="fa fa-robot"), " Claude AI")
        )
      )
    ),

    # ── inline JS: tab switching + password strength ─────────────────────────
    tags$script(HTML(sprintf("
(function() {
  // Tab switcher
  window.authSwitchTab = function(ns, tab) {
    var loginPanel    = document.getElementById(ns + 'panel_login');
    var registerPanel = document.getElementById(ns + 'panel_register');
    var loginTab      = document.getElementById(ns + 'tab_login');
    var registerTab   = document.getElementById(ns + 'tab_register');
    if (!loginPanel || !registerPanel) return;

    if (tab === 'login') {
      loginPanel.classList.remove('auth-panel--hidden');
      registerPanel.classList.add('auth-panel--hidden');
      loginTab.classList.add('auth-tab--active');
      registerTab.classList.remove('auth-tab--active');
    } else {
      registerPanel.classList.remove('auth-panel--hidden');
      loginPanel.classList.add('auth-panel--hidden');
      registerTab.classList.add('auth-tab--active');
      loginTab.classList.remove('auth-tab--active');
    }
  };

  // Password strength meter (runs after DOM ready)
  $(document).on('input', '#%s', function() {
    var pw   = $(this).val();
    var bar  = document.getElementById('%s');
    var lbl  = document.getElementById('%s');
    if (!bar || !lbl) return;
    var score = 0;
    if (pw.length >= 8)          score++;
    if (/[A-Z]/.test(pw))        score++;
    if (/[0-9]/.test(pw))        score++;
    if (/[^A-Za-z0-9]/.test(pw)) score++;
    var map = ['','weak','fair','good','strong'];
    var col = ['','#dc3545','#fd7e14','#ffc107','#28a745'];
    bar.style.width = (score * 25) + '%%';
    bar.style.background = col[score] || '#dee2e6';
    lbl.textContent = map[score] || '';
    lbl.style.color  = col[score] || '#999';
  });
})();
    ", ns("reg_password"), ns("pw_strength_bar"), ns("pw_strength_label"))))
  )
}

# ---------------------------------------------------------------------------
# Server
# ---------------------------------------------------------------------------

#' Authentication Page Server
#' Sets auth_state$logged_in, user_id, username, email, full_name, session_token
authPageServer <- function(id, auth_state) {
  moduleServer(id, function(input, output, session) {

    # ── local alert reactives ────────────────────────────────────────────────
    login_msg    <- reactiveVal(NULL)
    register_msg <- reactiveVal(NULL)

    output$login_alert <- renderUI({
      m <- login_msg(); if (is.null(m)) return(NULL)
      tags$div(class = paste0("auth-alert auth-alert--", m$type),
               tags$i(class = if(m$type=="error") "fa fa-exclamation-circle"
                              else if(m$type=="success") "fa fa-check-circle"
                              else "fa fa-info-circle"),
               " ", m$text)
    })

    output$register_alert <- renderUI({
      m <- register_msg(); if (is.null(m)) return(NULL)
      tags$div(class = paste0("auth-alert auth-alert--", m$type),
               tags$i(class = if(m$type=="error") "fa fa-exclamation-circle"
                              else if(m$type=="success") "fa fa-check-circle"
                              else "fa fa-info-circle"),
               " ", m$text)
    })

    # ── LOGIN ────────────────────────────────────────────────────────────────
    observeEvent(input$do_login, {
      login_input <- trimws(input$login_input %||% "")
      pw          <- input$login_password %||% ""

      if (nchar(login_input) == 0 || nchar(pw) == 0) {
        login_msg(list(type = "error", text = "Please fill in all fields."))
        return()
      }

      login_msg(list(type = "info", text = "Signing in…"))

      result <- tryCatch(
        authenticate_user(login_input, pw),
        error = function(e) list(success = FALSE, message = e$message)
      )

      if (isTRUE(result$success)) {
        login_msg(NULL)
        auth_state$logged_in          <- TRUE
        auth_state$user_id            <- result$user_id
        auth_state$username           <- result$username
        auth_state$email              <- result$email
        auth_state$full_name          <- result$full_name %||% result$username
        auth_state$session_token      <- result$token
        auth_state$current_project_id <- NULL
        auth_state$current_project    <- NULL
      } else {
        login_msg(list(type = "error", text = result$message))
      }
    })

    # ── REGISTER ─────────────────────────────────────────────────────────────
    observeEvent(input$do_register, {
      uname   <- trimws(input$reg_username %||% "")
      email   <- trimws(input$reg_email    %||% "")
      pw      <- input$reg_password %||% ""
      pw2     <- input$reg_confirm  %||% ""
      fname   <- trimws(input$reg_fullname %||% "")

      if (nchar(uname) == 0 || nchar(email) == 0 || nchar(pw) == 0) {
        register_msg(list(type = "error", text = "Please fill in all required fields."))
        return()
      }
      if (!identical(pw, pw2)) {
        register_msg(list(type = "error", text = "Passwords do not match."))
        return()
      }

      register_msg(list(type = "info", text = "Creating your account…"))

      result <- tryCatch(
        create_user(uname, email, pw, fname),
        error = function(e) list(success = FALSE, message = e$message)
      )

      if (isTRUE(result$success)) {
        register_msg(list(type = "success", text = "Account created! Signing you in…"))
        Sys.sleep(0.4)   # brief pause so message is visible

        # auto-login
        login_result <- tryCatch(
          authenticate_user(uname, pw),
          error = function(e) list(success = FALSE, message = e$message)
        )
        if (isTRUE(login_result$success)) {
          register_msg(NULL)
          auth_state$logged_in          <- TRUE
          auth_state$user_id            <- login_result$user_id
          auth_state$username           <- login_result$username
          auth_state$email              <- login_result$email
          auth_state$full_name          <- login_result$full_name %||% login_result$username
          auth_state$session_token      <- login_result$token
          auth_state$current_project_id <- NULL
          auth_state$current_project    <- NULL
        }
      } else {
        register_msg(list(type = "error", text = result$message))
      }
    })

  })
}
