# Track B — Security Hardening Implementation Plan

> **For agentic workers:** implement task-by-task; each task ends with a verification you can run.

**Goal:** Encrypt API keys at rest, enforce strong password hashing with transparent migration, add login rate-limiting/lockout, and provide a secure Change-Password flow — without breaking the current zero-config launch.

**Architecture:** New `modules/auth/auth_crypto.R` (AES-256-CBC via `openssl`, master key from `SCRNA_MASTER_KEY` env or an auto-generated `data/.keyfile`). `auth_db.R` encrypts on write / decrypts on read, migrates existing plaintext keys, enforces bcrypt with legacy-hash upgrade on login, and gates logins behind a `login_attempts` table. A Change-Password modal is added to the project dashboard.

**Tech Stack:** R, openssl, bcrypt, DBI/RSQLite, Shiny, testthat.

**Decisions:** Password reset = **change-password only** (email deferred). Master key = **env var + auto-gen keyfile** (zero-config).

---

## Task 1: `auth_crypto.R` — encrypt/decrypt secrets

**Files:** Create `modules/auth/auth_crypto.R`; Test `tests/testthat/test-auth-crypto.R`.

- [ ] **Step 1: Failing tests** — create `tests/testthat/test-auth-crypto.R`:

```r
Sys.setenv(SCRNA_MASTER_KEY = "unit-test-master-key")   # hermetic: no keyfile
source(file.path("..", "..", "modules", "auth", "auth_crypto.R"))

test_that("encrypt/decrypt round-trips", {
  s <- "sk-ant-api03-abcDEF123"
  enc <- encrypt_secret(s)
  expect_true(startsWith(enc, "enc:v1:"))
  expect_false(identical(enc, s))
  expect_equal(decrypt_secret(enc), s)
})

test_that("empty and NULL return empty string", {
  expect_equal(encrypt_secret(""), "")
  expect_equal(encrypt_secret(NULL), "")
  expect_equal(decrypt_secret(""), "")
  expect_equal(decrypt_secret(NULL), "")
})

test_that("legacy plaintext passes through decrypt unchanged", {
  expect_equal(decrypt_secret("sk-ant-plaintext-legacy"), "sk-ant-plaintext-legacy")
})

test_that("is_encrypted detects prefix", {
  expect_true(is_encrypted(encrypt_secret("x")))
  expect_false(is_encrypted("sk-ant-plain"))
  expect_false(is_encrypted(""))
})
```

- [ ] **Step 2: Run — expect FAIL** (`could not find function "encrypt_secret"`):
`R -e "testthat::test_file('tests/testthat/test-auth-crypto.R')"`

- [ ] **Step 3: Create `modules/auth/auth_crypto.R`:**

```r
# modules/auth/auth_crypto.R
# Symmetric encryption for secrets at rest (API keys). AES-256-CBC via openssl.
# Master key: SCRNA_MASTER_KEY env var (any string, hashed to 32 bytes) or an
# auto-generated data/.keyfile (32 random bytes, chmod 0600). Ciphertext format:
#   "enc:v1:" + base64( iv[16] || ciphertext )
# The prefix lets decrypt_secret() pass through legacy plaintext untouched.

if (!requireNamespace("openssl", quietly = TRUE))
  stop("auth_crypto.R requires the 'openssl' package.")

.crypto_key_cache <- new.env(parent = emptyenv())

get_master_key <- function() {
  if (!is.null(.crypto_key_cache$key)) return(.crypto_key_cache$key)
  env <- Sys.getenv("SCRNA_MASTER_KEY", "")
  key <- if (nzchar(env)) {
    as.raw(openssl::sha256(charToRaw(env)))         # 32 bytes
  } else {
    kf <- file.path(getwd(), "data", ".keyfile")
    if (!file.exists(kf)) {
      if (!dir.exists(dirname(kf))) dir.create(dirname(kf), recursive = TRUE)
      k <- openssl::rand_bytes(32)
      writeBin(k, kf)
      tryCatch(Sys.chmod(kf, mode = "0600"), error = function(e) NULL)
      k
    } else {
      readBin(kf, what = "raw", n = 32)
    }
  }
  .crypto_key_cache$key <- key
  key
}

is_encrypted <- function(x) {
  is.character(x) && length(x) == 1 && !is.na(x) && startsWith(x, "enc:v1:")
}

encrypt_secret <- function(plaintext) {
  if (is.null(plaintext) || length(plaintext) == 0) return("")
  plaintext <- as.character(plaintext)
  if (!nzchar(plaintext)) return("")
  if (is_encrypted(plaintext)) return(plaintext)   # already encrypted
  iv  <- openssl::rand_bytes(16)
  ct  <- openssl::aes_cbc_encrypt(charToRaw(plaintext), key = get_master_key(), iv = iv)
  paste0("enc:v1:", openssl::base64_encode(c(iv, as.raw(ct))))
}

decrypt_secret <- function(ciphertext) {
  if (is.null(ciphertext) || length(ciphertext) == 0) return("")
  ciphertext <- as.character(ciphertext)
  if (!nzchar(ciphertext)) return("")
  if (!is_encrypted(ciphertext)) return(ciphertext) # legacy plaintext
  raw <- openssl::base64_decode(sub("^enc:v1:", "", ciphertext))
  iv  <- raw[1:16]; ct <- raw[-(1:16)]
  tryCatch(rawToChar(openssl::aes_cbc_decrypt(ct, key = get_master_key(), iv = iv)),
           error = function(e) "")
}
```

- [ ] **Step 4: Run — expect PASS:** `R -e "testthat::test_file('tests/testthat/test-auth-crypto.R')"`

- [ ] **Step 5: Source before auth_db.R** — in `app.R`, change the auth block:

```r
safe_source("modules/auth/auth_crypto.R")
safe_source("modules/auth/auth_db.R")
```

- [ ] **Step 6: Checkpoint** — tests pass; `data/.keyfile` created with 0600 perms when no env var set.

---

## Task 2: Encrypt keys in `auth_db.R` (write/read) + migrate existing

**Files:** Modify `modules/auth/auth_db.R`.

- [ ] **Step 1: Encrypt on create** — in `create_project`, encrypt the key before INSERT. Replace the `dbExecute(con, "INSERT ...` param list value for `ai_api_key` by wrapping it. Concretely, just before the `tryCatch({` that runs the INSERT, add:

```r
  ai_api_key <- if (exists("encrypt_secret")) encrypt_secret(ai_api_key) else ai_api_key
```

- [ ] **Step 2: Encrypt on update** — in `update_project_metadata`, encrypt when `ai_api_key` is among the fields. Right after `params <- params[names(params) %in% allowed]`, add:

```r
  if (!is.null(params$ai_api_key) && exists("encrypt_secret"))
    params$ai_api_key <- encrypt_secret(params$ai_api_key)
```

- [ ] **Step 3: Decrypt on read** — in `get_project`, decrypt before returning. Replace the function body's return with:

```r
get_project <- function(project_id, user_id) {
  con <- .db_con(); on.exit(dbDisconnect(con))
  df <- dbGetQuery(con, "SELECT * FROM projects WHERE id=? AND user_id=?",
                   list(project_id, user_id))
  if (nrow(df) > 0 && "ai_api_key" %in% names(df) && exists("decrypt_secret"))
    df$ai_api_key <- vapply(df$ai_api_key, decrypt_secret, character(1), USE.NAMES = FALSE)
  df
}
```

> `get_user_projects` does NOT select `ai_api_key` (only card metadata), so no change needed there.

- [ ] **Step 4: One-time migration** — add a migration function and call it from `init_auth_db`. Add before `message("[Auth DB] Initialised...`:

```r
  # Encrypt any legacy plaintext API keys in place
  if (exists("encrypt_secret") && exists("is_encrypted")) {
    tryCatch({
      rows <- dbGetQuery(con, "SELECT id, ai_api_key FROM projects
                               WHERE ai_api_key IS NOT NULL AND ai_api_key <> ''")
      if (nrow(rows) > 0) {
        n <- 0L
        for (i in seq_len(nrow(rows))) {
          k <- rows$ai_api_key[i]
          if (!is_encrypted(k)) {
            dbExecute(con, "UPDATE projects SET ai_api_key=? WHERE id=?",
                      list(encrypt_secret(k), rows$id[i]))
            n <- n + 1L
          }
        }
        if (n > 0) message("[Auth DB] Encrypted ", n, " legacy API key(s).")
      }
    }, error = function(e) warning("[Auth DB] key migration: ", e$message))
  }
```

- [ ] **Step 5: Verification (R script)** — round-trip through the DB layer:

```r
R -e '
source("modules/ai_models.R"); source("modules/auth/auth_crypto.R"); source("modules/auth/auth_db.R")
init_auth_db()
u <- create_user(paste0("t",as.integer(Sys.time())), paste0("t",as.integer(Sys.time()),"@x.io"), "password12")
p <- create_project(u$user_id, "enc-test", ai_api_key="sk-ant-SECRET-123")
# raw DB row must be ciphertext:
con <- DBI::dbConnect(RSQLite::SQLite(), AUTH_DB_PATH)
raw <- DBI::dbGetQuery(con, "SELECT ai_api_key FROM projects WHERE id=?", list(p$project_id))$ai_api_key
DBI::dbDisconnect(con)
cat("stored starts enc:v1: ->", startsWith(raw, "enc:v1:"), "\n")
# get_project must decrypt:
cat("decrypted ok ->", identical(get_project(p$project_id, u$user_id)$ai_api_key, "sk-ant-SECRET-123"), "\n")
'
```
Expected: both `TRUE`.

- [ ] **Step 6: Checkpoint** — the create-modal hint "Encrypted at rest in the project database." is now truthful.

---

## Task 3: Enforce bcrypt + rate-limit/lockout in `auth_db.R`

**Files:** Modify `modules/auth/auth_db.R`.

- [ ] **Step 1: bcrypt-only hashing** — replace `hash_password`:

```r
hash_password <- function(password) {
  if (!requireNamespace("bcrypt", quietly = TRUE))
    stop("bcrypt package is required for secure password hashing.")
  bcrypt::hashpw(password)
}
```
(`verify_password` already handles both bcrypt `$2` and legacy `sha256:` — leave it.)

- [ ] **Step 2: login_attempts table** — in `init_auth_db`, after the `auth_sessions` table block, add:

```r
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS login_attempts (
      identifier   TEXT PRIMARY KEY,
      attempts     INTEGER DEFAULT 0,
      locked_until DATETIME,
      last_attempt DATETIME DEFAULT CURRENT_TIMESTAMP
    )")
```

- [ ] **Step 3: lockout + legacy-hash upgrade in `authenticate_user`** — replace the body from the `user <- dbGetQuery(...)` call through the token creation with:

```r
  # --- lockout check ---
  MAX_ATTEMPTS <- 5L; LOCK_MINUTES <- 15L
  la <- dbGetQuery(con, "SELECT attempts, locked_until FROM login_attempts WHERE identifier=?",
                   list(login_input))
  if (nrow(la) == 1 && !is.na(la$locked_until[1]) &&
      as.POSIXct(la$locked_until[1], tz = "UTC") > as.POSIXct(Sys.time(), tz = "UTC")) {
    return(list(success = FALSE,
                message = paste0("Too many failed attempts. Try again after ",
                                 format(as.POSIXct(la$locked_until[1]), "%H:%M"), ".")))
  }

  user <- dbGetQuery(con,
    "SELECT id, username, email, password_hash, full_name, is_active
     FROM users WHERE username=? OR email=?", list(login_input, login_input))

  record_fail <- function() {
    now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cur <- if (nrow(la) == 1) la$attempts[1] else 0L
    n <- cur + 1L
    lock <- if (n >= MAX_ATTEMPTS)
      format(Sys.time() + LOCK_MINUTES * 60, "%Y-%m-%d %H:%M:%S") else NA_character_
    dbExecute(con,
      "INSERT INTO login_attempts (identifier, attempts, locked_until, last_attempt)
       VALUES (?,?,?,?)
       ON CONFLICT(identifier) DO UPDATE SET attempts=?, locked_until=?, last_attempt=?",
      list(login_input, n, lock, now, n, lock, now))
  }

  if (nrow(user) == 0 || !isTRUE(as.logical(user$is_active[1]))) {
    record_fail(); return(list(success = FALSE, message = "Invalid username / password."))
  }
  if (!verify_password(password, user$password_hash[1])) {
    record_fail(); return(list(success = FALSE, message = "Invalid username / password."))
  }

  # success: clear attempts, upgrade legacy sha256 hash to bcrypt transparently
  dbExecute(con, "DELETE FROM login_attempts WHERE identifier=?", list(login_input))
  if (startsWith(user$password_hash[1], "sha256:")) {
    tryCatch(dbExecute(con, "UPDATE users SET password_hash=? WHERE id=?",
                       list(hash_password(password), user$id[1])),
             error = function(e) NULL)
  }
  dbExecute(con, "UPDATE users SET last_login=CURRENT_TIMESTAMP WHERE id=?", list(user$id[1]))
  token      <- .generate_token()
  expires_at <- format(Sys.time() + 7 * 86400, "%Y-%m-%d %H:%M:%S")
  dbExecute(con, "INSERT INTO auth_sessions (user_id, token, expires_at) VALUES (?,?,?)",
            list(user$id[1], token, expires_at))
  list(success = TRUE, message = "Login successful",
       user_id = user$id[1], username = user$username[1], email = user$email[1],
       full_name = user$full_name[1] %||% user$username[1], token = token)
```

- [ ] **Step 4: Verification** — after 5 wrong passwords, the 6th attempt (even with the right password) is rejected with a lockout message; correct login after unlock clears attempts. Test in the app or via a small R script that calls `authenticate_user` 6×.

- [ ] **Step 5: Checkpoint** — new accounts store `$2...` bcrypt hashes; a pre-existing `sha256:` account is re-hashed to bcrypt on next successful login.

---

## Task 4: Change-Password flow

**Files:** Modify `modules/auth/auth_db.R` (function) and `modules/auth/project_module.R` (UI + server).

- [ ] **Step 1: DB function** — add to `auth_db.R`:

```r
#' Change a user's password (verifies the current one).
change_user_password <- function(user_id, current_pw, new_pw) {
  if (nchar(new_pw) < 8)
    return(list(success = FALSE, message = "New password must be at least 8 characters."))
  con <- .db_con(); on.exit(dbDisconnect(con))
  r <- dbGetQuery(con, "SELECT password_hash FROM users WHERE id=?", list(user_id))
  if (nrow(r) == 0) return(list(success = FALSE, message = "User not found."))
  if (!verify_password(current_pw, r$password_hash[1]))
    return(list(success = FALSE, message = "Current password is incorrect."))
  dbExecute(con, "UPDATE users SET password_hash=? WHERE id=?",
            list(hash_password(new_pw), user_id))
  list(success = TRUE, message = "Password updated.")
}
```

- [ ] **Step 2: UI button** — in `projectDashboardUI`, in the header user actions (next to `logout_btn` at the `actionButton(ns("logout_btn"), ...)`), add before it:

```r
          actionButton(ns("change_pw_btn"), "Change Password",
                       icon  = icon("key"),
                       class = "btn-proj-logout"),
```

- [ ] **Step 3: Modal** — in `projectDashboardUI`, after the Edit modal `tags$div(id = ns("modal_edit"), ...)`, add:

```r
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
```

- [ ] **Step 4: Server** — in `projectDashboardServer`, before the logout observer, add:

```r
    pw_msg <- reactiveVal(NULL)
    output$pw_alert <- renderUI({
      m <- pw_msg(); if (is.null(m)) return(NULL)
      tags$div(class = paste0("auth-alert auth-alert--", m$type), m$text)
    })
    observeEvent(input$change_pw_btn, {
      pw_msg(NULL)
      updateTextInput(session, "pw_current", value = "")
      updateTextInput(session, "pw_new", value = "")
      updateTextInput(session, "pw_confirm", value = "")
      shinyjs::runjs(sprintf("projShowModal('%s')", ns("modal_pw")))
    })
    observeEvent(input$confirm_pw, {
      if (!identical(input$pw_new, input$pw_confirm)) {
        pw_msg(list(type = "error", text = "New passwords do not match.")); return()
      }
      res <- change_user_password(auth_state$user_id, input$pw_current %||% "", input$pw_new %||% "")
      if (isTRUE(res$success)) {
        shinyjs::runjs(sprintf("projHideModal('%s')", ns("modal_pw")))
        showNotification("Password updated.", type = "message", duration = 4)
      } else {
        pw_msg(list(type = "error", text = res$message))
      }
    })
```

- [ ] **Step 5: Checkpoint** — logged-in user changes password (wrong current → rejected; mismatch → rejected; success → can log in with new password).

---

## Task 5: Final verification

- [ ] **Full test suite:** `R -e "testthat::test_dir('tests/testthat')"` — all green (incl. new `test-auth-crypto.R`).
- [ ] **Parse-check:** `auth_crypto.R`, `auth_db.R`, `project_module.R`, `app.R`.
- [ ] **Headless boot:** launch the app; log confirms `[Auth DB] Initialised` and (if legacy keys existed) `Encrypted N legacy API key(s)`; `Listening on …`.
- [ ] **Manual:** register → create project with a key → confirm the raw `projects.ai_api_key` in `data/scrna4_auth.db` is `enc:v1:…`; open project → AI feature works (key decrypts); 5 bad logins → lockout; Change Password works.

---

## Self-review
- Coverage: encryption (T1-T2), bcrypt enforcement + legacy upgrade (T3), rate-limit/lockout (T3), change-password (T4), truthful UI claim (T2). Email reset intentionally deferred (per decision).
- Consistency: `encrypt_secret/decrypt_secret/is_encrypted/get_master_key/change_user_password` used with defined signatures.
