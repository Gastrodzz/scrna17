# =============================================================================
# auth_db.R — SQLite Database Layer for Authentication & Project Management
# scRNA-seq Platform
# =============================================================================
# Tables: users, auth_sessions, projects
# Password hashing: openssl SHA-256 + random salt (bcrypt if available)
# Tokens: openssl::rand_bytes or digest fallback
# =============================================================================

suppressPackageStartupMessages({
  for (.p in c("DBI", "RSQLite", "digest", "jsonlite")) {
    if (!requireNamespace(.p, quietly = TRUE))
      install.packages(.p, repos = "https://cran.r-project.org")
  }
  library(DBI)
  library(RSQLite)
  library(digest)
})

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
AUTH_DB_PATH  <- file.path(getwd(), "data", "scrna4_auth.db")
SESSIONS_DIR  <- file.path(getwd(), "data", "user_sessions")

# ---------------------------------------------------------------------------
# Database initialisation
# ---------------------------------------------------------------------------
init_auth_db <- function() {
  data_dir <- file.path(getwd(), "data")
  if (!dir.exists(data_dir))    dir.create(data_dir,    recursive = TRUE)
  if (!dir.exists(SESSIONS_DIR)) dir.create(SESSIONS_DIR, recursive = TRUE)

  con <- dbConnect(RSQLite::SQLite(), AUTH_DB_PATH)
  on.exit(dbDisconnect(con))

  dbExecute(con, "PRAGMA journal_mode = WAL")
  dbExecute(con, "PRAGMA foreign_keys = ON")

  # ---- users ----------------------------------------------------------------
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS users (
      id            INTEGER PRIMARY KEY AUTOINCREMENT,
      username      TEXT    UNIQUE NOT NULL,
      email         TEXT    UNIQUE NOT NULL,
      password_hash TEXT    NOT NULL,
      full_name     TEXT    DEFAULT '',
      created_at    DATETIME DEFAULT CURRENT_TIMESTAMP,
      last_login    DATETIME,
      is_active     INTEGER DEFAULT 1
    )
  ")

  # ---- auth_sessions --------------------------------------------------------
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS auth_sessions (
      id         INTEGER PRIMARY KEY AUTOINCREMENT,
      user_id    INTEGER NOT NULL,
      token      TEXT    UNIQUE NOT NULL,
      created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
      expires_at DATETIME NOT NULL,
      is_valid   INTEGER  DEFAULT 1,
      FOREIGN KEY (user_id) REFERENCES users(id) ON DELETE CASCADE
    )
  ")

  # ---- login_attempts (rate limiting / lockout) -----------------------------
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS login_attempts (
      identifier   TEXT PRIMARY KEY,
      attempts     INTEGER DEFAULT 0,
      locked_until DATETIME,
      last_attempt DATETIME DEFAULT CURRENT_TIMESTAMP
    )
  ")

  # ---- projects -------------------------------------------------------------
  dbExecute(con, "
    CREATE TABLE IF NOT EXISTS projects (
      id                   INTEGER PRIMARY KEY AUTOINCREMENT,
      user_id              INTEGER NOT NULL,
      name                 TEXT    NOT NULL,
      description          TEXT    DEFAULT '',
      species              TEXT    DEFAULT '',
      tissue               TEXT    DEFAULT '',
      conditions           TEXT    DEFAULT '',
      biological_question  TEXT    DEFAULT '',
      experiment_design    TEXT    DEFAULT '',
      project_background   TEXT    DEFAULT '',
      ai_api_key           TEXT    DEFAULT '',
      ai_model             TEXT    DEFAULT '',
      ai_custom_persona    TEXT    DEFAULT '',
      ai_project_summary   TEXT    DEFAULT '',
      rds_path             TEXT,
      status               TEXT    DEFAULT 'empty',
      analysis_flags       TEXT    DEFAULT '{}',
      created_at           DATETIME DEFAULT CURRENT_TIMESTAMP,
      last_modified        DATETIME DEFAULT CURRENT_TIMESTAMP,
      FOREIGN KEY (user_id) REFERENCES users(id) ON DELETE CASCADE
    )
  ")

  # Idempotent schema upgrade: add new columns on existing databases
  existing_cols <- tryCatch(
    dbGetQuery(con, "PRAGMA table_info(projects)")$name,
    error = function(e) character(0)
  )
  for (col in c("biological_question", "experiment_design", "project_background",
                "ai_api_key", "ai_model", "ai_custom_persona",
                "ai_project_summary")) {
    if (!(col %in% existing_cols)) {
      tryCatch(
        dbExecute(con, sprintf(
          "ALTER TABLE projects ADD COLUMN %s TEXT DEFAULT ''", col)),
        error = function(e) warning("[Auth DB] add column ", col, ": ", e$message)
      )
    }
  }

  # Encrypt any legacy plaintext API keys in place (one-time migration)
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

  message("[Auth DB] Initialised at: ", AUTH_DB_PATH)
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

.db_con <- function() dbConnect(RSQLite::SQLite(), AUTH_DB_PATH)

# Generate a random salt (hex)
.make_salt <- function(n = 32) {
  if (requireNamespace("openssl", quietly = TRUE)) {
    paste(as.character(openssl::rand_bytes(n)), collapse = "")
  } else {
    paste(sample(c(letters, LETTERS, 0:9), n * 2, replace = TRUE), collapse = "")
  }
}

# Hash password with bcrypt (required). Legacy "sha256:<salt>:<hash>" values are
# still accepted by verify_password() and upgraded to bcrypt on next login.
hash_password <- function(password) {
  if (!requireNamespace("bcrypt", quietly = TRUE))
    stop("bcrypt package is required for secure password hashing. Install with: install.packages('bcrypt')")
  bcrypt::hashpw(password)
}

# Verify password against stored hash
verify_password <- function(password, hash) {
  if (is.na(hash) || nchar(hash) == 0) return(FALSE)
  if (startsWith(hash, "$2")) {           # bcrypt
    if (requireNamespace("bcrypt", quietly = TRUE))
      return(bcrypt::checkpw(password, hash))
    return(FALSE)
  }
  if (startsWith(hash, "sha256:")) {      # our fallback
    parts <- strsplit(hash, ":", fixed = TRUE)[[1]]
    if (length(parts) != 3) return(FALSE)
    computed <- digest::digest(paste0(parts[2], password), algo = "sha256", serialize = FALSE)
    return(identical(computed, parts[3]))
  }
  FALSE
}

# Generate a secure random token (64-hex chars)
.generate_token <- function() {
  if (requireNamespace("openssl", quietly = TRUE)) {
    paste(as.character(openssl::rand_bytes(32)), collapse = "")
  } else {
    digest::digest(paste0(Sys.time(), runif(1), sample(.Machine$integer.max, 1)),
                   algo = "sha256", serialize = FALSE)
  }
}

# ---------------------------------------------------------------------------
# User management
# ---------------------------------------------------------------------------

#' Register a new user
#' @return list(success, message, user_id)
create_user <- function(username, email, password, full_name = "") {
  username <- tolower(trimws(username))
  email    <- tolower(trimws(email))

  if (nchar(username) < 3)
    return(list(success = FALSE, message = "Username must be at least 3 characters."))
  if (!grepl("^[a-zA-Z0-9._-]+$", username))
    return(list(success = FALSE, message = "Username may only contain letters, digits, ., _, -."))
  if (!grepl("^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$", email))
    return(list(success = FALSE, message = "Please enter a valid email address."))
  if (nchar(password) < 8)
    return(list(success = FALSE, message = "Password must be at least 8 characters."))

  con <- .db_con(); on.exit(dbDisconnect(con))

  if (nrow(dbGetQuery(con, "SELECT id FROM users WHERE username=?", list(username))) > 0)
    return(list(success = FALSE, message = "Username is already taken."))
  if (nrow(dbGetQuery(con, "SELECT id FROM users WHERE email=?", list(email))) > 0)
    return(list(success = FALSE, message = "Email is already registered."))

  tryCatch({
    dbExecute(con,
      "INSERT INTO users (username, email, password_hash, full_name) VALUES (?,?,?,?)",
      list(username, email, hash_password(password), trimws(full_name))
    )
    uid <- dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id
    list(success = TRUE, message = "Account created!", user_id = uid)
  }, error = function(e) list(success = FALSE, message = paste("Registration failed:", e$message)))
}

#' Authenticate a user (login)
#' @return list(success, message, user_id, username, email, full_name, token)
authenticate_user <- function(login_input, password) {
  login_input <- tolower(trimws(login_input))
  con <- .db_con(); on.exit(dbDisconnect(con))

  MAX_ATTEMPTS <- 5L; LOCK_MINUTES <- 15L

  # --- lockout check ---
  la <- dbGetQuery(con,
    "SELECT attempts, locked_until FROM login_attempts WHERE identifier=?",
    list(login_input))
  if (nrow(la) == 1 && !is.na(la$locked_until[1]) &&
      as.POSIXct(la$locked_until[1], tz = "UTC") > as.POSIXct(Sys.time(), tz = "UTC")) {
    return(list(success = FALSE,
                message = paste0("Too many failed attempts. Try again after ",
                                 format(as.POSIXct(la$locked_until[1]), "%H:%M"), ".")))
  }

  user <- dbGetQuery(con,
    "SELECT id, username, email, password_hash, full_name, is_active
     FROM users WHERE username=? OR email=?",
    list(login_input, login_input)
  )

  # Record a failed attempt; lock the identifier after MAX_ATTEMPTS
  record_fail <- function() {
    now <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    cur <- if (nrow(la) == 1) la$attempts[1] else 0L
    n   <- cur + 1L
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

  # success: clear failed attempts and upgrade any legacy sha256 hash to bcrypt
  dbExecute(con, "DELETE FROM login_attempts WHERE identifier=?", list(login_input))
  if (startsWith(user$password_hash[1], "sha256:")) {
    tryCatch(dbExecute(con, "UPDATE users SET password_hash=? WHERE id=?",
                       list(hash_password(password), user$id[1])),
             error = function(e) NULL)
  }
  dbExecute(con, "UPDATE users SET last_login=CURRENT_TIMESTAMP WHERE id=?", list(user$id[1]))

  token      <- .generate_token()
  expires_at <- format(Sys.time() + 7 * 86400, "%Y-%m-%d %H:%M:%S")
  dbExecute(con,
    "INSERT INTO auth_sessions (user_id, token, expires_at) VALUES (?,?,?)",
    list(user$id[1], token, expires_at)
  )

  list(
    success   = TRUE,
    message   = "Login successful",
    user_id   = user$id[1],
    username  = user$username[1],
    email     = user$email[1],
    full_name = user$full_name[1] %||% user$username[1],
    token     = token
  )
}

#' Change a user's password (verifies the current one first).
#' @return list(success, message)
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

#' Validate a session token
#' @return list(valid, user_id, username, email, full_name)
validate_session_token <- function(token) {
  if (is.null(token) || nchar(token) == 0) return(list(valid = FALSE))
  con <- .db_con(); on.exit(dbDisconnect(con))
  r <- dbGetQuery(con,
    "SELECT s.user_id, u.username, u.email, u.full_name
     FROM auth_sessions s JOIN users u ON s.user_id=u.id
     WHERE s.token=? AND s.is_valid=1
       AND s.expires_at > CURRENT_TIMESTAMP AND u.is_active=1",
    list(token)
  )
  if (nrow(r) == 0) return(list(valid = FALSE))
  list(valid = TRUE, user_id = r$user_id[1], username = r$username[1],
       email = r$email[1], full_name = r$full_name[1])
}

#' Invalidate (logout) a session token
invalidate_session_token <- function(token) {
  if (is.null(token)) return(invisible(FALSE))
  con <- .db_con(); on.exit(dbDisconnect(con))
  dbExecute(con, "UPDATE auth_sessions SET is_valid=0 WHERE token=?", list(token))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Project management
# ---------------------------------------------------------------------------

#' Create a new project
#' @return list(success, message, project_id)
create_project <- function(user_id, name, description = "",
                           species = "", tissue = "", conditions = "",
                           biological_question = "",
                           experiment_design = "",
                           project_background = "",
                           ai_api_key = "", ai_model = "",
                           ai_custom_persona = "",
                           ai_project_summary = "") {
  name <- trimws(name)
  if (nchar(name) == 0)
    return(list(success = FALSE, message = "Project name cannot be empty."))

  con <- .db_con(); on.exit(dbDisconnect(con))

  if (nrow(dbGetQuery(con, "SELECT id FROM projects WHERE user_id=? AND name=?",
                      list(user_id, name))) > 0)
    return(list(success = FALSE, message = "A project with that name already exists."))

  # Encrypt the API key at rest (no-op if the crypto layer is unavailable)
  ai_api_key <- if (exists("encrypt_secret")) encrypt_secret(ai_api_key) else ai_api_key

  tryCatch({
    dbExecute(con,
      "INSERT INTO projects (user_id, name, description,
                             species, tissue, conditions,
                             biological_question, experiment_design,
                             project_background,
                             ai_api_key, ai_model,
                             ai_custom_persona, ai_project_summary)
       VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
      list(user_id, name, description,
           species, tissue, conditions,
           biological_question, experiment_design, project_background,
           ai_api_key, ai_model, ai_custom_persona, ai_project_summary)
    )
    pid <- dbGetQuery(con, "SELECT last_insert_rowid() AS id")$id
    proj_dir <- file.path(SESSIONS_DIR, as.character(user_id), as.character(pid))
    if (!dir.exists(proj_dir)) dir.create(proj_dir, recursive = TRUE)
    list(success = TRUE, message = "Project created.", project_id = pid)
  }, error = function(e) list(success = FALSE, message = paste("Failed:", e$message)))
}

#' List all projects for a user (ordered by last_modified DESC)
get_user_projects <- function(user_id) {
  con <- .db_con(); on.exit(dbDisconnect(con))
  dbGetQuery(con,
    "SELECT id, name, description, species, tissue, conditions,
            status, rds_path, analysis_flags, created_at, last_modified
     FROM projects WHERE user_id=? ORDER BY last_modified DESC",
    list(user_id)
  )
}

#' Get a single project (with ownership check)
get_project <- function(project_id, user_id) {
  con <- .db_con(); on.exit(dbDisconnect(con))
  df <- dbGetQuery(con, "SELECT * FROM projects WHERE id=? AND user_id=?",
                   list(project_id, user_id))
  if (nrow(df) > 0 && "ai_api_key" %in% names(df) && exists("decrypt_secret"))
    df$ai_api_key <- vapply(df$ai_api_key, decrypt_secret, character(1), USE.NAMES = FALSE)
  df
}

#' Update project metadata (dynamic field set)
update_project_metadata <- function(project_id, user_id, ...) {
  allowed <- c("name","description","species","tissue","conditions",
               "biological_question","experiment_design","project_background",
               "ai_api_key","ai_model","ai_custom_persona","ai_project_summary",
               "status","rds_path","analysis_flags")
  params <- Filter(function(x) !is.null(x), list(...))
  params <- params[names(params) %in% allowed]
  if (length(params) == 0) return(invisible(FALSE))

  # Encrypt the API key at rest when it is being updated
  if (!is.null(params$ai_api_key) && exists("encrypt_secret"))
    params$ai_api_key <- encrypt_secret(params$ai_api_key)

  con <- .db_con(); on.exit(dbDisconnect(con))
  set_sql <- paste(paste(names(params), "=?"), collapse = ", ")
  sql <- paste0("UPDATE projects SET ", set_sql,
                ", last_modified=CURRENT_TIMESTAMP WHERE id=? AND user_id=?")
  tryCatch({
    dbExecute(con, sql, c(unname(params), list(project_id, user_id)))
    invisible(TRUE)
  }, error = function(e) { warning("update_project_metadata: ", e$message); invisible(FALSE) })
}

#' Save full analysis state to per-project RDS
save_project_state <- function(project_id, user_id, session_data) {
  rds_dir <- file.path(SESSIONS_DIR, as.character(user_id), as.character(project_id))
  if (!dir.exists(rds_dir)) dir.create(rds_dir, recursive = TRUE)
  rds_path <- file.path(rds_dir, "analysis.rds")

  tryCatch({
    saveRDS(session_data, rds_path, compress = "gzip")

    # Infer status from flags
    flags  <- session_data$analysis_flags %||% list()
    status <- if (isTRUE(flags$annotation_done))   "annotated"
         else if (isTRUE(flags$clustering_done))   "clustered"
         else if (isTRUE(flags$dimred_done))        "reduced"
         else if (isTRUE(flags$qc_done))            "qc_done"
         else if (isTRUE(flags$data_loaded))        "data_loaded"
         else "empty"

    flags_json <- tryCatch(
      jsonlite::toJSON(flags, auto_unbox = TRUE),
      error = function(e) "{}"
    )
    update_project_metadata(project_id, user_id,
      rds_path       = rds_path,
      status         = status,
      analysis_flags = as.character(flags_json)
    )
    list(success = TRUE, path = rds_path)
  }, error = function(e) list(success = FALSE, message = e$message))
}

#' Load analysis state from per-project RDS
load_project_state <- function(project_id, user_id) {
  proj <- get_project(project_id, user_id)
  if (nrow(proj) == 0) return(NULL)
  rds_path <- proj$rds_path[1]
  if (is.na(rds_path) || !file.exists(rds_path)) return(NULL)
  tryCatch(readRDS(rds_path), error = function(e) {
    warning("load_project_state: ", e$message); NULL
  })
}

#' Delete a project and its files
delete_project <- function(project_id, user_id) {
  if (nrow(get_project(project_id, user_id)) == 0)
    return(list(success = FALSE, message = "Project not found."))
  proj_dir <- file.path(SESSIONS_DIR, as.character(user_id), as.character(project_id))
  if (dir.exists(proj_dir)) unlink(proj_dir, recursive = TRUE)
  con <- .db_con(); on.exit(dbDisconnect(con))
  dbExecute(con, "DELETE FROM projects WHERE id=? AND user_id=?", list(project_id, user_id))
  list(success = TRUE, message = "Project deleted.")
}

# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------

#' Map status code → label + colour
get_status_info <- function(status) {
  switch(status %||% "empty",
    "empty"       = list(label = "New",         color = "#6c757d", pct = 0),
    "data_loaded" = list(label = "Data Loaded",  color = "#17a2b8", pct = 15),
    "qc_done"     = list(label = "QC Complete",  color = "#007bff", pct = 30),
    "reduced"     = list(label = "Dim Reduced",  color = "#6f42c1", pct = 50),
    "clustered"   = list(label = "Clustered",    color = "#fd7e14", pct = 70),
    "annotated"   = list(label = "Annotated",    color = "#28a745", pct = 90),
                    list(label = "In Progress",  color = "#ffc107", pct = 20)
  )
}
