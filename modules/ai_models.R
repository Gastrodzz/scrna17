# modules/ai_models.R
# Central registry for Claude models. Live-fetches from Anthropic's /v1/models
# endpoint when an API key is available, caches the result, and falls back to a
# curated static list otherwise. Single source of truth for every model dropdown
# and default-model resolution in the app.

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

# --- Curated fallback (newest -> oldest); used when no key / offline / API error.
# NOTE: this list is only shown before an API key is entered — once a key is
# present, list_claude_models() live-fetches the authoritative list from
# Anthropic's /v1/models endpoint. Keep these current so the pre-key view is sane.
AI_MODELS_FALLBACK <- c(
  "Claude Opus 4.8 (Most capable)" = "claude-opus-4-8",
  "Claude Sonnet 5 (Recommended)"  = "claude-sonnet-5",
  "Claude Haiku 4.5 (Fastest)"     = "claude-haiku-4-5-20251001"
)
AI_MODELS_DEFAULT <- "claude-sonnet-5"

# Families that support extended thinking (used to synthesize "-thinking" variants)
.AI_THINKING_FAMILIES <- c("opus", "sonnet")

# Module-level cache
.ai_models_cache <- new.env(parent = emptyenv())
.ai_models_cache$data <- NULL
.ai_models_cache$ts   <- 0

#' List available Claude models as a named vector (display_name -> id).
#' @param api_key Anthropic key; when supplied, live-fetches from /v1/models.
#' @param force Bypass the cache.
#' @param ttl Cache lifetime in seconds.
#' @param include_thinking Append synthesized "-thinking" display variants.
list_claude_models <- function(api_key = NULL, force = FALSE, ttl = 3600,
                               include_thinking = TRUE) {
  now <- as.numeric(Sys.time())
  if (!force && !is.null(.ai_models_cache$data) &&
      (now - .ai_models_cache$ts) < ttl) {
    return(.ai_models_cache$data)
  }

  live <- NULL
  if (!is.null(api_key) && nzchar(trimws(api_key)) &&
      requireNamespace("httr2", quietly = TRUE)) {
    live <- tryCatch({
      resp <- httr2::request("https://api.anthropic.com/v1/models") |>
        httr2::req_headers(
          "x-api-key"         = trimws(api_key),
          "anthropic-version" = "2023-06-01"
        ) |>
        httr2::req_url_query(limit = 100) |>
        httr2::req_timeout(20) |>
        httr2::req_perform()
      body <- httr2::resp_body_json(resp)
      data <- body$data
      if (is.null(data) || length(data) == 0) stop("empty model list")
      ids  <- vapply(data, function(m) m$id %||% NA_character_, character(1))
      disp <- vapply(data, function(m) m$display_name %||% m$id %||% NA_character_, character(1))
      keep <- !is.na(ids) & nzchar(ids) & grepl("^claude-", ids)
      ids <- ids[keep]; disp <- disp[keep]
      if (length(ids) == 0) stop("no claude models returned")
      stats::setNames(ids, disp)
    }, error = function(e) NULL)
  }

  base <- if (!is.null(live)) live else AI_MODELS_FALLBACK

  if (isTRUE(include_thinking)) {
    thinking <- character(0)
    for (i in seq_along(base)) {
      id <- base[[i]]; nm <- names(base)[i]
      if (any(vapply(.AI_THINKING_FAMILIES,
                     function(f) grepl(f, id, ignore.case = TRUE), logical(1)))) {
        thinking[paste0(nm, " — Extended Thinking")] <- paste0(id, "-thinking")
      }
    }
    base <- c(base, thinking)
  }

  # Cache live results for the full ttl; cache fallback only briefly so that a
  # later key entry re-triggers a live fetch soon.
  .ai_models_cache$data <- base
  .ai_models_cache$ts   <- if (!is.null(live)) now else (now - ttl + 60)
  base
}

#' Named choices vector for selectInput (display_name -> id).
ai_model_choices <- function(api_key = NULL, force = FALSE) {
  list_claude_models(api_key = api_key, force = force)
}

#' Resolve a sensible default model id (non-thinking sonnet preferred).
get_default_model <- function(api_key = NULL) {
  ch <- tryCatch(list_claude_models(api_key = api_key, include_thinking = FALSE),
                 error = function(e) AI_MODELS_FALLBACK)
  ids <- unname(ch)
  sonnet <- ids[grepl("sonnet", ids, ignore.case = TRUE) & !grepl("-thinking$", ids)]
  if (length(sonnet) > 0) return(sonnet[1])
  if (length(ids) > 0) return(ids[1])
  AI_MODELS_DEFAULT
}
