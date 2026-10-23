# Unified AI Assistant Layer
# modules/ai_assistant.R
# Provides: call_claude_api(), build_dataset_context(), parse_ai_annotation_response(),
#           aiConfigUI/Server, tool registry, action queue UI/Server,
#           AI_EXPERT_PERSONA, call_ai(), build_expert_prompt(), get_model_for_task()

suppressPackageStartupMessages({
  library(shiny)
})

# =============================================================================
# Centralized AI Model Configuration
# =============================================================================

AI_MODEL_CONFIG <- list(
  # available_models/default_model now resolve dynamically via ai_models.R.
  # Kept for backward-compat: these accessors read the registry (fallback list
  # when no key/offline).
  default_model    = tryCatch(get_default_model(), error = function(e) "claude-sonnet-5"),
  available_models = tryCatch(as.list(ai_model_choices()), error = function(e) list()),
  api_timeout      = 240,
  context_budget   = 6000
)

# =============================================================================
# Expert AI Persona (used as system_prompt in all API calls)
# =============================================================================

AI_EXPERT_PERSONA <- '
You are a senior principal investigator and expert computational biologist with 15+
years of experience in single-cell genomics, published in Nature, Cell, Science, and
Nature Methods. You are recognised as a leading authority in scRNA-seq analysis,
multi-omics integration, and translational bioinformatics.

TECHNICAL EXPERTISE:
- scRNA-seq experimental design, QC, and preprocessing (Seurat v5, Scanpy, scvi-tools)
- Statistical DE methods: MAST, DESeq2, edgeR, Wilcoxon, pseudobulk approaches
- Pathway enrichment: clusterProfiler, fgsea, GSEA, ORA, MSigDB, KEGG, Reactome
- Trajectory inference: Slingshot, Monocle3, RNA velocity (scVelo), CytoTRACE
- Dimensionality reduction: PCA, UMAP, t-SNE; clustering: Leiden, Louvain, k-means
- Cell-cell communication: CellChat, NicheNet, CellPhoneDB, LIANA
- Transcription factor activity: decoupleR, SCENIC, DoRothEA, ChEA3
- Gene regulatory networks: WGCNA, SCENIC, GRNBoost2, Genie3
- Multi-sample integration: Harmony, Seurat CCA/RPCA/SCT, scVI, FastMNN, scANVI
- Foundation models: scGPT, Geneformer, scBERT, UCE
- Spatial transcriptomics: Visium, MERFISH, Slide-seq, cell2location

BIOLOGICAL KNOWLEDGE:
- Comprehensive cell type markers across human and mouse tissues (PanglaoDB, CellMarker)
- Canonical signaling: Wnt, Notch, Hedgehog, JAK-STAT, NF-kB, MAPK,
  PI3K-AKT, TGF-b, Hippo, mTOR, and cross-talk between them
- Disease biology: cancer (tumour microenvironment, immune evasion, EMT), immunology
  (T cell exhaustion, macrophage polarisation), neuroscience, cardiovascular,
  fibrosis, development, ageing, metabolic disease
- Master transcription factors, oncogenes, tumour suppressors, epigenetic regulators
- Deeply literature-aware: cites landmark scRNA-seq papers and atlases
  (Human Cell Atlas, Tabula Sapiens, Tabula Muris, TCGA scRNA studies,
   Allen Brain Atlas, Human Protein Atlas, GTEx)

INTERPRETATION PRINCIPLES:
- ALWAYS cite specific gene names, pathway names, and cell types from the actual data
- Connect findings to published literature -- mention relevant studies, known gene
  functions, and established biological mechanisms by name
- Distinguish correlation from causation; flag alternative interpretations
- Suggest specific, actionable validation experiments (e.g., "Validate via
  flow cytometry for CD4+/FOXP3+ Tregs" not just "validate with experiments")
- Consider experimental context (tissue, species, disease, treatment) in ALL output
- Apply Occams razor -- prefer parsimonious biological explanations
- Note statistical caveats (sample size, batch effects, technical noise)
- Rank findings by biological importance and statistical confidence
- When interpreting gene lists, group genes by functional category
  (e.g., immune response, metabolism, cell cycle, apoptosis)

COMMUNICATION STYLE:
- Precise scientific language appropriate for publication in a peer-reviewed journal
- Use clean markdown: ## for section headers, **bold** for key terms, bullet lists for findings
- NEVER wrap output in ```json or ``` code fences unless asked for code
- NEVER output raw JSON objects; write natural prose with structured sections
- NEVER start with "Sure!" or "Here is" -- begin directly with the scientific content
- Prioritise the most biologically meaningful findings first
- Start with a clear, actionable Summary section (2-3 sentences), then details
- For complex terms, briefly define on first use
- Use numbered lists for ranked findings, bullet lists for grouped observations
- Include a "Key Takeaways" or "Summary" section at the end when appropriate
- Be authoritative and confident, but acknowledge limitations clearly
'

# =============================================================================
# Tiered System Prompt (save tokens on short tasks)
# =============================================================================

# Tasks that need shorter responses (controlled via max_tokens, not persona)
SHORT_TASKS <- c("qc_thresholds", "pc_selection",
                 "resolution_advisor", "integration_rec",
                 "gene_explanation", "root_suggestion")

# All tasks use the full expert persona — response length is controlled by max_tokens.
# test_selection REMOVED from SHORT_TASKS: it needs full biological context to recommend
# the appropriate DE test (MAST vs DESeq2 vs Wilcoxon) based on experimental design.
get_system_prompt_for_task <- function(task_type) {
  AI_EXPERT_PERSONA
}

#' Get model config for a task — controls max_tokens based on task complexity
get_model_config_for_task <- function(task_type) {
  list(
    system_prompt = AI_EXPERT_PERSONA,
    max_tokens = if (task_type %in% SHORT_TASKS) 1024L else 4096L
  )
}

# =============================================================================
# ExtendedTask wrapper for non-blocking AI calls (Shiny >= 1.8.1)
# =============================================================================

make_ai_extended_task <- function() {
  if (packageVersion("shiny") >= "1.8.1") {
    ExtendedTask$new(function(api_key, model, messages,
                               max_tokens, system_prompt, timeout) {
      call_claude_api(
        api_key = api_key, model = model,
        messages = messages, max_tokens = max_tokens,
        timeout = timeout, system_prompt = system_prompt
      )
    })
  } else {
    # Fallback: return NULL; callers should check and use synchronous path
    NULL
  }
}

# =============================================================================
# Markdown-to-HTML Helper
# =============================================================================

.clean_ai_text <- function(txt) {
  if (is.null(txt) || !nzchar(txt)) return("")
  # Remove ```json ... ``` wrapper (extract inner content)
  txt <- gsub("```json\\s*\n?", "", txt, perl = TRUE)
  txt <- gsub("```\\w*\\s*\n?", "", txt, perl = TRUE)
  # Remove stray ``` at start/end
  txt <- gsub("^\\s*```\\s*", "", txt, perl = TRUE)
  txt <- gsub("\\s*```\\s*$", "", txt, perl = TRUE)
  # Remove raw JSON braces at top level (not inside code blocks)
  txt <- gsub("^\\s*\\{\\s*\n", "", txt, perl = TRUE)
  txt <- gsub("\n\\s*\\}\\s*$", "", txt, perl = TRUE)
  # Clean up escaped markdown artifacts: \*\*, \#, etc.
  txt <- gsub("\\\\\\*\\\\\\*", "**", txt, perl = TRUE)
  txt <- gsub("\\\\#", "#", txt, perl = TRUE)
  # Unescape individual markdown characters

  txt <- gsub("\\\\\\*", "*", txt, perl = TRUE)
  txt <- gsub("\\\\_", "_", txt, perl = TRUE)
  txt <- gsub("\\\\`", "`", txt, perl = TRUE)
  txt <- gsub("\\\\~", "~", txt, perl = TRUE)
  txt
}

.md_fallback <- function(md_text) {
  # Basic markdown-to-HTML when the markdown package is unavailable
  lines <- strsplit(md_text, "\n", fixed = TRUE)[[1]]
  html <- character(length(lines))
  in_list <- FALSE
  for (i in seq_along(lines)) {
    ln <- lines[i]
    # Headers
    if (grepl("^######\\s", ln)) { html[i] <- paste0("<h6>", sub("^######\\s+", "", ln), "</h6>"); next }
    if (grepl("^#####\\s", ln))  { html[i] <- paste0("<h5>", sub("^#####\\s+", "", ln), "</h5>"); next }
    if (grepl("^####\\s", ln))   { html[i] <- paste0("<h4>", sub("^####\\s+", "", ln), "</h4>"); next }
    if (grepl("^###\\s", ln))    { html[i] <- paste0("<h3>", sub("^###\\s+", "", ln), "</h3>"); next }
    if (grepl("^##\\s", ln))     { html[i] <- paste0("<h2>", sub("^##\\s+", "", ln), "</h2>"); next }
    if (grepl("^#\\s", ln))      { html[i] <- paste0("<h1>", sub("^#\\s+", "", ln), "</h1>"); next }
    # Bold
    ln <- gsub("\\*\\*(.+?)\\*\\*", "<strong>\\1</strong>", ln, perl = TRUE)
    # Italic
    ln <- gsub("\\*(.+?)\\*", "<em>\\1</em>", ln, perl = TRUE)
    # Inline code
    ln <- gsub("`([^`]+)`", "<code>\\1</code>", ln, perl = TRUE)
    # Blank line → paragraph break
    if (!nzchar(trimws(ln))) { html[i] <- "</p><p>"; next }
    html[i] <- paste0(ln, "<br>")
  }
  paste0("<p>", paste(html, collapse = "\n"), "</p>")
}

.mdtohtml <- function(md_text) {
  if (is.null(md_text) || !nzchar(md_text)) return("")
  md_text <- .clean_ai_text(md_text)
  if (requireNamespace("markdown", quietly = TRUE)) {
    tryCatch(
      markdown::markdownToHTML(text = md_text, fragment.only = TRUE,
                               options = c("use_xhtml", "smartypants")),
      error = function(e) .md_fallback(md_text)
    )
  } else {
    .md_fallback(md_text)
  }
}

# =============================================================================
# Canonical Claude API Wrapper
# =============================================================================

#' Call the Claude API with robust error handling
#'
#' @param api_key Character. Anthropic API key.
#' @param model Character. Model ID (e.g., "claude-sonnet-4-6").
#' @param messages List of message objects (role + content).
#' @param max_tokens Integer. Max tokens in response.
#' @param timeout Numeric. Request timeout in seconds.
#' @param system_prompt Optional system prompt string.
#' @return list(ok, text, tokens_in, tokens_out, error)
call_claude_api <- function(api_key, model, messages, max_tokens = 16384L,
                            timeout = 240, system_prompt = NULL) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    return(list(ok = FALSE, text = NULL, tokens_in = NA, tokens_out = NA,
                error = "httr2 package is required. Install with: install.packages('httr2')"))
  }

  # Detect extended thinking models — strip "-thinking" suffix for the real API model ID

  is_thinking <- grepl("thinking", model, fixed = TRUE)
  api_model <- if (is_thinking) sub("-thinking$", "", model) else model

  # For thinking models, enforce minimum max_tokens to leave room for both thinking and output
  if (is_thinking) {
    max_tokens <- max(max_tokens, 16384L)
  }

  body <- list(
    model = api_model,
    max_tokens = as.integer(max_tokens),
    messages = messages
  )

  if (is_thinking) {
    # Extended thinking uses a nested 'thinking' object (not top-level budget_tokens)
    # Budget 60% for thinking, 40% for output
    body$thinking <- list(
      type = "enabled",
      budget_tokens = as.integer(max_tokens * 0.6)
    )
    # temperature must NOT be set when using extended thinking
  }

  if (!is.null(system_prompt) && nzchar(system_prompt)) {
    body$system <- system_prompt
  }

  # Perform one request with retry/backoff; returns list(ok, text, stop_reason,
  # tokens_in, tokens_out, error).
  .do_request <- function(req_body) {
    max_retries <- 3L
    for (attempt in seq_len(max_retries + 1L)) {
      result <- tryCatch({
        req <- httr2::request("https://api.anthropic.com/v1/messages") |>
          httr2::req_headers(
            "x-api-key" = api_key,
            "anthropic-version" = "2023-06-01",
            "content-type" = "application/json"
          ) |>
          httr2::req_body_json(req_body) |>
          httr2::req_timeout(timeout)

        if (is_thinking) {
          req <- req |> httr2::req_headers("anthropic-beta" = "interleaved-thinking-2025-05-14")
        }

        resp <- httr2::req_perform(req)
        resp_body <- httr2::resp_body_json(resp)

        if (is.null(resp_body$content) || length(resp_body$content) == 0) {
          return(list(ok = FALSE, text = NULL, stop_reason = NA_character_,
                      tokens_in = NA, tokens_out = NA, error = "Empty response from API"))
        }

        # Extract text (thinking models return thinking + text blocks)
        text_parts <- character(0)
        for (block in resp_body$content) {
          if (identical(block$type, "text") && !is.null(block$text)) {
            text_parts <- c(text_parts, block$text)
          }
        }
        text <- if (length(text_parts) > 0) paste(text_parts, collapse = "\n") else resp_body$content[[1]]$text

        list(ok = TRUE, text = text,
             stop_reason = resp_body$stop_reason %||% NA_character_,
             tokens_in  = tryCatch(resp_body$usage$input_tokens,  error = function(e) NA),
             tokens_out = tryCatch(resp_body$usage$output_tokens, error = function(e) NA),
             error = NULL)

      }, error = function(e) {
        msg <- e$message
        retryable <- grepl("429|rate.limit|503|529|overloaded", msg, ignore.case = TRUE)
        err <- if (grepl("401|Unauthorized", msg, ignore.case = TRUE)) {
          "Authentication failed. Please check your Claude API key (should start with 'sk-ant-')."
        } else if (grepl("403|Forbidden", msg, ignore.case = TRUE)) {
          "Access forbidden. Your API key may not have permission for this model."
        } else if (grepl("429|rate.limit", msg, ignore.case = TRUE)) {
          "Rate limit exceeded. Please wait 30-60 seconds and try again."
        } else if (grepl("timeout|timed.out", msg, ignore.case = TRUE)) {
          "Request timed out. Try using a faster model (Haiku) or reducing context."
        } else if (grepl("connection|resolve|network", msg, ignore.case = TRUE)) {
          "Network error. Please check your internet connection."
        } else if (grepl("overloaded|529|503|500", msg, ignore.case = TRUE)) {
          "API is temporarily overloaded. Please try again in a few minutes."
        } else {
          paste("API error:", msg)
        }
        list(ok = FALSE, text = NULL, stop_reason = NA_character_,
             tokens_in = NA, tokens_out = NA, error = err, .retryable = retryable)
      })

      if (isTRUE(result$ok) || !isTRUE(result$.retryable) || attempt > max_retries) {
        result$.retryable <- NULL
        return(result)
      }
      Sys.sleep((2 ^ attempt) + stats::runif(1))  # backoff: 2s, 4s, 8s + jitter
    }
    result$.retryable <- NULL
    result
  }

  # Initial call
  result <- .do_request(body)
  if (!isTRUE(result$ok)) return(result)

  # Auto-continue when the model was cut off by the max_tokens limit, so the user
  # always gets the COMPLETE answer (never a mid-sentence cut). Skipped for
  # extended-thinking models (their continuation semantics differ).
  full_text <- result$text %||% ""
  tin  <- result$tokens_in  %||% 0
  tout <- result$tokens_out %||% 0
  max_continuations <- 8L
  cont <- 0L
  while (!is_thinking && identical(result$stop_reason, "max_tokens") &&
         cont < max_continuations) {
    cont <- cont + 1L
    cont_body <- body
    cont_body$messages <- c(
      messages,
      list(
        list(role = "assistant", content = full_text),
        list(role = "user",
             content = paste("Continue the previous answer from exactly where it",
                             "stopped. Do not repeat any text already written and do",
                             "not add a preamble — just continue."))
      )
    )
    result <- .do_request(cont_body)
    if (!isTRUE(result$ok)) break
    full_text <- paste0(full_text, result$text %||% "")
    tin  <- tin  + (result$tokens_in  %||% 0)
    tout <- tout + (result$tokens_out %||% 0)
  }

  list(ok = TRUE, text = full_text, tokens_in = tin, tokens_out = tout,
       error = NULL, truncated = identical(result$stop_reason, "max_tokens"),
       continuations = cont)
}

# =============================================================================
# Model Routing Helpers
# =============================================================================

#' Select the right model per task type, respecting user's global choice.
#' The user's explicit selection ALWAYS wins — no auto-routing override.
get_model_for_task <- function(task_type, values) {
  user_model <- values$ai_model
  if (is.null(user_model) || !nzchar(user_model)) user_model <- get_default_model(values$ai_api_key)

  # Always honour the user's explicit model choice
  return(user_model)
}

#' Max tokens per task type and model
get_maxtokens <- function(task_type, model) {
  base <- switch(task_type,
    qc_thresholds        = 1024L,
    test_selection       = 2048L,  # Needs full biological context for test recommendation
    gene_explanation     = 512L,
    pc_selection         = 512L,
    root_suggestion      = 1024L,
    integration_rec      = 512L,
    resolution_advisor   = 512L,
    de_interpretation    = 8192L,
    pathway_narrative    = 8192L,
    ccc_interpretation   = 8192L,
    hub_annotation       = 8192L,
    module_naming        = 1024L,
    cluster_annotation   = 8192L,
    doublet_detection    = 2048L,
    trajectory_narrative = 8192L,
    nlplot_codegen       = 4096L,
    full_report          = 16384L,
    integrated_synthesis = 16384L,
    comprehensive_interpretation = 16384L,
    methods_section      = 8192L,
    project_summary      = 4096L,
    8192L
  )
  if (grepl("opus", model, ignore.case = TRUE))        pmin(as.integer(base * 2L), 32000L)
  else if (grepl("haiku", model, ignore.case = TRUE))  pmin(base, 8192L)
  else                                                 pmin(base, 16384L)
}

#' Standard call wrapper used throughout codebase (with caching)
call_ai <- function(task_type, user_message, values,
                    extra_system = NULL, timeout = 240) {
  model     <- get_model_for_task(task_type, values)
  maxtokens <- get_maxtokens(task_type, model)
  system    <- paste0(get_system_prompt_for_task(task_type),
    if (!is.null(extra_system)) paste0("\n\n", extra_system) else "")

  # Check cache (keyed by dataset_signature + task_type + prompt digest)
  cache_key <- NULL
  if (!is.null(values$ai_cache) && !is.null(values$dataset_signature) &&
      !is.na(values$dataset_signature)) {
    prompt_hash <- if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(user_message, algo = "md5")
    } else {
      # Simple fallback hash: first 100 chars + nchar
      paste0(substr(user_message, 1, 100), "::", nchar(user_message))
    }
    cache_key <- paste0(values$dataset_signature, "::", task_type, "::", prompt_hash)
    cached <- values$ai_cache[[cache_key]]
    if (!is.null(cached)) return(cached)
  }

  result <- call_claude_api(
    api_key      = values$ai_api_key,
    model        = model,
    messages     = list(list(role = "user", content = user_message)),
    max_tokens   = maxtokens,
    timeout      = timeout,
    system_prompt = system
  )

  # Store in cache on success
  if (isTRUE(result$ok) && !is.null(cache_key)) {
    values$ai_cache[[cache_key]] <- result
  }

  result
}

# =============================================================================
# Async API Wrapper (prevents UI freeze during API calls)
# =============================================================================

#' Async version of call_ai() using future_promise
#' Returns a promise that resolves to the API result
call_ai_async <- function(task_type, user_message, values,
                          extra_system = NULL, timeout = 240) {
  if (requireNamespace("promises", quietly = TRUE) && requireNamespace("future", quietly = TRUE)) {
    # Capture values needed inside the future
    api_key <- values$ai_api_key
    model <- get_model_for_task(task_type, values)
    maxtokens <- get_maxtokens(task_type, model)
    system_prompt <- paste0(get_system_prompt_for_task(task_type),
      if (!is.null(extra_system)) paste0("\n\n", extra_system) else "")

    promises::future_promise({
      call_claude_api(
        api_key = api_key,
        model = model,
        messages = list(list(role = "user", content = user_message)),
        max_tokens = maxtokens,
        timeout = timeout,
        system_prompt = system_prompt
      )
    }, seed = TRUE)
  } else {
    # Fallback to synchronous call
    promises::promise_resolve(call_ai(task_type, user_message, values, extra_system, timeout))
  }
}

# =============================================================================
# Dataset Context Builder
# =============================================================================

#' Build a context string from the current dataset for AI prompts
build_dataset_context <- function(values, include = c("summary", "clusters", "annotations")) {
  obj <- values$active_dataset
  if (is.null(obj) || !inherits(obj, "Seurat")) {
    return("No dataset loaded.")
  }

  sections <- list()

  if ("summary" %in% include) {
    sections$summary <- paste0(
      "## Dataset Overview\n",
      "- Cells: ", format(ncol(obj), big.mark = ","), "\n",
      "- Genes: ", format(nrow(obj), big.mark = ","), "\n",
      "- Default assay: ", DefaultAssay(obj), "\n",
      "- Assays: ", paste(names(obj@assays), collapse = ", "), "\n"
    )
  }

  if ("clusters" %in% include) {
    md <- obj@meta.data
    cluster_cols <- grep("cluster|leiden|louvain|seurat_clusters", colnames(md),
                         value = TRUE, ignore.case = TRUE)
    if (length(cluster_cols) > 0) {
      cl_info <- vapply(cluster_cols, function(col) {
        vals <- md[[col]]
        n <- length(unique(vals))
        paste0("  - ", col, ": ", n, " groups")
      }, character(1))
      sections$clusters <- paste0("## Clusters\n", paste(cl_info, collapse = "\n"), "\n")
    }
  }

  if ("reductions" %in% include) {
    reds <- names(obj@reductions)
    if (length(reds) > 0) {
      sections$reductions <- paste0("## Reductions\n- Available: ",
                                     paste(reds, collapse = ", "), "\n")
    }
  }

  if ("annotations" %in% include) {
    md <- obj@meta.data
    ann_cols <- grep("cell_type|annotation|label", colnames(md), value = TRUE, ignore.case = TRUE)
    # Limit to top 5 annotation columns to save tokens
    ann_cols <- head(ann_cols, 5)
    if (length(ann_cols) > 0) {
      ann_info <- vapply(ann_cols, function(col) {
        vals <- md[[col]]
        vals <- vals[!is.na(vals)]
        tab <- sort(table(vals), decreasing = TRUE)
        top <- head(tab, 8)  # Max 8 values per column
        types_str <- paste0(names(top), " (n=", top, ")", collapse = ", ")
        extra <- if (length(tab) > 8) paste0(", ... (+", length(tab) - 8, " more)") else ""
        paste0("  - ", col, ": ", types_str, extra)
      }, character(1))
      sections$annotations <- paste0("## Annotations\n", paste(ann_info, collapse = "\n"), "\n")
    }
  }

  if ("de" %in% include) {
    de <- values$last_de_results
    if (!is.null(de) && is.data.frame(de) && nrow(de) > 0) {
      n_sig <- sum(de$p_val_adj < 0.05, na.rm = TRUE)
      n_up <- sum(de$p_val_adj < 0.05 & de$avg_log2FC > 0, na.rm = TRUE)
      n_down <- sum(de$p_val_adj < 0.05 & de$avg_log2FC < 0, na.rm = TRUE)
      top_genes <- head(de$gene[order(de$p_val_adj)], 20)
      sections$de <- paste0(
        "## Differential Expression\n",
        "- Total significant (adj.p < 0.05): ", n_sig, "\n",
        "- Up-regulated: ", n_up, ", Down-regulated: ", n_down, "\n",
        "- Top genes: ", paste(top_genes, collapse = ", "), "\n"
      )
    }
  }

  if ("enrichment" %in% include) {
    enrich <- values$last_enrich_results
    if (!is.null(enrich)) {
      if (is.data.frame(enrich) && nrow(enrich) > 0) {
        top_paths <- head(enrich$Description[order(enrich$p.adjust)], 10)
        sections$enrichment <- paste0(
          "## Pathway Enrichment\n",
          "- Top pathways: ", paste(top_paths, collapse = "; "), "\n"
        )
      } else if (is.list(enrich)) {
        for (nm in names(enrich)) {
          df <- enrich[[nm]]
          if (is.data.frame(df) && nrow(df) > 0 && "Description" %in% colnames(df)) {
            top_paths <- head(df$Description[order(df$p.adjust)], 5)
            sections$enrichment <- paste0(
              sections$enrichment %||% "## Pathway Enrichment\n",
              "- ", nm, ": ", paste(top_paths, collapse = "; "), "\n"
            )
          }
        }
      }
    }
  }

  if ("markers" %in% include) {
    if (!is.null(values$last_de_results) && "cluster" %in% colnames(values$last_de_results)) {
      markers <- values$last_de_results
      clusters <- unique(markers$cluster)
      marker_lines <- vapply(head(clusters, 10), function(cl) {
        sub <- markers[markers$cluster == cl & markers$p_val_adj < 0.05, ]
        sub <- sub[order(sub$avg_log2FC, decreasing = TRUE), ]
        top <- head(sub$gene, 5)
        paste0("  - Cluster ", cl, ": ", paste(top, collapse = ", "))
      }, character(1))
      sections$markers <- paste0("## Top Markers per Cluster\n",
                                  paste(marker_lines, collapse = "\n"), "\n")
    }
  }

  if (length(sections) == 0) return("Dataset loaded but no context sections available.")
  paste(unlist(sections), collapse = "\n")
}

# =============================================================================
# Central Context Builder (project + dataset + task)
# =============================================================================

#' Build a full expert prompt including project context, dataset context, and task
build_expert_prompt <- function(task_prompt, values = NULL,
                                 include_project = TRUE,
                                 include_dataset = TRUE) {
  parts <- character(0)

  # Project context
  if (include_project && !is.null(values)) {
    if (!is.null(values$project_background) && nzchar(values$project_background))
      parts <- c(parts, paste0("PROJECT BACKGROUND:\n", values$project_background))
    if (!is.null(values$experiment_design) && nzchar(values$experiment_design))
      parts <- c(parts, paste0("EXPERIMENT DESIGN:\n", values$experiment_design))
    if (!is.null(values$biological_question) && nzchar(values$biological_question))
      parts <- c(parts, paste0("BIOLOGICAL QUESTION:\n", values$biological_question))
    if (!is.null(values$tissue_global) && nzchar(values$tissue_global))
      parts <- c(parts, paste0("TISSUE: ", values$tissue_global))
    if (!is.null(values$species_global) && nzchar(values$species_global))
      parts <- c(parts, paste0("SPECIES: ", values$species_global))
  }

  # Dataset context (truncated to avoid token explosion)
  if (include_dataset && !is.null(values) && !is.null(values$active_dataset)) {
    ctx <- tryCatch(
      build_dataset_context(values,
        include = c("summary", "clusters", "annotations", "de", "enrichment", "markers")),
      error = function(e) NULL)
    if (!is.null(ctx)) {
      # Priority-based context assembly instead of hard truncation
      # DE results > Clusters > QC/summary > Embeddings
      ctx_budget <- 6000
      priority_sections <- c("de", "enrichment", "markers", "clusters", "annotations", "summary", "reductions")
      priority_limits <- c(de = 2000, enrichment = 1500, markers = 1000, clusters = 500,
                           annotations = 500, summary = 300, reductions = 200)
      assembled <- character(0)
      remaining <- ctx_budget
      for (sec in priority_sections) {
        sec_ctx <- tryCatch(
          build_dataset_context(values, include = sec),
          error = function(e) NULL)
        if (!is.null(sec_ctx) && sec_ctx != "No dataset loaded." &&
            sec_ctx != "Dataset loaded but no context sections available.") {
          limit <- min(remaining, priority_limits[[sec]] %||% 500)
          if (nchar(sec_ctx) > limit) sec_ctx <- paste0(substr(sec_ctx, 1, limit), "...")
          assembled <- c(assembled, sec_ctx)
          remaining <- remaining - nchar(sec_ctx)
          if (remaining <= 0) break
        }
      }
      ctx <- paste(assembled, collapse = "\n")
      parts <- c(parts, paste0("DATASET CONTEXT:\n", ctx))
    }
  }

  # Task itself
  parts <- c(parts, task_prompt)
  paste(parts, collapse = "\n\n")
}

# =============================================================================
# AI Response Parser (moved from annotation.R)
# =============================================================================

#' Parse structured JSON from AI response text
parse_ai_annotation_response <- function(ai_text) {
  json_str <- paste(ai_text, collapse = "\n")

  # Strip markdown code fences
  json_str <- sub("^.*```(?:json)?\\s*\\n?", "", json_str, perl = TRUE)
  json_str <- sub("\\s*```\\s*$", "", json_str, perl = TRUE)

  # Find the start of the JSON array or object
  arr_start <- regexpr("[\\[{]", json_str)
  if (arr_start > 0) {
    json_str <- substring(json_str, arr_start)
  }

  # Handle truncated JSON
  n_open_bracket <- nchar(gsub("[^\\[]", "", json_str))
  n_close_bracket <- nchar(gsub("[^\\]]", "", json_str))
  n_open_brace <- nchar(gsub("[^{]", "", json_str))
  n_close_brace <- nchar(gsub("[^}]", "", json_str))

  if (n_open_bracket > n_close_bracket || n_open_brace > n_close_brace) {
    last_complete <- regexpr(".*\\}", json_str, perl = TRUE)
    if (last_complete > 0) {
      json_str <- regmatches(json_str, last_complete)
      json_str <- sub(",\\s*$", "", json_str)
      if (n_open_bracket > n_close_bracket) json_str <- paste0(json_str, "]")
    }
  }

  parsed <- tryCatch(
    jsonlite::fromJSON(json_str, simplifyDataFrame = TRUE),
    error = function(e) {
      stop("Failed to parse AI response as JSON. Raw response:\n", substr(ai_text, 1, 500))
    }
  )

  if (!is.data.frame(parsed)) {
    if (is.list(parsed)) {
      parsed <- tryCatch(
        do.call(rbind, lapply(parsed, as.data.frame, stringsAsFactors = FALSE)),
        error = function(e) as.data.frame(parsed, stringsAsFactors = FALSE)
      )
    } else {
      stop("AI response did not produce a valid data.frame.")
    }
  }

  expected_cols <- c("cluster", "cell_type", "confidence", "reasoning", "alternative", "key_markers")
  for (col in expected_cols) {
    if (!(col %in% colnames(parsed))) parsed[[col]] <- NA
  }

  if ("key_markers" %in% colnames(parsed) && is.list(parsed$key_markers)) {
    parsed$key_markers <- vapply(parsed$key_markers, function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      paste(x, collapse = ", ")
    }, character(1))
  }

  parsed$cluster <- as.character(parsed$cluster)
  if ("confidence" %in% colnames(parsed))
    parsed$confidence <- as.numeric(parsed$confidence)
  parsed
}

# =============================================================================
# AI Configuration UI/Server (Global API key + model)
# =============================================================================

aiConfigUI <- function(id) {
  ns <- NS(id)
  div(
    class = "ai-config-panel",
    style = "padding: 12px; margin: 10px 0; border: 1px solid #d9e2ef; border-radius: 8px; background: #f9fbff;",
    h5(style = "margin-top: 0;", "AI Assistant Configuration"),
    fluidRow(
      column(4,
        passwordInput(ns("api_key"), "Claude API Key",
                      placeholder = "sk-ant-api03-...", width = "100%"),
        helpText("Saved to this project only. Leave blank to keep the current key.",
                 style = "font-size: 11px; color: #888;")
      ),
      column(4,
        selectInput(ns("model"), "Default Model", width = "100%",
          choices = ai_model_choices(), selected = get_default_model()),
        conditionalPanel(
          condition = paste0("input['", ns("model"), "'].indexOf('opus') > -1"),
          helpText(
            icon("exclamation-triangle", style = "color:#e65100;"),
            "Opus models are ~5x more expensive than Sonnet. ",
            "Reserve for final reports and publication-ready text.",
            style = "font-size:11px; color:#e65100; margin-top:4px;"
          )
        ),
        conditionalPanel(
          condition = paste0("input['", ns("model"), "'].indexOf('thinking') > -1"),
          helpText(
            icon("lightbulb", style = "color:#1565c0;"),
            "Extended Thinking enables multi-step reasoning. ",
            "Best for complex biological interpretation.",
            style = "font-size:11px; color:#1565c0; margin-top:4px;"
          )
        )
      ),
      column(4,
        div(style = "margin-top: 25px;",
          uiOutput(ns("key_status"))
        ),
        br(),
        actionButton(ns("ai_test_connection"), "Test API Connection",
          icon = icon("plug"), class = "ai-btn btn-sm", width = "100%"),
        uiOutput(ns("connection_test_result")),
        br(),
        actionButton(ns("ai_suggest_actions"), "Suggest Next Actions",
          icon = icon("robot"), class = "btn-default btn-sm", width = "100%")
      )
    ),
    aiActionQueueUI(ns("actionQueue"))
  )
}

aiConfigServer <- function(id, values, auth_state = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Hydrate controls from the current project state when the panel appears
    observe({
      updateTextInput(session, "api_key", value = values$ai_api_key %||% "")
    })
    observeEvent(TRUE, {
      ch  <- tryCatch(ai_model_choices(api_key = values$ai_api_key),
                      error = function(e) ai_model_choices())
      sel <- values$ai_model %||% get_default_model(values$ai_api_key)
      updateSelectInput(session, "model", choices = ch,
                        selected = if (sel %in% ch) sel else get_default_model())
    }, once = TRUE)

    # Sync API key to global values (+ persist to the open project's DB row)
    observeEvent(input$api_key, {
      values$ai_api_key <- input$api_key
      if (!is.null(auth_state) && !is.null(auth_state$current_project_id) &&
          nzchar(trimws(input$api_key %||% ""))) {
        tryCatch(update_project_metadata(auth_state$current_project_id, auth_state$user_id,
                                         ai_api_key = trimws(input$api_key)),
                 error = function(e) NULL)
        # refresh the live model list against the new key
        ch <- tryCatch(ai_model_choices(api_key = input$api_key, force = TRUE),
                       error = function(e) NULL)
        if (!is.null(ch)) updateSelectInput(session, "model", choices = ch,
                                            selected = isolate(input$model))
      }
    }, ignoreNULL = FALSE)

    # Sync model to global values (+ persist to the open project's DB row)
    observeEvent(input$model, {
      values$ai_model <- input$model
      if (!is.null(auth_state) && !is.null(auth_state$current_project_id)) {
        tryCatch(update_project_metadata(auth_state$current_project_id, auth_state$user_id,
                                         ai_model = input$model),
                 error = function(e) NULL)
      }
    })

    # Status indicator with model badge
    output$key_status <- renderUI({
      key <- input$api_key
      model_id <- input$model
      model_badge <- if (grepl("opus", model_id)) {
        tags$span(class = "badge", style = "background:#7b1fa2; color:#fff; margin-left:6px;",
          if (grepl("thinking", model_id)) "Opus 4.6 Thinking" else "Opus 4.6")
      } else if (grepl("sonnet-4-6", model_id)) {
        tags$span(class = "badge", style = "background:#0d47a1; color:#fff; margin-left:6px;",
          if (grepl("thinking", model_id)) "Sonnet 4.6 Thinking" else "Sonnet 4.6")
      } else if (grepl("sonnet-4-5", model_id)) {
        tags$span(class = "badge", style = "background:#1565c0; color:#fff; margin-left:6px;", "Sonnet 4.5")
      } else if (grepl("haiku", model_id)) {
        tags$span(class = "badge", style = "background:#2e7d32; color:#fff; margin-left:6px;", "Haiku 4.5")
      } else {
        tags$span(class = "badge", style = "background:#555; color:#fff; margin-left:6px;",
          substr(model_id, 1, 20))
      }

      if (is.null(key) || !nzchar(trimws(key))) {
        tagList(
          tags$span(style = "color: #e65100; font-size: 12px;",
                    icon("exclamation-triangle"), " No API key set"),
          model_badge
        )
      } else if (grepl("^sk-ant-", key)) {
        tagList(
          tags$span(style = "color: #2e7d32; font-size: 12px;",
                    icon("check-circle"), " API key configured"),
          model_badge
        )
      } else {
        tagList(
          tags$span(style = "color: #ff9800; font-size: 12px;",
                    icon("question-circle"), " Key format unusual (expected sk-ant-...)"),
          model_badge
        )
      }
    })

    # Test API Connection button
    observeEvent(input$ai_test_connection, {
      api_key <- values$ai_api_key
      if (is.null(api_key) || !nzchar(trimws(api_key))) {
        output$connection_test_result <- renderUI({
          div(style = "margin-top: 8px; padding: 8px; background: #fff3e0; border-radius: 6px; font-size: 12px;",
            icon("exclamation-triangle", style = "color: #e65100;"),
            " Please enter your API key first."
          )
        })
        return()
      }
      output$connection_test_result <- renderUI({
        div(style = "margin-top: 8px; padding: 8px; background: #f3f4f6; border-radius: 6px; font-size: 12px;",
          icon("spinner", class = "fa-spin"), " Testing connection..."
        )
      })
      model <- values$ai_model %||% get_default_model(values$ai_api_key)
      start_time <- Sys.time()
      result <- tryCatch(
        call_claude_api(trimws(api_key), model,
          list(list(role = "user", content = "Reply with exactly: OK")),
          max_tokens = 16L, timeout = 30),
        error = function(e) list(ok = FALSE, error = e$message)
      )
      elapsed <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
      if (isTRUE(result$ok)) {
        output$connection_test_result <- renderUI({
          div(style = "margin-top: 8px; padding: 8px; background: #e8f5e9; border-radius: 6px; font-size: 12px;",
            icon("check-circle", style = "color: #2e7d32;"),
            tags$strong(" Connected "),
            tags$span(style = "color: #555;",
              paste0("Model: ", model, " | Response: ", elapsed, "s"))
          )
        })
      } else {
        output$connection_test_result <- renderUI({
          div(style = "margin-top: 8px; padding: 8px; background: #ffebee; border-radius: 6px; font-size: 12px;",
            icon("times-circle", style = "color: #c62828;"),
            tags$strong(" Connection failed: "),
            tags$span(style = "color: #c62828;", result$error)
          )
        })
      }
    })

    # Suggest Next Actions button
    observeEvent(input$ai_suggest_actions, {
      req(values$ai_api_key)
      task <- paste0(
        "TASK: Next Analysis Step Recommendations\n\n",
        "Current state:\n",
        "- Data loaded: ",     isTRUE(values$data_loaded), "\n",
        "- QC done: ",         isTRUE(values$qc_done), "\n",
        "- Normalised: ",      isTRUE(values$normalization_done), "\n",
        "- Clustered: ",       isTRUE(values$clustering_done), "\n",
        "- Annotated: ",       isTRUE(values$annotation_done), "\n",
        "- DE done: ",         !is.null(values$last_de_results), "\n",
        "- Enrichment done: ", !is.null(values$last_enrich_results), "\n\n",
        "Suggest next 3 analysis steps. Use tools if applicable:\n",
        "set_clustering_resolution, run_de_analysis, run_enrichment, annotate_clusters\n",
        "Return JSON array:\n",
        '[{"name":"...","description":"...","module":"...",',
        '"parameters":{},"execute_fn_name":"<tool or null>"}]'
      )
      result <- call_ai("resolution_advisor",
        build_expert_prompt(task, values, include_dataset = FALSE), values)

      if (result$ok) {
        parsed <- tryCatch(jsonlite::fromJSON(result$text, simplifyDataFrame = FALSE),
                           error = function(e) NULL)
        if (!is.null(parsed) && is.list(parsed)) {
          new_q <- lapply(parsed, function(a) {
            a$status <- "pending"; a$timestamp <- Sys.time(); a
          })
          values$ai_action_queue <- c(values$ai_action_queue %||% list(), new_q)
          showNotification(paste(length(new_q), "actions queued for approval"),
            type = "message", duration = 5)
        }
      }
    })

    # Action queue sub-module
    aiActionQueueServer("actionQueue", values)
  })
}

# =============================================================================
# Tool Registry (structured AI action definitions)
# =============================================================================

get_ai_tool_registry <- function() {
  list(
    list(
      name = "set_clustering_resolution",
      description = "Re-run clustering with a specific resolution parameter",
      module = "clustering",
      parameters = list(
        resolution = list(type = "numeric", min = 0.1, max = 3.0, default = 0.5)
      ),
      execute_fn_name = "ai_tool_set_resolution"
    ),
    list(
      name = "run_de_analysis",
      description = "Run differential expression between two groups",
      module = "downstream",
      parameters = list(
        group_col = list(type = "character", description = "Metadata column for grouping"),
        group1 = list(type = "character", description = "First group identity"),
        group2 = list(type = "character", description = "Second group identity"),
        test = list(type = "character", default = "wilcox",
                    choices = c("wilcox", "MAST", "t", "bimod"))
      ),
      execute_fn_name = "ai_tool_run_de"
    ),
    list(
      name = "run_enrichment",
      description = "Run pathway enrichment on DE results",
      module = "downstream",
      parameters = list(
        database = list(type = "character", default = "H",
                        choices = c("H", "KEGG", "REACTOME", "C5")),
        species = list(type = "character", default = "Homo sapiens")
      ),
      execute_fn_name = "ai_tool_run_enrichment"
    ),
    list(
      name = "annotate_clusters",
      description = "Run automatic cell type annotation",
      module = "annotation",
      parameters = list(
        method = list(type = "character", default = "enhanced",
                      choices = c("enhanced", "singler", "combined")),
        level = list(type = "character", default = "clusters",
                     choices = c("clusters", "cells"))
      ),
      execute_fn_name = "ai_tool_annotate"
    ),
    list(
      name = "generate_plot",
      description = "Generate a specific visualization",
      module = "visualization",
      parameters = list(
        plot_type = list(type = "character",
                         choices = c("dimplot", "featureplot", "violin", "dotplot", "heatmap")),
        genes = list(type = "character", description = "Comma-separated gene names"),
        group_by = list(type = "character", description = "Metadata column for grouping")
      ),
      execute_fn_name = "ai_tool_generate_plot"
    )
  )
}

# =============================================================================
# Action Queue UI/Server
# =============================================================================

aiActionQueueUI <- function(id) {
  ns <- NS(id)
  div(id = ns("action_queue_container"), style = "margin: 10px 0;",
    uiOutput(ns("action_queue_display"))
  )
}

aiActionQueueServer <- function(id, values) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      if (is.null(values$ai_action_queue)) values$ai_action_queue <- list()
    })

    # Render action queue
    output$action_queue_display <- renderUI({
      queue <- values$ai_action_queue
      if (is.null(queue) || length(queue) == 0) return(NULL)

      action_uis <- lapply(seq_along(queue), function(i) {
        action <- queue[[i]]
        status <- action$status %||% "pending"
        status_badge <- switch(status,
          "pending"  = tags$span(class = "badge", style = "background: #ff9800;", "Pending"),
          "approved" = tags$span(class = "badge", style = "background: #4caf50;", "Approved"),
          "rejected" = tags$span(class = "badge", style = "background: #f44336;", "Rejected"),
          "running"  = tags$span(class = "badge", style = "background: #2196f3;", "Running"),
          "done"     = tags$span(class = "badge", style = "background: #9e9e9e;", "Done"),
          "failed"   = tags$span(class = "badge", style = "background: #b71c1c;", "Failed"),
          tags$span(class = "badge", status)
        )
        params_str <- if (length(action$parameters) > 0) {
          paste(vapply(names(action$parameters), function(k)
            paste0(k, "=", action$parameters[[k]]), character(1)), collapse = ", ")
        } else ""

        div(style = "padding: 8px 12px; margin: 4px 0; border: 1px solid #e0e0e0; border-radius: 6px; background: white;",
          fluidRow(
            column(1, status_badge),
            column(3, tags$strong(action$name %||% "Unknown")),
            column(4, tags$span(style = "color: #666; font-size: 12px;", params_str)),
            column(2, tags$span(style = "color: #888; font-size: 11px;",
                                paste("Module:", action$module %||% "?"))),
            column(2,
              if (status == "pending") {
                div(
                  actionButton(ns(paste0("approve_", i)), "Approve",
                               class = "btn-success btn-xs", style = "margin-right: 4px;"),
                  actionButton(ns(paste0("reject_", i)), "Reject",
                               class = "btn-danger btn-xs")
                )
              }
            )
          )
        )
      })
      div(h5(style = "margin-top: 0;", "AI Action Queue"), do.call(tagList, action_uis))
    })

    # Observe approve/reject buttons
    observe({
      queue <- values$ai_action_queue
      if (is.null(queue) || length(queue) == 0) return()
      lapply(seq_along(queue), function(i) {
        observeEvent(input[[paste0("approve_", i)]], {
          q <- values$ai_action_queue
          if (i <= length(q)) { q[[i]]$status <- "approved"; values$ai_action_queue <- q }
        }, ignoreInit = TRUE, once = TRUE)
        observeEvent(input[[paste0("reject_", i)]], {
          q <- values$ai_action_queue
          if (i <= length(q)) { q[[i]]$status <- "rejected"; values$ai_action_queue <- q }
        }, ignoreInit = TRUE, once = TRUE)
      })
    })

    # Execute approved actions
    observe({
      queue <- values$ai_action_queue
      if (is.null(queue) || length(queue) == 0) return()
      lapply(seq_along(queue), function(i) {
        if (identical(queue[[i]]$status, "approved")) {
          q <- values$ai_action_queue
          q[[i]]$status <- "running"
          values$ai_action_queue <- q
          execute_tool(queue[[i]], values, session)
        }
      })
    })
  })
}

# =============================================================================
# Tool Executor
# =============================================================================

execute_tool <- function(action, values, session) {
  params <- action$parameters
  fn     <- action$execute_fn_name %||% ""

  result <- tryCatch(switch(fn,
    ai_tool_set_resolution = {
      req(!is.null(values$active_dataset))
      obj <- values$active_dataset
      resolution_use <- as.numeric(params$resolution)
      if (!is.finite(resolution_use) || resolution_use <= 0) {
        stop("Resolution must be a positive number.")
      }

      if (exists("run_standard_graph_clustering", mode = "function")) {
        last_run <- obj@misc$clustering_last_run %||% list()
        reduction_use <- last_run$reduction %||%
          if (exists("select_clustering_reduction", mode = "function")) {
            select_clustering_reduction(obj, requested = "pca", allow_2d = FALSE)
          } else if ("pca" %in% names(obj@reductions)) {
            "pca"
          } else {
            names(obj@reductions)[1]
          }
        if (is.null(reduction_use) || !reduction_use %in% names(obj@reductions)) {
          stop("No suitable reduction available for re-clustering.")
        }

        dims_use <- last_run$dims %||% min(30L, ncol(Seurat::Embeddings(obj, reduction_use)))
        k_use <- last_run$k_param %||% min(max(20L, floor(sqrt(ncol(obj)) / 2)), ncol(obj) - 1L)
        algo_use <- last_run$algorithm %||% "leiden"
        cluster_name_use <- last_run$cluster_name %||% "seurat_clusters"
        prune_use <- last_run$prune_snn %||% (1 / 15)
        graph_type_use <- last_run$graph_type %||% "snn"
        seed_use <- last_run$random_seed %||% 42L

        rerun <- run_standard_graph_clustering(
          obj,
          reduction = reduction_use,
          dims = dims_use,
          k_param = k_use,
          resolution = resolution_use,
          algorithm = algo_use,
          cluster_name = cluster_name_use,
          prune_snn = prune_use,
          graph_type = graph_type_use,
          graph_prefix = "ai_resolution",
          compute_quality = FALSE,
          random_seed = seed_use
        )
        obj <- rerun$object
        n_clusters <- length(unique(rerun$cluster_labels))
      } else {
        obj <- Seurat::FindClusters(
          obj,
          resolution = resolution_use,
          random.seed = 42L,
          verbose = FALSE
        )
        n_clusters <- length(unique(obj$seurat_clusters))
      }

      values$active_dataset  <- obj
      values$clustering_done <- TRUE
      list(ok = TRUE, msg = paste0("Re-clustered at resolution ",
        resolution_use, " -> ", n_clusters, " clusters"))
    },

    ai_tool_run_de = {
      req(!is.null(values$active_dataset))
      obj <- values$active_dataset
      Idents(obj) <- params$group_col
      markers <- Seurat::FindMarkers(obj, ident.1 = params$group1,
        ident.2 = if (params$group2 == "rest") NULL else params$group2,
        test.use = params$test %||% "wilcox", verbose = FALSE)
      markers$gene <- rownames(markers)
      values$last_de_results <- markers
      list(ok = TRUE, msg = paste0("DE: ", sum(markers$p_val_adj < 0.05, na.rm = TRUE), " sig genes"))
    },

    ai_tool_run_enrichment = {
      req(!is.null(values$last_de_results))
      de <- values$last_de_results
      sig <- de$gene[de$p_val_adj < 0.05 & !is.na(de$p_val_adj)]
      if (requireNamespace("clusterProfiler", quietly = TRUE)) {
        ora <- tryCatch(clusterProfiler::enrichGO(sig,
          OrgDb = if (params$species == "Homo sapiens") "org.Hs.eg.db" else "org.Mm.eg.db",
          keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH"), error = function(e) NULL)
        values$last_enrich_results <- list(ora = list(up = as.data.frame(ora)))
      }
      list(ok = TRUE, msg = "Enrichment complete")
    },

    ai_tool_annotate = {
      req(!is.null(values$active_dataset))
      if (exists("EnhancedCellTypeAnnotator") && exists("CELL_TYPE_MARKERS")) {
        obj <- values$active_dataset
        annotator <- EnhancedCellTypeAnnotator$new(CELL_TYPE_MARKERS)
        scores <- tryCatch(annotator$score_cell_types_enhanced(obj,
          assay = DefaultAssay(obj)), error = function(e) NULL)
        if (!is.null(scores)) {
          best <- colnames(scores$scores)[max.col(scores$scores, ties.method = "first")]
          obj@meta.data$cell_type_enhanced <- best
          values$active_dataset  <- obj
          values$annotation_done <- TRUE
        }
      }
      list(ok = TRUE, msg = "Annotation complete")
    },

    ai_tool_generate_plot = {
      values$nlplot_request <- params$plot_type
      list(ok = TRUE, msg = paste("Plot request sent:", params$plot_type))
    },

    list(ok = FALSE, msg = paste("Unknown tool:", fn))
  ), error = function(e) list(ok = FALSE, msg = e$message))

  # Update queue status
  q <- values$ai_action_queue
  idx <- which(vapply(q, function(x) identical(x$name, action$name), logical(1)))[1]
  if (!is.na(idx)) {
    q[[idx]]$status <- if (result$ok) "done" else "failed"
    q[[idx]]$result <- result$msg
    values$ai_action_queue <- q
  }
  if (exists("append_audit_log", mode = "function"))
    append_audit_log(values, module = "AI Queue", action = paste("executed:", fn),
      level = if (result$ok) "info" else "error",
      details = list(params = params, result = result$msg))
  showNotification(result$msg, type = if (result$ok) "message" else "error", duration = 6)
}

# =============================================================================
# AI Thread UI Helpers
# =============================================================================

#' Create a placeholder for an AI thread
ai_thread_ui <- function(ns, id) {
  uiOutput(ns(paste0("ai_thread_", id)))
}

#' Build a full thread message div
ai_thread_message <- function(content_html, title = "AI Assistant",
                               model = NULL, tokens_in = NA, tokens_out = NA,
                               apply_btn = NULL) {
  model_badge <- if (!is.null(model) && nzchar(model)) {
    badge_label <- if (grepl("opus", model, ignore.case = TRUE)) "OPUS"
                   else if (grepl("haiku", model, ignore.case = TRUE)) "HAIKU"
                   else "SONNET"
    tags$span(class = "ai-thread-model", badge_label)
  } else NULL

  token_content <- if (!is.na(tokens_in) && !is.na(tokens_out)) {
    tagList(icon("microchip"), paste0(tokens_in, " in / ", tokens_out, " out"))
  } else ""

  footer <- div(class = "ai-thread-footer",
    tags$span(class = "ai-thread-tokens", token_content),
    if (!is.null(apply_btn)) apply_btn else NULL
  )

  div(class = "ai-thread-container",
    # Gradient header bar
    div(class = "ai-thread-header-bar",
      div(class = "ai-thread-avatar", icon("robot")),
      tags$span(class = "ai-thread-sender", title),
      model_badge,
      tags$span(class = "ai-thread-time", format(Sys.time(), "%H:%M"))
    ),
    # Message body
    div(class = "ai-thread-message",
      div(class = "ai-thread-body",
        div(class = "ai-thread-content", HTML(content_html))
      ),
      footer
    )
  )
}

#' Build a loading spinner thread message
ai_thread_loading <- function(title = "AI is thinking...") {
  div(class = "ai-thread-loading",
    icon("spinner", class = "fa-spin"),
    tags$span(title)
  )
}

#' Build an error thread message
ai_thread_error <- function(error_msg) {
  div(class = "ai-thread-error",
    icon("exclamation-triangle"),
    div(style = "flex: 1;",
      tags$strong("Something went wrong"),
      tags$p(style = "margin: 4px 0 0 0; font-size: 13px; opacity: 0.85;", error_msg)
    )
  )
}

# =============================================================================
# Cross-Module AI Helper Button  — REMOVED
# =============================================================================
# The generic per-module "Ask AI" button was redundant with the main AI
# assistant and the purpose-built AI features (annotation, QC/DE advisors,
# downstream interpretation). These are kept as no-op stubs so any remaining
# call site is harmless: the button renders nothing and the server does nothing.

ai_helper_button_ui <- function(ns, module_name) NULL

ai_helper_server <- function(id, values, module_name, context_fn = NULL) invisible(NULL)
