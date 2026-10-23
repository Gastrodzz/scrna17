# modules/species_utils.R
# Shared species detection + gene-pattern resolution. One source of truth for
# mitochondrial/ribosomal regexes and species labels, used by QC, preprocessing,
# and annotation so patterns never drift between modules.

if (!exists("%||%")) `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

.SPECIES_MARKERS <- list(
  human = c("PTPRC","EPCAM","PECAM1","CD19","CD3D","CD4","CD8A","FOXP3",
            "COL1A1","KRT19","MS4A1","GATA3","NKG7","S100A4"),
  mouse = c("Ptprc","Epcam","Pecam1","Cd19","Cd3d","Cd4","Cd8a","Foxp3",
            "Col1a1","Krt19","Ms4a1","Gata3","Nkg7","S100a4")
)

.SPECIES_PATTERNS <- list(
  human     = list(mt = "^MT-",  ribo = "^RP[SL]\\d"),
  mouse     = list(mt = "^mt-",  ribo = "^Rp[sl]\\d"),
  rat       = list(mt = "^Mt-",  ribo = "^Rp[sl]\\d"),
  zebrafish = list(mt = "^mt-",  ribo = "^rp[sl]\\d")
)

#' Map a free-text / latin species label to the app's canonical key.
normalize_species <- function(x) {
  if (is.null(x) || !nzchar(trimws(x))) return("")
  s <- tolower(trimws(x))
  if (grepl("sapiens|human|homo", s))          return("human")
  if (grepl("musculus|mouse|\\bmus\\b", s))    return("mouse")
  if (grepl("norvegicus|\\brat\\b|rattus", s)) return("rat")
  if (grepl("rerio|zebrafish|danio", s))       return("zebrafish")
  if (s %in% c("human","mouse","rat","zebrafish")) return(s)
  ""
}

#' Detect species from gene identifiers.
#' @return list(species, mt_pattern, ribo_pattern, confidence, note)
detect_species <- function(genes) {
  genes <- as.character(genes)
  if (length(genes) == 0) {
    return(list(species = "unknown", mt_pattern = "^MT-", ribo_pattern = "^RP[SL]\\d",
                confidence = "none", note = "No genes supplied."))
  }
  if (mean(grepl("^ENS[A-Z]*G\\d", genes)) > 0.5) {
    return(list(species = "unknown", mt_pattern = "^MT-", ribo_pattern = "^RP[SL]\\d",
                confidence = "low",
                note = "Matrix uses Ensembl gene IDs; convert to symbols for reliable MT/ribo detection."))
  }
  n_mt_h <- sum(grepl("^MT-", genes)); n_mt_m <- sum(grepl("^mt-", genes))
  n_h <- sum(genes %in% .SPECIES_MARKERS$human)
  n_m <- sum(genes %in% .SPECIES_MARKERS$mouse)
  h_score <- n_mt_h * 2 + n_h
  m_score <- n_mt_m * 2 + n_m
  sp <- if (h_score == 0 && m_score == 0) "unknown"
        else if (h_score >= m_score) "human" else "mouse"
  conf <- if (sp == "unknown") "none"
          else if (abs(h_score - m_score) >= 6) "high"
          else if (abs(h_score - m_score) >= 2) "medium" else "low"
  pat <- .SPECIES_PATTERNS[[sp]] %||% .SPECIES_PATTERNS$human
  list(species = sp, mt_pattern = pat$mt, ribo_pattern = pat$ribo,
       confidence = conf, note = sprintf("human_score=%d, mouse_score=%d", h_score, m_score))
}

#' MT/ribo patterns for a known species key (free-text accepted).
species_patterns <- function(species_key) {
  .SPECIES_PATTERNS[[normalize_species(species_key)]] %||% .SPECIES_PATTERNS$human
}
