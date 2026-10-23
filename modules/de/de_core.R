# modules/de/de_core.R
# Pure helpers shared by de_methods.R, de_server.R. No Shiny imports.

#' Map a species name (Latin or lowercase) to its OrgDb package name.
#' @param species character(1). Examples: "Homo sapiens", "human", "MOUSE".
#' @return character(1). Package name like "org.Hs.eg.db".
orgdb_for_species <- function(species) {
  key <- .normalize_species_key(species)
  map <- c(
    human     = "org.Hs.eg.db",
    mouse     = "org.Mm.eg.db",
    rat       = "org.Rn.eg.db",
    zebrafish = "org.Dr.eg.db"
  )
  pkg <- map[key]
  if (is.na(pkg)) stop("Unsupported species: ", species,
                       ". Supported: human, mouse, rat, zebrafish.")
  unname(pkg)
}

#' Map a species name to its KEGG 3-letter organism code.
kegg_code_for_species <- function(species) {
  key <- .normalize_species_key(species)
  map <- c(human = "hsa", mouse = "mmu", rat = "rno", zebrafish = "dre")
  code <- map[key]
  if (is.na(code)) stop("No KEGG code for species: ", species)
  unname(code)
}

#' Internal: collapse Latin or English species name to a canonical lowercase key.
.normalize_species_key <- function(species) {
  if (is.null(species) || !nzchar(species)) return(NA_character_)
  s <- tolower(trimws(species))
  if (s %in% c("human", "homo sapiens", "h. sapiens", "hs"))      return("human")
  if (s %in% c("mouse", "mus musculus", "m. musculus", "mm"))     return("mouse")
  if (s %in% c("rat", "rattus norvegicus", "r. norvegicus", "rn")) return("rat")
  if (s %in% c("zebrafish", "danio rerio", "d. rerio", "dr"))     return("zebrafish")
  NA_character_
}

#' Pick the most plausible cell-type / cluster column from a Seurat object.
#' @return character(1) or NULL.
detect_cell_type_col <- function(obj) {
  candidates <- c("cell_type", "celltype", "cell.type", "annotation",
                  "seurat_clusters", "cluster")
  cols <- colnames(obj@meta.data)
  found <- intersect(candidates, cols)
  if (length(found) > 0) return(found[1])
  cat_cols <- cols[vapply(obj@meta.data,
                          function(x) is.factor(x) || is.character(x),
                          logical(1))]
  if (length(cat_cols) > 0) cat_cols[1] else NULL
}

#' Pick the most plausible sample / replicate column. Returns NULL if none found.
detect_sample_col <- function(obj) {
  candidates <- c("orig.ident", "sample", "sample_id", "donor",
                  "patient", "replicate", "batch")
  found <- intersect(candidates, colnames(obj@meta.data))
  if (length(found) > 0) found[1] else NULL
}
