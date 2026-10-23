# Extracted from test-gallery-integration.R:141

# prequel ----------------------------------------------------------------------
suppressPackageStartupMessages({
  library(shiny)
  source("../../modules/umap_gallery.R")
  source("../../modules/core_robustness.R")
  source("../../modules/de/de_core.R")
  source("../../modules/de/de_methods.R")
  source("../../modules/de/de_server.R")
  source("../../modules/pathway_network/pn_core.R")
  source("../../modules/pathway_network/pn_enrichment.R")
  source("../../modules/pathway_network/pn_network.R")
  source("../../modules/pathway_network/pn_server.R")
})
.fixture <- function() {
  set.seed(11)
  ifn <- c("STAT1","IRF1","JAK1","ISG15","MX1","OAS1","IFIT1","CXCL10","HLA-A",
           "B2M","TAP1","PSMB8","IRF7","IRF9","RSAD2","BST2","CCL5","CXCL11",
           "TAP2","PSMB9")
  emt <- c("VIM","SNAI1","ZEB1","TWIST1","FN1","CDH2","MMP2","TGFB1","COL1A1",
           "ACTA2","TAGLN","SPARC","CTGF","BMP2","SMAD2","SMAD3","LOX",
           "ITGB1","LOXL2","POSTN")
  mk <- function(g) {
    df <- data.frame(
      gene = g,
      avg_log2FC = runif(length(g), 1, 4),
      p_val      = 10^(-runif(length(g), 5, 30)),
      p_val_adj  = 10^(-runif(length(g), 4, 25)),
      pct.1      = runif(length(g), 0.2, 0.9),
      pct.2      = runif(length(g), 0.0, 0.3))
    attr(df, "de_audit") <- list(
      test = "wilcox", group_col = "genotype",
      ident1 = "ko", ident2 = "wt",
      n_cells_total = 200L, n_group1 = 100L, n_group2 = 100L,
      latent_vars = character(0),
      n_genes_tested = nrow(df),
      timestamp = Sys.time(),
      seurat_version = "5.4.0", r_version = R.version$version.string)
    df
  }
  de1 <- mk(ifn); de2 <- mk(emt)
  bg <- unique(c(ifn, emt, paste0("BG", 1:5000)))
  list(de1 = de1, de2 = de2, universe = bg,
       runs = list(
         "Epithelial: ko vs wt"  = list(data = de1, audit = attr(de1, "de_audit")),
         "Mesenchymal: ko vs wt" = list(data = de2, audit = attr(de2, "de_audit"))))
}

# test -------------------------------------------------------------------------
for (f in c("../../modules/de/de_server.R",
              "../../modules/pathway_network/pn_server.R")) {
    txt <- readLines(f, warn = FALSE)
    # Strip palette definition + comments + rgba tuples for transparent backgrounds
    bad <- grep('color = "#[0-9A-Fa-f]{6}"', txt, value = TRUE)
    # Should be zero — palette tokens are referenced as .DE_PAL$xxx / .PN_PAL$xxx
    expect_length(bad, 0)
  }
