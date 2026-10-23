# Cell Type Markers Database Module
# modules/cell_type_markers.R

# Comprehensive cell type marker database with scoring weights
# Weights range from 0.3 to 1.0, with highest weights for the most specific/canonical markers
# Negative weights indicate negative markers (markers that should NOT be expressed)
CELL_TYPE_MARKERS <- list(

  # ============================================================================
  # T-CELL LINEAGE
  # ============================================================================

  "T-cell (General)" = c(
    CD3D = 1.0, CD3E = 1.0, CD3G = 0.9, TRAC = 0.9, CD2 = 0.8,
    IL7R = 0.5, CD28 = 0.5, LCK = 0.4, ITK = 0.4, TRBC1 = 0.6
  ),

  "CD4+ Helper T-cell" = c(
    CD4 = 1.0, IL7R = 0.9, TCF7 = 0.8, LEF1 = 0.7, CCR7 = 0.7,
    CD3D = 0.6, CD3E = 0.6, SELL = 0.5, CD27 = 0.5, MAL = 0.4,
    CD40LG = 0.5, TRAC = 0.4
  ),

  "CD8+ Cytotoxic T-cell" = c(
    CD8A = 1.0, CD8B = 0.9, CD3D = 0.8, GZMA = 0.7, GZMB = 0.7,
    PRF1 = 0.6, NKG7 = 0.5, CST7 = 0.5, GZMK = 0.5, IFNG = 0.4,
    KLRG1 = 0.4, CD3E = 0.6
  ),

  "T-reg (Regulatory)" = c(
    FOXP3 = 1.0, IL2RA = 0.9, CTLA4 = 0.8, TIGIT = 0.7, IKZF2 = 0.7,
    CD4 = 0.5, TNFRSF18 = 0.6, TNFRSF4 = 0.5, BATF = 0.4, IL10 = 0.4,
    RTKN2 = 0.5, LRRC32 = 0.4
  ),

  "T-cell (Naive)" = c(
    CCR7 = 1.0, LEF1 = 0.9, SELL = 0.9, TCF7 = 0.8, CD27 = 0.7,
    CD3D = 0.6, CD3E = 0.6, IL7R = 0.5, NOSIP = 0.4, PIK3IP1 = 0.4,
    MAL = 0.5, LDHB = 0.3
  ),

  "T-cell (Memory)" = c(
    CD44 = 0.9, S100A4 = 0.8, ANXA1 = 0.7, IL7R = 0.7, GPR183 = 0.6,
    CD69 = 0.5, FOS = 0.4, JUN = 0.4, AQP3 = 0.5, KLRB1 = 0.4,
    CD3D = 0.5, TRAC = 0.4
  ),

  "T-cell (Exhausted)" = c(
    PDCD1 = 1.0, LAG3 = 0.9, HAVCR2 = 0.9, TIGIT = 0.8, CTLA4 = 0.7,
    TOX = 0.8, ENTPD1 = 0.7, CXCL13 = 0.6, LAYN = 0.5, CD3D = 0.4,
    ITGAE = 0.4, CD38 = 0.4
  ),

  "Th1 Cells" = c(
    TBX21 = 1.0, IFNG = 0.9, CXCR3 = 0.8, STAT4 = 0.7, TNF = 0.6,
    CCR5 = 0.6, IL12RB2 = 0.5, CD4 = 0.5, CD3D = 0.4, STAT1 = 0.4,
    IL2 = 0.4, EOMES = 0.3
  ),

  "Th2 Cells" = c(
    GATA3 = 1.0, IL4 = 0.9, IL13 = 0.8, IL5 = 0.8, CCR4 = 0.7,
    PTGDR2 = 0.7, STAT6 = 0.6, IL4R = 0.5, CD4 = 0.5, CD3D = 0.4,
    CRLF2 = 0.4, HPGDS = 0.5
  ),

  "Th17 Cells" = c(
    RORC = 1.0, IL17A = 0.9, IL17F = 0.8, CCR6 = 0.8, IL23R = 0.7,
    IL22 = 0.6, STAT3 = 0.5, CD4 = 0.5, CD3D = 0.4, AHR = 0.5,
    KLRB1 = 0.4, CCL20 = 0.4
  ),

  "Follicular Helper T Cells (Tfh)" = c(
    CXCR5 = 1.0, BCL6 = 0.9, PDCD1 = 0.8, ICOS = 0.8, CD200 = 0.7,
    SH2D1A = 0.7, MAF = 0.6, CXCL13 = 0.6, TOX2 = 0.5, CD4 = 0.5,
    CD3D = 0.4, IL21 = 0.5
  ),

  "Gamma-delta T-cell" = c(
    TRDC = 1.0, TRGC1 = 0.9, TRGC2 = 0.9, TRDV2 = 0.8, CD3D = 0.6,
    CD3E = 0.5, NKG7 = 0.5, KLRB1 = 0.4, GNLY = 0.4, GZMB = 0.4,
    IL17A = 0.3, ZBTB16 = 0.3
  ),

  "NKT Cells" = c(
    KLRB1 = 1.0, CD3D = 0.8, CD3E = 0.7, NCAM1 = 0.7, ZBTB16 = 0.8,
    NKG7 = 0.6, GNLY = 0.5, GZMB = 0.5, TRAV10 = 0.9, SLC7A5 = 0.4,
    IFNG = 0.4, CD1D = 0.3
  ),

  "MAIT cell" = c(
    SLC4A10 = 1.0, KLRB1 = 0.9, ZBTB16 = 0.7, RORC = 0.6, `TRAV1-2` = 0.9,
    CD3D = 0.5, CD3E = 0.5, NKG7 = 0.4, GNLY = 0.4, DPP4 = 0.5,
    IL18RAP = 0.5, NCR3 = 0.3
  ),

  # ============================================================================
  # B-CELL LINEAGE
  # ============================================================================

  "B-cell (General)" = c(
    MS4A1 = 1.0, CD19 = 1.0, CD79A = 0.9, CD79B = 0.8, PAX5 = 0.7,
    BANK1 = 0.6, CD22 = 0.6, BLK = 0.5, FCER2 = 0.5, TCL1A = 0.4,
    BLNK = 0.4, CD74 = 0.3
  ),

  "B-cell (Naive)" = c(
    MS4A1 = 1.0, IGHD = 1.0, IGHM = 0.9, FCER2 = 0.7, TCL1A = 0.7,
    CD19 = 0.6, CD79A = 0.6, IL4R = 0.5, PAX5 = 0.5, CD22 = 0.4,
    SELL = 0.4, BTK = 0.3
  ),

  "B-cell (Memory)" = c(
    MS4A1 = 1.0, CD27 = 0.9, IGHG1 = 0.8, IGHG2 = 0.7, CD80 = 0.6,
    AIM2 = 0.5, CD19 = 0.5, CD79A = 0.5, TNFRSF13B = 0.5, IGHG3 = 0.4,
    GPR183 = 0.4, IGHA1 = 0.4
  ),

  "Plasma cell" = c(
    SDC1 = 1.0, MZB1 = 0.9, XBP1 = 0.9, JCHAIN = 0.8, IGHG1 = 0.7,
    PRDM1 = 0.7, IRF4 = 0.6, TNFRSF17 = 0.6, CD38 = 0.6, DERL3 = 0.5,
    SSR4 = 0.4, FKBP11 = 0.4
  ),

  "Plasmablast" = c(
    CD38 = 1.0, SLAMF7 = 0.8, PRDM1 = 0.8, MZB1 = 0.7, XBP1 = 0.7,
    IRF4 = 0.7, TNFRSF17 = 0.6, JCHAIN = 0.6, MKI67 = 0.5, TOP2A = 0.4,
    SDC1 = 0.5, DERL3 = 0.4
  ),

  "Germinal Center B-cell" = c(
    BCL6 = 1.0, AICDA = 0.9, MME = 0.8, SUGCT = 0.6, RGS13 = 0.7,
    MARCKSL1 = 0.5, STMN1 = 0.5, MEF2B = 0.7, MKI67 = 0.4, ELL3 = 0.4,
    LRMP = 0.5, SERPINA9 = 0.6
  ),

  "Pre-B Cells" = c(
    VPREB1 = 1.0, IGLL1 = 0.9, CD79A = 0.8, CD19 = 0.7, RAG1 = 0.7,
    RAG2 = 0.7, DNTT = 0.6, CD34 = 0.5, IL7R = 0.5, PAX5 = 0.6,
    VPREB3 = 0.5, CD10 = 0.4
  ),

  "Pro-B Cells" = c(
    CD34 = 1.0, DNTT = 0.9, RAG1 = 0.8, RAG2 = 0.8, CD19 = 0.6,
    IL7R = 0.7, FLT3 = 0.6, PAX5 = 0.5, CD79A = 0.5, CD10 = 0.4,
    MME = 0.5, IGLL1 = 0.3
  ),

  "Transitional B Cells" = c(
    MS4A1 = 1.0, CD24 = 0.9, CD38 = 0.8, IGHM = 0.8, IGHD = 0.7,
    CD19 = 0.6, CD79A = 0.6, CD10 = 0.5, MME = 0.5, TCL1A = 0.6,
    CD9 = 0.4, FCER2 = 0.4
  ),

  # ============================================================================
  # NK CELLS AND ILCs
  # ============================================================================

  "NK cell" = c(
    NKG7 = 1.0, GNLY = 0.9, NCAM1 = 0.8, KLRF1 = 0.8, KLRD1 = 0.7,
    KLRK1 = 0.6, FCGR3A = 0.6, PRF1 = 0.5, GZMB = 0.5, GZMA = 0.4,
    NCR1 = 0.5, TYROBP = 0.4
  ),

  "Cytotoxic NK Cells" = c(
    FCGR3A = 1.0, FGFBP2 = 0.9, CX3CR1 = 0.8, GZMB = 0.8, PRF1 = 0.7,
    NKG7 = 0.7, GNLY = 0.6, KLRD1 = 0.6, GZMH = 0.5, SPON2 = 0.5,
    S1PR5 = 0.4, KLRF1 = 0.4
  ),

  "Regulatory NK Cells" = c(
    NCAM1 = 1.0, SELL = 0.9, GZMK = 0.8, XCL1 = 0.7, XCL2 = 0.7,
    GPR183 = 0.6, CD62L = 0.5, KLRC1 = 0.5, IL7R = 0.4, NKG7 = 0.4,
    CD27 = 0.5, KIT = 0.3
  ),

  "CD56-bright NK" = c(
    NCAM1 = 1.0, SELL = 0.9, GZMK = 0.8, XCL1 = 0.7, XCL2 = 0.7,
    IL7R = 0.5, CD27 = 0.5, KLRC1 = 0.5, GPR183 = 0.4, NKG7 = 0.5,
    GNLY = 0.4, CD62L = 0.4
  ),

  "CD56-dim NK" = c(
    FCGR3A = 1.0, FGFBP2 = 0.9, CX3CR1 = 0.8, GZMB = 0.8, PRF1 = 0.7,
    NKG7 = 0.7, GNLY = 0.6, KLRD1 = 0.6, S1PR5 = 0.5, SPON2 = 0.5,
    GZMH = 0.4, KLRF1 = 0.4
  ),

  "ILC1" = c(
    TBX21 = 1.0, IFNG = 0.9, IL12RB2 = 0.8, IL7R = 0.7, KLRB1 = 0.6,
    NCR1 = 0.5, EOMES = 0.4, CD127 = 0.5, CD161 = 0.5, TNF = 0.4,
    KIT = -0.6, RORC = -0.5
  ),

  "ILC2" = c(
    GATA3 = 1.0, IL1RL1 = 0.9, PTGDR2 = 0.9, HPGDS = 0.8, IL7R = 0.7,
    KLRB1 = 0.6, AREG = 0.5, IL13 = 0.5, IL4 = 0.4, IL5 = 0.4,
    KIT = 0.4, IL33 = 0.3
  ),

  "ILC3" = c(
    RORC = 1.0, IL23R = 0.9, IL22 = 0.8, KIT = 0.8, IL7R = 0.7,
    AHR = 0.7, NCR2 = 0.6, KLRB1 = 0.5, IL17A = 0.5, IL17F = 0.4,
    CCR6 = 0.5, LTA = 0.3
  ),

  # ============================================================================
  # MYELOID LINEAGE
  # ============================================================================

  "Myeloid (General)" = c(
    LYZ = 1.0, CD68 = 0.8, SPI1 = 0.7, AIF1 = 0.6, CSF1R = 0.5,
    CD14 = 0.5, TYROBP = 0.5, FCER1G = 0.4, FTL = 0.3, MNDA = 0.4
  ),

  "Monocyte" = c(
    CD14 = 1.0, VCAN = 0.9, FCN1 = 0.9, S100A9 = 0.8, S100A8 = 0.8,
    LYZ = 0.7, CD68 = 0.5, MNDA = 0.6, CTSS = 0.5, LST1 = 0.4,
    SERPINA1 = 0.4, MS4A6A = 0.3
  ),

  "Classical Monocytes (CD14++)" = c(
    CD14 = 1.0, LYZ = 0.9, S100A9 = 0.9, S100A8 = 0.8, VCAN = 0.8,
    FCN1 = 0.8, CD163 = 0.5, MNDA = 0.6, S100A12 = 0.7, CTSS = 0.5,
    IL1B = 0.4, CCR2 = 0.4
  ),

  "Non-classical Monocytes (CD16++)" = c(
    FCGR3A = 1.0, LST1 = 0.9, MS4A7 = 0.8, CDKN1C = 0.8, LILRB2 = 0.7,
    CSF1R = 0.6, CX3CR1 = 0.6, IFITM3 = 0.5, AIF1 = 0.5, HES4 = 0.5,
    TCF7L2 = 0.4, RHOC = 0.4
  ),

  "Intermediate Monocytes" = c(
    CD14 = 0.8, FCGR3A = 0.7, HLA_DPA1 = 0.8, HLA_DPB1 = 0.7, CCR5 = 0.7,
    CSF1R = 0.6, CD163 = 0.5, S100A8 = 0.5, S100A9 = 0.5, CLEC10A = 0.4,
    LILRB2 = 0.4, AIF1 = 0.3
  ),

  "Macrophage" = c(
    CD68 = 1.0, CD163 = 0.9, MRC1 = 0.8, APOE = 0.7, C1QA = 0.7,
    C1QB = 0.6, C1QC = 0.6, MSR1 = 0.5, MARCO = 0.5, SIGLEC1 = 0.4,
    MERTK = 0.5, AIF1 = 0.4
  ),

  "Macrophage M1" = c(
    CD80 = 0.9, CD86 = 0.8, TNF = 0.8, IL1B = 0.8, CCL3 = 0.7,
    CXCL9 = 0.7, CXCL10 = 0.7, CXCL11 = 0.6, NOS2 = 0.6, IDO1 = 0.5,
    STAT1 = 0.5, IRF5 = 0.4
  ),

  "Macrophage M2" = c(
    CD163 = 1.0, MRC1 = 0.9, ARG1 = 0.8, CCL18 = 0.7, CD206 = 0.7,
    TGFB1 = 0.5, MSR1 = 0.6, APOE = 0.6, C1QA = 0.5, MAF = 0.4,
    STAB1 = 0.4, FOLR2 = 0.5
  ),

  "Alveolar Macrophages" = c(
    MARCO = 1.0, FABP4 = 0.9, MCEMP1 = 0.8, RBP4 = 0.7, PPARG = 0.7,
    CD68 = 0.6, SIGLEC1 = 0.5, MRC1 = 0.5, INHBA = 0.5, SERPING1 = 0.4,
    AIF1 = 0.4, MSR1 = 0.4
  ),

  "Kupffer Cells" = c(
    CD68 = 1.0, CD163 = 0.9, MARCO = 0.8, TIMD4 = 0.8, VSIG4 = 0.7,
    C1QA = 0.7, C1QB = 0.6, CLEC4F = 0.9, AIF1 = 0.5, CD5L = 0.6,
    SLC40A1 = 0.4, FOLR2 = 0.5
  ),

  "Microglia" = c(
    AIF1 = 1.0, P2RY12 = 0.9, CX3CR1 = 0.8, CSF1R = 0.7, TMEM119 = 0.9,
    HEXB = 0.7, TREM2 = 0.6, GPR34 = 0.5, C1QA = 0.5, C1QB = 0.4,
    SALL1 = 0.6, SIGLECH = 0.5
  ),

  "Dendritic cell (DC)" = c(
    FCER1A = 0.9, CST3 = 0.8, HLA_DQA1 = 0.8, ITGAX = 0.7, HLA_DRA = 0.7,
    HLA_DRB1 = 0.6, CD74 = 0.5, HLA_DPA1 = 0.5, CPVL = 0.4, FLT3 = 0.5,
    ZBTB46 = 0.4, IRF4 = 0.3
  ),

  "Plasmacytoid DCs" = c(
    IL3RA = 1.0, TCF4 = 0.9, CLEC4C = 0.9, GZMB = 0.7, IRF7 = 0.8,
    LILRA4 = 0.8, JCHAIN = 0.6, PLD4 = 0.6, SERPINF1 = 0.5, ITM2C = 0.5,
    DERL3 = 0.4, SCT = 0.3
  ),

  "Conventional DC1" = c(
    CLEC9A = 1.0, XCR1 = 0.9, BATF3 = 0.8, IRF8 = 0.8, IDO1 = 0.6,
    CADM1 = 0.7, WDFY4 = 0.7, CLNK = 0.5, THBD = 0.6, HLA_DRA = 0.5,
    NAAA = 0.4, SNX22 = 0.4
  ),

  "Conventional DC2" = c(
    CD1C = 1.0, CLEC10A = 0.9, SIRPA = 0.8, FCER1A = 0.7, CD1E = 0.7,
    IRF4 = 0.6, ITGAX = 0.5, HLA_DQA1 = 0.5, PLD4 = 0.4, PKIB = 0.5,
    ADAM8 = 0.4, GSN = 0.3
  ),

  "Mast Cells" = c(
    KIT = 1.0, CPA3 = 0.9, TPSAB1 = 0.9, TPSB2 = 0.8, MS4A2 = 0.8,
    HPGDS = 0.7, HDC = 0.6, FCER1A = 0.6, GATA2 = 0.5, CLU = 0.4,
    CTSG = 0.4, RGS13 = 0.3
  ),

  "Basophils" = c(
    MS4A2 = 1.0, HDC = 0.9, GATA2 = 0.8, CLC = 0.8, CPA3 = 0.7,
    IL4 = 0.6, CD123 = 0.6, FCER1A = 0.7, CD200R3 = 0.5, ENPP3 = 0.5,
    IL1RL1 = 0.4, SLC45A3 = 0.3
  ),

  "Eosinophils" = c(
    RNASE2 = 1.0, RNASE3 = 0.9, CCR3 = 0.8, EPX = 0.8, CLC = 0.7,
    PRG2 = 0.7, SIGLEC8 = 0.8, IL5RA = 0.6, ALOX15 = 0.5, GATA1 = 0.5,
    OLIG2 = 0.4, CD69 = 0.3
  ),

  "Neutrophils" = c(
    S100A8 = 1.0, S100A9 = 1.0, FCGR3B = 0.9, CSF3R = 0.9, CXCR2 = 0.8,
    MMP9 = 0.7, ELANE = 0.6, CAMP = 0.6, LTF = 0.5, MMP8 = 0.5,
    CEACAM8 = 0.6, AZU1 = 0.4
  ),

  "Langerhans cell" = c(
    CD207 = 1.0, CD1A = 0.9, EPCAM = 0.7, FCER1A = 0.6, CD1C = 0.5,
    HLA_DRA = 0.5, ITGAX = 0.4, LANGERIN = 0.8, CD1E = 0.4, S100B = 0.3,
    AXL = 0.3, BIRC3 = 0.3
  ),

  # ============================================================================
  # HEMATOPOIETIC STEM AND PROGENITOR CELLS
  # ============================================================================

  "HSCs (Hematopoietic Stem Cells)" = c(
    CD34 = 1.0, KIT = 0.8, SPINK2 = 0.9, PROM1 = 0.7, THY1 = 0.7,
    CRHBP = 0.8, HLF = 0.7, AVP = 0.6, MLLT3 = 0.5, HOPX = 0.5,
    MECOM = 0.5, PROCR = 0.4
  ),

  "Common Myeloid Progenitors" = c(
    CD34 = 1.0, KIT = 0.9, CD38 = 0.7, IL3RA = 0.7, FLT3 = 0.6,
    GATA2 = 0.7, CEBPA = 0.5, MPO = 0.4, SPI1 = 0.4, GATA1 = 0.3,
    CSF2RA = 0.4, LMO2 = 0.3
  ),

  "Granulocyte-Monocyte Progenitors" = c(
    CD34 = 0.9, SPI1 = 0.8, CEBPA = 0.8, MPO = 0.7, CSF3R = 0.6,
    CSF2RA = 0.6, ELANE = 0.5, AZU1 = 0.5, CTSG = 0.4, PRTN3 = 0.4,
    LYZ = 0.4, KIT = 0.3
  ),

  "Megakaryocytes" = c(
    ITGA2B = 1.0, GP9 = 0.9, PF4 = 0.9, PPBP = 0.8, GP1BB = 0.8,
    TUBB1 = 0.7, TREML1 = 0.6, CMTM5 = 0.5, CLEC1B = 0.5, GATA1 = 0.5,
    FLI1 = 0.4, NFE2 = 0.4
  ),

  "Erythroid Progenitors" = c(
    GYPA = 1.0, KLF1 = 0.9, TFRC = 0.8, GATA1 = 0.8, HBB = 0.7,
    HBA1 = 0.7, EPOR = 0.6, CA1 = 0.5, AHSP = 0.6, ALAS2 = 0.6,
    TFR2 = 0.4, ANK1 = 0.4
  ),

  "Erythrocytes" = c(
    HBB = 1.0, HBA1 = 1.0, HBA2 = 1.0, ALAS2 = 0.8, SLC4A1 = 0.7,
    GYPA = 0.7, ANK1 = 0.6, CA1 = 0.6, SLC25A37 = 0.5, EPB42 = 0.5,
    AHSP = 0.5, SNCA = 0.3
  ),

  "Platelets" = c(
    PPBP = 1.0, PF4 = 1.0, GP9 = 0.9, ITGA2B = 0.8, GP1BB = 0.8,
    TUBB1 = 0.7, TREML1 = 0.6, CLEC1B = 0.5, CMTM5 = 0.5, SPARC = 0.4,
    GNG11 = 0.4, NRGN = 0.4
  ),

  # ============================================================================
  # EPITHELIAL LINEAGE
  # ============================================================================

  "Epithelial (General)" = c(
    EPCAM = 1.0, KRT8 = 0.9, KRT18 = 0.9, CDH1 = 0.8, KRT19 = 0.7,
    CLDN4 = 0.6, CLDN7 = 0.5, TJP1 = 0.5, OCLN = 0.4, MUC1 = 0.4,
    GRHL2 = 0.4, CRB3 = 0.3
  ),

  "AT1 Pneumocytes" = c(
    AGER = 1.0, PDPN = 0.9, AQP5 = 0.8, HOPX = 0.8, EMP2 = 0.7,
    CLIC5 = 0.7, CAV1 = 0.6, RTKN2 = 0.5, GPRC5A = 0.5, VEGFA = 0.4,
    ICAM1 = 0.3, SCNN1A = 0.3
  ),

  "AT2 Pneumocytes" = c(
    SFTPC = 1.0, SFTPA1 = 0.9, SFTPA2 = 0.9, SFTPB = 0.8, SFTPD = 0.7,
    ABCA3 = 0.8, SLC34A2 = 0.7, NAPSA = 0.6, LAMP3 = 0.5, PGC = 0.5,
    NKX2_1 = 0.5, ETV5 = 0.4
  ),

  "Club Cells" = c(
    SCGB1A1 = 1.0, SCGB3A2 = 0.9, CYP2F1 = 0.7, CYP2B6 = 0.6, BPIFB1 = 0.6,
    KRT15 = 0.5, MSLN = 0.4, MUC5B = 0.5, CP = 0.4, WFDC2 = 0.5,
    PIGR = 0.3, SAA1 = 0.3
  ),

  "Ciliated Cells" = c(
    FOXJ1 = 1.0, PIFO = 0.8, TMEM212 = 0.7, DNAH5 = 0.7, CFAP299 = 0.6,
    RSPH1 = 0.6, TPPP3 = 0.7, DNALI1 = 0.5, CAPS = 0.5, CCDC78 = 0.4,
    SNTN = 0.5, CCDC153 = 0.4
  ),

  "Goblet Cells" = c(
    MUC5AC = 1.0, MUC5B = 0.9, TFF3 = 0.8, SPDEF = 0.8, MUC2 = 0.7,
    AGR2 = 0.6, FCGBP = 0.5, CLCA1 = 0.5, ITLN1 = 0.5, ZG16 = 0.4,
    TFF1 = 0.4, ATOH1 = 0.3
  ),

  "Basal Epithelial Cells" = c(
    KRT5 = 1.0, KRT14 = 0.9, TP63 = 0.9, KRT17 = 0.7, KRT15 = 0.6,
    ITGA6 = 0.6, ITGB4 = 0.5, DST = 0.5, COL17A1 = 0.4, S100A2 = 0.5,
    LAMB3 = 0.4, LAMC2 = 0.3
  ),

  "Luminal Epithelial Cells" = c(
    KRT8 = 1.0, KRT18 = 0.9, KRT19 = 0.8, EPCAM = 0.8, CDH1 = 0.7,
    MUC1 = 0.6, CLDN3 = 0.5, CLDN4 = 0.5, ESR1 = 0.4, GATA3 = 0.4,
    FOXA1 = 0.4, AGR2 = 0.3
  ),

  "Myoepithelial Cells" = c(
    ACTA2 = 1.0, MYH11 = 0.8, KRT14 = 0.8, KRT5 = 0.7, TP63 = 0.7,
    MYLK = 0.6, CNN1 = 0.5, TAGLN = 0.5, OXTR = 0.4, CD10 = 0.4,
    ITGA6 = 0.4, SNAI2 = 0.3
  ),

  "Ductal Cells" = c(
    KRT19 = 1.0, KRT7 = 0.9, SOX9 = 0.8, CFTR = 0.7, SLC4A4 = 0.6,
    CA2 = 0.6, MUC1 = 0.5, SPP1 = 0.5, EPCAM = 0.4, HNF1B = 0.5,
    AQP1 = 0.4, TFF1 = 0.3
  ),

  "Acinar Cells" = c(
    PRSS1 = 1.0, CPA1 = 0.9, CTRB1 = 0.9, CELA3A = 0.8, AMY2A = 0.8,
    PNLIP = 0.7, CPA2 = 0.7, CELA2A = 0.6, PRSS2 = 0.6, REG1A = 0.5,
    CTRL = 0.4, CLPS = 0.4
  ),

  "Cholangiocytes" = c(
    KRT19 = 1.0, KRT7 = 0.9, SOX9 = 0.8, EPCAM = 0.7, HNF1B = 0.7,
    CFTR = 0.6, AQP1 = 0.6, TFF2 = 0.5, MUC6 = 0.5, SLC4A2 = 0.4,
    SPP1 = 0.4, CLDN1 = 0.3
  ),

  "Hepatocytes" = c(
    ALB = 1.0, TTR = 0.9, APOB = 0.8, CYP3A4 = 0.7, CYP2E1 = 0.7,
    HNF4A = 0.8, SERPINA1 = 0.6, TF = 0.6, ASGR1 = 0.6, TAT = 0.5,
    CPS1 = 0.5, PCK1 = 0.4, G6PC = 0.4
  ),

  "Urothelial Cells" = c(
    UPK1A = 1.0, UPK2 = 0.9, UPK3A = 0.9, KRT20 = 0.7, CLDN4 = 0.6,
    KRT18 = 0.6, KRT8 = 0.5, EPCAM = 0.5, PSCA = 0.5, KRT13 = 0.4,
    GATA3 = 0.4, PPARG = 0.3
  ),

  "Corneal Epithelial Cells" = c(
    KRT3 = 1.0, KRT12 = 0.9, KRT14 = 0.7, PAX6 = 0.8, TP63 = 0.6,
    KRT15 = 0.5, CLU = 0.5, ALDH3A1 = 0.6, DST = 0.4, CLDN1 = 0.4,
    KRT24 = 0.3, MUC16 = 0.3
  ),

  # ============================================================================
  # ENDOCRINE / ISLET CELLS
  # ============================================================================

  "Beta Cells (Pancreatic)" = c(
    INS = 1.0, IAPP = 0.9, MAFA = 0.8, PDX1 = 0.8, NKX6_1 = 0.7,
    NKX2_2 = 0.6, SLC30A8 = 0.7, HADH = 0.5, G6PC2 = 0.6, PCSK1 = 0.5,
    GCK = 0.4, ABCC8 = 0.4
  ),

  "Alpha Cells (Pancreatic)" = c(
    GCG = 1.0, ARX = 0.9, IRX2 = 0.8, TTR = 0.7, GC = 0.6,
    LOXL4 = 0.5, DPP4 = 0.5, FAP = 0.4, PCSK2 = 0.6, MAFB = 0.5,
    SLC7A2 = 0.4, CRYBA2 = 0.3
  ),

  "Delta Cells (Pancreatic)" = c(
    SST = 1.0, HHEX = 0.9, RBP4 = 0.7, LEPR = 0.6, PCSK1 = 0.5,
    GABRB3 = 0.6, PRG4 = 0.5, BCHE = 0.5, UCP2 = 0.4, SEC11C = 0.4,
    GHSR = 0.3, CD24 = 0.3
  ),

  "PP Cells (Pancreatic)" = c(
    PPY = 1.0, SERTM1 = 0.8, ETV1 = 0.7, ARX = 0.5, MEIS2 = 0.5,
    CARTPT = 0.4, SLITRK6 = 0.5, ENTPD3 = 0.4, ID4 = 0.4, FOLR1 = 0.3,
    ABCC9 = 0.3, SLC6A4 = 0.3
  ),

  # ============================================================================
  # NEURAL LINEAGE
  # ============================================================================

  "Neuron" = c(
    RBFOX3 = 1.0, SYN1 = 0.9, SNAP25 = 0.8, TUBB3 = 0.8, MAP2 = 0.7,
    NEFL = 0.6, NEFM = 0.6, SYP = 0.5, STMN2 = 0.6, ENO2 = 0.5,
    NRGN = 0.4, THY1 = 0.3
  ),

  "Excitatory Neurons" = c(
    SLC17A7 = 1.0, CAMK2A = 0.9, SATB2 = 0.8, NRGN = 0.7, SLC17A6 = 0.7,
    GRIA1 = 0.6, GRIA2 = 0.6, RBFOX3 = 0.5, SNAP25 = 0.5, SYN1 = 0.4,
    TBR1 = 0.5, NEUROD6 = 0.4
  ),

  "Inhibitory Neurons" = c(
    GAD1 = 1.0, GAD2 = 1.0, SLC32A1 = 0.9, DLX1 = 0.7, DLX2 = 0.7,
    DLX5 = 0.6, SST = 0.5, PVALB = 0.5, VIP = 0.5, ADARB2 = 0.4,
    LAMP5 = 0.4, RBFOX3 = 0.3
  ),

  "Dopaminergic Neurons" = c(
    TH = 1.0, SLC6A3 = 0.9, DDC = 0.8, NR4A2 = 0.8, PITX3 = 0.7,
    EN1 = 0.6, EN2 = 0.5, LMX1A = 0.6, FOXA2 = 0.5, SLC18A2 = 0.6,
    ALDH1A1 = 0.4, KCNJ6 = 0.4
  ),

  "GABAergic Neurons" = c(
    GAD1 = 1.0, GAD2 = 1.0, SLC32A1 = 0.9, ABAT = 0.7, SLC6A1 = 0.7,
    DLX1 = 0.6, DLX2 = 0.6, PVALB = 0.5, SST = 0.5, CALB1 = 0.4,
    CALB2 = 0.4, NPY = 0.3
  ),

  "Serotonergic Neurons" = c(
    TPH2 = 1.0, SLC6A4 = 0.9, FEV = 0.8, DDC = 0.7, GCHFR = 0.6,
    SLC18A2 = 0.6, HTR1A = 0.5, HTR1B = 0.5, MAOA = 0.4, SYN1 = 0.4,
    RBFOX3 = 0.3, SNAP25 = 0.3
  ),

  "Motor Neurons" = c(
    CHAT = 1.0, SLC5A7 = 0.9, ISL1 = 0.8, ISL2 = 0.7, MNX1 = 0.8,
    NEFL = 0.6, NEFM = 0.5, PRPH = 0.5, SLC18A3 = 0.6, RET = 0.4,
    LHX3 = 0.4, OLIG2 = 0.3
  ),

  "Astrocyte" = c(
    GFAP = 1.0, S100B = 0.9, AQP4 = 0.8, ALDH1L1 = 0.8, SLC1A3 = 0.7,
    SLC1A2 = 0.7, GJA1 = 0.6, GLUL = 0.6, SOX9 = 0.5, NDRG2 = 0.5,
    CLU = 0.4, FGFR3 = 0.4
  ),

  "Oligodendrocyte" = c(
    MBP = 1.0, MOG = 0.9, PLP1 = 0.9, OLIG2 = 0.8, MAG = 0.7,
    CNP = 0.7, CLDN11 = 0.6, GALC = 0.5, UGT8 = 0.5, MOBP = 0.6,
    TF = 0.4, FA2H = 0.4
  ),

  "OPCs (Oligodendrocyte Precursors)" = c(
    PDGFRA = 1.0, CSPG4 = 0.9, OLIG1 = 0.8, OLIG2 = 0.8, SOX10 = 0.7,
    NKX2_2 = 0.6, GPR17 = 0.6, PCDH15 = 0.5, VCAN = 0.5, PTPRZ1 = 0.5,
    BCAN = 0.4, NEU4 = 0.3
  ),

  "Ependymal Cells" = c(
    FOXJ1 = 1.0, CFAP126 = 0.8, PIFO = 0.7, TPPP3 = 0.7, DNAH12 = 0.6,
    CCDC153 = 0.6, TMEM212 = 0.6, S100B = 0.5, CD24 = 0.5, VIM = 0.4,
    SOX2 = 0.4, RARRES2 = 0.3
  ),

  "Schwann cell" = c(
    MPZ = 1.0, PMP22 = 0.9, SOX10 = 0.8, PLP1 = 0.6, MBP = 0.6,
    EGR2 = 0.7, PRX = 0.6, NCAM1 = 0.4, L1CAM = 0.5, CDH19 = 0.5,
    DHH = 0.4, NGFR = 0.4
  ),

  "Satellite Glial Cells" = c(
    FABP7 = 1.0, SOX10 = 0.8, SOX2 = 0.7, CDH19 = 0.7, S100B = 0.6,
    PLP1 = 0.5, GAP43 = 0.5, GLUL = 0.5, KCNJ10 = 0.4, GJB1 = 0.4,
    GDNF = 0.3, SCN7A = 0.3
  ),

  "Radial Glia" = c(
    VIM = 1.0, NES = 0.9, SOX2 = 0.8, HES1 = 0.8, PAX6 = 0.7,
    FABP7 = 0.7, GFAP = 0.5, TNC = 0.6, PTPRZ1 = 0.5, GLI3 = 0.4,
    EMX2 = 0.4, HES5 = 0.4
  ),

  "Neural Progenitor Cell" = c(
    SOX2 = 1.0, PAX6 = 0.9, NES = 0.8, HES1 = 0.7, HES5 = 0.6,
    ASCL1 = 0.6, DCX = 0.5, NEUROD1 = 0.4, DLL1 = 0.5, NOTCH1 = 0.4,
    SOX1 = 0.5, MKI67 = 0.3
  ),

  # ============================================================================
  # STROMAL / MESENCHYMAL LINEAGE
  # ============================================================================

  "Fibroblast" = c(
    COL1A1 = 1.0, DCN = 0.9, LUM = 0.8, PDGFRA = 0.7, COL1A2 = 0.8,
    COL3A1 = 0.7, VIM = 0.5, THY1 = 0.5, FAP = 0.4, FBLN1 = 0.5,
    COL6A1 = 0.4, MFAP5 = 0.4
  ),

  "Myofibroblast" = c(
    ACTA2 = 1.0, TAGLN = 0.9, MYL9 = 0.8, COL1A1 = 0.7, COL1A2 = 0.6,
    MYH11 = 0.5, CNN1 = 0.5, POSTN = 0.6, SPARC = 0.4, CTGF = 0.5,
    TPM2 = 0.4, CALD1 = 0.4
  ),

  "Cancer-Associated Fibroblasts (CAFs)" = c(
    FAP = 1.0, ACTA2 = 0.8, PDPN = 0.7, PDGFRA = 0.7, COL1A1 = 0.7,
    COL1A2 = 0.6, POSTN = 0.7, CXCL12 = 0.6, IL6 = 0.5, MMP2 = 0.5,
    COMP = 0.4, THY1 = 0.4, TNC = 0.5
  ),

  "Mesenchymal Stem Cell" = c(
    LEPR = 1.0, NT5E = 0.9, THY1 = 0.8, ENG = 0.7, VCAM1 = 0.6,
    ALCAM = 0.5, ITGB1 = 0.5, CD44 = 0.4, PDGFRB = 0.5, CXCL12 = 0.5,
    KITLG = 0.4, ANGPT1 = 0.3
  ),

  "Adipocyte" = c(
    ADIPOQ = 1.0, PLIN1 = 0.9, LEP = 0.8, FABP4 = 0.8, PPARG = 0.7,
    LPL = 0.7, CIDEC = 0.6, PCK1 = 0.5, LIPE = 0.5, PLIN4 = 0.4,
    CFD = 0.5, SLC2A4 = 0.3
  ),

  "Osteoblast" = c(
    RUNX2 = 1.0, SP7 = 0.9, BGLAP = 0.8, ALPL = 0.8, COL1A1 = 0.6,
    IBSP = 0.7, SPP1 = 0.6, DMP1 = 0.5, PHEX = 0.5, IFITM5 = 0.4,
    SOST = 0.3, BMP2 = 0.3
  ),

  "Osteoclast" = c(
    ACP5 = 1.0, CTSK = 0.9, MMP9 = 0.8, DCSTAMP = 0.7, OSCAR = 0.7,
    TNFRSF11A = 0.6, ATP6V0D2 = 0.6, CALCR = 0.5, CA2 = 0.4, ITGB3 = 0.4,
    NFATc1 = 0.5, TRAP = 0.3
  ),

  "Chondrocyte" = c(
    COL2A1 = 1.0, ACAN = 0.9, SOX9 = 0.9, COMP = 0.7, COL11A1 = 0.7,
    COL9A1 = 0.6, MATN3 = 0.5, GDF5 = 0.5, PRG4 = 0.5, HAPLN1 = 0.4,
    FRZB = 0.4, COL11A2 = 0.3
  ),

  "Tenocytes" = c(
    SCX = 1.0, TNMD = 0.9, MKX = 0.8, COL1A1 = 0.7, COL1A2 = 0.6,
    DCN = 0.6, FMOD = 0.6, TNC = 0.5, THBS4 = 0.5, COMP = 0.4,
    COL3A1 = 0.4, BGN = 0.3
  ),

  "Pericyte" = c(
    RGS5 = 1.0, PDGFRB = 0.9, ACTA2 = 0.7, CSPG4 = 0.7, MCAM = 0.6,
    NOTCH3 = 0.6, DES = 0.5, ABCC9 = 0.5, KCNJ8 = 0.5, DLC1 = 0.4,
    VTN = 0.4, HIGD1B = 0.3
  ),

  # ============================================================================
  # ENDOTHELIAL LINEAGE
  # ============================================================================

  "Endothelial (General)" = c(
    PECAM1 = 1.0, CDH5 = 1.0, VWF = 0.8, CLDN5 = 0.8, ERG = 0.7,
    KDR = 0.6, FLT1 = 0.6, EMCN = 0.5, ENG = 0.5, TIE1 = 0.4,
    TEK = 0.4, ESAM = 0.4
  ),

  "Vascular Endothelial Cells" = c(
    PECAM1 = 1.0, CDH5 = 0.9, VWF = 0.8, CLDN5 = 0.7, FLT1 = 0.7,
    KDR = 0.6, ERG = 0.6, ENG = 0.5, NOS3 = 0.5, TEK = 0.4,
    EMCN = 0.4, ICAM2 = 0.3
  ),

  "Lymphatic Endothelial Cells" = c(
    PROX1 = 1.0, PDPN = 0.9, LYVE1 = 0.9, CCL21 = 0.7, FLT4 = 0.8,
    MMRN1 = 0.6, TBX1 = 0.5, FOXC2 = 0.5, SOX18 = 0.4, NRP2 = 0.5,
    VEGFC = 0.4, CLDN5 = 0.3
  ),

  "Arterial Endothelial" = c(
    BMX = 1.0, HEY1 = 0.8, GJA5 = 0.8, EFNB2 = 0.7, GJA4 = 0.6,
    DLL4 = 0.6, SOX17 = 0.6, NOTCH4 = 0.5, CXCR4 = 0.4, CDH5 = 0.4,
    PECAM1 = 0.4, SEMA3G = 0.5
  ),

  "Venous Endothelial" = c(
    ACKR1 = 1.0, NR2F2 = 0.9, SELP = 0.7, EPHB4 = 0.6, VWF = 0.5,
    PLVAP = 0.5, CLU = 0.4, PECAM1 = 0.4, CDH5 = 0.4, SELE = 0.4,
    VCAM1 = 0.3, ICAM1 = 0.3
  ),

  "Capillary Endothelial" = c(
    CA4 = 1.0, RGCC = 0.8, PLVAP = 0.7, FABP4 = 0.6, SLC9A3R2 = 0.5,
    PECAM1 = 0.5, CDH5 = 0.5, CLDN5 = 0.4, VWF = 0.4, FCN3 = 0.4,
    SPARCL1 = 0.4, CD36 = 0.3
  ),

  # ============================================================================
  # MUSCLE LINEAGE
  # ============================================================================

  "Smooth Muscle Cell" = c(
    ACTA2 = 1.0, MYH11 = 1.0, TAGLN = 0.9, CNN1 = 0.8, DES = 0.7,
    MYLK = 0.7, MYL9 = 0.6, CALD1 = 0.6, LMOD1 = 0.5, SMTN = 0.5,
    ACTG2 = 0.5, TPM2 = 0.4
  ),

  "Skeletal Muscle Cell" = c(
    MYH1 = 1.0, DES = 0.9, ACTN2 = 0.8, TNNT1 = 0.7, MYH2 = 0.7,
    MYH7 = 0.6, TTN = 0.6, MYBPC1 = 0.5, NEB = 0.5, CKM = 0.5,
    ACTA1 = 0.5, TNNI2 = 0.4
  ),

  "Cardiac Muscle Cell" = c(
    TNNT2 = 1.0, MYH7 = 0.9, NPPA = 0.8, MYL2 = 0.7, MYH6 = 0.7,
    ACTC1 = 0.6, MB = 0.6, RYR2 = 0.5, KCNJ2 = 0.5, SCN5A = 0.4,
    GJA1 = 0.4, PLN = 0.4
  ),

  # ============================================================================
  # KIDNEY-SPECIFIC
  # ============================================================================

  "Proximal Tubular Cells" = c(
    SLC34A1 = 1.0, LRP2 = 0.9, CUBN = 0.8, SLC5A12 = 0.7, ALDOB = 0.7,
    SLC22A8 = 0.6, GATM = 0.6, SLC13A3 = 0.5, SLC5A2 = 0.5, ACSM2A = 0.5,
    HNF4A = 0.4, NAT8 = 0.4
  ),

  "Distal Tubular Cells" = c(
    SLC12A3 = 1.0, TRPM6 = 0.8, PVALB = 0.7, CALB1 = 0.7, SLC8A1 = 0.6,
    EGF = 0.6, CLDN8 = 0.5, UMOD = 0.5, ESRRB = 0.5, WNK1 = 0.4,
    KLHL3 = 0.4, FXYD2 = 0.4
  ),

  "Collecting Duct Principal Cells" = c(
    AQP2 = 1.0, AQP3 = 0.9, FXYD4 = 0.8, SCNN1B = 0.7, SCNN1G = 0.7,
    AQP4 = 0.6, AVPR2 = 0.6, SLC14A2 = 0.5, HCRT = 0.4, GATA3 = 0.4,
    CALB1 = 0.4, FOXI1 = -0.5
  ),

  "Collecting Duct Intercalated Cells" = c(
    SLC4A1 = 0.9, ATP6V1B1 = 1.0, ATP6V0D2 = 0.8, SLC26A4 = 0.7, CA2 = 0.7,
    FOXI1 = 0.8, DMRT2 = 0.5, AQP6 = 0.5, SLC4A9 = 0.5, KIT = 0.4,
    INSRR = 0.4, CLNK = 0.3
  ),

  "Loop of Henle Cells" = c(
    UMOD = 1.0, SLC12A1 = 0.9, CLDN16 = 0.8, CLDN19 = 0.7, CASR = 0.7,
    KCNJ1 = 0.6, SLC12A2 = 0.5, CLCNKB = 0.5, PTH1R = 0.4, IRXL1 = 0.4,
    EGF = 0.3, SLC44A5 = 0.3
  ),

  "Podocytes" = c(
    NPHS1 = 1.0, NPHS2 = 1.0, WT1 = 0.9, SYNPO = 0.8, PODXL = 0.7,
    PTPRO = 0.7, TCF21 = 0.6, MAFB = 0.5, VEGFA = 0.5, CLIC5 = 0.4,
    CDH3 = 0.4, NTNG1 = 0.3
  ),

  "Mesangial Cells" = c(
    PDGFRB = 1.0, GATA3 = 0.8, ITGA8 = 0.7, DES = 0.6, CNN1 = 0.5,
    ACTA2 = 0.5, REN = 0.4, TAGLN = 0.4, COL4A5 = 0.4, CSPG4 = 0.3,
    CFH = 0.3, COL4A4 = 0.3
  ),

  # ============================================================================
  # LIVER-SPECIFIC
  # ============================================================================

  "Liver Sinusoidal Endothelial Cells (LSECs)" = c(
    CLEC4G = 1.0, CLEC4M = 0.9, FCN2 = 0.8, STAB2 = 0.8, LYVE1 = 0.7,
    PECAM1 = 0.5, CDH5 = 0.5, DNASE1L3 = 0.6, OIT3 = 0.5, CLEC1B = 0.5,
    CD32B = 0.4, PLVAP = 0.4
  ),

  "Hepatic Stellate Cells" = c(
    PDGFRA = 0.9, COL1A1 = 0.8, RGS5 = 0.7, LRAT = 0.9, HGF = 0.6,
    DES = 0.5, ACTA2 = 0.5, COL3A1 = 0.5, HAND2 = 0.6, RELN = 0.5,
    CTGF = 0.4, RBP1 = 0.4
  ),

  "Stellate cell (Liver)" = c(
    PDGFRA = 0.9, COL1A1 = 0.8, RGS5 = 0.7, LRAT = 0.9, HGF = 0.6,
    DES = 0.5, ACTA2 = 0.5, COL3A1 = 0.5, HAND2 = 0.6, RELN = 0.5,
    CTGF = 0.4, RBP1 = 0.4
  ),

  # ============================================================================
  # LUNG-SPECIFIC
  # ============================================================================

  "Airway Basal Cells" = c(
    KRT5 = 1.0, KRT14 = 0.9, TP63 = 0.9, KRT17 = 0.7, KRT15 = 0.6,
    NGFR = 0.5, ITGA6 = 0.6, ITGB4 = 0.5, S100A2 = 0.5, DLK2 = 0.4,
    DAPL1 = 0.4, IL33 = 0.3
  ),

  "Airway Smooth Muscle" = c(
    ACTA2 = 1.0, MYH11 = 0.9, TAGLN = 0.8, DES = 0.7, CNN1 = 0.7,
    MYLK = 0.6, ACTG2 = 0.6, MYL9 = 0.5, CALD1 = 0.5, LMOD1 = 0.4,
    SMTN = 0.4, TPM2 = 0.3
  ),

  "Pulmonary Neuroendocrine Cells" = c(
    CHGA = 1.0, SYP = 0.9, ASCL1 = 0.8, GRP = 0.8, CALCA = 0.7,
    INSM1 = 0.7, NCAM1 = 0.6, SYN1 = 0.5, CHGB = 0.5, ENO2 = 0.5,
    DDC = 0.4, SCG2 = 0.3
  ),

  "Ionocytes" = c(
    FOXI1 = 1.0, ASCL3 = 0.9, CFTR = 0.8, ATP6V1C2 = 0.7, ATP6V0D2 = 0.6,
    BSND = 0.6, CLCNKB = 0.5, DMRT2 = 0.5, SLC26A4 = 0.4, FOXI2 = 0.4,
    INSRR = 0.4, KIT = 0.3
  ),

  # ============================================================================
  # HEART-SPECIFIC
  # ============================================================================

  "Atrial Cardiomyocytes" = c(
    NPPA = 1.0, NPPB = 0.8, MYH6 = 0.9, MYL7 = 0.8, TNNT2 = 0.7,
    GJA5 = 0.6, NR2F1 = 0.5, TBX5 = 0.5, HCN4 = 0.4, PITX2 = 0.5,
    KCNA5 = 0.4, SCN5A = 0.3
  ),

  "Ventricular Cardiomyocytes" = c(
    MYH7 = 1.0, MYL2 = 0.9, TNNT2 = 0.8, ACTC1 = 0.7, RYR2 = 0.7,
    MB = 0.6, IRX4 = 0.6, HAND1 = 0.5, HAND2 = 0.5, PLN = 0.5,
    HEY2 = 0.4, KCNJ2 = 0.4
  ),

  "Cardiac Fibroblasts" = c(
    COL1A1 = 1.0, DCN = 0.9, TCF21 = 0.8, PDGFRA = 0.7, COL3A1 = 0.7,
    LUM = 0.6, POSTN = 0.6, DDR2 = 0.5, FBLN1 = 0.5, GSN = 0.4,
    MFAP5 = 0.4, COL1A2 = 0.5
  ),

  "Endocardial Cells" = c(
    PECAM1 = 0.9, CDH5 = 0.8, NPR3 = 0.8, EMCN = 0.7, NFATC1 = 0.7,
    BMP2 = 0.5, EDN1 = 0.5, VWF = 0.5, TBX20 = 0.4, NOTCH1 = 0.4,
    CDH11 = 0.4, CYTL1 = 0.4
  ),

  # ============================================================================
  # GUT / INTESTINE-SPECIFIC
  # ============================================================================

  "Intestinal Stem Cells (LGR5+)" = c(
    LGR5 = 1.0, ASCL2 = 0.9, OLFM4 = 0.8, SOX9 = 0.7, EPHB2 = 0.7,
    SMOC2 = 0.6, RGMB = 0.5, CD44 = 0.5, MYC = 0.4, AXIN2 = 0.4,
    RNF43 = 0.4, ZNRF3 = 0.3
  ),

  "Transit Amplifying Cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, PCNA = 0.8, STMN1 = 0.7, HMGB2 = 0.6,
    TUBA1B = 0.5, UBE2C = 0.5, BIRC5 = 0.5, CENPF = 0.4, CDK1 = 0.4,
    TPX2 = 0.4, NUSAP1 = 0.3
  ),

  "Enterocytes" = c(
    APOA1 = 1.0, FABP1 = 0.9, VIL1 = 0.8, FABP2 = 0.8, SLC5A1 = 0.7,
    ALPI = 0.7, RBP2 = 0.6, SI = 0.6, ACE2 = 0.5, ANPEP = 0.5,
    SLC26A3 = 0.4, SLC15A1 = 0.4
  ),

  "Goblet Cells (Intestinal)" = c(
    MUC2 = 1.0, TFF3 = 0.9, CLCA1 = 0.8, SPDEF = 0.8, ZG16 = 0.7,
    FCGBP = 0.7, AGR2 = 0.6, ITLN1 = 0.6, ATOH1 = 0.5, REP15 = 0.4,
    KLK1 = 0.4, LRRC26 = 0.3
  ),

  "Paneth Cells" = c(
    LYZ = 1.0, DEFA5 = 0.9, DEFA6 = 0.9, MMP7 = 0.7, REG3A = 0.7,
    REG3G = 0.6, ITLN2 = 0.5, PLA2G2A = 0.5, WNT3 = 0.5, SOX9 = 0.4,
    EPHB3 = 0.4, CA7 = 0.3
  ),

  "Enteroendocrine Cells" = c(
    CHGA = 1.0, CHGB = 0.9, SYP = 0.8, NEUROD1 = 0.8, PAX4 = 0.6,
    INSM1 = 0.7, GCG = 0.5, GIP = 0.5, CCK = 0.5, SST = 0.4,
    TPH1 = 0.5, SCG5 = 0.4
  ),

  "M Cells" = c(
    GP2 = 1.0, SPIB = 0.9, CCL20 = 0.7, TNFAIP2 = 0.6, MARCKSL1 = 0.5,
    SPI1 = 0.4, ANXA5 = 0.5, SOX8 = 0.5, BTNL8 = 0.4, CCL9 = 0.4,
    NFKBIZ = 0.3, PGRP = 0.3
  ),

  "Tuft Cells" = c(
    DCLK1 = 1.0, TRPM5 = 0.9, POU2F3 = 0.9, HCK = 0.7, GFI1B = 0.6,
    SH2D6 = 0.6, ALOX5AP = 0.5, BMX = 0.5, AVIL = 0.5, GNG13 = 0.4,
    IL25 = 0.4, LRMP = 0.3
  ),

  # ============================================================================
  # SKIN-SPECIFIC
  # ============================================================================

  "Keratinocytes (Basal)" = c(
    KRT5 = 1.0, KRT14 = 0.9, TP63 = 0.9, COL17A1 = 0.7, ITGA6 = 0.7,
    ITGB4 = 0.6, DST = 0.5, KRT15 = 0.6, S100A2 = 0.5, LAMB3 = 0.4,
    LAMC2 = 0.4, LAMA3 = 0.3
  ),

  "Keratinocytes (Suprabasal)" = c(
    KRT1 = 1.0, KRT10 = 0.9, IVL = 0.8, FLG = 0.7, LOR = 0.7,
    SPRR1B = 0.6, SPRR2A = 0.5, DSG1 = 0.6, DSC1 = 0.5, TGM1 = 0.5,
    CDSN = 0.4, CASP14 = 0.3
  ),

  "Melanocytes" = c(
    PMEL = 1.0, MLANA = 0.9, TYR = 0.8, DCT = 0.7, TYRP1 = 0.7,
    MITF = 0.8, SOX10 = 0.6, S100B = 0.5, KIT = 0.5, EDNRB = 0.4,
    SLC45A2 = 0.4, CDH19 = 0.3
  ),

  "Merkel Cells" = c(
    ATOH1 = 1.0, KRT20 = 0.9, ISL1 = 0.8, SOX2 = 0.7, CHGA = 0.6,
    SYP = 0.5, NCAM1 = 0.5, PIEZO2 = 0.6, ENO2 = 0.4, KCNC1 = 0.4,
    CDH13 = 0.3, ASCL1 = 0.3
  ),

  "Langerhans Cells (Skin)" = c(
    CD207 = 1.0, CD1A = 0.9, EPCAM = 0.7, FCER1A = 0.6, CD1C = 0.5,
    HLA_DRA = 0.5, ITGAX = 0.4, CD1E = 0.5, S100B = 0.4, AXL = 0.3,
    BIRC3 = 0.3, TREM2 = 0.3
  ),

  "Dermal Papilla Cells" = c(
    SOX2 = 1.0, CRABP1 = 0.9, LEF1 = 0.8, WNT5A = 0.7, WIF1 = 0.7,
    VCAN = 0.6, CORIN = 0.6, BMP4 = 0.5, ALPL = 0.5, HHIP = 0.4,
    PTCH1 = 0.4, NOG = 0.3
  ),

  # ============================================================================
  # EYE / RETINA-SPECIFIC
  # ============================================================================

  "Rod Photoreceptors" = c(
    RHO = 1.0, NRL = 0.9, NR2E3 = 0.8, PDE6A = 0.7, PDE6B = 0.7,
    GNAT1 = 0.8, CNGA1 = 0.6, RCVRN = 0.5, ROM1 = 0.5, SAG = 0.5,
    RBP3 = 0.4, GNGT1 = 0.4
  ),

  "Cone Photoreceptors" = c(
    ARR3 = 1.0, OPN1SW = 0.8, OPN1MW = 0.8, OPN1LW = 0.7, GNAT2 = 0.8,
    PDE6C = 0.7, PDE6H = 0.6, CNGA3 = 0.6, GRK7 = 0.5, THRB = 0.5,
    GNGT2 = 0.4, RXRG = 0.4
  ),

  "Retinal Ganglion Cells" = c(
    POU4F2 = 1.0, RBPMS = 0.9, SLC17A6 = 0.8, ISL1 = 0.7, GAP43 = 0.7,
    NEFL = 0.6, NEFM = 0.5, THY1 = 0.5, SNCG = 0.6, BRN3B = 0.4,
    POU4F1 = 0.5, ATOH7 = 0.3
  ),

  "Muller Glia" = c(
    RLBP1 = 1.0, GLUL = 0.9, VIM = 0.8, CLU = 0.7, SOX9 = 0.7,
    SOX2 = 0.6, HES1 = 0.5, CA2 = 0.5, SLC1A3 = 0.6, DIO3 = 0.4,
    APOE = 0.4, CRALBP = 0.3
  ),

  "Retinal Pigment Epithelium" = c(
    RPE65 = 1.0, RLBP1 = 0.8, BEST1 = 0.8, LRAT = 0.7, TYR = 0.7,
    TYRP1 = 0.6, MITF = 0.6, TTR = 0.5, RDH5 = 0.6, MERTK = 0.5,
    RGR = 0.4, CRALBP = 0.4
  ),

  "Bipolar Cells (Retina)" = c(
    VSX2 = 1.0, VSX1 = 0.8, GRM6 = 0.7, OTX2 = 0.7, TRPM1 = 0.6,
    PRKCA = 0.6, CABP5 = 0.6, ISL1 = 0.5, LHFPL3 = 0.4, BHLHE23 = 0.5,
    SCGN = 0.4, TMEM215 = 0.3
  ),

  "Amacrine Cells" = c(
    GAD1 = 0.9, SLC6A9 = 0.8, PAX6 = 0.7, TFAP2A = 0.8, SLC32A1 = 0.7,
    CALB1 = 0.6, GAD2 = 0.6, CHAT = 0.5, TH = 0.4, NHLH2 = 0.5,
    MEIS2 = 0.4, SLC6A1 = 0.3
  ),

  # ============================================================================
  # DEVELOPMENT / REPRODUCTIVE
  # ============================================================================

  "Trophoblasts" = c(
    KRT7 = 1.0, TFAP2A = 0.8, GATA3 = 0.7, KRT19 = 0.7, EPCAM = 0.6,
    HLA_G = 0.8, TEAD4 = 0.6, CDX2 = 0.5, GCM1 = 0.5, TFAP2C = 0.5,
    TP63 = 0.4, ELF5 = 0.4
  ),

  "Syncytiotrophoblasts" = c(
    CGA = 1.0, CGB = 0.9, CSH1 = 0.8, CSH2 = 0.7, GH2 = 0.6,
    CYP19A1 = 0.7, SDC1 = 0.6, GCM1 = 0.6, ERVW_1 = 0.5, LGALS13 = 0.5,
    PSG1 = 0.5, HSD3B1 = 0.4
  ),

  "Extravillous Trophoblasts" = c(
    HLA_G = 1.0, ITGA5 = 0.8, ITGA1 = 0.8, MMP2 = 0.7, NOTCH2 = 0.6,
    FLT1 = 0.6, PAPPA2 = 0.6, ADAM12 = 0.5, PLAC8 = 0.5, ASCL2 = 0.5,
    KRT7 = 0.4, CDH5 = 0.3
  ),

  "Decidual Cells" = c(
    PRL = 1.0, IGFBP1 = 0.9, FOXO1 = 0.7, WNT4 = 0.6, BMP2 = 0.6,
    SST = 0.5, HAND2 = 0.5, VIM = 0.4, DKK1 = 0.5, MAOA = 0.5,
    LEFTY2 = 0.4, SCARA5 = 0.4
  ),

  "Sertoli Cells" = c(
    SOX9 = 1.0, AMH = 0.9, WT1 = 0.8, GATA4 = 0.7, CLU = 0.6,
    INHA = 0.6, INHBB = 0.5, FST = 0.5, CLDN11 = 0.5, DMRT1 = 0.5,
    DHH = 0.4, FSHR = 0.4
  ),

  "Leydig Cells" = c(
    CYP11A1 = 1.0, CYP17A1 = 0.9, HSD3B1 = 0.8, STAR = 0.8, INSL3 = 0.7,
    HSD17B3 = 0.7, NR5A1 = 0.6, LHCGR = 0.6, DLK1 = 0.5, IGF1 = 0.4,
    SULT2A1 = 0.4, TSPO = 0.3
  ),

  # ============================================================================
  # TUMOR SIGNATURES
  # ============================================================================

  "Proliferating Tumor Cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, PCNA = 0.8, BIRC5 = 0.8, CDK1 = 0.7,
    CCNB1 = 0.7, UBE2C = 0.6, CENPF = 0.6, NUSAP1 = 0.5, TPX2 = 0.5,
    STMN1 = 0.4, TYMS = 0.4
  ),

  "Hypoxic Tumor Cells" = c(
    HIF1A = 1.0, VEGFA = 0.9, CA9 = 0.9, SLC2A1 = 0.8, LDHA = 0.7,
    PGK1 = 0.7, BNIP3 = 0.6, ENO1 = 0.6, ADM = 0.5, NDRG1 = 0.5,
    LOX = 0.5, EGLN3 = 0.4
  ),

  "Mesenchymal Tumor Cells" = c(
    VIM = 1.0, CDH2 = 0.9, FN1 = 0.8, TWIST1 = 0.8, SNAI1 = 0.7,
    SNAI2 = 0.7, ZEB1 = 0.6, ZEB2 = 0.6, MMP2 = 0.5, SPARC = 0.5,
    ACTA2 = 0.4, ITGA5 = 0.4
  ),

  "Epithelial Tumor Cells" = c(
    EPCAM = 1.0, CDH1 = 0.9, KRT8 = 0.8, KRT18 = 0.8, KRT19 = 0.7,
    MUC1 = 0.6, CLDN4 = 0.6, CLDN7 = 0.5, ERBB2 = 0.5, TACSTD2 = 0.5,
    GRHL2 = 0.4, CD24 = 0.4
  ),

  "Cancer Stem Cells" = c(
    PROM1 = 1.0, ALDH1A1 = 0.9, CD44 = 0.8, SOX2 = 0.7, NANOG = 0.7,
    POU5F1 = 0.6, BMI1 = 0.6, ABCG2 = 0.6, NOTCH1 = 0.5, WNT1 = 0.4,
    LGR5 = 0.4, KLF4 = 0.4
  ),

  # ============================================================================
  # PROGENITOR / PRECURSOR CELLS
  # ============================================================================

  "Megakaryocyte-Erythroid Progenitors (MEPs)" = c(
    GATA1 = 1.0, KLF1 = 0.9, CD34 = 0.8, ITGA2B = 0.8, GYPA = 0.7,
    TAL1 = 0.7, GFI1B = 0.6, NFE2 = 0.6, ZFPM1 = 0.5, GATA2 = 0.5,
    LMO2 = 0.4, EPOR = 0.4
  ),

  "Common Lymphoid Progenitors (CLPs)" = c(
    CD34 = 1.0, IL7R = 0.9, FLT3 = 0.8, DNTT = 0.8, MME = 0.7,
    EBF1 = 0.6, IGLL1 = 0.6, RAG1 = 0.5, RAG2 = 0.5, PTPRC = 0.5,
    CXCR4 = 0.4, CD38 = 0.4
  ),

  "Pre-DCs / DC Precursors" = c(
    CLEC9A = 1.0, IRF8 = 0.9, BATF3 = 0.8, CD34 = 0.7, FLT3 = 0.7,
    IRF4 = 0.6, TCF4 = 0.6, ID2 = 0.5, ZBTB46 = 0.5, AXL = 0.5,
    CX3CR1 = 0.4, SIGLEC6 = 0.4
  ),

  # ============================================================================
  # TISSUE-SPECIFIC CELLS
  # ============================================================================

  "Mesothelial Cells" = c(
    MSLN = 1.0, WT1 = 0.9, KRT19 = 0.8, UPK3B = 0.8, CALB2 = 0.7,
    HP = 0.6, LRRN4 = 0.6, BNC1 = 0.5, PRG4 = 0.5, ITLN1 = 0.5,
    GPM6A = 0.4, PDPN = 0.4
  ),

  "Interstitial Cells of Cajal" = c(
    KIT = 1.0, ANO1 = 0.9, ETV1 = 0.8, PRKCQ = 0.7, CD34 = 0.6,
    SLC4A4 = 0.5, THBS4 = 0.5, GPR20 = 0.5, LRIG1 = 0.4, DES = -0.5
  ),

  # ============================================================================
  # CYCLING / PROLIFERATING STATES
  # ============================================================================

  "Cycling T-cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, CD3D = 0.8, CD3E = 0.8, STMN1 = 0.7,
    BIRC5 = 0.7, UBE2C = 0.6, CENPF = 0.6, TYMS = 0.5, TRAC = 0.5,
    CDK1 = 0.5, PCNA = 0.4
  ),

  "Cycling B-cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, MS4A1 = 0.8, CD79A = 0.8, STMN1 = 0.7,
    BIRC5 = 0.7, UBE2C = 0.6, CENPF = 0.6, TYMS = 0.5, CD79B = 0.5,
    CDK1 = 0.5, PCNA = 0.4
  ),

  "Cycling NK Cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, NKG7 = 0.8, GNLY = 0.8, STMN1 = 0.7,
    BIRC5 = 0.7, UBE2C = 0.6, CENPF = 0.6, TYMS = 0.5, KLRD1 = 0.5,
    CDK1 = 0.5, NCAM1 = 0.4
  ),

  # ============================================================================
  # LUNG / AIRWAY EPITHELIUM
  # ============================================================================

  "Alveolar Type 1 Cells (AT1)" = c(
    AGER = 1.0, PDPN = 0.9, CAV1 = 0.8, AQP5 = 0.8, CLIC5 = 0.7,
    EMP2 = 0.7, RTKN2 = 0.6, HOPX = 0.6, IGFBP2 = 0.5, AKAP5 = 0.5,
    SPOCK2 = 0.4, CEACAM6 = 0.3
  ),

  "Alveolar Type 2 Cells (AT2)" = c(
    SFTPC = 1.0, SFTPB = 0.9, ABCA3 = 0.8, SFTPA1 = 0.8, SFTPA2 = 0.7,
    SFTPD = 0.7, SLC34A2 = 0.6, NAPSA = 0.6, PGC = 0.5, LAMP3 = 0.5,
    NKX2_1 = 0.5, ETV5 = 0.4
  ),

  "Club Cells" = c(
    SCGB1A1 = 1.0, SCGB3A1 = 0.9, SCGB3A2 = 0.8, CYP2F1 = 0.7,
    HP = 0.6, MGP = 0.6, BPIFB1 = 0.5, LCN2 = 0.5, KRT15 = 0.4,
    SOX2 = 0.3, NOTCH3 = 0.3
  ),

  "Goblet Cells" = c(
    MUC5AC = 1.0, MUC5B = 0.9, SPDEF = 0.8, TFF3 = 0.7, AGR2 = 0.7,
    FCGBP = 0.6, MUC2 = 0.6, CLCA1 = 0.5, ITLN1 = 0.5, REP15 = 0.4,
    GOB5 = 0.4, ATOH1 = 0.3
  ),

  "Ionocytes" = c(
    FOXI1 = 1.0, CFTR = 0.9, ATP6V1C2 = 0.8, ATP6V0D2 = 0.7,
    ASCL3 = 0.7, CLCNKA = 0.6, TFCP2L1 = 0.6, BSND = 0.5,
    SLC26A4 = 0.5, DMRT2 = 0.4
  ),

  "Tuft Cells" = c(
    POU2F3 = 1.0, DCLK1 = 0.9, TRPM5 = 0.8, GFI1B = 0.7, SH2D6 = 0.7,
    LRMP = 0.6, IL25 = 0.6, BMX = 0.5, HCK = 0.5, PLCG2 = 0.5,
    AVIL = 0.4, GNG13 = 0.4
  ),

  # ============================================================================
  # SKIN / EPITHELIAL
  # ============================================================================

  "Keratinocytes (Basal)" = c(
    KRT5 = 1.0, KRT14 = 0.9, TP63 = 0.8, COL17A1 = 0.7, ITGA6 = 0.7,
    ITGB4 = 0.6, DST = 0.6, LAMB3 = 0.5, KRT15 = 0.5, MKI67 = 0.4,
    CDH3 = 0.4, S100A2 = 0.3
  ),

  "Keratinocytes (Suprabasal)" = c(
    KRT1 = 1.0, KRT10 = 0.9, IVL = 0.8, FLG = 0.7, LOR = 0.7,
    SPRR1A = 0.6, SPRR1B = 0.6, TGM1 = 0.5, DSP = 0.5, SBSN = 0.5,
    KRT2 = 0.4, CASP14 = 0.3
  ),

  # ============================================================================
  # NEURAL / GLIAL
  # ============================================================================

  "Schwann Cells" = c(
    MPZ = 1.0, PMP22 = 0.9, SOX10 = 0.8, S100B = 0.8, MBP = 0.7,
    PRX = 0.7, EGR2 = 0.6, NGFR = 0.6, NCAM1 = 0.5, CDH19 = 0.5,
    GAP43 = 0.4, L1CAM = 0.4
  ),

  "Satellite Glial Cells" = c(
    FABP7 = 1.0, CDH19 = 0.9, SOX2 = 0.7, S100B = 0.7, KCNJ10 = 0.6,
    GJA1 = 0.6, BCAN = 0.5, GLUL = 0.5, NTRK2 = 0.5, SOX10 = 0.4,
    PLP1 = 0.4, GFAP = 0.3
  ),

  # ============================================================================
  # KIDNEY
  # ============================================================================

  "Podocytes" = c(
    NPHS1 = 1.0, NPHS2 = 0.9, WT1 = 0.8, PODXL = 0.8, SYNPO = 0.7,
    CD2AP = 0.7, PTPRO = 0.6, TCF21 = 0.6, MAFB = 0.5, TJP1 = 0.5,
    VEGFA = 0.4, CLIC5 = 0.4
  ),

  # ============================================================================
  # LIVER
  # ============================================================================

  "Hepatocytes" = c(
    ALB = 1.0, APOB = 0.9, CYP3A4 = 0.8, HNF4A = 0.8, TF = 0.7,
    TTR = 0.7, SERPINA1 = 0.6, CYP2E1 = 0.6, FGA = 0.5, FGB = 0.5,
    ASGR1 = 0.5, PCK1 = 0.4
  ),

  "Cholangiocytes" = c(
    KRT19 = 1.0, KRT7 = 0.9, EPCAM = 0.8, SOX9 = 0.8, MUC1 = 0.7,
    CFTR = 0.6, SCTR = 0.6, SPP1 = 0.5, HNF1B = 0.5, CLDN4 = 0.5,
    ONECUT1 = 0.4, ANXA4 = 0.4
  ),

  # ============================================================================
  # PANCREAS
  # ============================================================================

  "Pancreatic Alpha Cells" = c(
    GCG = 1.0, ARX = 0.9, IRX2 = 0.8, TTR = 0.7, GC = 0.6,
    PCSK2 = 0.6, DPP4 = 0.5, MAFB = 0.5, LOXL4 = 0.5, TM4SF4 = 0.4,
    PLCE1 = 0.4, FAP = 0.3
  ),

  "Pancreatic Beta Cells" = c(
    INS = 1.0, IAPP = 0.9, MAFA = 0.8, NKX6_1 = 0.8, PDX1 = 0.7,
    SLC30A8 = 0.7, HADH = 0.6, PCSK1 = 0.6, DLK1 = 0.5, NPTX2 = 0.5,
    UCN3 = 0.4, GCK = 0.4
  ),

  "Pancreatic Delta Cells" = c(
    SST = 1.0, HHEX = 0.8, LEPR = 0.7, PCSK1 = 0.6, RBP4 = 0.6,
    GABRB3 = 0.5, UNC5B = 0.5, SEC11C = 0.5, CD24 = 0.4, BCHE = 0.4,
    PRG4 = 0.3, GHSR = 0.3
  ),

  # ============================================================================
  # CARDIAC
  # ============================================================================

  "Cardiomyocytes" = c(
    TNNT2 = 1.0, MYH7 = 0.9, MYH6 = 0.8, ACTC1 = 0.8, MYL2 = 0.7,
    MYL7 = 0.7, ACTN2 = 0.6, TTN = 0.6, RYR2 = 0.5, MYBPC3 = 0.5,
    NPPA = 0.5, NKX2_5 = 0.4
  ),

  # ============================================================================
  # ADDITIONAL IMMUNE STATES & SUBTYPES
  # ============================================================================

  "Tissue-Resident Memory T-cells (Trm)" = c(
    CD69 = 1.0, ITGAE = 0.9, ZNF683 = 0.9, CXCR6 = 0.8, ITGA1 = 0.7,
    CD3D = 0.6, CD3E = 0.6, RUNX3 = 0.5, BHLHE40 = 0.5, PRDM1 = 0.4,
    NR4A1 = 0.4, RGS1 = 0.5
  ),

  "Central Memory T-cells (Tcm)" = c(
    CCR7 = 1.0, SELL = 0.9, CD27 = 0.8, IL7R = 0.8, CD3D = 0.7,
    TCF7 = 0.7, CD28 = 0.6, LEF1 = 0.6, LTB = 0.5, S100A4 = 0.5,
    TRAC = 0.4, GPR183 = 0.4
  ),

  "Effector Memory T-cells (Tem)" = c(
    GZMK = 1.0, CCL5 = 0.9, NKG7 = 0.8, GZMA = 0.8, CD3D = 0.7,
    CST7 = 0.7, CXCR4 = 0.6, PRF1 = 0.5, KLRG1 = 0.5, CD44 = 0.5,
    EOMES = 0.4, IFNG = 0.4
  ),

  "M0 Macrophages" = c(
    CD68 = 1.0, CD163 = 0.6, CSF1R = 0.8, MERTK = 0.7, ADGRE1 = 0.6,
    AIF1 = 0.7, FCGR1A = 0.6, LYZ = 0.5, TYROBP = 0.5, C1QA = 0.4,
    FCER1G = 0.4, SPI1 = 0.4
  ),

  "Disease-Associated Macrophages (DAM)" = c(
    TREM2 = 1.0, SPP1 = 0.9, GPNMB = 0.8, LGALS3 = 0.8, FABP5 = 0.7,
    CD9 = 0.7, CTSD = 0.6, LPL = 0.6, CST3 = 0.5, FTL = 0.5,
    CTSB = 0.4, APOE = 0.5
  ),

  "Tumor-Associated Macrophages (TAM)" = c(
    CD163 = 1.0, MRC1 = 0.9, MSR1 = 0.8, CD68 = 0.8, MARCO = 0.7,
    TREM2 = 0.7, SPP1 = 0.6, APOE = 0.6, C1QA = 0.5, C1QB = 0.5,
    VEGFA = 0.4, CCL18 = 0.5
  ),

  "Tolerogenic Dendritic Cells" = c(
    IDO1 = 1.0, PDCD1LG2 = 0.9, CD274 = 0.8, IL10 = 0.7, ENTPD1 = 0.7,
    ITGAX = 0.6, HLA_DRA = 0.5, CCR7 = 0.4, SOCS1 = 0.5, TGFB1 = 0.4,
    LGALS9 = 0.4, FOXO3 = 0.3
  ),

  "Myeloid-Derived Suppressor Cells (MDSC)" = c(
    S100A8 = 1.0, S100A9 = 0.9, S100A12 = 0.8, ARG1 = 0.8, NOS2 = 0.7,
    CYBB = 0.7, CD33 = 0.6, ITGAM = 0.6, IL4R = 0.5, CSF3R = 0.5,
    CEACAM8 = 0.4, MMP9 = 0.4
  ),

  # ============================================================================
  # FIBROBLAST SUBTYPES
  # ============================================================================

  "Inflammatory Fibroblasts" = c(
    IL6 = 1.0, CXCL1 = 0.9, CXCL2 = 0.8, CCL2 = 0.8, IL1B = 0.7,
    CXCL8 = 0.7, COL1A1 = 0.5, PDGFRA = 0.5, FAP = 0.5, HAS2 = 0.5,
    POSTN = 0.4, MMP3 = 0.5
  ),

  "Reticular Fibroblasts" = c(
    PI16 = 1.0, MFAP5 = 0.9, SFRP2 = 0.8, PDGFRA = 0.7, DPP4 = 0.7,
    COL14A1 = 0.6, DCN = 0.6, LUM = 0.5, ELN = 0.5, CFD = 0.5,
    APOD = 0.4, GSN = 0.4
  ),

  "Adventitial Fibroblasts" = c(
    MFAP5 = 1.0, FBLN1 = 0.9, PDGFRA = 0.8, PI16 = 0.7, SFRP4 = 0.7,
    ELN = 0.6, MGP = 0.6, IGFBP6 = 0.5, COL6A3 = 0.5, DCN = 0.5,
    ADH1B = 0.5, CXCL14 = 0.4
  ),

  # ============================================================================
  # ADIPOCYTE SUBTYPES
  # ============================================================================

  "White Adipocytes" = c(
    LEP = 1.0, ADIPOQ = 0.9, FABP4 = 0.8, PLIN1 = 0.8, LPL = 0.7,
    PPARG = 0.7, CEBPA = 0.6, CIDEC = 0.6, CFD = 0.5, RBP4 = 0.5,
    FASN = 0.4, SLC2A4 = 0.4
  ),

  "Brown / Beige Adipocytes" = c(
    UCP1 = 1.0, CIDEA = 0.9, PRDM16 = 0.8, DIO2 = 0.8, PPARGC1A = 0.7,
    ELOVL3 = 0.7, TBX1 = 0.6, ADIPOQ = 0.5, FABP4 = 0.5, EBF2 = 0.5,
    COX7A1 = 0.5, TMEM26 = 0.4
  ),

  # ============================================================================
  # VASCULAR SUBTYPES
  # ============================================================================

  "Arterial Endothelial Cells" = c(
    GJA5 = 1.0, BMX = 0.9, EFNB2 = 0.8, DLL4 = 0.7, SOX17 = 0.7,
    HEY1 = 0.6, PECAM1 = 0.6, CDH5 = 0.5, SEMA3G = 0.5, CXCL12 = 0.5,
    MGP = 0.4, NOTCH4 = 0.4
  ),

  "Venous Endothelial Cells" = c(
    ACKR1 = 1.0, NR2F2 = 0.9, EPHB4 = 0.8, PECAM1 = 0.7, CDH5 = 0.7,
    VWF = 0.6, PLVAP = 0.6, SELP = 0.5, CLU = 0.5, VCAM1 = 0.5,
    EMCN = 0.4, SELE = 0.4
  ),

  "Tip Cells (Angiogenic)" = c(
    ESM1 = 1.0, APLN = 0.9, ANGPT2 = 0.8, PGF = 0.7, DLL4 = 0.7,
    CXCR4 = 0.6, PECAM1 = 0.5, KDR = 0.6, CDH5 = 0.5, NRP1 = 0.5,
    PDGFB = 0.4, MCAM = 0.4
  ),

  # ============================================================================
  # BONE MARROW NICHE
  # ============================================================================

  "BM Mesenchymal Stromal Cells" = c(
    CXCL12 = 1.0, KITL = 0.9, LEPR = 0.8, PDGFRA = 0.7, PDGFRB = 0.7,
    ANGPT1 = 0.6, VCAM1 = 0.6, THY1 = 0.5, FOXC1 = 0.5, EBF3 = 0.4,
    ADIPOQ = 0.3, ADRB3 = 0.3
  ),

  "Osteoclasts" = c(
    ACP5 = 1.0, CTSK = 0.9, MMP9 = 0.8, RANK = 0.7, DCSTAMP = 0.7,
    CALCR = 0.6, OCSTAMP = 0.6, TRAP = 0.5, ITGB3 = 0.5, NFATC1 = 0.4,
    CA2 = 0.4, OSCAR = 0.5
  ),

  "Osteoblasts" = c(
    BGLAP = 1.0, SP7 = 0.9, RUNX2 = 0.8, ALPL = 0.8, SPP1 = 0.7,
    COL1A1 = 0.6, IBSP = 0.6, DMP1 = 0.5, PHEX = 0.5, SOST = 0.4,
    BMP2 = 0.4, PTH1R = 0.4
  ),

  # ============================================================================
  # ADDITIONAL CYCLING STATES
  # ============================================================================

  "Cycling Myeloid Cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, CD14 = 0.7, LYZ = 0.7, STMN1 = 0.7,
    BIRC5 = 0.7, UBE2C = 0.6, CENPF = 0.6, TYMS = 0.5, CST3 = 0.5,
    CDK1 = 0.5, PCNA = 0.4
  ),

  "Cycling Epithelial Cells" = c(
    MKI67 = 1.0, TOP2A = 0.9, EPCAM = 0.8, KRT18 = 0.7, STMN1 = 0.7,
    BIRC5 = 0.7, UBE2C = 0.6, CENPF = 0.6, TYMS = 0.5, KRT8 = 0.5,
    CDK1 = 0.5, PCNA = 0.4
  ),

  # ============================================================================
  # STRESS / DAMAGE RESPONSE STATES
  # ============================================================================

  "Stress-Response Cells" = c(
    FOS = 1.0, JUN = 0.9, JUNB = 0.8, ATF3 = 0.8, EGR1 = 0.7,
    FOSB = 0.7, HSPA1A = 0.6, HSPA1B = 0.6, DUSP1 = 0.5, IER2 = 0.5,
    NR4A1 = 0.5, ZFP36 = 0.4
  ),

  "Interferon-Stimulated Cells" = c(
    ISG15 = 1.0, IFIT1 = 0.9, IFIT3 = 0.8, MX1 = 0.8, IFI6 = 0.7,
    OAS1 = 0.7, IFI44L = 0.6, RSAD2 = 0.6, STAT1 = 0.5, IRF7 = 0.5,
    IFITM1 = 0.5, BST2 = 0.4
  ),

  "Apoptotic Cells" = c(
    CASP3 = 1.0, BAX = 0.9, BCL2L11 = 0.8, CYCS = 0.7, PMAIP1 = 0.7,
    BBC3 = 0.6, BID = 0.6, BNIP3 = 0.5, DIABLO = 0.5, APAF1 = 0.4,
    TP53 = 0.4, FAS = 0.4
  )
)

# Pathway databases
PATHWAY_DATABASES <- list(
  Hallmark = list(
    HALLMARK_HYPOXIA = c("PGK1", "VEGFA", "ENO1", "LDHA", "ALDOA", "GAPDH", "CDKN1A", "TPI1", "P4HA1", "BNIP3"),
    HALLMARK_GLYCOLYSIS = c("PKM", "LDHA", "ALDOA", "ENO1", "GAPDH", "PGK1", "PGAM1", "TPI1", "GPI", "PFKL"),
    HALLMARK_OXIDATIVE_PHOSPHORYLATION = c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6C", "COX7A2", "COX7C", "NDUFA4", "NDUFB2", "NDUFB3"),
    HALLMARK_G2M_CHECKPOINT = c("CCNB1", "CCNB2", "CDK1", "TOP2A", "CENPA", "CENPE", "CENPF", "PLK1", "AURKB", "BUB1"),
    HALLMARK_E2F_TARGETS = c("MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "PCNA", "TYMS", "RRM1", "RRM2"),
    HALLMARK_MYC_TARGETS_V1 = c("NPM1", "RPLP0", "RPS2", "RPS3", "RPS6", "RPL3", "RPL4", "RPL5", "RPL6", "RPL7"),
    HALLMARK_MTORC1_SIGNALING = c("RPS6", "RPL3", "RPL4", "RPL5", "RPL6", "RPL7", "RPL8", "RPL9", "RPL10", "RPL11"),
    HALLMARK_INTERFERON_ALPHA_RESPONSE = c("ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "OAS3", "IFITM1", "IFITM3"),
    HALLMARK_INTERFERON_GAMMA_RESPONSE = c("STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "GBP1", "GBP2", "TAP1", "PSMB8", "PSMB9"),
    HALLMARK_INFLAMMATORY_RESPONSE = c("IL1B", "IL6", "TNF", "CXCL1", "CXCL2", "CXCL8", "CCL2", "CCL5", "PTGS2", "NFKB1"),
    HALLMARK_TNFA_SIGNALING_VIA_NFKB = c("NFKB1", "NFKB2", "REL", "RELA", "RELB", "TNF", "IL1B", "IL6", "CXCL1", "CXCL2"),
    HALLMARK_APOPTOSIS = c("CASP3", "CASP7", "CASP8", "CASP9", "FAS", "FADD", "BID", "BAX", "BAK1", "BCL2"),
    HALLMARK_P53_PATHWAY = c("CDKN1A", "TP53", "MDM2", "GADD45A", "GADD45B", "BAX", "BBC3", "PMAIP1", "FAS", "TNFRSF10B"),
    HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION = c("VIM", "CDH2", "FN1", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "MMP2", "MMP9", "TWIST1"),
    HALLMARK_ANGIOGENESIS = c("VEGFA", "VEGFB", "VEGFC", "FLT1", "KDR", "FLT4", "ANGPT1", "ANGPT2", "TEK", "TIE1")
  ),

  KEGG = list(
    Cell_Cycle = c("CCNA2", "CCNB1", "CCNB2", "CCND1", "CCNE1", "CDK1", "CDK2", "CDK4", "CDK6", "E2F1"),
    Apoptosis = c("CASP3", "CASP7", "CASP8", "CASP9", "CYCS", "APAF1", "BAX", "BCL2", "BID", "FAS"),
    JAK_STAT_Signaling = c("JAK1", "JAK2", "JAK3", "TYK2", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B"),
    MAPK_Signaling = c("MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14", "RAF1"),
    WNT_Signaling = c("WNT1", "WNT2", "WNT3", "WNT3A", "WNT5A", "CTNNB1", "APC", "AXIN1", "AXIN2", "GSK3B"),
    Notch_Signaling = c("NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", "HES1"),
    Hedgehog_Signaling = c("SHH", "IHH", "DHH", "PTCH1", "PTCH2", "SMO", "GLI1", "GLI2", "GLI3", "SUFU"),
    TGF_Beta_Signaling = c("TGFB1", "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "SMAD2", "SMAD3", "SMAD4", "SMAD7", "BMPR1A")
  )
)

# Ligand-Receptor pairs database
LIGAND_RECEPTOR_PAIRS <- list(
  # Immune checkpoints
  PD1_PDL1 = list(
    c("PDCD1", "CD274"),
    c("PDCD1", "PDCD1LG2")
  ),

  CTLA4 = list(
    c("CTLA4", "CD80"),
    c("CTLA4", "CD86")
  ),

  TIM3 = list(
    c("HAVCR2", "LGALS9"),
    c("HAVCR2", "CEACAM1")
  ),

  LAG3 = list(
    c("LAG3", "HLA-A"),
    c("LAG3", "HLA-B"),
    c("LAG3", "HLA-C")
  ),

  TIGIT = list(
    c("TIGIT", "PVR"),
    c("TIGIT", "PVRL2")
  ),

  # Cytokines and chemokines
  IL2 = list(
    c("IL2", "IL2RA"),
    c("IL2", "IL2RB"),
    c("IL2", "IL2RG")
  ),

  IL6 = list(
    c("IL6", "IL6R"),
    c("IL6", "IL6ST")
  ),

  IL10 = list(
    c("IL10", "IL10RA"),
    c("IL10", "IL10RB")
  ),

  TNF = list(
    c("TNF", "TNFRSF1A"),
    c("TNF", "TNFRSF1B")
  ),

  IFNG = list(
    c("IFNG", "IFNGR1"),
    c("IFNG", "IFNGR2")
  ),

  TGFB = list(
    c("TGFB1", "TGFBR1"),
    c("TGFB1", "TGFBR2"),
    c("TGFB2", "TGFBR1"),
    c("TGFB2", "TGFBR2")
  ),

  # Growth factors
  VEGF = list(
    c("VEGFA", "FLT1"),
    c("VEGFA", "KDR"),
    c("VEGFB", "FLT1"),
    c("VEGFC", "FLT4")
  ),

  EGF = list(
    c("EGF", "EGFR"),
    c("AREG", "EGFR"),
    c("EREG", "EGFR"),
    c("BTC", "EGFR")
  ),

  FGF = list(
    c("FGF1", "FGFR1"),
    c("FGF2", "FGFR1"),
    c("FGF2", "FGFR2")
  ),

  PDGF = list(
    c("PDGFA", "PDGFRA"),
    c("PDGFB", "PDGFRB")
  ),

  HGF = list(
    c("HGF", "MET")
  ),

  # Cell adhesion
  Selectins = list(
    c("SELE", "SELPLG"),
    c("SELP", "SELPLG"),
    c("SELL", "CD34")
  ),

  Integrins = list(
    c("ITGAL", "ICAM1"),
    c("ITGAM", "ICAM1"),
    c("ITGA4", "VCAM1")
  ),

  # Notch signaling
  Notch = list(
    c("DLL1", "NOTCH1"),
    c("DLL1", "NOTCH2"),
    c("DLL4", "NOTCH1"),
    c("JAG1", "NOTCH1"),
    c("JAG1", "NOTCH2")
  ),

  # WNT signaling
  WNT = list(
    c("WNT3A", "FZD1"),
    c("WNT3A", "FZD7"),
    c("WNT5A", "FZD5"),
    c("WNT5A", "ROR2")
  )
)

# Function to get all markers for a cell type
get_cell_type_markers <- function(cell_type) {
  if (cell_type %in% names(CELL_TYPE_MARKERS)) {
    return(CELL_TYPE_MARKERS[[cell_type]])
  } else {
    return(list())
  }
}

# Function to get all cell types
get_all_cell_types <- function() {
  return(names(CELL_TYPE_MARKERS))
}

# Function to search for markers across all cell types
search_markers <- function(gene_name) {
  results <- list()
  for (cell_type in names(CELL_TYPE_MARKERS)) {
    markers <- CELL_TYPE_MARKERS[[cell_type]]
    if (gene_name %in% names(markers)) {
      results[[cell_type]] <- markers[[gene_name]]
    }
  }
  return(results)
}

# =============================================================================
# Live Marker Database from PanglaoDB
# =============================================================================

#' Fetch marker genes from PanglaoDB (online TSV)
#' @param species Character: "Hs" (human), "Mm" (mouse), or "Mm Hs" (both)
#' @param tissue Character: optional tissue filter (e.g., "Brain", "Blood")
#' @param min_sensitivity Numeric: minimum sensitivity threshold (0-1)
#' @param min_specificity Numeric: minimum specificity threshold (0-1)
#' @return Named list of named numeric vectors (cell_type -> gene_name -> sensitivity)
#'         or NULL on failure
fetch_panglao_markers <- function(species = "Hs", tissue = NULL,
                                   min_sensitivity = 0.25,
                                   min_specificity = 0.50) {
  url <- "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz"
  df <- tryCatch({
    con <- gzcon(url(url, "rb"))
    on.exit(close(con), add = TRUE)
    read.table(con, sep = "\t", header = TRUE, stringsAsFactors = FALSE,
               quote = "")
  }, error = function(e) {
    message("PanglaoDB fetch failed: ", e$message, " - using built-in DB")
    return(NULL)
  })
  if (is.null(df)) return(NULL)

  # Filter by species
  df <- df[df$species %in% c(species, "Mm Hs"), ]

  # Filter by tissue if specified
  if (!is.null(tissue) && nzchar(tissue))
    df <- df[tolower(df$organ) == tolower(tissue), ]

  # Filter by sensitivity and specificity
  df <- df[!is.na(df$sensitivity) & df$sensitivity >= min_sensitivity, ]
  df <- df[!is.na(df$specificity) & df$specificity >= min_specificity, ]

  if (nrow(df) == 0) return(NULL)

  # Build marker database: list of named numeric vectors
  split(
    setNames(as.numeric(df$sensitivity), df$official.gene.symbol),
    df$cell.type
  )
}
