## =====================================================================
## Tumor-myeloid compartmentalization of polyamine metabolism in Ewing sarcoma
## Analysis code for: "Single-cell analysis of two Ewing sarcoma cohorts reveals
## tumor-myeloid compartmentalization of polyamine metabolism"
##
## Reproduces all results, figures and supplementary tables of the manuscript.
## Run from the directory that contains ES_project/ (see README.md for inputs).
## Time-consuming steps are cached in ES_project/02_results; delete a cached
## .rds file to recompute that step.
##
##   Part 1-2  Data and annotation              Methods 2.1-2.2, Fig. 1
##   Part 3    Patient-level module screen      Methods 2.3, Fig. 2
##   Part 4    Gene-level analysis              Methods 2.4, Fig. 3
##   Part 5    Reactome analysis                Methods 2.4, Fig. 4
##   Part 6    Exploratory bulk survival        Methods 2.6, Supp. Figs. 1-2
##   Part 7    Sensitivity analyses             Methods 2.5, Supp. Table 2
##   Part 8    Supplementary tables
##   Part 9    Exploratory ENPP1 analyses       Methods 2.7, Supp. Fig. 3
## =====================================================================

## =====================================================================
## Part 0. Setup: packages, paths, helper functions, gene sets
## =====================================================================

## ---- 0.1 Packages (installed if missing) ----
pkgs_cran <- c("Seurat", "harmony", "Matrix", "data.table", "ggplot2", "ggrepel",
               "patchwork", "survival", "openxlsx")
pkgs_bioc <- c("AnnotationDbi", "org.Hs.eg.db", "reactome.db",
               "hgu133plus2.db", "huex10sttranscriptcluster.db")
for (p in pkgs_cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in pkgs_bioc) if (!requireNamespace(p, quietly = TRUE))
  BiocManager::install(p, ask = FALSE, update = FALSE)
suppressPackageStartupMessages(
  for (p in c(pkgs_cran, pkgs_bioc)) library(p, character.only = TRUE))

## ---- 0.2 Paths (run from the directory that contains ES_project/) ----
proj     <- "ES_project"
if (!dir.exists(proj)) stop("ES_project/ not found in the working directory; set it with setwd()")
dir_data <- file.path(proj, "01_data")
dir_bulk <- file.path(dir_data, "bulk")
dir_res  <- file.path(proj, "02_results")
dir_fig  <- file.path(proj, "03_figures")
for (d in c(dir_data, dir_bulk, dir_res, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
cat("Working directory:", getwd(), "\n")

## ---- 0.3 Helper functions ----
## Load a cached result if present; otherwise compute and save it
cached <- function(file, expr, dir = dir_res) {
  f <- file.path(dir, file)
  if (file.exists(f)) { cat("Loading cached result:", file, "\n"); return(readRDS(f)) }
  x <- expr; saveRDS(x, f); x
}
## Save each figure as PDF, 300-dpi PNG and 300-dpi TIFF
save_fig <- function(p, name, w, h) {
  ggsave(file.path(dir_fig, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(dir_fig, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  ggsave(file.path(dir_fig, paste0(name, ".tiff")), p, width = w, height = h, dpi = 300,
         device = "tiff", compression = "lzw")                     # submission format (TIFF)
  cat("Saved figure:", name, "\n")
}
## z-score
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

col_t <- "#B8412F"   # tumor / tumor-high: red
col_m <- "#2F6DB5"   # myeloid / myeloid-high: blue

## ---- 0.4 Fourteen curated metabolic gene sets (Supplementary Table 1; current and former HGNC symbols, absent symbols ignored) ----
sets_full <- list(
 OXPHOS = c("NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA8","NDUFA9","NDUFA10","NDUFA11","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10","NDUFB11","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","SDHA","SDHB","SDHC","SDHD","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","UQCR10","CYC1","COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7B","COX7C","COX8A","ATP5A1","ATP5B","ATP5C1","ATP5D","ATP5E","ATP5F1","ATP5G1","ATP5G3","ATP5H","ATP5I","ATP5J","ATP5J2","ATP5L","ATP5O","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5PB","ATP5MC1","ATP5MC3","ATP5PD","ATP5ME","ATP5PF","ATP5MG","ATP5PO"),
 Glycolysis = c("SLC2A1","HK1","HK2","GPI","PFKL","PFKM","PFKP","ALDOA","ALDOC","TPI1","GAPDH","PGK1","PGAM1","ENO1","ENO2","PKM","LDHA","PDK1","SLC16A3"),
 TCA_cycle = c("CS","ACO2","IDH2","IDH3A","IDH3B","IDH3G","OGDH","DLST","DLD","SUCLG1","SUCLA2","SDHA","FH","MDH2","PDHA1","PDHB","DLAT"),
 FA_oxidation = c("CPT1A","CPT2","ACADM","ACADVL","ACADS","HADHA","HADHB","HADH","ECHS1","ACAA2","ACSL1","ACSL3","ACSL4","SLC25A20","ETFA","ETFB","PPARA","PPARGC1A"),
 FA_synthesis_uptake = c("FASN","ACACA","ACLY","SCD","ELOVL6","ELOVL1","SREBF1","FADS1","FADS2","CD36","FABP4","FABP5","LPL"),
 Cholesterol_synth = c("HMGCR","HMGCS1","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","DHCR7","DHCR24","SC5D","NSDHL","MSMO1","LDLR","INSIG1","SREBF2"),
 Glutamine_AA = c("GLS","GLS2","GLUD1","GOT1","GOT2","GPT2","SLC1A5","SLC7A5","SLC3A2","ASNS","PYCR1","ALDH18A1"),
 Serine_1C = c("PHGDH","PSAT1","PSPH","SHMT1","SHMT2","MTHFD1","MTHFD2","MTHFD1L","MTHFD2L","TYMS","DHFR","SLC1A4"),
 Redox_FerroDefense = c("GPX4","SLC7A11","GCLC","GCLM","GSS","GSR","TXNRD1","TXN","NFE2L2","SOD1","SOD2","PRDX1","PRDX3","NQO1","HMOX1","FTH1","FTL","AIFM2"),
 Nucleotide_synth = c("PPAT","GART","PFAS","PAICS","ATIC","IMPDH1","IMPDH2","GMPS","CAD","DHODH","UMPS","CTPS1","RRM1","RRM2","TYMS","ADSL","ADSS"),
 Pentose_phosphate = c("G6PD","PGD","TKT","TALDO1","RPIA","RPE","PGLS"),
 Polyamine = c("ODC1","SRM","SMS","AMD1","SAT1","PAOX","SMOX","AZIN1"),
 Hypoxia_HIF = c("HIF1A","EPAS1","VEGFA","CA9","BNIP3","BNIP3L","ENO1","PGK1","LDHA","SLC2A1","ADM","NDRG1","ERO1A","P4HA1","P4HA2"),
 Purinergic = c("ENPP1","NT5E","ENTPD1","ADORA2A","ADORA2B","CD38","NAMPT")
)
saveRDS(sets_full, file.path(dir_data, "metabolic_sets_full.rds"))
sc_mod <- names(sets_full)

## ---- 0.5 Map cell types to eight compartments ----
## Cohort A: cluster annotations from Part 1; glial and mast cells are not assigned (excluded)
broad1 <- c(Tumor_core = "Tumor", Tumor_ES030 = "Tumor", Tumor_ES030b = "Tumor",
            Tumor_ES026_016 = "Tumor", Tumor_neural = "Tumor", Tumor_cycling = "Tumor_cycling",
            Macrophage = "Myeloid", Monocyte = "Myeloid", DC = "Myeloid",
            T_cell = "T_NK", NK = "T_NK", B_plasma = "B", Endothelial = "Endothelial",
            Fibroblast = "CAF", Pericyte = "Pericyte")
## Cohort B: author-provided Cell.Type labels (no pericyte population)
broad2 <- c("Ewing sarcoma" = "Tumor", "Ewing sarcoma proliferating" = "Tumor_cycling",
            Myeloid = "Myeloid", T = "T_NK", B = "B", Endothelial = "Endothelial",
            "Matrix CAF" = "CAF", "Inflammatory CAF" = "CAF")
save(broad1, broad2, file = file.path(dir_data, "celltype_maps.RData"))

## Compute the 14 module scores per cell (columns named after the gene sets)
add_scores <- function(o, sets) {
  s <- lapply(sets, function(x) intersect(unique(x), rownames(o)))
  cat("Genes present per module:\n"); print(sapply(s, length))
  o <- AddModuleScore(o, features = s, name = "MS")
  idx <- match(paste0("MS", seq_along(s)), colnames(o@meta.data))
  colnames(o@meta.data)[idx] <- names(s)
  o
}


## =====================================================================
## Part 1. Cohort A: GSE243347 (CEL-seq2, 384-well plates)
## =====================================================================
fA <- file.path(dir_data, "obj_A_GSE243347.rds")
if (file.exists(fA)) {
  cat("Loading saved cohort A object\n")
  obj <- readRDS(fA)
} else {
  cat("Building cohort A from raw files (about 15-20 min)\n")
  cat("Select GSE243347_RAW.tar in the file dialog\n")
  tarfile <- file.choose()
  exdir <- file.path(dirname(tarfile), "raw")
  if (!dir.exists(exdir) || length(list.files(exdir)) == 0) untar(tarfile, exdir = exdir)
  files <- list.files(exdir, pattern = "GSM", full.names = TRUE)
  gsm <- as.numeric(sub("GSM(\\d+)_.*", "\\1", basename(files)))
  ## Plate (GSM) -> patient; MIXED and ES-048 plates are excluded
  pmap <- c(
   "7783783"="ES-024","7783784"="ES-024","7783785"="ES-026","7783786"="ES-026",
   "7783787"="ES-025","7783788"="ES-025","7783789"="ES-016","7783790"="ES-016",
   "7783791"="ES-027","7783792"="ES-027","7783793"="ES-030","7783794"="ES-030",
   "7783795"="ES-030","7783796"="ES-006","7783797"="ES-034","7783798"="ES-036",
   "7783799"="ES-006","7783800"="ES-039","7783801"="ES-042","7783802"="ES-042",
   "7783803"="ES-042","7783804"="ES-042","7783805"="ES-046","7783806"="ES-162",
   "7783807"="MIXED","7783808"="ES-048","7783809"="ES-049")
  lst <- list()
  for (i in seq_along(files)) {
    g <- as.character(gsm[i]); pt <- pmap[g]
    if (pt %in% c("MIXED", "ES-048")) next
    m <- read.table(files[i], header = TRUE, row.names = 1, sep = "\t",
                    check.names = FALSE, quote = "")
    m <- m[, colnames(m) != "UNK", drop = FALSE]          # drop unassigned barcodes
    sym <- sub(".*__", "", rownames(m))
    m <- rowsum(as.matrix(m), sym)
    m <- m[!grepl("^ERCC-", rownames(m)), , drop = FALSE] # drop ERCC spike-ins
    ok <- colSums(m) >= 1000 & colSums(m > 0) >= 500       # initial cell filter
    m <- as(m[, ok, drop = FALSE], "dgCMatrix")
    colnames(m) <- paste0("GSM", g, "_", colnames(m))
    lst[[g]] <- m
    cat(g, pt, ncol(m), "\n")
  }
  ## Each plate lists only genes detected on it, so the union of genes is used (an intersection would drop e.g. PTPRC)
  allg <- Reduce(union, lapply(lst, rownames))
  expand <- function(x) {
    s <- summary(x); idx <- match(rownames(x), allg)
    sparseMatrix(i = idx[s$i], j = s$j, x = s$x, dims = c(length(allg), ncol(x)),
                 dimnames = list(allg, colnames(x)))
  }
  counts <- do.call(cbind, lapply(lst, expand))
  meta <- data.frame(row.names = colnames(counts), gsm = sub("_.*", "", colnames(counts)))
  meta$patient <- unname(pmap[sub("GSM", "", meta$gsm)])
  obj <- CreateSeuratObject(counts, meta.data = meta, min.cells = 3)
  obj[["pct_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = pct_mt < 25)
  obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 3000) |>
         ScaleData() |> RunPCA(npcs = 30)
  obj <- RunHarmony(obj, group.by.vars = "patient")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20) |>
         FindNeighbors(reduction = "harmony", dims = 1:20) |>
         FindClusters(resolution = 0.6)
  cat("\nCluster x patient table (cluster 1 should be mostly ES-030; cluster 3 ES-026/ES-016):\n")
  print(table(obj$seurat_clusters, obj$patient))
  ## Annotation from canonical markers (PTPRC, CD3E, CD68, COL1A1, PECAM1, NKX2-2, STEAP1, etc.)
  ann <- c("0"="Tumor_core","1"="Tumor_ES030","2"="Endothelial","3"="Tumor_ES026_016",
           "4"="T_cell","5"="Macrophage","6"="Pericyte","7"="Tumor_cycling",
           "8"="Tumor_ES030b","9"="Fibroblast","10"="NK","11"="Tumor_neural",
           "12"="DC","13"="Monocyte","14"="B_plasma","15"="Glial","16"="Muscle","17"="Mast")
  obj$celltype <- unname(ann[as.character(obj$seurat_clusters)])
  obj <- subset(obj, subset = celltype != "Muscle")        # remove skeletal-muscle contamination
}
if (!all(sc_mod %in% colnames(obj@meta.data))) obj <- add_scores(obj, sets_full)
obj$broad <- unname(broad1[as.character(obj$celltype)])
saveRDS(obj, fA)


## =====================================================================
## Part 2. Cohort B: GSE277083 (10x Genomics, author cell-type labels)
## =====================================================================
fB <- file.path(dir_data, "obj_B_GSE277083.rds")
if (file.exists(fB)) {
  cat("Loading saved cohort B object\n")
  obj2 <- readRDS(fB)
} else {
  cat("Building cohort B: select the cell metadata file, then the counts matrix file (primary tumors)\n")
  f_meta <- file.choose()
  f_cnt  <- file.choose()
  meta2 <- as.data.frame(fread(f_meta)); rownames(meta2) <- meta2$Cell
  cnt <- fread(f_cnt)
  genes <- make.unique(cnt[[1]])
  mm0 <- as.matrix(cnt[, -1]); rownames(mm0) <- genes
  mm0 <- as(mm0, "dgCMatrix"); rm(cnt); gc()
  cat("Check: total UMI of first cell =", sum(mm0[, 1]), "; metadata nCount_RNA =", meta2$nCount_RNA[1], "\n")
  obj2 <- CreateSeuratObject(mm0, meta.data = meta2[colnames(mm0), ])
  obj2 <- NormalizeData(obj2)
}
if (!all(sc_mod %in% colnames(obj2@meta.data))) obj2 <- add_scores(obj2, sets_full)
obj2$broad <- unname(broad2[obj2$Cell.Type])
## Harmony + UMAP for visualization only (Fig. 1); annotation uses the author labels
if (!"umap" %in% Reductions(obj2)) {
  obj2 <- FindVariableFeatures(obj2, nfeatures = 3000) |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)
  obj2 <- RunHarmony(obj2, group.by.vars = "Sample2")
  obj2 <- RunUMAP(obj2, reduction = "harmony", dims = 1:20)
}
saveRDS(obj2, fB)

cat("\nCohort A:", ncol(obj), "cells,", length(unique(obj$patient)), "patients\n")
cat("Cohort B:", ncol(obj2), "cells,", length(unique(obj2$Sample2)), "patients\n")
cat("Cells per compartment:\n"); print(table(obj$broad, useNA = "ifany")); print(table(obj2$broad, useNA = "ifany"))

## ---- Fig. 1: UMAP of both cohorts ----
lv <- c("Tumor", "Tumor_cycling", "Myeloid", "T_NK", "B", "CAF", "Endothelial", "Pericyte")
obj$comp  <- factor(obj$broad,  levels = lv)
obj2$comp <- factor(obj2$broad, levels = lv)
pa <- DimPlot(obj,  group.by = "comp", label = TRUE, repel = TRUE) + ggtitle("Cohort A (GSE243347): compartment") + NoLegend()
pb <- DimPlot(obj2, group.by = "comp", label = TRUE, repel = TRUE) + ggtitle("Cohort B (GSE277083): compartment") + NoLegend()
pc <- DimPlot(obj,  group.by = "patient") + ggtitle("Cohort A: patient")
pd <- DimPlot(obj2, group.by = "Sample2") + ggtitle("Cohort B: patient")
save_fig((pa | pb) / (pc | pd), "Fig1_umap_both_cohorts", 14, 11)


## =====================================================================
## Part 3. Patient-level module screen: myeloid minus non-cycling tumor (Fig. 2)
## =====================================================================
## Mean score per patient and compartment (at least 10 cells), then myeloid minus tumor
pat_diff <- function(md, pat, a = "Myeloid", b = "Tumor", minn = 10) {
  d <- data.frame(pat = md[[pat]], broad = md$broad, md[, sc_mod, drop = FALSE])
  d <- d[!is.na(d$broad), ]
  n <- aggregate(list(n = d$pat), list(pat = d$pat, broad = d$broad), length)
  m <- merge(aggregate(d[, sc_mod], list(pat = d$pat, broad = d$broad), mean), n)
  m <- m[m$n >= minn, ]
  pats <- intersect(m$pat[m$broad == a], m$pat[m$broad == b])
  X <- m[m$broad == a, ]; X <- X[match(pats, X$pat), ]
  Y <- m[m$broad == b, ]; Y <- Y[match(pats, Y$pat), ]
  out <- as.matrix(X[, sc_mod] - Y[, sc_mod]); rownames(out) <- pats
  sdc <- sapply(sc_mod, function(s) sd(d[[s]]))   # cell-level SD, used only for plotting
  list(raw = out, z = sweep(out, 2, sdc, "/"))
}
pdA <- pat_diff(obj@meta.data,  "patient")
pdB <- pat_diff(obj2@meta.data, "Sample2")
cat("\nEvaluable patients: A =", nrow(pdA$raw), "(", paste(rownames(pdA$raw), collapse = ", "), "); B =",
    nrow(pdB$raw), "(", paste(rownames(pdB$raw), collapse = ", "), ")\n")

summ_pat <- function(r) data.frame(module = sc_mod, n_patients = nrow(r$raw),
                                   mean_diff = colMeans(r$raw),
                                   n_myeloid_higher = colSums(r$raw > 0), row.names = NULL)
resC <- list(A = summ_pat(pdA), B = summ_pat(pdB))
cat("\nPatient-level results (AddModuleScore):\n"); print(resC, digits = 3)
saveRDS(list(resC = resC, rA = pdA, rB = pdB), file.path(dir_res, "patient_level_results.rds"))

## ---- Fig. 2: patient x module heatmap (blue = higher in myeloid, red = higher in tumor) ----
ord <- names(sort(colMeans(rbind(pdA$z, pdB$z)), decreasing = TRUE))
lab <- sapply(ord, function(s) sprintf("%s  [A %d/%d, B %d/%d]", gsub("_", " ", s),
              sum(pdA$raw[, s] > 0), nrow(pdA$raw), sum(pdB$raw[, s] > 0), nrow(pdB$raw)))
to_long <- function(r, coh) { x <- as.data.frame(as.table(r$z)); colnames(x) <- c("patient", "pathway", "z"); x$cohort <- coh; x }
hd <- rbind(to_long(pdA, "Cohort A (GSE243347)"), to_long(pdB, "Cohort B (GSE277083)"))
hd$pathway <- factor(lab[as.character(hd$pathway)], levels = rev(lab))
hd$z <- pmax(pmin(hd$z, 1.5), -1.5)
p_fig2 <- ggplot(hd, aes(patient, pathway, fill = z)) + geom_tile(colour = "white") +
  facet_wrap(~ cohort, scales = "free_x") +
  scale_fill_gradient2(low = col_t, mid = "white", high = col_m, midpoint = 0,
                       name = "Myeloid \u2212 tumor\n(SD units)") +
  labs(x = NULL, y = NULL) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_fig2, "Fig2_patient_level_heatmap", 11, 4.8)


## =====================================================================
## Part 4. Gene-level analysis: 40 genes of the four modules, tumor vs myeloid (Fig. 3)
## =====================================================================
## The four modules higher in myeloid cells across patients in Part 3
sets4  <- sets_full[c("Redox_FerroDefense", "Pentose_phosphate", "Polyamine", "Purinergic")]
genes4 <- unique(unlist(sets4))
cat("\nCandidate genes:", length(genes4), "\n")

## Tumor (incl. cycling) vs myeloid; all cells, no pre-filter; p_val_adj is Seurat's Bonferroni correction
de_genes4 <- function(o) {
  cells <- colnames(o)[which(o$broad %in% c("Tumor", "Tumor_cycling", "Myeloid"))]
  o <- subset(o, cells = cells)
  Idents(o) <- factor(ifelse(o$broad == "Myeloid", "Myeloid", "Tumor"))
  FindMarkers(o, ident.1 = "Tumor", ident.2 = "Myeloid", features = intersect(genes4, rownames(o)),
              logfc.threshold = 0, min.pct = 0)
}
DE40 <- cached("DE_genes40_tumor_vs_myeloid.rds", list(A = de_genes4(obj), B = de_genes4(obj2)))
keep <- c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")
tt <- merge(data.frame(gene = rownames(DE40$A), DE40$A[, keep]),
            data.frame(gene = rownames(DE40$B), DE40$B[, keep]), by = "gene", suffixes = c("_A", "_B"))
tt$program <- sapply(tt$gene, function(g) paste(names(sets4)[sapply(sets4, function(s) g %in% s)], collapse = ";"))
ok_t <- with(tt, avg_log2FC_A >  0.5 & avg_log2FC_B >  0.5 & p_val_adj_A < 0.05 & p_val_adj_B < 0.05)
ok_m <- with(tt, avg_log2FC_A < -0.5 & avg_log2FC_B < -0.5 & p_val_adj_A < 0.05 & p_val_adj_B < 0.05)
tt$arm <- ifelse(ok_t, "Tumor arm", ifelse(ok_m, "Myeloid arm", "Other"))
write.csv(tt, file.path(dir_res, "tumor_vs_myeloid_program_genes.csv"), row.names = FALSE)

tumor_g   <- c("ODC1", "SRM", "SLC7A11", "AIFM2", "ENPP1")
myeloid_g <- c("HMOX1", "FTL", "FTH1", "ENTPD1", "NAMPT", "SAT1")
cat("Replicated tumor-high genes:", paste(tt$gene[ok_t], collapse = ", "), "\n")
cat("Replicated myeloid-high genes:", paste(tt$gene[ok_m], collapse = ", "), "\n")
if (!setequal(tt$gene[ok_t], tumor_g) || !setequal(tt$gene[ok_m], myeloid_g))
  warning("Replicated genes differ from the 11 reported in the paper")

## ---- Fig. 3A/B: dot plots of the 11 genes across compartments ----
dot <- function(o, ttl) {
  o <- subset(o, cells = colnames(o)[!is.na(o$broad)])
  o$broad <- factor(o$broad, levels = lv)
  DotPlot(o, features = c(tumor_g, myeloid_g), group.by = "broad", cols = c("grey90", "#333333")) +
    ggtitle(ttl) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())
}
save_fig(dot(obj,  "Cohort A (GSE243347)"), "Fig3A_dotplot_A", 8, 4.2)
save_fig(dot(obj2, "Cohort B (GSE277083)"), "Fig3B_dotplot_B", 8, 4.2)

## ---- Fig. 3C: log2FC in cohort A vs cohort B ----
p_fig3c <- ggplot(tt, aes(avg_log2FC_A, avg_log2FC_B, colour = arm)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
  geom_point(data = subset(tt, arm == "Other"), size = 1.4, alpha = .6) +
  geom_point(data = subset(tt, arm != "Other"), size = 2.6) +
  geom_text_repel(data = subset(tt, arm != "Other"), aes(label = gene), size = 3.2,
                  max.overlaps = Inf, box.padding = .5, show.legend = FALSE) +
  scale_colour_manual(values = c("Tumor arm" = col_t, "Myeloid arm" = col_m, "Other" = "grey70")) +
  scale_x_continuous(expand = expansion(mult = .12)) + scale_y_continuous(expand = expansion(mult = .12)) +
  labs(x = "log2FC tumor vs myeloid, cohort A", y = "log2FC tumor vs myeloid, cohort B", colour = NULL) +
  theme_bw()
save_fig(p_fig3c, "Fig3C_logFC_replication", 5.8, 5.2)


## =====================================================================
## Part 5. Genome-wide differential expression and Reactome over-representation (Fig. 4)
## =====================================================================
## ---- 5.1 Reactome human pathways -> Entrez IDs (as character) ----
pw <- stack(as.list(reactomeEXTID2PATHID)); colnames(pw) <- c("pid", "entrez")
pw$pid <- as.character(pw$pid); pw$entrez <- as.character(pw$entrez)
pw <- pw[grepl("^R-HSA", pw$pid), ]
pw_nm <- unlist(as.list(reactomePATHID2NAME)); pw_nm <- sub("^Homo sapiens: ", "", pw_nm)
pw_gene <- lapply(split(pw$entrez, pw$pid), unique)
pw_gene <- pw_gene[lengths(pw_gene) >= 10 & lengths(pw_gene) <= 300]
cat("\nReactome pathways with 10-300 genes:", length(pw_gene), "\n")

## ---- 5.2 Restrict to metabolic pathways by name keywords ----
kw <- paste(c("metabol","biosynthesis","catabol","degradation of","synthesis","glycolysis","gluconeogenesis","citric acid","TCA",
              "respiratory electron","oxidative phosphorylation","pentose phosphate","glutathione","glutamate","glutamine","serine",
              "amino acid","polyamine","arginine","purine","pyrimidine","nucleotide","iron","heme","haem","ferritin","transferrin",
              "fatty acid","lipid","cholesterol","sphingolipid","NAD","nicotinamide","vitamin","SLC-mediated","solute carrier",
              "Nucleobase","adenosine","ATP","hypoxia","HIF"), collapse = "|")
metab_ids <- names(pw_nm)[grepl(kw, pw_nm, ignore.case = TRUE) & grepl("^R-HSA", names(pw_nm))]
metab_ids <- intersect(metab_ids, names(pw_gene))
cat("Metabolic pathways retained:", length(metab_ids), "\n")

## ---- 5.3 Genome-wide DE: tumor (incl. cycling) vs myeloid ----
de_all <- function(o) {
  cells <- colnames(o)[which(o$broad %in% c("Tumor", "Tumor_cycling", "Myeloid"))]
  o <- subset(o, cells = cells)
  Idents(o) <- factor(ifelse(o$broad == "Myeloid", "Myeloid", "Tumor"))
  FindMarkers(o, ident.1 = "Tumor", ident.2 = "Myeloid", logfc.threshold = 0.5, min.pct = 0.1,
              max.cells.per.ident = 500, verbose = FALSE)
}
DEall <- cached("DE_tumor_vs_myeloid.rds", list(dA = de_all(obj), dB = de_all(obj2)))
dA <- DEall$dA; dB <- DEall$dB
de_set <- function(d, up) rownames(d)[d$p_val_adj < 0.05 & (if (up) d$avg_log2FC > 0.5 else d$avg_log2FC < -0.5)]
tumA <- de_set(dA, TRUE); tumB <- de_set(dB, TRUE); myeA <- de_set(dA, FALSE); myeB <- de_set(dB, FALSE)
cat("Tumor-high genes A/B:", length(tumA), length(tumB), "; myeloid-high genes A/B:", length(myeA), length(myeB), "\n")

## ---- 5.4 Hypergeometric test (background = all genes in the dataset), BH correction ----
ora1 <- function(genes, universe_sym) {
  uni <- na.omit(suppressMessages(mapIds(org.Hs.eg.db, universe_sym, "ENTREZID", "SYMBOL")))
  ent <- na.omit(suppressMessages(mapIds(org.Hs.eg.db, genes, "ENTREZID", "SYMBOL")))
  r <- do.call(rbind, lapply(metab_ids, function(p) {
    pg <- intersect(pw_gene[[p]], uni); k <- length(intersect(ent, pg))
    if (k < 2) return(NULL)
    data.frame(pid = p, k = k, size = length(pg),
               p = phyper(k - 1, length(pg), length(uni) - length(pg), length(ent), lower.tail = FALSE))
  }))
  r$q <- p.adjust(r$p, "BH"); r$pathway <- pw_nm[r$pid]; r
}
## Run in both cohorts and keep the weaker (larger) q
side <- function(gA, gB) {
  rA <- ora1(gA, rownames(obj)); rB <- ora1(gB, rownames(obj2))
  m <- merge(rA[, c("pid", "pathway", "k", "size", "q")], rB[, c("pid", "k", "q")], by = "pid", suffixes = c("_A", "_B"))
  m$q_max <- pmax(m$q_A, m$q_B); m[order(m$q_max), ]
}
tum_up <- side(tumA, tumB)
mye_up <- side(myeA, myeB)

## ---- 5.5 Exclude spurious metabolic signals ----
## (a) pathways in which >50% of overlapping genes are ribosomal protein genes (RPL/RPS)
rp_frac <- function(pid, hits) {
  sym <- suppressMessages(mapIds(org.Hs.eg.db, pw_gene[[pid]], "SYMBOL", "ENTREZID"))
  h <- intersect(sym, hits); if (!length(h)) return(NA)
  mean(grepl("^RP[LS][0-9PLX]|^RPLP", h))
}
## (b) nucleic-acid replication, repair and chromatin pathways (proliferation) and non-metabolic immune signaling pathways
excl <- "chromatin|RMTs|Excision|AP sites|Nucleotide-binding|IL10|NLR|DNA|non-coding RNA|Lagging|Gap-filling|Telomere"
tum_up$rp_frac <- sapply(tum_up$pid, rp_frac, hits = intersect(tumA, tumB))
mye_up$rp_frac <- sapply(mye_up$pid, rp_frac, hits = intersect(myeA, myeB))
flag <- function(df) {
  df$kept <- !is.na(df$rp_frac) & df$rp_frac <= 0.5 & !grepl(excl, df$pathway) & df$q_max < 0.05
  df
}
tum_up <- flag(tum_up); mye_up <- flag(mye_up)
cat("\nPathways excluded as ribosomal-protein driven:\n")
print(subset(tum_up, q_max < 0.05 & rp_frac > 0.5)[, c("pathway", "q_max", "rp_frac")], row.names = FALSE, digits = 2)
tum_f <- subset(tum_up, kept); mye_f <- subset(mye_up, kept)
cat("\nRetained tumor-side pathways:\n"); print(tum_f[, c("pathway", "k_A", "k_B", "q_max")], row.names = FALSE, digits = 2)
cat("\nRetained myeloid-side pathways:\n"); print(mye_f[, c("pathway", "k_A", "k_B", "q_max")], row.names = FALSE, digits = 2)
write.csv(tum_up, file.path(dir_res, "reactome_metab_tumor_side.csv"), row.names = FALSE)
write.csv(mye_up, file.path(dir_res, "reactome_metab_myeloid_side.csv"), row.names = FALSE)

## ---- Fig 4 ----
bar_df <- rbind(cbind(head(tum_f[, c("pathway", "q_max")], 8), side = "Tumor-high"),
                cbind(head(mye_f[, c("pathway", "q_max")], 8), side = "Myeloid-high"))
wrap <- function(s) paste(strwrap(s, 45), collapse = "\n")
bar_df$lab <- sapply(bar_df$pathway, wrap)
bar_df$lab <- factor(bar_df$lab, levels = rev(unique(bar_df$lab)))
p_fig4 <- ggplot(bar_df, aes(-log10(q_max), lab, fill = side)) + geom_col() +
  facet_grid(side ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("Tumor-high" = col_t, "Myeloid-high" = col_m)) +
  labs(x = "-log10(q), weaker of the two cohorts", y = NULL) + theme_bw() + theme(legend.position = "none")
save_fig(p_fig4, "Fig4_reactome_metabolic", 7.5, 4)


## =====================================================================
## Part 6. Exploratory bulk survival analysis (Supplementary Figs. 1-2)
## =====================================================================
## Requires the four GEO series matrix files in 01_data/bulk

## ---- 6.0 Download the series matrix files if missing (otherwise download manually from GEO) ----
options(timeout = 600)
for (id in c("GSE17674", "GSE17618", "GSE63155", "GSE63156")) {
  dest <- file.path(dir_bulk, paste0(id, "_series_matrix.txt.gz"))
  if (!file.exists(dest)) {
    url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
                   sub("\\d{3}$", "nnn", id), id, id)
    tryCatch(download.file(url, dest, mode = "wb", quiet = TRUE), error = function(e) NULL)
  }
  cat(id, ifelse(file.exists(dest), "available", "missing: download manually from GEO"), "\n")
}

## ---- 6.1 Read series matrix files: expression and sample characteristics ----
read_sm <- function(f) {
  L   <- readLines(gzfile(f))
  geo <- gsub("\"", "", strsplit(L[grepl("^!Sample_geo_accession", L)][1], "\t")[[1]][-1])
  ch  <- L[grepl("^!Sample_characteristics_ch1", L)]
  cl  <- sapply(ch, function(x) { v <- gsub("\"", "", strsplit(x, "\t")[[1]][-1]); sub("^[^:]*: ?", "", v) })
  nm  <- sapply(ch, function(x) sub(":.*", "", gsub("\"", "", strsplit(x, "\t")[[1]][2])))
  cl  <- as.data.frame(cl, stringsAsFactors = FALSE); colnames(cl) <- make.names(nm, unique = TRUE)
  rownames(cl) <- geo
  ex <- read.delim(gzfile(f), comment.char = "!", row.names = 1, check.names = FALSE)
  ex <- as.matrix(ex); ex <- ex[, geo, drop = FALSE]
  if (max(ex, na.rm = TRUE) > 100) ex <- log2(ex + 1)
  list(expr = ex, clin = cl)
}
## Probe set -> gene symbol; for genes with several probe sets, keep the one with the highest mean
to_gene <- function(ex, db) {
  sym <- suppressMessages(mapIds(db, rownames(ex), "SYMBOL", "PROBEID"))
  ok  <- !is.na(sym); ex <- ex[ok, , drop = FALSE]; sym <- sym[ok]
  o   <- order(rowMeans(ex), decreasing = TRUE)
  keep_p <- !duplicated(sym[o]); ex <- ex[o, ][keep_p, ]; rownames(ex) <- sym[o][keep_p]
  ex
}
Bk <- cached("bulk_list.rds", dir = dir_bulk, expr = {
  b <- list(GSE17674 = read_sm(file.path(dir_bulk, "GSE17674_series_matrix.txt.gz")),
            GSE17618 = read_sm(file.path(dir_bulk, "GSE17618_series_matrix.txt.gz")),
            GSE63155 = read_sm(file.path(dir_bulk, "GSE63155_series_matrix.txt.gz")),
            GSE63156 = read_sm(file.path(dir_bulk, "GSE63156_series_matrix.txt.gz")))
  for (id in c("GSE17674", "GSE17618")) b[[id]]$gene <- to_gene(b[[id]]$expr, hgu133plus2.db)
  for (id in c("GSE63155", "GSE63156")) b[[id]]$gene <- to_gene(b[[id]]$expr, huex10sttranscriptcluster.db)
  b
})
## Document the exclusion of GSE17618
try({
ov <- intersect(Bk$GSE17674$clin$sample.id, Bk$GSE17618$clin$sample.id)
cat("\nGSE17618 samples shared with GSE17674:", length(ov), "; unique to GSE17618:",
    sum(!Bk$GSE17618$clin$sample.id %in% Bk$GSE17674$clin$sample.id), "(diagnosis:",
    paste(unique(Bk$GSE17618$clin$diagnosis[!Bk$GSE17618$clin$sample.id %in% Bk$GSE17674$clin$sample.id]), collapse = ", "),
    ") -> GSE17618 excluded\n")
}, silent = TRUE)

## ---- 6.2 Survival data for the three cohorts ----
num  <- function(x) suppressWarnings(as.numeric(x))
pick <- function(df, pat) {
  m <- grep(pat, colnames(df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(m)) { warning("Column not found: ", pat); return(rep(NA, nrow(df))) }
  df[[m]]
}
c1 <- Bk$GSE17674; k1 <- c1$clin$diagnosis %in% c("Ewing", "PNET", "Askin")   # exclude normal muscle
c1$clin <- c1$clin[k1, ]; c1$gene <- c1$gene[, k1]
c2 <- Bk$GSE63155; c3 <- Bk$GSE63156
G <- list(GSE17674 = c1$gene, GSE63155 = c2$gene, GSE63156 = c3$gene)
D <- list(
  GSE17674 = data.frame(time = num(pick(c1$clin, "^ovs")),           event = as.integer(pick(c1$clin, "^status") == "Dead"),
                        age = num(pick(c1$clin, "^age"))),
  GSE63155 = data.frame(time = num(pick(c2$clin, "survival.time")),  event = as.integer(pick(c2$clin, "survival.status") == "Dead"),
                        age = num(pick(c2$clin, "^age"))),
  GSE63156 = data.frame(time = num(pick(c3$clin, "^os.time")),       event = num(pick(c3$clin, "^os\\.\\.0")),
                        age = num(pick(c3$clin, "^age"))))
for (n in names(D)) {
  rownames(D[[n]]) <- colnames(G[[n]])
  cat(n, ": tumors", nrow(D[[n]]), "; with survival time", sum(!is.na(D[[n]]$time)), "; events", sum(D[[n]]$event, na.rm = TRUE), "\n")
}
## Follow-up unit: taken from the GEO column name; if absent, inferred from the maximum (> 400 = days)
tcol <- c(GSE17674 = grep("^ovs", colnames(c1$clin), value = TRUE, ignore.case = TRUE)[1],
          GSE63155 = grep("survival.time", colnames(c2$clin), value = TRUE, ignore.case = TRUE)[1],
          GSE63156 = grep("^os.time", colnames(c3$clin), value = TRUE, ignore.case = TRUE)[1])
t_unit <- sapply(names(D), function(n) {
  nm <- tolower(tcol[n])
  if (grepl("day", nm)) "days" else if (grepl("month", nm)) "months" else if (grepl("year", nm)) "years"
  else if (max(D[[n]]$time, na.rm = TRUE) > 400) "days" else "months"
})
cat("\nSurvival time columns and units:\n"); print(data.frame(column_in_GEO = tcol, unit = t_unit,
    max_time = sapply(D, function(d) max(d$time, na.rm = TRUE))))
cat("(Units without day/month/year in the column name are inferred and should be checked on GEO)\n")

## ---- 6.3 Scores: mean of within-cohort z-scores ----
bscore <- function(g, genes) { genes <- intersect(genes, rownames(g)); colMeans(t(apply(g[genes, , drop = FALSE], 1, zs))) }
myeloid_sens <- c("HMOX1", "ENTPD1", "NAMPT", "SAT1")                       # without the broadly expressed FTL/FTH1
prolif       <- c("MKI67", "TOP2A", "PCNA", "MCM2", "CCNB1", "CDK1", "BIRC5")
myel_abund   <- c("CD68", "CD163", "CSF1R", "MS4A7", "LYZ")
for (n in names(D)) {
  g <- G[[n]]; d <- D[[n]]
  d$Tumor_arm   <- bscore(g, tumor_g)
  d$Myeloid_arm <- bscore(g, myeloid_g)
  d$Myeloid_sens <- bscore(g, myeloid_sens)
  d$Myeloid_minus_Tumor <- zs(d$Myeloid_arm) - zs(d$Tumor_arm)
  for (x in c(tumor_g, myeloid_g)) if (x %in% rownames(g)) d[[x]] <- zs(g[x, ])
  d$Polyamine_syn <- bscore(g, c("ODC1", "SRM"))        # post hoc
  d$ENPP1_AIFM2   <- bscore(g, c("ENPP1", "AIFM2"))     # post hoc
  d$Prolif        <- bscore(g, prolif)
  d$Myel_abund    <- bscore(g, myel_abund)
  D[[n]] <- d
}

## ---- 6.4 Cox regression (HR per 1 SD) and inverse-variance fixed-effect meta-analysis ----
cox_one <- function(d, v, covs = character(0)) {
  if (!v %in% colnames(d)) return(c(b = NA, se = NA))   # skip a cohort lacking the gene
  dd <- d; dd$x <- zs(dd[[v]]); for (cv in covs) dd[[cv]] <- zs(dd[[cv]])
  f <- as.formula(paste("Surv(time, event) ~", paste(c("x", covs), collapse = " + ")))
  s <- summary(coxph(f, data = dd))$coefficients
  c(b = s["x", 1], se = s["x", 3])
}
cox_by_cohort <- function(v, covs = character(0)) t(sapply(D, cox_one, v = v, covs = covs))
meta_cox <- function(v, covs = character(0)) {
  r <- cox_by_cohort(v, covs); ok <- !is.na(r[, "b"])
  w <- 1 / r[ok, "se"]^2; B <- sum(w * r[ok, "b"]) / sum(w); SE <- sqrt(1 / sum(w)); Q <- sum(w * (r[ok, "b"] - B)^2)
  data.frame(variable = v, adjusted_for = if (length(covs)) paste(covs, collapse = "+") else "none",
             HR_GSE17674 = exp(r["GSE17674", "b"]), HR_GSE63155 = exp(r["GSE63155", "b"]),
             HR_GSE63156 = exp(r["GSE63156", "b"]),
             HR_meta = exp(B), lo = exp(B - 1.96 * SE), hi = exp(B + 1.96 * SE),
             p_meta = 2 * pnorm(-abs(B / SE)), p_het = pchisq(Q, sum(ok) - 1, lower.tail = FALSE),
             same_direction = sum(sign(r[ok, "b"]) == sign(B)), row.names = NULL)
}

## (1) Fifteen pre-specified variables with BH correction
vars_pre <- c("Tumor_arm", "Myeloid_arm", "Myeloid_sens", "Myeloid_minus_Tumor", tumor_g, myeloid_g)
res_pre <- do.call(rbind, lapply(vars_pre, meta_cox))
res_pre$q_BH <- p.adjust(res_pre$p_meta, "BH")
cat("\nPre-specified variables:\n"); print(res_pre[order(res_pre$p_meta), ], row.names = FALSE, digits = 3)

## (2) Post hoc two-gene scores
res_post <- do.call(rbind, lapply(c("Polyamine_syn", "ENPP1_AIFM2"), meta_cox))
cat("\nPost hoc scores:\n"); print(res_post, row.names = FALSE, digits = 3)

## (3) Correlation with potential confounders
for (n in names(D)) with(D[[n]], cat(n,
  " Myeloid_arm~Myel_abund r =", round(cor(Myeloid_arm, Myel_abund, method = "spearman"), 2),
  " | ENPP1_AIFM2~Myel_abund r =", round(cor(ENPP1_AIFM2, Myel_abund, method = "spearman"), 2),
  " | Polyamine_syn~Prolif r =", round(cor(Polyamine_syn, Prolif, method = "spearman"), 2), "\n"))

## (4) Adjusted models (complete cases)
res_adj <- rbind(
  do.call(rbind, lapply(c("Myeloid_arm", "Myeloid_sens", "Polyamine_syn", "ENPP1_AIFM2", "ODC1", "SRM"), meta_cox, covs = "Myel_abund")),
  do.call(rbind, lapply(c("Polyamine_syn", "ENPP1_AIFM2", "ODC1", "SRM"), meta_cox, covs = "Prolif")),
  do.call(rbind, lapply(c("Polyamine_syn", "ENPP1_AIFM2"), meta_cox, covs = c("Prolif", "Myel_abund"))),
  do.call(rbind, lapply(c("Polyamine_syn", "ENPP1_AIFM2"), meta_cox, covs = c("Prolif", "Myel_abund", "age"))))
cat("\nAdjusted models:\n"); print(res_adj, row.names = FALSE, digits = 3)
write.csv(res_pre,  file.path(dir_res, "bulk_survival_prespecified.csv"), row.names = FALSE)
write.csv(res_post, file.path(dir_res, "bulk_survival_posthoc.csv"),      row.names = FALSE)
write.csv(res_adj,  file.path(dir_res, "bulk_survival_adjusted.csv"),     row.names = FALSE)
saveRDS(list(G = G, D = D), file.path(dir_bulk, "bulk_cohorts_scored.rds"))

## ---- Supplementary Fig. 1: forest plot (exploratory) ----
fdat <- do.call(rbind, lapply(c("Polyamine_syn", "ENPP1_AIFM2", "Myeloid_arm"), function(v) {
  r <- cox_by_cohort(v); w <- 1 / r[, "se"]^2; B <- sum(w * r[, "b"]) / sum(w); SE <- sqrt(1 / sum(w))
  data.frame(variable = v, cohort = c(rownames(r), "Meta (fixed)"), b = c(r[, "b"], B), se = c(r[, "se"], SE))
}))
fdat$HR <- exp(fdat$b); fdat$lo <- exp(fdat$b - 1.96 * fdat$se); fdat$hi <- exp(fdat$b + 1.96 * fdat$se)
fdat$cohort <- factor(fdat$cohort, levels = rev(c(names(D), "Meta (fixed)")))
## Panel order and titles as in the manuscript
fdat$variable <- factor(fdat$variable, levels = c("Polyamine_syn", "ENPP1_AIFM2", "Myeloid_arm"),
                        labels = c("Polyamine-synthesis score (ODC1, SRM)", "ENPP1/AIFM2 score",
                                   "Myeloid-arm score (HMOX1, FTL, FTH1, ENTPD1, NAMPT, SAT1)"))
p_fig5 <- ggplot(fdat, aes(HR, cohort, colour = cohort == "Meta (fixed)")) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi)) + facet_wrap(~ variable, ncol = 1) + scale_x_log10() +
  scale_colour_manual(values = c("grey30", col_t), guide = "none") +
  labs(x = "Hazard ratio per 1 SD (overall survival)", y = NULL) + theme_bw()
save_fig(p_fig5, "FigS1_bulk_forest", 5.5, 6)

## ---- Supplementary Fig. 2: Kaplan-Meier curves, median split of the polyamine-synthesis score ----
km_plot <- function() {
  par(mfrow = c(1, 3))
  for (n in names(D)) {
    d <- D[[n]]; d <- d[!is.na(d$time), ]
    d$grp <- ifelse(d$Polyamine_syn > median(d$Polyamine_syn), "High", "Low")
    p <- summary(coxph(Surv(time, event) ~ I(grp == "High"), data = d))$sctest["pvalue"]
    plot(survfit(Surv(time, event) ~ grp, data = d), col = c(col_t, col_m), lwd = 2,
         xlab = paste0("Time (", t_unit[n], ")"), ylab = "Overall survival",
         main = paste0(n, "  (n=", nrow(d), ")  log-rank p=", signif(p, 2)))
    legend("bottomleft", c("High", "Low"), col = c(col_t, col_m), lwd = 2, bty = "n")
  }
}
pdf(file.path(dir_fig, "FigS2_KM_polyamine_syn.pdf"), width = 12, height = 4); km_plot(); dev.off()
png(file.path(dir_fig, "FigS2_KM_polyamine_syn.png"), width = 12, height = 4, units = "in", res = 300); km_plot(); dev.off()
tiff(file.path(dir_fig, "FigS2_KM_polyamine_syn.tiff"), width = 12, height = 4, units = "in", res = 300, compression = "lzw"); km_plot(); dev.off()


## =====================================================================
## Part 7. Sensitivity analyses (Supplementary Table 2)
## =====================================================================
key_sets <- c("Redox_FerroDefense", "Polyamine", "Pentose_phosphate", "Purinergic")

## Mean log-normalized expression per patient and compartment, and cell-level SD per gene
pb_means <- function(o, pat_col, cells) {
  md <- o@meta.data[cells, ]; comp <- ifelse(md$broad == "Myeloid", "Myeloid", "Tumor")
  f <- factor(paste(md[[pat_col]], comp, sep = "|")); X <- GetAssayData(o, layer = "data")[, cells]
  Gm <- sparseMatrix(i = seq_along(f), j = as.integer(f), x = 1, dims = c(length(f), nlevels(f)),
                     dimnames = list(NULL, levels(f)))
  n <- colSums(Gm); M <- as.matrix(X %*% Gm); M <- sweep(M, 2, n, "/")
  X2 <- X; X2@x <- X2@x^2; mu <- rowMeans(X); sdv <- sqrt(pmax(rowMeans(X2) - mu^2, 1e-8))
  list(M = M, n = n, mu = mu, sd = sdv)
}
## Per gene and patient: (myeloid - tumor) / cell-level SD
pat_D <- function(P, minn = 10) {
  g <- colnames(P$M); pat <- sub("\\|.*", "", g); comp <- sub(".*\\|", "", g); ok <- P$n >= minn
  pats <- intersect(pat[ok & comp == "Myeloid"], pat[ok & comp == "Tumor"])
  Dm <- P$M[, paste0(pats, "|Myeloid"), drop = FALSE] - P$M[, paste0(pats, "|Tumor"), drop = FALSE]
  colnames(Dm) <- pats; list(pats = pats, Z = Dm / P$sd)
}
## (1) Direction of the 11 genes in each patient
gene_check <- function(R) {
  gg <- intersect(c(tumor_g, myeloid_g), rownames(R$Z))
  data.frame(gene = gg, expected = ifelse(gg %in% myeloid_g, "myeloid", "tumor"), n_patients = length(R$pats),
             n_myeloid_higher = sapply(gg, function(x) sum(R$Z[x, ] > 0)), row.names = NULL)
}
## (2) Alternative module score
set_check <- function(R) do.call(rbind, lapply(sc_mod, function(s) {
  g <- intersect(sets_full[[s]], rownames(R$Z)); v <- colMeans(R$Z[g, , drop = FALSE])
  data.frame(module = s, n_genes = length(g), mean_diff = mean(v), n_myeloid_higher = sum(v > 0), n_patients = length(v))
}))
## (3) Comparison with 2,000 size- and expression-matched random gene sets
null_check <- function(R, P, B = 2000) {
  gz <- rowMeans(R$Z); ok <- intersect(names(P$mu)[P$mu > 0.01], rownames(R$Z))
  bin <- cut(rank(P$mu[ok]), 10, labels = FALSE); names(bin) <- ok; pools <- split(ok, bin)
  do.call(rbind, lapply(key_sets, function(s) {
    g <- intersect(sets_full[[s]], ok); obs <- mean(gz[g])
    ex <- lapply(g, function(x) setdiff(pools[[as.character(bin[x])]], g))
    nul <- replicate(B, mean(gz[vapply(ex, function(p) p[sample.int(length(p), 1)], "")]))
    data.frame(module = s, n_genes = length(g), observed = obs, null_mean = mean(nul),
               p_empirical = (1 + sum(abs(nul) >= abs(obs))) / (B + 1))
  }))
}
## (4) Within each patient, five bins of detected genes; myeloid and tumor cells down-sampled to equal numbers per bin
match_cells <- function(o, pat_col, nb = 5) {
  md <- o@meta.data[which(o$broad %in% c("Tumor", "Myeloid")), ]
  unname(unlist(lapply(split(md, md[[pat_col]]), function(d) {
    br <- unique(quantile(d$nFeature_RNA, seq(0, 1, length.out = nb + 1)))
    d$bin <- cut(d$nFeature_RNA, br, include.lowest = TRUE)
    lapply(split(d, d$bin), function(x) {
      a <- rownames(x)[x$broad == "Myeloid"]; b <- rownames(x)[x$broad == "Tumor"]; k <- min(length(a), length(b))
      if (k == 0) character(0) else c(sample(a, k), sample(b, k))
    })
  })))
}
run_all <- function(o, pat_col, label) {
  cells <- colnames(o)[which(o$broad %in% c("Tumor", "Myeloid"))]
  P <- pb_means(o, pat_col, cells); R <- pat_D(P)
  cat("\n#####", label, ": evaluable patients =", length(R$pats), "(", paste(R$pats, collapse = ", "), ")\n")
  out <- list(gene = gene_check(R), set = set_check(R), null = null_check(R, P))
  cm <- match_cells(o, pat_col); Pm <- pb_means(o, pat_col, cm); Rm <- pat_D(Pm)
  cat("Cells after depth matching:", length(cm), "; evaluable patients =", length(Rm$pats), "\nMedian nFeature per compartment after matching:\n")
  print(tapply(o@meta.data[cm, "nFeature_RNA"], o@meta.data[cm, "broad"], median))
  out$gene_matched <- gene_check(Rm); out$set_matched <- set_check(Rm)
  out
}
ROB <- cached("robustness_results.rds", {
  set.seed(1)
  list(A = run_all(obj, "patient", "Cohort A"), B = run_all(obj2, "Sample2", "Cohort B"))
})
for (nm in names(ROB$A)) {
  cat("\n==== A:", nm, "====\n"); print(ROB$A[[nm]], row.names = FALSE, digits = 3)
  cat("\n==== B:", nm, "====\n"); print(ROB$B[[nm]], row.names = FALSE, digits = 3)
}


## =====================================================================
## Part 8. Supplementary tables (one Excel workbook, one sheet per table)
## =====================================================================
wb <- createWorkbook()
add_sheet <- function(name, df) { addWorksheet(wb, name); writeData(wb, name, df) }
add_sheet("S1_gene_sets", data.frame(
  module = names(sets_full), genes_listed = lengths(sets_full),
  present_cohort_A = sapply(sets_full, function(s) length(intersect(s, rownames(obj)))),
  present_cohort_B = sapply(sets_full, function(s) length(intersect(s, rownames(obj2)))),
  genes = sapply(sets_full, paste, collapse = ", "), row.names = NULL))
add_sheet("S2a_patient_modules", rbind(cbind(cohort = "A", resC$A), cbind(cohort = "B", resC$B)))
for (nm in names(ROB$A)) add_sheet(paste0("S2_", nm), rbind(cbind(cohort = "A", ROB$A[[nm]]), cbind(cohort = "B", ROB$B[[nm]])))
add_sheet("S3_reactome_tumor_side",   tum_up)
add_sheet("S3_reactome_myeloid_side", mye_up)
add_sheet("S4_gene_level_DE", tt)
add_sheet("S5_bulk_prespecified", res_pre)
add_sheet("S5_bulk_posthoc",      res_post)
add_sheet("S5_bulk_adjusted",     res_adj)
saveWorkbook(wb, file.path(dir_res, "Supplementary_Tables.xlsx"), overwrite = TRUE)
cat("\nSaved: 02_results/Supplementary_Tables.xlsx\n")


## =====================================================================
## Part 9. Exploratory analysis of the ENPP1 direction (three pre-specified tests; all results reported)
## ENPP1 is tumor-high and can generate immunosuppressive adenosine (cGAMP -> AMP -> adenosine),
## yet it was associated with longer survival in bulk cohorts.
##   E1 Is the adenosine pathway complete? AMP -> adenosine requires NT5E (CD73)
##   E2 Does ENPP1 index tumor-cell content? Correlation with tumor markers; HR adjusted for tumor content
##   E3 Does ENPP1 mark the EWSR1::FLI1-high state? (Franzetti et al. 2017; Aynaud et al. 2020)
## Uses objects created in Parts 0-6.
## =====================================================================
puri_g <- c("ENPP1", "NT5E", "ENTPD1", "CD38", "ADORA2A", "ADORA2B")   # adenosine generation and receptors
tum_mk <- c("CD99", "NKX2-2", "STEAP1")                                 # Ewing sarcoma tumor-cell markers
ews_up <- c("PRKCB", "LIPI", "CCND1", "NR0B1")                          # EWSR1::FLI1-activated targets (Aynaud et al. 2020)
ews_dn <- c("IGFBP3", "LOX", "VIM", "CXCL8", "IL8")                     # EWSR1::FLI1-repressed targets (IL8 = former symbol of CXCL8)

## ---- E1 Fraction of expressing cells and mean expression of adenosine-pathway genes per compartment ----
comp_expr <- function(o, genes, label) {
  o <- subset(o, cells = colnames(o)[!is.na(o$broad)])
  g <- intersect(genes, rownames(o)); X <- GetAssayData(o, layer = "data")[g, , drop = FALSE]
  cl <- split(colnames(o), o$broad)
  do.call(rbind, lapply(names(cl), function(k) data.frame(
    cohort = label, compartment = k, gene = g, n_cells = length(cl[[k]]),
    pct_expr = round(100 * Matrix::rowMeans(X[, cl[[k]], drop = FALSE] > 0), 1),
    mean_expr = round(Matrix::rowMeans(X[, cl[[k]], drop = FALSE]), 3), row.names = NULL)))
}
E1 <- rbind(comp_expr(obj, c(puri_g, tum_mk), "A"), comp_expr(obj2, c(puri_g, tum_mk), "B"))
for (coh in c("A", "B")) {
  cat("\n==== E1 cohort", coh, ": % of cells expressing each gene ====\n")
  print(xtabs(pct_expr ~ compartment + gene, subset(E1, cohort == coh)))
}
## Fraction of tumor cells co-expressing ENPP1 and NT5E
coexp <- function(o, label) {
  cells <- colnames(o)[which(o$broad %in% c("Tumor", "Tumor_cycling"))]
  X <- GetAssayData(o, layer = "data")[intersect(c("ENPP1", "NT5E"), rownames(o)), cells, drop = FALSE]
  data.frame(cohort = label, n_tumor_cells = length(cells),
             pct_ENPP1 = round(100 * mean(X["ENPP1", ] > 0), 1),
             pct_NT5E = if ("NT5E" %in% rownames(X)) round(100 * mean(X["NT5E", ] > 0), 1) else NA,
             pct_both = if ("NT5E" %in% rownames(X)) round(100 * mean(X["ENPP1", ] > 0 & X["NT5E", ] > 0), 1) else NA)
}
E1co <- rbind(coexp(obj, "A"), coexp(obj2, "B"))
cat("\n==== E1 ENPP1 / NT5E co-expression in tumor cells ====\n"); print(E1co, row.names = FALSE)

## Supplementary Fig. 3: adenosine-pathway genes, tumor markers and EWSR1::FLI1-activated targets across compartments
dot2 <- function(o, ttl) {
  o <- subset(o, cells = colnames(o)[!is.na(o$broad)]); o$broad <- factor(o$broad, levels = lv)
  DotPlot(o, features = intersect(c(puri_g, tum_mk, ews_up), rownames(o)), group.by = "broad",
          cols = c("grey90", "#333333")) + ggtitle(ttl) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank())
}
save_fig(dot2(obj, "Cohort A (GSE243347)") / dot2(obj2, "Cohort B (GSE277083)"), "FigS3_ENPP1_context", 9.5, 8.5)

## ---- E3 (single-cell): correlation of ENPP1 with EWSR1::FLI1 activity within tumor cells of each patient ----
ews_sc <- function(o, pat_col, label, minn = 30) {
  cells <- colnames(o)[which(o$broad %in% c("Tumor", "Tumor_cycling"))]
  gg <- intersect(c("ENPP1", ews_up, ews_dn), rownames(o))
  X <- as.matrix(GetAssayData(o, layer = "data")[gg, cells, drop = FALSE])
  arm <- function(g) {
    g <- intersect(g, rownames(X)); g <- g[apply(X[g, , drop = FALSE], 1, sd) > 0]
    cat(label, "genes used:", paste(g, collapse = ", "), "\n")
    colMeans(t(scale(t(X[g, , drop = FALSE]))))
  }
  up <- arm(ews_up); dn <- arm(ews_dn); act <- up - dn      # activity = activated minus repressed targets
  pat <- o@meta.data[cells, pat_col]
  do.call(rbind, lapply(split(seq_along(cells), pat), function(i) {
    if (length(i) < minn) return(NULL)
    e <- X["ENPP1", i]
    data.frame(cohort = label, patient = pat[i[1]], n_tumor_cells = length(i),
               pct_ENPP1 = round(100 * mean(e > 0), 1),
               rho_activity = suppressWarnings(cor(e, act[i], method = "spearman")),
               rho_up_arm   = suppressWarnings(cor(e, up[i],  method = "spearman")),
               rho_down_arm = suppressWarnings(cor(e, dn[i],  method = "spearman")))
  }))
}
E3sc <- rbind(ews_sc(obj, "patient", "A"), ews_sc(obj2, "Sample2", "B"))
cat("\n==== E3 single-cell: Spearman correlation of ENPP1 with EWSR1::FLI1 activity in tumor cells, per patient ====\n")
print(E3sc, row.names = FALSE, digits = 2)
for (coh in c("A", "B")) with(subset(E3sc, cohort == coh & !is.na(rho_activity)),
  cat("Cohort", coh, ": patients", length(rho_activity), "; median rho", round(median(rho_activity), 2),
      "; patients with rho > 0:", sum(rho_activity > 0), "\n"))

## ---- E2 + E3 (bulk) ----
for (n in names(D)) {
  D[[n]]$Tumor_content <- bscore(G[[n]], tum_mk)
  D[[n]]$EWS_up        <- bscore(G[[n]], ews_up)   # bulk: activated targets only (repressed targets are also expressed by stroma)
}
E2cor <- do.call(rbind, lapply(names(D), function(n) with(D[[n]], data.frame(
  cohort = n,
  ENPP1_vs_tumor_content      = cor(ENPP1, Tumor_content, method = "spearman"),
  ENPP1_vs_EWS_up             = cor(ENPP1, EWS_up,        method = "spearman"),
  ENPP1_vs_myeloid_abundance  = cor(ENPP1, Myel_abund,    method = "spearman"),
  ENPP1_vs_proliferation      = cor(ENPP1, Prolif,        method = "spearman"),
  tumor_content_vs_myeloid    = cor(Tumor_content, Myel_abund, method = "spearman")))))
cat("\n==== E2/E3 bulk: Spearman correlations ====\n"); print(E2cor, row.names = FALSE, digits = 2)
E23cox <- rbind(meta_cox("ENPP1"), meta_cox("ENPP1", "Tumor_content"), meta_cox("ENPP1", "EWS_up"),
                meta_cox("ENPP1_AIFM2"), meta_cox("ENPP1_AIFM2", "Tumor_content"), meta_cox("ENPP1_AIFM2", "EWS_up"))
cat("\n==== E2/E3 bulk: HR adjusted for tumor content or EWSR1::FLI1 activity (two coefficients per model) ====\n")
print(E23cox, row.names = FALSE, digits = 3)

## ---- Save ----
saveRDS(list(E1 = E1, E1_coexpression = E1co, E3_singlecell = E3sc, E2_bulk_cor = E2cor, E23_cox = E23cox),
        file.path(dir_res, "ENPP1_explanatory.rds"))
wb <- loadWorkbook(file.path(dir_res, "Supplementary_Tables.xlsx"))
for (sh in list(list("S6_ENPP1_E1_compartments", E1), list("S6_ENPP1_E1_coexpression", E1co),
                list("S6_ENPP1_E3_singlecell", E3sc), list("S6_ENPP1_bulk_correlation", E2cor),
                list("S6_ENPP1_bulk_cox", E23cox))) {
  if (sh[[1]] %in% names(wb)) removeWorksheet(wb, sh[[1]])
  addWorksheet(wb, sh[[1]]); writeData(wb, sh[[1]], sh[[2]])
}
saveWorkbook(wb, file.path(dir_res, "Supplementary_Tables.xlsx"), overwrite = TRUE)
cat("\nENPP1 analyses saved: 02_results/ENPP1_explanatory.rds, sheets S6_*, Supplementary Fig. 3\n")

writeLines(capture.output(sessionInfo()), file.path(dir_res, "sessionInfo.txt"))
cat("\nDone. Figures in ES_project/03_figures, tables in ES_project/02_results.\n")
