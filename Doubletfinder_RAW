# DOUBLET FINDER

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
})

# =========================
# PATHS
# =========================
sample_dirs <- list(
  adeno = "/media/user/8Tb/scRNAseq/adeno_wf/adeno/adeno.gene_raw_feature_bc_matrix",
  sham  = "/media/user/8Tb/scRNAseq/sham_wf/sham/sham.gene_raw_feature_bc_matrix"
)

out_dir <- "/media/user/8Tb/scRNAseq/seurat_analysis/doubletfinder_raw_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# =========================
# PARAMETRI
# =========================
min_cells <- 3
min_features <- 200
max_features <- 6000
max_mt <- 10
pcs_use <- 1:15
cluster_resolution <- 0.5
doublet_rate <- 0.075

# =========================
# COMPATIBILITÀ FUNZIONI DF
# =========================
has_paramSweep_v3 <- exists("paramSweep_v3")
has_paramSweep    <- exists("paramSweep")
has_df_v3         <- exists("doubletFinder_v3")
has_df            <- exists("doubletFinder")

if (!has_df_v3 && !has_df) {
  stop("Non trovo né doubletFinder_v3 né doubletFinder.")
}
if (!has_paramSweep_v3 && !has_paramSweep) {
  stop("Non trovo né paramSweep_v3 né paramSweep.")
}

run_paramSweep <- function(seu, PCs, sct = FALSE) {
  if (has_paramSweep_v3) {
    paramSweep_v3(seu, PCs = PCs, sct = sct)
  } else {
    paramSweep(seu, PCs = PCs, sct = sct)
  }
}

run_doubletFinder <- function(seu, PCs, pN, pK, nExp, reuse.pANN = NULL, sct = FALSE) {
  if (has_df_v3) {
    doubletFinder_v3(
      seu,
      PCs = PCs,
      pN = pN,
      pK = pK,
      nExp = nExp,
      reuse.pANN = reuse.pANN,
      sct = sct
    )
  } else {
    doubletFinder(
      seu,
      PCs = PCs,
      pN = pN,
      pK = pK,
      nExp = nExp,
      reuse.pANN = reuse.pANN,
      sct = sct
    )
  }
}

# =========================
# FUNZIONE PER CAMPIONE
# =========================
run_one_sample <- function(sample_name, sample_path) {
  cat("\n==============================\n")
  cat("PROCESSING:", sample_name, "\n")
  cat("PATH:", sample_path, "\n")
  cat("==============================\n")

  counts <- Read10X(data.dir = sample_path)

  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = min_cells,
    min.features = 0
  )

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

  cat("\nBefore QC\n")
  cat("Cells:", ncol(seu), "\n")
  cat("Genes:", nrow(seu), "\n")

  seu <- subset(
    seu,
    subset =
      nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      percent.mt <= max_mt
  )

  cat("\nAfter QC\n")
  cat("Cells:", ncol(seu), "\n")
  cat("Genes:", nrow(seu), "\n")

  cat("\nQC quantiles\n")
  cat("nFeature_RNA:\n")
  print(round(quantile(seu$nFeature_RNA, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)), 2))
  cat("nCount_RNA:\n")
  print(round(quantile(seu$nCount_RNA, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)), 2))
  cat("percent.mt:\n")
  print(round(quantile(seu$percent.mt, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1), na.rm = TRUE), 2))

  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = max(pcs_use), verbose = FALSE)

  seu <- FindNeighbors(seu, dims = pcs_use, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)

  cat("\nRunning paramSweep...\n")
  sweep.res.list <- run_paramSweep(seu, PCs = pcs_use, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  bcmvn$BCmetric <- as.numeric(as.character(bcmvn$BCmetric))
  best_pk_row <- bcmvn[which.max(bcmvn$BCmetric), , drop = FALSE]
  best_pK <- as.numeric(as.character(best_pk_row$pK))

  cat("Selected pK:", best_pK, "\n")

  nExp_poi <- round(doublet_rate * ncol(seu))
  cat("Initial nExp:", nExp_poi, "\n")

  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  cat("Homotypic proportion:", homotypic.prop, "\n")
  cat("Adjusted nExp:", nExp_poi.adj, "\n")

  seu <- run_doubletFinder(
    seu,
    PCs = pcs_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_poi,
    reuse.pANN = NULL,
    sct = FALSE
  )

  pANN_col <- paste("pANN", 0.25, best_pK, nExp_poi, sep = "_")

  seu <- run_doubletFinder(
    seu,
    PCs = pcs_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = pANN_col,
    sct = FALSE
  )

  class_col <- paste("DF.classifications", 0.25, best_pK, nExp_poi.adj, sep = "_")

  cat("\nFinal classification column:", class_col, "\n")
  print(table(seu@meta.data[[class_col]]))

  meta_file <- file.path(out_dir, paste0(sample_name, "_metadata_doubletfinder_raw.csv"))
  write.csv(seu@meta.data, meta_file, row.names = TRUE)

  rds_file <- file.path(out_dir, paste0(sample_name, "_seurat_doubletfinder_raw.rds"))
  saveRDS(seu, rds_file)

  singlet_cells <- rownames(seu@meta.data)[seu@meta.data[[class_col]] == "Singlet"]
  seu_singlets <- subset(seu, cells = singlet_cells)

  singlets_rds <- file.path(out_dir, paste0(sample_name, "_singlets_only_raw.rds"))
  saveRDS(seu_singlets, singlets_rds)

  singlet_meta_file <- file.path(out_dir, paste0(sample_name, "_singlets_metadata_raw.csv"))
  write.csv(seu_singlets@meta.data, singlet_meta_file, row.names = TRUE)

  cat("\nSaved files:\n")
  cat(meta_file, "\n")
  cat(rds_file, "\n")
  cat(singlets_rds, "\n")
  cat(singlet_meta_file, "\n")

  invisible(list(
    seu = seu,
    seu_singlets = seu_singlets,
    best_pK = best_pK,
    nExp_poi = nExp_poi,
    nExp_poi.adj = nExp_poi.adj,
    class_col = class_col
  ))
}

results <- lapply(names(sample_dirs), function(s) {
  run_one_sample(s, sample_dirs[[s]])
})
names(results) <- names(sample_dirs)

cat("\n==============================\n")
cat("ALL DONE\n")
cat("==============================\n")


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



