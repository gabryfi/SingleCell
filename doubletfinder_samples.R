suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
})

# =========================
# PATHS
# =========================
base_dir <- "/media/user/8Tb/scRNAseq/seurat_analysis/data_matrices"
out_dir  <- file.path(base_dir, "doubletfinder_output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sample_dirs <- list(
  adeno = file.path(base_dir, "adeno_data"),
  sham  = file.path(base_dir, "sham_data")
)

# =========================
# PARAMETRI SCELTI
# =========================
min_features <- 200
min_cells <- 3

# Scelta iniziale ragionevole
pcs_use <- 1:15

# Stima iniziale doublet rate
doublet_rate <- 0.075   # 7.5%

# Risoluzione clustering per stimare homotypic proportion
cluster_resolution <- 0.5

# =========================
# FUNZIONE PER UN CAMPIONE
# =========================
run_doubletfinder_one_sample <- function(sample_name, sample_path) {

  cat("\n==============================\n")
  cat("Processing sample:", sample_name, "\n")
  cat("Input path:", sample_path, "\n")
  cat("==============================\n")

  # 1. Leggi matrice 10X
  counts <- Read10X(data.dir = sample_path)

  # 2. Crea Seurat object
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = min_cells,
    min.features = min_features
  )

  # 3. QC di base
  # Per mouse, i geni mitocondriali iniziano spesso con "mt-"
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

  cat("Cells after object creation:", ncol(seu), "\n")
  cat("Genes:", nrow(seu), "\n")

  # 4. Filtro base
  # Questi cutoff sono prudenti e generici; puoi cambiarli dopo aver visto i plot QC
  seu <- subset(
    seu,
    subset = nFeature_RNA >= 200 &
             nFeature_RNA <= 7500 &
             percent.mt <= 15
  )

  cat("Cells after QC filtering:", ncol(seu), "\n")

  # 5. Preprocessing classico Seurat
  seu <- NormalizeData(seu, verbose = FALSE)
  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  seu <- ScaleData(seu, verbose = FALSE)
  seu <- RunPCA(seu, npcs = max(pcs_use), verbose = FALSE)

  # 6. Clustering preliminare per correzione omotipica
  seu <- FindNeighbors(seu, dims = pcs_use, verbose = FALSE)
  seu <- FindClusters(seu, resolution = cluster_resolution, verbose = FALSE)

  # 7. Sweep per scegliere pK
  cat("Running paramSweep...\n")
  sweep.res.list <- paramSweep_v3(seu, PCs = pcs_use, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  # Trova il pK col massimo BCmetric
  bcmvn$BCmetric <- as.numeric(as.character(bcmvn$BCmetric))
  best_pk_row <- bcmvn[which.max(bcmvn$BCmetric), , drop = FALSE]
  best_pK <- as.numeric(as.character(best_pk_row$pK))

  cat("Selected pK:", best_pK, "\n")

  # 8. Stima nExp
  nExp_poi <- round(doublet_rate * ncol(seu))
  cat("Initial nExp:", nExp_poi, "\n")

  # 9. Correzione omotipica
  homotypic.prop <- modelHomotypic(seu$seurat_clusters)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  cat("Homotypic proportion:", homotypic.prop, "\n")
  cat("Adjusted nExp:", nExp_poi.adj, "\n")

  # 10. Primo run DoubletFinder
  seu <- doubletFinder_v3(
    seu,
    PCs = pcs_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_poi,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  # Nome colonna pANN appena creata
  pANN_col <- paste("pANN", 0.25, best_pK, nExp_poi, sep = "_")

  # 11. Secondo run con nExp corretto
  seu <- doubletFinder_v3(
    seu,
    PCs = pcs_use,
    pN = 0.25,
    pK = best_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = pANN_col,
    sct = FALSE
  )

  # 12. Trova automaticamente la colonna classificazione finale
  class_col <- paste("DF.classifications", 0.25, best_pK, nExp_poi.adj, sep = "_")

  cat("Final classification column:", class_col, "\n")
  print(table(seu@meta.data[[class_col]]))

  # 13. Salva metadata completi
  meta_file <- file.path(out_dir, paste0(sample_name, "_metadata_doubletfinder.csv"))
  write.csv(seu@meta.data, meta_file, row.names = TRUE)

  # 14. Salva oggetto completo
  rds_file <- file.path(out_dir, paste0(sample_name, "_seurat_doubletfinder.rds"))
  saveRDS(seu, rds_file)

  # 15. Crea versione filtrata senza doublets
  seu_singlets <- subset(seu, subset = get(class_col) == "Singlet")

  singlets_rds <- file.path(out_dir, paste0(sample_name, "_singlets_only.rds"))
  saveRDS(seu_singlets, singlets_rds)

  singlet_meta_file <- file.path(out_dir, paste0(sample_name, "_singlets_metadata.csv"))
  write.csv(seu_singlets@meta.data, singlet_meta_file, row.names = TRUE)

  cat("Saved:\n")
  cat(" -", meta_file, "\n")
  cat(" -", rds_file, "\n")
  cat(" -", singlets_rds, "\n")
  cat(" -", singlet_meta_file, "\n")

  return(list(
    seu = seu,
    seu_singlets = seu_singlets,
    best_pK = best_pK,
    nExp_poi = nExp_poi,
    nExp_poi.adj = nExp_poi.adj,
    class_col = class_col
  ))
}

# =========================
# RUN SU TUTTI I CAMPIONI
# =========================
results <- lapply(names(sample_dirs), function(s) {
  run_doubletfinder_one_sample(s, sample_dirs[[s]])
})
names(results) <- names(sample_dirs)

cat("\nAll done.\n")
