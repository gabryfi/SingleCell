############################################################
# SoupX workflow con rho modificabile
#
# MODIFICA SOLO QUESTE 3 VARIABILI:
############################################################

RAW_DIR <- "adeno1_sup.gene_raw_feature_bc_matrix"
PROCESSED_DIR <- "adeno1_sup.gene_processed_feature_bc_matrix"

# Metti un numero per forzare rho manuale, es:

# 0.011 = come stima automatica precedente
# 0.02  = più pesante
# 0.03  = ancora più pesante

# MANUAL_RHO <- NA    = usa autoEstCont() senza forzare rho
MANUAL_RHO <- 0.02

############################################################
# Da qui in poi NON modificare
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

required_pkgs <- c("SoupX", "Seurat", "Matrix")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs)
}

suppressPackageStartupMessages({
  library(SoupX)
  library(Seurat)
  library(Matrix)
})

sample_name <- basename(getwd())

rho_label <- if (is.na(MANUAL_RHO)) {
  "AUTO_RHO"
} else {
  paste0("MANUAL_RHO_", gsub("\\.", "p", as.character(MANUAL_RHO)))
}

out_root <- paste0("SoupX_results_", sample_name, "_", rho_label)
out_matrix_dir <- file.path(out_root, "soupx_corrected_feature_bc_matrix")
qc_dir <- file.path(out_root, "QC")

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(out_matrix_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

check_10x_dir <- function(d) {
  missing <- required_files[!file.exists(file.path(d, required_files))]
  if (length(missing) > 0) {
    stop(
      "ERRORE: nella cartella ", d, " mancano: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

read_10x_matrix <- function(d) {
  m <- Seurat::Read10X(data.dir = d)

  if (is.list(m)) {
    if ("Gene Expression" %in% names(m)) {
      m <- m[["Gene Expression"]]
    } else {
      m <- m[[1]]
    }
  }

  m
}

read_barcodes <- function(d) {
  readLines(gzfile(file.path(d, "barcodes.tsv.gz")))
}

is_integer_sparse <- function(mat, n_check = 100000) {
  vals <- mat@x
  vals <- vals[seq_len(min(length(vals), n_check))]
  all(abs(vals - round(vals)) < 1e-8)
}

message("Controllo file input...")
check_10x_dir(RAW_DIR)
check_10x_dir(PROCESSED_DIR)

message("Cartella di lavoro:")
message(getwd())

message("RAW_DIR:")
message(RAW_DIR)

message("PROCESSED_DIR:")
message(PROCESSED_DIR)

message("MANUAL_RHO:")
print(MANUAL_RHO)

message("Output:")
message(out_root)

message("Lettura RAW matrix...")
tod <- read_10x_matrix(RAW_DIR)

message("Lettura barcode filtrati dalla PROCESSED...")
processed_barcodes <- read_barcodes(PROCESSED_DIR)

message("Dimensioni RAW/tod:")
print(dim(tod))

raw_is_integer <- is_integer_sparse(tod)

message("RAW contiene valori interi?")
print(raw_is_integer)

if (!raw_is_integer) {
  stop("ERRORE: la RAW matrix non contiene conteggi interi. Non procedo.", call. = FALSE)
}

missing_barcodes <- setdiff(processed_barcodes, colnames(tod))

if (length(missing_barcodes) > 0) {
  stop(
    "ERRORE: alcuni barcode della PROCESSED non sono presenti nella RAW. Esempi: ",
    paste(head(missing_barcodes, 10), collapse = ", "),
    call. = FALSE
  )
}

message("Costruzione toc corretta: RAW filtrata sui barcode della PROCESSED...")
toc <- tod[, processed_barcodes, drop = FALSE]

message("Dimensioni toc:")
print(dim(toc))

toc_is_integer <- is_integer_sparse(toc)

message("toc contiene valori interi?")
print(toc_is_integer)

if (!toc_is_integer) {
  stop("ERRORE: toc non è intera. Non procedo.", call. = FALSE)
}

input_summary <- data.frame(
  matrix = c("tod_raw_unfiltered", "toc_raw_filtered_barcodes"),
  genes = c(nrow(tod), nrow(toc)),
  cells = c(ncol(tod), ncol(toc)),
  total_counts = c(sum(tod), sum(toc)),
  nonzero_entries = c(length(tod@x), length(toc@x)),
  values_are_integer = c(raw_is_integer, toc_is_integer)
)

write.table(
  input_summary,
  file = file.path(qc_dir, "input_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Input summary:")
print(input_summary)

message("Clustering Seurat minimo sulla toc...")

seu <- CreateSeuratObject(
  counts = toc,
  project = sample_name,
  min.cells = 0,
  min.features = 0
)

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

clusters <- setNames(as.character(seu$seurat_clusters), colnames(seu))

write.table(
  as.data.frame(table(clusters)),
  file = file.path(qc_dir, "seurat_clusters_used_by_soupx.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Creazione SoupChannel...")
sc <- SoupChannel(tod = tod, toc = toc)
sc <- setClusters(sc, clusters[colnames(toc)])

message("Stima automatica contaminazione con autoEstCont(), usata come diagnosi...")
pdf(file.path(qc_dir, paste0(sample_name, ".soupx_autoEstCont.pdf")))
sc <- autoEstCont(sc)
dev.off()

auto_rho <- sc$fit$rhoEst

message("Rho automatico stimato da SoupX:")
print(auto_rho)

if (!is.null(sc$fit$markersUsed)) {
  write.table(
    sc$fit$markersUsed,
    file = file.path(qc_dir, "soupx_markers_used_autoEstCont.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

if (!is.na(MANUAL_RHO)) {
  message("Forzo rho manuale:")
  print(MANUAL_RHO)

  sc <- setContaminationFraction(
    sc,
    contFrac = MANUAL_RHO,
    forceAccept = TRUE
  )

  final_rho <- MANUAL_RHO
} else {
  message("Uso rho automatico stimato da autoEstCont().")
  final_rho <- auto_rho
}

message("Rho finale usato per adjustCounts:")
print(final_rho)

message("Correzione conteggi...")
corrected <- adjustCounts(sc, roundToInt = TRUE)

raw_total <- sum(toc)
corrected_total <- sum(corrected)
removed_total <- raw_total - corrected_total
removed_percent <- 100 * removed_total / raw_total

message("Raw/toc total counts: ", raw_total)
message("Corrected total counts: ", corrected_total)
message("Removed total counts: ", removed_total)
message("Removed percent: ", removed_percent)

if (corrected_total <= 0) {
  stop("ERRORE GRAVE: matrice corretta azzerata. Risultato non valido.", call. = FALSE)
}

if (removed_percent > 50) {
  warning(
    "ATTENZIONE: SoupX ha rimosso più del 50% dei conteggi totali: ",
    round(removed_percent, 2),
    "%. Controlla bene il risultato."
  )
}

message("QC prima/dopo...")

cell_qc <- data.frame(
  barcode = colnames(toc),
  raw_counts = as.numeric(Matrix::colSums(toc)),
  corrected_counts = as.numeric(Matrix::colSums(corrected))
)

cell_qc$removed_counts <- cell_qc$raw_counts - cell_qc$corrected_counts
cell_qc$removed_percent <- 100 * cell_qc$removed_counts / cell_qc$raw_counts

write.table(
  cell_qc,
  file = file.path(qc_dir, "cell_level_before_after_soupx.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

gene_raw <- Matrix::rowSums(toc)
gene_corrected <- Matrix::rowSums(corrected)

gene_qc <- data.frame(
  gene = rownames(toc),
  raw_counts = as.numeric(gene_raw),
  corrected_counts = as.numeric(gene_corrected),
  removed_counts = as.numeric(gene_raw - gene_corrected)
)

gene_qc$removed_percent_of_gene_raw <- 100 * gene_qc$removed_counts / gene_qc$raw_counts
gene_qc <- gene_qc[is.finite(gene_qc$removed_percent_of_gene_raw), ]
gene_qc <- gene_qc[order(gene_qc$removed_counts, decreasing = TRUE), ]

write.table(
  head(gene_qc, 300),
  file = file.path(qc_dir, "top300_removed_genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

marker_check <- c(
  "Snap25", "Syt1", "Tubb3", "Rbfox3",
  "Nefl", "Nefm", "Nefh",
  "Calca", "Tac1", "Trpv1", "Scn10a", "Scn11a",
  "Nppb", "Mrgprd", "Mrgprb5", "Sst",
  "Pvalb", "Th", "Fgf13", "Nrg3", "Stmn2",
  "Mbp", "Mpz", "Plp1", "Prx", "Mag", "Cnp", "Gldn", "Ugt8a",
  "Fabp7", "Apoe",
  "Ptprc", "Lyz2", "Cd74", "H2-Aa", "H2-Ab1", "C1qa", "C1qb", "C1qc",
  "Pecam1", "Kdr",
  "Dcn", "Col1a1",
  "Hbb-bs", "Hbb-bt", "Hba-a1", "Hba-a2", "Alas2", "Bpgm"
)

marker_qc <- gene_qc[gene_qc$gene %in% marker_check, ]

write.table(
  marker_qc,
  file = file.path(qc_dir, "biological_marker_before_after_soupx.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Creazione PDF QC...")

pdf(file.path(qc_dir, "soupx_QC_before_after.pdf"), width = 9, height = 7)

hist(
  cell_qc$removed_percent,
  breaks = 60,
  main = paste0("SoupX - % conteggi rimossi per cellula | rho=", final_rho),
  xlab = "% conteggi rimossi",
  ylab = "Numero cellule"
)

plot(
  cell_qc$raw_counts,
  cell_qc$corrected_counts,
  pch = 16,
  cex = 0.35,
  xlab = "Conteggi originali per cellula",
  ylab = "Conteggi corretti per cellula",
  main = paste0("Prima vs dopo SoupX | rho=", final_rho)
)
abline(0, 1, lwd = 2)

boxplot(
  cell_qc$removed_percent,
  horizontal = TRUE,
  main = "Distribuzione % rimossa per cellula",
  xlab = "% conteggi rimossi"
)

barplot(
  head(gene_qc$removed_counts, 30),
  names.arg = head(gene_qc$gene, 30),
  las = 2,
  cex.names = 0.55,
  main = "Top 30 geni rimossi da SoupX",
  ylab = "Conteggi rimossi"
)

dev.off()

message("Salvataggio risultati...")

saveRDS(
  sc,
  file = file.path(out_root, paste0(sample_name, ".soupx_channel.rds"))
)

saveRDS(
  corrected,
  file = file.path(out_root, paste0(sample_name, ".soupx_corrected_counts.rds"))
)

saveRDS(
  seu,
  file = file.path(out_root, paste0(sample_name, ".seurat_for_soupx_clusters.rds"))
)

Matrix::writeMM(
  corrected,
  file = file.path(out_matrix_dir, "matrix.mtx")
)

file.copy(
  from = file.path(PROCESSED_DIR, "barcodes.tsv.gz"),
  to = file.path(out_matrix_dir, "barcodes.tsv.gz"),
  overwrite = TRUE
)

file.copy(
  from = file.path(RAW_DIR, "features.tsv.gz"),
  to = file.path(out_matrix_dir, "features.tsv.gz"),
  overwrite = TRUE
)

system2("gzip", args = c("-f", file.path(out_matrix_dir, "matrix.mtx")))

capture.output({
  cat("Sample:", sample_name, "\n")
  cat("RAW_DIR:", RAW_DIR, "\n")
  cat("PROCESSED_DIR used only for barcodes:", PROCESSED_DIR, "\n")
  cat("MANUAL_RHO:", MANUAL_RHO, "\n")
  cat("auto_rho:", auto_rho, "\n")
  cat("final_rho_used:", final_rho, "\n")
  cat("Output dir:", out_root, "\n\n")

  cat("Input summary:\n")
  print(input_summary)

  cat("\nSoupX fit:\n")
  print(sc$fit)

  cat("\nRaw/toc total counts:", raw_total, "\n")
  cat("Corrected total counts:", corrected_total, "\n")
  cat("Removed total counts:", removed_total, "\n")
  cat("Removed percent:", removed_percent, "\n\n")

  cat("Cell removed percent summary:\n")
  print(summary(cell_qc$removed_percent))

  cat("\nNumber of cells removed 100%:\n")
  print(sum(cell_qc$removed_percent == 100))

  cat("\nBiological marker QC:\n")
  print(marker_qc)

  cat("\nTop 50 removed genes:\n")
  print(head(gene_qc, 50))

  cat("\nSession info:\n")
  print(sessionInfo())
}, file = file.path(qc_dir, "soupx_summary.txt"))

message("FINITO.")
message("Output salvato in:")
message(normalizePath(out_root))
message("Matrice corretta:")
message(normalizePath(out_matrix_dir))
message("Summary:")
message(file.path(qc_dir, "soupx_summary.txt"))
