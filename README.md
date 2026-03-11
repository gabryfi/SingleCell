# SingleCell
# 10X 3' GEM - Nanopore



@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
sample_dirs <- list(
  adeno = "/media/user/8Tb/scRNAseq/adeno_wf/adeno/adeno.gene_raw_feature_bc_matrix",
  sham  = "/media/user/8Tb/scRNAseq/sham_wf/sham/sham.gene_raw_feature_bc_matrix"
)

out_dir <- "/media/user/8Tb/scRNAseq/seurat_analysis/doubletfinder_raw_output"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
min_cells <- 3
min_features <- 200
max_features <- 6000
max_mt <- 10
pcs_use <- 1:15
cluster_resolution <- 0.5
doublet_rate <- 0.075
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# in_cells <- 3
Quando crei l’oggetto Seurat, un gene viene tenuto solo se è espresso in almeno 3 cellule.

# min_features <- 200
Una cellula deve avere almeno 200 geni rilevati per essere tenuta.

# max_features <- 6000
Una cellula con più di 6000 geni rilevati viene esclusa.
Perché
Spesso cellule con troppi geni sono sospette per:
doubletsl, multiplets, artefatti

# max_mt <- 10
Esclude cellule con più del 10% di reads mitocondriali.
Di solito un contenuto mitocondriale alto suggerisce:
cellula stressata, danneggiata o morente

# pcs_use <- 1:15
Userà le prime 15 componenti principali per:
neighbors, clustering, DoubletFinder

# cluster_resolution <- 0.5
Parametro del clustering Seurat.
Valori più alti → più cluster
Valori più bassi → meno cluster

# doublet_rate <- 0.075
Assume un tasso atteso di doublets del 7.5%. (10K cellule caricate = 6-8% doublet_rate)

####################################

