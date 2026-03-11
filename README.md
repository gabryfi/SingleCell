# SingleCell
# 10X 3' GEM - Nanopore

####################################

min_cells <- 3
min_features <- 200
max_features <- 6000
max_mt <- 10
pcs_use <- 1:15
cluster_resolution <- 0.5
doublet_rate <- 0.075

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

