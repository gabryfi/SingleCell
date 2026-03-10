#####################################################################
#### Run Scrublet con threshold manuale 0.3
import scrublet as scr
import scipy.io
import numpy as np
import os
import matplotlib.pyplot as plt
import pandas as pd
import gzip

# Directory di input (cartelle tipo BP1_skin, HC2_skin, ecc.)
base_dir = "/home/utenteunifi/data/03_projects_open/pemfigoide_singlecell/10924853/"

# Directory di output per i risultati
output_dir = os.path.join(base_dir, "scrublet_results")
os.makedirs(output_dir, exist_ok=True)

# Lista delle directories da processare (tutte le sottocartelle di base_dir)
directories = [
    os.path.join(base_dir, d) for d in os.listdir(base_dir)
    if os.path.isdir(os.path.join(base_dir, d))
]

# Threshold manuale per identificare doublets
threshold = 0.25

# Funzione per processare un singolo dataset
def process_dataset(input_dir, output_dir):
    print(f"Processing {input_dir}...")

    try:
        # Percorsi ai file
        matrix_path = os.path.join(input_dir, 'matrix.mtx.gz')
        features_path = os.path.join(input_dir, 'features.tsv.gz')

        if not os.path.exists(matrix_path) or not os.path.exists(features_path):
            print(f"Files missing in {input_dir}. Skipping...")
            return

        # Carica i dati
        with gzip.open(matrix_path, 'rb') as f:
            counts_matrix = scipy.io.mmread(f).T.tocsc()

        with gzip.open(features_path, 'rt') as f:
            features = np.array([line.strip().split('\t')[1] for line in f])

        print(f"Counts matrix shape: {counts_matrix.shape[0]} cells, {counts_matrix.shape[1]} genes")
        print(f"Number of features in feature list: {len(features)}")

        # Inizializza Scrublet
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.07)  # soglia per cellule caricate=10K

        # Calcola probabilità di doublet
        doublet_scores, _ = scrub.scrub_doublets()

        # Applica threshold manuale
        predicted_doublets = doublet_scores > threshold

        # Salva i risultati
        results_df = pd.DataFrame({
            'DoubletScore': doublet_scores,
            'PredictedDoublet': predicted_doublets
        })
        results_csv = os.path.join(output_dir, f"{os.path.basename(input_dir)}_doublet_scores.csv")
        results_df.to_csv(results_csv, index=False)
        print(f"Saved results to {results_csv}")

        # Genera e salva histogram
        plt.figure()
        scrub.plot_histogram()
        histogram_path = os.path.join(output_dir, f"{os.path.basename(input_dir)}_histogram.png")
        plt.savefig(histogram_path)
        plt.close()

        # Crea embedding UMAP e visualizza
        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
        embedding = scrub.embedding_

        plt.figure(figsize=(8, 6))
        plt.scatter(embedding[:, 0], embedding[:, 1], c=predicted_doublets, cmap='coolwarm', s=1)
        plt.colorbar(label='Doublet Prediction')
        plt.title(f'UMAP Embedding for {os.path.basename(input_dir)}')
        umap_path = os.path.join(output_dir, f"{os.path.basename(input_dir)}_umap.png")
        plt.savefig(umap_path)
        plt.close()

        print(f"Plots saved for {input_dir}")

    except Exception as e:
        print(f"Error processing {input_dir}: {e}")

# Ciclo su tutte le directories
for dir_name in directories:
    process_dataset(dir_name, output_dir)

print("Processing complete.")
