#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#####################################################################
# Scrublet indipendente su RAW gene matrices corrette
#
# PRINCIPI:
# - usa SOLO raw integer matrices
# - fa un prefiltro minimo per togliere barcode chiaramente non cellulari
# - NON impone un cutoff alto su nFeature/nCount prima di Scrublet
# - usa come chiamata principale quella AUTOMATICA di Scrublet
# - salva anche barcode e metriche QC
#####################################################################

import os
import gzip
import numpy as np
import pandas as pd
import scipy.io
import scrublet as scr
import matplotlib.pyplot as plt

# =========================================================
# INPUT: matrici RAW corrette
# =========================================================
datasets = {
    "adeno": "/media/user/8Tb/scRNAseq/adeno_wf/adeno/adeno.gene_raw_feature_bc_matrix",
    "sham":  "/media/user/8Tb/scRNAseq/sham_wf/sham/sham.gene_raw_feature_bc_matrix",
}

# =========================================================
# OUTPUT
# =========================================================
output_dir = "/media/user/8Tb/scRNAseq/seurat_analysis/scrublet_independent_raw_output"
os.makedirs(output_dir, exist_ok=True)

# =========================================================
# PARAMETRI SCRUBLET
# =========================================================
expected_doublet_rate = 0.07
min_counts_scrublet = 2
min_cells_scrublet = 3
min_gene_variability_pctl = 85
n_prin_comps = 15

# =========================================================
# PREFILTRO CELLULARE MINIMO (indipendente da DoubletFinder)
# Questo serve solo a togliere barcode troppo deboli / non-cellule.
# NON mettiamo cutoff alto su nFeature/nCount, per non eliminare doublets.
# =========================================================
min_features_cell = 200
min_counts_cell = 250
max_percent_mt = 10.0

# Se vuoi comunque salvare anche una chiamata manuale, lasciala qui.
# Non sarà la chiamata "principale".
manual_threshold = 0.25

# =========================================================
# FUNZIONI
# =========================================================
def read_mtx_gz(path):
    with gzip.open(path, "rb") as f:
        return scipy.io.mmread(f).tocsc()

def read_lines_gz(path):
    with gzip.open(path, "rt") as f:
        return [line.rstrip("\n") for line in f]

def load_10x_raw_matrix(input_dir):
    matrix_path = os.path.join(input_dir, "matrix.mtx.gz")
    features_path = os.path.join(input_dir, "features.tsv.gz")
    barcodes_path = os.path.join(input_dir, "barcodes.tsv.gz")

    missing = [p for p in [matrix_path, features_path, barcodes_path] if not os.path.exists(p)]
    if missing:
        raise FileNotFoundError(f"File mancanti in {input_dir}: {missing}")

    # matrix.mtx = genes x cells
    mat_gc = read_mtx_gz(matrix_path)

    features_raw = read_lines_gz(features_path)
    barcodes = np.array(read_lines_gz(barcodes_path))

    features_split = [x.split("\t") for x in features_raw]
    gene_ids = np.array([x[0] if len(x) > 0 else "" for x in features_split])
    gene_names = np.array([x[1] if len(x) > 1 else x[0] for x in features_split])

    return mat_gc, gene_ids, gene_names, barcodes

def compute_nfeatures(mat_gc):
    return np.asarray((mat_gc > 0).sum(axis=0)).ravel()

def compute_ncounts(mat_gc):
    return np.asarray(mat_gc.sum(axis=0)).ravel()

def compute_percent_mt(mat_gc, gene_names):
    mt_mask = np.array([g.lower().startswith("mt-") for g in gene_names])
    total_counts = np.asarray(mat_gc.sum(axis=0)).ravel()
    mt_counts = np.asarray(mat_gc[mt_mask, :].sum(axis=0)).ravel()

    percent_mt = np.zeros_like(total_counts, dtype=float)
    nz = total_counts > 0
    percent_mt[nz] = 100.0 * mt_counts[nz] / total_counts[nz]
    return percent_mt

def save_histogram(scrub, out_png):
    plt.figure()
    scrub.plot_histogram()
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

# =========================================================
# MAIN
# =========================================================
for sample_name, input_dir in datasets.items():
    print("\n====================================================")
    print(f"Processing sample: {sample_name}")
    print(f"Input dir: {input_dir}")
    print("====================================================")

    try:
        mat_gc, gene_ids, gene_names, barcodes = load_10x_raw_matrix(input_dir)

        vals = mat_gc.data
        all_integer = np.all(np.abs(vals - np.round(vals)) < 1e-8)

        print(f"Raw matrix shape (genes x cells): {mat_gc.shape[0]} x {mat_gc.shape[1]}")
        print(f"All non-zero values are integers: {all_integer}")

        if not all_integer:
            raise ValueError(
                f"{sample_name}: matrice non intera. Scrublet va eseguito su raw counts interi."
            )

        # -------------------------------------------------
        # QC minimo indipendente per selezionare vere cellule candidate
        # -------------------------------------------------
        nFeature_RNA = compute_nfeatures(mat_gc)
        nCount_RNA = compute_ncounts(mat_gc)
        percent_mt = compute_percent_mt(mat_gc, gene_names)

        qc_df = pd.DataFrame({
            "Barcode": barcodes,
            "nFeature_RNA": nFeature_RNA,
            "nCount_RNA": nCount_RNA,
            "percent_mt": percent_mt
        })

        pass_qc = (
            (qc_df["nFeature_RNA"] >= min_features_cell) &
            (qc_df["nCount_RNA"] >= min_counts_cell) &
            (qc_df["percent_mt"] <= max_percent_mt)
        )

        print(f"Cells before prefilter: {len(qc_df)}")
        print(f"Cells after prefilter:  {int(pass_qc.sum())}")
        print(f"Cells removed:          {int((~pass_qc).sum())}")

        keep_idx = np.where(pass_qc.values)[0]
        mat_gc_filt = mat_gc[:, keep_idx]
        qc_df_filt = qc_df.iloc[keep_idx].reset_index(drop=True)

        # Scrublet vuole cells x genes
        counts_matrix = mat_gc_filt.T.tocsc()

        print(f"Filtered matrix shape (cells x genes): {counts_matrix.shape[0]} x {counts_matrix.shape[1]}")

        # -------------------------------------------------
        # Scrublet
        # -------------------------------------------------
        scrub = scr.Scrublet(
            counts_matrix,
            expected_doublet_rate=expected_doublet_rate
        )

        doublet_scores, predicted_doublets_auto = scrub.scrub_doublets(
            min_counts=min_counts_scrublet,
            min_cells=min_cells_scrublet,
            min_gene_variability_pctl=min_gene_variability_pctl,
            n_prin_comps=n_prin_comps
        )

        auto_threshold = getattr(scrub, "threshold_", np.nan)
        predicted_doublets_manual = doublet_scores > manual_threshold

        # -------------------------------------------------
        # Salvataggi
        # -------------------------------------------------
        results_df = qc_df_filt.copy()
        results_df["DoubletScore"] = doublet_scores
        results_df["PredictedDoublet_Auto"] = predicted_doublets_auto
        results_df["PredictedDoublet_Manual"] = predicted_doublets_manual
        results_df["AutoThreshold"] = auto_threshold
        results_df["ManualThreshold"] = manual_threshold
        results_df["Sample"] = sample_name

        results_csv = os.path.join(output_dir, f"{sample_name}_scrublet_independent_results.csv")
        results_df.to_csv(results_csv, index=False)

        summary_df = pd.DataFrame([{
            "Sample": sample_name,
            "InputDir": input_dir,
            "CellsBeforePrefilter": int(len(qc_df)),
            "CellsAfterPrefilter": int(len(qc_df_filt)),
            "ExpectedDoubletRate": expected_doublet_rate,
            "AutoThreshold": auto_threshold,
            "ManualThreshold": manual_threshold,
            "AutoDoublets": int(np.sum(predicted_doublets_auto)),
            "AutoSinglets": int(np.sum(~predicted_doublets_auto)),
            "ManualDoublets": int(np.sum(predicted_doublets_manual)),
            "ManualSinglets": int(np.sum(~predicted_doublets_manual)),
            "MeanDoubletScore": float(np.mean(doublet_scores)),
            "MedianDoubletScore": float(np.median(doublet_scores))
        }])

        summary_csv = os.path.join(output_dir, f"{sample_name}_scrublet_independent_summary.csv")
        summary_df.to_csv(summary_csv, index=False)

        histogram_path = os.path.join(output_dir, f"{sample_name}_scrublet_independent_histogram.png")
        save_histogram(scrub, histogram_path)

        print("\nDone.")
        print(f"Saved: {results_csv}")
        print(f"Saved: {summary_csv}")
        print(f"Saved: {histogram_path}")

    except Exception as e:
        print(f"Error processing {sample_name}: {e}")

print("\nProcessing complete.")
