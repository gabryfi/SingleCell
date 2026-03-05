#!/usr/bin/env bash
set -euo pipefail

#############################################
# wf-single-cell wrapper (10x 3' v3 classic, Nanopore)
# Target: mouse GRCm39 (reference locale 10x-style)
#
# Requisiti:
#   - nextflow
#   - docker oppure apptainer/singularity (scegli con PROFILE)
#############################################

###############
# VARIABILI DA CONFIGURARE
###############

# Input: lascia come placeholder, come vuoi tu
# Può essere:
#  (i)   singolo FASTQ(.gz)
#  (ii)  directory con FASTQ(.gz)
#  (iii) directory con 1 livello di sottodirectory (una per libreria)
READS_INPUT="/path/to/fastq_or_folder_or_library_folders"

# Se READS_INPUT è (i) o (ii)
SAMPLE_NAME="lib_01"

# Se READS_INPUT è (iii) (4 librerie in sottocartelle), puoi (opzionale) passare un sample sheet
# vedi NOTE più sotto
SAMPLE_SHEET=""  # es: "/path/to/sample_sheet.csv"

# 10x kit
KIT="3prime:v3"

# Cellule attese (per libreria)
EXPECTED_CELLS=10000
ESTIMATE_CELL_COUNT="true"     # default true; lascia così

# Reference locale (obbligatoria per mouse)
# Deve essere formato 10x refdata:
#   refdata/fasta/genome.fa(.gz)
#   refdata/genes/genes.gtf(.gz)
REF_GENOME_DIR="/path/to/10x_refdata_mouse_GRCm39"

# Output
OUTDIR="./wf-single-cell_out"
WORKDIR="./wf-single-cell_work"
STORE_DIR="./wf-single-cell_store"

# Runtime
PROFILE="standard"   # tipico: standard (docker) oppure un profilo apptainer/singularity
RESUME="true"

# Pipeline behavior
FULL_LENGTH_ONLY="true"   # vuoi full-length only: ok
# min_read_qual: lasciamo default (qualità standard)

# JVM memory per Nextflow (non è la RAM usata dai processi bioinfo; è solo per Nextflow)
NXF_MEM_GB=16
export NXF_OPTS="-Xms2g -Xmx${NXF_MEM_GB}g"

###############
# CONTROLLI
###############
if [[ ! -e "$READS_INPUT" ]]; then
  echo "ERRORE: READS_INPUT non esiste: $READS_INPUT" >&2
  exit 1
fi
if [[ ! -d "$REF_GENOME_DIR" ]]; then
  echo "ERRORE: REF_GENOME_DIR non è una directory: $REF_GENOME_DIR" >&2
  exit 1
fi

mkdir -p "$OUTDIR" "$WORKDIR" "$STORE_DIR"

###############
# COMANDO NEXTFLOW
###############
NF_CMD=( nextflow run epi2me-labs/wf-single-cell
  --fastq "$READS_INPUT"
  --ref_genome_dir "$REF_GENOME_DIR"
  --kit "$KIT"
  --expected_cells "$EXPECTED_CELLS"
  --estimate_cell_count "$ESTIMATE_CELL_COUNT"
  --full_length_only "$FULL_LENGTH_ONLY"
  --out_dir "$OUTDIR"
  --store_dir "$STORE_DIR"
  -work-dir "$WORKDIR"
  -profile "$PROFILE"
)

# sample sheet vs sample name
if [[ -n "$SAMPLE_SHEET" ]]; then
  NF_CMD+=( --sample_sheet "$SAMPLE_SHEET" )
else
  NF_CMD+=( --sample "$SAMPLE_NAME" )
fi

# resume
if [[ "$RESUME" == "true" ]]; then
  NF_CMD+=( -resume )
fi

echo "Eseguo:"
printf ' %q' "${NF_CMD[@]}"
echo
"${NF_CMD[@]}"

echo "OK. Output in: $OUTDIR"
