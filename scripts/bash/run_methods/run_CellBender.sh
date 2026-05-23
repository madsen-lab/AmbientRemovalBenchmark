#!/bin/bash
# run_CellBender.sh
# Runs CellBender remove-background on all datasets.
# Run from the project root: bash scripts/bash/run_CellBender.sh
#
# Requires the cellbender conda environment (see Setup section).
# GPU (CUDA) is used for all runs.

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
source "$PROJECT_ROOT/scripts/bash/record_command_times.sh"

RAW="$PROJECT_ROOT/data/raw"
FILT="$PROJECT_ROOT/data/processed"
OUT="$PROJECT_ROOT/data/processed"


# ===========================================================================
# Setup: create conda environment and install CellBender (run once)
# ===========================================================================

# conda create -n cellbender python=3.7
# conda activate cellbender
# pip install cellbender


# ===========================================================================
# Create output directories
# ===========================================================================

mkdir -p "$OUT/low_complexity/CellBender"

mkdir -p "$OUT/medium_complexity/CellBender"
for ((i=483; i<=516; i++)); do
    mkdir -p "$OUT/medium_complexity/CellBender/SRR10751${i}"
done

mkdir -p "$OUT/high_complexity/CellBender"
for ((i=19; i<=34; i++)); do
    mkdir -p "$OUT/high_complexity/CellBender/SRR199629${i}"
done

mkdir -p "$OUT/smart-seq2/CellBender"
for i in AZ HP1502401 HP1504101T2D HP1504901 HP1506401 HP1507101 HP1508501T2D HP1509101 HP1525301T2D; do
    mkdir -p "$OUT/smart-seq2/CellBender/$i"
done

mkdir -p "$OUT/genotype/CellBender"
for i in rep1 rep2 rep3 nuc2 nuc3; do
    mkdir -p "$OUT/genotype/CellBender/$i"
done

mkdir -p "$OUT/synthetic/CellBender/normal/default_parameters"
mkdir -p "$OUT/synthetic/CellBender/normal/expected_cells_parameter"
mkdir -p "$OUT/synthetic/CellBender/wo_ambient/default_parameters"
mkdir -p "$OUT/synthetic/CellBender/wo_ambient/expected_cells_parameter"
for i in sample{1..6}; do
    mkdir -p "$OUT/synthetic/CellBender/normal/default_parameters/$i"
    mkdir -p "$OUT/synthetic/CellBender/normal/expected_cells_parameter/$i"
    mkdir -p "$OUT/synthetic/CellBender/wo_ambient/default_parameters/$i"
    mkdir -p "$OUT/synthetic/CellBender/wo_ambient/expected_cells_parameter/$i"
done

mkdir -p "$OUT/memory_use/cellbender"


# ===========================================================================
# Low complexity
# ===========================================================================

tic "CellBender: low complexity"
cellbender remove-background \
    --cuda \
    --input  "$RAW/low_complexity/raw_feature_bc_matrix/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5" \
    --output "$OUT/low_complexity/CellBender/low_complexity_post_CellBender.h5"
toc "CellBender: low complexity"


# ===========================================================================
# Medium complexity — human and mouse samples (SRR10751483–SRR10751516)
# Spike-in samples (SRR10751517–SRR10751518) are excluded; they are combined
# and filtered in R (see 04_preprocessing_and_quality_control.qmd).
# ===========================================================================

for ((i=483; i<=516; i++)); do
    SRR="SRR10751${i}"
    tic "CellBender: medium complexity $SRR"
    cellbender remove-background \
        --cuda \
        --input  "$RAW/medium_complexity/h5-files/${SRR}.h5" \
        --output "$OUT/medium_complexity/CellBender/$SRR/${SRR}_post_CellBender.h5"
    toc "CellBender: medium complexity $SRR"
done


# ===========================================================================
# High complexity
# ===========================================================================

for ((i=19; i<=34; i++)); do
    SRR="SRR199629${i}"
    tic "CellBender: high complexity $SRR"
    cellbender remove-background \
        --cuda \
        --input  "$RAW/high_complexity/h5-files/${SRR}.h5" \
        --output "$OUT/high_complexity/CellBender/$SRR/${SRR}_post_CellBender.h5"
    toc "CellBender: high complexity $SRR"
done


# ===========================================================================
# Genotype
# ===========================================================================

for i in rep1 rep2 rep3 nuc2 nuc3; do
    tic "CellBender: genotype $i"
    cellbender remove-background \
        --cuda \
        --input  "$RAW/genotype/h5-files/${i}_raw_feature_bc_matrix.h5" \
        --output "$OUT/genotype/CellBender/$i/${i}_post_CellBender.h5"
    toc "CellBender: genotype $i"
done


# ===========================================================================
# Smart-seq2
# ===========================================================================

for i in AZ HP1502401 HP1504101T2D HP1504901 HP1506401 HP1507101 HP1508501T2D HP1509101 HP1525301T2D; do
    tic "CellBender: smart-seq2 $i"
    cellbender remove-background \
        --cuda \
        --input  "$RAW/smart-seq2/h5-files/${i}.h5" \
        --output "$OUT/smart-seq2/CellBender/$i/${i}_post_CellBender.h5"
    toc "CellBender: smart-seq2 $i"
done


# ===========================================================================
# Synthetic — with ambient RNA
# ===========================================================================

# Default parameters
for i in sample{1..6}; do
    tic "CellBender: synthetic (default) $i"
    cellbender remove-background \
        --cuda \
        --input  "$RAW/synthetic/h5-files/normal/${i}.h5" \
        --output "$OUT/synthetic/CellBender/normal/default_parameters/$i/${i}_post_CellBender.h5"
    toc "CellBender: synthetic (default) $i"
done

# With --expected-cells parameter (used in the final analysis)
# Expected cells = 3240 based on the simulation parameters
for i in sample{1..6}; do
    tic "CellBender: synthetic (expected-cells) $i"
    cellbender remove-background \
        --cuda \
        --expected-cells 3240 \
        --input  "$RAW/synthetic/h5-files/normal/${i}.h5" \
        --output "$OUT/synthetic/CellBender/normal/expected_cells_parameter/$i/${i}_post_CellBender_with_expected_cells.h5"
    toc "CellBender: synthetic (expected-cells) $i"
done


# ===========================================================================
# Synthetic — without ambient RNA
# ===========================================================================

# Default parameters
for i in sample{1..6}; do
    tic "CellBender: synthetic wo_ambient (default) $i"
    cellbender remove-background \
        --cuda \
        --input  "$RAW/synthetic/h5-files/wo_ambient/${i}_wo_ambient.h5" \
        --output "$OUT/synthetic/CellBender/wo_ambient/default_parameters/$i/${i}_wo_ambient_post_CellBender.h5"
    toc "CellBender: synthetic wo_ambient (default) $i"
done

# With --expected-cells parameter (used in the final analysis)
for i in sample{1..6}; do
    tic "CellBender: synthetic wo_ambient (expected-cells) $i"
    cellbender remove-background \
        --cuda \
        --expected-cells 3240 \
        --input  "$RAW/synthetic/h5-files/wo_ambient/${i}_wo_ambient.h5" \
        --output "$OUT/synthetic/CellBender/wo_ambient/expected_cells_parameter/$i/${i}_wo_ambient_post_CellBender_with_expected_cells.h5"
    toc "CellBender: synthetic wo_ambient (expected-cells) $i"
done


# ===========================================================================
# Memory benchmarking (low complexity only)
# Runs a second time with /usr/bin/time -v to record peak RSS.
# Requires: sudo apt install time
# ===========================================================================

sudo apt-get install -y time

/usr/bin/time -v cellbender remove-background \
    --cuda \
    --input  "$RAW/low_complexity/raw_feature_bc_matrix/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5" \
    --output "$OUT/memory_use/cellbender/low_complexity_post_CellBender.h5" \
    2> "$PROJECT_ROOT/results/memory_use/python_based/CellBender_lc_mem.txt"
