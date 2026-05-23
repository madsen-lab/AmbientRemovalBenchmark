#!/bin/bash
# run_CellClear.sh
# Runs CellClear correct_expression on all datasets.
# Run from the project root: bash scripts/bash/run_CellClear.sh
#
# Requires the cellclear conda environment (see Setup section).
# Smart-seq2 was run interactively due to donor-specific parameter tuning;
# see the note in the Smart-seq2 section below.

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
source "$PROJECT_ROOT/scripts/bash/record_command_times.sh"

RAW="$PROJECT_ROOT/data/raw"
FILT="$PROJECT_ROOT/data/processed"
OUT="$PROJECT_ROOT/data/processed"


# ===========================================================================
# Setup: create conda environment and install CellClear (run once)
# ===========================================================================

# conda create -n cellclear
# conda activate cellclear
# pip install CellClear


# ===========================================================================
# Create output directories
# ===========================================================================

mkdir -p "$OUT/low_complexity/CellClear"

mkdir -p "$OUT/medium_complexity/CellClear"
for ((i=483; i<=516; i++)); do
    mkdir -p "$OUT/medium_complexity/CellClear/SRR10751${i}"
done

mkdir -p "$OUT/high_complexity/CellClear"
for ((i=19; i<=34; i++)); do
    mkdir -p "$OUT/high_complexity/CellClear/SRR199629${i}"
done

mkdir -p "$OUT/smart-seq2/CellClear"
for i in HP1502401 HP1504101T2D HP1504901 HP1506401 HP1507101 HP1508501T2D HP1509101 HP1525301T2D; do
    mkdir -p "$OUT/smart-seq2/CellClear/$i"
done

mkdir -p "$OUT/genotype/CellClear"
for i in rep1 rep2 rep3 nuc2 nuc3; do
    mkdir -p "$OUT/genotype/CellClear/$i"
done

# Synthetic — with ambient (two parameter sets)
mkdir -p "$OUT/synthetic/CellClear/normal/default_parameters"
mkdir -p "$OUT/synthetic/CellClear/normal/umi_parameters"
for i in sample{1..6}; do
    mkdir -p "$OUT/synthetic/CellClear/normal/default_parameters/$i"
    mkdir -p "$OUT/synthetic/CellClear/normal/umi_parameters/$i"
done

# Synthetic — without ambient (UMI-parameter run used in the final analysis)
mkdir -p "$OUT/synthetic/CellClear/wo_ambient/umi_parameters"
for i in sample{1..6}; do
    mkdir -p "$OUT/synthetic/CellClear/wo_ambient/umi_parameters/$i"
done

mkdir -p "$OUT/memory_use/cellclear"


# ===========================================================================
# Low complexity
# ===========================================================================

tic "CellClear: low complexity"
CellClear correct_expression \
    --filtered_mtx_path "$FILT/low_complexity/emptyDrops_filtered/h5-files/low_complexity_post_emptyDrops.h5" \
    --raw_mtx_path      "$RAW/low_complexity/raw_feature_bc_matrix/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5" \
    --prefix            low_complexity_post_CellClear \
    --output            "$OUT/low_complexity/CellClear"
toc "CellClear: low complexity"


# ===========================================================================
# Medium complexity — human and mouse samples
# ===========================================================================

for ((i=483; i<=516; i++)); do
    SRR="SRR10751${i}"
    tic "CellClear: medium complexity $SRR"
    CellClear correct_expression \
        --filtered_mtx_path "$FILT/medium_complexity/emptyDrops_filtered/h5-files/${SRR}_post_emptyDrops.h5" \
        --raw_mtx_path      "$RAW/medium_complexity/h5-files/${SRR}.h5" \
        --prefix            "${SRR}_post_CellClear" \
        --output            "$OUT/medium_complexity/CellClear/$SRR"
    toc "CellClear: medium complexity $SRR"
done


# ===========================================================================
# High complexity
# ===========================================================================

for ((i=19; i<=34; i++)); do
    SRR="SRR199629${i}"
    tic "CellClear: high complexity $SRR"
    CellClear correct_expression \
        --filtered_mtx_path "$FILT/high_complexity/emptyDrops_filtered/h5-files/${SRR}_post_emptyDrops.h5" \
        --raw_mtx_path      "$RAW/high_complexity/h5-files/${SRR}.h5" \
        --prefix            "${SRR}_post_CellClear" \
        --output            "$OUT/high_complexity/CellClear/$SRR"
    toc "CellClear: high complexity $SRR"
done


# ===========================================================================
# Genotype
# ===========================================================================

for i in rep1 rep2 rep3 nuc2 nuc3; do
    tic "CellClear: genotype $i"
    CellClear correct_expression \
        --filtered_mtx_path "$FILT/genotype/emptyDrops_filtered/h5-files/${i}_post_emptyDrops.h5" \
        --raw_mtx_path      "$RAW/genotype/h5-files/${i}_raw_feature_bc_matrix.h5" \
        --prefix            "${i}_post_CellClear" \
        --output            "$OUT/genotype/CellClear/$i"
    toc "CellClear: genotype $i"
done


# ===========================================================================
# Smart-seq2
# CellClear was run interactively for Smart-seq2 donors due to
# donor-specific tuning of --min_environ_umi and --max_environ_umi.
# The AZ donor was excluded as it lacked sufficient empty droplets.
# Commands for each donor are documented in manual_CellClear.sh.
# ===========================================================================


# ===========================================================================
# Synthetic — with ambient RNA
# ===========================================================================

# Default parameters (exploratory; UMI-parameter run used in the paper)
for i in sample{1..6}; do
    tic "CellClear: synthetic (default) $i"
    CellClear correct_expression \
        --filtered_mtx_path "$FILT/synthetic/filtered/h5-files/normal/${i}_filtered.h5" \
        --raw_mtx_path      "$RAW/synthetic/h5-files/normal/${i}.h5" \
        --prefix            "${i}_post_CellClear" \
        --output            "$OUT/synthetic/CellClear/normal/default_parameters/$i"
    toc "CellClear: synthetic (default) $i"
done

# With explicit UMI boundaries matching emptyDrops cutoffs (used in the paper)
for i in sample{1..6}; do
    tic "CellClear: synthetic (umi-parameters) $i"
    CellClear correct_expression \
        --filtered_mtx_path "$FILT/synthetic/filtered/h5-files/normal/${i}_filtered.h5" \
        --raw_mtx_path      "$RAW/synthetic/h5-files/normal/${i}.h5" \
        --min_environ_umi 0 \
        --max_environ_umi 1000 \
        --prefix            "${i}_post_CellClear_with_umi_parameters" \
        --output            "$OUT/synthetic/CellClear/normal/umi_parameters/$i"
    toc "CellClear: synthetic (umi-parameters) $i"
done


# ===========================================================================
# Synthetic — without ambient RNA
# ===========================================================================

# With explicit UMI boundaries (used in the paper)
for i in sample{1..6}; do
    tic "CellClear: synthetic wo_ambient (umi-parameters) $i"
    CellClear correct_expression \
        --filtered_mtx_path "$FILT/synthetic/filtered/h5-files/wo_ambient/${i}_wo_ambient_filtered.h5" \
        --raw_mtx_path      "$RAW/synthetic/h5-files/wo_ambient/${i}_wo_ambient.h5" \
        --min_environ_umi 0 \
        --max_environ_umi 500 \
        --prefix            "${i}_wo_ambient_post_CellClear_with_umi_parameters" \
        --output            "$OUT/synthetic/CellClear/wo_ambient/umi_parameters/$i"
    toc "CellClear: synthetic wo_ambient (umi-parameters) $i"
done


# ===========================================================================
# Memory benchmarking (low complexity only)
# Runs a second time with /usr/bin/time -v to record peak RSS.
# Requires: sudo apt install time
# ===========================================================================

sudo apt-get install -y time

/usr/bin/time -v CellClear correct_expression \
    --filtered_mtx_path "$FILT/low_complexity/emptyDrops_filtered/h5-files/low_complexity_post_emptyDrops.h5" \
    --raw_mtx_path      "$RAW/low_complexity/raw_feature_bc_matrix/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5" \
    --prefix            low_complexity_post_CellClear \
    --output            "$OUT/memory_use/cellclear" \
    2> "$PROJECT_ROOT/results/memory_use/python_based/CellClear_lc_mem.txt"
