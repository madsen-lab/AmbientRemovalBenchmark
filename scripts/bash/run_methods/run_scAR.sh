#!/bin/bash
# run_scAR.sh
# Runs scAR on all datasets using filtered H5 files as input.
# Run from the project root: bash scripts/bash/run_scAR.sh
#
# Requires the scar conda environment (see Setup section).
# GPU (CUDA) is used for all runs.
#
# Note: scAR was also tested with raw .pickle files as an alternative input
# format. Those runs are documented at the bottom of this script but were
# not used in the final analysis — the filtered H5 results were used.

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"
source "$PROJECT_ROOT/scripts/bash/record_command_times.sh"

RAW="$PROJECT_ROOT/data/raw"
FILT="$PROJECT_ROOT/data/processed"
OUT="$PROJECT_ROOT/data/processed"


# ===========================================================================
# Setup: create conda environment and install scAR (run once)
# ===========================================================================

# conda create -n scar
# conda activate scar
# conda install pytorch torchvision -c pytorch
# conda install bioconda::scar


# ===========================================================================
# Create output directories
# ===========================================================================

mkdir -p "$OUT/low_complexity/scAR"

mkdir -p "$OUT/medium_complexity/scAR"
for ((i=483; i<=516; i++)); do
    mkdir -p "$OUT/medium_complexity/scAR/SRR10751${i}"
done

mkdir -p "$OUT/high_complexity/scAR"
for ((i=19; i<=34; i++)); do
    mkdir -p "$OUT/high_complexity/scAR/SRR199629${i}"
done

mkdir -p "$OUT/smart-seq2/scAR"
for i in AZ HP1502401 HP1504101T2D HP1504901 HP1506401 HP1507101 HP1508501T2D HP1509101 HP1525301T2D; do
    mkdir -p "$OUT/smart-seq2/scAR/$i"
done

mkdir -p "$OUT/genotype/scAR"
for i in rep1 rep2 rep3 nuc2 nuc3; do
    mkdir -p "$OUT/genotype/scAR/$i"
done

mkdir -p "$OUT/synthetic/scAR"
for i in sample{1..6}; do
    mkdir -p "$OUT/synthetic/scAR/normal/$i"
    mkdir -p "$OUT/synthetic/scAR/wo_ambient/$i"
done

mkdir -p "$OUT/memory_use/scar"


# ===========================================================================
# Low complexity
# ===========================================================================

tic "scAR: low complexity"
scar "$FILT/low_complexity/emptyDrops_filtered/h5-files/low_complexity_post_emptyDrops.h5" \
    -ft "mRNA" -d "cuda" \
    -o  "$OUT/low_complexity/scAR/low_complexity_post_scAR.h5ad"
toc "scAR: low complexity"


# ===========================================================================
# Medium complexity — human and mouse samples
# ===========================================================================

for ((i=483; i<=516; i++)); do
    SRR="SRR10751${i}"
    tic "scAR: medium complexity $SRR"
    scar "$FILT/medium_complexity/emptyDrops_filtered/h5-files/${SRR}_post_emptyDrops.h5" \
        -ft "mRNA" -d "cuda" \
        -o  "$OUT/medium_complexity/scAR/$SRR/${SRR}_post_scAR.h5ad"
    toc "scAR: medium complexity $SRR"
done


# ===========================================================================
# High complexity
# ===========================================================================

for ((i=19; i<=34; i++)); do
    SRR="SRR199629${i}"
    tic "scAR: high complexity $SRR"
    scar "$FILT/high_complexity/emptyDrops_filtered/h5-files/${SRR}_post_emptyDrops.h5" \
        -ft "mRNA" -d "cuda" \
        -o  "$OUT/high_complexity/scAR/$SRR/${SRR}_post_scAR.h5ad"
    toc "scAR: high complexity $SRR"
done


# ===========================================================================
# Genotype
# ===========================================================================

for i in rep1 rep2 rep3 nuc2 nuc3; do
    tic "scAR: genotype $i"
    scar "$FILT/genotype/emptyDrops_filtered/h5-files/${i}_post_emptyDrops.h5" \
        -ft "mRNA" -d "cuda" \
        -o  "$OUT/genotype/scAR/$i/${i}_post_scAR.h5ad"
    toc "scAR: genotype $i"
done


# ===========================================================================
# Smart-seq2
# ===========================================================================

for i in AZ HP1502401 HP1504101T2D HP1504901 HP1506401 HP1507101 HP1508501T2D HP1509101 HP1525301T2D; do
    tic "scAR: smart-seq2 $i"
    scar "$FILT/smart-seq2/filtered/h5-files/${i}_post_filtering.h5" \
        -ft "mRNA" -d "cuda" \
        -o  "$OUT/smart-seq2/scAR/$i/${i}_post_scAR.h5ad"
    toc "scAR: smart-seq2 $i"
done


# ===========================================================================
# Synthetic — with ambient RNA
# ===========================================================================

for i in sample{1..6}; do
    tic "scAR: synthetic (normal) $i"
    scar "$FILT/synthetic/filtered/h5-files/${i}_filtered.h5" \
        -ft "mRNA" -d "cuda" \
        -o  "$OUT/synthetic/scAR/normal/$i/${i}_post_scAR.h5ad"
    toc "scAR: synthetic (normal) $i"
done


# ===========================================================================
# Synthetic — without ambient RNA
# ===========================================================================

for i in sample{1..6}; do
    tic "scAR: synthetic (wo_ambient) $i"
    scar "$FILT/synthetic/filtered/h5-files/wo_ambient/${i}_wo_ambient_filtered.h5" \
        -ft "mRNA" -d "cuda" \
        -o  "$OUT/synthetic/scAR/wo_ambient/$i/${i}_post_scAR.h5ad"
    toc "scAR: synthetic (wo_ambient) $i"
done


# ===========================================================================
# Memory benchmarking (low complexity only)
# Runs a second time with /usr/bin/time -v to record peak RSS.
# Requires: sudo apt install time
# ===========================================================================

sudo apt-get install -y time

/usr/bin/time -v scar \
    "$FILT/low_complexity/emptyDrops_filtered/h5-files/low_complexity_post_emptyDrops.h5" \
    -ft "mRNA" -d "cuda" \
    -o  "$OUT/memory_use/scar/low_complexity_post_scAR.h5ad" \
    2> "$PROJECT_ROOT/results/memory_use/python_based/scAR_lc_mem.txt"


# ===========================================================================
# Alternative: raw .pickle file input (exploratory, not used in the paper)
#
# scAR supports raw count matrices as .pickle files in addition to filtered
# H5 files. This section documents those runs for completeness. The filtered
# H5 results above were used in the final analysis.
#
# Conversion from H5 to pickle requires convert_h5_to_pickle.py (included
# in scripts/bash/).
# ===========================================================================

# python scripts/bash/convert_h5_to_pickle.py "$RAW/low_complexity/raw_feature_bc_matrix"
# python scripts/bash/convert_h5_to_pickle.py "$RAW/medium_complexity/h5-files"
# python scripts/bash/convert_h5_to_pickle.py "$RAW/high_complexity/h5-files"
# python scripts/bash/convert_h5_to_pickle.py "$RAW/smart-seq2/h5-files"
# python scripts/bash/convert_h5_to_pickle.py "$RAW/genotype/h5-files"
# python scripts/bash/convert_h5_to_pickle.py "$RAW/synthetic/h5-files"

# Low complexity
# scar "$RAW/low_complexity/raw_feature_bc_matrix/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.pickle" \
#     -ft "mRNA" -d "cuda" \
#     -o  "$OUT/low_complexity/scAR/low_complexity_post_scAR_from_raw.h5ad"

# Medium complexity
# for ((i=483; i<=516; i++)); do
#     SRR="SRR10751${i}"
#     scar "$RAW/medium_complexity/h5-files/${SRR}.pickle" \
#         -ft "mRNA" -d "cuda" \
#         -o  "$OUT/medium_complexity/scAR/$SRR/${SRR}_post_scAR_from_raw.h5ad"
# done

# High complexity
# for ((i=19; i<=34; i++)); do
#     SRR="SRR199629${i}"
#     scar "$RAW/high_complexity/h5-files/${SRR}.pickle" \
#         -ft "mRNA" -d "cuda" \
#         -o  "$OUT/high_complexity/scAR/$SRR/${SRR}_post_scAR_from_raw.h5ad"
# done

# Smart-seq2
# for i in AZ HP1502401 HP1504101T2D HP1504901 HP1506401 HP1507101 HP1508501T2D HP1509101 HP1525301T2D; do
#     scar "$RAW/smart-seq2/h5-files/${i}.pickle" \
#         -ft "mRNA" -d "cuda" \
#         -o  "$OUT/smart-seq2/scAR/$i/${i}_post_scAR_from_raw.h5ad"
# done

# Genotype
# for i in rep1 rep2 rep3 nuc2 nuc3; do
#     scar "$RAW/genotype/h5-files/${i}_raw_feature_bc_matrix.pickle" \
#         -ft "mRNA" -d "cuda" \
#         -o  "$OUT/genotype/scAR/$i/${i}_post_scAR_from_raw.h5ad"
# done

# Synthetic (normal)
# for i in sample{1..6}; do
#     scar "$RAW/synthetic/h5-files/${i}.pickle" \
#         -ft "mRNA" -d "cuda" \
#         -o  "$OUT/synthetic/scAR/normal/$i/${i}_post_scAR_from_raw.h5ad"
# done

# Synthetic (wo_ambient)
# for i in sample{1..6}; do
#     scar "$RAW/synthetic/h5-files/wo_ambient/${i}_wo_ambient.pickle" \
#         -ft "mRNA" -d "cuda" \
#         -o  "$OUT/synthetic/scAR/wo_ambient/$i/${i}_post_scAR_from_raw.h5ad"
# done
