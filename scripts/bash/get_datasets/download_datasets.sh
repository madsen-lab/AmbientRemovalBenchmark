#!/bin/bash
# download_datasets.sh
# Downloads and aligns all raw datasets used in the study.
# Run from the project root: bash scripts/bash/download_datasets.sh
#
# Requires a conda environment with sra-tools, STAR, fastqc, and cutadapt.
# See the "Setup" section below to create it.

PROJECT_ROOT="${PROJECT_ROOT:-$(pwd)}"

# ===========================================================================
# Setup: create conda environment (run once)
# ===========================================================================

# conda create -n sra-and-star
# conda activate sra-and-star
# conda install bioconda::sra-tools
# conda install bioconda::star
# conda install bioconda::fastqc
# conda install bioconda::cutadapt


# ===========================================================================
# Low complexity
# 10x Genomics 20k HGMM barnyard dataset — pre-aligned H5 file
# ===========================================================================

mkdir -p "$PROJECT_ROOT/data/raw/low_complexity/raw_feature_bc_matrix"
cd "$PROJECT_ROOT/data/raw/low_complexity/raw_feature_bc_matrix"

wget https://cf.10xgenomics.com/samples/cell-exp/6.1.0/20k_hgmm_3p_HT_nextgem_Chromium_X/20k_hgmm_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5


# ===========================================================================
# Medium complexity
# GSE147203 — human islets + Jurkat/32D spike-ins
# SRR10751483–SRR10751518 (19 human + 15 mouse + 2 shared spike-in samples)
#
# Some SRA files contain three reads (technical 8 nt, Read1 76 nt, Read2
# 98 nt); others contain only Read1 and Read2. The if-else block below
# detects which case applies and selects the correct read files.
# ===========================================================================

mkdir -p "$PROJECT_ROOT/data/raw/medium_complexity"
cd "$PROJECT_ROOT/data/raw/medium_complexity"

STAR_GENOME="$PROJECT_ROOT/data/raw/reference_sources/GRCh38-and-mm10-2020-A-10X-Download/star"
WHITELIST="$PROJECT_ROOT/data/raw/reference_sources/STAR-Solo-Whitelists/737K-august-2016.txt"

for ((i=483; i<=518; i++)); do
    SRR="SRR10751${i}"
    OUTDIR="${SRR}_STARsolo_Alignment"

    # Download and convert to FASTQ
    prefetch --max-size 500g --progress "$SRR"
    fasterq-dump "$SRR" --split-files --include-technical -e 30 \
        -O "$PROJECT_ROOT/data/raw/medium_complexity/$SRR"
    rm "$PROJECT_ROOT/data/raw/medium_complexity/$SRR/${SRR}.sra"

    # Select read files: three-read layout uses reads 3 (CB+UMI) and 2 (cDNA);
    # two-read layout uses reads 2 (CB+UMI) and 1 (cDNA)
    if [ -f "$PROJECT_ROOT/data/raw/medium_complexity/$SRR/${SRR}_3.fastq" ]; then
        READ1="$PROJECT_ROOT/data/raw/medium_complexity/$SRR/${SRR}_3.fastq"
        READ2="$PROJECT_ROOT/data/raw/medium_complexity/$SRR/${SRR}_2.fastq"
    else
        READ1="$PROJECT_ROOT/data/raw/medium_complexity/$SRR/${SRR}_2.fastq"
        READ2="$PROJECT_ROOT/data/raw/medium_complexity/$SRR/${SRR}_1.fastq"
    fi

    # Align with STARsolo
    mkdir -p "$OUTDIR"
    STAR \
        --genomeDir "$STAR_GENOME" \
        --readFilesIn "$READ1" "$READ2" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 \
        --soloBarcodeReadLength 0 --soloUMIlen 10 \
        --soloCellFilter None \
        --soloUMIdedup 1MM_CR \
        --soloFeatures Gene GeneFull GeneFull_Ex50pAS GeneFull_ExonOverIntron \
        --outFileNamePrefix "${OUTDIR}/${SRR}_" \
        --runThreadN 60

    # Clean up FASTQ files and large SAM output
    rm -r "$PROJECT_ROOT/data/raw/medium_complexity/$SRR"
    rm "${OUTDIR}/${SRR}_Aligned.out.sam"
done


# ===========================================================================
# High complexity
# GSE207393 — human islets xenografted into murine kidneys
# SRR19962919–SRR19962934 (16 samples)
#
# Based on FastQC quality assessment, 31 nt are trimmed from the 5'-end of
# Read 2 prior to alignment.
# ===========================================================================

mkdir -p "$PROJECT_ROOT/data/raw/high_complexity"
cd "$PROJECT_ROOT/data/raw/high_complexity"

for ((i=19; i<=34; i++)); do
    SRR="SRR199629${i}"
    OUTDIR="${SRR}_STARsolo_Alignment"

    # Download and convert to FASTQ
    prefetch --max-size 500g --progress "$SRR"
    fasterq-dump "$SRR" --split-files --include-technical -e 30 \
        -O "$PROJECT_ROOT/data/raw/high_complexity/$SRR"
    rm "$PROJECT_ROOT/data/raw/high_complexity/$SRR/${SRR}.sra"

    # Trim 31 nt from the 5'-end of Read 2 (CB+UMI read, identified via FastQC)
    cutadapt -u 31 \
        -o "$PROJECT_ROOT/data/raw/high_complexity/$SRR/${SRR}_3_trimmed.fastq" \
           "$PROJECT_ROOT/data/raw/high_complexity/$SRR/${SRR}_3.fastq" \
        --cores 60
    rm "$PROJECT_ROOT/data/raw/high_complexity/$SRR/${SRR}_3.fastq"

    # Align with STARsolo
    mkdir -p "$OUTDIR"
    STAR \
        --genomeDir "$STAR_GENOME" \
        --readFilesIn \
            "$PROJECT_ROOT/data/raw/high_complexity/$SRR/${SRR}_3_trimmed.fastq" \
            "$PROJECT_ROOT/data/raw/high_complexity/$SRR/${SRR}_2.fastq" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$WHITELIST" \
        --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 \
        --soloBarcodeReadLength 0 --soloUMIlen 10 \
        --soloCellFilter None \
        --soloUMIdedup 1MM_CR \
        --soloFeatures Gene GeneFull GeneFull_Ex50pAS GeneFull_ExonOverIntron \
        --outFileNamePrefix "${OUTDIR}/${SRR}_" \
        --runThreadN 60

    # Clean up
    rm -r "$PROJECT_ROOT/data/raw/high_complexity/$SRR"
    rm "${OUTDIR}/${SRR}_Aligned.out.sam"
done


# ===========================================================================
# Smart-seq2 (negative control)
# E-MTAB-5061 — human islets, plate-based scRNA-seq
# No alignment needed; count matrix downloaded directly from ArrayExpress.
# ===========================================================================

mkdir -p "$PROJECT_ROOT/data/raw/smart-seq2"
cd "$PROJECT_ROOT/data/raw/smart-seq2"

wget -O "pancreas_refseq_rpkms_counts_3514sc.txt" \
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/061/E-MTAB-5061/Files/pancreas_refseq_rpkms_counts_3514sc.txt"


# ===========================================================================
# Genotype (strain-mixing)
# GSE218853 — inbred mouse kidney cells (three strains)
# Pre-processed H5 files downloaded directly from GEO.
# ===========================================================================

mkdir -p "$PROJECT_ROOT/data/raw/genotype/h5-files"
cd "$PROJECT_ROOT/data/raw/genotype/h5-files"

# rep1 (GSM6757771)
wget --content-disposition \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6757771&format=file&file=GSM6757771%5Frep1%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# rep2 (GSM6757772)
wget --content-disposition \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6757772&format=file&file=GSM6757772%5Frep2%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# rep3 (GSM6757773)
wget --content-disposition \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6757773&format=file&file=GSM6757773%5Frep3%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# nuc2 (GSM6757774)
wget --content-disposition \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6757774&format=file&file=GSM6757774%5Fnuc2%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"

# nuc3 (GSM6757775)
wget --content-disposition \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM6757775&format=file&file=GSM6757775%5Fnuc3%5Fraw%5Ffeature%5Fbc%5Fmatrix%2Eh5"


# ===========================================================================
# Synthetic dataset
# GSE173634 — three breast cancer cell lines (MCF12A, HCC1500, HS578T)
#
# The cell lines are aligned individually using a 12-nt CB + 8-nt UMI
# chemistry. The command below is shown for HCC1500 (SRR14369363);
# the same command structure applies to MCF12A and HS578T with their
# respective SRR accessions. The simulated datasets themselves are generated
# in R (see 03_download_datasets.qmd).
# ===========================================================================

# Example alignment for HCC1500 (repeat for MCF12A and HS578T)
# STAR \
#     --genomeDir "$PROJECT_ROOT/data/raw/reference_sources/synthetic/STARindex" \
#     --readFilesIn HCC1500_SRR14369363_2.fastq.gz HCC1500_SRR14369363_1.fastq.gz \
#     --soloCBstart 1 --soloCBlen 12 \
#     --soloUMIstart 13 --soloUMIlen 8 \
#     --soloType CB_UMI_Simple \
#     --soloCBwhitelist None \
#     --outSAMtype BAM SortedByCoordinate \
#     --outSAMattributes NH HI nM AS CR UR CB UB \
#     --readFilesCommand zcat \
#     --runThreadN 30 \
#     --soloCellFilter None \
#     --soloUMIdedup 1MM_CR \
#     --soloUMIfiltering MultiGeneUMI_CR \
#     --outFileNamePrefix HCC1500
