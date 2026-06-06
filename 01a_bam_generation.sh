#!/bin/bash
# ==============================================================
# Script 01a: BAM File Generation (BWA-MEM Mapping)
# Project: rumex-landscape-phylogenomics
# Author:  Bikram Jnawali
#
# Description:
#   Maps paired-end cleaned reads to the reference genome using
#   BWA-MEM, sorts alignments with SAMtools, and indexes the
#   resulting BAM files. Produces one sorted, indexed BAM file
#   per sample, ready for ANGSD SNP calling (Script 01b).
#
# Input:
#   - Reference genome FASTA (BWA-indexed)
#   - Paired-end cleaned reads:
#       <CLEAN_READS_DIR>/<sample>/<sample>_1.clean.fq.gz
#       <CLEAN_READS_DIR>/<sample>/<sample>_2.clean.fq.gz
#   - Sample name list (one sample ID per line)
#
# Output:
#   - <BAM_DIR>/<sample>.sorted.bam
#   - <BAM_DIR>/<sample>.sorted.bam.bai
#
# Dependencies:
#   - BWA v0.7.17
#   - SAMtools v1.17
#   (install via: conda activate rumex-snp)
#
# Usage:
#   1. Edit the USER CONFIGURATION block below
#   2. Ensure the reference is indexed: bwa index <reference>
#   3. Run: bash scripts/01_snp_calling/01a_bam_generation.sh
# ==============================================================

set -euo pipefail   # Exit on error, undefined variable, or pipe failure

# ==============================================================
# USER CONFIGURATION — edit these paths for your system
# ==============================================================

# Reference genome (must be BWA-indexed: run `bwa index` once beforehand)
REF="data/raw/reference/genome.fasta"

# Directory containing cleaned reads, organized as one folder per sample
CLEAN_READS_DIR="data/raw/clean_reads"

# Output directory for sorted, indexed BAM files
BAM_DIR="results/snp_calling/bam"

# Sample name list (one ID per line, matching folder names in CLEAN_READS_DIR)
SAMPLE_LIST="data/sample_info/sample_names.txt"

# BWA-MEM parameters (match those used in the published analysis)
THREADS=50
SEED_LENGTH=19      # -k: minimum seed length
MISMATCH_PENALTY=3  # -B: mismatch penalty

# ==============================================================
# DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE CHANGING
# ==============================================================

# Verify required inputs
[ -f "${REF}" ]         || { echo "ERROR: Reference not found: ${REF}"; exit 1; }
[ -f "${REF}.bwt" ]     || { echo "ERROR: Reference not BWA-indexed. Run: bwa index ${REF}"; exit 1; }
[ -f "${SAMPLE_LIST}" ] || { echo "ERROR: Sample list not found: ${SAMPLE_LIST}"; exit 1; }
[ -d "${CLEAN_READS_DIR}" ] || { echo "ERROR: Reads directory not found: ${CLEAN_READS_DIR}"; exit 1; }

# Create output directory
mkdir -p "${BAM_DIR}"

# Log start
echo "=============================================="
echo "  BAM generation pipeline"
echo "=============================================="
echo "  Reference:    ${REF}"
echo "  Reads:        ${CLEAN_READS_DIR}"
echo "  Output:       ${BAM_DIR}"
echo "  Samples:      $(wc -l < "${SAMPLE_LIST}")"
echo "  Threads:      ${THREADS}"
echo "  Started:      $(date)"
echo "=============================================="

# Loop through samples
while read -r sample; do
    # Skip blank lines and comments
    [[ -z "${sample}" || "${sample}" =~ ^# ]] && continue

    R1="${CLEAN_READS_DIR}/${sample}/${sample}_1.clean.fq.gz"
    R2="${CLEAN_READS_DIR}/${sample}/${sample}_2.clean.fq.gz"
    OUT_BAM="${BAM_DIR}/${sample}.sorted.bam"

    # Skip if already done
    if [ -f "${OUT_BAM}.bai" ]; then
        echo "[$(date +%H:%M:%S)] Skipping ${sample} (already processed)"
        continue
    fi

    # Verify input reads exist
    [ -f "${R1}" ] || { echo "WARNING: Missing R1 for ${sample}, skipping"; continue; }
    [ -f "${R2}" ] || { echo "WARNING: Missing R2 for ${sample}, skipping"; continue; }

    echo "[$(date +%H:%M:%S)] Mapping ${sample}..."

    # Map, sort, and pipe to BAM
    bwa mem -t "${THREADS}" \
            -k "${SEED_LENGTH}" \
            -B "${MISMATCH_PENALTY}" \
            "${REF}" \
            "${R1}" "${R2}" \
        | samtools sort -@ "${THREADS}" -o "${OUT_BAM}" -

    # Index the sorted BAM
    samtools index "${OUT_BAM}"

    echo "[$(date +%H:%M:%S)] Done: ${OUT_BAM}"

done < "${SAMPLE_LIST}"

echo "=============================================="
echo "  All samples processed."
echo "  Finished: $(date)"
echo "=============================================="
