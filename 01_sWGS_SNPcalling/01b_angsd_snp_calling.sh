#!/bin/bash
# ==============================================================
# Script 01b: Nuclear SNP Calling via ANGSD
# Project: rumex-landscape-phylogenomics
# Author:  Bikram Jnawali
#
# Description:
#   Identifies nuclear SNPs across all samples using ANGSD with
#   the GATK genotype likelihood model (-GL 2). Produces:
#     - Minor allele frequency estimates (.mafs.gz)
#     - Posterior genotype probabilities (.geno.gz)
#     - Allele counts and depth summaries
#   Output forms the basis for all downstream population genetic
#   analyses (diversity, structure, GEA, demography, gene flow).
#
# Input:
#   - BAM list file: one path per line, pointing to sorted/indexed
#     BAM files produced by Script 01a
#   - Reference genome FASTA
#
# Output:
#   - <OUT_DIR>/<PREFIX>.mafs.gz
#   - <OUT_DIR>/<PREFIX>.geno.gz
#   - <OUT_DIR>/<PREFIX>.depthSample
#   - <OUT_DIR>/<PREFIX>.depthGlobal
#   - <OUT_DIR>/<PREFIX>.log
#
# Dependencies:
#   - ANGSD v0.935
#   (install via: conda activate rumex-snp)
#
# Usage:
#   1. Edit the USER CONFIGURATION block below
#   2. Generate a BAM list:
#        ls results/snp_calling/bam/*.sorted.bam > data/sample_info/bam_list.txt
#   3. Run: bash scripts/01_snp_calling/01b_angsd_snp_calling.sh
# ==============================================================

set -euo pipefail

# ==============================================================
# USER CONFIGURATION — edit these paths for your system
# ==============================================================

# BAM list — one absolute or relative path per line (output of Script 01a)
BAM_LIST="data/sample_info/bam_list.txt"

# Reference genome (same as used for mapping)
REF="data/raw/reference/genome.fasta"

# Output directory and file prefix
OUT_DIR="results/snp_calling/angsd"
PREFIX="rumex_snps_final"

# Threads
THREADS=20

# ----- ANGSD filtering & output parameters -----
GL=2                # Genotype likelihood model (2 = GATK)
DO_MAF=2            # MAF estimation method
DO_GENO=32          # Output binary posterior probabilities
DO_POST=1           # Posterior genotype probability
DO_MAJOR_MINOR=1    # Infer major/minor from genotype likelihoods
DO_COUNTS=1         # Output base counts
DO_DEPTH=1          # Output depth summaries

MIN_MAP_QUALITY=25
MIN_BASE_QUALITY=20
MIN_IND_PROP=0.5    # Minimum fraction of individuals with data (-minInd)
MIN_MAF=0.01        # Minimum minor allele frequency
SNP_PVAL=1e-6       # SNP-calling p-value threshold
MAX_DEPTH=1000      # Maximum total depth (filters duplicates / paralogs)

# ==============================================================
# DO NOT EDIT BELOW UNLESS YOU KNOW WHAT YOU ARE CHANGING
# ==============================================================

# Verify inputs
[ -f "${BAM_LIST}" ] || { echo "ERROR: BAM list not found: ${BAM_LIST}"; exit 1; }
[ -f "${REF}" ]      || { echo "ERROR: Reference not found: ${REF}"; exit 1; }

# Count samples and derive absolute minInd from proportion
N_SAMPLES=$(wc -l < "${BAM_LIST}")
MIN_IND=$(awk -v n="${N_SAMPLES}" -v p="${MIN_IND_PROP}" 'BEGIN { printf "%d", n * p }')

# Create output directory
mkdir -p "${OUT_DIR}"

# Log start
echo "=============================================="
echo "  ANGSD SNP calling"
echo "=============================================="
echo "  BAM list:        ${BAM_LIST}"
echo "  Samples:         ${N_SAMPLES}"
echo "  Reference:       ${REF}"
echo "  Output:          ${OUT_DIR}/${PREFIX}"
echo "  GL model:        ${GL} (GATK)"
echo "  Min mapping Q:   ${MIN_MAP_QUALITY}"
echo "  Min base Q:      ${MIN_BASE_QUALITY}"
echo "  Min individuals: ${MIN_IND} (${MIN_IND_PROP} of ${N_SAMPLES})"
echo "  Min MAF:         ${MIN_MAF}"
echo "  SNP p-value:     ${SNP_PVAL}"
echo "  Max depth:       ${MAX_DEPTH}"
echo "  Threads:         ${THREADS}"
echo "  Started:         $(date)"
echo "=============================================="

# ----- Run ANGSD -----
angsd \
    -bam "${BAM_LIST}" \
    -ref "${REF}" \
    -out "${OUT_DIR}/${PREFIX}" \
    -GL ${GL} \
    -doMaf ${DO_MAF} \
    -doGeno ${DO_GENO} \
    -doPost ${DO_POST} \
    -doMajorMinor ${DO_MAJOR_MINOR} \
    -doCounts ${DO_COUNTS} \
    -doDepth ${DO_DEPTH} \
    -minMapQ ${MIN_MAP_QUALITY} \
    -minQ ${MIN_BASE_QUALITY} \
    -minInd ${MIN_IND} \
    -minMaf ${MIN_MAF} \
    -SNP_pval ${SNP_PVAL} \
    -maxDepth ${MAX_DEPTH} \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -nThreads ${THREADS} \
    > "${OUT_DIR}/${PREFIX}.log" 2>&1

echo "=============================================="
echo "  ANGSD finished: $(date)"
echo "  Log:  ${OUT_DIR}/${PREFIX}.log"
echo "  Output files:"
ls -lh "${OUT_DIR}/${PREFIX}".* 2>/dev/null || true
echo "=============================================="
