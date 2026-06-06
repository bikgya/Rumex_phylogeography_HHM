#!/bin/bash
# ==============================================================
# Script 03a: ANGSD Hard Genotype Calling for GEA (RDA + LFMM2)
# Project: rumex-landscape-phylogenomics
#
# Calls hard genotypes (0/1/2) from the filtered BAM set for use
# in genotype-environment association (GEA) analyses. Hard calls
# are required by both RDA and LFMM2; the high posterior cutoff
# (0.95) yields conservative, high-confidence genotypes suitable
# for outlier-detection methods.
#
# Statistics computed
#   - Per-individual genotype calls (0 = homozygous major,
#     1 = heterozygous, 2 = homozygous minor, -1 = missing)
#   - Minor allele frequencies (.mafs.gz)
#   - Posterior genotype probabilities used for hard calls
#
# Input
#   - BAM list of GEA-filtered individuals (one path per line)
#   - Reference genome FASTA (same as Script 01b)
#
# Output (under OUTDIR)
#   gea_snps.geno.gz    hard genotype calls (rows = SNPs)
#   gea_snps.mafs.gz    SNP coordinates + MAFs
#   gea_snps.arg        ANGSD argument log
#
# Tools: ANGSD v0.935
# ==============================================================

set -euo pipefail

# ----- USER CONFIGURATION -----
BAMLIST="data/sample_info/gea_bamlist.txt"          # BAMs of GEA-filtered individuals
REF="data/raw/reference/genome.fasta"
OUTDIR="results/gea/genotypes"
OUTPREFIX="${OUTDIR}/gea_snps"

# ANGSD filtering
GL=2                # GATK genotype likelihood model
MIN_MAP_QUALITY=30  # NOTE: stricter than Script 01b (was 20) for GEA precision
MIN_BASE_QUALITY=20
MIN_MAF=0.05
SNP_PVAL=1e-6
POST_CUTOFF=0.95    # posterior threshold for hard genotype call
THREADS=80

# Min individuals = ~80% of sample size (set explicitly per dataset)
MIN_IND=121
# ------------------------------

mkdir -p "${OUTDIR}"
N_SAMPLES=$(wc -l < "${BAMLIST}")

echo "=============================================="
echo "  ANGSD hard genotype calling for GEA"
echo "=============================================="
echo "  Samples:         ${N_SAMPLES}"
echo "  Min individuals: ${MIN_IND}"
echo "  Min MAF:         ${MIN_MAF}"
echo "  Post cutoff:     ${POST_CUTOFF}"
echo "  Threads:         ${THREADS}"
echo "=============================================="

# NOTE: do NOT put comments after backslashes in line-continued commands -
# the comment breaks the continuation. Keep documentation in the block above.
angsd \
    -b "${BAMLIST}" \
    -ref "${REF}" \
    -out "${OUTPREFIX}" \
    -GL ${GL} \
    -doMajorMinor 1 \
    -doMaf 1 \
    -SNP_pval ${SNP_PVAL} \
    -doGeno 2 \
    -doPost 1 \
    -postCutoff ${POST_CUTOFF} \
    -minMapQ ${MIN_MAP_QUALITY} \
    -minQ ${MIN_BASE_QUALITY} \
    -minMaf ${MIN_MAF} \
    -minInd ${MIN_IND} \
    -doCounts 1 \
    -remove_bads 1 \
    -only_proper_pairs 1 \
    -nThreads ${THREADS}

echo "Done. Output files:"
ls -lh "${OUTPREFIX}".* 2>/dev/null || true
