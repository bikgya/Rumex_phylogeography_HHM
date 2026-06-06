#!/bin/bash
# ==============================================================
# Script 02: Non-coding SNP Extraction
# Project: rumex-landscape-phylogenomics
#
# Filters the SNP set from Script 01b to retain only non-coding
# SNPs (i.e., those outside annotated CDS/exon regions). The
# resulting TPED is used in neutral-marker analyses (population
# structure, demography, gene flow) where non-coding SNPs are
# preferred to minimize the influence of selection.
#
# Input:
#   - final_analysis.mafs.gz   (from Script 01b)
#   - final_analysis.tped      (from Script 01b, -doPlink 2)
#   - final_analysis.tfam      (from Script 01b)
#   - Reference GFF annotation
#
# Output:
#   - coding_regions.bed              (CDS + exon coordinates)
#   - non_coding_snps.bed             (SNPs outside coding regions)
#   - final_analysis_NONCODING.tped   (filtered TPED)
#   - final_analysis_NONCODING.tfam   (copy of original TFAM)
#
# Tools: bedtools, awk, grep
# ==============================================================

# ----- USER CONFIGURATION -----
WORKDIR="results/snp_calling/angsd"                    # output of Script 01b
GFF="data/raw/reference/R62326.gff"                    # reference annotation
OUTDIR="${WORKDIR}/Noncoding"
PREFIX="final_analysis"
# ------------------------------

mkdir -p "$OUTDIR"
cd "$WORKDIR"

echo "Step 1: Creating coding regions BED file..."
# Use the working pattern that we tested
grep -E "CDS|exon" "$GFF" | \
awk '$3 ~ /CDS|exon/ {print $1"\t"$4-1"\t"$5}' > "$OUTDIR/coding_regions.bed"
echo "Coding regions found: $(wc -l < "$OUTDIR/coding_regions.bed")"

echo "Step 2: Extracting SNP positions and filtering..."
zcat ${PREFIX}.mafs.gz | awk 'NR>1 {print $1"\t"$2-1"\t"$2}' > "$OUTDIR/all_snps.bed"
bedtools intersect -v -a "$OUTDIR/all_snps.bed" -b "$OUTDIR/coding_regions.bed" > "$OUTDIR/non_coding_snps.bed"

echo "Step 3: Creating filtered TPED..."
awk '{print $1"\t"$3}' "$OUTDIR/non_coding_snps.bed" > "$OUTDIR/non_coding_snps_positions.txt"

# Remove previous failed files
rm -f "$OUTDIR/${PREFIX}_NONCODING.tped" "$OUTDIR/${PREFIX}_NONCODING.tfam"

echo "Step 4: Filtering TPED file..."
awk '
BEGIN {OFS="\t"}
NR==FNR {noncoding[$1":"$2] = 1; next}
FNR==1 {print; next}
($1":"$4 in noncoding) {print}
' "$OUTDIR/non_coding_snps_positions.txt" ${PREFIX}.tped > "$OUTDIR/${PREFIX}_NONCODING.tped"

cp ${PREFIX}.tfam "$OUTDIR/${PREFIX}_NONCODING.tfam"

echo "=== FINAL RESULTS ==="
TOTAL_SNPS=$(wc -l < ${PREFIX}.tped)
NONCODING_SNPS=$(wc -l < "$OUTDIR/non_coding_snps_positions.txt")
FILTERED_SNPS=$(wc -l < "$OUTDIR/${PREFIX}_NONCODING.tped")
echo "Total SNPs in original: $((TOTAL_SNPS - 1))"
echo "Non-coding SNPs identified: $NONCODING_SNPS"
echo "SNPs in filtered TPED: $((FILTERED_SNPS - 1))"
echo "Coding SNPs removed: $((TOTAL_SNPS - FILTERED_SNPS))"
echo "Done! Non-coding SNP files created in: $OUTDIR"
