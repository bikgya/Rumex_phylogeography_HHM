#!/usr/bin/env bash
# ==============================================================
# Script 03: Genetic Diversity & Pairwise FST (ANGSD)
# Project: rumex-landscape-phylogenomics
#
# Computes per-population diversity statistics and pairwise FST
# across 162 Rumex hastatus samples (27 populations) using
# ANGSD's genotype-likelihood framework. Population labels are
# inferred automatically from sample IDs (alphabetic prefix
# before the first digit, e.g., RCX2 -> RCX, RDQA7 -> RDQA).
#
# Statistics computed
# --------------------------------------------------------------
#   Per population (genome-wide + 50 kb sliding windows):
#     - Watterson's theta (theta_W)  : population mutation rate
#     - Nucleotide diversity (theta_pi, pi)
#     - Tajima's D                    : neutrality test
#
#   Between populations:
#     - Pairwise F_ST (Hudson estimator, from joint 2D SFS)
#     - Symmetric F_ST matrix (CSV)
#
# Input
# --------------------------------------------------------------
#   - Sorted, indexed BAM files (from Script 01a)
#   - Reference genome FASTA (same as Script 01b)
#   - Ancestral sequence (reference is used here = folded SFS)
#
# Output (under ${OUTDIR})
# --------------------------------------------------------------
#   samples_and_pops.txt              sample -> population map
#   bamlists_per_pop/<pop>.bamlist.txt
#   angsd_out/<pop>.saf.idx           per-site allele freq likelihoods
#   angsd_out/<pop>.sfs               folded 1D SFS
#   diversity_results/<pop>.thetas.idx    per-site thetas
#   diversity_results/<pop>.theta.global.txt    genome-wide stats
#   diversity_results/<pop>.theta.windows.txt   sliding-window stats
#   fst_results/<p1>_<p2>.fst.stats.txt        pairwise FST output
#   fst_results/pairwise_fst_matrix.csv        symmetric FST matrix
#
# Tools: ANGSD v0.935 (angsd, realSFS, thetaStat)
# ==============================================================

set -euo pipefail

# ----- Activate ANGSD conda environment -----
if command -v conda >/dev/null 2>&1; then
    set +u
    eval "$(conda shell.bash hook)"
    conda activate ANGSD
    set -u
else
    echo "ERROR: conda not found. Load conda or add to PATH first."
    exit 1
fi

# ==============================================================
# USER CONFIGURATION
# ==============================================================
REF_FA="data/raw/reference/genome.fasta"        # Reference FASTA
ANC_FA="data/raw/reference/genome.fasta"        # Ancestral (= ref for folded SFS)
BAM_DIR="results/snp_calling/bam"               # Sorted BAMs from Script 01a
SAMPLE_LIST="data/sample_info/sample_names.txt" # One sample ID per line
OUTDIR="results/diversity"                      # All outputs go here

# ANGSD filtering parameters
NCORE=16
MINMAPQ=20
MINQ=20
MININD_FRAC=0.8       # min individuals = fraction of pop size
MINMAF=0.01
SNP_PVAL=1e-6
MAXDEPTH=1000

# Sliding-window settings (theta stats)
WINDOW=50000
STEP=10000

# Parallel FST jobs
MAX_FST_JOBS=10
# ==============================================================

mkdir -p "${OUTDIR}"
cd "${OUTDIR}"

# --------------------------------------------------------------
# STEP 1: Infer populations from sample IDs
#         (alphabetic prefix before first digit)
# --------------------------------------------------------------
echo "=== STEP 1: Inferring populations from sample names ==="
awk '{ match($0, /^[^0-9]+/); pop=substr($0,RSTART,RLENGTH); print $0"\t"pop }' \
    "../../${SAMPLE_LIST}" > samples_and_pops.txt

cut -f2 samples_and_pops.txt | sort -u > populations.txt
echo "Detected $(wc -l < populations.txt) populations:"
cat populations.txt
echo

# --------------------------------------------------------------
# STEP 2: Build one BAM list per population
# --------------------------------------------------------------
echo "=== STEP 2: Building per-population BAM lists ==="
mkdir -p bamlists_per_pop
> missing_bams.txt

while read -r sample pop; do
    bam="../../${BAM_DIR}/${sample}.sorted.bam"
    if [[ -f "$bam" ]]; then
        echo "$bam" >> "bamlists_per_pop/${pop}.bamlist.txt"
    else
        echo "$sample" >> missing_bams.txt
    fi
done < samples_and_pops.txt

[[ -s missing_bams.txt ]] && echo "WARNING: $(wc -l < missing_bams.txt) BAMs missing (see missing_bams.txt)"
echo

# --------------------------------------------------------------
# STEP 3: Per-population SAF, SFS, and theta statistics
# --------------------------------------------------------------
echo "=== STEP 3: Per-population diversity (SAF -> SFS -> theta) ==="
mkdir -p angsd_out diversity_results

for popfile in bamlists_per_pop/*.bamlist.txt; do
    pop=$(basename "${popfile}" .bamlist.txt)
    n_ind=$(wc -l < "${popfile}")
    min_ind=$(awk -v n="$n_ind" -v p="$MININD_FRAC" 'BEGIN { printf "%d", n*p + 0.5 }')

    echo "--- Population: ${pop} (n=${n_ind}, minInd=${min_ind}) ---"

    # 3a. Site allele frequency likelihoods (SAF)
    angsd -b "${popfile}" \
          -ref "../../${REF_FA}" -anc "../../${ANC_FA}" \
          -out "angsd_out/${pop}" \
          -GL 2 -doSaf 1 -doMajorMinor 1 -doMaf 1 \
          -minMapQ ${MINMAPQ} -minQ ${MINQ} -minInd ${min_ind} \
          -minMaf ${MINMAF} -SNP_pval ${SNP_PVAL} \
          -doCounts 1 -setMaxDepth ${MAXDEPTH} \
          -remove_bads 1 -only_proper_pairs 1 \
          -nThreads ${NCORE}

    # 3b. Folded SFS (reference used as ancestral)
    realSFS "angsd_out/${pop}.saf.idx" -fold 1 -P ${NCORE} > "angsd_out/${pop}.sfs"

    # 3c. Per-site thetas
    realSFS saf2theta "angsd_out/${pop}.saf.idx" \
        -sfs "angsd_out/${pop}.sfs" \
        -outname "diversity_results/${pop}" \
        -fold 1

    # 3d. Genome-wide stats (Watterson, pi, Tajima's D)
    thetaStat do_stat "diversity_results/${pop}.thetas.idx" \
        > "diversity_results/${pop}.theta.global.txt"

    # 3e. Sliding-window stats
    thetaStat do_stat "diversity_results/${pop}.thetas.idx" \
        -win ${WINDOW} -step ${STEP} \
        > "diversity_results/${pop}.theta.windows.txt"

    echo "Finished ${pop}"
    echo
done

# --------------------------------------------------------------
# STEP 4: Pairwise F_ST (parallel) and matrix
# --------------------------------------------------------------
echo "=== STEP 4: Pairwise F_ST ==="
mkdir -p fst_results
mapfile -t pops < populations.txt
cd angsd_out

fst_exists() {
    [[ -s "../fst_results/${1}_${2}.fst.stats.txt" ]] && return 0
    [[ -s "../fst_results/${2}_${1}.fst.stats.txt" ]] && return 0
    return 1
}

compute_fst() {
    local p1="$1" p2="$2"
    if fst_exists "$p1" "$p2"; then
        echo "SKIP: $p1 vs $p2 (already done)"
        return 0
    fi
    if [[ ! -s "${p1}.saf.idx" ]] || [[ ! -s "${p2}.saf.idx" ]]; then
        echo "MISSING SAF for $p1 or $p2 - skipping"
        return 1
    fi
    echo "START: $p1 vs $p2"
    realSFS "${p1}.saf.idx" "${p2}.saf.idx" -P ${THREADS} > "${p1}_${p2}.sfs"
    realSFS fst index "${p1}.saf.idx" "${p2}.saf.idx" \
        -sfs "${p1}_${p2}.sfs" -fstout "${p1}_${p2}"
    realSFS fst stats "${p1}_${p2}.fst.idx" \
        > "../fst_results/${p1}_${p2}.fst.stats.txt"
    echo "DONE:  $p1 vs $p2"
}

export -f compute_fst fst_exists
export THREADS=${NCORE}

# Launch all pairs (j > i) with a concurrency cap
running=0
for ((i=0; i<${#pops[@]}-1; i++)); do
    for ((j=i+1; j<${#pops[@]}; j++)); do
        compute_fst "${pops[i]}" "${pops[j]}" &
        running=$((running + 1))
        if (( running >= MAX_FST_JOBS )); then
            wait
            running=0
        fi
    done
done
wait
cd ..

# --------------------------------------------------------------
# STEP 5: Build symmetric F_ST matrix
# --------------------------------------------------------------
echo "=== STEP 5: Building pairwise F_ST matrix ==="
matrix="fst_results/pairwise_fst_matrix.csv"
{
    printf "POP"
    for p in "${pops[@]}"; do printf ",%s" "$p"; done
    printf "\n"

    for p1 in "${pops[@]}"; do
        printf "%s" "$p1"
        for p2 in "${pops[@]}"; do
            if [[ "$p1" == "$p2" ]]; then
                printf ",0"
            else
                f1="fst_results/${p1}_${p2}.fst.stats.txt"
                f2="fst_results/${p2}_${p1}.fst.stats.txt"
                stat_file=""
                [[ -s "$f1" ]] && stat_file="$f1"
                [[ -z "$stat_file" && -s "$f2" ]] && stat_file="$f2"
                if [[ -n "$stat_file" ]]; then
                    val=$(awk '/global Fst/ {print $3; exit}' "$stat_file" 2>/dev/null)
                    printf ",%s" "${val:-NA}"
                else
                    printf ",NA"
                fi
            fi
        done
        printf "\n"
    done
} > "${matrix}"

# --------------------------------------------------------------
# SUMMARY
# --------------------------------------------------------------
echo "=============================================="
echo " Diversity & FST analysis complete."
echo "=============================================="
echo "  Per-pop diversity : ${OUTDIR}/diversity_results/"
echo "  Pairwise FST stats: ${OUTDIR}/fst_results/"
echo "  FST matrix        : ${OUTDIR}/${matrix}"
echo "  Missing BAMs (if any): ${OUTDIR}/missing_bams.txt"
echo "=============================================="
