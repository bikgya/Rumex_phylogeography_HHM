# ==============================================================
# Script 03c: RDA + LFMM2 Intersection (Adaptive SNP Detection)
# Project: rumex-landscape-phylogenomics
#
# Identifies SNPs under environmental selection using a two-method
# consensus approach:
#   - RDA (Redundancy Analysis): multivariate ordination,
#     outliers at +/- 3 SD on constrained axes
#   - LFMM2 (Latent Factor Mixed Model, ridge regression):
#     univariate per-variable test with K latent factors to
#     control for population structure
# Adaptive SNPs = intersection (consensus) of both methods.
#
# Statistics computed
#   - RDA constrained axes + variance explained
#   - LFMM2 p-values per SNP per BIO variable
#   - Genomic inflation factor (lambda) per BIO
#   - q-values (BH or Storey's qvalue) for FDR control
#   - Consensus adaptive SNP set (RDA intersect LFMM)
#
# Input (all from Script 03b)
#   geno_ld.rds        LD-pruned imputed genotype matrix
#   env_scaled.rds     standardised current bioclim
#   meta.rds           sample metadata (sample_id, population, etc.)
#
# Output (under out_dir)
#   LFMM2_pvalues_<BIO>.csv         per-variable p-values
#   LFMM_pvalues_lfmm2.csv          combined p-value matrix
#   LFMM_qvalues.csv                BH/q-value matrix
#   LFMM_sig_counts_Nbios.txt       summary table
#   QQ_<BIO>.{pdf,png}              per-variable QQ plots
#   QQ_combined.{pdf,png}           overlaid QQ plot
#   adaptive_SNPs_RDA_LFMM.csv      consensus table with q-values
#   adaptive_snps_final.txt         consensus SNP IDs
#   neutral_snps_all.txt            all non-adaptive SNP IDs
#   neutral_snps_10k_random.txt     random subset for control analyses
#   LFMM_adaptive_set.RData         all R objects bundled
#
# Dependencies: LEA (Bioc), qvalue (Bioc), vegan, ggplot2
# ==============================================================

# ----- Load packages -----
required_pkgs <- c("LEA","qvalue","vegan","ggplot2")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Install missing packages: ",
                                   paste(missing_pkgs, collapse = ", "),
                                   "\nLEA + qvalue are Bioconductor packages.")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ==============================================================
# USER CONFIGURATION
# ==============================================================
# Inputs (from Script 03b)
prep_dir <- "results/gea/preprocessing"

# Output directory
out_dir <- "results/gea/rda_lfmm"

# Samples to drop (e.g., high missingness flagged in 03b); empty c() if none
remove_samples <- c()        # e.g. c("RXX9")

# LFMM2 settings
best_K       <- 2            # latent factors (chosen from snmf cross-entropy)
q_working    <- 0.01         # FDR threshold for intersection
outlier_sd   <- 3            # RDA outlier threshold (SD on constrained axes)
n_axes_lfmm  <- 2            # number of RDA axes to consider

# Bioclim variables to test (must match preprocessing output columns)
selected_vars <- c("bio_2",  "bio_3",  "bio_8",  "bio_9",
                   "bio_14", "bio_15", "bio_17", "bio_18", "bio_19")

# Plot colours
COL_HIGHLIGHT <- "#D7191C"
COL_CHR_A     <- "#4A4A4A"
COL_CHR_B     <- "#9E9E9E"
# ==============================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

cat("==========================================================\n")
cat("  RDA + LFMM2 intersection -- adaptive SNP detection\n")
cat("  Started:", format(Sys.time()), "\n")
cat("==========================================================\n\n")

# --------------------------------------------------------------
# STEP 1: Load preprocessing outputs
# --------------------------------------------------------------
cat("[1] Loading preprocessing data...\n")
geno_imputed <- readRDS(file.path("..", "..", "..", prep_dir, "geno_ld.rds"))
env_scaled   <- readRDS(file.path("..", "..", "..", prep_dir, "env_scaled.rds"))
meta         <- readRDS(file.path("..", "..", "..", prep_dir, "meta.rds"))

cat(sprintf("  Individuals: %d | SNPs: %d\n", nrow(geno_imputed), ncol(geno_imputed)))

# --------------------------------------------------------------
# STEP 2: Optionally drop high-missingness samples
# --------------------------------------------------------------
if (length(remove_samples) > 0) {
    cat(sprintf("\n[2] Removing samples: %s\n", paste(remove_samples, collapse = ", ")))
    id_col <- intersect(c("sample_id","SampleID","ind","ID"), colnames(meta))[1]
    keep_idx  <- !(meta[[id_col]] %in% remove_samples)
    geno_filt <- geno_imputed[keep_idx, ]
    meta_filt <- meta[keep_idx, ]
    env_filt  <- env_scaled[keep_idx, ]
} else {
    cat("\n[2] No samples flagged for removal.\n")
    geno_filt <- geno_imputed; meta_filt <- meta; env_filt <- env_scaled
}
cat(sprintf("  Final: %d individuals\n", nrow(geno_filt)))

# --------------------------------------------------------------
# STEP 3: Impute any residual NAs (mode per SNP)
# --------------------------------------------------------------
na_count <- sum(is.na(geno_filt))
if (na_count > 0) {
    cat(sprintf("\n[3] Imputing %d residual NAs by per-SNP mode...\n", na_count))
    for (j in which(apply(geno_filt, 2, anyNA))) {
        v <- geno_filt[, j]
        v[is.na(v)] <- as.numeric(names(sort(table(v), decreasing = TRUE))[1])
        geno_filt[, j] <- v
    }
}

# --------------------------------------------------------------
# STEP 4: RDA on filtered dataset; outliers at +/- outlier_sd SD
# --------------------------------------------------------------
cat(sprintf("\n[4] Running RDA (outliers at +/- %d SD on first %d axes)...\n",
            outlier_sd, n_axes_lfmm))
env_std_filt <- as.data.frame(env_filt)
t1 <- Sys.time()
rda_filt <- vegan::rda(geno_filt ~ ., data = env_std_filt)
cat(sprintf("  RDA finished in %.1f min\n",
            difftime(Sys.time(), t1, units = "mins")))

snp_sc <- vegan::scores(rda_filt, choices = 1:n_axes_lfmm,
                        display = "species", scaling = "sites")

adaptive_RDA <- c()
for (ax in 1:n_axes_lfmm) {
    sc <- snp_sc[, ax]
    outs <- names(which(sc > mean(sc) + outlier_sd * sd(sc) |
                        sc < mean(sc) - outlier_sd * sd(sc)))
    adaptive_RDA <- c(adaptive_RDA, outs)
}
adaptive_RDA <- unique(adaptive_RDA)
cat(sprintf("  RDA adaptive SNPs (+/- %d SD): %d\n", outlier_sd, length(adaptive_RDA)))

# --------------------------------------------------------------
# STEP 5: Build individual-level env matrix (pop means)
# --------------------------------------------------------------
cat("\n[5] Building individual-level env matrix (population means)...\n")
pop_col <- intersect(c("Population","population","Pop","pop"), colnames(meta_filt))[1]
meta_filt$Population <- meta_filt[[pop_col]]
env_vars <- colnames(env_filt)
pop_env_std <- aggregate(as.data.frame(env_filt),
                         by  = list(Population = meta_filt$Population),
                         FUN = mean)
rownames(pop_env_std) <- pop_env_std$Population
pop_env_std <- as.data.frame(scale(pop_env_std[, env_vars, drop = FALSE]))
env_indiv   <- pop_env_std[meta_filt$Population, env_vars]

# --------------------------------------------------------------
# STEP 6: Write LFMM2 input files
# --------------------------------------------------------------
cat("\n[6] Writing LFMM2 inputs...\n")
write.table(geno_filt, "geno.lfmm", row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = " ")
env_files <- setNames(character(length(selected_vars)), selected_vars)
for (bio in selected_vars) {
    ef <- paste0("env_", bio, ".env")
    write.table(env_indiv[, bio, drop = FALSE], ef,
                row.names = FALSE, col.names = FALSE, quote = FALSE, sep = " ")
    env_files[bio] <- ef
}

# --------------------------------------------------------------
# STEP 7: Run LFMM2 per BIO variable
# --------------------------------------------------------------
cat(sprintf("\n[7] Running LFMM2 (K = %d) for each BIO variable...\n", best_K))
n_snps <- ncol(geno_filt)
mean_pvalues <- matrix(NA, nrow = n_snps, ncol = length(selected_vars),
                       dimnames = list(colnames(geno_filt), selected_vars))
geno_mat <- as.matrix(geno_filt)

for (bio in selected_vars) {
    env_vec <- as.matrix(read.table(env_files[bio]))
    mod <- LEA::lfmm2(input = geno_mat, env = env_vec, K = best_K,
                      effect.sizes = TRUE)
    pv  <- LEA::lfmm2.test(object = mod, input = geno_mat, env = env_vec,
                           full = FALSE, genomic.control = TRUE)
    p <- pv$pvalues
    p[p == 0] <- .Machine$double.xmin   # fix underflow
    p[p > 1]  <- 1
    mean_pvalues[, bio] <- p

    cat(sprintf("  %s: min-P = %.2e | median = %.3f\n",
                bio, min(p, na.rm = TRUE), median(p, na.rm = TRUE)))
    write.csv(data.frame(pvalue = p),
              paste0("LFMM2_pvalues_", bio, ".csv"), row.names = FALSE)
}
write.csv(data.frame(SNP = rownames(mean_pvalues), mean_pvalues),
          "LFMM_pvalues_lfmm2.csv", row.names = FALSE)

# --------------------------------------------------------------
# STEP 8: Q-value correction (FDR)
# --------------------------------------------------------------
cat("\n[8] q-value correction (FDR)...\n")
qvalues <- matrix(NA, nrow = nrow(mean_pvalues), ncol = ncol(mean_pvalues),
                  dimnames = dimnames(mean_pvalues))
lambda_vec <- setNames(numeric(length(selected_vars)), selected_vars)

for (bio in selected_vars) {
    p_bio  <- mean_pvalues[, bio]
    lambda <- median(qchisq(1 - p_bio, df = 1), na.rm = TRUE) / qchisq(0.5, df = 1)
    lambda_vec[bio] <- lambda
    q <- tryCatch(qvalue::qvalue(p_bio)$qvalues,
                  error = function(e) p.adjust(p_bio, method = "BH"))
    qvalues[, bio] <- q
    cat(sprintf("  %s | lambda = %.3f | q < 0.05: %d | q < 0.01: %d\n",
                bio, lambda, sum(q < 0.05, na.rm = TRUE),
                sum(q < 0.01, na.rm = TRUE)))
}
write.csv(data.frame(SNP = rownames(qvalues), qvalues),
          "LFMM_qvalues.csv", row.names = FALSE)

count_tbl <- data.frame(
    BIO = selected_vars, lambda = round(lambda_vec, 3),
    count_q0.10 = sapply(selected_vars, function(b) sum(qvalues[,b] <= 0.10, na.rm = TRUE)),
    count_q0.05 = sapply(selected_vars, function(b) sum(qvalues[,b] <= 0.05, na.rm = TRUE)),
    count_q0.01 = sapply(selected_vars, function(b) sum(qvalues[,b] <= 0.01, na.rm = TRUE))
)
write.table(count_tbl, "LFMM_sig_counts_Nbios.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

adaptive_LFMM <- rownames(qvalues)[
    apply(qvalues, 1, function(x) any(x < q_working, na.rm = TRUE))]
cat(sprintf("  LFMM-significant SNPs (q < %.2f, any BIO): %d\n",
            q_working, length(adaptive_LFMM)))

# --------------------------------------------------------------
# STEP 9: SNP position table (parse LG / scaffold IDs)
# --------------------------------------------------------------
cat("\n[9] Building SNP position table...\n")
snp_ids <- colnames(geno_filt)
parse_snp <- function(id) {
    m <- regmatches(id, regexec("^LG([0-9]+)_([0-9]+)$", id))[[1]]
    if (length(m) == 3)
        list(contig = paste0("LG", m[2]), chr = as.integer(m[2]),
             pos = as.numeric(m[3]), type = "LG")
    else
        list(contig = id, chr = NA, pos = NA, type = "scaffold")
}
parsed     <- lapply(snp_ids, parse_snp)
contig_raw <- sapply(parsed, `[[`, "contig")
chr_raw    <- sapply(parsed, `[[`, "chr")
pos_raw    <- sapply(parsed, `[[`, "pos")
type_raw   <- sapply(parsed, `[[`, "type")
max_lg     <- max(chr_raw, na.rm = TRUE)
scaffold_ids <- unique(contig_raw[type_raw == "scaffold"])
scaf_map     <- setNames(seq(max_lg + 1, max_lg + length(scaffold_ids)), scaffold_ids)
chr_codes    <- chr_raw
chr_codes[type_raw == "scaffold"] <- scaf_map[contig_raw[type_raw == "scaffold"]]
pos_codes    <- pos_raw
pos_codes[type_raw == "scaffold"] <- seq_along(which(type_raw == "scaffold"))

snp_info <- data.frame(SNP = snp_ids, CHR = as.integer(chr_codes),
                       CHR_raw = contig_raw, BP = as.integer(pos_codes),
                       type = type_raw, SNP_index = seq_along(snp_ids),
                       stringsAsFactors = FALSE)

# --------------------------------------------------------------
# STEP 10: QQ plots (per BIO + combined)
# --------------------------------------------------------------
cat("\n[10] Creating QQ plots...\n")
make_qq <- function(bio_var, out_base) {
    p <- mean_pvalues[, bio_var]; p <- p[!is.na(p)]; n <- length(p)
    obs_all <- -log10(sort(p)); exp_all <- -log10(ppoints(n))
    lam <- round(lambda_vec[bio_var], 3)
    bonf_thresh <- -log10(0.05 / n)
    sig_mask <- obs_all > bonf_thresh
    set.seed(42)
    keep_nsig <- sample(which(!sig_mask),
                        size = min(sum(!sig_mask), round(0.10 * n)))
    keep_idx <- sort(c(which(sig_mask), keep_nsig))
    obs <- obs_all[keep_idx]; exp <- exp_all[keep_idx]

    do_plot <- function() {
        par(mar = c(5, 5, 3, 2), cex = 1)
        plot(exp, obs, pch = 20, cex = 0.4,
             col = ifelse(obs > bonf_thresh, COL_HIGHLIGHT, "grey60"),
             xlab = expression(Expected ~ -log[10](p)),
             ylab = expression(Observed ~ -log[10](p)),
             main = paste0(bio_var, " (lambda = ", lam, ")"),
             xlim = c(0, max(exp_all) * 1.05),
             ylim = c(0, max(obs_all) * 1.05))
        abline(0, 1, col = "#2C7BB6", lwd = 1.5)
    }
    cairo_pdf(paste0(out_base, ".pdf"), width = 5, height = 5,
              pointsize = 10, antialias = "none"); do_plot(); dev.off()
    png(paste0(out_base, ".png"), width = 5, height = 5, units = "in", res = 200)
    do_plot(); dev.off()
}
for (bio in selected_vars) make_qq(bio, paste0("QQ_", bio))

# Combined QQ
colors_bio <- colorRampPalette(c("#2C7BB6", "#D7191C"))(length(selected_vars))
do_combined_qq <- function() {
    par(mar = c(5, 5, 3, 2), cex = 1); first <- TRUE
    for (i in seq_along(selected_vars)) {
        bio <- selected_vars[i]
        p <- mean_pvalues[, bio]; p <- p[!is.na(p)]; n <- length(p)
        obs_all <- -log10(sort(p)); exp_all <- -log10(ppoints(n))
        top_idx <- which(obs_all >= quantile(obs_all, 0.95))
        set.seed(42 + i)
        keep_rest <- sample(which(obs_all < quantile(obs_all, 0.95)),
                            size = min(n - length(top_idx), round(0.05 * n)))
        keep_idx <- sort(c(top_idx, keep_rest))
        if (first) {
            plot(exp_all[keep_idx], obs_all[keep_idx], pch = 20, cex = 0.3,
                 col = colors_bio[i],
                 xlab = expression(Expected ~ -log[10](p)),
                 ylab = expression(Observed ~ -log[10](p)),
                 main = "QQ plots - all BIO variables",
                 xlim = c(0, max(exp_all) * 1.05),
                 ylim = c(0, max(obs_all) * 1.05))
            first <- FALSE
        } else {
            points(exp_all[keep_idx], obs_all[keep_idx], pch = 20, cex = 0.3, col = colors_bio[i])
        }
    }
    abline(0, 1, col = "black", lwd = 1.5, lty = 2)
    legend("topleft",
           legend = paste0(selected_vars, " (lambda=", round(lambda_vec, 2), ")"),
           col = colors_bio, pch = 20, pt.cex = 0.8, bty = "n", cex = 0.72)
}
cairo_pdf("QQ_combined.pdf", width = 6, height = 5, pointsize = 10, antialias = "none")
do_combined_qq(); dev.off()
png("QQ_combined.png", width = 6, height = 5, units = "in", res = 200)
do_combined_qq(); dev.off()

# --------------------------------------------------------------
# STEP 11: RDA-LFMM intersection (consensus adaptive SNPs)
# --------------------------------------------------------------
cat("\n[11] Building RDA intersect LFMM consensus set...\n")
adaptive_intersect <- intersect(adaptive_RDA, adaptive_LFMM)
cat(sprintf("  RDA adaptive: %d | LFMM adaptive: %d | Intersection: %d\n",
            length(adaptive_RDA), length(adaptive_LFMM), length(adaptive_intersect)))

intersect_table <- data.frame(
    SNP = adaptive_intersect,
    CHR_raw   = snp_info$CHR_raw[match(adaptive_intersect, snp_info$SNP)],
    LG        = snp_info$CHR[match(adaptive_intersect, snp_info$SNP)],
    BP        = snp_info$BP[match(adaptive_intersect, snp_info$SNP)],
    type      = snp_info$type[match(adaptive_intersect, snp_info$SNP)],
    SNP_index = snp_info$SNP_index[match(adaptive_intersect, snp_info$SNP)],
    stringsAsFactors = FALSE
)
for (bio in selected_vars)
    intersect_table[[paste0("q_", bio)]] <- qvalues[adaptive_intersect, bio]
intersect_table$min_q   <- apply(qvalues[adaptive_intersect, , drop = FALSE],
                                  1, min, na.rm = TRUE)
intersect_table$top_BIO <- selected_vars[apply(
    qvalues[adaptive_intersect, , drop = FALSE], 1, which.min)]
intersect_table <- intersect_table[order(intersect_table$min_q), ]
write.csv(intersect_table, "adaptive_SNPs_RDA_LFMM.csv", row.names = FALSE)

# Neutral SNP sets
neutral_final <- setdiff(colnames(geno_filt), adaptive_intersect)
set.seed(42)
neutral_subset_10k <- sample(neutral_final, min(10000, length(neutral_final)))
writeLines(adaptive_intersect, "adaptive_snps_final.txt")
writeLines(neutral_final,      "neutral_snps_all.txt")
writeLines(neutral_subset_10k, "neutral_snps_10k_random.txt")

# Per-LG and per-BIO summaries
lg_sum  <- as.data.frame(table(intersect_table$CHR_raw))
colnames(lg_sum) <- c("LG_or_scaffold", "N_adaptive_SNPs")
write.csv(lg_sum[order(-lg_sum$N_adaptive_SNPs), ], "adaptive_SNPs_per_LG.csv", row.names = FALSE)
bio_sum <- as.data.frame(table(intersect_table$top_BIO))
colnames(bio_sum) <- c("BIO_variable", "N_adaptive_SNPs")
write.csv(bio_sum[order(-bio_sum$N_adaptive_SNPs), ], "adaptive_SNPs_per_BIO.csv", row.names = FALSE)

# --------------------------------------------------------------
# STEP 12: Save R objects
# --------------------------------------------------------------
adaptive_snps_RDA_only  <- adaptive_RDA
adaptive_snps_LFMM_only <- adaptive_LFMM
adaptive_snps_intersect <- adaptive_intersect
save(geno_filt, meta_filt, env_filt, mean_pvalues, qvalues,
     best_K, count_tbl, selected_vars, snp_info, lambda_vec,
     adaptive_snps_RDA_only, adaptive_snps_LFMM_only, adaptive_snps_intersect,
     neutral_final, neutral_subset_10k, intersect_table,
     file = "LFMM_adaptive_set.RData")

cat("\n", strrep("=", 60), "\n", sep = "")
cat("  GEA ANALYSIS COMPLETE:", format(Sys.time()), "\n")
cat(sprintf("  RDA adaptive (+/- %d SD):    %d\n", outlier_sd, length(adaptive_RDA)))
cat(sprintf("  LFMM significant (q < %.2f): %d\n", q_working, length(adaptive_LFMM)))
cat(sprintf("  CONSENSUS adaptive SNPs:    %d\n", length(adaptive_intersect)))
cat(sprintf("  Neutral SNPs (all):         %d\n", length(neutral_final)))
cat(sprintf("  Neutral SNPs (10K subset):  %d\n", length(neutral_subset_10k)))
cat(strrep("=", 60), "\n")
