# ==============================================================
# Script 03b: GEA Preprocessing (ANGSD output -> RDS files)
# Project: rumex-landscape-phylogenomics
#
# Prepares the genotype matrix and environmental data for GEA
# analyses (Script 03c). Loads ANGSD hard-called genotypes,
# applies post-calling QC filters, extracts current and future
# bioclim values at sample locations, imputes residual missing
# genotypes with sNMF, and LD-prunes within each scaffold.
#
# Pipeline
#   1. Load ANGSD .geno.gz + .mafs.gz (from Script 03a)
#   2. Apply QC: per-site missingness <= 20%, MAF >= 0.05;
#      flag individuals with > 30% missing data
#   3. Extract current + future bioclim at occurrence coordinates,
#      standardise (z-score) and store transformation parameters
#   4. sNMF cross-entropy across K = 2:8 to choose best K,
#      then impute missing genotypes using mode method
#   5. LD-prune within scaffold (r^2 < 0.5, sliding window of 5)
#   6. Save genotype, env, and metadata matrices as .rds
#
# Input
#   - ANGSD hard genotypes (Script 03a):
#       gea_snps.geno.gz, gea_snps.mafs.gz
#   - Sample metadata CSV
#       (required columns: sample_id, population, latitude, longitude)
#   - WorldClim bioclim rasters (current and future)
#
# Output (under out_dir)
#   geno_ld.rds              LD-pruned imputed genotype matrix
#   snp_info_ld.rds          SNP coordinates after LD pruning
#   env_scaled.rds           scaled current bioclim (individuals)
#   env_future_scaled.rds    future bioclim, scaled w/ current centring
#   meta.rds                 sample metadata
#   bioclim_correlation.pdf  correlation plot of selected variables
#   snmf_cross_entropy.pdf   sNMF model-fit plot
#
# Dependencies: data.table, vegan, terra, LEA, lfmm, ggplot2,
#               dplyr, corrplot, geosphere
# ==============================================================

# ----- Load packages (assumes installed; see envs/) -----
required_pkgs <- c("data.table","vegan","terra","LEA","lfmm",
                   "ggplot2","dplyr","corrplot","geosphere")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Install missing packages: ",
                                   paste(missing_pkgs, collapse = ", "),
                                   "\nLEA + qvalue are Bioconductor packages.")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ==============================================================
# USER CONFIGURATION
# ==============================================================
# ANGSD output (from Script 03a)
angsd_geno <- "results/gea/genotypes/gea_snps.geno.gz"
angsd_mafs <- "results/gea/genotypes/gea_snps.mafs.gz"

# Sample metadata CSV (rows MUST match the ANGSD bamlist order)
# Required columns: sample_id, population, latitude, longitude
coords_file <- "data/sample_info/gea_pop_info.csv"

# WorldClim raster directories
current_bioclim_dir <- "data/raw/env/current"
future_bioclim_dir  <- "data/raw/env/future"

# Output directory
out_dir <- "results/gea/preprocessing"

# 9 selected bioclim variables (must match raster layer names)
selected_bios <- c("bio_2",  "bio_3",  "bio_8",  "bio_9",
                   "bio_14", "bio_15", "bio_17", "bio_18", "bio_19")

# Filters
site_missing_max <- 0.20    # drop sites with > 20% missing
ind_missing_warn <- 0.30    # warn (don't drop) individuals > 30% missing
maf_min          <- 0.05    # minimum MAF after recalculation

# sNMF imputation
snmf_K_range  <- 2:8
snmf_reps     <- 5
snmf_cpu      <- 8

# LD pruning
ld_r2_threshold <- 0.5
ld_window       <- 5        # compare against last N retained SNPs per scaffold
# ==============================================================

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------
# STEP 1: Load ANGSD genotypes + SNP info
# --------------------------------------------------------------
cat("[1] Loading ANGSD genotype data...\n")
geno_raw <- fread(angsd_geno, header = FALSE)
cat(sprintf("  Raw matrix: %d SNPs x %d individuals\n",
            nrow(geno_raw), ncol(geno_raw)))

# Transpose: rows = individuals, columns = SNPs
geno <- t(as.matrix(geno_raw))
rm(geno_raw); gc()
geno[geno == -1] <- NA   # ANGSD missing code

mafs <- fread(angsd_mafs, header = TRUE)
snp_info <- data.frame(
    chr   = mafs$chromo,  pos   = mafs$position,
    major = mafs$major,   minor = mafs$minor,
    maf   = mafs$knownEM,
    snp_id = paste0(mafs$chromo, "_", mafs$position),
    stringsAsFactors = FALSE
)
colnames(geno) <- snp_info$snp_id

meta <- read.csv(coords_file, stringsAsFactors = FALSE)
rownames(geno) <- meta$sample_id
stopifnot(nrow(geno) == nrow(meta))
cat(sprintf("  Loaded %d individuals from %d populations\n",
            nrow(meta), length(unique(meta$population))))

# --------------------------------------------------------------
# STEP 2: Post-calling QC filters
# --------------------------------------------------------------
cat("\n[2] Quality filters...\n")
keep_sites <- colMeans(is.na(geno)) <= site_missing_max
geno     <- geno[, keep_sites]
snp_info <- snp_info[keep_sites, ]
cat(sprintf("  After site missingness <= %.0f%%: %d SNPs\n",
            100 * site_missing_max, ncol(geno)))

ind_miss <- rowMeans(is.na(geno))
high_miss <- which(ind_miss > ind_missing_warn)
if (length(high_miss) > 0) {
    cat(sprintf("  WARNING: %d individuals > %.0f%% missing: %s\n",
                length(high_miss), 100 * ind_missing_warn,
                paste(meta$sample_id[high_miss], collapse = ", ")))
}

maf_recalc <- colMeans(geno, na.rm = TRUE) / 2
maf_recalc <- pmin(maf_recalc, 1 - maf_recalc)
keep_maf   <- maf_recalc >= maf_min
geno     <- geno[, keep_maf]
snp_info <- snp_info[keep_maf, ]
cat(sprintf("  After MAF >= %.2f: %d SNPs\n", maf_min, ncol(geno)))
cat(sprintf("  FINAL: %d individuals x %d SNPs\n", nrow(geno), ncol(geno)))

# --------------------------------------------------------------
# STEP 3: Extract bioclim variables
# --------------------------------------------------------------
cat("\n[3] Extracting bioclim variables...\n")
bio_files <- list.files(current_bioclim_dir, pattern = "\\.tif$", full.names = TRUE)
bio_stack <- rast(bio_files)
coords_sp <- vect(meta, geom = c("longitude", "latitude"), crs = "EPSG:4326")
env_current <- terra::extract(bio_stack, coords_sp)[, -1]
cat(sprintf("  Raster layers: %s\n", paste(names(env_current), collapse = ", ")))

env_sel <- env_current[, selected_bios]
if (any(is.na(env_sel)))
    cat(sprintf("  WARNING: %d individuals with NA bioclim values\n",
                sum(rowSums(is.na(env_sel)) > 0)))

# Standardise (z-score), store centre/scale for future climate
env_scaled <- scale(env_sel)
env_center <- attr(env_scaled, "scaled:center")
env_scale  <- attr(env_scaled, "scaled:scale")
env_scaled <- as.data.frame(env_scaled)
rownames(env_scaled) <- meta$sample_id

pdf(file.path(out_dir, "bioclim_correlation.pdf"), width = 8, height = 8)
corrplot(cor(env_scaled), method = "number", type = "lower", tl.cex = 0.8,
         title = "Correlation: 9 selected bioclim variables",
         mar = c(0, 0, 2, 0))
dev.off()

# Future bioclim (use CURRENT centre/scale - critical for genomic offset)
future_files <- list.files(future_bioclim_dir, pattern = "\\.tif$", full.names = TRUE)
bio_future   <- rast(future_files)
env_future_raw <- terra::extract(bio_future, coords_sp)[, -1]
env_future_sel <- env_future_raw[, selected_bios]
env_future_scaled <- as.data.frame(scale(env_future_sel,
                                         center = env_center, scale = env_scale))
rownames(env_future_scaled) <- meta$sample_id
cat("  Current + future bioclim extracted and standardised\n")

# --------------------------------------------------------------
# STEP 4: Impute missing genotypes with sNMF
# --------------------------------------------------------------
cat(sprintf("\n[4] sNMF imputation (K = %d:%d, %d reps)...\n",
            min(snmf_K_range), max(snmf_K_range), snmf_reps))
pct_miss <- round(100 * sum(is.na(geno)) / length(geno), 2)
cat(sprintf("  Missing before imputation: %.2f%%\n", pct_miss))

lea_dir <- file.path(out_dir, "LEA_temp")
dir.create(lea_dir, showWarnings = FALSE)
lfmm_file     <- file.path(lea_dir, "gea.lfmm")
geno_lea_file <- file.path(lea_dir, "gea.geno")
write.lfmm(geno, lfmm_file)
write.geno(geno, geno_lea_file)

snmf_proj <- snmf(geno_lea_file, K = snmf_K_range, repetitions = snmf_reps,
                  entropy = TRUE, project = "new", CPU = snmf_cpu)

pdf(file.path(out_dir, "snmf_cross_entropy.pdf"), width = 8, height = 5)
plot(snmf_proj, col = "steelblue", pch = 19, main = "sNMF Cross-Entropy")
dev.off()

ce_vals      <- sapply(snmf_K_range, function(k) min(cross.entropy(snmf_proj, K = k)))
best_K_snmf  <- snmf_K_range[which.min(ce_vals)]
cat(sprintf("  Best K for imputation: %d\n", best_K_snmf))

best_run     <- which.min(cross.entropy(snmf_proj, K = best_K_snmf))
imputed_file <- impute(snmf_proj, lfmm_file,
                       method = "mode", K = best_K_snmf, run = best_run)
geno_imp <- read.lfmm(imputed_file)
rownames(geno_imp) <- meta$sample_id
colnames(geno_imp) <- snp_info$snp_id
cat(sprintf("  Imputation complete. Remaining NAs: %d\n", sum(is.na(geno_imp))))

# --------------------------------------------------------------
# STEP 5: LD pruning (within scaffold)
# --------------------------------------------------------------
cat(sprintf("\n[5] LD pruning (r^2 < %.2f, sliding window of %d)...\n",
            ld_r2_threshold, ld_window))
scaffolds <- unique(snp_info$chr)
keep_ld <- c()
for (scf in scaffolds) {
    idx <- which(snp_info$chr == scf)
    if (length(idx) <= 1) { keep_ld <- c(keep_ld, idx); next }
    ord      <- idx[order(snp_info$pos[idx])]
    retained <- ord[1]
    for (i in 2:length(ord)) {
        recent <- tail(retained, ld_window)
        max_r2 <- max(sapply(recent, function(j)
            cor(geno_imp[, j], geno_imp[, ord[i]],
                use = "pairwise.complete.obs")^2), na.rm = TRUE)
        if (max_r2 < ld_r2_threshold) retained <- c(retained, ord[i])
    }
    keep_ld <- c(keep_ld, retained)
}
geno_ld     <- geno_imp[, keep_ld]
snp_info_ld <- snp_info[keep_ld, ]
cat(sprintf("  LD-pruned SNPs: %d (removed: %d)\n",
            ncol(geno_ld), ncol(geno_imp) - ncol(geno_ld)))

# --------------------------------------------------------------
# STEP 6: Save outputs
# --------------------------------------------------------------
cat("\n[6] Saving RDS files...\n")
saveRDS(geno_ld,           file.path(out_dir, "geno_ld.rds"))
saveRDS(snp_info_ld,       file.path(out_dir, "snp_info_ld.rds"))
saveRDS(env_scaled,        file.path(out_dir, "env_scaled.rds"))
saveRDS(env_future_scaled, file.path(out_dir, "env_future_scaled.rds"))
saveRDS(meta,              file.path(out_dir, "meta.rds"))

cat("\n", strrep("=", 60), "\n", sep = "")
cat("  Preprocessing complete:", format(Sys.time()), "\n")
cat(sprintf("  Individuals:          %d\n", nrow(geno_ld)))
cat(sprintf("  LD-pruned SNPs:       %d\n", ncol(geno_ld)))
cat(sprintf("  Bioclim variables:    %d\n", ncol(env_scaled)))
cat("  Next step: Script 03c (RDA + LFMM2 intersection)\n")
cat(strrep("=", 60), "\n")
