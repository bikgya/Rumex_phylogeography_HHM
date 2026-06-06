# ==============================================================
# Script 10: Ecological Niche Divergence (PCA-env, ecospat)
# Project: rumex-landscape-phylogenomics
#
# Tests whether the two phylogeographic lineages of
# Rumex hastatus - West Himalaya (WH) and Hengduan Mountains
# (HM) - occupy statistically distinct climatic niches.
#
# Pipeline (ecospat PCA-env framework, Broennimann et al. 2012)
#   1. Split occurrences by longitude (WH vs. HM)
#   2. Sample background environment from the study area
#   3. Calibrate PCA on combined background + occurrence env data
#   4. Build kernel density niches (z_WH, z_HM) in PCA space
#   5. Compute Schoener's D and Hellinger's I overlap metrics
#   6. Niche equivalency test  (H0: niches are equivalent)
#   7. Niche similarity tests  (both directions)
#   8. Generate publication figures + CSV summaries
#
# Input
#   - occ_csv : whole-region occurrence CSV (Longitude, Latitude)
#   - bio_dir : folder with WorldClim/bioclim GeoTIFFs
#               (bio2.tif, bio3.tif, ... - the VIF-selected set
#               from Script 08b)
#
# Output (under wd/)
#   Data    : niche_overlap_results.csv,
#             niche_background_env.csv,
#             niche_WH_scores.csv, niche_HM_scores.csv
#   Figures : Figure_PCAenv_niche_overlap.{pdf,png}
#             Figure_equivalency_test.pdf
#             Figure_similarity_test.pdf
#             Figure_PCA_biplot_loadings.pdf
#             Figure_niche_overlap_SINGLE.{pdf,png}
#
# Dependencies: ecospat, ade4, terra, ggplot2, dplyr, scales, gridExtra
# ==============================================================

# ----- Load packages -----
required_pkgs <- c("ecospat","ade4","terra","ggplot2","dplyr","scales","gridExtra")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Install missing packages first: ",
                                   paste(missing_pkgs, collapse = ", "))
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ==============================================================
# USER CONFIGURATION
# ==============================================================
# Working directory (all outputs saved here)
wd <- "results/niche_divergence"

# Occurrence CSV
occ_csv <- "data/occurrences/rumex_whole.csv"
lon_col <- "Longitude"
lat_col <- "Latitude"

# Bioclimatic raster directory (CURRENT climate)
bio_dir <- "data/raw/env/current"

# VIF-selected bioclim variables (from Script 08b); files named bioN.tif
bio_vars <- c(2, 3, 8, 9, 14, 15, 17, 18, 19)

# WH / HM split longitude
# West of split = WH (West Himalaya); East of split = HM (Hengduan)
split_lon <- 85.0

# Background extent (full SDM calibration area)
bg_lon_min <- 69.25;  bg_lon_max <- 104.10
bg_lat_min <- 26.22;  bg_lat_max <- 35.25

# Background samples
n_background  <- 10000

# Permutations for equivalency / similarity tests
#  - 99-499 is typical for exploratory work
#  - 1000+ recommended for final manuscript submission
n_permutations <- 499

# Grid resolution for ecospat density estimation
R_grid <- 100
# ==============================================================

dir.create(wd, recursive = TRUE, showWarnings = FALSE)
setwd(wd)

# --------------------------------------------------------------
# STEP 1: Load occurrences and split into WH / HM
# --------------------------------------------------------------
cat("\n-- STEP 1: Loading and splitting occurrence data --\n")
if (!file.exists(occ_csv)) stop("Occurrence CSV not found: ", occ_csv)

occ <- read.csv(occ_csv, stringsAsFactors = FALSE)
occ <- occ[!is.na(occ[[lon_col]]) & !is.na(occ[[lat_col]]), ]
cat(sprintf("  Total occurrences: %d\n", nrow(occ)))

occ_WH <- occ[occ[[lon_col]] <  split_lon, ]
occ_HM <- occ[occ[[lon_col]] >  split_lon, ]
cat(sprintf("  WH (longitude < %.1f): %d records\n", split_lon, nrow(occ_WH)))
cat(sprintf("  HM (longitude > %.1f): %d records\n", split_lon, nrow(occ_HM)))

if (nrow(occ_WH) < 5 | nrow(occ_HM) < 5)
    stop("Too few points in one group; check split_lon or occurrence data.")

xy_WH <- occ_WH[, c(lon_col, lat_col)]
xy_HM <- occ_HM[, c(lon_col, lat_col)]

# --------------------------------------------------------------
# STEP 2: Load bioclim rasters, sample background, extract env
# --------------------------------------------------------------
cat("\n-- STEP 2: Loading bioclim rasters and extracting environment --\n")
bio_files <- file.path(bio_dir, paste0("bio", bio_vars, ".tif"))
missing <- bio_files[!file.exists(bio_files)]
if (length(missing) > 0)
    stop("Missing raster files:\n", paste(" ", missing, collapse = "\n"))

bio_stack <- rast(bio_files)
names(bio_stack) <- paste0("bio", bio_vars)
bio_crop  <- crop(bio_stack, ext(bg_lon_min, bg_lon_max, bg_lat_min, bg_lat_max))

# Background sample
set.seed(42)
bg_pts <- spatSample(bio_crop, size = n_background, method = "random",
                     na.rm = TRUE, xy = TRUE)
bg_env <- bg_pts[, paste0("bio", bio_vars)]
cat(sprintf("  Background points: %d\n", nrow(bg_env)))

# Extract environment at occurrence points
extract_env <- function(xy_df) {
    pts <- vect(xy_df, geom = c(lon_col, lat_col),
                crs = "+proj=longlat +datum=WGS84")
    env <- extract(bio_crop, pts)[, -1]
    names(env) <- paste0("bio", bio_vars)
    env[complete.cases(env), ]
}
env_WH <- extract_env(xy_WH)
env_HM <- extract_env(xy_HM)
cat(sprintf("  WH env values: %d points\n", nrow(env_WH)))
cat(sprintf("  HM env values: %d points\n", nrow(env_HM)))

# --------------------------------------------------------------
# STEP 3: PCA-env (calibrated on background + occurrences)
# --------------------------------------------------------------
cat("\n-- STEP 3: PCA-env analysis --\n")
all_env <- rbind(bg_env, env_WH, env_HM)
pca_cal <- dudi.pca(all_env, center = TRUE, scale = TRUE, scannf = FALSE, nf = 2)

pc1_var <- round(pca_cal$eig[1] / sum(pca_cal$eig) * 100, 1)
pc2_var <- round(pca_cal$eig[2] / sum(pca_cal$eig) * 100, 1)
cat(sprintf("  PC1 = %.1f%%, PC2 = %.1f%%\n", pc1_var, pc2_var))

scores_bg <- pca_cal$li[1:nrow(bg_env), ]
scores_WH <- pca_cal$li[(nrow(bg_env) + 1):(nrow(bg_env) + nrow(env_WH)), ]
scores_HM <- pca_cal$li[(nrow(bg_env) + nrow(env_WH) + 1):nrow(all_env), ]

# Niche grids
z_WH <- ecospat.grid.clim.dyn(glob = scores_bg, glob1 = scores_bg,
                              sp = scores_WH, R = R_grid, th.sp = 0)
z_HM <- ecospat.grid.clim.dyn(glob = scores_bg, glob1 = scores_bg,
                              sp = scores_HM, R = R_grid, th.sp = 0)

# --------------------------------------------------------------
# STEP 4: Schoener's D / Hellinger's I overlap
# --------------------------------------------------------------
cat("\n-- STEP 4: Niche overlap (Schoener's D / Hellinger's I) --\n")
overlap <- ecospat.niche.overlap(z_WH, z_HM, cor = TRUE)
D_value <- overlap$D
I_value <- overlap$I
cat(sprintf("  D = %.4f, I = %.4f\n", D_value, I_value))
cat(sprintf("  Interpretation: %s\n",
            ifelse(D_value < 0.20, "Low - strong niche divergence",
            ifelse(D_value < 0.50, "Moderate - partial divergence",
                                   "High - niches largely similar"))))

# --------------------------------------------------------------
# STEP 5: Niche equivalency test
# --------------------------------------------------------------
cat(sprintf("\n-- STEP 5: Niche equivalency test (%d permutations) --\n", n_permutations))
equiv_test <- ecospat.niche.equivalency.test(
    z1 = z_WH, z2 = z_HM, rep = n_permutations,
    overlap.alternative = "lower", ncores = 1)
p_equiv_D <- equiv_test$p.D
p_equiv_I <- equiv_test$p.I
cat(sprintf("  p (D) = %.4f, p (I) = %.4f\n", p_equiv_D, p_equiv_I))
cat(sprintf("  Result: %s\n",
            ifelse(p_equiv_D < 0.05,
                   "SIGNIFICANT - niches are NOT equivalent",
                   "Not significant - cannot reject niche equivalency")))

# --------------------------------------------------------------
# STEP 6: Niche similarity tests (both directions)
# --------------------------------------------------------------
cat(sprintf("\n-- STEP 6: Niche similarity test (%d permutations, both directions) --\n", n_permutations))
sim_WH_to_HM <- ecospat.niche.similarity.test(z1 = z_WH, z2 = z_HM,
                                              rep = n_permutations,
                                              overlap.alternative = "lower",
                                              rand.type = 2, ncores = 1)
sim_HM_to_WH <- ecospat.niche.similarity.test(z1 = z_HM, z2 = z_WH,
                                              rep = n_permutations,
                                              overlap.alternative = "lower",
                                              rand.type = 2, ncores = 1)
p_sim_WH_HM_D <- sim_WH_to_HM$p.D
p_sim_HM_WH_D <- sim_HM_to_WH$p.D
cat(sprintf("  WH->HM (D): p = %.4f\n", p_sim_WH_HM_D))
cat(sprintf("  HM->WH (D): p = %.4f\n", p_sim_HM_WH_D))

# --------------------------------------------------------------
# STEP 7: Plots
# --------------------------------------------------------------
cat("\n-- STEP 7: Generating figures --\n")

# --- Figure 1: PCA-env niche overlap (two-panel) -- PDF ---
pdf("Figure_PCAenv_niche_overlap.pdf", width = 10, height = 5, onefile = FALSE)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0), family = "sans")
ecospat.plot.niche.dyn(
    z1 = z_WH, z2 = z_HM, quant = 0.10, interest = 1,
    title = "A  West Himalaya (WH) niche",
    name.axis1 = paste0("PC1 (", pc1_var, "% variance)"),
    name.axis2 = paste0("PC2 (", pc2_var, "% variance)"),
    col.unf = "#E8F5E9", col.exp = "#66BB6A", col.stab = "#1B5E20",
    col.bord.unf = "grey60", col.bord.exp = "#1B5E20", col.bord.stab = "#1B5E20"
)
ecospat.plot.niche.dyn(
    z1 = z_WH, z2 = z_HM, quant = 0.10, interest = 2,
    title = "B  Hengduan Mountains (HM) niche",
    name.axis1 = paste0("PC1 (", pc1_var, "% variance)"),
    name.axis2 = paste0("PC2 (", pc2_var, "% variance)"),
    col.unf = "#E3F2FD", col.exp = "#42A5F5", col.stab = "#0D47A1",
    col.bord.unf = "grey60", col.bord.exp = "#0D47A1", col.bord.stab = "#0D47A1"
)
mtext(bquote(italic("Rumex hastatus") ~
             "niche comparison | Schoener's D =" ~ .(round(D_value, 3))),
      outer = TRUE, cex = 1.1, font = 2)
dev.off()

# --- Figure 1: PNG version (300 DPI) ---
png("Figure_PCAenv_niche_overlap.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0), family = "sans")
ecospat.plot.niche.dyn(
    z1 = z_WH, z2 = z_HM, quant = 0.10, interest = 1,
    title = "A  West Himalaya (WH) niche",
    name.axis1 = paste0("PC1 (", pc1_var, "% variance)"),
    name.axis2 = paste0("PC2 (", pc2_var, "% variance)"),
    col.unf = "#E8F5E9", col.exp = "#66BB6A", col.stab = "#1B5E20",
    col.bord.unf = "grey60", col.bord.exp = "#1B5E20", col.bord.stab = "#1B5E20"
)
ecospat.plot.niche.dyn(
    z1 = z_WH, z2 = z_HM, quant = 0.10, interest = 2,
    title = "B  Hengduan Mountains (HM) niche",
    name.axis1 = paste0("PC1 (", pc1_var, "% variance)"),
    name.axis2 = paste0("PC2 (", pc2_var, "% variance)"),
    col.unf = "#E3F2FD", col.exp = "#42A5F5", col.stab = "#0D47A1",
    col.bord.unf = "grey60", col.bord.exp = "#0D47A1", col.bord.stab = "#0D47A1"
)
mtext(bquote(italic("Rumex hastatus") ~
             "niche comparison | Schoener's D =" ~ .(round(D_value, 3))),
      outer = TRUE, cex = 1.1, font = 2)
dev.off()

# --- Figure 2: Equivalency test histogram ---
pdf("Figure_equivalency_test.pdf", width = 6, height = 5, onefile = FALSE)
ecospat.plot.overlap.test(
    equiv_test, type = "D",
    title = paste0("Niche Equivalency Test (WH vs HM)\n",
                   "Observed D = ", round(D_value, 3),
                   " | p = ", round(p_equiv_D, 3))
)
dev.off()

# --- Figure 3: Similarity test histograms (both directions) ---
pdf("Figure_similarity_test.pdf", width = 10, height = 5, onefile = FALSE)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.5, 1.5), oma = c(0, 0, 2, 0))
ecospat.plot.overlap.test(sim_WH_to_HM, type = "D",
    title = paste0("Similarity: WH -> HM\np = ", round(p_sim_WH_HM_D, 3)))
ecospat.plot.overlap.test(sim_HM_to_WH, type = "D",
    title = paste0("Similarity: HM -> WH\np = ", round(p_sim_HM_WH_D, 3)))
mtext("Niche Similarity Tests | Rumex hastatus (WH vs HM)",
      outer = TRUE, cex = 1.1, font = 2)
dev.off()

# --- Figure 4: PCA biplot with variable loadings ---
pdf("Figure_PCA_biplot_loadings.pdf", width = 7, height = 6)
bg_df <- data.frame(scores_bg, group = "Background")
wh_df <- data.frame(scores_WH, group = "WH")
hm_df <- data.frame(scores_HM, group = "HM")
all_scores <- rbind(bg_df, wh_df, hm_df)

loadings <- as.data.frame(pca_cal$c1[, 1:2])
loadings$variable <- rownames(loadings)
names(loadings)[1:2] <- c("PC1", "PC2")
arrow_scale <- 2.5

plot(all_scores$Axis1, all_scores$Axis2,
     col = ifelse(all_scores$group == "WH", "#1B5E20",
                  ifelse(all_scores$group == "HM", "#0D47A1", "grey85")),
     pch = ifelse(all_scores$group == "Background", 1, 19),
     cex = ifelse(all_scores$group == "Background", 0.4, 1.2),
     xlab = paste0("PC1 (", pc1_var, "% variance)"),
     ylab = paste0("PC2 (", pc2_var, "% variance)"),
     main = "PCA-env Biplot - Environmental Drivers of Niche Separation",
     asp = 1)
arrows(0, 0, loadings$PC1 * arrow_scale, loadings$PC2 * arrow_scale,
       length = 0.08, col = "#C62828", lwd = 1.5)
text(loadings$PC1 * arrow_scale * 1.15, loadings$PC2 * arrow_scale * 1.15,
     labels = loadings$variable, col = "#C62828", cex = 0.85, font = 2)
legend("topright",
       legend = c("WH occurrences", "HM occurrences", "Background"),
       col = c("#1B5E20", "#0D47A1", "grey70"),
       pch = c(19, 19, 1), pt.cex = c(1.2, 1.2, 0.8), bty = "n", cex = 0.9)
abline(h = 0, v = 0, lty = 2, col = "grey60")
dev.off()

# --- Figure 5: Single overlaid niche plot ---
draw_overlaid_niches <- function() {
    xrange <- range(z_WH$x); yrange <- range(z_WH$y)
    plot(0, type = "n", xlim = xrange, ylim = yrange,
         xlab = paste0("PC1 (", pc1_var, "% variance)"),
         ylab = paste0("PC2 (", pc2_var, "% variance)"), main = "")

    wh_d <- z_WH$z.uncor; wh_d[wh_d == 0] <- NA
    image(z_WH$x, z_WH$y, wh_d,
          col = adjustcolor(colorRampPalette(c("#E8F5E9","#A5D6A7","#2E7D32"))(100), 0.85),
          add = TRUE)
    hm_d <- z_HM$z.uncor; hm_d[hm_d == 0] <- NA
    image(z_HM$x, z_HM$y, hm_d,
          col = adjustcolor(colorRampPalette(c("#E3F2FD","#90CAF9","#0D47A1"))(100), 0.70),
          add = TRUE)
    contour(z_WH$x, z_WH$y, z_WH$z.uncor,
            levels = quantile(z_WH$z.uncor[z_WH$z.uncor > 0], 0.10, na.rm = TRUE),
            col = "#1B5E20", lwd = 1.8, drawlabels = FALSE, add = TRUE)
    contour(z_HM$x, z_HM$y, z_HM$z.uncor,
            levels = quantile(z_HM$z.uncor[z_HM$z.uncor > 0], 0.10, na.rm = TRUE),
            col = "#0D47A1", lwd = 1.8, drawlabels = FALSE, add = TRUE)
    contour(z_WH$x, z_WH$y, z_WH$w, levels = 0.5,
            col = "#C0392B", lwd = 1.5, drawlabels = FALSE, add = TRUE)
    title(main = expression(paste(italic("Rumex hastatus"), " ecological niche comparison")),
          cex.main = 1.25, font.main = 1)
    legend("topright",
           legend = c("WH niche (West Himalaya)",
                      "HM niche (Hengduan Mts.)",
                      "Available climate space"),
           fill = c(adjustcolor("#2E7D32", 0.75), adjustcolor("#0D47A1", 0.70), NA),
           border = c("#1B5E20", "#0D47A1", "#C0392B"),
           lty = c(NA, NA, 1), col = c(NA, NA, "#C0392B"),
           bty = "n", cex = 0.9, pt.cex = 1.5)
    mtext(bquote("Schoener's" ~ italic(D) ~ "=" ~ .(round(D_value, 3)) ~
                 "  |  Equivalency" ~ italic(p) ~ "=" ~ .(round(p_equiv_D, 4))),
          side = 1, line = -1.8, adj = 0.02, cex = 0.82, col = "#555555", font = 3)
}
pdf("Figure_niche_overlap_SINGLE.pdf", width = 7, height = 6, onefile = FALSE)
par(mar = c(4.5, 4.5, 4, 2), family = "sans"); draw_overlaid_niches(); dev.off()
png("Figure_niche_overlap_SINGLE.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(4.5, 4.5, 4, 2), family = "sans"); draw_overlaid_niches(); dev.off()

# --------------------------------------------------------------
# STEP 8: Export results to CSV
# --------------------------------------------------------------
cat("\n-- STEP 8: Saving CSV outputs --\n")
results_df <- data.frame(
    Analysis = c("PCA PC1 variance (%)","PCA PC2 variance (%)",
                 "Schoener's D","Hellinger's I",
                 "Equivalency test p-value (D)","Equivalency test p-value (I)",
                 "Similarity test p-value WH->HM (D)",
                 "Similarity test p-value HM->WH (D)",
                 "WH occurrences used","HM occurrences used",
                 "Background points","Permutations","Split longitude"),
    Value = c(pc1_var, pc2_var,
              round(D_value, 4), round(I_value, 4),
              round(p_equiv_D, 4), round(p_equiv_I, 4),
              round(p_sim_WH_HM_D, 4), round(p_sim_HM_WH_D, 4),
              nrow(env_WH), nrow(env_HM), nrow(bg_env), n_permutations, split_lon),
    Interpretation = c("", "",
                       ifelse(D_value < 0.20, "Low - strong divergence",
                       ifelse(D_value < 0.50, "Moderate - partial divergence",
                              "High - conserved")), "",
                       ifelse(p_equiv_D < 0.05, "SIGNIFICANT - not equivalent", "Not significant"),
                       ifelse(p_equiv_I < 0.05, "SIGNIFICANT", "Not significant"),
                       ifelse(p_sim_WH_HM_D < 0.05, "SIGNIFICANT - divergence > random", "Not significant"),
                       ifelse(p_sim_HM_WH_D < 0.05, "SIGNIFICANT - divergence > random", "Not significant"),
                       "", "", "", "", "")
)
write.csv(results_df, "niche_overlap_results.csv", row.names = FALSE)
write.csv(data.frame(scores_bg, group = "Background"),
          "niche_background_env.csv", row.names = FALSE)
write.csv(data.frame(scores_WH, group = "WH"), "niche_WH_scores.csv", row.names = FALSE)
write.csv(data.frame(scores_HM, group = "HM"), "niche_HM_scores.csv", row.names = FALSE)

# --------------------------------------------------------------
# Final summary
# --------------------------------------------------------------
cat("\n", strrep("=", 65), "\n", sep = "")
cat("  NICHE DIVERGENCE ANALYSIS COMPLETE\n")
cat(strrep("=", 65), "\n")
cat(sprintf("  Schoener's D:        %.4f\n", D_value))
cat(sprintf("  Equivalency p (D):   %.4f  %s\n", p_equiv_D,
            ifelse(p_equiv_D < 0.05, "[SIGNIFICANT]", "[not significant]")))
cat(sprintf("  Similarity WH->HM:   %.4f  %s\n", p_sim_WH_HM_D,
            ifelse(p_sim_WH_HM_D < 0.05, "[SIGNIFICANT]", "[not significant]")))
cat(sprintf("  Similarity HM->WH:   %.4f  %s\n", p_sim_HM_WH_D,
            ifelse(p_sim_HM_WH_D < 0.05, "[SIGNIFICANT]", "[not significant]")))
cat("  Outputs in:", wd, "\n")
cat(strrep("=", 65), "\n")
