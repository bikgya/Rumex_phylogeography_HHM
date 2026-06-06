# ==============================================================
# Script 09: Dispersal Corridor Analysis (Waypoint-Guided LCP)
# Project: rumex-landscape-phylogenomics
#
# Tests for a Late Pleistocene (LGM) dispersal corridor between
# the West Himalaya (WH) and the Hengduan Mountains (HM) using
# the LGM ensemble suitability rasters produced by Script 08c.
#
# Pipeline
#   1. Load 3 LGM ensemble projections (CCSM4 / MIROC / MPI),
#      auto-rescale to 0-1, and build a multi-GCM mean
#   2. Extract TSS thresholds from the ensemble eval files;
#      build a binary consensus map (>= 2 of 3 GCMs agree)
#   3. Cluster occurrences along the WH-HM arc into waypoints
#      so the least-cost path follows the biologically
#      plausible southern Himalayan route (not a straight line
#      across the Tibetan Plateau)
#   4. Compute LCP through waypoints using gdistance
#      (also computes the direct LCP for comparison)
#   5. Extract suitability transect along the LCP, classify
#      each point as Corridor / Marginal / Barrier
#   6. Save figures (transect, LGM map with LCP, comparison)
#      and CSV summaries
#
# Input (from Script 08c)
#   - 3 LGM ensemble projection rasters (.tif), one per GCM
#   - Matching ensemble evaluation CSVs (eval_EM.csv)
#   - Whole-region occurrence CSV (Longitude, Latitude)
#
# Output (under out_dir)
#   Data    : LGM_suitability_mean.tif, LGM_binary_consensus.tif,
#             corridor_transect_along_LCP.csv,
#             corridor_connectivity_summary.csv
#   Figures : Figure_A_transect_along_LCP.{png,pdf}
#             Figure_B_LGM_map_with_LCP.{png,pdf}
#             Figure_C_comparison_straight_vs_LCP.{png,pdf}
#             Figure_composite_2panel.{png,pdf}  (A+B)
#             Figure_composite_3panel.{png,pdf}  (A+B+C)
#
# Dependencies: terra, gdistance, raster, geosphere, ggplot2,
#               dplyr, scales, tidyterra, sp, gridExtra, grid, cowplot
# ==============================================================

# ----- Load packages (assumes already installed; see envs/) -----
required_pkgs <- c("terra","gdistance","raster","geosphere","ggplot2",
                   "dplyr","scales","tidyterra","sp","gridExtra","grid","cowplot")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Install missing packages first: ",
                                   paste(missing_pkgs, collapse = ", "))
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ==============================================================
# USER CONFIGURATION
# ==============================================================
# Working directory (where outputs go)
wd <- "results/corridor"

# LGM ensemble projections from Script 08c (one per GCM)
raster_paths <- list(
    CCSM4 = "results/sdm/env/LGM_LGM_CCSM4_autocorrelated/biomod_output/LGM_CCSM4_EM/individual_projections/rumex_EMmeanByTSS_mergedData_mergedRun_mergedAlgo.tif",
    MIROC = "results/sdm/env/LGM_LGM_MIROC_autocorrelated/biomod_output/LGM_MIROC_EM/individual_projections/rumex_EMmeanByTSS_mergedData_mergedRun_mergedAlgo.tif",
    MPI   = "results/sdm/env/LGM_LGM_MPI_autocorrelated/biomod_output/LGM_MPI_EM/individual_projections/rumex_EMmeanByTSS_mergedData_mergedRun_mergedAlgo.tif"
)

# Ensemble eval CSVs from Script 08c (one per GCM)
eval_paths <- list(
    CCSM4 = "results/sdm/env/LGM_LGM_CCSM4_autocorrelated/biomod_output/rumex_LGM_CCSM4_eval_EM.csv",
    MIROC = "results/sdm/env/LGM_LGM_MIROC_autocorrelated/biomod_output/rumex_LGM_MIROC_eval_EM.csv",
    MPI   = "results/sdm/env/LGM_LGM_MPI_autocorrelated/biomod_output/rumex_LGM_MPI_eval_EM.csv"
)

# Whole-region occurrence CSV
occ_csv <- "data/occurrences/rumex_whole.csv"
lon_col <- "Longitude"
lat_col <- "Latitude"

# Analysis parameters
barrier_threshold   <- 0.30   # below = "Barrier"
corridor_threshold  <- 0.60   # above = "Corridor"; in-between = "Marginal"
consensus_min_gcms  <- 2      # >= N GCMs must agree for binary consensus

# LCP waypoint clustering (higher = more constrained to occurrence arc)
n_waypoint_clusters <- 8      # 6-10 is a reasonable range

# Number of evenly spaced points along the LCP for the transect
n_lcp_sample_points <- 200

# Study area bounding box (for the map figure)
study_lon_min <- 69.25;  study_lon_max <- 104.10
study_lat_min <- 26.22;  study_lat_max <- 35.25
# ==============================================================

out_dir <- file.path(wd, "outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n")

# --------------------------------------------------------------
# STEP 1: Load LGM rasters, inspect range, auto-rescale to 0-1
# --------------------------------------------------------------
cat("\n-- STEP 1: Loading LGM rasters --\n")
load_and_rescale <- function(path, gcm_name) {
    if (!file.exists(path)) stop("Raster not found for ", gcm_name, ": ", path)
    r <- rast(path)
    if (nlyr(r) > 1) r <- r[[1]]
    raw_max <- global(r, "max", na.rm = TRUE)[[1]]
    if (raw_max > 1) {
        r <- r / 1000
        cat(sprintf("  [%s] Rescaled to 0-1\n", gcm_name))
    }
    names(r) <- gcm_name
    r
}
suit_list <- mapply(load_and_rescale, raster_paths, names(raster_paths), SIMPLIFY = FALSE)

# Align to CCSM4 grid
ref <- suit_list[["CCSM4"]]
suit_list <- lapply(names(suit_list), function(nm) {
    r <- suit_list[[nm]]
    if (!compareGeom(ref, r, stopOnError = FALSE)) {
        r <- resample(r, ref, method = "bilinear"); names(r) <- nm
    }
    r
})
names(suit_list) <- c("CCSM4", "MIROC", "MPI")
suit_stack    <- rast(suit_list)
suit_lgm_mean <- mean(suit_stack, na.rm = TRUE)
names(suit_lgm_mean) <- "mean_suitability"

# --------------------------------------------------------------
# STEP 2: Extract TSS thresholds from eval CSVs
# --------------------------------------------------------------
cat("\n-- STEP 2: Extracting TSS thresholds --\n")
extract_tss_threshold <- function(eval_path, gcm_name) {
    if (!file.exists(eval_path)) {
        warning("[", gcm_name, "] eval file not found - using default 0.40")
        return(0.40)
    }
    tss_val <- tryCatch({
        eval_df <- read.csv(eval_path, stringsAsFactors = FALSE)
        metric_col <- grep("metric", names(eval_df), ignore.case = TRUE, value = TRUE)[1]
        cutoff_col <- grep("cutoff", names(eval_df), ignore.case = TRUE, value = TRUE)[1]
        if (is.na(metric_col) || is.na(cutoff_col))
            stop("Could not find 'metric' or 'cutoff' columns")
        tss_rows <- eval_df[grepl("TSS", eval_df[[metric_col]], ignore.case = TRUE), ]
        if (nrow(tss_rows) == 0) stop("No TSS rows found")
        thresh_raw <- mean(as.numeric(tss_rows[[cutoff_col]]), na.rm = TRUE)
        if (thresh_raw > 1) thresh_raw / 1000 else thresh_raw
    }, error = function(e) {
        warning("[", gcm_name, "] ", e$message, " - using 0.40")
        0.40
    })
    cat(sprintf("  [%s] TSS threshold = %.4f\n", gcm_name, tss_val))
    tss_val
}
tss_thresholds <- mapply(extract_tss_threshold, eval_paths, names(eval_paths))
tss_mean <- mean(tss_thresholds)

# --------------------------------------------------------------
# STEP 3: Binary consensus map (>= N GCMs agree)
# --------------------------------------------------------------
cat("\n-- STEP 3: Building binary consensus map --\n")
bin_ccsm4 <- ifel(suit_list[["CCSM4"]] >= tss_thresholds["CCSM4"], 1, 0)
bin_miroc <- ifel(suit_list[["MIROC"]] >= tss_thresholds["MIROC"], 1, 0)
bin_mpi   <- ifel(suit_list[["MPI"]]   >= tss_thresholds["MPI"],   1, 0)
bin_consensus <- ifel((bin_ccsm4 + bin_miroc + bin_mpi) >= consensus_min_gcms, 1, 0)
names(bin_consensus) <- "binary_consensus"

writeRaster(suit_lgm_mean, file.path(out_dir, "LGM_suitability_mean.tif"), overwrite = TRUE)
writeRaster(bin_consensus, file.path(out_dir, "LGM_binary_consensus.tif"), overwrite = TRUE)

# --------------------------------------------------------------
# STEP 4: Load occurrences and define endpoints
# --------------------------------------------------------------
cat("\n-- STEP 4: Loading occurrences and defining endpoints --\n")
if (!file.exists(occ_csv)) stop("CSV not found: ", occ_csv)
occ <- read.csv(occ_csv, stringsAsFactors = FALSE)
occ <- occ[!is.na(occ[[lon_col]]) & !is.na(occ[[lat_col]]), ]
cat(sprintf("  Occurrence records: %d\n", nrow(occ)))

wh_idx <- which.min(occ[[lon_col]])
hm_idx <- which.max(occ[[lon_col]])
wh_coords <- c(lon = occ[[lon_col]][wh_idx], lat = occ[[lat_col]][wh_idx])
hm_coords <- c(lon = occ[[lon_col]][hm_idx], lat = occ[[lat_col]][hm_idx])
total_dist_km <- distHaversine(c(wh_coords["lon"], wh_coords["lat"]),
                               c(hm_coords["lon"], hm_coords["lat"])) / 1000
cat(sprintf("  WH endpoint: %.4f E, %.4f N\n", wh_coords["lon"], wh_coords["lat"]))
cat(sprintf("  HM endpoint: %.4f E, %.4f N\n", hm_coords["lon"], hm_coords["lat"]))
cat(sprintf("  Straight-line distance: %.0f km\n", total_dist_km))

# --------------------------------------------------------------
# STEP 5: Build waypoints by k-means clustering on longitude
# --------------------------------------------------------------
cat("\n-- STEP 5: Building waypoints --\n")
set.seed(42)
occ_sorted <- occ[order(occ[[lon_col]]), ]
occ_coords <- data.frame(lon = occ_sorted[[lon_col]], lat = occ_sorted[[lat_col]])
km_fit <- kmeans(occ_coords$lon, centers = n_waypoint_clusters, nstart = 25)
occ_coords$cluster <- km_fit$cluster

waypoints <- occ_coords %>%
    group_by(cluster) %>%
    summarise(lon = median(lon), lat = median(lat), n_pts = n(), .groups = "drop") %>%
    arrange(lon)

all_waypoints <- rbind(
    data.frame(lon = wh_coords["lon"], lat = wh_coords["lat"]),
    data.frame(lon = waypoints$lon,    lat = waypoints$lat),
    data.frame(lon = hm_coords["lon"], lat = hm_coords["lat"])
)
rownames(all_waypoints) <- NULL
cat(sprintf("  Total waypoint sequence: %d (WH + %d clusters + HM)\n",
            nrow(all_waypoints), nrow(waypoints)))

# --------------------------------------------------------------
# STEP 6: Waypoint-guided least-cost path
# --------------------------------------------------------------
cat("\n-- STEP 6: Computing waypoint-guided LCP --\n")
suit_rl <- raster::raster(suit_lgm_mean)
tr <- gdistance::transition(suit_rl, transitionFunction = mean, directions = 8)
tr <- gdistance::geoCorrection(tr, type = "c", multpl = FALSE)

lcp_all_coords <- NULL
for (i in 1:(nrow(all_waypoints) - 1)) {
    from_sp <- sp::SpatialPoints(matrix(c(all_waypoints$lon[i], all_waypoints$lat[i]), nrow = 1),
                                 proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
    to_sp   <- sp::SpatialPoints(matrix(c(all_waypoints$lon[i+1], all_waypoints$lat[i+1]), nrow = 1),
                                 proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
    seg_lcp <- tryCatch(gdistance::shortestPath(tr, from_sp, to_sp, output = "SpatialLines"),
                        error = function(e) NULL)
    seg_coords <- if (!is.null(seg_lcp)) {
        sp::coordinates(seg_lcp)[[1]][[1]]
    } else {
        # Fallback: straight line for this segment
        matrix(c(all_waypoints$lon[i], all_waypoints$lat[i],
                 all_waypoints$lon[i+1], all_waypoints$lat[i+1]),
               ncol = 2, byrow = TRUE)
    }
    lcp_all_coords <- rbind(lcp_all_coords, seg_coords)
}

lcp_all_coords <- as.data.frame(lcp_all_coords)
names(lcp_all_coords) <- c("lon", "lat")
dup_idx <- c(FALSE, (diff(lcp_all_coords$lon) == 0) & (diff(lcp_all_coords$lat) == 0))
lcp_all_coords <- lcp_all_coords[!dup_idx, ]
cat(sprintf("  Full LCP: %d vertices\n", nrow(lcp_all_coords)))

# Direct LCP (no waypoints) for comparison
wh_sp <- sp::SpatialPoints(matrix(c(wh_coords["lon"], wh_coords["lat"]), nrow = 1),
                           proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
hm_sp <- sp::SpatialPoints(matrix(c(hm_coords["lon"], hm_coords["lat"]), nrow = 1),
                           proj4string = sp::CRS("+proj=longlat +datum=WGS84"))
lcp_direct <- gdistance::shortestPath(tr, wh_sp, hm_sp, output = "SpatialLines")
direct_coords <- as.data.frame(sp::coordinates(lcp_direct)[[1]][[1]])
names(direct_coords) <- c("lon", "lat")
cost_direct <- gdistance::costDistance(tr, wh_sp, hm_sp)[1, 1]

# --------------------------------------------------------------
# STEP 7: Extract suitability transect along the LCP
# --------------------------------------------------------------
cat("\n-- STEP 7: Extracting suitability along LCP --\n")
n_lcp_total <- nrow(lcp_all_coords)
sample_idx  <- if (n_lcp_total > n_lcp_sample_points)
    round(seq(1, n_lcp_total, length.out = n_lcp_sample_points)) else seq_len(n_lcp_total)
lcp_sample <- lcp_all_coords[sample_idx, ]

lcp_dists <- c(0, cumsum(
    distHaversine(as.matrix(lcp_sample[-nrow(lcp_sample), ]),
                  as.matrix(lcp_sample[-1, ])) / 1000))

lcp_sv <- vect(as.matrix(lcp_sample), crs = "+proj=longlat +datum=WGS84")
suit_mean_vals  <- extract(suit_lgm_mean,         lcp_sv)[, 2]
suit_ccsm4_vals <- extract(suit_list[["CCSM4"]], lcp_sv)[, 2]
suit_miroc_vals <- extract(suit_list[["MIROC"]], lcp_sv)[, 2]
suit_mpi_vals   <- extract(suit_list[["MPI"]],   lcp_sv)[, 2]
bin_vals        <- extract(bin_consensus,         lcp_sv)[, 2]

suit_mat <- cbind(suit_ccsm4_vals, suit_miroc_vals, suit_mpi_vals)
suit_sd  <- apply(suit_mat, 1, sd, na.rm = TRUE)
suit_lo  <- pmax(suit_mean_vals - suit_sd, 0)
suit_hi  <- pmin(suit_mean_vals + suit_sd, 1)

zone_class <- dplyr::case_when(
    suit_mean_vals <  barrier_threshold  ~ "Barrier",
    suit_mean_vals >= corridor_threshold ~ "Corridor",
    TRUE                                 ~ "Marginal"
)

transect_df <- data.frame(
    point_id = seq_along(sample_idx),
    longitude = lcp_sample$lon, latitude = lcp_sample$lat,
    distance_km = round(lcp_dists, 1),
    suit_mean = round(suit_mean_vals, 4),
    suit_CCSM4 = round(suit_ccsm4_vals, 4),
    suit_MIROC = round(suit_miroc_vals, 4),
    suit_MPI   = round(suit_mpi_vals, 4),
    suit_SD = round(suit_sd, 4),
    suit_lo = round(suit_lo, 4), suit_hi = round(suit_hi, 4),
    binary_consensus = bin_vals, habitat_zone = zone_class
)
total_lcp_dist_km <- max(transect_df$distance_km)
mean_lcp_suit <- mean(transect_df$suit_mean, na.rm = TRUE)
min_lcp_suit  <- min(transect_df$suit_mean, na.rm = TRUE)

write.csv(transect_df, file.path(out_dir, "corridor_transect_along_LCP.csv"), row.names = FALSE)

# --------------------------------------------------------------
# STEP 8: Connectivity summary CSV
# --------------------------------------------------------------
connectivity_df <- data.frame(
    metric = c("WH_longitude","WH_latitude","HM_longitude","HM_latitude",
               "straight_line_distance_km","LCP_distance_km","LCP_type",
               "n_waypoint_clusters","direct_LCP_cost_distance",
               "mean_suitability_along_LCP","min_suitability_along_LCP",
               "pct_corridor_points","pct_marginal_points","pct_barrier_points",
               "mean_TSS_threshold","corridor_assessment"),
    value = c(round(wh_coords["lon"], 4), round(wh_coords["lat"], 4),
              round(hm_coords["lon"], 4), round(hm_coords["lat"], 4),
              round(total_dist_km, 1), round(total_lcp_dist_km, 1),
              "waypoint-guided (occurrence clusters)", n_waypoint_clusters,
              round(cost_direct, 2), round(mean_lcp_suit, 4), round(min_lcp_suit, 4),
              round(100 * mean(zone_class == "Corridor", na.rm = TRUE), 1),
              round(100 * mean(zone_class == "Marginal", na.rm = TRUE), 1),
              round(100 * mean(zone_class == "Barrier",  na.rm = TRUE), 1),
              round(tss_mean, 4),
              ifelse(min_lcp_suit >= corridor_threshold, "CONTINUOUS CORRIDOR",
              ifelse(min_lcp_suit >= barrier_threshold, "MARGINAL CORRIDOR (stepping-stone)",
                                                       "BARRIER DETECTED")))
)
write.csv(connectivity_df, file.path(out_dir, "corridor_connectivity_summary.csv"), row.names = FALSE)

# --------------------------------------------------------------
# STEP 9: Figure A - Suitability transect along the LCP
# --------------------------------------------------------------
cat("\n-- STEP 9: Generating figures --\n")
zone_colors <- c("Corridor" = "#2E7D32", "Marginal" = "#F9A825", "Barrier" = "#C62828")

wp_on_transect <- sapply(seq_len(nrow(all_waypoints)), function(i) {
    d <- distHaversine(c(all_waypoints$lon[i], all_waypoints$lat[i]),
                       as.matrix(transect_df[, c("longitude", "latitude")]))
    transect_df$distance_km[which.min(d)]
})
wp_labels <- c("WH", paste0("W", seq_len(nrow(waypoints))), "HM")

p_transect <- ggplot(transect_df, aes(x = distance_km, y = suit_mean)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = barrier_threshold,
             fill = "#FFCDD2", alpha = 0.45) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = barrier_threshold, ymax = corridor_threshold,
             fill = "#FFF9C4", alpha = 0.45) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = corridor_threshold, ymax = 1,
             fill = "#C8E6C9", alpha = 0.45) +
    geom_ribbon(aes(ymin = suit_lo, ymax = suit_hi), fill = "#90CAF9", alpha = 0.40) +
    geom_hline(yintercept = barrier_threshold,  linetype = "dashed", color = "#C62828", linewidth = 0.8) +
    geom_hline(yintercept = corridor_threshold, linetype = "dashed", color = "#2E7D32", linewidth = 0.8) +
    geom_line(aes(y = suit_CCSM4), color = "#78909C", linewidth = 0.5, alpha = 0.6) +
    geom_line(aes(y = suit_MIROC), color = "#8D6E63", linewidth = 0.5, alpha = 0.6) +
    geom_line(aes(y = suit_MPI),   color = "#7B1FA2", linewidth = 0.5, alpha = 0.6) +
    geom_line(color = "#1565C0", linewidth = 1.2) +
    geom_point(aes(color = habitat_zone), size = 2.2, shape = 19) +
    scale_color_manual(values = zone_colors, name = "Habitat Zone") +
    geom_vline(xintercept = wp_on_transect, linetype = "dotted", color = "grey50",
               linewidth = 0.4, alpha = 0.6) +
    annotate("text", x = wp_on_transect, y = rep(0.98, length(wp_on_transect)),
             label = wp_labels, size = 2.5, angle = 45, hjust = 0,
             color = "#1A237E", fontface = "bold") +
    scale_x_continuous(breaks = seq(0, ceiling(max(transect_df$distance_km)/500)*500, 500),
                       expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(title = "LGM Habitat Suitability Along Waypoint-Guided Least-Cost Path",
         subtitle = paste0("Route through ", n_waypoint_clusters, " waypoints | LCP distance: ",
                           round(total_lcp_dist_km), " km | Mean TSS: ", sprintf("%.3f", tss_mean)),
         x = "Distance Along LCP from West Himalaya (km)",
         y = "Ensemble Habitat Suitability (0-1)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(color = "grey40", size = 9),
          legend.position = c(0.88, 0.85),
          legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
          panel.grid.minor = element_blank())

# --------------------------------------------------------------
# STEP 10: Figure B - LGM map with LCP overlay
# --------------------------------------------------------------
suit_df <- as.data.frame(suit_lgm_mean, xy = TRUE)
names(suit_df)[3] <- "suitability"
suit_df <- suit_df[!is.na(suit_df$suitability), ]

wp_df <- all_waypoints
wp_df$label <- c("WH", paste0("W", seq_len(nrow(waypoints))), "HM")
wp_df$type  <- c("endpoint", rep("waypoint", nrow(waypoints)), "endpoint")

p_map <- ggplot() +
    geom_raster(data = suit_df, aes(x = x, y = y, fill = suitability)) +
    scale_fill_gradientn(colours = c("#F5F5F5","#FFF176","#FFB300","#E65100","#1B5E20"),
                         values  = scales::rescale(c(0, 0.2, 0.4, 0.6, 1)),
                         limits  = c(0, 1), name = "Habitat\nSuitability",
                         guide   = guide_colorbar(barwidth = 1, barheight = 8)) +
    geom_point(data = occ, aes(x = .data[[lon_col]], y = .data[[lat_col]]),
               color = "cyan", shape = 3, size = 1.2, alpha = 0.7) +
    geom_path(data = direct_coords, aes(x = lon, y = lat),
              color = "grey80", linewidth = 0.8, linetype = "dashed") +
    geom_path(data = lcp_all_coords, aes(x = lon, y = lat), color = "white", linewidth = 1.4) +
    geom_path(data = lcp_all_coords, aes(x = lon, y = lat), color = "#1565C0", linewidth = 0.8) +
    geom_point(data = wp_df[wp_df$type == "waypoint", ], aes(x = lon, y = lat),
               color = "white", fill = "#FFA726", shape = 21, size = 3, stroke = 1.0) +
    geom_point(data = wp_df[wp_df$type == "endpoint", ], aes(x = lon, y = lat),
               color = "white", fill = "#1565C0", shape = 21, size = 4.5, stroke = 1.5) +
    geom_text(data = wp_df[wp_df$type == "endpoint", ], aes(x = lon, y = lat, label = label),
              color = "white", fontface = "bold", nudge_y = 0.8, size = 3.5) +
    geom_text(data = wp_df[wp_df$type == "waypoint", ], aes(x = lon, y = lat, label = label),
              color = "white", fontface = "bold", nudge_y = 0.5, size = 2.5) +
    coord_fixed(xlim = c(study_lon_min - 1, study_lon_max + 1),
                ylim = c(study_lat_min - 2, study_lat_max + 1)) +
    labs(title    = expression(paste("LGM Ensemble Habitat Suitability - ", italic("Rumex hastatus"))),
         subtitle = "Blue line = waypoint-guided LCP | Grey dashed = direct LCP | Orange dots = waypoints",
         x = "Longitude (E)", y = "Latitude (N)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(color = "grey40", size = 9),
          legend.position = "right")

# --------------------------------------------------------------
# STEP 11: Figure C - Straight-line vs LCP comparison
# --------------------------------------------------------------
n_straight <- 200
frac_seq   <- seq(0, 1, length.out = n_straight)
straight_pts <- do.call(rbind, lapply(frac_seq, function(f) {
    c(lon = wh_coords["lon"] + f * (hm_coords["lon"] - wh_coords["lon"]),
      lat = wh_coords["lat"] + f * (hm_coords["lat"] - wh_coords["lat"]))
}))
straight_dist_km <- c(0, cumsum(distHaversine(straight_pts[-nrow(straight_pts), , drop = FALSE],
                                              straight_pts[-1, , drop = FALSE]) / 1000))
straight_sv   <- vect(straight_pts, crs = "+proj=longlat +datum=WGS84")
straight_suit <- extract(suit_lgm_mean, straight_sv)[, 2]

compare_df <- rbind(
    data.frame(distance_km = straight_dist_km, suitability = straight_suit,
               method = "Straight line (WH->HM)"),
    data.frame(distance_km = transect_df$distance_km, suitability = transect_df$suit_mean,
               method = "Waypoint-guided LCP")
)

p_compare <- ggplot(compare_df, aes(x = distance_km, y = suitability, color = method)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = barrier_threshold,
             fill = "#FFCDD2", alpha = 0.35) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = barrier_threshold, ymax = corridor_threshold,
             fill = "#FFF9C4", alpha = 0.35) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = corridor_threshold, ymax = 1,
             fill = "#C8E6C9", alpha = 0.35) +
    geom_hline(yintercept = barrier_threshold,  linetype = "dashed", color = "#C62828", linewidth = 0.6) +
    geom_hline(yintercept = corridor_threshold, linetype = "dashed", color = "#2E7D32", linewidth = 0.6) +
    geom_line(linewidth = 1.1) +
    scale_color_manual(values = c("Straight line (WH->HM)" = "#B0BEC5",
                                  "Waypoint-guided LCP"     = "#1565C0"),
                       name = "Transect Method") +
    scale_x_continuous(breaks = seq(0, 5000, 500), expand = c(0.01, 0)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    labs(title = "Comparison: Straight-Line vs Waypoint-Guided LCP Transect",
         subtitle = paste0("Straight-line: ", round(total_dist_km), " km | LCP: ",
                           round(total_lcp_dist_km), " km"),
         x = "Distance from West Himalaya (km)", y = "Mean Habitat Suitability (0-1)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(color = "grey40", size = 9),
          legend.position = c(0.75, 0.85),
          legend.background = element_rect(fill = "white", color = "grey80", linewidth = 0.5),
          panel.grid.minor = element_blank())

# --------------------------------------------------------------
# STEP 12: Save figures (individual + composites)
# --------------------------------------------------------------
ggsave(file.path(out_dir, "Figure_A_transect_along_LCP.png"),
       p_transect, width = 22/2.54, height = 11/2.54, units = "in", dpi = 300)
ggsave(file.path(out_dir, "Figure_B_LGM_map_with_LCP.png"),
       p_map, width = 22/2.54, height = 11/2.54, units = "in", dpi = 300)
ggsave(file.path(out_dir, "Figure_C_comparison_straight_vs_LCP.png"),
       p_compare, width = 22/2.54, height = 11/2.54, units = "in", dpi = 300)
ggsave(file.path(out_dir, "Figure_A_transect_along_LCP.pdf"),
       p_transect, width = 20/2.54, height = 10/2.54, units = "in")
ggsave(file.path(out_dir, "Figure_B_LGM_map_with_LCP.pdf"),
       p_map, width = 20/2.54, height = 10/2.54, units = "in")
ggsave(file.path(out_dir, "Figure_C_comparison_straight_vs_LCP.pdf"),
       p_compare, width = 20/2.54, height = 10/2.54, units = "in")

# 2-panel composite (A+B)
p_transect_l <- p_transect + labs(tag = "A") + theme(plot.tag = element_text(face = "bold", size = 16))
p_map_l      <- p_map      + labs(tag = "B") + theme(plot.tag = element_text(face = "bold", size = 16))
p_compare_l  <- p_compare  + labs(tag = "C") + theme(plot.tag = element_text(face = "bold", size = 16))

composite_2 <- cowplot::plot_grid(p_transect_l, p_map_l, ncol = 1,
                                  rel_heights = c(0.48, 0.52), align = "v")
cowplot::ggsave2(file.path(out_dir, "Figure_composite_2panel.png"),
                 composite_2, width = 11, height = 11.5, dpi = 300)
cowplot::ggsave2(file.path(out_dir, "Figure_composite_2panel.pdf"),
                 composite_2, width = 11, height = 11.5)

# 3-panel composite (A+B+C)
composite_3 <- cowplot::plot_grid(p_transect_l, p_map_l, p_compare_l, ncol = 1,
                                  rel_heights = c(0.32, 0.38, 0.30), align = "v")
cowplot::ggsave2(file.path(out_dir, "Figure_composite_3panel.png"),
                 composite_3, width = 11, height = 16, dpi = 300)
cowplot::ggsave2(file.path(out_dir, "Figure_composite_3panel.pdf"),
                 composite_3, width = 11, height = 16)

# --------------------------------------------------------------
# Final summary
# --------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("  CORRIDOR ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")
cat(sprintf("  Straight-line:    %.0f km\n", total_dist_km))
cat(sprintf("  LCP distance:     %.0f km\n", total_lcp_dist_km))
cat(sprintf("  Min suitability:  %.4f\n",   min_lcp_suit))
cat(sprintf("  Mean suitability: %.4f\n",   mean_lcp_suit))
cat("\n  Interpretation:\n")
if (min_lcp_suit >= corridor_threshold) {
    cat("  CONTINUOUS CORRIDOR - supports WH->HM dispersal via southern Himalayan arc.\n")
} else if (min_lcp_suit >= barrier_threshold) {
    cat("  MARGINAL / STEPPING-STONE CORRIDOR (min suitability =",
        round(min_lcp_suit, 3), ")\n")
} else {
    cat("  BARRIER DETECTED on LCP (min suitability =",
        round(min_lcp_suit, 3), ")\n")
}
cat("  Outputs in:", out_dir, "\n")
cat(strrep("=", 70), "\n")
