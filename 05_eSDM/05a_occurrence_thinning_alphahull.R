# ==============================================================
# Script 08a: Occurrence Cleaning, Spatial Thinning, and Alpha Hull
# Project: rumex-landscape-phylogenomics
#
# Prepares occurrence data for ensemble SDM (Script 08c):
#   1. Reads cleaned occurrence records
#   2. Spatially thins records to a minimum inter-point distance
#      (default 10 km) to reduce sampling bias
#   3. Delineates an alpha hull around the thinned records to
#      define the study/background area for SDM training
#
# Input
#   - occ_csv  : CSV with at least "Longitude" and "Latitude" columns
#                (one record per row)
#
# Output
#   - thinned_csv  : CSV of spatially thinned records
#   - hull_shp     : Alpha-hull polygon shapefile (study area)
#
# Dependencies: spThin, alphahull, terra, sf, rangeBuilder, ggplot2, maps
# ==============================================================

# ----- USER CONFIGURATION -----
occ_csv        <- "data/occurrences/rumex_raw.csv"
thinned_csv    <- "results/sdm/occurrences/rumex_thinned.csv"
hull_shp       <- "results/sdm/study_area/rumex_alphahull.shp"
point_shp      <- "results/sdm/study_area/rumex_points.shp"
species_label  <- "rumex_hastatus"
thin_km        <- 10        # minimum distance (km) between points
thin_reps      <- 10        # spThin replicates
hull_buffer_m  <- 200000    # buffer for alpha hull (200 km)
plot_xlim      <- c(70, 105)
plot_ylim      <- c(20, 42)
# ------------------------------

# Libraries
library(spThin)
library(geosphere)
library(alphahull)
library(terra)
library(sf)
library(rangeBuilder)
library(ggplot2)
library(maps)

# Create output directories
dir.create(dirname(thinned_csv), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(hull_shp),    recursive = TRUE, showWarnings = FALSE)

# --------------------------------------------------------------
# STEP 1: Load and clean occurrence data
# --------------------------------------------------------------
data <- read.csv(occ_csv, encoding = "UTF-8")

# Coerce coordinate columns to numeric and drop incomplete rows
data$Longitude <- as.numeric(data$Longitude)
data$Latitude  <- as.numeric(data$Latitude)
data <- na.omit(data)
data <- unique(data[, c("Longitude", "Latitude")])

# OPTIONAL diagnostic: check minimum pairwise distance.
# If already > thin_km, you can skip spatial thinning entirely.
dist_matrix <- distm(data[, c("Longitude", "Latitude")])
dist_matrix[dist_matrix == 0] <- NA
min_dist_km <- min(dist_matrix, na.rm = TRUE) / 1000
message("Minimum pairwise distance (km): ", round(min_dist_km, 3))

# --------------------------------------------------------------
# STEP 2: Spatial thinning (spThin)
# --------------------------------------------------------------
thin_input <- data.frame(
    species            = species_label,
    decimalLatitude    = data$Latitude,
    decimalLongitude   = data$Longitude
)

thinned_list <- thin(
    loc.data                = thin_input,
    lat.col                 = "decimalLatitude",
    long.col                = "decimalLongitude",
    spec.col                = "species",
    thin.par                = thin_km,
    reps                    = thin_reps,
    locs.thinned.list.return = TRUE,
    write.files             = FALSE
)

# Select the replicate retaining the most records
best       <- which.max(sapply(thinned_list, nrow))
thinned_df <- thinned_list[[best]]
colnames(thinned_df) <- c("Longitude", "Latitude")

write.csv(thinned_df, thinned_csv, row.names = FALSE)
message("Thinned records: ", nrow(thinned_df), " (from ", nrow(data), ")")

# --------------------------------------------------------------
# STEP 3: Visual comparison (raw vs. thinned)
# --------------------------------------------------------------
world_map <- map_data("world")

p_raw <- ggplot() +
    geom_map(data = world_map, map = world_map,
             aes(map_id = region), fill = "gray90", color = "gray50") +
    geom_point(data = data, aes(x = Longitude, y = Latitude),
               color = "red", size = 1) +
    coord_cartesian(xlim = plot_xlim, ylim = plot_ylim) +
    ggtitle("Original Data") + theme_bw()

p_thinned <- ggplot() +
    geom_map(data = world_map, map = world_map,
             aes(map_id = region), fill = "gray90", color = "gray50") +
    geom_point(data = thinned_df, aes(x = Longitude, y = Latitude),
               color = "blue", size = 1) +
    coord_cartesian(xlim = plot_xlim, ylim = plot_ylim) +
    ggtitle(paste0("Thinned Data (", thin_km, " km)")) + theme_bw()

print(p_raw); print(p_thinned)

# --------------------------------------------------------------
# STEP 4: Alpha hull -> study area polygon
# --------------------------------------------------------------
point_2_alphahull_curve <- function(coords_df,
                                    long = "Longitude",
                                    lat  = "Latitude",
                                    crs  = "EPSG:4326",
                                    point_shp,
                                    hull_shp,
                                    buff = 200000) {
    pts <- vect(coords_df, geom = c(long, lat), crs = crs)
    writeVector(pts, point_shp, overwrite = TRUE)

    range <- getDynamicAlphaHull(coords_df,
                                 coordHeaders = c(long, lat),
                                 clipToCoast  = "no",
                                 fraction     = 1,
                                 partCount    = 1,
                                 buff         = buff)
    plot(range[[1]], col = transparentColor("dark green", 0.5), border = NA)
    points(coords_df[, c(long, lat)], cex = 0.5, pch = 3)

    st_write(range[[1]], hull_shp, delete_dsn = TRUE)
}

point_2_alphahull_curve(thinned_df,
                        point_shp = point_shp,
                        hull_shp  = hull_shp,
                        buff      = hull_buffer_m)

message("Done. Study area written to: ", hull_shp)
