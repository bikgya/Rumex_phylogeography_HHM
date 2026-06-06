# ==============================================================
# Script 08e: Range Size Change (Current vs. Past / Future)
# Project: rumex-landscape-phylogenomics
#
# Quantifies pixel-wise habitat change between the Current
# ensemble projection (Script 08c) and each Past (LGM) or
# Future scenario using biomod2's BIOMOD_RangeSize() function.
# For every comparison, the script writes:
#   - A GeoTIFF showing per-pixel change codes:
#         -2 = loss  (suitable now, not in scenario)
#         -1 = stable (suitable in both)
#          0 = never suitable
#          1 = gain  (not suitable now, suitable in scenario)
#   - The Compt.By.Models summary printed to console
#
# Input
#   - Binary ensemble projection rasters (.tif) from Script 08c,
#     one per scenario (Current + each Past/Future GCM)
#
# Output
#   - <out_dir>/RangeChange_<scenario>.tif    pixel-change raster
#
# Dependencies: biomod2, raster
# ==============================================================

# ----- USER CONFIGURATION -----
in_dir   <- "results/sdm/range_change/inputs"     # folder with binary .tif files
out_dir  <- "results/sdm/range_change/outputs"
current_raster <- "Rumexcurrent.tif"              # current scenario file

# Scenarios to compare against Current.
# Each entry: name = label used in output filename, file = raster filename in in_dir
scenarios <- list(
    list(name = "Future_CCSM4", file = "RumexfutCC4.tif"),
    list(name = "Future_MPI",   file = "RumexfutMPI.tif"),
    list(name = "LGM_CCSM4",    file = "LGM_CC4.tif"),
    list(name = "LGM_MIROC",    file = "LGMMIROC.tif"),
    list(name = "LGM_MPI",      file = "LGMMPI.tif")
)
# ------------------------------

library(biomod2)
library(raster)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load current projection once
current <- raster(file.path(in_dir, current_raster))
message("Current CRS: ", proj4string(current))

# Loop over scenarios
for (sc in scenarios) {
    cat("\n=========================================\n")
    cat("Comparing Current vs.", sc$name, "\n")
    cat("=========================================\n")

    future_path <- file.path(in_dir, sc$file)
    if (!file.exists(future_path)) {
        warning("Missing raster - skipping: ", future_path)
        next
    }

    future <- raster(future_path)
    message(sc$name, " CRS: ", proj4string(future))

    # Align CRS if they differ (BIOMOD_RangeSize requires matching grids)
    if (proj4string(current) != proj4string(future)) {
        message("  CRS mismatch -> reprojecting future to match current")
        future <- projectRaster(future, crs = projection(current))
    }

    # Compute range size change
    change <- BIOMOD_RangeSize(proj.current = current, proj.future = future)

    # Print summary
    print(change$Compt.By.Models)

    # Write the pixel-wise change raster
    out_tif <- file.path(out_dir, paste0("RangeChange_", sc$name, ".tif"))
    writeRaster(raster(change$Diff.By.Pixel),
                filename = out_tif,
                format   = "GTiff",
                overwrite = TRUE)
    cat("  Wrote:", out_tif, "\n")
}

cat("\n=========================================\n")
cat("All range-size comparisons complete.\n")
cat("Output rasters in:", out_dir, "\n")
cat("=========================================\n")
