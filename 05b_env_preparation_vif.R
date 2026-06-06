# ==============================================================
# Script 05b: Environmental Variable Preparation (VIF)
# Project: rumex-landscape-phylogenomics
#
# Prepares environmental rasters for ensemble SDM (Script 05c):
#   1. Crops/masks each climate scenario to the study area
#      (alpha-hull polygon from Script 08a)
#   2. Runs Variance Inflation Factor (VIF) selection on the
#      Current climate to remove multicollinear variables
#      (threshold = 10)
#   3. Subsets all scenarios (Current + Future GCMs + LGM GCMs)
#      to the VIF-selected variable set so that biomod2 can
#      train on Current and project consistently to Future/LGM.
#
# Input
#   - study_area_shp  : alpha-hull polygon from Script 08a
#   - <scenario>.csv  : one CSV per scenario listing absolute
#                       paths to that scenario's WorldClim/Soil
#                       rasters (one path per line)
#
# Output
#   - <out_root>/<scenario>/<var>.tif    masked, VIF-selected rasters
#
# Dependencies: raster, terra, sf, usdm, tidyverse
# ==============================================================

# ----- USER CONFIGURATION -----
study_area_shp <- "results/sdm/study_area/rumex_alphahull.shp"
out_root       <- "results/sdm/env"
vif_threshold  <- 10

# Path-CSV files (each lists raster paths for one scenario)
current_csv <- "data/raw/env/Current.csv"
future_csvs <- c(
    "Future_CCSM4" = "data/raw/env/Future_CCSM4.csv",
    "Future_MIROC" = "data/raw/env/Future_MIROC.csv",
    "Future_MPI"   = "data/raw/env/Future_MPI.csv"
)
lgm_csvs <- c(
    "LGM_CCSM4" = "data/raw/env/LGM_CCSM4.csv",
    "LGM_MIROC" = "data/raw/env/LGM_MIROC.csv",
    "LGM_MPI"   = "data/raw/env/LGM_MPI.csv"
)
# ------------------------------

library(raster)
library(tidyverse)
library(usdm)
library(stringr)

# Load study area
study_area <- shapefile(study_area_shp)

# Helper: load a scenario's rasters, crop & mask to study area
load_scenario <- function(csv_path, study_area) {
    rst_paths <- read_csv(csv_path, show_col_types = FALSE)[[1]]
    rst_stack <- stack(rst_paths)
    rst_masked <- mask(crop(rst_stack, study_area), study_area)
    names(rst_masked) <- str_replace(basename(rst_paths), ".tif$", "")
    rst_masked
}

# Helper: write a stack to a folder (one .tif per layer)
write_stack <- function(stk, out_dir) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    for (var in names(stk)) {
        writeRaster(stk[[var]],
                    filename  = file.path(out_dir, paste0(var, ".tif")),
                    overwrite = TRUE)
    }
}

# --------------------------------------------------------------
# STEP 1: Current climate — crop/mask + VIF selection
# --------------------------------------------------------------
message("=== Current: cropping and running VIF ===")
current_stack <- load_scenario(current_csv, study_area)

vif_df  <- na.omit(as.data.frame(current_stack))
vif_res <- vifstep(vif_df, th = vif_threshold)
selected_vars <- vif_res@results$Variables
message("VIF retained ", length(selected_vars), " variables: ",
        paste(selected_vars, collapse = ", "))

# Write VIF-selected current rasters
write_stack(current_stack[[selected_vars]],
            file.path(out_root, "Current_autocorrelated"))

# --------------------------------------------------------------
# STEP 2: Future and LGM scenarios — crop/mask + subset to VIF set
# --------------------------------------------------------------
process_other_scenarios <- function(scenario_csvs, label_prefix) {
    for (i in seq_along(scenario_csvs)) {
        scen_name <- names(scenario_csvs)[i]
        csv_path  <- scenario_csvs[i]
        message("=== ", scen_name, ": cropping and applying VIF subset ===")

        stk <- load_scenario(csv_path, study_area)

        # Match by variable name; warn if any VIF-selected variables are missing
        missing <- setdiff(selected_vars, names(stk))
        if (length(missing) > 0) {
            warning(scen_name, " is missing VIF-selected variables: ",
                    paste(missing, collapse = ", "))
        }
        keep <- intersect(selected_vars, names(stk))

        out_dir <- file.path(out_root, paste0(label_prefix, "_", scen_name, "_autocorrelated"))
        write_stack(stk[[keep]], out_dir)
    }
}

process_other_scenarios(future_csvs, "Future")
process_other_scenarios(lgm_csvs,    "LGM")

message("Done. VIF-selected variable set was applied to all scenarios.")
