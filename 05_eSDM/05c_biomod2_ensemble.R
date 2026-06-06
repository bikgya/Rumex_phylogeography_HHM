# ==============================================================
# Script 08c: Ensemble Species Distribution Modelling (biomod2)
# Project: rumex-landscape-phylogenomics
#
# Builds ensemble SDMs across multiple climate scenarios using
# biomod2. For each scenario (Current, Future GCMs, LGM GCMs):
#   1. Loads VIF-selected environmental rasters from Script 08b
#   2. Formats biomod2 data with disk-based pseudo-absences
#   3. Trains 10 algorithms (MAXENT, GAM, GLM, GBM, CTA, ANN,
#      SRE, FDA, MARS, RF) with block cross-validation
#   4. Builds an ROC/TSS-filtered ensemble model
#   5. Projects the ensemble onto the scenario's climate
#   6. Saves evaluations, variable importance, and projections
#
# Input
#   - occ_csv      : thinned occurrence CSV from Script 08a
#   - scenarios    : list of (name, path) pairs - each path is a
#                    folder of VIF-selected rasters from Script 08b
#
# Output (per scenario, under each scenario's biomod_output/)
#   - <species>_<scenario>_eval.csv         per-algorithm metrics
#   - <species>_<scenario>_imp.csv          variable importance
#   - <species>_<scenario>_eval_EM.csv      ensemble evaluation
#   - <species>_<scenario>_imp_EM.csv       ensemble var importance
#   - ROC_TSS_<scenario>.png                evaluation plot
#   - Projection_<scenario>.png             projection plot
#   - <scenario>_EM/ ... .tif               binary projection raster(s)
#
# Dependencies: biomod2, raster, dismo, tidyverse, tidyterra, ENMeval
# ==============================================================

# ----- USER CONFIGURATION -----
base_dir       <- "results/sdm"
occ_csv        <- "results/sdm/occurrences/rumex_thinned.csv"
species_label  <- "rumex"           # used in resp.name & output names
env_root       <- "results/sdm/env" # output of Script 08b

# Define scenarios (folder names match Script 08b output)
scenarios <- list(
    list(name = "Current",      path = file.path(env_root, "Current_autocorrelated")),
    list(name = "LGM_MPI",      path = file.path(env_root, "LGM_LGM_MPI_autocorrelated")),
    list(name = "LGM_CCSM4",    path = file.path(env_root, "LGM_LGM_CCSM4_autocorrelated")),
    list(name = "LGM_MIROC",    path = file.path(env_root, "LGM_LGM_MIROC_autocorrelated")),
    list(name = "Future_MPI",   path = file.path(env_root, "Future_Future_MPI_autocorrelated")),
    list(name = "Future_CCSM4", path = file.path(env_root, "Future_Future_CCSM4_autocorrelated")),
    list(name = "Future_MIROC", path = file.path(env_root, "Future_Future_MIROC_autocorrelated"))
)

# biomod2 settings
pa_nb_rep        <- 3
pa_nb_absences   <- 1000
pa_strategy      <- "disk"
pa_dist_min      <- 10000
cv_strategy      <- "block"
cv_nb_rep        <- 3
cv_perc          <- 0.7
algorithms       <- c("MAXENT","GAM","GLM","GBM","CTA","ANN","SRE","FDA","MARS","RF")
metric_eval      <- c("ROC","TSS")
ensemble_thresh  <- c(0.8, 0.6)   # ROC, TSS thresholds for ensemble inclusion
seed_val         <- 123
nb_cpu           <- 8
# ------------------------------

library(biomod2)
library(raster)
library(dismo)
library(tidyverse)
library(tidyterra)
library(ENMeval)
library(pryr)

message("Memory used at start: ", format(mem_used()))

# --------------------------------------------------------------
# Load and clean occurrence data
# --------------------------------------------------------------
occ_data <- read.csv(occ_csv)
occ_data$Longitude <- as.numeric(occ_data$Longitude)
occ_data$Latitude  <- as.numeric(occ_data$Latitude)
occ_data <- occ_data[!duplicated(occ_data[, c("Longitude","Latitude")]), ]
rownames(occ_data) <- NULL

# --------------------------------------------------------------
# Loop through scenarios
# --------------------------------------------------------------
for (i in seq_along(scenarios)) {
    scen <- scenarios[[i]]
    cat("\n========================================\n")
    cat("Processing scenario", i, "/", length(scenarios), ":", scen$name, "\n")
    cat("========================================\n")

    # Output directory inside the scenario's folder
    output_dir <- file.path(scen$path, "biomod_output")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    setwd(output_dir)

    # Load environmental rasters
    rst_files <- list.files(scen$path, pattern = "\\.tif$", full.names = TRUE)
    if (length(rst_files) == 0) {
        warning("No .tif files found in ", scen$path, " - skipping")
        setwd(base_dir); next
    }
    rst_env <- raster::stack(rst_files)
    cat("Loaded", length(rst_files), "raster layers\n")

    # Format biomod2 data with pseudo-absences
    myBiomodData <- BIOMOD_FormatingData(
        resp.name      = paste0(species_label, "_", scen$name),
        resp.var       = rep(1, nrow(occ_data)),
        expl.var       = rst_env,
        resp.xy        = occ_data[, c("Longitude","Latitude")],
        PA.nb.rep      = pa_nb_rep,
        PA.nb.absences = pa_nb_absences,
        PA.strategy    = pa_strategy,
        PA.dist.min    = pa_dist_min
    )

    # Modeling options
    myBiomodOptions <- bm_ModelingOptions(
        data.type = "binary",
        strategy  = "bigboss",
        bm.format = myBiomodData
    )

    # Train individual models
    myBiomodModelOut <- BIOMOD_Modeling(
        bm.format      = myBiomodData,
        bm.options     = myBiomodOptions,
        modeling.id    = paste0("AllModels_", scen$name),
        models         = algorithms,
        CV.strategy    = cv_strategy,
        CV.nb.rep      = cv_nb_rep,
        CV.perc        = cv_perc,
        var.import     = 3,
        metric.eval    = metric_eval,
        do.full.models = TRUE,
        nb.cpu         = nb_cpu,
        seed.val       = seed_val
    )

    # Save evaluations and variable importance
    write_csv(as.data.frame(get_evaluations(myBiomodModelOut)),
              paste0(species_label, "_", scen$name, "_eval.csv"))
    write_csv(as.data.frame(get_variables_importance(myBiomodModelOut)),
              paste0(species_label, "_", scen$name, "_imp.csv"))

    # ROC/TSS plot
    png(paste0("ROC_TSS_", scen$name, ".png"), width = 800, height = 600)
    bm_PlotEvalMean(bm.out = myBiomodModelOut, metric_eval)
    dev.off()

    # Build ensemble model (ROC >= 0.8 AND TSS >= 0.6)
    myBiomodEM <- BIOMOD_EnsembleModeling(
        bm.mod               = myBiomodModelOut,
        em.by                = "all",
        metric.select        = metric_eval,
        metric.select.thresh = ensemble_thresh,
        var.import           = 3,
        metric.eval          = metric_eval,
        nb.cpu               = nb_cpu
    )

    write_csv(as.data.frame(get_evaluations(myBiomodEM)),
              paste0(species_label, "_", scen$name, "_eval_EM.csv"))
    write_csv(as.data.frame(get_variables_importance(myBiomodEM)),
              paste0(species_label, "_", scen$name, "_imp_EM.csv"))

    # Project the ensemble
    myBiomodEMProj <- BIOMOD_EnsembleForecasting(
        bm.em         = myBiomodEM,
        proj.name     = paste0(scen$name, "_EM"),
        new.env       = rst_env,
        metric.binary = "all",
        nb.cpu        = nb_cpu,
        do.stack      = TRUE
    )

    png(paste0("Projection_", scen$name, ".png"), width = 1000, height = 800)
    plot(myBiomodEMProj)
    dev.off()

    # Report where projections landed
    proj_folder <- file.path(output_dir, paste0(scen$name, "_EM"))
    if (dir.exists(proj_folder)) {
        proj_files <- list.files(proj_folder, pattern = "\\.tif$", full.names = TRUE)
        cat("Projections saved (", length(proj_files), " files) in: ", proj_folder, "\n", sep = "")
    }

    # Free memory before the next scenario
    rm(myBiomodData, myBiomodOptions, myBiomodModelOut,
       myBiomodEM, myBiomodEMProj, rst_env, rst_files)
    gc()

    setwd(base_dir)
    cat("Finished:", scen$name, "\n")
}

cat("\n========================================\n")
cat("All scenarios completed successfully!\n")
cat("========================================\n")
