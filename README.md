# Rumex\_phylogeography\_HHM

Scripts for **Jnawali et al. (2026)** *Historical Landscapes and Climate Dynamics shape Phylogeography of Rumex hastatus in the Himalaya-Hengduan Mountains.* **Plant Diversity**, in press.

!\[License](https://img.shields.io/badge/License-MIT-yellow.svg)
!\[Status](https://img.shields.io/badge/Status-In\_Press-blue.svg)
!\[Plant Diversity](https://img.shields.io/badge/Journal-Plant\_Diversity-green.svg)

!\[Workflow](docs/workflow.svg)

!\[Study area](docs/HHM\_studyarea.svg)

## Overview

This repository contains the complete analytical pipeline for a phylogeographic study of *Rumex hastatus* D. Don (Polygonaceae) across the Himalaya–Hengduan Mountains (HHM). The pipeline integrates **shallow whole-genome sequencing (sWGS)** data with **species distribution modelling** and **ecological niche analyses** to understand how Pleistocene climate dynamics and topographic heterogeneity shaped the genetic structure and distribution of this widespread montane species.

### Genomic component (sWGS)

Population-level shallow whole-genome sequencing data (\~2–5× coverage) was generated for 162 *Rumex hastatus* individuals from 27 populations spanning the HHM. Because read depth at any single site is low under sWGS, hard-called genotypes are unreliable for many sites; we therefore work with **genotype likelihoods** throughout most of the genomic pipeline. Specifically, nuclear SNPs were identified using **ANGSD v0.935**, which estimates genotype likelihoods under the **GATK model (`-GL 2`)** designed to accommodate low-coverage data. Genotype likelihoods feed directly into the diversity, F<sub>ST</sub>, and SFS analyses (Steps 1–3) without requiring confident per-site calls. For the genotype-environment association analyses (Step 4), high-confidence hard genotypes were called separately using a strict posterior cutoff (0.95) on a filtered subset of individuals.

### Spatial component

The species distribution modelling, dispersal corridor analysis, and ecological niche divergence test (Steps 5–7) use occurrence records and WorldClim climate layers rather than the genomic data — these analyses complement the population genetics by reconstructing past, present, and future habitat availability across the HHM.

\---

## 1\. SNP\_calling — `01\_sWGS\_SNPcalling/`

(1) `01a\_bam\_generation.sh` — Maps paired-end cleaned reads to the reference genome using BWA-MEM, sorts and indexes with SAMtools, producing one BAM file per sample.

(2) `01b\_angsd\_snp\_calling.sh` — Calls nuclear SNPs across all samples using **ANGSD with the GATK genotype likelihood model (`-GL 2`)** to accommodate low coverage. Outputs MAF estimates, posterior genotypes, and PLINK-format files.

(3) `envs\_angsd\_setup.md` — Conda environment setup instructions for ANGSD v0.935, BWA, and SAMtools.

(4) `01\_snp\_calling\_README.md` — Folder-level documentation describing the SNP calling pipeline.

## 2\. Coding\_noncoding — `02\_sWGS\_noncoding/`

(1) `02\_noncoding\_snp\_extraction.sh` — Filters the SNP set from Step 1 into a non-coding subset using the reference GFF annotation and `bedtools intersect`. Non-coding SNPs are preferred for neutral analyses where the influence of selection should be minimised.

## 3\. Genetic\_diversity — `03\_sWGS\_geneticdiversity/`

(1) `03\_genetic\_diversity\_fst.sh` — Computes per-population diversity statistics (Watterson's θ, nucleotide diversity π, Tajima's D) and pairwise F<sub>ST</sub> using ANGSD (`realSFS`, `thetaStat`). Includes parallel F<sub>ST</sub> computation and a symmetric F<sub>ST</sub> matrix output.

## 4\. GEA — `04\_sWGS\_GEA/`

(1) `04a\_angsd\_hard\_genotyping.sh` — Calls hard genotypes (0/1/2) from filtered BAM files using ANGSD with a 0.95 posterior cutoff — required for RDA and LFMM2.

(2) `04b\_gea\_preprocessing.R` — Loads ANGSD hard genotypes, applies QC filters (missingness ≤ 20%, MAF ≥ 0.05), extracts bioclim values at sample locations, imputes missing genotypes via sNMF, and LD-prunes within scaffolds (r² < 0.5).

(3) `04c\_rda\_lfmm2\_intersection.R` — Identifies adaptive SNPs using two complementary methods: **Redundancy Analysis (RDA)** with ±3 SD outliers on constrained axes, and **Latent Factor Mixed Model (LFMM2)** with genomic-control-corrected p-values. The intersection of both methods defines the consensus adaptive SNP set.

## 5\. Ensemble\_SDM — `05\_eSDM/`

(1) `05a\_occurrence\_thinning\_alphahull.R` — Cleans and spatially thins occurrence records (`spThin`, 10 km), delineates an alpha-hull study area used as background for SDM training.

(2) `05b\_env\_preparation\_vif.R` — Crops/masks environmental rasters to the study area, runs Variance Inflation Factor (VIF, threshold = 10) selection on Current climate, and applies the selected variable subset consistently to Future and LGM scenarios.

(3) `05c\_biomod2\_ensemble.R` — Builds ensemble species distribution models using `biomod2` across Current, LGM (CCSM4, MIROC, MPI), and Future (CCSM4, MIROC, MPI) scenarios with 10 algorithms, block cross-validation, and ROC/TSS-filtered ensembling.

(4) `05d\_suitability\_visualization.R` — Generates publication-quality habitat suitability maps with three suitability categories (Low / Medium / High) overlaid on a regional basemap.

(5) `05e\_range\_size\_change.R` — Quantifies pixel-wise habitat change between the Current ensemble projection and each Past/Future scenario using `biomod2::BIOMOD\_RangeSize`.

(6) `05\_samplecsv` — Sample occurrence file used for SDM input (cleaned and thinned).

## 6\. Corridor\_analysis — `06\_Dispersal\_corridor/`

(1) `06\_dispersal\_corridor.R` — Tests for a Late Pleistocene dispersal corridor between West Himalaya (WH) and Hengduan Mountains (HM) using a waypoint-guided least-cost path (`gdistance`) on the LGM ensemble suitability surface. Outputs the LCP transect, suitability along the path, and habitat-zone classification (corridor / marginal / barrier).

## 7\. Niche\_divergence — `07\_Niche\_divergence/`

(1) `07\_niche\_divergence.R` — Tests whether WH and HM lineages occupy statistically distinct ecological niches using the PCA-env framework (`ecospat`). Computes Schoener's D, Hellinger's I, and runs niche equivalency and similarity tests with permutation.

\---

## Data sources

External datasets required to reproduce the pipeline (not redistributed here):

* **Bioclimatic variables (Current).** WorldClim v2.1, 19 bioclim variables at 2.5 arc-minute resolution. https://worldclim.org/data/worldclim21.html
* **Bioclimatic variables (LGM).** WorldClim v1.4 paleoclimate downscaling for the Last Glacial Maximum (\~22 kya), GCMs: **CCSM4, MIROC-ESM, MPI-ESM-P**, 2.5 arc-minute resolution. https://www.worldclim.org/data/v1.4/paleo1.4.html
* **Bioclimatic variables (Future, 2070).** WorldClim CMIP5 downscaling, GCMs: **CCSM4, MIROC5, MPI-ESM-LR**, 2.5 arc-minute resolution. https://www.worldclim.org/data/cmip5\_2.5m.html
* **Occurrence data.** Field collections and herbarium records (deposited at the Kunming Institute of Botany, KIB), supplemented by GBIF queries for *Rumex hastatus*. https://www.gbif.org

## Software and dependencies

|Tool|Version|Purpose|
|-|-|-|
|ANGSD|0.935|SNP calling, diversity, F<sub>ST</sub>, hard genotypes (sWGS)|
|BWA-MEM|0.7.17|Read mapping|
|SAMtools|1.17|BAM sorting and indexing|
|PLINK|1.9|File-format conversion|
|`biomod2` (R)|4.x|Ensemble SDM|
|`ecospat` (R)|4.x|Niche divergence (PCA-env)|
|`gdistance` (R)|1.x|Least-cost path|
|`LEA`, `vegan`, `lfmm` (R)|latest|GEA (LFMM2 + RDA)|

## Citation

If you use this code or data, please cite:

> Jnawali B., Sun H., Luo D. (2026). Historical Landscapes and Climate Dynamics shape Phylogeography of \*Rumex hastatus\* in the Himalaya–Hengduan Mountains. \*Plant Diversity\*, in press.

A `CITATION.cff` file is provided for automatic citation parsing by GitHub and Zenodo.

## License

The code in this repository is released under the [MIT License](LICENSE). External datasets retain their original licenses — see provider websites.

