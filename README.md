# Rumex_phylogeography_HHM

Scripts for **Jnawali et al. (2026)** *Historical Landscapes and Climate Dynamics shape Phylogeography of Rumex hastatus in the Himalaya-Hengduan Mountains.* **Plant Diversity**, in press.

This repository contains the complete analytical pipeline for a phylogeographic study of *Rumex hastatus* D. Don (Polygonaceae) across the Himalaya–Hengduan Mountains (HHM). The pipeline integrates population genomics, landscape genetics, species distribution modelling, and ecological niche analyses to understand how Pleistocene climate dynamics and topographic heterogeneity shaped the genetic and distributional patterns of this montane species.

## 1. SNP_calling

(1) `01a_bam_generation.sh` — Script to map paired-end cleaned reads to the reference genome using BWA-MEM, and sort/index alignments with SAMtools, producing one BAM file per sample.

(2) `01b_angsd_snp_calling.sh` — Script to call nuclear SNPs across all samples using ANGSD with the GATK genotype likelihood model (`-GL 2`), outputting MAF estimates, posterior genotypes, and PLINK-format files.

(3) `envs_angsd_setup.md` — Conda environment setup instructions for ANGSD v0.935, BWA, and SAMtools.

(4) `01_snp_calling_README.md` — Folder-level documentation describing the SNP calling pipeline.

## 2. Coding_noncoding

(1) `02_noncoding_snp_extraction.sh` — Script to filter the SNP set from Step 1 into a non-coding subset using the reference GFF annotation and `bedtools intersect`. Non-coding SNPs are used in downstream neutral analyses where the influence of selection should be minimised.

## 3. Genetic_diversity

(1) `03_genetic_diversity_fst.sh` — Script to compute per-population diversity statistics (Watterson's θ, nucleotide diversity π, Tajima's D) and pairwise F<sub>ST</sub> across all population pairs using ANGSD (`realSFS`, `thetaStat`). Includes parallel F<sub>ST</sub> computation and a symmetric F<sub>ST</sub> matrix.

## 4. GEA (Genotype-Environment Association)

(1) `04a_angsd_hard_genotyping.sh` — Script to call hard genotypes (0/1/2) from filtered BAM files using ANGSD with a 0.95 posterior cutoff, producing the input required by RDA and LFMM2.

(2) `04b_gea_preprocessing.R` — Script to load ANGSD hard genotypes, apply post-calling QC filters (missingness ≤ 20%, MAF ≥ 0.05), extract current and future bioclim values at sample locations, impute missing genotypes via sNMF, and perform LD pruning within scaffolds (r² < 0.5).

(3) `04c_rda_lfmm2_intersection.R` — Script to identify adaptive SNPs using two complementary methods: Redundancy Analysis (RDA) with ±3 SD outliers on constrained axes, and Latent Factor Mixed Model (LFMM2) with genomic-control-corrected p-values. The intersection of both methods defines the consensus adaptive SNP set.

## 5. Ensemble_SDM

(1) `05a_occurrence_thinning_alphahull.R` — Script to clean and spatially thin occurrence records (`spThin`, 10 km), then delineate an alpha-hull study area (`rangeBuilder::getDynamicAlphaHull`) used as background for SDM training.

(2) `05b_env_preparation_vif.R` — Script to crop/mask environmental rasters to the study area, run Variance Inflation Factor (VIF, threshold = 10) selection on Current climate, and apply the selected variable subset consistently to Future and LGM scenarios.

(3) `05c_biomod2_ensemble.R` — Script to build ensemble species distribution models using `biomod2` across Current, LGM (CCSM4, MIROC, MPI), and Future (CCSM4, MIROC, MPI) scenarios with 10 algorithms (GLM, GBM, RF, GAM, CTA, ANN, SRE, FDA, MARS, MAXENT), block cross-validation, and ROC/TSS-filtered ensembling.

(4) `05d_suitability_visualization.R` — Script to generate publication-quality habitat suitability maps with three suitability categories (Low / Medium / High) overlaid on a regional basemap with occurrence points.

(5) `05e_range_size_change.R` — Script to quantify pixel-wise habitat change between the Current ensemble projection and each Past/Future scenario using `biomod2::BIOMOD_RangeSize` (gain / loss / stable / never-suitable codes).

(6) `05_samplecsv.csv` — Sample occurrence file used for SDM input (cleaned and thinned).

## 6. Corridor_analysis

(1) `06_dispersal_corridor.R` — Script to test for a Late Pleistocene dispersal corridor between West Himalaya (WH) and Hengduan Mountains (HM) using waypoint-guided least-cost path (`gdistance`) analysis on the LGM ensemble suitability surface. Outputs include the LCP transect, suitability along the path, and habitat-zone classification (corridor / marginal / barrier).

## 7. Niche_divergence

(1) `07_niche_divergence.R` — Script to test whether WH and HM lineages occupy statistically distinct ecological niches using the PCA-env framework (`ecospat`). Computes Schoener's D, Hellinger's I, and runs niche equivalency and similarity tests (both directions) with permutation.

## Data sources

External datasets required to reproduce the pipeline (not redistributed in this repository):

- **Bioclimatic variables (Current).** WorldClim v2.1, 19 bioclim variables at 2.5 arc-minute resolution. https://worldclim.org/data/worldclim21.html
- **Bioclimatic variables (LGM).** WorldClim v1.4 paleoclimate downscaling for the Last Glacial Maximum (~22 kya), GCMs: **CCSM4, MIROC-ESM, MPI-ESM-P**, 2.5 arc-minute resolution. https://www.worldclim.org/data/v1.4/paleo1.4.html
- **Bioclimatic variables (Future, 2070).** WorldClim CMIP5 downscaling under RCP scenarios, GCMs: **CCSM4, MIROC5, MPI-ESM-LR**, 2.5 arc-minute resolution. https://www.worldclim.org/data/cmip5_2.5m.html
- **Occurrence data.** Field collections and herbarium records (deposited at the Kunming Institute of Botany), supplemented by GBIF queries for *Rumex hastatus*. https://www.gbif.org

## Software and dependencies

Core tools used:

| Tool | Version | Purpose |
|------|---------|---------|
| ANGSD | 0.935 | SNP calling, diversity, F<sub>ST</sub>, hard genotypes |
| BWA-MEM | 0.7.17 | Read mapping |
| SAMtools | 1.17 | BAM sorting and indexing |
| PLINK | 1.9 | File-format conversion |
| `biomod2` (R) | 4.x | Ensemble SDM |
| `ecospat` (R) | 4.x | Niche divergence (PCA-env) |
| `gdistance` (R) | 1.x | Least-cost path |
| LEA, vegan, lfmm (R) | latest | GEA (LFMM2 + RDA) |

## Citation

If you use this code or data, please cite:

> Jnawali B., Sun H., Luo D. (2026). Historical Landscapes and Climate Dynamics shape Phylogeography of *Rumex hastatus* in the Himalaya–Hengduan Mountains. *Plant Diversity*, in press.

A `CITATION.cff` file is provided for automatic citation parsing by GitHub and Zenodo.

## License

The code in this repository is released under the [MIT License](LICENSE).

External datasets retain their original licenses — see provider websites.
