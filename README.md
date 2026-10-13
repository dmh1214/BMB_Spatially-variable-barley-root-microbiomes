# Manuscript analysis code for Figures 2c, 3b and 5

This reviewer package contains the study-specific R workflows for spatial modulation of host genetic effects on microbial composition, local root-to-bulk host selection, within-genotype dispersion, SNP-associated biomarker occupancy, spatial likelihood-ratio tests, SNP-by-PCNM interactions, and biomarker stability. It implements established statistical methods using existing R packages; it does not introduce a new software tool or statistical algorithm.

The scope is deliberately limited to Figures 2c, 3b and 5 and the associated sensitivity analyses. This package is not an end-to-end reproduction of all manuscript results. Standard PCoA analyses, biomarker discovery, upstream sequencing processing and GEMMA GWAS are outside this package. GWAS-derived node and edge tables and previously selected PCNM axes are supplied as inputs.

## System requirements

R 4.5.1 on Windows 11 x64 was used for the validation reported in `VALIDATION.md`. Other operating systems have not been tested. See `tested_package_versions.csv` for the package versions used during testing. No specialist hardware or GPU is required. Figure 3b converts a large bacterial BIOM table to a dense matrix and temporarily allocates approximately 2.6 GiB for that conversion alone; allow several additional GB for intermediate objects. A desktop with 16 GB RAM is a reasonable recommendation, not a measured minimum.

Dependencies: tidyverse, qiime2R, vegan, ggeffects, lme4, lmerTest, patchwork, ggthemes, performance, broom, cowplot and ggsci. Dependencies and their version-specific APIs may affect numerical or graphical reproducibility.

## Installation

Unzip the package into a writable directory and install the listed dependencies in R:

```r
install.packages(c("tidyverse", "vegan", "ggeffects", "lme4", "lmerTest",
                   "patchwork", "ggthemes", "performance", "broom", "cowplot", "ggsci"))
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
```

For matching the tested environment, use the versions recorded in `tested_package_versions.csv`, rather than assuming the newest versions are equivalent. qiime2R is distributed through GitHub. Installation time depends on downloads and whether source compilation is necessary; budget approximately 10–30 minutes on a typical desktop. This is an estimate; installation was not timed because dependencies were already present.

## Demo

From the package directory:

```text
Rscript --vanilla run_selected.R demo
```

The demo uses the supplied real genotype-level summary table to regenerate the positive biomarker pool CV ratio plot. Expected outputs are `outputs/fig5d/positive_biomarker_pool_cv_ratio_spring.pdf` and `.png`. Values below one indicate a lower CV for pooled abundance than the median individual-biomarker CV. Allow approximately 5–30 seconds on a typical desktop; see `VALIDATION.md` for local execution results. This demo begins with an existing analysis summary, rather than recomputing host-selection strength from raw microbiome counts.

## Full selected workflows

```text
Rscript --vanilla run_selected.R fig2c
Rscript --vanilla run_selected.R fig3b
Rscript --vanilla run_selected.R fig5
Rscript --vanilla run_selected.R all
```

The Figure 5 driver runs occupancy counting before plotting and runs the SNP-by-PCNM model workflow through the crossover plotting script. Each module runs in its own R process so it does not depend on objects in an interactive workspace. Generated results and logs go to `outputs/`. Source data are read from `data/processed/`.

| Figure | Main scripts | Analysis |
|---|---|---|
| 2c | `scripts/fig2c_spatial_interaction.R` | Host genetic distance × raw XY/PCNM spatial metrics, top-ten selected PCNM axes, interaction confidence intervals and BH classification |
| 3b; associated outlier sensitivity and distance illustration | `scripts/fig3b_host_selection_dispersion.R` | Six nearest bulk samples, Bray-Curtis distances, genotype dispersion, mixed models and outlier sensitivity |
| 5a | `scripts/fig5a_occupancy_counts.R`, `scripts/fig5a_occupancy_plot.R`, `scripts/fig5a_network_weights.R` | Occupancy from rarefied counts and weighted SNP–biomarker edges |
| 5b | `scripts/fig5b_spatial_lrt.R` | Genotype random intercept, full/top-ten PCNM fixed effects, ML model comparison and marginal R² |
| 5c | `scripts/fig5c_snp_pcnm_models.R`, `scripts/fig5c_crossover_plot.R` | Multilocus SNP-by-PCNM models and allele-specific slopes |
| 5d | `scripts/fig5d_individual_stability.R`, `scripts/fig5d_pool_stability.R` | Individual biomarker CV, supplied pooled CV summaries, regression and Spearman tests |

The panel map was checked against the Nature Communications submission manuscript on 13 October 2026. A copy of its embedded Figure 5 is provided in `reference_results/` for comparison. Some scripts generate additional supporting plots. Final multipanel Illustrator/Cytoscape composition is manual; these scripts reproduce analysis tables and panel plots, not the final assembled artwork.

## Inputs and use with another dataset

Retain the supplied file names and schemas when replacing inputs. Count matrices, rCLR matrices, metadata and genotype matrices have different orientations; use the bundled files as schema examples. Sample IDs, genotype IDs, feature IDs, trial names, SNP names and PCNM names must agree across inputs. Do not substitute differently normalized tables and expect the published values to remain unchanged.

The rCLR matrices are exact column subsets of the original processed matrices, containing only biomarkers needed here. No transformation was recomputed after subsetting. PCNM selection and GWAS discovery are upstream inputs. Stability plotting reads supplied pooled host-selection/CV tables; regeneration of those upstream tables is not included. Individual biomarker CVs are recalculated from the supplied count matrices.

`reference_results/` contains previously exported quantitative results for comparison. `originals/` contains the source scripts before packaging adjustments. `provenance.csv` records the source file and packaged-file checksum. `CHANGELOG.md` explains the packaging adjustments. This version retains the original manuscript analysis and annotation behavior. In Figure 5d, rho is calculated with Spearman correlation, while `p_value` and the aggregate displayed P are from the linear-regression slope test. In Figure 5c, annotations are from the full multivariable model, while visual trend lines are separate within-allele linear regressions. These are the original workflows, not revised analyses. See `VALIDATION.md` for the exact reproduction limits.

## License and availability

This copy is supplied for editorial and peer-review evaluation. An open-source license for public redistribution has not been selected by the authors. Do not assume that the dependencies' licenses license the authors' analysis scripts. Before a public release, the authors should choose a license and add its full text. No public archive DOI or Code Ocean capsule is specified in this package.

## Figure 2c inputs

The added entry point uses the original `interaction_only_top10_q33.R` workflow. The three RData files under `data/processed/fig2c/` contain the original aligned metadata, sample IDs, genetic/microbial distance objects and PCNM scores required by that script. Large unrelated workspaces are omitted. The selected-model objects are represented by their original terms, preserving the top-ten axis order exactly; no distances or PCNM selection were recomputed during extraction. Pairwise PCNM context uses the original pair-mean PCNM-score calculation.

Run `Rscript --vanilla run_selected.R fig2c`. Results are written to `outputs/fig2c/fig2c_interaction_statistics.csv` and the `main_interaction_only_top10_q33_q67` PDF/PNG. The original script outputs four panels in a two-by-two arrangement; the manuscript places them in one row during final manual assembly. The suffix in the original filename is preserved.
