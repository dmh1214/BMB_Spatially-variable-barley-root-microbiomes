# Code accompanying the barley root microbiome manuscript

Manuscript: *Spatially variable barley root microbiomes do not imply weak host genetic control at field scale*.

This repository provides figure analysis and plotting code using established statistical methods. It covers main Figures 1–5, Supplementary Figures S2–S13, Extended Data Figure 1, and supporting table workflows. No dedicated script was located for the field map, experimental layout and weather artwork in Supplementary Figure S1. Final panel assembly, maps and network layout may involve manual editing.

## Scope of this submission: source code only

**Only source code and accompanying documentation are provided in this release.** Experimental datasets, sample metadata, fitted-workspace objects, archived numerical results, expected-output files and manuscript images are not included. Package/version inventories and figure/path indexes are documentation, not experimental data.

The code describes the study workflows, but it cannot reproduce the manuscript results or provide a runnable scientific demo from this repository alone. Execution requires the original input files and, for some plotting scripts, saved upstream results. `required_input_paths.csv` documents literal input paths; dynamically constructed paths are also present in the sources. This README does not specify a separate data-access arrangement; data availability should be described in the manuscript's Data Availability statement.

All analysis commands below are conditional on having the required external inputs. Listing modules is the only example here that can run without experimental data.

## Files

- `scripts/`: 36 figure modules, extracted from original R/R Markdown sources.
- `original_sources/`: unmodified source snapshots, including R Markdown notebooks, GEMMA commands and annotation utilities.
- `figure_code_index.csv` and `FIGURE_MAP.md`: figure numbers, module keys and original working directories.
- `source_provenance.csv`: original paths and SHA-256 checksums of the archived sources.
- `dependencies.csv`: packages explicitly loaded by the modules; transitive dependencies must also be installed.
- `run_figure.R`: optional launcher resolving external input locations and isolating common output operations.

## System requirements and versions

Local representative checks used R 4.5.1 on Windows 11 x64. Other operating systems have not been tested. `local_package_versions.csv` records available local dependency versions; inclusion in that file does not mean every module was executed. `selected_workflow_tested_versions.csv` records the earlier selected-workflow environment. The original notebooks do not specify a complete historical environment lockfile.

The R plotting workflows require no specialist hardware. At least 16 GB RAM is recommended for large microbiome matrices; this is not a measured minimum. Full feature selection and GWAS require more time and memory than panel plotting. Original GEMMA shell workflows require a Unix-compatible shell, GEMMA and PLINK, and appropriate genotype inputs. Those external executables and their versions were not validated in this release. Annotation utilities require Python and the genome/annotation inputs named in their sources; workbook assembly uses Node.js. These auxiliary workflows were archived, not newly executed.

## Installation

Download or clone the repository. In R, inspect the required packages:

```r
d <- read.csv("dependencies.csv")
packages <- sort(unique(c(d$package, "jsonlite", "tidyverse", "vegan")))
packages
install.packages(setdiff(packages, c("qiime2R", "linkET")))
install.packages("remotes")
remotes::install_github("jbisanz/qiime2R")
remotes::install_github("Hy4m/linkET")
```

Some packages may require compatible archived versions or additional build tools. Match the recorded versions when evaluating manuscript correspondence. Installation has not been timed; allow approximately 20–60 minutes on a typical desktop with binary packages and an internet connection, and longer if compilation is needed. This is an estimate, not a benchmark. Install GEMMA/PLINK separately only if running the upstream GWAS workflows.

## Demo availability

No real or simulated demo dataset is included. No executable data-based demo or expected-output reference is supplied. Demo runtime is therefore not applicable to this code-only release. The private-input timing below documents a local check and is not a bundled demo benchmark.

## Running with separately supplied inputs

The launcher accepts one module key at a time. Listing modules requires only base R and no data:

```text
Rscript --vanilla run_figure.R list
```

For actual analysis, provide the original project folder, retaining its input directory structure. Example in Windows PowerShell from this code folder:

```powershell
$env:BMB_DATA_ROOT = "C:/path/to/private/BMB_AUT_MS"
$env:BMB_OUTPUT_ROOT = "C:/path/to/private/figure_outputs"
Rscript --vanilla run_figure.R figS2_bacteria_taxonomic_composition --check
Rscript --vanilla run_figure.R figS2_bacteria_taxonomic_composition
```

The path check covers literal filenames only. When the required external inputs are present, that module generates `family_composition.pdf`, `notebook_plots.pdf` and `session_info.txt` in the module output folder. Local execution with dependencies already installed took about 10 seconds, including R startup; this is a private-input example, not a bundled demo. These generated files are not distributed with this release.

Other examples:

```text
Rscript --vanilla run_figure.R fig2c_spatial_interaction
Rscript --vanilla run_figure.R figS10_biomass
Rscript --vanilla run_figure.R fig5d_plot_archived_statistics
```

The optional `BMB_CURATED_ROOT` variable can point to the earlier selected-workflow input package. The launcher recognizes selected legacy paths through `input_aliases.json`. It changes common input/output paths and CSV BOM handling only, preserving analysis statements. Common output filenames are flattened into separate module directories. It is not a sandbox: namespace-qualified or external file operations may bypass its wrappers. Run on a private working copy of the original project. Upstream modules and saved results are not automatically connected into a single pipeline.

To run the original notebooks directly, use the original project inputs, original working directory recorded in the index, and the required original interactive dependencies. FigS10's executable module contains only the original spring/winter biomass rank-plot blocks (lines 1–57 and 254–286); the complete source, including unrelated field-map calls, is archived unchanged.

## Manuscript correspondence

Original statistical models, filters and analysis logic have been retained. Packaging changes include file-path adaptation, selected notebook-block extraction, explicit input loading in selected workflows, and minor file-reference or graphical-parameter corrections. Unmodified source snapshots are provided in `original_sources/`; the executable modules are not described as comment-only edits. See `VALIDATION.md` for what was checked. Newly archived sources should not be interpreted as a claim that every final manuscript panel has been reproduced or compared pixel by pixel.

Fig4 feature-selection runs are stochastic. To reconstruct the manuscript SHAP/LME/GWAS panels, use the original saved feature-selection and GEMMA result objects with the corresponding plotting modules. Re-running discovery can produce different selected features. The original shell output-directory spellings and original prevalence filters have been retained; input organization may require attention before a complete upstream rerun.

Fig5d has two distinct workflows: `fig5d_individual_stability` / `fig5d_pool_stability` retain original recomputation, while `fig5d_plot_archived_statistics` plots the original prepared data using archived manuscript annotation statistics. The latter requires original result CSVs externally and does not newly recompute the archived Spearman coefficients. Known numerical/annotation discrepancies from the earlier audit have not been silently corrected. Final multipanel composition and some manually placed labels are outside automated reproduction.

## Use with other data

Follow each source's input schema. Count/rCLR matrices, metadata, genotype distances, PCNM axes and GWAS result objects have different orientations and must have consistent sample, genotype, feature, trial and SNP identifiers. Replace inputs only after checking those assumptions; a different normalization or upstream feature set changes the scientific results. This release provides study-specific analysis code, not a general-purpose application.
