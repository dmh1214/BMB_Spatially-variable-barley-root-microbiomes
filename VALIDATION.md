# Validation scope

Checked on 13 October 2026 with R 4.5.1 on Windows 11 x64.

- 46 archived original/supporting source files match the project originals byte for byte (SHA-256 in `source_provenance.csv`).
- All 36 executable module bodies match the original R source or executable R Markdown chunks. FigS10 contains only the two original rank-plot blocks, without changing statements. Chunks explicitly marked `eval=FALSE` are retained in the archived notebook but excluded from executable modules.
- All 57 `.R` files in the release, including launcher and supporting archived scripts, passed R syntax parsing. This establishes syntax only.
- Six modules completed with external private project inputs: FigS2 bacterial composition, FigS3 fungal composition, FigS10 biomass rank plots, Fig2c interactions, FigS7 context regressions, and Extended Data Fig1 nonsignificant examples. Outputs were stored outside this code-only release. FigS2/S3 and the interaction workflows emitted package/plotting warnings; successful completion is not proof of pixel identity.
- FigS10's complete original source initially stopped at an unrelated `agronomymap` call with no function definition located. The executable module was limited to the manuscript rank plots and then completed for both seasons. The full original source was preserved.

See `module_validation.csv` for per-module scope. New full discovery, upstream GWAS, annotation utilities and all remaining notebook workflows were not executed in this release. No full comparison of newly included output panels with the final manuscript artwork was performed.

## Earlier selected-workflow audit

The previous selected package covered Fig2c, Fig3b and Fig5. That audit found matching Fig2c model statistics and historical standalone panel render; matching Fig3b coefficients/P values to stored precision; matching Fig5a occupancy counts, Fig5b LRT statistics, and Fig5c annotated model statistics. Several standalone plots matched archived renderings; final manuscript assembly was manual.

Fig5c retained the original distinction between multivariable annotated slopes and bivariate drawn regression lines. Fig5d retained original linear-model P values and the separate archived-label plotting mode; recalculated Spearman coefficients can differ because of numerical ranking among nearly tied sparse CV values. Some final artwork labels/legend details differ from automatically rounded source values. These points were documented previously and are not silently repaired here. The release does not claim all observed differences are floating-point errors.

Fig4 sources preserve their original filters and stochastic feature-selection routines. Code archival is not proof that a new discovery run will return the same biomarker counts. Original saved results are required for plotting the manuscript panels. Auxiliary GEMMA shell commands retain their original output-directory spellings, which may need reconciliation when preparing a separate rerun environment.

## Data and reproducibility limit

No measured inputs or numerical result tables are included. Listing modules works without data, but a referee needs separate access to the required inputs to test the scientific workflows. This is not a complete runnable reviewer dataset package.
