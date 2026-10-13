# Prepared for reviewer use; original analysis retained in originals/.
if (!exists("project_root", inherits=FALSE)) {
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), winslash="/", mustWork=TRUE)
}
module_output <- file.path(project_root, "outputs", "fig5a")
dir.create(module_output, recursive=TRUE, showWarnings=FALSE)
setwd(module_output)
dir.create("Network/k_tables", recursive=TRUE, showWarnings=FALSE)
on_error <- NULL

library(tidyverse)


edge <- read.csv(file.path(project_root, "data/processed/network/mgwas_sum_edge.csv"))

edge_df <- edge |>
  mutate(weight_p_raw = -log10(mGWAS_pw),
         weight_r_raw=lme_slope) |>
  group_by(shap_sign) |>
  mutate(
    weight_p_scaled = (weight_p_raw - min(weight_p_raw, na.rm = TRUE)) /
      (max(weight_p_raw, na.rm = TRUE) - min(weight_p_raw, na.rm = TRUE)),
    weight_r_scaled=(weight_r_raw - min(weight_r_raw, na.rm = TRUE)) /
      (max(weight_r_raw, na.rm = TRUE) - min(weight_r_raw, na.rm = TRUE))
  ) |>
  ungroup()

write.csv(edge_df,file.path(module_output, "Network", "edge_weighted.csv"))
