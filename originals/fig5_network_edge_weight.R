rm(list=ls())
library(tidyverse)


edge <- read.csv("Network/mgwas_sum_edge.csv")

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

write.csv(edge_df,"Network/edge_weighted.csv")
