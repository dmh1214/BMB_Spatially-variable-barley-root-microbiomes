# Figures: Fig3b; FigS8,S9
# Original source: 0_Submission/4_NatComm/Manuscript_Figure_Code/scripts/fig3b_host_selection_dispersion.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

# Prepared for reviewer use; original analysis retained in originals/.
if (!exists("project_root", inherits=FALSE)) {
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), winslash="/", mustWork=TRUE)
}
module_output <- file.path(project_root, "outputs", "fig3b")
dir.create(module_output, recursive=TRUE, showWarnings=FALSE)
setwd(module_output)
dir.create("Network/k_tables", recursive=TRUE, showWarnings=FALSE)
on_error <- NULL



library(qiime2R)
library(tidyverse)
library(vegan)
library(ggeffects)
library(lme4)
library(lmerTest)
library(patchwork)


bac <- read_q2biom(file.path(project_root, "data/processed/bulk/bulk_bacteria_feature_table.biom")) # sample as col, otu as row
colnames(bac) <- sub("^FL_", "FL", colnames(bac))

fun <- readRDS(file.path(project_root, "data/processed/bulk/fun_otu_count_with_bulk.rds")) |> t() |> as.data.frame()# make format consistant
colnames(fun) <- sub("^FL_", "FL", colnames(fun))
meta <- read.csv(file.path(project_root, "data/processed/bulk/meta_file.csv"))

feature_min_prevalence <- 0.05
feature_min_total_count <- 10


meta_aut <- meta |> 
  filter(Location=="Austria",Trial!="DivPan")

meta_aut <- meta |> 
  filter(Location == "Austria")


meta4bulk <- data.frame(sampleid = colnames(bac)) |>
  filter(str_detect(sampleid, "^bulk_?\\d+$")) |>
  mutate(
    bulk_num = as.numeric(str_extract(sampleid, "\\d+$")),
    Trial = case_when(
      bulk_num >= 1    & bulk_num <= 600  ~ "WinterPan",
      bulk_num >= 601  & bulk_num <= 1200 ~ "DivPan",
      bulk_num >= 1201 & bulk_num <= 1800 ~ "SpringPan",
      TRUE ~ NA_character_
    ),
    Plot = case_when(
      bulk_num >= 1    & bulk_num <= 600  ~ bulk_num,
      bulk_num >= 601  & bulk_num <= 1200 ~ bulk_num - 600,
      bulk_num >= 1201 & bulk_num <= 1800 ~ bulk_num - 1200,
      TRUE ~ NA_real_
    ),
    Entry = "Bulk_Soil"
  ) |>
  left_join(
    meta_aut |>
      select(Trial, Plot, Col, Row) |>
      distinct(),
    by = c("Trial", "Plot")
  ) |>
  mutate(Side = NA_character_) |>
  select(sampleid, Trial, Entry, Plot, Col, Row, Side)

meta4all <- rbind(meta_aut |> select(sampleid,Trial,Entry,Plot,Col,Row,Side), meta4bulk)


filter_low_abundance_features <- function(mat,
                                          min_prevalence = feature_min_prevalence,
                                          min_total_count = feature_min_total_count) {
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"

  prevalence <- rowMeans(mat > 0)
  total_count <- rowSums(mat)

  keep_features <- prevalence >= min_prevalence &
    total_count >= min_total_count

  mat[keep_features, , drop = FALSE]
}

calc_host_selection <- function(mat, meta, panel_name,
                                min_depth = 3000,
                                min_prevalence = feature_min_prevalence,
                                min_total_count = feature_min_total_count,
                                method = "bray",
                                x_col = "Col",
                                y_col = "Row",
                                k_nearest_bulk = 6) {
  
  meta_sub <- meta |>
    dplyr::filter(Trial == panel_name)
  
  keep_samples <- intersect(colnames(mat), meta_sub$sampleid)
  meta_sub <- meta_sub |>
    dplyr::filter(sampleid %in% keep_samples)
  
  mat_sub <- mat[, meta_sub$sampleid, drop = FALSE]
  
  # filter samples by sequencing depth
  mat_sub <- mat_sub[, colSums(mat_sub) > min_depth, drop = FALSE]
  meta_sub <- meta_sub |>
    dplyr::filter(sampleid %in% colnames(mat_sub))
  meta_sub <- meta_sub[match(colnames(mat_sub), meta_sub$sampleid), , drop = FALSE]
  
  # filter low-prevalence and low-abundance taxa
  n_features_before <- nrow(mat_sub)
  mat_sub <- filter_low_abundance_features(
    mat_sub,
    min_prevalence = min_prevalence,
    min_total_count = min_total_count
  )
  n_features_after <- nrow(mat_sub)
  mat_sub <- mat_sub[, colSums(mat_sub) > 0, drop = FALSE]
  
  meta_sub <- meta_sub |>
    dplyr::filter(sampleid %in% colnames(mat_sub))
  meta_sub <- meta_sub[match(colnames(mat_sub), meta_sub$sampleid), , drop = FALSE]
  
  if (ncol(mat_sub) == 0) stop("No samples remaining after filtering in panel: ", panel_name)
  if (nrow(mat_sub) == 0) stop("No taxa remaining after filtering in panel: ", panel_name)
  
  mat_sub <- as.matrix(mat_sub)
  storage.mode(mat_sub) <- "numeric"
  
  mat_rel <- vegan::decostand(mat_sub, method = "total", MARGIN = 2)
  dist_mat <- as.matrix(vegan::vegdist(t(mat_rel), method = method))
  
  bulk_ids  <- meta_sub$sampleid[meta_sub$Entry == "Bulk_Soil"]
  rhizo_ids <- meta_sub$sampleid[meta_sub$Entry != "Bulk_Soil"]
  
  if (length(bulk_ids) < k_nearest_bulk) {
    stop("Not enough bulk samples for k_nearest_bulk = ", k_nearest_bulk)
  }
  if (length(rhizo_ids) == 0) stop("No rhizosphere samples found in panel: ", panel_name)
  
  # global bulk baseline
  bulk_dist_mat <- dist_mat[bulk_ids, bulk_ids, drop = FALSE]
  bulk_baseline_global <- mean(bulk_dist_mat[upper.tri(bulk_dist_mat)], na.rm = TRUE)
  
  bulk_meta <- meta_sub |>
    dplyr::filter(sampleid %in% bulk_ids) |>
    dplyr::select(sampleid, !!rlang::sym(x_col), !!rlang::sym(y_col))
  
  rhizo_meta <- meta_sub |>
    dplyr::filter(sampleid %in% rhizo_ids) |>
    dplyr::select(sampleid, Entry, !!rlang::sym(x_col), !!rlang::sym(y_col))
  
  sample_level <- lapply(rhizo_ids, function(s) {
    one_rhizo <- rhizo_meta[rhizo_meta$sampleid == s, , drop = FALSE]
    
    # spatial distance to all bulk plots
    d_xy <- sqrt(
      (bulk_meta[[x_col]] - one_rhizo[[x_col]])^2 +
      (bulk_meta[[y_col]] - one_rhizo[[y_col]])^2
    )
    
    nearest_idx <- order(d_xy)[1:k_nearest_bulk]
    nearest_bulk_ids <- bulk_meta$sampleid[nearest_idx]
    
    # microbial distance from rhizosphere sample to bulk
    dist_global <- mean(dist_mat[s, bulk_ids], na.rm = TRUE)
    dist_local  <- mean(dist_mat[s, nearest_bulk_ids], na.rm = TRUE)
    
    # local bulk baseline: mean pairwise distance among nearest bulk samples
    bulk_local_mat <- dist_mat[nearest_bulk_ids, nearest_bulk_ids, drop = FALSE]
    bulk_baseline_local <- mean(bulk_local_mat[upper.tri(bulk_local_mat)], na.rm = TRUE)
    
    data.frame(
      sampleid = s,
      dist_to_bulk_global = dist_global,
      dist_to_bulk_local = dist_local,
      bulk_baseline_global = bulk_baseline_global,
      bulk_baseline_local = bulk_baseline_local,
      mean_xy_to_nearest_bulk = mean(d_xy[nearest_idx], na.rm = TRUE),
      nearest_bulk_n = k_nearest_bulk
    )
  }) |>
    bind_rows() |>
    mutate(
      dist_to_bulk_global_norm = dist_to_bulk_global / bulk_baseline_global,
      dist_to_bulk_local_norm  = dist_to_bulk_local / bulk_baseline_local
    ) |>
    left_join(meta_sub, by = "sampleid")
  
  spatial_entry <- sample_level |>
    group_by(Trial, Entry) |>
    summarise(
      mean_x = mean(.data[[x_col]], na.rm = TRUE),
      mean_y = mean(.data[[y_col]], na.rm = TRUE),
      spatial_dispersion = mean(
        sqrt((.data[[x_col]] - mean(.data[[x_col]], na.rm = TRUE))^2 +
             (.data[[y_col]] - mean(.data[[y_col]], na.rm = TRUE))^2),
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  entry_level <- sample_level |>
    group_by(Trial, Entry) |>
    summarise(
      host_selection_global = mean(dist_to_bulk_global, na.rm = TRUE),
      host_selection_local  = mean(dist_to_bulk_local, na.rm = TRUE),
      host_selection_global_norm = mean(dist_to_bulk_global_norm, na.rm = TRUE),
      host_selection_local_norm  = mean(dist_to_bulk_local_norm, na.rm = TRUE),
      mean_bulk_baseline_local = mean(bulk_baseline_local, na.rm = TRUE),
      sd_global = sd(dist_to_bulk_global, na.rm = TRUE),
      sd_local  = sd(dist_to_bulk_local, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) |>
    left_join(spatial_entry, by = c("Trial", "Entry"))
  
  list(
    sample_level = sample_level,
    entry_level = entry_level,
    filtered_mat = mat_sub,
    feature_filter = data.frame(
      Trial = panel_name,
      n_samples = ncol(mat_sub),
      n_features_before = n_features_before,
      n_features_after = n_features_after,
      min_prevalence = min_prevalence,
      min_total_count = min_total_count
    ),
    dist_mat = dist_mat
  )
}


res_bac_spring <- calc_host_selection(bac, meta4all, "SpringPan", k_nearest_bulk = 6)
res_bac_winter <- calc_host_selection(bac, meta4all, "WinterPan", k_nearest_bulk = 6)

res_fun_spring <- calc_host_selection(fun, meta4all, "SpringPan", k_nearest_bulk = 6)
res_fun_winter <- calc_host_selection(fun, meta4all, "WinterPan", k_nearest_bulk = 6)


host_selection_feature_filter_summary <- bind_rows(
  res_bac_spring$feature_filter |> mutate(Community = "Bacteria"),
  res_bac_winter$feature_filter |> mutate(Community = "Bacteria"),
  res_fun_spring$feature_filter |> mutate(Community = "Fungi"),
  res_fun_winter$feature_filter |> mutate(Community = "Fungi")
) |>
  select(Community, everything())

host_selection_feature_filter_summary

write.csv(
  host_selection_feature_filter_summary,
  "host_selection_feature_filter_summary.csv",
  row.names = FALSE
)


### helper function to add per-entry short red line
make_baseline_segments <- function(df) {
  df_sum <- df |>
    group_by(Entry) |>
    summarise(
      bulk_baseline_local = mean(bulk_baseline_local, na.rm = TRUE),
      .groups = "drop"
    )
  
  x_levels <- unique(df$Entry)
  df_sum$xnum <- match(df_sum$Entry, x_levels)
  df_sum
}

###bac spring###
m_ss_bac_spring <- aov(dist_to_bulk_global~Entry,data=res_bac_spring$sample_level)
m_ss_bac_spring |> summary()

res_bac_spring$sample_level |> 
ggplot(aes(Entry,dist_to_bulk_global))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_hline(yintercept = unique(res_bac_spring$sample_level$bulk_baseline_global),linetype="dashed")+
  theme_bw()

  
## bac winter ###
m_ss_bac_winter<- aov(dist_to_bulk_global~Entry,data=res_bac_winter$sample_level)
m_ss_bac_winter |> summary()

res_bac_winter$sample_level |> 
  ggplot(aes(Entry,dist_to_bulk_global))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_hline(yintercept = unique(res_bac_winter$sample_level$bulk_baseline_global),linetype="dashed")
theme_bw()

### fun spring ###
m_ss_fun_spring <- aov(dist_to_bulk_global~Entry,data=res_fun_spring$sample_level)
m_ss_fun_spring |> summary()

res_fun_spring$sample_level |> 
  ggplot(aes(Entry,dist_to_bulk_global))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_hline(yintercept = unique(res_fun_spring$sample_level$bulk_baseline),linetype="dashed")+
theme_bw()

### fun winter ###
m_ss_fun_winter<- aov(dist_to_bulk_global~Entry,data=res_fun_winter$sample_level)
m_ss_fun_winter |> summary()

res_fun_winter$sample_level |> 
  ggplot(aes(Entry,dist_to_bulk_global))+
  geom_boxplot(outlier.shape = NA)+
  geom_point()+
  geom_hline(yintercept = unique(res_fun_winter$sample_level$bulk_baseline),linetype="dashed")+
theme_bw()


### local####

###bac spring###
m_ss_bac_spring <- aov(dist_to_bulk_local ~ Entry, data = res_bac_spring$sample_level)
summary(m_ss_bac_spring)

seg_bac_spring <- make_baseline_segments(res_bac_spring$sample_level)

res_bac_spring$sample_level |>
  ggplot(aes(Entry, dist_to_bulk_local)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.15, height = 0)) +
  geom_segment(
    data = seg_bac_spring,
    aes(
      x = xnum - 0.28,
      xend = xnum + 0.28,
      y = bulk_baseline_local,
      yend = bulk_baseline_local
    ),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 1
  ) +
  theme_bw()


## bac winter ###
m_ss_bac_winter <- aov(dist_to_bulk_local ~ Entry, data = res_bac_winter$sample_level)
summary(m_ss_bac_winter)

seg_bac_winter <- make_baseline_segments(res_bac_winter$sample_level)

res_bac_winter$sample_level |>
  ggplot(aes(Entry, dist_to_bulk_local)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.15, height = 0)) +
  geom_segment(
    data = seg_bac_winter,
    aes(
      x = xnum - 0.28,
      xend = xnum + 0.28,
      y = bulk_baseline_local,
      yend = bulk_baseline_local
    ),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 1
  ) +
  theme_bw()


### fun spring ###
m_ss_fun_spring <- aov(dist_to_bulk_local ~ Entry, data = res_fun_spring$sample_level)
summary(m_ss_fun_spring)

seg_fun_spring <- make_baseline_segments(res_fun_spring$sample_level)

res_fun_spring$sample_level |>
  ggplot(aes(Entry, dist_to_bulk_local)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.15, height = 0)) +
  geom_segment(
    data = seg_fun_spring,
    aes(
      x = xnum - 0.28,
      xend = xnum + 0.28,
      y = bulk_baseline_local,
      yend = bulk_baseline_local
    ),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 1
  ) +
  theme_bw()

### fun winter ###
m_ss_fun_winter <- aov(dist_to_bulk_local ~ Entry, data = res_fun_winter$sample_level)
summary(m_ss_fun_winter)

seg_fun_winter <- make_baseline_segments(res_fun_winter$sample_level)

res_fun_winter$sample_level |>
  ggplot(aes(Entry, dist_to_bulk_local)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.15, height = 0)) +
  geom_segment(
    data = seg_fun_winter,
    aes(
      x = xnum - 0.28,
      xend = xnum + 0.28,
      y = bulk_baseline_local,
      yend = bulk_baseline_local
    ),
    inherit.aes = FALSE,
    color = "red",
    linewidth = 1
  ) +
  theme_bw()


calc_beta_dispersion <- function(mat, meta, panel_name,
                                 min_depth = 3000,
                                 min_prevalence = feature_min_prevalence,
                                 min_total_count = feature_min_total_count,
                                 method = "bray",
                                 min_reps = 2) {
  
  meta_sub <- meta |>
    dplyr::filter(Trial == panel_name)
  
  keep_samples <- intersect(colnames(mat), meta_sub$sampleid)
  meta_sub <- meta_sub |>
    dplyr::filter(sampleid %in% keep_samples)
  
  mat_sub <- mat[, meta_sub$sampleid, drop = FALSE]
  
  # filter samples by sequencing depth
  mat_sub <- mat_sub[, colSums(mat_sub) > min_depth, drop = FALSE]
  meta_sub <- meta_sub |>
    dplyr::filter(sampleid %in% colnames(mat_sub))
  meta_sub <- meta_sub[match(colnames(mat_sub), meta_sub$sampleid), , drop = FALSE]
  
  # keep rhizosphere only
  meta_rhizo <- meta_sub |>
    dplyr::filter(Entry != "Bulk_Soil")
  mat_rhizo <- mat_sub[, meta_rhizo$sampleid, drop = FALSE]
  meta_rhizo <- meta_rhizo[match(colnames(mat_rhizo), meta_rhizo$sampleid), , drop = FALSE]
  
  # filter low-prevalence and low-abundance taxa
  n_features_before <- nrow(mat_rhizo)
  mat_rhizo <- filter_low_abundance_features(
    mat_rhizo,
    min_prevalence = min_prevalence,
    min_total_count = min_total_count
  )
  n_features_after <- nrow(mat_rhizo)
  
  # remove zero-sum samples
  mat_rhizo <- mat_rhizo[, colSums(mat_rhizo) > 0, drop = FALSE]
  meta_rhizo <- meta_rhizo |>
    dplyr::filter(sampleid %in% colnames(mat_rhizo))
  meta_rhizo <- meta_rhizo[match(colnames(mat_rhizo), meta_rhizo$sampleid), , drop = FALSE]
  
  # keep entries with enough replicates
  entry_n <- table(meta_rhizo$Entry)
  keep_entries <- names(entry_n[entry_n >= min_reps])
  
  meta_rhizo <- meta_rhizo |>
    dplyr::filter(Entry %in% keep_entries)
  mat_rhizo <- mat_rhizo[, meta_rhizo$sampleid, drop = FALSE]
  meta_rhizo <- meta_rhizo[match(colnames(mat_rhizo), meta_rhizo$sampleid), , drop = FALSE]
  
  if (ncol(mat_rhizo) == 0) stop("No rhizosphere samples remaining in panel: ", panel_name)
  if (nrow(mat_rhizo) == 0) stop("No taxa remaining in panel: ", panel_name)
  
  mat_rhizo <- as.matrix(mat_rhizo)
  storage.mode(mat_rhizo) <- "numeric"
  
  mat_rel <- vegan::decostand(mat_rhizo, method = "total", MARGIN = 2)
  dist_rhizo <- vegan::vegdist(t(mat_rel), method = method)
  
  bd <- vegan::betadisper(dist_rhizo, group = meta_rhizo$Entry)
  
  sample_level <- data.frame(
    sampleid = names(bd$distances),
    beta_dispersion = as.numeric(bd$distances)
  ) |>
    left_join(meta_rhizo, by = "sampleid")
  
  entry_level <- sample_level |>
    group_by(Trial, Entry) |>
    summarise(
      beta_dispersion = mean(beta_dispersion, na.rm = TRUE),
      sd = sd(beta_dispersion, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    )
  
  list(
    sample_level = sample_level,
    entry_level = entry_level,
    filtered_mat = mat_rhizo,
    feature_filter = data.frame(
      Trial = panel_name,
      n_samples = ncol(mat_rhizo),
      n_features_before = n_features_before,
      n_features_after = n_features_after,
      min_prevalence = min_prevalence,
      min_total_count = min_total_count
    ),
    betadisper_obj = bd
  )
}


res_disp_bac_spring <- calc_beta_dispersion(bac, meta4all, "SpringPan")
res_disp_bac_winter <- calc_beta_dispersion(bac, meta4all, "WinterPan")

res_disp_fun_spring <- calc_beta_dispersion(fun, meta4all, "SpringPan")
res_disp_fun_winter <- calc_beta_dispersion(fun, meta4all, "WinterPan")


beta_dispersion_feature_filter_summary <- bind_rows(
  res_disp_bac_spring$feature_filter |> mutate(Community = "Bacteria"),
  res_disp_bac_winter$feature_filter |> mutate(Community = "Bacteria"),
  res_disp_fun_spring$feature_filter |> mutate(Community = "Fungi"),
  res_disp_fun_winter$feature_filter |> mutate(Community = "Fungi")
) |>
  select(Community, everything())

beta_dispersion_feature_filter_summary

write.csv(
  beta_dispersion_feature_filter_summary,
  "beta_dispersion_feature_filter_summary.csv",
  row.names = FALSE
)


df_cor_bac_spring <- res_bac_spring$sample_level |>
  select(sampleid, dist_to_bulk_local, dist_to_bulk_local_norm, dist_to_bulk_global, dist_to_bulk_global_norm) |>
  inner_join(
    res_disp_bac_spring$sample_level |>
      select(sampleid, beta_dispersion),
    by = "sampleid"
  )

df_cor_bac_winter <- res_bac_winter$sample_level |>
  select(sampleid, dist_to_bulk_local, dist_to_bulk_local_norm, dist_to_bulk_global, dist_to_bulk_global_norm) |>
  inner_join(
    res_disp_bac_winter$sample_level |>
      select(sampleid, beta_dispersion),
    by = "sampleid"
  )

df_cor_fun_spring <- res_fun_spring$sample_level |>
  select(sampleid, dist_to_bulk_local, dist_to_bulk_local_norm, dist_to_bulk_global, dist_to_bulk_global_norm) |>
  inner_join(
    res_disp_fun_spring$sample_level |>
      select(sampleid, beta_dispersion),
    by = "sampleid"
  )

df_cor_fun_winter <- res_fun_winter$sample_level |>
  select(sampleid, dist_to_bulk_local, dist_to_bulk_local_norm, dist_to_bulk_global, dist_to_bulk_global_norm) |>
  inner_join(
    res_disp_fun_winter$sample_level |>
      select(sampleid, beta_dispersion),
    by = "sampleid"
  )

cor.test(df_cor_bac_spring$dist_to_bulk_local, df_cor_bac_spring$beta_dispersion, method = "spearman")
cor.test(df_cor_bac_winter$dist_to_bulk_local, df_cor_bac_winter$beta_dispersion, method = "spearman")
cor.test(df_cor_fun_spring$dist_to_bulk_local, df_cor_fun_spring$beta_dispersion, method = "spearman")
cor.test(df_cor_fun_winter$dist_to_bulk_local, df_cor_fun_winter$beta_dispersion, method = "spearman")



cor.test(df_cor_bac_spring$dist_to_bulk_local_norm, df_cor_bac_spring$beta_dispersion, method = "spearman")
cor.test(df_cor_bac_winter$dist_to_bulk_local_norm, df_cor_bac_winter$beta_dispersion, method = "spearman")
cor.test(df_cor_fun_spring$dist_to_bulk_local_norm, df_cor_fun_spring$beta_dispersion, method = "spearman")
cor.test(df_cor_fun_winter$dist_to_bulk_local_norm, df_cor_fun_winter$beta_dispersion, method = "spearman")


cor.test(res_bac_spring$entry_level$host_selection_global,
         res_disp_bac_spring$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_bac_winter$entry_level$host_selection_global,
         res_disp_bac_winter$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_fun_spring$entry_level$host_selection_global,
         res_disp_fun_spring$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_fun_winter$entry_level$host_selection_global,
         res_disp_fun_winter$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_bac_spring$entry_level$host_selection_local,
         res_disp_bac_spring$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_bac_winter$entry_level$host_selection_local,
         res_disp_bac_winter$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_fun_spring$entry_level$host_selection_local,
         res_disp_fun_spring$entry_level$beta_dispersion,
         method = "spearman")

cor.test(res_fun_winter$entry_level$host_selection_local,
         res_disp_fun_winter$entry_level$beta_dispersion,
         method = "spearman")


make_lmm_df <- function(selection_res, dispersion_res, community, trial) {
  selection_res$sample_level |>
    select(sampleid, Entry, dist_to_bulk_global, dist_to_bulk_local,
           dist_to_bulk_global_norm, dist_to_bulk_local_norm) |>
    inner_join(
      dispersion_res$sample_level |>
        select(sampleid, beta_dispersion),
      by = "sampleid"
    ) |>
    mutate(
      Community = community,
      Trial = trial,
      panel = paste(Community, Trial, sep = "_"),
      Entry = factor(Entry)
    ) |>
    filter(
      is.finite(dist_to_bulk_local),
      is.finite(beta_dispersion),
      !is.na(Entry)
    )
}

df_lmm_raw <- bind_rows(
  make_lmm_df(res_bac_spring, res_disp_bac_spring, "Bacteria", "Spring"),
  make_lmm_df(res_bac_winter, res_disp_bac_winter, "Bacteria", "Winter"),
  make_lmm_df(res_fun_spring, res_disp_fun_spring, "Fungi", "Spring"),
  make_lmm_df(res_fun_winter, res_disp_fun_winter, "Fungi", "Winter")
) |>
  mutate(
    Community = factor(Community, levels = c("Bacteria", "Fungi")),
    Trial = factor(Trial, levels = c("Spring", "Winter")),
    panel = factor(panel)
  )

flag_panel_outliers_mahal <- function(df,
                                      vars = c("dist_to_bulk_local", "beta_dispersion"),
                                      cutoff_prob = 0.975) {
  df |>
    group_by(Community, Trial) |>
    group_modify(function(d, key) {
      x <- d |> select(all_of(vars))
      keep_complete <- complete.cases(x)

      d$mahal_dist <- NA_real_
      d$outlier_mahal <- FALSE

      if (sum(keep_complete) > length(vars) + 2) {
        x_complete <- x[keep_complete, , drop = FALSE]
        md <- mahalanobis(
          x_complete,
          center = colMeans(x_complete),
          cov = cov(x_complete)
        )
        cutoff <- qchisq(cutoff_prob, df = length(vars))

        d$mahal_dist[keep_complete] <- md
        d$outlier_mahal[keep_complete] <- md > cutoff
      }

      d
    }) |>
    ungroup()
}

df_lmm_raw_flagged <- df_lmm_raw |>
  flag_panel_outliers_mahal(
    vars = c("dist_to_bulk_local", "beta_dispersion"),
    cutoff_prob = 0.99
  )

df_lmm_raw_no_outlier <- df_lmm_raw_flagged |>
  filter(!outlier_mahal)

df_lmm_z <- df_lmm_raw |>
  group_by(Community, Trial) |>
  mutate(
    dist_to_bulk_local_z = as.numeric(scale(dist_to_bulk_local)),
    beta_dispersion_z = as.numeric(scale(beta_dispersion))
  ) |>
  ungroup()

df_lmm_z_no_outlier <- df_lmm_raw_no_outlier |>
  group_by(Community, Trial) |>
  mutate(
    dist_to_bulk_local_z = as.numeric(scale(dist_to_bulk_local)),
    beta_dispersion_z = as.numeric(scale(beta_dispersion))
  ) |>
  ungroup()

df_lmm_raw |>
  count(Community, Trial, Entry, name = "n_samples") |>
  arrange(Community, Trial, Entry)


outlier_summary <- df_lmm_raw_flagged |>
  count(Community, Trial, outlier_mahal, name = "n") |>
  arrange(Community, Trial, outlier_mahal)

outlier_samples <- df_lmm_raw_flagged |>
  filter(outlier_mahal) |>
  arrange(Community, Trial, desc(mahal_dist)) |>
  select(
    Community, Trial, sampleid, Entry,
    dist_to_bulk_local, beta_dispersion,
    mahal_dist
  )

outlier_summary
outlier_samples

write.csv(outlier_samples, "host_selection_vs_dispersion_flagged_outliers.csv", row.names = FALSE)


fit_lmm_by_panel <- function(df, response, predictor) {
  panel_data <- df |>
    group_by(Community, Trial) |>
    group_split()

  names(panel_data) <- map_chr(
    panel_data,
    ~ paste(unique(.x$Community), unique(.x$Trial), sep = "_")
  )

  map(panel_data, function(dat) {
    formula <- as.formula(paste(response, "~", predictor, "+ (1 | Entry)"))
    lmer(formula, data = dat, REML = TRUE)
  })
}

extract_lmm_slope <- function(models, predictor, scale_label) {
  imap_dfr(models, function(mod, panel_name) {
    coef_tab <- summary(mod)$coefficients
    slope <- coef_tab[predictor, , drop = FALSE]
    panel_info <- str_split_fixed(panel_name, "_", 2)

    data.frame(
      scale = scale_label,
      Community = panel_info[, 1],
      Trial = panel_info[, 2],
      predictor = predictor,
      estimate = slope[, "Estimate"],
      std_error = slope[, "Std. Error"],
      df = slope[, "df"],
      t_value = slope[, "t value"],
      p_value = slope[, "Pr(>|t|)"],
      n = nobs(mod),
      singular = isSingular(mod)
    )
  })
}

m_raw <- fit_lmm_by_panel(
  df = df_lmm_raw,
  response = "beta_dispersion",
  predictor = "dist_to_bulk_local"
)

m_raw_no_outlier <- fit_lmm_by_panel(
  df = df_lmm_raw_no_outlier,
  response = "beta_dispersion",
  predictor = "dist_to_bulk_local"
)

lmm_raw_summary <- extract_lmm_slope(
  models = m_raw,
  predictor = "dist_to_bulk_local",
  scale_label = "raw full data"
)

lmm_raw_no_outlier_summary <- extract_lmm_slope(
  models = m_raw_no_outlier,
  predictor = "dist_to_bulk_local",
  scale_label = "raw no outlier"
)

lmm_raw_summary
lmm_raw_no_outlier_summary

summary(m_raw$Bacteria_Spring)
summary(m_raw$Bacteria_Winter)
summary(m_raw$Fungi_Spring)
summary(m_raw$Fungi_Winter)


m_z <- fit_lmm_by_panel(
  df = df_lmm_z,
  response = "beta_dispersion_z",
  predictor = "dist_to_bulk_local_z"
)

m_z_no_outlier <- fit_lmm_by_panel(
  df = df_lmm_z_no_outlier,
  response = "beta_dispersion_z",
  predictor = "dist_to_bulk_local_z"
)

lmm_z_summary <- extract_lmm_slope(
  models = m_z,
  predictor = "dist_to_bulk_local_z",
  scale_label = "z full data"
)

lmm_z_no_outlier_summary <- extract_lmm_slope(
  models = m_z_no_outlier,
  predictor = "dist_to_bulk_local_z",
  scale_label = "z no outlier"
)

lmm_z_summary
lmm_z_no_outlier_summary

summary(m_z$Bacteria_Spring)
summary(m_z$Bacteria_Winter)
summary(m_z$Fungi_Spring)
summary(m_z$Fungi_Winter)


lmm_summary <- bind_rows(
  lmm_raw_summary,
  lmm_raw_no_outlier_summary,
  lmm_z_summary,
  lmm_z_no_outlier_summary
) |>
  arrange(Community, Trial, scale)

lmm_summary

write.csv(lmm_summary, "host_selection_vs_dispersion_lmm_outlier_sensitivity_summary.csv", row.names = FALSE)


make_lmm_prediction_df <- function(models, predictor, scale_label) {
  imap_dfr(models, function(mod, panel_name) {
    panel_info <- str_split_fixed(panel_name, "_", 2)
    ggpredict(mod, terms = predictor) |>
      as.data.frame() |>
      mutate(
        Community = panel_info[, 1],
        Trial = panel_info[, 2],
        scale = scale_label
      )
  })
}

pred_raw <- make_lmm_prediction_df(
  models = m_raw,
  predictor = "dist_to_bulk_local",
  scale_label = "raw full data"
)

pred_raw_no_outlier <- make_lmm_prediction_df(
  models = m_raw_no_outlier,
  predictor = "dist_to_bulk_local",
  scale_label = "raw no outlier"
)

p_lmm_raw <- ggplot(df_lmm_raw_flagged, aes(x = dist_to_bulk_local, y = beta_dispersion)) +
  geom_point(aes(colour = outlier_mahal), alpha = 0.45, size = 1.4) +
  geom_ribbon(
    data = pred_raw,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_raw,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  scale_colour_manual(
    values = c("FALSE" = "grey45", "TRUE" = "#D73027"),
    labels = c("FALSE" = "Kept", "TRUE" = "Flagged outlier"),
    name = NULL
  ) +
  facet_grid(Community ~ Trial, scales = "free") +
  theme_bw() +
  labs(
    title = "Raw-scale LMM: full data",
    x = "Local selection strength",
    y = "Distance to genotype centroid"
  )

p_lmm_raw

p_lmm_raw_no_outlier <- ggplot(df_lmm_raw_no_outlier, aes(x = dist_to_bulk_local, y = beta_dispersion)) +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_ribbon(
    data = pred_raw_no_outlier,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_raw_no_outlier,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  facet_grid(Community ~ Trial, scales = "free") +
  theme_bw() +
  labs(
    title = "Raw-scale LMM: flagged outliers removed",
    x = "Local selection strength",
    y = "Distance to genotype centroid"
  )

p_lmm_raw_no_outlier


pred_z <- make_lmm_prediction_df(
  models = m_z,
  predictor = "dist_to_bulk_local_z",
  scale_label = "z full data"
)

pred_z_no_outlier <- make_lmm_prediction_df(
  models = m_z_no_outlier,
  predictor = "dist_to_bulk_local_z",
  scale_label = "z no outlier"
)

p_lmm_z <- ggplot(
  df_lmm_z |>
    left_join(
      df_lmm_raw_flagged |> select(sampleid, Community, Trial, outlier_mahal),
      by = c("sampleid", "Community", "Trial")
    ),
  aes(x = dist_to_bulk_local_z, y = beta_dispersion_z)
) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(aes(colour = outlier_mahal), alpha = 0.45, size = 1.4) +
  geom_ribbon(
    data = pred_z,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_z,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  scale_colour_manual(
    values = c("FALSE" = "grey45", "TRUE" = "#D73027"),
    labels = c("FALSE" = "Kept", "TRUE" = "Flagged outlier"),
    name = NULL
  ) +
  facet_grid(Community ~ Trial) +
  theme_bw() +
  labs(
    title = "Within-panel standardized LMM: full data",
    x = "Local selection strength (z-score)",
    y = "Distance to genotype centroid (z-score)"
  )

p_lmm_z

p_lmm_z_no_outlier <- ggplot(df_lmm_z_no_outlier, aes(x = dist_to_bulk_local_z, y = beta_dispersion_z)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(alpha = 0.35, size = 1.2) +
  geom_ribbon(
    data = pred_z_no_outlier,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    inherit.aes = FALSE,
    alpha = 0.2
  ) +
  geom_line(
    data = pred_z_no_outlier,
    aes(x = x, y = predicted),
    inherit.aes = FALSE,
    linewidth = 1
  ) +
  facet_grid(Community ~ Trial) +
  theme_bw() +
  labs(
    title = "Within-panel standardized LMM: flagged outliers removed",
    x = "Local selection strength (z-score)",
    y = "Distance to genotype centroid (z-score)"
  )

p_lmm_z_no_outlier


p_data_filter_norm <- (p_lmm_raw / p_lmm_raw_no_outlier) | (p_lmm_z / p_lmm_z_no_outlier)
ggsave("FigS_host_selection_vs_dispersion_lmm_outlier_sensitivity.pdf", p_data_filter_norm, height = 12, width = 10)
ggsave("FigS_host_selection_vs_dispersion_lmm_outlier_sensitivity.png", p_data_filter_norm, height = 12, width = 10, dpi = 300)


library(ggthemes)
p_lmm_z_no_outlier <- p_lmm_z_no_outlier+
  theme_few()
p_lmm_z_no_outlier

ggsave("host_selection_vs_dispersion.pdf",p_lmm_z_no_outlier,height=8,width=8 )


make_selection_dispersion_explainer <- function(selection_res,
                                                dispersion_res,
                                                community = "Bacteria",
                                                season = "Winter",
                                                entry_to_highlight = NULL,
                                                k_nearest_bulk = 6) {
  dist_mat <- selection_res$dist_mat
  sample_ids <- rownames(dist_mat)

  rhizo_meta <- selection_res$sample_level |>
    select(
      sampleid, Entry, Trial, Col, Row,
      dist_to_bulk_local, bulk_baseline_local
    ) |>
    distinct()

  bulk_meta <- meta4all |>
    filter(sampleid %in% sample_ids, Entry == "Bulk_Soil") |>
    select(sampleid, Entry, Trial, Col, Row) |>
    distinct()

  if (is.null(entry_to_highlight)) {
    entry_to_highlight <- rhizo_meta |>
      count(Entry, sort = TRUE) |>
      slice(1) |>
      pull(Entry)
  }

  highlight_roots_meta <- rhizo_meta |>
    filter(Entry == entry_to_highlight)

  nearest_bulk_links <- highlight_roots_meta |>
    group_by(sampleid) |>
    group_modify(function(root_sample, key) {
      root_id <- key$sampleid[[1]]
      root_col <- root_sample$Col[[1]]
      root_row <- root_sample$Row[[1]]

        bulk_ranked <- bulk_meta |>
        mutate(
          d_xy = sqrt((Col - root_col)^2 + (Row - root_row)^2)
        ) |>
        arrange(d_xy) |>
        slice_head(n = k_nearest_bulk)

      bulk_ranked |>
        transmute(
          root_sampleid = root_id,
          bulk_sampleid = sampleid,
          x = root_col,
          y = root_row,
          xend = Col,
          yend = Row
        )
    }) |>
    ungroup()

  linked_bulk_ids <- unique(nearest_bulk_links$bulk_sampleid)

  p_xy <- ggplot() +
    geom_point(
      data = bulk_meta,
      aes(Col, Row),
      colour = "#9ECAE1",
      alpha = 0.35,
      size = 1.4
    ) +
    geom_segment(
      data = nearest_bulk_links,
      aes(x = x, y = y, xend = xend, yend = yend),
      colour = "#2166AC",
      alpha = 0.25,
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_point(
      data = bulk_meta |> filter(sampleid %in% linked_bulk_ids),
      aes(Col, Row),
      colour = "#2166AC",
      size = 2.2
    ) +
    geom_point(
      data = highlight_roots_meta,
      aes(Col, Row),
      colour = "#D73027",
      shape = 17,
      size = 2.2
    ) +
    coord_equal() +
    theme_bw() +
    labs(
      title = "Field XY map",
      subtitle = "Blue dashed lines connect each root to its 6 nearest local bulk soil plots",
      x = "Column",
      y = "Row"
    )

  plot_meta <- bind_rows(
    rhizo_meta |> select(sampleid, Entry, Trial),
    bulk_meta |> select(sampleid, Entry, Trial)
  ) |>
    right_join(data.frame(sampleid = sample_ids), by = "sampleid") |>
    mutate(
      Entry = if_else(is.na(Entry), "Bulk_Soil", Entry),
      Sample_type = if_else(Entry == "Bulk_Soil", "Bulk soil", "Root"),
      group = case_when(
        Entry == entry_to_highlight & Sample_type == "Root" ~ "Highlighted genotype",
        sampleid %in% linked_bulk_ids ~ "Linked local bulk soil",
        Sample_type == "Bulk soil" ~ "Other bulk soil",
        TRUE ~ "Other root"
      )
    )

  pcoa <- cmdscale(as.dist(dist_mat), eig = TRUE, k = 2)

  pcoa_df <- data.frame(
    sampleid = rownames(pcoa$points),
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2]
  ) |>
    left_join(plot_meta, by = "sampleid")

  highlight_roots_pcoa <- pcoa_df |>
    filter(
      Entry == entry_to_highlight,
      Sample_type == "Root"
    )

  root_coords <- highlight_roots_pcoa |>
    select(
      root_sampleid = sampleid,
      root_PCoA1 = PCoA1,
      root_PCoA2 = PCoA2
    )

  bulk_coords <- pcoa_df |>
    filter(Sample_type == "Bulk soil") |>
    select(
      bulk_sampleid = sampleid,
      bulk_PCoA1 = PCoA1,
      bulk_PCoA2 = PCoA2
    )

  root_bulk_segments_pcoa <- nearest_bulk_links |>
    select(root_sampleid, bulk_sampleid) |>
    distinct() |>
    inner_join(root_coords, by = "root_sampleid") |>
    inner_join(bulk_coords, by = "bulk_sampleid")

  centroid_df <- highlight_roots_pcoa |>
    summarise(
      PCoA1_centroid = mean(PCoA1, na.rm = TRUE),
      PCoA2_centroid = mean(PCoA2, na.rm = TRUE)
    )

  centroid_segments <- highlight_roots_pcoa |>
    mutate(
      PCoA1_centroid = centroid_df$PCoA1_centroid,
      PCoA2_centroid = centroid_df$PCoA2_centroid
    )

  eig <- pcoa$eig
  var_exp <- eig / sum(eig[eig > 0])

  p_pcoa_explain <- ggplot() +
    geom_point(
      data = pcoa_df |> filter(group == "Other root"),
      aes(PCoA1, PCoA2),
      colour = "grey70",
      alpha = 0.35,
      size = 1.2
    ) +
    geom_point(
      data = pcoa_df |> filter(group == "Other bulk soil"),
      aes(PCoA1, PCoA2),
      colour = "#9ECAE1",
      alpha = 0.35,
      size = 1.3
    ) +
    geom_segment(
      data = root_bulk_segments_pcoa,
      aes(
        x = root_PCoA1,
        y = root_PCoA2,
        xend = bulk_PCoA1,
        yend = bulk_PCoA2
      ),
      colour = "#2166AC",
      linewidth = 0.4,
      alpha = 0.45,
      linetype = "dashed"
    ) +
    geom_point(
      data = pcoa_df |> filter(group == "Linked local bulk soil"),
      aes(PCoA1, PCoA2),
      colour = "#2166AC",
      size = 2.1
    ) +
    geom_segment(
      data = centroid_segments,
      aes(
        x = PCoA1,
        y = PCoA2,
        xend = PCoA1_centroid,
        yend = PCoA2_centroid
      ),
      colour = "#D73027",
      shap=17,
      linewidth = 0.55,
      alpha = 0.65
    ) +
    geom_point(
      data = highlight_roots_pcoa,
      aes(PCoA1, PCoA2),
      colour = "#D73027",
      shape = 17,
      size = 2
    ) +
    geom_point(
      data = centroid_df,
      aes(PCoA1_centroid, PCoA2_centroid),
      colour = "#D73027",
      fill = "white",
      shape = 8,
      size = 3.4,
      stroke = 1
    ) +
    theme_bw() +
    labs(
      title = "Community PCoA",
      subtitle = "Blue dashed lines connect each highlighted root to its 6 nearest bulk soils; red lines show distance to genotype centroid",
      x = paste0("PCoA1 (", round(var_exp[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(var_exp[2] * 100, 1), "%)")
    )

  dispersion_values <- dispersion_res$sample_level |>
    filter(Entry == entry_to_highlight) |>
    summarise(value = mean(beta_dispersion, na.rm = TRUE)) |>
    mutate(metric = "Distance to genotype centroid")

  selection_values <- highlight_roots_meta |>
    summarise(value = mean(dist_to_bulk_local, na.rm = TRUE)) |>
    mutate(metric = "Local host selection strength")

  value_df <- bind_rows(selection_values, dispersion_values) |>
    mutate(metric = factor(
      metric,
      levels = c(
        "Local host selection strength",
        "Distance to genotype centroid"
      )
    ))

  p_values <- ggplot(value_df, aes(metric, value, fill = metric)) +
    geom_col(width = 0.65, colour = "grey30") +
    scale_fill_manual(values = c(
      "Local host selection strength" = "#2166AC",
      "Distance to genotype centroid" = "#D73027"
    )) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 35, hjust = 1)
    ) +
    labs(
      title = "Derived values",
      subtitle = "Mean values for highlighted genotype",
      x = NULL,
      y = "Community distance"
    )

  (p_xy | p_pcoa_explain) +
    plot_layout(widths = c(1, 1.15)) +
    plot_annotation(
      title = paste0(community, " ", season, ": ", entry_to_highlight),
      subtitle = "How local host selection strength and distance to genotype centroid were derived"
    )
}

p_explain_bac_spring <- make_selection_dispersion_explainer(
  selection_res = res_bac_spring,
  dispersion_res = res_disp_bac_spring,
  community = "Bacteria",
  season = "Spring",
  entry_to_highlight = NULL,
  k_nearest_bulk = 6
)

p_explain_bac_winter <- make_selection_dispersion_explainer(
  selection_res = res_bac_winter,
  dispersion_res = res_disp_bac_winter,
  community = "Bacteria",
  season = "Winter",
  entry_to_highlight = NULL,
  k_nearest_bulk = 6
)

p_explain_fun_spring <- make_selection_dispersion_explainer(
  selection_res = res_fun_spring,
  dispersion_res = res_disp_fun_spring,
  community = "Fungi",
  season = "Spring",
  entry_to_highlight = NULL,
  k_nearest_bulk = 6
)

p_explain_fun_winter <- make_selection_dispersion_explainer(
  selection_res = res_fun_winter,
  dispersion_res = res_disp_fun_winter,
  community = "Fungi",
  season = "Winter",
  entry_to_highlight = NULL,
  k_nearest_bulk = 6
)

p_explain_four_panels <-
  p_explain_bac_spring /
  p_explain_bac_winter /
  p_explain_fun_spring /
  p_explain_fun_winter

p_explain_four_panels # ALANDA highlighted as an example genotype in all panels since it had the most samples across all panels

ggsave("FigS_selection_dispersion_explainer.pdf", p_explain_four_panels, height = 16, width = 10)
ggsave("FigS_selection_dispersion_explainer.png", p_explain_four_panels, height = 16, width = 10, dpi = 300)
