# Figures: Fig3a
# Original source: BMB_bulk/root_bulk_soil_integrated_associations.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

rm(list = ls())

load("../pcnm_VPA.rdata")

library(qiime2R)
library(tidyverse)
library(vegan)
library(linkET)

set.seed(123)

feature_min_prevalence <- 0.01
feature_min_total_abundance <- 10

bac <- read_q2biom("AUT_Bac/bulk_merged_table/feature-table.biom")
colnames(bac) <- sub("^FL_", "FL", colnames(bac))

fun <- readRDS("AUT_Fun/fun_otu_count_with_bulk.rds") |>
  t() |>
  as.data.frame()
colnames(fun) <- sub("^FL_", "FL", colnames(fun))

meta <- read.csv("../Bac_Aut/0_raw_data/meta_file.csv")
meta_aut <- meta |>
  filter(Location == "Austria", Trial != "DivPan") |>
  mutate(Plot = as.character(Plot))

env <- read.csv("AUS_ENV_20241112.csv", check.names = TRUE) |>
  mutate(Plot = as.character(Plot))

pcnm_df$sampleid <- as.character(pcnm_df$sampleid)


make_bulk_meta <- function(sample_names, meta_aut) {
  data.frame(sampleid = sample_names) |>
    filter(str_detect(sampleid, "^bulk_?\\d+$")) |>
    mutate(
      bulk_num = as.numeric(str_extract(sampleid, "\\d+$")),
      Trial = case_when(
        bulk_num >= 1 & bulk_num <= 600 ~ "WinterPan",
        bulk_num >= 601 & bulk_num <= 1200 ~ "DivPan",
        bulk_num >= 1201 & bulk_num <= 1800 ~ "SpringPan",
        TRUE ~ NA_character_
      ),
      Plot = case_when(
        bulk_num >= 1 & bulk_num <= 600 ~ bulk_num,
        bulk_num >= 601 & bulk_num <= 1200 ~ bulk_num - 600,
        bulk_num >= 1201 & bulk_num <= 1800 ~ bulk_num - 1200,
        TRUE ~ NA_real_
      ),
      Plot = as.character(Plot),
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
}

meta4bulk_bac <- make_bulk_meta(colnames(bac), meta_aut)
meta4bulk_fun <- make_bulk_meta(colnames(fun), meta_aut)


aggregate_otu_to_plot <- function(otu, meta_plot, sample_col = "sampleid") {
  otu_sub <- otu[, colnames(otu) %in% meta_plot[[sample_col]], drop = FALSE]
  meta_plot <- meta_plot[
    match(colnames(otu_sub), meta_plot[[sample_col]]),
    ,
    drop = FALSE
  ]

  stopifnot(all(colnames(otu_sub) == meta_plot[[sample_col]]))

  colnames(otu_sub) <- as.character(meta_plot$Plot)

  if (any(table(colnames(otu_sub)) > 1)) {
    plot_ids <- unique(colnames(otu_sub))
    otu_plot <- sapply(plot_ids, function(p) {
      cols <- which(colnames(otu_sub) == p)
      rowMeans(as.matrix(otu_sub[, cols, drop = FALSE]))
    })
    otu_plot <- as.matrix(otu_plot)
    colnames(otu_plot) <- plot_ids
  } else {
    otu_plot <- as.matrix(otu_sub)
  }

  otu_plot
}

filter_features <- function(otu_plot,
                            min_prevalence = feature_min_prevalence,
                            min_total_abundance = feature_min_total_abundance) {
  prevalence <- rowMeans(otu_plot > 0)
  total_abundance <- rowSums(otu_plot)

  otu_plot[
    prevalence >= min_prevalence &
      total_abundance >= min_total_abundance,
    ,
    drop = FALSE
  ]
}

make_hellinger_plot_matrix <- function(otu_plot,
                                       min_prevalence = feature_min_prevalence,
                                       min_total_abundance = feature_min_total_abundance) {
  otu_plot <- filter_features(
    otu_plot,
    min_prevalence = min_prevalence,
    min_total_abundance = min_total_abundance
  )

  comm <- otu_plot |>
    t() |>
    as.data.frame()

  rownames(comm) <- colnames(otu_plot)

  decostand(comm, method = "hellinger")
}

mantel_to_row <- function(x, data_type, trial, association, distance) {
  data.frame(
    data_type = data_type,
    trial = trial,
    association = association,
    distance = distance,
    mantel_r = as.numeric(x$statistic),
    p_value = as.numeric(x$signif)
  )
}

get_selected_pcnm <- function(object_name) {
  if (exists(object_name, inherits = TRUE)) {
    get(object_name, inherits = TRUE)
  } else {
    NULL
  }
}


run_bulk_spatial_analysis <- function(comm_mat,
                                      bulk_meta,
                                      meta_aut,
                                      pcnm_df,
                                      trial,
                                      selected_pcnm_model = NULL,
                                      pcnm_plot_source_trial = trial) {
  bulk_meta_sub <- bulk_meta |>
    filter(Trial == trial, sampleid %in% colnames(comm_mat))

  bulk_plot <- aggregate_otu_to_plot(comm_mat, bulk_meta_sub)

  bulk_meta_plot <- bulk_meta_sub |>
    mutate(Plot = as.character(Plot)) |>
    distinct(Plot, .keep_all = TRUE)
  bulk_meta_plot <- bulk_meta_plot[
    match(colnames(bulk_plot), bulk_meta_plot$Plot),
    ,
    drop = FALSE
  ]

  stopifnot(all(colnames(bulk_plot) == bulk_meta_plot$Plot))

  bulk_hellinger <- make_hellinger_plot_matrix(bulk_plot)

  meta_xy <- bulk_meta_plot |>
    select(Plot, Row, Col) |>
    mutate(Row = as.numeric(Row), Col = as.numeric(Col))
  meta_xy <- meta_xy[
    match(rownames(bulk_hellinger), meta_xy$Plot),
    ,
    drop = FALSE
  ]

  stopifnot(all(rownames(bulk_hellinger) == meta_xy$Plot))

  pcnm_root <- pcnm_df |>
    inner_join(
      meta_aut |>
        filter(Trial == pcnm_plot_source_trial) |>
        select(sampleid, Plot),
      by = "sampleid"
    ) |>
    filter(Plot %in% bulk_meta_plot$Plot)

  pcnm_plot <- pcnm_root |>
    group_by(Plot) |>
    summarise(across(starts_with("PCNM"), mean), .groups = "drop") |>
    as.data.frame()

  rownames(pcnm_plot) <- pcnm_plot$Plot
  pcnm_plot$Plot <- NULL
  pcnm_plot <- pcnm_plot[
    match(rownames(bulk_hellinger), rownames(pcnm_plot)),
    ,
    drop = FALSE
  ]

  stopifnot(all(rownames(bulk_hellinger) == rownames(pcnm_plot)))

  pcnm_mat <- pcnm_plot[, grepl("^PCNM", colnames(pcnm_plot)), drop = FALSE]
  pcnm_mat <- pcnm_mat[, apply(pcnm_mat, 2, var, na.rm = TRUE) > 0, drop = FALSE]

  mod_xy <- rda(bulk_hellinger ~ Row + Col, data = meta_xy)
  mod_pcnm_full <- rda(bulk_hellinger ~ ., data = pcnm_mat)

  if (!is.null(selected_pcnm_model)) {
    sel_vars <- attr(terms(selected_pcnm_model), "term.labels")
    sel_vars <- intersect(sel_vars, colnames(pcnm_mat))
    pcnm_sel_mat <- pcnm_mat[, sel_vars, drop = FALSE]
    mod_pcnm_sel <- rda(bulk_hellinger ~ ., data = pcnm_sel_mat)
  } else {
    pcnm_sel_mat <- NULL
    mod_pcnm_sel <- NULL
  }

  pcnm_2 <- pcnm_mat[, 1:min(2, ncol(pcnm_mat)), drop = FALSE]
  pcnm_10 <- pcnm_mat[, 1:min(10, ncol(pcnm_mat)), drop = FALSE]
  mod_pcnm2 <- rda(bulk_hellinger ~ ., data = pcnm_2)
  mod_pcnm10 <- rda(bulk_hellinger ~ ., data = pcnm_10)

  comparison <- data.frame(
    model = c("XY", "PCNM_2", "PCNM_10", "PCNM_selected"),
    R2 = c(
      RsquareAdj(mod_xy)$r.squared,
      RsquareAdj(mod_pcnm2)$r.squared,
      RsquareAdj(mod_pcnm10)$r.squared,
      if (!is.null(mod_pcnm_sel)) RsquareAdj(mod_pcnm_sel)$r.squared else NA
    ),
    adj_R2 = c(
      RsquareAdj(mod_xy)$adj.r.squared,
      RsquareAdj(mod_pcnm2)$adj.r.squared,
      RsquareAdj(mod_pcnm10)$adj.r.squared,
      if (!is.null(mod_pcnm_sel)) RsquareAdj(mod_pcnm_sel)$adj.r.squared else NA
    )
  )

  list(
    trial = trial,
    bulk_meta_plot = bulk_meta_plot,
    bulk_hellinger = bulk_hellinger,
    meta_xy = meta_xy,
    pcnm_mat = pcnm_mat,
    pcnm_sel_mat = pcnm_sel_mat,
    mod_xy = mod_xy,
    mod_pcnm_full = mod_pcnm_full,
    mod_pcnm_sel = mod_pcnm_sel,
    comparison = comparison
  )
}


bulk_bac_spring <- run_bulk_spatial_analysis(
  comm_mat = bac,
  bulk_meta = meta4bulk_bac,
  meta_aut = meta_aut,
  pcnm_df = pcnm_df,
  trial = "SpringPan",
  selected_pcnm_model = get_selected_pcnm("pcnm_sel_spring_bac"),
  pcnm_plot_source_trial = "SpringPan"
)

bulk_bac_winter <- run_bulk_spatial_analysis(
  comm_mat = bac,
  bulk_meta = meta4bulk_bac,
  meta_aut = meta_aut,
  pcnm_df = pcnm_df,
  trial = "WinterPan",
  selected_pcnm_model = get_selected_pcnm("pcnm_sel_winter_bac"),
  pcnm_plot_source_trial = "SpringPan"
)

bulk_fun_spring <- run_bulk_spatial_analysis(
  comm_mat = fun,
  bulk_meta = meta4bulk_fun,
  meta_aut = meta_aut,
  pcnm_df = pcnm_df,
  trial = "SpringPan",
  selected_pcnm_model = get_selected_pcnm("pcnm_sel_spring_fun"),
  pcnm_plot_source_trial = "SpringPan"
)

bulk_fun_winter <- run_bulk_spatial_analysis(
  comm_mat = fun,
  bulk_meta = meta4bulk_fun,
  meta_aut = meta_aut,
  pcnm_df = pcnm_df,
  trial = "WinterPan",
  selected_pcnm_model = get_selected_pcnm("pcnm_sel_winter_fun"),
  pcnm_plot_source_trial = "SpringPan"
)

bulk_pcnm_summary <- bind_rows(
  bulk_bac_spring$comparison |> mutate(data_type = "Bacteria", trial = "SpringPan"),
  bulk_bac_winter$comparison |> mutate(data_type = "Bacteria", trial = "WinterPan"),
  bulk_fun_spring$comparison |> mutate(data_type = "Fungi", trial = "SpringPan"),
  bulk_fun_winter$comparison |> mutate(data_type = "Fungi", trial = "WinterPan")
) |>
  select(data_type, trial, model, R2, adj_R2)

bulk_pcnm_summary


run_root_bulk_mantel <- function(otu_root,
                                 otu_bulk,
                                 meta_root,
                                 meta_bulk,
                                 trial,
                                 data_type,
                                 min_prevalence = feature_min_prevalence,
                                 min_total_abundance = feature_min_total_abundance,
                                 permutations = 9999) {
  meta_root2 <- meta_root |>
    filter(
      Trial == trial,
      Plot %in% meta_bulk$Plot,
      sampleid %in% colnames(otu_root)
    )

  meta_bulk2 <- meta_bulk |>
    filter(
      Trial == trial,
      sampleid %in% colnames(otu_bulk)
    )

  common_plot <- intersect(unique(meta_root2$Plot), unique(meta_bulk2$Plot))
  meta_root2 <- meta_root2 |> filter(Plot %in% common_plot)
  meta_bulk2 <- meta_bulk2 |> filter(Plot %in% common_plot)

  root_plot <- aggregate_otu_to_plot(otu_root, meta_root2)
  bulk_plot <- aggregate_otu_to_plot(otu_bulk, meta_bulk2)

  common_plot <- intersect(colnames(root_plot), colnames(bulk_plot))
  root_plot <- root_plot[, common_plot, drop = FALSE]
  bulk_plot <- bulk_plot[, common_plot, drop = FALSE]

  stopifnot(identical(colnames(root_plot), colnames(bulk_plot)))

  root_plot_filtered <- filter_features(
    root_plot,
    min_prevalence = min_prevalence,
    min_total_abundance = min_total_abundance
  )
  bulk_plot_filtered <- filter_features(
    bulk_plot,
    min_prevalence = min_prevalence,
    min_total_abundance = min_total_abundance
  )

  root_hellinger <- make_hellinger_plot_matrix(
    root_plot_filtered,
    min_prevalence = 0,
    min_total_abundance = 0
  )
  bulk_hellinger <- make_hellinger_plot_matrix(
    bulk_plot_filtered,
    min_prevalence = 0,
    min_total_abundance = 0
  )

  bc_root <- vegdist(root_hellinger, method = "bray")
  bc_bulk <- vegdist(bulk_hellinger, method = "bray")

  jaccard_root <- root_plot_filtered |>
    t() |>
    vegdist(method = "jaccard", binary = TRUE)
  jaccard_bulk <- bulk_plot_filtered |>
    t() |>
    vegdist(method = "jaccard", binary = TRUE)

  stopifnot(identical(attr(bc_root, "Labels"), attr(bc_bulk, "Labels")))
  stopifnot(identical(attr(jaccard_root, "Labels"), attr(jaccard_bulk, "Labels")))

  mantel_bc <- mantel(bc_bulk, bc_root, method = "spearman", permutations = permutations)
  mantel_jaccard <- mantel(jaccard_bulk, jaccard_root, method = "spearman", permutations = permutations)

  summary <- bind_rows(
    mantel_to_row(mantel_bc, data_type, trial, "root microbiome - bulk microbiome", "Bray-Curtis"),
    mantel_to_row(mantel_jaccard, data_type, trial, "root microbiome - bulk microbiome", "Jaccard")
  )

  list(
    trial = trial,
    data_type = data_type,
    root_hellinger = root_hellinger,
    bulk_hellinger = bulk_hellinger,
    bc_root = bc_root,
    bc_bulk = bc_bulk,
    jaccard_root = jaccard_root,
    jaccard_bulk = jaccard_bulk,
    mantel_bc = mantel_bc,
    mantel_jaccard = mantel_jaccard,
    summary = summary
  )
}


root_bulk_bac_spring <- run_root_bulk_mantel(
  otu_root = bac,
  otu_bulk = bac,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_bac,
  trial = "SpringPan",
  data_type = "Bacteria"
)

root_bulk_bac_winter <- run_root_bulk_mantel(
  otu_root = bac,
  otu_bulk = bac,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_bac,
  trial = "WinterPan",
  data_type = "Bacteria"
)

root_bulk_fun_spring <- run_root_bulk_mantel(
  otu_root = fun,
  otu_bulk = fun,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_fun,
  trial = "SpringPan",
  data_type = "Fungi"
)

root_bulk_fun_winter <- run_root_bulk_mantel(
  otu_root = fun,
  otu_bulk = fun,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_fun,
  trial = "WinterPan",
  data_type = "Fungi"
)

root_bulk_summary <- bind_rows(
  root_bulk_bac_spring$summary,
  root_bulk_bac_winter$summary,
  root_bulk_fun_spring$summary,
  root_bulk_fun_winter$summary
 ) |>
  mutate(
    pair = case_when(
      data_type == "Bacteria" ~ "Root_Bacteria - Bulk_Bacteria",
      data_type == "Fungi" ~ "Root_Fungi - Bulk_Fungi",
      TRUE ~ association
    )
)

root_bulk_summary


root_bulk_correlation_plot <- root_bulk_summary |>
  mutate(
    pair = factor(
      pair,
      levels = c("Root_Bacteria - Bulk_Bacteria", "Root_Fungi - Bulk_Fungi")
    ),
    significance = case_when(
      p_value < 0.01 ~ "p < 0.01",
      p_value < 0.05 ~ "p < 0.05",
      TRUE ~ "p >= 0.05"
    )
  ) |>
  ggplot(aes(x = pair, y = mantel_r, colour = significance)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(size = 3) +
  facet_grid(distance ~ trial) +
  scale_colour_manual(values = c("p < 0.01" = "#D73027", "p < 0.05" = "#FC8D59", "p >= 0.05" = "#4575B4")) +
  labs(
    x = NULL,
    y = "Mantel's r",
    colour = "Mantel's p"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid.minor = element_blank()
  )

root_bulk_correlation_plot


prepare_env_matrix <- function(env, trial, common_plots, env_vars = NULL) {
  env_sub <- env |>
    filter(Trial == trial, Plot %in% common_plots) |>
    arrange(match(Plot, common_plots))

  if (is.null(env_vars)) {
    env_vars <- env_sub |>
      select(where(is.numeric)) |>
      select(-any_of(c("Plot0", "Plot"))) |>
      colnames()
  }

  env_mat <- env_sub |>
    select(all_of(env_vars)) |>
    select(where(is.numeric)) |>
    select(where(~ sd(., na.rm = TRUE) > 0))

  rownames(env_mat) <- env_sub$Plot
  env_mat
}

run_soil_mantel_qcorrplot <- function(env,
                                      trial,
                                      microbiome_mats,
                                      association,
                                      env_vars = NULL) {
  common_plots <- Reduce(
    intersect,
    c(lapply(microbiome_mats, rownames), list(env |> filter(Trial == trial) |> pull(Plot)))
  )

  microbe_mats2 <- lapply(microbiome_mats, function(x) {
    x[common_plots, , drop = FALSE]
  })

  env_var <- prepare_env_matrix(env, trial, common_plots, env_vars)

  microbe_mats2 <- lapply(microbe_mats2, function(x) {
    x[rownames(env_var), , drop = FALSE]
  })

  keep <- complete.cases(env_var)
  env_var <- env_var[keep, , drop = FALSE]
  microbe_mats2 <- lapply(microbe_mats2, function(x) {
    x[keep, , drop = FALSE]
  })

  stopifnot(all(vapply(microbe_mats2, function(x) {
    all(rownames(x) == rownames(env_var))
  }, logical(1))))

  spec_all <- bind_cols(lapply(names(microbe_mats2), function(nm) {
    x <- as.data.frame(microbe_mats2[[nm]])
    colnames(x) <- paste0(nm, "_", colnames(x))
    x
  }))

  group_sizes <- vapply(microbe_mats2, ncol, integer(1))
  group_ends <- cumsum(group_sizes)
  group_starts <- group_ends - group_sizes + 1
  spec_select <- map2(group_starts, group_ends, ~ .x:.y)
  names(spec_select) <- names(microbe_mats2)

  mantel_res <- mantel_test(
    spec = spec_all,
    env = env_var,
    spec_select = spec_select,
    spec_dist = "bray",
    env_dist = "euclidean",
    mantel_fun = "mantel"
  ) |>
    mutate(
      trial = trial,
      association = association,
      rd = cut(r, breaks = c(-Inf, 0.1, 0.3, Inf),
               labels = c("< 0.1", "0.1 - 0.3", ">= 0.3")),
      pd = cut(p, breaks = c(-Inf, 0.01, 0.05, Inf),
               labels = c("< 0.01", "0.01 - 0.05", ">= 0.05"))
    )

  mantel_summary <- mantel_res |>
    transmute(
      data_type = spec,
      trial = trial,
      association = association,
      soil_property = env,
      distance = "Bray-Curtis vs Euclidean",
      mantel_r = r,
      p_value = p
    )

  cor_env <- correlate(env_var, method = "spearman")

  p <- qcorrplot(cor_env, type = "lower", diag = FALSE) +
    geom_square() +
    geom_couple(
      aes(colour = pd, size = rd),
      data = mantel_res,
      curvature = nice_curvature()
    ) +
    scale_fill_gradient2(low = "#3B9AB2", mid = "white", high = "#F21A00") +
    scale_size_manual(values = c(0.5, 1, 2)) +
    scale_colour_manual(values = color_pal(3)) +
    guides(
      size = guide_legend(title = "Mantel's r"),
      colour = guide_legend(title = "Mantel's p"),
      fill = guide_colorbar(title = "Spearman's r")
    ) +
    ggtitle(paste(trial, association))

  list(
    env_var = env_var,
    microbiome_mats = microbe_mats2,
    mantel_res = mantel_res,
    mantel_summary = mantel_summary,
    cor_env = cor_env,
    plot = p
  )
}


root_bulk_soil_spring <- run_soil_mantel_qcorrplot(
  env = env,
  trial = "SpringPan",
  microbiome_mats = list(
    Bulk_Bacteria = bulk_bac_spring$bulk_hellinger,
    Bulk_Fungi = bulk_fun_spring$bulk_hellinger,
    Root_Bacteria = root_bulk_bac_spring$root_hellinger,
    Root_Fungi = root_bulk_fun_spring$root_hellinger
  ),
  association = "root and bulk microbiomes - soil properties"
)

root_bulk_soil_winter <- run_soil_mantel_qcorrplot(
  env = env,
  trial = "WinterPan",
  microbiome_mats = list(
    Bulk_Bacteria = bulk_bac_winter$bulk_hellinger,
    Bulk_Fungi = bulk_fun_winter$bulk_hellinger,
    Root_Bacteria = root_bulk_bac_winter$root_hellinger,
    Root_Fungi = root_bulk_fun_winter$root_hellinger
  ),
  association = "root and bulk microbiomes - soil properties"
)

root_bulk_soil_spring$plot
root_bulk_soil_spring$mantel_summary

root_bulk_soil_winter$plot
root_bulk_soil_winter$mantel_summary


soil_summary <- bind_rows(
  root_bulk_soil_spring$mantel_summary,
  root_bulk_soil_winter$mantel_summary
)

integrated_mantel_summary <- bind_rows(
  root_bulk_summary,
  soil_summary
) |>
  arrange(data_type, trial, association, distance, p_value)

integrated_mantel_summary


write.csv(
  bulk_pcnm_summary,
  "bulk_pcnm_model_summary.csv",
  row.names = FALSE
)

write.csv(
  integrated_mantel_summary,
  "root_bulk_soil_integrated_mantel_summary.csv",
  row.names = FALSE
)
