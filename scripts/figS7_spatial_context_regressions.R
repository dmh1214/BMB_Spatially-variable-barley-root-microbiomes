# Figures: FigS7
# Original source: 08_compare_xy_pcnm/supp_context_grouped_regressions_slope_top10_q33.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggsci)
  library(broom)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(patchwork)
})

base_theme <- theme_classic(base_size = 10.5) +
  theme(
    text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8.4, color = "black"),
    plot.title = element_text(size = 10.5, face = "plain", hjust = 0.5),
    plot.subtitle = element_text(size = 8.8, color = "grey30", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    legend.title = element_text(size = 9.5),
    legend.text = element_text(size = 8.5)
  )

root_dir <- "/Users/gaoqi/Library/CloudStorage/OneDrive-Aarhusuniversitet/Desktop/BMB_AUT_MS"
out_dir <- file.path(root_dir, "08_compare_xy_pcnm")

pcnm_env <- new.env(parent = emptyenv())
bac_env <- new.env(parent = emptyenv())
fun_env <- new.env(parent = emptyenv())

load(file.path(root_dir, "PCNM/PCNM.Rdata"), envir = pcnm_env)
load(file.path(root_dir, "Bac_Aut/snp.Rdata"), envir = bac_env)
load(file.path(root_dir, "Fun_Aut/snp.Rdata"), envir = fun_env)

extract_top10_pcnm <- function(model_obj) head(attr(terms(model_obj), "term.labels"), 10)

pair_mean_upper <- function(x) {
  m <- outer(x, x, "+") / 2
  m[upper.tri(m)]
}

build_panel_meta <- function(meta_df, sample_ids) {
  meta_panel <- meta_df %>%
    distinct(sampleid, .keep_all = TRUE) %>%
    filter(sampleid %in% sample_ids)
  meta_panel[match(sample_ids, meta_panel$sampleid), , drop = FALSE]
}

assign_tertiles <- function(x) {
  q1 <- as.numeric(quantile(x, 1/3, na.rm = TRUE))
  q2 <- as.numeric(quantile(x, 2/3, na.rm = TRUE))
  case_when(
    x < q1 ~ "Near",
    x > q2 ~ "Far",
    TRUE ~ "Mid"
  )
}

make_metric_pairs <- function(spec, metric_name) {
  meta_panel <- build_panel_meta(spec$meta, spec$sample_ids)

  if (metric_name == "XY distance") {
    spatial_metric <- as.numeric(stats::dist(as.matrix(meta_panel[, c("Row", "Col")]), method = "euclidean"))
  } else {
    meta_panel <- meta_panel %>% mutate(Plot = as.character(Plot))
    axis_scores <- pcnm_env$pcnm30 %>%
      mutate(Plot = as.character(Plot)) %>%
      dplyr::select(Plot, !!metric_name) %>%
      rename(axis_value = !!metric_name)
    spatial_metric <- meta_panel %>%
      left_join(axis_scores, by = "Plot") %>%
      pull(axis_value) %>%
      pair_mean_upper()
  }

  tibble(
    Prevosti = as.numeric(spec$d_prev),
    BrayCurtis = as.numeric(spec$d_comm),
    spatial_metric = spatial_metric
  ) %>%
    filter(Prevosti > 0, is.finite(Prevosti), is.finite(BrayCurtis), is.finite(spatial_metric)) %>%
    mutate(
      Prevosti_z = as.numeric(scale(Prevosti)),
      BrayCurtis_z = as.numeric(scale(BrayCurtis)),
      spatial_group = factor(assign_tertiles(spatial_metric), levels = c("Near", "Mid", "Far")),
      metric = metric_name,
      Panel = spec$Panel,
      Organism = spec$Organism
    )
}

fit_lines <- function(df) {
  group_models <- df %>%
    group_by(spatial_group) %>%
    group_map(~ lm(BrayCurtis_z ~ Prevosti_z, data = .x), .keep = TRUE)
  names(group_models) <- levels(df$spatial_group)
  global_model <- lm(BrayCurtis_z ~ Prevosti_z, data = df)
  xseq <- seq(min(df$Prevosti_z), max(df$Prevosti_z), length.out = 100)

  line_df <- bind_rows(lapply(names(group_models), function(g) {
    tibble(
      Prevosti_z = xseq,
      BrayCurtis_z = predict(group_models[[g]], newdata = tibble(Prevosti_z = xseq)),
      line_group = g
    )
  }))

  global_df <- tibble(
    Prevosti_z = xseq,
    BrayCurtis_z = predict(global_model, newdata = tibble(Prevosti_z = xseq)),
    line_group = "Overall"
  )

  list(lines = line_df, global = global_df)
}

summarise_slopes <- function(df, metric_name, panel_name, organism_name) {
  global_fit <- lm(BrayCurtis_z ~ Prevosti_z, data = df)
  df %>%
    group_by(spatial_group) %>%
    do(tidy(lm(BrayCurtis_z ~ Prevosti_z, data = .), conf.int = TRUE)) %>%
    filter(term == "Prevosti_z") %>%
    ungroup() %>%
    transmute(
      Panel = panel_name,
      Organism = organism_name,
      metric = metric_name,
      spatial_group,
      slope = estimate,
      conf.low,
      conf.high,
      global_slope = unname(coef(global_fit)[["Prevosti_z"]])
    )
}

panel_specs <- list(
  list(
    Panel = "SpringPan",
    Organism = "Bacteria",
    meta = bac_env$meta %>% filter(Trial == "SpringPan"),
    sample_ids = bac_env$sampleid_spring,
    d_prev = bac_env$D_prev_spring_sample,
    d_comm = bac_env$D_asv_spring,
    top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_spring_bac)
  ),
  list(
    Panel = "WinterPan",
    Organism = "Bacteria",
    meta = bac_env$meta %>% filter(Trial == "WinterPan"),
    sample_ids = bac_env$sampleid_winter,
    d_prev = bac_env$D_prev_winter_sample,
    d_comm = bac_env$D_asv_winter,
    top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_winter_bac)
  ),
  list(
    Panel = "SpringPan",
    Organism = "Fungi",
    meta = fun_env$meta %>% filter(Trial == "SpringPan"),
    sample_ids = fun_env$sampleid_spring,
    d_prev = fun_env$D_prev_spring_sample,
    d_comm = fun_env$D_otu_spring,
    top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_spring_fun)
  ),
  list(
    Panel = "WinterPan",
    Organism = "Fungi",
    meta = fun_env$meta %>% filter(Trial == "WinterPan"),
    sample_ids = fun_env$sampleid_winter,
    d_prev = fun_env$D_prev_winter_sample,
    d_comm = fun_env$D_otu_winter,
    top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_winter_fun)
  )
)

all_pairs <- map(panel_specs, function(spec) {
  metric_order <- c("XY distance", spec$top10)
  metric_pairs <- map(metric_order, ~ make_metric_pairs(spec, .x))
  names(metric_pairs) <- metric_order
  metric_pairs
})
names(all_pairs) <- map_chr(panel_specs, ~ paste(.x$Organism, .x$Panel, sep = "_"))

slope_tbl <- map2_dfr(panel_specs, names(all_pairs), function(spec, key) {
  metric_order <- c("XY distance", spec$top10)
  bind_rows(lapply(metric_order, function(metric_name) {
    summarise_slopes(
      all_pairs[[key]][[metric_name]],
      metric_name = metric_name,
      panel_name = spec$Panel,
      organism_name = spec$Organism
    )
  }))
}) %>%
  mutate(spatial_group = factor(spatial_group, levels = c("Near", "Mid", "Far")))

slope_xlim <- range(c(slope_tbl$conf.low, slope_tbl$conf.high), na.rm = TRUE)
line_cols <- c(Near = "#C98A8A", Mid = "#8F8F8F", Far = "#7FA6C9")

make_line_plot <- function(df, metric_name, show_x = FALSE, show_y = FALSE) {
  fits <- fit_lines(df)
  y_rng <- range(c(fits$lines$BrayCurtis_z, fits$global$BrayCurtis_z), na.rm = TRUE)
  y_pad <- diff(y_rng) * 0.14
  if (!is.finite(y_pad) || y_pad == 0) y_pad <- 0.12
  label_x <- -2.25
  label_y <- y_rng[2] + y_pad * 0.72

  ggplot() +
    geom_line(
      data = fits$global,
      aes(x = Prevosti_z, y = BrayCurtis_z),
      color = "grey40",
      linetype = "dashed",
      linewidth = 0.55,
      alpha = 0.68
    ) +
    geom_line(
      data = fits$lines,
      aes(x = Prevosti_z, y = BrayCurtis_z, color = line_group),
      linewidth = 0.95,
      alpha = 0.78
    ) +
    scale_color_manual(values = line_cols) +
    coord_cartesian(xlim = c(-2.4, 2.4), ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)) +
    annotate("text", x = label_x, y = label_y, label = metric_name, hjust = 0, vjust = 1, size = 3.45) +
    base_theme +
    labs(
      x = if (show_x) "Host genetic distance (z)" else NULL,
      y = if (show_y) "Community dissimilarity (z)" else NULL,
      color = "Context group"
    ) +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 9.4),
      axis.title.y = element_text(size = 10.8),
      axis.text.x = if (show_x) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.text.y = if (show_y) element_text(size = 8) else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      aspect.ratio = 0.53,
      plot.margin = margin(0.8, 1.1, 0.8, 1.1)
    )
}

make_slope_plot <- function(spec, show_x = FALSE) {
  panel_df <- slope_tbl %>%
    filter(Panel == spec$Panel, Organism == spec$Organism)
  y_levels <- c("XY distance", spec$top10)

  ggplot(panel_df, aes(x = slope, y = metric, color = spatial_group)) +
    geom_vline(aes(xintercept = global_slope), linetype = "dashed", color = "grey45", linewidth = 0.35, show.legend = FALSE) +
    geom_hline(
      yintercept = seq(1.5, length(y_levels) - 0.5, by = 1),
      linetype = "dashed",
      color = "grey88",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    geom_errorbar(
      aes(xmin = conf.low, xmax = conf.high),
      width = 0.16,
      linewidth = 0.45,
      position = position_dodge(width = 0.56),
      orientation = "y"
    ) +
    geom_point(size = 1.7, position = position_dodge(width = 0.56)) +
    scale_color_manual(values = line_cols) +
    scale_y_discrete(limits = rev(y_levels)) +
    coord_cartesian(xlim = slope_xlim) +
    base_theme +
    labs(
      title = "Slope",
      x = if (show_x) "Slope" else NULL,
      y = NULL,
      color = "Context group"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 9.6, hjust = 0.5),
      axis.text.y = element_text(size = 8.8),
      axis.title.x = element_text(size = 9.4),
      axis.text.x = if (show_x) element_text(size = 8) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.margin = margin(3, 2, 1, 1)
    )
}

row_blocks <- map(seq_along(panel_specs), function(i) {
  spec <- panel_specs[[i]]
  key <- paste(spec$Organism, spec$Panel, sep = "_")
  metric_order <- c("XY distance", spec$top10)

  left_column <- wrap_plots(
    map(seq_along(metric_order), function(j) {
      make_line_plot(
        all_pairs[[key]][[metric_order[[j]]]],
        metric_name = metric_order[[j]],
        show_x = i == length(panel_specs) && j == length(metric_order),
        show_y = j == ceiling(length(metric_order) / 2)
      )
    }),
    ncol = 1
  ) +
    plot_annotation(
      title = str_c(spec$Organism, "\n", spec$Panel),
      subtitle = "XY + top 10 PCNMs"
    ) &
    theme(
      plot.title = element_text(size = 11.5, hjust = 0.5),
      plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey35")
    )

  (left_column | make_slope_plot(spec, show_x = i == length(panel_specs))) +
    plot_layout(widths = c(1.25, 0.8))
})

combined_plot <- wrap_plots(row_blocks, ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(12, "pt"),
    legend.key.height = unit(8, "pt"),
    panel.spacing = unit(4.5, "pt"),
    plot.margin = margin(2, 2, 2, 2)
  )

ggsave(
  file.path(out_dir, "supp_context_grouped_regressions_slope_top10_q33_q67.pdf"),
  combined_plot,
  width = 14.6,
  height = 20.4,
  useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "supp_context_grouped_regressions_slope_top10_q33_q67.png"),
  combined_plot,
  width = 14.6,
  height = 20.4,
  dpi = 300
)

message("Finished supplementary grouped regressions + slope top10 figure in: ", out_dir)
