rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggthemes)
  library(ggsci)
  library(broom)
  library(purrr)
  library(patchwork)
})

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
      spatial_z = as.numeric(scale(spatial_metric))
    )
}

fit_interaction <- function(spec, metric_name) {
  make_metric_pairs(spec, metric_name) %>%
    { tidy(lm(BrayCurtis_z ~ Prevosti_z * spatial_z, data = .), conf.int = TRUE) } %>%
    filter(term == "Prevosti_z:spatial_z") %>%
    transmute(
      Panel = spec$Panel,
      Organism = spec$Organism,
      metric = metric_name,
      estimate,
      conf.low,
      conf.high,
      p.value
    )
}

panel_specs <- list(
  list(Panel = "SpringPan", Organism = "Bacteria", meta = bac_env$meta %>% filter(Trial == "SpringPan"),
       sample_ids = bac_env$sampleid_spring, d_prev = bac_env$D_prev_spring_sample, d_comm = bac_env$D_asv_spring,
       top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_spring_bac)),
  list(Panel = "WinterPan", Organism = "Bacteria", meta = bac_env$meta %>% filter(Trial == "WinterPan"),
       sample_ids = bac_env$sampleid_winter, d_prev = bac_env$D_prev_winter_sample, d_comm = bac_env$D_asv_winter,
       top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_winter_bac)),
  list(Panel = "SpringPan", Organism = "Fungi", meta = fun_env$meta %>% filter(Trial == "SpringPan"),
       sample_ids = fun_env$sampleid_spring, d_prev = fun_env$D_prev_spring_sample, d_comm = fun_env$D_otu_spring,
       top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_spring_fun)),
  list(Panel = "WinterPan", Organism = "Fungi", meta = fun_env$meta %>% filter(Trial == "WinterPan"),
       sample_ids = fun_env$sampleid_winter, d_prev = fun_env$D_prev_winter_sample, d_comm = fun_env$D_otu_winter,
       top10 = extract_top10_pcnm(pcnm_env$pcnm_sel_winter_fun))
)

interaction_tbl <- bind_rows(lapply(panel_specs, function(spec) {
  bind_rows(lapply(c("XY distance", spec$top10), function(metric_name) {
    fit_interaction(spec, metric_name)
  }))
})) %>%
  group_by(Panel, Organism) %>%
  mutate(
    p_adj = p.adjust(p.value, method = "BH"),
    direction = case_when(
      p_adj < 0.05 & estimate > 0 ~ "Positive",
      p_adj < 0.05 & estimate < 0 ~ "Negative",
      TRUE ~ "NS"
    ),
    direction = factor(direction, levels = c("Positive", "NS", "Negative"))
  ) %>%
  ungroup()

right_xlim <- range(c(interaction_tbl$conf.low, interaction_tbl$conf.high), na.rm = TRUE)
interaction_cols <- c(
  Positive = "#6B9EC9",
  NS = "grey70",
  Negative = "#C98A8A"
)

base_theme <- theme_classic(base_size = 9) +
  theme(
    text = element_text(color = "black"),
    axis.line = element_line(linewidth = 0.35, color = "black"),
    axis.ticks = element_line(linewidth = 0.35, color = "black"),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.2, color = "black"),
    plot.title = element_text(size = 10, hjust = 0.5),
    plot.subtitle = element_text(size = 8, hjust = 0.5, color = "grey35"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
    legend.title = element_text(size = 8.5),
    legend.text = element_text(size = 7.5)
  )

make_panel_plot <- function(spec, show_x = FALSE) {
  panel_df <- interaction_tbl %>% filter(Panel == spec$Panel, Organism == spec$Organism)
  y_levels <- c("XY distance", spec$top10)

  ggplot(panel_df, aes(x = estimate, y = metric, color = direction)) +
    annotate("rect", xmin = right_xlim[1], xmax = 0, ymin = -Inf, ymax = Inf, fill = "#F9EAEA", alpha = 0.28) +
    annotate("rect", xmin = 0, xmax = right_xlim[2], ymin = -Inf, ymax = Inf, fill = "#EAF3FA", alpha = 0.28) +
    geom_hline(
      yintercept = seq(1.5, length(y_levels) - 0.5, by = 1),
      linetype = "dashed", color = "grey90", linewidth = 0.25, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.36) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.18, linewidth = 0.58, orientation = "y") +
    geom_point(size = 2.1) +
    annotate("text", x = right_xlim[1] * 0.78, y = length(y_levels) + 0.45, label = "weaker across\ndissimilar contexts", size = 2.55, color = "#9E6F6F") +
    annotate("text", x = right_xlim[2] * 0.74, y = length(y_levels) + 0.45, label = "stronger across\ndissimilar contexts", size = 2.55, color = "#6B8FB2") +
    scale_color_manual(values = interaction_cols) +
    scale_y_discrete(limits = rev(y_levels)) +
    coord_cartesian(xlim = right_xlim, clip = "off") +
    base_theme +
    labs(
      title = paste(spec$Organism, spec$Panel, sep = "\n"),
      subtitle = "XY + top 10 PCNMs",
      x = if (show_x) "Interaction coefficient" else NULL,
      y = NULL,
      color = "Interaction"
    ) +
    theme(
      axis.text.x = if (show_x) element_text(size = 7.5) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    )
}

panel_plots <- lapply(seq_along(panel_specs), function(i) {
  make_panel_plot(panel_specs[[i]], show_x = i > 2)
})

combined_plot <- wrap_plots(panel_plots, ncol = 2, guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.width = unit(14, "pt"),
    legend.key.height = unit(9, "pt")
  )

ggsave(
  file.path(out_dir, "main_interaction_only_top10_q33_q67.pdf"),
  combined_plot,
  width = 11.2,
  height = 10.4,
  useDingbats = FALSE
)
ggsave(
  file.path(out_dir, "main_interaction_only_top10_q33_q67.png"),
  combined_plot,
  width = 11.2,
  height = 10.4,
  dpi = 300
)

message("Finished main interaction-only top10 figure in: ", out_dir)
