# Prepared for reviewer use; original analysis retained in originals/.
if (!exists("project_root", inherits=FALSE)) {
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), winslash="/", mustWork=TRUE)
}
module_output <- file.path(project_root, "outputs", "fig5c")
dir.create(module_output, recursive=TRUE, showWarnings=FALSE)
setwd(module_output)
dir.create("Network/k_tables", recursive=TRUE, showWarnings=FALSE)
on_error <- NULL

library(tidyverse)
library(cowplot)


source(file.path(script_dir, "fig5c_snp_pcnm_models.R"), echo = FALSE)

out_dir <- file.path(project_root, "outputs", "fig5c")
allele_levels <- c("A", "C", "G", "T")
allele_colors <- c(
  "A" = "#F8766D",
  "C" = "#7CAE00",
  "G" = "#00BFC4",
  "T" = "#C77CFF"
)

sig_terms <- term_effects %>%
  filter(
    term_type == "snp_pcnm_interaction",
    p.value < 0.05,
    !is.na(pcnm_variable),
    !is.na(snp_variable)
  )

get_simple_slopes <- function(row) {
  org <- row[["Organism"]]
  tr <- row[["Trial"]]
  bio <- row[["Biomarker"]]
  snp_var <- row[["snp_variable"]]
  pcnm_var <- row[["pcnm_variable"]]

  fit_row <- model_fits %>%
    filter(Organism == org, Trial == tr, Biomarker == bio)

  if (nrow(fit_row) != 1) return(tibble())

  model <- fit_row$model[[1]]
  df <- fit_row$model_data[[1]]
  if (is.null(model) || !snp_var %in% names(df) || !pcnm_var %in% names(df)) return(tibble())

  coefs <- coef(model)
  levs <- levels(df[[snp_var]])
  if (length(levs) != 2 || !pcnm_var %in% names(coefs)) return(tibble())

  int_name_1 <- paste0(snp_var, levs[2], ":", pcnm_var)
  int_name_2 <- paste0(pcnm_var, ":", snp_var, levs[2])
  int_coef <- if (int_name_1 %in% names(coefs)) {
    coefs[int_name_1]
  } else if (int_name_2 %in% names(coefs)) {
    coefs[int_name_2]
  } else {
    NA_real_
  }

  if (is.na(int_coef)) return(tibble())

  tibble(
    Organism = org,
    Trial = tr,
    Biomarker = bio,
    SNP = row[["snp_name"]],
    snp_variable = snp_var,
    PCNM = pcnm_var,
    interaction_p = row[["p.value"]],
    ref_allele = levs[1],
    alt_allele = levs[2],
    slope_ref = unname(coefs[pcnm_var]),
    slope_alt = unname(coefs[pcnm_var] + int_coef)
  )
}

slopes <- bind_rows(lapply(split(sig_terms, seq_len(nrow(sig_terms))), get_simple_slopes)) %>%
  mutate(
    pattern = if_else(slope_ref * slope_alt < 0, "crossover", "non_crossover")
  )

# Keep the strongest crossover examples per PCNM block
top_crossovers <- slopes %>%
  filter(pattern == "crossover") %>%
  arrange(interaction_p) %>%
  group_by(PCNM) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  group_by(PCNM) %>%
  mutate(panel_rank = row_number()) %>%
  ungroup() %>%
  mutate(panel_id = paste(PCNM, panel_rank, sep = "__"))

plot_dat <- if (nrow(top_crossovers) == 0) {
  tibble()
} else {
  model_fits %>%
    semi_join(
      top_crossovers %>% select(Organism, Trial, Biomarker),
      by = c("Organism", "Trial", "Biomarker")
    ) %>%
    select(Organism, Trial, Biomarker, model_data) %>%
    mutate(
      model_data = purrr::map(model_data, ~ .x %>% select(-any_of(c("Organism", "Trial", "Biomarker"))))
    ) %>%
    unnest(model_data) %>%
    inner_join(
      top_crossovers,
      by = c("Organism", "Trial", "Biomarker"),
      relationship = "many-to-many"
    ) %>%
    pivot_longer(
      cols = matches("^PCNM\\d+$"),
      names_to = "pcnm_current",
      values_to = "pcnm_value"
    ) %>%
    filter(pcnm_current == PCNM) %>%
    mutate(
      allele = case_when(
        snp_variable == "SNP1" ~ as.character(SNP1),
        snp_variable == "SNP2" ~ as.character(SNP2),
        TRUE ~ NA_character_
      ),
      allele = factor(allele, levels = allele_levels),
      biomarker_label = paste0(
        SNP, "\n",
        Organism, " | ", Trial, " | ", Biomarker
      ),
      interaction_label = paste0("interaction P = ", signif(interaction_p, 2)),
      slope_label = paste0(
        "slope(", ref_allele, ")=", signif(slope_ref, 2),
        ", slope(", alt_allele, ")=", signif(slope_alt, 2)
      ),
      biomarker_label = forcats::fct_inorder(biomarker_label)
    ) %>%
    filter(!is.na(allele), is.finite(pcnm_value))
}

make_panel_plot <- function(df_panel) {
  panel_label <- unique(df_panel$biomarker_label)
  p_label <- unique(df_panel$interaction_label)
  slope_label <- unique(df_panel$slope_label)
  label_df <- tibble(
    x = mean(range(df_panel$pcnm_value, na.rm = TRUE)),
    y = max(df_panel$abundance_rclr, na.rm = TRUE) * 0.92 +
      min(df_panel$abundance_rclr, na.rm = TRUE) * 0.08,
    label = p_label
  )
  slope_df <- tibble(
    x = mean(range(df_panel$pcnm_value, na.rm = TRUE)),
    y = min(df_panel$abundance_rclr, na.rm = TRUE) * 0.92 +
      max(df_panel$abundance_rclr, na.rm = TRUE) * 0.08,
    label = slope_label
  )
  ggplot(df_panel, aes(x = pcnm_value, y = abundance_rclr, colour = allele)) +
    geom_point(alpha = 0.32, size = 0.95) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      vjust = 0.2,
      size = 2.7,
      colour = "black"
    ) +
    geom_text(
      data = slope_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      vjust = 1.2,
      size = 2.45,
      colour = "black"
    ) +
    scale_colour_manual(values = allele_colors, breaks = allele_levels, drop = FALSE) +
    theme_bw(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 7, hjust = 0.5, lineheight = 0.95),
      axis.title = element_blank(),
      plot.margin = margin(4, 4, 4, 4)
    ) +
    labs(title = panel_label, colour = "Allele")
}

make_header_plot <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 4) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "grey90", colour = "grey50", linewidth = 0.6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

p <- if (nrow(plot_dat) == 0) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No crossover SNP-by-PCNM interactions found") +
    theme_void()
} else {
  pcnm_order <- top_crossovers %>%
    count(PCNM, sort = TRUE, name = "n_panels") %>%
    mutate(PCNM = factor(PCNM, levels = PCNM)) %>%
    pull(PCNM) %>%
    as.character()

  panel_meta <- plot_dat %>%
    distinct(panel_id, PCNM, biomarker_label) %>%
    mutate(PCNM = factor(PCNM, levels = pcnm_order)) %>%
    arrange(PCNM, panel_id)

  legend_plot <- ggplot(plot_dat, aes(x = pcnm_value, y = abundance_rclr, colour = allele)) +
    geom_point() +
    scale_colour_manual(values = allele_colors, breaks = allele_levels, drop = FALSE) +
    theme_void() +
    theme(legend.position = "top") +
    labs(colour = "Allele")

  shared_legend <- cowplot::get_legend(legend_plot)

  n_cols <- length(pcnm_order)
  max_rows <- panel_meta %>%
    count(PCNM) %>%
    summarise(max_n = max(n)) %>%
    pull(max_n)

  left_margin <- 0.06
  right_margin <- 0.02
  bottom_margin <- 0.06
  title_height <- 0.075
  legend_height <- 0.05
  header_height <- 0.04
  x_gap <- 0.008
  y_gap <- 0.008

  panel_width <- (1 - left_margin - right_margin - (n_cols - 1) * x_gap) / n_cols
  panel_height <- (1 - bottom_margin - title_height - legend_height - header_height - (max_rows - 1) * y_gap) / max_rows

  composed <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Crossover SNP-by-PCNM interactions grouped by PCNM context",
      x = left_margin,
      y = 0.99,
      hjust = 0,
      vjust = 1,
      fontface = "plain",
      size = 16
    ) +
    cowplot::draw_label(
      "Top crossover examples within each PCNM across bacteria and fungi",
      x = left_margin,
      y = 0.965,
      hjust = 0,
      vjust = 1,
      size = 11
    ) +
    cowplot::draw_plot(shared_legend, x = 0.32, y = 0.89, width = 0.36, height = legend_height)

  for (col_idx in seq_along(pcnm_order)) {
    pcnm_i <- pcnm_order[col_idx]
    x_left <- left_margin + (col_idx - 1) * (panel_width + x_gap)
    header_y <- bottom_margin + max_rows * panel_height + (max_rows - 1) * y_gap

    composed <- composed + cowplot::draw_plot(
      make_header_plot(pcnm_i),
      x = x_left,
      y = header_y,
      width = panel_width,
      height = header_height
    )

    panel_ids_i <- panel_meta %>%
      filter(PCNM == pcnm_i) %>%
      pull(panel_id)

    for (row_idx in seq_along(panel_ids_i)) {
      panel_id_i <- panel_ids_i[row_idx]
      y_bottom <- header_y - row_idx * panel_height - (row_idx - 1) * y_gap
      panel_df <- plot_dat %>% filter(panel_id == panel_id_i)

      composed <- composed + cowplot::draw_plot(
        make_panel_plot(panel_df),
        x = x_left,
        y = y_bottom,
        width = panel_width,
        height = panel_height
      )
    }
  }

  composed +
    cowplot::draw_label("PCNM value", x = 0.5, y = 0.02, size = 12) +
    cowplot::draw_label("rCLR biomarker abundance", x = 0.015, y = 0.42, angle = 90, size = 12)
}

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_by_organism_pcnm.png"),
  p,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_by_organism_pcnm.pdf"),
  p,
  width = 14,
  height = 11,
  bg = "white"
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_staircase_white.png"),
  p,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_staircase_white.pdf"),
  p,
  width = 14,
  height = 11,
  bg = "white"
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_staircase.png"),
  p,
  width = 14,
  height = 11,
  dpi = 300,
  bg = "white"
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_staircase.pdf"),
  p,
  width = 14,
  height = 11,
  bg = "white"
)

write.csv(
  top_crossovers,
  file.path(out_dir, "SNP_biomarker_PCNM_crossover_by_organism_pcnm.csv"),
  row.names = FALSE
)
