# Figures: ExtendedDataFig1
# Original source: 04_SNP_PCNM_interactions/plot_SNP_PCNM_nonsig_examples.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

library(tidyverse)
library(cowplot)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}

source(file.path(script_dir, "SNP_biomarker_PCNM_interaction_models_extracted.R"), echo = FALSE)

out_dir <- file.path(script_dir, "SNP_biomarker_PCNM_interaction_results")
allele_levels <- c("A", "C", "G", "T")
allele_colors <- c(
  "A" = "#F8766D",
  "C" = "#7CAE00",
  "G" = "#00BFC4",
  "T" = "#C77CFF"
)

nonsig_terms <- term_effects %>%
  filter(
    term_type == "snp_pcnm_interaction",
    p.value >= 0.05,
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

nonsig_slopes <- bind_rows(lapply(split(nonsig_terms, seq_len(nrow(nonsig_terms))), get_simple_slopes)) %>%
  mutate(
    slope_gap = abs(slope_alt - slope_ref),
    pcnm_num = readr::parse_number(PCNM)
  ) %>%
  arrange(pcnm_num, interaction_p, desc(slope_gap), Organism, Trial, Biomarker) %>%
  group_by(PCNM) %>%
  mutate(panel_rank = row_number()) %>%
  ungroup() %>%
  mutate(panel_id = paste(PCNM, panel_rank, sep = "__"))

plot_dat <- if (nrow(nonsig_slopes) == 0) {
  tibble()
} else {
  model_fits %>%
    semi_join(
      nonsig_slopes %>% select(Organism, Trial, Biomarker),
      by = c("Organism", "Trial", "Biomarker")
    ) %>%
    select(Organism, Trial, Biomarker, model_data) %>%
    mutate(
      model_data = purrr::map(model_data, ~ .x %>% select(-any_of(c("Organism", "Trial", "Biomarker"))))
    ) %>%
    unnest(model_data) %>%
    inner_join(
      nonsig_slopes,
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
        Organism, " | ", Biomarker, " | ", SNP, "\n",
        Trial, " | P=", signif(interaction_p, 2), "\n",
        "slope(", ref_allele, ")=", signif(slope_ref, 2),
        ", slope(", alt_allele, ")=", signif(slope_alt, 2)
      ),
      biomarker_label = forcats::fct_inorder(biomarker_label)
    ) %>%
    filter(!is.na(allele), is.finite(pcnm_value))
}

make_panel_plot <- function(df_panel) {
  panel_label <- unique(df_panel$biomarker_label)
  ggplot(df_panel, aes(x = pcnm_value, y = abundance_rclr, colour = allele)) +
    geom_point(alpha = 0.28, size = 0.9) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.78) +
    scale_colour_manual(values = allele_colors, breaks = allele_levels, drop = FALSE) +
    theme_bw(base_size = 8.5) +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = 6.4, hjust = 0.5, lineheight = 0.94),
      axis.title = element_blank(),
      axis.text = element_text(size = 6.5),
      plot.margin = margin(3, 3, 3, 3)
    ) +
    labs(title = panel_label, colour = "Allele")
}

make_header_plot <- function(label, n_panels) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.6, label = label, size = 4.1, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.23, label = paste0("n = ", n_panels), size = 3.1) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      plot.background = element_rect(fill = "grey92", colour = "grey55", linewidth = 0.6),
      plot.margin = margin(2, 2, 2, 2)
    )
}

p <- if (nrow(plot_dat) == 0) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No non-significant SNP-by-PCNM interactions found") +
    theme_void()
} else {
  pcnm_order_tbl <- nonsig_slopes %>%
    count(PCNM, pcnm_num, name = "n_panels") %>%
    arrange(pcnm_num, PCNM)

  pcnm_order <- pcnm_order_tbl$PCNM

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
  max_rows <- max(pcnm_order_tbl$n_panels)

  left_margin <- 0.038
  right_margin <- 0.012
  bottom_margin <- 0.04
  title_height <- 0.055
  subtitle_height <- 0.03
  legend_height <- 0.04
  header_height <- 0.03
  x_gap <- 0.0045
  y_gap <- 0.0045

  panel_width <- (1 - left_margin - right_margin - (n_cols - 1) * x_gap) / n_cols
  panel_height <- (1 - bottom_margin - title_height - subtitle_height - legend_height - header_height - (max_rows - 1) * y_gap) / max_rows

  composed <- cowplot::ggdraw() +
    cowplot::draw_label(
      "Non-significant SNP-by-PCNM interaction regressions grouped by PCNM",
      x = left_margin,
      y = 0.995,
      hjust = 0,
      vjust = 1,
      size = 16
    ) +
    cowplot::draw_label(
      "All panels with nominal interaction P >= 0.05, arranged with the same PCNM-centered logic as the main figure",
      x = left_margin,
      y = 0.968,
      hjust = 0,
      vjust = 1,
      size = 10.5
    ) +
    cowplot::draw_plot(shared_legend, x = 0.43, y = 0.915, width = 0.16, height = legend_height)

  for (col_idx in seq_along(pcnm_order)) {
    pcnm_i <- pcnm_order[col_idx]
    x_left <- left_margin + (col_idx - 1) * (panel_width + x_gap)
    header_y <- bottom_margin + max_rows * panel_height + (max_rows - 1) * y_gap
    n_panels_i <- pcnm_order_tbl$n_panels[pcnm_order_tbl$PCNM == pcnm_i]

    composed <- composed + cowplot::draw_plot(
      make_header_plot(pcnm_i, n_panels_i),
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
    cowplot::draw_label("PCNM value", x = 0.5, y = 0.014, size = 11) +
    cowplot::draw_label("rCLR biomarker abundance", x = 0.012, y = 0.42, angle = 90, size = 11)
}

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_nonsig_by_pcnm_columns_all.png"),
  p,
  width = 42,
  height = 32,
  units = "in",
  dpi = 300,
  limitsize = FALSE,
  bg = "white"
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_nonsig_by_pcnm_columns_all.pdf"),
  p,
  width = 42,
  height = 32,
  units = "in",
  limitsize = FALSE,
  bg = "white"
)

write.csv(
  nonsig_slopes %>%
    select(Organism, Trial, Biomarker, SNP, PCNM, interaction_p, ref_allele, alt_allele, slope_ref, slope_alt, slope_gap),
  file.path(out_dir, "SNP_biomarker_PCNM_nonsig_by_pcnm_columns_all.csv"),
  row.names = FALSE
)
