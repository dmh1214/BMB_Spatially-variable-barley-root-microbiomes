# Figures: Fig4e,f
# Original source: Fun_Aut/mGWAS/mGWAS_plot.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

## ============================================================
## Biomarker OTUs: SHAP + LME + GWAS summary (Spring & Winter)
## Final ordering:
##   - Keep original within-panel rank (variable_ranked levels)
##   - Final combined figure: SpringPan block first, WinterPan block second
## ============================================================

library(tidyverse)
library(readr)
library(viridis)
library(ggbeeswarm)
library(grid)
library(patchwork)
library(ggrastr)

rm(list = ls())

## ------------------------------------------------------------
## 0) User settings
## ------------------------------------------------------------
p_cutoff <- 5e-5
outdir_spring <- "Spring_gwas_res"
outdir_winter <- "winter_gwas_res"

## ------------------------------------------------------------
## 1) Helper: read GEMMA assoc file from one OTU folder
## ------------------------------------------------------------
read_one_assoc <- function(d) {
  OTU <- basename(d)
  f <- list.files(d, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  if (length(f) == 0) return(NULL)
  
  df <- read.table(f[1], header = TRUE, stringsAsFactors = FALSE)
  
  pos_col <- intersect(c("ps", "pos", "bp", "position"), colnames(df))[1]
  if (is.na(pos_col)) {
    stop("No SNP position column (ps/pos/bp/position) found in: ", f[1],
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "))
  }
  if (!all(c("chr", "p_wald") %in% colnames(df))) {
    stop("Missing 'chr' or 'p_wald' in ", f[1],
         "\nAvailable columns: ", paste(colnames(df), collapse = ", "))
  }
  
  df %>%
    transmute(
      OTU     = OTU,
      chr_raw = as.character(chr),
      chr     = str_replace_all(chr_raw, "[^0-9]", ""),
      chr     = as.integer(chr),
      ps      = as.numeric(.data[[pos_col]]),
      p_wald  = as.numeric(p_wald)
    ) %>%
    filter(!is.na(chr), chr != 0,
           !is.na(ps), ps > 0,
           !is.na(p_wald), p_wald > 0)
}

## ------------------------------------------------------------
## 2) Helper: read all assoc under a trial folder
## ------------------------------------------------------------
read_trial_assoc <- function(outdir) {
  OTU_dirs <- list.dirs(outdir, recursive = FALSE, full.names = TRUE) %>%
    keep(~ grepl("^OTU", basename(.x)))
  if (length(OTU_dirs) == 0) stop("No OTU_ folders found under ", outdir)
  
  all_assoc <- map_dfr(OTU_dirs, read_one_assoc)
  if (nrow(all_assoc) == 0) stop("No valid SNP rows read from assoc files under ", outdir)
  
  all_assoc
}

## ------------------------------------------------------------
## 3) Helper: OTU × chr summary for bubble plot
## ------------------------------------------------------------
summarise_otu_chr <- function(all_assoc, p_cutoff) {
  all_assoc %>%
    group_by(OTU, chr) %>%
    summarise(
      n_snps = n(),
      n_sig  = sum(p_wald < p_cutoff),
      min_p  = min(p_wald),
      .groups = "drop"
    )
}

## ------------------------------------------------------------
## 4) Load biomarkers (Spring + Winter)
## ------------------------------------------------------------
load("../feature_selection/spring_biomarker.Rdata")
positive_biomarkers <- shap_lme_f %>%
  filter(shap_cor > 0,
         slope > 0,
         direction_match == "Consistent") %>%
  pull(variable)
shap_lme <- shap_lme %>%
  filter(variable %in% positive_biomarkers)
shap_lme_f <- shap_lme_f %>%
  filter(variable %in% positive_biomarkers) %>%
  mutate(variable_ranked = droplevels(variable_ranked))
shap_long <- shap_long %>%
  filter(variable %in% positive_biomarkers)
shap_lme_spring   <- shap_lme
shap_lme_f_spring <- shap_lme_f
shap_long_spring  <- shap_long

load("../feature_selection/winter_biomarker.Rdata")
positive_biomarkers <- shap_lme_f %>%
  filter(shap_cor > 0,
         slope > 0,
         direction_match == "Consistent") %>%
  pull(variable)
shap_lme <- shap_lme %>%
  filter(variable %in% positive_biomarkers)
shap_lme_f <- shap_lme_f %>%
  filter(variable %in% positive_biomarkers) %>%
  mutate(variable_ranked = droplevels(variable_ranked))
shap_long <- shap_long %>%
  filter(variable %in% positive_biomarkers)
shap_lme_winter   <- shap_lme
shap_lme_f_winter <- shap_lme_f
shap_long_winter  <- shap_long

## Add Trial column
shap_lme_spring$Trial   <- "SpringPan"
shap_lme_winter$Trial   <- "WinterPan"
shap_lme_f_spring$Trial <- "SpringPan"
shap_lme_f_winter$Trial <- "WinterPan"
shap_long_spring$Trial  <- "SpringPan"
shap_long_winter$Trial  <- "WinterPan"

shap_lme_all   <- bind_rows(shap_lme_spring, shap_lme_winter)
shap_lme_f_all <- bind_rows(shap_lme_f_spring, shap_lme_f_winter)
shap_long_all  <- bind_rows(shap_long_spring, shap_long_winter)

## ------------------------------------------------------------
## 5) Read GWAS assoc (Spring + Winter)
## ------------------------------------------------------------
all_assoc_spring <- read_trial_assoc(outdir_spring) %>% mutate(Trial = "SpringPan")
all_assoc_winter <- read_trial_assoc(outdir_winter) %>% mutate(Trial = "WinterPan")
all_assoc_all <- bind_rows(all_assoc_spring, all_assoc_winter)

## ------------------------------------------------------------
## 6) OTU × chr summary and keep only biomarker OTUs per trial
## ------------------------------------------------------------
OTU_chr_spring <- summarise_otu_chr(all_assoc_spring, p_cutoff) %>%
  mutate(Trial = "SpringPan")
OTU_chr_winter <- summarise_otu_chr(all_assoc_winter, p_cutoff) %>%
  mutate(Trial = "WinterPan")

write_tsv(OTU_chr_spring, file.path(outdir_spring, "spring_OTU_chr_summary.tsv"))
write_tsv(OTU_chr_winter, file.path(outdir_winter, "winter_OTU_chr_summary.tsv"))

## Join with biomarker ranks (within each panel)
heat_spring <- OTU_chr_spring %>%
  rename(variable = OTU) %>%
  left_join(shap_lme_f_spring %>% select(variable, variable_ranked, Trial), by = c("variable","Trial")) %>%
  filter(!is.na(variable_ranked)) %>%
  mutate(
    chr      = factor(chr, levels = sort(unique(chr))),
    logP     = -log10(min_p),
    logP_cap = pmin(logP, 10)
  )

heat_winter <- OTU_chr_winter %>%
  rename(variable = OTU) %>%
  left_join(shap_lme_f_winter %>% select(variable, variable_ranked, Trial), by = c("variable","Trial")) %>%
  filter(!is.na(variable_ranked)) %>%
  mutate(
    chr      = factor(chr, levels = sort(unique(chr))),
    logP     = -log10(min_p),
    logP_cap = pmin(logP, 10)
  )

## Keep original within-panel rank order explicitly
heat_spring$variable_ranked <- factor(heat_spring$variable_ranked,
                                      levels = levels(shap_lme_f_spring$variable_ranked))
heat_winter$variable_ranked <- factor(heat_winter$variable_ranked,
                                      levels = levels(shap_lme_f_winter$variable_ranked))

df_all <- bind_rows(heat_spring, heat_winter)
rm_intersect <- intersect(
  heat_spring$variable %>% unique(),
  heat_winter$variable %>% unique()
)

df_all <- df_all %>%
  filter(!variable %in% rm_intersect)
## ------------------------------------------------------------
## 7) Filter OTUs: keep OTUs that have any significant SNPs within that trial
## ------------------------------------------------------------
df_all_sig <- df_all %>%
  group_by(Trial, variable) %>%
  filter(any(n_sig > 0)) %>%
  ungroup()

## ------------------------------------------------------------
## 8) Build final ordering (Spring block first, Winter block second)
##     - keep within-panel rank from variable_ranked levels
## ------------------------------------------------------------
spring_levels <- levels(shap_lme_f_spring$variable_ranked)
winter_levels <- levels(shap_lme_f_winter$variable_ranked)

## Only keep those that survived "has sig SNP" filter
spring_keep <- spring_levels[spring_levels %in% (df_all_sig %>% filter(Trial=="SpringPan") %>% pull(variable) %>% unique())]
winter_keep <- winter_levels[winter_levels %in% (df_all_sig %>% filter(Trial=="WinterPan") %>% pull(variable) %>% unique())]

row_levels <- c(
  paste("SpringPan", spring_keep, sep="__"),
  paste("WinterPan", winter_keep, sep="__")
)

## Labels for y-axis (you can replace this later with taxonomy labels)
row_to_label <- c(
  setNames(spring_keep, paste("SpringPan", spring_keep, sep="__")),
  setNames(winter_keep, paste("WinterPan", winter_keep, sep="__"))
)

## Add row_id to all plotting datasets
df_all_sig <- df_all_sig %>%
  mutate(
    row_id = paste(Trial, variable, sep="__"),
    row_id = factor(row_id, levels = row_levels),
    chr    = factor(chr, levels = as.character(1:7))
  )

## ------------------------------------------------------------
## 9) GWAS bubble plot (combined)
## ------------------------------------------------------------
p_gwas_bubble <- ggplot(df_all_sig, aes(x = chr, y = row_id)) +
  geom_point(
    data = filter(df_all_sig, n_sig > 0),
    aes(size = logP_cap, fill = factor(n_sig)),
    shape = 21, color = "grey20", alpha = 0.9
  ) +
  geom_point(
    data = filter(df_all_sig, n_sig == 0 & n_snps > 0),
    shape = 4, size = 1.4, stroke = 0.5, color = "grey50"
  ) +
  scale_x_discrete(name = "Chromosome", limits = as.character(1:7), drop = FALSE) +
  scale_y_discrete(limits = rev(row_levels), labels = row_to_label[rev(row_levels)]) +
  scale_fill_viridis(discrete = TRUE, name = "Number of significant SNPs") +
  scale_size_continuous(name = expression(-log[10](min~p[w]))) +
  labs(y = NULL, title = "Genetic control of biomarker OTUs (bubble)") +
  theme_bw(base_size = 12) +
  theme(axis.text.y = element_text(size = 6),
        panel.grid  = element_blank())

## ------------------------------------------------------------
## 10) SHAP plot (only OTUs retained in df_all_sig)
## ------------------------------------------------------------
df_fs <- shap_lme_all %>%
  filter(paste(Trial, variable, sep="__") %in% row_levels) %>%
  mutate(
    row_id = paste(Trial, variable, sep="__"),
    row_id = factor(row_id, levels = row_levels)
  )

shap_long_f <- shap_long_all %>%
  filter(paste(Trial, variable, sep="__") %in% row_levels) %>%
  mutate(
    row_id = paste(Trial, variable, sep="__"),
    row_id = factor(row_id, levels = row_levels)
  )

p_shap <- ggplot(df_fs, aes(x = 6 * mean_abs_shap, y = row_id)) +
  geom_bar(stat = "identity", fill = "grey", alpha = 0.3) +
  geom_quasirandom(
    data = shap_long_f,
    aes(x = shap_value, y = row_id, color = feature_value_clip),
    alpha = 0.6, size = 1.5, width = 0.2
  ) +
  geom_text(aes(label = round(mean_abs_shap, 3)), hjust = -0.1, size = 3.5) +
  scale_y_discrete(limits = rev(row_levels), labels = row_to_label[rev(row_levels)]) +
  labs(title = "SHAP based on XGBoost", y = NULL, x = "SHAP value") +
  theme_bw(base_size = 13) +
  theme(legend.position = "bottom")+
  scale_color_viridis_c(
    begin = 0.25, end = 0.85, option = "inferno",
    name = NULL,
    breaks = c(lower_clip, upper_clip),
    labels = c("Low", "High"),
    guide = guide_colorbar(
      title = NULL,
      label.position = "bottom",
      label.hjust = c(0, 1),
      barwidth = unit(80, "mm"),
      barheight = unit(2, "mm"),
      ticks = FALSE
    ))

## ------------------------------------------------------------
## 11) LME plot (same OTUs; hide y labels to avoid repetition)
## ------------------------------------------------------------
shap_lme_f2 <- shap_lme_f_all %>%
  filter(paste(Trial, variable, sep="__") %in% row_levels) %>%
  mutate(
    row_id = paste(Trial, variable, sep="__"),
    row_id = factor(row_id, levels = row_levels)
  )

p_lme <- ggplot(shap_lme_f2, aes(x = "OTU", y = row_id)) +
  geom_point(aes(size = logP, color = slope), alpha = 0.8, shape = 15) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "firebrick",
                        midpoint = 0, name = "Slope in LME") +
  scale_size_continuous(name = "-log10(p_adj)") +
  scale_y_discrete(limits = rev(row_levels), labels = NULL) +
  labs(x = NULL, y = NULL, title = "LME") +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.y  = element_blank(),
    plot.margin   = margin(5, 5, 5, 5),
    legend.position = "top"
  )

## ------------------------------------------------------------
## 12) Manhattan facets (trial-specific significant OTUs only)
##     - Grey for non-significant SNPs, red for significant SNPs
##     - Facet rows use row_id (Spring block then Winter block)
## ------------------------------------------------------------
sig_OTU_by_trial <- all_assoc_all %>%
  group_by(Trial, OTU) %>%
  summarise(has_sig = any(p_wald < p_cutoff), .groups = "drop")

assoc_sig_only <- all_assoc_all %>%
  left_join(sig_OTU_by_trial, by = c("Trial","OTU")) %>%
  filter(has_sig) %>%
  select(-has_sig)

chr_info_global <- all_assoc_all %>%
  group_by(chr) %>%
  summarise(chr_len = max(ps, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0),
         center = chr_start + chr_len/2)

man_dat <- assoc_sig_only %>%
  mutate(
    variable = OTU,
    row_id   = paste(Trial, variable, sep="__"),
    row_id   = factor(row_id, levels = row_levels)
  ) %>%
  left_join(chr_info_global %>% select(chr, chr_start), by = "chr") %>%
  mutate(
    x    = ps + chr_start,
    logp = -log10(p_wald),
    sig  = (p_wald < p_cutoff)
  ) %>%
  filter(row_id %in% row_levels)

p_gwas_manhattan <- ggplot(man_dat, aes(x = x, y = logp)) +
  ggrastr::geom_point_rast(
    data  = subset(man_dat, !sig),
    color = "grey75", size = 0.45, alpha = 0.6,
    raster.dpi = 300
  ) +
  ggrastr::geom_point_rast(
    data  = subset(man_dat, sig),
    aes(size = logp),
    color = "firebrick", alpha = 0.8,
    raster.dpi = 300
  ) +
  geom_hline(
    yintercept = -log10(p_cutoff),
    linetype = "dashed", linewidth = 0.35, color = "grey40"
  ) +
  facet_grid(
    rows = vars(row_id),
    drop = FALSE,
    labeller = labeller(row_id = row_to_label)
  ) +
  scale_x_continuous(
    breaks = chr_info_global$center,
    labels = chr_info_global$chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_size_continuous(range = c(0.8, 2.6), guide = "none") +
  labs(
    title = "Manhattan plots of biomarker OTUs (trial-specific genetic control)",
    x = "Chromosome",
    y = expression(-log[10](p[w]))
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.text.y = element_text(size = 7),
    strip.text.x = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

## ------------------------------------------------------------
## 13) Final combined figure (SHAP + LME + Manhattan)
## ------------------------------------------------------------
combined_plot <- p_shap + p_lme + p_gwas_manhattan +
  plot_layout(widths = c(5, 2, 8))

combined_plot

ggsave("fs_gwas.pdf", combined_plot, width = 12, height = 12, dpi = 300)

## Optional: save the bubble plot too
ggsave("fs_gwas_bubble.pdf", p_gwas_bubble, width = 7, height = 10, dpi = 300)
