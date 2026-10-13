

#setwd("./mGWAS")
library(tidyverse)
library(readr)
library(ggrastr)


rm(list=ls())
load("../feature_selection/spring_biomarker.Rdata")
positive_biomarkers <- shap_lme_f %>%
  filter(shap_cor > 0) %>%
  pull(variable)
shap_lme <- shap_lme %>%
  filter(variable %in% positive_biomarkers)
shap_lme_f <- shap_lme_f %>%
  filter(variable %in% positive_biomarkers) %>%
  mutate(variable_ranked = droplevels(variable_ranked))
shap_long <- shap_long %>%
  filter(variable %in% positive_biomarkers)
shap_lme_spring <- shap_lme
shap_lme_f_spring <- shap_lme_f
shap_lme_f$variable
shap_long_spring <- shap_long
## 1. Basic settings ----------------------------------------------------
# Folder containing asv_xxx subfolders with GEMMA assoc outputs
outdir <- "Spring_gwas_res"

# P-value cutoff to count "significant" SNPs
p_cutoff <- 5e-5

## 2. List asv folders ---------------------------------------------------
asv_dirs <- list.dirs(outdir, recursive = FALSE, full.names = TRUE) %>%
  keep(~ grepl("^ASV", basename(.x)))

if (length(asv_dirs) == 0) {
  stop("No asv_ folders found under ", outdir)
}

## 3. Function: read one assoc and clean chr column ----------------------
read_one_assoc <- function(d) {
  asv <- basename(d)
  f <- list.files(d, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  if (length(f) == 0) return(NULL)
  
  df <- read.table(f[1], header = TRUE, stringsAsFactors = FALSE)
  
  # find position column name (GEMMA often uses "ps")
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
      ASV     = asv,
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

## 4. Read all assoc files -----------------------------------------------
all_assoc <- map_dfr(asv_dirs, read_one_assoc)
all_assoc_spring <- all_assoc
if (nrow(all_assoc) == 0) {
  stop("No valid SNP rows read from assoc files.")
}

## 5. Summarise per ASV × chromosome -------------------------------------
asv_chr_summary <- all_assoc %>%
  group_by(ASV, chr) %>%
  summarise(
    n_snps = n(),
    n_sig  = sum(p_wald < p_cutoff),
    min_p  = min(p_wald),
    .groups = "drop"
  )

# Save summary table (optional but useful)
write_tsv(asv_chr_summary,
          file.path(outdir, "spring_ASV_chr_summary.tsv"))

## 6. Prepare data for bubble plot ---------------------------------------

# shap_lme_f: has columns `variable` (asv_ID) and `variable_ranked` (factor)
chr_levels_num <- sort(unique(asv_chr_summary$chr))  # 一般就是 1,2,3,4,5,6,7

# asv_chr_summary: has column `ASV` (same asv_ID)
heat_dat <- asv_chr_summary %>%
  dplyr::rename(variable = ASV) %>%  # 对齐列名，用 variable 便于 join
  dplyr::left_join(
    shap_lme_f %>%
      dplyr::select(variable, variable_ranked),
    by = "variable"
  ) %>%
  dplyr::filter(!is.na(variable_ranked)) %>%  # 只保留 biomarker ASVs
  dplyr::mutate(
    # 不加 H，只用数字，但强制按照 chr_levels_num 的顺序
    chr      = factor(chr, levels = chr_levels_num),
    logP     = -log10(min_p),
    logP_cap = pmin(logP, 10)
  )
# Order ASVs by biomarker ranks (same as shap_lme_f)
heat_dat$variable_ranked <- factor(
  heat_dat$variable_ranked,
  levels = levels(shap_lme_f$variable_ranked)  # 保留原来的 level 顺序
)
heat_dat_spring<- heat_dat
## 7. Bubble heatmap ------------------------------------------------------
library(viridis)
p_bubble_spring<- ggplot(
  heat_dat,
  aes(
    # 这里再 **显式** 把 x 设成 factor，并指定顺序
    x = factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7")),
    y = variable_ranked
  )
) +
  # 1) 有显著 SNP（n_sig > 0）的圆圈
  geom_point(
    data = filter(heat_dat, n_sig > 0),
    aes(size = logP_cap, fill = factor(n_sig)),
    shape = 21,
    color = "grey20",
    alpha = 0.9
  ) +
  # 2) 有 SNP 但没有显著（n_sig == 0）
  geom_point(
    data = filter(heat_dat, n_sig == 0 & n_snps > 0),
    shape  = 4,
    size   = 1.4,
    stroke = 0.5,
    color  = "grey50"
  ) +
  scale_fill_viridis(discrete = TRUE, name = "Number of significant SNPs") +
  scale_size_continuous(
    name = expression(-log[10](min~p[w]))
  ) +
  scale_y_discrete(
    labels = asv_to_label_ranked,
    limits = rev(levels(shap_lme_f$variable_ranked))
  ) +
  # **再上第二层保险：强制 x 轴顺序**
  scale_x_discrete(
    name   = "Chromosome",
    limits = c("1", "2", "3", "4", "5", "6", "7")
  ) +
  labs(
    y     = "",
    title = "Genetic control of biomarker ASVs (SpringPan)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid  = element_blank()
  )

# Save plots
ggsave(file.path(outdir, "spring_ASV_chr_bubble_plot.pdf"),
       p_bubble_spring, width = 7, height = 8)
ggsave(file.path(outdir, "spring_ASV_chr_bubble_plot.png"),
       p_bubble_spring, width = 7, height = 8, dpi = 300)

print(p_bubble_spring)



#winter ----
## -----------------------------------------------------------
## winterPan mGWAS summary: ASV × Chromosome bubble plot
## - Reads all asv_*/ *.assoc.txt under winter_gwas_res
## - For each ASV & chromosome:
##     * n_sig  = # SNPs with p_wald < p_cutoff
##     * min_p  = minimum p_wald
## - Plot:
##     color = -log10(min_p)
##     size  = n_sig
## -----------------------------------------------------------

load("../feature_selection/winter_biomarker.Rdata")
positive_biomarkers <- shap_lme_f %>%
  filter(shap_cor > 0) %>%
  pull(variable)
shap_lme <- shap_lme %>%
  filter(variable %in% positive_biomarkers)
shap_lme_f <- shap_lme_f %>%
  filter(variable %in% positive_biomarkers) %>%
  mutate(variable_ranked = droplevels(variable_ranked))
shap_long <- shap_long %>%
  filter(variable %in% positive_biomarkers)
shap_lme_winter <- shap_lme
shap_lme_f_winter <- shap_lme_f
shap_lme_f$variable
shap_long_winter <- shap_long
## 1. Basic settings ----------------------------------------------------
# Folder containing asv_xxx subfolders with GEMMA assoc outputs
outdir <- "winter_gwas_res"

# P-value cutoff to count "significant" SNPs
p_cutoff <- 5e-5

## 2. List asv folders ---------------------------------------------------
asv_dirs <- list.dirs(outdir, recursive = FALSE, full.names = TRUE) %>%
  keep(~ grepl("^ASV", basename(.x)))

if (length(asv_dirs) == 0) {
  stop("No asv_ folders found under ", outdir)
}

## 3. Function: read one assoc and clean chr column ----------------------
rread_one_assoc <- function(d) {
  asv <- basename(d)
  f <- list.files(d, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  if (length(f) == 0) return(NULL)
  
  df <- read.table(f[1], header = TRUE, stringsAsFactors = FALSE)
  
  # find position column name (GEMMA often uses "ps")
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
      ASV     = asv,
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

## 4. Read all assoc files -----------------------------------------------
all_assoc <- map_dfr(asv_dirs, read_one_assoc)
all_assoc_winter <- all_assoc
if (nrow(all_assoc) == 0) {
  stop("No valid SNP rows read from assoc files.")
}

## 5. Summarise per ASV × chromosome -------------------------------------
asv_chr_summary <- all_assoc %>%
  group_by(ASV, chr) %>%
  summarise(
    n_snps = n(),
    n_sig  = sum(p_wald < p_cutoff),
    min_p  = min(p_wald),
    .groups = "drop"
  )

# Save summary table (optional but useful)
write_tsv(asv_chr_summary,
          file.path(outdir, "winter_ASV_chr_summary.tsv"))

## 6. Prepare data for bubble plot ---------------------------------------

# shap_lme_f: has columns `variable` (asv_ID) and `variable_ranked` (factor)
chr_levels_num <- sort(unique(asv_chr_summary$chr))  # 一般就是 1,2,3,4,5,6,7

# asv_chr_summary: has column `ASV` (same asv_ID)
heat_dat <- asv_chr_summary %>%
  dplyr::rename(variable = ASV) %>%  # 对齐列名，用 variable 便于 join
  dplyr::left_join(
    shap_lme_f %>%
      dplyr::select(variable, variable_ranked),
    by = "variable"
  ) %>%
  dplyr::filter(!is.na(variable_ranked)) %>%  # 只保留 biomarker ASVs
  dplyr::mutate(
    # 不加 H，只用数字，但强制按照 chr_levels_num 的顺序
    chr      = factor(chr, levels = chr_levels_num),
    logP     = -log10(min_p),
    logP_cap = pmin(logP, 10)
  )
# Order ASVs by biomarker ranks (same as shap_lme_f)
heat_dat$variable_ranked <- factor(
  heat_dat$variable_ranked,
  levels = levels(shap_lme_f$variable_ranked)  # 保留原来的 level 顺序
)
heat_dat_winter <- heat_dat
## 7. Bubble heatmap ------------------------------------------------------
library(viridis)
p_bubble_winter <- ggplot(
  heat_dat,
  aes(
    # 这里再 **显式** 把 x 设成 factor，并指定顺序
    x = factor(chr, levels = c("1", "2", "3", "4", "5", "6", "7")),
    y = variable_ranked
  )
) +
  # 1) 有显著 SNP（n_sig > 0）的圆圈
  geom_point(
    data = filter(heat_dat, n_sig > 0),
    aes(size = logP_cap, fill = factor(n_sig)),
    shape = 21,
    color = "grey20",
    alpha = 0.9
  ) +
  # 2) 有 SNP 但没有显著（n_sig == 0）
  geom_point(
    data = filter(heat_dat, n_sig == 0 & n_snps > 0),
    shape  = 4,
    size   = 1.4,
    stroke = 0.5,
    color  = "grey50"
  ) +
  scale_fill_viridis(discrete = TRUE, name = "Number of significant SNPs") +
  scale_size_continuous(
    name = expression(-log[10](min~p[w]))
  ) +
  scale_y_discrete(
    labels = asv_to_label_ranked,
    limits = rev(levels(shap_lme_f$variable_ranked))
  ) +
  # **再上第二层保险：强制 x 轴顺序**
  scale_x_discrete(
    name   = "Chromosome",
    limits = c("1", "2", "3", "4", "5", "6", "7")
  ) +
  labs(
    y     = "",
    title = "Genetic control of biomarker ASVs (WinterPan)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid  = element_blank()
  )

ggsave(file.path(outdir, "winter_ASV_chr_bubble_plot.pdf"),
       p_bubble_winter, width = 7, height = 8)
ggsave(file.path(outdir, "winter_ASV_chr_bubble_plot.png"),
       p_bubble_winter, width = 7, height = 8, dpi = 300)

print(p_bubble_winter)



# library(patchwork)
# 
# combined_two <- p_bubble_spring + p_bubble_winter +
#   plot_layout(guides = "collect", widths = c(1, 1)) &
#   theme(legend.position = "bottom")
# 
# combined_two

heat_dat_spring$Trial <- "SpringPan"
heat_dat_winter$Trial <- "WinterPan"

df_all <- rbind(heat_dat_spring, heat_dat_winter)

## 1. 先过滤掉“两个 Trial 都没有显著 SNP 的 ASV”
df_all_filtered <- df_all %>%
  group_by(variable_ranked) %>%
  filter(any(n_sig > 0)) %>%
  ungroup() %>%
  mutate(
    variable_ranked = droplevels(variable_ranked)  # ★ 把没用到的 level 真的删掉
  )

## 2. 取过滤后的 level 顺序
var_levels <- levels(df_all_filtered$variable_ranked)

p <- ggplot(
  df_all_filtered,
  aes(x = chr, y = variable_ranked)
) +
  geom_point(
    data = filter(df_all_filtered, n_sig > 0),
    aes(size = logP_cap, fill = factor(n_sig)),
    shape = 21,
    color = "grey20",
    alpha = 0.9
  ) +
  geom_point(
    data = filter(df_all_filtered, n_sig == 0 & n_snps > 0),
    shape  = 4,
    size   = 1.4,
    stroke = 0.5,
    color  = "grey50"
  ) +
  scale_x_discrete(
    name   = "Chromosome",
    limits = c("1","2","3","4","5","6","7"),
    drop   = FALSE
  ) +
  scale_y_discrete(
    labels = asv_to_label_ranked[var_levels],  # 只给保留下来的 ASV 贴标签
    limits = rev(var_levels)                   # 只包含保留下来的 levels
  ) +
  scale_fill_viridis(discrete = TRUE, name = "Number of significant SNPs") +
  scale_size_continuous(
    name = expression(-log[10](min~p[w]))
  ) +
  labs(
    y     = NULL,
    title = "Genetic control of biomarker ASVs"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 6),
    panel.grid  = element_blank()
  )

p

shap_lme_spring$Trial <- "SpringPan"
shap_lme_winter$Trial <- "WinterPan"

shap_lme_all <- rbind(shap_lme_spring,shap_lme_winter)

shap_lme_f_spring$Trial <- "SpringPan"
shap_lme_f_winter$Trial <- "WinterPan"

shap_lme_f_all <- rbind(shap_lme_f_spring,shap_lme_f_winter)

res_df <- shap_lme_f_all |> left_join(df_all,join_by(variable,Trial)) |> 
  select(variable,mean_abs_shap,shap_cor,P_Value,p_adj,logP.x,sig_group,slope,shap_sign,slope_sign,direction_match,Trial,chr,n_snps,min_p,logP.y)

write.csv(res_df,"biomarker_out_put_summary.csv")

shap_long_spring$Trial <- "SpringPan"
shap_long_winter$Trial <- "WinterPan"


shap_long <- rbind(shap_long_spring,shap_long_winter)

# df_fs: SHAP summary，只保留在 GWAS 里有显著 SNP 的 ASV
df_fs <- shap_lme_all |>
  filter(variable %in% df_all_filtered$variable) |>
  mutate(
    variable_ranked = factor(variable, levels = var_levels)
  )



# shap_long_f: SHAP long 格式，同样只保留这些 ASV
shap_long_f <- shap_long |>
  filter(variable %in% df_all_filtered$variable) |>
  mutate(
    variable_ranked = factor(variable, levels = var_levels)
  )


library(ggplot2)
library(ggbeeswarm)   # geom_quasirandom 用到这个包
library(grid)         # unit() 用到

p_shap <- ggplot(df_fs, aes(x = 6 * mean_abs_shap, y = variable_ranked)) +
  geom_bar(stat = "identity", fill = "grey", alpha = 0.3) +
  geom_quasirandom(
    data = shap_long_f,   # ★ 用过滤后的 shap_long_f
    aes(x = shap_value, y = variable_ranked, color = feature_value_clip),
    alpha = 0.6, size = 1.5, width = 0.2
  ) +
  geom_text(
    aes(label = round(mean_abs_shap, 3)),
    hjust = -0.1, size = 3.5
  ) +
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
    )
  ) +
  scale_y_discrete(
    labels = ASV_to_label_ranked[var_levels],  # ★ 用同一个 ASV 列表
    limits = rev(var_levels)                   # ★ 和 GWAS 图完全一致
  ) +
  labs(
    title = "SHAP based on XGBoost",
    y = NULL,
    x = "SHAP value"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "bottom"
  )
p_shap 

shap_lme_f2 <- shap_lme_f_all |>
  filter(variable %in% df_all_filtered$variable) |>
  mutate(
    variable_ranked = factor(variable_ranked, levels = var_levels)
  )

p_lme <- ggplot(shap_lme_f2, aes(x = "ASV", y = variable_ranked)) +
  geom_point(aes(size = logP, color = slope), alpha = 0.8, shape = 15) +
  scale_color_gradient2(
    low = "blue", mid = "grey", high = "firebrick",
    midpoint = 0,
    name = "Slope in LME"
  ) +
  scale_size_continuous(name = "-log10(p_adj)") +
  scale_y_discrete(
    labels = NULL,              # 这里不要再显示 label，避免和两边重复
    limits = rev(var_levels)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "LME"
  ) +
  theme_bw(base_size = 13) +
  theme(
    axis.text.x   = element_blank(),
    axis.text.y   = element_blank(),
    axis.ticks.x  = element_blank(),
    axis.ticks.y  = element_blank(),
    plot.margin   = margin(5, 5, 5, 5),
    legend.position = "top"
  )

# p_gwas <- ggplot(
#   df_all_filtered,
#   aes(x = chr, y = variable_ranked)
# ) +
#   geom_point(
#     data = filter(df_all_filtered, n_sig > 0),
#     aes(size = logP_cap, fill = factor(n_sig)),
#     shape = 21,
#     color = "grey20",
#     alpha = 0.9
#   ) +
#   geom_point(
#     data = filter(df_all_filtered, n_sig == 0 & n_snps > 0),
#     shape  = 4,
#     size   = 1.4,
#     stroke = 0.5,
#     color  = "grey50"
#   ) +
#   scale_x_discrete(
#     name   = "Chromosome",
#     limits = c("1","2","3","4","5","6","7"),
#     drop   = FALSE
#   ) +
#   scale_y_discrete(
#     labels = ASV_to_label_ranked[var_levels],
#     limits = rev(var_levels)
#   ) +
#   scale_fill_viridis(discrete = TRUE, name = "Number of significant SNPs") +
#   scale_size_continuous(
#     name = expression(-log[10](min~p[w]))
#   ) +
#   labs(
#     y     = NULL,
#     title = "Genetic control of biomarker ASVs"
#   ) +
#   theme_bw(base_size = 13) +
#   theme(
#     axis.text.y = element_text(size = 6),
#     panel.grid  = element_blank(),
#     plot.margin = margin(5, 5, 5, 5)
#   )



make_manhattan_facets <- function(all_assoc_trial, trial_name,
                                  var_levels, ASV_to_label_ranked,
                                  p_cutoff = 5e-5) {
  
  # keep only selected ASVs
  dat <- all_assoc_trial %>%
    rename(variable = ASV) %>%
    filter(variable %in% var_levels) %>%
    mutate(
      Trial = trial_name,
      variable_ranked = factor(variable, levels = var_levels),
      chr = as.integer(chr)
    )
  
  # build cumulative genome position
  chr_info <- dat %>%
    group_by(chr) %>%
    summarise(chr_len = max(ps, na.rm = TRUE), .groups = "drop") %>%
    arrange(chr) %>%
    mutate(chr_start = lag(cumsum(chr_len), default = 0))
  
  dat2 <- dat %>%
    left_join(chr_info, by = "chr") %>%
    mutate(x = ps + chr_start,
           logp = -log10(p_wald))
  
  # x-axis tick position per chr
  axis_df <- chr_info %>%
    mutate(center = chr_start + chr_len / 2)
  
  ggplot(dat2, aes(x = x, y = logp)) +
    geom_point(aes(group = chr), size = 0.5, alpha = 0.6) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed") +
    facet_wrap(~ variable_ranked, ncol = 2, scales = "free_x",
               labeller = labeller(variable_ranked = ASV_to_label_ranked)) +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$chr,
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    labs(
      title = paste0("Manhattan (", trial_name, ")"),
      x = "Chromosome",
      y = expression(-log[10](p[w]))
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(size = 8)
    )
}
# var_levels 你已经有了：levels(df_all_filtered$variable_ranked)
var_levels <- levels(df_all_filtered$variable_ranked)
all_assoc_spring$Trial <- "SpringPan"
all_assoc_winter$Trial <- "WinterPan"

all_assoc_all <- bind_rows(all_assoc_spring, all_assoc_winter)


sig_asv_by_trial <- all_assoc_all %>%
  group_by(Trial, ASV) %>%
  summarise(has_sig = any(p_wald < p_cutoff), .groups = "drop")

# 只保留 Trial 内显著的 ASV 的 SNP 点
assoc_sig_only <- all_assoc_all %>%
  left_join(sig_asv_by_trial, by = c("Trial","ASV")) %>%
  filter(has_sig) %>%
  select(-has_sig)


# 统一的 chr length（跨两个 Trial）
chr_info_global <- all_assoc_all %>%
  group_by(chr) %>%
  summarise(chr_len = max(ps, na.rm = TRUE), .groups = "drop") %>%
  arrange(chr) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0),
         center = chr_start + chr_len/2)

# 变成 genome x
man_dat <- assoc_sig_only %>%
  filter(ASV %in% var_levels) %>%          # ★ 先限制
  mutate(
    variable = ASV,
    variable_ranked = factor(variable, levels = var_levels),
    Trial = factor(Trial, levels = c("SpringPan","WinterPan"))
  ) %>%
  left_join(chr_info_global %>% select(chr, chr_start), by = "chr") %>%
  mutate(x = ps + chr_start,
         logp = -log10(p_wald)) %>%
  filter(!is.na(variable_ranked))          # ★ 再保险一次

man_dat <- man_dat %>%
  mutate(
    sig = ifelse(p_wald < p_cutoff, "Significant", "Not significant"),
    sig = factor(sig, levels = c("Not significant", "Significant"))
  )

man_dat <- man_dat %>%
  mutate(
    sig  = p_wald < p_cutoff,
    logp = -log10(p_wald)
  )

p_gwas_new <- ggplot(man_dat, aes(x = x, y = logp)) +
  ## ① 不显著 SNP：灰色、固定小点（背景）
  ggrastr::geom_point_rast(
    data = subset(man_dat, !sig),
    color = "grey75",
    size  = 0.45,
    alpha = 0.6,
    raster.dpi = 300
  ) +
  
  ## ② 显著 SNP：红色，size = -log10(p)
  ggrastr::geom_point_rast(
    data = subset(man_dat, sig),
    aes(size = logp),
    color = "firebrick",
    alpha = 0.8,
    raster.dpi = 300
  ) +
  
  ## ③ 显著性阈值线
  geom_hline(
    yintercept = -log10(p_cutoff),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey40"
  ) +
  
  ## ④ facet：ASV × Trial（与你 SHAP / LME 完全对齐）
  facet_grid(
    rows = vars(variable_ranked),
    drop = FALSE,
    labeller = labeller(variable_ranked = ASV_to_label_ranked)
  ) +
  
  ## ⑤ 染色体刻度
  scale_x_continuous(
    breaks = chr_info_global$center,
    labels = chr_info_global$chr,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  ## ⑥ size 映射范围（很关键！）
  scale_size_continuous(
    range = c(0.8, 2.6),   # 可微调
    guide = "none"        # ★ 去掉 size 图例
  ) +
  
  labs(
    title = "Manhattan plots of biomarker ASVs (trial-specific genetic control)",
    x = "Chromosome",
    y = expression(-log[10](p[w]))
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",        # ★ 总保险
    strip.text.y = element_text(size = 7),
    strip.text.x = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

p_gwas_new 


library(patchwork)

# combined_plot <- p_shap + p_lme + p_gwas +
#   plot_layout(widths = c(9, 1.5, 6))
combined_plot <- p_shap + p_lme + p_gwas_new +
  plot_layout(widths = c(7, 1, 8))

combined_plot

plot_height <- max(8, length(var_levels) * 1.2)

ggsave("fs_gwas_pos.pdf", combined_plot,
       width = 12, height = 6, dpi = 300,
       useDingbats = FALSE)
ggsave("fs_gwas_pos.png", combined_plot,
       width = 12, height = plot_height, dpi = 300)
