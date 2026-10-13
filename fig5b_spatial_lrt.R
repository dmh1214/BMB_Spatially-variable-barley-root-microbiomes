# Figures: Fig5b; FigS13
# Original source: 0_Submission/4_NatComm/Manuscript_Figure_Code/scripts/fig5b_spatial_lrt.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

# Prepared for reviewer use; original analysis retained in originals/.
if (!exists("project_root", inherits=FALSE)) {
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), winslash="/", mustWork=TRUE)
}
module_output <- file.path(project_root, "outputs", "fig5b")
dir.create(module_output, recursive=TRUE, showWarnings=FALSE)
setwd(module_output)
dir.create("Network/k_tables", recursive=TRUE, showWarnings=FALSE)
on_error <- NULL


suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(performance)
  library(ggthemes)
})


bac_spring_rclr <- readRDS(file.path(project_root, "data/processed/bacteria/asv_spring_rclr.rds"))
bac_winter_rclr <- readRDS(file.path(project_root, "data/processed/bacteria/asv_winter_rclr.rds"))
fun_spring_rclr <- readRDS(file.path(project_root, "data/processed/fungi/otu_spring_rclr.rds"))
fun_winter_rclr <- readRDS(file.path(project_root, "data/processed/fungi/otu_winter_rclr.rds"))

meta_bac_spring <- readRDS(file.path(project_root, "data/processed/bacteria/meta_spring.rds"))
meta_bac_winter <- readRDS(file.path(project_root, "data/processed/bacteria/meta_winter.rds"))
meta_fun_spring <- readRDS(file.path(project_root, "data/processed/fungi/meta_spring.rds"))
meta_fun_winter <- readRDS(file.path(project_root, "data/processed/fungi/meta_winter.rds"))

load(file.path(project_root, "data/processed/derived/pcnm_release_objects.rdata"))

network_node <- read.csv(file.path(project_root, "data/processed/network/mgwas_sum_node.csv"), check.names = FALSE)


fix_pcnm_ids <- function(X, prefer_cols = c("sampleid", "sampleid.x", "sample_id", "SampleID", "ID")) {
  if (!is.null(rownames(X)) && any(grepl("[A-Za-z]", rownames(X)))) {
    return(X)
  }

  hit <- prefer_cols[prefer_cols %in% colnames(X)][1]
  if (is.na(hit)) {
    stop(
      "No usable sample id column found in PCNM object. Available cols: ",
      paste(colnames(X), collapse = ", ")
    )
  }

  X$sampleid <- as.character(X[[hit]])
  rownames(X) <- X$sampleid
  X
}

select_pcnm_set <- function(X_pcnm, sel_vec, mode = c("full", "top10"), K_max = 10) {
  mode <- match.arg(mode)
  sel_vec <- as.character(sel_vec)

  miss <- setdiff(sel_vec, colnames(X_pcnm))
  if (length(miss) > 0) {
    stop("Selected PCNM names not found in X_pcnm: ", paste(miss, collapse = ", "))
  }

  pcnm_all_sel <- X_pcnm[, sel_vec, drop = FALSE]

  if (mode == "full") {
    return(list(pcnm_mat = pcnm_all_sel, n_used = ncol(pcnm_all_sel), mode_label = "full"))
  }

  n_used <- min(K_max, ncol(pcnm_all_sel))
  list(
    pcnm_mat = pcnm_all_sel[, seq_len(n_used), drop = FALSE],
    n_used = n_used,
    mode_label = paste0("top", n_used)
  )
}

safe_r2m <- function(model) {
  tryCatch(
    suppressWarnings(performance::r2_nakagawa(model)$R2_marginal),
    error = function(e) NA_real_
  )
}

empty_result <- function(biomarker_id, n = NA_integer_, status) {
  tibble(
    Biomarker = biomarker_id,
    n = n,
    lrt_chisq = NA_real_,
    lrt_df = NA_real_,
    p_lrt = NA_real_,
    minus_log10_p_lrt = NA_real_,
    r2m_m0 = NA_real_,
    r2m_m1 = NA_real_,
    delta_r2m = NA_real_,
    AIC_m0 = NA_real_,
    AIC_m1 = NA_real_,
    delta_AIC = NA_real_,
    BIC_m0 = NA_real_,
    BIC_m1 = NA_real_,
    delta_BIC = NA_real_,
    status = status
  )
}

fit_one_biomarker_lrt <- function(biomarker_id, rclr_mat, meta_df, pcnm_df,
                                  geno_col = "Entry", min_non_na = 30) {
  if (!biomarker_id %in% colnames(rclr_mat)) {
    return(empty_result(biomarker_id, status = "biomarker_not_in_matrix"))
  }

  df <- meta_df |>
    mutate(.y = as.numeric(rclr_mat[, biomarker_id, drop = TRUE])) |>
    bind_cols(pcnm_df) |>
    filter(!is.na(.y), !is.na(.data[[geno_col]]))

  n <- nrow(df)
  if (n < min_non_na) {
    return(empty_result(biomarker_id, n = n, status = "too_few_samples"))
  }

  pcnm_terms <- paste(colnames(pcnm_df), collapse = " + ")
  f0 <- as.formula(paste0(".y ~ 1 + (1|", geno_col, ")"))
  f1 <- as.formula(paste0(".y ~ ", pcnm_terms, " + (1|", geno_col, ")"))

  m0 <- tryCatch(
    suppressWarnings(lmer(f0, data = df, REML = FALSE)),
    error = function(e) NULL
  )
  if (is.null(m0)) {
    return(empty_result(biomarker_id, n = n, status = "m0_failed"))
  }

  m1 <- tryCatch(
    suppressWarnings(lmer(f1, data = df, REML = FALSE)),
    error = function(e) NULL
  )
  if (is.null(m1)) {
    return(empty_result(biomarker_id, n = n, status = "m1_failed"))
  }

  lrt <- anova(m0, m1)
  p_lrt <- lrt$`Pr(>Chisq)`[2]
  r2m_m0 <- safe_r2m(m0)
  r2m_m1 <- safe_r2m(m1)

  tibble(
    Biomarker = biomarker_id,
    n = n,
    lrt_chisq = lrt$Chisq[2],
    lrt_df = lrt$`Chi Df`[2],
    p_lrt = p_lrt,
    minus_log10_p_lrt = -log10(p_lrt),
    r2m_m0 = r2m_m0,
    r2m_m1 = r2m_m1,
    delta_r2m = r2m_m1 - r2m_m0,
    AIC_m0 = AIC(m0),
    AIC_m1 = AIC(m1),
    delta_AIC = AIC(m0) - AIC(m1),  # >0 supports M1
    BIC_m0 = BIC(m0),
    BIC_m1 = BIC(m1),
    delta_BIC = BIC(m0) - BIC(m1),  # >0 supports M1
    status = "ok"
  )
}

run_panel <- function(panel_name, organism, trial, biomarker_ids,
                      rclr_mat, meta_df, X_pcnm, sel_pcnm,
                      geno_col = "Entry", sample_col = "sampleid",
                      min_non_na = 30, mode = c("full", "top10"),
                      K_max = 10) {
  mode <- match.arg(mode)

  if (sample_col %in% colnames(meta_df)) {
    rownames(meta_df) <- as.character(meta_df[[sample_col]])
  }
  if (!geno_col %in% colnames(meta_df)) {
    stop(panel_name, ": geno_col not found: ", geno_col)
  }
  if (is.null(rownames(rclr_mat)) || any(nchar(rownames(rclr_mat)) == 0)) {
    stop(panel_name, ": rclr_mat must have sample IDs as rownames.")
  }

  psel <- select_pcnm_set(X_pcnm, sel_pcnm, mode = mode, K_max = K_max)
  pcnm_sel <- psel$pcnm_mat

  common_ids <- Reduce(intersect, list(rownames(rclr_mat), rownames(pcnm_sel), rownames(meta_df)))
  if (length(common_ids) < min_non_na) {
    stop(panel_name, " [", mode, "]: too few common samples after alignment.")
  }

  rclr_sub <- rclr_mat[common_ids, , drop = FALSE]
  meta_sub <- meta_df[common_ids, , drop = FALSE]
  pcnm_sub <- as.data.frame(pcnm_sel[common_ids, , drop = FALSE])
  meta_sub[[geno_col]] <- as.factor(meta_sub[[geno_col]])

  biomarker_ids <- intersect(unique(as.character(biomarker_ids)), colnames(rclr_sub))
  if (length(biomarker_ids) == 0) {
    warning(panel_name, " [", mode, "]: no positive biomarkers found in rCLR matrix.")
    return(tibble())
  }

  map_dfr(
    biomarker_ids,
    \(id) fit_one_biomarker_lrt(
      biomarker_id = id,
      rclr_mat = rclr_sub,
      meta_df = meta_sub,
      pcnm_df = pcnm_sub,
      geno_col = geno_col,
      min_non_na = min_non_na
    )
  ) |>
    mutate(
      panel = panel_name,
      Organism = organism,
      Trial = trial,
      pcnm_mode = mode,
      pcnm_mode_label = psel$mode_label,
      n_pcnm_used = psel$n_used,
      n_positive_biomarkers_requested = length(biomarker_ids)
    )
}


positive_biomarkers <- network_node |>
  filter(
    group %in% c("Bacteria", "Fungi"),
    lme_slope == "Positive",
    Trial %in% c("SpringPan", "WinterPan")
  ) |>
  transmute(
    Biomarker = trimws(as.character(Id)),
    Organism = group,
    Trial = Trial
  ) |>
  distinct()

positive_biomarkers


make_id_list <- function(organism, trial, mat_cols) {
  positive_biomarkers |>
    filter(Organism == organism, Trial == trial, Biomarker %in% trimws(mat_cols)) |>
    pull(Biomarker)
}

bac_ids_spring <- make_id_list("Bacteria", "SpringPan", colnames(bac_spring_rclr))
bac_ids_winter <- make_id_list("Bacteria", "WinterPan", colnames(bac_winter_rclr))
fun_ids_spring <- make_id_list("Fungi", "SpringPan", colnames(fun_spring_rclr))
fun_ids_winter <- make_id_list("Fungi", "WinterPan", colnames(fun_winter_rclr))

tibble(
  panel = c("Bac_Spring", "Bac_Winter", "Fun_Spring", "Fun_Winter"),
  Organism = c("Bacteria", "Bacteria", "Fungi", "Fungi"),
  Trial = c("SpringPan", "WinterPan", "SpringPan", "WinterPan"),
  n_positive_biomarkers_in_matrix = c(
    length(bac_ids_spring),
    length(bac_ids_winter),
    length(fun_ids_spring),
    length(fun_ids_winter)
  )
)


X_pcnm_spring_bac <- fix_pcnm_ids(X_pcnm_spring_bac)
X_pcnm_winter_bac <- fix_pcnm_ids(X_pcnm_winter_bac)
X_pcnm_spring_fun <- fix_pcnm_ids(X_pcnm_spring_fun)
X_pcnm_winter_fun <- fix_pcnm_ids(X_pcnm_winter_fun)

panels <- list(
  Bac_Spring = list(
    organism = "Bacteria", trial = "SpringPan", biomarker_ids = bac_ids_spring,
    rclr_mat = bac_spring_rclr, meta_df = meta_bac_spring,
    X_pcnm = X_pcnm_spring_bac, sel_pcnm = sel_pcnm_spring_bac
  ),
  Bac_Winter = list(
    organism = "Bacteria", trial = "WinterPan", biomarker_ids = bac_ids_winter,
    rclr_mat = bac_winter_rclr, meta_df = meta_bac_winter,
    X_pcnm = X_pcnm_winter_bac, sel_pcnm = sel_pcnm_winter_bac
  ),
  Fun_Spring = list(
    organism = "Fungi", trial = "SpringPan", biomarker_ids = fun_ids_spring,
    rclr_mat = fun_spring_rclr, meta_df = meta_fun_spring,
    X_pcnm = X_pcnm_spring_fun, sel_pcnm = sel_pcnm_spring_fun
  ),
  Fun_Winter = list(
    organism = "Fungi", trial = "WinterPan", biomarker_ids = fun_ids_winter,
    rclr_mat = fun_winter_rclr, meta_df = meta_fun_winter,
    X_pcnm = X_pcnm_winter_fun, sel_pcnm = sel_pcnm_winter_fun
  )
)


all_results <- imap_dfr(
  panels,
  function(p, nm) {
    message("\n========== Running ", nm, " (FULL) ==========")
    res_full <- run_panel(
      panel_name = nm,
      organism = p$organism,
      trial = p$trial,
      biomarker_ids = p$biomarker_ids,
      rclr_mat = p$rclr_mat,
      meta_df = p$meta_df,
      X_pcnm = p$X_pcnm,
      sel_pcnm = p$sel_pcnm,
      mode = "full"
    )

    message("\n========== Running ", nm, " (TOP10) ==========")
    res_top10 <- run_panel(
      panel_name = nm,
      organism = p$organism,
      trial = p$trial,
      biomarker_ids = p$biomarker_ids,
      rclr_mat = p$rclr_mat,
      meta_df = p$meta_df,
      X_pcnm = p$X_pcnm,
      sel_pcnm = p$sel_pcnm,
      mode = "top10",
      K_max = 10
    )

    bind_rows(res_full, res_top10)
  }
) |>
  group_by(panel, pcnm_mode) |>
  mutate(
    p_adj = p.adjust(p_lrt, method = "BH"),
    sig_fdr05 = !is.na(p_adj) & p_adj < 0.05,
    spatial_impact = case_when(
      is.na(delta_r2m) ~ NA_character_,
      delta_r2m >= 0.10 ~ "large",
      delta_r2m >= 0.03 ~ "medium",
      delta_r2m >= 0.01 ~ "small",
      TRUE ~ "tiny"
    )
  ) |>
  ungroup() |>
  arrange(Organism, Trial, pcnm_mode, p_lrt)

all_results


write.csv(
  all_results,
  "Positive_biomarker_PCNM_LRT_all_panels_full_vs_top10.csv",
  row.names = FALSE
)

pcnm_top10_for_join <- all_results |>
  filter(pcnm_mode == "top10", status == "ok") |>
  transmute(
    Organism,
    Trial,
    Biomarker,
    panel,
    pcnm_mode,
    n,
    n_pcnm_used,
    pcnm_p_lrt = p_lrt,
    pcnm_p_adj = p_adj,
    pcnm_minus_log10_p = minus_log10_p_lrt,
    pcnm_delta_r2m = delta_r2m,
    pcnm_delta_AIC = delta_AIC,
    pcnm_delta_BIC = delta_BIC,
    pcnm_sig_fdr05 = sig_fdr05,
    pcnm_spatial_impact = spatial_impact
  ) |>
  arrange(Organism, Trial, Biomarker)

write.csv(
  pcnm_top10_for_join,
  "Positive_biomarker_PCNM_LRT_top10_for_interaction_join.csv",
  row.names = FALSE
)

summary_table <- all_results |>
  filter(status == "ok") |>
  group_by(Organism, Trial, panel, pcnm_mode) |>
  summarise(
    n_tested = n(),
    n_sig_fdr05 = sum(sig_fdr05, na.rm = TRUE),
    prop_sig_fdr05 = n_sig_fdr05 / n_tested,
    median_delta_r2m = median(delta_r2m, na.rm = TRUE),
    median_delta_AIC = median(delta_AIC, na.rm = TRUE),
    median_delta_BIC = median(delta_BIC, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_table, "Positive_biomarker_PCNM_LRT_summary.csv", row.names = FALSE)

top_hits <- all_results |>
  filter(status == "ok") |>
  group_by(Organism, Trial, panel, pcnm_mode) |>
  slice_min(order_by = p_lrt, n = 10, with_ties = FALSE) |>
  ungroup() |>
  select(
    Organism, Trial, panel, pcnm_mode, Biomarker, n,
    p_lrt, p_adj, delta_r2m, delta_AIC, delta_BIC, spatial_impact
  )

write.csv(top_hits, "Positive_biomarker_PCNM_LRT_top_hits.csv", row.names = FALSE)


summary_table
top_hits
pcnm_top10_for_join


p_delta_r2_top10 <- pcnm_top10_for_join |>
  mutate(
    Biomarker = fct_reorder(Biomarker, pcnm_delta_r2m),
    panel = factor(panel, levels = c("Bac_Spring", "Bac_Winter", "Fun_Spring", "Fun_Winter"))
  ) |>
  ggplot(aes(Biomarker, pcnm_delta_r2m, fill = Organism)) +
  geom_col(width = 0.75, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  labs(
    x = NULL,
    y = expression(Delta * " marginal " * R^2),
    title = "Spatial dependency of positive biomarkers",
    subtitle = "Top10 PCNM model vs Entry-only null model"
  ) +
  theme_few() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = "top"
  )

p_delta_r2_top10

ggsave(
  "Positive_biomarker_PCNM_LRT_top10_deltaR2.pdf",
  p_delta_r2_top10,
  width = 9,
  height = 5
)

ggsave(
  "Positive_biomarker_PCNM_LRT_top10_deltaR2.png",
  p_delta_r2_top10,
  width = 9,
  height = 5,
  dpi = 300
)


pcnm_full_vs_top10 <- all_results |>
  filter(status == "ok") |>
  select(Organism, Trial, Biomarker, pcnm_mode, delta_r2m, p_lrt) |>
  pivot_wider(
    names_from = pcnm_mode,
    values_from = c(delta_r2m, p_lrt)
  )

p_full_vs_top10 <- pcnm_full_vs_top10 |>
  ggplot(aes(delta_r2m_full, delta_r2m_top10, color = Organism, shape = Trial)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    x = expression(Delta * " marginal " * R^2 * " (full PCNM)"),
    y = expression(Delta * " marginal " * R^2 * " (top10 PCNM)"),
    title = "Top10 PCNM captures the dominant spatial signal"
  ) +
  theme_few() +
  theme(legend.position = "top")

p_full_vs_top10

ggsave(
  "Positive_biomarker_PCNM_LRT_full_vs_top10_deltaR2.png",
  p_full_vs_top10,
  width = 6.5,
  height = 5,
  dpi = 300
)

ggsave(
  "Positive_biomarker_PCNM_LRT_full_vs_top10_deltaR2.pdf",
  p_full_vs_top10,
  width = 6.5,
  height = 5
)
