args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
project_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

knitr::opts_knit$set(root.dir = project_root)
knitr::opts_chunk$set(error = FALSE)

rm(list = setdiff(ls(), c("args", "script_arg", "script_dir", "project_root")))

library(tidyverse)
library(broom)
library(patchwork)

out_dir <- file.path(script_dir, "SNP_biomarker_PCNM_interaction_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

network_node <- read.csv(file.path(project_root, "Network/mgwas_sum_node.csv"))
network_edge <- read.csv(file.path(project_root, "Network/mgwas_sum_edge.csv"), check.names = FALSE) |>
  as_tibble(.name_repair = "unique")

bac_spring_rclr <- readRDS(file.path(project_root, "Bac_Aut/1_processed_data/asv_spring_rclr.rds"))
bac_winter_rclr <- readRDS(file.path(project_root, "Bac_Aut/1_processed_data/asv_winter_rclr.rds"))
fun_spring_rclr <- readRDS(file.path(project_root, "Fun_Aut/1_processed_data/otu_spring_rclr.rds"))
fun_winter_rclr <- readRDS(file.path(project_root, "Fun_Aut/1_processed_data/otu_winter_rclr.rds"))

meta_bac_spring <- readRDS(file.path(project_root, "Bac_Aut/1_processed_data/meta_spring.rds"))
meta_bac_winter <- readRDS(file.path(project_root, "Bac_Aut/1_processed_data/meta_winter.rds"))
meta_fun_spring <- readRDS(file.path(project_root, "Fun_Aut/1_processed_data/meta_spring.rds"))
meta_fun_winter <- readRDS(file.path(project_root, "Fun_Aut/1_processed_data/meta_winter.rds"))

snp_bac <- readRDS(file.path(project_root, "Bac_Aut/1_processed_data/snp_geno_alleles_aut.rds"))
snp_fun <- readRDS(file.path(project_root, "Fun_Aut/1_processed_data/snp_geno_alleles_aut.rds"))

load(file.path(project_root, "data/processed/derived/pcnm_release_objects.rdata"))

clean_term_name <- function(x) {
  x |>
    str_replace_all("`", "") |>
    str_replace_all(":", "_X_") |>
    str_replace_all("[^A-Za-z0-9_]+", "_")
}

classify_model_term <- function(term) {
  case_when(
    term == "(Intercept)" ~ "intercept",
    str_detect(term, "^PCNM\\d+$") ~ "pcnm_main",
    str_detect(term, ":") & str_detect(term, "PCNM\\d+") ~ "snp_pcnm_interaction",
    str_detect(term, "^SNP\\d+") ~ "snp_main",
    TRUE ~ "other"
  )
}

get_pcnm_variable <- function(term) {
  str_extract(term, "PCNM\\d+")
}

get_snp_variable <- function(term) {
  str_extract(term, "SNP\\d+")
}

fix_pcnm_ids <- function(X, prefer_cols = c("sampleid", "sampleid.x", "sample_id", "SampleID", "ID")) {
  X <- as.data.frame(X)

  if ("sampleid" %in% colnames(X)) {
    X$sampleid <- as.character(X$sampleid)
    rownames(X) <- X$sampleid
    return(X)
  }

  if (!is.null(rownames(X)) && any(grepl("[A-Za-z]", rownames(X)))) {
    X$sampleid <- rownames(X)
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

select_top_pcnm <- function(X_pcnm, sel_vec, k_max = 10) {
  sel_vec <- as.character(sel_vec)
  sel_vec <- sel_vec[sel_vec %in% colnames(X_pcnm)]

  if (length(sel_vec) == 0) {
    return(tibble(sampleid = rownames(X_pcnm)))
  }

  keep <- head(sel_vec, k_max)
  X_tbl <- if ("sampleid" %in% colnames(X_pcnm)) {
    as_tibble(X_pcnm)
  } else {
    as_tibble(X_pcnm, rownames = "sampleid")
  }

  X_tbl |>
    select(sampleid, all_of(keep))
}

read_snp_genotypes <- function(snp_obj, organism) {
  if (!"Entry" %in% colnames(snp_obj)) {
    snp_obj <- snp_obj |>
      as.data.frame() |>
      tibble::rownames_to_column("Entry")
  } else {
    snp_obj <- as.data.frame(snp_obj)
  }

  snp_obj |>
    as_tibble(.name_repair = "unique") |>
    mutate(Entry_key = str_to_lower(Entry)) |>
    pivot_longer(
      cols = -c(Entry, Entry_key),
      names_to = "SNP",
      values_to = "allele_group"
    ) |>
    mutate(
      Organism = organism,
      allele_group = str_trim(as.character(allele_group)),
      allele_group = na_if(allele_group, ""),
      allele_group = if_else(
        str_to_upper(allele_group) %in% c("N", "NA", "N/A", "-", ".", "?"),
        NA_character_,
        allele_group
      )
    ) |>
    filter(!is.na(allele_group))
}

get_sig_biomarker_snp_pairs <- function(network_node,
                                        network_edge,
                                        p_cutoff = 5e-5,
                                        max_snps_per_biomarker = 2) {
  positive_biomarkers <- network_node |>
    filter(
      group %in% c("Bacteria", "Fungi"),
      lme_slope == "Positive"
    ) |>
    select(Id, group, Trial)

  network_edge |>
    filter(
      shap_sign == "Positive",
      as.numeric(mGWAS_pw) < p_cutoff
    ) |>
    transmute(
      Trial,
      Biomarker = Target,
      SNP = Source,
      p_wald = as.numeric(mGWAS_pw)
    ) |>
    inner_join(
      positive_biomarkers |>
        rename(Biomarker = Id, Organism = group),
      by = c("Biomarker", "Trial")
    ) |>
    group_by(Organism, Trial, Biomarker, SNP) |>
    summarise(p_wald = min(p_wald, na.rm = TRUE), .groups = "drop") |>
    group_by(Organism, Trial, Biomarker) |>
    arrange(p_wald, SNP, .by_group = TRUE) |>
    slice_head(n = max_snps_per_biomarker) |>
    ungroup()
}

extract_abundance_long <- function(rclr_mat, meta_df, organism, trial) {
  rclr_df <- as.data.frame(rclr_mat)

  if (!"sampleid" %in% colnames(rclr_df)) {
    rclr_df <- rclr_df |>
      tibble::rownames_to_column("sampleid")
  }

  rclr_df |>
    as_tibble(.name_repair = "minimal") |>
    pivot_longer(
      cols = -sampleid,
      names_to = "Biomarker",
      values_to = "abundance_rclr"
    ) |>
    inner_join(
      meta_df |>
        as_tibble() |>
        select(sampleid, Entry),
      by = "sampleid"
    ) |>
    mutate(
      Organism = organism,
      Trial = trial,
      Entry_key = str_to_lower(Entry)
    )
}

build_biomarker_model_data <- function(pair_df, abundance_df, snp_df, pcnm_df) {
  snps <- pair_df |>
    arrange(p_wald, SNP) |>
    pull(SNP)

  base_df <- abundance_df |>
    filter(
      Organism == pair_df$Organism[[1]],
      Trial == pair_df$Trial[[1]],
      Biomarker == pair_df$Biomarker[[1]]
    ) |>
    distinct(sampleid, .keep_all = TRUE) |>
    left_join(
      pcnm_df |>
        filter(
          Organism == pair_df$Organism[[1]],
          Trial == pair_df$Trial[[1]]
        ) |>
        select(-Organism, -Trial),
      by = "sampleid"
    ) |>
    select(sampleid, Entry, Entry_key, Organism, Trial, Biomarker, abundance_rclr, everything())

  for (i in seq_along(snps)) {
    snp_i <- snps[[i]]
    snp_col <- paste0("SNP", i)

    allele_df <- snp_df |>
      filter(
        Organism == pair_df$Organism[[1]],
        SNP == snp_i
      ) |>
      select(Entry_key, allele_group) |>
      rename(!!snp_col := allele_group)

    base_df <- base_df |>
      left_join(allele_df, by = "Entry_key")
  }

  attr(base_df, "snps") <- snps
  panel_pcnm_cols <- names(base_df)[
    str_detect(names(base_df), "^PCNM\\d+$") &
      purrr::map_lgl(names(base_df), ~ any(!is.na(base_df[[.x]])))
  ]
  attr(base_df, "pcnm_cols") <- panel_pcnm_cols
  base_df
}

make_model_formula <- function(snp_cols, pcnm_cols) {
  main_terms <- c(snp_cols, pcnm_cols)
  interaction_terms <- expand.grid(snp_cols, pcnm_cols) |>
    transmute(term = paste(Var1, Var2, sep = ":")) |>
    pull(term)

  as.formula(
    paste(
      "abundance_rclr ~",
      paste(c(main_terms, interaction_terms), collapse = " + ")
    )
  )
}

make_reduced_formula <- function(snp_cols, pcnm_cols) {
  as.formula(
    paste(
      "abundance_rclr ~",
      paste(c(snp_cols, pcnm_cols), collapse = " + ")
    )
  )
}

fit_one_biomarker_pcnm <- function(pair_df,
                                   abundance_df,
                                   snp_df,
                                   pcnm_df,
                                   min_n = 30,
                                   min_per_allele = 5) {
  model_df <- build_biomarker_model_data(pair_df, abundance_df, snp_df, pcnm_df)
  snps <- attr(model_df, "snps")
  snp_cols <- paste0("SNP", seq_along(snps))
  pcnm_cols <- attr(model_df, "pcnm_cols")

  if (length(pcnm_cols) == 0) {
    return(list(
      summary = tibble(
        Organism = pair_df$Organism[[1]],
        Trial = pair_df$Trial[[1]],
        Biomarker = pair_df$Biomarker[[1]],
        SNPs = paste(snps, collapse = " / "),
        PCNMs = "",
        n_snps = length(snps),
        n_pcnm = 0L,
        n = 0L,
        interaction_p = NA_real_,
        interaction_p_adj = NA_real_,
        delta_adj_r2 = NA_real_,
        full_adj_r2 = NA_real_,
        status = "no_selected_pcnm"
      ),
      terms = tibble(),
      data = model_df,
      model = NULL
    ))
  }

  model_df <- model_df |>
    filter(
      is.finite(abundance_rclr),
      if_all(all_of(pcnm_cols), is.finite)
    )

  missing_snps <- snp_cols[
    purrr::map_lgl(snp_cols, ~ all(is.na(model_df[[.x]]) | model_df[[.x]] == ""))
  ]

  if (length(missing_snps) > 0) {
    missing_snp_names <- snps[match(missing_snps, snp_cols)]
    return(list(
      summary = tibble(
        Organism = pair_df$Organism[[1]],
        Trial = pair_df$Trial[[1]],
        Biomarker = pair_df$Biomarker[[1]],
        SNPs = paste(snps, collapse = " / "),
        PCNMs = paste(pcnm_cols, collapse = " / "),
        n_snps = length(snps),
        n_pcnm = length(pcnm_cols),
        n = 0L,
        interaction_p = NA_real_,
        interaction_p_adj = NA_real_,
        delta_adj_r2 = NA_real_,
        full_adj_r2 = NA_real_,
        status = paste0("snp_not_in_genotype_table: ", paste(missing_snp_names, collapse = " / "))
      ),
      terms = tibble(),
      data = model_df,
      model = NULL
    ))
  }

  model_df <- model_df |>
    filter(if_all(all_of(snp_cols), ~ !is.na(.x) & .x != "")) |>
    mutate(across(all_of(snp_cols), factor))

  if (nrow(model_df) < min_n) {
    return(list(
      summary = tibble(
        Organism = pair_df$Organism[[1]],
        Trial = pair_df$Trial[[1]],
        Biomarker = pair_df$Biomarker[[1]],
        SNPs = paste(snps, collapse = " / "),
        PCNMs = paste(pcnm_cols, collapse = " / "),
        n_snps = length(snps),
        n_pcnm = length(pcnm_cols),
        n = nrow(model_df),
        interaction_p = NA_real_,
        interaction_p_adj = NA_real_,
        delta_adj_r2 = NA_real_,
        full_adj_r2 = NA_real_,
        status = "too_few_samples"
      ),
      terms = tibble(),
      data = model_df,
      model = NULL
    ))
  }

  allele_ok <- purrr::map_lgl(snp_cols, function(snp_col) {
    counts <- table(model_df[[snp_col]])
    nlevels(model_df[[snp_col]]) >= 2 && all(counts >= min_per_allele)
  })

  if (!all(allele_ok)) {
    bad_snps <- snps[!allele_ok]
    return(list(
      summary = tibble(
        Organism = pair_df$Organism[[1]],
        Trial = pair_df$Trial[[1]],
        Biomarker = pair_df$Biomarker[[1]],
        SNPs = paste(snps, collapse = " / "),
        PCNMs = paste(pcnm_cols, collapse = " / "),
        n_snps = length(snps),
        n_pcnm = length(pcnm_cols),
        n = nrow(model_df),
        interaction_p = NA_real_,
        interaction_p_adj = NA_real_,
        delta_adj_r2 = NA_real_,
        full_adj_r2 = NA_real_,
        status = paste0("single_or_rare_allele: ", paste(bad_snps, collapse = " / "))
      ),
      terms = tibble(),
      data = model_df,
      model = NULL
    ))
  }

  reduced_formula <- make_reduced_formula(snp_cols, pcnm_cols)
  full_formula <- make_model_formula(snp_cols, pcnm_cols)

  reduced_model <- tryCatch(lm(reduced_formula, data = model_df), error = function(e) NULL)
  full_model <- tryCatch(lm(full_formula, data = model_df), error = function(e) NULL)

  if (is.null(reduced_model) || is.null(full_model)) {
    return(list(
      summary = tibble(
        Organism = pair_df$Organism[[1]],
        Trial = pair_df$Trial[[1]],
        Biomarker = pair_df$Biomarker[[1]],
        SNPs = paste(snps, collapse = " / "),
        PCNMs = paste(pcnm_cols, collapse = " / "),
        n_snps = length(snps),
        n_pcnm = length(pcnm_cols),
        n = nrow(model_df),
        interaction_p = NA_real_,
        interaction_p_adj = NA_real_,
        delta_adj_r2 = NA_real_,
        full_adj_r2 = NA_real_,
        status = "model_failed"
      ),
      terms = tibble(),
      data = model_df,
      model = NULL
    ))
  }

  model_cmp <- anova(reduced_model, full_model)
  snp_lookup <- setNames(snps, snp_cols)

  term_table <- broom::tidy(full_model) |>
    mutate(
      Organism = pair_df$Organism[[1]],
      Trial = pair_df$Trial[[1]],
      Biomarker = pair_df$Biomarker[[1]],
      SNPs = paste(snps, collapse = " / "),
      PCNMs = paste(pcnm_cols, collapse = " / "),
      term_clean = clean_term_name(term),
      term_type = classify_model_term(term),
      pcnm_variable = get_pcnm_variable(term),
      snp_variable = get_snp_variable(term),
      snp_name = unname(snp_lookup[snp_variable])
    ) |>
    select(
      Organism, Trial, Biomarker, SNPs, PCNMs,
      term, term_clean, term_type, pcnm_variable, snp_variable, snp_name,
      estimate, std.error, statistic, p.value
    )

  summary_row <- tibble(
    Organism = pair_df$Organism[[1]],
    Trial = pair_df$Trial[[1]],
    Biomarker = pair_df$Biomarker[[1]],
    SNPs = paste(snps, collapse = " / "),
    PCNMs = paste(pcnm_cols, collapse = " / "),
    n_snps = length(snps),
    n_pcnm = length(pcnm_cols),
    n = nrow(model_df),
    interaction_p = model_cmp$`Pr(>F)`[2],
    interaction_p_adj = NA_real_,
    delta_adj_r2 = summary(full_model)$adj.r.squared - summary(reduced_model)$adj.r.squared,
    full_adj_r2 = summary(full_model)$adj.r.squared,
    status = "ok"
  )

  list(
    summary = summary_row,
    terms = term_table,
    data = model_df,
    model = full_model
  )
}

pcnm_long <- bind_rows(
  select_top_pcnm(fix_pcnm_ids(X_pcnm_spring_bac), sel_pcnm_spring_bac, 10) |> mutate(Organism = "Bacteria", Trial = "SpringPan"),
  select_top_pcnm(fix_pcnm_ids(X_pcnm_winter_bac), sel_pcnm_winter_bac, 10) |> mutate(Organism = "Bacteria", Trial = "WinterPan"),
  select_top_pcnm(fix_pcnm_ids(X_pcnm_spring_fun), sel_pcnm_spring_fun, 10) |> mutate(Organism = "Fungi", Trial = "SpringPan"),
  select_top_pcnm(fix_pcnm_ids(X_pcnm_winter_fun), sel_pcnm_winter_fun, 10) |> mutate(Organism = "Fungi", Trial = "WinterPan")
)

sig_pairs <- get_sig_biomarker_snp_pairs(network_node, network_edge)

snp_long <- bind_rows(
  read_snp_genotypes(snp_bac, "Bacteria"),
  read_snp_genotypes(snp_fun, "Fungi")
)

abundance_long <- bind_rows(
  extract_abundance_long(bac_spring_rclr, meta_bac_spring, "Bacteria", "SpringPan"),
  extract_abundance_long(bac_winter_rclr, meta_bac_winter, "Bacteria", "WinterPan"),
  extract_abundance_long(fun_spring_rclr, meta_fun_spring, "Fungi", "SpringPan"),
  extract_abundance_long(fun_winter_rclr, meta_fun_winter, "Fungi", "WinterPan")
)

biomarker_models <- sig_pairs |>
  group_by(Organism, Trial, Biomarker) |>
  nest(pair_data = c(SNP, p_wald)) |>
  ungroup()

model_fits <- biomarker_models |>
  mutate(
    fit = purrr::pmap(
      list(Organism, Trial, Biomarker, pair_data),
      function(Organism, Trial, Biomarker, pair_data) {
        pair_df <- pair_data |>
          mutate(
            Organism = Organism,
            Trial = Trial,
            Biomarker = Biomarker
          ) |>
          select(Organism, Trial, Biomarker, SNP, p_wald)

        fit_one_biomarker_pcnm(pair_df, abundance_long, snp_long, pcnm_long)
      }
    ),
    summary = purrr::map(fit, "summary"),
    terms = purrr::map(fit, "terms"),
    model_data = purrr::map(fit, "data"),
    model = purrr::map(fit, "model")
  )

interaction_results <- model_fits |>
  select(summary) |>
  unnest(summary) |>
  group_by(Organism, Trial) |>
  mutate(interaction_p_adj = p.adjust(interaction_p, method = "BH")) |>
  ungroup() |>
  arrange(interaction_p)

term_effects_raw <- model_fits |>
  select(terms) |>
  unnest(terms)

term_effects <- if (nrow(term_effects_raw) == 0) {
  tibble(
    Organism = character(),
    Trial = character(),
    Biomarker = character(),
    SNPs = character(),
    PCNMs = character(),
    term = character(),
    term_clean = character(),
    term_type = character(),
    pcnm_variable = character(),
    snp_variable = character(),
    snp_name = character(),
    estimate = numeric(),
    std.error = numeric(),
    statistic = numeric(),
    p.value = numeric(),
    p_adj_by_term_type = numeric()
  )
} else {
  term_effects_raw |>
    group_by(Organism, Trial, term_type) |>
    mutate(p_adj_by_term_type = p.adjust(p.value, method = "BH")) |>
    ungroup()
}

interaction_summary <- interaction_results |>
  group_by(Organism, Trial) |>
  summarise(
    n_models = n(),
    n_tested = sum(status == "ok", na.rm = TRUE),
    n_nominal_p05 = sum(interaction_p < 0.05, na.rm = TRUE),
    n_fdr_p05 = sum(interaction_p_adj < 0.05, na.rm = TRUE),
    median_delta_adj_r2 = median(delta_adj_r2, na.rm = TRUE),
    .groups = "drop"
  )

term_summary <- if (nrow(term_effects) == 0) {
  tibble(
    Organism = character(),
    Trial = character(),
    n_terms = numeric(),
    n_nominal_p05 = numeric(),
    n_fdr_p05 = numeric()
  )
} else {
  term_effects |>
    filter(term_type == "snp_pcnm_interaction") |>
    group_by(Organism, Trial) |>
    summarise(
      n_terms = n(),
      n_nominal_p05 = sum(p.value < 0.05, na.rm = TRUE),
      n_fdr_p05 = sum(p_adj_by_term_type < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
}

interaction_summary
term_summary

p_term_summary <- if (nrow(term_summary) == 0) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No estimable SNP-by-PCNM interaction terms") +
    theme_void() +
    labs(title = "SNP-by-PCNM interaction terms with nominal P < 0.05")
} else {
  term_summary |>
    ggplot(aes(x = Trial, y = n_nominal_p05, fill = Organism)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    theme_bw(base_size = 12) +
    labs(
      title = "SNP-by-PCNM interaction terms with nominal P < 0.05",
      x = NULL,
      y = "Number of terms",
      fill = NULL
    )
}

p_term_summary

top_model_hits <- interaction_results |>
  filter(status == "ok", is.finite(interaction_p)) |>
  group_by(Organism, Trial) |>
  slice_min(interaction_p, n = 3, with_ties = FALSE) |>
  ungroup()

top_term_hits <- term_effects |>
  filter(
    term_type == "snp_pcnm_interaction",
    is.finite(p.value)
  ) |>
  arrange(p.value)

top_model_hits
top_term_hits |> print(n = 50)

p_top_terms <- if (nrow(top_term_hits) == 0) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No SNP-by-PCNM interaction terms available") +
    theme_void() +
    labs(title = "Top individual SNP-by-PCNM interaction terms")
} else {
  top_term_hits |>
    slice_head(n = 20) |>
    mutate(
      label = paste(Organism, Trial, Biomarker, term, sep = " | "),
      label = fct_reorder(label, -log10(p.value))
    ) |>
    ggplot(aes(x = -log10(p.value), y = label, fill = pcnm_variable)) +
    geom_col(width = 0.7) +
    theme_bw(base_size = 10) +
    labs(
      title = "Top individual SNP-by-PCNM interaction terms",
      x = "-log10(P)",
      y = NULL,
      fill = "PCNM"
    )
}

p_top_terms

if (nrow(top_term_hits) == 0) {
  p_interaction_effect_size <- ggplot() +
    annotate("text", x = 1, y = 1, label = "No SNP-by-PCNM interaction terms available") +
    theme_void() +
    labs(title = "Effect sizes of top SNP-by-PCNM interaction terms")
} else {
  plot_interaction_terms <- top_term_hits |>
    slice_head(n = 20) |>
    mutate(
      conf_low = estimate - 1.96 * std.error,
      conf_high = estimate + 1.96 * std.error,
      significant = p.value < 0.05,
      label = paste(Organism, Trial, Biomarker, snp_name, pcnm_variable, sep = " | "),
      label = fct_reorder(label, estimate)
    )

  p_interaction_effect_size <- ggplot(
    plot_interaction_terms,
    aes(x = estimate, y = label, colour = pcnm_variable, alpha = significant)
  ) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.45) +
    geom_errorbarh(
      aes(xmin = conf_low, xmax = conf_high),
      height = 0.18,
      linewidth = 0.65
    ) +
    geom_point(size = 2.4) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.35)) +
    theme_bw(base_size = 10) +
    labs(
      title = "Effect sizes of top SNP-by-PCNM interaction terms",
      subtitle = "Points show model coefficients; intervals are estimate +/- 1.96 SE",
      x = "Interaction estimate",
      y = NULL,
      colour = "PCNM axis",
      alpha = "P < 0.05"
    )
}

p_interaction_effect_size

top_slope_terms <- top_term_hits |>
  filter(
    !is.na(pcnm_variable),
    !is.na(snp_variable),
    p.value < 0.05
  ) |>
  arrange(p.value) |>
  select(
    Organism, Trial, Biomarker, term, term_type,
    pcnm_variable, snp_variable, snp_name, p.value
  )

top_slope_dat <- if (nrow(top_slope_terms) == 0) {
  tibble()
} else {
  model_fits |>
    semi_join(
      top_slope_terms |> select(Organism, Trial, Biomarker),
      by = c("Organism", "Trial", "Biomarker")
    ) |>
    select(Organism, Trial, Biomarker, model_data) |>
    mutate(
      model_data = purrr::map(model_data, ~ .x |> select(-any_of(c("Organism", "Trial", "Biomarker"))))
    ) |>
    unnest(model_data) |>
    inner_join(
      top_slope_terms,
      by = c("Organism", "Trial", "Biomarker"),
      relationship = "many-to-many"
    ) |>
    pivot_longer(
      cols = matches("^PCNM\\d+$"),
      names_to = "pcnm_current",
      values_to = "x_value"
    ) |>
    filter(pcnm_current == pcnm_variable) |>
    mutate(
      allele = case_when(
        snp_variable == "SNP1" ~ as.character(SNP1),
        snp_variable == "SNP2" ~ as.character(SNP2),
        TRUE ~ NA_character_
      ),
      panel = paste0(
        Organism, " ", Trial, "\n",
        Biomarker, " | ", snp_name, " | ", pcnm_variable,
        "\nP = ", signif(p.value, 2)
      )
    ) |>
    filter(!is.na(allele), is.finite(x_value))
}

p_interaction_slope_top <- if (nrow(top_slope_dat) == 0) {
  ggplot() +
    annotate("text", x = 1, y = 1, label = "No significant SNP-by-PCNM slopes available") +
    theme_void() +
    labs(title = "Allele-specific PCNM slopes")
} else {
  ggplot(
    top_slope_dat,
    aes(x = x_value, y = abundance_rclr, colour = allele)
  ) +
    geom_point(alpha = 0.42, size = 1.25) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
    facet_wrap(~ panel, scales = "free", ncol = 2) +
    theme_bw(base_size = 10) +
    theme(strip.text = element_text(size = 8)) +
    labs(
      title = "Significant allele-specific PCNM slopes reveal context-dependent recruitment",
      x = "PCNM value",
      y = "rCLR biomarker abundance",
      colour = "Allele"
    )
}

p_interaction_slope_top

write.csv(
  interaction_results,
  file.path(out_dir, "SNP_biomarker_PCNM_multisnp_model_summary.csv"),
  row.names = FALSE
)

write.csv(
  term_effects,
  file.path(out_dir, "SNP_biomarker_PCNM_multisnp_term_effects.csv"),
  row.names = FALSE
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_multisnp_term_summary.png"),
  p_term_summary,
  width = 10,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_multisnp_top_terms.png"),
  p_top_terms,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_multisnp_interaction_effect_sizes.png"),
  p_interaction_effect_size,
  width = 12,
  height = 8,
  dpi = 300
)

ggsave(
  file.path(out_dir, "SNP_biomarker_PCNM_multisnp_interaction_slopes.png"),
  p_interaction_slope_top,
  width = 12,
  height = 10,
  dpi = 300
)
