library(readxl)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  getwd()
}
base_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)

clean_csv <- function(path) {
  d <- read.csv(path, check.names = FALSE)
  empty_names <- names(d) == ""
  if (any(empty_names)) {
    names(d)[empty_names] <- paste0("blank_", seq_len(sum(empty_names)))
  }
  d
}

taxon_lookup <- bind_rows(
  clean_csv(file.path(base_dir, "Bac_Aut/feature_selection/bac_spring_output.csv")) %>%
    transmute(Biomarker = variable_ranked, Organism = "Bacteria", Trial = "SpringPan", Taxon = Genus),
  clean_csv(file.path(base_dir, "Bac_Aut/feature_selection/bac_winter_output.csv")) %>%
    transmute(Biomarker = variable_ranked, Organism = "Bacteria", Trial = "WinterPan", Taxon = Genus),
  clean_csv(file.path(base_dir, "Fun_Aut/feature_selection/fun_spring_output.csv")) %>%
    transmute(Biomarker = variable_ranked, Organism = "Fungi", Trial = "SpringPan", Taxon = genus),
  clean_csv(file.path(base_dir, "Fun_Aut/feature_selection/fun_winter_output.csv")) %>%
    transmute(Biomarker = variable_ranked, Organism = "Fungi", Trial = "WinterPan", Taxon = genus)
)

feature_stats <- bind_rows(
  clean_csv(file.path(base_dir, "Bac_Aut/feature_selection/bac_spring_output.csv")) %>%
    transmute(
      Biomarker = variable_ranked, Organism = "Bacteria", Trial = "SpringPan", Taxon = Genus,
      Mean_abs_SHAP = mean_abs_shap,
      SHAP_cor = NA_real_,
      SHAP_sign = NA_character_,
      LME_slope = slope,
      LME_sign = slope_sign,
      LME_p = P_Value,
      LME_FDR = p_adj,
      Direction_match = NA_character_
    ),
  clean_csv(file.path(base_dir, "Bac_Aut/feature_selection/bac_winter_output.csv")) %>%
    transmute(
      Biomarker = variable_ranked, Organism = "Bacteria", Trial = "WinterPan", Taxon = Genus,
      Mean_abs_SHAP = mean_abs_shap,
      SHAP_cor = NA_real_,
      SHAP_sign = NA_character_,
      LME_slope = slope,
      LME_sign = slope_sign,
      LME_p = P_Value,
      LME_FDR = p_adj,
      Direction_match = NA_character_
    ),
  clean_csv(file.path(base_dir, "Fun_Aut/feature_selection/fun_spring_output.csv")) %>%
    transmute(
      Biomarker = variable_ranked, Organism = "Fungi", Trial = "SpringPan", Taxon = genus,
      Mean_abs_SHAP = mean_abs_shap,
      SHAP_cor = shap_cor,
      SHAP_sign = shap_sign,
      LME_slope = slope,
      LME_sign = slope_sign,
      LME_p = P_Value,
      LME_FDR = p_adj,
      Direction_match = direction_match
    ),
  clean_csv(file.path(base_dir, "Fun_Aut/feature_selection/fun_winter_output.csv")) %>%
    transmute(
      Biomarker = variable_ranked, Organism = "Fungi", Trial = "WinterPan", Taxon = genus,
      Mean_abs_SHAP = mean_abs_shap,
      SHAP_cor = shap_cor,
      SHAP_sign = shap_sign,
      LME_slope = slope,
      LME_sign = slope_sign,
      LME_p = P_Value,
      LME_FDR = p_adj,
      Direction_match = direction_match
    )
)

biomarker_stats <- bind_rows(
  clean_csv(file.path(base_dir, "Bac_Aut/mGWAS/biomarker_out_put_summary.csv")) %>%
    mutate(Organism = "Bacteria", Biomarker = variable),
  clean_csv(file.path(base_dir, "Fun_Aut/mGWAS/biomarker_out_put_summary.csv")) %>%
    mutate(Organism = "Fungi", Biomarker = variable)
) %>%
  group_by(Biomarker, Organism, Trial) %>%
  summarise(
    Mean_abs_SHAP = first(mean_abs_shap),
    SHAP_cor = first(shap_cor),
    SHAP_sign = first(shap_sign),
    LME_slope = first(slope),
    LME_sign = first(slope_sign),
    LME_p = first(P_Value),
    LME_FDR = first(p_adj),
    Direction_match = first(direction_match),
    .groups = "drop"
  ) %>%
  left_join(taxon_lookup, by = c("Biomarker", "Organism", "Trial"))

network_hits <- bind_rows(
  read_excel(file.path(base_dir, "Network/mgwas_sum.xlsx"), sheet = "bac") %>%
    mutate(Organism = "Bacteria", Biomarker = asv),
  read_excel(file.path(base_dir, "Network/mgwas_sum.xlsx"), sheet = "fun") %>%
    mutate(Organism = "Fungi", Biomarker = asv)
) %>%
  group_by(Organism, Trial, shap_sign, Biomarker) %>%
  summarise(
    In_network_summary = "Yes",
    Network_sig_SNP_n = n_distinct(rs),
    Network_lead_SNP = rs[which.min(p_wald)][1],
    Network_lead_p = min(p_wald, na.rm = TRUE),
    Network_SNP_list = paste(sort(unique(rs)), collapse = "; "),
    .groups = "drop"
  )

chr_summary <- bind_rows(
  read.delim(
    file.path(base_dir, "Bac_Aut/mGWAS/Spring_gwas_res/spring_ASV_chr_summary.tsv"),
    check.names = FALSE
  ) %>% mutate(Organism = "Bacteria", Trial = "SpringPan", Biomarker = ASV),
  read.delim(
    file.path(base_dir, "Bac_Aut/mGWAS/Winter_gwas_res/winter_ASV_chr_summary.tsv"),
    check.names = FALSE
  ) %>% mutate(Organism = "Bacteria", Trial = "WinterPan", Biomarker = ASV),
  read.delim(
    file.path(base_dir, "Fun_Aut/mGWAS/Spring_gwas_res/spring_OTU_chr_summary.tsv"),
    check.names = FALSE
  ) %>% mutate(Organism = "Fungi", Trial = "SpringPan", Biomarker = OTU),
  read.delim(
    file.path(base_dir, "Fun_Aut/mGWAS/Winter_gwas_res/winter_OTU_chr_summary.tsv"),
    check.names = FALSE
  ) %>% mutate(Organism = "Fungi", Trial = "WinterPan", Biomarker = OTU)
) %>%
  group_by(Organism, Trial, Biomarker) %>%
  summarise(
    Chr_tested = n_distinct(chr),
    Total_SNPs_tested = sum(n_snps, na.rm = TRUE),
    Sig_chr_n = sum(n_sig > 0, na.rm = TRUE),
    Sig_SNP_n = sum(n_sig, na.rm = TRUE),
    Min_p_chr = chr[which.min(min_p)][1],
    Min_chr_p = min(min_p, na.rm = TRUE),
    Sig_chr_list = paste(chr[n_sig > 0], collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    Sig_chr_list = ifelse(Sig_chr_list == "", "", Sig_chr_list),
    mGWAS_sig = ifelse(Sig_SNP_n > 0, "Yes", "No"),
    mGWAS_tested = "Yes"
  )

output <- feature_stats %>%
  left_join(
    biomarker_stats %>%
      rename(
        Mean_abs_SHAP_mgwas = Mean_abs_SHAP,
        SHAP_cor_mgwas = SHAP_cor,
        SHAP_sign_mgwas = SHAP_sign,
        LME_slope_mgwas = LME_slope,
        LME_sign_mgwas = LME_sign,
        LME_p_mgwas = LME_p,
        LME_FDR_mgwas = LME_FDR,
        Direction_match_mgwas = Direction_match
      ),
    by = c("Biomarker", "Organism", "Trial", "Taxon")
  ) %>%
  mutate(
    Mean_abs_SHAP = coalesce(Mean_abs_SHAP_mgwas, Mean_abs_SHAP),
    SHAP_cor = coalesce(SHAP_cor_mgwas, SHAP_cor),
    SHAP_sign = coalesce(SHAP_sign_mgwas, SHAP_sign),
    LME_slope = coalesce(LME_slope_mgwas, LME_slope),
    LME_sign = coalesce(LME_sign_mgwas, LME_sign),
    LME_p = coalesce(LME_p_mgwas, LME_p),
    LME_FDR = coalesce(LME_FDR_mgwas, LME_FDR),
    Direction_match = coalesce(Direction_match_mgwas, Direction_match)
  ) %>%
  left_join(chr_summary, by = c("Organism", "Trial", "Biomarker")) %>%
  left_join(
    network_hits,
    by = c("Organism", "Trial", "Biomarker", "SHAP_sign" = "shap_sign")
  ) %>%
  mutate(
    mGWAS_tested = replace_na(mGWAS_tested, "No"),
    mGWAS_sig = replace_na(mGWAS_sig, "No"),
    Chr_tested = replace_na(Chr_tested, 0L),
    Total_SNPs_tested = replace_na(Total_SNPs_tested, 0L),
    Sig_chr_n = replace_na(Sig_chr_n, 0L),
    Sig_SNP_n = replace_na(Sig_SNP_n, 0L),
    Sig_chr_list = replace_na(Sig_chr_list, ""),
    In_network_summary = replace_na(In_network_summary, "No"),
    Network_sig_SNP_n = replace_na(Network_sig_SNP_n, 0L),
    Network_lead_SNP = replace_na(Network_lead_SNP, ""),
    Network_SNP_list = replace_na(Network_SNP_list, ""),
    SHAP_sign = replace_na(SHAP_sign, ""),
    Direction_match = replace_na(Direction_match, "")
  ) %>%
  arrange(Organism, Trial, desc(SHAP_sign), Biomarker) %>%
  select(
    Biomarker, Organism, Trial, Taxon,
    Mean_abs_SHAP, SHAP_cor, SHAP_sign,
    LME_slope, LME_sign, LME_p, LME_FDR, Direction_match,
    mGWAS_tested, mGWAS_sig,
    Chr_tested, Total_SNPs_tested, Sig_chr_n, Sig_SNP_n,
    Min_p_chr, Min_chr_p, Sig_chr_list,
    In_network_summary, Network_sig_SNP_n, Network_lead_SNP, Network_lead_p, Network_SNP_list
  )

write.csv(
  output,
  file.path(base_dir, "supplementary_biomarker_mgwas_summary.csv"),
  row.names = FALSE,
  na = ""
)

message("Wrote ", nrow(output), " rows to supplementary_biomarker_mgwas_summary.csv")
