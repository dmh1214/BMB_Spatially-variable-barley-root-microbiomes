# Prepared for reviewer use; original analysis retained in originals/.
if (!exists("project_root", inherits=FALSE)) {
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), winslash="/", mustWork=TRUE)
}
module_output <- file.path(project_root, "outputs", "fig5d")
dir.create(module_output, recursive=TRUE, showWarnings=FALSE)
setwd(module_output)
dir.create("Network/k_tables", recursive=TRUE, showWarnings=FALSE)
on_error <- NULL

library(tidyverse)
library(vegan)

out_dir <- module_output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_png <- file.path(out_dir, "positive_biomarker_cv_vs_host_selection_individuals.png")
out_pdf <- file.path(out_dir, "positive_biomarker_cv_vs_host_selection_individuals.pdf")
out_csv <- file.path(out_dir, "positive_biomarker_cv_vs_host_selection_individuals_slopes.csv")

# Aggregate positive-pool CV already matched to host selection
agg_dat <- read.csv(
  file.path(project_root, "data", "processed", "derived", "positive_biomarker_cv_vs_host_selection_spring_entry_level.csv"),
  check.names = FALSE
) %>%
  as_tibble() %>%
  filter(Organism %in% c("Bacteria", "Fungi"), Trial == "SpringPan") %>%
  transmute(
    Organism,
    Trial,
    Entry,
    host_selection_local_norm,
    aggregate_cv = cv_abund
  )

biomarkers <- read.csv(
  file.path(project_root, "data", "processed", "network", "mgwas_sum_node.csv"),
  check.names = FALSE
) %>%
  as_tibble() %>%
  filter(group %in% c("Bacteria", "Fungi"))

meta <- read.csv(
  file.path(project_root, "data", "processed", "bulk", "meta_file.csv"),
  check.names = FALSE
) %>%
  as_tibble() %>%
  filter(Location == "Austria", Trial == "SpringPan") %>%
  select(sampleid, Entry, Trial)

bac_spring <- readRDS(
  file.path(project_root, "data", "processed", "bacteria", "asv_spring.rds")
)
fun_spring <- readRDS(
  file.path(project_root, "data", "processed", "fungi", "otu_spring.rds")
)

bac_spring_ra <- decostand(t(bac_spring), "total") * 100
fun_spring_ra <- decostand(t(fun_spring), "total") * 100

get_entry_level_biomarker_cv <- function(ra_mat, biomarker_df, meta_df, organism_label) {
  ids <- intersect(biomarker_df$Id, colnames(ra_mat))
  if (length(ids) == 0) return(tibble())

  as.data.frame(ra_mat[, ids, drop = FALSE]) %>%
    rownames_to_column("sampleid") %>%
    pivot_longer(-sampleid, names_to = "Biomarker", values_to = "Abundance") %>%
    left_join(meta_df, by = "sampleid") %>%
    group_by(Organism = organism_label, Trial, Entry, Biomarker) %>%
    summarise(
      mean_abundance = mean(Abundance, na.rm = TRUE),
      sd_abundance = sd(Abundance, na.rm = TRUE),
      cv_abundance = ifelse(mean_abundance == 0, NA_real_, sd_abundance / mean_abundance),
      n = sum(!is.na(Abundance)),
      .groups = "drop"
    )
}

indiv_cv <- bind_rows(
  get_entry_level_biomarker_cv(
    bac_spring_ra,
    biomarkers %>% filter(group == "Bacteria", Trial == "SpringPan", lme_slope == "Positive"),
    meta,
    "Bacteria"
  ),
  get_entry_level_biomarker_cv(
    fun_spring_ra,
    biomarkers %>% filter(group == "Fungi", Trial == "SpringPan", lme_slope == "Positive"),
    meta,
    "Fungi"
  )
) %>%
  inner_join(
    agg_dat %>% select(Organism, Trial, Entry, host_selection_local_norm),
    by = c("Organism", "Trial", "Entry")
  ) %>%
  filter(is.finite(cv_abundance), is.finite(host_selection_local_norm))

# Pick a readable subset: strongest biomarker-level CV~host-selection associations per panel
indiv_slopes <- indiv_cv %>%
  group_by(Organism, Trial, Biomarker) %>%
  group_modify(~ {
    d <- .x
    if (nrow(d) < 10) {
      return(tibble(slope = NA_real_, p_value = NA_real_, rho = NA_real_))
    }
    fit <- lm(cv_abundance ~ host_selection_local_norm, data = d)
    sp <- suppressWarnings(cor.test(d$host_selection_local_norm, d$cv_abundance, method = "spearman"))
    tibble(
      slope = unname(coef(fit)[["host_selection_local_norm"]]),
      p_value = summary(fit)$coefficients["host_selection_local_norm", "Pr(>|t|)"],
      rho = unname(sp$estimate)
    )
  }) %>%
  ungroup()

top_biomarkers <- indiv_slopes %>%
  filter(is.finite(p_value)) %>%
  group_by(Organism, Trial) %>%
  arrange(p_value, desc(abs(rho)), .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  mutate(
    biomarker_label = paste0(
      Biomarker,
      " (rho=",
      sprintf("%.2f", rho),
      ")"
    )
  )

plot_indiv <- indiv_cv %>%
  semi_join(top_biomarkers, by = c("Organism", "Trial", "Biomarker")) %>%
  left_join(
    top_biomarkers %>% select(Organism, Trial, Biomarker, biomarker_label),
    by = c("Organism", "Trial", "Biomarker")
  ) %>%
  mutate(
    panel = factor(
      paste(Organism, Trial, sep = " | "),
      levels = c("Bacteria | SpringPan", "Fungi | SpringPan")
    )
  )

agg_plot <- agg_dat %>%
  mutate(
    panel = factor(
      paste(Organism, Trial, sep = " | "),
      levels = c("Bacteria | SpringPan", "Fungi | SpringPan")
    )
  )

agg_stats <- agg_plot %>%
  group_by(panel) %>%
  group_modify(~ {
    fit <- lm(aggregate_cv ~ host_selection_local_norm, data = .x)
    sp <- suppressWarnings(cor.test(.x$host_selection_local_norm, .x$aggregate_cv, method = "spearman"))
    tibble(
      slope = unname(coef(fit)[["host_selection_local_norm"]]),
      p_value = summary(fit)$coefficients["host_selection_local_norm", "Pr(>|t|)"],
      rho = unname(sp$estimate)
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "Aggregate: rho = ", sprintf("%.2f", rho),
      "\nP = ", format.pval(p_value, digits = 2, eps = 1e-3)
    )
  )

write.csv(top_biomarkers, out_csv, row.names = FALSE)

p <- ggplot() +
  geom_point(
    data = plot_indiv,
    aes(
      x = host_selection_local_norm,
      y = cv_abundance,
      colour = biomarker_label
    ),
    alpha = 0.65,
    size = 1.9
  ) +
  geom_smooth(
    data = plot_indiv,
    aes(
      x = host_selection_local_norm,
      y = cv_abundance,
      colour = biomarker_label
    ),
    method = "lm",
    se = FALSE,
    linewidth = 0.9
  ) +
  geom_point(
    data = agg_plot,
    aes(
      x = host_selection_local_norm,
      y = aggregate_cv
    ),
    colour = "black",
    fill = "black",
    alpha = 0.85,
    size = 2.4,
    shape = 16
  ) +
  geom_smooth(
    data = agg_plot,
    aes(
      x = host_selection_local_norm,
      y = aggregate_cv
    ),
    method = "lm",
    se = TRUE,
    linewidth = 1.2,
    colour = "black",
    fill = "grey70"
  ) +
  facet_wrap(~ panel, scales = "free") +
  geom_text(
    data = agg_stats,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.04,
    vjust = 1.2,
    size = 3.8
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Host selection acts more coherently on positive biomarker pools than on individual biomarkers",
    subtitle = "Coloured points and lines show single positive biomarkers; black points and line show the aggregate pool",
    x = "Mean local host selection strength\n(root-to-local-bulk distance / local bulk baseline)",
    y = "CV of relative abundance",
    colour = "Positive biomarker"
  )

ggsave(out_png, p, width = 10.5, height = 5.2, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 10.5, height = 5.2, bg = "white")

write.csv(indiv_slopes, file.path(module_output, "all_individual_tests.csv"), row.names=FALSE)
write.csv(agg_stats, file.path(module_output, "aggregate_tests.csv"), row.names=FALSE)
