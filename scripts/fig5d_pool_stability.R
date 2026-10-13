# Figures: Fig5
# Original source: 0_Submission/4_NatComm/Manuscript_Figure_Code/scripts/fig5d_pool_stability.R
# Analysis statements retained; input/output locations are handled by run_figure.R.

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


infile <- file.path(project_root, "data", "processed", "derived", "biomarker_pool_stability_entry_compare.csv")
out_dir <- module_output
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_png <- file.path(out_dir, "positive_biomarker_pool_cv_ratio_spring.png")
out_pdf <- file.path(out_dir, "positive_biomarker_pool_cv_ratio_spring.pdf")

dat <- read.csv(infile, check.names = FALSE) %>%
  as_tibble() %>%
  filter(
    pool == "Positive",
    paste(Organism, Trial, sep = " | ") %in% c(
      "Bacteria | SpringPan",
      "Fungi | SpringPan"
    ),
    is.finite(cv_ratio)
  ) %>%
  mutate(
    panel = factor(
      paste(Organism, Trial, sep = " | "),
      levels = c("Bacteria | SpringPan", "Fungi | SpringPan")
    )
  )

p <- ggplot(dat, aes(x = panel, y = cv_ratio, color = panel)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.6) +
  geom_jitter(width = 0.14, height = 0, alpha = 0.7, size = 1.6) +
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.42,
    linewidth = 0.5,
    color = "black"
  ) +
  scale_color_manual(
    values = c(
      "Bacteria | SpringPan" = "#F08A7E",
      "Fungi | SpringPan" = "#1FBFCF"
    )
  ) +
  coord_cartesian(ylim = c(0, max(1.6, quantile(dat$cv_ratio, 0.98, na.rm = TRUE)))) +
  labs(
    title = "Relative stability of positive biomarker pools",
    subtitle = "Each point is one genotype; values below 1 indicate lower CV for the aggregate pool",
    x = NULL,
    y = "Aggregate CV / median individual CV",
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 15, hjust = 1)
  )

ggsave(out_png, p, width = 6.2, height = 4.8, dpi = 300, bg = "white")
ggsave(out_pdf, p, width = 6.2, height = 4.8, bg = "white")
