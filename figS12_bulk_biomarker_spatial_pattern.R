# Figures: FigS12
# Original source: BMB_bulk/spatial_pattern_SNP_biomarkers.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

project_dir <- "C:/Users/au747146/OneDrive - Aarhus universitet/Desktop/BMB_AUT_MS"
knitr::opts_knit$set(root.dir = project_dir)
setwd(project_dir)


library(tidyverse)
library(grid)


node <- read.csv("Network/mgwas_sum_node.csv", check.names = FALSE)

biomarkers <- node |>
  filter(group %in% c("Bacteria", "Fungi")) |>
  select(Id, group, lme_slope, Trial) |>
  rename(
    Biomarker = Id,
    Organism  = group,
    Direction = lme_slope,
    Panel     = Trial
  ) |>
  filter(Direction == "Positive")


library(qiime2R)

# Raw community matrices (OTU/ASV x sample, both contain bulk samples)
bac_raw <- read_q2biom("BMB_bulk/AUT_Bac/bulk_merged_table/feature-table.biom")
colnames(bac_raw) <- sub("^FL_", "FL", colnames(bac_raw))

fun_raw <- readRDS("BMB_bulk/AUT_Fun/fun_otu_count_with_bulk.rds") |>
  t() |> as.data.frame()
colnames(fun_raw) <- sub("^FL_", "FL", colnames(fun_raw))

# Taxonomy: maps QIIME2 hash FeatureIDs -> ASV IDs
tax <- readRDS("Bac_Aut/0_raw_data/taxonomy.RDS")

# Full metadata with spatial coordinates
meta <- read.csv("Bac_Aut/0_raw_data/meta_file.csv")
meta_aut <- meta |> filter(Location == "Austria", Trial != "DivPan")


plot_coord <- meta_aut |>
  filter(Trial != "DivPan") |>
  mutate(
    Plot = as.numeric(Plot),
    Col  = as.numeric(Col),
    Row  = as.numeric(Row)
  ) |>
  group_by(Trial, Plot) |>
  summarise(
    Col = first(na.omit(Col)),
    Row = first(na.omit(Row)),
    n_coord = n_distinct(paste(Col, Row)),
    .groups = "drop"
  )

plot_coord |> filter(n_coord > 1)


make_bulk_meta <- function(sample_names, plot_coord) {
  tibble(sampleid = sample_names) |>
    filter(str_detect(sampleid, "^bulk_?\\d+$")) |>
    mutate(
      bulk_num = as.numeric(str_extract(sampleid, "\\d+$")),
      Trial = case_when(
        bulk_num >= 1    & bulk_num <= 600  ~ "WinterPan",
        bulk_num >= 601  & bulk_num <= 1200 ~ "DivPan",
        bulk_num >= 1201 & bulk_num <= 1800 ~ "SpringPan",
        TRUE ~ NA_character_
      ),
      Plot = case_when(
        bulk_num >= 1    & bulk_num <= 600  ~ bulk_num,
        bulk_num >= 601  & bulk_num <= 1200 ~ bulk_num - 600,
        bulk_num >= 1201 & bulk_num <= 1800 ~ bulk_num - 1200,
        TRUE ~ NA_real_
      )
    ) |>
    filter(Trial != "DivPan") |>
    left_join(plot_coord, by = c("Trial", "Plot")) |>
    filter(!is.na(Col), !is.na(Row)) |>
    filter((Col - Row) %% 20 %in% c(0, 10))
}


meta_bulk_bac <- make_bulk_meta(colnames(bac_raw), plot_coord)
meta_bulk_fun <- make_bulk_meta(colnames(fun_raw), plot_coord)

bind_rows(
  meta_bulk_bac |> mutate(dataset = "Bacteria"),
  meta_bulk_fun |> mutate(dataset = "Fungi")
) |>
  filter(Col < 1 | Col > 20 | Row < 1 | Row > 40 | is.na(Col) | is.na(Row)) |>
  arrange(dataset, Trial, Plot)


bac_biomarkers <- biomarkers |> filter(Organism == "Bacteria") |> pull(Biomarker)
fun_biomarkers <- biomarkers |> filter(Organism == "Fungi")   |> pull(Biomarker)

# Bacteria: map hash FeatureIDs to ASV IDs, then extract bulk samples
bac_id_map <- tax |>
  filter(ASVid %in% bac_biomarkers) |>
  select(FeatureID, ASVid)

bac_bulk_counts <- as.matrix(bac_raw[bac_id_map$FeatureID, meta_bulk_bac$sampleid])
bac_col_totals  <- Matrix::colSums(bac_raw[, meta_bulk_bac$sampleid])
bac_bulk_ra     <- sweep(bac_bulk_counts, 2, bac_col_totals, "/")
rownames(bac_bulk_ra) <- bac_id_map$ASVid[match(rownames(bac_bulk_ra), bac_id_map$FeatureID)]

# Fungi: OTUs already named, extract bulk samples
fun_all_bulk    <- as.matrix(fun_raw[, meta_bulk_fun$sampleid])
fun_col_totals  <- colSums(fun_all_bulk)
fun_bulk_ra     <- sweep(fun_all_bulk[fun_biomarkers, , drop = FALSE], 2, fun_col_totals, "/")

cat("Bac biomarker RA matrix:", dim(bac_bulk_ra), "\n")
cat("Fun biomarker RA matrix:", dim(fun_bulk_ra), "\n")


to_long <- function(ra_mat, meta_bulk) {
  as.data.frame(ra_mat) |>
    rownames_to_column("Biomarker") |>
    pivot_longer(-Biomarker, names_to = "sampleid", values_to = "RA") |>
    left_join(meta_bulk |> select(sampleid, Trial, Col, Row), by = "sampleid") |>
    group_by(Biomarker, Trial, Col, Row) |>
    summarise(mean_RA = mean(RA, na.rm = TRUE), .groups = "drop")
}

plot_df <- bind_rows(
  to_long(bac_bulk_ra, meta_bulk_bac),
  to_long(fun_bulk_ra, meta_bulk_fun)
) |>
  left_join(biomarkers, by = "Biomarker")

max_ra <- max(plot_df$mean_RA, na.rm = TRUE)
if (!is.finite(max_ra)) max_ra <- 1
ra_limits <- c(0, max_ra)

head(plot_df)


plot_field_heatmap <- function(df, biomarker_id) {
  meta_bm <- biomarkers |> filter(Biomarker == biomarker_id)
  df_sub  <- df |> filter(Biomarker == biomarker_id, Trial == meta_bm$Panel)

  # Absent plots (RA = 0) get NA fill so they render as a distinct grey
  df_sub <- df_sub |> mutate(fill_val = ifelse(mean_RA > 0, mean_RA, NA_real_))

  ggplot(df_sub, aes(x = Col, y = Row)) +
    geom_tile(aes(fill = fill_val), color = "grey40", linewidth = 0.3) +
    geom_text(
      data = df_sub |> filter(mean_RA == 0),
      aes(label = "x"),
      color = "grey50", size = 2.5
    ) +
    scale_fill_viridis_c(
      name = "RA",
      option = "magma",
      na.value = "grey92",
      limits = ra_limits
    ) +
    scale_x_continuous(breaks = seq(2, 20, by = 2)) +
    scale_y_continuous(breaks = seq(5, 30, by = 5)) +
    coord_fixed() +
    labs(
      title    = biomarker_id,
      subtitle = paste0(meta_bm$Organism, "  |  ", meta_bm$Panel,
                        "  |  Direction: ", meta_bm$Direction),
      x = "Col", y = "Row"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank())
}

# Generate one plot per biomarker
heatmap_list <- setNames(
  lapply(biomarkers$Biomarker, plot_field_heatmap, df = plot_df),
  biomarkers$Biomarker
)


library(patchwork)

# Compact plot for use inside combined figures
plot_field_heatmap_compact <- function(df, biomarker_id) {
  meta_bm <- biomarkers |> filter(Biomarker == biomarker_id)
  df_sub  <- df |> filter(Biomarker == biomarker_id, Trial == meta_bm$Panel) |>
    mutate(fill_val = ifelse(mean_RA > 0, mean_RA, NA_real_))

  dir_symbol <- "(+)"
  title_col  <- "#D55E00"

  ggplot(df_sub, aes(x = Col, y = Row)) +
    geom_tile(aes(fill = fill_val), color = "grey40", linewidth = 0.2) +
    geom_text(
      data = df_sub |> filter(mean_RA == 0),
      aes(label = "x"), color = "grey60", size = 1.8
    ) +
    scale_fill_viridis_c(
      name = "RA", option = "magma", na.value = "grey92",
      limits = ra_limits
    ) +
    scale_x_continuous(breaks = seq(4, 20, by = 4)) +
    scale_y_continuous(breaks = seq(10, 30, by = 10)) +
    labs(title = paste(biomarker_id, dir_symbol), x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid        = element_blank(),
      plot.title        = element_text(size = 8, face = "bold", hjust = 0.5,
                                       color = title_col),
      legend.key.height = unit(0.4, "cm"),
      legend.title      = element_text(size = 7),
      legend.text       = element_text(size = 6)
    )
}

make_combined <- function(organism, panel, ncol = NULL, collect_guides = TRUE) {
  bm     <- biomarkers |> filter(Organism == organism, Panel == panel)
  bm_pos <- bm |> pull(Biomarker)

  if (is.null(ncol)) ncol <- min(5, nrow(bm))

  make_section <- function(ids, label) {
    plots <- lapply(ids, plot_field_heatmap_compact, df = plot_df)
    section <- wrap_plots(plots, ncol = ncol) +
      plot_annotation(
        subtitle = label,
        theme = theme(plot.subtitle = element_text(size = 10, face = "italic",
                                                   color = "grey30"))
      )

    if (collect_guides) {
      section <- section + plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
    }

    section
  }

  combined <- if (length(bm_pos) > 0) make_section(bm_pos, "Positive") else plot_spacer()

  combined + plot_annotation(
    title = paste(organism, "|", panel),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )
}


combined_all <- wrap_plots(
  make_combined("Bacteria", "SpringPan", collect_guides = FALSE),
  make_combined("Fungi", "SpringPan", ncol = 5, collect_guides = FALSE),
  make_combined("Fungi", "WinterPan", collect_guides = FALSE),
  ncol = 1,
  guides = "collect"
) &
  theme(legend.position = "bottom")

combined_all

ggsave(
  filename = "BMB_bulk/spatial_pattern_positive_biomarkers.png",
  plot = combined_all,
  width = 8,
  height = 8,
  dpi = 300
)
ggsave(
  filename = "BMB_bulk/spatial_pattern_positive_biomarkers.pdf",
  plot = combined_all,
  width = 8,
  height = 8,
)
