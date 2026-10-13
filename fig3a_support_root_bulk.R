# Figures: Fig3a
# Original source: BMB_bulk/comm_cor_root_bulk.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

rm(list=ls())


library(qiime2R)
library(tidyverse)
library(vegan)


bac <- read_q2biom("AUT_Bac/bulk_merged_table/feature-table.biom") # sample as col, otu as row
colnames(bac) <- sub("^FL_", "FL", colnames(bac))

fun <- readRDS("AUT_Fun/fun_otu_count_with_bulk.rds") |> t() |> as.data.frame()# make format consistant
colnames(fun) <- sub("^FL_", "FL", colnames(fun))
meta <- read.csv("../Bac_Aut/0_raw_data/meta_file.csv")


meta_aut <- meta |> 
  filter(Location=="Austria",Trial!="DivPan")


meta4bulk_bac <- data.frame(sampleid = colnames(bac)) |>
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
    ),
    Entry = "Bulk_Soil"
  ) |>
  left_join(
    meta_aut |>
      select(Trial, Plot, Col, Row) |>
      distinct(),
    by = c("Trial", "Plot")
  ) |>
  mutate(Side = NA_character_) |>
  select(sampleid, Trial, Entry, Plot, Col, Row, Side)

meta4all_bac <- rbind(meta_aut |> select(sampleid,Trial,Entry,Plot,Col,Row,Side), meta4bulk_bac)


meta4bulk_bac_spring <- meta4bulk_bac |> filter(Trial=="SpringPan")
meta4bulk_bac_winter <- meta4bulk_bac |> filter(Trial=="WinterPan")

# bac spring
bac_bulk_spring <- bac[,meta4bulk_bac_spring$sampleid]
bac_bulk_spring <- bac_bulk_spring[rowSums(bac_bulk_spring)>2,]
dim(bac_bulk_spring)
#bac winter
bac_bulk_winter <- bac[,meta4bulk_bac_winter$sampleid]
bac_bulk_winter <- bac_bulk_winter[rowSums(bac_bulk_winter)>2,]
dim(bac_bulk_winter)

all(meta4bulk_bac_spring$sampleid==colnames(bac_bulk_spring))

eu_spatial_bulk_spring <- vegdist(meta4bulk_bac_spring |> select(Row,Col),"eu")
eu_spatial_bulk_winter<- vegdist(meta4bulk_bac_winter |> select(Row,Col),"eu")

bc_bac_bulk_spring <- bac_bulk_spring |> t() |> 
  decostand("hellinger") |> 
  vegdist("bray")


bc_bac_bulk_winter <- bac_bulk_winter |> t() |> 
  decostand("hellinger") |> 
  vegdist("bray")


mantel(eu_spatial_bulk_spring ,bc_bac_bulk_spring,method="spearman")
mantel(eu_spatial_bulk_winter ,bc_bac_bulk_winter,method="spearman")


meta4bulk_fun <- data.frame(sampleid = colnames(fun)) |>
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
    ),
    Entry = "Bulk_Soil"
  ) |>
  left_join(
    meta_aut |>
      select(Trial, Plot, Col, Row) |>
      distinct(),
    by = c("Trial", "Plot")
  ) |>
  mutate(Side = NA_character_) |>
  select(sampleid, Trial, Entry, Plot, Col, Row, Side)


meta4bulk_fun_spring <- meta4bulk_fun |> filter(Trial=="SpringPan")
meta4bulk_fun_winter <- meta4bulk_fun |> filter(Trial=="WinterPan")
# fun spring
fun_bulk_spring <- fun |> select(meta4bulk_fun_spring$sampleid)
fun_bulk_spring <- fun_bulk_spring[rowSums(fun_bulk_spring)>2,]
dim(fun_bulk_spring)
#fun winter
fun_bulk_winter <- fun[,meta4bulk_fun_winter$sampleid]
fun_bulk_winter <- fun_bulk_winter[rowSums(fun_bulk_winter)>2,]
dim(fun_bulk_winter)

all(meta4bulk_fun_spring$sampleid==colnames(fun_bulk_spring))

eu_spatial_bulk_spring <- vegdist(meta4bulk_fun_spring |> select(Row,Col),"eu")
eu_spatial_bulk_winter<- vegdist(meta4bulk_fun_winter |> select(Row,Col),"eu")

bc_fun_bulk_spring <- fun_bulk_spring |> t() |> 
  decostand("hellinger") |> 
  vegdist("bray")


bc_fun_bulk_winter <- fun_bulk_winter |> t() |> 
  decostand("hellinger") |> 
  vegdist("bray")


mantel(eu_spatial_bulk_spring ,bc_fun_bulk_spring,method="spearman")
mantel(eu_spatial_bulk_winter ,bc_fun_bulk_winter,method="spearman")


# =====================================================
# Generic pipeline: bulk-root Mantel at plot level
# =====================================================

run_bulk_root_mantel <- function(
    otu_root,
    otu_bulk,
    meta_root,
    meta_bulk,
    trial_name,
    root_filter_threshold = 10,
    permutations = 9999
) {
  
  # -----------------------------
  # Root metadata: keep samples in selected trial
  # -----------------------------
  meta_root2 <- meta_root |>
    filter(
      Trial == trial_name,
      Plot %in% meta_bulk$Plot,
      sampleid %in% colnames(otu_root)
    )
  
  # subset root OTU table
  otu_root2 <- otu_root[, colnames(otu_root) %in% meta_root2$sampleid, drop = FALSE]
  
  # reorder metadata to match root OTU table
  meta_root2 <- meta_root2[
    match(colnames(otu_root2), meta_root2$sampleid),
  ]
  
  stopifnot(all(colnames(otu_root2) == meta_root2$sampleid))
  
  # -----------------------------
  # Bulk metadata: keep samples present in bulk OTU table
  # -----------------------------
  meta_bulk2 <- meta_bulk |>
    filter(sampleid %in% colnames(otu_bulk))
  
  meta_bulk2 <- meta_bulk2[
    match(colnames(otu_bulk), meta_bulk2$sampleid),
  ]
  
  keep_bulk <- !is.na(meta_bulk2$sampleid)
  otu_bulk2 <- otu_bulk[, keep_bulk, drop = FALSE]
  meta_bulk2 <- meta_bulk2[keep_bulk, ]
  
  stopifnot(all(colnames(otu_bulk2) == meta_bulk2$sampleid))
  
  # -----------------------------
  # Keep only common plots
  # -----------------------------
  common_plot <- intersect(unique(meta_root2$Plot), unique(meta_bulk2$Plot))
  
  meta_root2 <- meta_root2 |>
    filter(Plot %in% common_plot)
  
  meta_bulk2 <- meta_bulk2 |>
    filter(Plot %in% common_plot)
  
  otu_root2 <- otu_root2[, meta_root2$sampleid, drop = FALSE]
  otu_bulk2 <- otu_bulk2[, meta_bulk2$sampleid, drop = FALSE]
  
  # reorder again after filtering
  meta_root2 <- meta_root2[
    match(colnames(otu_root2), meta_root2$sampleid),
  ]
  
  meta_bulk2 <- meta_bulk2[
    match(colnames(otu_bulk2), meta_bulk2$sampleid),
  ]
  
  # -----------------------------
  # Aggregate root to plot level
  # if 2 samples -> mean of 2
  # if 1 sample  -> itself
  # -----------------------------
  plot_order <- meta_bulk2$Plot
  
  otu_root_plot_mean <- sapply(plot_order, function(p) {
    cols <- which(meta_root2$Plot == p)
    rowMeans(otu_root2[, cols, drop = FALSE])
  })
  
  colnames(otu_root_plot_mean) <- plot_order
  
  # rename bulk columns to Plot
  colnames(otu_bulk2) <- meta_bulk2$Plot
  
  # make sure root and bulk are aligned
  stopifnot(
    identical(colnames(otu_root_plot_mean), colnames(otu_bulk2))
  )
  
  # -----------------------------
  # Bray-Curtis distance
  # -----------------------------
  bc_root <- otu_root_plot_mean[rowSums(otu_root_plot_mean) > root_filter_threshold, , drop = FALSE] |>
    t() |>
    decostand("hellinger") |>
    vegdist(method = "bray")
  
  bc_bulk <- otu_bulk2 |>
    t() |>
    decostand("hellinger") |>
    vegdist(method = "bray")
  
  stopifnot(
    identical(attr(bc_root, "Labels"), attr(bc_bulk, "Labels"))
  )
  
  mantel_bc <- mantel(
    bc_bulk,
    bc_root,
    permutations = permutations,
    method = "spearman"
  )
  
  # -----------------------------
  # Jaccard distance
  # -----------------------------
  jaccard_root <- otu_root_plot_mean[rowSums(otu_root_plot_mean) > root_filter_threshold, , drop = FALSE] |>
    t() |>
    vegdist(method = "jaccard", binary = TRUE)
  
  jaccard_bulk <- otu_bulk2 |>
    t() |>
    vegdist(method = "jaccard", binary = TRUE)
  
  stopifnot(
    identical(attr(jaccard_root, "Labels"), attr(jaccard_bulk, "Labels"))
  )
  
  mantel_jaccard <- mantel(
    jaccard_bulk,
    jaccard_root,
    permutations = permutations,
    method = "spearman"
  )
  
  # -----------------------------
  # Return everything
  # -----------------------------
  return(list(
    meta_root = meta_root2,
    meta_bulk = meta_bulk2,
    otu_root_plot_mean = otu_root_plot_mean,
    otu_bulk_plot = otu_bulk2,
    bc_root = bc_root,
    bc_bulk = bc_bulk,
    jaccard_root = jaccard_root,
    jaccard_bulk = jaccard_bulk,
    mantel_bc = mantel_bc,
    mantel_jaccard = mantel_jaccard
  ))
}


# bacteria spring
res_bac_spring <- run_bulk_root_mantel(
  otu_root = bac,
  otu_bulk = bac_bulk_spring,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_bac_spring,
  trial_name = "SpringPan"
)

# bacteria winter
res_bac_winter <- run_bulk_root_mantel(
  otu_root = bac,
  otu_bulk = bac_bulk_winter,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_bac_winter,
  trial_name = "WinterPan"
)

# fungi spring
res_fun_spring <- run_bulk_root_mantel(
  otu_root = fun,
  otu_bulk = fun_bulk_spring,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_fun_spring,
  trial_name = "SpringPan"
)

# fungi winter
res_fun_winter <- run_bulk_root_mantel(
  otu_root = fun,
  otu_bulk = fun_bulk_winter,
  meta_root = meta_aut,
  meta_bulk = meta4bulk_fun_winter,
  trial_name = "WinterPan"
)


res_bac_spring$mantel_bc


res_bac_spring$mantel_jaccard


res_bac_winter$mantel_bc


res_bac_winter$mantel_jaccard


res_fun_spring$mantel_bc


res_fun_spring$mantel_jaccard


res_fun_winter$mantel_bc


res_fun_winter$mantel_jaccard
