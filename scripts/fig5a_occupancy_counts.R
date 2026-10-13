# Prepared for reviewer use; original analysis retained in originals/.
if (!exists("project_root", inherits=FALSE)) {
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash="/", mustWork=TRUE)) else getwd()
project_root <- normalizePath(file.path(script_dir, ".."), winslash="/", mustWork=TRUE)
}
module_output <- file.path(project_root, "outputs", "fig5a")
dir.create(module_output, recursive=TRUE, showWarnings=FALSE)
setwd(module_output)
dir.create("Network/k_tables", recursive=TRUE, showWarnings=FALSE)
on_error <- NULL

library(tidyverse)

# ============================================================
# Inputs
# ============================================================
edge <- read.csv(file.path(project_root, "data/processed/network/mgwas_sum_edge.csv"), stringsAsFactors = FALSE)

asv <- readRDS(file.path(project_root, "data/processed/bacteria/asv_rarefied.RDS"))  # Sample x ASV
otu <- readRDS(file.path(project_root, "data/processed/fungi/otu_rarefied.RDS"))  # Sample x OTU

meta_spring_bac <- readRDS(file.path(project_root, "data/processed/bacteria/meta_spring.rds"))
meta_winter_bac <- readRDS(file.path(project_root, "data/processed/bacteria/meta_winter.rds"))
meta_spring_fun <- readRDS(file.path(project_root, "data/processed/fungi/meta_spring.rds"))
meta_winter_fun <- readRDS(file.path(project_root, "data/processed/fungi/meta_winter.rds"))

out_dir <- "Network/k_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Sanity checks
# ============================================================
stopifnot("sampleid" %in% colnames(meta_spring_bac))
stopifnot("Entry"    %in% colnames(meta_spring_bac))

# Ensure character (avoid factor issues)
meta_spring_bac$sampleid <- as.character(meta_spring_bac$sampleid)
meta_winter_bac$sampleid <- as.character(meta_winter_bac$sampleid)
meta_spring_fun$sampleid <- as.character(meta_spring_fun$sampleid)
meta_winter_fun$sampleid <- as.character(meta_winter_fun$sampleid)

# ============================================================
# Extract targets
# ============================================================
bac_spring_targets <- edge |> filter(Trial=="SpringPan", grepl("^ASV", Target)) |> distinct(Target) |> pull(Target)
bac_winter_targets <- edge |> filter(Trial=="WinterPan", grepl("^ASV", Target)) |> distinct(Target) |> pull(Target)

fun_spring_targets <- edge |> filter(Trial=="SpringPan", grepl("^OTU", Target)) |> distinct(Target) |> pull(Target)
fun_winter_targets <- edge |> filter(Trial=="WinterPan", grepl("^OTU", Target)) |> distinct(Target) |> pull(Target)

# ============================================================
# Helper: compute k=1..6 Entry counts (no renaming, no transpose)
# ============================================================
calc_k <- function(count_mat, meta_df, targets, label="") {
  
  cat("\n---", label, "---\n")
  
  # Which targets actually exist
  targets_found <- intersect(targets, colnames(count_mat))
  cat("Targets total:", length(targets), "\n")
  cat("Targets found:", length(targets_found), "\n")
  if (length(targets_found) == 0) return(NULL)
  
  # Sample overlap
  sample_keep <- intersect(rownames(count_mat), meta_df$sampleid)
  cat("Samples in count table:", nrow(count_mat), "\n")
  cat("Samples in meta:", nrow(meta_df), "\n")
  cat("Samples overlapped:", length(sample_keep), "\n")
  
  # Presence / absence
  pa <- (count_mat[sample_keep, targets_found, drop = FALSE] > 0)
  
  # k per Entry
  k_by_entry <- pa |>
    as.data.frame() |>
    rownames_to_column("sampleid") |>
    pivot_longer(-sampleid, names_to="Target", values_to="present") |>
    left_join(meta_df |> select(sampleid, Entry), by="sampleid") |>
    filter(!is.na(Entry)) |>
    group_by(Target, Entry) |>
    summarise(k = sum(present), .groups="drop")
  
  # Count Entries for k=1..6
  k_counts <- k_by_entry |>
    count(Target, k, name="n_entry") |>
    filter(k %in% 1:6) |>
    tidyr::complete(Target, k=1:6, fill=list(n_entry=0)) |>
    arrange(Target, k)
  
  list(k_by_entry = k_by_entry, k_counts = k_counts)
}

# ============================================================
# Run
# ============================================================
res_bac_spring <- calc_k(asv, meta_spring_bac, bac_spring_targets, "Spring bacteria")
res_bac_winter <- calc_k(asv, meta_winter_bac, bac_winter_targets, "Winter bacteria")

res_fun_spring <- calc_k(otu, meta_spring_fun, fun_spring_targets, "Spring fungi")
res_fun_winter <- calc_k(otu, meta_winter_fun, fun_winter_targets, "Winter fungi")

# ============================================================
# Save
# ============================================================
if (!is.null(res_bac_spring)) {
  write.csv(res_bac_spring$k_counts, file.path(out_dir,"bac_spring_k1to6_entry_counts.csv"), row.names=FALSE)
  write.csv(res_bac_spring$k_by_entry, file.path(out_dir,"bac_spring_k_by_entry.csv"), row.names=FALSE)
}

if (!is.null(res_bac_winter)) {
  write.csv(res_bac_winter$k_counts, file.path(out_dir,"bac_winter_k1to6_entry_counts.csv"), row.names=FALSE)
  write.csv(res_bac_winter$k_by_entry, file.path(out_dir,"bac_winter_k_by_entry.csv"), row.names=FALSE)
}

if (!is.null(res_fun_spring)) {
  write.csv(res_fun_spring$k_counts, file.path(out_dir,"fun_spring_k1to6_entry_counts.csv"), row.names=FALSE)
  write.csv(res_fun_spring$k_by_entry, file.path(out_dir,"fun_spring_k_by_entry.csv"), row.names=FALSE)
}

if (!is.null(res_fun_winter)) {
  write.csv(res_fun_winter$k_counts, file.path(out_dir,"fun_winter_k1to6_entry_counts.csv"), row.names=FALSE)
  write.csv(res_fun_winter$k_by_entry, file.path(out_dir,"fun_winter_k_by_entry.csv"), row.names=FALSE)
}

# ============================================================
# Inspect
# ============================================================
cat("\n=== Spring bacteria k counts ===\n")
print(res_bac_spring$k_counts)
