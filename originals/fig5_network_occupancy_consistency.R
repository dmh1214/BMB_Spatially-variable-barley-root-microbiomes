#rm(list=ls())
library(tidyverse)

# Read k-count tables
bac_spring <- read.csv("Network/k_tables/bac_spring_k1to6_entry_counts.csv")
bac_winter <- read.csv("Network/k_tables/bac_winter_k1to6_entry_counts.csv")
fun_spring <- read.csv("Network/k_tables/fun_spring_k1to6_entry_counts.csv")
fun_winter <- read.csv("Network/k_tables/fun_winter_k1to6_entry_counts.csv")

# ============================================================
# Inputs
# ============================================================
edge <- read.csv("Network/mgwas_sum_edge.csv", stringsAsFactors = FALSE)

asv <- readRDS("Bac_Aut/1_processed_data/asv_rarefied.RDS")  # Sample x ASV
otu <- readRDS("Fun_Aut/1_processed_data/otu_rarefied.RDS")  # Sample x OTU


out_dir <- "Network/k_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Sanity checks
# ============================================================
stopifnot("sampleid" %in% colnames(meta_bac_spring))
stopifnot("Entry"    %in% colnames(meta_bac_spring))

# Ensure character (avoid factor issues)
meta_bac_spring$sampleid <- as.character(meta_bac_spring$sampleid)
meta_bac_winter$sampleid <- as.character(meta_bac_winter$sampleid)
meta_fun_spring$sampleid <- as.character(meta_fun_spring$sampleid)
meta_fun_winter$sampleid <- as.character(meta_fun_winter$sampleid)

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
calc_k <- function(count_mat, meta_df, targets, label = "") {
  
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
  
  # Restrict meta to overlapped samples and valid Entry
  meta_keep <- meta_df |>
    dplyr::filter(sampleid %in% sample_keep) |>
    dplyr::select(sampleid, Entry) |>
    dplyr::filter(!is.na(Entry))
  
  # How many replicates actually exist per Entry (after overlap)
  entry_nrep <- meta_keep |>
    dplyr::count(Entry, name = "n_rep_entry")
  
  # Presence / absence
  pa <- (count_mat[sample_keep, targets_found, drop = FALSE] > 0)
  
  # k per Entry + attach n_rep_entry
  k_by_entry <- pa |>
    as.data.frame() |>
    tibble::rownames_to_column("sampleid") |>
    tidyr::pivot_longer(-sampleid, names_to = "Target", values_to = "present") |>
    dplyr::left_join(meta_keep, by = "sampleid") |>
    dplyr::group_by(Target, Entry) |>
    dplyr::summarise(k = sum(present), .groups = "drop") |>
    dplyr::left_join(entry_nrep, by = "Entry") |>
    dplyr::mutate(
      n_rep_entry = as.integer(n_rep_entry),
      k_prop = k / n_rep_entry
    )
  
  # Count Entries for k=1..6 (still works; missing-rep Entries simply can't contribute to high k)
  k_counts <- k_by_entry |>
    dplyr::count(Target, k, name = "n_entry") |>
    dplyr::filter(k %in% 1:6) |>
    tidyr::complete(Target, k = 1:6, fill = list(n_entry = 0)) |>
    dplyr::arrange(Target, k)
  
  # Optional diagnostics: how many Entries have <6 reps?
  diag_entry_rep <- entry_nrep |>
    dplyr::count(n_rep_entry, name = "n_entry") |>
    dplyr::arrange(n_rep_entry)
  
  list(
    k_by_entry = k_by_entry,
    k_counts = k_counts,
    entry_nrep = entry_nrep,
    diag_entry_rep = diag_entry_rep
  )
}


# ============================================================
# Run
# ============================================================
res_bac_spring <- calc_k(asv, meta_bac_spring, bac_spring_targets, "Spring bacteria")
res_bac_winter <- calc_k(asv, meta_bac_winter, bac_winter_targets, "Winter bacteria")

res_fun_spring <- calc_k(otu, meta_fun_spring, fun_spring_targets, "Spring fungi")
res_fun_winter <- calc_k(otu, meta_fun_winter, fun_winter_targets, "Winter fungi")

k_ge_counts <- function(k_by_entry, k_max = 6) {
  # For each Target, count number of Entries with k >= threshold
  map_dfr(1:k_max, function(th) {
    k_by_entry |>
      group_by(Target) |>
      summarise(n_entry = sum(k >= th), .groups = "drop") |>
      mutate(k = th)
  }) |>
    arrange(Target, k)
}

bac_spring_ge <- k_ge_counts(res_bac_spring$k_by_entry)
bac_winter_ge <- k_ge_counts(res_bac_winter$k_by_entry)
fun_spring_ge <- k_ge_counts(res_fun_spring$k_by_entry)
fun_winter_ge <- k_ge_counts(res_fun_winter$k_by_entry)

write.csv(bac_spring_ge, "Network/k_tables/bac_spring_kGE_entry_counts.csv", row.names = FALSE)
write.csv(bac_winter_ge, "Network/k_tables/bac_winter_kGE_entry_counts.csv", row.names = FALSE)
write.csv(fun_spring_ge, "Network/k_tables/fun_spring_kGE_entry_counts.csv", row.names = FALSE)
write.csv(fun_winter_ge, "Network/k_tables/fun_winter_kGE_entry_counts.csv", row.names = FALSE)


plot_k_ge_line <- function(df, title_text) {
  ggplot(df, aes(x = k, y = n_entry, color = Target, group = Target)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = 1:6) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(
      x = "k (minimum number of replicates detected, ≥k)",
      y = "Number of entries",
      title = title_text
    )
}

p_bac_spring <- plot_k_ge_line(bac_spring_ge, "SpringPan – Bacteria (ASVs)")
p_bac_winter <- plot_k_ge_line(bac_winter_ge, "WinterPan – Bacteria (ASVs)")
p_fun_spring <- plot_k_ge_line(fun_spring_ge, "SpringPan – Fungi (OTUs)")
p_fun_winter <- plot_k_ge_line(fun_winter_ge, "WinterPan – Fungi (OTUs)")

# bac spring ----
## stack barplot----
# =========================
# 1) Build exact-k count table from k_by_entry
# =========================
# exact-k counts per Target
# 1) exact-k counts per Target
df_exact <- res_bac_spring$k_by_entry |>
  count(Target, k, name = "n_entry") |>
  filter(k %in% 1:6) |>
  tidyr::complete(Target, k = 1:6, fill = list(n_entry = 0))

# 2) set your desired Target order (no reversing)
target_levels <- c("ASV16927","ASV137553","ASV68333","ASV26933","ASV42347","ASV25606")

df_plot <- df_exact |>
  mutate(
    Target = factor(Target, levels = rev(target_levels)),
    k = factor(k, levels = 1:6)
  )

# Optional: reverse display order (uncomment if you want the list upside-down)
# df_plot <- df_plot |> mutate(Target = fct_rev(Target))

# 3) normal horizontal stacked bar plot
p <- ggplot(df_plot, aes(x = n_entry, y = Target, fill = k)) +
  geom_col(width = 0.8) +
  scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_y_discrete(position = "right") +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    
    # show Target axis on the right only
    axis.line.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left  = element_blank(),
    axis.line.y.right  = element_line(),
    axis.ticks.y.right = element_line(),
    axis.text.y.right  = element_text(hjust = 0),
    
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Number of entries",
    y = "ASV",
    fill = "Exact k\n(# of replicates)",
    title = "Occupancy consistency of bacterial biomarkers (SpringPan)"
  )

p

ggsave("Network/Occupancy consistency bac_spring.pdf",p,height=6,width = 6)

### pieplot ----

# df_plot should already have:
# Target (factor, ordered), k (1–6), n_entry

# 1) Prepare proportions
df_pie <- df_plot |>
  group_by(Target) |>
  mutate(prop = n_entry / sum(n_entry)) |>
  ungroup()

# 2) Make long-row pie plot
p_pies <- ggplot(df_pie, aes(x = 1, y = prop, fill = k)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  facet_wrap(~ Target, nrow = 1) +   # <<< one long row
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

p_pies
ggsave("Network/Occupancy_pies_bac_spring.pdf",p_pies,height=2,width=10)

# bac winter ----
## stack barplot----
# =========================
# 1) Build exact-k count table from k_by_entry
# =========================
# exact-k counts per Target
# 1) exact-k counts per Target
df_exact <- res_bac_winter$k_by_entry |>
  count(Target, k, name = "n_entry") |>
  filter(k %in% 1:6) |>
  tidyr::complete(Target, k = 1:6, fill = list(n_entry = 0))

# 2) set your desired Target order (no reversing)


df_plot <- df_exact |>
  mutate(
    Target = factor(Target),
    k = factor(k, levels = 1:6)
  )

# Optional: reverse display order (uncomment if you want the list upside-down)
# df_plot <- df_plot |> mutate(Target = fct_rev(Target))

# 3) normal horizontal stacked bar plot
p <- ggplot(df_plot, aes(x = n_entry, y = Target, fill = k)) +
  geom_col(width = 0.8) +
  scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_y_discrete(position = "right") +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    
    # show Target axis on the right only
    axis.line.y.left  = element_blank(),
    axis.ticks.y.left = element_blank(),
    axis.text.y.left  = element_blank(),
    axis.line.y.right  = element_line(),
    axis.ticks.y.right = element_line(),
    axis.text.y.right  = element_text(hjust = 0),
    
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Number of entries",
    y = "ASV",
    fill = "Exact k\n(# of replicates)",
    title = "Occupancy consistency of bacterial biomarkers (winterPan)"
  )

p

ggsave("Network/Occupancy consistency bac_winter.pdf",p,height=6,width = 6)

### pieplot ----

# df_plot should already have:
# Target (factor, ordered), k (1–6), n_entry

# 1) Prepare proportions
df_pie <- df_plot |>
  group_by(Target) |>
  mutate(prop = n_entry / sum(n_entry)) |>
  ungroup()

# 2) Make long-row pie plot
p_pies <- ggplot(df_pie, aes(x = 1, y = prop, fill = k)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  facet_wrap(~ Target, nrow = 1) +   # <<< one long row
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    strip.background = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

p_pies
ggsave("Network/Occupancy_pies_bac_winter.pdf",p_pies,height=2,width=10)

# fun spring ----
## stack barplot----
# =========================
# 1) Build exact-k count table from k_by_entry
# =========================
# exact-k counts per Target
# 1) exact-k counts per Target
df_exact <- res_fun_spring$k_by_entry |>
  count(Target, k, name = "n_entry") |>
  filter(k %in% 1:6) |>
  tidyr::complete(Target, k = 1:6, fill = list(n_entry = 0))

# 2) set your desired Target order (no reversing)

df_plot <- df_exact |>
  mutate(
    Target = factor(Target,levels = c('OTU_1939','OTU_1233','OTU_4115','OTU_1836','OTU_3119','OTU_4309','OTU_4022','OTU_294','OTU_4182','OTU_4188','OTU_3575','OTU_863','OTU_4088','OTU_4190','OTU_4803','OTU_4263','OTU_3588','OTU_3665','OTU_4625','OTU_3978','OTU_3662','OTU_739')),
    k = factor(k, levels = 1:6)
  )

# Optional: reverse display order (uncomment if you want the list upside-down)
# df_plot <- df_plot |> mutate(Target = fct_rev(Target))

# 3) normal horizontal stacked bar plot
p <- ggplot(df_plot, aes(x = n_entry, y = Target, fill = k)) +
  geom_col(width = 0.8) +
  scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # show Target axis on the right only
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Number of entries",
    y = "",
    fill = "Exact k\n(# of replicates)",
    title = "Occupancy consistency of funterial biomarkers (SpringPan)"
  )

p

ggsave("Network/Occupancy consistency fun_spring.pdf",p,height=6,width = 6)

### pieplot ----

# df_plot should already have:
# Target (factor, ordered), k (1–6), n_entry

# 1) Prepare proportions
df_pie <- df_plot |>
  group_by(Target) |>
  mutate(prop = n_entry / sum(n_entry)) |>
  ungroup()

# 2) Make long-row pie plot
p_pies <- ggplot(df_pie, aes(x = 1, y = prop, fill = k)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  facet_wrap(~ Target, nrow = 1) +   # <<< one long row
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    strip.funkground = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    plot.funkground = element_rect(fill = "transparent", color = NA),
    panel.funkground = element_rect(fill = "transparent", color = NA)
  )

p_pies
ggsave("Network/Occupancy_pies_fun_spring.pdf",p_pies,height=2,width=10)

# fun winter ----
## stack barplot----
# =========================
# 1) Build exact-k count table from k_by_entry
# =========================
# exact-k counts per Target
# 1) exact-k counts per Target
df_exact <- res_fun_winter$k_by_entry |>
  count(Target, k, name = "n_entry") |>
  filter(k %in% 1:6) |>
  tidyr::complete(Target, k = 1:6, fill = list(n_entry = 0))

# 2) set your desired Target order (no reversing)


df_plot <- df_exact |>
  mutate(
    Target = factor(Target,levels=c("OTU_1067","OTU_2720","OTU_773","OTU_1241","OTU_2962")),
    k = factor(k, levels = 1:6)
  )

# Optional: reverse display order (uncomment if you want the list upside-down)
# df_plot <- df_plot |> mutate(Target = fct_rev(Target))

# 3) normal horizontal stacked bar plot
p <- ggplot(df_plot, aes(x = n_entry, y = Target, fill = k)) +
  geom_col(width = 0.8) +
  scale_x_continuous(limits = c(0, NA), breaks = scales::pretty_breaks()) +
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    # show Target axis on the right only
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Number of entries",
    y = "",
    fill = "Exact k\n(# of replicates)",
    title = "Occupancy consistency of funterial biomarkers (winterPan)"
  )

p

ggsave("Network/Occupancy consistency fun_winter.pdf",p,height=6,width = 6)

### pieplot ----

# df_plot should already have:
# Target (factor, ordered), k (1–6), n_entry

# 1) Prepare proportions
df_pie <- df_plot |>
  group_by(Target) |>
  mutate(prop = n_entry / sum(n_entry)) |>
  ungroup()

# 2) Make long-row pie plot
p_pies <- ggplot(df_pie, aes(x = 1, y = prop, fill = k)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  facet_wrap(~ Target, nrow = 1) +   # <<< one long row
  scale_fill_brewer(palette = "YlOrRd", direction = 1) +
  theme_void() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 9, face = "bold"),
    strip.funkground = element_blank(),
    panel.spacing.x = unit(0.8, "lines"),
    plot.funkground = element_rect(fill = "transparent", color = NA),
    panel.funkground = element_rect(fill = "transparent", color = NA)
  )

p_pies
ggsave("Network/Occupancy_pies_fun_winter.pdf",p_pies,height=2,width=10)




