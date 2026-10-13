# Figures: Fig4f
# Original source: Fun_Aut/mGWAS_format.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

rm(list=ls())

snp_hapmap <- read.delim("0_raw_data/BMB_Geno_Hapmap v4_2024_05_31.txt")
meta <- readRDS("1_processed_data/meta4snp.rds")
microbe <- readRDS("mGWAS/biomarkers_rclr_winter_fun.rds")


microbe_entry <- microbe[,colnames(microbe)%in%meta$sampleid] |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(select(meta,sampleid,Entry_lower)) |> 
  group_by(Entry_lower) |> 
  summarise(across(where(is.numeric),mean)) |> 
  column_to_rownames("Entry_lower") |> 
  as.data.frame() |> 
  t()

microbe_entry 
colnames(microbe_entry) <- tolower(colnames(microbe_entry))
colnames(snp_hapmap)    <- tolower(colnames(snp_hapmap))


common_ids <- intersect(colnames(microbe_entry), colnames(snp_hapmap))
length(common_ids)
common_ids[1:10]


anno_n <- 11  
hapmap_sample_cols <- colnames(snp_hapmap)[(anno_n + 1):ncol(snp_hapmap)]
keep_samples <- hapmap_sample_cols[hapmap_sample_cols %in% common_ids]

snp_hapmap_sub <- snp_hapmap[, c(colnames(snp_hapmap)[1:anno_n], keep_samples)]

microbe_entry_f <- microbe_entry[, common_ids, drop = FALSE]

microbe_entry_reordered <- microbe_entry_f[, keep_samples]

all(colnames(microbe_entry_reordered) == keep_samples)




write.table(
  snp_hapmap_sub,
  file = "mGWAS/fun_winter_snp_hapmap.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  t(microbe_entry_reordered),
  file = "mGWAS/fun_winter_biomarker.phe",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
) 

write.table(
  rownames(microbe_entry_reordered),
  file = "mGWAS/fun_winter_biomarker.phe.id",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
) 
