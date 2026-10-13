# Figures: Fig2a; FigS4,S5,S6a; TableS1
# Original source: Bac_Aut/SNP.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  cache = T
)


rm(list=ls())
library(tidyverse)
library(ape)


snp_alle <- readRDS("1_processed_data/snp_geno_alleles_aut.rds")


meta <- readRDS("1_processed_data/meta4snp.rds")


snp_alle <- snp_alle[which(rownames(snp_alle)!="abivoyager"),] |> as.matrix()
meta <- meta |> filter(Entry_lower!="abivoyager")


library(adegenet)
library(poppr)
## make a group vector for spring, winter and check lines

spring_entry <- filter(meta,Trial=="SpringPan")$Entry_lower
winter_entry <- filter(meta,Trial=="WinterPan")$Entry_lower

## shared check lines

shared_check_entry <- intersect(spring_entry, winter_entry)
shared_check_entry 

meta <- meta |> 
  mutate(
    snp_group = case_when(
      Entry_lower %in% spring_entry & !Entry_lower %in% shared_check_entry ~ "Spring Barley",
      Entry_lower %in% winter_entry & !Entry_lower %in% shared_check_entry ~ "Winter Barley",
      Entry_lower %in% shared_check_entry ~ "Shared Check lines",
      TRUE ~ "Unknown"   # optional: handle entries in none of the lists
    )
  )




# ==========================================================
# 0) INPUT
# ----------------------------------------------------------
# 1) Convert to genind (treat your per-locus nucleotide as haploid call)
# ----------------------------------------------------------
gi <- df2genind(snp_alle,
                ploidy = 1,       # your calls are A/C/G/T per locus (haploid-like encoding)
                type   = "codom") # codominant SNPs

# Optional: remove loci that are monomorphic (only 1 allele observed)
keep_loci <- locNames(gi)[nAll(gi) > 1]
gi <- gi[, loc=keep_loci, drop = TRUE]

# already done: remove loci with too much missing data (e.g., >20% missing)
# (poppr::info_table shows missingness; poppr::missingno filters)
#gi <- missingno(gi, type = "loci", cutoff = 0.20, quiet = TRUE)



cat("Samples:", nInd(gi), "  Loci after filtering:", nLoc(gi), "\n")


# ==========================================================
# 2) Prevosti’s absolute genetic distance
# ----------------------------------------------------------
D_prev <- prevosti.dist(gi)   # as.dist object


# Save distance if needed
# write.csv(as.matrix(D_prev), "prevosti_distance.csv", quote = FALSE)

# ==========================================================
# 3) PCoA on Prevosti distance
# ----------------------------------------------------------
pco <- ape::pcoa(D_prev)  # returns eigenvalues & vectors
scores <- as.data.frame(pco$vectors[, 1:2, drop = FALSE])
colnames(scores) <- c("PCoA1", "PCoA2")
scores$Entry_lower <- rownames(scores)

# Percent variance explained
var_pct <- round(100 * pco$values$Relative_eig[1:2], 2)


pc_df <- scores |> left_join(meta,multiple = "first")

library(ggthemes)
library(ggrepel)
# 4) 画图（可把 color 改成你的 metadata 列）
p_snp_pcoa_all <- 
ggplot(pc_df, aes(PCoA1, PCoA2,color=snp_group)) +
  geom_hline(yintercept=0,color="grey20",linetype="dashed")+
  geom_vline(xintercept=0,color="grey20",linetype="dashed")+
    geom_point(size = 2, alpha = 0.7) +   # add points for clarity
  geom_text_repel(aes(label = Entry), size = 2, max.overlaps = Inf) +

  theme_few() +
  labs(x = paste0("PCo1 (", var_pct[1], "%)"),
       y = paste0("PCo2 (", var_pct[2], "%)"))+
    scale_color_manual(values=c("grey2","#458a74","#cc5b45"))+
  theme(legend.position = "top",legend.position.inside = c(0.1,0.1))

p_snp_pcoa_all
ggsave("SNP/all_entry_snp_pcoa.pdf",p_snp_pcoa_all,height=12,width=12,dpi=300)


entry_list_spring <- filter(meta,snp_group=="Spring Barley")$Entry_lower |> 
  unique()
entry_list_winter <- filter(meta,snp_group=="Winter Barley")$Entry_lower  |>
  unique()


snp_alle_spring <- snp_alle[entry_list_spring,]
snp_alle_winter <- snp_alle[entry_list_winter,]


gi_spring <- df2genind(snp_alle_spring,
                ploidy = 1,       # your calls are A/C/G/T per locus (haploid-like encoding)
                type   = "codom") # codominant SNPs

# Optional: remove loci that are monomorphic (only 1 allele observed)
keep_loci <- locNames(gi_spring)[nAll(gi_spring) > 1]
gi_spring <- gi_spring[, loc=keep_loci, drop = TRUE]

# already done: remove loci with too much missing data (e.g., >20% missing)
# (poppr::info_table shows missingness; poppr::missingno filters)
gi_spring <- missingno(gi_spring, type = "loci", cutoff = 0.20, quiet = TRUE)



cat("Samples:", nInd(gi_spring), "  Loci after filtering:", nLoc(gi_spring), "\n")


# ==========================================================
# 2) Prevosti’s absolute genetic distance
# ----------------------------------------------------------
D_prev_spring <- prevosti.dist(gi_spring)   # as.dist object


gi_winter <- df2genind(snp_alle_winter,
                ploidy = 1,       # your calls are A/C/G/T per locus (haploid-like encoding)
                type   = "codom") # codominant SNPs

# Optional: remove loci that are monomorphic (only 1 allele observed)
keep_loci <- locNames(gi_winter)[nAll(gi_winter) > 1]
gi_winter <- gi_winter[, loc=keep_loci, drop = TRUE]

# already done: remove loci with too much missing data (e.g., >20% missing)
# (poppr::info_table shows missingness; poppr::missingno filters)
gi_winter <- missingno(gi_winter, type = "loci", cutoff = 0.20, quiet = TRUE)



cat("Samples:", nInd(gi_winter), "  Loci after filtering:", nLoc(gi_winter), "\n")


# ==========================================================
# 2) Prevosti’s absolute genetic distance
# ----------------------------------------------------------
D_prev_winter <- prevosti.dist(gi_winter)   # as.dist object


# Save distance if needed
saveRDS(D_prev_spring, "snp/prevosti_distance_spring.rds")
saveRDS(D_prev_winter, "snp/prevosti_distance_winter.rds")

# ==========================================================
# 3) PCoA on Prevosti distance
# ----------------------------------------------------------
pco_spring <- ape::pcoa(D_prev_spring)  # returns eigenvalues & vectors
scores_spring <- as.data.frame(pco_spring$vectors[, 1:2, drop = FALSE])
colnames(scores_spring) <- c("PCoA1", "PCoA2")
scores_spring$Entry_lower <- rownames(scores_spring)

# Percent variance explained
var_pct_spring <- round(100 * pco_spring$values$Relative_eig[1:2], 2)

pc_df_spring <- scores_spring |> left_join(meta,multiple = "first")

pco_winter <- ape::pcoa(D_prev_winter)  # returns eigenvalues & vectors
scores_winter <- as.data.frame(pco_winter$vectors[, 1:2, drop = FALSE])
colnames(scores_winter) <- c("PCoA1", "PCoA2")
scores_winter$Entry_lower <- rownames(scores_winter)

# Percent variance explained
var_pct_winter <- round(100 * pco_winter$values$Relative_eig[1:2], 2)

pc_df_winter <- scores_winter |> left_join(meta,multiple = "first")



library(ggthemes)
library(ggrepel)
# 4) plot
p_snp_pcoa_spring <- 
ggplot(pc_df_spring, aes(PCoA1, PCoA2,color=Group)) +
  geom_hline(yintercept=0,color="grey20",linetype="dashed")+
  geom_vline(xintercept=0,color="grey20",linetype="dashed")+
    geom_point(size = 2, alpha = 0.7) +   # add points for clarity
  #geom_text_repel(aes(label = Entry), size = 2, max.overlaps = Inf) +

  theme_few() +
  labs(x = paste0("PCo1 (", var_pct_spring[1], "%)"),
       y = paste0("PCo2 (", var_pct_spring[2], "%)"))+
    #scale_color_manual(values=c("grey2","#458a74","#cc5b45"))+
  theme(legend.position = "top",legend.position.inside = c(0.1,0.1))

p_snp_pcoa_spring
ggsave("SNP/spring_entry_snp_pcoa.pdf",p_snp_pcoa_all,height=12,width=12,dpi=300)


p_snp_pcoa_winter <- 
ggplot(pc_df_winter, aes(PCoA1, PCoA2,color=Group)) +
  geom_hline(yintercept=0,color="grey20",linetype="dashed")+
  geom_vline(xintercept=0,color="grey20",linetype="dashed")+
    geom_point(size = 2, alpha = 0.7) +   # add points for clarity
  #geom_text_repel(aes(label = Entry), size = 2, max.overlaps = Inf) +

  theme_few() +
  labs(x = paste0("PCo1 (", var_pct_winter[1], "%)"),
       y = paste0("PCo2 (", var_pct_winter[2], "%)"))+
    #scale_color_manual(values=c("grey2","#458a74","#cc5b45"))+
  theme(legend.position = "top",legend.position.inside = c(0.1,0.1))

p_snp_pcoa_winter


library(vegan)

sampleid_spring <- filter(meta,Entry_lower%in%entry_list_spring)$sampleid
sampleid_winter <- filter(meta,Entry_lower%in%entry_list_winter)$sampleid
 
asv_spring <- readRDS("1_processed_data/asv_spring.rds")
asv_winter <- readRDS("1_processed_data/asv_winter.rds")


asv_spring4snp <- asv_spring[,sampleid_spring]
asv_spring4snp_tss <- asv_spring4snp |> t() |> decostand("total")

asv_winter4snp <- asv_winter[,sampleid_winter]
asv_winter4snp_tss <- asv_winter4snp |> t() |> decostand("total")

# calculate the average relative abundance of ASVs across entries
asv_entry_spring4snp <- asv_spring4snp_tss |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(meta) |> 
  group_by(Entry_lower) |> 
  summarize(across(where(is.numeric),mean)) |>
  as.data.frame()

asv_entry_spring4snp |> head(5,5)
asv_entry_spring4snp <- asv_entry_spring4snp |> column_to_rownames("Entry_lower")

#### make sure the snp and asv table have the same order of entries ####
table(rownames(asv_entry_spring4snp)==rownames(snp_alle_spring))
asv_entry_spring4snp <- asv_entry_spring4snp[rownames(snp_alle_spring),]
table(rownames(asv_entry_spring4snp)==rownames(snp_alle_spring))
#----

#### calculate bray-curtis dist based on hellering transfromed data
D_asv_entry_spring <- asv_entry_spring4snp |> 
  decostand("hel") |> 
  vegdist("bray")

#----

####  winter  ####
asv_entry_winter4snp <- asv_winter4snp_tss |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(meta) |> 
  group_by(Entry_lower) |> 
  summarize(across(where(is.numeric),mean)) |>
  as.data.frame()

asv_entry_winter4snp |> head(5,5)
asv_entry_winter4snp <- asv_entry_winter4snp |> column_to_rownames("Entry_lower")

#### make sure the snp and asv table have the same order of entries ####
table(rownames(asv_entry_winter4snp)==rownames(snp_alle_winter))
asv_entry_winter4snp <- asv_entry_winter4snp[rownames(snp_alle_winter),]
table(rownames(asv_entry_winter4snp)==rownames(snp_alle_winter))
#----

#### calculate bray-curtis dist based on hellering transfromed data
D_asv_entry_winter <- asv_entry_winter4snp |> 
  decostand("hel") |> 
  vegdist("bray")


#----
saveRDS(D_asv_entry_spring,"1_processed_data/D_asv_entry_spring.rds")
saveRDS(D_asv_entry_winter,"1_processed_data/D_asv_entry_winter.rds")


saveRDS(D_prev_spring,"1_processed_data/D_prev_spring.rds")
saveRDS(D_prev_winter,"1_processed_data/D_prev_winter.rds")


mantel_spring <- mantel(D_asv_entry_spring,D_prev_spring,method="spearman",permutations = 9999)
mantel_winter <- mantel(D_asv_entry_winter,D_prev_winter,method="spearman",permutations = 9999)

mantel_spring
mantel_winter

## plot

df_mantel_entry_spring <- data.frame(
  Prevosti = as.numeric(D_prev_spring),
  BrayCurtis = as.numeric(D_asv_entry_spring))


p_mantel_entry_spring <- df_mantel_entry_spring |> 
  ggplot(aes(Prevosti,BrayCurtis))+
  geom_point(size=2,alpha=0.5)+
  geom_smooth(method="lm",color="green4")+
  theme_bw()+
  labs(x="Prevosti's genetic distance between genotypes",y="Bray Curtis's distance between genotypes",title="Bacteria SpringPan")

p_mantel_entry_spring


df_mantel_entry_winter <- data.frame(
  Prevosti = as.numeric(D_prev_winter),
  BrayCurtis = as.numeric(D_asv_entry_winter))

p_mantel_entry_winter <- df_mantel_entry_winter |> 
  ggplot(aes(Prevosti,BrayCurtis))+
  geom_point(size=2,alpha=0.5)+
  geom_smooth(method="lm",color="green4")+
  theme_bw()+
  labs(x="Prevosti's genetic distance between genotypes",y="Bray Curtis's distance between genotypes",title="Bacteria WinterPan")

p_mantel_entry_winter

ggsave("SNP/mantel_entry_spring_snp_bray.pdf",p_mantel_entry_spring, height = 6,width=6, dpi=300)
ggsave("SNP/mantel_entry_winter_snp_bray.pdf",p_mantel_entry_winter, height = 6,width=6, dpi=300)


## prensence absence?

D_asv_entry_spring_jaccard <- vegdist(asv_entry_spring4snp,binary = T,method="jaccard")

mantel_spring_jaccard <- mantel(D_asv_entry_spring_jaccard,D_prev_spring,method="spearman")
mantel_spring_jaccard


D_asv_entry_winter_jaccard <- vegdist(asv_entry_winter4snp,binary = T,method="jaccard")

mantel_winter_jaccard <- mantel(D_asv_entry_winter_jaccard,D_prev_winter,method="spearman")
mantel_winter_jaccard


snp_alle_sample <- readRDS("1_processed_data/snp_geno_alleles_aut_dup4sample.rds")


snp_sample_spring <- snp_alle_sample[sampleid_spring,]

snp_sample_winter <- snp_alle_sample[sampleid_winter,]


gi_spring_sample <- df2genind(snp_sample_spring,
                ploidy = 1,       # your calls are A/C/G/T per locus (haploid-like encoding)
                type   = "codom") # codominant SNPs

# Optional: remove loci that are monomorphic (only 1 allele observed)
keep_loci <- locNames(gi_spring_sample)[nAll(gi_spring_sample) > 1]
gi_spring_sample <- gi_spring_sample[, loc=keep_loci, drop = TRUE]

# already done: remove loci with too much missing data (e.g., >20% missing)
# (poppr::info_table shows missingness; poppr::missingno filters)
gi_spring_sample <- missingno(gi_spring_sample, type = "loci", cutoff = 0.20, quiet = TRUE)

D_prev_spring_sample <- prevosti.dist(gi_spring_sample)


### winter ###
gi_winter_sample <- df2genind(snp_sample_winter,
                ploidy = 1,       # your calls are A/C/G/T per locus (haploid-like encoding)
                type   = "codom") # codominant SNPs

# Optional: remove loci that are monomorphic (only 1 allele observed)
keep_loci <- locNames(gi_winter_sample)[nAll(gi_winter_sample) > 1]
gi_winter_sample <- gi_winter_sample[, loc=keep_loci, drop = TRUE]

# already done: remove loci with too much missing data (e.g., >20% missing)
# (poppr::info_table shows missingness; poppr::missingno filters)
gi_winter_sample <- missingno(gi_winter_sample, type = "loci", cutoff = 0.20, quiet = TRUE)

D_prev_winter_sample <- prevosti.dist(gi_winter_sample)



D_asv_spring <- vegdist(asv_spring4snp_tss) |> 
  decostand("hel") |> 
  vegdist("bray")


D_asv_winter <- vegdist(asv_winter4snp_tss) |> 
  decostand("hel") |> 
  vegdist("bray")

save.image("snp.Rdata")


mantel_sample_spring <- mantel(D_prev_spring_sample,D_asv_spring,method="spearman",parallel = 4)
mantel_sample_winter <- mantel(D_prev_winter_sample,D_asv_winter,method="spearman",parallel = 4)


meta_4snp_spring <- data.frame(sampleid=rownames(asv_spring4snp_tss)) |> 
  left_join(meta)

meta_4snp_winter <- data.frame(sampleid=rownames(asv_winter4snp_tss)) |> 
  left_join(meta)

adonis_spring_gxe <- adonis2(D_asv_spring~Entry+Col+Row+Entry*Col+Entry*Row,data=meta_4snp_spring,by="terms")
adonis_spring_gxe

adonis_winter_gxe <- adonis2(D_asv_winter~Entry+Col+Row+Entry*Col+Entry*Row,data=meta_4snp_winter,by="terms")
 adonis_winter_gxe


# 


# Prepare data for visualization

D_spatialxy_spring <- dist(meta_4snp_spring[,c("Col",'Row')])
D_spatialxy_winter <- dist(meta_4snp_winter[,c("Col",'Row')])

df_mantel_spring <- data.frame(
  Prevosti = as.numeric(D_prev_spring_sample),
  BrayCurtis = as.numeric(D_asv_spring),
  spatialxy=as.numeric(D_spatialxy_spring)
)

mantel_sample_spring_xy <- mantel(D_asv_spring,D_spatialxy_spring,method="spearman" ) 

p_mantel_spring <- 
df_mantel_spring |> ggplot()+
  geom_point(aes( scale(Prevosti),BrayCurtis),alpha=0.1)+
  geom_smooth(aes( scale(Prevosti),BrayCurtis),method="lm",se=F,color="green4")+
  geom_point(aes(scale(spatialxy),BrayCurtis),alpha=0.1)+
   geom_smooth(aes(scale(spatialxy),BrayCurtis),method="lm",se=F,color="brown4")+
 theme_few()
df_mantel_winter <- data.frame(
  Prevosti = as.numeric(D_prev_winter_sample),
  BrayCurtis = as.numeric(D_asv_winter),
  spatialxy=as.numeric(D_spatialxy_winter)
)

mantel_sample_winter_xy <- mantel(D_asv_winter,D_spatialxy_winter,method="spearman" ) 

p_mantel_winter <- 
  df_mantel_winter |> ggplot()+
  geom_point(aes( scale(Prevosti),BrayCurtis),alpha=0.1)+
  geom_smooth(aes( scale(Prevosti),BrayCurtis),method="lm",se=F,color="green4")+
  geom_point(aes(scale(spatialxy),BrayCurtis),alpha=0.1)+
  geom_smooth(aes(scale(spatialxy),BrayCurtis),method="lm",se=F,color="brown4")+
  theme_few()


p_mantel_winter


library(adespatial)


### 1. Hellinger-transform the asv/ASV table
asv_spring4snp_tss_f <- asv_spring4snp_tss[,colSums(asv_spring4snp_tss>0)>0.05*ncol(asv_spring4snp_tss)]
asv_spring4snp_tss_f |> dim()
asv_spring4snp_tss |> dim()
asv_spring_hel <- asv_spring4snp_tss_f |> decostand("hellinger")


### 2. Extract all PCNM eigenvectors
# pcnm_res_spring is obtained from pcnm(dist(coords))
pcnm_res_spring <- pcnm(D_spatialxy_spring)
pcnm_df_spring_all <- as.data.frame(pcnm_res_spring$vectors)

## Keep only PCNM with positive eigenvalues (spatially positive structures)
pos_idx <- which(pcnm_res_spring$values > 0)
pcnm_df_spring <- pcnm_df_spring_all[, pos_idx, drop = FALSE]

### 3. Fit a full RDA model using all PCNM vectors
# This provides the adjusted R² threshold for forward selection.
mod_full_spring <- rda(asv_spring_hel ~ ., data = pcnm_df_spring)
vif.cca(mod_full_spring)

global_adjR2_spring <- RsquareAdj(mod_full_spring)$adj.r.squared
global_adjR2_spring

### 4. Forward selection of significant PCNM vectors
# This step selects only those spatial variables that significantly
# explain community structure, preventing overfitting.
fs_spring <- forward.sel(
  Y = asv_spring_hel,
  X = pcnm_df_spring[,1:100],
  adjR2thresh = global_adjR2_spring,
  nperm = 999,
  alpha = 0.01,                    
  R2more = FALSE 
)

fs_spring

### 5. Extract selected PCNM vectors
selected_pcnm_names <- fs_spring$order
pcnm_sel_spring <- pcnm_df_spring[, selected_pcnm_names, drop = FALSE]
pcnm_sel_spring


### copare pcnm and xy

mod_xy <- capscale(asv_spring_hel ~ Row + Col, data = meta_4snp_spring)
RsquareAdj(mod_xy)                   # 看 R2 和 adj.R2
anova(mod_xy, permutations = 999)    # 看是否显著

mod_pcnm <- capscale(asv_spring_hel ~ ., data = pcnm_sel_spring)

RsquareAdj(mod_pcnm)
anova(mod_pcnm, permutations = 999)

### 6 db-RDA using selected PCNM + genotype PCs
head(fs_spring)
tail(fs_spring)
plot(fs_spring$AdjR2Cum, type = "b")
plot(fs_spring$R2, type = "b")


snp_pcoa_spring_pcos <- pco_spring$vectors |> as.data.frame() |> rownames_to_column("Entry_lower")
table(meta_4snp_spring$Entry_lower%in%snp_pcoa_spring_pcos$Entry_lower)
snp_pcoa_spring <- meta_4snp_spring |> select("sampleid","Entry_lower") |> 
  left_join(snp_pcoa_spring_pcos) |> select(-Entry_lower) |> column_to_rownames("sampleid")

fs_spring_snp <- forward.sel(
  Y = asv_spring_hel,
  X = snp_pcoa_spring,
  adjR2thresh = global_adjR2_spring,
  nperm = 999,
  alpha = 0.01,
  R2more = FALSE
)
fs_spring_snp 
sel_snp_spring <- fs_spring_snp$variables

table(rownames(snp_pcoa_spring)==meta_4snp_spring$sampleid)


X_snp_spring   <- snp_pcoa_spring[, sel_snp_spring, drop = FALSE]
colnames(X_snp_spring) <- gsub("^Axis\\.", "SNP", colnames(X_snp_spring))
X_space_spring <- pcnm_sel_spring
X_all_spring <- cbind(X_snp_spring,X_space_spring)


Y <- D_asv_spring 

## 1) only SNP
mod_snp   <- capscale(Y ~ ., data = X_snp_spring)

## 2) only Space (PCNM)
mod_space <- capscale(Y ~ ., data = X_space_spring)

## 3) SNP + Space 
mod_all   <- capscale(Y ~ ., data = X_all_spring)


RsquareAdj(mod_snp)$adj.r.squared
RsquareAdj(mod_space)$adj.r.squared
RsquareAdj(mod_all)$adj.r.squared

anova(mod_snp,   permutations = 999)
anova(mod_space, permutations = 999)
anova(mod_all,   permutations = 999)

anova_snp_space_terms <-  anova(mod_all, by = "terms", permutations = 999)  # check the significance of each SNP-PCoA axis or PCNMs
anova_snp_space_terms

# dbRDA_spring <- capscale(D_bc_spring~.,data=X_snp_spring)
# RsquareAdj(dbRDA_spring)$adj.r.squared
# anova.cca(dbRDA_spring,by="terms")


vp_spring <- varpart(D_asv_spring,
                     X_snp_spring,
                     X_space_spring)

vp_spring

plot(vp_spring, Xnames = c("SNP", "Space"))

vpa_spring_df <- data.frame(
  component = c("SNP", "Space", "Shared", "Unexplained"),
  R2 = c(
    vp_spring$part$indfract$Adj.R.square[1],  # SNP pure
    vp_spring$part$indfract$Adj.R.square[2],  # Space pure
    vp_spring$part$indfract$Adj.R.square[3],  # Shared
    1 - sum(vp_spring$part$indfract$Adj.R.square[1:3]) # Unexplained
  )
)
vpa_spring_df$component <- factor(vpa_spring_df$component, 
                           levels = c("SNP", "Space", "Shared", "Unexplained"))

vpa_spring_df


library("rdacca.hp")
iv_list_spring <- list(
  Space = X_space_spring,
  Geno  = X_snp_spring
)

## 4. Run rdacca.hp on Hellinger-transformed community
res_hp_spring <- rdacca.hp(
  dv       = D_asv_spring,      # response matrix
  iv       = iv_list_spring,    # list of matrices = groups (Env/Soil/Space/Geno)
  method   = "dbRDA",      # RDA / dbRDA / CCA
  type     = "adjR2",    # use adjusted R2
  var.part = TRUE        # also output full variation partitioning table
)

## Check the object structure
res_hp_spring
res_hp_spring$Total_explained_variation

plot(res_hp_spring, plot.perc = TRUE)


## Importance of each PCNM axis
res_pcnm_spring <- rdacca.hp(
  dv     = D_asv_spring,
  iv     = X_space_spring,    # all PCNM axes as individual predictors
  method = "dbRDA",
  type   = "adjR2"
)

pcnm_hp <- res_pcnm$Hier.part
pcnm_hp[order(pcnm_hp[,"I"], decreasing = TRUE), ]   # rank PCNM axes


library(corrplot)

X_all_spring <- cbind(X_space_spring, X_snp_spring)  # only selected axes
# 1. 先算整体相关矩阵
cor_mat_spring <- cor(X_all_spring, method = "spearman")

# 2. 找出 PCNM 和 SNP
pcnm_cols_spring <- grepl("^PCNM", colnames(cor_mat_spring))
snp_cols_spring  <- grepl("^SNP",  colnames(cor_mat_spring))

# 3. 子矩阵：相关系数
cor_mat_sub_spring <- cor_mat_spring[pcnm_cols_spring, snp_cols_spring]

# 4. 对整矩阵算 p 值
p_all_spring <- cor.mtest(cor_mat_spring)$p

# 5. 提取同样位置的子矩阵：p 值
p_mat_sub_spring <- p_all_spring[pcnm_cols_spring, snp_cols_spring]

# 6. 画图
corrplot(
  cor_mat_sub_spring,
  method   = "square",
  tl.cex   = 0.7,
  number.cex = 0.7,
  is.corr  = F,
  p.mat    = p_mat_sub_spring,
  sig.level = 0.05
)


### 1. Hellinger-transform the asv/ASV table
asv_winter4snp_tss_f <- asv_winter4snp_tss[,colSums(asv_winter4snp_tss>0)>0.05*ncol(asv_winter4snp_tss)]
asv_winter4snp_tss_f |> dim()
asv_winter4snp_tss |> dim()
asv_winter_hel <- asv_winter4snp_tss_f |> decostand("hellinger")


### 2. Extract all PCNM eigenvectors
# pcnm_res_winter is obtained from pcnm(dist(coords))
pcnm_res_winter <- pcnm(D_spatialxy_winter)
pcnm_df_winter_all <- as.data.frame(pcnm_res_winter$vectors)

## Keep only PCNM with positive eigenvalues (spatially positive structures)
pos_idx <- which(pcnm_res_winter$values > 0)
pcnm_df_winter <- pcnm_df_winter_all[, pos_idx, drop = FALSE]

### 3. Fit a full RDA model using all PCNM vectors
# This provides the adjusted R² threshold for forward selection.
mod_full_winter <- rda(asv_winter_hel ~ ., data = pcnm_df_winter)
vif.cca(mod_full_winter)

global_adjR2_winter <- RsquareAdj(mod_full_winter)$adj.r.squared
global_adjR2_winter

### 4. Forward selection of significant PCNM vectors
# This step selects only those spatial variables that significantly
# explain community structure, preventing overfitting.
fs_winter <- forward.sel(
  Y = asv_winter_hel,
  X = pcnm_df_winter[,1:100],
  adjR2thresh = global_adjR2_winter,
  nperm = 999,
  alpha = 0.01,                    
  R2more = FALSE 
)

fs_winter

### 5. Extract selected PCNM vectors
selected_pcnm_names <- fs_winter$order
pcnm_sel_winter <- pcnm_df_winter[, selected_pcnm_names, drop = FALSE]
pcnm_sel_winter


### compare pcnm and xy

mod_xy <- capscale(asv_winter_hel ~ Row + Col, data = meta_4snp_winter)
RsquareAdj(mod_xy)                   # 看 R2 和 adj.R2
anova(mod_xy, permutations = 999)    # 看是否显著

mod_pcnm <- capscale(asv_winter_hel ~ ., data = pcnm_sel_winter)

RsquareAdj(mod_pcnm)
anova(mod_pcnm, permutations = 999)

### 6 db-RDA using selected PCNM + genotype PCs
head(fs_winter)
tail(fs_winter)
plot(fs_winter$AdjR2Cum, type = "b")
plot(fs_winter$R2, type = "b")


snp_pcoa_winter_pcos <- pco_winter$vectors |> as.data.frame() |> rownames_to_column("Entry_lower")
table(meta_4snp_winter$Entry_lower%in%snp_pcoa_winter_pcos$Entry_lower)
snp_pcoa_winter <- meta_4snp_winter |> select("sampleid","Entry_lower") |> 
  left_join(snp_pcoa_winter_pcos) |> select(-Entry_lower) |> column_to_rownames("sampleid")

fs_winter_snp <- forward.sel(
  Y = asv_winter_hel,
  X = snp_pcoa_winter,
  adjR2thresh = global_adjR2_winter,
  nperm = 999,
  alpha = 0.01,
  R2more = FALSE
)
fs_winter_snp 
sel_snp_winter <- fs_winter_snp$variables

table(rownames(snp_pcoa_df)==meta_4snp_winter$sampleid)


X_snp_winter   <- snp_pcoa_df[, sel_snp_winter, drop = FALSE]
colnames(X_snp_winter) <- gsub("^Axis\\.", "SNP", colnames(X_snp_winter))
X_space_winter <- pcnm_sel_winter
X_all_winter <- cbind(X_snp_winter,X_space_winter)


Y <- D_asv_winter 

## 1) only SNP
mod_snp   <- capscale(Y ~ ., data = X_snp_winter)

## 2) only Space (PCNM)
mod_space <- capscale(Y ~ ., data = X_space_winter)

## 3) SNP + Space 
mod_all   <- capscale(Y ~ ., data = X_all_winter)


RsquareAdj(mod_snp)$adj.r.squared
RsquareAdj(mod_space)$adj.r.squared
RsquareAdj(mod_all)$adj.r.squared

anova(mod_snp,   permutations = 999)
anova(mod_space, permutations = 999)
anova(mod_all,   permutations = 999)

anova_snp_space_terms <-  anova(mod_all, by = "terms", permutations = 999)  # check the significance of each SNP-PCoA axis or PCNMs
anova_snp_space_terms



# dbRDA_winter <- capscale(D_bc_winter~.,data=X_snp_winter)
# RsquareAdj(dbRDA_winter)$adj.r.squared
# anova.cca(dbRDA_winter,by="terms")


vp_winter <- varpart(D_asv_winter,
                     X_snp_winter,
                     X_space_winter)

vp_winter

plot(vp_winter, Xnames = c("SNP", "Space"))

vpa_winter_df <- data.frame(
  component = c("SNP", "Space", "Shared", "Unexplained"),
  R2 = c(
    vp_winter$part$indfract$Adj.R.square[1],  # SNP pure
    vp_winter$part$indfract$Adj.R.square[2],  # Space pure
    vp_winter$part$indfract$Adj.R.square[3],  # Shared
    1 - sum(vp_winter$part$indfract$Adj.R.square[1:3]) # Unexplained
  )
)
vpa_winter_df$component <- factor(vpa_winter_df$component, 
                           levels = c("SNP", "Space", "Shared", "Unexplained"))

vpa_winter_df


library("rdacca.hp")
iv_list_winter <- list(
  Space = X_space_winter,
  Geno  = X_snp_winter
)

## 4. Run rdacca.hp on Hellinger-transformed community
res_hp_winter <- rdacca.hp(
  dv       = D_asv_winter,      # response matrix
  iv       = iv_list_winter,    # list of matrices = groups (Env/Soil/Space/Geno)
  method   = "dbRDA",      # RDA / dbRDA / CCA
  type     = "adjR2",    # use adjusted R2
  var.part = TRUE        # also output full variation partitioning table
)

## Check the object structure
res_hp_winter
res_hp_winter$Total_explained_variation


library(corrplot)

X_all_winter <- cbind(X_space_winter, X_snp_winter)  # only selected axes
# 1. 先算整体相关矩阵
cor_mat_winter <- cor(X_all_winter, method = "spearman")

# 2. 找出 PCNM 和 SNP
pcnm_cols_winter <- grepl("^PCNM", colnames(cor_mat_winter))
snp_cols_winter  <- grepl("^SNP",  colnames(cor_mat_winter))

# 3. 子矩阵：相关系数
cor_mat_sub_winter <- cor_mat_winter[pcnm_cols_winter, snp_cols_winter]

# 4. 对整矩阵算 p 值
p_all_winter <- cor.mtest(cor_mat_winter)$p

# 5. 提取同样位置的子矩阵：p 值
p_mat_sub_winter <- p_all_winter[pcnm_cols_winter, snp_cols_winter]

corrplot(
  cor_mat_sub_winter,
  method   = "square",
  tl.cex   = 0.7,
  number.cex = 0.7,
  is.corr  = F,
  p.mat    = p_mat_sub_winter,
  sig.level = 0.05
)


env <- read.csv("0_raw_data/AUS_ENV_20241112.csv")
meta_4snp_spring

rownames(X_space_spring) <- rownames(X_snp_spring)

rownames(X_space_winter) <- rownames(X_snp_winter)

x_space_spring_plot <- X_space_spring |> 
  rownames_to_column("sampleid") |> 
  left_join(select(meta_4snp_spring,sampleid,Plot)) |> 
  group_by(Plot) |> 
  summarise(across(where(is.numeric), mean)) |> 
  as.data.frame()

x_space_winter_plot <- X_space_winter |> 
  rownames_to_column("sampleid") |> 
  left_join(select(meta_4snp_winter,sampleid,Plot)) |> 
  group_by(Plot) |> 
  summarise(across(where(is.numeric), mean)) |> 
  as.data.frame()

env_spring <- env |> filter(Trial=="SpringPan")
env_winter <- env |> filter(Trial=="WinterPan")
env_spring$Plot <- env_spring$Plot-1200
env_winter$Plot
spring_comm_plot <- intersect(env_spring$Plot,x_space_spring_plot$Plot)
winter_comm_plot <- intersect(env_winter$Plot,x_space_winter_plot$Plot)

env_pncm_spring_d <- data.frame(Plot=spring_comm_plot) |> 
  left_join(env_spring) |> 
  left_join(x_space_spring_plot) |> 
  column_to_rownames("Plot") |> 
  select(-Trial)

env_pncm_spring_num <- env_pncm_spring_d |>
  select(where(is.numeric))   # 只保留 numeric 列

cor_mat_env_spring <- cor(env_pncm_spring_d,method="spearman", use    = "pairwise.complete.obs")

pcnm_plot_cols_spring <- grep("^PCNM", colnames(cor_mat_env_spring), value = TRUE)

env_plot_cols_spring  <- setdiff(colnames(cor_mat_env_spring),pcnm_plot_cols_spring)


cor_mat_env_sub_spring <- cor_mat_env_spring[pcnm_plot_cols_spring, env_plot_cols_spring]




p_env_spring <- cor.mtest(cor_mat_env_spring,method="spearman" )$p

# 5. 提取同样位置的子矩阵：p 值
p_env_sub_spring <- p_env_spring[pcnm_plot_cols_spring, env_plot_cols_spring]

corrplot(
  cor_mat_env_sub_spring,
  method   = "square",
  tl.cex   = 0.7,
  number.cex = 0.7,
  is.corr  = TRUE,
  p.mat    = p_env_sub_spring,
  sig.level = 0.05,
  insig = "blank"
)
cor_mat_env_sub_spring |> hist()
cor_mat_env_sub_spring |> max()
  cor_mat_env_sub_spring |> min()



env_pncm_winter_d <- data.frame(Plot=winter_comm_plot) |> 
  left_join(env_winter) |> 
  left_join(x_space_winter_plot) |> 
  column_to_rownames("Plot") |> 
  select(-Trial)

env_pncm_winter_num <- env_pncm_winter_d |>
  select(where(is.numeric))   # 只保留 numeric 列

cor_mat_env_winter <- cor(env_pncm_winter_d,method="spearman", use    = "pairwise.complete.obs")

pcnm_plot_cols_winter <- grep("^PCNM", colnames(cor_mat_env_winter), value = TRUE)

env_plot_cols_winter  <- setdiff(colnames(cor_mat_env_winter),pcnm_plot_cols_winter)


cor_mat_env_sub_winter <- cor_mat_env_winter[pcnm_plot_cols_winter, env_plot_cols_winter]




p_env_winter <- cor.mtest(cor_mat_env_winter,method="spearman" )$p

# 5. 提取同样位置的子矩阵：p 值
p_env_sub_winter <- p_env_winter[pcnm_plot_cols_winter, env_plot_cols_winter]

corrplot(
  cor_mat_env_sub_winter,
  method   = "square",
  tl.cex   = 0.7,
  number.cex = 0.7,
  is.corr  = TRUE,
  p.mat    = p_env_sub_winter,
  sig.level = 0.05,
  insig = "blank"
)

cor_mat_env_sub_winter |> hist()
cor_mat_env_sub_winter |> max()
cor_mat_env_sub_winter |> min()


library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidygraph)

## ---- 1. Prepare edge table from correlation & p-value matrices ----------
# Convert correlation matrix into long format
df_cor_spring <- as.data.frame(cor_mat_env_sub_spring) |>
  tibble::rownames_to_column("PCNM") |>
  pivot_longer(-PCNM, names_to = "ENV", values_to = "cor")

# Convert p-value matrix into long format
df_p_spring <- as.data.frame(p_env_sub_spring) |>
  tibble::rownames_to_column("PCNM") |>
  pivot_longer(-PCNM, names_to = "ENV", values_to = "p")

# Merge correlation and p-values, and apply significance thresholds
edges_spring <- df_cor_spring |>
  left_join(df_p_spring, by = c("PCNM","ENV")) |>
  filter(p < 0.05)   # Thresholds can be adjusted


## ---- 2. Build node table including manual x/y coordinates ---------------
# Extract unique node names from both PCNM and ENV sides
nodes_spring <- unique(c(edges_spring$PCNM, edges_spring$ENV))

# Create node metadata and assign coordinates
node_df_spring <- tibble(
  name = nodes_spring,
  type = ifelse(grepl("^PCNM", nodes_spring), "PCNM", "ENV")
) |>
  group_by(type) |>
  arrange(type, name, .by_group = TRUE) |>
  mutate(
    # PCNM nodes placed on the left (x = 0), ENV on the right (x = 1)
    x = ifelse(type == "PCNM", 0, 1),

    # Vertical order within each category
    y = row_number()
  ) |>
  ungroup()


## ---- 3. Construct graph object ------------------------------------------
# Create igraph object with node metadata (including x/y)
g_spring <- graph_from_data_frame(
  edges_spring,
  directed = FALSE,
  vertices = node_df_spring
)

# Convert to tbl_graph for ggraph compatibility
g_tbl_spring <- as_tbl_graph(g_spring)


## ---- 4. Plot bipartite correlation network -------------------------------
p_spring <- ggraph(g_tbl_spring, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(
      edge_alpha = abs(cor),
      edge_width = abs(cor),
      edge_color = cor
    ),
    show.legend = TRUE
  ) +
  scale_edge_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Correlation"
  ) +
  scale_edge_width(range = c(0.3, 1.5), guide = "none") +
  scale_edge_alpha(range = c(0.3, 1.0), guide = "none") +

  geom_node_point(aes(color = type), size = 3) +
  scale_color_manual(values = c(PCNM = "#1f78b4", ENV = "#33a02c"),
                     name = "") +

  geom_node_text(
    aes(
      label = name,
      hjust = ifelse(type == "PCNM", 1.1, -0.1)
    ),
    size = 3
  ) +

  xlim(-0.2, 1.2) +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )
p_spring




library(dplyr)
library(tidyr)
library(circlize)
## spring
## 1. Get unique PCNM and ENV names -------------------------------
pcnm_nodes <- sort(unique(edges_spring$PCNM))
env_nodes  <- sort(unique(edges_spring$ENV))

## 2. Combined order around the circle: first all PCNM, then all ENV
all_nodes_order <- c(pcnm_nodes, env_nodes)

## 3. Define grid colors for sectors (nodes) ----------------------
# PCNM nodes: blue, ENV nodes: green
grid_col <- c(
  setNames(rep("#1f78b4", length(pcnm_nodes)), pcnm_nodes),
  setNames(rep("#33a02c", length(env_nodes)),  env_nodes)
)

## 4. Define link color mapping based on correlation --------------
# Negative: blue, zero: white, positive: red
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),
  colors = c("blue", "white", "red")
)


## 5. Clear any previous circlize plots --------------------------
circos.clear()

## 6. Prepare link colors for each edge (row) --------------------
# One color per correlation value
link_cols <- col_fun(edges_spring$cor)

## 7. Draw chord diagram (keep labels, no ticks) -----------------
chordDiagram(
  x = edges_spring |> dplyr::select(PCNM, ENV, cor),
  order    = all_nodes_order,        # sector order on the circle
  grid.col = grid_col,               # sector colors
  transparency = 0.3,                # link transparency
  col = link_cols,                   # link colors based on correlation
  annotationTrack = c("name", "grid")  # <- keep names + grid, drop axis ticks
  
)


library(dplyr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(tidygraph)

## ---- 1. Prepare edge table from correlation & p-value matrices ----------
# Convert correlation matrix into long format
df_cor_winter <- as.data.frame(cor_mat_env_sub_winter) |>
  tibble::rownames_to_column("PCNM") |>
  pivot_longer(-PCNM, names_to = "ENV", values_to = "cor")

# Convert p-value matrix into long format
df_p_winter <- as.data.frame(p_env_sub_winter) |>
  tibble::rownames_to_column("PCNM") |>
  pivot_longer(-PCNM, names_to = "ENV", values_to = "p")

# Merge correlation and p-values, and apply significance thresholds
edges_winter <- df_cor_winter |>
  left_join(df_p_winter, by = c("PCNM","ENV")) |>
  filter(p < 0.05)   # Thresholds can be adjusted


## ---- 2. Build node table including manual x/y coordinates ---------------
# Extract unique node names from both PCNM and ENV sides
nodes_winter <- unique(c(edges_winter$PCNM, edges_winter$ENV))

# Create node metadata and assign coordinates
node_df_winter <- tibble(
  name = nodes_winter,
  type = ifelse(grepl("^PCNM", nodes_winter), "PCNM", "ENV")
) |>
  group_by(type) |>
  arrange(type, name, .by_group = TRUE) |>
  mutate(
    # PCNM nodes placed on the left (x = 0), ENV on the right (x = 1)
    x = ifelse(type == "PCNM", 0, 1),

    # Vertical order within each category
    y = row_number()
  ) |>
  ungroup()


## ---- 3. Construct graph object ------------------------------------------
# Create igraph object with node metadata (including x/y)
g_winter <- graph_from_data_frame(
  edges_winter,
  directed = FALSE,
  vertices = node_df_winter
)

# Convert to tbl_graph for ggraph compatibility
g_tbl_winter <- as_tbl_graph(g_winter)


## ---- 4. Plot bipartite correlation network -------------------------------
p_winter <- ggraph(g_tbl_winter, layout = "manual", x = x, y = y) +
  geom_edge_link(aes(
      edge_alpha = abs(cor),
      edge_width = abs(cor),
      edge_color = cor
    ),
    show.legend = TRUE
  ) +
  scale_edge_color_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0, name = "Correlation"
  ) +
  scale_edge_width(range = c(0.3, 1.5), guide = "none") +
  scale_edge_alpha(range = c(0.3, 1.0), guide = "none") +

  geom_node_point(aes(color = type), size = 3) +
  scale_color_manual(values = c(PCNM = "#1f78b4", ENV = "#33a02c"),
                     name = "") +

  geom_node_text(
    aes(
      label = name,
      hjust = ifelse(type == "PCNM", 1.1, -0.1)
    ),
    size = 3
  ) +

  xlim(-0.2, 1.2) +
  theme(
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )
p_winter




library(dplyr)
library(tidyr)
library(circlize)
## winter
## 1. Get unique PCNM and ENV names -------------------------------
pcnm_nodes <- sort(unique(edges_winter$PCNM))
env_nodes  <- sort(unique(edges_winter$ENV))

## 2. Combined order around the circle: first all PCNM, then all ENV
all_nodes_order <- c(pcnm_nodes, env_nodes)

## 3. Define grid colors for sectors (nodes) ----------------------
# PCNM nodes: blue, ENV nodes: green
grid_col <- c(
  setNames(rep("#1f78b4", length(pcnm_nodes)), pcnm_nodes),
  setNames(rep("#33a02c", length(env_nodes)),  env_nodes)
)

## 4. Define link color mapping based on correlation --------------
# Negative: blue, zero: white, positive: red
col_fun <- colorRamp2(
  breaks = c(-1, 0, 1),
  colors = c("blue", "white", "red")
)


## 5. Clear any previous circlize plots --------------------------
circos.clear()

## 6. Prepare link colors for each edge (row) --------------------
# One color per correlation value
link_cols <- col_fun(edges_winter$cor)

## 7. Draw chord diagram (keep labels, no ticks) -----------------
chordDiagram(
  x = edges_winter |> dplyr::select(PCNM, ENV, cor),
  order    = all_nodes_order,        # sector order on the circle
  grid.col = grid_col,               # sector colors
  transparency = 0.3,                # link transparency
  col = link_cols,                   # link colors based on correlation
  annotationTrack = c("name", "grid")  # <- keep names + grid, drop axis ticks
  
)


save.image("SNP/snp.Rdata")
