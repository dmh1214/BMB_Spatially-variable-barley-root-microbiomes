# Figures: Fig2b; FigS6b
# Original source: 10_PCNM_VPA/pcnm_vpa.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

rm(list=ls())
## package ----
library(adespatial)
library(tidyverse)
library(vegan)
library(ggcorrplot)
library(corrplot)
## data----
meta_raw <- read.csv("Bac_Aut/0_raw_data/meta_file.csv")
meta <- meta_raw |> filter(Location=="Austria",Trial!="DivPan")
meta_spring <- meta |> filter(Trial=="SpringPan")
meta_winter <- meta |> filter(Trial=="WinterPan")


spatialXY <- meta_spring |> select(sampleid,Col,Row,Side) |> column_to_rownames("sampleid") # spring and winter panels have the same layout of plots

## Plots were splited into 2 subplots from row, so adjust the row with 0.25 and 0.75 (N-0.25, S-0.75)

spatialXY <- spatialXY |> 
  mutate(Row=case_when(Side=="N"~Row-0.25,
                       Side=="S"~ Row-0.75)) |> 
  select(Row,Col) |> 
  as.matrix()


D_spatialXY <- dist(spatialXY, method = "euclidean")

PCNMs <-pcnm(D_spatialXY)


pcnm_df_all <- as.data.frame(PCNMs$vectors)
pcnm_df_all <- pcnm_df_all[, PCNMs$values > 0, drop = FALSE]


env <- read.csv("Bac_Aut/0_raw_data/AUS_ENV_20241112.csv")

env_spring <- env |> filter(Trial=="SpringPan") 
env_winter <- env |> filter(Trial=="WinterPan")
# correct the naming 
#(BOKU ranked Plots in Spring and Winter panles together, so -1200 for springpan)
env_spring$Plot <- env_spring$Plot-1200
env_winter$Plot


# correlation between env and pcnms-----
pcnm30 <- pcnm_df_all[,1:30] |> 
  rownames_to_column("sampleid") |> 
  left_join(select(meta,sampleid,Plot)) |> 
  select(-sampleid)

env_spring_pcnm <- env_spring |> 
  left_join(pcnm30,multiple = "first") |>
  select(where(is.numeric)) |> 
  select(-Plot)

env_spring_pcnm |> dim()

# cor mat
cor_mat_env_spring <- cor(env_spring_pcnm, method = "spearman",use    = "pairwise.complete.obs")
pcnm_plot_cols_spring <- grep("^PCNM", colnames(cor_mat_env_spring ), value = TRUE)
env_plot_cols_spring  <- setdiff(colnames(cor_mat_env_spring ),pcnm_plot_cols_spring)
cor_mat_env_sub_spring <- cor_mat_env_spring[pcnm_plot_cols_spring, env_plot_cols_spring]



# p mat
p_env_spring <- cor.mtest(env_spring_pcnm,method="spearman")$p
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


## winter

env_winter_pcnm <- env_winter |> 
  left_join(pcnm30,multiple = "first") |>
  select(where(is.numeric)) |> 
  select(-Plot)

# cor mat
cor_mat_env_winter <- cor(env_winter_pcnm, method = "spearman",use    = "pairwise.complete.obs")
pcnm_plot_cols_winter <- grep("^PCNM", colnames(cor_mat_env_winter ), value = TRUE)
env_plot_cols_winter  <- setdiff(colnames(cor_mat_env_winter ),pcnm_plot_cols_winter)
cor_mat_env_sub_winter <- cor_mat_env_winter[pcnm_plot_cols_winter, env_plot_cols_winter]



# p mat
p_env_winter <- cor.mtest(env_winter_pcnm,method="spearman")$p
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


bac_spring <- readRDS("Bac_Aut/1_processed_data/asv_spring.rds")
bac_winter <- readRDS("Bac_Aut/1_processed_data/asv_winter.rds")

fun_spring <- readRDS("Fun_Aut/1_processed_data/otu_spring.rds")
fun_winter <- readRDS("Fun_Aut/1_processed_data/otu_winter.rds")


bac_spring <- bac_spring[rowSums(bac_spring>0)>0.05*ncol(bac_spring),]
bac_winter<- bac_winter[rowSums(bac_winter>0)>0.05*ncol(bac_winter),]

fun_spring <- fun_spring[rowSums(fun_spring>0)>0.05*ncol(fun_spring),]
fun_winter<- fun_winter[rowSums(fun_winter>0)>0.05*ncol(fun_winter),]


meta4snp_spring_bac <- readRDS("Bac_Aut/1_processed_data/meta4snp_spring.rds") |> filter(Panel=="SpringPan") # only use unique lines for each panel, make sure comparing within spring/winter barley
meta4snp_winter_bac <- readRDS("Bac_Aut/1_processed_data/meta4snp_winter.rds")|> filter(Panel=="WinterPan")
meta4snp_spring_fun <- readRDS("Fun_Aut/1_processed_data/meta4snp_spring.rds")|> filter(Panel=="SpringPan")
meta4snp_winter_fun <- readRDS("Fun_Aut/1_processed_data/meta4snp_winter.rds")|> filter(Panel=="WinterPan")


meta4snp_spring_bac |> head()

bac_spring4snp <- bac_spring |> select(meta4snp_spring_bac$sampleid)
bac_spring4snp |> dim()
bac_spring4snp_hel <- bac_spring4snp |> 
  t() |>
  decostand("total") |> 
  decostand("hel") 

bac_winter4snp <- bac_winter |> select(meta4snp_winter_bac$sampleid)
bac_winter4snp |> dim()
bac_winter4snp_hel <- bac_winter4snp |> 
  t() |>
  decostand("total") |> 
  decostand("hel") 

fun_spring4snp <- fun_spring |>
  as.data.frame() |>
  select(meta4snp_spring_fun$sampleid)

fun_spring4snp |> dim()
fun_spring4snp_hel <- fun_spring4snp |> 
  t() |>
  decostand("total") |> 
  decostand("hel") 


fun_winter4snp <- fun_winter |>   as.data.frame() |>select(meta4snp_winter_fun$sampleid)

fun_winter4snp |> dim()
fun_winter4snp_hel <- fun_winter4snp |> 
  t() |>
  decostand("total") |> 
  decostand("hel") 

# genetic distance data based on 50k SNPs
D_snp_bac_spring <- readRDS("Bac_Aut/SNP/prevosti_distance_spring.rds")
D_snp_bac_winter <- readRDS("Bac_Aut/SNP/prevosti_distance_winter.rds")
D_snp_fun_spring <- readRDS("Fun_Aut/SNP/prevosti_distance_spring.rds")
D_snp_fun_winter <- readRDS("Fun_Aut/SNP/prevosti_distance_winter.rds")


pcnm_df <- pcnm_df_all |> rownames_to_column("sampleid")

pcnm_df4rda_spring_bac <- meta4snp_spring_bac |> 
  select(sampleid) |>
  left_join(pcnm_df) 

pcnm_df4rda_winter_bac <- meta4snp_winter_bac |> 
  select(Plot,Side,sampleid) |>
  left_join(select(meta_spring,sampleid,Plot,Side),join_by(Plot,Side)) |>
  left_join(pcnm_df,join_by(sampleid.y==sampleid)) |> 
  select(-Plot,-Side,-sampleid.y)

pcnm_df4rda_spring_fun <- meta4snp_spring_fun |> select(sampleid) |> 
  left_join(pcnm_df)


pcnm_df4rda_winter_fun <- meta4snp_winter_fun |> 
  select(Plot,Side,sampleid) |>
  left_join(select(meta_spring,sampleid,Plot,Side),join_by(Plot,Side)) |>
  left_join(pcnm_df,join_by(sampleid.y==sampleid)) |> 
 select(-Plot,-Side,-sampleid.y)


table(pcnm_df4rda_spring_bac$sampleid==rownames(bac_spring4snp_hel)) #check if they have same order
table(pcnm_df4rda_winter_bac$sampleid==rownames(bac_winter4snp_hel))


comm <- bac_spring4snp_hel
dim(comm)
rda_pcnm_bac_spring <- capscale(comm ~ ., data = pcnm_df4rda_spring_bac[,2:101], distance = "bray")

RsquareAdj(rda_pcnm_bac_spring)$adj.r.squared


pcnm_sel_spring_bac <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=pcnm_df4rda_spring_bac[,2:101]),
  scope = formula(rda_pcnm_bac_spring),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


pcoa_snp_bac_spring <- cmdscale(D_snp_bac_spring,k=100,eig=T)
sum(pcoa_snp_bac_spring$eig[1:100])/sum(pcoa_snp_bac_spring$eig)

pcoa_snp_bac_spring_points <- pcoa_snp_bac_spring$points |> 
  as.data.frame() |> 
  rownames_to_column("Entry_lower") |> 
  left_join(meta4snp_spring_bac,join_by(Entry_lower))

snp_df4rda_bac_spring <- meta4snp_spring_bac |> select(sampleid,Plot) |> 
  left_join(pcoa_snp_bac_spring_points)

snp_df4rda_bac_spring[is.na(snp_df4rda_bac_spring$V1),]$Entry_lower |> unique() # check lines should be removed from the analysi
snp_df4rda_bac_spring_clean <- snp_df4rda_bac_spring[!is.na(snp_df4rda_bac_spring$V1),]


comm <- bac_spring4snp_hel[snp_df4rda_bac_spring_clean$sampleid, ]

table(rownames(comm)==snp_df4rda_bac_spring_clean$sampleid)

rda_snp_bac_spring <- capscale(comm ~ ., data =snp_df4rda_bac_spring_clean[,4:103], distance = "bray")

snp_sel_spring_bac <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=snp_df4rda_bac_spring_clean[,4:103]),
  scope = formula(rda_snp_bac_spring),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


sel_pcnm_spring_bac <- attr(terms(pcnm_sel_spring_bac), "term.labels")
sel_snp_spring_bac  <- attr(terms(snp_sel_spring_bac),  "term.labels")

id_comm <- rownames(bac_spring4snp_hel)

id_pcnm <- pcnm_df4rda_spring_bac$sampleid
id_snp  <- snp_df4rda_bac_spring_clean$sampleid

id_keep <- Reduce(intersect, list(id_comm, id_pcnm, id_snp))
length(id_keep)
comm_bac_spring <- bac_spring4snp_hel[id_keep, , drop = FALSE]
stopifnot(all(rownames(comm_bac_spring) == id_keep))


tmp_snp <- snp_df4rda_bac_spring |>
  filter(sampleid %in% id_keep) |>
  arrange(match(sampleid, id_keep))

tmp_pcnm <- pcnm_df4rda_spring_bac |>
  filter(sampleid %in% id_keep) |>
  arrange(match(sampleid, id_keep))

stopifnot(all(tmp_snp$sampleid  == rownames(comm_bac_spring)))
stopifnot(all(tmp_pcnm$sampleid == rownames(comm_bac_spring)))


X_snp_spring_bac  <- tmp_snp  #|> select(all_of(sel_snp_spring_bac))
X_pcnm_spring_bac <- tmp_pcnm #|> select(all_of(sel_pcnm_spring_bac))

vpa_bac_spring <- varpart(vegdist(comm_bac_spring,method="bray"), X_pcnm_spring_bac[,sel_pcnm_spring_bac], X_snp_spring_bac[,sel_snp_spring_bac ])
plot(vpa_bac_spring)
vpa_bac_spring 


comm <- bac_winter4snp_hel
dim(comm)
rda_pcnm_bac_winter <- capscale(comm ~ ., data = pcnm_df4rda_winter_bac[,2:101], distance = "bray")

RsquareAdj(rda_pcnm_bac_winter)$adj.r.squared


pcnm_sel_winter_bac <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=pcnm_df4rda_winter_bac[,2:101]),
  scope = formula(rda_pcnm_bac_winter),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


pcoa_snp_bac_winter <- cmdscale(D_snp_bac_winter,k=100,eig=T)
sum(pcoa_snp_bac_winter$eig[1:100])/sum(pcoa_snp_bac_winter$eig)

pcoa_snp_bac_winter_points <- pcoa_snp_bac_winter$points |> 
  as.data.frame() |> 
  rownames_to_column("Entry_lower") |> 
  left_join(meta4snp_winter_bac,join_by(Entry_lower))

snp_df4rda_bac_winter <- meta4snp_winter_bac |> select(sampleid,Plot) |> 
  left_join(pcoa_snp_bac_winter_points)

snp_df4rda_bac_winter[is.na(snp_df4rda_bac_winter$V1),]$Entry_lower |> unique() # check lines should be removed from the analysi
snp_df4rda_bac_winter_clean <- snp_df4rda_bac_winter[!is.na(snp_df4rda_bac_winter$V1),]


comm <- bac_winter4snp_hel[snp_df4rda_bac_winter_clean$sampleid, ]

table(rownames(comm)==snp_df4rda_bac_winter_clean$sampleid)

rda_snp_bac_winter <- capscale(comm ~ ., data =snp_df4rda_bac_winter_clean[,4:103], distance = "bray")

snp_sel_winter_bac <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=snp_df4rda_bac_winter_clean[,4:103]),
  scope = formula(rda_snp_bac_winter),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


sel_pcnm_winter_bac <- attr(terms(pcnm_sel_winter_bac), "term.labels")
sel_snp_winter_bac  <- attr(terms(snp_sel_winter_bac),  "term.labels")

id_comm <- rownames(bac_winter4snp_hel)

id_pcnm <- pcnm_df4rda_winter_bac$sampleid
id_snp  <- snp_df4rda_bac_winter_clean$sampleid

id_keep <- Reduce(intersect, list(id_comm, id_pcnm, id_snp))
length(id_keep)
comm_bac_winter <- bac_winter4snp_hel[id_keep, , drop = FALSE]
stopifnot(all(rownames(comm_bac_winter) == id_keep))


tmp_snp <- snp_df4rda_bac_winter |>
  filter(sampleid %in% id_keep) |>
  arrange(match(sampleid, id_keep))

tmp_pcnm <- pcnm_df4rda_winter_bac |>
  filter(sampleid.x %in% id_keep) |>
  arrange(match(sampleid.x, id_keep))

stopifnot(all(tmp_snp$sampleid  == rownames(comm_bac_winter)))
stopifnot(all(tmp_pcnm$sampleid == rownames(comm_bac_winter)))


X_snp_winter_bac  <- tmp_snp  #|> select(all_of(sel_snp_winter_bac))
X_pcnm_winter_bac <- tmp_pcnm #|> select(all_of(sel_pcnm_winter_bac))

vpa_bac_winter <- varpart(vegdist(comm_bac_winter,method="bray"), X_pcnm_winter_bac[,sel_pcnm_winter_bac], X_snp_winter_bac[,sel_snp_winter_bac])
plot(vpa_bac_winter)
vpa_bac_winter 


comm <- fun_spring4snp_hel
dim(comm)
rda_pcnm_fun_spring <- capscale(comm ~ ., data = pcnm_df4rda_spring_fun[,2:101], distance = "bray")

RsquareAdj(rda_pcnm_fun_spring)$adj.r.squared


pcnm_sel_spring_fun <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=pcnm_df4rda_spring_fun[,2:101]),
  scope = formula(rda_pcnm_fun_spring),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


pcoa_snp_fun_spring <- cmdscale(D_snp_fun_spring,k=100,eig=T)
sum(pcoa_snp_fun_spring$eig[1:100])/sum(pcoa_snp_fun_spring$eig)

pcoa_snp_fun_spring_points <- pcoa_snp_fun_spring$points |> 
  as.data.frame() |> 
  rownames_to_column("Entry_lower") |> 
  left_join(meta4snp_spring_fun,join_by(Entry_lower))

snp_df4rda_fun_spring <- meta4snp_spring_fun |> select(sampleid,Plot) |> 
  left_join(pcoa_snp_fun_spring_points)

snp_df4rda_fun_spring[is.na(snp_df4rda_fun_spring$V1),]$Entry_lower |> unique() # check lines should be removed from the analysi
snp_df4rda_fun_spring_clean <- snp_df4rda_fun_spring[!is.na(snp_df4rda_fun_spring$V1),]


comm <- fun_spring4snp_hel[snp_df4rda_fun_spring_clean$sampleid, ]

table(rownames(comm)==snp_df4rda_fun_spring_clean$sampleid)

rda_snp_fun_spring <- capscale(comm ~ ., data =snp_df4rda_fun_spring_clean[,4:103], distance = "bray")

snp_sel_spring_fun <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=snp_df4rda_fun_spring_clean[,4:103]),
  scope = formula(rda_snp_fun_spring),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


sel_pcnm_spring_fun <- attr(terms(pcnm_sel_spring_fun), "term.labels")
sel_snp_spring_fun  <- attr(terms(snp_sel_spring_fun),  "term.labels")

id_comm <- rownames(fun_spring4snp_hel)

id_pcnm <- pcnm_df4rda_spring_fun$sampleid
id_snp  <- snp_df4rda_fun_spring_clean$sampleid

id_keep <- Reduce(intersect, list(id_comm, id_pcnm, id_snp))
length(id_keep)
comm_fun_spring <- fun_spring4snp_hel[id_keep, , drop = FALSE]
stopifnot(all(rownames(comm_fun_spring) == id_keep))


tmp_snp <- snp_df4rda_fun_spring |>
  filter(sampleid %in% id_keep) |>
  arrange(match(sampleid, id_keep))

tmp_pcnm <- pcnm_df4rda_spring_fun |>
  filter(sampleid %in% id_keep) |>
  arrange(match(sampleid, id_keep))

stopifnot(all(tmp_snp$sampleid  == rownames(comm_fun_spring)))
stopifnot(all(tmp_pcnm$sampleid == rownames(comm_fun_spring)))


X_snp_spring_fun  <- tmp_snp  #|> select(all_of(sel_snp_spring_fun))
X_pcnm_spring_fun <- tmp_pcnm #|> select(all_of(sel_pcnm_spring_fun))

vpa_fun_spring <- varpart(vegdist(comm_fun_spring,method="bray"), X_pcnm_spring_fun[,sel_pcnm_spring_fun], X_snp_spring_fun[,sel_snp_spring_fun ])
plot(vpa_fun_spring)
vpa_fun_spring 


comm <- fun_winter4snp_hel
dim(comm)
rda_pcnm_fun_winter <- capscale(comm ~ ., data = pcnm_df4rda_winter_fun[,2:101], distance = "bray")

RsquareAdj(rda_pcnm_fun_winter)$adj.r.squared


pcnm_sel_winter_fun <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=pcnm_df4rda_winter_fun[,2:101]),
  scope = formula(rda_pcnm_fun_winter),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


pcoa_snp_fun_winter <- cmdscale(D_snp_fun_winter,k=100,eig=T)
sum(pcoa_snp_fun_winter$eig[1:100])/sum(pcoa_snp_fun_winter$eig)

pcoa_snp_fun_winter_points <- pcoa_snp_fun_winter$points |> 
  as.data.frame() |> 
  rownames_to_column("Entry_lower") |> 
  left_join(meta4snp_winter_fun,join_by(Entry_lower))

snp_df4rda_fun_winter <- meta4snp_winter_fun |> select(sampleid,Plot) |> 
  left_join(pcoa_snp_fun_winter_points)

snp_df4rda_fun_winter[is.na(snp_df4rda_fun_winter$V1),]$Entry_lower |> unique() # check lines should be removed from the analysi
snp_df4rda_fun_winter_clean <- snp_df4rda_fun_winter[!is.na(snp_df4rda_fun_winter$V1),]


comm <- fun_winter4snp_hel[snp_df4rda_fun_winter_clean$sampleid, ]

table(rownames(comm)==snp_df4rda_fun_winter_clean$sampleid)

rda_snp_fun_winter <- capscale(comm ~ ., data =snp_df4rda_fun_winter_clean[,4:103], distance = "bray")

snp_sel_winter_fun <- ordiR2step(
  capscale(comm ~ 1, distance = "bray",data=snp_df4rda_fun_winter_clean[,4:103]),
  scope = formula(rda_snp_fun_winter),
  direction = "forward",
  R2scope = TRUE,
  permutations = 999
)


sel_pcnm_winter_fun <- attr(terms(pcnm_sel_winter_fun), "term.labels")
sel_snp_winter_fun  <- attr(terms(snp_sel_winter_fun),  "term.labels")

id_comm <- rownames(fun_winter4snp_hel)

id_pcnm <- pcnm_df4rda_winter_fun$sampleid
id_snp  <- snp_df4rda_fun_winter_clean$sampleid

id_keep <- Reduce(intersect, list(id_comm, id_pcnm, id_snp))
length(id_keep)
comm_fun_winter <- fun_winter4snp_hel[id_keep, , drop = FALSE]
stopifnot(all(rownames(comm_fun_winter) == id_keep))


tmp_snp <- snp_df4rda_fun_winter |>
  filter(sampleid %in% id_keep) |>
  arrange(match(sampleid, id_keep))

tmp_pcnm <- pcnm_df4rda_winter_fun |>
  filter(sampleid.x %in% id_keep) |>
  arrange(match(sampleid.x, id_keep))

stopifnot(all(tmp_snp$sampleid  == rownames(comm_fun_winter)))
stopifnot(all(tmp_pcnm$sampleid == rownames(comm_fun_winter)))


X_snp_winter_fun  <- tmp_snp  #|> select(all_of(sel_snp_winter_fun))
X_pcnm_winter_fun <- tmp_pcnm #|> select(all_of(sel_pcnm_winter_fun))

vpa_fun_winter <- varpart(vegdist(comm_fun_winter,method="bray"), X_pcnm_winter_fun[,sel_pcnm_winter_fun], X_snp_winter_fun[,sel_snp_winter_fun])
plot(vpa_fun_winter)
vpa_fun_winter 


library(rdacca.hp)


dv <- vegdist(comm_bac_spring,method="bray")
iv <- list(pcnm= X_pcnm_spring_bac[,sel_pcnm_spring_bac],snp=X_snp_spring_bac[,sel_snp_spring_bac])
rdacca_bac_spring <- rdacca.hp(dv,iv,method = "dbRDA",type="adjR2")
rdacca_bac_spring


dv <- vegdist(comm_bac_winter,method="bray")
iv <- list(pcnm= X_pcnm_winter_bac[,sel_pcnm_winter_bac],snp=X_snp_winter_bac[,sel_snp_winter_bac])
rdacca_bac_winter <- rdacca.hp(dv,iv,method = "dbRDA",type="adjR2")
rdacca_bac_winter


dv <- vegdist(comm_fun_spring,method="bray")
iv <- list(pcnm= X_pcnm_spring_fun[,sel_pcnm_spring_fun],snp=X_snp_spring_fun[,sel_snp_spring_fun])
rdacca_fun_spring <- rdacca.hp(dv,iv,method = "dbRDA",type="adjR2")
rdacca_fun_spring


dv <- vegdist(comm_fun_winter,method="bray")
iv <- list(pcnm= X_pcnm_winter_fun[,sel_pcnm_winter_fun],snp=X_snp_winter_fun[,sel_snp_winter_fun])
rdacca_fun_winter <- rdacca.hp(dv,iv,method = "dbRDA",var.par=T,type = "adjR2")
rdacca_fun_winter


# ---- helper: convert ordiR2step $anova table to plotting df ----
prep_cum_df <- function(anova_tbl, group_name, p_col = "Pr(>F)", r2_col = "R2.adj",
                        drop_all = TRUE, p_floor = 1e-6) {

  df <- as.data.frame(anova_tbl)
  df$term <- rownames(df)

  # Drop the "<All variables>" row (full model scope) for the cumulative curve
  if (drop_all) {
    df <- df[df$term != "<All variables>", , drop = FALSE]
  }

  # Ensure numeric
  df[[r2_col]] <- as.numeric(df[[r2_col]])
  df[[p_col]]  <- as.numeric(df[[p_col]])

  # Handle p=0 or NA safely (avoid Inf on -log10)
  # p_floor can be set to 1/(nperm+1) if you know permutations; 1e-6 is a conservative default.
  df$p_plot <- pmax(df[[p_col]], p_floor, na.rm = TRUE)
  df$mlog10p <- -log10(df$p_plot)

  df$n_axes <- seq_len(nrow(df))
  df$group  <- group_name

  df %>% select(group, term, n_axes, R2.adj = all_of(r2_col), p = all_of(p_col), mlog10p)
}

# ---- Example usage ----
# Replace these with your objects:
# pcnm_sel_spring_bac, snp_sel_spring_bac (or *_fun etc.)
df_pcnm <- prep_cum_df(pcnm_sel_spring_bac$anova, "PCNM", p_floor = 1e-6)
df_snp  <- prep_cum_df(snp_sel_spring_bac$anova,  "SNP-PCoA", p_floor = 1e-6)

cum_df <- bind_rows(df_pcnm, df_snp)
# PCNM full model adjR²
full_pcnm <- as.numeric(
  pcnm_sel_spring_bac$anova["<All variables>", "R2.adj"]
)

# SNP-PCoA full model adjR²
full_snp <- as.numeric(
  snp_sel_spring_bac$anova["<All variables>", "R2.adj"]
)

full_df <- data.frame(
  group = c("PCNM", "SNP-PCoA"),
  full_R2 = c(full_pcnm, full_snp)
)

# ---- Plot: two curves + point size = -log10(p) ----
p_cum <- ggplot(cum_df, aes(x = n_axes, y = R2.adj, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = mlog10p), alpha = 0.85) +

  geom_hline(
    data = full_df,
    aes(yintercept = full_R2, color = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +

  scale_size_continuous(name = expression(-log[10](p)),range = c(1, 3) ) +
  labs(
    x = "Number of selected axes (forward selection order)",
    y = "Cumulative adjusted R²",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "inside",legend.position.inside = c(0.9,0.1))+
  scale_color_manual(values=c("#894b00","#2c5600"))
ggsave("cum_curve_dbrda_bac_spring.pdf",p_cum,height = 4,width=4)


# ---- helper: convert ordiR2step $anova table to plotting df ----
prep_cum_df <- function(anova_tbl, group_name, p_col = "Pr(>F)", r2_col = "R2.adj",
                        drop_all = TRUE, p_floor = 1e-6) {
  
  df <- as.data.frame(anova_tbl)
  df$term <- rownames(df)
  
  # Drop the "<All variables>" row (full model scope) for the cumulative curve
  if (drop_all) {
    df <- df[df$term != "<All variables>", , drop = FALSE]
  }
  
  # Ensure numeric
  df[[r2_col]] <- as.numeric(df[[r2_col]])
  df[[p_col]]  <- as.numeric(df[[p_col]])
  
  # Handle p=0 or NA safely (avoid Inf on -log10)
  # p_floor can be set to 1/(nperm+1) if you know permutations; 1e-6 is a conservative default.
  df$p_plot <- pmax(df[[p_col]], p_floor, na.rm = TRUE)
  df$mlog10p <- -log10(df$p_plot)
  
  df$n_axes <- seq_len(nrow(df))
  df$group  <- group_name
  
  df %>% select(group, term, n_axes, R2.adj = all_of(r2_col), p = all_of(p_col), mlog10p)
}

# ---- Example usage ----
# Replace these with your objects:
# pcnm_sel_winter_bac, snp_sel_winter_bac (or *_fun etc.)
df_pcnm <- prep_cum_df(pcnm_sel_winter_bac$anova, "PCNM", p_floor = 1e-6)
df_snp  <- prep_cum_df(snp_sel_winter_bac$anova,  "SNP-PCoA", p_floor = 1e-6)

cum_df <- bind_rows(df_pcnm, df_snp)
# PCNM full model adjR²
full_pcnm <- as.numeric(
  pcnm_sel_winter_bac$anova["<All variables>", "R2.adj"]
)

# SNP-PCoA full model adjR²
full_snp <- as.numeric(
  snp_sel_winter_bac$anova["<All variables>", "R2.adj"]
)

full_df <- data.frame(
  group = c("PCNM", "SNP-PCoA"),
  full_R2 = c(full_pcnm, full_snp)
)

# ---- Plot: two curves + point size = -log10(p) ----
p_cum <- ggplot(cum_df, aes(x = n_axes, y = R2.adj, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = mlog10p), alpha = 0.85) +
  
  geom_hline(
    data = full_df,
    aes(yintercept = full_R2, color = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_size_continuous(name = expression(-log[10](p)),range = c(1, 3) ) +
  labs(
    x = "Number of selected axes (forward selection order)",
    y = "Cumulative adjusted R²",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",legend.position.inside = c(0.9,0.1))+
  scale_color_manual(values=c("#894b00","#2c5600"))
ggsave("cum_curve_dbrda_bac_winter.pdf",p_cum,height = 4,width=4)


# ---- helper: convert ordiR2step $anova table to plotting df ----
prep_cum_df <- function(anova_tbl, group_name, p_col = "Pr(>F)", r2_col = "R2.adj",
                        drop_all = TRUE, p_floor = 1e-6) {
  
  df <- as.data.frame(anova_tbl)
  df$term <- rownames(df)
  
  # Drop the "<All variables>" row (full model scope) for the cumulative curve
  if (drop_all) {
    df <- df[df$term != "<All variables>", , drop = FALSE]
  }
  
  # Ensure numeric
  df[[r2_col]] <- as.numeric(df[[r2_col]])
  df[[p_col]]  <- as.numeric(df[[p_col]])
  
  # Handle p=0 or NA safely (avoid Inf on -log10)
  # p_floor can be set to 1/(nperm+1) if you know permutations; 1e-6 is a conservative default.
  df$p_plot <- pmax(df[[p_col]], p_floor, na.rm = TRUE)
  df$mlog10p <- -log10(df$p_plot)
  
  df$n_axes <- seq_len(nrow(df))
  df$group  <- group_name
  
  df %>% select(group, term, n_axes, R2.adj = all_of(r2_col), p = all_of(p_col), mlog10p)
}

# ---- Example usage ----
# Replace these with your objects:
# pcnm_sel_spring_fun, snp_sel_spring_fun (or *_fun etc.)
df_pcnm <- prep_cum_df(pcnm_sel_spring_fun$anova, "PCNM", p_floor = 1e-6)
df_snp  <- prep_cum_df(snp_sel_spring_fun$anova,  "SNP-PCoA", p_floor = 1e-6)

cum_df <- bind_rows(df_pcnm, df_snp)
# PCNM full model adjR²
full_pcnm <- as.numeric(
  pcnm_sel_spring_fun$anova["<All variables>", "R2.adj"]
)

# SNP-PCoA full model adjR²
full_snp <- as.numeric(
  snp_sel_spring_fun$anova["<All variables>", "R2.adj"]
)

full_df <- data.frame(
  group = c("PCNM", "SNP-PCoA"),
  full_R2 = c(full_pcnm, full_snp)
)

# ---- Plot: two curves + point size = -log10(p) ----
p_cum <- ggplot(cum_df, aes(x = n_axes, y = R2.adj, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = mlog10p), alpha = 0.85) +
  
  geom_hline(
    data = full_df,
    aes(yintercept = full_R2, color = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_size_continuous(name = expression(-log[10](p)),range = c(1, 3) ) +
  labs(
    x = "Number of selected axes (forward selection order)",
    y = "Cumulative adjusted R²",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",legend.position.inside = c(0.9,0.1))+
  scale_color_manual(values=c("#894b00","#2c5600"))
ggsave("cum_curve_dbrda_fun_spring.pdf",p_cum,height = 4,width=4)


# ---- helper: convert ordiR2step $anova table to plotting df ----
prep_cum_df <- function(anova_tbl, group_name, p_col = "Pr(>F)", r2_col = "R2.adj",
                        drop_all = TRUE, p_floor = 1e-6) {
  
  df <- as.data.frame(anova_tbl)
  df$term <- rownames(df)
  
  # Drop the "<All variables>" row (full model scope) for the cumulative curve
  if (drop_all) {
    df <- df[df$term != "<All variables>", , drop = FALSE]
  }
  
  # Ensure numeric
  df[[r2_col]] <- as.numeric(df[[r2_col]])
  df[[p_col]]  <- as.numeric(df[[p_col]])
  
  # Handle p=0 or NA safely (avoid Inf on -log10)
  # p_floor can be set to 1/(nperm+1) if you know permutations; 1e-6 is a conservative default.
  df$p_plot <- pmax(df[[p_col]], p_floor, na.rm = TRUE)
  df$mlog10p <- -log10(df$p_plot)
  
  df$n_axes <- seq_len(nrow(df))
  df$group  <- group_name
  
  df %>% select(group, term, n_axes, R2.adj = all_of(r2_col), p = all_of(p_col), mlog10p)
}

# ---- Example usage ----
# Replace these with your objects:
# pcnm_sel_winter_fun, snp_sel_winter_fun (or *_fun etc.)
df_pcnm <- prep_cum_df(pcnm_sel_winter_fun$anova, "PCNM", p_floor = 1e-6)
df_snp  <- prep_cum_df(snp_sel_winter_fun$anova,  "SNP-PCoA", p_floor = 1e-6)

cum_df <- bind_rows(df_pcnm, df_snp)
# PCNM full model adjR²
full_pcnm <- as.numeric(
  pcnm_sel_winter_fun$anova["<All variables>", "R2.adj"]
)

# SNP-PCoA full model adjR²
full_snp <- as.numeric(
  snp_sel_winter_fun$anova["<All variables>", "R2.adj"]
)

full_df <- data.frame(
  group = c("PCNM", "SNP-PCoA"),
  full_R2 = c(full_pcnm, full_snp)
)

# ---- Plot: two curves + point size = -log10(p) ----
p_cum <- ggplot(cum_df, aes(x = n_axes, y = R2.adj, color = group)) +
  geom_line(linewidth = 1) +
  geom_point(aes(size = mlog10p), alpha = 0.85) +
  
  geom_hline(
    data = full_df,
    aes(yintercept = full_R2, color = group),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_size_continuous(name = expression(-log[10](p)),range = c(1, 3) ) +
  labs(
    x = "Number of selected axes (forward selection order)",
    y = "Cumulative adjusted R²",
    color = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none",legend.position.inside = c(0.9,0.1))+
  scale_color_manual(values=c("#894b00","#2c5600"))
ggsave("cum_curve_dbrda_fun_winter.pdf",p_cum,height = 4,width=4)


save.image("pcnm_VPA.rdata")
