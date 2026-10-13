# Figures: Fig1a,c; FigS11a,b
# Original source: Bac_Aut/beta.Rmd
# Analysis statements retained; input/output locations are handled by run_figure.R.

library(tidyverse)
library(vegan)
library(ggthemes)
library(ggpubr)


#raw data
rm(list=ls())
meta_all <- readRDS("1_processed_data/meta_aut.RDS")
meta_spring <- readRDS("1_processed_data/meta_spring.rds")
meta_winter <- readRDS("1_processed_data/meta_winter.rds")

asv_all_raw_c <- readRDS("1_processed_data/asv_count_filtered.RDS")
asv_spring_raw_c <- readRDS("1_processed_data/asv_spring.rds")
asv_winter_raw_c <- readRDS("1_processed_data/asv_winter.rds")

asv_all_rclr <- readRDS("1_processed_data/asv_rclr.RDS")
asv_spring_rclr <- readRDS("1_processed_data/asv_spring_rclr.rds")
asv_winter_rclr <- readRDS("1_processed_data/asv_winter_rclr.rds")

if (!dir.exists("beta")) {
  dir.create("beta")
} else {
  print("dir already exist")
}


table(meta_all$sampleid==colnames(asv_all_raw_c))
table(meta_all$sampleid==rownames(asv_all_rclr))

table(meta_spring$sampleid==colnames(asv_spring_raw_c))
table(meta_spring$sampleid==rownames(asv_spring_rclr))

table(meta_winter$sampleid==colnames(asv_winter_raw_c))
table(meta_winter$sampleid==rownames(asv_winter_rclr))


meta_all <- 
meta_all |> 
  mutate(check=if_else(Group=="Check","Check","Unique"))

pca_all_AD <- rda(asv_all_rclr)

pca_all_AD_eigvals <- summary(pca_all_AD)$cont$importance

pca_all_AD_loadings <- data.frame(
  PC1 = round(pca_all_AD_eigvals[2, 1] * 100, 2),
  PC2 = round(pca_all_AD_eigvals[2, 2] * 100, 2)
  )

pca_all_AD_scores <-  scores(pca_all_AD, display = "sites", scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(meta_all)


pca_all_AD_scores |> dim()
p_pca_all <- 
ggplot(pca_all_AD_scores, aes(x = PC1, y = PC2, fill=Trial,shape=check)) +
  geom_point(size=2,alpha=0.5) +
  theme_bw() +
  labs(
    x = paste0("PC1 (", pca_all_AD_loadings$PC1, "%)"),
    y = paste0("PC2 (", pca_all_AD_loadings$PC2, "%)"),
    title = "PCA based on Aitchison distance (n=2323)"
  )+
  scale_shape_manual(values=c(23,21))+
  theme(legend.position = "inside",legend.position.inside = c(0.85,0.3))

ggsave("PCA_Aitchison_all.pdf",p_pca_all,height = 6,width=6,dpi=300)


meta_all_checkLine <- meta_all |> filter(check=="Check")
asv_raw_count_check <- asv_all_raw_c[,meta_all_checkLine$sampleid]
asv_raw_count_check <- asv_raw_count_check [rowSums(asv_raw_count_check>0)>5,]


X <- t(as.matrix(asv_raw_count_check))

rn <- rownames(X)  # sample names
cn <- colnames(X)  # ASV ids

asv_check_rclr <- vegan::decostand(X, method = "rclr")
# re-attach names
rownames(asv_check_rclr) <- rn
colnames(asv_check_rclr) <- cn




pca_check_AD <- rda(asv_check_rclr)

pca_check_AD_eigvals <- summary(pca_check_AD)$cont$importance

pca_check_AD_loadings <- data.frame(
  PC1 = round(pca_check_AD_eigvals[2, 1] * 100, 2),
  PC2 = round(pca_check_AD_eigvals[2, 2] * 100, 2)
  )

pca_check_AD_scores <-  scores(pca_check_AD, display = "sites", scaling = 2) |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(meta_all_checkLine)

# === 计算每个组的几何中心（group centroid） ===
# Compute geometric centers (centroids) for each group (Entry)
centroids <- pca_check_AD_scores %>%
  group_by(Entry,Trial) %>%
  summarise(
    PC1_mean = mean(PC1, na.rm = TRUE),  # 每组PC1均值 / mean PC1 for each group
    PC2_mean = mean(PC2, na.rm = TRUE)   # 每组PC2均值 / mean PC2 for each group
  )

pca_check_AD_scores <- pca_check_AD_scores |> left_join(centroids)

# === 绘制PCA散点图并添加连线 ===
# Draw PCA scatter plot and add centroid connection lines
p_pca_check <- ggplot(pca_check_AD_scores, aes(x = PC1, y = PC2, shape=Trial,fill=Entry))+
  geom_point(aes(color=Entry),size = 2, alpha = 0.5) +  # 绘制样本点 / plot sample points
  
  # 添加各组几何中心之间的连线
  # Add connection lines between group centroids
 geom_segment( aes(x = PC1, y = PC2,
               xend = PC1_mean, yend = PC2_mean, color = Entry),
                      alpha = 0.6, show.legend = FALSE)+
  
  
  # 图形主题设置
  # Set theme and labels
  theme_bw() +
  labs(
    x = paste0("PC1 (", pca_check_AD_loadings$PC1, "%)"),   # X轴标签 / X-axis label
    y = paste0("PC2 (", pca_check_AD_loadings$PC2, "%)"),   # Y轴标签 / Y-axis label
    title = "PCA based on Aitchison distance (Check Lines)"       # 图标题 / Plot title
  ) +
  theme(
    legend.position = "inside",                # 将图例放入图内 / place legend inside plot
    legend.position.inside = c(0.85, 0.5)      # 图例位置坐标 / legend coordinates inside plot
  )+
  scale_shape_manual(values=c(22,25))


# === 显示图形 / Show the plot ===
p_pca_check

adonis.check <- adonis2(asv_check_rclr~Trial*Entry,data=meta_all_checkLine,method="eu",by="terms")
adonis.check

D_aitchison_check <- vegdist(asv_check_rclr,"eu")

cap.check.ait <- capscale(D_aitchison_check~Entry*Trial,data=meta_all_checkLine,add=T)


summary(cap.check.ait)
anova.cca(cap.check.ait, permutations = 999, by = "terms")

ggsave("PCA_Aitchison_checkLine.pdf",p_pca_check,height = 6,width=6,dpi=300)


D_bc_check <- asv_raw_count_check |>
  t() |> 
  decostand("hellinger") |> 
  vegdist(method="bray")

adonis.check.bc <- adonis2(D_bc_check~ Entry*Trial,data=meta_all_checkLine,by="terms")
adonis.check.bc


pcoa_check_bray <- cmdscale(D_bc_check,k=2,eig=T)

pca_check_bray_p_data <- pcoa_check_bray$points |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(meta_all_checkLine)

pca_check_bray_loadings <- c(round((pcoa_check_bray$eig[1]*100/sum(pcoa_check_bray$eig)),2),round((pcoa_check_bray$eig[2]*100/sum(pcoa_check_bray$eig)),2))


# === 计算每个组的几何中心（group centroid） ===
# Compute geometric centers (centroids) for each group (Entry)
centroids <- pca_check_bray_p_data %>%
  group_by(Entry,Trial) %>%
  summarise(
    V1_mean = mean(V1, na.rm = TRUE),  # 每组PC1均值 / mean PC1 for each group
    V2_mean = mean(V2, na.rm = TRUE)   # 每组PC2均值 / mean PC2 for each group
  )

pca_check_bray_p_data <-pca_check_bray_p_data |> left_join(centroids)

# === 绘制PCA散点图并添加连线 ===
# Draw PCA scatter plot and add centroid connection lines
p_pca_check <- ggplot(pca_check_bray_p_data, aes(x = V1, y = V2, shape=Trial,fill=Entry))+
  geom_point(aes(color=Entry),size = 2, alpha = 0.5) +  # 绘制样本点 / plot sample points
  
  # 添加各组几何中心之间的连线
  # Add connection lines between group centroids
 geom_segment( aes(x = V1, y = V2,
               xend = V1_mean, yend = V2_mean, color = Entry),
                      alpha = 0.6, show.legend = FALSE)+
  
  
  # 图形主题设置
  # Set theme and labels
  theme_bw() +
  labs(
    x = paste0("PCo1 (", pca_check_bray_loadings[1], "%)"),   # X轴标签 / X-axis label
    y = paste0("PCo2 (", pca_check_bray_loadings[2], "%)"),   # Y轴标签 / Y-axis label
    title = "PCoA of Check Lines"       # 图标题 / Plot title
  ) +
  theme(
    legend.position = "top",                # 将图例放入图内 / place legend inside plot
  )+
  scale_shape_manual(values=c(22,25))+
  theme(legend.position="none")


# === 显示图形 / Show the plot ===
p_pca_check
ggsave("beta/pcoa_checkLine.pdf",p_pca_check,height = 6,width = 6,dpi=300)


asv_all_raw_c_f <- asv_all_raw_c[rowSums(asv_all_raw_c>0)>50,]
D_bc_all <- asv_all_raw_c_f |>
  t() |> 
  decostand("hellinger") |> 
  vegdist(method="bray")
saveRDS(D_bc_all,"1_processed_data/D_bray_all_samples")

#adonis.all.bc <- adonis2(D_bc_all~ Entry*Trial,data=meta_all_allLine,by="terms")
#adonis.all.bc


pcoa_all_bray <- cmdscale(D_bc_all,k=2,eig=T)

pca_all_bray_p_data <- pcoa_all_bray$points |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(meta_all)

pca_all_bray_loadings <- c(round((pcoa_all_bray$eig[1]*100/sum(pcoa_all_bray$eig)),2),round((pcoa_all_bray$eig[2]*100/sum(pcoa_all_bray$eig)),2))


# === 计算每个组的几何中心（group centroid） ===
# Compute geometric centers (centroids) for each group (Entry)
centroids <- pca_all_bray_p_data %>%
  group_by(Entry,Trial) %>%
  summarise(
    V1_mean = mean(V1, na.rm = TRUE),  # 每组PC1均值 / mean PC1 for each group
    V2_mean = mean(V2, na.rm = TRUE)   # 每组PC2均值 / mean PC2 for each group
  )

pca_all_bray_p_data <-pca_all_bray_p_data |> left_join(centroids)

# === 绘制PCA散点图并添加连线 ===
# Draw PCA scatter plot and add centroid connection lines
p_pca_all <- ggplot(pca_all_bray_p_data, aes(x = V1, y = V2, fill=Trial))+
  geom_hline(yintercept=0,color="grey20",linetype="dashed")+
  geom_vline(xintercept=0,color="grey20",linetype="dashed")+
  geom_point(color="grey20",size = 2, alpha = 0.5,aes(shape=check)) +  
  theme_few() +
  labs(
    x = paste0("PCo1 (", pca_all_bray_loadings[1], "%)"),
    y = paste0("PCo2 (", pca_all_bray_loadings[2], "%)"))+    
  theme(
    legend.position = "inside",legend.position.inside = c(0.1,0.85)                
  )+
  scale_shape_manual(values=c(22,21))+
  scale_fill_manual(values=c("#458a74","#cc5b45"))



# === 显示图形 / Show the plot ===
p_pca_all


ggsave("beta/pcoa_bray_all.pdf",p_pca_all,height=6,width=6,dpi=300)


### permanove of panel

adoni2_panel <- adonis2(D_bc_all~Trial, data=pca_all_bray_p_data)
adoni2_panel


asv_spring <- readRDS("1_processed_data/asv_spring.rds")

table(colnames(asv_spring)==meta_spring$sampleid)

D_bc_spring <- asv_spring |>
  t() |> 
  decostand("hellinger") |> 
  vegdist(method="bray")

saveRDS(D_bc_spring,"1_processed_data/D_bray_spring.rds")

pcoa_spring_bray <- cmdscale(D_bc_spring,k=20,eig=T)

field_spring <- readRDS("1_processed_data/field_spring.rds")

field_spring_s <- field_spring |> 
  select(Plot,Biomass.t_ha.,GY.t_ha.,S,N,mean_SN) |> 
   pivot_longer(
    cols = c(S, N),
    names_to = "Side",     # 新列名，表示元素类型（S 或 N）
    values_to = "agb_side"       # 新列名，表示对应的数值
  )
  
spring_meta_field <- meta_spring |> left_join(field_spring_s)

pca_spring_bray_p_data <- pcoa_spring_bray$points |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(spring_meta_field)

pca_spring_bray_loadings <- round((pcoa_spring_bray$eig*100/sum(pcoa_spring_bray$eig)),2)

# === 计算每个组的几何中心（group centroid） ===
# Compute geometric centers (centroids) for each group (Entry)
centroids <- pca_spring_bray_p_data %>%
  group_by(Entry,Trial) %>%
  summarise(
    V1_mean = mean(V1, na.rm = TRUE),  # 每组PC1均值 / mean PC1 for each group
    V2_mean = mean(V2, na.rm = TRUE)   # 每组PC2均值 / mean PC2 for each group
  )

pca_spring_bray_p_data <-pca_spring_bray_p_data |> left_join(centroids)




p_pca_spring <- ggplot(pca_spring_bray_p_data, aes(x = V1, y = V2,fill=mean_SN))+
  geom_hline(yintercept=0,color="grey20",linetype="dashed")+
  geom_vline(xintercept=0,color="grey20",linetype="dashed")+
  geom_point(aes(size=mean_SN),color="grey20", alpha = 0.5,shape=21) +  
  theme_few() +
  labs(
    x = paste0("PCo1 (", pca_spring_bray_loadings[1], "%)"),
    y = paste0("PCo2 (", pca_spring_bray_loadings[2], "%)"))+    
  theme(
    legend.position = "inside",legend.position.inside = c(0.85,0.85)                
  )+
   scale_fill_gradient(low="#6a8ec9", high="#f5a216")
p_pca_spring 

ggsave("beta/pcoa_biomass_spring.pdf",p_pca_spring,height=6,width = 6, dpi=300 )

pcor_springV1 <-  
  ggplot(pca_spring_bray_p_data, aes(x = V1, y = mean_SN))+
  geom_point(alpha=0.6)+
  geom_smooth(method="lm",color="grey3")+
  stat_cor(method="pearson")+
 labs(
  x = "PCo1",
  y= "Barley biomass (g/5 plants)")+
  theme_few()

pcor_springV1

ggsave("beta/cor_pco1_biomass_spring.pdf",pcor_springV1,height=2,width = 6, dpi=300 ) 

pcor_springV2 <-  
  ggplot(pca_spring_bray_p_data, aes(x = V2, y = mean_SN))+
  geom_point(alpha=0.6)+
  geom_smooth(method="lm",color="grey3")+
  stat_cor(method="pearson")+
 labs(
  x = "PCo2",
  y= "Barley biomass (g/5 plants)")+
  theme_few()

pcor_springV2

ggsave("beta/cor_pco2_biomass_spring.pdf",pcor_springV2,height=2,width = 6, dpi=300 ) 


field.spring <- readRDS("1_processed_data/field_spring.rds")
field.spring1 <- field.spring |> filter(!is.na(GY.t_ha.))
field.spring1 |> dim()

field.spring1$green_biomass <- field.spring1$Biomass.t_ha.-field.spring1$GY.t_ha.

cor.test(~green_biomass+mean_SN,data=field.spring1)
cor.test(~green_biomass+GY.t_ha.,data=field.spring1)
cor.test(~mean_SN+GY.t_ha.,data=field.spring1)

## correlated, but not strong. 
### keep it simple, just continue with above_ground biomass

#total biomass
## check sample with NAs
is.na(pca_spring_bray_p_data$Biomass.t_ha.) |> table()


adonis_spring_tBiomass <- adonis2(D_bc_spring~Biomass.t_ha.,data=pca_spring_bray_p_data,by="terms")
adonis_spring_tBiomass

# above ground green biomass
## check sample with NAs
is.na(pca_spring_bray_p_data$mean_SN) |> table()

## remove the NA sample and re calculate the D
NAsample <- pca_spring_bray_p_data[is.na(pca_spring_bray_p_data$mean_SN),]$sampleid

asv_spring_agb <- select(asv_spring,-all_of(NAsample))
asv_spring_agb |> ncol();pca_spring_bray_p_data |> nrow()

D_bc_spring_agb <- asv_spring_agb |> 
  t() |> 
  decostand("hellinger") |> 
  vegdist("bray")

pca_spring_bray_p_dat_agb <- pca_spring_bray_p_data |> filter(!is.na(mean_SN))
table(colnames(asv_spring_agb)==pca_spring_bray_p_dat_agb$sampleid)


adonis_spring_agb <- adonis2(D_bc_spring_agb~mean_SN,data=pca_spring_bray_p_dat_agb,by="terms")
adonis_spring_agb


asv_winter <- readRDS("1_processed_data/asv_winter.rds")

table(colnames(asv_winter)==meta_winter$sampleid)

D_bc_winter <- asv_winter |>
  t() |> 
  decostand("hellinger") |> 
  vegdist(method="bray")

saveRDS(D_bc_winter,"1_processed_data/D_bray_winter.rds")

pcoa_winter_bray <- cmdscale(D_bc_winter,k=20,eig=T)

field_winter <- readRDS("1_processed_data/field_winter.rds")

field_winter_s <- field_winter |> 
  select(Plot,Biomass.t_ha.,GY.t_ha.,S,N,mean_SN) |> 
   pivot_longer(
    cols = c(S, N),
    names_to = "Side",     # 新列名，表示元素类型（S 或 N）
    values_to = "agb_side"       # 新列名，表示对应的数值
  )
  
winter_meta_field <- meta_winter |> left_join(field_winter_s)

pca_winter_bray_p_data <- pcoa_winter_bray$points |> 
  as.data.frame() |> 
  rownames_to_column("sampleid") |> 
  left_join(winter_meta_field)

pca_winter_bray_loadings <- round((pcoa_winter_bray$eig*100/sum(pcoa_winter_bray$eig)),2)

# === 计算每个组的几何中心（group centroid） ===
# Compute geometric centers (centroids) for each group (Entry)
centroids <- pca_winter_bray_p_data %>%
  group_by(Entry,Trial) %>%
  summarise(
    V1_mean = mean(V1, na.rm = TRUE),  # 每组PC1均值 / mean PC1 for each group
    V2_mean = mean(V2, na.rm = TRUE)   # 每组PC2均值 / mean PC2 for each group
  )

pca_winter_bray_p_data <-pca_winter_bray_p_data |> left_join(centroids)




p_pca_winter <- ggplot(pca_winter_bray_p_data, aes(x = V1, y = V2,fill=mean_SN))+
  geom_hline(yintercept=0,color="grey20",linetype="dashed")+
  geom_vline(xintercept=0,color="grey20",linetype="dashed")+
  geom_point(aes(size=mean_SN),color="grey20", alpha = 0.5,shape=21) +  
  theme_few() +
  labs(
    x = paste0("PCo1 (", pca_winter_bray_loadings[1], "%)"),
    y = paste0("PCo2 (", pca_winter_bray_loadings[2], "%)"))+    
  theme(
    legend.position = "inside",legend.position.inside = c(0.85,0.85)                
  )+
   scale_fill_gradient(low="#6a8ec9", high="#f5a216")
p_pca_winter 

ggsave("beta/pcoa_biomass_winter.pdf",p_pca_winter,height=6,width = 6, dpi=300 )

cor.test(pca_winter_bray_p_data$V1,pca_winter_bray_p_data$mean_SN)
cor.test(pca_winter_bray_p_data$V2,pca_winter_bray_p_data$mean_SN)

pcor_winterV1 <-  
  ggplot(pca_winter_bray_p_data, aes(x = V1, y = mean_SN))+
  geom_point(alpha=0.6)+
  geom_smooth(method="lm",color="grey3")+
  stat_cor(method="pearson")+
 labs(
  x = "PCo1",
  y= "Barley biomass (g/5 plants)")+
  theme_few()

pcor_winterV1

ggsave("beta/cor_pco1_biomass_winter.pdf",pcor_winterV1,height=2,width = 6, dpi=300 ) 

pcor_winterV2 <-  
  ggplot(pca_winter_bray_p_data, aes(x = V2, y = mean_SN))+
  geom_point(alpha=0.6)+
  geom_smooth(method="lm",color="grey3")+
  stat_cor(method="pearson")+
 labs(
  x = "PCo2",
  y= "Barley biomass (g/5 plants)")+
  theme_few()

pcor_winterV2

ggsave("beta/cor_pco2_biomass_winter.pdf",pcor_winterV2,height=2,width = 6, dpi=300 ) 


field.winter1 <- field_winter |> filter(!is.na(GY.t_ha.))
field.winter1 |> dim()

field.winter1$green_biomass <- field.winter1$Biomass.t_ha.-field.winter1$GY.t_ha.

cor.test(~green_biomass+mean_SN,data=field.winter1)
cor.test(~green_biomass+GY.t_ha.,data=field.winter1)
cor.test(~mean_SN+GY.t_ha.,data=field.winter1)

## correlated, but not strong. 
### keep it simple, just continue with above_ground biomass

#total biomass
## check sample with NAs
is.na(pca_winter_bray_p_data$Biomass.t_ha.) |> table()


adonis_winter_tBiomass <- adonis2(D_bc_winter~Biomass.t_ha.,data=pca_winter_bray_p_data,by="terms")
adonis_winter_tBiomass

adonis_winter_gy <- adonis2(D_bc_winter~GY.t_ha.,data=pca_winter_bray_p_data,by="terms")
adonis_winter_gy


# above ground green biomass
## check sample with NAs
is.na(pca_winter_bray_p_data$mean_SN) |> table()

## remove the NA sample and re calculate the D
NAsample <- pca_winter_bray_p_data[is.na(pca_winter_bray_p_data$mean_SN),]$sampleid

asv_winter_agb <- select(asv_winter,-all_of(NAsample))
asv_winter_agb |> ncol();pca_winter_bray_p_data |> nrow()

D_bc_winter_agb <- asv_winter_agb |> 
  t() |> 
  decostand("hellinger") |> 
  vegdist("bray")

pca_winter_bray_p_dat_agb <- pca_winter_bray_p_data |> filter(!is.na(mean_SN))
table(colnames(asv_winter_agb)==pca_winter_bray_p_dat_agb$sampleid)


adonis_winter_agb <- adonis2(D_bc_winter_agb~mean_SN,data=pca_winter_bray_p_dat_agb,by="terms")
adonis_winter_agb


beta_disper_spring <- betadisper(D_bc_spring,group=pca_spring_bray_p_data$Entry)

beta_disper_dat_spring <- beta_disper_spring$distances |> 
  as.data.frame() |> 
  rownames_to_column("sampleid")
colnames(beta_disper_dat_spring)[2] <- "betadisper"

beta_disper_dat_spring <- 
beta_disper_dat_spring |> 
  left_join(pca_spring_bray_p_data) |> 
  group_by(Entry) |> 
  summarise(tBiomass=mean(Biomass.t_ha.),
            tGY=mean(GY.t_ha.),
            tAGB=mean(mean_SN),
            sdAGB=sd(mean_SN),
            mean_disper=mean(betadisper)
            ) |> 
  as.data.frame()

shapiro.test(beta_disper_dat_spring$mean_disper)
shapiro.test(beta_disper_dat_spring$tAGB)

cor.test(~tAGB+mean_disper,data=beta_disper_dat_spring)
cor.test(~sdAGB+mean_disper,data=beta_disper_dat_spring)

p_disper_cor_biomass_spring <-beta_disper_dat_spring |> 
  ggplot(aes(mean_disper,tAGB))+
  geom_smooth(method="lm")+
  stat_cor()+
  geom_point()+
  theme_bw()+
  labs(x="Variation in bacterial community composition (beta dispersion) per genotype",
       y="Average biomass of barley per genotype")
p_disper_cor_biomass_spring 


## winter ##
beta_disper_winter <- betadisper(D_bc_winter,group=pca_winter_bray_p_data$Entry)

beta_disper_dat_winter <- beta_disper_winter$distances |> 
  as.data.frame() |> 
  rownames_to_column("sampleid")
colnames(beta_disper_dat_winter)[2] <- "betadisper"

beta_disper_dat_winter <- 
beta_disper_dat_winter |> 
  left_join(pca_winter_bray_p_data) |> 
  group_by(Entry) |> 
  summarise(tBiomass=mean(Biomass.t_ha.),
            tGY=mean(GY.t_ha.),
            tAGB=mean(mean_SN),
            mean_disper=mean(betadisper)
            ) |> 
  as.data.frame()

shapiro.test(beta_disper_dat_winter$mean_disper)
shapiro.test(beta_disper_dat_winter$tAGB)

cor.test(~tAGB+mean_disper,data=beta_disper_dat_winter)
cor.test(~tBiomass+mean_disper,data=beta_disper_dat_winter)


p_disper_cor_biomass_winter <-beta_disper_dat_winter |> 
  ggplot(aes(mean_disper,tAGB))+
  geom_smooth(method="lm")+
  stat_cor()+
  geom_point()+
  theme_bw()+
  labs(x="Variation in bacterial community composition (beta dispersion) per genotype",
       y="Average biomass of barley per genotype")
p_disper_cor_biomass_winter 


library(patchwork)

p_disper_cor_biomass_spring +p_disper_cor_biomass_winter


save.image("beta/beta.Rdata")
