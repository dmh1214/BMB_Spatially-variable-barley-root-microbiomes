##agro heatmap visulize biomass in the fields.

rm(list=ls())
# data----
library(tidyverse)
agb <- read.csv("Fun_Aut/0_raw_data/2023_aboveground_biomass_AUT.csv")
meta <- read.csv("Fun_Aut/0_raw_data/meta_file.csv")

## merge 

agb_spring <- agb |> 
  filter(Trial=="SpringPan") |> 
  pivot_longer(cols=-c(Trial,plot0,Plot) ,names_to = "Side",values_to = "agb")
agb_winter <- agb |> 
  filter(Trial=="WinterPan")|> 
  pivot_longer(cols=-c(Trial,plot0,Plot) ,names_to = "Side",values_to = "agb")

meta_spring <- meta |> filter(Location=="Austria",
                              Trial=="SpringPan")


meta_winter <- meta |> filter(Location=="Austria",
                              Trial=="WinterPan")


df_spring <- agb_spring |> 
  left_join(meta_spring) |> 
  filter(!is.na(agb))
df_spring <- 
df_spring |> mutate(
  Column=case_when(Side=="S"~Col-0.5,
                   Side=="N"~Col)
)







df_spring <- df_spring |>
  mutate(Entry = reorder(Entry, agb, FUN = mean))

p_biomass_entry_spring <- 
ggplot(df_spring, aes(x=Entry,y=agb) )+
         geom_point(alpha=0.5)+
  theme_bw()+  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )+stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "green4",
    linewidth = 0.5
  )+
  ylim(0,100)
ggsave("biomass_rank_spring.pdf",p_biomass_entry_spring,height = 6,width=15,dpi=300)

## anova

mod_spring_aov <- aov(agb~Entry,data=df_spring)
mod_spring_aov |> summary()


df_spring1 <- df_spring |> 
  group_by(Entry) |> 
  summarise(mean_agb=mean(agb),sd_agb=sd(agb))

library(ggpubr)
library(ggrepel)

# 生物量 top10
top_mean  <- df_spring1 |> slice_max(mean_agb, n = 10) |> pull(Entry)

# 方差 top10
top_sd    <- df_spring1 |> slice_max(sd_agb,   n = 10) |> pull(Entry)

# 交集
top_both  <- intersect(top_mean, top_sd)

# 为每个 Entry 赋分类别
df_spring1 <- df_spring1 |>
  mutate(group = case_when(
    Entry %in% top_both ~ "High Mean & SD",     # 交集 → 红色
    Entry %in% top_mean ~ "High Mean(top10)",               # 生物量 top10 → 绿色
    Entry %in% top_sd   ~ "High SD(top10)",                 # 方差 top10 → 黄色
    TRUE ~ "Other"
  ))

cols <- c(
  "High Mean(top10)" = "green4",
  "High SD(top10)"   = "gold",
  "High Mean & SD" = "red4",
  "Other" = "grey80"
)


p_corr_mean_sd_spring <- 
  ggplot(df_spring1, aes(mean_agb, sd_agb, color = group)) +
  geom_point(size = 2.5) +
  
  geom_text_repel(
    data = df_spring1 |> filter(group != "Other"),
    aes(label = Entry, color = group),
    size = 4,
    max.overlaps = Inf
  ) +
  
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  
  # 只用整体数据算相关性
  stat_cor(
    aes(x = mean_agb, y = sd_agb),
    inherit.aes = FALSE,
    label.x =  min(df_spring1$mean_agb),
    label.y =  max(df_spring1$sd_agb)
  ) +
  
  scale_color_manual(values = cols) +
  theme_bw(base_size = 14)+
  labs(x="Mean biomass (g)",y="SD Biomass (g)")+
  theme(legend.position="inside",
        legend.position.inside = c(0.8,0.2))
  

p_corr_mean_sd_spring
ggsave("Corr_bioamss_mean_sd_spring.pdf",p_corr_mean_sd_spring,height = 6,width=6,dpi=300)



df_spring2 <- df_spring |> 
  group_by(Entry,Col,Row,Plot) |> 
  summarise(mean_biomass=mean(agb))

library(RColorBrewer)

ggplot(df_spring2, aes(x=Row, y=Col, fill=mean_biomass)) +
  geom_tile(color="white") +
  geom_text(aes(label=Entry), color="black", size=3) +
  scale_fill_gradientn(
    colours=colorRampPalette(brewer.pal(9, "Blues"))(100),
     name="Biomass", na.value="white") +
  labs(x="Row", y="Column") +
  theme_classic(base_size=15) +
  theme(legend.position="right",
        legend.title= element_text(size=12),
        legend.key.size=unit(0.5,'cm'),
        legend.key=element_rect(color=alpha("white",.05), fill=alpha("white",.05)),
        legend.text=element_text(size=12),
        legend.background=element_rect(fill=alpha("white",.05)),
        panel.grid=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        axis.line=element_line(linewidth=0.5, colour="black"))




p <- agronomymap(
  data       = df_spring,
  map        = c("Row", "Column"),   # 这里仍然用 Column, 能看见 N/S 差异
  variable   = "agb",
  palette    = "YlGn",
  label_title = "Biomass (g/5 plants)",
  point_fill = "grey25",
  add_border = T  # 先关掉自带的 border，方便我们自己加“大格”
)

p

ggsave("biomass_agronomymap_spring.pdf",p,height = 8,width=10,dpi=300)


# 只保留需要的列（可选）
df_spring_clean <- df_spring |>
  mutate(
    Column = case_when(
      Side == "N" ~ Col + 0.25,
      Side == "S" ~ Col - 0.25
    )
  ) |>
  select(Row, Col, Side, Column, agb)


# 每个大格子：一个 Row + 一个 Col
# 大格子坐标
grid_big <- df_spring_clean |>
  dplyr::distinct(Row, Col)

# N/S 分割线坐标
ns_lines <- grid_big |>
  dplyr::mutate(y = Col)
library(ggthemes)

p_biomass_ns <-
  ggplot() +
  # ① 大格子边框：黑色实线
  geom_rect(
    data = grid_big,
    aes(
      xmin = Row - 0.5,
      xmax = Row + 0.5,
      ymin = Col - 0.5,
      ymax = Col + 0.5
    ),
    fill  = NA,
    color = "black",
    linewidth = 0.7
  ) +
  
  # ② 中间 N/S 分割线：灰色虚线
  geom_segment(
    data = ns_lines,
    aes(
      x    = Row - 0.5,
      xend = Row + 0.5,
      y    = y,
      yend = y
    ),
    color    = "grey60",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  # ③ N/S tile，上下半格，不要边框线
  geom_tile(
    data = df_spring_clean,
    aes(
      x    = Row,
      y    = Column,
      fill = agb
    ),
    width  = 0.9,
    height = 0.45,
    color  = NA      # 关键：去掉这里的灰色边框
  ) +
  
  scale_fill_distiller(
    palette   = "YlGn",
    direction = 1,
    name      = "Biomass",
    na.value  = "white"
  ) +
  coord_equal() +
  theme_few(base_size = 14)+
  labs(x="Row",y="Column",title="SpringPan")

p_biomass_ns


ggsave("biomass_ns_spring.pdf",p_biomass_ns,width=10,height = 6,dpi=300)



# winter -------
df_winter <- agb_winter |> 
  left_join(meta_winter) |> 
  filter(!is.na(agb))
df_winter <- 
  df_winter |> mutate(
    Column=case_when(Side=="S"~Col-0.5,
                     Side=="N"~Col)
  )







df_winter <- df_winter |>
  mutate(Entry = reorder(Entry, agb, FUN = mean))

p_biomass_entry_winter <- 
  ggplot(df_winter, aes(x=Entry,y=agb) )+
  geom_point(alpha=0.5)+
  theme_bw()+  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )+stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "green4",
    linewidth = 0.5
  )+
  ylim(0,100)
ggsave("biomass_rank_winter.pdf",p_biomass_entry_winter,height =6,width=15,dpi=300)


mod_winter_aov <- aov(agb~Entry,data=df_winter)
mod_winter_aov |> summary()

df_winter1 <- df_winter |> 
  group_by(Entry) |> 
  summarise(mean_agb=mean(agb),sd_agb=sd(agb))

library(ggpubr)
library(ggrepel)

# 生物量 top10
top_mean  <- df_winter1 |> slice_max(mean_agb, n = 10) |> pull(Entry)

# 方差 top10
top_sd    <- df_winter1 |> slice_max(sd_agb,   n = 10) |> pull(Entry)

# 交集
top_both  <- intersect(top_mean, top_sd)

# 为每个 Entry 赋分类别
df_winter1 <- df_winter1 |>
  mutate(group = case_when(
    Entry %in% top_both ~ "High Mean & SD",     # 交集 → 红色
    Entry %in% top_mean ~ "High Mean(top10)",               # 生物量 top10 → 绿色
    Entry %in% top_sd   ~ "High SD(top10)",                 # 方差 top10 → 黄色
    TRUE ~ "Other"
  ))

cols <- c(
  "High Mean(top10)" = "green4",
  "High SD(top10)"   = "gold",
  "High Mean & SD" = "red4",
  "Other" = "grey80"
)


p_corr_mean_sd_winter <- 
  ggplot(df_winter1, aes(mean_agb, sd_agb, color = group)) +
  geom_point(size = 2.5) +
  
  geom_text_repel(
    data = df_winter1 |> filter(group != "Other"),
    aes(label = Entry, color = group),
    size = 4,
    max.overlaps = Inf
  ) +
  
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  
  # 只用整体数据算相关性
  stat_cor(
    aes(x = mean_agb, y = sd_agb),
    inherit.aes = FALSE,
    label.x =  min(df_winter1$mean_agb),
    label.y =  max(df_winter1$sd_agb)
  ) +
  
  scale_color_manual(values = cols) +
  theme_bw(base_size = 14)+
  labs(x="Mean biomass (g)",y="SD Biomass (g)")+
  theme(legend.position="inside",
        legend.position.inside = c(0.8,0.2))


p_corr_mean_sd_winter
ggsave("Corr_bioamss_mean_sd_winter.pdf",p_corr_mean_sd_winter,height = 6,width=6,dpi=300)



df_winter2 <- df_winter |> 
  group_by(Entry,Col,Row,Plot) |> 
  summarise(mean_biomass=mean(agb))

library(RColorBrewer)

ggplot(df_winter2, aes(x=Row, y=Col, fill=mean_biomass)) +
  geom_tile(color="white") +
  geom_text(aes(label=Entry), color="black", size=3) +
  scale_fill_gradientn(
    colours=colorRampPalette(brewer.pal(9, "Blues"))(100),
    name="Biomass", na.value="white") +
  labs(x="Row", y="Column") +
  theme_classic(base_size=15) +
  theme(legend.position="right",
        legend.title= element_text(size=12),
        legend.key.size=unit(0.5,'cm'),
        legend.key=element_rect(color=alpha("white",.05), fill=alpha("white",.05)),
        legend.text=element_text(size=12),
        legend.background=element_rect(fill=alpha("white",.05)),
        panel.grid=element_blank(),
        panel.border=element_rect(color="black", fill=NA, linewidth=0.5),
        axis.line=element_line(linewidth=0.5, colour="black"))




p <- agronomymap(
  data       = df_winter,
  map        = c("Row", "Column"),   # 这里仍然用 Column, 能看见 N/S 差异
  variable   = "agb",
  palette    = "YlGn",
  label_title = "Biomass (g/5 plants)",
  point_fill = "grey25",
  add_border = T  # 先关掉自带的 border，方便我们自己加“大格”
)

p
ggsave("biomass_agronomymap_winter.pdf",p,width = 10,height=8,dpi=300)



# 只保留需要的列（可选）
df_winter_clean <- df_winter |>
  mutate(
    Column = case_when(
      Side == "N" ~ Col + 0.25,
      Side == "S" ~ Col - 0.25
    )
  ) |>
  select(Row, Col, Side, Column, agb)


# 每个大格子：一个 Row + 一个 Col
# 大格子坐标
grid_big <- df_winter_clean |>
  dplyr::distinct(Row, Col)

# N/S 分割线坐标
ns_lines <- grid_big |>
  dplyr::mutate(y = Col)
library(ggthemes)

p_biomass_ns <-
  ggplot() +
  # ① 大格子边框：黑色实线
  geom_rect(
    data = grid_big,
    aes(
      xmin = Row - 0.5,
      xmax = Row + 0.5,
      ymin = Col - 0.5,
      ymax = Col + 0.5
    ),
    fill  = NA,
    color = "black",
    linewidth = 0.7
  ) +
  
  # ② 中间 N/S 分割线：灰色虚线
  geom_segment(
    data = ns_lines,
    aes(
      x    = Row - 0.5,
      xend = Row + 0.5,
      y    = y,
      yend = y
    ),
    color    = "grey60",
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  # ③ N/S tile，上下半格，不要边框线
  geom_tile(
    data = df_winter_clean,
    aes(
      x    = Row,
      y    = Column,
      fill = agb
    ),
    width  = 0.9,
    height = 0.45,
    color  = NA      # 关键：去掉这里的灰色边框
  ) +
  
  scale_fill_distiller(
    palette   = "YlGn",
    direction = 1,
    name      = "Biomass",
    na.value  = "white"
  ) +
  coord_equal() +
  theme_few(base_size = 14)+
  labs(x="Row",y="Column",title="winterPan")

p_biomass_ns


ggsave("biomass_ns_winter.pdf",p_biomass_ns,width=10,height = 6,dpi=300)

