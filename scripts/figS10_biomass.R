# FigS10: exact blocks from Biomass_heatmap.R, lines 1-57 and 254-286.
# The unrelated agronomymap sections remain in the archived original source.
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
