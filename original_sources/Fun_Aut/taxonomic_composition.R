## ============================================================================
## Taxonomic composition at family level for Spring and Winter panels
## 春季与冬季群体在科水平上的群落组成可视化（按基因型汇总）
## 要求：data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABIAAAASCAYAAABWzo5XAAAAbElEQVR4Xs2RQQrAMAgEfZgf7W9LAguybljJpR3wEse5JOL3ZObDb4x1loDhHbBOFU6i2Ddnw2KNiXcdAXygJlwE8OFVBHDgKrLgSInN4WMe9iXiqIVsTMjH7z/GhNTEibOxQswcYIWYOR/zAjBJfiXh3jZ6AAAAAElFTkSuQmCC
## - Spring 和 Winter 各自取 top15 family，然后取并集（最多 30 个）
## - 其余 family 归为 Others（灰色）
## - 两个季节使用相同的 family 顺序和颜色
## - 堆积条中：丰度最大的 family 在最底层，Others 在最顶层
## ============================================================================

## 清空环境 / Clear workspace -------------------------------------------------
rm(list = ls())

## 载入数据 / Load data -------------------------------------------------------
otu_all    <- readRDS("1_processed_data/otu_count_filtered.RDS")
otu_spring <- readRDS("1_processed_data/otu_spring.rds")
otu_winter <- readRDS("1_processed_data/otu_winter.rds")

meta_spring <- readRDS("1_processed_data/meta_spring.rds")
meta_winter <- readRDS("1_processed_data/meta_winter.rds")

taxa <- readRDS("0_raw_data/Fungal_taxonomy.rds")  # 应包含 otuid, order, family 等列
colnames(taxa)

## 加载 R 包 / Load required packages -----------------------------------------
library(vegan)       # decostand for TSS
library(tidyverse)   # dplyr, tidyr, ggplot2, etc.
library(ggthemes)    # theme_few
library(RColorBrewer)
library(patchwork)   # optional, for combining plots / 可选，用于拼图

## ----------------------------------------------------------------------------
## 1. 将 otu 计数转换为相对丰度（TSS）
##    Transform otu counts to relative abundance (Total Sum Scaling)
## ----------------------------------------------------------------------------
## 注意：decostand("total") 假定“行 = 样本，列 = OTU/otu”
## Note : decostand("total") assumes rows = samples, cols = OTUs/otus

spring_tss <- otu_spring |>
  t() |>
  decostand("total")

winter_tss <- otu_winter |>
  t() |>
  decostand("total")

## ----------------------------------------------------------------------------
## 2. 定义函数：聚合到 family_clean 水平
##    Define helper: aggregate otus to family_clean
##
## 处理逻辑 / Logic:
## - 如果 family 有明确名称且不是 "f_uncultured"：直接使用该科名
## - 如果 family == "f_uncultured" 且有 order：使用 order_uncultured 拆分不同目
## - 其余情况：归为 "Unclassified"
## ----------------------------------------------------------------------------

aggregate_to_family_clean <- function(tss_mat, taxa_tbl) {
  tss_mat |>
    t() |>
    as.data.frame() |>
    rownames_to_column("otuid") |>
    left_join(taxa_tbl |> select(otuid, order, family), by = "otuid") |>
    mutate(
      family_clean = case_when(
        !is.na(family) & family != "Unclassified" ~ family,
        family == "Unclassified" & !is.na(order)  ~ paste0(order, "_Unclassified"),
        TRUE ~ "Unclassified"
      )
    ) |>
    group_by(family_clean) |>
    summarise(across(where(is.numeric), sum), .groups = "drop")
}

## 对 Spring 和 Winter 进行科水平聚合 / Aggregate Spring & Winter to family_clean
spring_family <- aggregate_to_family_clean(spring_tss, taxa)
winter_family <- aggregate_to_family_clean(winter_tss, taxa)

## ----------------------------------------------------------------------------
## 3. 各自取 top15 family，然后取并集
##    Take top15 families in Spring and Winter separately, then union
## ----------------------------------------------------------------------------

tax_level <- "family_clean"
topN <- 15

## Spring: 计算每个 family 的总丰度，并取 top15
## Spring: total abundance per family and select top15
family_totals_spring <- spring_family |>
  filter(.data[[tax_level]] != "Unclassified") |>  
  mutate(total = rowSums(across(where(is.numeric)))) |>
  arrange(desc(total))

top_spring <- head(family_totals_spring[[tax_level]], topN)

## Winter: 同样操作 / Same for Winter
family_totals_winter <- winter_family |>
  filter(.data[[tax_level]] != "Unclassified") |>  
  mutate(total = rowSums(across(where(is.numeric)))) |>
  arrange(desc(total))

top_winter <- head(family_totals_winter[[tax_level]], topN)

## Spring 和 Winter 的 top15 的并集（最多 30 个）
## Union of Spring/Winter top15 (up to 30 families)
union_top <- union(top_spring, top_winter)

## ----------------------------------------------------------------------------
## 4. 展开为长表，并根据 union_top 标记 Taxon_plot（union_top vs Others）
##    Long format and label Taxon_plot (union_top vs Others)
## ----------------------------------------------------------------------------

spring_family_long <- spring_family |>
  pivot_longer(
    cols      = -all_of(tax_level),
    names_to  = "sampleid",
    values_to = "Abundance"
  ) |>
  mutate(
    Taxon_plot = if_else(
      !!sym(tax_level) %in% union_top,
      !!sym(tax_level),
      "Others"
    )
  ) |>
  group_by(sampleid, Taxon_plot) |>
  summarise(Abundance = sum(Abundance), .groups = "drop") |>
  left_join(meta_spring, by = "sampleid")

winter_family_long <- winter_family |>
  pivot_longer(
    cols      = -all_of(tax_level),
    names_to  = "sampleid",
    values_to = "Abundance"
  ) |>
  mutate(
    Taxon_plot = if_else(
      !!sym(tax_level) %in% union_top,
      !!sym(tax_level),
      "Others"
    )
  ) |>
  group_by(sampleid, Taxon_plot) |>
  summarise(Abundance = sum(Abundance), .groups = "drop") |>
  left_join(meta_winter, by = "sampleid")

## ----------------------------------------------------------------------------
## 5. 按 Entry 求平均相对丰度（一个基因型一条柱子）
##    Average relative abundance per Entry (one bar per genotype)
## ----------------------------------------------------------------------------

spring_family_entry <- spring_family_long |>
  group_by(Entry, Taxon_plot) |>
  summarise(mean_Abundance = mean(Abundance), .groups = "drop")

winter_family_entry <- winter_family_long |>
  group_by(Entry, Taxon_plot) |>
  summarise(mean_Abundance = mean(Abundance), .groups = "drop")

## ----------------------------------------------------------------------------
## 6. 为两个季节统一堆积顺序和颜色
##    Define stacking order and color palette shared by Spring & Winter
##
## 要求 / Requirements:
## - Taxon_plot 包含：union_top 中的 family + "Others"
## - 堆积顺序：最大的 family 在最底层；Others 在最顶层
##   (在 ggplot 中，factor 的最后一个 level 会画在最底层)
## ----------------------------------------------------------------------------

# 合并 spring & winter，用于计算整体丰度排序
# Combine Spring & Winter to compute global ordering
combined_entry <- bind_rows(
  spring_family_entry |> mutate(Season = "Spring"),
  winter_family_entry |> mutate(Season = "Winter")
)

family_order <- combined_entry |>
  filter(Taxon_plot != "Others") |>
  group_by(Taxon_plot) |>
  summarise(total = sum(mean_Abundance), .groups = "drop") |>
  arrange(total) |>          # 升序：small → big，这样最大的在最后一个 level
  pull(Taxon_plot)

## 注意：family_order 只包含 union_top 中实际出现的 family
## Note: family_order contains all real families in union_top (excluding Others)

# factor levels：Others 放在最前（堆积条最上），其余 family 按丰度从小到大
# With ggplot stacking, first level is on top; last level is at bottom.
final_levels <- c("Others", family_order)

spring_family_entry$Taxon_plot <- factor(
  spring_family_entry$Taxon_plot,
  levels = final_levels
)

winter_family_entry$Taxon_plot <- factor(
  winter_family_entry$Taxon_plot,
  levels = final_levels
)

## 定义颜色：union_top 使用调色板，Others 固定灰色
## Define colors: union_top families use palettes, Others = grey

n_taxa <- length(family_order)  # 实际有颜色的 family 数量（不含 Others）

cols_main <- c(
  brewer.pal(min(n_taxa, 12), "Paired"),
  brewer.pal(max(n_taxa - 12, 8), "Set3")
) |> head(n_taxa)

names(cols_main) <- family_order
cols <- c(cols_main, Others = "gray80")

## ----------------------------------------------------------------------------
## 7. 绘图：Spring & Winter 分别作图，共用同一图例和颜色
##    Plot Spring and Winter separately, sharing legend and colors
## ----------------------------------------------------------------------------

p_spring <- ggplot(spring_family_entry,
                   aes(x = Entry, y = 100 * mean_Abundance, fill = Taxon_plot)) +
  geom_col(width = 1) +
  ylab("Relative abundance (%)") +
  xlab("Genotypes (Spring panel)") +
  theme_few() +
  scale_fill_manual(values = cols, drop = FALSE, name = "family",guide = guide_legend(ncol = 10)) +
  theme(
    legend.position = "top",           # ← 图例放上面 / legend at top
    legend.direction = "horizontal",   # ← 横向排列 / horizontal legend
    legend.box = "horizontal",         # ← 所有条目在同一行  
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )

p_winter <- ggplot(winter_family_entry,
                   aes(x = Entry, y = 100 * mean_Abundance, fill = Taxon_plot)) +
  geom_col(width = 1) +
  ylab("Relative abundance (%)") +
  xlab("Genotypes (Winter panel)") +
  theme_few() +
  scale_fill_manual(values = cols, drop = FALSE, name = "family",guide = guide_legend(ncol = 10)) +
  theme(
    legend.position = "top",           # ← 图例放上面
    legend.direction = "horizontal",
    legend.box = "horizontal",
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank()
  )


## 如果需要在一张图中展示两个 panel，可以使用 patchwork 合并并共享图例
## If you want a combined figure for submission, collect guides with patchwork:
p_family <- (p_spring + p_winter) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "top")
p_family

ggsave("family_composition.pdf",p_family,width=10,height=6,dpi=300)
