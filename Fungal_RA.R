setwd("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity")
## Load the packages
library(readr)
library(reshape2)
library(dplyr)
library(GUniFrac)
library(vegan)
library(stringr)
library(tidyr)
library(wilkoxmisc)
library(fossil)
library(OTUtable)
library(plyr)

# 加载必要的包
library(tidyverse)
library(vegan)  # 用于计算物种累积曲线
library(ggplot2)


# 定义颜色
season_colors <- c(
  "Spring" = "#979797",
  "Summer" = "#e4dfc3",
  "Autumn" = "#d6c65c",
  "Winter" = "#ebb17c"
)

group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

############### OTU table
table <- read_tsv('/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/feature_table.tsv', show_col_types = FALSE)

# 识别RNA列
rna_cols <- grep("^R", colnames(table), value = TRUE)

# 提取RNA数据
RNA_ASV_OTU <- table %>% select(ASV, all_of(rna_cols))

# 提取DNA数据
DNA_ASV_OTU <- table %>% select(-all_of(rna_cols))


write.table(
  DNA_ASV_OTU,
  file = "fungi_DNA_ASV_OTU.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  RNA_ASV_OTU,
  file = "fungi_RNA_ASV_OTU.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
otu_data <- DNA_ASV_OTU

otu_data <- otu_data %>% gather(Sample, Count, 2:(ncol(otu_data))) %>% filter(Count > 0)
names(otu_data)[1] <- c("ASV")
otu_data <- ddply(otu_data, .(Sample, ASV), summarise, count = sum(Count))
## Convert count to relative abundance and add column (require dplyr)
otu_data <- ddply(otu_data, .(Sample), mutate, RelativeAbundance = (count * 100) / sum(count))


write_tsv(otu_data,"Fungi_ASV_abundance.tidy.txt")


library(wilkoxmisc)
library(reshape2)
library(dplyr) 
library(ggplot2)

OTU <- read_tsv("Fungi_ASV_abundance.tidy.txt")


library(dplyr)
library(readr)

# 读取分类表
tax_data <- read_csv(
  '/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/clean_tax.csv'
)

# 再把 tax 信息合并进 OTU
OTU <- OTU %>%
  left_join(
    tax_data %>%
      select(
        ASV,
        Phylum,
        Class,
        Order,
        Family,
        Genus
      ),
    by = "ASV"
  )
write_tsv(OTU,"Fungi_ASV_abundance_taxonomy.tidy.txt")
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}

library(dplyr)
library(stringr)
###########Class_level
top21_class <- OTU %>%
  filter(
    !is.na(Class),
    str_to_lower(Class) != "nan"
  ) %>%
  group_by(Class) %>%
  summarise(
    TotalAbundance = sum(RelativeAbundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 21) %>%
  pull(Class)




OTUTable <- OTU %>%
  mutate(
    Class_final = case_when(
      !is.na(Class) &
        str_to_lower(Class) != "nan" &
        Class %in% top21_class ~ Class,
      TRUE ~ "Minor/Unclassified"
    )
  ) %>%
  group_by(Sample, Class = Class_final) %>%
  summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )




# 读取 metadata（先读，后面要用 Season）

sam_data <- read_csv(
  '/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/metadata_fungi.csv'
)

# 将 Season 合并进 OTU（按 Sample）
# 将 Season 合并进 OTU（按 Sample）
OTUTable <- OTUTable %>%
  left_join(
    sam_data %>% select(Sample, Season,Group),
    by = c("Sample" = "Sample")
  )



library(dplyr)
library(tidyr)
library(ggplot2)

# 1. 过滤掉 Minor/Unclassified 的行
OTUTable_filtered <- OTUTable %>%
  filter(Class != "Minor/Unclassified")  # 去掉 Minor/Unclassified

# 2. 指定季节顺序
season_order <- c("Autumn", "Winter", "Spring", "Summer")

# 3. 过滤掉 NA 季节，并确保只有这四个季节
OTUTable_filtered <- OTUTable_filtered %>%
  filter(Season %in% season_order) %>%
  mutate(Season = factor(Season, levels = season_order))

# 4. 计算每个 Group 和 Season 下每个 Class 的 RelativeAbundance 平均值
df_avg <- OTUTable_filtered %>%
  group_by(Group, Season, Class) %>%
  summarise(average_RA = mean(RelativeAbundance, na.rm = TRUE), .groups = "drop")

# 5. 获取所有可能的 Class 和 Group 的组合（只包含指定的季节）
all_combinations <- expand.grid(
  Group = unique(df_avg$Group),
  Season = factor(season_order, levels = season_order),
  Class = unique(df_avg$Class)
)

# 6. 将原始数据框与所有组合进行合并，确保每个 Class 在所有 Group 和 Season 下都有值
df_avg_complete <- left_join(all_combinations, df_avg, by = c("Group", "Season", "Class")) %>%
  mutate(average_RA = ifelse(is.na(average_RA), 0, average_RA))  # 用 0 填补缺失值

# 7. 绘制热图
# 使用自定义渐变颜色的热图
HeaMap1 <- ggplot(df_avg_complete, aes(x = Season, y = Class, fill = average_RA)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(. ~ Group, scales = "free") +
  
  # 方案C: 自定义渐变颜色 - 红黄蓝经典配色
  scale_fill_gradientn(
    colours = c("#4575B4", "#91BFDB", "#E0F3F8", "#FFFFBF", "#FEE090", "#FC8D59", "#D73027"),
    name = "Relative Abundance",
    na.value = "grey90",
    limits = c(0, max(df_avg_complete$average_RA, na.rm = TRUE))
  ) +
  
  labs(x = "Season", y = "Class") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    
    # === Facet strip：极简 & 对齐 ===
    strip.background = element_rect(
      fill = "grey95",
      colour = NA
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(size = 8),
    
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.4, "cm")
  )


# 显示图形
print(HeaMap1)
write_tsv(df_avg_complete,"Fungi__Top20_RA.txt")

############################Genus_Top
###########Genus_level

########### Genus_level
top21_genus <- OTU %>%
  filter(
    !is.na(Genus),
    str_to_lower(Genus) != "nan"
  ) %>%
  group_by(Genus) %>%
  summarise(
    TotalAbundance = sum(RelativeAbundance),
    .groups = "drop"
  ) %>%
  arrange(desc(TotalAbundance)) %>%
  slice_head(n = 21) %>%
  pull(Genus)
OTUTable_genus <- OTU %>%
  mutate(
    Genus_final = case_when(
      !is.na(Genus) &
        str_to_lower(Genus) != "nan" &
        Genus %in% top21_genus ~ Genus,
      TRUE ~ "Minor/Unclassified"
    )
  ) %>%
  group_by(Sample, Genus = Genus_final) %>%
  summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )
OTUTable_genus <- OTUTable_genus %>%
  left_join(
    sam_data %>% select(Sample, Season, Group),
    by = "Sample"
  )
OTUTable_genus_filtered <- OTUTable_genus %>%
  filter(Genus != "Minor/Unclassified")

season_order <- c("Autumn", "Winter", "Spring", "Summer")

OTUTable_genus_filtered <- OTUTable_genus_filtered %>%
  filter(Season %in% season_order) %>%
  mutate(Season = factor(Season, levels = season_order))
df_avg_genus <- OTUTable_genus_filtered %>%
  group_by(Group, Season, Genus) %>%
  summarise(
    average_RA = mean(RelativeAbundance, na.rm = TRUE),
    .groups = "drop"
  )
all_combinations_genus <- expand.grid(
  Group = unique(df_avg_genus$Group),
  Season = factor(season_order, levels = season_order),
  Genus = unique(df_avg_genus$Genus)
)

df_avg_genus_complete <- left_join(
  all_combinations_genus,
  df_avg_genus,
  by = c("Group", "Season", "Genus")
) %>%
  mutate(average_RA = ifelse(is.na(average_RA), 0, average_RA))
HeatMap_Genus <- ggplot(
  df_avg_genus_complete,
  aes(x = Season, y = Genus, fill = average_RA)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(. ~ Group, scales = "free") +
  
  scale_fill_gradientn(
    colours = c(
      "#4575B4", "#91BFDB", "#E0F3F8",
      "#FFFFBF",
      "#FEE090", "#FC8D59", "#D73027"
    ),
    name = "Relative Abundance",
    limits = c(0, max(df_avg_genus_complete$average_RA, na.rm = TRUE)),
    na.value = "grey90"
  ) +
  
  labs(x = "Season", y = "Genus") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text = element_text(face = "bold", size = 10),
    
    axis.text.y = element_text(size = 7),
    
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(1.2, "cm"),
    legend.key.width = unit(0.4, "cm")
  )

print(HeatMap_Genus)


write_tsv(df_avg_genus_complete,"Fungi_Genus_Top20_RA.txt")

library(dplyr)
library(readr)

genus_group_summary <- df_avg_genus_complete %>%
  group_by(Group, Genus) %>%
  summarise(
    mean_RA = mean(average_RA, na.rm = TRUE),
    sd_RA   = sd(average_RA, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Group, desc(mean_RA))

write_tsv(
  genus_group_summary,
  "Genus_average_RA_by_Group_mean_sd.tsv"
)


########################## Venn_Group

OTU_Ven <- read_tsv("Fungi_ASV_abundance_taxonomy.tidy.txt")
# 读取 metadata（先读，后面要用 Season）

# 将 Season 合并进 OTU（按 Sample）
OTU_Ven <- OTU_Ven %>%
  left_join(
    sam_data %>% select(Sample, Season,Group),
    by = c("Sample" = "Sample")
  )
library(dplyr)
library(ggvenn)
library(tidyr)
library(dplyr)
library(ggvenn)
library(tidyr)
library(stringr)

# 定义颜色
group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

# 1. 先过滤掉 Genus 中的 NA / nan（关键）
genus_data <- OTU_Ven %>%
  filter(
    Group %in% c("Underground", "Aboveground", "Elevator", "Outdoor"),
    !is.na(Genus),
    str_to_lower(Genus) != "nan"
  ) %>%
  select(Group, Genus) %>%
  distinct()   # presence / absence，只保留一次

# 2. 按 Group 生成 Genus 列表
venn_data <- genus_data %>%
  split(.$Group) %>%
  lapply(`[[`, "Genus")

# 3. 固定 Group 顺序（防止 ggvenn 自动重排）
venn_data <- venn_data[c(
  "Underground",
  "Aboveground",
  "Elevator",
  "Outdoor"
)]

# 4. 画 Venn 图
ggvenn(
  venn_data,
  fill_color = group_colors,
  stroke_color = NA,
  set_name_color = "black",
  show_percentage = TRUE
)

total_genus <- genus_data %>%
  pull(Genus) %>%
  unique() %>%
  length()

total_genus

############################Venn-Season
OTU_Ven <- read_tsv("Fungi_ASV_abundance_taxonomy.tidy.txt")
# 读取 metadata（先读，后面要用 Season）

# 将 Season 合并进 OTU（按 Sample）
OTU_Ven <- OTU_Ven %>%
  left_join(
    sam_data %>% select(Sample, Season,Group),
    by = c("Sample" = "Sample")
  )
library(dplyr)
library(ggvenn)
library(tidyr)
library(dplyr)
library(ggvenn)
library(tidyr)
library(stringr)

# 定义颜色
season_colors <- c(
  "Spring" = "#979797",
  "Summer" = "#e4dfc3",
  "Autumn" = "#d6c65c",
  "Winter" = "#ebb17c"
)

# 1. 先过滤掉 Genus 中的 NA / nan（关键）
genus_data <- OTU_Ven %>%
  filter(
    Season %in% c("Autumn", "Winter", "Spring", "Summer"),
    !is.na(Genus),
    str_to_lower(Genus) != "nan"
  ) %>%
  select(Season, Genus) %>%
  distinct()   # presence / absence，只保留一次

# 2. 按 Group 生成 Genus 列表
venn_data <- genus_data %>%
  split(.$Season) %>%
  lapply(`[[`, "Genus")

# 3. 固定 Group 顺序（防止 ggvenn 自动重排）
venn_data <- venn_data[c(
  "Autumn", "Winter", "Spring", "Summer"
)]

# 4. 画 Venn 图
ggvenn(
  venn_data,
  fill_color = season_colors,
  stroke_color = NA,
  set_name_color = "black",
  show_percentage = TRUE
)



###########################3
library(dplyr)
library(stringr)

# 1. 再次确保 Genus 干净（防御性写法）
genus_data_clean <- genus_data %>%
  filter(
    !is.na(Genus),
    str_to_lower(Genus) != "nan"
  )

# 2. 按 Group 拆分 Genus 列表
genus_list <- genus_data_clean %>%
  split(.$Group) %>%
  lapply(function(df) unique(df$Genus))

# 3. 求四个 Group 的交集
common_genus <- Reduce(intersect, genus_list)

common_genus_df <- tibble(Genus = common_genus)

# 4. 处理病原体表（关键：统一格式）
pathogen_genus <- pathogen_genus_bsl2_3 %>%
  filter(!is.na(Genus)) %>%
  mutate(
    Genus = str_trim(Genus),          # 去空格
    Genus = str_replace(Genus, "_.*$", "")  # 去掉类似 Alternaria_sp
  ) %>%
  distinct(Genus)

# 5. 同样处理你的 OTU Genus（保证能 match）
common_genus_df <- common_genus_df %>%
  mutate(
    Genus_clean = str_trim(Genus),
    Genus_clean = str_replace(Genus_clean, "_.*$", "")
  )

# 6. 匹配并标注
result_table <- common_genus_df %>%
  left_join(
    pathogen_genus %>% mutate(Pathogen = "Pathogen"),
    by = c("Genus_clean" = "Genus")
  ) %>%
  mutate(
    Pathogen = ifelse(is.na(Pathogen), "NotPathogen", Pathogen)
  ) %>%
  select(Genus, Pathogen)

# 7. 查看结果
print(result_table)

# 8. 可选：保存
write.csv(result_table, "fungi_common_genus_pathogen_annotation.csv", row.names = FALSE)
############################Pathogen

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 读取病原体数据
pathogen_data <- read_csv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/pathegon_fungiTaxa_list.csv")

# 筛选出 Biosafety level 为 BSL-2 和 BSL-3 的行
pathogen_bsl2_3 <- pathogen_data %>%
  dplyr::filter(BiosafetyLevel %in% c("BSL-2", "BSL-3"))

# 提取 Pathogen name 中的 Genus 名称（前部分）
pathogen_bsl2_3 <- pathogen_bsl2_3 %>%
  dplyr::mutate(Genus = word(PathogenName, 1))  # 提取物种名称中的 Genus 部分

# 保存为 tsv 文件
write_tsv(pathogen_bsl2_3, "Fungi_Pathogen_Genus_BSL2-3.tsv")


# 读取Tax_ASV数据
Tax_ASV <- read_tsv("Fungi_ASV_abundance_taxonomy.tidy.txt")

# 读取 Pathogen_Genus_BSL2-3 数据
pathogen_genus_bsl2_3 <- read_tsv("Fungi_Pathogen_Genus_BSL2-3.tsv")

# 过滤 Genus 存在于 Pathogen_Genus_BSL2-3 的行
Pathogen_ASV <- Tax_ASV %>%
  filter(Genus %in% pathogen_genus_bsl2_3$Genus)

# 保存生成的 table
write_tsv(Pathogen_ASV, "Fungi_DNA_Pathogen_ASV_table.tsv")



library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(readr)


library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)

# 读取 Pathogen ASV 数据
Pathogen_ASV <- read_tsv("Fungi_DNA_Pathogen_ASV_table.tsv")

# 将 Season 合并进 Pathogen_ASV 数据（按 Sample）
Pathogen_ASV <- Pathogen_ASV %>%
  left_join(sam_data %>% select(Sample, Season, Group), by = c("Sample" = "Sample"))

# 过滤四个大类和四个季节数据
valid_groups <- c("Underground", "Aboveground", "Elevator", "Outdoor")
valid_seasons <- c("Spring", "Summer", "Autumn", "Winter")

Pathogen_ASV_filtered <- Pathogen_ASV %>%
  filter(Group %in% valid_groups & Season %in% valid_seasons)

write_tsv(Pathogen_ASV_filtered, "Fungi_DNA_Pathogen_Relative_Abundance_Metadata.tsv")



library(dplyr)
library(vegan)

alpha_pathogen <- Pathogen_ASV_filtered %>%
  filter(!is.na(Genus)) %>%
  group_by(Sample, Group, Season) %>%
  summarise(
    Richness = n_distinct(ASV),
    Shannon  = diversity(RelativeAbundance, index = "shannon"),
    Evenness = ifelse(Richness > 1, Shannon / log(Richness), NA),
    .groups = "drop"
  )
library(tidyr)

alpha_long <- alpha_pathogen %>%
  pivot_longer(
    cols = c(Richness, Shannon, Evenness),
    names_to = "alpha",
    values_to = "value"
  ) %>%
  mutate(
    Season = factor(Season, levels = c("Autumn", "Winter", "Spring", "Summer")),
    Group  = factor(Group,  levels = c("Underground", "Aboveground", "Elevator", "Outdoor"))
  ) %>%
  filter(!is.na(value))
write_tsv(alpha_long, "Fungi_DNA_Pathogen_AlphaDiversity.tsv")
library(rstatix)

kw_group <- alpha_long %>%
  group_by(alpha) %>%
  kruskal_test(value ~ Group)

pw_group <- alpha_long %>%
  group_by(alpha) %>%
  wilcox_test(
    value ~ Group,
    p.adjust.method = "BH"
  )
kw_season <- alpha_long %>%
  group_by(alpha) %>%
  kruskal_test(value ~ Season)

pw_season <- alpha_long %>%
  group_by(alpha) %>%
  wilcox_test(
    value ~ Season,
    p.adjust.method = "BH"
  )
write_tsv(kw_group,  "Pathogen_Alpha_Kruskal_Group.tsv")
write_tsv(pw_group,  "Pathogen_Alpha_Wilcoxon_Group_pairwise.tsv")
write_tsv(kw_season, "Pathogen_Alpha_Kruskal_Season.tsv")
write_tsv(pw_season, "Pathogen_Alpha_Wilcoxon_Season_pairwise.tsv")
library(ggplot2)
library(ggpubr)
group_pairs <- list(
  c("Underground", "Aboveground"),
  c("Underground", "Elevator"),
  c("Underground", "Outdoor"),
  c("Aboveground", "Elevator"),
  c("Aboveground", "Outdoor"),
  c("Elevator", "Outdoor")
)
Pahtoge_alpha_group <-ggplot(alpha_long, aes(x = Group, y = value)) +
  geom_boxplot(aes(fill = Group), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Group), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  facet_wrap(~alpha, scales = "free_y", nrow = 1) +
  
  # ① Kruskal–Wallis（整体）
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y.npc = 0.95,
    size = 4
  ) +
  
  # ② Wilcoxon（两两）
  stat_compare_means(
    comparisons = group_pairs,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    step.increase = 0.08,
    size = 3
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Alpha diversity",
    title = "Alpha diversity of airborne pathogenic fungi across spatial environments"
  )

season_pairs <- list(
  c("Spring", "Summer"),
  c("Summer", "Autumn"),
  c("Autumn", "Winter"),
  c("Winter", "Spring")
)
Pahtoge_alpha_season <-ggplot(alpha_long, aes(x = Season, y = value)) +
  geom_boxplot(aes(fill = Season), outlier.shape = NA, width = 0.6) +
  geom_jitter(aes(color = Season), width = 0.15, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = season_colors) +
  scale_color_manual(values = season_colors) +
  facet_wrap(~alpha, scales = "free_y", nrow = 1) +
  
  # ① Kruskal–Wallis
  stat_compare_means(
    method = "kruskal.test",
    label = "p.format",
    label.y.npc = 0.95,
    size = 4
  ) +
  
  # ② Wilcoxon pairwise
  stat_compare_means(
    comparisons = season_pairs,
    method = "wilcox.test",
    p.adjust.method = "BH",
    label = "p.signif",
    hide.ns = TRUE,
    step.increase = 0.08,
    size = 3
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Alpha diversity",
    title = "Seasonal variation in alpha diversity of airborne pathogenic fungi"
  )








###################PathogenVenn
library(dplyr)
library(ggvenn)
library(stringr)
library(tidyr)
Pathogen_ASV_filtered <- read_tsv("Fungi_DNA_Pathogen_Relative_Abundance_Metadata.tsv")
#----------------------------------------------------
# 1. 定义 Season 颜色（沿用你之前的）
#----------------------------------------------------
season_colors <- c(
  "Spring" = "#979797",
  "Summer" = "#e4dfc3",
  "Autumn" = "#d6c65c",
  "Winter" = "#ebb17c"
)

#----------------------------------------------------
# 2. 过滤 & 提取 Genus（presence / absence）
#----------------------------------------------------
genus_data <- Pathogen_ASV_filtered %>%
  filter(
    Season %in% c("Autumn", "Winter", "Spring", "Summer"),
    !is.na(Genus),
    str_to_lower(Genus) != "nan"
  ) %>%
  select(Season, Genus) %>%
  distinct()

#----------------------------------------------------
# 3. 按 Season 生成 Genus list（ggvenn 输入）
#----------------------------------------------------
venn_data <- genus_data %>%
  split(.$Season) %>%
  lapply(`[[`, "Genus")

#----------------------------------------------------
# 4. 固定 Season 顺序（防止 ggvenn 自动打乱）
#----------------------------------------------------
venn_data <- venn_data[c(
  "Autumn", "Winter", "Spring", "Summer"
)]

#----------------------------------------------------
# 5. 绘制 Venn 图
#----------------------------------------------------
ggvenn(
  venn_data,
  fill_color = season_colors,
  stroke_color = NA,
  set_name_color = "black",
  set_name_size = 4,
  text_size = 4,
  show_percentage = TRUE
)


library(dplyr)
library(ggvenn)
library(stringr)
library(readr)

# 颜色
group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

# 过滤 & 提取 genus presence
genus_data_group <- Pathogen_ASV_filtered %>%
  filter(
    Group %in% c("Underground", "Aboveground", "Elevator", "Outdoor"),
    !is.na(Genus),
    str_to_lower(Genus) != "nan"
  ) %>%
  select(Group, Genus) %>%
  distinct()

# 按 Group 生成 list
venn_data_group <- genus_data_group %>%
  split(.$Group) %>%
  lapply(`[[`, "Genus")

# 固定顺序（非常重要）
venn_data_group <- venn_data_group[c(
  "Underground", "Aboveground", "Elevator", "Outdoor"
)]

# 画 Venn
ggvenn(
  venn_data_group,
  fill_color = group_colors,
  stroke_color = NA,
  set_name_color = "black",
  set_name_size = 4,
  text_size = 4,
  show_percentage = TRUE
)
library(dplyr)
library(readr)



library(dplyr)
library(tidyr)
library(stringr)
library(readr)

venn_overlap_table <- Pathogen_ASV_filtered %>%
  filter(
    !is.na(Genus),
    str_to_lower(Genus) != "nan",
    Group %in% c("Underground", "Aboveground", "Elevator", "Outdoor")
  ) %>%
  distinct(Genus, Group) %>%
  group_by(Genus) %>%
  summarise(
    overlap = paste(sort(unique(Group)), collapse = ";"),
    .groups = "drop"
  ) %>%
  arrange(desc(str_count(overlap, ";")), overlap, Genus)

write_tsv(
  venn_overlap_table,
  "Pathogen_Genus_overlap_by_Group.tsv"
)


#######################

# 计算每个大类和季节的病原体相对丰度（RA）的平均值
pathogen_stats <- Pathogen_ASV_filtered %>%
  group_by(Group, Season) %>%
  summarise(
    avg_RA = mean(RelativeAbundance, na.rm = TRUE),
    .groups = "drop"
  )

# 计算每个大类的总平均值，用于绘制虚线
group_avg <- pathogen_stats %>%
  group_by(Group) %>%
  summarise(group_avg = mean(avg_RA))

# 查看生成的表格
print(pathogen_stats)

# 保存为 CSV 文件
write_tsv(pathogen_stats, "Fungi_DNA_Pathogen_Relative_Abundance_Stats_filtered.tsv")


# 定义颜色
# 定义颜色
season_colors <- c(
  "Spring" = "#979797",
  "Summer" = "#e4dfc3",
  "Autumn" = "#d6c65c",
  "Winter" = "#ebb17c"
)

group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)


# 计算每个bar上的数字标签位置
pathogen_stats <- pathogen_stats %>%
  mutate(label_y = avg_RA + max(pathogen_stats$avg_RA) * 0.02)

# 设置 Season 顺序：Autumn → Winter → Spring → Summer
pathogen_stats <- pathogen_stats %>%
  mutate(
    Season = factor(
      Season,
      levels = c("Autumn", "Winter", "Spring", "Summer")
    )
  )


# 创建一个合并的图例数据框
legend_data <- data.frame(
  type = c(rep("Season", 4), rep("Group Avg", 4)),
  name = c(names(season_colors), names(group_colors)),
  color = c(season_colors, group_colors)
)

# 绘制柱状图
p_bar <- ggplot() +
  # 首先添加柱状图
  geom_bar(
    data = pathogen_stats,
    aes(x = Group, y = avg_RA, fill = Season),
    stat = "identity", 
    position = position_dodge(width = 0.8), 
    width = 0.7
  ) +
  
  # 为每个柱状图添加数字标签
  geom_text(
    data = pathogen_stats,
    aes(x = Group, y = label_y, label = sprintf("%.2f", avg_RA), group = Season),
    position = position_dodge(width = 0.8),
    size = 3.5,
    fontface = "bold",
    vjust = 0
  ) +
  # 添加每个大类的平均值虚线（使用group_colors）
  geom_hline(
    data = group_avg, 
    aes(yintercept = group_avg, color = Group),
    linetype = "dashed", 
    linewidth = 1
  ) +
  
  # 设置填充和颜色比例
  scale_fill_manual(
    name = "Season",
    values = season_colors,
    guide = guide_legend(order = 1)
  ) +
  scale_color_manual(
    name = "Group Average",
    values = group_colors,
    guide = guide_legend(
      order = 2,
      override.aes = list(linetype = "dashed", fill = NA)
    )
  ) +
  
  labs(x = "Group", y = "Average Relative Abundance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.spacing.y = unit(0.5, "cm"),
    legend.key = element_rect(fill = "white", color = NA)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # 为顶部留出更多空间

# 显示柱状图
print(p_bar)
#################Barplot_V2
season_avg <- pathogen_stats %>%
  group_by(Season) %>%
  summarise(
    season_avg = mean(avg_RA, na.rm = TRUE),
    .groups = "drop"
  )
pathogen_stats <- pathogen_stats %>%
  mutate(
    Season = factor(
      Season,
      levels = c("Autumn", "Winter", "Spring", "Summer")
    ),
    Season_num = as.numeric(Season),
    label_y = avg_RA + max(avg_RA, na.rm = TRUE) * 0.03
  )

season_avg <- season_avg %>%
  mutate(
    Season = factor(
      Season,
      levels = c("Autumn", "Winter", "Spring", "Summer")
    ),
    Season_num = as.numeric(Season)
  )
p_bar_season <- ggplot(pathogen_stats) +
  
  # 柱状图（每个 Season 下 4 个 Group）
  geom_bar(
    aes(x = Season, y = avg_RA, fill = Group),
    stat = "identity",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  # 数值标签
  geom_text(
    aes(
      x = Season,
      y = label_y,
      label = sprintf("%.2f", avg_RA),
      group = Group
    ),
    position = position_dodge(width = 0.8),
    size = 3.5,
    fontface = "bold",
    vjust = 0
  ) +
  
  # ⭐ 每个 Season 的平均虚线（关键修改）
  geom_segment(
    data = season_avg,
    aes(
      x = Season_num - 0.45,
      xend = Season_num + 0.45,
      y = season_avg,
      yend = season_avg
    ),
    inherit.aes = FALSE,
    linetype = "dashed",
    linewidth = 1,
    color = "black"
  ) +
  
  # 颜色
  scale_fill_manual(
    name = "Group",
    values = group_colors
  ) +
  
  # 坐标 & 主题
  labs(
    x = "Season",
    y = "Average Relative Abundance"
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    legend.position = "top",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_rect(fill = "white", color = NA)
  )

print(p_bar_season)

###################################HeatMap_Pathogen
###################################### HeatMap of Pathogen
Pathogen_ASV_metadata <- read_tsv("Fungi_DNA_Pathogen_Relative_Abundance_Metadata.tsv")
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}
#Pathogen_ASV_metadata <- collapse_taxon_table(Pathogen_ASV_metadata, n = 16, Rank = "Genus")
Pathogen_genus <- Pathogen_ASV_metadata %>%
  dplyr::filter(!is.na(Genus), Genus != "Unclassified") %>%
  dplyr::group_by(Sample, Genus, Season, Group) %>%
  dplyr::summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )
top16 <- Pathogen_genus %>%
  group_by(Genus) %>%
  summarise(
    MeanRA = mean(RelativeAbundance),
    .groups = "drop"
  ) %>%
  arrange(desc(MeanRA)) %>%
  slice_head(n = 15) %>%
  pull(Genus)
Pathogen_genus_top <- Pathogen_genus %>%
  mutate(
    Genus = ifelse(Genus %in% top16, Genus, "Minor/Unclassified")
  ) %>%
  group_by(Sample, Genus, Season, Group) %>%
  summarise(
    RelativeAbundance = sum(RelativeAbundance),
    .groups = "drop"
  )
Pathogen_ASV_metadata <- Pathogen_ASV_metadata %>%
  left_join(sam_data %>% select(Sample, Season, Group), by = c("Sample" = "Sample"))

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(readr)

# 假设您的数据已经读取并存储在 Pathogen_ASV_metadata 中
# 如果没有，请先读取数据：
# Pathogen_ASV_metadata <- read_tsv("your_file.tsv")

library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(readr)

# 1. 过滤掉 Minor/Unclassified
Pathogen_ASV_filtered <- Pathogen_ASV_metadata %>%
  filter(Genus != "Minor/Unclassified")

# 2. 筛选四个大类和四个季节
valid_groups <- c("Underground", "Aboveground", "Elevator", "Outdoor")
season_order <- c("Autumn", "Winter", "Spring", "Summer")

Pathogen_ASV_filtered <- Pathogen_ASV_filtered %>%
  filter(Group %in% valid_groups, Season %in% season_order) %>%
  mutate(
    Season = factor(Season, levels = season_order),
    Group = factor(Group, levels = valid_groups)
  )

# 3. 计算每个Genus在所有样本中的总相对丰度，选择Top 20
genus_rank <- Pathogen_ASV_filtered %>%
  group_by(Genus) %>%
  summarise(total_RA = mean(RelativeAbundance, na.rm = TRUE)) %>%
  arrange(desc(total_RA)) %>%
  slice_head(n = 15)

# 4. 只保留Top 20的Genus
Pathogen_top20 <- Pathogen_ASV_filtered %>%
  filter(Genus %in% genus_rank$Genus) %>%
  mutate(Genus = factor(Genus, levels = rev(genus_rank$Genus)))

# 5. 计算每个Group、Season和Genus的平均相对丰度
df_avg <- Pathogen_top20 %>%
  group_by(Group, Season, Genus) %>%
  summarise(average_RA = mean(RelativeAbundance, na.rm = TRUE), .groups = "drop")

# 6. 获取所有可能的组合
all_combinations <- expand.grid(
  Group = factor(valid_groups, levels = valid_groups),
  Season = factor(season_order, levels = season_order),
  Genus = unique(df_avg$Genus)
)

# 7. 合并数据，填补缺失值为0
df_avg_complete <- left_join(all_combinations, df_avg, 
                             by = c("Group", "Season", "Genus")) %>%
  mutate(average_RA = ifelse(is.na(average_RA), 0, average_RA))

# 8. 计算每个Genus在所有Group中的最大丰度，用于排序
genus_order <- df_avg_complete %>%
  group_by(Genus) %>%
  summarise(max_RA = max(average_RA)) %>%
  arrange(desc(max_RA))

df_avg_complete <- df_avg_complete %>%
  mutate(Genus = factor(Genus, levels = rev(genus_order$Genus)))

# 9. 分析数据分布，选择合适的配色
summary(df_avg_complete$average_RA)
write_tsv(df_avg_genus_complete,"Fungi_Pathogen_Genus_Top10_RA.txt")
# 显示数据分布
cat("数据分布统计:\n")
cat("最小值:", min(df_avg_complete$average_RA), "\n")
cat("中位数:", median(df_avg_complete$average_RA), "\n")
cat("平均值:", mean(df_avg_complete$average_RA), "\n")
cat("最大值:", max(df_avg_complete$average_RA), "\n")

# 查看百分位数
quantiles <- quantile(df_avg_complete$average_RA, probs = seq(0, 1, 0.1))
print(quantiles)

# 10. 绘制热图 - 使用经典的Set3配色方案

max_RA <- max(df_avg_complete$average_RA, na.rm = TRUE)

Heatmap_Pathogen_Set3 <- ggplot(
  df_avg_complete,
  aes(x = Season, y = Genus, fill = average_RA)
) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  
  # ✅ 正确的分段 plasma 映射
  scale_fill_viridis_c(
    option = "plasma",
    direction = 1,
    name = "Relative\nAbundance",
    na.value = "grey90",
    limits = c(0, max_RA),
    values = c(0, 2 / max_RA, 1)
  ) +
  
  labs(
    x = "Season",
    y = "Top 20 Pathogen Genera",
    title = "Average Relative Abundance of Top 20 Pathogen Genera"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5
    ),
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      color = "grey50"
    ),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm")
  )

print(Heatmap_Pathogen_Set3)

####################3
#########Log2Version
# 生成 log2 转换后的平均丰度
df_avg_complete <- df_avg_complete %>%
  mutate(log2_RA = log2(average_RA + 1e-6))  # 避免 log2(0)

# 绘制 log2 热图
Heatmap_Pathogen_log2 <- ggplot(df_avg_complete, 
                                aes(x = Season, y = Genus, fill = log2_RA)) +
  geom_tile(color = "white", linewidth = 0.3, width = 1, height = 1) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  
  # 使用 viridis 的连续色彩
  scale_fill_viridis_c(
    option = "plasma",
    direction = 1,
    name = "log2(Relative\nAbundance + 1e-6)",
    na.value = "grey90",
    limits = c(min(df_avg_complete$log2_RA), max(df_avg_complete$log2_RA))
  ) +
  labs(
    x = "Season",
    y = "Top 20 Pathogen Genera",
    title = "Log2-transformed Average Relative Abundance of Top 20 Pathogen Genera"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1,
      vjust = 1,
      size = 10,
      face = "bold"
    ),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    plot.title = element_text(
      size = 14, 
      face = "bold", 
      hjust = 0.5,
      margin = margin(b = 15)
    ),
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "grey90",
      color = "grey50",
      linewidth = 0.5
    ),
    strip.text = element_text(
      face = "bold",
      size = 11,
      margin = margin(5, 0, 5, 0)
    ),
    panel.spacing = unit(0.5, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.5, "cm")
  )

print(Heatmap_Pathogen_log2)

##########
##################################################
# 方法：计算每个Group的每个季节中，所有Genus的平均相对丰度（整体平均）
group_season_avg <- Pathogen_ASV_filtered %>%
  group_by(Group, Season) %>%
  summarise(
    # 这个季节中所有样本的所有ASV的平均相对丰度
    Avg_RelativeAbundance = mean(RelativeAbundance, na.rm = TRUE),
    # 标准差
    SD = sd(RelativeAbundance, na.rm = TRUE),
    # 该组内的样本数量
    Sample_Count = n_distinct(Sample),
    # 该组内的Genus数量
    Genus_Count = n_distinct(Genus),
    # 该组内的ASV总数
    ASV_Count = n(),
    .groups = "drop"
  ) %>%
  arrange(Group, Season)

# 查看结果
print(group_season_avg)

# 保存结果
write_tsv(group_season_avg, "DNA_Pathogen_Genus_Group_Season_Overall_Avg.tsv")