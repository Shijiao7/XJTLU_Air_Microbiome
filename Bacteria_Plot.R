setwd("/Users/candiceqi/Desktop/lab/Air_Bac")
library(tidyverse)
library(data.table)
library(devtools)
library(wilkoxmisc)
library(reshape2) 
library(readr)
library(plyr)
library(dplyr)
library(tidyr)

library(tidyverse)

# 读取数据
otu_file <- "otu_data.tsv"
table <- read_tsv(otu_file, show_col_types = FALSE)

# 识别RNA列
rna_cols <- grep("^R", colnames(table), value = TRUE)

# 提取RNA数据
RNA_ASV_OTU <- table %>% select(OTU_ASV, all_of(rna_cols))

# 提取DNA数据
DNA_ASV_OTU <- table %>% select(-all_of(rna_cols))


write.table(
  DNA_ASV_OTU,
  file = "DNA_ASV_OTU.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  RNA_ASV_OTU,
  file = "RNA_ASV_OTU.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

###############


otu_data <- DNA_ASV_OTU

otu_data <- otu_data %>% gather(Sample, Count, 2:(ncol(otu_data))) %>% filter(Count > 0)
names(otu_data)[1] <- c("ASV")
otu_data <- ddply(otu_data, .(Sample, ASV), summarise, count = sum(Count))
## Convert count to relative abundance and add column (require dplyr)
otu_data <- ddply(otu_data, .(Sample), mutate, RelativeAbundance = (count * 100) / sum(count))


write_tsv(otu_data,"DNA_ASV_abundance.tidy.txt")

library(wilkoxmisc)
library(reshape2)
library(dplyr) 
library(ggplot2)

OTU <- read_tsv("DNA_ASV_abundance.tidy.txt")


library(dplyr)
library(readr)

# 读取分类表
tax_data <- read_tsv(
  "tax_data.tsv",
  show_col_types = FALSE
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
write_tsv(OTU,"DNA_ASV_abundance_taxonomy.tidy.txt")
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}
OTUTable <- collapse_taxon_table(OTU, n = 21, Rank = "Class")

# 读取 metadata（先读，后面要用 Season）
sam_data <- read_tsv(
  "sam_data.txt",
  show_col_types = FALSE
)

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




###########################################Venn-Plot
OTU_Ven <- read_tsv("DNA_ASV_abundance_taxonomy.tidy.txt")
# 读取 metadata（先读，后面要用 Season）
sam_data <- read_tsv(
  "sam_data.txt",
  show_col_types = FALSE
)

# 将 Season 合并进 OTU（按 Sample）
OTU_Ven <- OTU_Ven %>%
  left_join(
    sam_data %>% select(Sample, Season,Group),
    by = c("Sample" = "Sample")
  )
library(dplyr)
library(ggvenn)
library(tidyr)
group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

library(dplyr)
library(ggvenn)
library(tidyr)

# 1. 筛选出需要的 Group 和 Genus
genus_data <- OTU_Ven %>%
  filter(Group %in% c("Underground", "Aboveground", "Elevator", "Outdoor")) %>%
  select(Sample, Group, Genus) %>%
  distinct()  # 保证每个 Genus 只出现一次

# 2. 获取每个 Group 中出现的 Genus
genus_by_group <- genus_data %>%
  group_by(Group) %>%
  summarise(Genus = list(Genus)) %>%
  ungroup()

venn_data <- genus_data %>%
  split(.$Group) %>%
  lapply(`[[`, "Genus")

venn_data <- venn_data[c("Underground", "Aboveground", "Elevator", "Outdoor")]

ggvenn(
  venn_data,
  fill_color = group_colors,
  stroke_color = NA,
  set_name_color = "black",
  show_percentage = TRUE
)



Pathogen_ASV <- read_tsv("DNA_Pathogen_ASV_table.tsv")

library(dplyr)

# 1. 去掉 NA 的 Genus（很重要）
genus_data_clean <- genus_data %>%
  filter(!is.na(Genus))

# 2. 找到四个 Group 的交集 Genus
genus_list <- genus_data_clean %>%
  split(.$Group) %>%
  lapply(function(df) unique(df$Genus))

common_genus <- Reduce(intersect, genus_list)

# 转为数据框
common_genus_df <- tibble(Genus = common_genus)

# 3. 提取病原体 Genus 列（去重）
pathogen_genus <- pathogen_genus_bsl2_3 %>%
  filter(!is.na(Genus)) %>%
  distinct(Genus)

# 4. 匹配并标注
result_table <- common_genus_df %>%
  left_join(
    pathogen_genus %>% mutate(Pathogen = "Pathogen"),
    by = "Genus"
  ) %>%
  mutate(
    Pathogen = ifelse(is.na(Pathogen), "NotPathogen", Pathogen)
  )

# 5. 查看结果
print(result_table)

# 如需保存
write.csv(result_table, "common_genus_pathogen_annotation.csv", row.names = FALSE)
########################################## Part-Pathogen
Tax_ASV <- read_tsv("DNA_ASV_abundance_taxonomy.tidy.txt")
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

# 读取病原体数据
pathogen_data <- read_csv("/Users/candiceqi/Desktop/tmp/pathegon_bacteriaTaxa_list.csv")

# 筛选出 Biosafety level 为 BSL-2 和 BSL-3 的行
pathogen_bsl2_3 <- pathogen_data %>%
  dplyr::filter(BiosafetyLevel %in% c("BSL-2", "BSL-3"))

# 提取 Pathogen name 中的 Genus 名称（前部分）
pathogen_bsl2_3 <- pathogen_bsl2_3 %>%
  dplyr::mutate(Genus = word(PathogenName, 1))  # 提取物种名称中的 Genus 部分

# 保存为 tsv 文件
write_tsv(pathogen_bsl2_3, "Pathogen_Genus_BSL2-3.tsv")


# 读取Tax_ASV数据
Tax_ASV <- read_tsv("DNA_ASV_abundance_taxonomy.tidy.txt")

# 读取 Pathogen_Genus_BSL2-3 数据
pathogen_genus_bsl2_3 <- read_tsv("Pathogen_Genus_BSL2-3.tsv")

# 过滤 Genus 存在于 Pathogen_Genus_BSL2-3 的行
Pathogen_ASV <- Tax_ASV %>%
  filter(Genus %in% pathogen_genus_bsl2_3$Genus)

# 保存生成的 table
write_tsv(Pathogen_ASV, "DNA_Pathogen_ASV_table.tsv")
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
Pathogen_ASV <- read_tsv("DNA_Pathogen_ASV_table.tsv")

# 将 Season 合并进 Pathogen_ASV 数据（按 Sample）
Pathogen_ASV <- Pathogen_ASV %>%
  left_join(sam_data %>% select(Sample, Season, Group), by = c("Sample" = "Sample"))

# 过滤四个大类和四个季节数据
valid_groups <- c("Underground", "Aboveground", "Elevator", "Outdoor")
valid_seasons <- c("Spring", "Summer", "Autumn", "Winter")

Pathogen_ASV_filtered <- Pathogen_ASV %>%
  filter(Group %in% valid_groups & Season %in% valid_seasons)

write_tsv(Pathogen_ASV_filtered, "DNA_Pathogen_Relative_Abundance_Metadata.tsv")


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
write_tsv(pathogen_stats, "DNA_Pathogen_Relative_Abundance_Stats_filtered.tsv")
#pathogen_stats <- read_tsv("DNA_Pathogen_Relative_Abundance_Stats_filtered.tsv")
# 定义颜色
season_colors <- c(
  "Autumn" = "#d6c65c",
  "Winter" = "#ebb17c",
  "Spring" = "#979797",
  "Summer" = "#e4dfc3"
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

##############################SumPathoge_RA

###################################### HeatMap of Pathogen
Pathogen_ASV_metadata <- read_tsv("DNA_Pathogen_Relative_Abundance_Metadata.tsv")
if ("package:plyr" %in% search()) {
  detach("package:plyr", unload = TRUE)
}
Pathogen_ASV_metadata <- collapse_taxon_table(Pathogen_ASV_metadata, n = 16, Rank = "Genus")

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

# 显示数据分布
cat("数据分布统计:\n")
cat("最小值:", min(df_avg_complete$average_RA), "\n")
cat("中位数:", median(df_avg_complete$average_RA), "\n")
cat("平均值:", mean(df_avg_complete$average_RA), "\n")
cat("最大值:", max(df_avg_complete$average_RA), "\n")

# 查看百分位数
quantiles <- quantile(df_avg_complete$average_RA, probs = seq(0, 1, 0.1))
print(quantiles)
write_tsv(df_avg_complete, "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/DNA_Pathogen_Average_RA_heatmapData.tsv")
# 10. 绘制热图 - 使用经典的Set3配色方案
df_avg_complete <- read_tsv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/DNA_Pathogen_Average_RA_heatmapData_Top16.tsv")
#######33


library(RColorBrewer)
library(viridis)

Heatmap_Pathogen_Set3 <- ggplot(df_avg_complete, 
                                aes(x = Season, y = Genus, fill = average_RA)) +
  geom_tile(color = "white", linewidth = 0.3, width = 1, height = 1) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  
  # plasma 默认：小值=紫色，大值=黄色（direction = 1）
  scale_fill_viridis(
    option = "plasma",  
    direction = 1,       # 改为 1（默认），小值紫色，大值黄色
    name = "Relative\nAbundance",
    na.value = "grey90",
    limits = c(0, max(df_avg_complete$average_RA))
  ) +
  
  labs(
    x = "Season",
    y = "Top 20 Pathogen Genera",
    title = "Average Relative Abundance of Top 20 Pathogen Genera (Plasma配色)"
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
#######3
# 方案1: 使用RColorBrewer的Set3配色（离散色）
Heatmap_Pathogen_Set3 <- ggplot(df_avg_complete, 
                                aes(x = Season, y = Genus, fill = average_RA)) +
  geom_tile(color = "white", linewidth = 0.3, width = 1, height = 1) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  
  # 方案A: 使用RColorBrewer的Set3配色（离散色）
  scale_fill_distiller(
    palette = "Plasma",  # 或 "RdYlBu", "RdYlGn"
    direction = -1,        # 反转颜色方向
    name = "Relative\nAbundance",
    na.value = "grey90",
    limits = c(0, max(df_avg_complete$average_RA))
  ) +
  
  labs(
    x = "Season",
    y = "Top 20 Pathogen Genera",
    title = "Average Relative Abundance of Top 20 Pathogen Genera (Set3配色)"
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

print(Heatmap_Pathogen_Set3)

##########
# 1. 确保 dplyr 已加载
library(dplyr)

# 2. 直接定义函数（不依赖 wilkoxmisc 包）
collapse_taxon_table <- function(
    TaxonTable,
    Rank = "Phylum",
    n = 8,
    UnclassifiedTerms = c("")
) {
  
  # Prepare taxon table
  message("Preparing taxon table...")
  
  TaxonTable <- TaxonTable %>%
    rename("Rank" = all_of(Rank)) %>%
    select(Sample, Rank, RelativeAbundance) %>%
    group_by(Sample, Rank) %>%
    summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    mutate(Rank = ifelse(Rank %in% UnclassifiedTerms, "Unclassified", Rank)) %>%
    ungroup()
  
  # Identify top taxa
  message("Identifying top taxa...")
  nSamples <- TaxonTable %>% pull(Sample) %>% unique() %>% length()
  
  TopTaxa <- TaxonTable %>%
    group_by(Rank) %>%
    summarise(MeanRelativeAbundance = sum(RelativeAbundance) / nSamples, .groups = "drop") %>%
    arrange(desc(MeanRelativeAbundance)) %>%
    filter(Rank != "Unclassified") %>%
    slice_head(n = n - 1) %>%
    pull(Rank)
  
  # Collapse the table
  message("Collapsing table...")
  TaxonTable <- TaxonTable %>%
    mutate(Rank = ifelse(Rank %in% TopTaxa, Rank, "Minor/Unclassified")) %>%
    group_by(Sample, Rank) %>%
    summarise(RelativeAbundance = sum(RelativeAbundance), .groups = "drop") %>%
    ungroup()
  
  # Rename the Rank column back to original name
  names(TaxonTable) <- c("Sample", Rank, "RelativeAbundance")
  
  # Return
  return(TaxonTable)
}

# 3. 运行你的分析
Pathogen_ASV_metadata <- collapse_taxon_table(Pathogen_ASV_metadata, n = 21, Rank = "Genus")

# 4. 查看结果
head(Pathogen_ASV_metadata)
#######33
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