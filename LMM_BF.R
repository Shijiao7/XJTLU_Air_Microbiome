# 加载必要的包
library(dplyr)
library(ggplot2)
library(lubridate)
library(tidyr)
library(patchwork)
library(scales)
setwd("/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/LMM")
# 设置输出目录
outDir <- "./AlphaDiversity_Analysis/"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# 定义分组颜色
group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

# ==================== 1. 加载数据 ====================
# 加载元数据
metadata <- read.csv("/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/combined_meta.csv", 
                     stringsAsFactors = FALSE)

# 加载细菌alpha多样性数据
bacteria_alpha <- read.delim("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/Bacteria_DNA_alpha_diversity_rarefied.txt", 
                             check.names = FALSE)

# 加载真菌alpha多样性数据
fungi_alpha <- read.delim("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/Fungi_alpha_diversity_rarefied.txt", 
                          check.names = FALSE)

# ==================== 2. 数据预处理 ====================
# 处理日期
metadata$Date_full <- as.Date(paste0("20", substr(metadata$Date, 1, 2), "-",
                                     substr(metadata$Date, 3, 4), "-",
                                     substr(metadata$Date, 5, 6)))

# 计算从起始日期开始的天数（用于排序）
base_date <- as.Date("2023-09-01")
metadata$Days <- as.integer(metadata$Date_full - base_date)

# 转换温度和湿度为数值
metadata$Temperature_num <- as.numeric(metadata$Temperature)
metadata$Humidness_num <- as.numeric(gsub("%", "", metadata$Humidness))

# ==================== 3. 合并细菌数据 ====================
bacteria_data <- metadata %>%
  inner_join(bacteria_alpha, by = "Sample") %>%
  mutate(
    log_richness = log(richness),
    log_shannon = log(Shannon),
    Type = "Bacteria"
  )

# ==================== 4. 合并真菌数据 ====================
fungi_data <- metadata %>%
  inner_join(fungi_alpha, by = "Sample") %>%
  mutate(
    log_richness = log(richness),
    log_shannon = log(Shannon),
    Type = "Fungi"
  )

# ==================== 5. 定义统一的日期断点函数 ====================
# 获取数据的日期范围
date_range <- range(c(bacteria_data$Date_full, fungi_data$Date_full), na.rm = TRUE)

# 创建更详细的日期断点（每2周或每个月）
# 方法1：每2周一个断点
date_breaks_2weeks <- seq(from = floor_date(date_range[1], "month"), 
                          to = ceiling_date(date_range[2], "month"), 
                          by = "2 weeks")

# 方法2：每个月一个断点
date_breaks_month <- seq(from = floor_date(date_range[1], "month"), 
                         to = ceiling_date(date_range[2], "month"), 
                         by = "month")

# 方法3：每2个月一个断点（如果数据时间跨度长）
date_breaks_2months <- seq(from = floor_date(date_range[1], "month"), 
                           to = ceiling_date(date_range[2], "month"), 
                           by = "2 months")

# 创建日期标签格式函数
date_labels_format <- function(x) {
  format(x, "%Y-%m-%d")
}

month_labels_format <- function(x) {
  format(x, "%Y-%m")
}

# 选择最适合的断点间隔（根据数据范围自动选择）
date_span <- as.numeric(diff(date_range))
if(date_span <= 180) {
  # 半年以内的数据，每2周一个标签
  date_breaks <- date_breaks_2weeks
  date_labels <- date_labels_format
  x_axis_angle <- 45
  x_axis_hjust <- 1
} else if(date_span <= 365) {
  # 一年以内的数据，每个月一个标签
  date_breaks <- date_breaks_month
  date_labels <- month_labels_format
  x_axis_angle <- 45
  x_axis_hjust <- 1
} else {
  # 一年以上的数据，每2个月一个标签
  date_breaks <- date_breaks_2months
  date_labels <- month_labels_format
  x_axis_angle <- 45
  x_axis_hjust <- 1
}

# ==================== 6. 细菌Richness分析 ====================
# 6.1 细菌Richness - 所有点在一张图上，按Group分组，每个Group一条拟合线
p_bacteria_richness <- ggplot(bacteria_data, aes(x = Date_full, y = log_richness, color = Group)) +
  # 添加抖动点
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  # 添加平滑线（每个Group一条线）
  geom_smooth(method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Bacteria - Richness over Time",
       subtitle = "Each group has its own loess smoothing curve with 95% CI",
       x = "Date", y = "Log Richness",
       color = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 6.2 细菌Richness - 分面子图（按Group分面）
p_bacteria_richness_facet <- ggplot(bacteria_data, aes(x = Date_full, y = log_richness, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Bacteria - Richness by Group",
       x = "Date", y = "Log Richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

# ==================== 7. 细菌Shannon分析 ====================
# 7.1 细菌Shannon - 所有点在一张图上，按Group分组，每个Group一条拟合线
p_bacteria_shannon <- ggplot(bacteria_data, aes(x = Date_full, y = log_shannon, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Bacteria - Shannon over Time",
       subtitle = "Each group has its own loess smoothing curve with 95% CI",
       x = "Date", y = "Log Shannon",
       color = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 7.2 细菌Shannon - 分面子图
p_bacteria_shannon_facet <- ggplot(bacteria_data, aes(x = Date_full, y = log_shannon, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Bacteria - Shannon by Group",
       x = "Date", y = "Log Shannon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

# ==================== 8. 真菌Richness分析 ====================
# 8.1 真菌Richness - 所有点在一张图上，按Group分组，每个Group一条拟合线
p_fungi_richness <- ggplot(fungi_data, aes(x = Date_full, y = log_richness, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Fungi - Richness over Time",
       subtitle = "Each group has its own loess smoothing curve with 95% CI",
       x = "Date", y = "Log Richness",
       color = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 8.2 真菌Richness - 分面子图
p_fungi_richness_facet <- ggplot(fungi_data, aes(x = Date_full, y = log_richness, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Fungi - Richness by Group",
       x = "Date", y = "Log Richness") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

# ==================== 9. 真菌Shannon分析 ====================
# 9.1 真菌Shannon - 所有点在一张图上，按Group分组，每个Group一条拟合线
p_fungi_shannon <- ggplot(fungi_data, aes(x = Date_full, y = log_shannon, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Fungi - Shannon over Time",
       subtitle = "Each group has its own loess smoothing curve with 95% CI",
       x = "Date", y = "Log Shannon",
       color = "Group") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

# 9.2 真菌Shannon - 分面子图
p_fungi_shannon_facet <- ggplot(fungi_data, aes(x = Date_full, y = log_shannon, color = Group)) +
  geom_point(alpha = 0.6, size = 2, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_color_manual(values = group_colors) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Fungi - Shannon by Group",
       x = "Date", y = "Log Shannon") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        legend.position = "none")

# ==================== 10. 创建组合图 ====================
# 10.1 细菌Richness组合图（总图+分面图并排）
bacteria_richness_combined <- (p_bacteria_richness | p_bacteria_richness_facet) +
  plot_annotation(
    title = "Bacteria Richness Analysis",
    subtitle = "Left: All groups together | Right: Faceted by group",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# 10.2 细菌Shannon组合图
bacteria_shannon_combined <- (p_bacteria_shannon | p_bacteria_shannon_facet) +
  plot_annotation(
    title = "Bacteria Shannon Analysis",
    subtitle = "Left: All groups together | Right: Faceted by group",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# 10.3 真菌Richness组合图
fungi_richness_combined <- (p_fungi_richness | p_fungi_richness_facet) +
  plot_annotation(
    title = "Fungi Richness Analysis",
    subtitle = "Left: All groups together | Right: Faceted by group",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# 10.4 真菌Shannon组合图
fungi_shannon_combined <- (p_fungi_shannon | p_fungi_shannon_facet) +
  plot_annotation(
    title = "Fungi Shannon Analysis",
    subtitle = "Left: All groups together | Right: Faceted by group",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# 10.5 细菌所有指标组合图（上下排列）
bacteria_all <- (p_bacteria_richness / p_bacteria_shannon) +
  plot_annotation(
    title = "Bacteria Alpha Diversity Analysis",
    subtitle = "Richness (top) and Shannon (bottom) with loess smoothing",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# 10.6 真菌所有指标组合图
fungi_all <- (p_fungi_richness / p_fungi_shannon) +
  plot_annotation(
    title = "Fungi Alpha Diversity Analysis",
    subtitle = "Richness (top) and Shannon (bottom) with loess smoothing",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# ==================== 11. 细菌和真菌对比图 ====================
# 合并数据用于对比
combined_data <- bind_rows(
  bacteria_data %>% select(Date_full, log_richness, log_shannon, Group, Type),
  fungi_data %>% select(Date_full, log_richness, log_shannon, Group, Type)
)

# 11.1 Richness对比 - 所有点在一张图上，按Type分组
p_compare_richness <- ggplot(combined_data, aes(x = Date_full, y = log_richness, color = Type)) +
  geom_point(alpha = 0.5, size = 1.5, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = c("Bacteria" = "#2c7fb8", "Fungi" = "#d95f02")) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Richness Comparison: Bacteria vs Fungi",
       subtitle = "Loess smoothing with 95% confidence intervals",
       x = "Date", y = "Log Richness",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

# 11.2 Shannon对比
p_compare_shannon <- ggplot(combined_data, aes(x = Date_full, y = log_shannon, color = Type)) +
  geom_point(alpha = 0.5, size = 1.5, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1.2, alpha = 0.3) +
  scale_color_manual(values = c("Bacteria" = "#2c7fb8", "Fungi" = "#d95f02")) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Shannon Comparison: Bacteria vs Fungi",
       subtitle = "Loess smoothing with 95% confidence intervals",
       x = "Date", y = "Log Shannon",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 10),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        legend.position = "right",
        legend.title = element_text(face = "bold"))

# 11.3 对比组合图
comparison_combined <- (p_compare_richness | p_compare_shannon) +
  plot_annotation(
    title = "Bacteria vs Fungi Alpha Diversity Comparison",
    subtitle = "Richness (left) and Shannon (right) with loess smoothing",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
                  plot.subtitle = element_text(hjust = 0.5, size = 10))
  )

# 11.4 Richness按Group分面对比
p_compare_richness_facet <- ggplot(combined_data, aes(x = Date_full, y = log_richness, color = Type)) +
  geom_point(alpha = 0.5, size = 1.5, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Bacteria" = "#2c7fb8", "Fungi" = "#d95f02")) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Richness Comparison by Group",
       x = "Date", y = "Log Richness",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold", size = 12))

# 11.5 Shannon按Group分面对比
p_compare_shannon_facet <- ggplot(combined_data, aes(x = Date_full, y = log_shannon, color = Type)) +
  geom_point(alpha = 0.5, size = 1.5, position = position_jitter(width = 2, height = 0.02)) +
  geom_smooth(method = "loess", se = TRUE, size = 1, alpha = 0.3) +
  facet_wrap(~Group, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Bacteria" = "#2c7fb8", "Fungi" = "#d95f02")) +
  scale_x_date(
    breaks = date_breaks,
    labels = date_labels,
    limits = date_range
  ) +
  labs(title = "Shannon Comparison by Group",
       x = "Date", y = "Log Shannon",
       color = "Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = x_axis_angle, hjust = x_axis_hjust, size = 8),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = NA),
        strip.text = element_text(face = "bold", size = 12))

# ==================== 12. 保存所有PDF文件 ====================
# 12.1 保存细菌分析PDF
pdf(paste0(outDir, "Bacteria_Alpha_Diversity.pdf"), width = 14, height = 10)
print(bacteria_richness_combined)
print(bacteria_shannon_combined)
print(bacteria_all)
dev.off()
cat("✓ 已保存: Bacteria_Alpha_Diversity.pdf\n")

# 12.2 保存真菌分析PDF
pdf(paste0(outDir, "Fungi_Alpha_Diversity.pdf"), width = 14, height = 10)
print(fungi_richness_combined)
print(fungi_shannon_combined)
print(fungi_all)
dev.off()
cat("✓ 已保存: Fungi_Alpha_Diversity.pdf\n")

# 12.3 保存对比分析PDF
pdf(paste0(outDir, "Bacteria_vs_Fungi_Comparison.pdf"), width = 14, height = 10)
print(comparison_combined)
print(p_compare_richness_facet)
print(p_compare_shannon_facet)
dev.off()
cat("✓ 已保存: Bacteria_vs_Fungi_Comparison.pdf\n")

# 12.4 保存单独的图表PDF（每个图表单独一页）
pdf(paste0(outDir, "Individual_Plots.pdf"), width = 10, height = 8)
print(p_bacteria_richness)
print(p_bacteria_richness_facet)
print(p_bacteria_shannon)
print(p_bacteria_shannon_facet)
print(p_fungi_richness)
print(p_fungi_richness_facet)
print(p_fungi_shannon)
print(p_fungi_shannon_facet)
print(p_compare_richness)
print(p_compare_shannon)
print(p_compare_richness_facet)
print(p_compare_shannon_facet)
dev.off()
cat("✓ 已保存: Individual_Plots.pdf\n")

# ==================== 13. 保存处理后的数据 ====================
write.csv(bacteria_data, paste0(outDir, "Bacteria_Alpha_Data.csv"), row.names = FALSE)
write.csv(fungi_data, paste0(outDir, "Fungi_Alpha_Data.csv"), row.names = FALSE)
write.csv(combined_data, paste0(outDir, "Combined_Alpha_Data.csv"), row.names = FALSE)
cat("✓ 已保存: CSV数据文件\n")

# ==================== 14. 输出摘要信息 ====================
cat("\n========== 分析完成！==========\n")
cat("所有结果保存在：", outDir, "\n\n")

cat("X轴时间设置：\n")
cat("- 数据时间范围：", format(date_range[1], "%Y-%m-%d"), "至", format(date_range[2], "%Y-%m-%d"), "\n")
cat("- 时间跨度：", date_span, "天\n")
if(date_span <= 180) {
  cat("- 刻度间隔：每2周\n")
  cat("- 标签格式：YYYY-MM-DD\n")
} else if(date_span <= 365) {
  cat("- 刻度间隔：每月\n")
  cat("- 标签格式：YYYY-MM\n")
} else {
  cat("- 刻度间隔：每2个月\n")
  cat("- 标签格式：YYYY-MM\n")
}
cat("\n")

cat("生成的文件包括：\n")
cat("1. Bacteria_Alpha_Diversity.pdf - 细菌分析（含组合图）\n")
cat("   - Richness总图+分面图\n")
cat("   - Shannon总图+分面图\n")
cat("   - Richness和Shannon上下组合图\n\n")

cat("2. Fungi_Alpha_Diversity.pdf - 真菌分析（含组合图）\n")
cat("   - Richness总图+分面图\n")
cat("   - Shannon总图+分面图\n")
cat("   - Richness和Shannon上下组合图\n\n")

cat("3. Bacteria_vs_Fungi_Comparison.pdf - 对比分析\n")
cat("   - Richness和Shannon对比（总图）\n")
cat("   - 按Group分面的Richness对比\n")
cat("   - 按Group分面的Shannon对比\n\n")

cat("4. Individual_Plots.pdf - 所有单独的图表\n\n")

cat("5. CSV数据文件：\n")
cat("   - Bacteria_Alpha_Data.csv\n")
cat("   - Fungi_Alpha_Data.csv\n")
cat("   - Combined_Alpha_Data.csv\n\n")

cat("图表说明：\n")
cat("- 所有点都按Group颜色显示\n")
cat("- 每条线都是每个Group独立的loess平滑曲线\n")
cat("- 灰色阴影区域表示95%置信区间\n")
cat("- 总图：所有Group的点在一张图上，每个Group有独立的拟合线\n")
cat("- 分面图：每个Group单独一个子图\n")
cat("- X轴现在显示更详细的时间信息（根据数据范围自动调整间隔）\n")

#######################LMM
library(mgcv)

# 确保是factor
bacteria_data$Group <- as.factor(bacteria_data$Group)
fungi_data$Group <- as.factor(fungi_data$Group)

# ==================== 15.1 细菌 ====================

gam_bac_rich <- gam(log_richness ~ s(Days, by = Group) + Group +
                      Temperature_num + Humidness_num,
                    data = bacteria_data,
                    method = "REML")

gam_bac_shan <- gam(log_shannon ~ s(Days, by = Group) + Group +
                      Temperature_num + Humidness_num,
                    data = bacteria_data,
                    method = "REML")

summary(gam_bac_rich)
summary(gam_bac_shan)
# ==================== 15.2 真菌 ====================

gam_fun_rich <- gam(log_richness ~ s(Days, by = Group) + Group +
                      Temperature_num + Humidness_num,
                    data = fungi_data,
                    method = "REML")

gam_fun_shan <- gam(log_shannon ~ s(Days, by = Group) + Group +
                      Temperature_num + Humidness_num,
                    data = fungi_data,
                    method = "REML")

summary(gam_fun_rich)
summary(gam_fun_shan)

# ==================== 15.3 合并模型 ====================

combined_data <- bind_rows(
  bacteria_data,
  fungi_data
)

combined_data$Type <- as.factor(combined_data$Type)

gam_compare <- gam(log_richness ~ 
                     s(Days, by = interaction(Type, Group)) +
                     Type + Group +
                     Temperature_num + Humidness_num,
                   data = combined_data,
                   method = "REML")

summary(gam_compare)

newdata_bac <- expand.grid(
  Days = seq(min(bacteria_data$Days), max(bacteria_data$Days), length.out = 200),
  Group = levels(bacteria_data$Group),
  Temperature_num = mean(bacteria_data$Temperature_num, na.rm = TRUE),
  Humidness_num = mean(bacteria_data$Humidness_num, na.rm = TRUE)
)

newdata_bac$pred <- predict(gam_bac_rich, newdata = newdata_bac)

p_bacteria_richness_gam <- ggplot(bacteria_data, aes(x = Days, y = log_richness, color = Group)) +
  geom_point(alpha = 0.4) +
  geom_line(data = newdata_bac, aes(y = pred), size = 1.3) +
  scale_color_manual(values = group_colors) +
  labs(title = "Bacteria Richness (GAM controlled for Temp & Humidity)") +
  theme_minimal()

newdata_fun <- expand.grid(
  Days = seq(min(fungi_data$Days), max(fungi_data$Days), length.out = 200),
  Group = levels(fungi_data$Group),
  Temperature_num = mean(fungi_data$Temperature_num, na.rm = TRUE),
  Humidness_num = mean(fungi_data$Humidness_num, na.rm = TRUE)
)

newdata_fun$pred <- predict(gam_fun_rich, newdata = newdata_fun)

p_fungi_richness_gam <- ggplot(fungi_data, aes(x = Days, y = log_richness, color = Group)) +
  geom_point(alpha = 0.4) +
  geom_line(data = newdata_fun, aes(y = pred), size = 1.3) +
  scale_color_manual(values = group_colors) +
  theme_minimal()

newdata_comb <- expand.grid(
  Days = seq(min(combined_data$Days), max(combined_data$Days), length.out = 200),
  Group = levels(combined_data$Group),
  Type = levels(combined_data$Type),
  Temperature_num = mean(combined_data$Temperature_num, na.rm = TRUE),
  Humidness_num = mean(combined_data$Humidness_num, na.rm = TRUE)
)

newdata_comb$pred <- predict(gam_compare, newdata = newdata_comb)
newdata_comb$Date_full <- base_date + newdata_comb$Days
p_compare_gam <- ggplot(combined_data, aes(x = Date_full, y = log_richness, color = Type)) +
  geom_point(alpha = 0.3) +
  geom_line(data = newdata_comb, aes(x = Date_full, y = pred), size = 1.2) +
  facet_wrap(~Group) +
  scale_color_manual(values = c("Bacteria" = "#2c7fb8", "Fungi" = "#d95f02")) +
  scale_x_date(date_breaks = "1 month", date_labels = "%Y-%m") +
  theme_minimal()
p_compare_gam <- ggplot(combined_data, 
                        aes(x = Date_full, y = log_richness, color = Type)) +
  
  # 🔹 原始点（弱化）
  geom_point(alpha = 0.25, size = 1.2) +
  
  # 🔹 GAM / 拟合曲线（强调）
  geom_line(data = newdata_comb,
            aes(x = Date_full, y = pred, color = Type),
            size = 1.2) +
  
  # 🔹 分面（控制排版）
  facet_wrap(~ Group, ncol = 4) +
  
  # 🔹 配色（保持你原本很好的一组）
  scale_color_manual(values = c("Bacteria" = "#2c7fb8",
                                "Fungi" = "#d95f02")) +
  
  # 🔹 时间轴
  scale_x_date(date_breaks = "1 month",
               date_labels = "%Y-%m",
               expand = expansion(mult = c(0.01, 0.02))) +
  
  # 🔹 坐标轴标签
  labs(x = "Sampling time",
       y = "Log richness",
       color = NULL) +
  
  # 🔹 主题（顶刊风格）
  theme_classic(base_size = 14) +
  theme(
    # 分面标题
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(face = "bold", size = 12),
    
    # 面板边框
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    
    # X轴竖排（关键）
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    
    # Y轴更清晰
    axis.text.y = element_text(color = "black"),
    
    # legend
    legend.position = "top",
    legend.text = element_text(size = 11),
    
    # 去掉多余网格
    panel.grid = element_blank()
  )

p_compare_gam

ggsave("GAM_Richness_Time.pdf",
       p_compare_gam,
       width = 16, height = 6, dpi = 300)



gam.check(gam_bac_rich)
gam.check(gam_fun_rich)
gam.check(gam_compare)
########################################GAM_V3
library(mgcv)
library(dplyr)
library(ggplot2)
library(gratia)      # 推荐：更好的GAM可视化和诊断

# ── 数据合并 ──────────────────────────────────────────────
bacteria_data$Type <- "Bacteria"
fungi_data$Type    <- "Fungi"
combined_data <- bind_rows(bacteria_data, fungi_data)

combined_data$Group    <- as.factor(combined_data$Group)
combined_data$Type     <- as.factor(combined_data$Type)
combined_data$Days     <- as.numeric(combined_data$Date_full)
combined_data$TypeGroup <- as.factor(interaction(combined_data$Type, combined_data$Group))

# ── GAM模型（修正k值）────────────────────────────────────
# ⚠️ 关键修正：k从10提升至15，避免基维度不足
gam_compare <- gam(
  log_richness ~
    s(Days, by = TypeGroup, k = 15) +   # k提升
    s(Temperature_num, k = 8) +          # 适当增加
    s(Humidness_num, k = 8) +
    Type * Group,
  data   = combined_data,
  method = "REML"
)

summary(gam_compare)
gam.check(gam_compare)

# ── 预测数据 ──────────────────────────────────────────────
newdata <- expand.grid(
  Date_full = seq(min(combined_data$Date_full),
                  max(combined_data$Date_full),
                  by = "3 days"),
  Group = levels(combined_data$Group),
  Type  = levels(combined_data$Type)
) %>%
  mutate(
    Days      = as.numeric(Date_full),
    TypeGroup = factor(interaction(Type, Group),
                       levels = levels(combined_data$TypeGroup)),
    Temperature_num = mean(combined_data$Temperature_num, na.rm = TRUE),
    Humidness_num   = mean(combined_data$Humidness_num,   na.rm = TRUE)
  )

pred <- predict(gam_compare, newdata = newdata, se.fit = TRUE)
newdata$fit   <- pred$fit
newdata$se    <- pred$se.fit
newdata$upper <- newdata$fit + 1.96 * newdata$se
newdata$lower <- newdata$fit - 1.96 * newdata$se

# ── 配色（顶刊风格）──────────────────────────────────────
col_bacteria <- "#2c7fb8"
col_fungi    <- "#d95f02"

# ── 绘图（四个正方形子图）────────────────────────────────
# 关键：coord_fixed() 或 aspect.ratio = 1 控制正方形
p_compare_richness_facet <- ggplot() +
  
  geom_point(
    data = combined_data,
    aes(x = Date_full, y = log_richness, color = Type),
    alpha = 0.45, size = 1.2,
    position = position_jitter(width = 2, height = 0.02),
    shape = 16
  ) +
  
  geom_ribbon(
    data = newdata,
    aes(x = Date_full, ymin = lower, ymax = upper, fill = Type),
    alpha = 0.18
  ) +
  
  geom_line(
    data = newdata,
    aes(x = Date_full, y = fit, color = Type),
    linewidth = 1.0
  ) +
  
  facet_wrap(~ Group, nrow = 1) +
  
  scale_color_manual(
    values = c("Bacteria" = col_bacteria, "Fungi" = col_fungi),
    labels = c("Bacteria", "Fungi")
  ) +
  scale_fill_manual(
    values = c("Bacteria" = col_bacteria, "Fungi" = col_fungi),
    labels = c("Bacteria", "Fungi")
  ) +
  
  scale_x_date(
    breaks = seq(min(combined_data$Date_full),
                 max(combined_data$Date_full),
                 by = "2 months"),
    labels = function(x) format(x, "%b\n%Y"),
    expand = c(0.02, 0)
  ) +
  
  scale_y_continuous(expand = c(0.05, 0)) +
  
  labs(
    x     = "Sampling date",
    y     = expression(Log[10]~"(observed richness)"),
    color = "Microbial type",
    fill  = "Microbial type"
  ) +
  
  # ✅ 正方形子图的核心设置
  coord_fixed(ratio = diff(range(as.numeric(combined_data$Date_full))) /
                diff(range(combined_data$log_richness, na.rm = TRUE)) * 0.25) +
  
  theme_classic(base_size = 11) +
  theme(
    # 子图标签
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 11, margin = margin(b = 4)),
    strip.placement  = "outside",
    
    # 坐标轴
    axis.text.x  = element_text(size = 8, color = "black"),
    axis.text.y  = element_text(size = 9, color = "black"),
    axis.title   = element_text(size = 10),
    axis.line    = element_line(linewidth = 0.4),
    axis.ticks   = element_line(linewidth = 0.4),
    
    # 图例
    legend.position  = "top",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.4, "cm"),
    legend.spacing.x = unit(0.3, "cm"),
    
    # 面板
    panel.spacing  = unit(0.8, "cm"),
    plot.margin    = margin(8, 8, 8, 8),
    
    # 正方形的另一种方案（如果coord_fixed效果不理想）
    aspect.ratio = 1
  )

p_compare_richness_facet

# 保存（顶刊推荐300dpi）
ggsave("richness_GAM_facet.pdf",
       plot   = p_compare_richness_facet,
       width  = 30, height = 21,
       units  = "cm", dpi = 300)
ggsave("richness_GAM_facet.png",
       plot   = p_compare_richness_facet,
       width  = 16, height = 5,
       units  = "cm", dpi = 300, bg = "white")

#######################3GAM_V4
############################################################
# 0. 加载包
############################################################
library(mgcv)
library(dplyr)
library(ggplot2)
library(gratia)
library(vegan)
library(variancePartition)

############################################################
# 1. 数据准备（关键修正）
############################################################

# 添加类型标签
bacteria_data$Type <- "Bacteria"
fungi_data$Type    <- "Fungi"

combined_data <- bind_rows(bacteria_data, fungi_data)

# 因子化
combined_data$Group <- as.factor(combined_data$Group)
combined_data$Type  <- as.factor(combined_data$Type)

# ✅ 正确时间变量（关键修正）
combined_data$Days <- as.numeric(
  combined_data$Date_full - min(combined_data$Date_full)
)

# 组合因子
combined_data$TypeGroup <- interaction(combined_data$Type, combined_data$Group)

# ✅ 标准化环境变量（提升模型稳定性）
combined_data$Temperature_num <- scale(combined_data$Temperature_num)
combined_data$Humidness_num   <- scale(combined_data$Humidness_num)

############################################################
# 2. GAM模型（最终论文版本）
############################################################

gam_compare <- gam(
  log_richness ~
    s(Days, by = TypeGroup, k = 15) +
    s(Temperature_num, k = 8) +
    s(Humidness_num, k = 8) +
    Type * Group,
  data   = combined_data,
  method = "REML"
)

############################################################
# 3. 模型诊断 & 统计输出
############################################################

summary_gam <- summary(gam_compare)
print(summary_gam)

# 关键统计结果提取
cat("\n===== GAM Smooth Terms =====\n")
print(summary_gam$s.table)

cat("\n===== GAM Parametric Terms =====\n")
print(summary_gam$p.table)

# 模型检验
gam.check(gam_compare)

# ANOVA（增强论文说服力）
anova_gam <- anova(gam_compare, test = "F")
print(anova_gam)

############################################################
# 4. 方差分解（核心新增）
############################################################


# 建立线性模型
lm_model <- lm(
  log_richness ~ Group + Type + Days + Temperature_num + Humidness_num,
  data = combined_data
)

# 查看方差贡献（ANOVA分解）
anova_res <- anova(lm_model)
print(anova_res)

# 计算每个变量解释的方差比例
var_explained <- anova_res$`Sum Sq` / sum(anova_res$`Sum Sq`)
var_table <- data.frame(
  Factor = rownames(anova_res),
  Variance_Explained = round(var_explained * 100, 2)
)

cat("\n===== Variance Explained (%) =====\n")
print(var_table)

############################################################
# 5. 预测数据（用于画图）
############################################################

newdata <- expand.grid(
  Date_full = seq(min(combined_data$Date_full),
                  max(combined_data$Date_full),
                  by = "3 days"),
  Group = levels(combined_data$Group),
  Type  = levels(combined_data$Type)
) %>%
  mutate(
    Days      = as.numeric(Date_full - min(combined_data$Date_full)),
    TypeGroup = interaction(Type, Group),
    Temperature_num = 0,  # 标准化后均值=0
    Humidness_num   = 0
  )

pred <- predict(gam_compare, newdata = newdata, se.fit = TRUE)

newdata$fit   <- pred$fit
newdata$se    <- pred$se.fit
newdata$upper <- newdata$fit + 1.96 * newdata$se
newdata$lower <- newdata$fit - 1.96 * newdata$se

############################################################
# 6. 绘图（投稿级）
############################################################

col_bacteria <- "#2c7fb8"
col_fungi    <- "#d95f02"

p_compare <- ggplot() +
  
  geom_point(
    data = combined_data,
    aes(x = Date_full, y = log_richness, color = Type),
    alpha = 0.35, size = 1.2,
    position = position_jitter(width = 2, height = 0.02)
  ) +
  
  geom_ribbon(
    data = newdata,
    aes(x = Date_full, ymin = lower, ymax = upper, fill = Type),
    alpha = 0.18
  ) +
  
  geom_line(
    data = newdata,
    aes(x = Date_full, y = fit, color = Type),
    linewidth = 1.1
  ) +
  
  facet_wrap(~ Group, nrow = 1) +
  
  scale_color_manual(values = c("Bacteria" = col_bacteria,
                                "Fungi" = col_fungi)) +
  scale_fill_manual(values = c("Bacteria" = col_bacteria,
                               "Fungi" = col_fungi)) +
  
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b\n%Y"
  ) +
  
  labs(
    x = "Sampling date",
    y = expression(Log[10]~"richness"),
    color = NULL,
    fill  = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

print(p_compare)

ggsave("Final_GAM_Richness.pdf",
       p_compare,
       width = 16, height = 5, dpi = 300)

############################################################
# 7. 自动生成论文结果（关键！）
############################################################

cat("\n========== RESULTS SUMMARY ==========\n")

# 提取R²
cat("Adjusted R-sq:",
    round(summary_gam$r.sq, 3), "\n")

# 提取显著性smooth
sig_smooth <- summary_gam$s.table[
  summary_gam$s.table[,4] < 0.05, ]

cat("\nSignificant smooth terms:\n")
print(sig_smooth)

cat("\nInterpretation:\n")
cat("- Temporal dynamics are significant if s(Days) p < 0.05\n")
cat("- Temperature / Humidity effects evaluated via smooth terms\n")
cat("- Group differences via parametric terms\n")

############################################################
# END
############################################################



############################################################
# GAM_V5 — 细菌与真菌分别建模
# 解释变量：空间(Group) + 时间(Days) + 温度 + 湿度
############################################################

############################################################
# 0. 加载包
############################################################
library(mgcv)
library(dplyr)
library(ggplot2)
library(gratia)
library(vegan)
library(patchwork)   # 拼图用

############################################################
# 1. 数据准备
############################################################

bacteria_data$Type <- "Bacteria"
fungi_data$Type    <- "Fungi"

combined_data <- bind_rows(bacteria_data, fungi_data)

combined_data$Group <- as.factor(combined_data$Group)
combined_data$Type  <- as.factor(combined_data$Type)

# 时间变量（天数，从第一天起算）
combined_data$Days <- as.numeric(
  combined_data$Date_full - min(combined_data$Date_full)
)

# 标准化环境变量（提升模型数值稳定性，但保留实际变量）
combined_data$Temperature_num <- scale(combined_data$Temperature_num)
combined_data$Humidness_num   <- scale(combined_data$Humidness_num)

# 分别提取细菌和真菌数据集
bac_data <- combined_data %>% filter(Type == "Bacteria")
fun_data <- combined_data %>% filter(Type == "Fungi")

############################################################
# 2. 分别建立GAM模型
# 结构：空间(Group，参数项) + 时间(Days，by=Group平滑项)
#        + 温度(平滑项) + 湿度(平滑项)
############################################################

# ---- 细菌 GAM ----
gam_bac <- gam(
  log_richness ~
    Group +                              # 空间主效应（参数项）
    s(Days, by = Group, k = 15) +        # 各空间独立时间趋势
    s(Temperature_num, k = 8) +          # 温度非线性效应
    s(Humidness_num, k = 8),             # 湿度非线性效应
  data   = bac_data,
  method = "REML"
)

# ---- 真菌 GAM ----
gam_fun <- gam(
  log_richness ~
    Group +
    s(Days, by = Group, k = 15) +
    s(Temperature_num, k = 8) +
    s(Humidness_num, k = 8),
  data   = fun_data,
  method = "REML"
)

############################################################
# 3. 模型诊断 & 统计输出
############################################################

cat("\n\n========================================\n")
cat("===== BACTERIA GAM Summary =====\n")
cat("========================================\n")
summary_bac <- summary(gam_bac)
print(summary_bac)
gam.check(gam_bac)

cat("\n\n========================================\n")
cat("===== FUNGI GAM Summary =====\n")
cat("========================================\n")
summary_fun <- summary(gam_fun)
print(summary_fun)
gam.check(gam_fun)

############################################################
# 4. ANOVA（参数项显著性）
############################################################

cat("\n===== BACTERIA ANOVA =====\n")
print(anova(gam_bac))

cat("\n===== FUNGI ANOVA =====\n")
print(anova(gam_fun))

############################################################
# 5. 方差分解（线性近似，分别对细菌/真菌）
# 解释变量：Group + Days + Temperature + Humidity
# 不含 Type（细菌和真菌已分开）
############################################################

# ---- 细菌方差分解 ----
lm_bac <- lm(
  log_richness ~ Group + Days + Temperature_num + Humidness_num,
  data = bac_data
)
anova_bac <- anova(lm_bac)
var_bac <- data.frame(
  Factor = rownames(anova_bac),
  SS     = anova_bac$`Sum Sq`,
  Variance_Explained_pct = round(
    anova_bac$`Sum Sq` / sum(anova_bac$`Sum Sq`) * 100, 2)
)

cat("\n===== Bacteria — Variance Explained (%) =====\n")
print(var_bac)

# ---- 真菌方差分解 ----
lm_fun <- lm(
  log_richness ~ Group + Days + Temperature_num + Humidness_num,
  data = fun_data
)
anova_fun <- anova(lm_fun)
var_fun <- data.frame(
  Factor = rownames(anova_fun),
  SS     = anova_fun$`Sum Sq`,
  Variance_Explained_pct = round(
    anova_fun$`Sum Sq` / sum(anova_fun$`Sum Sq`) * 100, 2)
)

cat("\n===== Fungi — Variance Explained (%) =====\n")
print(var_fun)

############################################################
# 6. 预测数据（用于画图，温度和湿度固定在均值=0用于可视化）
############################################################

make_newdata <- function(data) {
  expand.grid(
    Date_full = seq(min(data$Date_full),
                    max(data$Date_full),
                    by = "3 days"),
    Group = levels(data$Group)
  ) %>%
    mutate(
      Days            = as.numeric(Date_full - min(data$Date_full)),
      Temperature_num = 0,   # 标准化均值=0，画图用
      Humidness_num   = 0
    )
}

newdata_bac <- make_newdata(bac_data)
newdata_fun <- make_newdata(fun_data)

pred_bac <- predict(gam_bac, newdata = newdata_bac, se.fit = TRUE)
pred_fun <- predict(gam_fun, newdata = newdata_fun, se.fit = TRUE)

newdata_bac <- newdata_bac %>%
  mutate(fit   = pred_bac$fit,
         se    = pred_bac$se.fit,
         upper = fit + 1.96 * se,
         lower = fit - 1.96 * se,
         Type  = "Bacteria")

newdata_fun <- newdata_fun %>%
  mutate(fit   = pred_fun$fit,
         se    = pred_fun$se.fit,
         upper = fit + 1.96 * se,
         lower = fit - 1.96 * se,
         Type  = "Fungi")

newdata_all <- bind_rows(newdata_bac, newdata_fun)

############################################################
# 7. 绘图（投稿级，细菌/真菌分面 × 空间分面）
############################################################

col_bacteria <- "#2c7fb8"
col_fungi    <- "#d95f02"

p_final <- ggplot() +
  
  geom_point(
    data = combined_data,
    aes(x = Date_full, y = log_richness, color = Type),
    alpha = 0.35, size = 1.2,
    position = position_jitter(width = 2, height = 0.02)
  ) +
  
  geom_ribbon(
    data = newdata_all,
    aes(x = Date_full, ymin = lower, ymax = upper, fill = Type),
    alpha = 0.18
  ) +
  
  geom_line(
    data = newdata_all,
    aes(x = Date_full, y = fit, color = Type),
    linewidth = 1.1
  ) +
  
  facet_grid(Type ~ Group) +          # 行=Type，列=空间，更清晰
  
  scale_color_manual(values = c("Bacteria" = col_bacteria,
                                "Fungi"    = col_fungi)) +
  scale_fill_manual(values  = c("Bacteria" = col_bacteria,
                                "Fungi"    = col_fungi)) +
  
  scale_x_date(
    date_breaks = "2 months",
    date_labels = "%b\n%Y"
  ) +
  
  labs(
    x     = "Sampling date",
    y     = expression(Log[10]~"richness"),
    color = NULL,
    fill  = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    strip.text       = element_text(face = "bold"),
    legend.position  = "top",
    axis.text.x      = element_text(angle = 90, vjust = 0.5)
  )

print(p_final)

ggsave("Final_GAM_Richness_V5.pdf",
       p_final,
       width = 16, height = 7, dpi = 300)

############################################################
# 8. 结果摘要输出
############################################################

cat("\n========== RESULTS SUMMARY ==========\n")

cat("\n--- Bacteria ---\n")
cat("Adjusted R-sq:", round(summary_bac$r.sq, 3), "\n")
cat("Deviance explained:", round(summary_bac$dev.expl * 100, 1), "%\n")

cat("\n--- Fungi ---\n")
cat("Adjusted R-sq:", round(summary_fun$r.sq, 3), "\n")
cat("Deviance explained:", round(summary_fun$dev.expl * 100, 1), "%\n")

cat("\n--- Significant smooth terms (p < 0.05) ---\n")
cat("Bacteria:\n")
print(summary_bac$s.table[summary_bac$s.table[, 4] < 0.05, ])

cat("Fungi:\n")
print(summary_fun$s.table[summary_fun$s.table[, 4] < 0.05, ])

cat("\n--- Variance Explained Summary ---\n")
cat("Bacteria:\n"); print(var_bac)
cat("Fungi:\n");    print(var_fun)


####################3LMM
# ============================================================
# VERSION 1: LMM + Fourier Terms
# 傅里叶项线性混合效应模型
# ============================================================

library(dplyr)
library(ggplot2)
library(lubridate)
library(lme4)
library(lmerTest)   # 提供 summary 中的 p 值
library(patchwork)

# ── 输出目录 ────────────────────────────────────────────────
outDir <- "./AlphaDiversity_LMM_Fourier/"
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# ── 配色 ────────────────────────────────────────────────────
col_bac  <- "#2c7fb8"
col_fun  <- "#d95f02"
col_vals <- c("Bacteria" = col_bac, "Fungi" = col_fun)

# ============================================================
# 1. 数据准备
# ============================================================

# 合并细菌数据
bacteria_data <- metadata %>%
  inner_join(bacteria_alpha, by = "Sample") %>%
  mutate(
    log_richness = log(richness),
    log_shannon  = log(Shannon),
    Type         = "Bacteria"
  )

# 合并真菌数据
fungi_data <- metadata %>%
  inner_join(fungi_alpha, by = "Sample") %>%
  mutate(
    log_richness = log(richness),
    log_shannon  = log(Shannon),
    Type         = "Fungi"
  )

# 合并数据
combined_data <- bind_rows(bacteria_data, fungi_data) %>%
  mutate(
    Group    = as.factor(Group),
    Type     = as.factor(Type),
    # 从数据最早日期起算的天数（更稳定）
    Days     = as.numeric(Date_full - min(Date_full, na.rm = TRUE)),
    # 标准化环境变量
    Temp_sc  = as.numeric(scale(Temperature_num)),
    Humid_sc = as.numeric(scale(Humidness_num)),
    # 年内天数（用于季节性）
    DOY      = as.integer(format(Date_full, "%j"))
  )

# 分别提取
bac_data <- combined_data %>% filter(Type == "Bacteria")
fun_data <- combined_data %>% filter(Type == "Fungi")

# ============================================================
# 2. 傅里叶项 LMM 模型
# 周期 = 365天  → ω = 2π/365
# 包含年周期（ω）和半年周期（2ω）两组傅里叶项
# ============================================================

omega <- 2 * pi / 365

# 添加傅里叶项到数据
add_fourier <- function(df) {
  df %>% mutate(
    sin1 = sin(omega * DOY),
    cos1 = cos(omega * DOY),
    sin2 = sin(2 * omega * DOY),
    cos2 = cos(2 * omega * DOY)
  )
}

bac_data <- add_fourier(bac_data)
fun_data <- add_fourier(fun_data)
combined_data <- add_fourier(combined_data)

# ---- 细菌 LMM（随机效应：采样地点/SampleSite 或 Group）
# Temperature 和 Humidity 作为随机截距
lmm_bac <- lmer(
  log_richness ~
    sin1 + cos1 + sin2 + cos2 +   # 年周期 + 半年周期
    Group +                         # 空间主效应
    Temp_sc + Humid_sc +            # 环境固定效应
    (1 | Group),                    # 随机截距（组内重复）
  data    = bac_data,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

# ---- 真菌 LMM
lmm_fun <- lmer(
  log_richness ~
    sin1 + cos1 + sin2 + cos2 +
    Group +
    Temp_sc + Humid_sc +
    (1 | Group),
  data    = fun_data,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

# ---- 合并模型（Type 作为固定效应）
lmm_combined <- lmer(
  log_richness ~
    sin1 + cos1 + sin2 + cos2 +
    Type + Group + Type:Group +
    Temp_sc + Humid_sc +
    (1 | Group),
  data    = combined_data,
  REML    = TRUE,
  control = lmerControl(optimizer = "bobyqa")
)

# ── 统计输出 ─────────────────────────────────────────────
cat("\n========== LMM Bacteria Summary ==========\n")
print(summary(lmm_bac))

cat("\n========== LMM Fungi Summary ==========\n")
print(summary(lmm_fun))

cat("\n========== LMM Combined Summary ==========\n")
print(summary(lmm_combined))

# ANOVA 显著性
cat("\n--- Bacteria ANOVA ---\n")
print(anova(lmm_bac))
cat("\n--- Fungi ANOVA ---\n")
print(anova(lmm_fun))
cat("\n--- Combined ANOVA ---\n")
print(anova(lmm_combined))

# ── 计算 R² (marginal & conditional) ────────────────────
library(MuMIn)
cat("\n--- R² (Bacteria) ---\n"); print(r.squaredGLMM(lmm_bac))
cat("\n--- R² (Fungi) ---\n");    print(r.squaredGLMM(lmm_fun))
cat("\n--- R² (Combined) ---\n"); print(r.squaredGLMM(lmm_combined))

# ============================================================
# 3. 生成预测数据（用于画图）
# ============================================================

# 生成平滑DOY序列用于季节性预测
make_pred_lmm <- function(model, data, type_label) {
  # 用实际采样日期序列预测（展示真实时间轴）
  date_seq <- seq(min(data$Date_full), max(data$Date_full), by = "3 days")
  doy_seq  <- as.integer(format(date_seq, "%j"))
  
  newdf <- expand.grid(
    Group = levels(data$Group)
  ) %>%
    crossing(
      tibble(
        Date_full = date_seq,
        DOY       = doy_seq
      )
    ) %>%
    mutate(
      Days     = as.numeric(Date_full - min(data$Date_full)),
      sin1     = sin(omega * DOY),
      cos1     = cos(omega * DOY),
      sin2     = sin(2 * omega * DOY),
      cos2     = cos(2 * omega * DOY),
      Temp_sc  = 0,   # 固定在均值
      Humid_sc = 0,
      Type     = type_label
    )
  
  # re.form = NA：只用固定效应预测（适合画总体趋势）
  pred <- predict(model, newdata = newdf, re.form = NA, se.fit = FALSE)
  # lme4 不直接给 se，用 bootMer 或手动计算；这里用近似方法
  newdf$fit <- pred
  # 近似 SE（用模型残差SD）
  resid_sd   <- sigma(model)
  newdf$se   <- resid_sd / sqrt(10)   # 近似值
  newdf$upper <- newdf$fit + 1.96 * newdf$se
  newdf$lower <- newdf$fit - 1.96 * newdf$se
  newdf
}

# 使用 tidyr::crossing
library(tidyr)

pred_bac_lmm <- make_pred_lmm(lmm_bac, bac_data, "Bacteria")
pred_fun_lmm <- make_pred_lmm(lmm_fun, fun_data, "Fungi")
pred_all_lmm <- bind_rows(pred_bac_lmm, pred_fun_lmm) %>%
  mutate(Type = as.factor(Type),
         Group = as.factor(Group))

# ============================================================
# 4. 提取傅里叶项振幅（用于图内统计标注）
# ============================================================

extract_fourier_sig <- function(model, type_label) {
  coef_df <- as.data.frame(summary(model)$coefficients)
  coef_df$term <- rownames(coef_df)
  fourier_terms <- coef_df %>%
    filter(term %in% c("sin1","cos1","sin2","cos2")) %>%
    mutate(
      sig = case_when(
        `Pr(>|t|)` < 0.001 ~ "***",
        `Pr(>|t|)` < 0.01  ~ "**",
        `Pr(>|t|)` < 0.05  ~ "*",
        `Pr(>|t|)` < 0.10  ~ ".",
        TRUE               ~ "ns"
      )
    )
  # 年周期振幅
  A1 <- sqrt(coef_df$Estimate[coef_df$term=="sin1"]^2 +
               coef_df$Estimate[coef_df$term=="cos1"]^2)
  A2 <- sqrt(coef_df$Estimate[coef_df$term=="sin2"]^2 +
               coef_df$Estimate[coef_df$term=="cos2"]^2)
  cat(type_label, "— Annual amplitude:", round(A1,3),
      " | Semi-annual amplitude:", round(A2,3), "\n")
  fourier_terms
}

cat("\n--- Fourier term significance ---\n")
ft_bac <- extract_fourier_sig(lmm_bac, "Bacteria")
ft_fun <- extract_fourier_sig(lmm_fun, "Fungi")
print(ft_bac); print(ft_fun)

# ============================================================
# 5. 画图（4个分面，细菌+真菌两条线）
# ============================================================

# 提取 R² 用于图注
r2_bac <- round(r.squaredGLMM(lmm_bac)[1, "R2m"], 3)
r2_fun <- round(r.squaredGLMM(lmm_fun)[1, "R2m"], 3)

# 构建图内标注（放在每个分面右下角）
annot_lmm <- data.frame(
  Group = levels(combined_data$Group)
) %>%
  mutate(
    label_bac = paste0("Bac R²m=", r2_bac),
    label_fun = paste0("Fun R²m=", r2_fun),
    x_pos     = max(combined_data$Date_full) - 20,
    y_bac     = max(combined_data$log_richness, na.rm=TRUE) * 0.97,
    y_fun     = max(combined_data$log_richness, na.rm=TRUE) * 0.90
  )

p_lmm <- ggplot() +
  
  # 原始散点
  geom_point(
    data = combined_data,
    aes(x = Date_full, y = log_richness, color = Type),
    alpha = 0.30, size = 0.9, shape = 16,
    position = position_jitter(width = 2, height = 0.02)
  ) +
  
  # 置信带
  geom_ribbon(
    data = pred_all_lmm,
    aes(x = Date_full, ymin = lower, ymax = upper, fill = Type),
    alpha = 0.15
  ) +
  
  # 拟合曲线
  geom_line(
    data = pred_all_lmm,
    aes(x = Date_full, y = fit, color = Type),
    linewidth = 1.0
  ) +
  
  facet_wrap(~ Group, nrow = 1) +
  
  scale_color_manual(values = col_vals) +
  scale_fill_manual(values  = col_vals) +
  
  scale_x_date(
    breaks = seq(min(combined_data$Date_full),
                 max(combined_data$Date_full),
                 by = "2 months"),
    labels = function(x) format(x, "%b\n%Y"),
    expand = c(0.02, 0)
  ) +
  
  labs(
    x     = "Sampling date",
    y     = expression(Log[10]~"(observed richness)"),
    color = "Microbial type",
    fill  = "Microbial type",
    title = "LMM + Fourier terms (annual + semi-annual seasonality)"
  ) +
  
  theme_classic(base_size = 10) +
  theme(
    aspect.ratio     = 1,
    plot.title       = element_text(size = 10, hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text.x      = element_text(size = 7.5, color = "black"),
    axis.text.y      = element_text(size = 8.5, color = "black"),
    axis.title       = element_text(size = 9.5),
    axis.line        = element_line(linewidth = 0.35),
    axis.ticks       = element_line(linewidth = 0.35),
    legend.position  = "top",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 9),
    legend.key.size  = unit(0.35, "cm"),
    panel.spacing    = unit(0.7, "cm"),
    plot.margin      = margin(6, 6, 6, 6)
  )

p_lmm

ggsave(paste0(outDir, "Fig_LMM_Fourier_Richness.pdf"),
       p_lmm, width = 22, height = 7, units = "cm", dpi = 300)
ggsave(paste0(outDir, "Fig_LMM_Fourier_Richness.png"),
       p_lmm, width = 22, height = 7, units = "cm", dpi = 300, bg = "white")

# ============================================================
# 6. 季节性曲线图（补充）：DOY 1–365 的预测
# ============================================================

make_doy_pred <- function(model, type_label) {
  doy_df <- expand.grid(
    DOY   = 1:365,
    Group = levels(combined_data$Group)
  ) %>%
    mutate(
      sin1     = sin(omega * DOY),
      cos1     = cos(omega * DOY),
      sin2     = sin(2 * omega * DOY),
      cos2     = cos(2 * omega * DOY),
      Days     = DOY,
      Temp_sc  = 0,
      Humid_sc = 0,
      Type     = type_label,
      Date_full = as.Date(DOY - 1, origin = "2023-01-01")
    )
  doy_df$fit <- predict(model, newdata = doy_df, re.form = NA)
  doy_df
}

doy_bac <- make_doy_pred(lmm_bac, "Bacteria")
doy_fun <- make_doy_pred(lmm_fun, "Fungi")
doy_all <- bind_rows(doy_bac, doy_fun) %>%
  mutate(Type = as.factor(Type), Group = as.factor(Group))

p_doy_lmm <- ggplot(doy_all, aes(x = DOY, y = fit, color = Type)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ Group, nrow = 1) +
  scale_color_manual(values = col_vals) +
  scale_x_continuous(
    breaks = c(1, 91, 182, 274, 365),
    labels = c("Jan", "Apr", "Jul", "Oct", "Dec")
  ) +
  labs(
    x     = "Day of year",
    y     = expression("Predicted"~Log[10]~"(richness)"),
    color = "Microbial type",
    title = "LMM Fourier — estimated seasonal pattern (annual + semi-annual)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    aspect.ratio     = 1,
    plot.title       = element_text(size = 9, hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = 8, color = "black"),
    legend.position  = "top",
    panel.spacing    = unit(0.7, "cm")
  )

p_doy_lmm

ggsave(paste0(outDir, "Fig_LMM_Seasonal_DOY.pdf"),
       p_doy_lmm, width = 22, height = 7, units = "cm", dpi = 300)

cat("\n✓ Version 1 (LMM + Fourier) complete.\n")


############################
###########################
###################################GAMM
# ============================================================
# GAMM 完整版（顶刊投稿）
# 模型：Group固定效应 + 全局平滑 + 组间偏差平滑
# 输出：四张独立图（细菌/真菌 × 时间趋势/季节性）
# ============================================================

# ── 0. 依赖包 ──────────────────────────────────────────────
library(mgcv)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(stringr)

col_bac  <- "#2c7fb8"
col_fun  <- "#d95f02"

# ============================================================
# 1. 数据整理
# ============================================================

metadata <- metadata %>%
  mutate(
    Date_full = as.Date(as.character(Date_full)),
    Days      = as.numeric(Date_full - min(Date_full, na.rm = TRUE)),
    DOY       = as.integer(format(Date_full, "%j")),
    Temp_sc   = as.numeric(scale(Temperature_num)),
    Humid_sc  = as.numeric(scale(Humidness_num)),
    Group     = factor(Group)
  )

days_origin <- min(metadata$Date_full, na.rm = TRUE)
cat("Days 起点：", as.character(days_origin), "\n")

data_bac <- bacteria_alpha %>%
  inner_join(metadata, by = "Sample") %>%
  filter(!is.na(richness), richness > 0) %>%
  mutate(log_richness = log(richness))

data_fun <- fungi_alpha %>%
  inner_join(metadata, by = "Sample") %>%
  filter(!is.na(richness), richness > 0) %>%
  mutate(log_richness = log(richness))

all_groups <- union(levels(data_bac$Group), levels(data_fun$Group))
data_bac$Group <- factor(data_bac$Group, levels = all_groups)
data_fun$Group <- factor(data_fun$Group, levels = all_groups)
group_levels   <- levels(data_bac$Group)

cat("[细菌] n =", nrow(data_bac), "\n")
cat("[真菌] n =", nrow(data_fun), "\n")

# ============================================================
# 2. 建模：细菌
# ============================================================

cat("\n===== 拟合细菌 GAMM =====\n")

model_bac <- gam(
  log_richness ~
    Group +                                          # 空间类型主效应（固定参数项）
    Temp_sc +                                        # 温度（线性）
    Humid_sc +                                       # 湿度（线性）
    s(Days, k = 10) +                                # 全局长期趋势
    s(Days, by = Group, k = 6, m = 1) +              # 各组趋势偏差
    s(DOY, bs = "cc", k = 12) +                      # 全局季节性
    s(DOY, by = Group, bs = "cc", k = 8, m = 1),    # 各组季节偏差
  data   = data_bac,
  method = "REML",
  family = gaussian()
)

smry_bac <- summary(model_bac)
cat("\n[细菌] 模型摘要：\n"); print(smry_bac)
cat("\n[细菌] gam.check：\n"); gam.check(model_bac)

# ============================================================
# 3. 建模：真菌
# ============================================================

cat("\n===== 拟合真菌 GAMM =====\n")

model_fun <- gam(
  log_richness ~
    Group +
    Temp_sc +
    Humid_sc +
    s(Days, k = 10) +
    s(Days, by = Group, k = 6, m = 1) +
    s(DOY, bs = "cc", k = 12) +
    s(DOY, by = Group, bs = "cc", k = 8, m = 1),
  data   = data_fun,
  method = "REML",
  family = gaussian()
)

smry_fun <- summary(model_fun)
cat("\n[真菌] 模型摘要：\n"); print(smry_fun)
cat("\n[真菌] gam.check：\n"); gam.check(model_fun)

# ============================================================
# 4. 方差分解（量化各因素贡献，支撑"空间类型为主导"的结论）
# ============================================================

compute_var_decomp_new <- function(full_model, data) {
  r2_full <- summary(full_model)$r.sq
  
  # 去掉 Group 主效应
  m_no_group <- update(full_model,
                       . ~ Temp_sc + Humid_sc +
                         s(Days, k = 10) +
                         s(Days, by = Group, k = 6, m = 1) +
                         s(DOY, bs = "cc", k = 12) +
                         s(DOY, by = Group, bs = "cc", k = 8, m = 1))
  
  # 去掉全局趋势
  m_no_trend <- update(full_model,
                       . ~ Group + Temp_sc + Humid_sc +
                         s(DOY, bs = "cc", k = 12) +
                         s(DOY, by = Group, bs = "cc", k = 8, m = 1))
  
  # 去掉全局季节
  m_no_seas <- update(full_model,
                      . ~ Group + Temp_sc + Humid_sc +
                        s(Days, k = 10) +
                        s(Days, by = Group, k = 6, m = 1))
  
  # 去掉组间偏差（趋势+季节）
  m_no_dev <- update(full_model,
                     . ~ Group + Temp_sc + Humid_sc +
                       s(Days, k = 10) +
                       s(DOY, bs = "cc", k = 12))
  
  # 去掉温湿度
  m_no_env <- update(full_model,
                     . ~ Group +
                       s(Days, k = 10) +
                       s(Days, by = Group, k = 6, m = 1) +
                       s(DOY, bs = "cc", k = 12) +
                       s(DOY, by = Group, bs = "cc", k = 8, m = 1))
  
  data.frame(
    Component = c(
      "Spatial type (Group)",
      "Global long-term trend",
      "Global seasonal cycle",
      "Environment-specific deviation",
      "Temperature + Humidity",
      "Residual/Unexplained"
    ),
    Delta_R2 = round(c(
      r2_full - summary(m_no_group)$r.sq,
      r2_full - summary(m_no_trend)$r.sq,
      r2_full - summary(m_no_seas)$r.sq,
      r2_full - summary(m_no_dev)$r.sq,
      r2_full - summary(m_no_env)$r.sq,
      1 - r2_full
    ), 3)
  )
}

cat("\n===== 细菌 方差分解 =====\n")
vd_bac <- compute_var_decomp_new(model_bac, data_bac)
print(vd_bac)

cat("\n===== 真菌 方差分解 =====\n")
vd_fun <- compute_var_decomp_new(model_fun, data_fun)
print(vd_fun)

# ============================================================
# 5. 统计摘要输出（用于 Results 撰写）
# ============================================================

print_model_summary <- function(smry, domain, n) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("  ", domain, " GAMM 统计摘要\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")
  
  cat("\n[PARAMETRIC TERMS]\n")
  print(round(smry$p.table, 4))
  
  cat("\n[SMOOTH TERMS]\n")
  print(round(smry$s.table, 4))
  
  cat("\nAdj. R2 =", round(smry$r.sq, 3))
  cat("  |  Dev. explained =", round(smry$dev.expl * 100, 1), "%")
  cat("  |  n =", n, "\n")
}

print_model_summary(smry_bac, "细菌", nrow(data_bac))
print_model_summary(smry_fun, "真菌", nrow(data_fun))

# ============================================================
# 6. 生成预测数据
# ============================================================

date_seq <- seq(min(metadata$Date_full, na.rm = TRUE),
                max(metadata$Date_full, na.rm = TRUE),
                by = "3 days")

base_nd <- expand.grid(
  Date_full = date_seq,
  Group     = group_levels
) %>%
  mutate(
    Days     = as.numeric(Date_full - days_origin),
    DOY      = as.integer(format(Date_full, "%j")),
    Temp_sc  = 0,
    Humid_sc = 0
  )

# 细菌时间预测
p_bac <- predict(model_bac, newdata = base_nd, se.fit = TRUE)
pred_bac_time <- base_nd %>%
  mutate(fit = p_bac$fit, se = p_bac$se.fit,
         upper = fit + 1.96*se, lower = fit - 1.96*se)

# 真菌时间预测
p_fun <- predict(model_fun, newdata = base_nd, se.fit = TRUE)
pred_fun_time <- base_nd %>%
  mutate(fit = p_fun$fit, se = p_fun$se.fit,
         upper = fit + 1.96*se, lower = fit - 1.96*se)

# DOY 预测网格
doy_nd_bac <- expand.grid(DOY = 1:365, Group = group_levels) %>%
  mutate(Days = median(data_bac$Days), Temp_sc = 0, Humid_sc = 0,
         Date_full = as.Date(DOY - 1, origin = "2023-01-01"))

doy_nd_fun <- expand.grid(DOY = 1:365, Group = group_levels) %>%
  mutate(Days = median(data_fun$Days), Temp_sc = 0, Humid_sc = 0,
         Date_full = as.Date(DOY - 1, origin = "2023-01-01"))

# 细菌季节预测
p_bac_doy <- predict(model_bac, newdata = doy_nd_bac, se.fit = TRUE)
pred_bac_doy <- doy_nd_bac %>%
  mutate(fit = p_bac_doy$fit, se = p_bac_doy$se.fit,
         upper = fit + 1.96*se, lower = fit - 1.96*se)

# 真菌季节预测
p_fun_doy <- predict(model_fun, newdata = doy_nd_fun, se.fit = TRUE)
pred_fun_doy <- doy_nd_fun %>%
  mutate(fit = p_fun_doy$fit, se = p_fun_doy$se.fit,
         upper = fit + 1.96*se, lower = fit - 1.96*se)

# ============================================================
# 7. 构建图内统计标注函数
# ============================================================

# 提取参数项标注（空间类型效应，用于时间图）
make_param_annot <- function(smry, group_levels, raw_data, x_date_end) {
  pt <- as.data.frame(smry$p.table) %>%
    rownames_to_column("term") %>%
    filter(grepl("^Group", term) | term == "(Intercept)") %>%
    mutate(
      Group = case_when(
        term == "(Intercept)" ~ group_levels[1],
        TRUE ~ str_remove(term, "^Group")
      ),
      sig = case_when(
        `Pr(>|t|)` < 0.001 ~ "***",
        `Pr(>|t|)` < 0.01  ~ "**",
        `Pr(>|t|)` < 0.05  ~ "*",
        `Pr(>|t|)` < 0.10  ~ ".",
        TRUE                ~ "ns"
      ),
      # 截距（参照组）无需标注 beta，其他组标注偏差量
      param_lbl = ifelse(
        term == "(Intercept)",
        paste0("Ref. group\nμ = ", round(Estimate, 2)),
        paste0("β = ", round(Estimate, 2), sig)
      )
    ) %>%
    select(Group, param_lbl)
  
  y_pos <- raw_data %>%
    group_by(Group) %>%
    summarise(y_max = quantile(log_richness, 0.98, na.rm = TRUE), .groups = "drop")
  
  pt %>%
    mutate(Group = factor(Group, levels = group_levels)) %>%
    left_join(y_pos, by = "Group") %>%
    mutate(x_pos = x_date_end - 15)
}

# 提取平滑项标注（时间趋势，用于时间图）
make_smooth_annot_time <- function(smry, group_levels, raw_data, x_date_end) {
  st <- as.data.frame(smry$s.table) %>%
    rownames_to_column("term") %>%
    mutate(
      sig = case_when(
        `p-value` < 0.001 ~ "***",
        `p-value` < 0.01  ~ "**",
        `p-value` < 0.05  ~ "*",
        `p-value` < 0.10  ~ ".",
        TRUE               ~ "ns"
      ),
      edf_r = round(edf, 2)
    )
  
  # 全局趋势（所有分面一样）
  g_trend_row <- st %>% filter(term == "s(Days)")
  g_trend_lbl <- if (nrow(g_trend_row) > 0)
    paste0("Global trend: edf=", g_trend_row$edf_r[1], g_trend_row$sig[1])
  else "Global trend: ns"
  
  # 各组趋势偏差
  dev_rows <- st %>%
    filter(str_detect(term, "s\\(Days\\):Group")) %>%
    mutate(Group = str_extract(term, paste(group_levels, collapse = "|")),
           dev_lbl = ifelse(`p-value` < 0.10,
                            paste0("Δtrend: edf=", edf_r, sig),
                            ""))
  
  y_pos <- raw_data %>%
    group_by(Group) %>%
    summarise(y_max = quantile(log_richness, 0.98, na.rm = TRUE), .groups = "drop")
  
  data.frame(Group = group_levels) %>%
    left_join(dev_rows %>% select(Group, dev_lbl), by = "Group") %>%
    mutate(
      dev_lbl  = replace_na(dev_lbl, ""),
      full_lbl = ifelse(dev_lbl == "",
                        g_trend_lbl,
                        paste0(g_trend_lbl, "\n", dev_lbl)),
      Group    = factor(Group, levels = group_levels)
    ) %>%
    left_join(y_pos, by = "Group") %>%
    mutate(x_pos = x_date_end - 15,
           y_pos = y_max * 0.88)
}

# 提取季节性标注（用于 DOY 图）
make_smooth_annot_seas <- function(smry, group_levels, raw_data) {
  st <- as.data.frame(smry$s.table) %>%
    rownames_to_column("term") %>%
    mutate(
      sig = case_when(
        `p-value` < 0.001 ~ "***",
        `p-value` < 0.01  ~ "**",
        `p-value` < 0.05  ~ "*",
        `p-value` < 0.10  ~ ".",
        TRUE               ~ "ns"
      ),
      edf_r = round(edf, 2)
    )
  
  g_seas_row <- st %>% filter(term == "s(DOY)")
  g_seas_lbl <- if (nrow(g_seas_row) > 0)
    paste0("Global seas.: edf=", g_seas_row$edf_r[1], g_seas_row$sig[1])
  else "Global seas.: ns"
  
  dev_rows <- st %>%
    filter(str_detect(term, "s\\(DOY\\):Group")) %>%
    mutate(Group = str_extract(term, paste(group_levels, collapse = "|")),
           dev_lbl = ifelse(`p-value` < 0.10,
                            paste0("Δseas.: edf=", edf_r, sig),
                            ""))
  
  y_pos <- raw_data %>%
    group_by(Group) %>%
    summarise(y_max = quantile(log_richness, 0.98, na.rm = TRUE), .groups = "drop")
  
  data.frame(Group = group_levels) %>%
    left_join(dev_rows %>% select(Group, dev_lbl), by = "Group") %>%
    mutate(
      dev_lbl  = replace_na(dev_lbl, ""),
      full_lbl = ifelse(dev_lbl == "",
                        g_seas_lbl,
                        paste0(g_seas_lbl, "\n", dev_lbl)),
      Group    = factor(Group, levels = group_levels)
    ) %>%
    left_join(y_pos, by = "Group") %>%
    mutate(x_pos = 355, y_pos = y_max * 0.88)
}

# ============================================================
# 8. 生成标注数据
# ============================================================

x_end <- max(metadata$Date_full, na.rm = TRUE)

# 细菌标注
annot_bac_param <- make_param_annot(smry_bac, group_levels, data_bac, x_end)
annot_bac_time  <- make_smooth_annot_time(smry_bac, group_levels, data_bac, x_end)
annot_bac_seas  <- make_smooth_annot_seas(smry_bac, group_levels, data_bac)

# 真菌标注
annot_fun_param <- make_param_annot(smry_fun, group_levels, data_fun, x_end)
annot_fun_time  <- make_smooth_annot_time(smry_fun, group_levels, data_fun, x_end)
annot_fun_seas  <- make_smooth_annot_seas(smry_fun, group_levels, data_fun)

# ============================================================
# 9. 统一主题
# ============================================================

my_theme <- function(base = 10) {
  theme_classic(base_size = base) +
    theme(
      aspect.ratio     = 1,
      strip.background = element_rect(fill = "grey96", color = "grey60",
                                      linewidth = 0.4),
      strip.text       = element_text(face = "bold", size = base + 1),
      axis.text.x      = element_text(size = base - 2, color = "black"),
      axis.text.y      = element_text(size = base - 1, color = "black"),
      axis.title       = element_text(size = base, face = "bold"),
      axis.line        = element_line(linewidth = 0.35),
      axis.ticks       = element_line(linewidth = 0.35),
      legend.position  = "none",
      panel.spacing    = unit(0.7, "cm"),
      plot.title       = element_text(face = "bold", size = base + 2,
                                      hjust = 0.5, margin = margin(b = 6)),
      plot.subtitle    = element_text(size = base - 1, hjust = 0.5,
                                      color = "grey40", margin = margin(b = 4)),
      plot.margin      = margin(8, 10, 8, 8)
    )
}

date_breaks <- seq(min(metadata$Date_full, na.rm = TRUE),
                   max(metadata$Date_full, na.rm = TRUE),
                   by = "2 months")
doy_breaks  <- c(1, 60, 121, 182, 244, 305, 365)
doy_labels  <- c("Jan","Mar","May","Jul","Sep","Nov","Dec")

# ============================================================
# 图1：细菌 — 时间趋势
# ============================================================

p1 <- ggplot() +
  # 原始散点
  geom_point(
    data = data_bac,
    aes(x = Date_full, y = log_richness),
    color = col_bac, alpha = 0.30, size = 0.85, shape = 16,
    position = position_jitter(width = 2, height = 0.02)
  ) +
  # 95% CI
  geom_ribbon(
    data = pred_bac_time,
    aes(x = Date_full, ymin = lower, ymax = upper),
    fill = col_bac, alpha = 0.18
  ) +
  # 拟合曲线
  geom_line(
    data = pred_bac_time,
    aes(x = Date_full, y = fit),
    color = col_bac, linewidth = 0.95
  ) +
  # 空间类型效应标注（β值）
  geom_text(
    data = annot_bac_param,
    aes(x = x_pos, y = y_max * 1.00, label = param_lbl),
    color = col_bac, size = 2.2, hjust = 1, lineheight = 1.05
  ) +
  # 时间趋势平滑统计标注
  geom_text(
    data = annot_bac_time,
    aes(x = x_pos, y = y_pos, label = full_lbl),
    color = "grey30", size = 1.9, hjust = 1, lineheight = 1.05
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_x_date(
    breaks = date_breaks,
    labels = function(x) format(x, "%b\n%Y"),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    name   = expression(Log[e]~"(bacterial richness)"),
    expand = c(0.12, 0)
  ) +
  labs(
    x        = "Sampling date",
    title    = "Bacteria — Temporal Dynamics",
    subtitle = paste0("Adj. R² = ", round(smry_bac$r.sq, 3),
                      "  |  Dev. expl. = ",
                      round(smry_bac$dev.expl * 100, 1), "%",
                      "  |  n = ", nrow(data_bac))
  ) +
  my_theme()

p1
ggsave("Fig1_Bacteria_Temporal.pdf",  p1, width = 24, height = 8, units = "cm", dpi = 300)
ggsave("Fig1_Bacteria_Temporal.png",  p1, width = 24, height = 8, units = "cm", dpi = 300, bg = "white")

# ============================================================
# 图2：细菌 — 季节性
# ============================================================

p2 <- ggplot() +
  geom_ribbon(
    data = pred_bac_doy,
    aes(x = DOY, ymin = lower, ymax = upper),
    fill = col_bac, alpha = 0.18
  ) +
  geom_line(
    data = pred_bac_doy,
    aes(x = DOY, y = fit),
    color = col_bac, linewidth = 0.95
  ) +
  geom_text(
    data = annot_bac_seas,
    aes(x = x_pos, y = y_pos, label = full_lbl),
    color = "grey30", size = 1.9, hjust = 1, lineheight = 1.05
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_x_continuous(
    breaks = doy_breaks, labels = doy_labels,
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    name   = expression("Predicted"~Log[e]~"(bacterial richness)"),
    expand = c(0.12, 0)
  ) +
  labs(
    x        = "Day of year",
    title    = "Bacteria — Seasonal Dynamics",
    subtitle = paste0("Cyclic cubic spline on DOY  |  n = ", nrow(data_bac))
  ) +
  my_theme()

p2
ggsave("Fig2_Bacteria_Seasonal.pdf", p2, width = 24, height = 8, units = "cm", dpi = 300)
ggsave("Fig2_Bacteria_Seasonal.png", p2, width = 24, height = 8, units = "cm", dpi = 300, bg = "white")

# ============================================================
# 图3：真菌 — 时间趋势
# ============================================================

p3 <- ggplot() +
  geom_point(
    data = data_fun,
    aes(x = Date_full, y = log_richness),
    color = col_fun, alpha = 0.30, size = 0.85, shape = 16,
    position = position_jitter(width = 2, height = 0.02)
  ) +
  geom_ribbon(
    data = pred_fun_time,
    aes(x = Date_full, ymin = lower, ymax = upper),
    fill = col_fun, alpha = 0.18
  ) +
  geom_line(
    data = pred_fun_time,
    aes(x = Date_full, y = fit),
    color = col_fun, linewidth = 0.95
  ) +
  geom_text(
    data = annot_fun_param,
    aes(x = x_pos, y = y_max * 1.00, label = param_lbl),
    color = col_fun, size = 2.2, hjust = 1, lineheight = 1.05
  ) +
  geom_text(
    data = annot_fun_time,
    aes(x = x_pos, y = y_pos, label = full_lbl),
    color = "grey30", size = 1.9, hjust = 1, lineheight = 1.05
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_x_date(
    breaks = date_breaks,
    labels = function(x) format(x, "%b\n%Y"),
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    name   = expression(Log[e]~"(fungal richness)"),
    expand = c(0.12, 0)
  ) +
  labs(
    x        = "Sampling date",
    title    = "Fungi — Temporal Dynamics",
    subtitle = paste0("Adj. R² = ", round(smry_fun$r.sq, 3),
                      "  |  Dev. expl. = ",
                      round(smry_fun$dev.expl * 100, 1), "%",
                      "  |  n = ", nrow(data_fun))
  ) +
  my_theme()

p3
ggsave("Fig3_Fungi_Temporal.pdf",  p3, width = 24, height = 8, units = "cm", dpi = 300)
ggsave("Fig3_Fungi_Temporal.png",  p3, width = 24, height = 8, units = "cm", dpi = 300, bg = "white")

# ============================================================
# 图4：真菌 — 季节性
# ============================================================

p4 <- ggplot() +
  geom_ribbon(
    data = pred_fun_doy,
    aes(x = DOY, ymin = lower, ymax = upper),
    fill = col_fun, alpha = 0.18
  ) +
  geom_line(
    data = pred_fun_doy,
    aes(x = DOY, y = fit),
    color = col_fun, linewidth = 0.95
  ) +
  geom_text(
    data = annot_fun_seas,
    aes(x = x_pos, y = y_pos, label = full_lbl),
    color = "grey30", size = 1.9, hjust = 1, lineheight = 1.05
  ) +
  facet_wrap(~ Group, nrow = 1) +
  scale_x_continuous(
    breaks = doy_breaks, labels = doy_labels,
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    name   = expression("Predicted"~Log[e]~"(fungal richness)"),
    expand = c(0.12, 0)
  ) +
  labs(
    x        = "Day of year",
    title    = "Fungi — Seasonal Dynamics",
    subtitle = paste0("Cyclic cubic spline on DOY  |  n = ", nrow(data_fun))
  ) +
  my_theme()

p4
ggsave("Fig4_Fungi_Seasonal.pdf", p4, width = 24, height = 8, units = "cm", dpi = 300)
ggsave("Fig4_Fungi_Seasonal.png", p4, width = 24, height = 8, units = "cm", dpi = 300, bg = "white")

# ============================================================
# 10. 方差分解结果整合输出
# ============================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("  方差分解汇总（用于 Results 撰写）\n")
cat(strrep("=", 60), "\n\n", sep = "")

cat("[细菌]\n")
print(vd_bac)
cat("\n[真菌]\n")
print(vd_fun)