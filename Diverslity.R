# =========================================
# 载入必要的包
# =========================================
library(ggplot2)
library(dplyr)
library(ggpubr)
library(readr)
setwd("/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro")
# =========================
# 1. 加载包
# =========================
library(tidyverse)
library(ggpubr)
library(rstatix)

# =========================
# 2. 读取数据
# =========================
bacteria <- read.delim("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/Bacteria_DNA_alpha_diversity_rarefied.txt")
fungi <- read.delim("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/Fungi_alpha_diversity_rarefied.txt")
metadata <- read.csv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/metadata_fungi.csv")

# =========================
# 3. 合并
# =========================
bacteria <- left_join(bacteria, metadata, by = "Sample")
fungi <- left_join(fungi, metadata, by = "Sample")

# =========================
# 4. 颜色
# =========================
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

season_order <- c("Autumn", "Winter", "Spring", "Summer")
group_order <- c("Underground", "Aboveground", "Elevator", "Outdoor")

# =========================
# 5. 输出目录
# =========================
outdir <- "/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/Diversity/"
dir.create(outdir, showWarnings = FALSE)

# =========================
# 6. 统计函数
# =========================
run_stats <- function(df, xvar, yvar, prefix){
  
  df <- df %>% drop_na(.data[[xvar]], .data[[yvar]])
  
  # Kruskal
  kw <- df %>%
    kruskal_test(as.formula(paste(yvar, "~", xvar)))
  
  write.csv(kw,
            paste0(outdir, prefix, "_", yvar, "_", xvar, "_Kruskal.csv"),
            row.names = FALSE)
  
  # Wilcox with BH adjustment
  pw <- df %>%
    pairwise_wilcox_test(as.formula(paste(yvar, "~", xvar)),
                         p.adjust.method = "BH")
  
  write.csv(pw,
            paste0(outdir, prefix, "_", yvar, "_", xvar, "_Wilcox.csv"),
            row.names = FALSE)
}

# =========================
# 7. 生成 comparisons（修复Season顺序）
# =========================
get_comparisons <- function(df, xvar){
  
  if(xvar == "Season"){
    groups <- season_order
    groups <- groups[groups %in% unique(df[[xvar]])]
  } else if(xvar == "Group"){
    groups <- group_order
    groups <- groups[groups %in% unique(df[[xvar]])]
  } else {
    groups <- unique(df[[xvar]])
  }
  
  combn(groups, 2, simplify = FALSE)
}

# =========================
# 8. 作图函数（顶刊风格，修复统计学显示）
# =========================
plot_alpha <- function(df, xvar, yvar, colors, title){
  
  df <- df %>% drop_na(.data[[xvar]], .data[[yvar]])
  
  # 强制因子顺序
  if(xvar == "Season"){
    df[[xvar]] <- factor(df[[xvar]], levels = season_order)
  } else if(xvar == "Group"){
    df[[xvar]] <- factor(df[[xvar]], levels = group_order)
  }
  
  comparisons <- get_comparisons(df, xvar)
  
  ymax <- max(df[[yvar]], na.rm = TRUE)
  y_range <- diff(range(df[[yvar]], na.rm = TRUE))
  
  # 使用自定义标注函数，只显示星号
  p <- ggplot(df, aes_string(x = xvar, y = yvar, fill = xvar)) +
    
    geom_boxplot(width = 0.65,
                 outlier.shape = NA,
                 color = "black",
                 size = 0.6) +
    
    geom_jitter(aes_string(color = xvar),
                width = 0.12,
                size = 1.8,
                alpha = 0.7) +
    
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    # 总体Kruskal检验
    annotate("text",
             x = 1.5,
             y = ymax + y_range * 0.08,
             label = paste("Kruskal p =", format(kruskal.test(as.formula(paste(yvar, "~", xvar)), data = df)$p.value, digits = 3)),
             size = 4,
             hjust = 0) +
    
    # 两两Wilcox检验，只显示星号
    stat_compare_means(comparisons = comparisons,
                       method = "wilcox.test",
                       label = "p.signif",
                       hide.ns = TRUE,
                       step.increase = 0.08,
                       size = 6,
                       bracket.size = 0.5) +
    
    labs(title = title,
         x = NULL,
         y = yvar) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_text(size = 11),
      legend.position = "none",
      axis.line = element_line(size = 0.8)
    )
  
  return(p)
}

# =========================
# 9. 批量画图（细菌+真菌）- A4尺寸
# =========================
alpha_vars <- c("richness", "Shannon", "evenness")

datasets <- list(
  Bacteria = bacteria,
  Fungi = fungi
)

for (name in names(datasets)) {
  
  df <- datasets[[name]]
  
  for (var in alpha_vars) {
    
    # ===== Group =====
    run_stats(df, "Group", var, name)
    
    p1 <- plot_alpha(df, "Group", var, group_colors,
                     paste(name, var, "by Group"))
    
    ggsave(paste0(outdir, name, "_", var, "_Group.pdf"),
           p1, width = 7, height = 6, dpi = 300)  # A4尺寸
    
    # ===== Season =====
    run_stats(df, "Season", var, name)
    
    p2 <- plot_alpha(df, "Season", var, season_colors,
                     paste(name, var, "by Season"))
    
    ggsave(paste0(outdir, name, "_", var, "_Season.pdf"),
           p2, width = 7, height = 6, dpi = 300)  # A4尺寸
  }
}


# =========================
# 10. 相关性分析（标注R和p值）
# =========================
merged <- inner_join(
  bacteria %>% select(Sample, richness, Shannon, evenness, Group),
  fungi %>% select(Sample, richness, Shannon, evenness),
  by = "Sample",
  suffix = c("_bac", "_fun")
)

# 计算相关性并保存
calculate_and_save_correlations <- function(df, var) {
  
  results <- data.frame()
  
  # 总体相关性
  overall_cor <- cor.test(df[[paste0(var, "_bac")]], 
                          df[[paste0(var, "_fun")]], 
                          method = "spearman")
  
  overall_result <- data.frame(
    Group = "Overall",
    Variable = var,
    Correlation = overall_cor$estimate,
    P_value = overall_cor$p.value,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(results, overall_result)
  
  # 每个Group的相关性
  for(g in unique(df$Group)) {
    sub_df <- df %>% filter(Group == g)
    if(nrow(sub_df) >= 3) {  # 至少需要3个点才能计算相关性
      cor_test <- cor.test(sub_df[[paste0(var, "_bac")]], 
                           sub_df[[paste0(var, "_fun")]], 
                           method = "spearman")
      
      group_result <- data.frame(
        Group = g,
        Variable = var,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value,
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, group_result)
    } else {
      group_result <- data.frame(
        Group = g,
        Variable = var,
        Correlation = NA,
        P_value = NA,
        stringsAsFactors = FALSE
      )
      results <- rbind(results, group_result)
    }
  }
  
  write.csv(results, 
            paste0(outdir, "Correlation_", var, "_statistics.csv"),
            row.names = FALSE)
  
  return(results)
}

# 格式化p值的函数
format_p_value <- function(p) {
  if(p < 0.001) return("p < 0.001")
  if(p < 0.01) return(paste("p =", format(p, digits = 3)))
  if(p < 0.05) return(paste("p =", format(p, digits = 3)))
  return(paste("p =", format(p, digits = 3)))
}

# 修改后的作图函数（标注R和p值）
plot_corr <- function(df, var){
  
  # 计算总体和分组的相关系数
  cor_stats <- calculate_and_save_correlations(df, var)
  
  # 总体相关性标注（显示ρ和p值）
  overall_cor <- cor_stats %>% filter(Group == "Overall")
  
  # 添加星号标记
  sig_symbol <- ifelse(overall_cor$P_value < 0.001, "***",
                       ifelse(overall_cor$P_value < 0.01, "**",
                              ifelse(overall_cor$P_value < 0.05, "*", "")))
  
  cor_label <- paste0("Overall: ρ = ", round(overall_cor$Correlation, 3),
                      sig_symbol, "\n",
                      format_p_value(overall_cor$P_value))
  
  p <- ggplot(df,
              aes_string(x = paste0(var, "_bac"),
                         y = paste0(var, "_fun"),
                         color = "Group")) +
    
    geom_point(size = 2.5, alpha = 0.8) +
    
    geom_smooth(method = "lm",
                se = TRUE,
                color = "black",
                linewidth = 0.8) +
    
    annotate("text",
             x = min(df[[paste0(var, "_bac")]], na.rm = TRUE),
             y = max(df[[paste0(var, "_fun")]], na.rm = TRUE),
             label = cor_label,
             size = 4,
             hjust = 0,
             vjust = 1,
             fontface = "italic") +
    
    scale_color_manual(values = group_colors) +
    
    labs(
      x = paste("Bacteria", var),
      y = paste("Fungi", var),
      title = paste("Correlation:", var)
    ) +
    
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.title = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 9)
    )
  
  return(p)
}

# 生成每个Group单独的相关性图（标注R和p值）
plot_corr_by_group <- function(df, var, group_name){
  
  sub_df <- df %>% filter(Group == group_name)
  
  if(nrow(sub_df) >= 3) {
    cor_test <- cor.test(sub_df[[paste0(var, "_bac")]], 
                         sub_df[[paste0(var, "_fun")]], 
                         method = "spearman")
    
    # 添加星号标记
    sig_symbol <- ifelse(cor_test$p.value < 0.001, "***",
                         ifelse(cor_test$p.value < 0.01, "**",
                                ifelse(cor_test$p.value < 0.05, "*", "")))
    
    cor_label <- paste0("ρ = ", round(cor_test$estimate, 3),
                        sig_symbol, "\n",
                        format_p_value(cor_test$p.value))
    
    p <- ggplot(sub_df,
                aes_string(x = paste0(var, "_bac"),
                           y = paste0(var, "_fun"))) +
      geom_point(size = 3, alpha = 0.8, color = group_colors[group_name]) +
      geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
      annotate("text",
               x = min(sub_df[[paste0(var, "_bac")]], na.rm = TRUE),
               y = max(sub_df[[paste0(var, "_fun")]], na.rm = TRUE),
               label = cor_label,
               size = 4,
               hjust = 0,
               vjust = 1,
               fontface = "italic") +
      labs(
        x = paste("Bacteria", var),
        y = paste("Fungi", var),
        title = paste("Correlation:", var, "-", group_name)
      ) +
      theme_classic(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
      )
    
    return(p)
  } else {
    return(NULL)
  }
}

# 生成总体相关性图（A4尺寸）
for (var in alpha_vars) {
  
  p <- plot_corr(merged, var)
  
  ggsave(paste0(outdir, "Correlation_", var, ".pdf"),
         p, width = 7, height = 6)  # A4尺寸
}

# 生成每个Group单独的相关性图（A4尺寸）
for (var in alpha_vars) {
  for(g in unique(merged$Group)) {
    p <- plot_corr_by_group(merged, var, g)
    if(!is.null(p)) {
      ggsave(paste0(outdir, "Correlation_", var, "_", g, ".pdf"),
             p, width = 8.27, height = 11.69, units = "in")
    }
  }
}
# =========================================
# PCoA 绘图函数
# =========================================


season_shapes <- c(
  "Spring" = 4,  
  "Summer" = 10,  
  "Autumn" = 23,  # 菱形
  "Winter" = 18   # 三角
)

group_shapes <- c(
  "Underground" = 16,  # 实心方
  "Aboveground" = 12,  # 实心圆
  "Elevator"    = 6,  # 实心三角
  "Outdoor"     = 8   # 实心菱形
)


plot_pcoa <- function(pcoa_file, var1="PCo1", var2="PCo2", title_prefix="PCoA") {
  
  PCoA <- read.delim(pcoa_file, header=TRUE)
  
  # 提取解释度（如果列名包含百分比可以进一步优化）
  var1_perc <- round(var(PCoA[[var1]]) / 
                       (var(PCoA[[var1]]) + var(PCoA[[var2]])) * 100, 1)
  var2_perc <- round(var(PCoA[[var2]]) / 
                       (var(PCoA[[var1]]) + var(PCoA[[var2]])) * 100, 1)
  
  #========================
  # 1️⃣ Season主导（颜色=Season，形状=Group）
  #========================
  p1 <- ggplot(PCoA, aes(x = .data[[var1]], y = .data[[var2]])) +
    
    geom_point(aes(color = Season, shape = Group),
               size = 3, alpha = 0.9) +
    
    stat_ellipse(aes(color = Season),
                 type = "t", linewidth = 1) +
    
    scale_color_manual(values = season_colors) +
    scale_shape_manual(values = group_shapes) +
    
    theme_classic(base_size = 14) +
    
    labs(
      x = paste0(var1," (", var1_perc, "%)"),
      y = paste0(var2," (", var2_perc, "%)"),
      title = paste0(title_prefix," – Season (color) + Group (shape)")
    )
  
  ggsave(paste0(title_prefix,"_SeasonColor_GroupShape.pdf"),
         p1, width = 7, height = 6, dpi = 300)
  
  
  #========================
  # 2️⃣ Group主导（颜色=Group，形状=Season）
  #========================
  p2 <- ggplot(PCoA, aes(x = .data[[var1]], y = .data[[var2]])) +
    
    geom_point(aes(color = Group, shape = Season),
               size = 3, alpha = 0.9) +
    
    stat_ellipse(aes(color = Group),
                 type = "t", linewidth = 1) +
    
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = season_shapes) +
    
    theme_classic(base_size = 14) +
    
    labs(
      x = paste0(var1," (", var1_perc, "%)"),
      y = paste0(var2," (", var2_perc, "%)"),
      title = paste0(title_prefix," – Group (color) + Season (shape)")
    )
  
  ggsave(paste0(title_prefix,"_GroupColor_SeasonShape.pdf"),
         p2, width = 7, height = 6, dpi = 300)
  
  return(PCoA)
}
# 示例：绘制细菌 PCoA
PCoA_unwei <- plot_pcoa("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/PCoA_unwei_unifrac_PCoA.tsv", title_prefix="PCoA_Unweighted")
PCoA_wei <- plot_pcoa("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/PCoA_wei_unifrac_PCoA.tsv", title_prefix="PCoA_Weighted")

############################3Pathogen——AlpahDiversity
library(tidyverse)
library(rstatix)
library(ggpubr)

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

fungi_raw <- read_tsv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/Fungi_DNA_Pathogen_AlphaDiversity.tsv")

fungi_shannon <- fungi_raw %>%
  filter(alpha == "Shannon") %>%
  rename(Shannon = value) %>%
  mutate(
    Group = factor(Group, levels = c("Underground","Aboveground","Elevator","Outdoor")),
    Season = factor(Season, levels = c("Spring","Summer","Autumn","Winter"))
  )
# 读取细菌 Shannon
bac_shannon <- read_tsv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/Pathogen_shannon.txt")
sam_data <- read_csv(
  '/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/metadata_fungi.csv'
)
# 读取 metadata
meta_bac<- sam_data

# ---- 关键：构建 Group 和 Season ----
meta_bac <- meta_bac %>%
  select(Sample, Group, Season)

# 合并
bac_df <- bac_shannon %>%
  left_join(meta_bac, by = "Sample") %>%
  filter(!is.na(Shannon)) %>%
  mutate(
    Group = factor(Group, levels = c("Underground","Aboveground","Elevator","Outdoor")),
    Season = factor(Season, levels = c("Spring","Summer","Autumn","Winter"))
  )

run_stats <- function(df, xvar){
  list(
    kw = df %>% kruskal_test(as.formula(paste("Shannon ~", xvar))),
    pw = df %>% wilcox_test(as.formula(paste("Shannon ~", xvar)), p.adjust.method = "BH")
  )
}

# 真菌
fungi_group_stats  <- run_stats(fungi_shannon, "Group")
fungi_season_stats <- run_stats(fungi_shannon, "Season")

# 细菌
bac_group_stats  <- run_stats(bac_df, "Group")
bac_season_stats <- run_stats(bac_df, "Season")

plot_alpha <- function(df, xvar, colors, title){
  
  pairs <- if(xvar == "Group"){
    list(
      c("Underground","Aboveground"),
      c("Underground","Elevator"),
      c("Underground","Outdoor"),
      c("Aboveground","Elevator"),
      c("Aboveground","Outdoor"),
      c("Elevator","Outdoor")
    )
  } else {
    list(
      c("Spring","Summer"),
      c("Summer","Autumn"),
      c("Autumn","Winter"),
      c("Winter","Spring")
    )
  }
  
  ggplot(df, aes_string(x = xvar, y = "Shannon")) +
    geom_boxplot(aes_string(fill = xvar), outlier.shape = NA, width = 0.6) +
    geom_jitter(aes_string(color = xvar), width = 0.15, size = 1.5, alpha = 0.6) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    
    stat_compare_means(
      method = "kruskal.test",
      label = "p.format",
      label.y.npc = 0.95,
      size = 4
    ) +
    
    stat_compare_means(
      comparisons = pairs,
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
      y = "Shannon diversity",
      title = title
    )
}

# ---- 真菌 ----
p_fungi_group <- plot_alpha(
  fungi_shannon,
  "Group",
  group_colors,
  "Fungal Shannon diversity across environments"
)

p_fungi_season <- plot_alpha(
  fungi_shannon,
  "Season",
  season_colors,
  "Seasonal variation of fungal Shannon diversity"
)

# ---- 细菌 ----
p_bac_group <- plot_alpha(
  bac_df,
  "Group",
  group_colors,
  "Bacterial Shannon diversity across environments"
)

p_bac_season <- plot_alpha(
  bac_df,
  "Season",
  season_colors,
  "Seasonal variation of bacterial Shannon diversity"
)

ggsave("Pathongen_Fungi_Shannon_Group.pdf", p_fungi_group, width = 7, height = 6)
ggsave("Pathongen_Fungi_Shannon_Season.pdf", p_fungi_season, width = 7, height = 6)

ggsave("Pathongen_Bacteria_Shannon_Group.pdf", p_bac_group, width = 7, height = 6)
ggsave("Pathongen_Bacteria_Shannon_Season.pdf", p_bac_season, width = 7, height = 6)

##############################Pathogen_PCOA
library(tidyverse)
library(vegan)
library(ggplot2)
library(ggthemes)
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

season_shapes <- c(
  "Spring" = 4,
  "Summer" = 10,
  "Autumn" = 23,
  "Winter" = 18
)

group_shapes <- c(
  "Underground" = 16,
  "Aboveground" = 12,
  "Elevator"    = 6,
  "Outdoor"     = 8
)



run_pcoa <- function(df){
  
  # ✅ 1. 强制 count 为 numeric
  df <- df %>%
    mutate(count = suppressWarnings(as.numeric(count)))
  
  # NA → 0（非常关键）
  df$count[is.na(df$count)] <- 0
  
  # ✅ 2. pivot
  df_cast <- df %>%
    pivot_wider(
      names_from = ASV,
      values_from = count,
      values_fill = 0
    )
  
  # 转 data.frame
  df_cast <- as.data.frame(df_cast)
  
  # 行名
  rownames(df_cast) <- df_cast$Sample
  
  # 去掉 Sample
  df_cast <- df_cast[ , -1]
  
  # ✅ 3. 强制 numeric（双保险）
  df_cast[] <- lapply(df_cast, function(x){
    x <- suppressWarnings(as.numeric(x))
    x[is.na(x)] <- 0
    return(x)
  })
  
  # ✅ 4. 去掉全0列（防止 colSums 报错）
  keep_cols <- colSums(df_cast) > 0
  
  # 如果全是0（极端情况保护）
  if(sum(keep_cols) == 0){
    stop("All ASVs are zero after filtering. Check your data.")
  }
  
  df_cast <- df_cast[ , keep_cols, drop = FALSE]
  
  # ✅ 5. Bray-Curtis
  dist <- vegdist(df_cast, method = "bray")
  
  # ✅ 6. PCoA
  pcoa <- cmdscale(dist, k = 2, eig = TRUE)
  
  pcoa_df <- data.frame(
    Sample = rownames(pcoa$points),
    PCoA1 = pcoa$points[,1],
    PCoA2 = pcoa$points[,2]
  )
  
  eig <- pcoa$eig
  var1 <- round(eig[1]/sum(eig)*100,2)
  var2 <- round(eig[2]/sum(eig)*100,2)
  
  list(pcoa_df = pcoa_df, dist = dist, var1 = var1, var2 = var2)
}


run_beta_stats <- function(dist, meta, factor_name, prefix){
  
  permanova <- adonis2(dist ~ meta[[factor_name]], permutations = 999)
  
  write.csv(as.data.frame(permanova),
            paste0(prefix, "_Pathogen_PERMANOVA.csv"))
  
  disp <- betadisper(dist, meta[[factor_name]])
  disp_test <- permutest(disp, permutations = 999)
  
  write.csv(as.data.frame(disp_test$tab),
            paste0(prefix, "_Pathogen_PERMDISP.csv"))
}

plot_pcoa <- function(df, var1, var2, color_var, shape_var,
                      color_values, shape_values, title){
  
  ggplot(df, aes(PCoA1, PCoA2)) +
    geom_point(aes_string(color = color_var, shape = shape_var),
               size = 3, alpha = 0.8) +
    stat_ellipse(aes_string(color = color_var),
                 type = "t", linewidth = 1) +
    scale_color_manual(values = color_values) +
    scale_shape_manual(values = shape_values) +
    labs(
      x = paste0("PCoA1 (", var1, "%)"),
      y = paste0("PCoA2 (", var2, "%)"),
      title = title
    ) +
    theme_classic(base_size = 14)
}

plot_pcoa <- function(df, var1, var2, color_var, shape_var,
                      color_values, shape_values, title){
  
  ggplot(df, aes(PCoA1, PCoA2)) +
    geom_point(aes_string(color = color_var, shape = shape_var),
               size = 3, alpha = 0.8) +
    stat_ellipse(aes_string(color = color_var),
                 type = "t", linewidth = 1) +
    scale_color_manual(values = color_values) +
    scale_shape_manual(values = shape_values) +
    labs(
      x = paste0("PCoA1 (", var1, "%)"),
      y = paste0("PCoA2 (", var2, "%)"),
      title = title
    ) +
    theme_classic(base_size = 14)
}

fungi <- read_tsv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/Fungi_DNA_Pathogen_Relative_Abundance_Metadata.tsv")

fungi2 <- fungi %>%
  select(Sample, ASV, count, Group, Season) %>%
  distinct()

fungi_res <- run_pcoa(fungi2)

fungi_pcoa <- fungi_res$pcoa_df %>%
  left_join(fungi2 %>% distinct(Sample, Group, Season), by = "Sample") %>%
  mutate(
    Group = factor(Group, levels = c("Underground","Aboveground","Elevator","Outdoor")),
    Season = factor(Season, levels = c("Spring","Summer","Autumn","Winter"))
  )
p_fungi_season <- plot_pcoa(
  fungi_pcoa,
  fungi_res$var1,
  fungi_res$var2,
  "Season",
  "Group",
  season_colors,
  group_shapes,
  "Fungal pathogen PCoA (colored by Season)"
)

ggsave("Pathogen_Fungi_PCoA_Season.pdf", p_fungi_season, width = 7, height = 6)

run_beta_stats(fungi_res$dist, fungi_pcoa, "Season", "Fungi_Season")

p_fungi_group <- plot_pcoa(
  fungi_pcoa,
  fungi_res$var1,
  fungi_res$var2,
  "Group",
  "Season",
  group_colors,
  season_shapes,
  "Fungal pathogen PCoA (colored by Group)"
)

ggsave("Pathogen_Fungi_PCoA_Group.pdf", p_fungi_group, width = 7, height = 6)

run_beta_stats(fungi_res$dist, fungi_pcoa, "Group", "Fungi_Group")

bac <- read_tsv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/Pathogen_ASV_filtered.tidy.tsv")

meta <- sam_data %>%
  select(Sample, Group, Season)

bac2 <- bac %>%
  left_join(meta, by = "Sample") %>%
  select(Sample, ASV, count = Count, Group, Season)


bac_res <- run_pcoa(bac2)

bac_pcoa <- bac_res$pcoa_df %>%
  left_join(meta, by = "Sample") %>%
  mutate(
    Group = factor(Group, levels = c("Underground","Aboveground","Elevator","Outdoor")),
    Season = factor(Season, levels = c("Spring","Summer","Autumn","Winter"))
  )

p_bac_season <- plot_pcoa(
  bac_pcoa,
  bac_res$var1,
  bac_res$var2,
  "Season",
  "Group",
  season_colors,
  group_shapes,
  "Bacterial pathogen PCoA (colored by Season)"
)

ggsave("Pathogen_Bacteria_PCoA_Season.pdf", p_bac_season, width = 7, height = 6)

run_beta_stats(bac_res$dist, bac_pcoa, "Season", "Bacteria_Season")

p_bac_group <- plot_pcoa(
  bac_pcoa,
  bac_res$var1,
  bac_res$var2,
  "Group",
  "Season",
  group_colors,
  season_shapes,
  "Bacterial pathogen PCoA (colored by Group)"
)

ggsave("Pathogen_Bacteria_PCoA_Group.pdf", p_bac_group, width = 7, height = 6)

run_beta_stats(bac_res$dist, bac_pcoa, "Group", "Bacteria_Group")
 ####################################3
#########################################Correlation with Temp
library(tidyverse)

#---------------------------
# 1. 读取数据
#---------------------------
bac <- read.table("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/Pathogen_richness.txt",
                  header = TRUE, sep = "\t")

fungi <- read.table("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/Fungi_DNA_Pathogen_AlphaDiversity.tsv",
                    header = TRUE, sep = "\t")

meta <- read.csv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/metadata_fungi.csv")

#---------------------------
# 2. 数据整理
#---------------------------

# 细菌
bac_df <- bac %>%
  rename(Richness = richness) %>%
  mutate(Type = "Bacteria") %>%
  left_join(meta, by = "Sample") %>%
  select(-Type.y) %>%          # 删除metadata里的Type
  rename(Type = Type.x)        # 保留我们自己的Type

# 真菌
fungi_df <- fungi %>%
  filter(alpha == "Richness") %>%
  select(Sample, value) %>%
  rename(Richness = value) %>%
  mutate(Type = "Fungi") %>%
  left_join(meta, by = "Sample") %>%
  select(-Type.y) %>%
  rename(Type = Type.x)

# 合并
all_df <- bind_rows(bac_df, fungi_df)

#---------------------------
# 3. 计算 Spearman 相关性
#---------------------------

# Temperature
cor_temp <- all_df %>%
  group_by(Type) %>%
  summarise(
    R = cor(Temperature, Richness, method = "spearman", use = "complete.obs"),
    p = cor.test(Temperature, Richness, method = "spearman")$p.value
  )

# Humidity
cor_hum <- all_df %>%
  group_by(Type) %>%
  summarise(
    R = cor(Humidness, Richness, method = "spearman", use = "complete.obs"),
    p = cor.test(Humidness, Richness, method = "spearman")$p.value
  )

#---------------------------
# 4. 画图函数
#---------------------------

plot_func <- function(data, xvar, cor_df, xlab) {
  
  ggplot(data, aes_string(x = xvar, y = "Richness",
                          color = "Type", shape = "Type")) +
    
    geom_point(size = 3, alpha = 0.8) +
    
    geom_smooth(method = "loess", se = TRUE, color = "black") +
    
    scale_color_manual(values = c("Bacteria" = "#2c7fb8",
                                  "Fungi" = "#d95f02")) +
    
    scale_shape_manual(values = c("Bacteria" = 15,
                                  "Fungi" = 17)) +
    
    theme_bw() +
    
    labs(x = xlab, y = "Richness") +
    
    theme(
      legend.title = element_blank(),
      text = element_text(size = 14)
    ) +
    
    annotate("text",
             x = min(data[[xvar]], na.rm = TRUE),
             y = max(data$Richness, na.rm = TRUE),
             hjust = 0,
             label = paste0("Bacteria: R = ",
                            round(cor_df$R[cor_df$Type == "Bacteria"], 2),
                            ", p = ",
                            signif(cor_df$p[cor_df$Type == "Bacteria"], 2))) +
    
    annotate("text",
             x = min(data[[xvar]], na.rm = TRUE),
             y = max(data$Richness, na.rm = TRUE) * 0.9,
             hjust = 0,
             label = paste0("Fungi: R = ",
                            round(cor_df$R[cor_df$Type == "Fungi"], 2),
                            ", p = ",
                            signif(cor_df$p[cor_df$Type == "Fungi"], 2)))
}

#---------------------------
# 5. 作图
#---------------------------

p_temp <- plot_func(all_df, "Temperature", cor_temp, "Temperature (°C)")
p_hum  <- plot_func(all_df, "Humidness", cor_hum, "Humidity")

p_temp
p_hum

#---------------------------
# 6. 保存（可选）
#---------------------------

ggsave("Temperature_vs_Richness.png", p_temp, width = 6, height = 5, dpi = 300)
ggsave("Humidity_vs_Richness.png", p_hum, width = 6, height = 5, dpi = 300)

#####################Correlation Temp V2
library(tidyverse)

#---------------------------
# 1. 读取数据
#---------------------------
bac <- read.table("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/Pathogen_richness.txt",
                  header = TRUE, sep = "\t")

fungi <- read.table("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/Fungi_DNA_Pathogen_AlphaDiversity.tsv",
                    header = TRUE, sep = "\t")

meta <- read.csv("/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/metadata_fungi.csv")

#---------------------------
# 2. 数据整理
#---------------------------

# 去掉 metadata 里的 Type（避免冲突）
meta2 <- meta %>% select(-Type)

# 细菌
bac_df <- bac %>%
  rename(Richness = richness) %>%
  mutate(Type = "Bacteria") %>%
  left_join(meta2, by = "Sample")

# 真菌
fungi_df <- fungi %>%
  filter(alpha == "Richness") %>%
  select(Sample, value) %>%
  rename(Richness = value) %>%
  mutate(Type = "Fungi") %>%
  left_join(meta2, by = "Sample")

#---------------------------
# 3. 转换成长格式（方便分面）
#---------------------------
bac_long <- bac_df %>%
  pivot_longer(cols = c(Temperature, Humidness),
               names_to = "Variable",
               values_to = "Env")

fungi_long <- fungi_df %>%
  pivot_longer(cols = c(Temperature, Humidness),
               names_to = "Variable",
               values_to = "Env")

#---------------------------
# 4. 计算 Spearman（分别）
#---------------------------

calc_cor <- function(df) {
  df %>%
    group_by(Variable) %>%
    summarise(
      R = cor(Env, Richness, method = "spearman", use = "complete.obs"),
      p = cor.test(Env, Richness, method = "spearman")$p.value
    )
}

cor_bac <- calc_cor(bac_long)
cor_fun <- calc_cor(fungi_long)

#---------------------------
# 5. 画图函数（单一Type）
#---------------------------

plot_single <- function(df, cor_df, point_shape, point_color, title_name) {
  
  ggplot(df, aes(x = Env, y = Richness)) +
    
    geom_point(shape = point_shape, color = point_color, size = 2, alpha = 0.8) +
    
    geom_smooth(method = "loess", color = "black", fill = "grey80") +
    
    facet_wrap(~Variable, scales = "free_x",
               labeller = labeller(
                 Variable = c(Temperature = "Temperature (°C)",
                              Humidness = "Humidity")
               )) +
    
    theme_bw() +
    
    labs(y = "Richness", x = NULL, title = title_name) +
    
    theme(
      text = element_text(size = 14),
      strip.text = element_text(face = "bold")
    ) +
    
    # 添加相关性
    geom_text(
      data = cor_df,
      aes(x = -Inf, y = Inf,
          label = paste0("R = ", round(R, 2),
                         "\np = ", signif(p, 2))),
      hjust = -0.1, vjust = 1.2,
      inherit.aes = FALSE
    )
}

#---------------------------
# 6. 作图
#---------------------------

p_bac <- plot_single(bac_long, cor_bac,
                     point_shape = 16,
                     point_color = "#2c7fb8",
                     title_name = "Bacteria")

p_fun <- plot_single(fungi_long, cor_fun,
                     point_shape = 17,
                     point_color = "#d95f02",
                     title_name = "Fungi")
#p_fun <- plot_single(fungi_long, cor_fun,
#                     point_shape = 17,
#                     point_color = "#d95f02",
#                     title_name = "Fungi") +
#  
#  scale_y_continuous(limits = c(0, NA))
# 显示
p_bac
p_fun

#---------------------------
# 7. 保存
#---------------------------
ggsave("/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/Bacteria_env_correlation.pdf", p_bac, width = 8, height = 4, dpi = 300)
ggsave("/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/Fungi_env_correlation.pdf", p_fun, width = 8, height = 4, dpi = 300)
#####################
####################MetaPlot

# 加载包
library(dplyr)
library(ggplot2)
library(lubridate)

# 自定义分组颜色（顺序与期望的四个组一致）
group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

# 数据清洗与过滤
df_filtered <- metadata %>%
  filter(
    Notes != "Drop",                      # 移除 Notes 为 Drop 的行
    !grepl("^R", Sample),                # 保留 Sample 首字母不是 R 的行
    Group %in% c("Underground", "Aboveground", "Elevator", "Outdoor")  # 仅保留这四个组
  ) %>%
  mutate(
    Date = as.Date(paste0("20", Date), format = "%Y%m%d"),
    # 将 Group 转换为因子，按指定顺序排列（确保颜色映射正确）
    Group = factor(Group, levels = names(group_colors))
  ) %>%
  arrange(Date)

# 检查过滤后的组别构成（可选）
# table(df_filtered$Group)

# ===================== 温度趋势图 =====================
p_temp <- ggplot(df_filtered, aes(x = Date, y = Temperature, color = Group, fill = Group)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, size = 1.2) +
  geom_point(alpha = 0.5, size = 1.8, stroke = 0.2) +
  facet_wrap(~ Group, nrow = 1, scales = "free_x") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(
    title = "Temperature trends over time",
    x = "Date",
    y = "Temperature (°C)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# ===================== 湿度趋势图 =====================
p_humid <- ggplot(df_filtered, aes(x = Date, y = Humidness, color = Group, fill = Group)) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2, size = 1.2) +
  geom_point(alpha = 0.5, size = 1.8, stroke = 0.2) +
  facet_wrap(~ Group, nrow = 1, scales = "free_x") +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(
    title = "Humidity trends over time",
    x = "Date",
    y = "Relative humidity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# 显示图形
print(p_temp)
print(p_humid)

# 保存为高分辨率图片（PDF 矢量格式，适合投稿）
ggsave("Temperature_trend.pdf", p_temp, width = 16, height = 5, dpi = 300, device = "pdf")
ggsave("Humidity_trend.pdf", p_humid, width = 16, height = 5, dpi = 300, device = "pdf")

################V2
# 加载包
library(dplyr)
library(ggplot2)
library(lubridate)

# 自定义分组颜色
group_colors <- c(
  "Underground" = "#4c628c",
  "Aboveground" = "#f0c396",
  "Elevator"    = "#ac9ec2",
  "Outdoor"     = "#dbb1bb"
)

# 数据清洗与过滤
df_filtered <- metadata %>%
  filter(
    Notes != "Drop",
    !grepl("^R", Sample),
    Group %in% c("Underground", "Aboveground", "Elevator", "Outdoor")
  ) %>%
  mutate(
    Date = as.Date(paste0("20", Date), format = "%Y%m%d"),
    Group = factor(Group, levels = names(group_colors))
  ) %>%
  arrange(Date)

# ===================== 温度趋势图（折线，无平滑） =====================
p_temp_raw <- ggplot(df_filtered, aes(x = Date, y = Temperature, color = Group, group = 1)) +
  geom_line(size = 1.0) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ Group, nrow = 1, scales = "free_x") +
  scale_color_manual(values = group_colors) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(
    title = "Temperature trends over time (raw data)",
    x = "Date",
    y = "Temperature (°C)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# ===================== 湿度趋势图（折线，无平滑） =====================
p_humid_raw <- ggplot(df_filtered, aes(x = Date, y = Humidness, color = Group, group = 1)) +
  geom_line(size = 1.0) +
  geom_point(size = 2, alpha = 0.7) +
  facet_wrap(~ Group, nrow = 1, scales = "free_x") +
  scale_color_manual(values = group_colors) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%b %d") +
  labs(
    title = "Humidity trends over time (raw data)",
    x = "Date",
    y = "Relative humidity"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.title = element_text(size = 12),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.spacing = unit(1, "lines")
  )

# 显示图形
print(p_temp_raw)
print(p_humid_raw)

# 保存图片
ggsave("Temperature_trend_raw.pdf", p_temp_raw, width = 16, height = 5, dpi = 300, device = "pdf")
ggsave("Humidity_trend_raw.pdf", p_humid_raw, width = 16, height = 5, dpi = 300, device = "pdf")



############Check infornation
########Meta Info Merge
library(dplyr)

# ============================================================
# 1. 读取两个 meta 表
# ============================================================
meta_bac <- read.csv(
  "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/metadata_v2.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)

meta_fungi <- read.csv(
  "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/feature table/metadata_fungi.csv",
  stringsAsFactors = FALSE, check.names = FALSE
)

# 清理 Sample 列空格
meta_bac$Sample   <- trimws(meta_bac$Sample)
meta_fungi$Sample <- trimws(meta_fungi$Sample)

cat("细菌 meta 行数:", nrow(meta_bac),   "| 列数:", ncol(meta_bac), "\n")
cat("真菌 meta 行数:", nrow(meta_fungi), "| 列数:", ncol(meta_fungi), "\n")

# ============================================================
# 2. 统一列名映射（真菌列名 → 细菌列名）
# ============================================================
# 真菌列名        细菌列名
col_map <- c(
  "Indoor_Outdoor" = "Indoor.Outdoor",
  "Air_out"        = "Air.out",
  "DNA_RNA"        = "DNARNA"
)

# 将真菌 meta 中需要重命名的列改为细菌 meta 的列名
for (old_name in names(col_map)) {
  new_name <- col_map[old_name]
  if (old_name %in% colnames(meta_fungi)) {
    colnames(meta_fungi)[colnames(meta_fungi) == old_name] <- new_name
  }
}

cat("\n真菌 meta 列名（重命名后）:\n")
print(colnames(meta_fungi))
cat("\n细菌 meta 列名:\n")
print(colnames(meta_bac))

# ============================================================
# 3. 找出真菌 meta 独有的样本（不在细菌 meta 中的）
# ============================================================
samples_in_bac   <- meta_bac$Sample
samples_in_fungi <- meta_fungi$Sample

only_in_fungi <- setdiff(samples_in_fungi, samples_in_bac)
both           <- intersect(samples_in_fungi, samples_in_bac)

cat("\n两表共有样本数:", length(both), "\n")
cat("仅在真菌 meta 中的样本数 (需追加):", length(only_in_fungi), "\n")
if (length(only_in_fungi) > 0) print(only_in_fungi)

# ============================================================
# 4. 从真菌 meta 中提取需要追加的行，对齐列到细菌 meta
# ============================================================
fungi_to_append <- meta_fungi %>%
  filter(Sample %in% only_in_fungi)

# 以细菌 meta 的列为基准：
#   - 真菌有的列 → 保留值
#   - 真菌没有的列 → 填 NA
all_cols <- colnames(meta_bac)

fungi_aligned <- fungi_to_append %>%
  # 只保留与细菌 meta 列名相同的列
  select(any_of(all_cols)) %>%
  # 补充细菌 meta 中有但真菌没有的列（填 NA）
  bind_cols(
    setNames(
      as.data.frame(
        matrix(NA, nrow = nrow(fungi_to_append),
               ncol = sum(!all_cols %in% colnames(fungi_to_append)))
      ),
      all_cols[!all_cols %in% colnames(fungi_to_append)]
    )
  ) %>%
  # 按细菌 meta 列顺序排列
  select(all_of(all_cols))

# ============================================================
# 5. 合并：细菌 meta（主体）+ 真菌独有行
# ============================================================
meta_combined <- bind_rows(meta_bac, fungi_aligned)

# 按 Sample 排序，方便后续查阅
meta_combined <- meta_combined %>% arrange(Sample)

cat("\n合并后总行数:", nrow(meta_combined), "\n")
cat("合并后总列数:", ncol(meta_combined), "\n")

# ============================================================
# 6. 保存到本地
# ============================================================
output_path <- "/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro//combined_meta.csv"
write.csv(meta_combined, output_path, row.names = FALSE)
cat("\n✓ 合并后的 meta 已保存至:", output_path, "\n")

# ============================================================
# 7. 简单诊断：检查关键列的填充情况
# ============================================================
cat("\n===== 各列 NA 数量统计 =====\n")
na_counts <- colSums(is.na(meta_combined))
print(na_counts[na_counts > 0])   # 只显示有缺失的列


##########
# 加载必要的包
library(dplyr)
library(tibble)

# 定义文件路径
bacteria_dna_file <- "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/DNA_ASV_OTU.txt"
fungi_dna_file <- "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/fungi_DNA_ASV_OTU.txt"
bacteria_rna_file <- "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/RNA_ASV_OTU.txt"

# 辅助函数：读取 OTU 表并自动处理行名/列名
read_otu_table <- function(file_path) {
  # 读取数据，第一列为行名（ASV ID），第一行为列名（样本名）
  otu <- read.delim(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  return(otu)
}

# 读取三个表格
bacteria_dna <- read_otu_table(bacteria_dna_file)
fungi_dna <- read_otu_table(fungi_dna_file)
bacteria_rna <- read_otu_table(bacteria_rna_file)

# 正确统计：样本数 = 列数，ASV 数 = 行数
bact_dna_samples <- ncol(bacteria_dna)
bact_dna_asvs <- nrow(bacteria_dna)

fungi_dna_samples <- ncol(fungi_dna)
fungi_dna_asvs <- nrow(fungi_dna)

bact_rna_samples <- ncol(bacteria_rna)
bact_rna_asvs <- nrow(bacteria_rna)

# 输出基本信息
cat("===== 细菌 DNA OTU 表 =====\n")
cat("样本数 (Samples):", bact_dna_samples, "\n")
cat("ASV 数目 (ASVs):", bact_dna_asvs, "\n\n")

cat("===== 真菌 DNA OTU 表 =====\n")
cat("样本数 (Samples):", fungi_dna_samples, "\n")
cat("ASV 数目 (ASVs):", fungi_dna_asvs, "\n\n")

cat("===== 细菌 RNA OTU 表 =====\n")
cat("样本数 (Samples):", bact_rna_samples, "\n")
cat("ASV 数目 (ASVs):", bact_rna_asvs, "\n\n")

# 提取样本名（列名）并清理可能的隐藏字符
bact_dna_sample_names <- colnames(bacteria_dna) %>% trimws()  # 去除首尾空格
fungi_dna_sample_names <- colnames(fungi_dna) %>% trimws()

# 检查是否存在非打印字符（可选）
# bact_dna_sample_names <- gsub("[^ -~]", "", bact_dna_sample_names)  # 移除控制字符

# 计算配对样本数
paired_samples <- intersect(bact_dna_sample_names, fungi_dna_sample_names)
paired_count <- length(paired_samples)

cat("===== DNA 样本配对情况 =====\n")
cat("细菌 DNA 样本数:", bact_dna_samples, "\n")
cat("真菌 DNA 样本数:", fungi_dna_samples, "\n")
cat("共同样本名数量 (Paired):", paired_count, "\n")

if (paired_count > 0) {
  cat("前 20 个共同样本名示例:\n")
  print(head(paired_samples, 20))
} else {
  # 如果为 0，尝试找出不匹配的原因
  cat("\n警告：未找到共同样本名！可能原因：\n")
  cat("1. 样本名大小写不一致\n")
  cat("2. 存在不可见字符（如换行符、回车符）\n")
  cat("3. 真菌表或细菌表的列名并非样本名（请检查文件格式）\n")
  cat("\n检查细菌 DNA 表前 5 个列名：\n")
  print(head(bact_dna_sample_names, 5))
  cat("\n检查真菌 DNA 表前 5 个列名：\n")
  print(head(fungi_dna_sample_names, 5))
  
  # 尝试模糊匹配：看是否有包含关系（例如 "230922_10A" vs "230922_10A\n"）
  # 直接比较字符串差异
  common_start <- intersect(substr(bact_dna_sample_names, 1, 10), substr(fungi_dna_sample_names, 1, 10))
  if (length(common_start) > 0) {
    cat("\n可能存在后缀或前缀差异。例如，前10字符匹配的有", length(common_start), "个。\n")
  }
}

##############3
library(dplyr)

# ============================================================
# 1. 读取 OTU 表
# ============================================================
read_otu_table <- function(file_path) {
  read.delim(file_path, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
}

bacteria_dna_file <- "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/DNA_ASV_OTU.txt"
fungi_dna_file    <- "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Fungi/Air_fungal_study/fungi/1alpha_beta_diversity/fungi_DNA_ASV_OTU.txt"
bacteria_rna_file <- "/Users/shijiaoqi/Desktop/PHD_Materials/Air_Bacteria/Air_Bac/RNA_ASV_OTU.txt"

bacteria_dna <- read_otu_table(bacteria_dna_file)
fungi_dna    <- read_otu_table(fungi_dna_file)
bacteria_rna <- read_otu_table(bacteria_rna_file)

# ============================================================
# 2. 提取并清理各 OTU 表的样本名（列名）
# ============================================================
bac_dna_samples  <- trimws(colnames(bacteria_dna))
fungi_dna_samples <- trimws(colnames(fungi_dna))
bac_rna_samples  <- trimws(colnames(bacteria_rna))

# ============================================================
# 3. 计算 Paired / Unpaired 样本名
# ============================================================
paired_samples      <- intersect(bac_dna_samples, fungi_dna_samples)
bac_unpaired        <- setdiff(bac_dna_samples,   fungi_dna_samples)   # 细菌有，真菌没有
fungi_unpaired      <- setdiff(fungi_dna_samples, bac_dna_samples)     # 真菌有，细菌没有

cat("Paired samples (Bac DNA ∩ Fungi DNA):", length(paired_samples), "\n")
cat("Bac DNA only (unpaired):", length(bac_unpaired), "\n")
cat("Fungi DNA only (unpaired):", length(fungi_unpaired), "\n")
cat("Bac RNA samples:", length(bac_rna_samples), "\n")

# ============================================================
# 4. 读取 meta 表（请修改路径）
# ============================================================
meta_file <- "/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/combined_meta.csv"   # ← 按实际路径修改
meta <- read.csv(meta_file, stringsAsFactors = FALSE, check.names = FALSE)
meta$Sample <- trimws(meta$Sample)   # 去除首尾空格，与 OTU 列名保持一致

# ============================================================
# 5. 提取 meta 子集并导出 CSV
# ============================================================
output_dir <- "/Users/shijiaoqi/Desktop/PHD_Materials/AirMicro/"   # ← 输出目录，按需修改

# (A) Paired 样本 meta
paired_meta <- meta %>% filter(Sample %in% paired_samples)
write.csv(paired_meta, file.path(output_dir, "pairedSample_meta.csv"), row.names = FALSE)
cat("pairedSample_meta.csv 已保存，共", nrow(paired_meta), "行\n")

# (B) 细菌 DNA Unpaired 样本 meta
bac_unpaired_meta <- meta %>% filter(Sample %in% bac_unpaired)
write.csv(bac_unpaired_meta, file.path(output_dir, "Bac_unpairedSample_meta.csv"), row.names = FALSE)
cat("Bac_unpairedSample_meta.csv 已保存，共", nrow(bac_unpaired_meta), "行\n")

# (C) 真菌 DNA Unpaired 样本 meta
fungi_unpaired_meta <- meta %>% filter(Sample %in% fungi_unpaired)
write.csv(fungi_unpaired_meta, file.path(output_dir, "Fungi_unpairedSample_meta.csv"), row.names = FALSE)
cat("Fungi_unpairedSample_meta.csv 已保存，共", nrow(fungi_unpaired_meta), "行\n")

# (D) 细菌 RNA 样本 meta
bac_rna_meta <- meta %>% filter(Sample %in% bac_rna_samples)
write.csv(bac_rna_meta, file.path(output_dir, "Bac_RNA_sample_meta.csv"), row.names = FALSE)
cat("Bac_RNA_sample_meta.csv 已保存，共", nrow(bac_rna_meta), "行\n")

# ============================================================
# 6. 诊断：OTU 中有但 meta 里没有的样本名
# ============================================================
check_missing <- function(label, samples) {
  missing <- setdiff(samples, meta$Sample)
  if (length(missing) > 0) {
    cat("\n⚠ [", label, "] 以下样本名在 meta 中未找到 (共", length(missing), "):\n")
    print(missing)
  } else {
    cat("\n✓ [", label, "] 所有样本名均在 meta 中找到\n")
  }
}

check_missing("Paired samples",    paired_samples)
check_missing("Bac DNA unpaired",  bac_unpaired)
check_missing("Fungi DNA unpaired",fungi_unpaired)
check_missing("Bac RNA samples",   bac_rna_samples)