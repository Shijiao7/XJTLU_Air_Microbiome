#microtable 对象meco，根据meco$sample_table$Group指等于group1_name和group2_name，取两个组，
#找出两个组共有的asv，做微生物DNA和RNA相关性分析，最后用spearman判断相关性可信度

# 加载必要的包
library(microeco)
library(ggplot2)
library(ggpubr)

##############################################################################>
# rna活性分析封装成函数  ----
##############################################################################>
f_analyze_rna_activity <-function(outDir,a_meco,group1_name,group2_name,taxa,main_core,main_core2) {
  # 假设你的microtable对象名为meco
  #a_meco<-meco
  # 假设你要比较的两个组名为group1_name和group2_name
  #group1_name <- "Sum_DNA_sam"  # 替换为实际的组名
  #group2_name <- "Sum_RNA"  # 替换为实际的组名
  #group1_name <- "Win_DNA_p"  # 替换为实际的组名
  #group2_name <- "Win_RNA_p"  # 替换为实际的组名
  #taxa="Species"  
  #taxa="Genus"
  
  
  # 将meco对象中的OTU表转换为相对丰度
  # 转换方法1
  tmp <- trans_norm$new(dataset = a_meco)
  meco_rel <- tmp$norm(method = "TSS")  #TSS cls   rarefy
  save_df_as_tsv(meco_rel$otu_table,"o2.tsv")
  
  # 转换方法2
  # meco_rel<-clone(a_meco)
  # otu_table<-meco_rel$otu_table
  # otu_table_rel <- apply(otu_table, 2, function(x) x / sum(x) )
  # meco_rel$otu_table <-  otu_table_rel
  # save_df_as_tsv(meco_rel$otu_table,"o2.tsv")
  
  #view(meco_rel$otu_table)
  
  ######################################################################>
  # Part1.取物种表 形成“asvID+指定物种级别”的数据框
  tax_table <- meco_rel$tax_table %>% tibble::rownames_to_column("ID")
  tax_table_taxa <- as.data.frame( tax_table[,c("ID",taxa)] ) # 仅仅保留ID和Genus列
  #数据框tax_table_taxa第2列，取_后面的替换，如为空填英文的未知
  library(stringr)
  tax_table_taxa <- tax_table_taxa %>%
    mutate(
      across(2, ~ {
        # 先提取__后面的内容
        extracted <- str_extract(.x, "(?<=__).+")
        # 判断条件：
        # 1. 如果原字符串包含__但提取结果为NA（即__后无内容），设为"unknown"
        # 2. 如果原字符串不包含__，设为"unknown"
        # 3. 否则使用提取结果
        case_when(
          grepl("__", .x) & is.na(extracted) ~ "unknown",
          !grepl("__", .x) ~ "unknown",
          TRUE ~ extracted
        )
      })
    )
  
  ######################################################################>
  # Part2.
  # 2-1. 提取两个组的样本
  sample_subset_dna <- meco_rel$sample_table[meco_rel$sample_table$Group %in% c(group1_name), ]
  sample_subset_rna <- meco_rel$sample_table[meco_rel$sample_table$Group %in% c(group2_name), ]
  # 2-2. 提取这两个组的ASV表(DNA和RNA)
  asv_table<-t(meco_rel$otu_table)
  #view(asv_table)
  dna_subset <- asv_table[rownames(sample_subset_dna), ]
  rna_subset <- asv_table[rownames(sample_subset_rna), ]  
  #转置并将asv关联Genus
  dna_subset_t<-as.data.frame( t(dna_subset) )
  rna_subset_t<-as.data.frame( t(rna_subset) )
  
  # 2-3.处理DNA和物种关联及取丰度的均值的log10+1
  dna_subset_tmp <- dna_subset_t %>% tibble::rownames_to_column("ID") # 将行名转换为显式列（left_join 需要列名）
  dna_subset_tmp <- left_join(dna_subset_tmp, tax_table_taxa, by = "ID") # 按 ID 列左连接
  dna_subset_tmp <- dna_subset_tmp %>% column_to_rownames(var = "ID") #将ID转成行名
  dna_subset_tmp$dna_row_means <- dna_subset_tmp %>% 
    dplyr::select(-all_of(taxa)) %>%  # 使用 all_of() 安全地排除列
    rowMeans(na.rm = TRUE)
  dna_subset_tmp$dna_mean_log10 <- log10(dna_subset_tmp$dna_row_means*100 + 1) #新增log10（平均值*100+1）列
  #20250730注释
  #dna_subset_tmp <- dna_subset_tmp %>% 
  #  filter(dna_mean_log10 != 0)  # 保留 log(row_mean+1) 不等于 0 的行
  
  # 20250730新增：将宽表转换为长表,并新增log10字段将每个样本的值*100+1后取log10
  dna_subset_tmp_long <- dna_subset_tmp %>%
    tibble::rownames_to_column("asv") %>%
    pivot_longer(cols = -c(asv, taxa,dna_row_means,dna_mean_log10), 
                 names_to = "sample_name", values_to = "dna") %>%
    filter(!near(dna, 0)) #将dna为零的行删除
  dna_subset_tmp_long$dna_log10<-  log10(dna_subset_tmp_long$dna*100+1) #新增log10（平均值*100+1）列
  
  # 2-4.处理RNA和物种关联及取丰度的均值
  rna_subset_tmp <- rna_subset_t %>% tibble::rownames_to_column("ID") # 将行名转换为显式列（left_join 需要列名）
  rna_subset_tmp <- left_join(rna_subset_tmp, tax_table_taxa, by = "ID") # 按 ID 列左连接
  rna_subset_tmp <- rna_subset_tmp %>% column_to_rownames(var = "ID") #将ID转成行名
  rna_subset_tmp$rna_row_means <- rna_subset_tmp %>% 
    dplyr::select(-all_of(taxa)) %>%  # 使用 all_of() 安全地排除列
    rowMeans(na.rm = TRUE)
  rna_subset_tmp$rna_mean_log10 <- log10(rna_subset_tmp$rna_row_means*100 + 1) #新增log10（平均值*100+1）列
  #20250730注释
  #rna_subset_tmp <- rna_subset_tmp %>% 
  # filter(rna_mean_log10 != 0)  # 保留 log(row_mean+1) 不等于 0 的行
  
  # 20250730新增：将宽表转换为长表,并新增log10字段将每个样本的值*100+1后取log10
  rna_subset_tmp_long <- rna_subset_tmp %>%
    tibble::rownames_to_column("asv") %>%
    pivot_longer(cols = -c(asv, taxa,rna_row_means,rna_mean_log10), 
                 names_to = "sample_name", values_to = "rna") %>%
    mutate(sample_name = str_remove(sample_name, "^R")) %>% #删除sample_name列开头的R字母，便于后续合并数据进行关联
    filter(!near(rna, 0)) #将rna为零的行删除
  rna_subset_tmp_long$rna_log10<-  log10(rna_subset_tmp_long$rna*100+1) #新增log10（平均值*100+1）列
  
  
  # 20250730新增： 合并数据框（仅保留共有的 asv + taxa + sample_name 组合）
  merged_dna_rna_long <- inner_join(
    dna_subset_tmp_long,
    rna_subset_tmp_long,
    by = c("asv", taxa,"sample_name")
  )
  merged_dna_rna_long<- merged_dna_rna_long %>%
    relocate(sample_name, .before = 1)  # 将sample_name放在第1列
  merged_dna_rna_long$rna_div_dna_log10<- log10( (merged_dna_rna_long$rna/merged_dna_rna_long$dna)+1 ) # 新增rna_div_dna_log10列=log10(rna除dna*100 + 1)
  merged_dna_rna_long$rna_div_dna<-merged_dna_rna_long$rna/merged_dna_rna_long$dna
  merged_dna_rna_long$rna_div_dna_log2<- log2( (merged_dna_rna_long$rna/merged_dna_rna_long$dna)+1 ) # 新增rna_div_dna_log10列=log10(rna除dna*100 + 1)
  # 
  # 2-5.找出共同存在asv
  common_rows <- intersect(rownames(dna_subset_tmp), rownames(rna_subset_tmp)) #找出共有行名
  #根据共有行名提取数据
  dna_common <- dna_subset_tmp[common_rows, ]
  rna_common <- rna_subset_tmp[common_rows, ]
  
  # 2-6.dna_common 和 rna_common 按行名进行合并 , 创建数据框用于相关性分析
  # 将行名转换为显式列（dplyr 需要列名才能合并）
  dna_common_tmp <- dna_common %>% tibble::rownames_to_column("asv_id")
  rna_common_tmp <- rna_common %>% tibble::rownames_to_column("asv_id")
  # 按 asv_id 合并（全连接）
  cor_data <- full_join(dna_common_tmp, rna_common_tmp, by = "asv_id")
  # 恢复行名
  cor_data <- cor_data %>% tibble::column_to_rownames("asv_id")
  cor_data$taxa<-cor_data[[paste0(taxa,".x")]] #新生成taxa变量的列
  
  #剔除两个丰度都为零的dna_mean_log10和rna_mean_log10
  cor_data <- cor_data %>%
    filter(!(dna_mean_log10 == 0 & rna_mean_log10 == 0))
  cor_data$rna_div_dna<-cor_data$rna_row_means/cor_data$dna_row_means
  cor_data$rna_div_dna_log10<-log10(cor_data$rna_div_dna+1)
  cor_data$rna_div_dna_log2<-log2(cor_data$rna_div_dna+1)
  cor_data2<-cor_data
  cor_data3<-cor_data
  pathogenDB<-read_tsv("data/pathogen_b.tsv",show_col_types = FALSE)
  if(taxa=='Species'){ #Species
    pathogens<-paste0("",pathogenDB$Genus,"_",pathogenDB$Species)
    cor_data3$is_pathogen<-cor_data3$Species.y %in% pathogens
    merged_dna_rna_long$is_pathogen<-merged_dna_rna_long$Species %in% pathogens
  } else {  #Genus
    pathogens<-paste0("",pathogenDB$Genus) #,"_",pathogenDB$Species)
    cor_data3$is_pathogen<-cor_data3$Genus.y %in% pathogens
    merged_dna_rna_long$is_pathogen<-merged_dna_rna_long$Genus %in% pathogens
  }
  
  merged_dna_rna_long$season<-season
  save_df_as_tsv(merged_dna_rna_long,
                 paste0(outDir,group1_name,"-",group2_name,"_",taxa,"_merged_dna_rna_long",".tsv") )
  merged_dna_rna_long<-read_tsv(paste0(outDir,group1_name,"-",group2_name,"_",taxa,"_merged_dna_rna_long",".tsv"))
  ###############################################################################
  ###############################################################################
  #save File and Fig 1
  ###############################################################################
  ###############################################################################
  #从 cor_data 数据框中：
  # a.按 taxa 列分组，计算 rna_mean_log10 的总和。
  # b.选取总和最大的 12 个 taxa，其余 taxa 替换为 "Minor/Unclassified"。
  # 1. 计算每个 taxa 的 rna_mean_log10 总和
  tmp_core<-cor_data %>% dplyr::select (taxa,dna_mean_log10,rna_mean_log10)
  tmp_core$Season<-substr(group1_name, 1, 3)
  tmp_core$Type<-"all"
  main_core$main_RNA<-rbind(main_core$main_RNA,tmp_core)
  taxa_sums <- cor_data %>%
    filter(taxa != "unknown") %>%  # 剔除 "unknown"
    group_by(taxa) %>%
    summarise(total_rna = sum(rna_mean_log10, na.rm = TRUE)) %>%
    arrange(desc(total_rna))  # 按总和降序排序
  # 2. 提取前 12 个 taxa 名称
  top12_taxa <- head(taxa_sums$taxa, 12)
  # 3. 在原始数据中，将非前 12 的 taxa 替换为 "Minor/Unclassified"
  cor_data <- cor_data %>%
    mutate(taxa = if_else(taxa %in% top12_taxa, taxa, "Minor/Unclassified"))
  save_df_as_tsv(cor_data %>% dplyr::select (taxa,dna_mean_log10,rna_mean_log10),
                 paste0(outDir,group1_name,"-",group2_name,"_",taxa,".tsv") )
  
  
  ############################################>
  # Part3. 计算DNA和RNA的相关性
  #        绘制散点图并计算Spearman相关性
  
  # 3.1 计算Spearman相关系数及其p值
  cor_test_result <- cor.test(cor_data$dna_mean_log10, cor_data$rna_mean_log10, method = "spearman")
  # 输出结果
  cat("Spearman correlation coefficient:", cor_test_result$estimate, "\n")
  cat("p-value:", cor_test_result$p.value, "\n")
  # 判断相关性是否显著
  if(cor_test_result$p.value < 0.05) {
    cat("The correlation is statistically significant (p < 0.05).\n")
  } else {
    cat("The correlation is not statistically significant (p >= 0.05).\n")
  }
  
  #   3.2 绘制散点图
  # 将 'Minor/UnGenusified' 放在最后，其他 taxa 按照原来的顺序排列
  # 明确指定要排除的模式（不区分大小写）
  unique_taxa<- unique(cor_data$taxa)
  other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
  cor_data$taxa <- factor(
    as.character(cor_data$taxa),
    levels = c(other_levels,"Minor/Unclassified"
    )
  )
  levels(cor_data$taxa)
  # 获取颜色调色板，并去掉灰色
  taxa_colors <- RColorBrewer::brewer.pal(12, "Set3")
  taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色
  
  # 创建颜色映射，保持 "Minor/Unclassified" 为灰色
  taxa_color_map <- c("Minor/Unclassified" = "gray", 
                      setNames(taxa_colors, setdiff(unique(cor_data$taxa), "Minor/Unclassified")))
  
  p <- ggplot(cor_data, aes(x = dna_mean_log10, y = rna_mean_log10,color=taxa)) +
    geom_point(alpha = 1)  +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
    # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
    # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
    geom_vline(xintercept = log10(0.005 + 1), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    # 新增水平虚线（RNA转录活性的参考线）
    geom_hline(yintercept = log10(0.005 + 1),  # 使用与垂直线相同的阈值
               linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste0("",group1_name," - ",group2_name," RNA Amount (All)"),
         x = "log10(relative aboundance + 1)",
         y = "log10(relative transcription  + 1)", 
         color = "taxa" )+  # 设置图例的标题为 "taxa"
    scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                       name = taxa ) +  # 使用手动颜色映射
    theme_bw()
  x_pos <- min(cor_data$dna_mean_log10, na.rm = TRUE)
  y_pos <- max(cor_data$rna_mean_log10, na.rm = TRUE)
  p <- p + annotate("text", 
                    x = x_pos, y = y_pos,
                    label = paste("p =", round(cor_test_result$estimate, 4),
                                  "\nP =", format(cor_test_result$p.value, scientific = TRUE, digits = 4)),
                    hjust = -0.2,  # 更大幅度左移
                    vjust = 2,     # 更大幅度下移
                    size = 3, color = "black")
  p <- p + coord_cartesian(clip = "off")
  # 显示图形
  print(p)
  
  ###############################################################################  
  ###############################################################################
  #save File and Fig 2
  ###############################################################################
  ###############################################################################
  taxa_sums <- cor_data3 %>%
    filter(is_pathogen != "FALSE") %>%  # 剔除 "FALSE"
    group_by(taxa) %>%
    summarise(total_rna = sum(rna_mean_log10, na.rm = TRUE)) %>%
    arrange(desc(total_rna))  # 按总和降序排序
  top10_taxa <- head(taxa_sums$taxa, 10)
  # 3. 在原始数据中，将非前 12 的 taxa 替换为 "Minor/Unclassified"
  cor_data3 <- cor_data3 %>%
    mutate(taxa = if_else(taxa %in% top10_taxa, taxa, "Minor/Unclassified"))
  save_df_as_tsv(cor_data3 %>% dplyr::select (taxa,dna_mean_log10,rna_mean_log10),
                 paste0(outDir,group1_name,"-",group2_name,"_",taxa,"_pathogen.tsv") )
  
  ##########################################
  # Part3. 计算DNA和RNA的相关性
  #        绘制散点图并计算Spearman相关性
  
  # 3.1 计算Spearman相关系数及其p值
  cor_test_result <- cor.test(cor_data3$dna_mean_log10, cor_data3$rna_mean_log10, method = "spearman")
  # 输出结果
  cat("Spearman correlation coefficient:", cor_test_result$estimate, "\n")
  cat("p-value:", cor_test_result$p.value, "\n")
  # 判断相关性是否显著
  if(cor_test_result$p.value < 0.05) {
    cat("The correlation is statistically significant (p < 0.05).\n")
  } else {
    cat("The correlation is not statistically significant (p >= 0.05).\n")
  }
  
  #   3.2 绘制散点图
  # 将 'Minor/UnGenusified' 放在最后，其他 taxa 按照原来的顺序排列
  # 明确指定要排除的模式（不区分大小写）
  unique_taxa<- unique(cor_data3$taxa)
  other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
  cor_data3$taxa <- factor(
    as.character(cor_data3$taxa),
    levels = c(other_levels,"Minor/Unclassified"
    )
  )
  levels(cor_data3$taxa)
  # 获取颜色调色板，并去掉灰色
  taxa_colors <- RColorBrewer::brewer.pal(10, "Set3")
  taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色
  
  # 创建颜色映射，保持 "Minor/Unclassified" 为灰色
  taxa_color_map <- c("Minor/Unclassified" = "gray", 
                      setNames(taxa_colors, setdiff(unique(cor_data3$taxa), "Minor/Unclassified")))
  
  p <- ggplot(cor_data3, aes(x = dna_mean_log10, y = rna_mean_log10,color=taxa)) +
    geom_point(alpha = 1)  +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
    # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
    # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
    geom_vline(xintercept = log10(0.005 + 1), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    # 新增水平虚线（RNA转录活性的参考线）
    geom_hline(yintercept = log10(0.005 + 1),  # 使用与垂直线相同的阈值
               linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste0("",group1_name," - ",group2_name," RNA Amount (pathogen)"),
         x = "log10(relative aboundance + 1)",
         y = "log10(relative transcription + 1)", 
         color = "taxa" )+  # 设置图例的标题为 "taxa"
    scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                       name = taxa ) +  # 使用手动颜色映射
    theme_bw()
  x_pos <- min(cor_data3$dna_mean_log10, na.rm = TRUE)
  y_pos <- max(cor_data3$rna_mean_log10, na.rm = TRUE)
  p <- p + annotate("text", 
                    x = x_pos, y = y_pos,
                    label = paste("p =", round(cor_test_result$estimate, 4),
                                  "\nP =", format(cor_test_result$p.value, scientific = TRUE, digits = 4)),
                    hjust = -0.2,  # 更大幅度左移
                    vjust = 2,     # 更大幅度下移
                    size = 3, color = "black")
  p <- p + coord_cartesian(clip = "off")
  # 显示图形
  print(p)
  
  
  ###############################################################################
  #save File and Fig 3
  ###############################################################################
  cor_data<-cor_data2
  tmp_core<-cor_data %>% dplyr::select (taxa,dna_mean_log10,rna_div_dna_log10)
  tmp_core$Season<-substr(group1_name, 1, 3)
  main_core$main_RNA_div_DNA<-rbind(main_core$main_RNA_div_DNA,tmp_core)
  taxa_sums <- cor_data %>%
    filter(taxa != "unknown") %>%  # 剔除 "unknown"
    group_by(taxa) %>%
    summarise(total_rna = sum(rna_div_dna_log10, na.rm = TRUE)) %>%
    arrange(desc(total_rna))  # 按总和降序排序
  # 2. 提取前 12 个 taxa 名称
  top12_taxa <- head(taxa_sums$taxa, 12)
  cor_data <- cor_data %>%
    mutate(taxa = if_else(taxa %in% top12_taxa, taxa, "Minor/Unclassified"))
  
  save_df_as_tsv(cor_data %>% dplyr::select (taxa,dna_mean_log10,rna_div_dna_log10,rna_div_dna_log2,rna_div_dna),
                 paste0(outDir,group1_name,"-",group2_name,"_",taxa,"_RNAdivDNA.tsv") )
  
  
  # 3.1 计算Spearman相关系数及其p值
  cor_test_result <- cor.test(cor_data$dna_mean_log10, cor_data$rna_div_dna_log10, method = "spearman")
  # 输出结果
  cat("Spearman correlation coefficient:", cor_test_result$estimate, "\n")
  cat("p-value:", cor_test_result$p.value, "\n")
  # 判断相关性是否显著
  if(cor_test_result$p.value < 0.05) {
    cat("The correlation is statistically significant (p < 0.05).\n")
  } else {
    cat("The correlation is not statistically significant (p >= 0.05).\n")
  }
  
  #   3.2 绘制散点图
  # 将 'Minor/UnGenusified' 放在最后，其他 taxa 按照原来的顺序排列
  # 明确指定要排除的模式（不区分大小写）
  unique_taxa<- unique(cor_data$taxa)
  other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
  cor_data$taxa <- factor(
    as.character(cor_data$taxa),
    levels = c(other_levels,"Minor/Unclassified"
    )
  )
  levels(cor_data$taxa)
  # 获取颜色调色板，并去掉灰色
  taxa_colors <- RColorBrewer::brewer.pal(12, "Set3")
  taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色
  
  # 创建颜色映射，保持 "Minor/Unclassified" 为灰色
  taxa_color_map <- c("Minor/Unclassified" = "gray", 
                      setNames(taxa_colors, setdiff(unique(cor_data$taxa), "Minor/Unclassified")))
  
  p <- ggplot(cor_data, aes(x = dna_mean_log10, y = rna_div_dna_log10,color=taxa)) +
    geom_point(alpha = 1)  +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
    # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
    # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
    geom_vline(xintercept = log10(0.005 + 1), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    # 新增水平虚线（RNA转录活性的参考线）
    geom_hline(yintercept = log10(2), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste0("",group1_name," - ",group2_name," RNA/DNA  Activity (All)"),
         x = "log10(relative aboundance + 1)",
         y = "log10(relative transcription/aboundance + 1)", 
         color = "taxa" )+  # 设置图例的标题为 "taxa"
    scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                       name = taxa ) +  # 使用手动颜色映射
    theme_bw()
  p <- p + annotate("text", 
                    x = x_pos, y = y_pos,
                    label = paste("p =", round(cor_test_result$estimate, 4),
                                  "\nP =", format(cor_test_result$p.value, scientific = TRUE, digits = 4)),
                    hjust = -0.2,  # 更大幅度左移
                    vjust = 2,     # 更大幅度下移
                    size = 3, color = "black")
  p <- p + coord_cartesian(clip = "off")
  # 显示图形
  print(p)
  
  
  ###############################################################################
  #save File and Fig 4
  ###############################################################################
  cor_data4<-cor_data3
  taxa_sums <- cor_data4 %>%
    filter(is_pathogen != "FALSE") %>%  # 剔除 "FALSE"
    group_by(taxa) %>%
    summarise(total_rna = sum(rna_div_dna_log10, na.rm = TRUE)) %>%
    arrange(desc(total_rna))  # 按总和降序排序
  # 2. 提取前 12 个 taxa 名称
  top10_taxa <- head(taxa_sums$taxa, 10)
  cor_data4 <- cor_data4 %>%
    mutate(taxa = if_else(taxa %in% top10_taxa, taxa, "Minor/Unclassified"))
  save_df_as_tsv(cor_data4 %>% dplyr::select (taxa,dna_mean_log10,rna_div_dna_log10,rna_div_dna_log2,rna_div_dna),
                 paste0(outDir,group1_name,"-",group2_name,"_",taxa,"RNAdivDNA.tsv") )
  # 3.1 计算Spearman相关系数及其p值
  cor_test_result <- cor.test(cor_data3$dna_mean_log10, cor_data3$rna_div_dna_log10, method = "spearman")
  # 输出结果
  cat("Spearman correlation coefficient:", cor_test_result$estimate, "\n")
  cat("p-value:", cor_test_result$p.value, "\n")
  # 判断相关性是否显著
  if(cor_test_result$p.value < 0.05) {
    cat("The correlation is statistically significant (p < 0.05).\n")
  } else {
    cat("The correlation is not statistically significant (p >= 0.05).\n")
  }
  
  #   3.2 绘制散点图
  # 将 'Minor/UnGenusified' 放在最后，其他 taxa 按照原来的顺序排列
  # 明确指定要排除的模式（不区分大小写）
  unique_taxa<- unique(cor_data4$taxa)
  other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
  cor_data4$taxa <- factor(
    as.character(cor_data4$taxa),
    levels = c(other_levels,"Minor/Unclassified"
    )
  )
  levels(cor_data4$taxa)
  # 获取颜色调色板，并去掉灰色
  taxa_colors <- RColorBrewer::brewer.pal(10, "Set3")
  taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色
  
  # 创建颜色映射，保持 "Minor/Unclassified" 为灰色
  taxa_color_map <- c("Minor/Unclassified" = "gray", 
                      setNames(taxa_colors, setdiff(unique(cor_data4$taxa), "Minor/Unclassified")))
  
  p <- ggplot(cor_data4, aes(x = dna_mean_log10, y = rna_div_dna_log10,color=taxa)) +
    geom_point(alpha = 1)  +
    geom_smooth(method = "loess", se = TRUE, color = "blue") +
    # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
    # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
    # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
    geom_vline(xintercept = log10(0.005 + 1), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    # 新增水平虚线（RNA转录活性的参考线）
    geom_hline(yintercept = log10(2), 
               linetype = "dashed", color = "black", linewidth = 0.5) +
    scale_x_log10() +
    scale_y_log10() +
    labs(title = paste0("",group1_name," - ",group2_name," RNA/DNA Activity (Pathogen)"),
         x = "log10(relative aboundance + 1)",
         y = "log10(relative transcription/aboundance + 1)", 
         color = "taxa" )+  # 设置图例的标题为 "taxa"
    scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                       name = taxa ) +  # 使用手动颜色映射
    theme_bw()
  p <- p + annotate("text", 
                    x = x_pos, y = y_pos,
                    label = paste("p =", round(cor_test_result$estimate, 4),
                                  "\nP =", format(cor_test_result$p.value, scientific = TRUE, digits = 4)),
                    hjust = -0.2,  # 更大幅度左移
                    vjust = 2,     # 更大幅度下移
                    size = 3, color = "black")
  p <- p + coord_cartesian(clip = "off")
  # 显示图形
  print(p)
  mains<-list(main_core=main_core,main_core2=merged_dna_rna_long)
  #return(main_core)
  
  
  
  
  
  
  return(mains)
}







################################################################################>
# 主程序开始 ----
################################################################################>
#group1_name="aAut_p"
#group2_name="bWin_p"
#不需要了 Groups<-f_get_sampleFilterGroupNameByOrder()  ##增加序号a b c  
Groups<-getSampleArg$sampleFilterGroupName  #取分组

# 按照顺序生成两两配对，如分组是：A、B、C、D，则生成两对 A对B、C对D
pairs <- list()
for (i in seq(1, length(Groups), by = 2)) {
  # 如果剩余元素不足2个则跳过
  if (i + 1 <= length(Groups)) {
    pair <- c(Groups[i], Groups[i + 1])
    # 将配对结果添加到列表
    pairs[[paste("Pair", (i + 1) / 2)]] <- pair
  }
}
#############################################################################>
# 生成组图1：计算的是两组样本（group1_name 和 group2_name）的相对丰度的对数转换值
#         （log10(Relative Abundance + 1)），并直接比较这两组数据的相关性。
#taxa<-"Species" 
taxa<-"Genus" 
outDir<-paste0("result/",getSampleArg$microbeType,"/",getSampleArg$ver,"/",
               getSampleArg$samplesName,"/","02Plot/20DNA2RNA/",taxa,"/")
dir.create(outDir,recursive = TRUE,showWarnings = FALSE)
print(paste0("输出目录为：",outDir))
pdf(paste0(outDir,f_microbeType(),"_",getSampleArg$samplesName,"_rna_activity(",taxa,").pdf"), 
    width = 8, height = 5, pointsize = 300/72)
main_RNA<-data.frame()
main_RNA_div_DNA<-data.frame()
main_core<-list(main_RNA = main_RNA,
                main_RNA_div_DNA = main_RNA_div_DNA)
main_core2<-list()
for (i in 1:length(pairs)) {
  pair <- pairs[[i]]
  group1_name <- pair[1]
  group2_name <- pair[2]
  print(paste("处理函数1，对以下两组对比:",group1_name," 与 ", group2_name ) )
  if(i==1){
    season<-"Winter"
  }else{
    season<-"Plum_Rain"
  }
  tryCatch( {
    mains<-f_analyze_rna_activity(outDir,meco,group1_name,group2_name,taxa,main_core) #最后一个参数是物种
    main_core<-mains$main_core
    main_core2[[i]]<-mains$main_core2
  }, error = function(e) {
    message(paste(group1_name," 与 ", group2_name,"对比，通过f_spearman_plot生成对比图时，错误发生了:"))
    print(e)
  })
}
graphics.off()  # 关闭所有打开的图形设备（包括PDF、PNG等）


core_win<-main_core2[[1]]
cor_test_result_win <- cor.test(core_win$dna, core_win$rna, method = "spearman")
core_win_major<-subset(core_win, dna > 0.005)
cor_test_result_win_major <- cor.test(core_win_major$dna, core_win_major$rna, method = "spearman")
cor_test_result_win2 <- cor.test(core_win$dna, core_win$rna_div_dna, method = "spearman")
cor_test_result_win_major2 <- cor.test(core_win_major$dna, core_win_major$rna_div_dna, method = "spearman")


core_sum<-main_core2[[2]]
cor_test_result_sum <- cor.test(core_sum$dna, core_sum$rna, method = "spearman")
core_sum_major<-subset(core_sum, dna > 0.005)
cor_test_result_sum_major <- cor.test(core_sum_major$dna, core_sum_major$rna, method = "spearman")
cor_test_result_sum2 <- cor.test(core_sum$dna, core_sum$rna_div_dna, method = "spearman")
cor_test_result_sum_major2 <- cor.test(core_sum_major$dna, core_sum_major$rna_div_dna, method = "spearman")

pdf(paste0(outDir,taxa,"main.pdf"), 
    width = 12,   # 宽度（默认 7 英寸）
    height = 6)
core_RNA<-main_core$main_RNA
taxa_sums <- core_RNA %>%
  filter(taxa != "unknown") %>%  # 剔除 "unknown"
  group_by(taxa) %>%
  summarise(total_rna = sum(rna_mean_log10, na.rm = TRUE)) %>%
  arrange(desc(total_rna))  # 按总和降序排序
# 2. 提取前 12 个 taxa 名称
top12_taxa <- head(taxa_sums$taxa, 12)
core_RNA <- core_RNA %>%
  mutate(taxa = if_else(taxa %in% top12_taxa, taxa, "Minor/Unclassified"))
unique_taxa<- unique(core_RNA$taxa)
other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
core_RNA$taxa <- factor(
  as.character(core_RNA$taxa),
  levels = c(other_levels,"Minor/Unclassified"
  )
)
levels(core_RNA$taxa)
# 获取颜色调色板，并去掉灰色
taxa_colors <- RColorBrewer::brewer.pal(12, "Set3")
taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色

# 创建颜色映射，保持 "Minor/Unclassified" 为灰色
taxa_color_map <- c("Minor/Unclassified" = "gray", 
                    setNames(taxa_colors, setdiff(unique(core_RNA$taxa), "Minor/Unclassified")))
# facet_grid(facets = Season ~ .)+
#   facet_grid(Season ~ .)+
#   facet_grid(rows = vars(Season))+
#   facet_grid(rows = "Season")+
#   facet_grid(facets = ~ Type)+
#   facet_grid(. ~ Type)+
#   facet_grid(cols = vars(Type))+
p <- ggplot(core_RNA, aes(x = dna_mean_log10, y = rna_mean_log10,color=taxa)) +
  facet_grid(. ~ Season)+
  geom_point(alpha = 1)  +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
  # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
  # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
  geom_vline(xintercept = log10(0.005 + 1), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  # 新增水平虚线（RNA转录活性的参考线）
  geom_hline(yintercept = log10(0.005 + 1),  # 使用与垂直线相同的阈值
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = paste0("",group1_name," - ",group2_name," RNA Amount (All)"),
       x = "log10(relative aboundance + 1)",
       y = "log10(relative transcription  + 1)", 
       color = "taxa" )+  # 设置图例的标题为 "taxa"
  scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                     name = taxa ) +  # 使用手动颜色映射
  theme_bw()
x_pos <- min(core_RNA$dna_mean_log10, na.rm = TRUE)
y_pos <- max(core_RNA$rna_mean_log10, na.rm = TRUE)
p <- p + annotate("text", 
                  x = x_pos, y = y_pos,
                  label = paste("",""),
                  hjust = -0.2,  # 更大幅度左移
                  vjust = 2,     # 更大幅度下移
                  size = 3, color = "black")
p <- p + coord_cartesian(clip = "off")
# 显示图形
print(p)
graphics.off() 



pdf(paste0(outDir,taxa,"main2.pdf"), 
    width = 12,   # 宽度（默认 7 英寸）
    height = 6)
core_RNA_div_DNA<-main_core$main_RNA_div_DNA
taxa_sums <- core_RNA_div_DNA %>%
  filter(taxa != "unknown") %>%  # 剔除 "unknown"
  group_by(taxa) %>%
  summarise(total_rna = sum(rna_div_dna_log10, na.rm = TRUE)) %>%
  arrange(desc(total_rna))  # 按总和降序排序
# 2. 提取前 12 个 taxa 名称
top12_taxa <- head(taxa_sums$taxa, 12)
core_RNA_div_DNA <- core_RNA_div_DNA %>%
  mutate(taxa = if_else(taxa %in% top12_taxa, taxa, "Minor/Unclassified"))
unique_taxa<- unique(core_RNA_div_DNA$taxa)
other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
core_RNA_div_DNA$taxa <- factor(
  as.character(core_RNA_div_DNA$taxa),
  levels = c(other_levels,"Minor/Unclassified"
  )
)
levels(core_RNA_div_DNA$taxa)
# 获取颜色调色板，并去掉灰色
taxa_colors <- RColorBrewer::brewer.pal(12, "Set3")
taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色

# 创建颜色映射，保持 "Minor/Unclassified" 为灰色
taxa_color_map <- c("Minor/Unclassified" = "gray", 
                    setNames(taxa_colors, setdiff(unique(core_RNA_div_DNA$taxa), "Minor/Unclassified")))
# facet_grid(facets = Season ~ .)+
#   facet_grid(Season ~ .)+
#   facet_grid(rows = vars(Season))+
#   facet_grid(rows = "Season")+
#   facet_grid(facets = ~ Type)+
#   facet_grid(. ~ Type)+
#   facet_grid(cols = vars(Type))+
p <- ggplot(core_RNA_div_DNA, aes(x = dna_mean_log10, y = rna_div_dna_log10,color=taxa)) +
  facet_grid(. ~ Season)+
  geom_point(alpha = 1)  +
  geom_smooth(method = "loess", se = TRUE, color = "blue") +
  # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
  # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
  # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
  geom_vline(xintercept = log10(0.005 + 1), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  # 新增水平虚线（RNA转录活性的参考线）
  geom_hline(yintercept = log10(0.005 + 1),  # 使用与垂直线相同的阈值
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = paste0("",group1_name," - ",group2_name," RNA_div_DNA (All)"),
       x = "log10(relative aboundance + 1)",
       y = "log10(relative transcription  + 1)", 
       color = "taxa" )+  # 设置图例的标题为 "taxa"
  scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                     name = taxa ) +  # 使用手动颜色映射
  theme_bw()
x_pos <- min(core_RNA_div_DNA$dna_mean_log10, na.rm = TRUE)
y_pos <- max(core_RNA_div_DNA$rna_div_dna_log10, na.rm = TRUE)
p <- p + annotate("text", 
                  x = x_pos, y = y_pos,
                  label = paste("",""),
                  hjust = -0.2,  # 更大幅度左移
                  vjust = 2,     # 更大幅度下移
                  size = 3, color = "black")
p <- p + coord_cartesian(clip = "off")
# 显示图形
print(p)
graphics.off() 


pdf(paste0(outDir,taxa,"new_main.pdf"), 
    width = 12,   # 宽度（默认 7 英寸）
    height = 6)
main_core3<-rbind(main_core2[[1]],main_core2[[2]])
main_core3$season <- factor(main_core3$season, 
                            levels = c("Winter","Plum_Rain"))

taxa_sums <- main_core3 %>%
  filter(Genus != "unknown") %>%
  group_by(Genus) %>%
  summarise(total_rna = sum(rna, na.rm = TRUE)) %>%
  arrange(desc(total_rna))
top12_taxa <- head(taxa_sums$Genus, 12)
main_core3 <- main_core3 %>%
  mutate(Genus = if_else(Genus %in% top12_taxa, Genus, "Minor/Unclassified"))

unique_taxa<- unique(main_core3$Genus)
other_levels <- sort(unique_taxa[!grepl("minor/unclassified", unique_taxa, ignore.case = TRUE)])
main_core3$Genus <- factor(
  as.character(main_core3$Genus),
  levels = c(other_levels,"Minor/Unclassified"
  )
)
levels(main_core3$Genus)
# 获取颜色调色板，并去掉灰色
taxa_colors <- RColorBrewer::brewer.pal(12, "Set3")
taxa_colors <- taxa_colors[taxa_colors != "gray"]  # 排除灰色

# 创建颜色映射，保持 "Minor/Unclassified" 为灰色
taxa_color_map <- c("Minor/Unclassified" = "gray", 
                    setNames(taxa_colors, setdiff(unique(main_core3$Genus), "Minor/Unclassified")))

seasons <- unique(main_core3$season)
main_core3_major<-subset(main_core3,dna>0.005)
cor_summary <- data.frame(season = character(), rho = numeric(), p_value = numeric())

for (s in seasons) {
  subset_data <- subset(main_core3, season == s)
  test_result <- cor.test(subset_data$dna, subset_data$rna, method = "spearman")
  cor_summary <- rbind(cor_summary, data.frame(
    season = s,
    rho = test_result$estimate,
    p_value = test_result$p.value
  ))
}
row.names(cor_summary)<-cor_summary$season

cor_summary2 <- data.frame(season = character(), rho = numeric(), p_value = numeric())

for (s in seasons) {
  subset_data <- subset(main_core3_major, season == s)
  test_result <- cor.test(subset_data$dna, subset_data$rna, method = "spearman")
  cor_summary2 <- rbind(cor_summary2, data.frame(
    season = s,
    rho_major = test_result$estimate,
    p_value_major = test_result$p.value
  ))
}
row.names(cor_summary2)<-cor_summary2$season
cor_summary2 <- subset(cor_summary2, select = -which(names(cor_summary2) == "season"))
cor_summary<-cbind(cor_summary,cor_summary2)


p<-ggplot(main_core3, aes(x = dna_log10, y = rna_log10,color=Genus)) +
  facet_grid(. ~ season)+
  geom_point(alpha = 1)  +
  geom_smooth(method = "loess", se = TRUE, color = "blue")
x_pos <- min(main_core3$dna_log10, na.rm = TRUE)
y_pos <- max(main_core3$rna_log10, na.rm = TRUE)
p <- p + geom_text(
  data = cor_summary,
  aes(x = x_pos, y = y_pos, label = paste("r =", round(rho, 4),"\nP =", signif(p_value, 3))),
  hjust = -0.2, vjust = 2, size = 2, color = "black",
  inherit.aes = FALSE  # 避免继承color=Genus的映射
)+  
  # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
  # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
  # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
  #geom_vline(xintercept = log10(1.5), 
  #           linetype = "dashed", color = "black", linewidth = 0.5) +
  # 新增水平虚线（RNA转录活性的参考线）
  #geom_hline(yintercept = log10(1.5), 
  #          linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RNA vs DNA",
       x = "log10(relative aboundance% + 1)",
       y = "log10(relative transcription%+ 1)", 
       color = "taxa" )+  # 设置图例的标题为 "taxa"
  scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                     name = taxa ) +  # 使用手动颜色映射
  theme_bw()
p <- p + coord_cartesian(clip = "off")
# 显示图形
print(p)


cor_summary_div <- data.frame(season = character(), rho = numeric(), p_value = numeric())

for (s in seasons) {
  subset_data <- subset(main_core3, season == s)
  test_result <- cor.test(subset_data$dna, subset_data$rna_div_dna, method = "spearman")
  cor_summary_div <- rbind(cor_summary_div, data.frame(
    season = s,
    rho = test_result$estimate,
    p_value = test_result$p.value
  ))
}
row.names(cor_summary_div)<-cor_summary_div$season

cor_summary_div2 <- data.frame(season = character(), rho = numeric(), p_value = numeric())

for (s in seasons) {
  subset_data <- subset(main_core3_major, season == s)
  test_result <- cor.test(subset_data$dna, subset_data$rna_div_dna, method = "spearman")
  cor_summary_div2 <- rbind(cor_summary_div2, data.frame(
    season = s,
    rho_major = test_result$estimate,
    p_value_major = test_result$p.value
  ))
}
row.names(cor_summary_div2)<-cor_summary_div2$season
cor_summary_div2 <- subset(cor_summary_div2, select = -which(names(cor_summary_div2) == "season"))
cor_summary_div<-cbind(cor_summary_div,cor_summary_div2)
p<-ggplot(main_core3, aes(x = dna_log10, y = rna_div_dna_log10,color=Genus)) +
  facet_grid(. ~ season)+
  geom_point(alpha = 1)  +
  geom_smooth(method = "loess", se = TRUE, color = "blue") 
x_pos <- 0.8
y_pos <- max(main_core3$rna_div_dna_log10, na.rm = TRUE)
p <- p + geom_text(
  data = cor_summary_div,
  aes(x = x_pos, y = y_pos, label = paste("r =", round(rho, 4), "\nP =", signif(p_value, 3))),
  hjust = -0.2, vjust = 2, size = 2, color = "black",
  inherit.aes = FALSE  # 避免继承color=Genus的映射
)+ 
  # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
  # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
  # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
  #geom_vline(xintercept = log10(1.5), 
  #           linetype = "dashed", color = "black", linewidth = 0.5) +
  # 新增水平虚线（RNA转录活性的参考线）
  geom_hline(yintercept = log10(2), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RNA/DNA vs DNA",
       x = "log10(relative aboundance% + 1)",
       y = "log10(relative transcription/abundance+ 1)", 
       color = "taxa" )+  # 设置图例的标题为 "taxa"
  scale_color_manual(values = c("gray" = "gray", taxa_color_map),
                     name = taxa ) +  # 使用手动颜色映射
  theme_bw()
p <- p + coord_cartesian(clip = "off")
# 显示图形
print(p)
graphics.off() 

pdf(paste0(outDir,taxa,"_RNA_div_DNA2.pdf"))
Win_r_d<-read_tsv(paste0(outDir,"Win_DNA_sam-Win_RNA_",taxa,"RNAdivDNA.tsv"))
Win_r_d$win_rna_div_dna<-Win_r_d$rna_div_dna_log2
Sum_r_d<-read_tsv(paste0(outDir,"Sum_DNA_sam-Sum_RNA_",taxa,"RNAdivDNA.tsv"))
Sum_r_d$sum_rna_div_dna<-Sum_r_d$rna_div_dna_log2
W_R<-merge(Win_r_d,Sum_r_d,by="id")
W_R<-W_R%>%
  filter(!(dna_mean_log10.x == 0 | win_rna_div_dna == 0 | dna_mean_log10.y == 0 | sum_rna_div_dna == 0))
W_R$taxa<-W_R$taxa.x

common_taxa <- unique(W_R$taxa[W_R$taxa != "no_pathogen"])    # Get unique taxa excluding "no_pathogen"
ordered_levels <- c(sort(setdiff(common_taxa, "Minor/Unclassified")), 
                    "Minor/Unclassified")
#将Minor/Unclassified 设置为灰色
main_taxa <- setdiff(common_taxa, "Minor/Unclassified")
taxa_colors <- get_sci_palette(main_taxa)
# Create the color mapping
taxa_color_map <- c("no_pathogen" = "gray", 
                    setNames(taxa_colors, main_taxa),
                    "Minor/Unclassified" = "grey50")
p <- ggplot(W_R, aes(x = win_rna_div_dna, y = sum_rna_div_dna,
                     color= factor(taxa, levels = ordered_levels)  )) +
  geom_point(alpha = 1)  +
  # 添加垂直线：log10(0.005 + 1) ≈ -2.3（若横轴是log10(relative abundance +1)）
  # 垂直虚线（位于横坐标 log₁₀(相对丰度 +1) ≈ -0.3 处，对应原始相对丰度0.5）
  # 用于划分微生物的高丰度群（≥0.5%) 和低丰度群（<0.5%）
  geom_vline(xintercept = log2(2), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  # 新增水平虚线（RNA转录活性的参考线）
  geom_hline(yintercept = log2(2), 
             linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "Winter and Summer transcription",
       x = "log2(Winter relative transcription/aboundance + 1)",
       y = "log2(Summer relative transcription/aboundance + 1)", 
       color = "taxa" )+  # 设置图例的标题为 "taxa"
  scale_color_manual(values = taxa_color_map ,
                     breaks = ordered_levels,  # 控制图例顺序
                     name = "Microbial Taxa") +  # 使用手动颜色映射
  theme_bw()
p <- p + coord_cartesian(clip = "off")
# 显示图形
print(p)
graphics.off()  # 关闭所有打开的图形设备（包括PDF、PNG等)

