# R 语言
# 1、读下面目录（result/bacteria/20241019/15 group_by_floor/02Plot/06sloan/all）下面的所有子目录，子目录名称为“Group1-Group2”，将读到的写进
# result/net_bacteria2fungi/15 group_by_floor/sloan目录下面Excel文件（migration summarize.xlsx），
# Excel文件子表为“bacteria”，写入头两列"source", "sink"，“15 group_by_floor"是输入目录
# 2、读下面目录（result/fungi/20241020/15 group_by_floor/02Plot/06sloan/all）下面的所有子目录，子目录名称为“Group1-Group2”，将读到的写进
# result/net_bacteria2fungi/15 group_by_floor/sloan目录下面Excel文件（migration summarize.xlsx），
# Excel文件子表为“bacteria”，写入头两列"source", "sink"，“15 group_by_floor"是输入目录
# 3、再读子目录下面的“residence_stats.tsv"文件，如没有则"migration_rate", "R_square"为空字符串，如找到文件，
#    将文件的m列和Rsqr列分别放到"migration_rate"和"R_square"列
# 4、再读子目录下面的“residence_result.tsv“文件，如没有则"Above"、"Neutral"、”Below"为空字符串，如找到文件，
#    按文件Partition列统计"Above"、"Neutral"、”Below"出现次数，统计后分别放到"Above"、"Neutral"、”Below"列
# 5、再弄一页画拼图，每行3个，用每行"Above"、"Neutral"、”Below"各列占比画，Above绿色、Neutral暗蓝色、Below红色，
#    标题为组名，封装成函数，输入参数是数据框，存储文件名是sheet_name_饼图.pdf

# 加载必要的包
library(openxlsx)
library(stringr)
library(ggplot2)
library(gridExtra)

#######################################################>
#Part1  画饼图 ----
draw_pie_charts <- function(result_df, output_pdf) {
  # 确保输出目录存在
  output_dir <- dirname(output_pdf)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 开启PDF设备（调整PDF尺寸适应小饼图）
  pdf(output_pdf, width = 10, height = ceiling(nrow(result_df)/5) )  # 增加每行显示数量
  
  # 创建绘图列表
  plot_list <- list()
  
  # 为每一行数据创建饼图
  for (i in 1:nrow(result_df)) {
    # 提取数据
    row_data <- result_df[i, ]
    group_name <- paste(row_data$source, row_data$sink, sep = "-")
    
    # 准备饼图数据并设置因子水平保证顺序
    pie_data <- data.frame(
      category = factor(c("Above", "Neutral", "Below"), 
                        levels = c("Above", "Neutral", "Below")),
      value = as.numeric(c(row_data$Above, row_data$Neutral, row_data$Below)),
      color = c("#4DAF4A", "#377EB8", "#E41A1C")  # 绿色、蓝色、红色
    )
    
    # 计算百分比
    total <- sum(pie_data$value)
    if (total > 0) {
      pie_data$percent <- paste0(round(pie_data$value/total * 100, 0), "%")
    } else {
      pie_data$percent <- "0%"
    }
    
    # 创建实心饼图（无空心）
    p <- ggplot(pie_data, aes(x = 1, y = value, fill = category)) +
      geom_col(color = NA, width = 1) +  # 设置width=1确保完全填充
      coord_polar("y", start = 0) +
      scale_fill_manual(values = setNames(pie_data$color, levels(pie_data$category))) +
      geom_text(
        aes(label = ifelse(value > 0, percent, "")),  # 只显示非零值标签
        position = position_stack(vjust = 0.5),
        size = 2.5,  # 缩小标签字号
        color = "black"  # 白色标签提高对比度
      ) +
      labs(title = group_name) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 8),  # 缩小标题
        legend.position = "right",
        legend.title = element_text(size = 8),  # 图例标题大小
        legend.text = element_text(size = 6),  # 图例文字大小
        plot.margin = unit(c(0, 0, 0, 0), "cm"),  # 完全去除边距
        panel.spacing = unit(0, "cm")  # 去除面板间距
      )
    
    plot_list[[i]] <- p
  }
  
  # 每行显示5个饼图（可根据需要调整ncol值）
  grid.arrange(
    grobs = plot_list, 
    ncol = 5,  # 增加每行饼图数量
    padding = unit(0, "cm")  # 去除网格间距
  )
  
  # 关闭PDF设备
  dev.off()
}


#######################################################>
#Part2  读取目录中各分组文件拼接成excel表（各组扩散率、R平方 ----
#       和样本在图中各区域个数（Above, Neutral, Below），
#       最后调用Part1函数画饼图
process_directories <- function(srcType,exclude_dirs,input_dir, output_file, sheet_name) {
  # 获取所有子目录
  subdirs <- list.dirs(input_dir, full.names = TRUE, recursive = FALSE)
  
  # 筛选出符合"X-Y"模式的目录（其中X和Y可以是任意字符）并排除exclude_dirs指定的目录
  group_dirs <- subdirs[grepl("[^-]+-[^-]+$", basename(subdirs)) & 
                          !(basename(subdirs) %in% exclude_dirs)]
  # 初始化结果数据框
  result_df <- data.frame(
    source = character(),
    sink = character(),
    migration_rate = character(),
    R_square = character(),
    Above = character(),
    Neutral = character(),
    Below = character(),
    stringsAsFactors = FALSE
  )
  
  # 遍历每个符合条件的目录
  for (dir in group_dirs) {
    # 提取source和sink
    dir_name <- basename(dir)
    pair <- strsplit(dir_name, "-")[[1]]
    source_group <- pair[1]
    sink_group <- pair[2]
    
    # 初始化migration_rate和R_square
    migration_rate <- ""
    R_square <- ""
    
    # 检查residence_stats.tsv文件是否存在
    tsv_file1 <- file.path(dir, "residence_stats.tsv")
    if (file.exists(tsv_file1)) {
      # 读取tsv文件
      tsv_data1 <- tryCatch(
        {
          read.delim(tsv_file1, header = TRUE, stringsAsFactors = FALSE)
        },
        error = function(e) NULL
      )
      
      # 如果成功读取且包含所需列
      if (!is.null(tsv_data1)) {
        if ("m" %in% colnames(tsv_data1)) {
          migration_rate <- as.character(tsv_data1$m[1])  # 取第一行的值
        }
        if ("Rsqr" %in% colnames(tsv_data1)) {
          R_square <- as.character(tsv_data1$Rsqr[1])  # 取第一行的值
        }
      }
    }
    
    # 初始化Above, Neutral, Below计数
    Above <- "0"
    Neutral <- "0"
    Below <- "0"
    Pa_Above <- "0"
    Pa_Neutral <- "0"
    Pa_Below <- "0"
    # 检查residence_result.tsv文件是否存在
    tsv_file2 <- file.path(dir, "residence_result.tsv")
    if (file.exists(tsv_file2)) {
      # 读取tsv文件
      tsv_data2 <- tryCatch(
        {
          read.delim(tsv_file2, header = TRUE, stringsAsFactors = FALSE)
        },
        error = function(e) NULL
      )
    tsv_data_pa<- tsv_data2 %>%filter(isPathogen == TRUE)
    
      
      # 如果成功读取且包含Partition列
      if (!is.null(tsv_data2) && "Partition" %in% colnames(tsv_data2)) {
        # 统计各分类的出现次数
        counts <- table(tsv_data2$Partition)

        if ("Above" %in% names(counts)) Above <- as.character(counts["Above"])
        if ("Neutral" %in% names(counts)) Neutral <- as.character(counts["Neutral"])
        if ("Below" %in% names(counts)) Below <- as.character(counts["Below"])
      }
      if (!is.null(tsv_data_pa) && "Partition" %in% colnames(tsv_data_pa)) {
      counts_pa<-table(tsv_data_pa$Partition)
      if ("Above" %in% names(counts_pa)) Pa_Above <- as.character(counts_pa["Above"])
      if ("Neutral" %in% names(counts_pa)) Pa_Neutral <- as.character(counts_pa["Neutral"])
      if ("Below" %in% names(counts_pa)) Pa_Below <- as.character(counts_pa["Below"])
      
    }
    }
    
    

    
    # 添加到结果数据框
    result_df <- rbind(result_df, data.frame(
      source = source_group,
      sink = sink_group,
      migration_rate = migration_rate,
      R_square = R_square,
      Above = Above,
      Neutral = Neutral,
      Below = Below,
      Pa_Above=Pa_Above,
      Pa_Neutral=Pa_Neutral,
      Pa_Below=Pa_Below,
      stringsAsFactors = FALSE
    ))
  }
  
  # 检查输出文件是否存在
  if (file.exists(output_file)) {
    # 如果文件存在，加载现有工作簿
    wb <- loadWorkbook(output_file)
    
    # 检查工作表是否存在，存在则删除
    if (sheet_name %in% names(wb)) {
      removeWorksheet(wb, sheet_name)
    }
    
    # 添加新工作表
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = result_df, startCol = 1, startRow = 1)
    
  } else {
    # 如果文件不存在，创建新工作簿
    wb <- createWorkbook()
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = result_df, startCol = 1, startRow = 1)
    
    # 确保输出目录存在
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  # 保存工作簿
  saveWorkbook(wb, output_file, overwrite = TRUE)

  # 绘制饼图并保存为PDF
  pie_chart_pdf <- file.path(dirname(output_file), 
                             paste0(srcType,"_",sheet_name, "_饼图.pdf"))
  draw_pie_charts(result_df, pie_chart_pdf)
  
}


##############################################################>
# 主程序1：调用上面函数进行处理 ----
# 安全关闭所有图形设备（即使当前没有打开的图形设备也不会报错）
graphics.off()
#开始处理
#############################################################>
# >要排除的目录名列表（完整名称）------
exclude_dirs <- c("elevator-Outdoor", "Outdoor-elevator",
                  "upper-Outdoor", "under-Outdoor",
                  "不需要的目录")
############################################################>

outDir<-paste0("result/net_bacteria2fungi/",getSampleArg$samplesName,"/sloan/")
#for (srcType in c("all","pathogen")) {
srcType<-"all"
  output_excel_file <- paste0(outDir,"migration summarize ",srcType,".xlsx")
  
  # 处理细菌数据
  bacteria_input_dir <- paste0("result/bacteria/20241019/",getSampleArg$samplesName,
                               "/02Plot/06sloan/",srcType)
  process_directories(srcType,exclude_dirs,bacteria_input_dir, output_excel_file, "bacteria")
  # 处理真菌数据
  fungi_input_dir <- paste0("result/fungi/20241020/",getSampleArg$samplesName,
                            "/02Plot/06sloan/",srcType)
  
  fungi_output_file <- "result/net_bacteria2fungi/15 group_by_floor/sloan/migration summarize.xlsx"
  process_directories(srcType,exclude_dirs,fungi_input_dir, output_excel_file, "fungi")
#}


#######################################################>
#Part3  画饼图2补充病原菌数量 ----
draw_pie_charts_add_pathogenCnt <- function(result_df, output_pdf) {
  # 确保输出目录存在
  output_dir <- dirname(output_pdf)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # 开启PDF设备（调整PDF尺寸适应小饼图）
  pdf(output_pdf, width = 10, height = ceiling(nrow(result_df)/5) )  # 增加每行显示数量
  
  # 创建绘图列表
  plot_list <- list()
  
  # 为每一行数据创建饼图
  for (i in 1:nrow(result_df)) {
    # 提取数据
    row_data <- result_df[i, ]
    group_name <- paste(row_data$source, row_data$sink, sep = "-")
    
    # 准备饼图数据并设置因子水平保证顺序
    pie_data <- data.frame(
      category = factor(c("Above", "Neutral", "Below"), 
                        levels = c("Above", "Neutral", "Below")),
      value = as.numeric(c(row_data$Above, row_data$Neutral, row_data$Below)),
      pathogenvalue = as.numeric(c(row_data$Pa_Above, row_data$Pa_Neutral, 
                               row_data$Pa_Below)),
      color = c("#fbb4ae", "#b3cde3", "#ccebc5")  # 绿色、蓝色、红色
    )
    
    # 计算百分比
    total <- sum(pie_data$value)
    if (total > 0) {
      pie_data$percent <- paste0(round(pie_data$value/total * 100, 0), "%",
                                "\n","(",pie_data$pathogenvalue,")")
    } else {
      pie_data$percent <- "0%"
    }
    
    # 创建实心饼图（无空心）
    p <- ggplot(pie_data, aes(x = 1, y = value, fill = category)) +
      geom_col(color = NA, width = 1) +  # 设置width=1确保完全填充
      coord_polar("y", start = 0) +
      scale_fill_manual(values = setNames(pie_data$color, levels(pie_data$category))) +
      geom_text(
        aes(label = ifelse(value > 0, percent, "")),  # 只显示非零值标签
        position = position_stack(vjust = 0.5),
        size = 1.5,  # 缩小标签字号
        color = "black"  # 白色标签提高对比度
      ) +
      labs(title = group_name) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 8),  # 缩小标题
        legend.position = "right",
        legend.title = element_text(size = 8),  # 图例标题大小
        legend.text = element_text(size = 6),  # 图例文字大小
        plot.margin = unit(c(0, 0, 0, 0), "cm"),  # 完全去除边距
        panel.spacing = unit(0, "cm")  # 去除面板间距
      )
    
    plot_list[[i]] <- p
  }
  
  # 每行显示4个饼图（可根据需要调整ncol值）
  grid.arrange(
    grobs = plot_list, 
    ncol = 4,  # 增加每行饼图数量
    padding = unit(0, "cm")  # 去除网格间距
  )
  
  # 关闭PDF设备
  dev.off()
}




###################################################>
#   主程序2：饼图上补充病原菌数量 ----
###################################################>

#需求1、在outDir变量指定的目录下，读取excel文件migration summarize all.xlsx和migration summarize pathogen.xlsx到数据框，每个excel有两个子表（bacteria和fungi）
#需求2、读入的将migration summarize pathogen.xlsx对应子表进行归并
library(readxl)
# 1. 设置 outDir 路径（替换为你的实际路径）
#outDir <- "your_directory_path_here/"

# 2. 读取所有 Excel 子表
# bacteria
all_bacteria <- read_excel(paste0(outDir, "migration summarize all.xlsx"), sheet = "bacteria")
#pathogen_bacteria <- read_excel(paste0(outDir, "migration summarize pathogen.xlsx"), sheet = "bacteria")

# fungi
all_fungi <- read_excel(paste0(outDir, "migration summarize all.xlsx"), sheet = "fungi")
#pathogen_fungi <- read_excel(paste0(outDir, "migration summarize pathogen.xlsx"), sheet = "fungi")

# 3. 处理 bacteria：添加列名后缀并按列合并
#pathogen_bacteria_renamed <- pathogen_bacteria %>% 
#  rename_with(~ paste0(., "_pathogen"))
#merged_bacteria <- bind_cols(all_bacteria, pathogen_bacteria_renamed)

# 4. 处理 fungi：完全相同的逻辑
#pathogen_fungi_renamed <- pathogen_fungi %>% 
#  rename_with(~ paste0(., "_pathogen"))
#merged_fungi <- bind_cols(all_fungi, pathogen_fungi_renamed)

# 5. 查看结果
#cat("Merged bacteria columns:\n")
#print(colnames(merged_bacteria))

#cat("\nMerged fungi columns:\n")
#print(colnames(merged_fungi))




# 绘制饼图并保存为PDF
pie_chart_pdf <- file.path(outDir, 
                      paste0("merged_bacteria" , "_饼图_all上增加病原菌数量_new.pdf"))
#draw_pie_charts_add_pathogenCnt(merged_bacteria, pie_chart_pdf)
draw_pie_charts_add_pathogenCnt(all_bacteria, pie_chart_pdf)  
pie_chart_pdf <- file.path(outDir, 
                           paste0("merged_fungi" , "_饼图_all上增加病原菌数量.pdf"))
draw_pie_charts_add_pathogenCnt(all_fungi, pie_chart_pdf)
#draw_pie_charts_add_pathogenCnt(merged_fungi, pie_chart_pdf)

