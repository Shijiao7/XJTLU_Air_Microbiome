
library(tidygraph)
library(igraph)
library(readxl)
library(patchwork)
library(purrr)
library(dplyr)
#install.packages("ggrepel")
library("ggrepel")
library(ggraph)

# 安全关闭所有图形设备（即使当前没有打开的图形设备也不会报错）
graphics.off()

#加载showtext包，没有就先安装
if (!require(showtext)) install.packages("showtext")
library(showtext)
# 自动加载系统字体
showtext_auto()

# 设置输出目录
outDir <- paste0("result/net_bacteria2fungi/", getSampleArg$samplesName, "/sloan/")
print(paste0("输出目录:",outDir))
for (srcType in (c("all","pathogen"))){
  fileName <- paste0(outDir, "migration summarize ",srcType,".xlsx")
  print(paste("处理输入文件：",fileName))
  # 获取所有工作表名称
  sheet_names <- excel_sheets(fileName)
  create_network_plot <- function(sheet_name) {
    tryCatch({
      # 在读取数据后立即过滤(过滤掉 NA 边​)
      diffusion_data <- read_excel(fileName, sheet = sheet_name) %>% 
        dplyr:: select(1:4) %>%  
        setNames(c("source", "sink", "migration_rate", "R_square")) %>% 
        mutate(across(3:4, as.numeric)) %>%  # Convert cols 3-4 to numeric
        filter(!is.na(migration_rate)) %>%  # 关键：移除 migration_rate 为 NA 的行
        filter(!if_all(everything(), is.na)) %>% 
        filter(!is.na(source), !is.na(sink))
      
      # 检查是否有有效数据
      if (nrow(diffusion_data) == 0) {
        message(paste("工作表", sheet_name, "无有效数据"))
        return(NULL)
      }
      
      # 计算最大R²值
      max_r_square <- max(diffusion_data$R_square, na.rm = TRUE)
      min_r_square <- min(diffusion_data$R_square, na.rm = TRUE)
      # 转换为图形对象
      graph <- as_tbl_graph(diffusion_data, directed = TRUE)
      
      # 计算节点大小
      max_label_length <- max(nchar(c(diffusion_data$source, diffusion_data$sink)))
      fixed_node_size <- 6 + max_label_length * 0.8  # 调整缩放系数
      create_fixed_layout_matrix <- function(graph) {
        node_names <- igraph::V(graph)$name
        n_nodes <- length(node_names)
        
        # 创建基础布局矩阵
        layout_matrix <- matrix(runif(n_nodes * 2, -1, 1), ncol = 2)
        
        # 定义关键节点固定位置
        key_positions <- list(
          "Underground" = c(0, -1),    # 下方
          "Aboveground" = c(0, 1),      # 上方
          "Outdoor" = c(-1, 0),         # 左方
          "Elevator" = c(1, 0)          # 右方
        )
        
        # 应用关键节点位置
        for (node_name in names(key_positions)) {
          if (node_name %in% node_names) {
            node_index <- which(node_names == node_name)
            layout_matrix[node_index, ] <- key_positions[[node_name]]
          }
        }
        
        return(layout_matrix)
      }
      
      # 创建图形 - 使用矩阵格式
      set.seed(123)
      layout_matrix <- create_fixed_layout_matrix(graph)
      
      # 正确的使用方法：直接提供矩阵，不指定参数名
      p <- ggraph(graph, layout = layout_matrix) +
        geom_edge_fan(
          aes(width = migration_rate, color = R_square,
              label = sprintf("%.4f", migration_rate),label_pos = 0.4,angle_calc = 'along'),
          alpha = 0.9,
          spread = 1.5,  # 控制多条连线间的间距（值越大间距越大）
          arrow = arrow(type = "closed", length = unit(3, "mm")),
          end_cap = circle(6, "mm"),start_cap = circle(6, "mm")
        ) +
        geom_node_point(
          size = fixed_node_size,color = "lightblue",alpha = 0.8
        ) +
        geom_node_text(
          aes(label = name),size = 2.5,color = "black",vjust = 0.5,hjust = 0.5
        ) +
        scale_edge_width_continuous(
          name = "Migration Rate",range = c(0.3, 2),breaks = scales::pretty_breaks(5),
          guide = "none"  # 关键：隐藏图例
        ) +
        scale_edge_color_gradient(
          name = expression(R^2),low = "#FF6B6B",high = "#4ECDC4",
          limits = c(0, 1)  # 颜色标尺
        ) +
        labs(title = sheet_name) +
        theme_graph(base_size = 8) +
        theme(
          legend.position = "right",
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          legend.key.size = unit(0.3, "cm")
        )
      return(p)
    }, error = function(e) {
      message(paste("处理工作表", sheet_name, "时出错:", e$message))
      return(NULL)
    })
  }
  
  # 安全处理所有工作表
  plot_list <- list()
  for (sheet in sheet_names) {
    plot <- create_network_plot(sheet)
    if (!is.null(plot)) {
      plot_list[[sheet]] <- plot
    }
  }
  if (length(plot_list) == 0) {
    stop("没有可用的有效数据来创建网络图")
  }
  
  # 计算排列方式（每行最多3列）
  n_plots <- length(plot_list)
  n_cols <- min(3, n_plots)
  n_rows <- ceiling(n_plots / n_cols)
  
  # 合并图形
  combined_plot <- wrap_plots(plot_list, ncol = n_cols) +
    plot_annotation(
      title = paste0( "Microbial Migration Networks (",srcType,")" ),
      theme = theme(legend.position = "right",  # 将图例放在右侧
        plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
    )
  
  # 动态调整输出尺寸
  base_width <- 8 # 单图基础宽度
  base_height <- 6   # 单图基础高度
  output_width <- base_width * n_cols + 1  # 额外空间给标题
  output_height <- base_height * n_rows + 0.5
  
  # 保存合并图形
  # 根据系统选择是否设置family参数
  pdf_name<-paste0(outDir,srcType, "_sloan_combined_networks.pdf")
  if(Sys.info()["sysname"] == "Windows") {
    cairo_pdf(pdf_name,width = output_width,height = output_height,family="Microsoft YaHei")
  } else {
    pdf(pdf_name, width = output_width, height = output_height, family = "Courier")
  }
  print(combined_plot)
  dev.off()
  
  # 单独保存每个图形（仅保存成功的）
  walk2(plot_list, sheet_names[seq_along(plot_list)], ~ {
    pdf_name<-paste0(outDir, srcType , "_sloan_network_", .y, ".pdf")
    if(Sys.info()["sysname"] == "Windows") {
      cairo_pdf(pdf_name,width = base_width,height = base_height,family="Microsoft YaHei")
    } else {
      pdf(pdf_name,width = base_width, height = base_height, family = "Courier")
    }
    # 获取当前标题
    current_title <- .x$labels$title
    #补充原先标题，增加Microbial Migration Networks
    .x<- .x + ggtitle(paste0("Microbial Migration Networks (", current_title,")")) 
    print( .x )
    dev.off()
  })
}  # end for

