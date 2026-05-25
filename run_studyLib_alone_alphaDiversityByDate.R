library(dplyr)
library(ggplot2)
library(lubridate)
library(lme4)
library(grid)  
ZSR_LIB_PROG="run_studyLib20250706.R"
source(ZSR_LIB_PROG)          #加载自定义函数库
########  01数据准备
meta<-meco$sample_table
# 使用 mutate 添加将alpha多样性添加到meta文件
meta <- meta %>% mutate(alpha_diversity_Observed = meco$alpha_diversity$Observed)
meta <- meta %>% mutate(alpha_diversity_Chao1 = meco$alpha_diversity$Chao1)
meta <- meta %>% mutate(alpha_diversity_Shannon = meco$alpha_diversity$Shannon)

##############################20240105新增代码
outDir<-paste0(getSampleArg$outDir,getSampleArg$samplesName,"/02Plot/10alpha+LMM/")
dir.create(outDir,recursive = TRUE,showWarnings = FALSE)

fileName<-"alpha_diversity"
save_df_as_tsv(meta ,paste0(outDir,fileName,".tsv") )
# 对 alpha_diversity 进行log变换
meta$log_alpha_diversity <- log(meta$alpha_diversity_Observed )

# 温度转换成数字
meta$Temperature2  <- as.integer(meta$Temperature)
# 湿度转换成数字
meta$Humidness2  <- as.integer(sub("%$", "", meta$Humidness))
# 采样点转换成数字
meta$SampleSiteNum  <- as.integer(substr(meta$SampleSite, 1, nchar(meta$SampleSite) - 1))
# 计算每个样本非零 OTU 数量
# 对 OTU 表中的每列（每个样本）应用函数，统计非零值的数量
non_zero_counts <- colSums(meco$otu_table != 0) 
# 将非零 OTU 数量添加到 meta 数据框的 asvCounts 列
meta$asvCounts <- non_zero_counts[match(meta$SampleName, names(non_zero_counts))]

###Date为字符串列，格式为yymmdd，计算从230901开始的天数，放进列days
# 转换 Date 列为日期格式
meta$Date2 <- as.Date(meta$Date, format = "%y%m%d")
# 基准日期 "230901" 对应的日期
base_date <- as.Date("2023-09-01")
# 计算每个日期与基准日期的天数差异
meta$Days <- as.integer(meta$Date2 - base_date)

meta_filtered <- meta[meta$Days < 120, ]
meta_filtered2 <- meta[meta$Days > 120, ]

# 将 Date 列从字符型转换为日期格式
meta$Date_full <- ymd(paste0("20", substr(meta$Date, 1, 2), "-",
                                     substr(meta$Date, 3, 4), "-",
                                     substr(meta$Date, 5, 6)))
meta$Date_full <-meta$Date
#查看结果
#View(meta)
meta_filtered$Date_full <- ymd(paste0("20", substr(meta_filtered$Date, 1, 2), "-",
                             substr(meta_filtered$Date, 3, 4), "-",
                             substr(meta_filtered$Date, 5, 6)))
meta_filtered2$Date_full <- ymd(paste0("20", substr(meta_filtered2$Date, 1, 2), "-",
                                      substr(meta_filtered2$Date, 3, 4), "-",
                                      substr(meta_filtered2$Date, 5, 6)))
###########################################################
# start part1  
# 1.1  画所有分组的曲线
pdf(file=paste0(outDir,"alphaDiversityByDate.pdf"), width = 8, height = 6, pointsize = 300/72)
# p<-ggplot(meta, aes(x = Date_full, y = log_alpha_diversity)) +
#   geom_jitter(width = 0.01, height = 0.01, color = "black", size = 1) +  # 添加抖动效果
#   geom_smooth(aes(group = 1),method = "loess", color = "black", size = 0.5, se = TRUE) +  # 添加置信区间（se = TRUE）
#   labs(x = "", y = "Alpha Diversity(log)", title =
#          paste0("Alpha Diversity over Time with Prediction and Confidence Interval(all)")) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(size = 8, hjust = 0.5),  # 标题字体大小
#     axis.text.x = element_text(size = 12,angle = 90),  # X轴刻度标签
#     legend.text = element_text(size = 16),          # 图例文字大小（可取消注释斜体）
#     legend.title = element_text(size = 16, face = "bold"),  # 图例标题加粗
#     axis.text.y = element_text(size = 12,  hjust = 1),          # Y轴刻度标签
#     axis.title.x = element_text(size = 16, face = "bold"),  # X轴标题加粗
#     axis.title.y = element_text(size = 16, face = "bold"),   # Y轴标题加粗
#     plot.margin = unit(c(10, 100, 10, 10), "pt"),  # 右侧预留 100pt 空白
#     # 增加坐标轴线设置
#     #axis.line = element_line(color = "black", size = 0.5),  # 坐标轴线
#     #panel.grid.major = element_line(color = "grey90", size = 0.2),  # 主网格线
#     #panel.grid.minor = element_blank()  # 移除次网格线
#   )
p <- ggplot(meta, aes(x = Date_full, y = log_alpha_diversity)) +
  # 通过 aes() 映射颜色到分组变量（这里用固定值 "all" 作为分组）
  geom_jitter(aes(color = "all"), width = 0.01, height = 0.01, size = 1) +  
  geom_smooth(aes(group = 1), method = "loess", color = "black", size = 0.5, se = TRUE) +
  labs(
    x = "", 
    y = "Alpha Diversity(log)", 
    title = "Alpha Diversity over Time with Prediction and Confidence Interval(all)",
    color = "Group"  # 设置图例标题
  ) +
  # 手动设置颜色和图例标签
  scale_color_manual(
    values = c("all" = "black"),  # 点的颜色设为黑色
    labels = c("all")             # 图例标签设为 "all"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 8, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 90),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 12, hjust = 1),
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    plot.margin = unit(c(10, 30, 10, 10), "pt")
  )
print(p)
# 1.2  画单个分组的曲线
split_group<-split(meta,meta$Group)  #按Group字段拆分成多个list
for (i_group in seq_along(split_group)) {
  group=split_group[[i_group]]
  group$Date_full <- as.factor(group$Date_full) #转成因子
  p<-ggplot(group, aes(x = Date_full, y = log_alpha_diversity)) +
    geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 1) +  # 添加抖动效果
    geom_smooth(aes(group = 1),method = "loess", color = "black", size = 0.5, se = TRUE) +  # 添加置信区间（se = TRUE）
    labs(x = "", y = "Alpha Diversity(log)", title = 
           paste0("Alpha Diversity over Time with Prediction and Confidence Interval(",names(split_group)[i_group],")")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 8, hjust = 0.5),  # 标题字体大小
      axis.text.x = element_text(size = 12,angle = 90),  # X轴刻度标签
      legend.text = element_text(size = 16),          # 图例文字大小（可取消注释斜体）
      legend.title = element_text(size = 16, face = "bold"),  # 图例标题加粗
      axis.text.y = element_text(size = 12),          # Y轴刻度标签
      axis.title.x = element_text(size = 16, face = "bold"),  # X轴标题加粗
      axis.title.y = element_text(size = 16, face = "bold"),   # Y轴标题加粗
      # 增加坐标轴线设置
      #axis.line = element_line(color = "black", size = 0.5),  # 坐标轴线
      #panel.grid.major = element_line(color = "grey90", size = 0.2),  # 主网格线
      #panel.grid.minor = element_blank()  # 移除次网格线
    )
  
  print(p)
}
dev.off()
# part1  end
#########################################################

#########################################################
# start part2  划各组曲线
################ 02建立LMM模型
#https://zhuanlan.zhihu.com/p/649029208 R使用线性混合效应模型一些常见的结果指标
#https://zhuanlan.zhihu.com/p/461859366 在R中使用线性混合效应模型的实现Linear Mixed-Effects Modeling
#https://www.jindouyun.cn/document/industry/details/184541 R统计绘图-线性混合效应模型详解
## linear mixed models - 建立公式
#2/365=0.005479
# re<-"log_alpha_diversity ~ sin(0.005479*pi*Days)+cos(0.005479*pi*Days)+"
# fe<-"(1|Temperature2)+(1|Humidness2)"
# formula1<-paste0(re,fe)
# l1<-lmer(formula = formula1,data = meta)
formulas <- list(
  log_alpha_diversity ~ sin(0.005479*pi*Days)+cos(0.005479*pi*Days)+(1|Temperature2)+(1|Humidness2),
  log_alpha_diversity ~ sin(0.010958*pi*Days)+cos(0.010958*pi*Days)+(1|Temperature2)+(1|Humidness2)
  #(temp.mean+temp.mean2+lat:(sin(2*pi*yday/365)+cos(2*pi*yday/365)))
)
## linear mixed models - 根据公式建立模型及计算残差
#2.1 根据分组+公式生成模型
models=list()
split_group2=list()

for (i_formulas in c(1:length(formulas))){
  global_model<-lmer( formulas[[i_formulas]], meta)## linear mixed models - reference values from older code
  global_predictions <- predict(global_model, newdata = meta, re.form = NULL, se.fit = TRUE)
  global_predict<-data.frame(global_predictions$fit)
  global_predict$global_se<-global_predictions$se.fit
  global_predict$SampleName<-rownames(global_predict)

  
  for (i_group in seq_along(split_group)) {
    group=split_group[[i_group]]
    group <- merge(group, global_predict, by = "SampleName")
    print(paste0("根据公式生成模型:",i))
    model<-lmer( formulas[[i_formulas]], group )## linear mixed models - reference values from older code
    models[paste0(names(split_group)[i_group],"_formula",i_formulas)]<-model  #存储模型到list
    # 计算拟合值和残差，绘制残差 vs 拟合值图
    #plot(fitted(model), residuals(model),main= deparse(formula(model)) )
    #abline(h = 0, col = "red")
    
    ###生成画图数据
    #print(model)
    # 生成模型的预测值
    predictions <- predict(model, newdata = group, re.form = NULL, se.fit = TRUE)
    # 提取预测值和标准误差
    group$predicted <- predictions$fit
    group$se <- predictions$se.fit
    # 计算95%置信区间（假设正态分布）
    group$global_predictions<-group$global_predictions.fit
    group$lower <- group$predicted - 1.96 * group$se
    group$upper <- group$predicted + 1.96 * group$se
    group$global_lower <- group$predicted - 1.96 * group$se
    group$global_upper <- group$predicted + 1.96 * group$se
    # 绘制图形
    group$All_Samples<-"all_sample"
    #将meta文件保存
    #fileName<-paste0("Meta ",i)#,"-",gsub(".*~\\s*", "", deparse(formula(model))),"")
    #save_df_as_tsv(meta ,paste0(outDir,fileName,".tsv") )
    
    # 扩展第一级（确保长度足够）
    while (length(split_group2) < i_formulas) {
      split_group2 <- c(split_group2, list(NULL))
    }
    # 初始化第二级（如果尚未存在）
    if (is.null(split_group2[[i_formulas]])) {
      split_group2[[i_formulas]] <- list()
    }
    # 添加对象到指定位置
    split_group2[[i_formulas]][[i_group]] <- group
  }
  
  
  ###生成画图数据
  #print(model)
  # 生成模型的预测值
  
  
}



#2.2 模型计算
# 创建一个文本文件并写入 ANOVA 结果
sink(paste0(outDir,"models_anova_results.txt"))
for (i in seq_along(models)) {
    # 执行 ANOVA 并将结果打印到文本文件
    cat("\n============ Model",names(models)[i],"->",deparse(formula(models[[i]])),":\n")
    #查看模型准确性
    print(anova(models[[i]]))
    cat("\n")
}
# 关闭文件连接
sink()

#####  03 画alpha曲线及增加拟合曲线
pdf(file=paste0(outDir,"alphaDiversityByDate2.pdf"), width = 8, height = 6, pointsize = 300/72)
tryCatch({
  for(i in seq_along(split_group2) ) {
    split_group_i <- do.call(rbind, split_group2[[i]]) #取第i个
    split_group_i$Date_full <- as.factor(split_group_i$Date_full) #转成因子
    p <- ggplot(split_group_i, aes(x = Date_full, y = log_alpha_diversity, color = Group, group = Group)) +
      geom_jitter(width = 0.01, height = 0.01, size = 1) +  # 添加抖动效果
      geom_smooth(aes(y = predicted), method = "loess", size = 1, 
                  fill = "grey", alpha = 0.2, se = TRUE, level = 0.8) +  # 使用 loess 平滑并显示灰色置信区间
      geom_smooth(aes(y = global_predictions, group = 1),method = "loess", size = 1, color = "black", linetype = "solid") +
      labs(x = "", y = "Alpha Diversity(log)", 
           title = paste0("Alpha Diversity over Date [ ",
                          paste(unlist(deparse(formulas[[i]])),collapse = "\n")," ]") ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 8, hjust = 0.5),  # 标题字体大小
        axis.text.x = element_text(size = 12,angle = 90),  # X轴刻度标签
        legend.text = element_text(size = 16),          # 图例文字大小（可取消注释斜体）
        legend.title = element_text(size = 16, face = "bold"),  # 图例标题加粗
        axis.text.y = element_text(size = 14),          # Y轴刻度标签
        axis.title.x = element_text(size = 16, face = "bold"),  # X轴标题加粗
        axis.title.y = element_text(size = 16, face = "bold")   # Y轴标题加粗
        )
    p<-p+ scale_color_manual( #统一分组颜色
      values = get_sci_palette(getSampleArg$sampleFilterGroupName)
    )
    
    
    print(p)
    
  }
}, error = function(e) {
  message("出错: ", e$message)
}, finally = {
  dev.off()
})

  # fileName<-paste0("Meta ",i)#,"-",gsub(".*~\\s*", "", deparse(formula(model))),"")
  # #fileName<-substr(gsub("\\*","_m_",gsub("\\|", "｜", fileName)),1,55)
  # #print(paste0(outDir,fileName) )
  # save_df_as_tsv(split_group ,paste0(outDir,fileName,".tsv") )
  # #############################start 原先代码
  # # 绘制散点图
  # # 使用 geom_jitter() 绘制抖动点图，并添加 Loess 曲线（非线性拟合）以及95%置信区间
  # p<-ggplot(meta, aes(x = Date_full, y = log_alpha_diversity)) +
  #   geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 1) +  # 添加抖动效果
  #   geom_smooth(method = "loess", color = "red", size = 1, se = TRUE) +  # 添加置信区间（se = TRUE）
  #   labs(x = "Date", y = "Alpha Diversity(log)", title = "Alpha Diversity over Time with Prediction and Confidence Interval") +
  #   theme_minimal() 
  # print(p)

#############################End 原先代码
stop("end")
# part1  end
###############################################################


for (i in c(1:length(formulas))){
  print(i)
  models[i] <- lmer( formulas[[i]], meta )## linear mixed models - reference values from older code
  #print(models[[i]])
  # 计算拟合值和残差，绘制残差 vs 拟合值图
  plot(fitted(models[[i]]), residuals(models[[i]]),main= deparse(formula(models[[i]])) )
  abline(h = 0, col = "red")
}

for (i in c(1:length(formulas))){
  print(i)
  models2[i] <- lmer( formulas[[i]], meta_filtered )## linear mixed models - reference values from older code
  #print(models[[i]])
  # 计算拟合值和残差，绘制残差 vs 拟合值图
  plot(fitted(models2[[i]]), residuals(models2[[i]]),main= deparse(formula(models2[[i]])) )
  abline(h = 0, col = "red")
}

for (i in c(1:length(formulas))){
  print(i)
  models3[i] <- lmer( formulas[[i]], meta_filtered2 )## linear mixed models - reference values from older code
  #print(models[[i]])
  # 计算拟合值和残差，绘制残差 vs 拟合值图
  plot(fitted(models3[[i]]), residuals(models3[[i]]),main= deparse(formula(models3[[i]])) )
  abline(h = 0, col = "red")
}

dev.off()

# 创建一个文本文件并写入 ANOVA 结果
sink(paste0(outDir,"anova_results.txt"))
for (i in c(1:length(formulas)) ) {
  # 执行 ANOVA 并将结果打印到文本文件
  cat("\n============ Model",i,"->",deparse(formula(models[[i]])),":\n")
  #查看模型准确性
  print(anova(models[[i]]))
  cat("\n")
}
for (i in c(1:length(formulas)) ) {
  # 执行 ANOVA 并将结果打印到文本文件
  cat("\n============ Model",i,"->",deparse(formula(models2[[i]])),":\n")
  #查看模型准确性
  print(anova(models2[[i]]))
  cat("\n")
}


#所有模型anova比较
cat("\n============ all",":\n")
print( do.call(anova, models) )
# 关闭文件连接
sink()

#####  03 画alpha曲线及增加拟合曲线
pdf(file=paste0(outDir,"alphaDiversityByDate2.pdf"), width = 16, height = 8, pointsize = 300/72)
i=1
for ( model in  models){
  #print(model)
  # 生成模型的预测值
  predictions <- predict(model, newdata = meta, re.form = NULL, se.fit = TRUE)
  # 提取预测值和标准误差
  meta$predicted <- predictions$fit
  meta$se <- predictions$se.fit
  # 计算95%置信区间（假设正态分布）
  meta$lower <- meta$predicted - 1.96 * meta$se
  meta$upper <- meta$predicted + 1.96 * meta$se
  # 绘制图形
  p<-ggplot(meta, aes(x = Date_full, y = log_alpha_diversity)) +
    geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 1) +  # 添加抖动效果
    #geom_line(aes(y = predicted), color = "red", size = 1) +  # 绘制拟合曲线
    #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +  # 添加95%置信区间
    geom_smooth(aes(y = predicted), method = "loess", color = "red", size = 1, 
                fill = "grey", alpha = 0.2, se = TRUE, level=0.8) +  # 使用 loess 平滑并显示灰色置信区间
    labs(x = "Date", y = "Alpha Diversity(log)", 
         title = paste0("Alpha Diversity over Date [ ", deparse(formula(model))," ]") ) +
    theme_minimal()
  print(p)
  #将meta文件保存
  fileName<-paste0("Meta ",i)#,"-",gsub(".*~\\s*", "", deparse(formula(model))),"")
  #fileName<-substr(gsub("\\*","_m_",gsub("\\|", "｜", fileName)),1,55)
  
  #print(paste0(outDir,fileName) )
  save_df_as_tsv(meta ,paste0(outDir,fileName,".tsv") )
  i<-i+1
}

dev.off()


pdf(file=paste0(outDir,"filtered1_alphaDiversityByDate2.pdf"), width = 16, height = 8, pointsize = 300/72)
i=1
for ( model in  models2){
  #print(model)
  # 生成模型的预测值
  predictions <- predict(model, newdata = meta_filtered, re.form = NULL, se.fit = TRUE)
  # 提取预测值和标准误差
  meta_filtered$predicted <- predictions$fit
  meta_filtered$se <- predictions$se.fit
  # 计算95%置信区间（假设正态分布）
  meta_filtered$lower <- meta_filtered$predicted - 1.96 * meta_filtered$se
  meta_filtered$upper <- meta_filtered$predicted + 1.96 * meta_filtered$se
  # 绘制图形
  p<-ggplot(meta_filtered, aes(x = Date_full, y = log_alpha_diversity)) +
    geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 1) +  # 添加抖动效果
    #geom_line(aes(y = predicted), color = "red", size = 1) +  # 绘制拟合曲线
    #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +  # 添加95%置信区间
    geom_smooth(aes(y = predicted), method = "loess", color = "red", size = 1, 
                fill = "grey", alpha = 0.2, se = TRUE, level=0.8) +  # 使用 loess 平滑并显示灰色置信区间
    labs(x = "Date", y = "Alpha Diversity(log)", 
         title = paste0("Alpha Diversity over Date [ ", deparse(formula(model))," ]") ) +
    theme_minimal()
  print(p)
  #将meta文件保存
  fileName<-paste0("Meta ",i)#,"-",gsub(".*~\\s*", "", deparse(formula(model))),"")
  #fileName<-substr(gsub("\\*","_m_",gsub("\\|", "｜", fileName)),1,55)
  
  #print(paste0(outDir,fileName) )
  save_df_as_tsv(meta ,paste0(outDir,fileName,".tsv") )
  i<-i+1
}
dev.off()
pdf(file=paste0(outDir,"filtered2_alphaDiversityByDate2.pdf"), width = 16, height = 8, pointsize = 300/72)
i=1
for ( model in  models3){
  #print(model)
  # 生成模型的预测值
  predictions <- predict(model, newdata = meta_filtered2, re.form = NULL, se.fit = TRUE)
  # 提取预测值和标准误差
  meta_filtered2$predicted <- predictions$fit
  meta_filtered2$se <- predictions$se.fit
  # 计算95%置信区间（假设正态分布）
  meta_filtered2$lower <- meta_filtered2$predicted - 1.96 * meta_filtered$se
  meta_filtered2$upper <- meta_filtered2$predicted + 1.96 * meta_filtered$se
  # 绘制图形
  p<-ggplot(meta_filtered2, aes(x = Date_full, y = log_alpha_diversity)) +
    geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 1) +  # 添加抖动效果
    #geom_line(aes(y = predicted), color = "red", size = 1) +  # 绘制拟合曲线
    #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "red", alpha = 0.2) +  # 添加95%置信区间
    geom_smooth(aes(y = predicted), method = "loess", color = "red", size = 1, 
                fill = "grey", alpha = 0.2, se = TRUE, level=0.8) +  # 使用 loess 平滑并显示灰色置信区间
    labs(x = "Date", y = "Alpha Diversity(log)", 
         title = paste0("Alpha Diversity over Date [ ", deparse(formula(model))," ]") ) +
    theme_minimal()
  print(p)
  #将meta文件保存
  fileName<-paste0("Meta ",i)#,"-",gsub(".*~\\s*", "", deparse(formula(model))),"")
  #fileName<-substr(gsub("\\*","_m_",gsub("\\|", "｜", fileName)),1,55)
  
  #print(paste0(outDir,fileName) )
  save_df_as_tsv(meta ,paste0(outDir,fileName,".tsv") )
  i<-i+1
}
dev.off()
# ggplot(meco$sample_table, aes(x = Date_full, y = log_alpha_diversity ,fill=Group)) +
#   geom_jitter(width = 0.01, height = 0.01, color = "blue", size = 1) +  # 添加抖动效果
#   geom_smooth(method = "loess", color = "red", size = 1, se = TRUE) +  # 添加置信区间（se = TRUE）
#   labs(x = "Date", y = "Alpha Diversity(log)", title = "Alpha Diversity over Time with Prediction and Confidence Interval") +
#   theme_minimal() 
# 

