#############################>
# Sloan模型研究
############################>



####################################>
# >中性模型neutral model from Sloan ------------------------------------------------
####################################>
#Adam Burns - 2/10/2015
#aburns2@uoregon.edu
#From Burns et al. Contribution of neutral processes to the assembly of the gut microbial communities changes over host development
#Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity when stats=FALSE. For use in R.
#spp: A community table for communities of interest with local communities/samples as rows and taxa as columns. All samples must be rarefied to the same depth.
#pool: A community table for defining source community (optional; Default=NULL).
#taxon: A table listing the taxonomic calls for each otu, with OTU ids as row names and taxonomic classifications as columns.
#If stats=TRUE the function will return fitting statistics.
#If stats=FALSE the function will return a table of observed and predicted values for each otu.
sncm.fit <- function(spp, pool=TRUE, stats=FALSE, taxon=NULL){
  require(minpack.lm)
  require(Hmisc)
  require(stats4)

  options(warn=-1)
  tryCatch({
    #Calculate the number of individuals per community
    N <- mean(apply(spp, 1, sum))
    
    
    #Calculate the average relative abundance of each taxa across communities
    if(is.null(pool)){
      p.m <- apply(spp, 2, mean)
      p.m <- p.m[p.m != 0]
      p <- p.m/N
    } else {
      p.m <- apply(pool, 2, mean)
      p.m <- p.m[p.m != 0]
      p <- p.m/N
    }
    
    #Calculate the occurrence frequency of each taxa across communities
    spp.bi <- 1*(spp>0)
    freq <- apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]
    
    #Combine
    C <- merge(p, freq, by=0)
    C <- C[order(C[,2]),]
    C <- as.data.frame(C)
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
    p <- C.0[,2]
    freq <- C.0[,3]
    names(p) <- C.0[,1]
    names(freq) <- C.0[,1]
    
    #Calculate the limit of detection
    d = 1/N
    
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
    m.ci <- confint(m.fit, 'm', level=0.95)
    
    ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
    sncm.LL <- function(m, sigma){
      R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
      R = dnorm(R, 0, sigma)
      -sum(log(R))
    }
    m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
    
    ##Calculate Akaike's Information Criterion (AIC)
    aic.fit <- AIC(m.mle, k=2)
    bic.fit <- BIC(m.mle)
    
    ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
    
    pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for binomial model
    bino.LL <- function(mu, sigma){
      R = freq - pbinom(d, N, p, lower.tail=FALSE)
      R = dnorm(R, mu, sigma)
      -sum(log(R))
    }
    bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    
    aic.bino <- AIC(bino.mle, k=2)
    bic.bino <- BIC(bino.mle)
    
    ##Goodness of fit for binomial model
    bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
    Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
    
    bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma){
      R = freq - ppois(d, N*p, lower.tail=FALSE)
      R = dnorm(R, mu, sigma)
      -sum(log(R))
    }
    pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    
    aic.pois <- AIC(pois.mle, k=2)
    bic.pois <- BIC(pois.mle)
    
    ##Goodness of fit for Poisson model
    pois.pred <- ppois(d, N*p, lower.tail=FALSE)
    Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
    
    pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Results
    if(stats==TRUE){
      fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
      fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
      return(fitstats)
    } else {
      A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
      A <- as.data.frame(A)
      colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
      if(is.null(taxon)){
        B <- A[order(A[,1]),]
      } else {
        B <- merge(A, taxon, by=0, all=TRUE)
        row.names(B) <- B[,1]
        B <- B[,-1]
        B <- B[order(B[,1]),]
      }
      return(B)
    }
    }
  ,error = function(e)  {
    # 打印错误信息
    message("函数sncm.fit执行错误，错误信息：", conditionMessage(e))
  })
}

#################################################################>
# 函数：f_Phyloseq_splitByGroup
# >将PS对象按Group拆成PS对象列表 --------------------------------------------------
# 将Phyloseq对象按照样本meta文件的Group字段拆分成Phyloseq对象列表
#################################################################>
f_Phyloseq_splitByGroup<-function(pys) {
  #pys<-ps_pathogen_fungi
  # 获取所有唯一的分组名
  all_groups <- unique(sample_data(pys)$Group)
  print(all_groups)
  # 根据每个分组生成不同的phyloseq对象
  i<-1
  list_of_physeq<-list()
  for (byGroup in all_groups) {
    #print(byGroup)
    #将byGroup放置为全局变量，因为subset_samples参数需要全局变量
    assign( "byGroup",byGroup, env=globalenv() )
    subset_physeq <- subset_samples(pys, Group== byGroup)
    list_of_physeq[i]<-subset_physeq
    i<-i+1
  }
  return(list_of_physeq)  # 返回生成的phyloseq对象列表
}
###############################################>
#f_mecoAbund_getTaxaName(taxa,taxaOrder)
# 参数1：taxa meco生成丰度的组合字符串，如k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__acnes
# 参数2：taxaOrder ，需要取参数1的第几个进行查找
# 返回：跟进参数2指定物种类别，将参数1对应物种类别返回，如没有返回NA   "unclassified"
##############################################>
f_mecoAbund_getTaxa<-function(taxa,taxaOrder){
  ############################>
  # 01 取taxa第taxaOrder名称
  # 使用 strsplit 函数按照|进行分割
  parts <- unlist(strsplit(taxa, "\\|"))
  # 取第taxaOrder个值
  taxa1 <- parts[taxaOrder]
  # 使用 strsplit 函数按照 __ 进行分割，并取出每个部分的最后一个值
  part2s <- unlist(strsplit(taxa1, "__"))
  # 检查是否存在 __ 分隔符
  if (length(part2s) > 1) {
    # 如果存在，取__分隔符后面的字符串
    return( tail(part2s, 1) )
  } else {
    return(NA) #没有注释到物种，返回不在
  }
  return(NA)
}

###############################################>
#f_isCommonTaxa(taxa,taxaOrder,common_taxas)
# 参数1：taxa meco生成丰度的组合字符串，如k__Bacteria|p__Actinobacteria|c__Actinobacteria|o__Actinomycetales|f__Propionibacteriaceae|g__Propionibacterium|s__acnes
# 参数2：taxaOrder ，需要取参数1的第几个进行查找
# 参数3：common_taxas，需要查找的物种
# 返回：如果在参数3找到，就返回找到的公共物种，如找不到返回"unclassified"
##############################################>
f_isInCommonTaxa<-function(taxa,taxaOrder,common_taxas){
  ############################>
  # 01 取taxa第taxaOrder名称
  # 使用 strsplit 函数按照|进行分割
  #print(taxa)
  parts <- unlist(strsplit(taxa, "\\|"))
  # 取第taxaOrder个值
  taxa1 <- parts[taxaOrder]
  # 使用 strsplit 函数按照 __ 进行分割，并取出每个部分的最后一个值
  part2s <- unlist(strsplit(taxa1, "__"))
  # 检查是否存在 __ 分隔符
  if (length(part2s) > 1) {
    # 如果存在，取__分隔符后面的字符串
    taxa2 <- tail(part2s, 1)
    # 使用 %in% 运算符检查 taxa2 是否在 common_taxas 中
    if (taxa2 %in% common_taxas) {
      return(taxa2)
    } else {
      return("unclassified")
    }
  } else {
    # 如果不存在，设置为 NA 或者其他你想要的默认值
    #taxa2 <- NA
    return("unclassified") #没有注释到物种，返回不在
  }
  return("unclassified")
}

#################################################################>
# 函数：f_microeco_otuSumByTaxa
# 将Phyloseq对象otu表按照taxa汇总
#################################################################>
f_microeco_otuSumByTaxa<-function(pys_src,pys_dst,taxaOrder,abund_dir){
  #pys_src=ps_src
  #pys_dst=ps_dst
  #abund_dir="result/sloan/meco_abund"
  #taxaOrder=7
  
  #################################>
  # 01 pys转换成microeco对象
  #################################>
  print("01 ps对象转meco对象...")
  meco_src <- phyloseq2meco(pys_src)
  meco_dst <- phyloseq2meco(pys_dst)
  
  #################################>
  # 02 meco对象处理成丰度...
  #################################>
  print("02 meco对象处理成丰度...")
  meco_src$cal_abund(rel = FALSE)
  meco_src$save_abund(dirpath = paste0(abund_dir,"_src"))
  meco_dst$cal_abund(rel = FALSE)
  meco_dst$save_abund(dirpath = paste0(abund_dir,"_dst"))

  #################################>
  #03 源和目的取共同的物种...
  #################################>
  print("03 源和目的取共同的物种...")
  df_src<-meco_src$taxa_abund[[taxaOrder]]
  df_dst<-meco_dst$taxa_abund[[taxaOrder]]
  # 获取数据框的行名（行索引）
  row_names <- rownames(df_src)
  # 使用 sapply 函数将每个行名作为参数传递给函数
  taxas_src <- sapply(row_names, function(row_name) f_mecoAbund_getTaxa(row_name,taxaOrder))
  # 获取数据框的行名（行索引）
  row_names <- rownames(df_dst)
  # 使用 sapply 函数将每个行名作为参数传递给函数 
  taxas_dst <- sapply(row_names, function(row_name) f_mecoAbund_getTaxa(row_name,taxaOrder))
  # 找出两个meco对象中共有的taxa信息
  common_taxas <- intersect(taxas_src, taxas_dst)
  head(common_taxas)
  #去掉NA
  library(purrr)
  common_taxas <- discard(common_taxas,is.na)
  print(paste0("共有一致的物种数：",length(common_taxas) ) )
  #print(" 查看《筛选共有类别的taxa》 ")
  #print(common_taxas)
  
  
  #################################>
  #04 源和目的按照共同的物种汇总...
  #################################>
  #初始化返回列表
  df_summary_ret<-list()
  for (i in 1:2) {
    if (i==1) {
      meco_tmp=meco_src
    } else {
      meco_tmp=meco_dst
    }
    #按照参数取物种分类
    df<-meco_tmp$taxa_abund[[taxaOrder]]
    # 获取数据框的行名（行索引）
    row_names <- rownames(df)
    # 提取竖线分隔的第taxaOrder个部分作为新的列,如在公共物种找到就返回物种，否则返回unclassified
    taxa_column <- sapply(row_names, function(row_name) f_isInCommonTaxa(row_name,taxaOrder,common_taxas))
    # 将新列插入到数据框的第一列位置
    df <- cbind(taxa_column = taxa_column, df)
    #按照第一列taxa_column进行汇总
    df_summary <- df %>%
      group_by(taxa_column) %>%
      summarise_all(sum) %>%
      ungroup()

    # 将第一列作为行名
    rownames(df_summary) <- df_summary$taxa_column
    #转置
    df_summary<-t(df_summary)
    #转换成数据框
    df_summary <-as.data.frame(df_summary)
    # 删除行名为 taxa_column 的行
    df_summary <- df_summary[!rownames(df_summary) %in% c("taxa_column"), ]
    #将汇总结果存储进列表
    df_summary_ret <- append(df_summary_ret, list(df_summary))
  }
  return(df_summary_ret)
}


#################################################################>
# 函数：f_gen_sloan
# >生成sloan扩散模型 ----
#################################################################>
f_gen_sloan<-function(sloanDir,srcType,f_microbeType,formatted_pair,spp,pool,
                      pathogens,dst_name,src_name) {
  tryCatch({
    #sloanDir<-paste0(sloanDir,"/")
    #通过参数传递
    #pool <- read_tsv(paste0(sloanDir,"src.tsv"))
    #spp <- read_tsv(paste0(sloanDir,"dst.tsv"))

    stats.otu <- sncm.fit(spp,pool,stats=FALSE)
    stats.otu$isPathogen<-rownames(stats.otu) %in% pathogens  #找到是病原菌的
    write.table(stats.otu, paste0(sloanDir,"residence.tsv"), sep="\t", col.names=NA)
    
    stats.all <- sncm.fit(spp,pool,stats=TRUE)
    write_tsv(stats.all, paste0(sloanDir,"residence_stats.tsv"))
    
    
    library(readr)
    library(dplyr)
    library(ggplot2)
    
    place <- read_tsv( paste0(sloanDir,"residence.tsv"))
    names(place)[1] <- c("Genus")
    place <- place %>% mutate(Partition=ifelse(freq < pred.lwr, "Below", "Neutral"))
    place <- place %>% mutate(Partition=ifelse(freq > pred.upr, "Above", Partition))
    write_tsv(place, paste0(sloanDir,"residence_Partition.tsv"))
    
    spp <- read_tsv(paste0(sloanDir,"dst.tsv"))
    #spp$unclassified <- 0
    N <- mean(apply(pool, 1, sum))
    p.m <- apply(pool, 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/N
    p <- as.data.frame(p)
    colnames(p) <- c("mean_abundance")
    ############>
    # 将行名存储到一个新列中，并命名为 "Genus"
    p$Genus <- rownames(p)
    # 移除默认的行名
    rownames(p) <- NULL
    # 调整列的顺序，确保 "Genus" 列在第一列
    p <- p[, c("Genus", names(p)[!names(p) %in% "Genus"])]
    ############>
    write.table(p, paste0(sloanDir,"residence_mean_abundance.tsv"), sep="\t", row.names=FALSE,col.names = TRUE)
    #Plese change the colnames in abundance manually
    #merge mean_abundance with sloean prediciton results
    abundance <- read_tsv(paste0(sloanDir,"residence_mean_abundance.tsv"))
    names(abundance)[1] <- c("Genus")
    merge <- left_join(place, abundance)
    # 按照 mean_abundance 排序
    #merge <- merge[order(merge$freq.pred), ]
    write_tsv(merge, paste0(sloanDir,"residence_result.tsv"))
    
    #Plot the Sloan neutral model for all individuals
    
    gg_title<- paste0("   ",srcType," ",f_microbeType(),formatted_pair,"")
    cbPalette <- c("#e41a1c", "#377eb8", "#4daf4a")
    Table <- read_tsv(paste0(sloanDir,"residence_result.tsv")) %>% filter(!Genus == "unclassified")
    Plot <- ggplot(Table,aes(x=log(mean_abundance))) 
    Plot <- Plot + geom_point(aes(y=freq, colour=Partition,shape = isPathogen),alpha=0.8
                              , size = ifelse(Table$isPathogen == "TRUE", 2, 1)) 
    Plot <- Plot + scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 4)) # Customizing shapes: TRUE -> 叉形 (X), FALSE -> circle
    Plot <- Plot + labs(
                title = gg_title,  # 图形标题
                shape = "Pathogen associated",  # 为 shape 添加图例标题
                colour = "Partition"   # 为 color 添加图例标题
                  )
    Plot <- Plot + guides(size = "none") # 去除大小（size）的图例
    #Plot <- Plot + theme(legend.title = element_blank()) #隐藏图例标题
    Plot <- Plot + geom_line(aes(y=freq.pred))
    #Plot <- Plot + geom_line(aes(y=pred.lwr), linetype="dotted")
    #Plot <- Plot + geom_line(aes(y=pred.upr), linetype="dotted")
    Plot <- Plot + geom_ribbon(aes(ymin=pred.lwr, ymax=pred.upr),fill="#969696",alpha=0.5)
    #Plot <- Plot + scale_color_brewer(palette = "Dark2")
    Plot <- Plot + scale_color_manual(values = cbPalette)
    Plot <- Plot + theme_classic()
    Plot <- Plot + theme(plot.title = element_text(hjust = 0.5, face="bold"))
    #Plot <- Plot + ggtitle(gg_title)
    Plot <- Plot + xlab(paste0("log (mean relative abundance of " ,src_name, " Genus)")) +ylab(paste0("Occurrence frequency of ",dst_name," Genus"))
    Plot <- Plot + annotate(geom = "text", x = -Inf,y = Inf, hjust = 0, vjust = 1, 
                      label = sprintf("m: %s\nRsqr: %s",
                                      formatC(stats.all[1,"m"], format = "f", digits = 4), 
                                      formatC(stats.all[1,"Rsqr"], format = "f", digits = 4)), 
                      color = "black", size = 5,fontface = "italic")
    ggsave(paste0(sloanDir,srcType," ",f_microbeType(),"",formatted_pair,"_residence_result.pdf"), height=7, width=7, unit="in",dpi=1200)
    ggsave(paste0(sloanDir,srcType," ",f_microbeType(),"",formatted_pair,"_residence_result_sloan_neutral_model_all_individuals.pdf"), height=7, width=8, unit="in", dpi=1200)
  },error=function(e){
    # 打印错误信息
    message("函数f_gen_sloan执行错误，错误信息：", conditionMessage(e))
  })
}

#################################################################>
# 函数：f_gen_sloan2
# >生成sloan扩散模型2 废弃 ---- 
# 来源:https://blog.csdn.net/weixin_43367441/article/details/123682669
#################################################################>
f_gen_sloan2<-function(sloanDir,spp) {
  ##Fits the neutral model from Sloan et al. 2006 to an OTU table and returns several fitting statistics. Alternatively, will return predicted occurrence frequencies for each OTU based on their abundance in the metacommunity
  #Install the following packages if they haven't been availabled in your computer yet 
  library(Hmisc)
  library(minpack.lm)
  library(stats4)
  #using Non-linear least squares (NLS) to calculate R2:
  #spp: A community table with taxa as rows and samples as columns
  #通过参数传递 spp<-read.csv('spp.txt',head=T,stringsAsFactors=F,row.names=1,sep = "\t")
  #spp<-t(spp)
  N <- mean(apply(spp, 1, sum))
  p.m <- apply(spp, 2, mean)
  p.m <- p.m[p.m != 0]
  p <- p.m/N
  spp.bi <- 1*(spp>0)
  freq <- apply(spp.bi, 2, mean)
  freq <- freq[freq != 0]
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  p <- C.0[,2]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  d = 1/N
  ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
  m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
  m.fit #get the m value
  m.ci <- confint(m.fit, 'm', level=0.95)
  freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
  pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  Rsqr# get the R2 value
  
  #Optional: write 3 files: p.csv, freq.csv and freq.pred.csv
  if (!file.exists(paste0(sloanDir,"neutr_model"))) {
    dir.create(paste0(sloanDir,"neutr_model"), recursive = TRUE)  # recursive = TRUE 表示创建所有不存在的父目录
  }
  write.csv(p, file = paste0(sloanDir,"neutr_model/p.csv"))
  write.csv(freq, file = paste0(sloanDir,"neutr_model/freq.csv"))
  write.csv(freq.pred, file = paste0(sloanDir,"neutr_model/freq.pred.csv"))
  
  #Drawing the figure using grid package:
  #p is the mean relative abundance
  #freq is occurrence frequency
  #freq.pred is predicted occurrence frequency
  bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
  inter.col<-rep('black',nrow(bacnlsALL))
  inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#define the color of below points
  inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#define the color of up points
  library(grid)
  # 将绘制的图形保存为 PDF 文件
  pdf(paste0(sloanDir,"neutr_model/plot_output.pdf"))
  
  grid.newpage()
  pushViewport(viewport(h=0.6,w=0.6))
  pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
  grid.rect()
  grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=20,gp=gpar(col=inter.col,cex=0.7))
  grid.yaxis()
  grid.xaxis()
  grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='blue',lwd=2),default='native')
  
  grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='blue',lwd=2,lty=2),default='native') 
  grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='blue',lwd=2,lty=2),default='native')  
  grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
  grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
  #grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
  #grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
  #grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
  #grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
  #grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
  draw.text <- function(just, i, j) {
    grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
    #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
    #          gp=gpar(col="grey", fontsize=8))
  }
  x <- unit(1:4/5, "npc")
  y <- unit(1:4/5, "npc")
  draw.text(c("centre", "bottom"), 4, 1)
  
  dev.off()  # 关闭 PDF 设备
}

#############################################################################>
#函数结束，下面是程序开始
#############################################################################>

# >按分组生成sloans扩散模型 ----
f_gen_sloans<-function(getSampleArg,ps,ps_pathogen){
  f_printMsg("=================================")
  f_printMsg("06：sloan分析 ")
  f_printMsg("=================================")
  #读入病原菌数据库，生成病原菌列表Genus_Genus
  if (f_microbeType()=="Bac_") {
    pathogenDB<-read_tsv("data/pathogen_b2.tsv",show_col_types = FALSE)
    #pathogens<-paste0("",pathogenDB$Genus,"_",pathogenDB$Genus)
    pathogens<-paste0("",pathogenDB$Genus)
  } else {
    pathogenDB<-read_tsv("data/pathogen_f.tsv",show_col_types = FALSE)
    #pathogens<-paste0("",pathogenDB$Genus,"_",pathogenDB$Genus)
    pathogens<-paste0("",pathogenDB$Genus)
  }
  srcTypes <- c("all","pathogen")
  for (srcType in srcTypes) {
    ##################################################################>
    # 00 设置需要处理的ps对象 
    ##################################################################>
    if (srcType=="all") {
      ps_sameDepth<-ps
    } else {
      ps_sameDepth<-ps_pathogen
    }
    ##################################################################>
    # 01 按照分组生成phyloseq对象列表
    ##################################################################>
    # 01-1 将otu统一深度All samples must be rarefied to the same
    # 获取 OTU 表达矩阵
    otu_mat <- otu_table(ps_sameDepth)
    # 计算每个样本的 OTU 丰度和
    sums_per_sample <- apply(otu_mat, 2, sum)
    if (srcType != "all") { #病原菌可能有丰度太小的样本，要剔除太小的样本
      # 过滤丰度小于20的样本
        #  找到丰度和小于 20 的样本的样本名
      samples_to_remove <- names(sums_per_sample)[sums_per_sample < 20]  # 获取丰度和小于 50 的样本名
        #  从 phyloseq 对象中去除这些样本
      ps_sameDepth <- prune_samples(!sample_names(ps_sameDepth) %in% samples_to_remove, ps_sameDepth)
      
      # 再次计算过滤后每个样本的 OTU 丰度和
      sums_per_sample <- apply(ps_sameDepth@otu_table, 2, sum)
    }
    # 找出这些和的最小值
    min_sum <- min(sums_per_sample)
    # 使用指定深度进行稀释
    print(paste0("otu统一深度为：",min_sum))
    ps_sameDepth <- rarefy_even_depth(ps_sameDepth, sample.size = min_sum , rngseed=1)
    
    #根据分组情况，俩俩对比研究扩散情况
    formatted_pairs <- f_group_pairs()
    for (formatted_pair in formatted_pairs) {
      print(paste0("进行分组生成Sloan：",srcType,'->',formatted_pair))
      #根据 formatted_pair 变量，它包含一个以 "GroupA-GroupB" 格式表示的两个组名字符串，
      #我们可以根据这些组名从phyloseq 对象中筛选出对应的样本，
      #生成两个新的 phyloseq 对象：ps_src 和 ps_dst
      # Step 1: 从 formatted_pair 中提取两个组名
      #   formatted_pair <- "GroupA-GroupB"  # 示例字符串，实际情况可能来自变量
      group_names <- strsplit(formatted_pair, "-")[[1]]  # 将字符串分割为两个组名
      src_name <- group_names[1]  # 获取 GroupA 名称
      dst_name <- group_names[2]  # 获取 GroupB 名称
      # Step 2: 获取样本数据并提取组信息
      sample_data_df <- sample_data(ps_sameDepth)  # 获取样本数据
      # Step 3: 根据组名筛选样本
      src_samples <- sample_names(ps_sameDepth)[sample_data_df$Group == src_name]
      dst_samples <- sample_names(ps_sameDepth)[sample_data_df$Group == dst_name]
      # Step 4: 创建两个新的 phyloseq 对象
      ps_src <- prune_samples(src_samples, ps_sameDepth)  # GroupA 的样本
      ps_dst <- prune_samples(dst_samples, ps_sameDepth)  # GroupB 的样本
      
      # 原先代码，废弃
      # # 01-2 分组
      # ps_groups <- f_Phyloseq_splitByGroup(ps_sameDepth)
      # #查看“分组生成的phyloseq对象列表”
      # ps_groups[[1]]@sam_data$Group
      # ps_groups[[2]]@sam_data$Group
      # #ps_groups[[3]]@sam_data$Group
      # # 01-3 设置研究的两个ps对象（源和目的）
      # ps_src<-ps_groups[[1]] #比较的源ps，如地下
      # ps_dst<-ps_groups[[2]]  #比较的目的ps，如电梯
      
      
      ##################################################################>
      # 03将两个phyloseq对象，OTU表按照物种进行汇总，并保存文件
      #################################################################>
      sloanDir=paste0(getSampleArg$outDir,getSampleArg$samplesName,"/02Plot/06sloan/",
                      srcType,"/",formatted_pair,"/")
      dir.create(sloanDir, recursive = TRUE,showWarnings = FALSE)  

      ####02-1 处理src
      #按照物种汇总
      #ps_sumByTaxa<-f_microeco_otuSumByTaxa(ps_src,ps_dst,taxaOrder=7,
      #                        abund_dir=paste0(sloanDir,"/meco_abund")) #"Genus" 7
      ps_sumByTaxa<-f_microeco_otuSumByTaxa(ps_src,ps_dst,taxaOrder=6,
                                            abund_dir=paste0(sloanDir,"/meco_abund")) #"Genus" 7
      # 将数据框保存为 TSV 文件，不包含行名
      write.table(ps_sumByTaxa[[1]], file = paste0(sloanDir,"/src.tsv"), sep = "\t", row.names = FALSE)
      write.table(ps_sumByTaxa[[2]], file = paste0(sloanDir,"/dst.tsv"), sep = "\t", row.names = FALSE)
      
      ##################################################################>
      # 04将步骤03生成文件 调用函数：f_sloan 进行扩散研究
      ##################################################################>
      #pool <- NULL
      pool <- read_tsv(paste0(sloanDir,"/src.tsv"))
      spp <-  read_tsv(paste0(sloanDir,"/dst.tsv"))
      f_gen_sloan(sloanDir,srcType,f_microbeType(),formatted_pair,spp,pool,pathogens,
                  dst_name,src_name)
      #f_sloan2(sloanDir,spp)
      f_printMsg("函数f_gen_sloans执行完成（",srcType,"-",formatted_pair,")")
    }
  }
}
# 测试函数
#f_gen_sloans(getSampleArg,ps,ps_pathogen)
