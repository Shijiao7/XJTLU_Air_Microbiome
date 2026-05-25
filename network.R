####################################>
# >单域网络分析 ggClusterNet ####
####################################>

get_sci_palette <- function(unique_values, reverse = FALSE) {
  n <- length(unique_values)
  # 检测传入unique_values是否是季节数据开头（Aut, Win, Spr, Sum），如是reverse为True，取后面四个颜色
  if (n == 4 && all(startsWith(as.character(unique_values), c("Aut", "Win", "Spr", "Sum")))) {
    reverse <- TRUE
  }
  # 定义完整的颜色组（不截取）
  full_palette <- if (n <= 8) {
    c(
      "#3366CC",  # 强蓝色 (1)
      "#EEAA33",  # 橙黄色 (2)
      "#CC3311",  # 饱和红 (3)
      "#117733",  # 深绿色 (4)
      "#EE7733",  # 橙红色 (7)
      "#66CC99",  # 薄荷绿 (6)
      "#0099BB",  # 青蓝色 (8)
      "#994FCC"   # 紫罗兰色 (5)
    )
  } else if (n <= 12) {
    brewer.pal(12, "Set3")
  } else if (n <= 20) {
    viridis(n)
  } else {
    warning("Number of categories exceeds 20 - consider merging minor categories")
    rainbow(n)
  }
  
  # 根据是否倒序选择颜色
  if (reverse) {
    # 从后往前选取n个颜色
    palette <- rev(full_palette)[1:n]
  } else {
    # 从前往后选取n个颜色
    palette <- full_palette[1:n]
  }
  
  # 确保颜色数量正确（当n > length(full_palette)时）
  palette <- palette[1:n]
  
  return(setNames(palette, unique_values))
}

f_ggClusterNet<-function(getSampleArg,ps,ps_pathogen,N,taxFill,displaySpeciesDegrees,displaySpeciesDegreesPathogen){
  # 网络分析主函数#-------->
  # https://mp.weixin.qq.com/s/aDnmmcNxoGRKasTP1U1Jyw
  #library(ggClusterNet)
  #library(phyloseq)
  #library(tidyverse)
  #library(igraph)
  #library(dplyr)
  #outDir<-"./result_big_1000/"
  
  f_printMsg("==========================================")
  f_printMsg("08-1：ggClusterNet 单网络分析 ")
  f_printMsg("==========================================")
  srcTypes <- c("all","pathogen")
  for (srcType in srcTypes) {
    ##################################################################>
    # 01 设置需要处理的ps对象和输出目录 
    ##################################################################>
    if (srcType=="all") {
      ps_net<-f_ps_addGroupOrder(ps)
    } else {
      ps_net<-f_ps_addGroupOrder(ps_pathogen)
    }
    #户外样本冬季只有两个：231222_8B 240229_8B，样本不足，为避免影响网络，剔除240229_8B
    # 剔除特定样本
    ps_net <- phyloseq::subset_samples(ps_net, !SampleName %in% c("231222_8B", "240229_8B"))
    #增加一列表示是否病原菌
    ps_net<-f_ps_tax_addIsPathogen(ps_net)  #增加一列表示是否病原菌
    #ps_net的tax表增加一列 only_Species取Species列_后面的字符串，如没有_就取Species列  
    #如果isPathogen列是TRUE，则在后面加*号
    tax_table(ps_net) <- tax_table(ps_net) %>% 
      as.data.frame() %>%
      mutate(
        Genus = paste0(
          #ifelse(str_detect(Species, "_"), str_extract(Species, "[^_]+$"), Species),
          ifelse(str_detect(Genus, "_"), str_extract(Genus, "[^_]+$"), Genus),
          ifelse(isPathogen == TRUE, "*", "")
        )
      ) %>%
      as.matrix()
    #view(ps_net@tax_table)
    netDir=paste0(getSampleArg$outDir,getSampleArg$samplesName,"/02Plot/08ggClusterNet/",
                  srcType,"/")
    dir.create(netDir, recursive = TRUE,showWarnings = FALSE) 
    ##################################################################>
    # 02 调用ggclusterNet绘制网络图 
    ##################################################################>
    
    ##################################################################>
    #Part1 将PS的各组分别生成网络模块图 
    ##################################################################>
    if(1==12) { #调试控制用，不执行Prat1，改成1!=1
      tryCatch( 
        {
          #assign("ps_net_debug",ps_net, envir = globalenv()) #调试用
          # Part1.1 拆成各分组，按照季节展示
          try(f_net_generate_model_by_group("Season",ps_net,"_All Group",netDir)) #按照Season分组
          group_levels <- levels(ps_net@sam_data$Group)
          i=1
          for (for_group in group_levels) {
            # 针对每个组的操作，例如：
            print(paste("生成组:", for_group," 的各组网络模块图"))
            # 进一步操作，如子集数据：
            assign("for_group",for_group, envir = globalenv()) #调试用
            subset_data <- phyloseq::subset_samples(ps_net, Group == for_group)
            # ...其他处理...
            try(f_net_generate_model_by_group("Season",subset_data,
                                              paste0(LETTERS[i],"-",for_group),netDir)) #按照Season分组
            i=i+1
          }
          # Part1.2 拆成各季节，按照分组展示
          try(f_net_generate_model_by_group("Group",ps_net,"_All Season",netDir)) #按照Group分组
          season_levels <- levels(ps_net@sam_data$Season)
          i=1
          for (for_season in season_levels) {
            #for_season="Autumn"
            #netDir="."
            # 针对每个季节的操作，例如：
            print(paste("生成季节:", for_season," 的各组网络模块图"))
            # 进一步操作，如子集数据：
            print(paste0("进一步操作，如子集数据 for_season1:",for_season))
            assign("for_season",for_season, envir = globalenv()) #调试用
            subset_data <- phyloseq::subset_samples(ps_net, Season == for_season)
            # ...其他处理...
            print(paste0("...其他处理..."))
            print(paste0("...其他处理... for_season2:",for_season))
            try(f_net_generate_model_by_group("Group",
                                              subset_data,paste0( LETTERS[i],"-",for_season),
                                              netDir)) #按照Group分组
            i=i+1
          }
        },error=function(e){
          f_printMsg("执行函数f_ggClusterNet中f_net_generate_model_by_group出错：")
          print(e)
        }
      )
    } #if(1==1) { #调试控制用，不执行Prat1，改成1!=1
    ##################################################################>
    #End Part1 将PS的各组分别生成网络模块图 
    ##################################################################>
    
    ##################################################################>
    #Part2 敲除节点处理
    ##################################################################>
    if(1==12) { #调试控制用，不执行Prat2，改成1!=1
      outFileName<-paste0(netDir,srcType," ",getSampleArg$samplesName," attack.pdf")
      igraph_list<-f_gen_igraph_by_group(ps_net,N=N,r.threshold = 0.2, 
                                         p.threshold = 0.05,method="pearson")
      f_attack_ver2(igraph_list,outFileName)
      if (srcType=="all") {
        # 获取分类信息
        tax_table_ps <- tax_table(ps_net)
        
        # 选择 'isPathogen' 为 TRUE 的行
        pathogen_ASVs <- rownames(tax_table_ps)[tax_table_ps[, "isPathogen"] == TRUE]
        # 从phyloseq对象中删除这些ASV
        ps_clean <- prune_taxa(!rownames(tax_table_ps) %in% pathogen_ASVs, ps_net)
        igraph_list<-f_gen_igraph_by_group(ps_clean,N=N,r.threshold = 0.2, 
                                           p.threshold = 0.05,method="pearson")
        outFileName2<-paste0(netDir,srcType," ",getSampleArg$samplesName," attack_noPathogen.pdf")
        f_attack_ver2(igraph_list,outFileName2)
      }
    } #if(1==1) { #调试控制用，不执行Prat2，改成1!=1
    ##################################################################>
    #End Part2 敲除节点处理
    ##################################################################>    
    
    ##################################################################>  
    #Part 3 按物种间相互作用
    ##################################################################>
    if(1==1) { #调试控制用，不执行Prat3，改成1!=1
      tryCatch({
        tab.r = network.pip_modify(
          ps = ps_net,N = N,# ra = 0.05,
          big = TRUE,select_layout = FALSE,layout_net = "model_maptree2",
          r.threshold = 0.6,p.threshold = 0.05,maxnode = 2,method = "spearman",
          label = FALSE,lab = "Species",group = "Group",fill = taxFill,size = "igraph.degree",
          zipi = TRUE,ram.net = TRUE,clu_method = "cluster_fast_greedy",
          step = 100,R=10,ncpus = 2)
        if (is.null(tab.r)) {
          next
        }
        #  建议保存一下输出结果为R对象，方便之后不进行相关矩阵的运算，节约时间
        saveRDS(tab.r,paste0(netDir,srcType," network.pip.sparcc.rds"))
        #tab.r<-readRDS("result\\bacteria\\20241019\\15 group_by_floor\\02Plot\\08ggClusterNet\\all\\all network.pip.sparcc.rds")
        
        # 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
        dat = tab.r[[2]]
        cortab = dat$net.cor.matrix$cortab
        saveRDS(cortab,paste0(netDir,srcType," cor.matrix.all.group.rds"))
        
        #将绘制结果保存到pdf文件
        plot = tab.r[[1]]
        p1 = plot[[1]]
        #调整画板的循序
        p1<-p1+ facet_wrap(. ~ paste0(Group," (",str_extract(label, "nodes: \\d+"),")"), 
                           scales = "free_y")
        
        #增加标题      
        p1<-p1+ggtitle(paste("Network (",netDir,")"))+
          theme(plot.title = element_text(size = 6,hjust = 0))
        # 对网络图的点进行处理，点的大小按照度的大小，度的大小超过displaySpeciesDegrees则显示物种分类的名称
        #取出点数据
        tmp_node<-tab.r[[2]]$net.cor.matrix$node
        tmp_node$isPathogen<-tmp_node$Genus %in% pathogenDB$Genus
        #tmp_node$lab<-tmp_node[["Species"]]
        #画点
        # 使用数据中的唯一值动态生成调色板
        #unique_values <- unique(tmp_node[[taxFill]])
        #color_palette <- setNames(RColorBrewer::brewer.pal(length(unique_values), "Set2"), unique_values)
        # 为特定值（如 "+" 和 "-"）手动指定颜色
        #color_palette["+"] <- "#FF0000"  # 设置 "+" 的颜色为番茄色
        #color_palette["-"] <- "#0000FF"  # 设置 "-" 的颜色为钢蓝色
        
        # 获取唯一值
        #unique_values <- unique(tmp_node[[taxFill]])
        # 获取taxFill列的唯一值并按degree降序排列
        # 首先替换原列增加度显示内容
        # 将字符串"NA"转换为真正的NA
        tmp_node[[taxFill]] <- ifelse(tmp_node[[taxFill]] == "NA", 
                                      NA_character_, 
                                      tmp_node[[taxFill]])
        # 替换原列并过滤NA值
        tmp_node[[taxFill]] <- ifelse(is.na(tmp_node[[taxFill]]), 
                                      NA,  # 保持NA不变，后续过滤
                                      paste0(tmp_node[[taxFill]], " (", tmp_node$igraph.degree, ")"))
        
        # 创建数据框并过滤NA记录
        unique_df <- data.frame(
          taxon = as.character(tmp_node[[taxFill]]),
          degree = tmp_node$igraph.degree,
          isPathogen=tmp_node$isPathogen,
          stringsAsFactors = FALSE
        ) %>%
          filter(!is.na(taxon)) %>%  # 剔除taxon为NA的记录
          distinct(taxon, .keep_all = TRUE) %>%
          arrange(desc(degree))%>%
          filter(isPathogen == TRUE)
        
        # 按degree排序后的unique_values
        sorted_unique <- unique_df$taxon
        
        # # 如果唯一值超过 8，使用 Set3 调色板
        # if (length(unique_values) <= 12) {
        #   color_palette <- setNames(RColorBrewer::brewer.pal(12, "Set3"), unique_values)
        # } else {
        #   # 如果唯一值超过 12，使用更多颜色或者手动指定颜色
        #   color_palette <- setNames(RColorBrewer::brewer.pal(12, "Set3"), unique_values[1:12])
        # }
        
        if (length(sorted_unique) <= 12) {
          color_palette <- setNames(RColorBrewer::brewer.pal(12, "Paired"), sorted_unique)
        } else {
          # 如果唯一值超过 12，使用更多颜色或者手动指定颜色
          color_palette <- setNames(RColorBrewer::brewer.pal(12, "Paired"), sorted_unique[1:12])
        }
        patho<-c("TRUE","FALSE","+","-")
        color_palette2 <- setNames(RColorBrewer::brewer.pal(4, "Paired"), patho)
        color_palette2["with_star"] <- "purple"  # 设置 "+" 的颜色为番茄色
        color_palette2["no_star"] <- "black"  # 设置 "-" 的颜色为钢蓝色
        
        
        
        # 为特定值（如 "+" 和 "-"）手动指定颜色
        color_palette["+"] <- "gray80"  # 设置 "+" 的颜色为番茄色
        color_palette["-"] <- "#FF0000"  # 设置 "-" 的颜色为钢蓝色
        color_palette2["+"] <- "gray80"  # 设置 "+" 的颜色为番茄色
        color_palette2["-"] <- "#FF0000"  # 设置 "-" 的颜色为钢蓝色
        
        
        color_palette_isPathogen <- c(
          "with_star" = "blue",    # 含星号的标签为红色
          "no_star" = "gray40",   # 不含星号为灰色
          "+"       = "yellow",  #连线+颜色
          "-"       = "red"      #连线-颜色
        )
        # nc.DCA <- nc.attack(p1)
        # png("connectivity_cheek.png",width=4,height=4,units="in",res=1200)
        # points(seq(0,0.8,len=length(nc.DCA)),nc.DCA,type='l')
        # dev.off()
        # 处理标签 显示物种名称和指定颜色
        #tmp_node$lab<-tmp_node[["Species"]]
        tmp_node$lab<-tmp_node[[taxFill]]
        if (srcType=="bacteria") {
          labelDispDegrees<-displaySpeciesDegrees
        } else {
          labelDispDegrees<-displaySpeciesDegreesPathogen
        }
        #生成填充变量fill_group
        tmp_node <- tmp_node %>%
          mutate(
            Fill_group = ifelse(grepl("\\*", lab), "with_star", "no_star")
          )
        #f_debug_save_var("tmp_node")
        #tmp_node<-f_debug_load_var("tmp_node")
        #画图
        p1 <- p1 + 
          geom_point(aes(X1, X2, 
                         #fill = isPathogen, 
                         fill = Fill_group,  # 与标签颜色逻辑一致
                         size = igraph.degree/4), 
                     pch = 21, data = tmp_node, 
                     color = "gray60", alpha = 0.8) + # 边框颜色固定
          scale_fill_manual(
            values = color_palette_isPathogen,  # 与标签颜色一致
            labels = c("with_star" = "* Pathogen", "no_star" = "Normal")
            ,guide = "none"  # 隐藏颜色图例（避免重复）
          ) 
        #scale_size_continuous(range = c(2, 8)) 
        
        # 标签图层
        p1<-p1+geom_text_repel(aes(X1, X2,
                                   label = ifelse(igraph.degree >= labelDispDegrees, lab, ""),
                                   size = igraph.degree/10,
                                   #color = !!sym(taxFill),
                                   color = ifelse(grepl("\\*", lab), "with_star", "no_star"),  # 检测星号
                                   fontface = ifelse(isPathogen == "TRUE", "bold.italic", "plain")),  # 如果isPathogen为TRUE则设置为粗体
                               data = tmp_node
        ) +
          # 设置标记颜色
          scale_color_manual(
            values = color_palette2,  # 与填充色一致
            labels = c("with_star" = "* Pathogen", "no_star" = "Normal")
            #,guide = "none"  # 隐藏颜色图例（避免重复）
          ) +
          # end 3. 统一颜色标度
          # 统一大小标度
          #scale_size_continuous(range = c(1, 5))
          scale_size_continuous(
            range = c(1, 8),  # 合并节点和标签的大小范围
            name = "Degree"   # 图例标题
          )
        # p1 <- p1 +
        #   geom_label_repel(aes(X1, X2,
        #                        label = ifelse(igraph.degree >= displaySpeciesDegrees, lab, ""),
        #                        size = igraph.degree,
        #                        color = !!sym(taxFill),
        #                        fontface = ifelse(isPathogen == "TRUE", "italic", "plain")),
        #                    data = tmp_node,
        #                    fill = ifelse(tmp_node$isPathogen == "TRUE", "gray80", NA),  # 仅为 isPathogen == TRUE 的标签设置灰色背景
        #                    label.padding = 0.3,  # 标签文字与背景之间的内边距
        #                    box.padding = 0.5) +  # 标签背景的内边距
        #   scale_color_manual(values = color_palette) +  # 设置标签的颜色
        #   scale_size_continuous(range = c(1, 5))  # 设置标签的大小范围
        
        # p1 <- p1 +
        #   geom_text_repel(aes(X1, X2,
        #                       label = ifelse(igraph.degree >= displaySpeciesDegrees, lab, ""),
        #                       size = igraph.degree,
        #                       color = !!sym(taxFill),
        #                       fontface = ifelse(isPathogen == "TRUE", "italic", "plain")),
        #                   data = tmp_node,
        #                   box.padding = 0.3,  # 设置标签与背景的内边距
        #                   label.padding = 0.3,  # 设置标签文字与背景的内边距
        #                   fill = ifelse(tmp_node$isPathogen == "TRUE", "gray80", NA),  # 只有isPathogen为TRUE时有灰色背景
        #                   label.size = 0,  # 设置label.size为0去除边框
        #                   label.r = 0) +  # 去除背景圆角的效果
        #   scale_color_manual(values = color_palette) +  # 设置标签的颜色
        #   scale_size_continuous(range = c(1, 5))  # 设置标签的大小范围
        
        
        ggsave(paste0(netDir,srcType," ",getSampleArg$samplesName," plot.network.pdf"),p1,width = 12,height = 5)
        ggsave(paste0(netDir,srcType," ",getSampleArg$samplesName," plot.network2.pdf"),p1,width = 16,height = 10)
        
        p2 = plot[[2]]
        p2<-p2+ggtitle(paste("zipi (",netDir,")"))+
          theme(plot.title = element_text(size = 6,hjust = 0))
        ggsave(paste0(netDir,srcType," ",getSampleArg$samplesName," plot.zipi.pdf"),p2,width = 12,height = 5)
        p3 = plot[[3]]
        p3<-p3+ggtitle(paste("Random (",netDir,")"))+
          theme(plot.title = element_text(size = 6,hjust = 0))
        ggsave(paste0(netDir,srcType," ",getSampleArg$samplesName," plot.random.pdf"),p3,width = 12,height = 5)
        
        #将运算节点和边保存
        dat = tab.r[[2]]
        node = dat$net.cor.matrix$node
        # 指定要检查的列
        cols_to_check <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genus_Species")
        # 将指定列中的 NA 替换为空字符
        node[cols_to_check] <- lapply(node[cols_to_check], function(x) ifelse(is.na(x), "", x))
        # 将指定列中的字符串 "NA" 替换为空字符
        node[cols_to_check] <- lapply(node[cols_to_check], function(x) ifelse(x == "NA", "", x))
        
        # 使用 split() 按照 Group 列拆分数据框
        split_df <- split(node, node$Group)
        # 循环保存每个子数据框为 CSV 文件
        for(group in names(split_df)) {
          FileName <- paste(netDir,srcType," ",getSampleArg$samplesName,"_",group,"_node","",".csv", sep = "")
          write_csv(split_df[[group]],FileName)
        }
        
        edge = dat$net.cor.matrix$edge
        # 创建 Gephi 边文件的数据框
        gephi_edges <- data.frame(
          Source = edge$OTU_1,  # 起始节点
          Target = edge$OTU_2,  # 目标节点
          Weight = edge$weight,  # 权重
          Correlation = edge$cor,  # 相关性（可选）
          Group = edge$Group  # 相关性（可选）
        )
        # 使用 split() 按照 Group 列拆分数据框
        split_df <- split(gephi_edges, gephi_edges$Group)
        # 循环保存每个子数据框为 CSV 文件
        for(group in names(split_df)) {
          FileName <- paste(netDir,srcType," ",getSampleArg$samplesName,"_",group,"_edge","",".csv", sep = "")
          write_csv(split_df[[group]],FileName)
        }
      },error=function(e){
        f_printMsg("执行函数f_ggClusterNet出错：")
        print(e)
      })  
      #head(edge)
      #head(node)
    } #if(1==1) { #调试控制用，不执行Prat3，改成1!=1
    ##################################################################>  
    #End Part 3 按物种间相互作用
    ##################################################################>
  }  #End   for (srcType in srcTypes)
}


# >将PS的各组分别生成网络模块图 --------------------------------------------------------
f_net_generate_model_by_group<-function(groupBy,arg_ps,seasonName,netDir) {
  if (groupBy=="Group") {
    groupNames<-getSampleArg$sampleFilterGroupName
  } else { #按照季节来
    groupNames<- c("Autumn", "Winter", "Spring", "Summer")
  }
  pdf(paste0(netDir,"/network_by_",groupBy,"_model_",seasonName,".pdf"),width = 10,height = 8)
  p_lists<- list()
  for (groupName in groupNames){ #循环将各组进行网络模块图生成
    print(paste0("分组规则：",groupBy,",循环将各组进行网络模块图生成,现在生成的组为:",groupName))
    if (groupBy=="Group") {
      a_ps<-f_filter_ps_by_group(arg_ps,groupName,group_var="Group_original")#取groupName现成ps对象
    } else { #按季节来
      a_ps<-f_filter_ps_by_group(arg_ps,groupName,group_var="Season")
    }
    if ( is.null(a_ps) ) {
      next
    }
    # cor_Big_micro2 增加了标准化方法和p值矫正方法
    result = cor_Big_micro2(ps = a_ps,
                            N = 800,
                            r.threshold=0.85,
                            p.threshold=0.05,
                            method = "pearson",
                            scale = FALSE
    )
    
    #--提取相关矩阵
    cor = result[[1]]
    #dim(cor)
    #> [1] 1000 1000
    
    # model_igraph2
    result2 <- model_igraph2(cor = cor,
                             method = "cluster_fast_greedy",
                             seed = 12
    )
    node = result2[[1]]
    #dim(node)
    dat = result2[[2]]
    #head(dat)
    
    tem = data.frame(mod = dat$model,col = dat$color) %>%  
      dplyr::distinct( mod, .keep_all = TRUE)  
    col = tem$col
    names(col) = tem$mod
    
    #---node节点注释#-----------
    otu_table = as.data.frame(t(vegan_otu(a_ps)))
    tax_table = as.data.frame(vegan_tax(a_ps))
    nodes = nodeadd(plotcord =node,otu_table = otu_table,tax_table = tax_table)
    #head(nodes)
    #-----计算边#--------
    edge = edgeBuild(cor = cor,node = node)
    colnames(edge)[8] = "cor"
    #head(edge)
    tem2 = dat %>% 
      dplyr::select(OTU,model,color) %>%
      dplyr::right_join(edge,by =c("OTU" = "OTU_1" ) ) %>%
      dplyr::rename(OTU_1 = OTU,model1 = model,color1 = color)
    #head(tem2)
    
    tem3 = dat %>% 
      dplyr::select(OTU,model,color) %>%
      dplyr::right_join(edge,by =c("OTU" = "OTU_2" ) ) %>%
      dplyr::rename(OTU_2 = OTU,model2 = model,color2 = color)
    #head(tem3)
    
    tem4 = tem2 %>%inner_join(tem3)
    #head(tem4)
    
    edge2 = tem4 %>% mutate(color = ifelse(model1 == model2,as.character(model1),"across"),
                            manual = ifelse(model1 == model2,as.character(color1),"#C1C1C1")
    )
    
    
    col_edge = edge2 %>% dplyr::distinct(color, .keep_all = TRUE)  %>% 
      dplyr::select(color,manual)
    col0 = col_edge$manual
    names(col0) = col_edge$color
    
    library(ggnewscale)
    
    p1 <- ggplot() + geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2,color = color),
                                  data = edge2, size = 1) +
      scale_colour_manual(values = col0) 
    # ggsave("./cs1.pdf",p1,width = 16,height = 14)
    p2 = p1 +
      new_scale_color() +
      geom_point(aes(X1, X2,color =model), data = dat,size = 4) +
      scale_colour_manual(values = col) +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      theme(panel.background = element_blank()) +
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme(legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white",  colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())+
      ggtitle(paste0(seasonName," network by model (",groupName,")"))
    print(p2)
    p_title<-paste0(seasonName," network by model (",groupName,")")
    if(!(grepl("Outdoor", p_title) & grepl("Winter", p_title) )) {  #判断户外及冬天，样本不够，剔除不展示
      p_lists <- append(p_lists, list(p2)) #将批（ggplot加入列表）
    }
  }
  dev.off()
  # 使用patchwork画组合
  print("将得到的网络模块图进行组合")
  #print(p_lists)
  try(
    {
      library(patchwork)
      pdf(paste0(netDir,"/network_by_",groupBy,"_model_",seasonName,"_merge.pdf"),width = 20,height = 16)
      p3<-wrap_plots(p_lists, nrow = 2)
      print(p3)
      dev.off()
      if (!is.null(dev.list())) dev.off()
    })
}
#f_net_generate_model_by_group(ps,"result\\bacteria\\20241019\\15 group_by_floor\\02Plot\\08ggClusterNet\\all")


####################################>
# >细菌和真菌跨域网络 ggClusterNet_BacAndFung ####
# displaySpeciesDegrees 为节点度大于等于多少显示物种的Species名称
####################################>
f_ggClusterNet_BacAndFung<-function(srcType,samplesName,ver_bac,ver_fung,
                                    N,r.threshold,p.threshold,displaySpeciesDegrees,skipGroup,r.threshold2,displaySpeciesDegreesPathogen){
  
  
  natcon <- function(ig) {
    N <- vcount(ig)
    adj <- get.adjacency(ig)
    evals <- eigen(adj)$value
    nc <- log(mean(exp(evals)))
    nc / (N - log(N))
  }
  nc.attack <- function(ig) {
    hubord <- order(rank(igraph::degree(ig)), decreasing=TRUE)
    sapply(1:round(vcount(ig)*.95), function(i) {
      ind <- hubord[1:i]
      tmp <- delete_vertices(ig, V(ig)$name[ind])
      natcon(tmp)
    })
  }
  #srcType="all";N=100;r.threshold=0.8;p.threshold=0.05;skipGroup=NULL;samplesName="Aut_p ^ Win_p ^ Spr_p ^ Sum_p"
  f_printMsg("==========================================")
  f_printMsg("08-2：ggClusterNet 细菌和真菌跨域网络分析(",srcType,")")
  f_printMsg("==========================================")
  
  # 1、建立输出目录
  if (srcType=="all") {
    outDir<-paste0("result/net_bacteria2fungi/",samplesName,"/ggClusterNet/")
    dir.create(outDir,recursive = T,showWarnings=FALSE)
  } else {
    outDir<-paste0("result/net_bacteria2fungi/",samplesName,"/ggClusterNet/_pathogen")
    dir.create(outDir,recursive = T,showWarnings=FALSE)
  }
  outDir1=paste0(outDir,"01细菌和真菌内部及跨域网络/")
  dir.create(outDir1,recursive = T,showWarnings=FALSE)
  
  
  f_printMsg("细菌和真菌跨域网络分析输出目录:",outDir)
  if (srcType=="all") {
    #读入已生成好的细菌ps对象
    rdsFile_bac=paste0( "result/bacteria/",ver_bac,"/",samplesName,"/01samples/","/ps/Phyloseq-",samplesName,"-",ver_bac,".rds" )
    ps_bac<-readRDS(rdsFile_bac)
    ps_bac <- f_sam_convert_factor(ps_bac)
    #读入已生成好的真菌ps对象
    rdsFile_fung=paste0( "result/fungi/",ver_fung,"/",samplesName,"/01samples/","/ps/Phyloseq-",samplesName,"-",ver_fung,".rds" )
    ps_fung<-readRDS(rdsFile_fung)
    ps_fung <- f_sam_convert_factor(ps_fung)
    # 细菌和真菌 处理成共同有的样本对象
    ps1names<-sample_names(ps_bac)
    ps2names<-sample_names(ps_fung)
    common_samples<-intersect(ps1names,ps2names)
    ps_bac_common<-prune_samples(common_samples,ps_bac)
    ps_fung_common<-prune_samples(common_samples,ps_fung)
  } else {
    #读入已生成好的细菌对应病原菌ps对象
    rdsFile_bac_pathogen=paste0( "result/bacteria/",ver_bac,"/",samplesName,"/01samples/","/ps/Phyloseq_pathogen-",samplesName,"-",ver_bac,".rds" )
    ps_bac_pathogen<-readRDS(rdsFile_bac_pathogen)
    ps_bac_pathogen <- f_sam_convert_factor(ps_bac_pathogen)
    #读入已生成好的真菌对应病原菌ps对象
    rdsFile_fung_pathogen=paste0( "result/fungi/",ver_fung,"/",samplesName,"/01samples/","/ps/Phyloseq_pathogen-",samplesName,"-",ver_fung,".rds" )
    ps_fung_pathogen<-readRDS(rdsFile_fung_pathogen)
    ps_fung_pathogen <- f_sam_convert_factor(ps_fung_pathogen)
    # 细菌和真菌 处理成共同有的样本对象
    ps1names<-sample_names(ps_bac_pathogen)
    ps2names<-sample_names(ps_fung_pathogen)
    common_samples<-intersect(ps1names,ps2names)
    ps_bac_common<-prune_samples(common_samples,ps_bac_pathogen)
    ps_fung_common<-prune_samples(common_samples,ps_fung_pathogen)
  }
  #  去掉 ps_bac_common 中丰度为 0 的特征
  ps_bac_common_pruned <- prune_taxa(taxa_sums(ps_bac_common) > 0, ps_bac_common)
  #  去掉 ps_fung_common 中丰度为 0 的特征
  ps_fung_common_pruned <- prune_taxa(taxa_sums(ps_fung_common) > 0, ps_fung_common)
  
  # 4、合并细菌和真菌ps对象
  ps.merge <- ggClusterNet::merge16S_ITS(ps16s = ps_bac_common,
                                         psITS = ps_fung_common,
                                         N16s = N,
                                         NITS = N)
  ps.merge<-f_ps_addGroupOrder(ps.merge) #给分组列前面增加序号
  #户外样本冬季只有两个：231222_8B 240229_8B，样本不足，为避免影响网络，剔除240229_8B
  # 剔除特定样本
  ps.merge <- phyloseq::subset_samples(ps.merge, !SampleName %in% c("231222_8B", "240229_8B"))
  # 查看当前的列名
  colnames(tax_table(ps.merge))
  
  # 假设"filed"是现有的列名，将其改为"Pathogenicity"
  # 先找到"filed"列的位置
  current_names <- colnames(tax_table(ps.merge))
  col_index <- which(current_names == "filed")
  
  if(length(col_index) > 0) {
    # 修改列名
    colnames(tax_table(ps.merge))[col_index] <- "Pathogenicity"
    cat("已将第", col_index, "列 'filed' 改名为 'Pathogenicity'\n")
  } else {
    cat("警告：未找到名为 'filed' 的列。当前列名为：", paste(current_names, collapse = ", "), "\n")
  }
  
  # 验证修改
  colnames(tax_table(ps.merge))
  ###################################################################>
  ## Part 1 网络稳定性评估，不需要可以不执行
  ###################################################################>
  if (1!=1){ #网络稳定性评估  如不需要改成 1!=1
    outFileName<-paste0(outDir1,srcType," ",getSampleArg$samplesName,"attack.pdf")
    igraph_list<-f_gen_igraph_by_group(ps.merge,N=200,r.threshold = r.threshold , 
                                       p.threshold = p.threshold,method="pearson")
    f_attack_ver2(igraph_list,outFileName)
  }
  ###################################################################>
  ##End  Part 1 网络稳定性评估，不需要可以不执行
  ###################################################################>
  
  ###################################################################>
  # Start 增加是否病原菌列和生成Pathogenicity 列
  # 合并病原菌数据库
  pathogenDB_b <- pathogenDB <- f_load_pathogenDB("bacteria")
  #pathogenDB_b<-as.data.frame(pathogenDB_b$Genus_Species)
  pathogenDB_b<-as.data.frame(pathogenDB_b$Genus)
  #colnames(pathogenDB_b) <- "Genus_Species"
  colnames(pathogenDB_b) <- "Genus"
  pathogenDB_f <- pathogenDB <- f_load_pathogenDB("fungi")
  #pathogenDB_f<-as.data.frame(pathogenDB_f$Genus_Species)
  pathogenDB_f<-as.data.frame(pathogenDB_f$Genus)
  #colnames(pathogenDB_f) <- "Genus_Species"
  colnames(pathogenDB_f) <- "Genus"
  pathogenDB<-rbind(pathogenDB_b,pathogenDB_f)
  # 创建 ispathogen 列，检查 Species 是否出现在 pathogenDB$Genus_Species 中
  tax_df <- as.data.frame(ps.merge@tax_table)
  #tax_df$ispathogen <- ifelse(tax_df$Species %in% pathogenDB$Genus_Species, TRUE, FALSE)
  tax_df$ispathogen <- ifelse(tax_df$Genus %in% pathogenDB$Genus, TRUE, FALSE)
  tax_df$Pathogenicity  <- ifelse(tax_df$ispathogen == TRUE, 
                                  paste0(tax_df$Pathogenicity , "_pathogen"), 
                                  tax_df$Pathogenicity )
  # 将更新的 tax_df 再次赋值给 ps.merge@tax_table
  ps.merge@tax_table <- tax_table(as.matrix(tax_df))
  #view(ps.merge@tax_table)
  
  #查看生成的ps.merge
  f_printMsg("合并后的对象如下：")
  print(ps.merge)
  # End 增加是否病原菌列
  ###################################################################>
  
  ##############################################################>
  ##################################################################>
  #  Part 2 调用 跨域网络交互网络图 
  ##################################################################>
  
  dir.create(outDir1,recursive = T,showWarnings = FALSE)
  # 重新设置 group 和 category 的显示顺序
  #ps.merge@sam_data$Group <- factor( ps.merge@sam_data$Group , 
  #                            levels = sort( unique(sample_data(ps.merge)$Group,decreasing =TRUE ) ))  # 按配置循序排序
  
  taxFill="Pathogenicity "
  tab.r = network.pip_modify(
    ps = ps.merge,N = N,# ra = 0.05,
    big = TRUE,select_layout = FALSE,layout_net = "model_maptree2",
    r.threshold = r.threshold2,p.threshold = p.threshold,maxnode = 2,method = "spearman",
    label = FALSE,lab = "Species",group = "Group",fill = taxFill,size = "igraph.degree",
    zipi = TRUE,ram.net = TRUE,clu_method = "cluster_fast_greedy",
    step = 100,R=10,ncpus = 2,attack=FALSE
  )
  
  
  #  建议保存一下输出结果为R对象，方便之后不进行相关矩阵的运算，节约时间
  saveRDS(tab.r,paste0(outDir1,srcType," network.pip.sparcc.rds"))
  #tab.r<-readRDS(paste0(outDir1,srcType," network.pip.sparcc.rds"))
  # 大型相关矩阵跑出来不容易，建议保存，方便各种网络性质的计算
  dat = tab.r[[2]]
  cortab = dat$net.cor.matrix$cortab
  saveRDS(cortab,paste0(outDir1,srcType," cor.matrix.all.group.rds"))
  
  #将绘制结果保存到pdf文件
  plot = tab.r[[1]]
  p1 = plot[[1]]
  #p1<-p1+ggtitle(paste("Network (",outDir1,")"))+
  #  theme(plot.title = element_text(size = 6,hjust = 0))
  library(stringr)
  #调整画板的循序
  p1<-p1+ facet_wrap(. ~ paste0(Group," (",str_extract(label, "nodes: \\d+"),")"), scales = "free_y")
  
  clean_facet_label <- function(labels) {
    sapply(labels, function(label) {
      # 去掉前2个字符
      if(nchar(label) > 2) {
        cleaned <- substr(label, 3, nchar(label))
      } else {
        cleaned <- label
      }
      # 首字母大写
      paste0(toupper(substr(cleaned, 1, 1)), substr(cleaned, 2, nchar(cleaned)))
    })
  }
  
  # 使用labeller参数
  p1 <- p1 + 
    facet_wrap(
      . ~ paste0(Group, " (", str_extract(label, "nodes: \\d+"), ")"), 
      scales = "free_y",
      labeller = as_labeller(clean_facet_label)
    )
  
  res = Robustness.Random.removal(ps = ps.merge, Top = N, r.threshold= r.threshold, p.threshold=p.threshold, method = "spearman" )
  
  res2= Robustness.Targeted.removal(ps = ps, Top = N, degree = TRUE, zipi = FALSE, r.threshold= r.threshold, p.threshold=p.threshold, method = "spearman")
  
  ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName," robon.pdf"),res[[1]],width = 12,height = 5)
  ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName," robon2.pdf"),res2[[1]],width = 12,height = 5)
  
  # 对网络图的点进行处理，点的大小按照度的大小，度的大小超过displaySpeciesDegrees则显示物种分类的名称
  #取出点数据
  tmp_node<-tab.r[[2]]$net.cor.matrix$node
  tmp_node$lab<-tmp_node[["Genus"]]
  
  # 在绘图前，在数据中创建一个新变量用于图例
  tmp_node <- tmp_node %>%
    mutate(Pathogenicity = case_when(
      Pathogenicity == "bac" ~ "Bacteria",
      Pathogenicity == "bac_pathogen" ~ "Bacterial pathogen", 
      Pathogenicity == "fun_pathogen" ~ "Fungal pathogen",
      Pathogenicity == "fun" ~ "Fungi",
      TRUE ~ as.character(Pathogenicity)  # 保留其他值
    ))
  
  # 将Pathogen_Type转换为因子，控制图例顺序
  tmp_node$Pathogenicity <- factor(tmp_node$Pathogenicity, 
                                   levels = c("Bacteria", "Bacterial pathogen", 
                                              "Fungal pathogen", "Fungi"))
  tmp_node$Degree<-tmp_node$igraph.degree
  
  
  #画点
  p1<-p1+ geom_point(aes(X1, X2, fill = Pathogenicity , size = Degree), 
                     pch = 21, data = tmp_node, color = "gray40",alpha = 0.5) +
    scale_size_continuous(range = c(2, 8))
  #调整填充颜色
  p1<-p1+ scale_fill_manual(values = c("Bacteria" = "#FF6B6F",  # 粉色
                                       "Bacterial pathogen" = "#FF0000",  # 红色
                                       "Fungal pathogen" = "#0000FF",  # 深蓝
                                       "Fungi" = "#8D98F6"))   # 浅蓝
  #处理标签 显示物种名称和指定颜色
  p2 <- p1 + 
    geom_text_repel(
      aes(
        X1, X2,
        label = ifelse(
          igraph.degree >= displaySpeciesDegreesPathogen & 
            gr("pathogen", Pathogenicity),
          lab, ""
        ),
        # 注意：这里我们不再映射color，而是通过其他方式
      ),
      family = "Arial",
      fontface = "italic",
      data = tmp_node,
      size = 3,
      max.overlaps = 20
    ) +
    # 重要：移除所有额外的scale_color_manual
    # 只保留必要的主题设置
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "right",
      legend.box = "vertical"
    )
  
  # 如果需要文本颜色，可以在aes外部设置
  p1 <- p1 + 
    geom_text_repel(
      aes(
        X1, X2,
        label = ifelse(
          igraph.degree >= displaySpeciesDegreesPathogen & 
            gr("pathogen", Pathogenicity),
          lab, ""
        )
      ),
      family = "Arial",
      fontface = "italic",
      data = tmp_node,
      size = 3,
      max.overlaps = 20,
      # 根据Pathogenicity设置颜色
      color = ifelse(tmp_node$Pathogenicity == "Bacteria", "#FF6B6F",
                     ifelse(tmp_node$Pathogenicity == "Bacterial pathogen", "#FF0000",
                            ifelse(tmp_node$Pathogenicity == "Fungal pathogen", "#0000FF",
                                   "#8D98F6")))
    )  +
    # 坐标轴和主题
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "right",
      legend.box = "vertical"
    )
  
  
  #保存图形成文件
  ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName," plot.network3.pdf"),p1,width = 12,height = 5)
  ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName," plot.network4.pdf"),p1,width = 16,height = 10)
  p2 = plot[[2]]
  p2<-p2+ggtitle(paste("zipi (",outDir1,")"))+
    theme(plot.title = element_text(size = 6,hjust = 0))
  ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName," plot.zipi.pdf"),p2,width = 12,height = 5)
  p3 = plot[[3]]
  p3<-p3+ggtitle(paste("Random (",outDir1,")"))+
    theme(plot.title = element_text(size = 6,hjust = 0))
  ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName," plot.random.pdf"),p3,width = 12,height = 5)
  
  #将运算节点和边保存
  dat = tab.r[[2]]
  node = dat$net.cor.matrix$node
  # 指定要检查的列
  cols_to_check <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Genus_Species")
  # 将指定列中的 NA 替换为空字符
  node[cols_to_check] <- lapply(node[cols_to_check], function(x) ifelse(is.na(x), "", x))
  # 将指定列中的字符串 "NA" 替换为空字符
  node[cols_to_check] <- lapply(node[cols_to_check], function(x) ifelse(x == "NA", "", x))
  
  # 使用 split() 按照 Group 列拆分数据框
  split_df <- split(node, node$Group)
  # 循环保存每个子数据框为 CSV 文件
  for(group in names(split_df)) {
    FileName <- paste(outDir1,srcType," ",getSampleArg$samplesName,"_",group,"_node","",".csv", sep = "")
    write_csv(split_df[[group]],FileName)
  }
  
  
  edge = dat$net.cor.matrix$edge
  
  node_lab<-as.data.frame(tmp_node$elements,tmp_node$group)
  split_nodes<-split(node_lab,rownames(node_lab))
  split_edges<-split(edge,edge$group)
  for(f in 1:length(split_nodes)){
    temp_n<-data.frame(split_nodes[f])
    temp_n$count<-NA
    temp_e<-data.frame(split_edges[[f]])
    for(i in 1:length(rownames(temp_n))){
      k=0
      for(j in 1:length(rownames(temp_e))){
        if(temp_e$OTU_1[j]==node_lab$`tmp_node$elements`[i]){
          k=k+1
        }
        if(temp_e$OTU_2[j]==node_lab$`tmp_node$elements`[i]){
          k=k+1
        }
      }
      temp_n$count[i]=k
    }
    my_df <- temp_n[temp_n$count != 0, ]
    count_freq <- table(my_df$count)
    count_freq_df <- as.data.frame(count_freq)
    group_name<-rownames(temp_n)[1]
    count_plot<-ggplot(count_freq_df, aes(x = Var1, y = Freq)) +
      geom_line() +                 # 绘制折线
      geom_point() +                # 在折线图上添加数据点
      labs(title = "Frequency of Count Values", # 添加标题
           x = "Count Value",                    # x轴标签
           y = "Frequency")                                              # 点的形状
    ggsave(paste0(outDir1,srcType," ",getSampleArg$samplesName,"_",group_name," degree_count.pdf"),count_plot,width = 12, height = 8)
  }
  
  
  # 创建 Gephi 边文件的数据框
  gephi_edges <- data.frame(
    Source = edge$OTU_1,  # 起始节点
    Target = edge$OTU_2,  # 目标节点
    Weight = edge$weight,  # 权重
    Correlation = edge$cor,  # 相关性（可选）
    Group = edge$Group  # 相关性（可选）
  )
  # 使用 split() 按照 Group 列拆分数据框
  split_df <- split(gephi_edges, gephi_edges$Group)
  # 循环保存每个子数据框为 CSV 文件
  for(group in names(split_df)) {
    FileName <- paste(outDir1,srcType," ",getSampleArg$samplesName,"_",group,"_edge","",".csv", sep = "")
    write_csv(split_df[[group]],FileName)
  }
  # End  Part 2 调用 跨域网络交互网络图 
  ##################################################################>
  
  
  ##################################################################>
  #  Part 3 调用 仅仅显示跨域网络交互图 , 根据需要选择执行
  ##################################################################>
  if (1!=1) { #如不需要改成 1!=1
    #跳过组
    if ( !is.null(skipGroup) ) {
      f_printMsg("剔除不能实现的组名")
      # 获取样本的 Group 信息
      group_info <- sample_data(ps.merge)$Group
      # 选择不在 skipGroup 中的样本（去除 A 和 B 组）
      ps.merge <- prune_samples(!(group_info %in% skipGroup), ps.merge)
    }
    
    #建立存放目录
    outDir2=paste0(outDir,"02仅仅细菌和真菌之间跨域网络/")
    dir.create(outDir2,recursive = T,showWarnings = FALSE)
    
    # 细菌和真菌OTU网络-域网络-二分网络#-------
    # 仅仅关注细菌和真菌之间的相关，不关注细菌内部和真菌内部相关
    #细菌和真菌跨域网络
    #view(ps.merge@otu_table)
    #view(ps.merge@tax_table)
    result <- corBionetwork(ps = ps.merge,
                            N = 0,
                            r.threshold = r.threshold, # 相关阈值
                            p.threshold = p.threshold,
                            group = "Group",
                            # env = data1, # 环境指标表格
                            # envGroup = Gru,# 环境因子分组文件表格
                            # layout = "fruchtermanreingold",
                            path = outDir2,# 结果文件存储路径
                            fill = "filed2", #Phylum", # 出图点填充颜色用什么值
                            size = "igraph.degree", # 出图点大小用什么数据
                            scale = TRUE, # 是否要进行相对丰度标准化
                            bio = TRUE, # 是否做二分网络
                            zipi = F, # 是否计算ZIPI
                            step = 100, # 随机网络抽样的次数
                            width = 12,
                            label = TRUE,
                            height = 10,
                            big = TRUE,
                            select_layout = TRUE,
                            layout_net = "model_maptree2",
                            clu_method = "cluster_fast_greedy")
  }
  # End  Part 3 调用 仅仅显示跨域网络交互图 , 根据需要选择执行
  ##################################################################>
  
}



f_ggClusterNet(getSampleArg,ps,ps_pathogen,N=200,taxFill="Genus", #, #"Species",  #Phylum","Top_Species"
               #displaySpeciesDegrees=1000,displaySpeciesDegreesPathogen=600)    
               displaySpeciesDegrees=20,displaySpeciesDegreesPathogen=10)  