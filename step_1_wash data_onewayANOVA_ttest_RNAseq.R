rm(list = ls())
library(ggplot2)
library(openxlsx)
###############################input data 
dir_path <- "C:\\Users\\xut2\\Desktop\\RNA_seq data\\"
dir_path_name <- dir(dir_path,pattern = "*.",full.names = T)
data_1 <- read.table(grep("NAseq_genes_rpkm_20180929.gct",dir_path_name,value = T),header = T,skip = 2,stringsAsFactors = F)
data_11 <- data_1
data_11[, -c(1, 2)] <- sapply(data_11[, -c(1, 2)], as.numeric)
###############################remove outliers using the capping method
data_12 <- data.frame(apply(data_1[,-c(1,2)], 2, max))
colnames(data_12) <- "max"
p1 <- ggplot(data_12) +
  stat_density(aes(x = max),size = 0.6,
               alpha=0.5, bw = 5,  geom="line",position="identity", color = "blue") 
p <- ggplot_build(p1)
p$data[[1]]$x[which.max(p$data[[1]]$y)] #[1] 10324.26
p1 + geom_vline(xintercept = p$data[[1]]$x[which.max(p$data[[1]]$y)])
data_11[, -c(1, 2)][data_11[, -c(1, 2)] >= 10324.26] <- 10324.26
#############################log2 transformation 
data_list <- list()
for (i in 1:nrow(data_11)) {
  print(i)
    data_2 <- data_11[i, ]
    data_21 <- log2(as.numeric(data_2[, -c(1,2)]))
    data_21[data_21 == -Inf] <- -max(abs(data_21[data_21 != -Inf]))
    data_2[, -c(1,2)] <- data_21
    data_list[[i]] <- data_2
}
data_3 <- do.call("rbind", data_list)
data_ori <- data_3[apply(data_3[,-c(1,2)], 1, sum) != Inf, ]
write.csv(data_ori, paste0(dir_path,Sys.Date(),"-","rna_seq_cap_log2.csv"),row.names = FALSE,na = "")
################################one-way ANOVA
data_31 <- data_ori
row.names(data_31) <- data_31$Name
data_31$Name <- NULL
data_31$Description <- NULL
data_32 <- data.frame(table(sub(".*?_", "", colnames(data_31))))
colnames(data_32) <- c("Cancer cell type","Number of cell lines")
data_32 <- data_32[order(data_32$`Number of cell lines`, decreasing = T),]
row.names(data_32) <- 1:nrow(data_32)
write.csv(data_32, paste0(dir_path,Sys.Date(),"-","Cancer cell type_STAT.csv"),row.names = T)
type_1 <- as.character(data_32$`Cancer cell type`[data_32$`Number of cell lines` > 20])
colnames(data_31)[!sub(".*?_", "", colnames(data_31)) %in% type_1] <- "other"
data_31[nrow(data_31)+1, ] <- sub(".*?_", "", colnames(data_31))
data_3 <- data.frame(t(data_31))
data_list <- list()
for (i in 1:length(grep("EN",colnames(data_3)))) {
  data_4 <- data_3[,c(grep("EN",colnames(data_3))[i], grep("X55415",colnames(data_3)))]
  data_4[,1] <- as.numeric(data_4[,1])
  data_4$X55415 <- as.factor(data_4$X55415)
  gene_1 <- colnames(data_4)[1]
  colnames(data_4)[1] <- "Gene"
  oneway <- oneway.test(Gene ~ X55415, data = data_4, var.equal = TRUE)
  data_5 <- data.frame(gene_1, oneway$p.value)
  data_list[[i]] <- data_5 
}
data_stat <- do.call("rbind",data_list)
dir.create(paste0(dir_path,"anova_step_1\\"))
write.csv(data_stat,paste0(dir_path, "anova_step_1\\", Sys.Date(),"-rna_seq_cap_log2_each_gene_anova.csv"),row.names = F)
################################t-test
data_33 <- data_ori
row.names(data_33) <- data_33$Name
data_33$Name <- NULL
data_33$Description <- NULL
data_33 <- data_33[apply(data_33, 1, var) != 0,]
cell_cat <- read.xlsx(grep("Cell line annotations.xlsx",dir_path_name,value = T), sheet = 1)
cell_cat_1 <- cell_cat[,c("CCLE_ID", "type_refined")]
row.names(cell_cat_1) <- cell_cat_1$CCLE_ID
cell_cat_2 <- data.frame(t(data.frame(t(cell_cat_1))))
cell_cat_2$CCLE_ID <- row.names(cell_cat_2)
cell_1 <- read.csv(grep("Cancer cell type.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
cell_2 <- cell_1[cell_1$Freq >= 15,]
num_i <- seq(0,8,1)
seq_y <-  c(1e-10*10^num_i,0.05)
data_2 <- data_33
dir.create(paste0(dir_path, "2_ttest\\"))
dir_path <- paste0(dir_path, "2_ttest\\")
for (i in 1:length(cell_2$Cancer.cell.type)) {
  data_33$pvalue <- "pvalue"
  dir.create(paste0(dir_path,cell_2$Cancer.cell.type[i],"-1"))
  out_dir <- paste0(dir_path,cell_2$Cancer.cell.type[i],"-1\\")
  type_1 <- cell_cat_2[cell_cat_2$type_refined == cell_2$Cancer.cell.type[i],]$CCLE_ID
  type_1 <- as.character(na.omit(type_1))
  data_special <- data_2[,colnames(data_2) %in% type_1]
  data_other <- data_2[,!colnames(data_2) %in% type_1]
  ####################################ttest
  for (j in 1:nrow(data_2)) {
    tryCatch({
      data_33$pvalue[j] <- t.test(as.numeric(data_special[j,]), 
                                  as.numeric(data_other[j,]))$p.value
    }, error = function(e) {
      cat("ERROR :",conditionMessage(e), "\n")
      cat("ERROR :", conditionMessage(e),"---",i,"---",gsub("\\:","-",Sys.time()),file = "error.txt", append = TRUE, "\n")
    })
  }
  data_33_fdr <- data_33
  data_33_fdr$fdr <- p.adjust(data_33_fdr$pvalue, method = "fdr", n = length(data_33_fdr$pvalue))
  write.csv(data_33_fdr[,c("pvalue","fdr")], paste0(out_dir,Sys.Date(),"-",cell_2$Cancer.cell.type[i],"_pvalue_FDR_rnaseq_only_log2.csv"),row.names = T)
  ###################################################
  method_1 <- "complete"
  data_list <- list()
  for (k in 1:length(seq_y)) {
    data_mad_1 <- data_33_fdr[as.numeric(data_33_fdr$fdr) < seq_y[k],-c(ncol(data_33_fdr),ncol(data_33_fdr)-1)]
    data_out <- data_33_fdr[as.numeric(data_33_fdr$fdr) < seq_y[k],c(ncol(data_33_fdr),ncol(data_33_fdr)-1)]
    data_list[[k]] <- data.frame(nrow(data_out))
    data_list[[k]]$source <- seq_y[k]
    data_3 <- data.frame(t(data_mad_1))
    ##########################################cluster_method and label
    data_d <- dist(scale(data_3))
    hc1 <- hclust(data_d,method_1[1])
    fam <- data.frame(row.names(data_3),row.names(data_3))
    fam[] <- sapply(fam[], as.character)
    colnames(fam) <- c("label","family")
    fam$family[!fam$family %in% type_1] <- "Non_specical"
    fam$family[fam$family %in% type_1] <- "Specical"
    dend <- as.dendrogram(hc1)
    dendata <- dendro_data(dend,type = "rectangle")
    labs <- dendata$labels
    labs <- merge(labs, fam, by = "label")
    data_plot <- segment(dendata)
    data_plot$label <- "ALL"
    data_plot_1 <- data_plot[data_plot$yend == 0,]
    label_2 <- labs[grep("Specical",labs$family),]
    colnames(label_2) <- gsub("y","yend",colnames(label_2))
    data_plot_1$label <- NULL
    data_plot_2 <- merge(data_plot_1,label_2,by = c("x", "yend"),all.x = T)
    data_plot_2$familyend[is.na(data_plot_2$familyend)] <- "Non_specical"
    data_plot_2$label <- NULL
    colnames(data_plot_2) <- gsub("familyend", "label",colnames(data_plot_2))
    data_plot_2$y <- 0
    data_plot_2$yend <- -5
    data_p2 <- rbind(data_plot,data_plot_2)
    p <- ggplot(data_p2) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend, colour = label),size = 0.01) + 
      coord_flip() + scale_y_reverse() + 
      scale_colour_manual(values = c("#E5E5E5","darkgreen","red"))+
      coord_polar(theta="x")
    
    p0 <- p + theme(legend.position = "")+ 
      theme(axis.line =element_blank(),
            axis.ticks = element_blank(),
            axis.text=element_blank(),
            axis.title =element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank()) +
      labs(x = "", y = "", title = "")
    
    ggsave(filename = paste0(Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-ttest_circle_log2.pdf"), 
           plot = p0, device = "pdf", path = out_dir,
           width = 4, height = 4, units = "cm",
           dpi = 300, limitsize = TRUE)
    write.csv(data_p2,paste0(out_dir,Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-data_for_plot_rnaseq_only_log2.csv"),row.names = F)
    write.csv(labs,paste0(out_dir,Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-data_for_plot_labs_position_rnaseq_only_log2.csv"),row.names = F)
    
    p1 <- ggplot(data_p2) +
      geom_segment(aes(x=x, y=y, xend=xend, yend=yend, colour = label),size = 0.01) + 
      coord_flip() + scale_y_reverse() + 
      scale_colour_manual(values = c("#E5E5E5","darkgreen","red"))
    
    p2 <- p1 + theme(legend.position = "")+ 
      theme(axis.line =element_blank(),
            axis.ticks = element_blank(),
            axis.text=element_blank(),
            axis.title =element_blank(),
            panel.background=element_rect(fill="white"),
            panel.grid=element_blank()) +
      labs(x = "", y = "", title = "")
    
    ggsave(filename = paste0(Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-ttest_rectangle_log2.pdf"), 
           plot = p2, device = "pdf", path = out_dir,
           width = 4, height = 4, units = "cm",
           dpi = 300, limitsize = TRUE)
    
  }
  data_stat <- do.call("rbind",data_list)
  colnames(data_stat) <- c("Number of gene","FDR")
  write.csv(data_stat,paste0(out_dir,Sys.Date(),"-",method_1[1],"-",cell_2$Cancer.cell.type[i],"-data_for_gene_stat_rnaseq_only_log2.csv"),row.names = F)
}



