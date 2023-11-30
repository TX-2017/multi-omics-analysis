rm(list = ls())
library(ggdendro)
library(ggplot2)
################# ##############input data 
dir_path <- "C:\\Users\\xut2\\Desktop\\step_2_proteomics\\"
dir_path_name <- dir(dir_path,pattern = "*.",full.names = T )
dir_path_name
library(openxlsx)
getSheetNames(grep("1-s2.0-S0092867419313856-mmc2.xlsx",dir_path_name,value = T))
data_1 <- read.xlsx(grep("1-s2.0-S0092867419313856-mmc2.xlsx",dir_path_name,value = T),sheet = 2)
data_t <- data_1
data_1 <- data_1[,-grep("Column",colnames(data_1))]
data_1 <- data_1[,-grep("Peptides",colnames(data_1))]
protein_cell <- colnames(data_1)[-c(1:6)]
protein_cell <- gsub("TenPx.*", "", protein_cell)
protein_cell <- gsub("_$","",protein_cell)
######################inout cancer category
cell_cat <- read.xlsx(grep("Cell line annotations.xlsx",dir_path_name,value = T), sheet = 1)
data_2 <- cell_cat[cell_cat$CCLE_ID %in% intersect(unique(cell_cat$CCLE_ID),as.character(protein_cell)),]
write.csv(data_2, paste0(dir_path,Sys.Date(),"-","375_cells_overlapping_anno.csv"),row.names = FALSE)
data_stat <- data.frame(table(data_2$type_refined))
data_stat <- data_stat[order(data_stat$Freq, decreasing = T),]
colnames(data_stat)[1] <- "Cancer cell type"
write.csv(data_stat, paste0(dir_path,Sys.Date(),"-","Cancer cell type_proteomics.csv"),row.names = FALSE)
###################################################t-test
data_t <- data_t[,-grep("Column",colnames(data_t))]
colnames(data_t)
data_t <- data_t[,-grep("Peptides",colnames(data_t))]
data_t$Protein_Id <- NULL
data_t$Description <- NULL
data_t$Group_ID <- data_t$Uniprot <- data_t$Uniprot_Acc <- NULL
#dim(data_t) #[1] 12755   379
data_t <- data_t[rowSums(is.na(data_t[,-1])) < ncol(data_t)-1, ]
#head(data_t)
data_t <- aggregate(data_t[,-1], list(Gene_Symbol = data_t$Gene_Symbol), FUN = mean, na.rm=TRUE, na.action=NULL)
dim(data_t) #[1] 12196   379
#View(data_t)
row.names(data_t) <- data_t$Gene_Symbol
data_t$Gene_Symbol <- NULL
length(unique(data_t$Gene_Symbol)) #[1] 12196
#data_t <- data_t[apply(data_t, 1, var) != 0,]
dim(data_t) #[1] 12196   378
colnames(data_t) <- gsub("_TenPx.*", "", colnames(data_t))
######################inout cancer category
cell_cat <- read.xlsx(grep("Cell line annotations.xlsx",dir_path_name,value = T), sheet = 1)
cell_cat_1 <- cell_cat[,c("CCLE_ID", "type_refined")]
cell_cat_2 <- cell_cat_1
cell_1 <- data_stat
head(cell_1)
cell_2 <- cell_1[cell_1$Freq >= 10,]
###############################################################################
num_i <- seq(6,8,1)
seq_y <-  c(1e-10*10^num_i,seq(0.01,0.05,0.01),seq(0.1,1,0.1))
length(seq_y)
data_t <- na.omit(data_t)
data_2 <- data_t
#data_2$pvalue <- NULL
dim(data_2) #[1] 12196   378
dim(data_t) #[1] 12196   378
for (i in 1:length(cell_2$Cancer.cell.type)) {
  #prepare data
  #i = 10
  data_t$pvalue <- "pvalue"
  dir.create(paste0(dir_path,cell_2$Cancer.cell.type[i]))
  out_dir <- paste0(dir_path,cell_2$Cancer.cell.type[i],"\\")
  #View(cell_cat_2)
  type_1 <- cell_cat_2[cell_cat_2$type_refined == cell_2$Cancer.cell.type[i],]$CCLE_ID
  type_1 <- as.character(na.omit(type_1))
  data_special <- data_2[,colnames(data_2) %in% type_1]
  dim(data_special) #[1] 12196    14
  data_other <- data_2[,!colnames(data_2) %in% type_1]
  dim(data_other) #[1] 12196   364
  intersect(colnames(data_special),colnames(data_other))
  #dim(data_2) #[1] 56202   128
  ####################################ttest
  for (j in 1:nrow(data_2)) {
    print(j)
    #j =3
    tryCatch ({
      data_t$pvalue[j] <- t.test(as.numeric(na.omit(as.numeric(data_special[j,]))),as.numeric(na.omit(as.numeric(data_other[j,]))))$p.value
    }, error = function(e) {
      cat("ERROR :",conditionMessage(e), "\n")
    })
  }
  #str(data_t$pvalue)
  data_t <- data_t[data_t$pvalue != "pvalue",]
  data_t$pvalue <- as.numeric(data_t$pvalue)
  data_t <- data_t[!is.na(data_t$pvalue),]
  dim(data_t) #[1] 10813   379
  #View(data_t[1:10,])
  data_t_fdr <- data_t
  data_t_fdr$fdr <- p.adjust(data_t_fdr$pvalue, method = "fdr", n = length(data_t_fdr$pvalue))
  #View(head(data_t_fdr))
  write.csv(data_t_fdr[,c("pvalue","fdr")], paste0(out_dir,Sys.Date(),"-",cell_2$Cancer.cell.type[i],"_pvalue_FDR.csv"),row.names = T)
  method_1 <- "complete"
  data_list <- list()
  for (k in 1:length(seq_y)) {
    #k = 1
    tryCatch ({
      print(k)
      data_mad_1 <- data_t_fdr[as.numeric(data_t_fdr$fdr) < seq_y[k],-c(ncol(data_t_fdr),ncol(data_t_fdr)-1)]
      dim(data_mad_1) #[1] 1124 1019
      #View(head(data_mad_1))
      data_out <- data_t_fdr[as.numeric(data_t_fdr$fdr) < seq_y[k],c(ncol(data_t_fdr),ncol(data_t_fdr)-1)]
      #write.csv(data_out, paste0(out_dir,Sys.Date(),"-",method_1[1],"-",seq_y[k],"-data_for_gene_list.csv"),row.names = T)
      data_list[[k]] <- data.frame(nrow(data_out))
      data_list[[k]]$source <- seq_y[k]
      data_3 <- data.frame(t(data_mad_1))
      dim(data_3) #[1] 378   6
      #str(data_3)
      #View(data_3)
      data_3 <- na.omit(data_3)
      ##########################################cluster_method and label
      data_d <- dist(scale(data_3))
      hc1 <- hclust(data_d,method_1[1])
      fam <- data.frame(row.names(data_3),row.names(data_3))
      str(fam)
      #View(fam)
      fam[] <- sapply(fam[], as.character)
      #dim(fam) #[1] 1019    2
      #head(fam)
      colnames(fam) <- c("label","family")
      fam$family[!fam$family %in% type_1] <- "Non_specical"
      fam$family[fam$family %in% type_1] <- "Specical"
      ###############################
      dend <- as.dendrogram(hc1)
      dendata <- dendro_data(dend,type = "rectangle")
      #add column family to labs dataframe (look like VLOOKUP in excel)
      labs <- dendata$labels
      labs <- merge(labs, fam, by = "label")
      ###########################################plot
      data_plot <- segment(dendata)
      data_plot$label <- "ALL"
      ###########################################
      data_plot_1 <- data_plot[data_plot$yend == 0,]
      dim(data_plot_1) #[1] 1019    5
      label_2 <- labs[grep("Specical",labs$family),]
      colnames(label_2) <- gsub("y","yend",colnames(label_2))
      data_plot_1$label <- NULL
      data_plot_2 <- merge(data_plot_1,label_2,by = c("x", "yend"),all.x = T)
      #View(data_plot_2)
      data_plot_2$familyend[is.na(data_plot_2$familyend)] <- "Non_specical"
      data_plot_2$label <- NULL
      colnames(data_plot_2) <- gsub("familyend", "label",colnames(data_plot_2))
      data_plot_2$y <- 0
      data_plot_2$yend <- -5
      #################################
      head(data_plot,2);head(data_plot_2,2)
      data_p2 <- rbind(data_plot,data_plot_2)
      #################################
      head(data_plot)
      head(data_p2)
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
      
      ggsave(filename = paste0(Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-ttest_circle.pdf"), 
             plot = p0, device = "pdf", path = out_dir,
             width = 4, height = 4, units = "cm",
             dpi = 300, limitsize = TRUE)
      write.csv(data_p2,paste0(out_dir,Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-data_for_plot.csv"),row.names = F)
      write.csv(labs,paste0(out_dir,Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-data_for_plot_labs_position.csv"),row.names = F)
      
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
      
      ggsave(filename = paste0(Sys.Date(),"-",method_1[1],"-",seq_y[k],"-",cell_2$Cancer.cell.type[i],"-ttest_rectangle.pdf"), 
             plot = p2, device = "pdf", path = out_dir,
             width = 4, height = 4, units = "cm",
             dpi = 300, limitsize = TRUE)
    }, error = function(e) {
      cat("ERROR :",conditionMessage(e), "\n")
    })
  }
  
  data_stat <- do.call("rbind",data_list)
  colnames(data_stat) <- c("Number of gene","FDR")
  write.csv(data_stat,paste0(out_dir,Sys.Date(),"-",method_1[1],"-",cell_2$Cancer.cell.type[i],"-data_for_gene_stat.csv"),row.names = F)
}
########################t-test

