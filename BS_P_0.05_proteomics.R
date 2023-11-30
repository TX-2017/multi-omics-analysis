rm(list = ls())
###############################input data 
dir_path <- "/Users/xut2/Desktop/pathway_add_bootstrap/pro_pathway/"
cancer_list <- list.dirs(path = dir_path, full.names = F, recursive = F)
dir_path_fisher <- dir(paste0(dir_path),pattern = ".*.csv",full.names = T, recursive = T)
dir_path_fisher #
###############################bath input data 
dir.create(paste0(dir_path,"4_final_results/"))
output_1 <- paste0(dir_path,"4_final_results/")
stat_list <- list()
for (i in 1:length(cancer_list)) {
  #i =1
  data_1 <- read.csv(grep(paste0(cancer_list[i],"-pathway_original"),dir_path_fisher,value = T,fixed = T),header = T,stringsAsFactors = F)
  #dim(data_1)
  #head(data_1)
  data_1 <- data_1[data_1$pvalue < 0.05, ]
  colnames(data_1)[c(1,3)] <- paste0(colnames(data_1)[c(1,3)], "_ori")
  data_2 <- read.csv(grep(paste0(cancer_list[i], "-pathway_add_bootstrap"),dir_path_fisher,value = T,fixed = T), header = T,stringsAsFactors = F)
  colnames(data_2)[1:2] <- c("pvalue-boot","Pathway")
  head(data_2,2 )
  head(data_1, 2)
  #colnames(data_1);colnames(data_2)
  print(c(dim(data_1),dim(data_2)))
  data_12 <- merge(data_1, data_2, by = "Pathway", all = T)
  print(head(data_12,6))
  #print(i)
  #data_12 <- data_12[which(data_12$`pvalue-ori` < 0.05),]
  #dat_0.05 <- length(unique(data_12$endpoint)) #[1] 43
  data_12$difference <- data_12$`pvalue-boot`- data_12$pvalue_ori
  data_12$label[data_12$difference < 0] <- 1
  data_12$label[data_12$difference >= 0] <- 0
  #data_13 <- data_12[which(data_12$difference > 0),]
  #data_13$label_1 <- 1
  data_13 <- data_12
  head(data_13)
  str(data_13)
  #table(data_13$label)
  data_14 <- aggregate(data_13["label"],by = list(Pathway = data_13$Pathway,
                                                         pvalue_ori = data_13$pvalue_ori,
                                                         cancer_type = data_13$cancer_type,
                                                         fdr_ori = data_13$fdr_ori),
                                                         FUN = sum, na.rm = T)
  head(data_14)
  #dim(data_14)
  #unique(data_13$source)
  #data_14$label_1 <- 10 - data_14$label_1
  data_15 <- data_14[which(data_14$label < 50),]
  #View(data_15)
  #dim(data_15)
  #table(data_15$label_1)
  data_15 <- data_15[data_15$fdr_ori < 0.05, ]
  write.csv(data_15,paste0(output_1,Sys.Date(),"-BS_filter_",cancer_list[i],".csv"),row.names = F)
  #write.csv(data_13,paste0(output_1,Sys.Date(),"-detail-",unique(dir_path_name_1)[i]),row.names = F)
  stat_list[[i]] <- data.frame(cancer_list[i],length(unique(data_15$Pathway)))
}
data_stat <- do.call("rbind",stat_list)
colnames(data_stat) <- c("source","bs_pvalue_0.05_num")
write.csv(data_stat,paste0(output_1,Sys.Date(),"-stat_p_0.05_BS_pro.csv"),row.names = F)
#View(data_15)


