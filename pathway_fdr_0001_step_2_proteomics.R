rm(list = ls())
###############################input data 
dir_path <- "C:\\Users\\liyix\\OneDrive\\Desktop\\2022\\2022_CCLE\\20220903_pathway_fisher\\20220903_pathway_fisher\\pro_pathway\\"
file_list <- list.files(dir_path, pattern = ".*.pathway_without_gene_add_fdr.csv",full.names = T, recursive = T)
file_list
##input data
all_data_stat <- data_fdr_0001 <- all_data <- list()
for (i in 1:length(file_list)) {
  #i = 1
  data_1 <- read.csv(file_list[i],header = T,stringsAsFactors = F)
  dim(data_1)
  print(dim(data_1))
  data_stat <- data.frame(table(cut(data_1$fdr, breaks = c(0,0.001,0.01,0.05,1),include.lowest = T)))
  data_stat$cancer_category <- gsub("\\/.*","",gsub(".*\\\\/","",file_list[i]))
  all_data_stat[[i]] <-  data_stat
  #print(data_stat)
  #View(data_1)
  all_data[[i]] <- data_1
  data_2 <- data_1[data_1$fdr < 0.05, ]
  print(dim(data_2)) #[1] 101   8
  data_fdr_0001 [[i]] <- data_2
}
data_stat_all <- do.call("rbind", all_data_stat)
data_all <- do.call("rbind", data_fdr_0001)
length(unique(data_all$cancer_type)) #[1] 23
length(unique(data_all$pathway)) #[1] 383
max(data_all$fdr) #[1] 0.000990652
length(data_all$pathway) #[1] 1227
data_3 <- do.call("rbind", all_data)
write.csv(data_stat_all, paste0(dir_path,Sys.Date(),"-","pathway_stat.csv"),row.names = FALSE)
write.csv(data_all, paste0(dir_path,Sys.Date(),"-","pathway_fdr_005.csv"),row.names = FALSE)
write.csv(data_3, paste0(dir_path,Sys.Date(),"-","pathway_all.csv"),row.names = FALSE)

