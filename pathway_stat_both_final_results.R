rm(list = ls())
###############################input data 
dir_path <- "C:\\Users\\xut2\\Desktop\\step_3_rnaseq_pro_pathway\\"
dir_path_name <- list.files(pattern = "4_final_results",dir_path,full.names = T, recursive = T, include.dirs = T)
list_name_1 <- list.files(pattern = ".*csv", dir_path_name[1], full.names = F, recursive = F)
list_name_2 <- gsub(".*filter_","", list_name_1)
list_name_3 <- gsub(".csv","",list_name_2)[1:c(length(list_name_2)-1)]
cancer_type <- list.files(pattern = ".*csv", dir_path_name, full.names = T, recursive = F)
##################################
data_list_1 <- data_list_2 <- data_list_rna <- data_list_pro <- list()
for (i in 1:length(list_name_3)) {
  data_files <- grep(list_name_3[i],cancer_type,value = T)
  data_rna <- read.csv(grep("RNA",data_files,value = T), header = T, stringsAsFactors = F)
  data_pro <- read.csv(grep("PRO",data_files,value = T), header = T, stringsAsFactors = F)
  colnames(data_rna)[-1] <- paste0(colnames(data_rna)[-1], "_RNA")
  colnames(data_pro)[-1] <- paste0(colnames(data_pro)[-1], "_PRO")
  data_merge <- merge(data_rna, data_pro, by = "Pathway")
  data_stat <- data.frame(rna = nrow(data_rna), pro = nrow(data_pro), both = nrow(data_merge))
  data_stat$cancer_type <- list_name_3[i]
  data_list_rna[[i]] <- data_rna
  data_list_pro[[i]] <- data_pro
  data_list_1[[i]] <- data_merge
  data_list_2[[i]] <- data_stat
}
data_stat_all <- do.call("rbind", data_list_2)
data_merge_pathway <- do.call("rbind", data_list_1)
write.csv(data_stat_all, paste0(dir_path,Sys.Date(),"-","pathway_stat_both.csv"),row.names = FALSE,na = "")
write.csv(data_merge_pathway, paste0(dir_path,Sys.Date(),"-","pathway_detail_merge.csv"),row.names = FALSE,na = "")
##########################################
data_rna_pathway <- do.call("rbind", data_list_rna)
data_pro_pathway <- do.call("rbind", data_list_pro)
write.csv(data_rna_pathway, paste0(dir_path,Sys.Date(),"-","pathway_rna_boot_all.csv"),row.names = FALSE,na = "")
write.csv(data_pro_pathway, paste0(dir_path,Sys.Date(),"-","pathway_pro_boot_all.csv"),row.names = FALSE,na = "")
