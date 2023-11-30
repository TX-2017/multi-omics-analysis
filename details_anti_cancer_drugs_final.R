rm(list = ls())
###############################input data 
dir_path <- "C:\\Users\\xut2\\Desktop\\step_5_pathway_drugs\\"
dir_path_name <- list.files(pattern = ".*csv",dir_path,full.names = T, recursive = T)
####################################rna_pro_gene_list___________1
data_rna <- read.csv(grep("2023-11-02-rna_list_details.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
data_rna$Name <- NULL
data_pro <- read.csv(grep("2022-08-26-anova_filter_pro_detail.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
data_pro$SYMBOL <- data_pro$pvalue <- data_pro$GENENAME <- NULL
data_rna <- data_rna[, c("cancer_type", "entrez_id")]
colnames(data_rna)[2] <- "ENTREZID"
rna_pro <- rbind(data_rna, data_pro)
rna_pro <- unique(rna_pro)
#########################################pathway_gene_rna_pro_____________2
pathway_ori <- read.csv(grep("2020-09-15-pathway_remove_biop.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
pathway_rna_pro <- read.csv(grep("2023-11-08-pathway_detail_merge.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
pathway_rna_pro <- pathway_rna_pro[, c("Pathway", "cancer_type_RNA")]
colnames(pathway_rna_pro)[1] <- "pathway"
pathway_gene_rna_pro <- merge(pathway_rna_pro, pathway_ori, by = "pathway")
#############################################cancer_type_pathwa_gene______3
colnames(pathway_gene_rna_pro)[2:3] <- c("cancer_type","ENTREZID")
data_1 <- merge(pathway_gene_rna_pro,rna_pro, by =  c("cancer_type","ENTREZID"))
################################################multiple_pathways_for_gene_________4
data_1 <- unique(data_1)
data_2 <- data.frame(table(data_1$ENTREZID, data_1$cancer_type))
colnames(data_2) <- c("GeneID","cancer_type","number of involved pathways")
data_3 <- data_2[data_2$`number of involved pathways` > 1, ]
data_stat <- data.frame(table(data_3$cancer_type))
colnames(data_stat) <- c("cancer_type","number of involved multiple pathways")
################################cancer_target_drugs_based_on genes number of involved multiple pathways______5
data_drug <- read.csv(grep("2021-04-08-compound_gene_add_ann.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
colnames(data_3)[1] <- "ENTREZID"
data_4 <- merge(data_3, data_drug, by = "ENTREZID")
data_stat_drug <- data.frame(table(data_4$cancer_type))
colnames(data_stat_drug) <- c("cancer_type","number of involved drugs")
#######################################################add__mapping id_________8
data_mapping <- read.csv(grep("2021-09-01-sample_name_mapping.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
data_4 <- merge(data_mapping, data_4, by = "Sample.Name")
#############################################final_drugs_stat______6
data_list_g <- list()
for (i in 1:length(unique(data_4$cancer_type))) {
  data_m <- data_4[data_4$cancer_type == unique(data_4$cancer_type)[i], ]
  data_list_1 <- list()
  for (j in 1:10) {
    data_m1 <- data_m[data_m$`number of involved pathways` >= j, ]
    data_s1 <- data.frame(length(unique(data_m1$Mapping.ID)))
    data_s1$cancer_type <- unique(data_4$cancer_type)[i]
    data_s1$number_target_gene <- j
    data_list_1[[j]] <- data_s1
  }
  data_list_g[[i]] <- do.call("rbind", data_list_1)
}
data_stat <- do.call("rbind", data_list_g)
colnames(data_stat)[1] <- "number_of_compound"
library(tidyr)
data_spread <- spread(data_stat, key = number_target_gene, value = number_of_compound, -1)
write.csv(data_spread, paste0(dir_path,Sys.Date(),"-","stat_number_drugs.csv"),row.names = FALSE)
#########################################cancer_drug_list_stat_______________7
data_stat_1 <- data_stat[data_stat$number_of_compound < 100, ]
df2 <- lapply(split(data_stat_1, data_stat_1$cancer_type),
              function(x) x[which.max(x$number_of_compound),])
df1 <- do.call('rbind', df2)
#######################################
data_list_g <- list()
for (i in 1:length(unique(data_4$cancer_type))) {
  data_m <- data_4[data_4$cancer_type == unique(data_4$cancer_type)[i], ]
  df3 <- df1[df1$cancer_type == unique(data_4$cancer_type)[i], ]
  #dim(data_m)
  data_m1 <- data_m[data_m$`number of involved pathways` >= df3$number_target_gene, ]
  data_list_g[[i]] <- data_m1
}
data_drugs <- do.call("rbind", data_list_g)
write.csv(data_drugs, paste0(dir_path,Sys.Date(),"-","details_anti_cancer_drugs_final.csv"),row.names = FALSE)
data_final <- merge(data_1, data_drugs, by = c("cancer_type", "ENTREZID"))
write.csv(data_final, paste0(dir_path,Sys.Date(),"-","details_anti_cancer_drugs_final.csv"),row.names = FALSE)