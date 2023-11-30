rm(list = ls())
###############################input data 
dir_path <- "C:\\Users\\liyix\\OneDrive\\Desktop\\2022\\2022_CCLE\\20220803_protein\\2_anova\\"
dir_path_name1 <- list.files(pattern ="2022-08-26-pro_anova_pvalue.csv",dir_path,full.names = T, recursive = T)
dir_path_name1
data_avova <- read.csv(dir_path_name1,header = T,stringsAsFactors = F)
dim(data_avova) #[1] 5118    2
#table(data_avova$cancer_type)
max(data_avova$oneway.p.value) # 0.04999903
data_avova <- data_avova[data_avova$oneway.p.value < 0.01, ]
dim(data_avova) #[1] 4890    2
############################################################
dir_path_name <- list.files(pattern = "2020-11-08-gene_for_pathway_add_id_0.05",dir_path,full.names = T, recursive = T)
dir_path_name
data_1 <- read.csv(dir_path_name[1],header = T,stringsAsFactors = F)
head(data_1,2)
dim(data_1) #[1] 23022     5
max(data_1$pvalue)
data_2 <- data_1[data_1$SYMBOL %in% data_avova$gene_1, ]
dim(data_2) #[1] 22015     5
table(data_2$cancer_type)
head(data_2)
##########################################################
data_3 <- data.frame(table(data_2$cancer_type))
dir_path_name2 <-list.files(pattern ="2021-01-24-pro_rna_gene_count.csv",dir_path,full.names = T, recursive = T)
dir_path_name2
data_type <- read.csv(dir_path_name2,header = T,stringsAsFactors = F)
data_type$cancer_type
data_3$Var1 <- as.character(data_3$Var1)
data_3$Var1 <- tolower(data_3$Var1)
tolower(data_type$cancer_type)[!tolower(data_type$cancer_type) %in% data_3$Var1]
############################
data_4 <- data_3[data_3$Var1 %in% tolower(data_type$cancer_type), ]
dim(data_4) #[1] 16  2
write.csv(data_4, paste0(dir_path,Sys.Date(),"-anova_filter_pro_stat.csv"),row.names = FALSE,na = "")
#################################
table(data_2$cancer_type)
data_5 <- data_2[tolower(data_2$cancer_type) %in% tolower(data_type$cancer_type), ]
unique(data_5$cancer_type)
write.csv(data_5, paste0(dir_path,Sys.Date(),"-anova_filter_pro_detail.csv"),row.names = FALSE,na = "")
table(data_5$cancer_type)

