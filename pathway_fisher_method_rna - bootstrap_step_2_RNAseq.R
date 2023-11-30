rm(list = ls())
library(venn)
library(ggplot2)
library(ggpolypath)
library(tidyverse)
###########################fda data
dir_path <- "C:\\Users\\xut2\\Desktop\\step_3_pathway_analysis\\"
dir_path_files <- list.files(pattern = ".*.csv",dir_path,recursive = T, all.files = F,full.names = T)
dir_path_files
data_1 <- read.csv(dir_path_files[grep("2023-11-02-rna_list_details.csv",dir_path_files)],header = T,stringsAsFactors = F)
data_1 <- data_1[,c("cancer_type", "entrez_id")]
data_gene <- data_1
#####################################################
pathway <- read.csv(dir_path_files[grep("2020-09-15-pathway_remove_biop.csv",dir_path_files)],header = T,stringsAsFactors = F)
pathway <- pathway[,c(1,2)]
pathway <- unique(pathway)
pathway$label <- 1
pathway_1 <- spread(pathway, key = pathway, value = label, fill = 0)
pathway <- pathway_1
for (k in 1:length(unique(data_1$cancer_type))) {
data_2 <- data_gene[data_gene$cancer_type == unique(data_gene$cancer_type)[k], c("entrez_id","cancer_type")]
colnames(data_2) <- c("GeneID","cancer_type")
dir_path_1 <- paste0(dir_path,unique(data_gene$cancer_type)[k],"\\")
out_dir <- dir_path_1 
###########################################
data_all_1 <- data_all_2 <- list()
dat_1 <- data_2[which(data_2$cancer_type == unique(data_2$cancer_type)[1]),]
dat_5 <- unique(na.omit(dat_1))
data_3 <- merge(dat_5, pathway,by = "GeneID", all.y = T)
data_3$cancer_type[!is.na(data_3$cancer_type)] <- 1
data_3$cancer_type[is.na(data_3$cancer_type)] <- 0
data_3$GeneID <- NULL
data_3 <- data_3[, c(2:ncol(data_3),1)]
#####################################p value 0.05
pvalue <- read.csv(dir_path_files[grep(unique(data_gene$cancer_type)[k],dir_path_files)][2],header = T,stringsAsFactors = F)
pvalue <- pvalue[pvalue$pvalue <0.05, ]
#dim(pvalue) #[1] 190   4
head(pvalue)
#pvalue$Pathway %in% colnames(data_3)
#colnames(data_3)[ncol(data_3)]
data_3 <- data_3[, c(pvalue$Pathway, colnames(data_3)[ncol(data_3)])]
#dim(data_3) #[1] 11729   191
#########################################################
#View(head(data_3))
object_file_lisT_select <- list()
#k = 100 ##########################times of boottrap
for (j in 1:1000) {
  #j = 1
  print(paste0(k,"-",j))
  set.seed(2020*j)
  data_4 <- data_3
  data_4[, ncol(data_4)] <- sample(data_4[,ncol(data_4)],replace = F, size = nrow(data_4))
  #dim(data_321);table(data_321$endpoint)
  a <- data_4
  b <- NULL
  system.time(
  for (mm in 1:(ncol(a)-1)) {
    dd <- spread(data.frame(table(a[,mm], a[,ncol(a)])),key = Var1, value = Freq)
    b <- c(b,fisher.test(dd[,-1], alternative = "two.sided")$p.value)
  })
  e <- data.frame(b,colnames(a)[1:(ncol(a)-1)])
  colnames(e) <- c("pvalue",colnames(a)[ncol(a)])
  e$label <- j
  object_file_lisT_select[[j]] <- e
}
data_pvalue <- do.call("rbind",object_file_lisT_select)
print(head(data_pvalue))
write.csv(data_pvalue,paste0(dir_path_1,Sys.Date(),"-",unique(data_2$cancer_type)[1],"-pathway_add_bootstrap",".csv"),row.names = F)
}

