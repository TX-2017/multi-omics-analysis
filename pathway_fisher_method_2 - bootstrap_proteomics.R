rm(list = ls())
library(tidyverse)
###########################fda data
dir_path <- "/Users/xut2/Desktop/20220903_pathway_add_bootstrap/pro_pathway/"
dir_path_files <- list.files(pattern = ".*.csv",dir_path,recursive = T, all.files = F,full.names = T)
dir_path_files
data_0 <- read.csv(dir_path_files[grep("2022-08-26-anova_filter_pro_detail.csv",dir_path_files)],header = T,stringsAsFactors = F)
data_1 <- data_0[,c("cancer_type", "ENTREZID")]
colnames(data_0) <- gsub("ENTREZID", "entrez_id", colnames(data_0))
data_gene <- data_0
pathway <- read.csv(dir_path_files[grep("2020-09-15-pathway_remove_biop.csv",dir_path_files)],header = T,stringsAsFactors = F)
pathway <- pathway[,c(1,2)]
pathway <- unique(pathway)
pathway$label <- 1
pathway_1 <- spread(pathway, key = pathway, value = label, fill = 0)
dim(pathway_1) #[1] 11729  2271
pathway <- pathway_1
library(venn)
#dir.create(paste0(dir_path,"1.0_sig_pathway/"))
library(ggplot2)
for (k in 1:length(unique(data_gene$cancer_type))) {
  data_2 <- data_gene[data_gene$cancer_type == unique(data_gene$cancer_type)[k], c("entrez_id","cancer_type")]
  colnames(data_2) <- c("GeneID","cancer_type")
  dir_path_1 <- paste0(dir_path,unique(data_gene$cancer_type)[k],"/")
  ##############################################################pre n3 and n4
  #head(data_2);colnames(data_2) #[1] "GeneID"      "liver_cancer"
  out_dir <- dir_path_1 
###########################################
  data_all_1 <- data_all_2 <- list()
  dat_1 <- data_2[which(data_2$cancer_type == unique(data_2$cancer_type)[1]),]
  #dim(dat_1) #[1] 409   2
  #View(dat_1)
  dat_5 <- unique(na.omit(dat_1))
  #head(dat_1)
  data_3 <- merge(dat_5, pathway,by = "GeneID", all.y = T)
  #dim(data_3) #[1] 11729  2272
  data_3$cancer_type[!is.na(data_3$cancer_type)] <- 1
  data_3$cancer_type[is.na(data_3$cancer_type)] <- 0
  data_3$GeneID <- NULL
  data_3 <- data_3[, c(2:ncol(data_3),1)]
  dim(data_3) #[1] 11729  2271
  #####################################p value 0.05
  pvalue <- read.csv(dir_path_files[grep(unique(data_gene$cancer_type)[k],dir_path_files)],header = T,stringsAsFactors = F)
  #dim(pvalue) #[1] 2270    4
  #head(pvalue)
  pvalue <- pvalue[pvalue$pvalue <0.05, ]
  data_3 <- data_3[, c(pvalue$Pathway, colnames(data_3)[ncol(data_3)])]
  #########################################################
    #View(data_3)
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
      #dim(a)
      b <- NULL
      #memory.limit(999999999)
      ####################
        for (mm in 1:(ncol(a)-1)) {

          dd <- spread(data.frame(table(a[,mm], a[,ncol(a)])),key = Var1, value = Freq)
          b <- c(b,fisher.test(dd[,-1], alternative = "two.sided")$p.value)
          #print(b)
        }
      ###################
      e <- data.frame(b,colnames(a)[1:(ncol(a)-1)])
      colnames(e) <- c("pvalue",colnames(a)[ncol(a)])
      e$label <- j
      object_file_lisT_select[[j]] <- e
    }
    data_pvalue <- do.call("rbind",object_file_lisT_select)
     write.csv(data_pvalue,paste0(dir_path_1,Sys.Date(),"-",unique(data_2$cancer_type)[1],"-pathway_add_bootstrap",".csv"),row.names = F)
  
}
