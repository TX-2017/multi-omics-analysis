rm(list = ls())
library(tidyverse)
###########################fda data
dir_path <- "C:\\Users\\xut2\\Desktop\\step_3_pathway_analysis_pro_seq\\"
dir_path_files <- list.files(pattern = ".*.csv",dir_path,recursive = T, all.files = F,full.names = T)
dir_path_files
data_0 <- read.csv(dir_path_files[grep("2022-08-26-anova_filter_pro_detail.csv",dir_path_files)],header = T,stringsAsFactors = F)
dim(data_0) #[1] 22015     5
head(data_0)
data_1 <- data_0[,c("cancer_type", "ENTREZID")]
colnames(data_0) <- gsub("ENTREZID", "entrez_id", colnames(data_0))
data_gene <- data_0
dim(data_gene) #[1] 12921     2
head(data_gene)
table(data_gene$cancer_type)
#1.2 input pathway_all
pathway <- read.csv(dir_path_files[grep("2020-09-15-pathway_remove_biop.csv",dir_path_files)],header = T,stringsAsFactors = F)
dim(pathway) #[1] 108639      2
grep("bioplanet",pathway$pathway) #integer(0)
length(unique(pathway$pathway)) #[1] 2270
head(pathway)
pathway <- pathway[,c(1,2)]
pathway <- unique(pathway)
pathway$label <- 1
pathway_1 <- spread(pathway, key = pathway, value = label, fill = 0)
dim(pathway_1) #[1] 11729  2271
pathway <- pathway_1
library(venn)
library(ggplot2)
for (k in 1:length(unique(data_gene$cancer_type))) {
  #k= 16
  print(k)
  #k= 1
  data_2 <- data_gene[data_gene$cancer_type == unique(data_gene$cancer_type)[k], c("entrez_id","cancer_type")]
 
  colnames(data_2) <- c("GeneID","cancer_type")
  #head(data_2)
  dir.create(paste0(dir_path,unique(data_gene$cancer_type)[k],"/"))
  dir_path_1 <- paste0(dir_path,unique(data_gene$cancer_type)[k],"/")
  ##############################################################pre n3 and n4
  #head(data_2);colnames(data_2) #[1] "GeneID"      "liver_cancer"
  out_dir <- dir_path_1 
  venn::venn(list("Specific genes" = unique(na.omit(data_2$GeneID)),
                  "Pathway related genes" = unique(na.omit(pathway$GeneID))),
             color = c("darkblue", "darkred"),box = F, size = 0.5,
             ggplot = T) +
    annotate(geom="text", x= 500, y= 900, label= data_2$cancer_type[1],
             color="black") 
  
  ggsave(filename = paste0(Sys.Date(),"-",data_2$cancer_type[1],".tif"), 
         plot = last_plot(), 
         device = "tiff", path = out_dir,
         scale = 1, compression = "lzw", width = 8, height = 8, units = "cm",
         dpi = 300, limitsize = TRUE)
###########################################
  data_all_1 <- data_all_2 <- list()
  dat_1 <- data_2[which(data_2$cancer_type == unique(data_2$cancer_type)[1]),]
  dat_5 <- unique(na.omit(dat_1))
  data_3 <- merge(dat_5, pathway,by = "GeneID", all.y = T)
  #dim(data_3) #[1] 11729  2272
  data_3$cancer_type[!is.na(data_3$cancer_type)] <- 1
  data_3$cancer_type[is.na(data_3$cancer_type)] <- 0
  data_3$GeneID <- NULL
  data_3 <- data_3[, c(2:ncol(data_3),1)]
  #View(data_3)
  a <- data_3
  str(data_3)
  ######################################################
  b <- NULL
  for (mm in 1:(ncol(a)-1)) {
    #mm=1
    #dim(a)
    #print(mm)
    tox_taget <- a[which(a[,mm] == 1 & a[,ncol(a)] == 1), ]
    tox_non_taget <- a[which(a[,mm] == 0 & a[,ncol(a)] == 1), ]
    non_tox_target <- a[which(a[,mm] == 1 & a[,ncol(a)] == 0), ]
    non_tox_non_taget <- a[which(a[,mm] == 0 & a[,ncol(a)] == 0), ]
    #dim(tox_taget) #[1]  76 367
    #dim(tox_non_taget) #[1] 132 367
    #dim(non_tox_target) #[1] 164 367
    #dim(non_tox_non_taget) #[1] 380 367
    #nrow(tox_taget)
    ##################################fisher
    tox_target <- matrix(c(nrow(tox_taget), nrow(non_tox_target), nrow(tox_non_taget), nrow(non_tox_non_taget)),nrow = 2,
                         dimnames = list(tox = c("tox", "non_tox"),gene = c("tar", "nontar")))
    tox_target_pvalue <- fisher.test(tox_target, alternative = "two.sided")$p.value
    b <- c(b,tox_target_pvalue)
  }
  e <- data.frame(b,colnames(a)[1:(ncol(a)-1)])
  colnames(e) <- c("pvalue","Pathway")
  ######################
  e2 <- e
  #head(e2)
  #head(f)
  f <- e2[order(e2$pvalue),]
  #min(f$pvalue)
  #head(f)
  f$fdr <- p.adjust(f$pvalue,method = "fdr")
  #head(f)
  f$cancer_type <-unique(data_2$cancer_type)
  write.csv(f,paste0(dir_path_1,Sys.Date(),"-",unique(data_2$cancer_type)[1],"-pathway_original",".csv"),row.names = F)
  
}
