rm(list = ls())
###############################input data 
dir_path <- "/Users/xut2/Desktop/2_anova_protein/"
dir_path_name <- dir(dir_path,pattern = "*.",full.names = T)
dir_path_name
library(openxlsx)
getSheetNames(grep("1-s2.0-S0092867419313856-mmc2.xlsx",dir_path_name,value = T))
data_1 <- read.xlsx(grep("1-s2.0-S0092867419313856-mmc2.xlsx",dir_path_name,value = T),sheet = 2)
colnames(data_1)[427] #[1] "Column1"
data_1 <- data_1[,-grep("Column",colnames(data_1))]
colnames(data_1)
data_1 <- data_1[,-grep("Peptides",colnames(data_1))]
colnames(data_1)
#protein_cell <- colnames(data_1)[-c(1:6)]
#protein_cell <- gsub("TenPx.*", "", protein_cell)
#length(unique(protein_cell)) #[1] 375
#protein_cell[duplicated(protein_cell) | duplicated(protein_cell,fromLast = T)]
#View(protein_cell)
#protein_cell <- gsub("_$","",protein_cell)
#dim(data_1) #[1] 12755   384
#length(unique(data_1$Gene_Symbol)) #[1] 12197
#View(data_1[duplicated(data_1$Gene_Symbol),])
colnames(data_1)[1:10]
data_1$Protein_Id <- NULL
data_1$Description <- NULL
data_1$Group_ID <- data_1$Uniprot <- data_1$Uniprot_Acc <- NULL
#dim(data_1) #[1] 12755   379
data_1 <- data_1[rowSums(is.na(data_1[,-1])) < ncol(data_1)-1, ]
max(rowSums(is.na(data_1[,-1]))) #[1] 369
dim(data_1) #[1] 12755   379
colnames(data_1)[1]
str(data_1)
data_1 <- aggregate(data_1[,-1], list(Gene_Symbol = data_1$Gene_Symbol), FUN = mean, na.rm=TRUE, na.action=NULL)
dim(data_1) #[1] 12196   379
#View(data_1)
row.names(data_1) <- data_1$Gene_Symbol
data_1$Gene_Symbol <- NULL
length(unique(data_1$Gene_Symbol)) #[1] 12196
#data_1 <- data_1[apply(data_1, 1, var) != 0,]
dim(data_1) #[1] 12196   378
colnames(data_1) <- gsub("_TenPx.*", "", colnames(data_1))
#View(data_1)
data_1 <-na.omit(data_1)
data_2 <- data.frame(t(data_1))
dim(data_2) #[1]  378 5118
#View(data_2)
data_2$cancer <- rownames(data_2)
######################inout cancer category
cell_cat <- read.xlsx(grep("Cell line annotations.xlsx",dir_path_name,value = T), sheet = 1)
#View(cell_cat)
#View(cell_cat)
#View(table(cell_cat$type_refined))
colnames(cell_cat)
cell_cat_1 <- cell_cat[,c("CCLE_ID", "type_refined")]
#head(cell_cat_1)
#row.names(cell_cat_1) <- cell_cat_1$CCLE_ID
#cell_cat_2 <- data.frame(t(data.frame(t(cell_cat_1))))
#cell_cat_2$CCLE_ID <- row.names(cell_cat_2)
cell_cat_2 <- cell_cat_1
#######
head(cell_cat_2)
cell_1 <- read.csv(grep("2020-10-05-Cancer cell type_proteomics.csv",dir_path_name,value = T),header = T,stringsAsFactors = F)
head(cell_1)
cell_2 <- cell_1[cell_1$Freq >= 10,]
head(cell_2)
intersect(cell_2$Cancer.cell.type, unique(cell_cat_2$type_refined))
intersect(cell_cat_2$CCLE_ID, unique(colnames(data_1))) #375
setdiff(colnames(data_1),cell_cat_2$CCLE_ID) #character(0)
data_2$cancer[!data_2$cancer %in% cell_cat_2$CCLE_ID]
data_2$cancer <- gsub("^X","",data_2$cancer)
sum(data_2$cancer %in% cell_cat_2$CCLE_ID) #375
colnames(cell_cat_2)[1] <- "cancer"
data_3 <- merge(data_2, cell_cat_2, by  = "cancer")
dim(data_3) #[1]  375 5120
head(cell_cat_2) #cancer          type_refined
table(data_3$type_refined)
data_3$type_refined[!data_3$type_refined %in% cell_2$Cancer.cell.type] <- "other"
write.csv(data_3, paste0(dir_path,Sys.Date(),"-protein_for_anova.csv"),row.names = F)
data_3$cancer <- NULL
dim(data_3) #[1]  375 5119
data_3$type_refined
###############################################################################
data_list <- list()
number_pro <- grep("type_refined",colnames(data_3))-1
number_pro  #5118
for (i in 1:number_pro) {
  #i =1 
  print(i)
  data_4 <- data_3[,c(i, grep("type_refined",colnames(data_3)))]
  data_4[,1] <- as.numeric(data_4[,1])
  data_4$type_refined <- as.factor(data_4$type_refined)
  gene_1 <- colnames(data_4)[1]
  colnames(data_4)[1] <- "Gene"
  oneway <- oneway.test(Gene ~ type_refined, data = data_4, var.equal = TRUE)
  data_5 <- data.frame(gene_1, oneway$p.value)
  data_list[[i]] <- data_5 
}
data_stat <- do.call("rbind",data_list)
head(data_stat)
dim(data_stat)
table(cut(data_stat$oneway.p.value,breaks = seq(0,1,0.05)))
#colnames(data_stat) <- c("Number of gene","FDR")
write.csv(data_stat,paste0(dir_path,Sys.Date(),"-pro_anova_pvalue.csv"),row.names = F)
