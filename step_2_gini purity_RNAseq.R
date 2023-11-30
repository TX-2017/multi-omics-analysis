rm(list = ls())
library(mltools)
################# ##############input data 
dir_path <- "C:\\Users\\xut2\\Desktop\\RNA_seq data\\2_ttest\\"
dirs_list <- list.dirs(dir_path,full.names = T, recursive = F)
data_list_all <- data_list_all_type <- list()
for (i in 1:length(dirs_list)) {
  files_list <- dir(dirs_list[i],pattern = "*.data_for_plot_labs_position",full.names = T)
  data_list <- daat_type <- list()
  for (j in 1:length(files_list)) {
    data_1 <- read.csv(files_list[j],header = T,stringsAsFactors = F)
    data_2 <- data_1[data_1$family == "Specical", ]
    daat_type[[j]] <-  data.frame(nrow(data_2),gsub("-data_for_plot_labs_position.csv","",gsub(".*complete-","",files_list[j])))
    ifelse (length(seq(min(data_2$x),max(data_2$x),1)) - length(sort(data_2$x)) > length(sort(data_2$x)), 
            gini_purity <- gini_impurity(seq(min(data_2$x),max(data_2$x),1) %in% sort(data_2$x)),
            gini_purity <- 1 - gini_impurity(seq(min(data_2$x),max(data_2$x),1) %in% sort(data_2$x)))
    
    data_3 <- data.frame(gini_purity, data_2$label[1],gsub("-data_for_plot_labs_position.csv","",gsub(".*complete-","",files_list[j])))
    colnames(data_3) <- c("Gini purity", "Cancer type", "ttest with p value")
    data_list[[j]] <-  data_3
  }
  data_4 <- do.call("rbind",data_list)
  data_list_all[[i]] <- data_4
  data_list_all_type[[i]] <- do.call("rbind",daat_type)
  
}
data_5 <- do.call("rbind", data_list_all)
data_5_type <- do.call("rbind", data_list_all_type)
colnames(data_5_type) <- c("number of cell lines", "Cancer type")
data_5$pvale <- data_5$`ttest with p value`
data_5$pvale <- sub("-[[:alpha:]].*", "", data_5$pvale)
data_5$`Cancer type_1` <- sub(".*[[:digit:]]-", "", data_5$`ttest with p value`)
data_5$`Cancer type_1` <- sub("-data.*", "", data_5$`Cancer type_1`)
data_5$pvale <- as.numeric(data_5$pvale)
data_6 <- data_5
data_6 <- data_5[order(data_5$`Cancer type_1`, -data_5$pvale),]
write.csv(data_6, paste0(dir_path,Sys.Date(),"-discribution_each_subject_Gini_purity_rnaseq_only_log2.csv"),row.names = F)
write.csv(data_5_type, paste0(dir_path,Sys.Date(),"-discribution_each_subject_Gini_purity_stat_rnaseq_only_log2.csv"),row.names = F)
#########################
library(ggplot2)
head(data_6)
ggplot(data=data_6 , aes(x=pvale, y=`Gini purity`, group=`Cancer type_1`)) +
  geom_line(linetype="dashed", size=0.5,color= "black")+
  geom_point(color= "blue", size=1.5) +
  labs(title="RNA-SEQ",x ="FDR", y = "Gini purity")+
  theme(legend.position = "", 
        legend.background = element_rect(fill="gray95", size=.5, linetype="dotted"),
        panel.background = element_rect(fill = NA, colour = "black",size = 1,linetype = "solid"),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x =  element_line(colour = "grey90",linetype="dotted"),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(color="black", size= 10, face="bold.italic", angle=0),
        strip.text.y = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.2, "lines"),
        axis.title = element_text(color="black", size=10, face="bold"),
        axis.text.y = element_text(color="black", size=10, face="plain"),
        plot.title = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(color="black", size=10, face="plain", angle = 90, hjust = 1,vjust = 0.3),
        #plot.title = element_text(hjust = 0.5,color="red", size= 10, face="bold.italic"),
        axis.ticks.length = unit(.15, "cm"),
        axis.ticks.x= element_line(size = 1),
        axis.ticks.y= element_line(size = 1)) + 
  scale_y_continuous(expand = c(0.01,0.01),limits = c(-0.1,1.1),breaks = seq(0,1,0.2)) +
  facet_wrap( `Cancer type_1` ~., ncol=4,scales = "free_y" )+
  geom_hline(yintercept= 1, colour="gray50", linetype="dashed") +
  geom_point(data=data_6[data_6$`Gini purity` == "1",],color="red",size=2) +
  scale_x_log10(breaks = unique(data_6$pvale))
#output pic
ggsave(filename = paste0(Sys.Date(),"-discribution_each_subject_Gini purity_rnaseq_only_log2.tif"), plot = last_plot(), 
       device = "tiff", path = dir_path,
       scale = 1, width = 20, height = 18, units = "cm",dpi = 300, 
       limitsize = TRUE, compression = "lzw")


