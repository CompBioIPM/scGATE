rm(list = ls())

library(readr)
library(ROCR)
################################################################################
setwd(paste0("D:\\scGATE\\scGATE\\R"))
load("base_GRNs.RData")
################################################################################
tissue_of_interest    <- "Spleen-10X_P4_7"

# load data
tissue         <- strsplit(tissue_of_interest, "-")[[1]][1]
target_names   <- unlist(read.table(paste0("D:\\scRNA logic\\mouse spleen dataset\\BEELINE\\data_cor\\", tissue_of_interest, "\\","tg_filtered.txt")))
tf_names       <- unlist(read.table(paste0("D:\\scRNA logic\\mouse spleen dataset\\BEELINE\\data_cor\\", tissue_of_interest, "\\","tf_filtered.txt")))

data           <- as.data.frame(read.csv(paste0("D:\\scRNA logic\\mouse spleen dataset\\BEELINE\\data_cor\\" , tissue_of_interest, "\\", "ExpressionData.csv") , header = TRUE))
gene_names     <- data[ ,1]
data           <- t(data[ ,2:ncol(data)])
colnames(data) <- gene_names

# make sure tfs and target genes are present in the data
tf_filtered    <- intersect(tf_names, colnames(data))
tg_filtered    <- intersect(target_names, colnames(data))

# sort expression profile matrix with tfs in the first columns and targets in the next columns
data           <- data[ ,c(tf_filtered, tg_filtered)]

# remove cells with low read counts
data           <- data[apply(data > 0, 1, sum) > 0.1*ncol(data), ]

# data normalization
data           <- t(apply(data, 1, function(x) x/sum(x)))*100
data           <- na.omit(data)
q975           <- quantile(as.numeric(data), probs = 0.975 , na.rm = FALSE)
data           <- data/q975
data           <- replace(data, data>=1, 0.9999)
data[1:10,1:10]
################################################################################
# load base grns for each tissue
if(tissue=="Spleen"){
  base_grn     <- atac_candidates_peak[[1]]
}else if(tissue=="Lung"){
  base_grn     <- atac_candidates_peak[[2]]
}else if(tissue=="Liver"){
  base_grn     <- atac_candidates_peak[[3]]
}else if(tissue=="Kidney"){
  base_grn     <- atac_candidates_peak[[4]]
}else if(tissue=="Heart_and_Aorta"){
  base_grn     <- atac_candidates_peak[[5]]
}
############################################################################ identify TF-gene network
source("scGATE_edge.R")
abs_cor        <- 0
k_act          <- 0.7
h_act          <- 7
res            <-  scGATE_edge(data = data, number_of_em_iterations = 3, max_num_regulators = 3, k_act = k_act, h_act = h_act, tf_names=tf_filtered, tg_names=tg_filtered, abs_cor, base_grn)
##############################################################################
source("label_calc_beeline.R")
TF_TARGET             <- as.data.frame(read.csv(paste0("D:\\scRNA logic\\mouse spleen dataset\\BEELINE\\data_cor\\" , tissue_of_interest, "\\", "refNetwork.csv"), header = TRUE))
TF_TARGET             <- TF_TARGET[TF_TARGET$Gene1 != TF_TARGET$Gene2, ]

ranked_edge_list      <- res$ranked_edge_list
ranked_edge_list$label<- 0

ranked_edge_list      <- label_calc(tf_numbers=length(tf_filtered), tg_numbers=length(tg_filtered), TF_TARGET, ranked_edge_list)
ranked_list           <- ranked_edge_list$BF_score
labels                <- ranked_edge_list$label
pred                  <- prediction(ranked_list, labels)
auc                   <- performance(pred, measure = "auc")@y.values[[1]]
##############################################################################
ranked_edge_list$label<- as.numeric(ranked_edge_list$label)
edge_numbers          <- sum(ranked_edge_list$label)
EPR                   <- sum(ranked_edge_list[1:edge_numbers, 4])
##############################################################################
print(c(tissue_of_interest, auc, EPR, EPR/edge_numbers))




