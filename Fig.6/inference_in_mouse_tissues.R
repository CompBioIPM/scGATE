
rm(list = ls())
library(scGATE)
################################################################################
tissue_of_interest    <- "Spleen-10X_P4_7"

# load data
load("base_GRNs_mouse_tissues.RData")
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

tissue         <- strsplit(tissue_of_interest, "-")[[1]][1]
target_names   <- unlist(read.table(paste0("\\", tissue_of_interest, "\\","tg_filtered.txt")))
tf_names       <- unlist(read.table(paste0("\\", tissue_of_interest, "\\","tf_filtered.txt")))

data           <- as.data.frame(read.csv(paste0("\\" , tissue_of_interest, "\\", "ExpressionData.csv") , header = TRUE))
gene_names     <- data[ ,1]
data           <- t(data[ ,2:ncol(data)])
colnames(data) <- gene_names

# make sure tfs and target genes are present in the expression data
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
################################################################################ identify TF-gene network
abs_cor        <- 0
k_act          <- 0.7
h_act          <- 7
res            <-  scGATE_edge(data = data, number_of_em_iterations = 3, max_num_regulators = 3, k_act = k_act, h_act = h_act, tf_names=tf_filtered, tg_names=tg_filtered, base_grn, abs_cor)
print(res$ranked_edge_list)
################################################################################
