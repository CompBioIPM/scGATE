# Context specific network inference in mouse tissue scRNA-seq datasets
# 1. Please refer to the Jupyter notebook for instructions on how to perform scATAC-seq analysis to derive the candidate TF lists (base GRNs) in *.parquet file format.
# 2. Load scGATE package and data (base GRN and scRNA-seq data and TF list) in example_data folder 

rm(list=ls())
library(scGATE)
# Load base GRN derived from external hints
candidate_tf_target <- as.data.frame(read_parquet("/example_data/Cusanovich2018_Spleen_peak_base_GRN_dataframe.parquet"))
candidate_tf_target <- read_base_GRN(candidate_tf_target)

# Load scRNA-seq data
data           <- as.data.frame(read.csv(paste0("/example_data/Tabula_Muris2018_Spleen-10X_P4_7_ExpressionData.csv") , header = TRUE))
gene_names     <- data[ ,1]
data           <- t(data[ ,2:ncol(data)])
colnames(data) <- gene_names

head(data[ , 1:10])

# Load TF list
# This step is optional
tf_names       <- unlist(read.table("/example_data/Tabula_Muris2018_Spleen-10X_P4_7_tf_lists.txt"))

# 3. scRNA-seq data preprocessing (library size normalization, quantile normalization technique to fit the scRNA-seq data within the (0,1) interval) 
data           <- scRNA_seq_preprocessing(data = data, library_size_normalization = "True", tf_list = tf_names)

# 4. Run scGATE_edge() function
ranked_edge_list <-  scGATE_edge(data = data, base_GRN = candidate_tf_target, h_act = 7)
print(head(ranked_edge_list))
