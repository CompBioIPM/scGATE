# 1. Load scGATE package and data
rm(list=ls())
library(scGATE)

data           <- as.data.frame(read.csv(/Fig.4/net1_data/net1-3000-1/ExpressionData.csv", header = TRUE))
gene_names     <- data[ ,1]
data           <- t(data[ ,2:ncol(data)])
colnames(data) <- gene_names
head(data[ , 1:10])

# Load TF list
# This step is optional
tf_names       <- paste0("T", 1:15)
print(head(tf_names))

# 2. data preprocessing 
# For scGATE simulated data, library size normalization is not performed. 
# However, the simulated data is only re-scaled using the quantile normalization technique to fit within the (0,1) interval.
data           <- scRNA_seq_preprocessing(data = data, library_size_normalization = "False", tf_list = tf_names)

# 3. Run scGATE_edge() function
ranked_edge_list <-  scGATE_edge(data = data, base_GRN = NA, h_act = 2, number_of_em_iterations = 20, abs_cor = 0.05)
print(head(ranked_edge_list))
