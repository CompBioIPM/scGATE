# Context-specific network and logic gate inference in synthetic toggle switch
# 1. Please refer to the Jupyter notebook for instructions on how to perform Louvain clustering on the cells in the BoolODE simulated data.
# 2. Retrieve the data from Cluster I of cells, which was obtained in the previous step.
# Load scGATE package and data in example_data folder
 
rm(list = ls())
library(scGATE)
data         <- as.matrix(read.csv(paste0("/example_data/ClusterI.csv"))[ ,2:15])

# 3. data preprocessing 
# For scGATE simulated data, library size normalization is not performed. 
# However, the simulated data is only re-scaled using the quantile normalization technique to fit the data within the (0,1) interval.
data         <- scRNA_seq_preprocessing(data = data, library_size_normalization = "False")

# 4. Remove genes with low variability (scGATE operates on highly variable genes per context).
# This step is optional
data$n_counts<- data$n_counts[ , which(sqrt(apply(data$n_counts,2,var))> 0.20)]

# 5. Run scGATE_logic() function
# Please note that the likelihood values can be affected by the Louvain clustering results.
gates        <- scGATE_logic(data = data, top_gates = 1, run_mode = "fast")
print(head(gates))
