# Context-specific network and logic gate inference in the mouse haematopoiesis scRNA-seq data
# 1. Please refer to the Jupyter notebook for instructions on how to perform Louvain clustering on the cells in the mouse haematopoiesis scRNA-seq dataset.
# 2. Retrieve the data from Megakaryocyte cells (Cluster 11).
# Load scGATE package and data in example_data folder

rm(list = ls())
library(scGATE)
data         <- as.data.frame(read.csv("/example_data/subset_counts_cluster_11.csv" , header = TRUE))

# select genes involved in the MegE differentiation
gene_list   <- c("Gata1", "Fli1", "Klf1", "Spi1", "Zfpm1", "Tal1", "Gata2")
data        <- data[  , gene_list]
data        <- na.omit(data)
print(head(data))

# Load base GRN
base_GRN    <- read.csv(file = "/example_data/base_grn_mouse_blood_cell_differentiation_toggle_switch.csv")

# 3. data preprocessing 
# The dataset underwent library size normalization in Jupyter Notebook. To fit the scRNA-seq data within the (0,1) interval, we applied quantile normalization as a technique to rescale the data.
data        <- scRNA_seq_preprocessing(data = data, library_size_normalization = "False")

# 4. Run scGATE_logic() function
gates       <- scGATE_logic(data = data, base_GRN = base_GRN, number_of_em_iterations = 10, top_gates = 1, run_mode = "slow")
print(head(gates))
