
rm(list = ls())
library(scGATE)
################################################################################
# load data
gene_list   <- c( "Gata1", "Fli1", "Klf1", "Spi1", "Zfpm1", "Tal1", "Gata2")

# MegE progenitor cells
data_7      <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_7.csv") , header = TRUE))

# Meg. cells (Cluster 11)
data_11     <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_11.csv") , header = TRUE))

# Early and Late Erythroid cells
data_2      <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_2.csv") , header = TRUE))
data_3      <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_3.csv") , header = TRUE))
data_6      <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_6.csv") , header = TRUE))

data_MegE   <- rbind(data_2, data_3, data_6, data_11)
data_MegE   <- data_MegE[ ,2:ncol(data_MegE)]

# select genes involved in the MegE differentiation
data        <- data_MegE[  , gene_list]
data        <- na.omit(data)
################################################################################
# data normalization
q975        <- quantile(as.numeric(as.matrix(data)), probs = 0.975 , na.rm = FALSE)
data        <- data/q975
data        <- replace(data, data>=1, 0.9999)
################################################################################
# base_grn definition
base_grn    <- rbind(c("Gata1", "Klf1"), c("Fli1", "Klf1"))
base_grn    <- rbind(base_grn, c("Gata1", "Fli1"), c("Klf1", "Fli1"))
base_grn    <- rbind(base_grn, c("Gata1", "Zfpm1"))
base_grn    <- rbind(base_grn, c("Gata1", "Tal1"), c("Spi1", "Tal1"))
base_grn    <- rbind(base_grn, c("Zfpm1", "Gata2"), c("Spi1", "Gata2"))
base_grn    <- rbind(base_grn, c("Gata1", "Gata2"))
################################################################################
h_set       <- c(1, 1.25, 1.75, 2.5, 5, 7) # the range of h values in the Hill climbing function
gates       <- scGATE_logic(data = data, number_of_em_iterations = 10, max_num_regulators = 3, h_set = h_set, base_grn = base_grn, abs_cor = 0, top_gates = 20, grid_par = "complex")
print(gates)
################################################################################



