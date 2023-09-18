
rm(list = ls())
library(scGATE)
################################################################################
# load data
gene_list <- c( "Gata1", "Fli1", "Klf1", "Spi1", "Zfpm1", "Tal1", "Gata2")

# MegE progenitor cells
data_7    <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_7.csv") , header = TRUE))

# Meg. cells (Cluster 11)
data_11   <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_11.csv") , header = TRUE))

# Early and Late Erythroid cells
data_2    <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_2.csv") , header = TRUE))
data_3    <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_3.csv") , header = TRUE))
data_6    <- as.data.frame(read.csv(paste0("\\write_fig3\\", "subset_counts_cluster_6.csv") , header = TRUE))
################################################################################
for(cluster in c(2,3,7,11)){
  if(cluster==2){
    data  <- data_2
  }else if(cluster==3){
    data  <- data_3
  }else if(cluster==7){
    data  <- data_7
  }else if(cluster==11){
    data  <- data_11
  }
  data    <- data[ ,2:ncol(data)]
  # select genes involved in the MegE differentiation
  data    <- data[  , gene_list]
  data    <- na.omit(data)
  ################################################################################
  # data normalization
  q975    <- quantile(as.numeric(as.matrix(data)), probs = 0.975 , na.rm = FALSE)
  data    <- data/q975
  data    <- replace(data, data>=1, 0.9999)
  ################################################################################
  # base_grn definition
  base_grn<- rbind( c("Fli1", "Gata1"), c("Spi1", "Gata1"), c("Gata2", "Gata1"))
  ################################################################################
  h_set   <- c(1, 1.25, 1.75, 2.5, 5, 7) # the range of h values in the Hill climbing function
  gates   <- scGATE_logic(data = data, number_of_em_iterations = 10, max_num_regulators = 3, h_set = h_set, base_grn = base_grn, abs_cor = 0, top_gates = 20, grid_par = "complex")
  print(gates)
} # end for cluster
################################################################################
















