
rm(list = ls())
library(scGATE)
################################################################################
# load data
ClusterI  <- as.data.frame(read.csv(paste0("\\write_fig2\\ClusterI.csv"))[ ,2:15])
ClusterII <- as.data.frame(read.csv(paste0("\\write_fig2\\ClusterII.csv"))[ ,2:15])
data      <- ClusterII
################################################################################
# data normalization
q975      <- quantile(c(as.matrix(data)), probs = 0.975 , na.rm = FALSE)
data      <- data/q975
data      <- replace(data, data>=1, 0.9999)
################################################################################
# remove genes with low variability
data      <- data[ , which(sqrt(apply(data,2,var))> 0.20)]
################################################################################
# base grn definition
gene_list <- colnames(data)
base_grn  <- c()
for(i in 1:length(gene_list)){
  for(j in 1:length(gene_list)){
    if(i!=j){
      base_grn  <- rbind(base_grn, c(gene_list[j], gene_list[i]))
    }
  }
}
################################################################################
# Please note that the likelihood values can be affected by the tSNE clustering results and the number of EM iterations.

h_set     <- c(1, 1.25, 1.75, 2.5, 5, 7) # the range of h values in the Hill climbing function
gates     <- scGATE_logic(data = data, number_of_em_iterations = 5, max_num_regulators = 3, h_set = h_set, base_grn = base_grn, abs_cor = 0.05, top_gates = 1, grid_par = "simple")
print(gates)
################################################################################































