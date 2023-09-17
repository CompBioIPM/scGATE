label_calc      <- function(tf_numbers, tg_numbers, TF_TARGET, ranked_edge_list){
  # tf_numbers  <- 4
  # tg_numbers  <- 200
  for(gene_id in 1:tg_numbers){
    # gene_id =1
    # k       =1
    input_edge_tmp    <- TF_TARGET[TF_TARGET$Gene2==colnames(data)[tf_numbers+gene_id], ]
    if(length(input_edge_tmp)>0){
      for(k in 1:nrow(input_edge_tmp)){
        string_vec    <- input_edge_tmp[k, ]
        matching_rows <- which(apply(ranked_edge_list[,c(2,1)], 1, function(x) all(x %in% string_vec)))
        ranked_edge_list[matching_rows, 4] <- 1
        #ranked_edge_list[matching_rows,  ]
      }
    }
  } # end for gene_id
  return(ranked_edge_list)
}

# colnames(data)[1:6]
# input_edge_tmp    <- TF_TARGET[TF_TARGET$V2=="Alcam", ]
# input_edge_tmp