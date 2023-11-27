## scGATE: single-cell gene regulatory gate inference
We have developed scGATE (single-cell gene regulatory gate) as a logic-based model for deciphering tissue- and cell-type-specific gene regulatory networks from single-cell RNA sequencing (scRNA-seq). While previous efforts have focused on reconstructing directed transcription factor (TF) to target gene networks, logic-based models enable the exploration of more complex combinatorial relationships between regulators. In particular, Boolean logic models can capture higher-order TF interactions, represented through AND, OR, and XOR logic gates. Our novel approach infers TF-gene networks from scRNA-seq data while simultaneously elucidating the underlying Boolean logic that combines TF activities.
The methodology integrates external genomic data, such as single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq) and motif analysis, to narrow down candidate TFs for each target gene. TF-gene links are further refined using scRNA-seq data, and logic rules are derived. To enhance statistical power, scGATE focuses on reconstructing context-specific regulatory networks, including tissue- or cell-type-specific networks. This study presents an integrative framework for deducing cell-type-specific gene regulation, moving beyond TF-gene pairs to capture the complex logic operations underlying combinatorial control.

<br>
Malekpour, S.A., Haghverdi, L., Sadeghi, M., Single-cell multi-omics analysis identifies context specific gene regulatory gates and mechanisms. Briefings In Bioinformatics, 2023.

##### If you find our study useful and relevant to your research, please kindly cite us. Your citation means a lot to us and helps acknowledge our contributions.
<br>
<br>

![Fig1](https://github.com/CompBioIPM/scGATE/assets/47293318/4be29239-0cd7-4871-aa01-18268ba7bb7d)

## Step 1. scGATE installation
The scGATE codes are written in R version 4.1.3 and have been tested in both Windows and Linux environments. 

### Installation
1. Download the compiled package file `scGATE_0.1.0.tar.gz` from this GitHub page.
2. Install the scGATE package by running the following command in R:
   
   ```R
   install.packages("path/to/scGATE_0.1.0.tar.gz", repos = NULL, type = "source")
   ```
<br>

### Dependencies  
Please ensure that you have the following packages installed:

```R
install.packages("VGAM")  
install.packages("truncnorm")
install.packages("arrow")
```
These commands will install the VGAM, truncnorm, and arrow packages, which are required for running scGATE.  

To load the packages, use the following commands:  
```R
library(scGATE)  
library(VGAM)  
library(truncnorm)  
library(arrow) 
```

<br>

## Step 2. Prepare input files
### Preprocessing base GRN generated from external hints
To summarize information in the base GRN file in ".parquet" format, previously generated using external hints like scATAC-seq and TF binding motif analyses, you can use the `read_base_GRN()` function from the scGATE package.

```R
# Read and summarize base GRN file
candidate_tf_target <- as.data.frame(read_parquet("res_Buenrostro2018_base_GRN_dataframe.parquet"))
candidate_tf_target <- read_base_GRN(candidate_tf_target)
```

### Preprocessing scRNA-seq count data
To preprocess raw scRNA-seq data, including steps such as normalization and rescaling, you can use the scRNA_seq_preprocessing() function from the scGATE package.

```R
# Preprocess scRNA-seq count data
normalized_counts <- scRNA_seq_preprocessing(data = data_scRNA_seq, library_size_normalization = "True", tf_list = NA)
```

Parameter Descriptions  
data: The scRNA-seq raw data matrix with cells in rows and genes in columns.  
library_size_normalization: A flag indicating whether library size normalization should be performed. The default value is "True". Set it to "True" if you want to perform library size normalization.  
tf_list: A list of transcription factors (TFs) to consider. The default value is NA, which means all columns in the data matrix will be considered as TFs.  

<br>

## Step 3. Run scGATE 
scGATE provides two functions for TF-target network inference: `scGATE_edge()` and `scGATE_gate()`. These functions infer the TF-target network without and with predicted Boolean logic gates in the output, respectively.
<br>
### TF-Target Network Inference (edge mode)
To infer the TF-target network without logic gates in the output, you can use the `scGATE_edge()` function.

```R
# Infer TF-target network without logic gates in the output
res <- scGATE_edge(n_counts = normalized_counts, base_GRN = candidate_tf_target, k_act = 0.7, h_act = 10, number_of_em_iterations = 3, max_num_regulators = 3, abs_cor = 0)
print(res$ranked_edge_list)
```

Parameter Descriptions  
n_counts: The normalized and scaled scRNA-seq count data from the previous preprocessing step.  
base_GRN: The TF-target gene network inferred from previous steps using external hints. Leave it empty if no base GRN is available.  
k_act and h_act: Hill function parameters used in the inference process.  
number_of_em_iterations: The number of iterations in the expectation-maximization (EM) algorithm.  
max_num_regulators: The maximum number of TFs in a Boolean logic gate. In the main manuscript, a value of 3 is used.  
abs_cor: This parameter varies in the (0, 1) interval and further removes edges with low absolute Pearson correlations between TFs and their targets. A value of 0 indicates no filtration based on correlations.  

<br>

## Example usage of scGATE 
### Context-specific network and logic gate inference in synthetic toggle switch 

```R
# 1. Please refer to the Jupyter notebook for instructions on how to perform Louvain clustering on the cells in the BoolODE simulated data.
# 2. Retrieve the data from Cluster I of cells, which was obtained in the previous step.
# Load scGATE package and data in example_data folder
 
rm(list = ls())
library(scGATE)
data         <- as.matrix(read.csv(paste0("D:\\scGATE_files\\example_data\\ClusterI.csv"))[ ,2:15])
print(head(data))
             gA       gB         gC        gC1        gC2         gD        gD1        gD2       gE        gE1      gE2         gF        gF1        gF2
[1,] 0.02764677 2.028944 0.01688577 0.01946526 0.02380772 0.01852824 0.02069895 0.02093184 1.932168 0.06889533 1.824497 0.04963150 0.05794413 0.04217521
[2,] 0.02643986 2.027956 0.01730882 0.01963009 0.02459190 0.01834800 0.01909300 0.02050692 1.965628 0.06075294 1.829349 0.04624483 0.04405836 0.02271849
[3,] 0.02593749 2.033592 0.01729984 0.01796269 0.02422372 0.01760521 0.01906345 0.02125864 1.976888 0.04514580 1.817391 0.03175314 0.03738465 0.02141319
[4,] 0.02595885 2.019971 0.01756862 0.01787157 0.02412755 0.01791644 0.02112435 0.02106185 1.980759 0.03720293 1.836962 0.02033092 0.03677651 0.01974638
[5,] 0.02629885 2.015461 0.01753645 0.01921252 0.02491328 0.01909465 0.02101008 0.02132906 1.986872 0.03738554 1.837226 0.01999704 0.03683252 0.01955996
[6,] 0.02640293 2.009388 0.01748322 0.02028304 0.02449585 0.01990073 0.02096272 0.02065111 1.982286 0.03776144 1.834406 0.02226978 0.03658167 0.01937109
```

```R
# 3. data preprocessing 
# For scGATE simulated data, library size normalization is not performed. 
# However, the simulated data is only re-scaled using the quantile normalization technique to fit within the (0,1) interval.
data         <- scRNA_seq_preprocessing(data = data, library_size_normalization = "False")
```

```R
# 4. Remove genes with low variability (scGATE operates on highly variable genes per context).
# This step is optional
data$n_counts<- data$n_counts[ , which(sqrt(apply(data$n_counts,2,var))> 0.20)]
```

```R
# 5. Run scGATE_logic() function
# Please note that the likelihood values can be affected by the Louvain clustering results.
gates        <- scGATE_logic(data = data, number_of_em_iterations = 5, max_num_regulators = 3, abs_cor = 0.05, top_gates = 1, run_mode = "fast")
print(head(gates))
   gene_name -log10 L0 -log10 L1 log10 BF logic_gate
  gene_name -log10 L0 -log10 L1 log10 BF logic_gate
1        gE     173.9   -278.87   452.76        ~gF
2       gE1      59.2    -275.1   334.30    gE.~gE2
3       gE2     41.52   -277.04   318.56    gE.~gE1
4        gF    170.38   -288.46   458.84        ~gE
5       gF1     80.36   -262.54   342.90    gF.~gF2
6       gF2      67.6   -264.85   332.45    gF.~gF1
```


<br>

## Datasets
Raw scRNA-seq data for cell-type specific logic gate inference in the mouse haematopoiesis (Fig.3) is available at https://doi.org/10.5281/zenodo.8353409.
Base grns reconstructed with scATAC-seq data and TF binding site motifs (Fig.6) are also available at https://doi.org/10.5281/zenodo.8353409.



<br>
