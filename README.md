## scGATE: single-cell gene regulatory gate inference
We have developed scGATE (single-cell gene regulatory gate) as a logic-based model for deciphering tissue- and cell-type-specific gene regulatory networks from single-cell RNA sequencing (scRNA-seq). While previous efforts have focused on reconstructing directed transcription factor (TF) to target gene networks, logic-based models enable the exploration of more complex combinatorial relationships between regulators. In particular, Boolean logic models can capture higher-order TF interactions, represented through AND, OR, and XOR logic gates. Our novel approach infers TF-gene networks from scRNA-seq data while simultaneously elucidating the underlying Boolean logic that combines TF activities.
The methodology integrates external genomic data, such as single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq) and motif analysis, to narrow down candidate TFs for each target gene. TF-gene links are further refined using scRNA-seq data, and logic rules are derived. To enhance statistical power, scGATE focuses on reconstructing context specific regulatory networks, including tissue- or cell-type-specific networks. This study presents an integrative framework for deducing cell-type-specific gene regulation, moving beyond TF-gene pairs to capture the complex logic operations underlying combinatorial control.

<br>
Malekpour, S.A., Haghverdi, L., Sadeghi, M., Single-cell multi-omics analysis identifies context specific gene regulatory gates and mechanisms. Briefings In Bioinformatics, 2023.

##### If you find our study useful and relevant to your research, please kindly cite us. Your citation means a lot to us and helps acknowledge our contributions.
<br>
<br>

![Fig1](https://github.com/CompBioIPM/scGATE/assets/47293318/0d596d39-ec44-4e23-93d3-e850ed8727d1)


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
install.packages("doParallel")
install.packages("foreach")
install.packages("doSNOW")
```
In order to run scGATE with parallel computing, the packages doParallel, foreach, and doSNOW need to be installed.

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
candidate_tf_target <- as.data.frame(read_parquet("Buenrostro2018_base_GRN_dataframe.parquet"))
candidate_tf_target <- read_base_GRN(candidate_tf_target)
```

### Preprocessing scRNA-seq count data
To preprocess raw scRNA-seq data, including steps such as normalization and rescaling, you can use the scRNA_seq_preprocessing() function from the scGATE package.

```R
# Preprocess scRNA-seq count data
normalized_counts <- scRNA_seq_preprocessing(data = data_scRNA_seq, library_size_normalization = "True", tf_list = NA)
```

Parameter Descriptions  
```bash
# data                       The scRNA-seq raw data matrix with cells in rows and genes in columns.

# library_size_normalization A flag indicating whether library size normalization should be performed.
#                            The default value is "True".
#                            Set it to "False" if you don't want to perform library size normalization.
  
# tf_list                    A list of transcription factors (TFs) to consider.
#                            The default value is NA, which means all columns in the data matrix will be considered as TFs.  
```

<br>

## Step 3. Run scGATE 
scGATE provides two functions for TF-target network inference: `scGATE_gate()` and `scGATE_edge()`. These functions infer the TF-target network with and without predicted Boolean logic gates in the output, respectively.
The scGATE_gate() function in the scGATE package is more suitable for small networks or when the base gene regulatory network (GRN) is available from external sources such as scATAC-seq and TF motif data.  


### TF-Target Network Inference (gate mode)
To infer the TF-target network with logic gates in the output, you can use the `scGATE_gate()` function.

```R
# Infer TF-target network without logic gates in the output
gates <- scGATE_logic(data = data, base_GRN = NA, h_set = NA, number_of_em_iterations = NA, max_num_regulators = NA, abs_cor = NA, top_gates = NA, run_mode = NA)
print(head(gates))
```
Parameter Descriptions   
```bash
# data                    A gene expression matrix with normalized counts within the (0,1) interval,
#                         where samples are represented as rows and genes as columns.
#                         The gene expression matrix should have been preprocessed using the scRNA_seq_preprocessing() function.

# base_GRN                Base TF-gene interaction network derived from external hints
#                         (e.g., scATAC-seq data and TF binding site motifs on DNA).
#                         Leave it empty if no base GRN is available.

# h_set                   The range of possible values for the "h" parameter in the Hill climbing function.
 
# number_of_em_iterations The number of iterations in the expectation-maximization (EM) algorithm.
#                         The default value is 3.

# max_num_regulators      Maximum number of TFs in a logic gate that can regulate the target gene profile.
#                         The default value is 3.

# abs_cor                 This parameter varies in the (0, 1) interval and further removes edges with low absolute Pearson correlations between TFs and their targets.
#                         A (default) value of 0 indicates no filtration based on correlations.
  
# top_gates               The number of top Boolean logic gates to be reported for each target gene, based on Bayes Factor.
#                         The default value is 1.
  
# run_mode                Use "simple" for a faster algorithm run and "complex" for more precise results that take more time.
#                         The argument is relevant to the possible complexities in the hill function parameter space for regulatory TFs and target genes.
#                         The default value is "simple".

# weight_threshold        The output form scGATE will present the logic combination or partition that yields a certain percentage of the target gene,
#                         specifically when it is above the weight_threshold.
#                         The default value is 0.05.

# num_cores               Specify the number of parallel workers (adjust according to your system).
```

<br>

### TF-Target Network Inference (edge mode)
To infer the TF-target network without logic gates in the output, you can use the `scGATE_edge()` function.

```R
# Infer TF-target network without logic gates in the output
edges <- scGATE_edge(data = data, base_GRN = candidate_tf_target, h_act = NA, number_of_em_iterations = NA, max_num_regulators = NA, abs_cor = NA)
print(head(edges))
```

Parameter Descriptions 
```bash
# data                    A gene expression matrix with normalized counts within the (0,1) interval,
#                         where samples are represented as rows and genes as columns.
#                         The gene expression matrix should have been preprocessed using the scRNA_seq_preprocessing() function.
  
# base_GRN                Base TF-gene interactions derived from external hints
#                         (e.g., scATAC-seq data and TF binding site motifs on DNA).
#                         Leave it empty if no base GRN is available.

# h_act                   Parameter of the Hill climbing function.
#                         It is the hill coefficient that represents the cooperativity or sigmoidicity of the TF regulatory response.
#                         The default value is 7.
  
# number_of_em_iterations The number of iterations in the expectation-maximization (EM) algorithm.
#                         The default value is 3.
 
# max_num_regulators      Maximum number of TFs in a logic gate that can regulate the target gene profile.
#                         The default value is 3.

# abs_cor                 This parameter varies in the (0, 1) interval and further removes edges with low absolute Pearson correlations between TFs and their targets.
#                         A (default) value of 0 indicates no filtration based on correlations.

# num_cores               Specify the number of parallel workers (adjust according to your system) 
```
<br>

## Example usage of scGATE 
### I. Context specific network and logic gate inference in synthetic toggle switch 

```R
# 1. Please refer to the Jupyter notebook for instructions on how to perform Louvain clustering on the cells in the BoolODE simulated data.
# 2. Retrieve the data from Cluster I of cells, which was obtained in the previous step.
# Load scGATE package and data in example_data folder
 
rm(list = ls())
library(scGATE)
data <- as.matrix(read.csv("/example_data/ClusterI.csv")[ ,2:15])
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
# However, the simulated data is only re-scaled using the quantile normalization technique to fit the data within the (0,1) interval.
data <- scRNA_seq_preprocessing(data = data, library_size_normalization = "False")
```

```R
# 4. Remove genes with low variability (scGATE operates on highly variable genes per context).
# This step is optional
data$n_counts <- data$n_counts[ , which(sqrt(apply(data$n_counts,2,var))> 0.20)]
```

```R
# 5. Run scGATE_logic() function
# Please note that the likelihood values can be affected by the Louvain clustering results.
gates <- scGATE_logic(data = data, top_gates = 1, run_mode = "simple")
print(head(gates))
  gene_name -log10 L0 -log10 L1 log10 BF logic_gate
1        gE     173.9   -268.57   442.47        ~gF
2       gE1     51.85   -234.65   286.50    gE.~gE2
3       gE2     38.43   -235.48   273.91    gE.~gE1
4        gF    170.38   -278.57   448.95        ~gE
5       gF1     80.36   -215.32   295.68    gF.~gF2
6       gF2      67.6   -217.88   285.48    gF.~gF1
```

<br>

### II. Context specific network and logic gate inference in the mouse haematopoiesis scRNA-seq data 

```R
# 1. Please refer to the Jupyter notebook for instructions on how to perform Louvain clustering on the cells in the mouse haematopoiesis scRNA-seq dataset.
# 2. Retrieve the data from Megakaryocyte cells (Cluster 11).
# Load scGATE package and data in example_data folder

rm(list = ls())
library(scGATE)
data <- as.data.frame(read.csv("/example_data/subset_counts_cluster_11.csv" , header = TRUE))

# select genes involved in the MegE differentiation
gene_list <- c("Gata1", "Fli1", "Klf1", "Spi1", "Zfpm1", "Tal1", "Gata2")
data <- data[  , gene_list]
data <- na.omit(data)
print(head(data))
      Gata1      Fli1     Klf1      Spi1     Zfpm1      Tal1     Gata2
1 0.6931472 1.0986123 0.000000 0.6931472 0.0000000 0.6931472 0.0000000
2 0.0000000 1.3862944 0.000000 0.0000000 0.0000000 0.6931472 1.0986123
3 0.6931472 1.6094380 0.000000 0.0000000 0.0000000 0.0000000 0.6931472
4 0.0000000 0.0000000 1.098612 0.0000000 0.6931472 0.0000000 1.6094380
5 0.0000000 0.0000000 0.000000 0.0000000 0.6931472 0.6931472 1.3862944
6 0.0000000 0.6931472 0.000000 0.0000000 0.6931472 1.0986123 0.0000000

# Load base GRN
base_GRN <- read.csv("/example_data/base_grn_mouse_blood_cell_differentiation_toggle_switch.csv")
```

```R
# 3. data preprocessing 
# The dataset underwent library size normalization in Jupyter Notebook. To fit the scRNA-seq data within the (0,1) interval, we applied quantile normalization as a technique to rescale the data.
data <- scRNA_seq_preprocessing(data = data, library_size_normalization = "False")
```

```R
# 4. Run scGATE_logic() function
gates <- scGATE_logic(data = data, base_GRN = base_GRN, number_of_em_iterations = 10, top_gates = 1, run_mode = "complex")
print(head(gates))


# To effectively derive Boolean rules from extensive scRNA-seq datasets containing over 10 TFs, we recommend employing scGATE with the following configuration.
gates <- scGATE_logic(data = data, base_GRN = base_GRN, number_of_em_iterations = 10, max_num_regulators = 2, top_gates = 50, run_mode = "simple")
print(gates)

or you may use,
h_set <- c(1.25, 2.25)
gates <- scGATE_logic(data = data, base_GRN = base_GRN, h_set = h_set, number_of_em_iterations = 10, max_num_regulators = 2, top_gates = 50, run_mode = "complex")
print(gates)
```

<br>

### III. Context specific network inference in mouse tissue scRNA-seq datasets

```R
# 1. Please refer to the Jupyter notebook for instructions on how to perform scATAC-seq analysis to derive the candidate TF lists (base GRNs) in *.parquet file format.
# 2. Load scGATE package and data (base GRN and scRNA-seq data and TF list) in example_data folder 

rm(list=ls())
library(scGATE)
# Load base GRN derived from external hints
candidate_tf_target <- as.data.frame(read_parquet("/example_data/Cusanovich2018_Spleen_peak_base_GRN_dataframe.parquet"))
candidate_tf_target <- read_base_GRN(candidate_tf_target)

# Load scRNA-seq data
data           <- as.data.frame(read.csv("/example_data/Tabula_Muris2018_Spleen-10X_P4_7_ExpressionData.csv" , header = TRUE))
gene_names     <- data[ ,1]
data           <- t(data[ ,2:ncol(data)])
colnames(data) <- gene_names

head(data[ , 1:10])
                   Batf Stat5b Ctcf H2-Eb1 AW112010 Ly6d Rplp0 Id2 Dok2 Gimap3
AAACCTGAGAAGGACA.1    0      0    0     18        0    0    10   0    0      0
AAACCTGAGCTAAGAT.1    0      0    1      0       19    0     5   1    1      1
AAACCTGCAACAACCT.1    0      0    0     22        0    5    12   0    0      2
AAACCTGCAGCCAATT.1    0      0    0     14        1    5    21   0    0      1
AAACCTGCAGCTCCGA.1    0      0    1     30        1    2    64   0    0      0
AAACCTGTCAGGTAAA.1    0      0    0     23        3    8    24   0    0      0


# Load TF list
# This step is optional
tf_names       <- unlist(read.table("/example_data/Tabula_Muris2018_Spleen-10X_P4_7_tf_lists.txt"))
print(head(tf_names))
      V1       V2       V3 
  "Batf" "Stat5b"   "Ctcf"
```

```R
# 3. scRNA-seq data preprocessing (library size normalization, quantile normalization technique to fit the scRNA-seq data within the (0,1) interval) 
data           <- scRNA_seq_preprocessing(data = data, library_size_normalization = "True", tf_list = tf_names)
```

```R
# 4. Run scGATE_edge() function
ranked_edge_list <-  scGATE_edge(data = data, base_GRN = candidate_tf_target, h_act = 7)
print(head(ranked_edge_list))
    from    to BF_score
1   Ctcf Rps19 2013.587
2   Batf Rps19 2012.551
3 Stat5b Rplp0 1850.334
4   Ctcf Rplp0 1849.896
5   Ctcf Rpl36 1649.263
6   Ctcf Eif5a 1559.044
```

<br>

### IV. Context specific network inference in human haematopoiesis scRNA-seq dataset

```R
# 1. Please refer to the Jupyter notebook for instructions on how to perform scATAC-seq analysis to derive the candidate TF lists (base GRNs) in *.parquet file format.
# 2. Load scGATE package and data (base GRN and scRNA-seq data and TF list) in example_data folder 

rm(list=ls())
library(scGATE)
# Load base GRN derived from external hints
candidate_tf_target <- as.data.frame(read_parquet("/example_data/Buenrostro2018_base_GRN_dataframe.parquet"))
candidate_tf_target <- read_base_GRN(candidate_tf_target)

# Load scRNA-seq data
data           <- as.data.frame(read.csv("/example_data/Buenrostro2018_ExpressionData.csv" , header = TRUE))
gene_names     <- data[ ,1]
data           <- t(data[ ,2:ncol(data)])
colnames(data) <- gene_names

head(data[ , 1:10])
      IRF8 FOS MAFF SPI1 JUNB SPIB IRF7 TFDP1 GATA1 RAD21
hsc_1    0   2    0    0    2    0    0     0     0     1
hsc_2    0   6    7    0    3    0    0     0     0     1
hsc_3    0   2    0    0    5    0    0     0     0     2
hsc_4    0   6    0    0    1    0    0     1     0     1
hsc_5    0   1    5    2    1    0    0     0     0     0
hsc_6    0   3    0    0    1    0    0     0     0     0

# Load TF list
# This step is optional
tf_names       <- unlist(read.table("/example_data/Buenrostro2018_tf_lists.txt"))
print(head(tf_names))
    V1     V2     V3     V4     V5     V6 
"IRF8"  "FOS" "MAFF" "SPI1" "JUNB" "SPIB" 
```

```R
# 3. scRNA-seq data preprocessing (library size normalization, quantile normalization technique to fit the scRNA-seq data within the (0,1) interval)
data           <- scRNA_seq_preprocessing(data = data, library_size_normalization = "True", tf_list = tf_names)
```

```R
# 4. Run scGATE_edge() function
ranked_edge_list <-  scGATE_edge(data = data, base_GRN = candidate_tf_target, h_act = 7)
print(head(ranked_edge_list))
     from     to BF_score
1    E2F1 MALAT1 13415.34
2 BHLHE40 MALAT1 13415.34
3   TFDP1 MALAT1 13415.32
4    NFE2 MALAT1 13414.98
5    IRF8 MALAT1 13414.26
6   RUNX2   PTMA 11592.68
```

<br>

## Datasets
Raw scRNA-seq data for cell-type specific logic gate inference in the mouse haematopoiesis dataset (Fig.3) is also available at Zenodo https://doi.org/10.5281/zenodo.8353409.  
The base GRNs reconstructed with scATAC-seq and TF binding site motif, in mouse tissue and human haematopoiesis datasets, together with other intermediate and processed files are available at Zenodo https://doi.org/10.5281/zenodo.8353409.



<br>
