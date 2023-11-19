## scGATE: single-cell gene regulatory gate inference
We have developed scGATE (single-cell gene regulatory gate) as a logic-based model for deciphering tissue- and cell-type-specific gene regulatory networks from single-cell RNA sequencing (scRNA-seq). While previous efforts have focused on reconstructing directed transcription factor (TF) to target gene networks, logic-based models enable the exploration of more complex combinatorial relationships between regulators. In particular, Boolean logic models can capture higher-order TF interactions, represented through AND, OR, and XOR logic gates. Our novel approach infers TF-gene networks from scRNA-seq data while simultaneously elucidating the underlying Boolean logic that combines TF activities.
The methodology integrates external genomic data, such as single-cell assay for transposase-accessible chromatin sequencing (scATAC-seq) and motif analysis, to narrow down candidate TFs for each target gene. TF-gene links are further refined using scRNA-seq data, and logic rules are derived. To enhance statistical power, scGATE focuses on reconstructing context-specific regulatory networks, including tissue- or cell type-specific networks. This study presents an integrative framework for deducing cell-type-specific gene regulation, moving beyond TF-gene pairs to capture the complex logic operations underlying combinatorial control.

<br>
Malekpour, S.A., Haghverdi, L., Sadeghi, M., Single-cell multi-omics analysis identifies context specific gene regulatory gates and mechanisms. Briefings In Bioinformatics, 2023.

##### If you find our study useful and relevant to your research, please kindly cite us. Your citation means a lot to us and helps acknowledge our contributions.
<br>
<br>

![Fig1](https://github.com/CompBioIPM/scGATE/assets/47293318/4be29239-0cd7-4871-aa01-18268ba7bb7d)

## Step 1. scGATE installation



<br>

## Step 2. Prepare input files



<br>

## Step 3. Run scGATE 



<br>

## Datasets
Raw scRNA-seq data for cell-type specific logic gate inference in the mouse haematopoiesis (Fig.3) is available at https://doi.org/10.5281/zenodo.8353409.
Base grns reconstructed with scATAC-seq data and TF binding site motifs (Fig.6) are also available at https://doi.org/10.5281/zenodo.8353409.



<br>
