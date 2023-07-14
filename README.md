# viral_coinfections

The human respiratory tract hosts a community of viruses that cause epidemics and pandemics. 
Respiratory Syncytial Virus (RSV) and Influenza are two significant contributors to respiratory infections. 
We sought to characterize the mechanism of RSV-IAV interactions by analyzing bulk transcriptomics data from experimental conditions of air-liquid interface cell cultures representing the human respiratory tract.
 

![image](https://github.com/marthaluka/viral_coinfections/figures/workflow.png)

Figure 1. Wet-lab and dry-lab workflows. 
(A) The laboratory workflow: A549 cells were infected with RSV, IAV or both viruses and incubated for 17 hours before fixation, and RNA sequencing using the Illumina NextSeq 500. 
(B) Bioinformatic analysis which consisted of quality control, mapping of reads and feature counts. The count table was subsequently imported into R and analysed using a suite of packages.


![image](https://github.com/marthaluka/viral_coinfections/figures/figure2.png)

Figure 2. Transcriptional perturbations by single and double infections in lung alveolar cell cultures. 
(A) Principal component analysis of replicate transcriptomes. Each point represents a sample and points are colored by infection status. 
(B) Number of differentially expressed genes following virus infection at 17 hours post infection. 
(C) Mean expression levels for each gene per infection status, with ten most significant genes labelled. Upregulated genes are colored orange, downregulated genes blue and insignificant genes grey. 

![image](https://github.com/marthaluka/viral_coinfections/figures/figure3.png)

Figure 3. Expression levels of common/shared genes. 
(A) Pairwise log2foldchange (L2FC) values comparisons across the different infection statuses. Each point represents a gene, and lighter colour indicates a higher L2FC value while darker colour indicates lower L2FC values. Genes dysregulated in an opposite direction are labelled. 
(B) A heatmap generated from transcriptomic data using commonly dysregulated genes across virus infection conditions. The intensity of the color in each cell represents the expression level of the gene in the corresponding sample, with red indicating higher relative expression and blue lower relative expression.


![image](https://github.com/marthaluka/viral_coinfections/figures/figure4.png)

Figure 4. Top 30 enriched GO biological processes enriched in (A) RSV, (B) IAV, or (C) Coinfected status when compared with mock infection at 17 hpi. (D) Top 30 enriched GO biological processes enriched of the 200 shared/ common DEGs across the infections.  


![image](https://github.com/marthaluka/viral_coinfections/figures/figure5.png)
Figure 5. Top 5 Gene concept networks of differentially expressed genes. In (A)RSV, (B) IAV and (C) coinfected cells. Red: up-regulated genes, Blue: down-regulated genes ; Size: number of genes per node. 