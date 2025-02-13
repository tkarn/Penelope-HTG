## Penelope-HTG / Fig4

## Supplementary Resource to:  

Denkert C et al. 2025, Cancer Cell 43, 1–16, February 10, 2025 
https://doi.org/10.1016/j.ccell.2025.01.002

Dynamics of molecular heterogeneity in high-risk luminal breast cancer - From intrinsic to adaptive subtyping.

PMID: 39933898

************************************************************

## This subdirectory contains the code and data to generate Figure 4 from the paper:

### Relationship of gene clusters and survival:
### Scatter plots of 1/HR for biopsy vs resect samples for all 335 genes colored according different classifications of genes 

1. [*HRscat_BiopResec_335genes_Diff-Colors.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/HRscat_BiopResec_335genes_Diff-Colors.pdf):  An R Markdown file of the analysis comparing the prognostic value of individual genes and their association with gene clusters from unsupervised clustering of Penelope-B samples. Data are presented as scatter plots of (inverted) Hazard Ratios for each of the 335 genes measured in preTx- and postTx biopsies, and dots are colored according unsupervised clustering and functional annotations.
2. [*HRscat_BiopResec_335genes_Diff-Colors.R*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/HRscat_BiopResec_335genes_Diff-Colors.R):  The R-script that generates this R-Markdown file.

#### Input files:
[*Inverse_HR_biopsy_resec_335genes.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/Inverse_HR_biopsy_resec_335genes.txt): (Inverted) hazard ratios for DFS for all 335 genes from the Penelope-B clusters, based on pre-treatment and post-treatment samples.

#### Output files:
[*1_HRscat_BiopResec_335genes_col-mainclust_labs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/1_HRscat_BiopResec_335genes_col-mainclust_labs.svg)

[*2_HRscat_BiopResec_335genes_col-mainclust_Nolabs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/2_HRscat_BiopResec_335genes_col-mainclust_Nolabs.svg)

[*3_HRscat_BiopResec_335genes_col-GeneClass_labs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/3_HRscat_BiopResec_335genes_col-GeneClass_labs.svg)

[*4_HRscat_BiopResec_335genes_col-GeneClass_Nolabs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/4_HRscat_BiopResec_335genes_col-GeneClass_Nolabs.svg)

[*5_HRscat_BiopResec_335genes_col-subclust_labs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/5_HRscat_BiopResec_335genes_col-subclust_labs.svg)

[*6_HRscat_BiopResec_335genes_col-subclust_Nolabs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/6_HRscat_BiopResec_335genes_col-subclust_Nolabs.svg)

[*Fig4.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig4/Fig4.pdf): Aligned figure


************************************************************

